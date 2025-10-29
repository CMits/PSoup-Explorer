# pages/02_psoup_gxe.py

import streamlit as st
import pandas as pd

from psoup_gxe import PSoupSimulatorGxE
from analysis import (
    plot_sustained_growth_curve_interactive,
    plot_delta_vs_WT_curve_interactive,
)
from scenarios import (
    MODIFIER_FIELDS,
    parse_modifiers_from_text_inputs,
    run_and_package,
    bulk_generate,
)
from plots import (
    make_multi_sustained_growth_plot,
    make_multi_delta_plot,
    make_scaling_fig,
)

st.title("🍬 PSoup G×E Explorer (sucrose × SL interaction)")
st.markdown("""
This page uses **PSoupSimulatorGxE**, which encodes sucrose relieving a
strigolactone (SL) brake **only if** the SL pathway is intact.

Biology we're capturing:
- Wild type (mSL_pathway_intact ≈ 1):
  - SL brake is active at low sucrose.
  - High sucrose turns that brake off → big jump in branching.
- max2 / smxl6_7_8-like mutant (mSL_pathway_intact ≈ 0):
  - SL brake basically gone even at low sucrose.
  - So branching starts high already, and adding more sucrose gives a smaller boost.
  - Slope flattens ⇒ true genotype × sucrose interaction.

This is what field people and branching papers actually report.
""")

# ---------------------------------------------------------
# Session init
# ---------------------------------------------------------
if "results_accum_gxe" not in st.session_state:
    st.session_state["results_accum_gxe"] = []

# ---------------------------------------------------------
# Sidebar: manual scenario
# ---------------------------------------------------------
st.sidebar.header("New Scenario (G×E PSoup)")
st.sidebar.markdown("""
Blank = 1.0 (WT).
To mimic max2 / smxl6_7_8, set `mSL_pathway_intact` close to 0.
""")

scenario_name = st.sidebar.text_input("Scenario name", value="my_scenario_gxe")
raw_modifiers = {}
for field in MODIFIER_FIELDS:
    raw_modifiers[field] = st.sidebar.text_input(field, value="", key=f"gxe_{field}")

if st.sidebar.button("Run scenario (G×E)"):
    sim = PSoupSimulatorGxE()
    mods_parsed = parse_modifiers_from_text_inputs(raw_modifiers)
    row = run_and_package(sim, scenario_name, mods_parsed)
    st.session_state["results_accum_gxe"].append(row)
    st.sidebar.success(f"Scenario '{scenario_name}' added (G×E).")

# ---------------------------------------------------------
# Sidebar: bulk scenario generation
# ---------------------------------------------------------
st.sidebar.markdown("---")
st.sidebar.header("Bulk scenario generation (G×E PSoup)")

bulk_mode = st.sidebar.selectbox(
    "Mode",
    [
        "all_modifiers_change",
        "sparse_modifiers_change",
        "mutant_single_knockout",
        "mutant_double_knockout",
    ],
    help="Knocking out `mSL_pathway_intact` in this model should flatten sucrose response."
)

num_random = st.sidebar.number_input(
    "How many scenarios?",
    min_value=1, max_value=1000, value=50, step=10,
    key="gxe_num_random"
)
min_val = st.sidebar.number_input(
    "Min modifier value (random modes)",
    min_value=0.0, max_value=10.0, value=0.0, step=0.1,
    key="gxe_min_val"
)
max_val = st.sidebar.number_input(
    "Max modifier value (random modes)",
    min_value=0.0, max_value=10.0, value=2.0, step=0.1,
    key="gxe_max_val"
)
p_change = st.sidebar.slider(
    "Probability modifier changes (sparse mode)",
    min_value=0.0, max_value=1.0, value=0.3, step=0.05,
    key="gxe_p_change"
)

if st.sidebar.button("Generate bulk (G×E)"):
    new_rows = bulk_generate(
        sim_class=PSoupSimulatorGxE,
        mode=bulk_mode,
        n=num_random,
        min_val=min_val,
        max_val=max_val,
        p_change=p_change,
        existing_count=len(st.session_state["results_accum_gxe"])
    )
    st.session_state["results_accum_gxe"].extend(new_rows)
    st.sidebar.success(f"Generated {num_random} '{bulk_mode}' scenarios (G×E).")

# ---------------------------------------------------------
# Sidebar: reset
# ---------------------------------------------------------
st.sidebar.markdown("---")
if st.sidebar.button("🔁 Reset G×E scenarios"):
    st.session_state["results_accum_gxe"] = []
    st.sidebar.warning("All G×E scenarios cleared.")
    st.stop()

# ---------------------------------------------------------
# Scenario table
# ---------------------------------------------------------
st.subheader("All simulated scenarios (G×E PSoup)")
if len(st.session_state["results_accum_gxe"]) == 0:
    st.info("No scenarios yet. Add one in the sidebar.")
    st.stop()

df_all = pd.DataFrame(st.session_state["results_accum_gxe"])
st.dataframe(df_all, use_container_width=True)

# ---------------------------------------------------------
# Per-scenario / overlay
# ---------------------------------------------------------
st.subheader("Single-scenario detail (G×E PSoup)")

show_all = st.checkbox(
    "Show ALL G×E scenarios overlaid",
    value=False,
    key="gxe_show_all"
)

scenario_options = [row["genotype"] for row in st.session_state["results_accum_gxe"]]
selected_name = st.selectbox(
    "Inspect scenario",
    options=scenario_options,
    key="gxe_selected"
)

sel_row = next(
    row for row in st.session_state["results_accum_gxe"]
    if row["genotype"] == selected_name
)

col_a, col_b = st.columns(2)

with col_a:
    st.markdown("**Sustained growth vs sucrose (G×E)**")
    if show_all:
        fig_all = make_multi_sustained_growth_plot(st.session_state["results_accum_gxe"])
        st.plotly_chart(fig_all, use_container_width=True)
    else:
        fig_one = plot_sustained_growth_curve_interactive(sel_row)
        st.plotly_chart(fig_one, use_container_width=True)

with col_b:
    st.markdown("**Mutant − WT across sucrose (G×E)**")
    if show_all:
        fig_all2 = make_multi_delta_plot(st.session_state["results_accum_gxe"])
        st.plotly_chart(fig_all2, use_container_width=True)
    else:
        fig_one2 = plot_delta_vs_WT_curve_interactive(sel_row)
        st.plotly_chart(fig_one2, use_container_width=True)

# ---------------------------------------------------------
# Global scaling / linearity
# ---------------------------------------------------------
st.subheader("Global scaling across all G×E scenarios")

df_scaling = pd.DataFrame([
    {
        "genotype": r["genotype"],
        "sg_low":  r["sg_low"],
        "sg_base": r["sg_base"],
        "sg_high": r["sg_high"]
    }
    for r in st.session_state["results_accum_gxe"]
])

col_c, col_d, col_e = st.columns(3)

with col_c:
    fig_lb, footer_lb = make_scaling_fig(
        df_scaling,
        xcol="sg_low",
        ycol="sg_base",
        xtitle="Sustained_growth @ low sucrose (0.5×)",
        ytitle="Sustained_growth @ base sucrose (1.0×)",
        title="low → base scaling (G×E)",
        fit_label_template="fit: base={slope:.3f}·low+{intercept:.3f} (R²={r2:.4f})"
    )
    st.plotly_chart(fig_lb, use_container_width=True)
    st.markdown(footer_lb)

with col_d:
    fig_bh, footer_bh = make_scaling_fig(
        df_scaling,
        xcol="sg_base",
        ycol="sg_high",
        xtitle="Sustained_growth @ base sucrose (1.0×)",
        ytitle="Sustained_growth @ high sucrose (1.5×)",
        title="base → high scaling (G×E)",
        fit_label_template="fit: high={slope:.3f}·base+{intercept:.3f} (R²={r2:.4f})"
    )
    st.plotly_chart(fig_bh, use_container_width=True)
    st.markdown(footer_bh)

with col_e:
    fig_lh, footer_lh = make_scaling_fig(
        df_scaling,
        xcol="sg_low",
        ycol="sg_high",
        xtitle="Sustained_growth @ low sucrose (0.5×)",
        ytitle="Sustained_growth @ high sucrose (1.5×)",
        title="low → high scaling (G×E)",
        fit_label_template="fit: high={slope:.3f}·low+{intercept:.3f} (R²={r2:.4f})"
    )
    st.plotly_chart(fig_lh, use_container_width=True)
    st.markdown(footer_lh)

st.markdown("""
If SL-intact mutants (`mSL_pathway_intact ≈ 0`) start out high at low sucrose
and then *do not* increase as steeply with sucrose, they will peel away from
the single global line and R² will drop.

That break in linearity is exactly the sucrose-by-genotype (G×E) behaviour
branching biologists report in max2 / smxl6_7_8 mutants.
""")
