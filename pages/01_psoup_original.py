# pages/01_psoup_original.py

import streamlit as st
import pandas as pd

from psoup import PSoupSimulator
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

st.title("🌱 PSoup Explorer (original model)")
st.markdown("""
This page uses the **original PSoup equations** (no explicit sucrose × SL interaction).

What you'll usually see:
- All genotypes scale almost perfectly with sucrose.
- The low→base, base→high, and low→high plots all have R² ≈ 1.
- Even knockouts and double knockouts sit on the same line.

Interpretation:
> Sucrose is acting like a single global gain, not an environment that
> different genotypes interpret differently. That means **no true G×E**.
""")

# ---------------------------------------------------------
# Session init for this page
# ---------------------------------------------------------
if "results_accum_psoup" not in st.session_state:
    st.session_state["results_accum_psoup"] = []

# ---------------------------------------------------------
# Sidebar: Manual scenario
# ---------------------------------------------------------
st.sidebar.header("New Scenario (original PSoup)")
st.sidebar.markdown("Blank = 1.0 (WT).")

scenario_name = st.sidebar.text_input("Scenario name", value="my_scenario_orig")
raw_modifiers = {}
for field in MODIFIER_FIELDS:
    raw_modifiers[field] = st.sidebar.text_input(field, value="", key=f"orig_{field}")

if st.sidebar.button("Run scenario (original)"):
    sim = PSoupSimulator()
    mods_parsed = parse_modifiers_from_text_inputs(raw_modifiers)
    row = run_and_package(sim, scenario_name, mods_parsed)
    st.session_state["results_accum_psoup"].append(row)
    st.sidebar.success(f"Scenario '{scenario_name}' added (original).")

# ---------------------------------------------------------
# Sidebar: Bulk scenarios
# ---------------------------------------------------------
st.sidebar.markdown("---")
st.sidebar.header("Bulk scenario generation (original PSoup)")

bulk_mode = st.sidebar.selectbox(
    "Mode",
    [
        "all_modifiers_change",
        "sparse_modifiers_change",
        "mutant_single_knockout",
        "mutant_double_knockout",
    ],
    help="mutant_* knocks out (sets to 0) 1 or 2 modifiers."
)

num_random = st.sidebar.number_input(
    "How many scenarios?",
    min_value=1, max_value=1000, value=50, step=10,
    key="orig_num_random"
)
min_val = st.sidebar.number_input(
    "Min modifier value (random modes)",
    min_value=0.0, max_value=10.0, value=0.0, step=0.1,
    key="orig_min_val"
)
max_val = st.sidebar.number_input(
    "Max modifier value (random modes)",
    min_value=0.0, max_value=10.0, value=2.0, step=0.1,
    key="orig_max_val"
)
p_change = st.sidebar.slider(
    "Probability modifier changes (sparse mode)",
    min_value=0.0, max_value=1.0, value=0.3, step=0.05,
    key="orig_p_change"
)

if st.sidebar.button("Generate bulk (original)"):
    new_rows = bulk_generate(
        sim_class=PSoupSimulator,
        mode=bulk_mode,
        n=num_random,
        min_val=min_val,
        max_val=max_val,
        p_change=p_change,
        existing_count=len(st.session_state["results_accum_psoup"])
    )
    st.session_state["results_accum_psoup"].extend(new_rows)
    st.sidebar.success(f"Generated {num_random} '{bulk_mode}' scenarios (original).")

# ---------------------------------------------------------
# Sidebar: reset
# ---------------------------------------------------------
st.sidebar.markdown("---")
if st.sidebar.button("🔁 Reset original scenarios"):
    st.session_state["results_accum_psoup"] = []
    st.sidebar.warning("All original-model scenarios cleared.")
    st.stop()

# ---------------------------------------------------------
# Table of scenarios
# ---------------------------------------------------------
st.subheader("All simulated scenarios (original PSoup)")
if len(st.session_state["results_accum_psoup"]) == 0:
    st.info("No scenarios yet. Add one in the sidebar.")
    st.stop()

df_all = pd.DataFrame(st.session_state["results_accum_psoup"])
st.dataframe(df_all, use_container_width=True)

# ---------------------------------------------------------
# Per-scenario / overlay view
# ---------------------------------------------------------
st.subheader("Single-scenario detail (original PSoup)")

show_all = st.checkbox(
    "Show ALL original-model scenarios overlaid",
    value=False,
    key="orig_show_all"
)

scenario_options = [row["genotype"] for row in st.session_state["results_accum_psoup"]]
selected_name = st.selectbox(
    "Inspect scenario",
    options=scenario_options,
    key="orig_selected"
)

sel_row = next(
    row for row in st.session_state["results_accum_psoup"]
    if row["genotype"] == selected_name
)

col_a, col_b = st.columns(2)

with col_a:
    st.markdown("**Sustained growth vs sucrose**")
    if show_all:
        fig_all = make_multi_sustained_growth_plot(st.session_state["results_accum_psoup"])
        st.plotly_chart(fig_all, use_container_width=True)
    else:
        fig_one = plot_sustained_growth_curve_interactive(sel_row)
        st.plotly_chart(fig_one, use_container_width=True)

with col_b:
    st.markdown("**Mutant − WT across sucrose**")
    if show_all:
        fig_all2 = make_multi_delta_plot(st.session_state["results_accum_psoup"])
        st.plotly_chart(fig_all2, use_container_width=True)
    else:
        fig_one2 = plot_delta_vs_WT_curve_interactive(sel_row)
        st.plotly_chart(fig_one2, use_container_width=True)

# ---------------------------------------------------------
# Global scaling / linearity
# ---------------------------------------------------------
st.subheader("Global scaling across all scenarios (original PSoup)")

df_scaling = pd.DataFrame([
    {
        "genotype": r["genotype"],
        "sg_low":  r["sg_low"],
        "sg_base": r["sg_base"],
        "sg_high": r["sg_high"]
    }
    for r in st.session_state["results_accum_psoup"]
])

col_c, col_d, col_e = st.columns(3)

with col_c:
    fig_lb, footer_lb = make_scaling_fig(
        df_scaling,
        xcol="sg_low",
        ycol="sg_base",
        xtitle="Sustained_growth @ low sucrose (0.5×)",
        ytitle="Sustained_growth @ base sucrose (1.0×)",
        title="low → base scaling (original)",
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
        title="base → high scaling (original)",
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
        title="low → high scaling (original)",
        fit_label_template="fit: high={slope:.3f}·low+{intercept:.3f} (R²={r2:.4f})"
    )
    st.plotly_chart(fig_lh, use_container_width=True)
    st.markdown(footer_lh)

st.markdown("""
If R² stays ~1 and all points lie on the same line in all three panels,
then sucrose is acting like a uniform gain and **the model cannot express G×E**.
""")
