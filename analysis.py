# analysis.py
import numpy as np
import plotly.graph_objects as go
from sklearn.linear_model import LinearRegression
import pandas as pd

# sucrose scaling levels we test
SUCROSE_LEVELS = {
    "low": 0.5,
    "base": 1.0,
    "high": 1.5,
}

def run_scenario(sim, modifiers, sucrose_scale):
    """
    sim: PSoupSimulator or PSoupSimulatorGxE (already constructed)
    modifiers: dict of gene modifier multipliers (e.g. mMAX2_n = 0)
    sucrose_scale: float (0.5, 1.0, 1.5)
    returns sustained growth steady-state
    """
    # start from default WT gene values for THIS simulator
    gene_vals = sim.default_gene_vals.copy()

    # apply user's genotype tweaks
    if modifiers:
        gene_vals.update(modifiers)

    # apply the sucrose level
    gene_vals["mSUC_n"] = gene_vals["mSUC_n"] * sucrose_scale

    # run steady-state and get sustained growth
    sg_val = sim.sustained_growth(gene_vals=gene_vals)
    return sg_val

def simulate_genotype(sim, scenario_name, modifiers):
    """
    Run one genotype across low/base/high sucrose.
    Compute deltas vs WT at each sucrose level.
    Return a dict of results.
    """
    sg_low  = run_scenario(sim, modifiers, SUCROSE_LEVELS["low"])
    sg_base = run_scenario(sim, modifiers, SUCROSE_LEVELS["base"])
    sg_high = run_scenario(sim, modifiers, SUCROSE_LEVELS["high"])

    # WT baselines at each sucrose, using SAME simulator class but WT modifiers (all 1.0)
    wt_mods = {}  # empty means all default 1.0
    sg_low_wt  = run_scenario(sim, wt_mods, SUCROSE_LEVELS["low"])
    sg_base_wt = run_scenario(sim, wt_mods, SUCROSE_LEVELS["base"])
    sg_high_wt = run_scenario(sim, wt_mods, SUCROSE_LEVELS["high"])

    # deltas mutant - WT
    delta_low  = sg_low  - sg_low_wt
    delta_base = sg_base - sg_base_wt
    delta_high = sg_high - sg_high_wt

    return {
        "genotype": scenario_name,
        "sg_low": sg_low,
        "sg_base": sg_base,
        "sg_high": sg_high,
        "delta_low_vs_WT":  delta_low,
        "delta_base_vs_WT": delta_base,
        "delta_high_vs_WT": delta_high,
    }

# ---------- Plot helpers (same idea we discussed before, with epsilon clamp) ----------

EPS = 1e-3

def _clean_delta(val):
    if val is None:
        return val
    try:
        if abs(float(val)) < EPS:
            return 0.0
    except Exception:
        pass
    return val

def plot_sustained_growth_curve_interactive(row):
    """
    row contains sg_low, sg_base, sg_high for one genotype.
    """
    suc_x = [SUCROSE_LEVELS["low"], SUCROSE_LEVELS["base"], SUCROSE_LEVELS["high"]]
    sg_y = [row["sg_low"], row["sg_base"], row["sg_high"]]

    fig = go.Figure()
    fig.add_trace(
        go.Scatter(
            x=suc_x,
            y=sg_y,
            mode="lines+markers+text",
            text=["low","base","high"],
            textposition="top center",
            marker=dict(size=10),
            hovertemplate=(
                f"<b>{row['genotype']}</b><br>"
                "Sucrose x%{x:.2f}<br>"
                "Sustained_growth=%{y:.4f}<extra></extra>"
            ),
            name=row["genotype"]
        )
    )

    fig.update_layout(
        title="Sustained growth vs sucrose",
        xaxis_title="Sucrose multiplier (mSUC_n)",
        yaxis_title="Sustained_growth (steady state)",
        margin=dict(l=40,r=40,t=60,b=40),
        legend_title="Scenario"
    )
    return fig

def plot_delta_vs_WT_curve_interactive(row):
    """
    row contains delta_low_vs_WT, delta_base_vs_WT, delta_high_vs_WT.
    We'll clamp tiny numeric noise to 0 so WT doesn't look 'different from itself'.
    """
    suc_x = [SUCROSE_LEVELS["low"], SUCROSE_LEVELS["base"], SUCROSE_LEVELS["high"]]
    delta_raw = [
        row["delta_low_vs_WT"],
        row["delta_base_vs_WT"],
        row["delta_high_vs_WT"],
    ]
    delta_y = [_clean_delta(v) for v in delta_raw]

    fig = go.Figure()
    fig.add_trace(
        go.Scatter(
            x=suc_x,
            y=delta_y,
            mode="lines+markers+text",
            text=["low","base","high"],
            textposition="top center",
            marker=dict(size=10),
            hovertemplate=(
                f"<b>{row['genotype']}</b><br>"
                "Sucrose x%{x:.2f}<br>"
                "Δ vs WT=%{y:.4f}<extra></extra>"
            ),
            name=row["genotype"]
        )
    )

    fig.add_hline(y=0, line_dash="dash", line_color="gray", opacity=0.5)

    fig.update_layout(
        title="Δ(sustained growth) vs WT across sucrose",
        xaxis_title="Sucrose multiplier (mSUC_n)",
        yaxis_title="Mutant - WT sustained_growth",
        margin=dict(l=40,r=40,t=60,b=40),
        legend_title="Scenario"
    )
    return fig

# ---------- Global scaling plots (low→base, base→high, low→high) ----------

def fit_scaling(df_scaling: pd.DataFrame, xcol: str, ycol: str):
    """
    Fit y = a*x + b and return slope, intercept, R^2, and line coords for plotting.
    df_scaling has sg_low, sg_base, sg_high per genotype.
    """
    X = df_scaling[[xcol]].values
    y = df_scaling[ycol].values
    model = LinearRegression().fit(X, y)

    slope = model.coef_[0]
    intercept = model.intercept_
    r2 = model.score(X, y)

    # line for plotting
    xline = np.linspace(X.min(), X.max(), 100)
    yline = model.predict(xline.reshape(-1, 1))

    line_df = pd.DataFrame({xcol: xline, ycol: yline})
    return slope, intercept, r2, line_df

def make_scaling_fig(df_scaling, xcol, ycol,
                     xtitle, ytitle, title, fit_label_template):
    """
    Make scatter + regression line for one of:
    - low → base
    - base → high
    - low → high
    """
    slope, intercept, r2, line_df = fit_scaling(df_scaling, xcol, ycol)

    fig = go.Figure()

    # add each genotype as its own marker with text label
    for _, r in df_scaling.iterrows():
        fig.add_trace(
            go.Scatter(
                x=[r[xcol]],
                y=[r[ycol]],
                mode="markers+text",
                text=[r["genotype"]],
                textposition="top center",
                marker=dict(size=8),
                hovertemplate=(
                    "<b>%{text}</b><br>"
                    f"{xcol}=%{{x:.4f}}<br>"
                    f"{ycol}=%{{y:.4f}}<extra></extra>"
                ),
                name=r["genotype"]
            )
        )

    # regression line
    fig.add_trace(
        go.Scatter(
            x=line_df[xcol],
            y=line_df[ycol],
            mode="lines",
            line=dict(width=2),
            hoverinfo="skip",
            name="fit"
        )
    )

    fig.update_layout(
        title=title,
        xaxis_title=xtitle,
        yaxis_title=ytitle,
        legend_title="Scenario",
        margin=dict(l=40,r=40,t=60,b=40),
    )

    footer = fit_label_template.format(
        slope=slope,
        intercept=intercept,
        r2=r2
    )

    return fig, footer
