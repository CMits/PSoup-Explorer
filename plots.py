# plots.py

import plotly.graph_objects as go
import pandas as pd
from analysis import (
    fit_scaling,
    SUCROSE_LEVELS,
)

def make_multi_sustained_growth_plot(all_rows):
    suc_x = [SUCROSE_LEVELS["low"], SUCROSE_LEVELS["base"], SUCROSE_LEVELS["high"]]
    fig = go.Figure()
    for row in all_rows:
        sg_y = [row["sg_low"], row["sg_base"], row["sg_high"]]
        fig.add_trace(
            go.Scatter(
                x=suc_x,
                y=sg_y,
                mode="lines+markers+text",
                text=["low", "base", "high"],
                textposition="top center",
                marker=dict(size=10),
                hovertemplate=(
                    f"<b>{row['genotype']}</b><br>"
                    "Sucrose x%{x:.2f}<br>"
                    "Sustained_growth=%{y:.3f}<extra></extra>"
                ),
                name=row["genotype"]
            )
        )
    fig.update_layout(
        title="Sustained growth vs sucrose (all scenarios)",
        xaxis_title="Sucrose multiplier (mSUC_n)",
        yaxis_title="Sustained_growth",
        margin=dict(l=40, r=40, t=60, b=40),
        legend_title="Scenario"
    )
    return fig


def make_multi_delta_plot(all_rows):
    suc_x = [SUCROSE_LEVELS["low"], SUCROSE_LEVELS["base"], SUCROSE_LEVELS["high"]]
    fig = go.Figure()
    for row in all_rows:
        delta_y = [
            row["delta_low_vs_WT"],
            row["delta_base_vs_WT"],
            row["delta_high_vs_WT"]
        ]
        fig.add_trace(
            go.Scatter(
                x=suc_x,
                y=delta_y,
                mode="lines+markers+text",
                text=["low", "base", "high"],
                textposition="top center",
                marker=dict(size=10),
                hovertemplate=(
                    f"<b>{row['genotype']}</b><br>"
                    "Sucrose x%{x:.2f}<br>"
                    "Δ vs WT=%{y:.3f}<extra></extra>"
                ),
                name=row["genotype"]
            )
        )
    fig.add_hline(y=0, line_dash="dash", line_color="gray", opacity=0.5)
    fig.update_layout(
        title="Δ(sustained growth) vs WT across sucrose (all scenarios)",
        xaxis_title="Sucrose multiplier (mSUC_n)",
        yaxis_title="Mutant - WT sustained_growth",
        margin=dict(l=40, r=40, t=60, b=40),
        legend_title="Scenario"
    )
    return fig


def make_scaling_fig(df_scaling: pd.DataFrame, xcol: str, ycol: str,
                     xtitle: str, ytitle: str, title: str,
                     fit_label_template: str):
    """
    Generic scatter+fit builder for scaling relationships.
    fit_label_template should include {slope}, {intercept}, {r2}.
    """
    slope, intercept, r2, line_df = fit_scaling(df_scaling, xcol, ycol)

    fig = go.Figure()

    # one trace per scenario for color separation
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
                    f"{xcol}=%{{x:.3f}}<br>"
                    f"{ycol}=%{{y:.3f}}<extra></extra>"
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
            name=fit_label_template.format(
                slope=slope, intercept=intercept, r2=r2
            )
        )
    )

    fig.update_layout(
        title=title,
        xaxis_title=xtitle,
        yaxis_title=ytitle,
        legend_title="Scenario",
        margin=dict(l=40, r=40, t=60, b=40),
    )

    footer_text = f"**{ycol} ≈ {slope:.3f}·{xcol} + {intercept:.3f} (R²={r2:.4f})**"
    return fig, footer_text
