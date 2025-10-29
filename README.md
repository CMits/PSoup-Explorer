# PSoup G×Sucrose Explorer

This Streamlit app lets you explore how different modifier combinations
(genotype-like scenarios) affect sustained growth across sucrose environments
in the PSoup branching network model.

## How it works

- The PSoup equations are taken directly from the exported model.
- A "scenario" is just a set of modifier values (like `mMAX3_n = 0`, `mBRC1_2_n = 0.2`, etc.).
- For each scenario we:
  - run steady state at sucrose 0.5×, 1.0×, 1.5× (implemented by scaling `mSUC_n`)
  - record `mSustained_growth`
  - compare to WT (reference values provided)
- We then:
  - plot sustained growth vs sucrose,
  - plot mutant − WT vs sucrose,
  - plot all scenarios together in low→base and base→high scaling space
    and re-fit a global line (slope/intercept/R²).

This tells you if sucrose is acting as a global linear gain (points on one line, R²≈1)
or if you've found genotype×sucrose interaction (points off-line, R² drops).

## Usage

1. Install deps:
   ```bash
   pip install -r requirements.txt
