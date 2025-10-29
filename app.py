# app.py
import streamlit as st

st.set_page_config(page_title="PSoup Home", layout="wide")

st.title("🌿 PSoup Modeling Sandbox")
st.markdown("""
Welcome. Pick a page from the sidebar:

### 1. PSoup Explorer (original model)
- Uses the unmodified PSoup equations.
- You can mutate pathway modifiers, generate random mutants, etc.
- We plot sustained growth at 3 sucrose levels (low/base/high).
- We check linearity: low→base, base→high, low→high.
- Empirically this model gives almost perfect straight lines for *all* genotypes,
  meaning sucrose just scales everything and there's basically no G×E.

### 2. PSoup G×E Explorer (modified model)
- Uses a biologically informed extension.
- We add an explicit **sucrose × strigolactone (SL) brake** interaction:
  sucrose relieves SL-mediated bud repression, *but only if the SL pathway is intact*.
- Now 'max2 / smxl6_7_8'-like genotypes (SL pathway broken) can:
  - branch more even at low sucrose,
  - and respond less strongly to extra sucrose,
  - so the slope vs sucrose flattens.
- This should create visible genotype × sucrose interaction (true G×E), which is
  what people actually report in branching / tillering biology.

Use this app to:
- stress-test PSoup as a predictor for breeding,
- show reviewers why the naive algebra cannot produce G×E,
- and show how a minimal biologically grounded change *can* produce G×E.
""")

st.markdown("---")
st.markdown("⬅ Use the left sidebar to switch pages.")
