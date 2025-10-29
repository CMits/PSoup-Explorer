# scenarios.py
import random
from typing import Dict, List, Type
from analysis import simulate_genotype

# all tunable gene multipliers, including the new G×E control dial
MODIFIER_FIELDS = [
    "mSL_n",
    "mBRC1_2_n",
    "mCK_n",
    "mAux_n",
    "mABA_n",
    "mPIN_1_3_4_7_n",
    "mAux_Bud_n",
    "mSUC_n",
    "mGA_n",
    "mDecap_signal_n",
    "mCXE_n",
    "mSMXL6_7_8_n",
    "mFHY3_FAR1_n",
    "mSPL9_15_n",
    "mMAX2_n",
    "mD14_n",
    "mD27_n",
    "mMAX3_n",
    "mMAX4_n",
    "mMAX1_n",
    "mLBO_n",
    "mCLAMT_n",
    "mPhyB_n",
    "mPIFs_n",
    "mHB21_40_53_n",
    "mNCED3_n",
    # new dial usable only in G×E, but harmless in original:
    "mSL_pathway_intact",
]

def parse_modifiers_from_text_inputs(raw_mods: Dict[str, str]) -> Dict[str, float]:
    out = {}
    for k, v in raw_mods.items():
        v_clean = v.strip()
        if v_clean == "":
            out[k] = 1.0
        else:
            try:
                out[k] = float(v_clean)
            except ValueError:
                out[k] = 1.0
    return out

def _random_all_modifiers(min_val: float, max_val: float) -> Dict[str, float]:
    return {f: random.uniform(min_val, max_val) for f in MODIFIER_FIELDS}

def _random_sparse_modifiers(min_val: float, max_val: float, p_change: float) -> Dict[str, float]:
    mods = {}
    for f in MODIFIER_FIELDS:
        if random.random() < p_change:
            mods[f] = random.uniform(min_val, max_val)
        else:
            mods[f] = 1.0
    return mods

def _random_mutant_single() -> Dict[str, float]:
    mods = {f: 1.0 for f in MODIFIER_FIELDS}
    knocked = random.choice(MODIFIER_FIELDS)
    mods[knocked] = 0.0
    return mods

def _random_mutant_double() -> Dict[str, float]:
    mods = {f: 1.0 for f in MODIFIER_FIELDS}
    knocked_two = random.sample(MODIFIER_FIELDS, k=2)
    for k in knocked_two:
        mods[k] = 0.0
    return mods

def run_and_package(sim,
                    name: str,
                    modifiers: Dict[str, float]) -> Dict[str, float]:
    """
    sim: simulator instance (PSoupSimulator or PSoupSimulatorGxE)
    name: scenario label
    modifiers: dict of per-gene multipliers
    """
    result_row = simulate_genotype(sim, name, modifiers)
    # attach param_ columns for transparency/debug
    for k, v in modifiers.items():
        result_row[f"param_{k}"] = v
    return result_row

def bulk_generate(sim_class: Type,
                  mode: str,
                  n: int,
                  min_val: float,
                  max_val: float,
                  p_change: float,
                  existing_count: int) -> List[Dict[str, float]]:
    """
    sim_class: class, not instance (PSoupSimulator or PSoupSimulatorGxE)
    mode: one of
          - "all_modifiers_change"
          - "sparse_modifiers_change"
          - "mutant_single_knockout"
          - "mutant_double_knockout"
    """
    rows = []
    for i in range(n):
        if mode == "all_modifiers_change":
            mods = _random_all_modifiers(min_val, max_val)
        elif mode == "sparse_modifiers_change":
            mods = _random_sparse_modifiers(min_val, max_val, p_change)
        elif mode == "mutant_single_knockout":
            mods = _random_mutant_single()
        elif mode == "mutant_double_knockout":
            mods = _random_mutant_double()
        else:
            mods = {f: 1.0 for f in MODIFIER_FIELDS}

        sim = sim_class()
        auto_name = f"rand_{mode}_{existing_count + i}"
        row = run_and_package(sim, auto_name, mods)
        rows.append(row)

    return rows
