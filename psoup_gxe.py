# psoup_gxe.py
from typing import Dict

class PSoupSimulatorGxE:
    """
    Same structure as PSoupSimulator, but Bud_release encodes sucrose × SL interaction.
    We also expose mSL_pathway_intact in gene_vals so genotypes like max2/smxl6_7_8
    can flatten their sucrose response.
    """

    def __init__(self, tmax: int = 100, threshold: float = 1e-4):
        self.TMAX = tmax
        self.THRESHOLD = threshold

        # Gene multipliers
        self.default_gene_vals = {
            'mSL_n': 1.0,
            'mBRC1_2_n': 1.0,
            'mCK_n': 1.0,
            'mAux_n': 1.0,
            'mABA_n': 1.0,
            'mPIN_1_3_4_7_n': 1.0,
            'mAux_Bud_n': 1.0,
            'mSUC_n': 1.0,
            'mGA_n': 1.0,
            'mDecap_signal_n': 1.0,
            'mCXE_n': 1.0,
            'mSMXL6_7_8_n': 1.0,
            'mFHY3_FAR1_n': 1.0,
            'mSPL9_15_n': 1.0,
            'mMAX2_n': 1.0,
            'mD14_n': 1.0,
            'mD27_n': 1.0,
            'mMAX3_n': 1.0,
            'mMAX4_n': 1.0,
            'mMAX1_n': 1.0,
            'mLBO_n': 1.0,
            'mCLAMT_n': 1.0,
            'mPhyB_n': 1.0,
            'mPIFs_n': 1.0,
            'mHB21_40_53_n': 1.0,
            'mNCED3_n': 1.0,

            # NEW dial:
            'mSL_pathway_intact': 1.0,
        }

        # Node states init
        self.default_node_vals = {
            'mSL': 1.0,
            'mBRC1_2': 1.0,
            'mCK': 1.0,
            'mAux_Shoot': 1.0,
            'mABA': 1.0,
            'mBud_release': 1.0,
            'mPAT_Bud': 1.0,
            'mAux_Bud': 1.0,
            'mSUC': 1.0,
            'mSustained_growth': 1.0,
            'mGA': 1.0,
            'mDecap_signal': 1.0,
            'mCXE': 1.0,
            'mSMXL6_7_8': 1.0,
            'mFHY3_FAR1': 1.0,
            'mSPL9_15': 1.0,
            'mSMXL6_7_8_Transcription': 1.0,
            'mLow_R_FR': 1.0,
            'mMAX2': 1.0,
            'mD14': 1.0,
            'mPerception_SL': 1.0,
            'mD27': 1.0,
            'mMAX3': 1.0,
            'mMAX4': 1.0,
            'mMAX1': 1.0,
            'mLBO': 1.0,
            'mCLAMT': 1.0,
            'mPhyB': 1.0,
            'mHigh_R_FR': 1.0,
            'mPIFs': 1.0,
            'mHB21_40_53': 1.0,
            'mNCED3': 1.0,
        }

        self.default_exo_vals = {k: 0.0 for k in self.default_node_vals.keys()}

    def _sucrose_relief_function(self, sucrose_val: float) -> float:
        # saturating relief of SL repression by sucrose
        K = 1.0
        return sucrose_val / (K + sucrose_val)

    def step_next(self,
                  gene_vals: Dict[str, float],
                  node_vals: Dict[str, float],
                  exo_vals: Dict[str, float]) -> Dict[str, float]:

        old = node_vals.copy()
        new = {}

        # same updates as PSoupSimulator EXCEPT Bud_release
        new['mSL'] = ((2/(1 + old['mCXE'])) * (gene_vals['mSL_n'])) * (old['mLBO']) + exo_vals['mSL']

        new['mBRC1_2'] = (
            (2 * ((old['mPIFs'] + old['mSPL9_15'])/2)
             / (1 + (old['mCK'] + old['mSUC'] + old['mPhyB'])/3))
            * (gene_vals['mBRC1_2_n'])
            + exo_vals['mBRC1_2']
        )

        new['mCK'] = (
            (2 * (old['mSUC'])
             / (1 + (old['mAux_Shoot'] + old['mPerception_SL'])/2))
            * (gene_vals['mCK_n'])
            + exo_vals['mCK']
        )

        new['mAux_Shoot'] = (
            (2/(1 + old['mDecap_signal'])) * (gene_vals['mAux_n'])
            + exo_vals['mAux_Shoot']
        )

        new['mABA'] = (old['mNCED3']) * (gene_vals['mABA_n']) + exo_vals['mABA']

        # -------- MODIFIED Bud_release (G×E with sucrose-specific bonus) --------
        suc = old['mSUC']
        ck  = old['mCK']

        # (1) Carbon baseline drive
        base_drive = (suc + ck) / 2.0

        # (2) Sucrose hormonal bonus, only if SL pathway intact
        intact = gene_vals.get('mSL_pathway_intact', 1.0)  # WT=1, SL mutant=0
        suc_relief = suc / (1.0 + suc)  # 0→1 saturating curve with sucrose
        bonus_strength = 0.5            # tune this number if effect too weak/too strong
        hormonal_bonus = bonus_strength * intact * suc_relief
        growth_drive = base_drive * (1.0 + hormonal_bonus)

        # (3) Brakes
        ABA_brake = (old['mABA'] + old['mBRC1_2']) / 2.0

        SL_brake_base = old['mPerception_SL']
        SL_brake_effective = SL_brake_base * intact * (1.0 - suc_relief)
        # WT (intact=1): big SL brake at low sucrose, almost none at high
        # SL mutant (intact=0): basically zero SL brake always

        # (4) Final Bud_release
        new['mBud_release'] = (
            2.0 * growth_drive /
            (1.0 + ABA_brake + SL_brake_effective)
        ) + exo_vals['mBud_release']
        # ------------------------------------------------------------------------


        new['mPAT_Bud'] = (
            (2 * ((old['mAux_Bud'] + old['mCK'] + old['mDecap_signal'])/3)
             / (1 + old['mAux_Shoot']))
            * (gene_vals['mPIN_1_3_4_7_n'])
            + exo_vals['mPAT_Bud']
        )

        new['mAux_Bud'] = (gene_vals['mAux_Bud_n']) + exo_vals['mAux_Bud']

        new['mSUC'] = (old['mDecap_signal']) * (gene_vals['mSUC_n']) + exo_vals['mSUC']

        new['mSustained_growth'] = (
            (2 * (min(old['mGA'], old['mPAT_Bud']))/(1 + old['mSUC']))
            * (old['mBud_release'])
            + exo_vals['mSustained_growth']
        )

        new['mGA'] = (old['mAux_Bud']) * (gene_vals['mGA_n']) + exo_vals['mGA']

        new['mDecap_signal'] = (gene_vals['mDecap_signal_n']) + exo_vals['mDecap_signal']

        new['mCXE'] = (gene_vals['mCXE_n']) + exo_vals['mCXE']

        new['mSMXL6_7_8'] = (
            (2 * ((old['mSMXL6_7_8_Transcription'] + old['mFHY3_FAR1'])/2)
             / (1 + old['mPerception_SL']))
            * (gene_vals['mSMXL6_7_8_n'])
            + exo_vals['mSMXL6_7_8']
        )

        new['mFHY3_FAR1'] = (
            (2 * (old['mHigh_R_FR'])/(1 + old['mLow_R_FR']))
            * (gene_vals['mFHY3_FAR1_n'])
            + exo_vals['mFHY3_FAR1']
        )

        new['mSPL9_15'] = (
            (2/(1 + (old['mFHY3_FAR1'] + old['mSMXL6_7_8'])/2))
            * (gene_vals['mSPL9_15_n'])
            + exo_vals['mSPL9_15']
        )

        new['mSMXL6_7_8_Transcription'] = (
            2/(1 + old['mSMXL6_7_8']) + exo_vals['mSMXL6_7_8_Transcription']
        )

        new['mLow_R_FR'] = 1 + exo_vals['mLow_R_FR']

        new['mMAX2'] = (
            (2/(1 + (old['mSUC'] + old['mPhyB'])/2))
            * (gene_vals['mMAX2_n'])
            + exo_vals['mMAX2']
        )

        new['mD14'] = (gene_vals['mD14_n']) + exo_vals['mD14']

        new['mPerception_SL'] = (
            min(old['mD14'], min(old['mMAX2'], old['mSL'])) + exo_vals['mPerception_SL']
        )

        new['mD27'] = (gene_vals['mD27_n']) + exo_vals['mD27']

        new['mMAX3'] = (
            ((2 * (old['mAux_Shoot'])/(1 + old['mPerception_SL']))
             * (gene_vals['mMAX3_n']))
            * (old['mD27'])
            + exo_vals['mMAX3']
        )

        new['mMAX4'] = (
            ((2 * (old['mAux_Shoot'])/(1 + (old['mPerception_SL'] + old['mPhyB'])/2))
             * (gene_vals['mMAX4_n']))
            * (old['mMAX3'])
            + exo_vals['mMAX4']
        )

        new['mMAX1'] = (
            ((2/(1 + old['mPerception_SL'])) * (gene_vals['mMAX1_n']))
            * (old['mMAX4'])
            + exo_vals['mMAX1']
        )

        new['mLBO'] = ((gene_vals['mLBO_n'])) * (old['mCLAMT']) + exo_vals['mLBO']

        new['mCLAMT'] = ((gene_vals['mCLAMT_n'])) * (old['mMAX1']) + exo_vals['mCLAMT']

        new['mPhyB'] = (
            (2 * (old['mHigh_R_FR'])/(1 + old['mLow_R_FR']))
            * (gene_vals['mPhyB_n'])
            + exo_vals['mPhyB']
        )

        new['mHigh_R_FR'] = 1 + exo_vals['mHigh_R_FR']

        new['mPIFs'] = (
            (2/(1 + old['mHigh_R_FR']))
            * (gene_vals['mPIFs_n'])
            + exo_vals['mPIFs']
        )

        new['mHB21_40_53'] = (
            min(old['mBRC1_2'], old['mLow_R_FR'])
            * (gene_vals['mHB21_40_53_n'])
            + exo_vals['mHB21_40_53']
        )

        new['mNCED3'] = (
            (old['mHB21_40_53'])
            * (gene_vals['mNCED3_n'])
            + exo_vals['mNCED3']
        )

        return new

    def _converged(self, old_vals: Dict[str, float], new_vals: Dict[str, float]) -> bool:
        for k in old_vals:
            if abs(old_vals[k] - new_vals[k]) > self.THRESHOLD:
                return False
        return True

    def simulate_to_steady_state(self,
                                 gene_vals: Dict[str, float] = None,
                                 node_vals: Dict[str, float] = None,
                                 exo_vals: Dict[str, float] = None,
                                 verbose: bool = False) -> Dict[str, float]:

        if gene_vals is None:
            gene_vals = self.default_gene_vals.copy()
        if node_vals is None:
            node_vals = self.default_node_vals.copy()
        if exo_vals is None:
            exo_vals = self.default_exo_vals.copy()

        current = node_vals.copy()

        for t in range(self.TMAX):
            prev = current
            current = self.step_next(gene_vals, current, exo_vals)
            if t > 0 and self._converged(prev, current):
                if verbose:
                    print(f"GxE PSoup converged at t={t}")
                break

        return current

    def get_sustained_growth(self,
                             gene_vals: Dict[str, float] = None,
                             node_vals: Dict[str, float] = None,
                             exo_vals: Dict[str, float] = None) -> float:
        final_vals = self.simulate_to_steady_state(
            gene_vals=gene_vals,
            node_vals=node_vals,
            exo_vals=exo_vals
        )
        return final_vals['mSustained_growth']

    def sustained_growth(self,
                         gene_vals: Dict[str, float] = None,
                         node_vals: Dict[str, float] = None,
                         exo_vals: Dict[str, float] = None) -> float:
        """Compatibility with analysis.run_scenario()"""
        return self.get_sustained_growth(
            gene_vals=gene_vals,
            node_vals=node_vals,
            exo_vals=exo_vals
        )
