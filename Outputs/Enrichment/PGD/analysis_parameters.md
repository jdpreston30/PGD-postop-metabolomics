# Mummichog Analysis Parameters

**Analysis Date:** 2025-12-02

**Database:** hsa_mfn

**MetaboAnalystR 'Set' Function Outputs:**
- SetPeakFormat: mpr
- SetPeakEnrichMethod: mummichog (v2)
- SetMummichogPvalFromPercent: 0.1 (top 10% of peaks)

**Instrument Parameters (UpdateInstrumentParameters):**
- instrumentOpt: 5
- msModeOpt: mixed
- force_primary_ion: yes
- rt_frac: 0.02

**Analysis Parameters:**
- Peak filtering method: Top 10% of peaks (dynamic)
- Peak filtering threshold (rounded): 0.05
- Peak filtering threshold (precise): 0.0744428039061391
- Peaks analyzed: 848 out of 16283
- Pathways analyzed: 73
- Significant pathways (p < 0.05): 3
- Pathway p-values range: 0.002308 to 0.9954
- Pathway FDR: Not calculated (using raw p-values)
- Pathway enrichment FDR threshold: 0.05 (fixed)
- Minimum pathway size: 3
- Background permutations: 100

**Input Data:**
- Number of features: 16283
- Output directory: /Users/jdp2019/Desktop/PGD-postop-metabolomics/Outputs/Enrichment/PGD

