# Results snapshot (Blocks 0–7)

## Data

- Samples used: **10159**; Proteins: **2919**; Metabolites: **249**.

## Pairwise protein–metabolite associations

- Tests: **726831**. Significant (Pearson, FDR≤0.05): **498829** (of which |r|≥0.1: **161128**). Spearman sig: **519028**. Overlap: **482975**.

## Sex-stratified associations (age-adjusted)

- Female sig edges: **169605**; Male sig edges: **147521**; Edges with significant F–M difference: **8361**.

## Canonical correlation analysis (age/sex residualized)

- First canonical correlation: **0.974** (perm p = **0.004975**). Subsequent components also significant.

## Network (FDR≤0.01, |r|≥0.20)

- Edges: **27229**; Nodes: **1456**; Modules: **5**; Modularity: **0.264**.

## Files included in this bundle

- `key_numbers.csv` (this table)

- `block5_assoc_summary.csv`, `block5a_unadjusted_assoc_summary.csv`

- `block5c_sex_stratified_summary.csv`, `block5c_edges_top20k_by_pdiff.csv`

- `block6_cca_summary.csv`, top loadings under `results/block6_toploadings/`

- `block7_network_summary.csv`, `block7_module_summary.csv`

- Optional figures: main network PNG/SVG if present

