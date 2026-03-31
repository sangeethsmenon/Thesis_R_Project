# Block5b Results Summary (protein-first modules + enrichment)

Generated: 2026-01-05 02:37:06

## Inputs checked

All expected inputs were found.

## Protein–protein network (PP)

Thresholds (reported): q ≤ 0.01, |r| ≥ 0.3.
PP edges retained: 380,024.
PP nodes (unique proteins in PP edges): 2,369.

## Protein modules (Louvain on PP graph)

Number of modules: 13.
Proteins assigned to modules: 2,369.
Module size (min / median / max): 2 / 3 / 1,099.

Top modules by size (see `results/block5b_results_summary_tables/top_modules_by_size.csv` for the table).

## PM edges mapped onto protein modules

PM edges with module annotation: 5,188.
Unique metabolites represented: 245.
Unique proteins represented: 257.
Unique modules represented: 3.

## Metabolite → module enrichment (Fisher)

Significance cutoff (reported): q ≤ 0.05.
Significant metabolite–module pairs: 232.
Unique significant metabolites: 232.
Unique significant modules: 3.

Top signals table: `results/block5b_results_summary_tables/top_metabolite_module_enrichment.csv`.

## ID mapping (SYMBOL → ENTREZ)

Symbols in mapping table: 2,369.
Symbols mapped to an ENTREZ ID: 2,341.
Mapping rate: 98.8%.

## GO enrichment (clusterProfiler)

Ontology (reported): BP.
Modules with ≥1 GO term (in results table): 5.
Total GO terms reported (all modules combined): 2,732.
Top GO terms per module table: `results/block5b_results_summary_tables/top_GO_terms_by_module.csv`.

## KEGG enrichment (clusterProfiler)

Modules with ≥1 KEGG pathway (in results table): 5.
Total KEGG pathways reported (all modules combined): 189.
Top KEGG pathways per module table: `results/block5b_results_summary_tables/top_KEGG_pathways_by_module.csv`.

## Integrative summary table

Rows in `block5b_module_integrative_summary.csv`: 13.
Distinct modules in integrative summary: 13.

## Notes for thesis writing

Use this file for the Results chapter narrative (counts + headline findings), and move full enrichment tables and full module tables to the Appendix. The small tables in `results/block5b_results_summary_tables/` are already in ‘top-only’ form for main-text reporting.
