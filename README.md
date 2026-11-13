# rscripts

Standalone R scripts for post-processing LC-MS proteomics search results.

| Script | Purpose |
|---|---|
| `maxQuantReport.R` | Summary report from MaxQuant `txt/` output |
| `analyzeMaxQuant_MBR_noimputation.R` | MaxQuant LFQ analysis with match-between-runs, no imputation |
| `analyzeMaxQuant_MBR_withimputation.R` | Same, with missing-value imputation |
| `plotMqQcSummary.R` | QC plots from MaxQuant `summary.txt` (see also [mqqc](https://github.com/heuselm/mqqc)) |
| `diannReport.R` | Summary report from DIA-NN output |
| `compareDiannResultSets.R` | Compare identifications/quantities between DIA-NN result sets |

Scripts are self-contained; open in RStudio (`rscripts.Rproj`), adjust the input paths at the top,
and run. For DIA-NN-based differential testing see [DiffTestR](https://github.com/heuselm/DiffTestR).

## License

MIT — see [LICENSE](LICENSE).

---
<sub>README optimized with Claude Opus 5.5.</sub>
