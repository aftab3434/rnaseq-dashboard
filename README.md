# RNA-seq Results Dashboard (V2)

An app-like, interactive dashboard for exploring the results of the companion
analysis project **[rnaseq-differential-expression](../rnaseq-differential-expression)** (V1).

**Open it:** double-click `index.html` - it runs entirely in the browser, no
install, no internet needed. Works on any laptop or phone.

## What is inside

| Tab | Contents |
|---|---|
| Mouse mammary: lactation (GSE60450) | Interactive volcano plot, PCA, top-50 heatmap, pathway enrichment, searchable table of all 1,356 significant genes |
| Human airway: dexamethasone (airway) | Interactive volcano plot, PCA, searchable table of all 1,151 significant genes |
| How to read this | Plain-language guide and honest caveats |

Everything is bundled into the single HTML file (plots, data, the Plotly
library), which is why the file is ~5 MB.

## What it is not

- It does **not** recompute the analysis and does **not** accept uploads. It
  presents the exact results produced by the V1 pipeline (DESeq2 statistics).
- The mouse dataset has 2 replicates per group: treat those results as a
  demonstration of the analysis workflow, not strong new biology. The airway
  dataset (4 per group) is the statistically stronger example.

## Rebuild from the V1 results (optional, for the author)

```bash
pip install pandas numpy
python build_dashboard.py --v1 ../rnaseq-differential-expression --out index.html
```

`plotly.min.js` (MIT, Plotly Inc.) and `data/airway_symbols.json` (gene symbol
map, cached from mygene.info) are committed so the build works offline.

Data sources: GSE60450 (Fu et al. 2015, *Nature Cell Biology*) and GSE52778
(Himes et al. 2014) via recount3 (SRP033351, counts converted with the recount3
`compute_read_counts` method).

## License

MIT. Plotly.js is MIT (Plotly, Inc.). Please cite the original dataset papers if reused.
