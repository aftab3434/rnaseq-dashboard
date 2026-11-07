#!/usr/bin/env python3
"""Build the V2 RNA-seq results dashboard: a single self-contained index.html.

Reads the RESULTS committed by the V1 pipeline (rnaseq-differential-expression)
and packages them into one interactive HTML file that works offline in any
browser - no install, no server. It does NOT recompute anything and does not
accept arbitrary uploads; it presents the bundled V1 results.

Usage:
    python build_dashboard.py [--v1 ../rnaseq-differential-expression] [--out index.html]
Requires: pandas, numpy, and plotly.min.js next to this script (downloaded by
01_fetch_assets or already committed).
"""
import argparse, base64, json
from pathlib import Path
import numpy as np
import pandas as pd

rng = np.random.default_rng(42)
PADJ, LFC = 0.05, 1.0
CLOUD_N = 4000  # downsampled non-significant background points per volcano

def volcano_payload(df, gene_col, symbol_map=None):
    df = df.dropna(subset=["padj"]).copy()
    df["neglog"] = -np.log10(df["padj"])
    sig = (df.padj < PADJ) & (df.log2FoldChange.abs() >= LFC)
    keep = sig | (df.padj < 1e-20)
    rest = df[~keep]
    if len(rest) > CLOUD_N:
        rest = rest.iloc[rng.choice(len(rest), CLOUD_N, replace=False)]
    out = {}
    for name, sub, color in [("not significant", rest, "#b6bec9"),
                             ("significant", df[keep & sig], "#C44E52")]:
        genes = sub[gene_col].astype(str)
        if symbol_map:
            genes = [symbol_map.get(g.split(".")[0], g.split(".")[0]) for g in genes]
        out[name] = {
            "x": sub.log2FoldChange.round(4).tolist(),
            "y": sub.neglog.round(4).tolist(),
            "gene": list(genes),
            "color": color,
        }
    return out

def table_payload(df, gene_col, symbol_map=None, n=None):
    d = df.sort_values("padj")
    if n: d = d.head(n)
    rows = []
    for _, r in d.iterrows():
        g = str(r[gene_col])
        sym = symbol_map.get(g.split(".")[0], "") if symbol_map else str(r.get("symbol", ""))
        rows.append([sym or g.split(".")[0], round(float(r.log2FoldChange), 3),
                     round(float(r.baseMean), 1), f"{r.padj:.3e}"])
    return rows

def b64img(path):
    return "data:image/png;base64," + base64.b64encode(Path(path).read_bytes()).decode()

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--v1", default="../rnaseq-differential-expression")
    ap.add_argument("--out", default="index.html")
    a = ap.parse_args()
    v1 = Path(a.v1)
    here = Path(__file__).parent

    plotly_js = (here / "plotly.min.js").read_text()
    airway_symbols = json.loads((here / "data" / "airway_symbols.json").read_text())

    mouse_full = pd.read_csv(v1 / "results/tables/de_results_lactating_vs_pregnant_full.csv")
    mouse_sig  = pd.read_csv(v1 / "results/tables/de_genes_significant.csv")
    airway_full = pd.read_csv(v1 / "results/airway/airway_dex_vs_untreated_full.csv")
    airway_sig  = pd.read_csv(v1 / "results/airway/airway_significant_genes.csv")

    data = {
        "mouse": {
            "volcano": volcano_payload(mouse_full, "symbol"),
            "table": table_payload(mouse_sig, "symbol"),
            "metrics": {
                "Genes tested": f"{mouse_full.pvalue.notna().sum():,}",
                "With adjusted p-values": f"{mouse_full.padj.notna().sum():,}",
                "Significant genes": f"{len(mouse_sig):,}",
                "Up in lactation": f"{int((mouse_sig.log2FoldChange>0).sum()):,}",
                "Down in lactation": f"{int((mouse_sig.log2FoldChange<0).sum()):,}",
            },
        },
        "airway": {
            "volcano": volcano_payload(airway_full, "gene_id", airway_symbols),
            "table": table_payload(airway_sig, "gene_id", airway_symbols),
            "metrics": {
                "Genes tested": f"{airway_full.pvalue.notna().sum():,}",
                "With adjusted p-values": f"{airway_full.padj.notna().sum():,}",
                "Significant genes": f"{len(airway_sig):,}",
                "Up with dexamethasone": f"{int((airway_sig.log2FoldChange>0).sum()):,}",
                "Down with dexamethasone": f"{int((airway_sig.log2FoldChange<0).sum()):,}",
            },
        },
    }
    imgs = {
        "mouse_pca": b64img(v1 / "results/figures/pca_plot.png"),
        "mouse_heatmap": b64img(v1 / "results/figures/heatmap_top50.png"),
        "mouse_enrich": b64img(v1 / "results/figures/enrichment_dotplot.png"),
        "airway_pca": b64img(v1 / "results/airway/pca_plot.png"),
    }
    html = TEMPLATE.replace("__DATA__", json.dumps(data)).replace("__IMGS__", json.dumps(imgs)).replace("__PLOTLY__", plotly_js)
    Path(a.out).write_text(html)
    print(f"wrote {a.out} ({len(html)/1e6:.1f} MB)")

TEMPLATE = r"""<!DOCTYPE html>
<html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>RNA-seq Results Dashboard</title>
<style>
:root{--bg:#0f141b;--card:#182030;--ink:#e8edf4;--muted:#93a1b3;--accent:#4da3ff;--red:#e0605e;--blue:#5e9ae0}
*{box-sizing:border-box}body{margin:0;background:var(--bg);color:var(--ink);font-family:-apple-system,'Segoe UI',Roboto,Helvetica,Arial,sans-serif}
header{padding:28px 32px 12px}h1{margin:0 0 4px;font-size:26px}header p{margin:0;color:var(--muted);max-width:900px}
.tabs{display:flex;gap:8px;padding:16px 32px 0;flex-wrap:wrap}
.tab{padding:10px 18px;border-radius:10px 10px 0 0;background:#141c29;color:var(--muted);cursor:pointer;border:1px solid #232f42;border-bottom:none;font-size:15px}
.tab.active{background:var(--card);color:var(--ink)}
.panel{display:none;padding:24px 32px 40px;background:var(--card);border-top:1px solid #232f42}
.panel.active{display:block}
.cards{display:flex;gap:14px;flex-wrap:wrap;margin-bottom:18px}
.card{background:#101826;border:1px solid #232f42;border-radius:12px;padding:14px 20px;min-width:170px}
.card .v{font-size:24px;font-weight:650}.card .k{color:var(--muted);font-size:13px;margin-top:2px}
.grid{display:grid;grid-template-columns:1fr;gap:22px}@media(min-width:1100px){.grid{grid-template-columns:1fr 1fr}}
.box{background:#101826;border:1px solid #232f42;border-radius:12px;padding:16px;overflow-x:auto}
.box h3{margin:0 0 8px;font-size:16px}.box img{width:100%;border-radius:8px;background:#fff}
.full{grid-column:1/-1}
input[type=search]{width:100%;max-width:420px;padding:9px 12px;border-radius:8px;border:1px solid #2c3a52;background:#0d1420;color:var(--ink);margin-bottom:10px}
table{border-collapse:collapse;width:100%;font-size:14px}
th,td{padding:7px 10px;text-align:left;border-bottom:1px solid #22304a}
th{cursor:pointer;color:var(--accent);user-select:none;position:sticky;top:0;background:#101826}
.wrap{max-height:430px;overflow:auto;border-radius:8px}
.note{color:var(--muted);font-size:13px;margin-top:6px}
h2{font-size:20px}li{margin:6px 0;max-width:880px}.pill{display:inline-block;background:#101c2e;border:1px solid #2c3a52;border-radius:999px;padding:2px 10px;font-size:12px;color:var(--muted);margin-left:8px}
</style></head><body>
<header><h1>RNA-seq Differential Expression Dashboard</h1>
<p>Interactive view of the results produced by the companion analysis project
(<b>rnaseq-differential-expression</b>). Two real public datasets, analyzed with DESeq2
statistics. Everything is bundled in this one file - it works offline.</p></header>
<div class="tabs">
  <div class="tab active" data-t="mouse">Mouse mammary: lactation (GSE60450)</div>
  <div class="tab" data-t="airway">Human airway: dexamethasone (airway)</div>
  <div class="tab" data-t="guide">How to read this</div>
</div>

<div class="panel active" id="p-mouse">
  <div class="cards" id="m-cards"></div>
  <div class="note" style="margin:-8px 0 14px">"Significant" throughout = adjusted p-value (FDR) &lt; 0.05 AND |log2 fold change| &ge; 1. "Genes tested" = genes DESeq2 computed a p-value for; "With adjusted p-values" = retained after DESeq2 independent filtering.</div>
  <div class="grid">
    <div class="box full"><h3>Volcano plot - every gene tested (hover for gene names, drag to zoom)</h3><div id="m-volcano" style="height:520px"></div>
      <div class="note">Each dot is one gene. Right of the right dashed line + above the horizontal line = significantly higher in lactation; left side = lower. Grey cloud is a 4,000-point sample of non-significant genes.</div></div>
    <div class="box"><h3>PCA - samples group by stage</h3><img id="m-pca" alt="PCA plot"></div>
    <div class="box"><h3>Top 50 genes - expression pattern across samples</h3><img id="m-heat" alt="heatmap"></div>
    <div class="box full"><h3>Pathway enrichment - what the changed genes do</h3><img id="m-enrich" alt="enrichment" style="background:#fff"></div>
    <div class="box full"><h3>All 1,356 significant genes</h3>
      <input type="search" id="m-search" placeholder="Search a gene name, e.g. Csn1s2b">
      <div class="wrap"><table id="m-table"></table></div></div>
  </div>
</div>

<div class="panel" id="p-airway">
  <div class="cards" id="a-cards"></div>
  <div class="note" style="margin:-8px 0 14px">"Significant" throughout = adjusted p-value (FDR) &lt; 0.05 AND |log2 fold change| &ge; 1. "Genes tested" = genes DESeq2 computed a p-value for; "With adjusted p-values" = retained after DESeq2 independent filtering.</div>
  <div class="grid">
    <div class="box full"><h3>Volcano plot - dexamethasone vs untreated (4 vs 4)</h3><div id="a-volcano" style="height:520px"></div>
      <div class="note">The community-standard DESeq2 benchmark dataset: 4 cell lines, each measured with and without the drug (a paired design, not 8 independent donors). Same statistics as the mouse tab.</div></div>
    <div class="box"><h3>PCA - treatment separates on PC1</h3><img id="a-pca" alt="airway PCA"></div>
    <div class="box full"><h3>All 1,123 significant genes</h3>
      <input type="search" id="a-search" placeholder="Search a gene name, e.g. DUSP1">
      <div class="wrap"><table id="a-table"></table></div></div>
  </div>
</div>

<div class="panel" id="p-guide">
  <h2>How to read this dashboard</h2>
  <ul>
    <li><b>What this is:</b> the finished results of an RNA-seq data-analysis project. RNA-seq measures how active every gene is in a sample. This project compares gene activity between two biological states and finds which genes change.</li>
    <li><b>Dataset 1 (mouse mammary):</b> compares breast cells during milk production (lactation) vs late pregnancy. Milk-protein genes switch on; cell-division genes switch off.</li>
    <li><b>Dataset 2 (human airway):</b> the standard teaching dataset for this method - 4 human airway cell lines, each measured with and without a steroid drug (dexamethasone). A paired design: 4 independent cell lines, not 8 independent people.</li>
    <li><b>Volcano plot:</b> one dot per gene. Far left/right = big change; high up = strong statistical evidence the change is not random noise. Strong evidence is not the same as proven biological importance. Hover any dot to see the gene.</li>
    <li><b>PCA plot:</b> each point is one sample. Samples from the same condition clustering together is what we expect when the condition drives large differences - reassuring, but not by itself proof the data is error-free.</li>
    <li><b>Tables:</b> every gene that passed the significance threshold (adjusted p-value &lt; 0.05 and at least 2-fold change). Use the search box.</li>
    <li><b>Honest caveats:</b> the mouse dataset has only 2 samples per group, so its results are a demonstration of the analysis skill set, not strong new biology. The airway dataset (4 per group) is the statistically stronger example. This dashboard presents bundled results from the analysis pipeline; it does not recompute anything or accept new uploads.</li>
  </ul>
  <p class="note">Built from the rnaseq-differential-expression project (V1). Datasets: GSE60450 (Fu et al. 2015) and GSE52778 via recount3 (Himes et al. 2014). Statistics: PyDESeq2.</p>
</div>

<script>__PLOTLY__</script>
<script>
const DATA = __DATA__;
const IMGS = __IMGS__;
function cards(id, m){document.getElementById(id).innerHTML=Object.entries(m).map(([k,v])=>`<div class="card"><div class="v">${v}</div><div class="k">${k}</div></div>`).join('');}
function volcano(id, v, xlab){
  const traces=[["not significant",10],["significant",14]].map(([nm,sz])=>({
    x:v[nm].x,y:v[nm].y,text:v[nm].gene,mode:'markers',type:'scattergl',name:nm,
    marker:{size:v[nm].gene.length>3000?4:6,color:v[nm].color,opacity:nm==='not significant'?0.35:0.65},
    hovertemplate:'%{text}<br>log2FC %{x:.2f}<br>-log10 padj %{y:.1f}<extra>'+nm+'</extra>'}));
  Plotly.newPlot(id,traces,{paper_bgcolor:'#101826',plot_bgcolor:'#101826',
    font:{color:'#e8edf4'},margin:{t:20,b:60,l:70,r:20},
    xaxis:{title:xlab,gridcolor:'#22304a',zeroline:false},
    yaxis:{title:'-log10 adjusted p-value',gridcolor:'#22304a'},
    legend:{orientation:'h',y:1.1},
    shapes:[{type:'line',x0:1,x1:1,y0:0,y1:1,yref:'paper',line:{dash:'dash',color:'#93a1b3'}},
            {type:'line',x0:-1,x1:-1,y0:0,y1:1,yref:'paper',line:{dash:'dash',color:'#93a1b3'}},
            {type:'line',x0:0,x1:1,xref:'paper',y0:-Math.log10(0.05),y1:-Math.log10(0.05),line:{dash:'dash',color:'#93a1b3'}}]},
    {responsive:true,displaylogo:false});
}
function table(id, searchId, rows){
  const cols=['Gene','log2 fold change','Average expression','Adjusted p-value'];
  let sortCol=0,asc=true,cur=rows;
  function render(){
    const t=document.getElementById(id);
    t.innerHTML='<tr>'+cols.map((c,i)=>`<th data-i="${i}">${c}${sortCol===i?(asc?' ▲':' ▼'):''}</th>`).join('')+'</tr>'+
      cur.map(r=>'<tr>'+r.map((c,i)=>`<td>${i===0?'<i>'+c+'</i>':c}</td>`).join('')+'</tr>').join('');
    t.querySelectorAll('th').forEach(th=>th.onclick=()=>{const i=+th.dataset.i;
      if(sortCol===i)asc=!asc;else{sortCol=i;asc=true;}
      cur=[...cur].sort((a,b)=>{const x=a[i],y=b[i];const nx=parseFloat(x),ny=parseFloat(y);
        const c=(isNaN(nx)||isNaN(ny))?String(x).localeCompare(y):nx-ny;return asc?c:-c;});
      render();});
  }
  document.getElementById(searchId).addEventListener('input',e=>{
    const q=e.target.value.toLowerCase();
    cur=rows.filter(r=>r[0].toLowerCase().includes(q));render();});
  render();
}
document.querySelectorAll('.tab').forEach(t=>t.onclick=()=>{
  document.querySelectorAll('.tab').forEach(x=>x.classList.remove('active'));
  document.querySelectorAll('.panel').forEach(x=>x.classList.remove('active'));
  t.classList.add('active');document.getElementById('p-'+t.dataset.t).classList.add('active');
  window.dispatchEvent(new Event('resize'));});
cards('m-cards',DATA.mouse.metrics);cards('a-cards',DATA.airway.metrics);
volcano('m-volcano',DATA.mouse.volcano,'log2 fold change (lactating vs pregnant)');
volcano('a-volcano',DATA.airway.volcano,'log2 fold change (dexamethasone vs untreated)');
document.getElementById('m-pca').src=IMGS.mouse_pca;
document.getElementById('m-heat').src=IMGS.mouse_heatmap;
document.getElementById('m-enrich').src=IMGS.mouse_enrich;
document.getElementById('a-pca').src=IMGS.airway_pca;
table('m-table','m-search',DATA.mouse.table);
table('a-table','a-search',DATA.airway.table);
</script></body></html>"""

if __name__ == "__main__":
    main()
