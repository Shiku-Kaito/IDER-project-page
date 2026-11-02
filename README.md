# IDER — NeurIPS 2026 Project Page

Project website for **[Preserving DEG Rankings for Gene Discovery in Histology-Based Spatial Gene Expression Prediction](https://arxiv.org/abs/2609.33928)** (NeurIPS 2026).

The original landing page has been extended with an **interactive DEG metric tutorial** based on the evaluation code in [Shiku-Kaito/IDER](https://github.com/Shiku-Kaito/IDER):

- One-vs-rest Mann–Whitney U-scores per gene and cluster (average ranks for ties)
- Spearman correlation between the predicted and measured U-score vectors, averaged across clusters
- Binary-relevance nDCG@k, averaged across clusters, including sklearn-style averaging for ties in predicted U-scores
- A simple top-k overlap visualization for interpretation (additional to the paper's metrics)

The demo runs entirely in the visitor's browser using JavaScript and **does not need any backend, GPU, PyTorch, or data download**. It accepts toy data with adjustable noise and seed, or user-pasted expression matrices and cluster IDs. Custom data is never uploaded.

For original Python code, visit the [metric notebook](https://github.com/Shiku-Kaito/IDER/blob/main/tutorial/01_tutorial_metrics.ipynb) or [`culc_DEG_metrics.py`](https://github.com/Shiku-Kaito/IDER/blob/main/script/culc_DEG_metrics.py).

## Files

```
index.html                    Landing page, interactive demo, BibTeX
assets/css/style.css          Original styling + demo components
assets/js/main.js             Figure display and BibTeX copy
assets/js/metric-demo.js      Standalone, dependency-free metric calculations + demo UI
assets/images/ider_concept.png
LICENSE
```

## Preview locally

```bash
python3 -m http.server 8000
```

Open `http://localhost:8000/`. No build step is required.

## Publish on GitHub Pages

1. Copy the **contents** of this directory to the root of the `Shiku-Kaito/IDER-project-page` repository and push.
2. In **Settings → Pages**, set **Deploy from a branch**, branch `main`, folder `/ (root)` (if your plan and repository visibility permit it).
3. The expected project-page URL is `https://shiku-kaito.github.io/IDER-project-page/` (once Pages is enabled and deployment succeeds).

Note: GitHub Pages visibility and whether Pages can publish from a private repository depend on the account plan. Do not assume a private source repository makes the website private.

## Metric details

This demo ports the numeric definitions of `DEG_ranking` and `DEG_metric_func` into browser-native JavaScript, rather than running those PyTorch and scikit-learn modules in the browser. The evaluated quantities are **not** DEG p-values. The U-statistic uses midranks for expression ties; Spearman uses average ranks; nDCG@k uses binary relevance for measured top-k genes, and mean gains for tied predicted scores. In cases of *tied measured U-scores at the top-k cutoff*, the selected relevant genes may depend on the sorting algorithm's tie order; the demo deterministically breaks such ties by gene index. Avoid treating top-k membership under exact ties as biologically meaningful.

The demo is educational and uses synthetic data by default. Actual paper results require the full repository and data pipeline.

## Validate the metric implementation

If you have Node.js, run:

```bash
node tests/check-metrics.js
```

## Attribution and caveats

Paper: Kaito Shiku, Kazuya Nishimura, Yasuhiro Kojima, Ryoma Bise (2026).

Check that all figures and other assets may be redistributed before public release. This project does not bundle the HER2ST dataset, CONCH model weights, or their licenses.
