# IDER — NeurIPS 2026 Project Page

Project page for **[Preserving DEG Rankings for Gene Discovery in Histology-Based Spatial Gene Expression Prediction](https://arxiv.org/abs/2609.33928)**, accepted at **NeurIPS 2026**.

This is a lightweight, static GitHub Pages website. It includes the research concept figure, abstract, paper and code links, citation, and a direct link to the **[Metric Tutorial notebook](https://github.com/Shiku-Kaito/IDER/blob/main/tutorial/01_tutorial_metrics.ipynb)** in the code repository.

The project page does **not** contain a browser-based interactive demo. The original Python metric implementation remains in the [IDER repository](https://github.com/Shiku-Kaito/IDER/blob/main/script/culc_DEG_metrics.py).

## Site structure

```text
index.html                      Project page (overview, abstract, citation)
assets/css/style.css            Layout and subtle section backgrounds
assets/js/main.js               Concept figure fallback and BibTeX copy
assets/images/ider_concept.png  Research concept figure
LICENSE
```

## Preview locally

Open `index.html` directly in a browser, or run:

```bash
python3 -m http.server 8000
```

Visit `http://localhost:8000/`. No build process or dependencies are required.

## Publish with GitHub Pages

Copy the **contents** of `IDER-project-page/` (not the surrounding folder) into the root of the project-page repository. In the repository's **Settings → Pages**, choose **Deploy from a branch**, `main`, `/ (root)`.

Expected URL once enabled: `https://shiku-kaito.github.io/IDER-project-page/` (adjust for your actual repository name).

## Notes

Check that all figures and assets can be redistributed. This website does not include the HER2ST dataset or CONCH weights; those resources have separate licensing and access requirements.
