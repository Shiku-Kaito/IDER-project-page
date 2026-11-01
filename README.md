# IDER — NeurIPS 2026 Project Page

Static HTML/CSS/JS project page for *Preserving DEG Rankings for Gene Discovery in Histology-Based Spatial Gene Expression Prediction*.

## Add the research concept figure

Save your figure as:

`assets/images/ider_concept.png`

PNG is recommended. The page automatically displays the figure when the file exists, and shows a placeholder when it is absent.

## Preview locally (without publishing GitHub Pages)

```bash
cd IDER-project-page
python3 -m http.server 8000
```

Visit http://localhost:8000/.

## Publish later

Upload the contents of this directory to the **root** of your private GitHub repository `IDER`. When ready to publish, enable GitHub Pages from Settings → Pages (main branch, /(root)); note that GitHub Pages websites are public even if the repository is private on typical personal plans.

## Before public release

- The BibTeX presently retains the `Anonymous` author field and citation key as provided. Replace with the final author list if appropriate for camera-ready.
- Add a real GitHub code link when the code repository is public.
- Check that paper figures and other assets can be redistributed.
