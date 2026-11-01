const conceptImage = document.getElementById('concept-image');
const conceptPlaceholder = document.getElementById('concept-placeholder');
conceptImage.addEventListener('load', () => { conceptImage.hidden = false; conceptPlaceholder.hidden = true; });
conceptImage.addEventListener('error', () => { conceptImage.hidden = true; conceptPlaceholder.hidden = false; });
if (conceptImage.complete && conceptImage.naturalWidth > 0) { conceptImage.hidden = false; conceptPlaceholder.hidden = true; }
document.getElementById('copy-bib').addEventListener('click', async (event) => {
  const text = document.getElementById('bibtex').textContent;
  try { await navigator.clipboard.writeText(text); event.currentTarget.textContent = 'Copied!'; }
  catch { event.currentTarget.textContent = 'Select text to copy'; }
  setTimeout(() => { event.currentTarget.textContent = 'Copy'; }, 1800);
});
