/*
 * IDER interactive metric demo. Pure JavaScript, no server or dependencies.
 * Port of: https://github.com/Shiku-Kaito/IDER/blob/main/script/culc_DEG_metrics.py
 *
 * DEG_ranking: one-vs-rest Mann-Whitney U via midranks for tied expression.
 * DEG_metric_func: cluster-average Spearman on U-score vectors and binary
 * nDCG@k, using mean gain for ties in predicted scores (as in sklearn).
 */
(function (root) {
  'use strict';

  function sortedIndicesDescending(values) {
    return values.map((_, i) => i).sort((a, b) => values[b] - values[a] || a - b);
  }

  // Average ranks, 1-based, for equal values (used for Spearman and U).
  function averageRanks(values) {
    const indices = values.map((_, i) => i).sort((a, b) => values[a] - values[b] || a - b);
    const ranks = new Array(values.length);
    for (let i = 0; i < indices.length;) {
      let j = i + 1;
      while (j < indices.length && values[indices[j]] === values[indices[i]]) j++;
      const rank = (i + 1 + j) / 2;
      for (let p = i; p < j; p++) ranks[indices[p]] = rank;
      i = j;
    }
    return ranks;
  }

  function pearson(x, y) {
    const n = x.length;
    const mx = x.reduce((a, b) => a + b, 0) / n;
    const my = y.reduce((a, b) => a + b, 0) / n;
    let num = 0, vx = 0, vy = 0;
    for (let i = 0; i < n; i++) {
      const dx = x[i] - mx, dy = y[i] - my;
      num += dx * dy; vx += dx * dx; vy += dy * dy;
    }
    if (vx === 0 || vy === 0) return NaN;
    return num / Math.sqrt(vx * vy);
  }

  function spearman(x, y) {
    return pearson(averageRanks(x), averageRanks(y));
  }

  // U is the count of (a>b) pairs plus 0.5 per tie, implemented as rank-sum.
  function mannWhitneyU(matrix, clusters) {
    const n = matrix.length, genes = matrix[0].length;
    const groupIds = [...new Set(clusters)].sort((a, b) => String(a).localeCompare(String(b), undefined, {numeric: true}));
    const output = {};
    for (const cid of groupIds) {
      const mask = clusters.map(c => c === cid);
      const na = mask.filter(Boolean).length;
      if (na === 0 || na === n) continue;
      const u = new Array(genes);
      for (let g = 0; g < genes; g++) {
        const values = matrix.map(row => row[g]);
        const ranks = averageRanks(values);
        let rankSum = 0;
        for (let i = 0; i < n; i++) if (mask[i]) rankSum += ranks[i];
        u[g] = rankSum - na * (na + 1) / 2;
      }
      output[cid] = u;
    }
    return output;
  }

  // Binary relevance: genes in the ground-truth top-k have gain=1.
  // sklearn.ndcg_score's default behavior averages relevance within tied
  // predicted score groups, not arbitrarily breaking a predicted tie.
  function binaryNdcg(pred, gt, k) {
    const n = gt.length;
    if (!Number.isInteger(k) || k < 1 || k > n) throw new Error(`k must be between 1 and ${n}.`);
    const relevant = new Array(n).fill(0);
    const groundOrder = sortedIndicesDescending(gt);
    for (let i = 0; i < k; i++) relevant[groundOrder[i]] = 1;
    const predictedOrder = sortedIndicesDescending(pred);
    let dcg = 0;
    for (let i = 0; i < n && i < k;) {
      let j = i + 1;
      while (j < n && pred[predictedOrder[j]] === pred[predictedOrder[i]]) j++;
      let totalGain = 0;
      for (let p = i; p < j; p++) totalGain += relevant[predictedOrder[p]];
      const meanGain = totalGain / (j - i);
      for (let p = i; p < j && p < k; p++) dcg += meanGain / Math.log2(p + 2);
      i = j;
    }
    let idealDcg = 0;
    for (let i = 0; i < k; i++) idealDcg += 1 / Math.log2(i + 2);
    return dcg / idealDcg;
  }

  function validateData(gt, pred, clusters) {
    if (!Array.isArray(gt) || !gt.length || !Array.isArray(pred) || !Array.isArray(clusters)) {
      throw new Error('Enter both expression matrices and cluster IDs.');
    }
    if (gt.length !== pred.length || gt.length !== clusters.length) {
      throw new Error('Measured, predicted, and cluster inputs must contain the same number of spots (rows).');
    }
    if (gt.length > 400) throw new Error('This browser demo supports up to 400 spots.');
    const g = gt[0].length;
    if (g < 2 || g > 300 || gt.length * g > 50000) throw new Error('Use between 2 and 300 genes, with at most 50,000 values in each matrix.');
    for (const matrix of [gt, pred]) {
      for (const row of matrix) {
        if (!Array.isArray(row) || row.length !== g || row.some(v => !Number.isFinite(v))) {
          throw new Error('Every expression row must have the same number of genes and contain only finite numbers.');
        }
      }
    }
    if (new Set(clusters).size < 2) throw new Error('At least two different cluster IDs are required.');
    if (clusters.some(id => !String(id).trim())) throw new Error('Cluster IDs cannot be empty.');
  }

  function evaluate(gt, pred, clusters, k) {
    validateData(gt, pred, clusters);
    if (!Number.isInteger(k) || k < 1 || k > gt[0].length) throw new Error('Top-k must be within the gene count.');
    const measuredU = mannWhitneyU(gt, clusters);
    const predictedU = mannWhitneyU(pred, clusters);
    const ids = Object.keys(measuredU);
    let totalSpearman = 0, totalNdcg = 0, totalOverlap = 0;
    for (const id of ids) {
      const s = spearman(predictedU[id], measuredU[id]);
      if (!Number.isFinite(s)) throw new Error('Spearman is undefined because all U-scores are equal for at least one cluster. Try more genes or varied data.');
      totalSpearman += s;
      totalNdcg += binaryNdcg(predictedU[id], measuredU[id], k);
      const gtTop = new Set(sortedIndicesDescending(measuredU[id]).slice(0, k));
      totalOverlap += sortedIndicesDescending(predictedU[id]).slice(0, k).filter(g => gtTop.has(g)).length / k;
    }
    return {
      measuredU, predictedU, ids,
      spearman: totalSpearman / ids.length,
      ndcg: totalNdcg / ids.length,
      overlap: totalOverlap / ids.length,
      spots: gt.length,
      genes: gt[0].length,
      clusters: ids.length,
      k
    };
  }

  function parseMatrix(text, label) {
    const lines = text.trim().split(/\r?\n/).map(l => l.trim()).filter(Boolean);
    if (!lines.length) throw new Error(`Please enter ${label}.`);
    return lines.map((line, idx) => {
      const cells = line.split(/[,;\s]+/).filter(Boolean);
      if (!cells.length || cells.some(x => !x.trim() || !Number.isFinite(Number(x)))) {
        throw new Error(`Invalid numeric value in ${label}, row ${idx + 1}.`);
      }
      return cells.map(Number);
    });
  }

  function parseClusterIds(text) {
    return text.trim().split(/[,;\s]+/).filter(Boolean);
  }

  // Seeded uniform and normal RNG to make example comparisons reproducible.
  function randomNormal(seed) {
    let state = (seed >>> 0) || 0x6d2b79f5;
    function uniform() {
      state ^= state << 13; state ^= state >>> 17; state ^= state << 5;
      return ((state >>> 0) + 0.5) / 4294967296;
    }
    let spare = null;
    return () => {
      if (spare !== null) { const val = spare; spare = null; return val; }
      const radius = Math.sqrt(-2 * Math.log(uniform()));
      const theta = 2 * Math.PI * uniform();
      spare = radius * Math.sin(theta);
      return radius * Math.cos(theta);
    };
  }

  function makeToyData(noise, seed) {
    const normal = randomNormal(seed), n = 90, genes = 30;
    const clusters = Array.from({length: n}, (_, i) => String(i % 3));
    const gt = [], pred = [];
    for (let i = 0; i < n; i++) {
      const gRow = [], pRow = [];
      for (let g = 0; g < genes; g++) {
        const base = normal() + (g >= Number(clusters[i]) * 5 && g < (Number(clusters[i]) + 1) * 5 ? 2 : 0);
        gRow.push(base); pRow.push(base + noise * normal());
      }
      gt.push(gRow); pred.push(pRow);
    }
    return {gt, pred, clusters};
  }

  // Public, testable numeric implementation.
  const API = {averageRanks, mannWhitneyU, spearman, binaryNdcg, evaluate, makeToyData, parseMatrix};
  if (typeof module !== 'undefined' && module.exports) module.exports = API;
  root.IDERMetrics = API;

  if (typeof document === 'undefined') return;
  const byId = id => document.getElementById(id);
  if (!byId('demo')) return;

  const tabToy = byId('tab-toy'), tabCustom = byId('tab-custom');
  const panelToy = byId('panel-toy'), panelCustom = byId('panel-custom');
  const noiseSlider = byId('noise-slider'), seedInput = byId('seed-input');
  const kSelect = byId('topk-select'), clusterSelect = byId('cluster-select');
  const errorBox = byId('demo-error');
  const customInitial = ['custom-clusters', 'custom-gt', 'custom-pred'].map(id => byId(id).value);
  let mode = 'toy', data = null, result = null;

  function showError(message) { errorBox.textContent = message; errorBox.hidden = !message; }
  function selectMode(newMode) {
    mode = newMode;
    const isToy = mode === 'toy';
    tabToy.classList.toggle('active', isToy);
    tabCustom.classList.toggle('active', !isToy);
    tabToy.setAttribute('aria-selected', String(isToy));
    tabCustom.setAttribute('aria-selected', String(!isToy));
    panelToy.hidden = !isToy;
    panelCustom.hidden = isToy;
    kSelect.value = isToy ? '10' : '2';
    showError('');
    run();
  }

  function updateKOptions(count) {
    const options = count >= 20 ? [5, 10, 20] : count >= 10 ? [3, 5, 10] : count >= 5 ? [1, 2, 3] : Array.from({length: count}, (_, i) => i + 1);
    const selected = Number(kSelect.value);
    kSelect.replaceChildren(...options.map(v => {
      const item = document.createElement('option'); item.value = String(v); item.textContent = String(v); return item;
    }));
    kSelect.value = options.includes(selected) ? String(selected) : String(options[Math.min(1, options.length - 1)]);
  }

  function readData() {
    if (mode === 'toy') {
      const noise = Number(noiseSlider.value), seed = Number(seedInput.value);
      if (!Number.isInteger(seed) || seed < 0 || seed > 999999) throw new Error('Random seed must be an integer from 0 to 999999.');
      return makeToyData(noise, seed);
    }
    return {
      gt: parseMatrix(byId('custom-gt').value, 'measured expression'),
      pred: parseMatrix(byId('custom-pred').value, 'predicted expression'),
      clusters: parseClusterIds(byId('custom-clusters').value)
    };
  }

  function run() {
    showError('');
    try {
      const newData = readData();
      validateData(newData.gt, newData.pred, newData.clusters);
      const previousGenes = data?.gt[0].length;
      if (previousGenes !== newData.gt[0].length || (mode === 'custom' && Number(kSelect.value) > newData.gt[0].length)) {
        updateKOptions(newData.gt[0].length);
      }
      data = newData;
      result = evaluate(data.gt, data.pred, data.clusters, Number(kSelect.value));
      render();
    } catch (err) {
      showError(err.message);
      byId('demo-results').hidden = true;
    }
  }

  function render() {
    byId('demo-results').hidden = false;
    const r = result;
    const format = x => x.toFixed(3);
    byId('demo-summary').textContent = `${mode === 'toy' ? 'Synthetic example' : 'Your data'} · ${r.spots} spots · ${r.genes} genes · ${r.clusters} clusters`;
    byId('metric-spearman').textContent = format(r.spearman);
    byId('metric-ndcg').textContent = format(r.ndcg);
    byId('metric-overlap').textContent = `${Math.round(r.overlap * 100)}%`;
    byId('metric-ndcg-name').textContent = `Binary nDCG@${r.k}`;
    byId('metric-overlap-name').textContent = `Top-${r.k} overlap`;
    byId('overlap-detail').textContent = 'Average fraction shared';
    const previousCluster = clusterSelect.value;
    clusterSelect.replaceChildren(...r.ids.map(cid => {
      const option = document.createElement('option'); option.value = cid; option.textContent = cid; return option;
    }));
    if (r.ids.includes(previousCluster)) clusterSelect.value = previousCluster;
    renderRankings();
  }

  function renderRankings() {
    if (!result) return;
    const cid = clusterSelect.value;
    const r = result, n = Math.min(7, r.genes);
    const measuredTop = sortedIndicesDescending(r.measuredU[cid]).slice(0, n);
    const predictedTop = sortedIndicesDescending(r.predictedU[cid]).slice(0, n);
    const gtSet = new Set(sortedIndicesDescending(r.measuredU[cid]).slice(0, r.k));
    const containers = [byId('measured-genes'), byId('predicted-genes')];
    const rows = [[measuredTop, r.measuredU[cid]], [predictedTop, r.predictedU[cid]]];
    for (let col = 0; col < 2; col++) {
      const [top, scores] = rows[col];
      const maxU = Math.max(...scores, 1);
      containers[col].replaceChildren(...top.map((gene, index) => {
        const item = document.createElement('div'); item.className = 'gene-row';
        const rank = document.createElement('span'); rank.className = 'gene-rank'; rank.textContent = String(index + 1).padStart(2, '0');
        const geneLabel = document.createElement('span'); geneLabel.className = 'gene-label';
        geneLabel.textContent = `Gene ${String(gene + 1).padStart(2, '0')}`;
        if (col === 1 && gtSet.has(gene)) { geneLabel.classList.add('matched'); geneLabel.title = `In the measured top ${r.k}`; }
        const bar = document.createElement('div'); bar.className = 'gene-track';
        const fill = document.createElement('div'); fill.className = `gene-fill ${col === 0 ? 'measured' : 'predicted'}`;
        fill.style.width = `${Math.max(0, scores[gene] / maxU * 100)}%`; bar.append(fill);
        const score = document.createElement('span'); score.className = 'gene-score'; score.textContent = scores[gene].toFixed(1);
        item.append(rank, geneLabel, bar, score); return item;
      }));
    }
  }

  tabToy.addEventListener('click', () => selectMode('toy'));
  tabCustom.addEventListener('click', () => selectMode('custom'));
  noiseSlider.addEventListener('input', () => { byId('noise-label').textContent = Number(noiseSlider.value).toFixed(2); });
  noiseSlider.addEventListener('change', run);
  byId('toy-run').addEventListener('click', run);
  seedInput.addEventListener('change', run);
  byId('custom-run').addEventListener('click', run);
  byId('custom-reset').addEventListener('click', () => {
    ['custom-clusters', 'custom-gt', 'custom-pred'].forEach((id, i) => { byId(id).value = customInitial[i]; });
    run();
  });
  kSelect.addEventListener('change', () => {
    if (!data) return;
    try { result = evaluate(data.gt, data.pred, data.clusters, Number(kSelect.value)); showError(''); render(); }
    catch (err) { showError(err.message); }
  });
  clusterSelect.addEventListener('change', renderRankings);

  // Populate the demonstration on first load, so users see real outputs immediately.
  run();
})(typeof window !== 'undefined' ? window : globalThis);
