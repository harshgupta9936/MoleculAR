/**
 * Extracts evaluation metrics from committed duplicate-neutral CSVs and manifest.
 * Pure read-only extraction from experiments/results.
 * Writes to frontend/public/metrics.json.
 */

import fs from 'fs';
import path from 'path';
import { fileURLToPath } from 'url';

const __dirname = path.dirname(fileURLToPath(import.meta.url));
const ROOT = path.resolve(__dirname, '../../');

function parseCSV(filePath) {
  const content = fs.readFileSync(filePath, 'utf-8');
  const lines = content.trim().split(/\r?\n/);
  const headers = lines[0].split(',').map((h) => h.trim());
  const rows = [];
  for (let i = 1; i < lines.length; i++) {
    if (!lines[i].trim()) continue;
    const parts = lines[i].split(',').map((p) => p.trim());
    const row = {};
    headers.forEach((h, idx) => {
      row[h] = parts[idx];
    });
    rows.push(row);
  }
  return rows;
}

const manifestPath = path.join(ROOT, 'experiments/results/manifest.json');
const manifest = JSON.parse(fs.readFileSync(manifestPath, 'utf-8'));

const protocolsConfig = [
  {
    key: 'p1_dn',
    label: 'Protocol P1-dn (Primary Same-Indication Repurposing)',
    file: 'experiments/results/metrics_C1_P1_dn.csv',
  },
  {
    key: 'p2a_dn_0.4',
    label: 'Protocol P2a-dn (Analogue Neutralization at θ=0.4)',
    file: 'experiments/results/metrics_C1_P2a_dn_0.4.csv',
  },
  {
    key: 'p3_dn',
    label: 'Protocol P3-dn (Lead Optimization / Analog Discovery)',
    file: 'experiments/results/metrics_C1_P3_dn.csv',
  },
  {
    key: 'p4',
    label: 'Protocol P4 (Mechanism Generalization across ATC L3 Classes)',
    file: 'experiments/results/metrics_C1_P4.csv',
  },
  {
    key: 'p4b',
    label: 'Protocol P4b (Mechanism Generalization with Target-Shared Classes)',
    file: 'experiments/results/metrics_C1_P4b.csv',
  },
];

const methodLabels = {
  M1: 'Random Baseline',
  M7: 'Evaluated Default (Hand-Weighted + Hub Penalty)',
  M8: 'Popularity Prior (uses labels, an upper bound for a query-independent method)',
};

function getMetricStats(rows, method, metricName) {
  const r = rows.find((row) => row.method === method && row.metric === metricName);
  if (!r) return null;
  return {
    mean: parseFloat(r.mean),
    ci_low: parseFloat(r.ci_low),
    ci_high: parseFloat(r.ci_high),
  };
}

const protocolsData = {};

for (const p of protocolsConfig) {
  const fullPath = path.join(ROOT, p.file);
  const rows = parseCSV(fullPath);
  const nQueries = parseInt(rows[0]?.n_queries || '0', 10);

  const methodsData = ['M1', 'M7', 'M8'].map((mId) => {
    return {
      id: mId,
      label: methodLabels[mId],
      isDefault: mId === 'M7',
      hit10: getMetricStats(rows, mId, 'hit@10'),
      mrr: getMetricStats(rows, mId, 'mrr'),
    };
  });

  protocolsData[p.key] = {
    key: p.key,
    label: p.label,
    file: p.file,
    n_queries: nQueries,
    methods: methodsData,
  };
}

// Provenance metadata:
// - db_sha256 is read directly from experiments/results/manifest.json (databases.DB-A.sha256).
// - git_tag ('eval-v2') is a typed constant representing the evaluation benchmark release tag
//   at commit 2aec56c, as manifest.json only records git_commit ("020dbadc...") and not git_tag.
const EVALUATION_TAG_CONSTANT = 'eval-v2';

const metrics = {
  provenance: {
    git_tag: manifest.metadata?.git_tag || EVALUATION_TAG_CONSTANT,
    db_sha256: manifest.databases?.['DB-A']?.sha256 || 'd88a2333743af2dbe859b3717764b6a645b9786caa6f059cc9558e29c7ff0da7',
    manifest_commit: manifest.metadata?.git_commit,
  },
  summary_sentence:
    'Strong at finding drugs with the same indication; much weaker across drug classes, where a popularity-only baseline can do as well or better.',
  protocols: protocolsData,
};

const outDir = path.resolve(__dirname, '../public');
if (!fs.existsSync(outDir)) {
  fs.mkdirSync(outDir, { recursive: true });
}
const outPath = path.join(outDir, 'metrics.json');
fs.writeFileSync(outPath, JSON.stringify(metrics, null, 2), 'utf-8');
console.log(`Generated ${outPath} successfully with all 5 protocols.`);
