/**
 * CSV Export utility.
 * Exports candidate results to a downloadable CSV string.
 * Hard rule: The "Known Indication" column is excluded. The word "indication" must never appear.
 */

export function generateCandidateCSV({ mode, query, queryEnvelope: _queryEnvelope, candidates }) {
  const timestamp = new Date().toISOString();
  const dateFormatted = timestamp.replace(/:/g, '-').slice(0, 19);

  const lines = [];

  // Metadata comments at top
  lines.push(`# MoleculAR Research Tool - Candidate Export`);
  lines.push(`# Mode: ${mode === 'drug' ? 'Drug Search' : 'Target Search'}`);
  lines.push(`# Query: ${query}`);
  lines.push(`# Timestamp: ${timestamp}`);
  lines.push(`# Note: Rankings show similarity, not proven treatments. Not medical advice.`);
  lines.push('');

  // CSV Headers - NO indication column
  const headers = [
    'Rank',
    'Candidate',
    'Score',
    'Structure_Similarity',
    'Target_Similarity',
    'Side_Effect_Similarity',
    'Pre_Penalty_Sum',
    'Hubness_Percentile',
    'Has_Target_Data',
    'Has_Side_Effect_Data',
    'Salt_Or_Stereoisomer_Form',
    'Shared_Targets',
    'Shared_Side_Effects',
  ];
  lines.push(headers.join(','));

  // Rows
  candidates.forEach((cand, idx) => {
    const s = cand.struct_sim ?? 0;
    const t = cand.target_sim ?? 0;
    const e = cand.side_effect_sim ?? 0;
    const prePenalty = 0.4 * s + 0.4 * t + 0.2 * e;

    const row = [
      idx + 1,
      `"${(cand.candidate || '').replace(/"/g, '""')}"`,
      Number(cand.score).toFixed(4),
      Number(cand.struct_sim).toFixed(4),
      Number(cand.target_sim).toFixed(4),
      Number(cand.side_effect_sim).toFixed(4),
      prePenalty.toFixed(4),
      cand.hubness_percentile !== null && cand.hubness_percentile !== undefined
        ? Number(cand.hubness_percentile).toFixed(2)
        : '',
      cand.has_targets ? 'Yes' : 'No',
      cand.has_side_effects ? 'Yes' : 'No',
      cand.is_duplicate_ingredient ? 'Yes' : 'No',
      `"${(cand.shared_targets || []).join('; ').replace(/"/g, '""')}"`,
      `"${(cand.shared_side_effects || []).join('; ').replace(/"/g, '""')}"`,
    ];
    lines.push(row.join(','));
  });

  return {
    filename: `molecular_export_${mode}_${query}_${dateFormatted}.csv`,
    content: lines.join('\n'),
  };
}

export function downloadCSV(filename, csvContent) {
  const blob = new Blob([csvContent], { type: 'text/csv;charset=utf-8;' });
  const url = URL.createObjectURL(blob);
  const link = document.createElement('a');
  link.setAttribute('href', url);
  link.setAttribute('download', filename);
  document.body.appendChild(link);
  link.click();
  document.body.removeChild(link);
  URL.revokeObjectURL(url);
}
