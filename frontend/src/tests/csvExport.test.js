import { describe, it, expect } from 'vitest';
import { generateCandidateCSV } from '../utils/csvExport';

describe('csvExport', () => {
  const sampleCandidates = [
    {
      candidate: 'VARDENAFIL',
      score: 0.4888,
      struct_sim: 0.82,
      target_sim: 0.95,
      side_effect_sim: 0.41,
      hubness_percentile: 57.48,
      has_targets: true,
      has_side_effects: true,
      shared_targets: ['PDE5A'],
      shared_side_effects: ['headache', 'flushing'],
      indication: 'Erectile dysfunction', // Present in DB/API data
    },
    {
      candidate: 'TADALAFIL',
      score: 0.4651,
      struct_sim: 0.71,
      target_sim: 0.91,
      side_effect_sim: 0.38,
      hubness_percentile: 33.17,
      has_targets: true,
      has_side_effects: true,
      shared_targets: ['PDE5A'],
      shared_side_effects: ['headache'],
      indication: 'Pulmonary hypertension',
    },
  ];

  it('generates valid CSV with required columns', () => {
    const { filename, content } = generateCandidateCSV({
      mode: 'drug',
      query: 'SILDENAFIL',
      candidates: sampleCandidates,
    });

    expect(filename).toContain('molecular_export_drug_SILDENAFIL_');
    expect(content).toContain('Rank,Candidate,Score');
    expect(content).toContain('Salt_Or_Stereoisomer_Form');
    expect(content).toContain('"VARDENAFIL",0.4888');
    expect(content).toContain('"TADALAFIL",0.4651');
  });

  it('MUST NOT contain the word or column "indication" anywhere in output', () => {
    const { content } = generateCandidateCSV({
      mode: 'drug',
      query: 'SILDENAFIL',
      candidates: sampleCandidates,
    });

    // Check case-insensitive occurrence of "indication"
    const hasIndication = /indication/i.test(content);
    expect(hasIndication).toBe(false);

    // Also verify actual indication values are stripped
    expect(content).not.toContain('Erectile dysfunction');
    expect(content).not.toContain('Pulmonary hypertension');
  });
});
