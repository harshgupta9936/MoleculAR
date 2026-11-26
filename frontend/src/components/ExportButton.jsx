import { useState } from 'react';
import { generateCandidateCSV, downloadCSV } from '../utils/csvExport';

export function ExportButton({ mode, query, queryEnvelope, candidates }) {
  const [downloading, setDownloading] = useState(false);

  if (!candidates || candidates.length === 0) {
    return null;
  }

  const handleExport = () => {
    setDownloading(true);
    try {
      const { filename, content } = generateCandidateCSV({
        mode,
        query,
        queryEnvelope,
        candidates,
      });
      downloadCSV(filename, content);
    } catch (err) {
      console.error('Export failed:', err);
    } finally {
      setTimeout(() => setDownloading(false), 400);
    }
  };

  return (
    <button
      className="btn-export"
      onClick={handleExport}
      disabled={downloading}
      title="Download candidates as CSV (excludes indication data)"
      aria-label="Export results as CSV"
    >
      <svg
        className="btn-icon"
        viewBox="0 0 16 16"
        width="14"
        height="14"
        fill="currentColor"
        aria-hidden="true"
      >
        <path d="M.5 9.9a.5.5 0 0 1 .5.5v2.5a1 1 0 0 0 1 1h12a1 1 0 0 0 1-1v-2.5a.5.5 0 0 1 1 0v2.5a2 2 0 0 1-2 2H2a2 2 0 0 1-2-2v-2.5a.5.5 0 0 1 .5-.5z" />
        <path d="M7.646 11.854a.5.5 0 0 0 .708 0l3-3a.5.5 0 0 0-.708-.708L8.5 10.293V1.5a.5.5 0 0 0-1 0v8.793L5.354 8.146a.5.5 0 1 0-.708.708l3 3z" />
      </svg>
      {downloading ? 'Exporting...' : 'Export CSV'}
    </button>
  );
}
