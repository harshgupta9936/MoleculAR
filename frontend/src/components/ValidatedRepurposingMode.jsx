import { useState } from 'react';

export default function RepurposingMode() {
  const [query, setQuery] = useState('');
  const [results, setResults] = useState([]);
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState(null);

  const handleSearch = async (e) => {
    e.preventDefault();
    if (!query.trim()) return;

    setLoading(true);
    setError(null);
    setResults([]);

    try {
      const response = await fetch('http://localhost:8000/api/rank', {
        method: 'POST',
        headers: {
          'Content-Type': 'application/json',
        },
        body: JSON.stringify({ drug_name: query, top_n: 10 }),
      });

      if (!response.ok) {
        if (response.status === 404) {
          throw new Error('Drug not found in dataset. Please check spelling.');
        }
        throw new Error('Server error occurred.');
      }

      const data = await response.json();
      setResults(data.results);
    } catch (err) {
      if (err.name === 'TypeError' || err.message === 'Failed to fetch') {
        setError('Connection Error: Unable to reach the backend ranking engine.');
      } else {
        setError(err.message);
      }
    } finally {
      setLoading(false);
    }
  };

  return (
    <div className="mode-container repurposing-mode">
      <div className="badge-container">
        <span className="badge repurposing-badge">✓ Phase 4.4 Hand-Weighted</span>
      </div>
      
      <p className="description">
        Enter an approved drug to find structurally and mechanistically similar candidates for repurposing.
      </p>

      <form onSubmit={handleSearch} className="search-form">
        <input 
          type="text" 
          value={query} 
          onChange={(e) => setQuery(e.target.value)} 
          placeholder="e.g., SILDENAFIL"
          className="search-input"
        />
        <button type="submit" disabled={loading} className="search-btn">
          {loading ? 'Ranking...' : 'Search Candidates'}
        </button>
      </form>

      {error && <div className="error-banner">{error}</div>}

      <div className="results-list">
        {results.map((result, idx) => (
          <div key={idx} className="candidate-card">
            <div className="card-header">
              <h2>#{idx + 1} {result.candidate}</h2>
              <span className="score">Score: {result.score.toFixed(4)}</span>
            </div>
            <p className="indication">Current Indication: {result.indication}</p>
            
            <div className="features">
              {result.shared_targets.length > 0 && (
                <div className="feature-group">
                  <strong>Shared Targets:</strong>
                  <div className="tags">
                    {result.shared_targets.map(t => <span key={t} className="tag target-tag">{t}</span>)}
                  </div>
                </div>
              )}

              {result.shared_side_effects.length > 0 ? (
                <div className="feature-group">
                  <strong>Shared Specific Side Effects:</strong>
                  <div className="tags">
                    {result.shared_side_effects.map(se => <span key={se} className="tag se-tag">{se}</span>)}
                  </div>
                </div>
              ) : (
                <div className="feature-group">
                  <em>Driven primarily by structural and phenotypic similarity</em>
                </div>
              )}
            </div>
          </div>
        ))}
      </div>
    </div>
  );
}
