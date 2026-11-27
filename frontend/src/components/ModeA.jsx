import { useState } from 'react'
import { API_BASE_URL } from '../utils/api'

function ModeA() {
  const [query, setQuery] = useState('')
  const [topN, setTopN] = useState(10)
  const [results, setResults] = useState(null)
  const [error, setError] = useState(null)
  const [isLoading, setIsLoading] = useState(false)

  const handleSearch = async (e) => {
    e.preventDefault()
    if (!query.trim()) return

    setIsLoading(true)
    setError(null)
    setResults(null)

    try {
      const response = await fetch(`${API_BASE_URL}/api/rank`, {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({ drug_name: query.trim(), top_n: Number(topN) })
      })

      if (!response.ok) {
        if (response.status === 404) {
          throw new Error(`Drug '${query.toUpperCase()}' not found in the database. Please check spelling or refer to the supported drug list.`)
        }
        const errData = await response.json().catch(() => ({}))
        throw new Error(errData.detail || `Server error: ${response.status}`)
      }

      const data = await response.json()
      setResults(data.results)
    } catch (err) {
      if (err.message === 'Failed to fetch') {
        setError('Cannot connect to MoleculAR backend. Ensure the FastAPI server is running on port 8000.')
      } else {
        setError(err.message)
      }
    } finally {
      setIsLoading(false)
    }
  }

  return (
    <div className="panel">
      <h2>Drug Search</h2>
      <p style={{ color: 'var(--text-secondary)', marginBottom: '1.5rem' }}>
        Rank approved drugs by similarity of structure, protein targets and side effects.
      </p>

      <form onSubmit={handleSearch} style={{ display: 'flex', gap: '1rem', alignItems: 'flex-end', marginBottom: '2rem' }}>
        <div className="form-group" style={{ margin: 0, flex: 1 }}>
          <label htmlFor="drugInput">Query Drug Name</label>
          <input 
            id="drugInput"
            type="text" 
            className="form-control" 
            value={query}
            onChange={(e) => setQuery(e.target.value)}
            placeholder="e.g., SILDENAFIL"
            disabled={isLoading}
          />
        </div>
        <div className="form-group" style={{ margin: 0, width: '120px' }}>
          <label htmlFor="topNInput">Top N</label>
          <input 
            id="topNInput"
            type="number" 
            min="1"
            max="100"
            className="form-control" 
            value={topN}
            onChange={(e) => setTopN(e.target.value)}
            disabled={isLoading}
          />
        </div>
        <button type="submit" className="btn-primary" disabled={isLoading || !query.trim()}>
          {isLoading ? 'Searching...' : 'Search'}
        </button>
      </form>

      {error && (
        <div className="alert alert-danger">
          <strong>Error: </strong> {error}
        </div>
      )}

      {results && (
        <div>
          <h3>Results for {query.toUpperCase()}</h3>
          <table className="data-table">
            <thead>
              <tr>
                <th style={{ width: '50px' }}>Rank</th>
                <th>Candidate Drug</th>
                <th>Total Score</th>
                <th>Struct / Target / SE</th>
                {/* <th>Known Indication</th> */}
              </tr>
            </thead>
            <tbody>
              {results.map((res, idx) => (
                <tr key={idx}>
                  <td>#{idx + 1}</td>
                  <td style={{ fontWeight: 500 }}>{res.candidate}</td>
                  <td>
                    <span className="score-badge">{res.score.toFixed(4)}</span>
                  </td>
                  <td style={{ color: 'var(--text-secondary)' }}>
                    {res.struct_sim.toFixed(3)} / {res.target_sim.toFixed(3)} / {res.side_effect_sim.toFixed(3)}
                  </td>
                  {/* <td>{res.indication}</td> */}
                </tr>
              ))}
            </tbody>
          </table>
        </div>
      )}
    </div>
  )
}

export default ModeA
