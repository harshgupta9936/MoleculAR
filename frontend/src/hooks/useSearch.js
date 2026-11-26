import { useState, useEffect, useCallback } from 'react';

const CACHE = new Map();
const RECENT_KEY = 'molecular_recent_searches';

function getStoredRecent() {
  try {
    const saved = localStorage.getItem(RECENT_KEY);
    return saved ? JSON.parse(saved) : [];
  } catch {
    return [];
  }
}

function saveRecent(item) {
  try {
    const current = getStoredRecent();
    const filtered = current.filter(
      (r) => !(r.mode === item.mode && r.query.toUpperCase() === item.query.toUpperCase())
    );
    const updated = [item, ...filtered].slice(0, 10);
    localStorage.setItem(RECENT_KEY, JSON.stringify(updated));
    return updated;
  } catch {
    return [];
  }
}

export function useSearch(mode, query, topN = 50) {
  const [data, setData] = useState(null);
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState(null);
  const [recentSearches, setRecentSearches] = useState(getStoredRecent);
  const [retryCount, setRetryCount] = useState(0);

  useEffect(() => {
    const trimmed = (query || '').trim();
    if (!trimmed) {
      queueMicrotask(() => {
        setData(null);
        setError(null);
        setLoading(false);
      });
      return;
    }

    const cacheKey = `${mode}:${trimmed.toUpperCase()}:${topN}`;
    if (CACHE.has(cacheKey)) {
      const cached = CACHE.get(cacheKey);
      queueMicrotask(() => {
        setData(cached);
        setError(null);
        setLoading(false);
      });
      return;
    }

    let isMounted = true;
    const controller = new AbortController();

    // Defer state update to next microtask to prevent cascading render warning
    queueMicrotask(() => {
      if (isMounted) {
        setLoading(true);
        setError(null);
      }
    });

    const endpoint = mode === 'drug' ? '/api/rank' : '/api/target-search';
    const body =
      mode === 'drug'
        ? { drug_name: trimmed, top_n: topN }
        : { target_id: trimmed, top_n: topN };

    fetch(`http://127.0.0.1:8000${endpoint}`, {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify(body),
      signal: controller.signal,
    })
      .then(async (resp) => {
        if (!resp.ok) {
          const errData = await resp.json().catch(() => ({}));
          throw new Error(errData.detail || `Server returned ${resp.status}: ${resp.statusText}`);
        }
        return resp.json();
      })
      .then((result) => {
        if (!isMounted) return;
        CACHE.set(cacheKey, result);
        setData(result);
        setError(null);

        const updatedRecent = saveRecent({
          mode,
          query: trimmed.toUpperCase(),
          timestamp: Date.now(),
        });
        setRecentSearches(updatedRecent);
      })
      .catch((err) => {
        if (!isMounted || err.name === 'AbortError') return;
        setError(err.message || 'Network error occurred while fetching candidates.');
        setData(null);
      })
      .finally(() => {
        if (isMounted) {
          setLoading(false);
        }
      });

    return () => {
      isMounted = false;
      controller.abort();
    };
  }, [mode, query, topN, retryCount]);

  const clearRecent = useCallback(() => {
    try {
      localStorage.removeItem(RECENT_KEY);
      setRecentSearches([]);
    } catch {}
  }, []);

  const retry = useCallback(() => {
    const trimmed = (query || '').trim();
    if (trimmed) {
      CACHE.delete(`${mode}:${trimmed.toUpperCase()}:${topN}`);
      setRetryCount((c) => c + 1);
    }
  }, [mode, query, topN]);

  return {
    data,
    loading,
    error,
    recentSearches,
    clearRecent,
    retry,
  };
}
