import { useState, useEffect, useCallback } from 'react';

/**
 * Hook to manage URL search parameters without react-router-dom.
 * Supports:
 *   mode: 'drug' | 'target' (default: 'drug')
 *   q: search query string
 *   topN: number of candidates to request from API (default: 50 for smooth client-side pagination)
 */
export function useUrlState() {
  const readParams = useCallback(() => {
    const params = new URLSearchParams(window.location.search);
    const mode = params.get('mode') === 'target' ? 'target' : 'drug';
    const q = params.get('q') || '';
    const topN = parseInt(params.get('n') || '50', 10);
    return { mode, q, topN: isNaN(topN) ? 50 : topN };
  }, []);

  const [urlState, setUrlState] = useState(readParams);

  // Sync state when user clicks Back / Forward
  useEffect(() => {
    const handlePopState = () => {
      setUrlState(readParams());
    };
    window.addEventListener('popstate', handlePopState);
    return () => window.removeEventListener('popstate', handlePopState);
  }, [readParams]);

  const updateUrl = useCallback((nextState, replace = false) => {
    setUrlState((prev) => {
      const merged = typeof nextState === 'function' ? nextState(prev) : { ...prev, ...nextState };
      const params = new URLSearchParams();

      if (merged.mode === 'target') {
        params.set('mode', 'target');
      } else {
        params.set('mode', 'drug');
      }

      if (merged.q) {
        params.set('q', merged.q);
      }

      if (merged.topN && merged.topN !== 50) {
        params.set('n', merged.topN.toString());
      }

      const queryString = params.toString();
      const newUrl = queryString ? `?${queryString}` : window.location.pathname;

      if (replace) {
        window.history.replaceState(null, '', newUrl);
      } else {
        window.history.pushState(null, '', newUrl);
      }

      return merged;
    });
  }, []);

  return [urlState, updateUrl];
}
