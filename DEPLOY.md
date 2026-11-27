# Deploying MoleculAR

MoleculAR runs as a decoupled FastAPI backend (`server.py`, hosted on Render) and React/Vite frontend (`frontend/`, hosted on Vercel). The scored SQLite database (`drugrepo.db`, 17.94 MB, expected SHA-256 `D88A2333743AF2DBE859B3717764B6A645B9786CAA6F059CC9558E29C7FF0DA7`) is not stored in Git because it integrates SIDER 4.1 (`CC BY-NC-SA 4.0`).

---

## 1. Host `drugrepo.db` Privately

Choose one private host so the database is never exposed in public source control:

### Option A: Private Hugging Face Dataset Repository
1. Create a **private** dataset repository on Hugging Face (for example `<hf-user>/molecular-db`).
2. Upload `drugrepo.db` (`SHA-256: D88A2333743AF2DBE859B3717764B6A645B9786CAA6F059CC9558E29C7FF0DA7`) to the root of that dataset repository.
3. Create a fine-grained **Read** user access token in Hugging Face settings (`https://huggingface.co/settings/tokens`) scoped to that private dataset.
4. Use these Render environment values:
   - `DB_URL`: `https://huggingface.co/datasets/<hf-user>/molecular-db/resolve/main/drugrepo.db`
   - `DB_URL_TOKEN`: `<your-hf-read-token>` (`server.py` sends it as `Authorization: Bearer <DB_URL_TOKEN>`).

### Option B: Private Object Storage (Cloudflare R2 / AWS S3 / Backblaze B2) or Private GitHub Release Asset
- **Pre-signed URL (R2 / S3 / B2):** Upload `drugrepo.db` to a private bucket and generate a long-lived HTTPS pre-signed download URL for `DB_URL` (leave `DB_URL_TOKEN` unset).
- **Private GitHub Release asset:** Upload `drugrepo.db` as a release asset in a separate **private** repository (such as your private archive repo), create a fine-grained PAT with **Contents: Read-only** permission scoped to that single repo, set `DB_URL` to `https://api.github.com/repos/<owner>/<private-repo>/releases/assets/<asset-id>`, and set `DB_URL_TOKEN` to that PAT.

At startup, if `drugrepo.db` is missing, `server.py` downloads `DB_URL`, checks its SHA-256 against `D88A2333743AF2DBE859B3717764B6A645B9786CAA6F059CC9558E29C7FF0DA7`, and deletes the file and exits non-zero on any mismatch.

---

## 2. Deploy Backend to Render

1. Push updates to `master` and connect the repository in Render as a **Blueprint** from `render.yaml` (or create a Python **Web Service** on the `free` plan — local peak RSS is ~165 MiB, within Render's 512 MB limit). Note that Render's free web service spins down after inactivity and wakes in about a minute on the next request.
2. Confirm settings (pre-filled by `render.yaml`):
   - **Runtime:** Python (`PYTHON_VERSION=3.14.5`; Render supports `3.14.x` natively and all packages in `requirements.txt` ship `cp314` `manylinux_2_28_x86_64` wheels on PyPI)
   - **Build Command:** `pip install -r requirements.txt`
   - **Start Command:** `uvicorn server:app --host 0.0.0.0 --port $PORT`
   - **Health Check Path:** `/api/health`
3. Set secrets/environment variables in Render:
   - `DB_URL` — private HTTPS download URL for `drugrepo.db`
   - `DB_URL_TOKEN` — bearer token if required by `DB_URL` (leave unset for pre-signed URLs)
   - `CORS_ORIGINS` — temporary placeholder until Vercel finishes, then set to `https://<your-vercel-app>.vercel.app` (comma-separated if multiple origins)
4. After deploy completes, confirm `GET https://<your-render-service>.onrender.com/api/health` returns `200 OK` with `"sha256_verified": true` and expected row counts (`drugs: 2267`, `targets: 7068`, `side_effects: 119227`, `crosswalk: 2267`).

---

## 3. Deploy Frontend to Vercel

1. Import the repository in Vercel on branch `master`.
2. Configure **Project Settings**:
   - **Root Directory:** `frontend`
   - **Framework Preset:** `Vite`
   - **Build Command:** `npm run build`
   - **Output Directory:** `dist`
3. Add **Environment Variable**:
   - `VITE_API_BASE_URL` = `https://<your-render-service>.onrender.com` (no trailing slash)
4. Deploy and copy your live `https://<your-vercel-app>.vercel.app` URL.
5. Set `CORS_ORIGINS=https://<your-vercel-app>.vercel.app` on Render and trigger a Render restart/redeploy so browser requests from Vercel succeed.
6. **Final check:** Open the deployed dashboard in your browser, open Developer Tools -> **Network** tab, run a search, and confirm `/api/rank`, `/api/target-search`, and `/api/structure` requests go to your `https://<your-render-service>.onrender.com` Render URL and not `127.0.0.1`.

---

## 4. Local Development (No New Environment Variables Required)

- When running locally with `drugrepo.db` present at the repository root:
  - `server.py` skips `DB_URL` and opens `drugrepo.db` read-only (`?mode=ro`).
  - `CORS_ORIGINS` defaults to `http://localhost:5173`, `http://127.0.0.1:5173`, `http://localhost:3000`, `http://127.0.0.1:3000`, `http://localhost:8000`, and `http://127.0.0.1:8000`.
  - `VITE_API_BASE_URL` defaults to `http://127.0.0.1:8000`.
