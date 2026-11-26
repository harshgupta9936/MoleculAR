"""
experiments/verify_restore.py — Verifies R5a through R5g after git stash apply + git rm --cached drugrepo.db.
Saves complete raw output to experiments/logs/restore_verify.log.
"""

import hashlib
import os
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
BACKUP_ROOT = Path(os.environ.get("MOLECULAR_BACKUP_ROOT", str(ROOT.parent / "MoleculAR_BACKUP_20261110")))
LOG_PATH = ROOT / "experiments" / "logs" / "restore_verify.log"
STATUS_BEFORE = ROOT / "experiments" / "logs" / "status_before_restore.txt"


def run_cmd(cmd):
    p = subprocess.run(cmd, cwd=ROOT, capture_output=True, text=True, encoding="utf-8", errors="replace")
    return p.returncode, p.stdout.strip(), p.stderr.strip()


def main():
    lines = []
    def log(msg=""):
        print(msg)
        lines.append(str(msg))

    log("=== R5 RESTORE VERIFICATION ===")

    # Parse stash show --name-status
    _, stash_ns, _ = run_cmd(["git", "stash", "show", "--name-status", "stash@{0}"])
    m_files = []
    d_files = []
    for line in stash_ns.splitlines():
        parts = line.split("\t", 1)
        if len(parts) == 2:
            st, path = parts[0].strip(), parts[1].strip()
            if st == "M":
                m_files.append(path)
            elif st == "D":
                d_files.append(path)

    # R5a: M files hash check
    log("\n--- R5a: Modified (M) files hash verification ---")
    log(f"{'File':<30} {'WorkingTree_Hash':<42} {'Stash_Blob_Hash':<42} {'Match'}")
    log("-" * 122)
    all_m_match = True
    for f in m_files:
        _, wt_h, _ = run_cmd(["git", "hash-object", f])
        _, ls_out, _ = run_cmd(["git", "ls-tree", "stash@{0}", f])
        st_h = ls_out.split()[2] if ls_out and len(ls_out.split()) >= 3 else "MISSING"
        match = (wt_h == st_h)
        if not match:
            all_m_match = False
        log(f"{f:<30} {wt_h:<42} {st_h:<42} {match}")
    log(f"All M files match stash: {all_m_match} ({len(m_files)}/{len(m_files)})")

    # R5b: D files absence check
    log("\n--- R5b: Deleted (D) files absence verification ---")
    log(f"{'File':<45} {'Exists_On_Disk'}")
    log("-" * 60)
    all_d_absent = True
    for f in d_files:
        exists = (ROOT / f).exists()
        if exists:
            all_d_absent = False
        log(f"{f:<45} {exists}")
    log(f"All D files absent on disk: {all_d_absent} ({len(d_files)}/{len(d_files)})")

    # R5c: SHA-256 of drugrepo.db
    log("\n--- R5c: Root drugrepo.db SHA-256 ---")
    db_sha = hashlib.sha256((ROOT / "drugrepo.db").read_bytes()).hexdigest().upper()
    expected_sha = "D88A2333743AF2DBE859B3717764B6A645B9786CAA6F059CC9558E29C7FF0DA7"
    log(f"drugrepo.db SHA-256: {db_sha}")
    log(f"Expected SHA-256:    {expected_sha}")
    log(f"Match: {db_sha == expected_sha}")

    # R5d: git stash list & git branch -a
    log("\n--- R5d: Stash and Rescue Branch Check ---")
    _, slist, _ = run_cmd(["git", "stash", "list"])
    _, blist, _ = run_cmd(["git", "branch", "-a"])
    log("$ git stash list\n" + slist)
    log("\n$ git branch -a\n" + blist)

    # R5e: Untracked list diff vs status_before_restore.txt
    log("\n--- R5e: Untracked Files Diff vs status_before_restore.txt ---")
    raw_bytes = STATUS_BEFORE.read_bytes()
    before_text = raw_bytes.decode("utf-16") if raw_bytes.startswith((b"\xff\xfe", b"\xfe\xff")) else raw_bytes.decode("utf-8")
    before_untracked = []
    for l in before_text.splitlines():
        if l.startswith("?? "):
            before_untracked.append(l[3:].strip())
    _, now_status, _ = run_cmd(["git", "status", "--short"])
    now_untracked = []
    for l in now_status.splitlines():
        if l.startswith("?? "):
            now_untracked.append(l[3:].strip())

    before_set = set(before_untracked)
    now_set = set(now_untracked)
    in_before_not_now = sorted(list(before_set - now_set))
    in_now_not_before = sorted(list(now_set - before_set))
    log(f"Untracked count before restore: {len(before_untracked)}")
    log(f"Untracked count after restore:  {len(now_untracked)}")
    log(f"In before_restore but NOT now: {in_before_not_now}")
    log(f"In now but NOT before_restore: {in_now_not_before}")
    if in_before_not_now:
        log("Checking .gitignore for items in before_restore but not now:")
        gitignore_txt = (ROOT / ".gitignore").read_text(encoding="utf-8", errors="ignore")
        for line in gitignore_txt.splitlines():
            log(f"  .gitignore: {line}")
        for item in in_before_not_now:
            clean_item = item.strip('"').rstrip('/')
            log(f"  Test-Path on disk '{clean_item}': {(ROOT / clean_item).exists()}")

    # R5f: Content check of restored pipeline.py, server.py, requirements.txt
    log("\n--- R5f: Content Check of Restored Files ---")
    pipe_lines = (ROOT / "pipeline.py").read_text(encoding="utf-8", errors="ignore").splitlines()
    log("pipeline.py keying lines:")
    for idx, l in enumerate(pipe_lines, 1):
        if any(k in l for k in ["SELECT ", "drug_features[", "chembl_id", "generic_base_name"]):
            log(f"  L{idx}: {l}")

    srv_lines = (ROOT / "server.py").read_text(encoding="utf-8", errors="ignore").splitlines()
    log("\nserver.py CORS / allow_origins lines (L11-24):")
    for idx in range(10, min(24, len(srv_lines))):
        log(f"  L{idx+1}: {srv_lines[idx]}")

    req_txt = (ROOT / "requirements.txt").read_text(encoding="utf-8", errors="ignore")
    fastapi_lines = [l for l in req_txt.splitlines() if l.lower().startswith("fastapi")]
    xgboost_lines = [l for l in req_txt.splitlines() if l.lower().startswith("xgboost")]
    log(f"\nrequirements.txt contains fastapi: {fastapi_lines}")
    log(f"requirements.txt contains xgboost: {xgboost_lines}")

    # Check data_pipeline/*.py for hardcoded C:\Users paths as well
    dp_hits = 0
    for pyf in (ROOT / "data_pipeline").glob("*.py"):
        if "C:\\Users" in pyf.read_text(encoding="utf-8", errors="ignore"):
            dp_hits += 1
    log(f"data_pipeline/*.py hardcoded C:\\Users path hits: {dp_hits}")

    # R5g: Backup folder MoleculAR_BACKUP_20261110 check + source comparison
    log("\n--- R5g: Backup Folder & Source File Count / Bytes Check ---")
    def scan_dir(p: Path):
        count = 0
        total = 0
        for root_d, dirs, files in os.walk(p):
            for fn in files:
                fp = Path(root_d) / fn
                count += 1
                total += fp.stat().st_size
        return count, total

    bak_cnt, bak_bytes = scan_dir(BACKUP_ROOT)
    src_cnt, src_bytes = scan_dir(ROOT)
    s_set = {p.relative_to(ROOT).as_posix() for p in ROOT.rglob("*") if p.is_file()}
    b_set = {p.relative_to(BACKUP_ROOT).as_posix() for p in BACKUP_ROOT.rglob("*") if p.is_file()}
    log(f"Backup MoleculAR_BACKUP_20261110: files={bak_cnt}, bytes={bak_bytes}")
    log(f"Backup matches 26734 & 31547551986: files_match={bak_cnt == 26734}, bytes_match={bak_bytes == 31547551986}")
    log(f"Current Source MoleculAR:         files={src_cnt}, bytes={src_bytes}")
    log(f"Delta Source vs Backup:           files={src_cnt - bak_cnt:+d}, bytes={src_bytes - bak_bytes:+d}")
    log(f"  Files in Backup only (24 deleted by stash apply): {len(b_set - s_set)}")
    log(f"  Files in Source only (new logs/scripts/.git refs since backup): {sorted(s_set - b_set)}")

    LOG_PATH.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"\nSaved {LOG_PATH}")


if __name__ == "__main__":
    main()
