# New PC Migration Checklist — Pandoravirus Project

**Created:** 2026-05-04
**Scope:** Everything needed to bring the Pandoravirus project (1.0 + 2.0 paper repo + monorepo subprojects used by the paper) up on a new Windows machine, in submission-ready state.

For broader project portfolio migration (DDI-Risk, KRAS, GEL, etc.), this is just one section of the larger move. Reference your `project-folder-organization.md` (canonical: `C:\Users\Eddie\Code Projects\Admin\folder-organization\PROJECT-STRUCTURE.md`) for the full layout.

---

## 1. Software to install (in order)

| Tool | Version on this machine | Why | How to install |
|---|---|---|---|
| Git | (current) | Source control | https://git-scm.com/download/win |
| Git Bash (MINGW64) | (bundled with Git) | Bash shell | included with Git installer |
| Python | 3.11.9 | All analysis scripts | https://www.python.org/downloads/release/python-3119/ — match the **3.11.x** line for psycopg2 / scipy compatibility |
| PostgreSQL | 15 (currently at `D:\Program Files\PostgreSQL\15\`) | Database | https://www.postgresql.org/download/windows/ — pick **15.x** to keep dump-restore compatibility |
| Pandoc | 3.9.0.2 (currently at `C:\Users\Eddie\AppData\Local\Pandoc\`) | Word doc generation from markdown | https://pandoc.org/installing.html |
| GitHub CLI (`gh`) | 2.87.3 | Auth, PRs, releases | https://cli.github.com/ |
| BLAST+ | (locally installed) | Local BLAST runs (used in Pandoravirus paper validation) | https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ |
| VS Code | (current) | Editor | https://code.visualstudio.com/ |

**VS Code extensions already installed (per memory):** Code Spell Checker, ESLint, Prettier, GitLens, Tailwind CSS IntelliSense (`bradlc.vscode-tailwindcss`), Thunder Client, Live Share. Settings Sync is on for the GitHub account `mathamagician`, so most of this should auto-restore on the new PC after sign-in.

**git-filter-repo** (used during the credential cleanup) is `pip install git-filter-repo` if ever needed again. Not needed for normal workflow.

---

## 2. GitHub clones to make on new PC

Both repos are now on clean public history (force-pushed 2026-05-04). On the new PC, **re-clone fresh** rather than copying the old machine's clone — that's how the cleaned history lands authoritatively:

```bash
mkdir -p "C:\Users\Eddie\Code Projects"
cd "C:\Users\Eddie\Code Projects"

# Monorepo (contains Pandoravirus 1.0 plus 15 other research subprojects)
git clone https://github.com/mathamagician/genomics-research.git genomics

# Standalone v2.0 paper repo (cloned alongside, not inside, the monorepo —
# the monorepo's .gitignore already excludes Pandoravirus/Pandoravirus 2.0/)
cd "genomics/Pandoravirus"
git clone https://github.com/mathamagician/pandoravirus-denovo-genebirth.git "Pandoravirus 2.0"
```

**Important:** the monorepo and the v2.0 paper repo are *separate* GitHub repositories. The monorepo has the Pandoravirus 1.0 work + 15 other research projects; the standalone v2.0 paper repo is the one being submitted to GBE.

### Disk-size note for the monorepo

The local monorepo working tree is currently **~52 GB**, but most of that is local-only research data that's gitignored and not in the GitHub clone. Breakdown:

| Subproject | Local disk | In GitHub? |
|---|---|---|
| `ddi-risk/` | 37 GB | partial — `data/`, `logs/` gitignored |
| `Reconstruct HODDI/` | 14 GB | partial — `data/`, `raw_data/`, `logs/` gitignored |
| `genomic-enrichment-layer/` | 774 MB | yes (mostly) |
| `project-exodus/` | 485 MB | yes (mostly) |
| `Pandoravirus/` | 119 MB | yes |
| `mimivirus/` | 8 MB | yes |
| `ncldv_pipeline/` | 7 MB | yes |
| `.git/` (monorepo history) | 167 MB | (the repo itself) |

**The fresh `git clone` of the monorepo will likely be ~1-2 GB**, not 52 GB. The big sub-directories (DDI-Risk data, HODDI raw data) live on disk locally, not in git. Those need to migrate via copy of `D:\Research Data Library\` (per section 6) or per-project data restore — not via git clone.

---

## 3. PostgreSQL database migration

The Pandoravirus paper depends on a local Postgres database named `pandoravirus` (22 tables, ~4M rows: position-indexed sequences, gene annotations, intergenic regions, codon usage, k-mer frequencies, derived analytical tables).

### On the OLD PC — dump the database

```powershell
& "D:\Program Files\PostgreSQL\15\bin\pg_dump.exe" -U postgres -d pandoravirus -F c -f "C:\Users\Eddie\Code Projects\genomics\Pandoravirus\pandoravirus_db.dump"
```
(Will prompt for password `postgres_pizza_87`.) Custom format (`-F c`) preserves indexes, sequences, ownership. Expect a few hundred MB to a couple GB depending on table sizes.

Copy `pandoravirus_db.dump` to the new PC (USB, network share, OneDrive, whatever).

### On the NEW PC — restore the database

After installing PostgreSQL 15 and starting the service:

```powershell
# First create the database
& "C:\Program Files\PostgreSQL\15\bin\psql.exe" -U postgres -c "CREATE DATABASE pandoravirus;"

# Then restore from the dump
& "C:\Program Files\PostgreSQL\15\bin\pg_restore.exe" -U postgres -d pandoravirus "C:\path\to\pandoravirus_db.dump"
```

**Note:** the new PC's Postgres password should be set during install. Pick a new strong password — don't reuse `postgres_pizza_87`. Update the `.env` files in step 4 to match.

---

## 4. Local config files (`.env`) — recreate on new PC

These files are **gitignored**, so they won't be in the cloned repos. You must recreate them after cloning. Templates are committed as `.env.example` in each subproject.

### Files to create on new PC

| Location | Source template | Contents |
|---|---|---|
| `Pandoravirus/Pandoravirus 1.0/.env` | (no template; create from scratch) | `DB_HOST=localhost`<br>`DB_PORT=5432`<br>`DB_NAME=pandoravirus`<br>`DB_USER=postgres`<br>`DB_PASSWORD=<your-new-password>` |
| `mimivirus/.env` | `mimivirus/.env.example` | Same five vars |
| `ncldv_pipeline/.env` | `ncldv_pipeline/.env.example` | Same five vars |
| `Pandoravirus 2.0/.env` | `Pandoravirus 2.0/.env.example` | Same five vars |

For each `.env.example`, just `cp .env.example .env` and edit the password.

---

## 5. Python dependencies

For each subproject that uses Python:

### Pandoravirus 2.0 (paper repo)
```bash
cd "C:\Users\Eddie\Code Projects\genomics\Pandoravirus\Pandoravirus 2.0"
python -m venv .venv
.venv\Scripts\activate     # PowerShell: .venv\Scripts\Activate.ps1
pip install -r requirements.txt
```

### Pandoravirus 1.0 / monorepo subprojects
The monorepo doesn't have a single root `requirements.txt`. Subprojects with their own dependencies (mimivirus/, ncldv_pipeline/, etc.) — install per directory. Core packages used across the Pandoravirus work:

```bash
pip install psycopg2 biopython numpy scipy scikit-learn pandas matplotlib seaborn pyrodigal
```

Versions on this PC (for reproducibility): Python 3.11.9, psycopg2 2.9.11, biopython 1.86, scipy 1.17.1.

---

## 6. Research Data Library (Tier-1 datasets)

Per your `data-architecture.md` memory, the Tier-1 data lives at `D:\Research Data Library\` (~250 GB across 86 datasets). This is **separate from the GitHub repos** — it doesn't migrate via git.

**Migrate by:** external drive or network copy of the entire `D:\Research Data Library\` tree. If the new PC's data drive is also `D:\`, paths in scripts will work without modification. If the data drive letter differs, search-and-replace `D:\Research Data Library\` across scripts on the new PC.

For Pandoravirus specifically, the relevant raw data is:
- `D:\Research Data Library\raw\genomics\` — viral genomes (FASTA), BAM, VCF, etc.
- `D:\Research Data Library\raw\ddi\rxnorm\` — RxNorm (used by DDI project, not Pandoravirus)

Most Pandoravirus working data lives **inside the repos**: `Pandoravirus 2.0/data/` (FASTA), `results/` (JSON), and the Postgres database itself.

---

## 7. Authentication / accounts (manual reconnects on new PC)

| Service | Setup needed |
|---|---|
| GitHub | Sign in to VS Code; `gh auth login`; SSH keys regenerated or copied if you use SSH (currently `https://` cloning works fine via `gh` token) |
| GitHub CLI (`gh`) | `gh auth login` interactively, then accept the browser auth flow |
| NCBI Entrez | Email is set in scripts (`Eddie.Bradford@gmail.com`); the api_key is `None` in code (3 req/s default). If you use a higher rate, the API key would be in `.env` or env var. |
| UMLS / UTS | API key currently in `D:\Research Data Library\... config\.env` and `genomics/faers-v3-pipeline/config/.env` (per memory). Move with the data library. |
| Zenodo | Already set up; v1.3 release cut today will auto-mint v1.3 DOI under concept `10.5281/zenodo.19046141`. No new-PC action needed. |

---

## 8. Pre-submission state of the Pandoravirus paper (current as of 2026-05-04)

This is what should look identical on the new PC after you re-clone and migrate the database:

- v2.0 paper repo HEAD: `d41e827` ("Update Zenodo DOI references to concept DOI")
- Monorepo HEAD: `787440d` ("Pandoravirus cleanup: scrub embedded credentials from cross-NCLDV pipeline scripts")
- v1.3 GitHub release cut at https://github.com/mathamagician/pandoravirus-denovo-genebirth/releases/tag/v1.3
- Three Word docs ready to submit: [Pandoravirus_combo_paper_draft6_GBE.docx](../Pandoravirus_combo_paper_draft6_GBE.docx), [supplementary_materials_GBE.docx](../supplementary_materials_GBE.docx), [cover_letter_gbe_GBE.docx](../cover_letter_gbe_GBE.docx)
- `.zenodo.json`, `statistics_audit.json`, all script outputs intact
- Postgres password rotated to `postgres_pizza_87` (will need to be set fresh on new PC during Postgres install)

---

## 9. Cleanup items for new PC (not blockers)

These were flagged during the Pandoravirus cleanup as out-of-scope and left for context-appropriate follow-up. Address whenever you next work on each project on the new PC:

- **`curable-api/.env.example`** — was hidden by over-broad `.env.*` gitignore (now fixed). Untracked in monorepo. Commit when you next work on curable-api.
- **`giant-virus-atlas/`** — large untracked block of new pipeline code (scripts/01–08, diamond binaries, etc.). Decide: keep in monorepo or extract to its own repo. Diamond `.exe` binaries should likely be gitignored regardless.
- **`genomic-enrichment-layer/`** — untracked changes per memory at v0.7.1 BETA. Confirm whether these should land in the GEL standalone repo or stay in the monorepo.

These are project hygiene, not security.

---

## 10. Deferred follow-ups (already in queue, not migration-blocking)

- **Old Zenodo versions (v1.0 / v1.1 / v1.2)** still contain the rotated-and-invalidated `pandora2026` credential. Either mark them "superseded" via the Zenodo dashboard description, or request curator deletion via Zenodo support citing security rotation. The new v1.3 (in progress) is clean.
- **GBE submission** — three Word docs are ready in [docs/](../). Cover letter, manuscript, supplementary. Do NOT include `REVIEWER_RESPONSE_2026-04-26.md` per prior decision.
- **Pandoravirus Paper 2** (heterologous reporter assay) — see prior conversation; deferred to post-publication.

---

## 11. Quick verification once new PC is set up

After completing steps 1–7, run:

```bash
cd "C:\Users\Eddie\Code Projects\genomics\Pandoravirus\Pandoravirus 2.0"
.venv\Scripts\activate
python -c "import os; from scripts.phase1b_corrected_analysis import get_conn; c = get_conn(); print('DB OK:', c is not None)"
```

Should print `DB OK: True`. If you get `KeyError: 'DB_PASSWORD'`, your `.env` isn't being loaded — check the file exists and the venv is activated.

To regenerate Word docs:
```bash
cd "C:\Users\Eddie\Code Projects\genomics\Pandoravirus\Pandoravirus 2.0\docs"
pandoc Pandoravirus_combo_paper_draft6_GBE.md -o Pandoravirus_combo_paper_draft6_GBE.docx
```

Should produce a 45 KB .docx with no errors.

---

*End of migration checklist.*
