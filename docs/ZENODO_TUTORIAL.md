# Zenodo DOIs explained — a tutorial for this repo

> Written 2026-04-25 to clarify the four DOIs currently associated with this work and what to cite where.

---

## Your four DOIs

When you connect a GitHub repo to Zenodo, Zenodo creates a **two-tier DOI system**:

| Tier | DOI | Type | Resolves to |
| --- | --- | --- | --- |
| Concept | `10.5281/zenodo.19046141` | umbrella | always the latest version |
| Version | `10.5281/zenodo.19046142` | v1.0 (2026-03-16) | snapshot of v1.0 |
| Version | `10.5281/zenodo.19347286` | v1.1 (2026-03-31) | snapshot of v1.1 |
| Version | `10.5281/zenodo.19669876` | v1.2 (2026-04-20) | snapshot of v1.2 |

**Note on the concept DOI:** Zenodo reserves the number *before* the first version DOI as the concept DOI. So your concept is `19046141` and v1.0 is `19046142` (one higher). This is how Zenodo always works — you can't see the concept DOI on the version pages directly, but it's exposed via the API.

---

## What each tier is for

### Concept DOI (`19046141`) — the "moving target"
Resolves to whatever the most recent version is. Today that's v1.2. Next time you release v1.3, the same concept DOI quietly switches to point at v1.3 without you doing anything.

**Citing the concept DOI says:** "Cite the current version of this work, whatever it happens to be."

### Version DOIs — frozen snapshots
Each release gets its own permanent DOI. These never change. v1.0 will always be v1.0; v1.1 will always be v1.1.

**Citing a version DOI says:** "Cite this exact frozen state of the work, so the reader can reproduce my results against the same files I had."

---

## Which one should the manuscript cite?

There's a real choice here, and both options are defensible:

### Option A — concept DOI in the manuscript (recommended for most papers)
- **Pros:** Citation always resolves to the latest version. If you fix a typo in a script after publication and re-release, readers automatically get the corrected version. No reader will ever land on a stale v1.0 archive.
- **Cons:** The exact bytes a future reader downloads may differ from what the original reviewers saw. For most research that's fine; the *interpretation* doesn't depend on which trailing-whitespace patch was current.

### Option B — version DOI in the manuscript (for reproducibility-critical work)
- **Pros:** Reviewers can verify your numbers against the exact same files you had. Readers years from now can reproduce against the historical state.
- **Cons:** If you find a bug post-publication, the manuscript still points at the buggy version. You'd need to publish an erratum to update the citation.

**For your paper, my recommendation is Option B (cite a specific version DOI)**, because:
1. The work explicitly markets `statistics_audit.json` as a per-claim reproducibility artifact. That promise is only honored if the reader can pull the exact file you audited against.
2. *GBE* readers and reviewers who care enough to check the repository care about reproducibility. They'll appreciate landing on the frozen artifact.
3. Concept DOIs are good for software libraries where users want the latest fixes. Research papers are different — you're documenting a specific moment.

If you ever publish a v2 of the paper itself (a major revision with new analyses), the new paper would cite the new Zenodo version. The old paper would still link to the old Zenodo version, which is exactly right.

---

## So which version DOI?

You currently have v1.0, v1.1, v1.2. **Use the most recent version that contains all the data, scripts, and figures referenced in the manuscript.**

- v1.0 (2026-03-16) — too old; pre-dates Draft 6 numbers
- v1.1 (2026-03-31) — has the exhaustive *P. dulcis* BLASTp results and final Draft 6 manuscript
- v1.2 (2026-04-20) — adds peer-review fixes, technical walkthrough, Figure 4 layout fix
- (a future v1.3) — would also include the BLAST 478-cluster validation and Scope of Generalization table from today's revision

**Today's recommendation:** Cut a v1.3 release after the senior reviewer fixes are committed, and cite v1.3's DOI in the paper. That way the paper's "frozen state" matches the reviewer-incorporated state.

---

## What about the README?

The README is a different artifact than the paper. The README is the front door to the GitHub repo and should always reflect the current state. Two reasonable approaches:

- **Cite concept DOI in README** (good — auto-updates the badge to point to the latest snapshot)
- **Cite the latest version DOI in README** (also good — matches what the paper cites)

I lean toward citing the concept DOI in the README. The README is for someone who arrived via GitHub and just wants the latest archive; the paper is for someone who arrived via citation and needs a specific frozen state.

---

## What's currently inconsistent

| Location | DOI cited | What you might want |
| --- | --- | --- |
| README badge (line 3) | v1.0 (`19046142`) | Either the concept (`19046141`) or whichever version you cut next |
| Paper Data Availability | v1.1 (`19347286`) | The latest version DOI you commit to (e.g. v1.3 once cut) |
| Cover letter | v1.1 (`19347286`) | Same as the paper |

---

## Suggested next steps (when you're ready)

1. Commit today's reviewer fixes (manuscript edits, README updates, BLAST 478 results, supp tables, cover letter, statistics_audit additions)
2. Push to GitHub
3. Cut a v1.3 release on GitHub — Zenodo will mint a new version DOI automatically
4. Once the new DOI is live (~5 minutes), update:
   - **README badge** → concept DOI `10.5281/zenodo.19046141`
   - **Paper Data Availability** → new v1.3 version DOI
   - **Cover letter** → same v1.3 version DOI
5. Regenerate Word docs
6. Final commit ("update DOIs to v1.3")

This produces a clean, reviewer-friendly state where:
- The paper points at a frozen archive that matches what reviewers will see
- The README points at "always-latest"
- Anyone who lands on v1.0 or v1.1 will see Zenodo's "newer version available" banner directing them forward

---

## TL;DR

- **4 DOIs is normal.** 1 concept (umbrella) + 3 version DOIs (one per release).
- **Concept DOI = always-latest.** Version DOIs = frozen snapshots.
- **For your paper, cite a specific version DOI** so reviewers see the exact files you used.
- **Cut a v1.3 after today's edits**, then cite v1.3 in the paper.
- The README can use either; I lean toward the concept DOI for "current state" framing.
