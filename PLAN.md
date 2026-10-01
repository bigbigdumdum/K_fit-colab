Copyright (C) 2026 Mukundan S

# K_fit-colab: implementation plan

Colab front end (repo 2). All parsing, checking, fitting and writing is done
by the `K_fit` engine. This repo contains only user-interface code.

## Rules

- No scientific logic here. If a cell needs to compute something, the
  computation belongs in `K_fit`.
- UI with `ipywidgets` (pre-installed in Colab; works with 7.7 and 8.x).
  Results are downloaded with `google.colab.files`.
- Names are standard and describe what the code does. Every function has a
  docstring.
- Nothing is committed without the owner's permission.
- License: GPL-3.0-or-later. Each source file starts with
  `# Copyright (C) 2026 Mukundan S` and
  `# SPDX-License-Identifier: GPL-3.0-or-later`.

## Layout

```
K_fit-colab/
├── README.md, PLAN.md, LICENSE, pyproject.toml
├── K_fit.ipynb          # the notebook: one section per stage
├── k_fit_colab/
│   ├── __init__.py
│   ├── session.py       # NotebookSession: shared state and logic (no widgets)
│   ├── colors.py        # residue / element colour schemes and legend
│   └── widgets.py       # one build_* function per stage
└── tests/
    └── test_session.py  # session logic + widget smoke test (button clicks)
```

## Notebook stages

| Stage | Cell | Code |
|---|---|---|
| 1 Install | clone or pull both repos, pip install | notebook cell |
| 2 Input | rows (upload / PDB ID / path), Add row, SOdc, RMSDc, Load | `build_input_panel` |
| 3 Input checks | counts, hetero residues, waters | `show_input_summary` |
| 4 Atom selection | model/chain, residue and atom dropdowns; typed matomids; pairs CSV; copy button; Fit all models; Check pairs | `build_selection_panel` |
| 5 Fitting | Run fitting, statistics table, warnings | `build_fit_panel` |
| 6 Download | zip download | `download_results` |

## Milestones

| | Milestone | Status |
|---|---|---|
| C1 | Stage 1 and Stage 2 widgets, upload and PDB ID | done |
| C2 | Stage 3 summary | done |
| C3 | Selection methods, copies, pairs CSV, pair check | done |
| C4 | Fitting, results table, download | done |
| C5 | End-to-end run in Colab; **Open in Colab** badge in the README | open: needs the GitHub URLs |

13 tests pass locally (including the GFP 1EMA / 6L26 test by PDB ID) with ipywidgets 7.7.1 and 8.1.9.

## Open items

- The install cell assumes `https://github.com/bigbigdumdum/K_fit.git` and
  `https://github.com/bigbigdumdum/K_fit-colab.git`. Change `K_FIT_REPO` and
  `K_FIT_COLAB_REPO` if they differ.
- `ipywidgets.FileUpload` has not been tested in a live Colab runtime yet. The
  text field also accepts a file path, as a fallback.
