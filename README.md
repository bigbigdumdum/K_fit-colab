Copyright (C) 2026 Mukundan S

# K_fit-colab

A Google Colab notebook that superimposes protein structures (PDB or mmCIF)
onto a reference structure by Kearsley fitting of atoms you choose. The
calculations are done by the [K_fit](https://github.com/bigbigdumdum/K_fit) engine.

**[Open the notebook in Colab](https://colab.research.google.com/github/bigbigdumdum/K_fit-colab/blob/main/K_fit.ipynb)**

## How to use

Run the notebook cells from top to bottom:

1. **Install:** sets up K_fit and its dependencies.
2. **Input:** upload files or enter PDB IDs, one row per structure. The first
   row is the reference; the others are fitted onto it. Optionally set the
   cutoffs SOdc and RMSDc (default 1.5 Å each).
3. **Input checks:** shows the atoms, residues, chains and hetero residues
   of each structure.
4. **Atom selection:** for each structure pick model/chain, residue and atom
   and click **Add atom** (**Add N, CA, C** adds a residue's backbone). Atoms
   are paired in order: the i-th atom of each target is fitted onto the i-th
   atom of the reference. Click **Check pairs** when done.
5. **Fitting:** shows RMSD, structure overlap and FIT / NO FIT for each
   structure.
6. **Download:** a zip with the fitted structures, a report and CSV files.

A fit is reported as FIT only if the RMSD is at most RMSDc and every selected
atom lies within SOdc of its reference atom.

## Atom IDs

Each atom has a **matomid**, `<model number>-<atom id>`, for example `1-245`.
You can type matomids instead of picking atoms, or upload a **pairs CSV**:
the first column is the reference, each further column one target, with the
structure names as headers and one atom pair per row.

```csv
ref.pdb,model_a.pdb,1ABC.cif
1-10,1-12,1-305
1-18,1-20,1-313
1-25,1-27,1-320
```

## More options

- **Fit all models:** fits every model of a multi-model file (e.g. NMR)
  separately. With a single input, all models are fitted onto the model you
  selected atoms from.
- **Add another atom set from this structure:** fit two parts of the same
  file onto each other.
- **Copy matomids to:** reuse one structure's selection for another.

## Output files

| File | Content |
|---|---|
| `<name>_fit.pdb` / `.cif` | Fitted structure, same format as the input (not written if no fit passed) |
| `report.txt` | Inputs with SHA-256 checksums, cutoffs, and for each fit the atom count, RMSD, structure overlap, FIT / NO FIT and the transformation |
| `<name>__pairs.csv` | The exact atom pairs used |
| `transforms.csv` | Rotation R and translation t for each fit (x' = R·x + t) |

## Running locally

The notebook also runs in Jupyter. Clone both repositories side by side,
then use either a Python virtual environment:

```bash
git clone https://github.com/bigbigdumdum/K_fit.git
git clone https://github.com/bigbigdumdum/K_fit-colab.git
python3 -m venv .venv && source .venv/bin/activate
pip install ./K_fit ./K_fit-colab jupyterlab
jupyter lab K_fit-colab/K_fit.ipynb
```

or a conda environment:

```bash
conda create -n k_fit -c conda-forge python=3.12 numpy numba biopython ipywidgets jupyterlab
conda activate k_fit
pip install ./K_fit ./K_fit-colab
jupyter lab K_fit-colab/K_fit.ipynb
```

Results are written to `k_fit_work/` in the folder Jupyter was started from.

## License

GNU General Public License v3.0 or later. See [LICENSE](LICENSE).
The K_fit engine is licensed separately under the MIT License.
