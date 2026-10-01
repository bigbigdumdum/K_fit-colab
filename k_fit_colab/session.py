# Copyright (C) 2026 Mukundan S
# SPDX-License-Identifier: GPL-3.0-or-later
"""Notebook state and logic, kept apart from the widgets so it can be tested.

``NotebookSession`` holds everything the stages share: the loaded
structures (index 0 is the reference), the atom selections, the parameters
and the last result. The widgets module only reads and writes this object.
"""

from __future__ import annotations

import csv
import io
import os
import re
from dataclasses import dataclass

from K_fit.checker import check_atom_pairs, read_pairs_csv
from K_fit.fetch import PDB_ID_PATTERN, fetch_pdb_entry
from K_fit.fitter import DEFAULT_RMSD_CUTOFF, DEFAULT_SO_CUTOFF
from K_fit.parser import load_structure, make_unique_name
from K_fit.pipeline import FitJob, run_superposition, self_fit_job

def default_work_dir() -> str:
    """Return /content/k_fit in Colab, else ./k_fit_work in the current folder.

    /content only exists (and is writable) in a Colab runtime.
    """
    if os.path.isdir("/content") and os.access("/content", os.W_OK):
        return "/content/k_fit"
    return os.path.join(os.getcwd(), "k_fit_work")


class SessionError(Exception):
    """Raised for user errors (missing input, empty selection, ...)."""


@dataclass
class InputSpec:
    """One input row: an uploaded file and/or the text field.

    ``text`` may be a PDB ID (downloaded from RCSB) or the path of a file
    already on disk (for example uploaded through the Colab file browser).
    An uploaded file takes priority over the text.
    """

    upload_name: str | None = None
    upload_bytes: bytes | None = None
    text: str = ""

    def is_empty(self) -> bool:
        """True when the row has neither an upload nor text."""
        return not self.upload_bytes and not self.text.strip()


class NotebookSession:
    """Shared state of one notebook run."""

    def __init__(self, work_dir: str | None = None):
        """Create the session; inputs go to work_dir/inputs, outputs to work_dir/outputs.

        ``work_dir`` defaults to ``default_work_dir()``.
        """
        work_dir = work_dir or default_work_dir()
        self.work_dir = work_dir
        self.input_dir = os.path.join(work_dir, "inputs")
        self.output_dir = os.path.join(work_dir, "outputs")
        self.structures = []        # [reference, target, target, ...]
        self.selections = {}        # structure name -> ordered list of matomids
        self.so_cutoff = DEFAULT_SO_CUTOFF
        self.rmsd_cutoff = DEFAULT_RMSD_CUTOFF
        self.fit_all_models = False
        self.result = None

    # Structures

    @property
    def reference(self):
        """The reference structure (first input)."""
        if not self.structures:
            raise SessionError("no structures loaded yet (run Stage 2)")
        return self.structures[0]

    def names(self) -> list:
        """Names of all loaded structures, reference first."""
        return [s.name for s in self.structures]

    def get(self, name: str):
        """Return the loaded structure with this name."""
        for structure in self.structures:
            if structure.name == name:
                return structure
        raise SessionError(f"structure {name} is not loaded")

    def load_inputs(self, specs: list) -> list:
        """Load every non-empty input row; return one message per structure.

        Replaces any previously loaded structures. Duplicate file names are
        renamed ``<stem>_<n><ext>``. Uploaded files are saved to
        ``input_dir/<row number>/<file name>`` so that the original file name is kept.
        """
        specs = [spec for spec in specs if not spec.is_empty()]
        if not specs:
            raise SessionError("provide at least one structure (upload or PDB ID)")
        self.structures, self.selections, self.result = [], {}, None
        messages = []
        for row, spec in enumerate(specs, start=1):
            path, source = self._resolve_input(row, spec)
            name = make_unique_name(os.path.basename(path), self.names())
            structure = load_structure(path, source=source, name=name)
            self.structures.append(structure)
            role = "reference" if row == 1 else "target"
            renamed = f" (renamed from {structure.file_name})" if name != structure.file_name else ""
            messages.append(f"{role}: {name}{renamed}")
        return messages

    def _resolve_input(self, row: int, spec: InputSpec) -> tuple:
        """Return (path, source) for one input row, saving or downloading as needed."""
        row_dir = os.path.join(self.input_dir, str(row))
        if spec.upload_bytes:
            os.makedirs(row_dir, exist_ok=True)
            path = os.path.join(row_dir, os.path.basename(spec.upload_name or "upload.pdb"))
            with open(path, "wb") as handle:
                handle.write(spec.upload_bytes)
            return path, "upload"
        text = spec.text.strip()
        for candidate in (text, os.path.join(self.input_dir, text), os.path.join("/content", text)):
            if os.path.isfile(candidate):
                return candidate, "upload"
        if PDB_ID_PATTERN.match(text):
            return fetch_pdb_entry(text, row_dir), "pdb_id"
        raise SessionError(f"row {row}: '{text}' is neither a file nor a 4-character PDB ID")

    def add_copy(self, name: str):
        """Add a copy of a loaded structure as a new target and return it.

        Used when two atom sets come from the same input file. The copy is
        named like a duplicate input (``<stem>_<n><ext>``).
        """
        copy = self.get(name).copy(make_unique_name(self.get(name).name, self.names()))
        self.structures.append(copy)
        return copy

    # Selections

    def set_selection(self, name: str, matomids: list) -> None:
        """Store the ordered matomids selected for one structure."""
        self.get(name)
        self.selections[name] = list(matomids)

    def apply_pairs_csv(self, text: str) -> list:
        """Replace all selections with those of a pairs CSV; return messages.

        The first column must be the reference. A header like ``<stem>_<n><ext>``
        that is not loaded but whose base file is loaded creates a copy, just
        like a duplicate input.
        """
        rows = [row for row in csv.reader(io.StringIO(text)) if any(c.strip() for c in row)]
        header = [cell.strip() for cell in rows[0]] if rows else []
        messages = []
        for name in header[1:]:
            if name in self.names():
                continue
            base = _base_name(name)
            if base in self.names():
                copy = self.add_copy(base)
                if copy.name != name:
                    raise SessionError(f"pairs CSV column {name}: expected the next copy to be "
                                       f"named {copy.name}; number copies consecutively")
                messages.append(f"created copy {name} of {base}")
        reference_name, selections = read_pairs_csv(text, self.names())
        if reference_name != self.reference.name:
            raise SessionError(f"pairs CSV first column is {reference_name}, but the reference "
                               f"(first input) is {self.reference.name}")
        self.selections = {}
        for target_name, (ref_ids, tgt_ids) in selections.items():
            self.selections[reference_name] = ref_ids
            self.selections[target_name] = tgt_ids
        messages.append(f"pairs CSV: {len(selections)} target(s), "
                        f"{len(self.selections[reference_name])} pair(s) each")
        return messages

    # Jobs

    def build_jobs(self) -> list:
        """Turn the selections into FitJobs (see CLAUDE.md Stage 3 rules).

        * The reference selection is required.
        * Every target with a selection becomes one job; a loaded target without a
          selection is an error, so nothing is skipped silently.
        * With no target selections, "fit all models" fits every model of a copy
          of the reference onto the selected reference model; otherwise at
          least two atom sets are required.
        """
        reference = self.reference
        ref_ids = self.selections.get(reference.name, [])
        if not ref_ids:
            raise SessionError(f"no atoms selected for the reference {reference.name}")
        jobs = []
        for target in self.structures[1:]:
            tgt_ids = self.selections.get(target.name, [])
            if not tgt_ids:
                raise SessionError(f"no atoms selected for target {target.name}; select atoms "
                                   "or reload without it")
            jobs.append(FitJob(target=target, reference_ids=ref_ids, target_ids=tgt_ids))
        if not jobs:
            if not self.fit_all_models:
                raise SessionError("at least two sets of atoms are needed: add a target, or use "
                                   "'Add another atom set from this structure', or tick "
                                   "'Fit all models'")
            jobs.append(self_fit_job(reference, ref_ids, self.names()))
        return jobs

    def check(self) -> list:
        """Return the PairCheckResults of all jobs (used by the Check button)."""
        return [check_atom_pairs(self.reference, job.target, job.reference_ids, job.target_ids)
                for job in self.build_jobs()]

    def run(self):
        """Run the fitting and writing; store and return the RunResult."""
        self.result = run_superposition(
            self.reference, self.build_jobs(), self.output_dir,
            so_cutoff=self.so_cutoff, rmsd_cutoff=self.rmsd_cutoff,
            fit_all_models=self.fit_all_models)
        return self.result


def _base_name(name: str) -> str:
    """Return the original name of a duplicate: "model_2.pdb" -> "model.pdb"."""
    stem, ext = os.path.splitext(name)
    return re.sub(r"_\d+$", "", stem) + ext
