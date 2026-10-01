# Copyright (C) 2026 Mukundan S
# SPDX-License-Identifier: GPL-3.0-or-later
"""ipywidgets panels for the notebook stages.

Each ``build_*`` function returns a widget to ``display`` in one notebook
cell. The panels only read and write the shared ``NotebookSession``; all
logic lives in ``session.py`` and in the K_fit engine.

Stage 2  build_input_panel     input rows, SOdc / RMSDc, Load button
Stage 3  show_input_summary    atom / residue / chain counts, hetero residues
Stage 4  build_selection_panel coloured atom selectors, selection table, pairs CSV
Stage 5  build_fit_panel       Run button and statistics table
Stage 6  download_results      zip download (google.colab.files)

To change the layout of a stage, edit its ``build_*`` function only.
"""

from __future__ import annotations

import itertools
import traceback
from html import escape

import ipywidgets as w
from IPython.display import display

from K_fit.checker import parse_matomid_list
from K_fit.errors import KFitError

from . import colors
from .session import InputSpec, NotebookSession, SessionError

STRUCTURE_FILE_TYPES = ".pdb,.ent,.cif,.mmcif"


# Small helpers

def read_upload(upload: w.FileUpload) -> tuple:
    """Return (file name, bytes) of the first file in a FileUpload, or (None, None).

    Handles both ipywidgets 7 (dict value) and ipywidgets 8 (tuple value).
    """
    value = upload.value
    if not value:
        return None, None
    if isinstance(value, dict):                                   # ipywidgets 7
        name, item = next(iter(value.items()))
        return name, bytes(item["content"])
    item = value[0]                                               # ipywidgets 8
    return item["name"], bytes(item["content"])


def run_safely(output: w.Output, action) -> None:
    """Run ``action()`` and print user errors plainly (tracebacks for bugs) into ``output``."""
    with output:
        output.clear_output()
        try:
            action()
        except (SessionError, KFitError, ValueError, KeyError, OSError) as err:
            print(f"ERROR: {err}")
        except Exception:                                         # noqa: BLE001
            traceback.print_exc()


def atom_option_label(atom) -> str:
    """Return the atom list label, e.g. ``CA [2]`` or ``CB1 altloc A (occ 0.85) [1020]``.

    Shows the atom name and the atom id from the file (not the matomid).
    """
    text = atom.atom_name
    if atom.altloc:
        text += f" altloc {atom.altloc} (occ {atom.occupancy:.2f})"
    return f"{text} [{atom.atom_id}]"


def residue_option_label(residue_key: tuple) -> str:
    """Return the residue list label, e.g. ``ALA 42`` or ``HOH 501A``."""
    _chain, number, icode, name = residue_key
    return f"{name} {number}{icode.strip()}"


# Stage 2: input

def build_input_panel(session: NotebookSession) -> w.Widget:
    """Return the Stage 2 panel: input rows, cutoffs and the Load button.

    Each row has a file upload and a text field for a PDB ID (or the path of
    a file uploaded through the Colab file browser). The first row is the
    reference. "Add row" adds another input.
    """
    rows = []
    rows_box = w.VBox()

    def add_row(_=None):
        role = "Reference" if not rows else f"Target {len(rows)}"
        upload = w.FileUpload(accept=STRUCTURE_FILE_TYPES, multiple=False,
                              description="Upload")
        text = w.Text(placeholder="PDB ID (e.g. 1UBQ) or file path")
        rows.append((upload, text))
        rows_box.children = rows_box.children + (
            w.HBox([w.Label(role, layout=w.Layout(width="90px")), upload, text]),)

    add_row()
    add_button = w.Button(description="Add row", icon="plus")
    add_button.on_click(add_row)

    so_field = w.BoundedFloatText(value=session.so_cutoff, min=0.01, max=100.0, step=0.1,
                                  description="SOdc (A)")
    rmsd_field = w.BoundedFloatText(value=session.rmsd_cutoff, min=0.01, max=100.0, step=0.1,
                                    description="RMSDc (A)")
    load_button = w.Button(description="Load structures", button_style="primary")
    output = w.Output()

    def load(_):
        def action():
            session.so_cutoff = so_field.value
            session.rmsd_cutoff = rmsd_field.value
            specs = []
            for upload, text in rows:
                name, data = read_upload(upload)
                specs.append(InputSpec(upload_name=name, upload_bytes=data, text=text.value))
            for message in session.load_inputs(specs):
                print("Loaded", message)
            print(f"SOdc = {session.so_cutoff} A, RMSDc = {session.rmsd_cutoff} A")
            print("Run the Stage 3 cell next.")
        run_safely(output, action)

    load_button.on_click(load)
    return w.VBox([rows_box, add_button, w.HBox([so_field, rmsd_field]), load_button, output])


# Stage 3: input checks

def show_input_summary(session: NotebookSession) -> None:
    """Print the atom, residue and chain counts and hetero residues of each structure."""
    for index, structure in enumerate(session.structures):
        role = "Reference" if index == 0 else "Target"
        print(f"[{role}] {structure.name}  ({structure.file_format}, {structure.source}, "
              f"SHA-256 {structure.sha256[:16]}...)")
        print(structure.summary().as_text())
        print()
    if not session.structures:
        print("No structures loaded. Run Stage 2 first.")


# Stage 4: atom selection

_selector_ids = itertools.count(1)


def selection_table_html(structure, matomids: list) -> str:
    """Return the selected atoms as an HTML table (#, matomid, model, chain, residue, atom).

    Residue names use residue colours and atom names use element colours
    (see colors.py). matomids not found in the structure are shown in red.
    """
    if not matomids:
        return "<i>No atoms selected yet.</i>"
    cell = "padding: 1px 8px; border-bottom: 1px solid #ddd;"
    head = "<tr>" + "".join(f"<th style='{cell} text-align: left'>{h}</th>"
                            for h in ("#", "matomid", "model", "chain", "residue", "atom")) + "</tr>"
    rows = []
    for i, matomid in enumerate(matomids, start=1):
        atom = structure.atoms.get(matomid)
        if atom is None:
            rows.append(f"<tr><td style='{cell}'>{i}</td><td style='{cell}'>{escape(matomid)}</td>"
                        f"<td colspan='4' style='{cell} color: #E00000'>not found in "
                        f"{escape(structure.name)}</td></tr>")
            continue
        residue = f"{atom.residue_name} {atom.residue_number}{atom.insertion_code.strip()}"
        atom_text = atom.atom_name + (f" ({atom.altloc})" if atom.altloc else "")
        rows.append(
            f"<tr><td style='{cell}'>{i}</td><td style='{cell}'>{matomid}</td>"
            f"<td style='{cell}'>{atom.model_number}</td><td style='{cell}'>{escape(atom.chain_id)}</td>"
            f"<td style='{cell} {colors.residue_style(atom.residue_name)} font-weight: 600'>"
            f"{escape(residue)}</td>"
            f"<td style='{cell} {colors.element_style(atom.element)} font-weight: 600'>"
            f"{escape(atom_text)}</td></tr>")
    return ("<div style='max-height: 320px; overflow-y: auto'>"
            f"<table style='border-collapse: collapse; font-size: 13px'>{head}{''.join(rows)}"
            "</table></div>")


class AtomSelector:
    """Atom selection widgets for one structure (Stage 4).

    Model/chain dropdown -> residue list -> atom list. The residue list is
    coloured by residue type (non-standard residues highlighted yellow), the
    atom list by element. Lists are used instead of closed dropdowns because
    browsers do not reliably show colours inside a closed dropdown.

    The selection is stored in ``session.selections[structure.name]`` and
    shown as a table that grows with each addition. matomids can also be
    typed in and added, or pasted from another structure with "Copy matomids to".
    Attributes are public so tests can drive the widgets.
    """

    BACKBONE_ATOMS = ("N", "CA", "C")     # atoms added by "Add N, CA, C", in this order

    def __init__(self, session: NotebookSession, structure, role: str,
                 on_duplicate=None, on_copy_to=None):
        """Build the widgets.

        ``on_duplicate(name)`` is called by "Add another atom set from this
        structure"; ``on_copy_to(source_name, target_name)`` by the Copy button.
        """
        self.session = session
        self.structure = structure
        uid = next(_selector_ids)
        self.residue_class = f"kfit-res-{uid}"
        self.atom_class = f"kfit-atom-{uid}"

        self.chain_dd = w.Dropdown(
            options=[(f"Model {m} / Chain {c}", (m, c))
                     for m in structure.models() for c in structure.chains(m)],
            description="Model/chain")
        self.residue_list = w.Select(rows=8, description="Residue",
                                     layout=w.Layout(width="260px"))
        self.atom_list = w.Select(rows=8, description="Atom", layout=w.Layout(width="320px"))
        self.residue_list.add_class(self.residue_class)
        self.atom_list.add_class(self.atom_class)
        self.residue_style = w.HTML()          # <style> blocks colouring the two lists
        self.atom_style = w.HTML()

        self.add_atom_button = w.Button(description="Add atom", icon="plus")
        self.add_backbone_button = w.Button(description="Add N, CA, C",
                                            tooltip="Add the backbone atoms N, CA, C of the residue")
        self.remove_last_button = w.Button(description="Remove last")
        self.clear_button = w.Button(description="Clear")
        self.id_entry = w.Text(placeholder="type matomids, e.g. 1-2 1-10 1-18",
                               layout=w.Layout(width="320px"))
        self.add_typed_button = w.Button(description="Add typed")
        self.message = w.HTML()
        self.table = w.HTML()
        self.copy_target_dd = w.Dropdown(description="Copy matomids to",
                                         style={"description_width": "initial"})
        self.copy_to_button = w.Button(description="Copy")

        self.chain_dd.observe(self.refresh_residues, names="value")
        self.residue_list.observe(self.refresh_atoms, names="value")
        self.add_atom_button.on_click(
            lambda _: self.append([self.atom_list.value] if self.atom_list.value else []))
        self.add_backbone_button.on_click(self.add_backbone)
        self.remove_last_button.on_click(lambda _: self.set_ids(self.selected_ids()[:-1]))
        self.clear_button.on_click(lambda _: self.set_ids([]))
        self.add_typed_button.on_click(self.add_typed)

        buttons = [self.add_atom_button, self.add_backbone_button,
                   self.remove_last_button, self.clear_button]
        if on_duplicate is not None:
            self.duplicate_button = w.Button(
                description="Add another atom set from this structure",
                layout=w.Layout(width="auto"))
            self.duplicate_button.on_click(lambda _: on_duplicate(structure.name))
            buttons.append(self.duplicate_button)
        if on_copy_to is not None:
            self.copy_to_button.on_click(
                lambda _: self.copy_target_dd.value
                and on_copy_to(structure.name, self.copy_target_dd.value))

        self.refresh_residues()
        self.refresh_table()
        self.widget = w.VBox([
            w.HTML(f"<b>{escape(role)}: {escape(structure.name)}</b>"),
            self.chain_dd,
            w.HBox([self.residue_list, self.atom_list]),
            self.residue_style, self.atom_style,
            w.HBox(buttons),
            w.HBox([self.id_entry, self.add_typed_button]),
            self.message,
            self.table,
            w.HBox([self.copy_target_dd, self.copy_to_button]),
        ], layout=w.Layout(border="1px solid #ccc", padding="6px", margin="4px 0"))

    # list contents

    @staticmethod
    def _set_options(select: w.Select, options: list) -> None:
        """Replace the options and select the first (ipywidgets 8 leaves none selected)."""
        select.options = options
        if options and select.value is None:
            select.index = 0

    def refresh_residues(self, _=None) -> None:
        """Fill the residue list for the chosen model/chain and colour it."""
        model, chain = self.chain_dd.value
        residues = self.structure.residues(model, chain)
        self.residue_style.value = colors.list_style_html(
            self.residue_class, [colors.residue_style(r[3]) for r in residues])
        self._set_options(self.residue_list, [(residue_option_label(r), r) for r in residues])
        self.refresh_atoms()

    def refresh_atoms(self, _=None) -> None:
        """Fill the atom list for the chosen residue and colour it by element."""
        model, _chain = self.chain_dd.value
        residue = self.residue_list.value
        atoms = self.structure.atoms_in_residue(model, residue) if residue else []
        self.atom_style.value = colors.list_style_html(
            self.atom_class, [colors.element_style(a.element) for a in atoms])
        self._set_options(self.atom_list, [(atom_option_label(a), a.matomid) for a in atoms])

    # selection

    def selected_ids(self) -> list:
        """Return the matomids selected so far, in fitting order."""
        return list(self.session.selections.get(self.structure.name, []))

    def set_ids(self, matomids: list) -> None:
        """Replace the selection and redraw the table."""
        self.session.selections[self.structure.name] = list(matomids)
        self.refresh_table()

    def append(self, matomids: list) -> None:
        """Add matomids to the end of the selection, skipping ones already selected."""
        current = self.selected_ids()
        self.set_ids(current + [m for m in matomids if m not in current])

    def add_backbone(self, _=None) -> None:
        """Add N, CA and C of the chosen residue (highest-occupancy altloc of each)."""
        model, _chain = self.chain_dd.value
        residue = self.residue_list.value
        atoms = self.structure.atoms_in_residue(model, residue) if residue else []
        chosen, missing = [], []
        for name in self.BACKBONE_ATOMS:
            candidates = [a for a in atoms if a.atom_name == name]
            if candidates:
                chosen.append(max(candidates, key=lambda a: a.occupancy).matomid)
            else:
                missing.append(name)
        self.append(chosen)
        self.message.value = (f"<span style='color: #E00000'>{residue_option_label(residue)} has no "
                              f"{', '.join(missing)}</span>" if residue and missing else "")

    def set_copy_targets(self, names: list) -> None:
        """Set the structures offered by "Copy matomids to" (all except this one)."""
        options = [n for n in names if n != self.structure.name]
        self.copy_target_dd.options = options
        if options and self.copy_target_dd.value is None:   # ipywidgets 8 selects nothing
            self.copy_target_dd.index = 0

    def paste_ids(self, matomids: list, source_name: str) -> None:
        """Put matomids copied from another structure into the entry field.

        They are not added yet: the user checks them and clicks "Add typed".
        """
        self.id_entry.value = " ".join(matomids)
        self.message.value = (f"Pasted {len(matomids)} matomid(s) from {escape(source_name)}. "
                              "Check them and click <b>Add typed</b>.")

    def add_typed(self, _=None) -> None:
        """Add the matomids typed in the entry field; report unknown ones."""
        typed = parse_matomid_list(self.id_entry.value)
        unknown = [m for m in typed if m not in self.structure.atoms]
        self.append([m for m in typed if m in self.structure.atoms])
        self.message.value = (f"<span style='color: #E00000'>Not in {escape(self.structure.name)}: "
                              f"{escape(', '.join(unknown))}</span>" if unknown else "")
        if not unknown:
            self.id_entry.value = ""

    def refresh_table(self) -> None:
        """Redraw the selection table from the session."""
        self.table.value = selection_table_html(self.structure, self.selected_ids())


def build_selection_panel(session: NotebookSession) -> w.Widget:
    """Return the Stage 4 panel: one selector per structure, pairs CSV, checks.

    The i-th selected atom of each target is fitted onto the i-th atom of the
    reference. A pairs CSV replaces all manual selections. The AtomSelectors
    are available as ``panel.selectors`` (dict: structure name -> AtomSelector).
    """
    selectors = {}
    selectors_box = w.VBox()
    output = w.Output()

    def add_selector(structure, role):
        selector = AtomSelector(session, structure, role, on_duplicate=make_copy,
                                on_copy_to=copy_ids)
        selectors[structure.name] = selector
        selectors_box.children = selectors_box.children + (selector.widget,)
        for each in selectors.values():
            each.set_copy_targets(list(selectors))

    def copy_ids(source_name, target_name):
        selectors[target_name].paste_ids(selectors[source_name].selected_ids(), source_name)

    def make_copy(name):
        def action():
            copy = session.add_copy(name)
            add_selector(copy, "Target (copy)")
            print(f"Added {copy.name}, a copy of {name}")
        run_safely(output, action)

    for index, structure in enumerate(session.structures):
        add_selector(structure, "Reference" if index == 0 else "Target")

    fit_all = w.Checkbox(value=session.fit_all_models, description="Fit all models",
                         indent=False)
    fit_all.observe(lambda change: setattr(session, "fit_all_models", change["new"]),
                    names="value")

    csv_upload = w.FileUpload(accept=".csv", multiple=False, description="Pairs CSV")
    csv_button = w.Button(description="Apply pairs CSV")

    def apply_csv(_):
        def action():
            _name, data = read_upload(csv_upload)
            if not data:
                raise SessionError("upload a pairs CSV first")
            for message in session.apply_pairs_csv(data.decode("utf-8-sig")):
                print(message)
            for structure in session.structures:
                if structure.name not in selectors:
                    add_selector(structure, "Target (copy)")
                selectors[structure.name].refresh_table()
        run_safely(output, action)

    csv_button.on_click(apply_csv)

    check_button = w.Button(description="Check pairs", button_style="primary")

    def check(_):
        def action():
            results = session.check()
            for result in results:
                print(result.as_text())
                print()
            if all(r.ok for r in results):
                print("All selections are valid. Run the Stage 5 cell next.")
        run_safely(output, action)

    check_button.on_click(check)
    help_text = w.HTML(
        "Select atoms in fitting order: the i-th atom of each target is fitted onto the "
        "i-th atom of the reference. Atoms are listed as <i>name [atom id]</i>. "
        "<i>Add N, CA, C</i> adds the backbone atoms of the chosen residue. "
        "<i>Copy matomids to</i> pastes this structure's matomids into another structure's "
        "entry field. "
        "<i>Fit all models</i> fits every model of each target with the same selection; "
        "with a single input, every model is fitted onto the selected model."
        + colors.legend_html())
    panel = w.VBox([help_text, fit_all, w.HBox([csv_upload, csv_button]), selectors_box,
                    check_button, output])
    panel.selectors = selectors
    return panel


# Stage 5: fitting

def build_fit_panel(session: NotebookSession) -> w.Widget:
    """Return the Stage 5 panel: Run button and the superposition statistics."""
    run_button = w.Button(description="Run fitting", button_style="primary")
    output = w.Output()

    def run(_):
        def action():
            result = session.run()
            print(result.summary_text())
            warnings = [f"{c.target_name}: {msg}" for c in result.checks for msg in c.warnings]
            if warnings:
                print("\nWarnings:")
                print("\n".join(f"  {msg}" for msg in warnings))
            print(f"\nOutputs in {session.output_dir}. Run the Stage 6 cell to download.")
        run_safely(output, action)

    run_button.on_click(run)
    return w.VBox([run_button, output])


# Stage 6: download

def download_results(session: NotebookSession) -> None:
    """Download the results zip in Colab, or print its path elsewhere."""
    if session.result is None or not session.result.zip_path:
        print("No results yet. Run Stage 5 first.")
        return
    try:
        from google.colab import files                           # only available in Colab
    except ImportError:
        print(f"Results zip: {session.result.zip_path}")
        return
    files.download(session.result.zip_path)


def show(widget: w.Widget) -> None:
    """Display a panel (small wrapper so notebook cells stay one line)."""
    display(widget)
