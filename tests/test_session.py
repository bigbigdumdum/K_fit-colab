# Copyright (C) 2026 Mukundan S
# SPDX-License-Identifier: GPL-3.0-or-later
"""Tests for k_fit_colab.session and a smoke test of the widgets.

Structure files are taken from the K_fit repository's tests/data folder
(set K_FIT_TEST_DATA to override). Tests are skipped if it is missing.
"""

import json
import os

import pytest

from k_fit_colab.session import InputSpec, NotebookSession, SessionError

HERE = os.path.dirname(__file__)
DATA_DIR = os.environ.get(
    "K_FIT_TEST_DATA", os.path.join(HERE, "..", "..", "K_fit", "tests", "data"))
pytestmark = pytest.mark.skipif(not os.path.isdir(DATA_DIR),
                                reason="K_fit test data not found")


def data_bytes(name):
    with open(os.path.join(DATA_DIR, name), "rb") as handle:
        return handle.read()


def upload(name):
    return InputSpec(upload_name=name, upload_bytes=data_bytes(name))


def ca_ids(structure, model=1):
    return [a.matomid for a in structure.atoms_in_model(model) if a.atom_name == "CA"]


@pytest.fixture
def session(tmp_path):
    return NotebookSession(work_dir=str(tmp_path))


def test_duplicate_names_are_renamed(session):
    messages = session.load_inputs([upload("1UBQ.pdb"), InputSpec(), upload("1UBQ.pdb")])
    assert session.names() == ["1UBQ.pdb", "1UBQ_1.pdb"]
    assert "renamed from 1UBQ.pdb" in messages[1]
    assert session.structures[1].file_name == "1UBQ.pdb"


def test_text_field_accepts_path(session):
    session.load_inputs([InputSpec(text=os.path.join(DATA_DIR, "1UBQ.cif"))])
    assert session.reference.source == "upload" and session.reference.file_format == "cif"


def test_bad_text_and_empty_rows(session):
    with pytest.raises(SessionError):
        session.load_inputs([InputSpec(), InputSpec(text="  ")])
    with pytest.raises(SessionError, match="neither a file nor"):
        session.load_inputs([InputSpec(text="not-a-file")])


def test_two_atom_sets_from_one_file(session):
    session.load_inputs([upload("1L2Y_3models.pdb")])
    ref = session.reference
    with pytest.raises(SessionError, match="at least two sets"):
        session.set_selection(ref.name, ca_ids(ref)) or session.build_jobs()
    copy = session.add_copy(ref.name)
    assert copy.name == "1L2Y_3models_1.pdb"
    session.set_selection(copy.name, ca_ids(copy, 2))
    session.so_cutoff = session.rmsd_cutoff = 5.0
    result = session.run()
    assert [f.target_model for f in result.fits] == [2] and result.fits[0].passed
    assert os.path.isfile(result.zip_path)


def test_fit_all_models_with_single_input(session):
    session.load_inputs([upload("1L2Y_3models.pdb")])
    session.set_selection(session.reference.name, ca_ids(session.reference))
    session.fit_all_models = True
    result = session.run()
    assert [f.target_model for f in result.fits] == [1, 2, 3]


def test_target_without_selection_is_error(session):
    session.load_inputs([upload("1UBQ.pdb"), upload("1UBQ.cif")])
    session.set_selection("1UBQ.pdb", ca_ids(session.reference))
    with pytest.raises(SessionError, match="no atoms selected for target"):
        session.build_jobs()


def test_pairs_csv_creates_copies(session):
    session.load_inputs([upload("1UBQ.pdb")])
    ids = ca_ids(session.reference)[:5]
    text = "1UBQ.pdb,1UBQ_1.pdb\n" + "".join(f"{i},{i}\n" for i in ids)
    messages = session.apply_pairs_csv(text)
    assert "created copy 1UBQ_1.pdb of 1UBQ.pdb" in messages
    assert session.selections["1UBQ_1.pdb"] == ids
    checks = session.check()
    assert checks[0].ok and len(checks[0].pairs) == 5


def test_pairs_csv_must_start_with_reference(session):
    session.load_inputs([upload("1UBQ.pdb"), upload("1UBQ.cif")])
    with pytest.raises(SessionError, match="first column"):
        session.apply_pairs_csv("1UBQ.cif,1UBQ.pdb\n1-2,1-2\n")


def test_widgets_smoke(session):
    """Build every panel and drive it with button clicks, as a user would."""
    widgets = pytest.importorskip("k_fit_colab.widgets")
    session.load_inputs([upload("1UBQ.pdb")])
    widgets.build_input_panel(session)

    panel = widgets.build_selection_panel(session)
    selector = panel.selectors["1UBQ.pdb"]
    assert selector.chain_dd.options[0][0] == "Model 1 / Chain A"
    assert selector.residue_list.options[0][0] == "MET 1"
    assert selector.atom_list.options[1][0] == "CA [2]"       # atom name [atom id]
    selector.add_backbone_button.click()                       # MET 1: N, CA, C
    assert [session.reference.atoms[m].atom_name for m in session.selections["1UBQ.pdb"]] \
        == ["N", "CA", "C"]
    assert selector.table.value.count("<tr>") == 1 + 3         # header + one row per atom
    selector.residue_list.value = selector.residue_list.options[1][1]   # GLN 2
    selector.add_atom_button.click()                           # first atom: N of GLN 2
    assert len(session.selections["1UBQ.pdb"]) == 4
    selector.remove_last_button.click()
    assert len(session.selections["1UBQ.pdb"]) == 3
    selector.id_entry.value = "1-8 1-99999"
    selector.add_typed_button.click()
    assert session.selections["1UBQ.pdb"][-1] == "1-8" and "1-99999" in selector.message.value

    selector.duplicate_button.click()                          # second atom set, same file
    copy = panel.selectors["1UBQ_1.pdb"]
    assert selector.copy_target_dd.options == ("1UBQ_1.pdb",)
    selector.copy_to_button.click()                            # paste into the copy's entry
    assert copy.id_entry.value.split() == selector.selected_ids()
    copy.add_typed_button.click()
    assert copy.selected_ids() == selector.selected_ids()
    session.so_cutoff = session.rmsd_cutoff = 1.5
    fit_panel = widgets.build_fit_panel(session)
    fit_panel.children[0].click()
    assert session.result is not None and session.result.fits[0].passed


def test_colouring(session):
    """Residue list, atom list and table use residue / element colours."""
    widgets = pytest.importorskip("k_fit_colab.widgets")
    from k_fit_colab import colors
    session.load_inputs([upload("1UBQ.pdb")])
    selector = widgets.build_selection_panel(session).selectors["1UBQ.pdb"]
    # MET 1 is the first residue (sulfur colour); waters are listed last.
    assert colors.RESIDUE_COLORS["MET"] in selector.residue_style.value
    assert colors.WATER_COLOR in selector.residue_style.value
    assert colors.ELEMENT_COLORS["N"] in selector.atom_style.value
    assert colors.residue_style("CRO") == colors.NONSTANDARD_STYLE
    assert colors.residue_style("ALA") != colors.NONSTANDARD_STYLE
    selector.add_atom_button.click()
    assert colors.RESIDUE_COLORS["MET"] in selector.table.value
    assert "non-standard" in colors.legend_html()


def test_gfp_by_pdb_id(session):
    """CLAUDE.md test case through the notebook logic: PDB IDs typed in two rows,
    6L26 fitted onto 1EMA by the 15 chromophore (CRO 66) carbons."""
    try:
        session.load_inputs([InputSpec(text="1EMA"), InputSpec(text="6l26")])
    except OSError as err:                     # FetchError or network problems
        pytest.skip(f"cannot download: {err}")
    except Exception as err:                   # noqa: BLE001
        if "download failed" in str(err):
            pytest.skip(str(err))
        raise
    assert session.names() == ["1EMA.cif", "6L26.cif"]
    assert all(s.source == "pdb_id" for s in session.structures)
    for structure in session.structures:
        carbons = {}
        for atom in structure.atoms_in_residue(1, ("A", 66, " ", "CRO")):
            if atom.element == "C" and (atom.atom_name not in carbons
                                        or atom.occupancy > carbons[atom.atom_name].occupancy):
                carbons[atom.atom_name] = atom
        order = [a.atom_name for a in session.reference.atoms_in_residue(1, ("A", 66, " ", "CRO"))
                 if a.element == "C" and not a.altloc]
        session.set_selection(structure.name, [carbons[name].matomid for name in order])
    result = session.run()
    fit = result.fits[0]
    assert fit.passed and fit.n_atoms == 15 and abs(fit.rmsd - 0.124) < 0.002
    assert "input mode: PDB ID" in open(result.report_path).read()


def test_notebook_is_valid_json():
    path = os.path.join(HERE, "..", "K_fit.ipynb")
    with open(path, encoding="utf-8") as handle:
        notebook = json.load(handle)
    sources = "".join("".join(c["source"]) for c in notebook["cells"])
    for call in ("build_input_panel", "show_input_summary", "build_selection_panel",
                 "build_fit_panel", "download_results"):
        assert call in sources


def test_add_backbone_altlocs_and_missing(session):
    """'Add N, CA, C' takes the highest-occupancy altloc and reports missing atoms."""
    widgets = pytest.importorskip("k_fit_colab.widgets")
    session.load_inputs([upload("1EJG_res1-8.pdb")])
    selector = widgets.build_selection_panel(session).selectors["1EJG_res1-8.pdb"]
    selector.add_backbone_button.click()                       # THR 1, altlocs A 0.82 / B 0.18
    atoms = [session.reference.atoms[m] for m in selector.selected_ids()]
    # N and CA have altlocs A (0.82) / B (0.18); C has none.
    assert [(a.atom_name, a.altloc) for a in atoms] == [("N", "A"), ("CA", "A"), ("C", "")]

    session.load_inputs([upload("1UBQ.pdb")])                  # a water has no backbone
    selector = widgets.build_selection_panel(session).selectors["1UBQ.pdb"]
    water = next(v for _label, v in selector.residue_list.options if v[3] == "HOH")
    selector.residue_list.value = water
    selector.add_backbone_button.click()
    assert selector.selected_ids() == [] and "has no N, CA, C" in selector.message.value
