# Copyright (C) 2026 Mukundan S
# SPDX-License-Identifier: GPL-3.0-or-later
"""Colour schemes for the atom selection lists and table.

Residues: RasMol "amino" colours (by residue type), darkened where needed so
they are readable as text on a white background. Non-standard residues get a
yellow background with black text. Waters are light blue.

Atoms: CPK / Jmol element colours, darkened where needed for readability.

To change a colour, edit the dictionaries below; ``legend_html`` follows
automatically.
"""

from __future__ import annotations

from html import escape

# Residue type groups (RasMol amino colours, text-friendly variants).
RESIDUE_GROUPS = [
    ("acidic", ("ASP", "GLU"), "#E60A0A"),
    ("basic", ("LYS", "ARG"), "#145AFF"),
    ("histidine", ("HIS",), "#6262C2"),
    ("hydroxyl", ("SER", "THR"), "#E08600"),
    ("amide", ("ASN", "GLN"), "#00A0A0"),
    ("sulfur", ("CYS", "MET"), "#A89800"),
    ("aromatic", ("PHE", "TYR"), "#3232AA"),
    ("tryptophan", ("TRP",), "#B45AB4"),
    ("aliphatic", ("LEU", "VAL", "ILE"), "#0F820F"),
    ("small", ("ALA", "GLY"), "#707070"),
    ("proline", ("PRO",), "#C07050"),
]
RESIDUE_COLORS = {name: color for _group, names, color in RESIDUE_GROUPS for name in names}

# Nucleotides count as standard residues (shown in dark grey).
NUCLEOTIDES = {"A", "C", "G", "U", "I", "DA", "DC", "DG", "DT", "DI"}
NUCLEOTIDE_COLOR = "#404040"

WATER_NAMES = {"HOH", "WAT", "H2O", "DOD", "D2O"}
WATER_COLOR = "#4A9BD6"

NONSTANDARD_STYLE = "background-color: #FFFF00; color: #000000;"

# Element colours (CPK / Jmol, text-friendly variants).
ELEMENT_COLORS = {
    "C": "#505050",
    "N": "#2040D0",
    "O": "#E00000",
    "S": "#B09000",
    "P": "#E07000",
    "H": "#909090",
    "D": "#909090",
    "SE": "#C07800",
    "F": "#1FA01F", "CL": "#1FA01F", "BR": "#A52A2A", "I": "#7A0099",
    "FE": "#C05000", "ZN": "#6A6AA0", "MG": "#2A9A00", "CA": "#3D8C3D",
    "MN": "#9C7AC7", "CU": "#C88033", "NA": "#8A3AC0", "K": "#8A3AC0",
}
OTHER_ELEMENT_COLOR = "#A0208F"


def is_standard_residue(name: str) -> bool:
    """True for the 20 amino acids, nucleotides and water."""
    return name in RESIDUE_COLORS or name in NUCLEOTIDES or name in WATER_NAMES


def residue_style(name: str) -> str:
    """Return the CSS for a residue name (text colour, or yellow highlight)."""
    if name in RESIDUE_COLORS:
        return f"color: {RESIDUE_COLORS[name]};"
    if name in NUCLEOTIDES:
        return f"color: {NUCLEOTIDE_COLOR};"
    if name in WATER_NAMES:
        return f"color: {WATER_COLOR};"
    return NONSTANDARD_STYLE


def element_style(element: str) -> str:
    """Return the CSS text colour for an element symbol."""
    return f"color: {ELEMENT_COLORS.get(element.upper(), OTHER_ELEMENT_COLOR)};"


def list_style_html(css_class: str, styles: list) -> str:
    """Return a <style> block colouring the options of a list widget.

    ``css_class`` is the class added to the ipywidgets Select; ``styles`` holds
    one CSS string per option, in option order.
    """
    rules = [f".{css_class} option:nth-child({i}) {{ {style} font-weight: 600; }}"
             for i, style in enumerate(styles, start=1)]
    return "<style>\n" + "\n".join(rules) + "\n</style>"


def _chip(text: str, style: str) -> str:
    """Return one legend entry."""
    return f'<span style="{style} font-weight: 600; padding: 0 3px;">{escape(text)}</span>'


def legend_html() -> str:
    """Return a short HTML legend for residue and element colours."""
    residues = " ".join(_chip("/".join(n.title() for n in names), f"color: {color};")
                        for _group, names, color in RESIDUE_GROUPS)
    elements = " ".join(_chip(e, f"color: {ELEMENT_COLORS[e]};")
                        for e in ("C", "N", "O", "S", "P", "H"))
    return (f"<div style='line-height: 1.8'><b>Residues:</b> {residues} "
            f"{_chip('HOH', f'color: {WATER_COLOR};')} "
            f"{_chip('non-standard', NONSTANDARD_STYLE)}<br>"
            f"<b>Atoms:</b> {elements} {_chip('other', element_style('XX'))}</div>")
