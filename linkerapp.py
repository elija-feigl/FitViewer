#!/usr/bin/env python
# -*- coding: utf-8 -*-

import ipywidgets as widgets
import MDAnalysis as mda

import mrcfile
from pathlib import Path
from IPython.display import display

from core.linker import Linker
from core.project import Project
from core.pdbCorrection import PDB_Corr, Logic


__authors__ = ["Elija Feigl"]
__version__ = "1.0"
__license__ = "GPL-3.0"

__descr__ = """
    VIEWERTOOL:
    This module defines the classes used to define the connectivity and
    geometry of a DNA structure.

    A DNA structure consists of a number of scaffold and staple strands
    (DNA origami), or oligo strands alone, bound together to form a designed
    geometric shape.
    """
__status__ = "Development"
__maintainer__ = "Elija Feigl"
__email__ = "elija.feigl@tum.de"


class Viewer(object):
    def __init__(self, folder, files):
        name = files[0].stem
        self.project = Project(input=Path(folder),
                               output=Path(folder) / "analysis",
                               name=name,
                               files=files,
                               )
        self.linker = Linker(self.project)
        self.selection_scaffold = None
        self.selection_staples = None

    def select_by_helixandbase(self, helices, bases):

        self._parse_selection(bases, helices)
        link = self.linker.create_linkage()

        atoms_scaffold = mda.AtomGroup([], self.linker.fit.u)
        atoms_staple = mda.AtomGroup([], self.linker.fit.u)

        for base in self.selection_scaffold:
            atoms_scaffold += self.linker.fit.u.residues[
                link.DidFid[link.DhpsDid[(base.h, base.p, True)]]].atoms

        for base in self.selection_staples:
            atoms_staple += self.linker.fit.u.residues[link.DidFid[
                link.DhpsDid[(base.h, base.p, False)]]].atoms

        atoms_selection = atoms_scaffold + atoms_staple

        return atoms_selection, atoms_scaffold, atoms_staple, link.Dcolor

    def _parse_selection(self, base_selection, helix_selection):
        # get helices from design
        helices = []
        for idx in helix_selection:
            helix = self.linker.design.design.structure_helices_map[idx]
            helices.append(helix)
        # get bases from helices
        self.selection_scaffold = []
        self.selection_staples = []
        for helix in helices:
            selection_scaffold_add = [
                base for base in helix.scaffold_bases
                if base.p in base_selection and base.num_deletions == 0
                and base.num_insertions == 0]
            selection_staples_add = [
                base for base in helix.staple_bases
                if base.p in base_selection and base.num_deletions == 0
                and base.num_insertions == 0]

            self.selection_scaffold += selection_scaffold_add
            self.selection_staples += selection_staples_add

        return

    def make_sliders(self):

        layout_w = widgets.Layout(width='900px', height='70px')
        layout_h = widgets.Layout(width='100px')
        layout_box = widgets.Layout(width='1000px', border='1px solid black')
        layout_button = widgets.Layout(width='30px', height='30px')
        design = self.linker.design.design

        base_p = [base.p for base in self.linker.design.allbases]
        maxb = max(base_p)
        minb = min(base_p)
        row = [coor[0] for coor, _ in
               design.structure_helices_coord_map.items()]
        minr = min(row)
        maxr = max(row)
        col = [coor[1] for coor, _ in
               design.structure_helices_coord_map.items()]
        minc = min(col)
        maxc = max(col)

        lattice = ("square" if design.lattice_type == 0
                   else "honeycomb")

        helix_buttons = []
        if lattice == "square":
            for r in range(minr, maxr + 1):
                row = []
                for c in range(minc, maxc + 1):
                    try:
                        id = str(
                            design.structure_helices_coord_map[(r, c)].id)
                    except KeyError:
                        id = ""

                    but = widgets.ToggleButton(
                        description=id, layout=layout_button)
                    row.append(but)
                helix_buttons.append(row)
        else:  # TODO: -low- improve
            for r in range(minr, maxr + 1):
                row = []
                for c in range(minc, maxc + 1):
                    try:
                        id = str(
                            design.structure_helices_coord_map[(r, c)].id)
                    except KeyError:
                        id = ""

                    but = widgets.ToggleButton(
                        description=id, layout=layout_button)
                    row.append(but)
                helix_buttons.append(row)

        HBoxes = []
        for row in helix_buttons:
            HBoxes.append(widgets.HBox(row))

        buttons_h = widgets.VBox(HBoxes)

        slider_b = widgets.IntRangeSlider(
            value=[minb, maxb],
            min=minb,
            max=maxb,
            step=1,
            description='base:',
            disabled=False,
            continuous_update=True,
            orientation='horizontal',
            readout=True,
            readout_format='d',
            layout=layout_w
        )

        slider_c = widgets.IntSlider(
            value=4,
            min=1,
            max=10,
            step=1,
            description='context [Angstr.]:',
            disabled=False,
            continuous_update=True,
            orientation='vertical',
            readout=True,
            readout_format='d',
            layout=layout_h
        )

        sliders_w = widgets.VBox([buttons_h, slider_b])
        sliders = widgets.Box([sliders_w, slider_c], layout=layout_box)
        display(sliders)

        return (helix_buttons, slider_b, slider_c)

    def eval_sliders(self, helix_buttons, slider_b, slider_c):

        selection_bases = range(slider_b.lower, slider_b.upper + 1)

        selection_helices = []
        flat_list = [item for sublist in helix_buttons for item in sublist]
        for button in flat_list:
            if button.value:
                try:
                    selection_helices.append(int(button.description))
                except ValueError:
                    pass
        with mrcfile.open(self.project.files.mrc, mode='r+') as mrc:
            voxel_size = mrc.voxel_size.x
        context = int(slider_c.value / voxel_size) + 1

        return (selection_helices, selection_bases), context

    def writepdb(self, atoms_selection, single_frame=True, frame=-1):
        inp = self.project.input
        path_in = inp / "{}-tmp.pdb".format(self.project.name)
        path_corr = inp / "{}-tmp_chim.pdb".format(self.project.name)

        with mda.Writer(path_in, multiframe=True,
                        n_atoms=atoms_selection.n_atoms) as W:
            if single_frame:
                self.linker.fit.u.trajectory[frame]
                W.write(atoms_selection)
            else:
                for _ in self.linker.fit.u.trajectory:
                    W.write(atoms_selection)
        logic = Logic(keep_header=False,
                      remove_H=False,
                      )
        pdb_Corr = PDB_Corr(reverse=False)
        with open(path_in, "r") as file_init:
            newFile = pdb_Corr.correct_pdb(pdb_file=file_init,
                                           logic=logic,
                                           )
        with open(path_corr, "w") as file_corr:
            file_corr.write(newFile)

        return
