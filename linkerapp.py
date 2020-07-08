#!/usr/bin/env python
# -*- coding: utf-8 -*-

import ipywidgets as widgets
import MDAnalysis as mda

import mrcfile
from pathlib import Path
from IPython.display import display

from core.linker import Linker
from core.project import Project, Files
from core.pdbCorrection import PDB_Corr, Logic

from linkerappUtil import (_create_voxel_mask, _mrc_cut_minimalbox,
                           _get_mrc_properties)


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
    def __init__(self, folder, files: Files):
        name = files[0].stem
        self.project = Project(folder=Path(folder),
                               name=name,
                               files=files,
                               )
        self.linker = Linker(self.project)
        self.link = self.linker.create_linkage()
        self.Hid2H = self.linker.design.design.structure_helices_map
        self.Hcoor2H = self.linker.design.design.structure_helices_coord_map
        self.selection_scaffold = None
        self.selection_staples = None

    def select_by_helixandbase(self, helices, bases):

        self._parse_selection(bases, helices)
        atoms_scaffold = mda.AtomGroup([], self.link.u)
        atoms_staple = mda.AtomGroup([], self.link.u)

        for base in self.selection_scaffold:
            atoms_scaffold += self.link.u.residues[
                self.link.DidFid[
                    self.link.DhpsDid[(base.h, base.p, True)]
                ]
            ].atoms

        for base in self.selection_staples:
            atoms_staple += self.link.u.residues[
                self.link.DidFid[
                    self.link.DhpsDid[(base.h, base.p, False)]
                ]
            ].atoms

        atoms_selection = atoms_scaffold + atoms_staple

        return atoms_selection, self.link.Dcolor

    def _parse_selection(self, base_selection, helix_selection):
        # get helices from design
        helices = [self.Hid2H[idx] for idx in helix_selection]

        # get bases from helices
        self.selection_scaffold = []
        self.selection_staples = []
        # TODO: test without deletions
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

    def make_sliders(self):
        layout_w = widgets.Layout(width='900px', height='70px')
        layout_h = widgets.Layout(width='100px')
        layout_box = widgets.Layout(width='1000px', border='1px solid black')
        layout_button = widgets.Layout(width='30px', height='30px')

        base_p = [base.p for base in self.linker.design.allbases]
        minb, maxb = min(base_p), max(base_p)
        row = [coor[0] for coor in self.Hcoor2H.keys()]
        minr, maxr = min(row), max(row)
        col = [coor[1] for coor in self.Hcoor2H.keys()]
        minc, maxc = min(col), max(col)

        lattice = ("square" if self.linker.design.design.lattice_type == 0
                   else "honeycomb")

        helix_buttons = []
        if lattice in ["square", "honeycomb"]:  # TODO improve honeycomb
            for r in range(minr, maxr + 1):
                row = []
                for c in range(minc, maxc + 1):
                    h = self.Hcoor2H.get((r, c), None)
                    h_id = "" if h is None else str(h.id)
                    but = widgets.ToggleButton(
                        description=h_id,
                        layout=layout_button,
                        disabled=(h is None),
                    )
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

        flat_list = [item for sublist in helix_buttons for item in sublist]
        selection_helices = [int(b.description) for b in flat_list
                             if (not b.disabled and b.value)
                             ]

        with mrcfile.open(self.project.files.mrc, mode='r+') as mrc:
            voxel_size = mrc.voxel_size.x
        context = int(slider_c.value / voxel_size) + 1

        return (selection_helices, selection_bases), context

    def writepdb(self, atoms_selection, name, single_frame=True, frame=-1, chimeraX=True):
        inp = self.project.folder
        path_in = inp / "{}.pdb".format(name)
        path_corr = inp / "{}_chim.pdb".format(name)

        with mda.Writer(path_in, multiframe=True,
                        n_atoms=atoms_selection.n_atoms) as W:
            if single_frame:
                self.link.u.trajectory[frame]
                W.write(atoms_selection)
            else:
                for _ in self.link.u.trajectory:
                    W.write(atoms_selection)

        if chimeraX:
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

    def writemrc(self, atoms_selection, name: str, context=4, cut_box=True):
        if not len(atoms_selection):
            print("EXIT - no atoms in this selection")
            return

        inp = self.project.folder
        path_out = inp / "{}.mrc".format(name)
        u = atoms_selection.universe
        u.trajectory[-1]

        with mrcfile.open(self.project.files.mrc) as mrc:
            voxel_size, grid, origin, full_data = _get_mrc_properties(mrc)

        data_mask = _create_voxel_mask(
            atoms_selection, grid, origin, voxel_size, context)
        data = full_data * data_mask

        if cut_box:
            data, origin, voxel_size = _mrc_cut_minimalbox(data, origin, voxel_size)

        with mrcfile.new(path_out, overwrite=True) as mrc_out:
            mrc_out.set_data(data.transpose())
            mrc_out._set_voxel_size(*voxel_size)
            mrc_out.header["origin"] = tuple(origin)
