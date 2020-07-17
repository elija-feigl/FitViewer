#!/usr/bin/env python
# -*- coding: utf-8 -*-

import ipywidgets as widgets
import MDAnalysis as mda

import mrcfile
from pathlib import Path
from typing import List
from IPython.display import display

from core.linker import Linker
from core.project import Project, Files
from core.pdbCorrection import PDB_Corr, Logic

from linkerappUtil import _create_voxel_mask, _mrc_cutbox, _get_mrc_properties


__authors__ = ["Elija Feigl"]
__version__ = "0.4"
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
    def __init__(self, folder: str, files: Files):
        name = files[0].stem
        self.project = Project(folder=Path(folder),
                               name=name,
                               files=files,
                               )
        self.linker = Linker(self.project)
        self.link = self.linker.create_linkage()
        self.Hid2H = self.linker.design.design.structure_helices_map
        self.Hcoor2H = self.linker.design.design.structure_helices_coord_map

    def select_by_helixandbase(self, helices: List, bases: List):
        def _DhpsFid(h, p, s) -> int:
            return self.link.DidFid[
                self.link.DhpsDid[(h, p, s)]
            ]
        u = self.link.u
        stsc_bases = self._parse_selection(bases, helices)
        atoms = mda.AtomGroup([], u)
        for idx, bases in enumerate(stsc_bases):
            for base in bases:
                Fid = _DhpsFid(base.h, base.p, bool(idx))
                atoms += u.residues[Fid].atoms
        return atoms, self.link.Dcolor

    def _parse_selection(self, base_pos: List, helix_ids: List) -> List:
        helices = [self.Hid2H[idx] for idx in helix_ids]
        scaffold = [
            b for h in helices for b in h.scaffold_bases if b.p in base_pos
        ]
        staples = [
            b for h in helices for b in h.staple_bases if b.p in base_pos
        ]
        return [staples, scaffold]

    def select_widget(self):
        def _button(r, c, lattice):
            h = self.Hcoor2H.get((r, c), None)
            h_id = "" if h is None else str(h.id)
            if lattice in ["square", "honeycomb"]:
                return widgets.ToggleButton(
                    description=h_id,
                    layout=widgets.Layout(width='30px', height='30px'),
                    disabled=(h is None),
                )
            else:
                raise TypeError

        def _minmax(alist: List, asrange=False):
            if asrange:
                return range(min(alist), max(alist) + 1)
            else:
                return min(alist), max(alist)

        minb, maxb = _minmax([base.p for base in self.linker.design.allbases])
        r_range = _minmax([coor[0] for coor in self.Hcoor2H.keys()], True)
        c_range = _minmax([coor[1] for coor in self.Hcoor2H.keys()], True)

        lattice = ("square" if self.linker.design.design.lattice_type == 0
                   else "honeycomb")
        helixButtons = [
            [_button(r, c, lattice) for c in c_range] for r in r_range
        ]

        buttonsBox = widgets.VBox([widgets.HBox(row) for row in helixButtons])
        baseSlider = widgets.IntRangeSlider(
            value=[minb, maxb], min=minb, max=maxb, step=1,
            description='base:',
            disabled=False, continuous_update=True, orientation='horizontal',
            readout=True, readout_format='d',
            layout=widgets.Layout(width='900px', height='70px')
        )
        contextSlider = widgets.IntSlider(
            value=4, min=1, max=10, step=1,
            description='context [Angstr.]:',
            disabled=False, continuous_update=True, orientation='vertical',
            readout=True, readout_format='d',
            layout=widgets.Layout(width='100px')
        )

        vBox = widgets.VBox([buttonsBox, baseSlider])
        mainWidget = widgets.Box(
            [vBox, contextSlider],
            layout=widgets.Layout(width='1000px', border='1px solid black')
        )
        display(mainWidget)
        return (helixButtons, baseSlider, contextSlider)

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

    def writepdb(self, atoms, name, singleframe=True, frame=-1, chimeraX=True):
        import warnings
        warnings.filterwarnings('ignore')

        inp = self.project.folder
        path_in = inp / "{}.pdb".format(name)
        path_corr = inp / "{}_chim.pdb".format(name)

        with mda.Writer(path_in, multiframe=True,
                        n_atoms=atoms.n_atoms) as W:
            if singleframe:
                self.link.u.trajectory[frame]
                W.write(atoms)
            else:
                for _ in self.link.u.trajectory:
                    W.write(atoms)

        if chimeraX:
            logic = Logic(keep_header=False,
                          remove_H=False,
                          )
            pdb_Corr = PDB_Corr()
            with open(path_in, "r") as file_init:
                newFile = pdb_Corr.correct_pdb(pdb_file=file_init,
                                               logic=logic,
                                               )
            with open(path_corr, "w") as file_corr:
                file_corr.write(newFile)

    def writemrc(self, atomsXX, name: str, context=4, cut_box=True):
        if not len(atomsXX):
            raise ValueError("no atoms in this selection")
        path_out = self.project.folder / "{}.mrc".format(name)
        u = atomsXX.universe
        u.trajectory[-1]
        with mrcfile.open(self.project.files.mrc) as mrc:
            voxel, grid, origin, full_data = _get_mrc_properties(mrc)

        data_mask = _create_voxel_mask(atomsXX, grid, origin, voxel, context)
        data = full_data * data_mask
        if cut_box:
            data, origin, voxel = _mrc_cutbox(data, origin, voxel)

        with mrcfile.new(path_out, overwrite=True) as mrc_out:
            mrc_out.set_data(data.transpose())
            mrc_out._set_voxel_size(*voxel)
            mrc_out.header["origin"] = tuple(origin)
