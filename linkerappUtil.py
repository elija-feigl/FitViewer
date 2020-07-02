#!/usr/bin/env python
# -*- coding: utf-8 -*-

import mrcfile
import numpy as np
import ipywidgets as widgets
import nglview as nv
import os

from PyQt5.QtWidgets import QFileDialog
from pathlib import Path
from collections import namedtuple


__authors__ = ["Elija Feigl"]
__version__ = "1.0"
__license__ = "GPL-3.0"

__descr__ = "VIEWERTOOL: utility scripts for viewer-tool"
__status__ = "Development"
__maintainer__ = "Elija Feigl"
__email__ = "elija.feigl@tum.de"


Files = namedtuple("Files", ["json", "psf", "coor", "mrc", "seq"])


def mrc(atoms_selection, path, context=4, cut_box=True):

    if not len(atoms_selection):
        print("EXIT - no atoms in this selection")
        return

    u = atoms_selection.universe
    u.trajectory[-1]
    with mrcfile.open(str(path) + ".mrc", mode='r+') as mrc:
        m_o = np.array(mrc.header["origin"])
        m_origin = np.array([m_o["x"], m_o["y"], m_o["z"]])
        m_c = np.array(mrc.header["cella"])
        m_cell = np.array([m_c["x"], m_c["y"], m_c["z"]])
        m_grid = np.array(
            [mrc.header["nx"], mrc.header["ny"], mrc.header["nz"]])
        m_spacing = m_cell / m_grid
        m_data = mrc.data.transpose()

    data_mask = np.zeros(m_grid, dtype=np.float32)
    voxel_size = m_cell / m_grid
    v_context = np.full(3, context / voxel_size).astype(int) + 1
    v_context_x, v_context_y, v_context_z = v_context

    grid_positions = np.rint(
        ((atoms_selection.positions - m_origin) / m_spacing)).astype(int)
    for pos in grid_positions:
        x, y, z = pos[0], pos[1], pos[2]  # fast to slow axis
        data_mask[x - v_context_x: x + v_context_x,
                  y - v_context_y: y + v_context_y,
                  z - v_context_z: z + v_context_z
                  ] = 1.
    data = m_data * data_mask

    if cut_box:
        # get rid of zero-padding
        idx_data = np.nonzero(data)

        x_min, x_max = np.min(idx_data[0]), np.max(idx_data[0])
        y_min, y_max = np.min(idx_data[1]), np.max(idx_data[1])
        z_min, z_max = np.min(idx_data[2]), np.max(idx_data[2])

        xyz_diff = max(x_max - x_min, y_max - y_min, z_max - z_min)
        x_pad = 0.5 * (xyz_diff + x_min - x_max)
        y_pad = 0.5 * (xyz_diff + y_min - y_max)
        z_pad = 0.5 * (xyz_diff + z_min - z_max)
        x_low = int(x_pad) if (x_pad % 1.) == 0. else int(x_pad) + 1
        y_low = int(y_pad) if (y_pad % 1.) == 0. else int(y_pad) + 1
        z_low = int(z_pad) if (z_pad % 1.) == 0. else int(z_pad) + 1

        data_small = np.zeros(
            (xyz_diff, xyz_diff, xyz_diff), dtype=np.float32)
        data_small[x_low: -int(x_pad) or None, y_low: -int(y_pad) or None,
                   z_low: -int(z_pad) or None] = data[x_min: x_max,
                                                      y_min: y_max,
                                                      z_min: z_max]

        grid = np.shape(data_small)
        origin = (m_origin + ((x_min - x_low) * m_spacing[0], (y_min - y_low) *
                              m_spacing[1], (z_min - z_low) * m_spacing[2]))
        cell = grid * m_spacing
        spacing = (cell / grid)
    else:
        data_small = data
        spacing = m_spacing
        origin = m_origin

    with mrcfile.new(str(path) + "-tmp.mrc", overwrite=True) as mrc_out:
        mrc_out.set_data(data_small.transpose())
        mrc_out._set_voxel_size(*spacing)
        mrc_out.header["origin"] = tuple(origin)

    return


def fit_widget(atoms_selection, atoms_scaffold, atoms_staple, path,
               opacity_map=0.9, isolevel=5., backbone=True, color_dict=None):

    if not len(atoms_selection):
        print("NOTHING TO DISPLAY - no atoms in this selection")
        return

    if backbone:
        scaffold_view = atoms_scaffold.select_atoms(
            "name C1' C2' C3' C4' C5' O1P O2P O3' O4' O5' P")
        staple_view = atoms_staple.select_atoms(
            "name C1' C2' C3' C4' C5' O1P O2P O3' O4' O5' P")
    else:
        scaffold_view = atoms_scaffold
        staple_view = atoms_staple

    # TODO: -low only last frame
    view = nv.show_mdanalysis(scaffold_view)

    view.clear(component=0)
    view.add_representation("spacefill", component=0,
                            radius=0.6, color="blue")

    view.add_component(str(path) + "-tmp.mrc")
    view.clear(component=1)
    view.add_representation("surface", component=1, color='grey',
                            wireframe=True, opacity=opacity_map,
                            isolevel=isolevel)
    if color_dict is None:  # no staple coloring
        view.add_trajectory(staple_view)
        view.clear(component=2)
        view.add_representation("spacefill", component=2,
                                radius=0.6, color="red")
    else:
        segnames = [seg.segid for seg in staple_view.segments]
        segindices = [seg.segindex for seg in staple_view.segments]
        for i, (segid, segindex) in enumerate(zip(segnames, segindices)):
            s_view = staple_view.select_atoms("segid " + str(segid))
            view.add_trajectory(s_view)
            view.clear(component=(i + 2))
            view.add_representation("spacefill", component=(
                i + 2), radius=0.6, color=color_dict[segindex])

    return view


class FileBrowser(object):
    """
    https://gist.github.com/DrDub/6efba6e522302e43d055
    """

    def __init__(self):
        self.path = os.getcwd()
        self._update_files()

    def _update_files(self):
        self.files = list()
        self.dirs = list()
        if(os.path.isdir(self.path)):
            for f in os.listdir(self.path):
                ff = self.path + "/" + f
                if os.path.isdir(ff):
                    self.dirs.append(f)
                else:
                    self.files.append(f)

    def widget(self):
        box = widgets.VBox()
        self._update(box)
        return box

    def _update(self, box):

        layout_button = widgets.Layout(width='800px', height='30px',
                                       text_align="start",
                                       border="1px solid black")

        def on_click(b):
            if b.description == '..':
                self.path = os.path.split(self.path)[0]
            else:
                self.path = self.path + "/" + b.description
            self._update_files()
            self._update(box)

        buttons = []
        if self.files:
            button = widgets.Button(description='..',
                                    background_color='#d0d0ff',
                                    layout=layout_button)
            button.on_click(on_click)
            buttons.append(button)
        for f in self.dirs:
            button = widgets.Button(description=f,
                                    background_color='#d0d0ff',
                                    layout=layout_button)
            button.on_click(on_click)
            buttons.append(button)

        box.children = tuple(
            [widgets.HTML("<h3>%s</h3>" % (self.path,))] + buttons)
        return


def gui_fileDialog(wd='./', filter="all files (*)"):
    """Select a file via a dialog and return the file name."""
    fnames = QFileDialog.getOpenFileName(None, "Select file...", wd, filter=filter)
    return Path(fnames[0])


def gui_folderDialog(wd='./'):
    """Select a file via a dialog and return the file name."""
    return QFileDialog.getExistingDirectory(None, "Choose folder...", wd)
