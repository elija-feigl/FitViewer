#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import nglview as nv

from operator import attrgetter

from MDAnalysis.core.groups import AtomGroup

__authors__ = ["Elija Feigl"]
__version__ = "1.0.0"
__license__ = "GPL-3.0"

__descr__ = "VIEWERTOOL: utility scripts for viewer-tool"
__status__ = "Development"
__maintainer__ = "Elija Feigl"
__email__ = "elija.feigl@tum.de"


def fit_widget(atoms_selection, mrc,
               opacity_map=0.9, isolevel=5., backbone=True, color_dict=None):

    if not len(atoms_selection):
        print("NOTHING TO DISPLAY - no atoms in this selection")
        return

    u = atoms_selection.universe
    u.trajectory[-1]
    scaffold = max(
        atoms_selection.segments, key=attrgetter("residues.n_residues")
    )
    atoms_scaffold = scaffold.atoms & atoms_selection
    atoms_staple = atoms_selection - atoms_scaffold

    if backbone:
        scaffold_view = atoms_scaffold.select_atoms(
            "name C1' C2' C3' C4' C5' O1P O2P O3' O4' O5' P")
        staple_view = atoms_staple.select_atoms(
            "name C1' C2' C3' C4' C5' O1P O2P O3' O4' O5' P")
    else:
        scaffold_view = atoms_scaffold
        staple_view = atoms_staple

    view = nv.show_mdanalysis(scaffold_view)

    view.clear(component=0)
    view.add_representation("spacefill", component=0,
                            radius=0.6, color="blue")

    view.add_component(str(mrc))
    view.clear(component=1)
    view.add_representation("surface", component=1, color='grey',
                            wireframe=True, opacity=opacity_map,
                            isolevel=isolevel)
    if color_dict is None:
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


def _create_voxel_mask(atoms, grid, origin, voxel_size, context):
    data_mask = np.zeros(grid, dtype=np.float32)
    v_context = np.full(3, context / voxel_size).astype(int) + 1
    v_context_x, v_context_y, v_context_z = v_context
    print(len(atoms))
    grid_positions = np.rint(
        ((atoms.positions - origin) / voxel_size)).astype(int)
    print(len(grid_positions))
    for pos in grid_positions:
        x, y, z = pos[0], pos[1], pos[2]  # fast to slow axis
        data_mask[x - v_context_x: x + v_context_x,
                  y - v_context_y: y + v_context_y,
                  z - v_context_z: z + v_context_z
                  ] = 1.
    return data_mask


def _get_mrc_properties(mrc):
    o = np.array(mrc.header["origin"])
    origin = np.array([o["x"], o["y"], o["z"]])
    c = np.array(mrc.header["cella"])
    cell = np.array([c["x"], c["y"], c["z"]])
    grid = np.array(
        [mrc.header["nx"], mrc.header["ny"], mrc.header["nz"]])
    data = mrc.data.transpose()
    voxel_size = cell / grid

    return (voxel_size, grid, origin, data)


def _mrc_cutbox(data, m_origin, m_spacing):
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

    # cumpute new origin
    origin = (m_origin + ((x_min - x_low) * m_spacing[0],
                          (y_min - y_low) * m_spacing[1],
                          (z_min - z_low) * m_spacing[2]
                          )
              )
    grid = np.shape(data_small)
    cell = grid * m_spacing
    spacing = (cell / grid)
    return data_small, origin, spacing
