#!/usr/bin/env python
# -*- coding: utf-8 -*-
import attr
from typing import Dict, Tuple, Set

from core.project import Project
from core.fit import Fit
from core.design import Design
from core.linkage import Linkage


__authors__ = ["Elija Feigl"]
""" VIEWERTOOL:
    Linker class, responsible for cross-referencing base IDs of cadnano design
    and atomic model
"""


@attr.s
class Linker(object):
    project: Project = attr.ib()
    Fbp: Dict[int, int] = dict()
    DidFid: Dict[int, int] = dict()
    DhpsDid: Dict[Tuple[int, int, bool], int] = dict()
    Fnicks: Dict[int, int] = dict()

    def __attrs_post_init__(self) -> None:
        self.fit: Fit = Fit(self.project)
        self.design: Design = Design(self.project)
        self.Dhp_skips: Set[Tuple[int, int]] = self.design.Dhp_skips

    def create_linkage(self) -> Linkage:
        """ invoke _link_scaffold, _link_staples, _link_bp to compute mapping
            of every base design-id to fit-id as well as the basepair mapping.
            basepairs are mapped from scaffold to staple, unique (invertable).
            updates linker attributes corresponding to the respective mapping
            and returns them.
        """
        self._link()
        self._identify_bp()
        self.link = Linkage(
            Fbp=self.Fbp,
            DidFid=self.DidFid,
            DhpsDid=self.DhpsDid,
            Dcolor=self.Dcolor,
            u=self.fit.u,
            Dhp_skips=self.Dhp_skips
        )
        return self.link

    def _link(self) -> Tuple[Dict[int, int],
                             Dict[Tuple[int, int, bool], int],
                             ]:
        def link_scaffold() -> Tuple[Dict[int, int],
                                     Dict[Tuple[int, int, bool], int],
                                     ]:
            """ collect position in scaffold (0-x) by comparing index in list
                of scaffold_design positions
            -------
                Returns
                -------
                DidFid
                    design-id -> fit-id
                DhpsDid
                    helix-number, base-position, is_scaffold -> design-id
            """
            Dscaffold = self.design.scaffold
            Did = [base.id for base in Dscaffold]
            Dhp = [(base.h, base.p, True) for base in Dscaffold]
            Fid_local = [Did.index(base.id) for base in Dscaffold]
            Fid_global = self.fit.scaffold.residues[Fid_local].resindices

            DidFid = dict(zip(Did, Fid_global))
            DhpsDid = dict(zip(Dhp, Did))
            return (DidFid, DhpsDid)

        def link_staples() -> Tuple[Dict[int, int],
                                    Dict[Tuple[int, int, bool], int],
                                    Dict[int, int],
                                    ]:
            """same procedure as scaffold for each
            -------
            Returns
                -------
                DidFid
                    design-id -> fit-id
                DhpsDid
                    helix-number, base-position, is_scaffold -> design-id
                dict color
                    fit-segment-id -> color
            """
            def get_resid(segindex: int, resindex_local: int) -> int:
                segment = self.fit.staples[segindex]
                return segment.residues[resindex_local].resindex

            DidFid: Dict[int, int] = {}
            DhpsDid: Dict[Tuple[int, int, bool], int] = {}
            color: Dict[int, int] = {}

            for i, staple in enumerate(self.design.staples):
                seg_id = self.design.stapleorder[i]

                Did = [base.id for base in staple]
                Dhp = [(base.h, base.p, False) for base in staple]

                Fid_local = [Did.index(base.id)for base in staple]
                Fid_global = [get_resid(seg_id, resid) for resid in Fid_local]

                icolor = self.design.design.strands[staple[0].strand].icolor
                segidxforcolor = self.fit.staples[seg_id].segindex
                color[segidxforcolor] = icolor

                DidFid_add = dict(zip(Did, Fid_global))
                DhpsDid_add = dict(zip(Dhp, Did))
                DidFid = {**DidFid, **DidFid_add}
                DhpsDid = {**DhpsDid, **DhpsDid_add}

            return (DidFid, DhpsDid, color)

        DidFid_sc, DhpsDid_sc = link_scaffold()
        DidFid_st, DhpsDid_st, self.Dcolor = link_staples()

        self.DidFid = {**DidFid_sc, **DidFid_st}
        self.DhpsDid = {**DhpsDid_sc, **DhpsDid_st}
        return (self.DidFid, self.DhpsDid)

    def _identify_bp(self) -> Dict[int, int]:
        """ link basepairs by mapping indices according to json (cadnano).
            basepairs are mapped from scaffold to staple, unique (invertable).
        -------
         Returns
            -------
            self.Fbp
                fit-id -> fit-id
        """
        self.Fbp = {
            self.DidFid[base.id]: self.DidFid[base.across.id]
            for base in self.design.scaffold
            if base.across is not None
        }
        return self.Fbp
