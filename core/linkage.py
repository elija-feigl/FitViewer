#!/usr/bin/env python
# -*- coding: utf-8 -*-
import pickle
import attr

from typing import Dict, Tuple, List, Set
import MDAnalysis as mda

from project import Project
from crossover import Crossover, CrossoverPicklable
from basepair import BasePair

__authors__ = ["Elija Feigl"]
""" VIEWERTOOL:
    Linkage class stores the translation from cadnano base-indexing to
    namd-indexing of bases.

    COMMENTS:
"""


@attr.s(auto_attribs=True)
class Linkage(object):
    Fbp: Dict[int, int] = {}
    DidFid: Dict[int, int] = {}
    DhpsDid: Dict[Tuple[int, int, bool], int] = {}
    Dcolor: Dict[int, int] = {}
    Dhp_skips: Set[Tuple[int, int]] = set()
    u: "mda.universe" = None

    def __attrs_post_init__(self) -> None:
        self._reverse()

    def _reverse(self) -> None:
        def reverse_d(dict: dict) -> dict:
            return {v: k for k, v in iter(dict.items())}

        self.FidDid = reverse_d(self.DidFid)
        self.DidDhps = reverse_d(self.DhpsDid)
        self.Fbp_rev = reverse_d(self.Fbp)
        self.Fbp_full = {**self.Fbp, **self.Fbp_rev}
