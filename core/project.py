#!/usr/bin/env python
# -*- coding: utf-8 -*-

import attr
from pathlib import Path
from collections import namedtuple


__authors__ = ["Elija Feigl"]
""" VIEWERTOOL:
    Context Class to store user input and pass it between various scripts.
"""


Files = namedtuple("Files", ["json", "psf", "coor", "mrc", "seq"])


@attr.s(slots=True)
class Project(object):

    folder: Path = attr.ib()
    name: str = attr.ib()
    files: "Files" = attr.ib()
    mrdna: bool = attr.ib()
