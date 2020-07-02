#!/usr/bin/env python
# -*- coding: utf-8 -*-

import attr
from pathlib import Path
from typing import NamedTuple

__authors__ = ["Elija Feigl"]
""" VIEWERTOOL:
    Context Class to store user input and pass it between various scripts.
"""


@attr.s(slots=True)
class Project(object):

    input: Path = attr.ib()
    output: Path = attr.ib()
    name: str = attr.ib()
    files: NamedTuple = attr.ib()
