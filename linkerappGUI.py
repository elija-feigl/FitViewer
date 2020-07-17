#!/usr/bin/env python
# -*- coding: utf-8 -*-

from PyQt5.QtWidgets import QFileDialog
from pathlib import Path


def fileDialog(wd='./', filter="all files (*)"):
    """Select a file via a dialog and return the file name."""
    fnames = QFileDialog.getOpenFileName(
        None, "Select file...", wd, filter=filter
    )
    return Path(fnames[0])


def folderDialog(wd='./'):
    """Select a file via a dialog and return the file name."""
    return QFileDialog.getExistingDirectory(None, "Choose folder...", wd)
