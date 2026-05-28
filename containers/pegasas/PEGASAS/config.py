# -*- coding: UTF-8 -*-
"""General configuration for the PEGASAS package."""

from pkg_resources import resource_filename
import os
import sys

CURRENT_VERSION = "v1.1.1-py3"


def update_progress(progress):
    barLength = 20
    status = ""
    if isinstance(progress, int):
        progress = float(progress)
    if not isinstance(progress, float):
        progress = 0
        status = "error: progress var must be float\r\n"
    if progress < 0:
        progress = 0
        status = "Halt...\r\n"
    if progress >= 1:
        progress = 1
        status = "Done...\r\n"
    block = int(round(barLength * progress))
    text = "\rPercent: [{0}] {1:.1f}% {2}".format(
        "#" * block + "-" * (barLength - block), progress * 100, status
    )
    sys.stdout.write(text)
    sys.stdout.flush()


def file_len(fin):
    with open(fin) as f:
        for i, _ in enumerate(f):
            pass
    return i + 1


HALLMARKS50   = resource_filename("PEGASAS.data", "hallmarks50.gmt.txt")
MAT_REORDER   = resource_filename("PEGASAS", "prepareGeneMatrixOrdered.py")
MAT_GENERATE  = resource_filename("PEGASAS", "generateMatrixbySample.py")
MAT_CORR      = resource_filename("PEGASAS", "cor_matrix_direct_perm.R")
GO_PLOT       = resource_filename("PEGASAS", "GO_plot.R")
GO_PLOT_LIB   = resource_filename("PEGASAS", "GO_enrichr_plot.R")
