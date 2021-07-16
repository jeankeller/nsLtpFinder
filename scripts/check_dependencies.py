#!/usr/bin/python3
# -*- coding: utf-8 -*-
# @Time    : 12/07/2021 17:54
# @Author  : Jean Keller; University of Toulouse III
# @Email   : kellerjeanphd@gmail.com
# @File    : check_dependencies.py
# @Software: PyCharm

import sys
import subprocess
import pkg_resources


def check_installed_modules():
    """
    Function to test if all Python modules are correctly installed and accessible
    :return: nothing, raise error if module not found
    """
    installed_modules = {pkg.key for pkg in pkg_resources.working_set}
    for module in ["pandas", "biopython", "argparse"]:
        sys.stdout.write(f"Testing if {module} is installed and accessible...  ")
        sys.stdout.flush()
        if module in installed_modules:
            sys.stdout.write("ok\n")
            sys.stdout.flush()
        else:
            raise ModuleNotFoundError(f"{module} is not installed or not accessible!")


def check_installed_programs():
    """
    Function to test if external programs are correctly installed and accessible
    By default, assumes that programs are in the path (global or user-specific)
    :return: nothing, raise error if program not found
    """
    programs = ["hmmsearch", "signalp5", "meme"]
    for p in programs:
        sys.stdout.write(f"Testing if {p} is installed and accessible...  ")
        sys.stdout.flush()
        test_prog_cl = subprocess.call(["which", p])
        if test_prog_cl == 0:
            pass
        else:
            raise ResourceWarning(f"Program {p} not found in path")
            # raise FileNotFoundError(f"Program {p} not found in path")
