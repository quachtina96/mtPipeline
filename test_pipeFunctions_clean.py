import shlex
import subprocess as sp
import os
import numpy
import sys
import pipeFunctions as pf

path = "/gpfs/home/quacht/debug/ID18_Mother"
depths = []

pf.clean(path)