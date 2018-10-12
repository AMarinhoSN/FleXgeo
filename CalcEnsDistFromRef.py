#!/usr/bin/env python
# -*- coding: utf-8 -*-

import FXgeoToolBox
from sys import argv
import pickle
import textwrap
'''
This scripts calculate the euclidean distance in the curvature and torsion
space between all conformations of the ensemble and an arbitrary reference
conformation. The script also generates a heatmap plot for visualization.

The default behaviour is calculate against the first conformation on the
ensemble, but the user can specify another by changing -cref_idx.

If you want to use an external conformation as a reference (such as a
crystallographic structure), just provide the .csv output from FleXgeo for the
structure you need using the -ext_ref.
'''

description = "This scripts calculate the euclidean distance in the curvature \
and torsion space between all conformations of the ensemble and an arbitrary \
reference conformation. The script also generates a heatmap plot for \
visualization. \n\n The default behaviour is use the first conformation on the ensemble as the \
reference, but the user can specify another conformation via '-cref_idx'. \
\n\n If you want to use an external conformation as a reference (such as a \
crystallographic structure), just provide the .csv output from FleXgeo for the \
structure you need via '-ext_ref'."

wraptxt = textwrap.fill(description, width=69, initial_indent='    ',
                        subsequent_indent = '  ' )

final_desc = textwrap.indent(wraptxt, '|')
print("|---------------------------------------------------------------------|")
print("|          >> CALCulate ENSemble DISTance From a Reference << ")
print(final_desc)
print("|---------------------------------------------------------------------|")
print(": USAGE: python CalcEnsDistFromRef.py  -in=XgeoObjFilname.p [OPTIONS]")
print("|---------------------------------------------------------------------|")
print("| [OPTIONS]")
print("|     -refidx=[int]        : the reference conformation index")
print("|                            (default=0)")
print("|     -ext_ref=[_xgeo.csv] : external reference state xgeo csv output")
print("|---------------------------------------------------------------------|")


print(": INPUT: ")

# 0 - Process input
pkl_flnm = None
refidx = None
ext_ref = None

for arg in argv:
    if arg == argv[0]:
        continue

    if arg.startswith('-in='):
        pkl_flnm = arg[4:len(arg)]
        print(":  -in = {}".format(pkl_flnm))

    if arg.startswith("-refidx="):
        try:
            refidx = int(arg[8:len(arg)])
        except(ValueError):
            print("! ERROR: You must provide a valid integer.")
            exit(1)
        print(":  -refidx = {}".format(refidx))

    if arg.startswith("-ext_ref="):
        ext_ref = arg[9:len(arg)]
        print(":  -ext_ref = {}".format(ext_ref))

# default behaviour is assume the first internal conf will be used
if refidx == None and ext_ref == None:
    refidx = 0

if refidx != None and ext_ref != None:
    print("! ERROR: This script only use one conformation as a reference. Choose")
    print("!       wisely.")
    exit(1)
print("|---------------------------------------------------------------------|")
# 1 Load and calculate from reference conf
print("@ Loading ensemble data from {}".format(pkl_flnm))
XgeoDataObj = pickle.load(open(pkl_flnm, "rb"))

if refidx != None:
    print("@ Calculate KT euclidean distance from reference conformation")
    XgeoDataObj.calcD2ConfperRes(cref_idx=refidx, showplot=False)
    print(":: DONE ::")

if ext_ref != None:
    print("@ loading data from {}".format(ext_ref))
    Ext_Refx = FXgeoToolBox.FXGeoData(ext_ref)
    print("@ Calculate KT euclidean distance from external reference conformation")
    XgeoDataObj.calcD2ConfperRes(ext_ref=Ext_Refx, showplot=False)
    print(":: DONE ::")
