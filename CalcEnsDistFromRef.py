import FXgeoToolBox
from sys import argv
import pickle

'''
This scripts calculate the euclidean distance in the curvature and torsion
space between all conformations of the ensemble and an arbitrary reference
conformation. The script also generates a heatmap plot for visualization.

The default behaviour is calculate against the first conformation on the
ensemble, but the user can specify another by changing -cref_idx.
'''

print(" >> CALCulate ENSemble DISTance From a Reference << ")
print(" USAGE: python CalcEnsDistFromRef.py  -in=XgeoObjFilname.p [OPTIONS]")
print(" [OPTIONS]")
print("        -refidx=[int] : the reference conformation index (default=0)")
print(" INPUT: ")
pkl_flnm = None
refidx = 0
for arg in argv:
    if arg == argv[0]:
        continue
    if arg.startswith('-in='):
        pkl_flnm = arg[4:len(arg)]
        print("   -in = {}".format(pkl_flnm))
    if arg.startswith("-refidx="):
        try:
            refidx = int(arg[8:len(arg)])
        except(ValueError):
            print("ERROR: You must provide a valid integer.")
            exit(1)
        print("   -refidx = {}".format(refidx))

print("@ Load data")
XgeoDataObj = pickle.load(open(pkl_flnm, "rb"))

print("@ Calculate KT euclidean distance from reference conformation")
XgeoDataObj.calcD2ConfperRes(cref_idx=refidx, showplot=False)
print(":: DONE ::")
