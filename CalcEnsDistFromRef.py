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
print(" USAGE: python CalcEnsDistFromRef.py  ")

print("@ Load data")
pkl_flnm = argv[1]
XgeoDataObj = pickle.load(open(pkl_flnm, "rb"))

print("@ Calculate KT euclidean distance from reference conformation")
XgeoDataObj.calcD2ConfperRes(cref_idx=0, showplot=False)
print(":: DONE ::")
