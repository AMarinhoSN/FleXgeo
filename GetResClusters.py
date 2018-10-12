#!/usr/bin/env python
# -*- coding: utf-8 -*-
import FXgeoToolBox
import argparse
from os import getcwd
#from sys import argv

'''
< GetResCluster >
This script is designed to clustering residues based on curvature and torsion
distribution using HDBSCAN.
'''

# define arguments
parser = argparse.ArgumentParser(description='''Cluster conformations based on
residues distribution of curvature and torsion using HDBSCAN.''')

parser.add_argument("in_csv", type=str, help='''csv file with FleXgeo data''')

parser.add_argument("-res",default='ALL', type=str, help='''specify residue to
be clustered. (default: 'ALL')''')

parser.add_argument("-out_path",default=getcwd(), type=str, help='''specify
dir to write output files. (default: working dir)''')

parser.add_argument("-min_pcluster", default=.05, type=float,
help= '''Set the minimum conformations percentage a cluster must have
(default: .05)''')

parser.parse_args()
parser.print_help()
args = parser.parse_args()
print("---")

# 1 Load CSV input data
Csv_filenm = args.in_csv
res = args.res
min_pcent_clstr = args.min_pcluster
out_path = args.out_path

print("@ Load CSV data of ", Csv_filenm)
CSVdata = FXgeoToolBox.FXGeoData(Csv_filenm)

# 2 Run HDBScan for a given residue

print("@ Run HDBScan")
# 2.1 - setting min size preference as at least 5% of total confs
MinSizeClstr = int(CSVdata.nconfs * min_pcent_clstr)
print(" > Min Size Cluster set as ", MinSizeClstr, " conformations")

if res != 'ALL':
    try:
        res = int(res)
        err_msg = '''residue {} is not present at ensemble input'''.format(res)
        assert(res in CSVdata.reslist), err_msg
        print("--> @ ", res)
        CSVdata.doHDBScanClstForRes(res, min_cluster_size=MinSizeClstr,
                                    lang="EN", savePDF=True, showplot=True,
                                    AllowSingle=False)
    except(ValueError):
        print('!!ERROR {} is not a valid integer.'.format(res))
        exit(1)
if res == 'ALL':
    for res_n in CSVdata.reslist:
        print("--> @ ", res_n)
        CSVdata.doHDBScanClstForRes(res_n, min_cluster_size=MinSizeClstr,
                                    showplot=False, lang="EN", out_dir=out_path)

print("@ Writing StrClstr.clst")
if out_path[-1] != '/':
    out_path = out_path+'/'
CSVdata.writeStrucClstFiles(clstflnm = out_path+"StrClstr.clst")

print("@ Writing summary on StrClstr_sum.csv...")
Sum_outFile = open(out_path+"StrClstr_sum.csv", "w")
print(" > Residues with more than a single cluster:")
for res in CSVdata.lclustersPerRes:
    Sum_outFile.write(str(res[0])+","+str(res[1].nclusters)+"\n")
    if res[1].nclusters >= 2:
        print(" {} | n_clusters = {}".format(res[0], res[1].nclusters))

Sum_outFile.close()
print(":: DONE ::")
