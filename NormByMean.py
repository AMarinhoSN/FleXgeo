import numpy as np
import pandas as pd
import argparse
import matplotlib.pyplot as plt
from os import getcwd

import seaborn as sns; sns.set(color_codes=True)

def normalization(array):
    '''normalization by keeping mean value '''
    x_average = np.mean(array)
    x_std_dev = np.std(array)
    new = [(x - x_average)/x_std_dev for x in array]
    return new

# define arguments
parser = argparse.ArgumentParser(description='''Generate normalized values of
FleXgeo descriptors by keeping the same mean values observed in the original
data.''')

parser.add_argument("in_csv", type=str, help='''csv file with FleXgeo data''')

parser.add_argument("-out_dir",default=getcwd(), type=str, help='''specify
dir to write output files. (default: working dir)''')

parser.add_argument("-out_sfx", default="DiffGeo_NORM_mean", type=str,
help= '''Suffix of csv output (default: DiffGeo_NORM_mean)''')

parser.parse_args()
parser.print_help()
args = parser.parse_args()

DIR = args.out_dir
CSV_IN = args.in_csv
OUT_SFX = args.out_sfx

print("@ loading {} data".format(CSV_IN))
df = pd.read_csv(CSV_IN,
                 sep=',',
                 decimal='.',
                 header=None,
                 names=['conf','res_num', 'curvature', 'torsion', 'arc', 'writhing',
                        'x', 'y', 'z', 'res'])

print("@ computing normalized values")
curv = df['curvature'].tolist()
tors = df['torsion'].tolist()
wris = df['writhing'].tolist()
arcs = df['arc'].tolist()
idx = df.index.tolist()

curvm = normalization(curv)
torsm = normalization(tors)
wrism = normalization(wris)
arcsm = normalization(arcs)

print("@ writing {}".format(OUT_SFX+".csv"))
dct = {'c_norm':curvm, 't_norm':torsm, 'w_norm':wrism, 'a_norm':arcsm}
cols_to_add_df = pd.DataFrame(dct)
new_df = pd.concat([df, cols_to_add_df], axis=1)
new_df.to_csv(DIR+OUT_SFX+".csv")
print(":: DONE ::")
