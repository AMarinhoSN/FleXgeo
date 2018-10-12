import FXgeoToolBox
from sys import argv
import pickle
import argparse

'''
This scripts generate plots for FleXgeo results,
'''
# Define arguments
parser = argparse.ArgumentParser(description='''This scripts generate plots for
FleXgeo results.''')

parser.add_argument("in_xgeo", type=str, help='''"_xgeo.csv" fleXgeo data file.''')
parser.parse_args()
parser.print_help()
args = parser.parse_args()
print("------")

# 1 Load CSV input data
Csv_filenm = args.in_xgeo
print("@ Load CSV data of ", Csv_filenm)
CSVdata = FXgeoToolBox.FXGeoData(Csv_filenm)

print("@ Save CSVdata as pickle file 'XgeoObj.p'")
# \--\ save state \--\
pickle.dump(CSVdata, open("XgeoObj.p", "wb"))
# \-------------------\

# 2 - Plot
print("@ Plot conformations curvature and torsion per residue")
CSVdata.plotConfsKTData(showplot=False, REF=None, REFlabel=None, lang="EN",
                        savePDF=False)
print("@ Generate violin plots of FleXgeo values")
CSVdata.plotDataViolin(showplot=False, scale="count", inner="stick",
                       linewidth=.1, bw=.1, lang="EN")
