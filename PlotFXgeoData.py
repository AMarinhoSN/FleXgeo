import FXgeoToolBox
from sys import argv
'''
This scripts generate plots for FleXgeo results,
'''

print(" >> Plot Xgeo Data << ")
print(" USAGE: python PlotXgeoData.py DiffGeo_raw.csv ")

# 1 Load CSV input data
Csv_filenm = argv[1]
print("@ Load CSV data of ", Csv_filenm)
CSVdata = FXgeoToolBox.FXGeoData(Csv_filenm)
# 2 - Plot curvature and
print("@ Plot conformations curvature and torsion per residue")
CSVdata.plotConfsKTData(showplot=False, REF=None, REFlabel=None, lang="EN",
                        savePDF=False)
