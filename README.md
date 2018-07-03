# FleXgeo

The FleXgeo is a software package designed for protein conformational ensemble analyses based on a differential geometry representation of protein backbones. The package is composed of a binary of the core program which calculates the differential geometry descriptors (more details bellow) and a set of scripts designed for the analyses of the results.

## How to run it?
There are two stages of FleXgeo applications, the differential geometry descriptors calculation and the analyses part. A website with tutorials will be provided in the future, but for now check the **quick and dirt** guide:
- **1 . Calculate Differential Geometry**
	>$ /path/to/FleXgeo_bin -pdb=ensemble.pdb [options]

FleXgeo accetps the following arguments:

|    OPTIONS       | DESCRIPTION               | DEFAULT                 |
|----------------|-----------------------------|--------------------------|
|`-pdb=[filename.pdb]`|Set input .pdb filename|User must provide  |
|`-ncpus=[int]`      |Set the number of cpus FleXgeo will use | All cpus available|
| `-isSingle`        |Indicate if input pdb is a single conformation pdb | FALSE|
|`-outprfx=[prefix]` | Set the output files prefix to be used | 'Diffgeo_' |

- **2. Analyses**
	* Plot xgeo data
	>$ python3.5 /path/to/PlotFXgeoData.py DiffGeo_raw.csv

	* Calculate distance between all conformations and a reference conformation on the ensemble
	>$ python3.5 /path/to/CalcEnsDistFromRef.py -in=XgeoObjFilname.p

	* [TO ADD] Calculate distance between all conformations and an external reference conformation
	* [TO ADD] Calculate dMax
	* [TO ADD] Clustering by HDBScan

## What are the contents of the output files?
FleXgeo outputs 5 .csvs files:

[ADICIONAR DESCRIÇAO]

## How to cite?
If you use results produced by FleXgeo on your amazing publication, please cite:
[adicionar paper]

## Why just binaries files for the core of FleXgeo?
Unfortunatly, we use Numerical Recipies on our code and we are not allowed to distribute the source code. We plan to rewrite those part of the code in the future, but this is not on our top priorities right now. If you have some trouble on running a binary file on your machine, feel free to contact Antonio at amarinho@cent.uw.edu.pl and we can try to provide an specific binary for you.
