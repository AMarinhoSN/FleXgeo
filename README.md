# FleXgeo

The FleXgeo is a software package designed for protein conformational ensemble analyses based on a differential geometry representation of protein backbones. The package is composed of a binary of the core program which calculates the differential geometry descriptors (more details bellow) and a set of python scripts designed for the analyses of the results. Currently,  there are ready to use scripts to  :
* **cluster protein conformations**, via global clustering solution per residue based on its curvature and torsion distribution.
*  **quantify protein residues flexibility**, via  the computation of dmax
*  **compare protein conformations to a reference structure**,  via the computation of euclidean distances on the curvature and torsion space.

An object oriented solution to  is also provided for users that need a better modularity and include FleXgeo data on their own analyses pipeline.  

FleXgeo code was written by PhD. Antonio Marinho da Silva Neto and PhD. Rinaldo Wander Montalvao.

## INSTALLATION

- **1 . Download FleXgeo from GitHub**
	```bash
	$ git clone https://github.com/AMarinhoSN/FleXgeo.git
	```
- **2 . Install requirements**
	You can will a need a [Python 3](https://www.python.org/download/releases/3.0/)  and [pip3](https://pip.pypa.io/en/stable/installing/) to install the python libraries. You can install directly on your default python enviroment
	```bash
	$ pip3 install cython
	$ pip3 install -r /path/to/FleXgeo/requirements.txt
	```
	Or you can create an virtual enviroment just for FleXgeo by
	```bash
	$ pip3 install virtualenv
	$ cd virtual_envs_location
	$ virtualenv flexgeo_env
	$ source flexgeo_env/bin/activate
	$ cd FLEXGEO_LOCATION
	$ pip3 install -r requirements.txt
	```
  You also gonna need [Lua](https://www.lua.org/start.html) interpreter to run the part of the dmax analyses. On Linux, you can install it by:
	```bash
	$ sudo apt-get install luajit
	```


## How to run it?
There are two stages of FleXgeo applications, the differential geometry descriptors calculation and the analyses part. A website with tutorials will be provided in the future, but for now check the **quick and dirt** guide:
- **1 . Calculate Differential Geometry**
	```bash
	$ /path/to/FleXgeo/bin/FleXgeo_[MY_OS] -pdb=ensemble.pdb [options]
	```
FleXgeo accepts the following arguments:

|    OPTIONS       | DESCRIPTION               | DEFAULT                 |
|----------------|-----------------------------|--------------------------|
|`-pdb=[filename.pdb]`|Set input .pdb filename|User must provide  |
|`-ncpus=[int]`      |Set the number of cpus FleXgeo will use | All cpus available|
| `-isSingle`        |Indicate if input pdb is a single conformation pdb | FALSE|
|`-outprfx=[prefix]` | Set the output files prefix to be used | 'Diffgeo_' |

 - **2. Analyses**
	* Plot xgeo data
	```bash
	$ python3 /path/to/PlotFXgeoData.py DiffGeo_xgeo.csv
	```
			usage: PlotFXgeoData.py [-h] in_xgeo
			This scripts generate plots for FleXgeo results.
			positional arguments:
			  in_xgeo   'xgeo.csv' FleXgeo data file.
			optional arguments:
			  -h, --help  show this help message and exit

	* Calculate distance between all conformations and a reference conformation on the ensemble
	```bash
	$ python3.5 /path/to/CalcEnsDistFromRef.py -in=XgeoObjFilname.p
	```

	* Calculate distance between all conformations and an external reference conformation
	```bash
	$ python3.5 /path/to/CalcEnsDistFromRef.py -in=XgeoObjFilname.p -ext_ref=/path/to/ref_xgeo.csv
	```
	* Calculate residues Max Euclidean distance observed (dMax)
		* Eliminate extreme bins outliers using a percentage treshold [default = 1% of the total conformations on the ensemble]
	```bash
	$ luajit /path/to/FleXgeo/HistProc.lua
	```
	NOTE: This lua script will use the "Diffgeo_stat.lua" file generated by FleXgeo as input and output a "DiffGeoSpec.ssv". This ".ssv" contains the lenght of the interval to be considered on each dimension for dMax.

	* Compute dMax for each residue
	```bash
	$ python3.5 /path/to/FleXgeo/ComputeResDMax.py
	```
			This script use "DiffGeoSpec.ssv" as input to compute the dMax and outputs:
			1) "maxdt.csv" : the dMax values per residues
			2) "aa_list.csv" : the most flexible residues detected (more details bellow)
			3) "z-score.csv": the dMax z-score used on most flexible residues detection

	In addition to the computation of dMax, this script also run a Haar Wavelet transformation of the Z-score of residues dmax in order to automatically identify the most flexible residues.

	* Clustering conformations
	    * Run HDBSCAN
	```bash
	>$ python3.5 /path/to/GetResClusters.py DiffGeo_xgeo.csv
	```

		USAGE: GetResClusters.py [-h] [-res RES] [-out_path OUT_PATH]
		[-min_pcluster MIN_PCLUSTER] in_csv
		positional arguments:
			in_csv                csv file with FleXgeo data
		optional arguments:
			-h, --help            show this help message and exit
			-res RES              specify residue to be clustered. (default: 'ALL')
			-out_path OUT_PATH    specify dir to write output files. (default: working dir)
			-min_pcluster MIN_PCLUSTER set the minimum conformations percentage a cluster  must have (default: .05)
	* Write cluster pdbs
	```bash
	$ python3.5 /path/to/WriteClustersPDB.py cluster.clstr source.pdb res
	```

		USAGE: WriteClustersPDB.py [-h] [-out_path OUT_PATH] in_clstr src_pdb res
		Write pdb files of cluster from '.clst'.

		positional arguments:
			in_clstr            clstr file.
			src_pdb             source pdb file.
			res                 specify residue to write clusters pdbs.

		optional arguments:
			-h, --help          show this help message and exit
			-out_path OUT_PATH  specify dir to write output files. (default: working dir)

## What are the FleXgeo output files content?
FleXgeo outputs 5 .csv files:

	 1. **Diffgeo_xgeo.csv** :  contains the calculated differential geometry descriptors of the input ensemble.
	 2.  Diffgeo_NORM.csv : Same values of "Diffgeo_xgeo,csv" but normalized to [0,1].
	 3.  Diffgeo_MEAN.csv: Mean values of FleXgeo descriptors per residue
	 4.  Diffgeo_STD.csv: Standard deviation of FleXgeo descriptors per residue
	 5.  Diffgeo_VAR.csv: Variation of FleXgeo descriptors per residue
	 6.  DiffgeoStat.lua : The histogram bins position and value of each descriptor for each residue.

## How to cite?
The FleXgeo paper is currently under consideration for publication and will be added here as soon as possible.

## Why just binaries files for the core of FleXgeo?
Unfortunately, we use Numerical Recipes on our code and we are not allowed to distribute the source code. We plan to rewrite those part of the code in the future, but this is not on our top priorities right now. If you have some trouble on running a binary file on your machine, feel free to contact Antonio at amarinho@cent.uw.edu.pl and we can try to provide an specific binary for you.
