import argparse
from os import getcwd

'''
This script write pdbs based on '.clstr' files
'''

# Define arguments
parser = argparse.ArgumentParser(description='''Write pdb files of cluster
from '.clst'.''')

parser.add_argument("in_clstr", type=str, help='''clstr file.''')

parser.add_argument("src_pdb", type=str, help='''source pdb file.''')

parser.add_argument("res", type=int, help='''specify residue to write clusters
 pdbs.''')

parser.add_argument("-out_path",default=getcwd(), type=str, help='''specify
dir to write output files. (default: working dir)''')

parser.parse_args()
parser.print_help()
args = parser.parse_args()

# load .clstr and .pdb
clstrFile = open(args.in_clstr, "r")
pdbFile = open(args.src_pdb,"r")
resIn = args.res
out_path = args.out_path

print("@ reading data from {}".format(clstrFile))
thisRes = False

# Get cluster confs info
for line in clstrFile:

    # check residue
    if line.startswith("> "):
        curr_res = int(line[2:-1])
        if curr_res == resIn:
            thisRes = True
            continue
        else:
            thisRes = False
            continue

    # if res, get confs data
    if thisRes is True and line.startswith(">>"):
        # Cluster Info
        ddot_idx = line.index(":")
        k = int(line[2:ddot_idx])
        kconfs_str = line[ddot_idx+1:-1].split(",")
        kconfs = [int(c) for c in kconfs_str]
        print("k =", k)
        pdbnm = out_path+"Res_"+str(resIn)+"Clstr_"+str(k)+".pdb"
        kpdb_out = open(pdbnm, "w")
        thisModel = False
        print(":: Writing {}".format(pdbnm))
        ### FIX CODE BUG ####
        # after first iteration on pdbFile, there is not anything left to read
        # We need to restart the file pointer to the beginning
        pdbFile.seek(0)
        #####################
        for line_pdb in pdbFile:
            if line_pdb.startswith("MODEL"):
                curr_model = int(line_pdb[5:-1])
                if curr_model in kconfs:
                    thisModel = True
                    kpdb_out.write(line_pdb)
                    continue
                else:
                    thisModel = False
                    continue

            if thisModel is True:
                if line_pdb.startswith("ATOM"):
                    kpdb_out.write(line_pdb)
                    continue
                if line_pdb.startswith("ENDMDL") or line_pdb.startswith("END") or line_pdb.startswith("TER"):
                    kpdb_out.write(line_pdb)
                    continue
        kpdb_out.close()
        continue

print(":: DONE ::")
