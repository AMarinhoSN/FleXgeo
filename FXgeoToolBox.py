#!/usr/bin/env python
# -*- coding: utf-8 -*-
import sys
import os
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import seaborn as sns
import sklearn.cluster as cluster
import hdbscan
import time

# dep numpy, matplotlib, seaborn, sklearn, hdbscan

'''
                            >> FXgeoToolBox <<

Here we define the FleXgeo Tool Box, which aims to provide classes, methods and
functions for protein conformational ensemble analyses based on FleXgeo output.

[ADD GENERAL DISCRIPTION]
- FUNCTIONS

- CLASSES
  -- conf : provide data representation and methods related to operations with
            a single conformation
  --

'''

## FUNCTIONS ##
def loadCSV(cvs_file_name):
    """ Load data from .csv files """
    data = np.loadtxt(cvs_file_name,
                      delimiter=',',
                      dtype={'names': ('conf_i',
                                       'res_n',
                                      #'res_name',
                                       'res_curvature',
                                       'res_torsion',
                                       'res_arc_len',
                                       'res_writhing',
                                       'ca_x',
                                       'ca_y',
                                       'ca_z'
                                       ),
                            'formats':
                                      ('i4',
                                       'i4',
                                      #'S3',
                                       'f4',
                                       'f4',
                                       'f4',
                                       'f4',
                                       'f4',
                                       'f4',
                                       'f4')
                            }
                      )
    return data

### CLASSES ####

class Conf:
    '''
    Xgeo conformation class
    '''
    def __init__(self, i, nres):

        # Define the conformation atributes
        self.conf_i = i                # conformation number
        self.data = np.zeros((nres,5)) # Xgeo representation matrix
        self.nres = len(self.data)     # number of residues


    def calcD2KTperResOf(self, ConfObj_j):
        '''
        Calculate KT Euclidean distance matrix of conf i and conf j
        '''
        D2_ijMtx = np.zeros((self.nres,2))

        for ni, res_ni in enumerate(self.data):

            res_nj = ConfObj_j.data[ni]
            d2_ij = np.sqrt(#pow(res_ni[0]-res_nj[0],2) +
                            pow(res_ni[1]-res_nj[1],2) + # curvature
                            pow(res_ni[2]-res_nj[2],2))  # torsion

            D2_ijMtx[ni][0] = res_ni[0]
            D2_ijMtx[ni][1] = d2_ij

        return D2_ijMtx

class FXGeoData:
    '''
    FleXgeo ensemble representation class.

    USAGE::
       >>> import FXgeoToolBox
       >>> MyEnsemble = FXGeoData("Diffgeo_raw.csv")
       >>> MyEnsemble.plotResData(7)

    :param csv_filenm: the .csv output from FleXgeo file name.
    :rtype: class
    '''

    def __init__(self, csv_filenm):
        ''' '''
        # 1 - get data set name and data matrix
        self.name = csv_filenm.replace('.csv', '')
        self.data = loadCSV(csv_filenm)

        # 2 - get a list of conformations indexes on .csv
        try:
            self.confslist = []

            for i, conf_i in enumerate(set(self.data['conf_i'])):
                self.confslist.append(conf_i)

        except(TypeError):
            # If only one conf on data
            self.confslist = self.data['conf_i']

        # 3 - get a list of residues indexes on .csv
        try:
            self.reslist = []
            for n, res_n in enumerate(set(self.data['res_n'])):
                self.reslist.append(res_n)

        except(TypeError):
            self.reslist = self.data['res_n']

        # is it unique conformation data?
        if len(self.confslist) > 1:
            self.isUniqueConf = False
        if len(self.confslist) == 1:
            self.isUniqueConf = True

        #|---| general data slices for code readability conveniece |-----------|
        self.nconfs = len(self.confslist)
        self.nres = len(self.reslist)
        self.allres = self.data['res_n']
        self.allarcl = self.data['res_arc_len']
        self.allcurv = self.data['res_curvature']
        self.alltors = self.data['res_torsion']
        self.allwrit = self.data['res_writhing']
        self.allca_x = self.data['ca_x']
        self.allca_y = self.data['ca_y']
        self.allca_z = self.data['ca_z']
        #|---------------------------------------------------------------------|

        #4 -  Get confs dictionaries
        self.confsObjs = []

        isFirst = True
        prev_conf = None
        firstconf = self.data[0]['conf_i']
        lastconf = self.data[-1]['conf_i']
        last_idx = len(self.data) - 1

        for idx, point in enumerate(self.data):
            cur_conf = point['conf_i']
            ##| DEBUG |###
            #print "c =", cur_conf, ' p=', prev_conf
            ##############

            # if is the first conf to be read, than create a Conf object and get
            # data point
            if isFirst == True:
                ConfObj = Conf(cur_conf, self.nres)
                n = self.getindexOfRes(point['res_n'])

                ConfObj.data[n][0] = point['res_n']
                ConfObj.data[n][1] = point['res_curvature']
                ConfObj.data[n][2] = point['res_torsion']
                ConfObj.data[n][3] = point['res_arc_len']
                ConfObj.data[n][4] = point['res_writhing']
                prev_conf = cur_conf
                isFirst = False
                continue

            # if last data point, store it
            if idx == last_idx:
                self.confsObjs.append([prev_conf, ConfObj])
                continue

            # if data for a new conf, store conf index and ConfObj, create a new
            # one for the new conf

            if cur_conf != prev_conf:
                # store prev conf data
                self.confsObjs.append([prev_conf, ConfObj])
                # reset confObj to store the new data
                ConfObj = Conf(cur_conf, self.nres)
                n = self.getindexOfRes(point['res_n'])
                ConfObj.data[n][0] = point['res_n']
                ConfObj.data[n][1] = point['res_curvature']
                ConfObj.data[n][2] = point['res_torsion']
                ConfObj.data[n][3] = point['res_arc_len']
                ConfObj.data[n][4] = point['res_writhing']
                prev_conf = cur_conf
                continue

            # if still in the same conf data as previous conformation, store the
            # new data on ConfObj
            if cur_conf == prev_conf:
                n = self.getindexOfRes(point['res_n'])
                ConfObj.data[n][0] = point['res_n']
                ConfObj.data[n][1] = point['res_curvature']
                ConfObj.data[n][2] = point['res_torsion']
                ConfObj.data[n][3] = point['res_arc_len']
                ConfObj.data[n][4] = point['res_writhing']
                continue
        #----------------------------------------------------------------------|
        #--| Clustering Atributes |--------------------------------------------|
        self.lclustersPerRes = []

        #--| StructVar Matrices |----------------------------------------------|
        self.d2perResSerialMTX = []
        self.d2VarperResMTX = []
        self.d2fromrefSerialMTX = []

    #---| Methods to get res and confs index |---------------------------------|
    # These functions are needed because set objects has no atribute index
    # and are mainly used for code convenience.
    # TODO: For the final version, substitute the set function for a more
    #      decent solution. Maybe write your own set function, this problem
    #      belongs to Antonio of the future.
    ############################################################################
    def getindexOfRes(self, res):
        ''' Get index of a residues in the reslist class atribute '''
        n = None

        for n_idx, res_n in enumerate(self.reslist):
            if res_n == res:
                n = n_idx
                break
        return n

    def getResOfIndex(self, n):
        ''' Get the residues of index n in the reslist class atribute '''
        res = None

        for n_idx, res_n in enumerate(self.reslist):
            if n_idx == n:
                res = res_n
                break
        return res

    def getindexOfConf(conf):
        ''' Get index of a conf number i in the conflist class atribute '''
        i = None

        for i_idx, conf_i in enumerate(self.confslist):
            if conf_i == conf:
                i = i_idx
                break
        return i

    def getConfOfIndex(self, i):
        ''' Get conf of index i in the conflist class atribute '''
        conf = None

        for i_idx, conf_i in enumerate(self.confslist):
            if i_idx == i:
                conf = conf_i
                break
        return conf

    #|-------------------------------------------------------------------------|

    #|----| Methods for CSV data analyses |--------------------------|
    def WriteClstrMultiModelPDB(self,pdb_in, out_name,out_pdbdir):
        '''
         Write a multimodel PDB with conformations. This is usefull for cluster
         objects.
        '''
        working_dir = os.getcwd()
        PDBFile = open(pdb_in, 'r')
        PDBout = open(out_pdbdir+'/'+out_name+'.pdb', 'w')

        curr_model = 0

        for line in PDBFile:
            if line.startswith("MODEL"):
                curr_model = int(line[5:])

            if curr_model in self.confslist:
                PDBout.write(line)
                continue

    def WriteMultiModelPartialPDB(self,pdb_in, out_name,out_pdbdir):
        ''' Write PDB from Global Clusters, not all residues are on a cluster.'''
        working_dir = os.getcwd()
        PDBFile = open(pdb_in, 'r')
        PDBout = open(out_pdbdir+'/'+out_name+'.pdb', 'w')

        curr_model = 0
        curr_res = 0

        ## DEBUG ###
        #print self.name
        #print self.confslist
        #print self.reslist
        ############

        for line in PDBFile:
            if line.startswith("MODEL"):
                curr_model = int(line[5:])
                if curr_model in self.confslist:
                    PDBout.write(line)

            if line.startswith("ATOM"):
                curr_res = int(line[23:28])

            if curr_model in self.confslist and curr_res in self.reslist:
                PDBout.write(line)
                continue
            if line.startswith("ENDMDL"):
                PDBout.write(line)

    def plotConfsKTData(self, showplot=False, REF=None, REFlabel=None,
                        lang="EN", savePDF=False, Outflnm="KTrawData"):
        '''
        Plot conformations Curvature and Torsion values per residues as line
        plots.

        This function generate the plot and save it as .png file by default.
        Is possible to save as PDF, include an reference conformation on the
        plot (check parameters description). English is the default language for
        plot x and y labels, but PT-BR is supported also.

        .. note::
        If you want to plot a reference structure against your data,
        set REF as FXGeoData Object

        '''

        # If a reference conformation is provided, plot this data first
        if REF is not None:

            ref_ri_data = REF.confsObjs[0][1].data[1:-1, 0]
            ref_ki_data = REF.confsObjs[0][1].data[1:-1, 1]
            ref_ti_data = REF.confsObjs[0][1].data[1:-1, 2]
            plt.plot(ref_ri_data, ref_ki_data, color='orange', linestyle='-',
                     markerfacecolor='blue', lw=3.0, label=REFlabel)
            plt.plot(ref_ri_data, ref_ti_data, color='orange', linestyle='-',
                     markerfacecolor='blue', lw=3.0)  # label='Torsion')

        c = 0
        for conf_i in self.confsObjs:

            i = conf_i[0]
            ri_data = conf_i[1].data[1:-1, 0]
            ki_data = conf_i[1].data[1:-1, 1]
            ti_data = conf_i[1].data[1:-1, 2]
            c +=1

            if c == 1:
                if lang == "EN":
                    plt.plot(ri_data, ki_data, color='blue', linestyle='solid',
                             markerfacecolor='blue', lw=0.2, label='Curvature')
                    plt.plot(ri_data, ti_data, color='red', linestyle='solid',
                             markerfacecolor='blue', lw=0.2, label='Torsion')
                if lang == "PT":
                    plt.plot(ri_data, ki_data, color='blue', linestyle='solid',
                             markerfacecolor='blue', lw=0.2, label='Curvatura')
                    plt.plot(ri_data, ti_data, color='red', linestyle='solid',
                             markerfacecolor='blue', lw=0.2, label='Torção')
            else:
                plt.plot(ri_data, ki_data, color='blue', linestyle='solid',
                         alpha=0.3, markerfacecolor='blue', lw=0.1)
                plt.plot(ri_data, ti_data, color='red', linestyle='solid',
                         alpha=0.3, markerfacecolor='blue', lw=0.1)
        plt.legend()
        if lang == "PT":
            plt.title("Valores de Curvatura e Torção")
            plt.xlabel("Resíduo")
        if lang == "EN":
            plt.title("Curvature and Torsion values")
            plt.xlabel("Residue")

        plt.xlim(self.reslist[0], self.reslist[-1])
        plt.savefig(Outflnm+".png", dpi=600)
        if savePDF is True:
            plt.savefig("KTrawData.pdf", dpi=300)
        if showplot is True:
            plt.show()
        plt.close()
        return None


    def plotResData(self, res):
        ''' Plot data for individual residues.'''

        ca_x = []
        ca_y = []
        ca_z = []
        k = []
        t = []
        a = []
        w = []

        for each in self.data:
            if each['res_n'] == res:
                ca_x.append(each['ca_x'])
                ca_y.append(each['ca_y'])
                ca_z.append(each['ca_z'])
                k.append(each['res_curvature'])
                t.append(each['res_torsion'])
                a.append(each['res_arc_len'])

        Plot3D(ca_x,ca_y,ca_z,self.name+'_xyz')
        Plot3D(k,t,n,self.name+'_ktn')
        return None
