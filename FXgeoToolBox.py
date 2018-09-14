#!/usr/bin/env python
# -*- coding: utf-8 -*-
import sys
import os
import numpy as np
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import seaborn as sns
sns.set()
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

# TODO add dmax script, Comparison script, add clustering script

# FUNCTIONS #
def loadCSV(cvs_file_name):
    """ Load data from .csv files """

    data = np.loadtxt(cvs_file_name,
                      delimiter=',',
                      dtype={'names': ('conf_i',
                                       'res_n',
                                      # 'res_name',
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
                                      # 'S3',
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

# CLASSES #

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
        # ---------------------------------------------------------------------|
        # -| Clustering Atributes |--------------------------------------------|
        self.lclustersPerRes = []

        # -| StructVar Matrices |----------------------------------------------|
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

    # -------------------------------------------------------------------------|

    # ----| Methods for CSV data analyses |--------------------------|
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


    def plotDataViolin(self, showplot=False, scale="count", inner="stick",
                       linewidth=.1, bw=.1, lang="EN", xticks=10, xUpLim=None,
                       xLoLim=None, Tor_yLoLim=None,  Tor_yUpLim=None,
                       Cur_yLoLim=None,  Cur_yUpLim=None):

        ''' Plot violin graphs of xgeo descriptors values'''
        if lang == "EN":
            tor_label = "Torsion"
            curv_label = "Curvature"
            res_label = "Residue"
        if lang == "PT":
            tor_label = "Torção"
            curv_label = "Curvatura"
            res_label = "Resíduo"

        # plot with 2 axes
        plt.figure(1)
        # Torsion
        plt.subplot(211)
        ax = sns.violinplot(x=self.data['res_n'], y=self.data['res_torsion'],
                            hue=None, data=None, split=True, inner=inner,
                            bw=bw, scale=scale, linewidth=linewidth)
        plt.xticks(range(1, self.nres, xticks),
                   range(self.data['res_n'][0], self.data['res_n'][-1], xticks))
        if Tor_yLoLim is not None:
            plt.ylim(ymin=Tor_yLoLim)
        if Tor_yUpLim is not None:
            plt.ylim(ymax=Tor_yUpLim)
        if xLoLim is not None:
            plt.xlim(xmin=xLoLim)
        if xUpLim is not None:
            plt.xlim(xmax=xUpLim)

        # plt.xlabel("Residues")
        plt.ylabel(tor_label)
        plt.grid(True)

        # Curvature
        plt.subplot(212)
        ax = sns.violinplot(x=self.data['res_n'], y=self.data['res_curvature'],
                            hue=None, data=None, split=True, inner=inner,
                            bw=bw, scale=scale, linewidth=linewidth)
        plt.xlabel(res_label)
        plt.ylabel(curv_label)
        plt.xticks(range(1, self.nres, xticks),
                   range(self.data['res_n'][0], self.data['res_n'][-1], xticks))
        plt.grid(True)

        if Cur_yLoLim is not None:
            plt.ylim(ymin=Cur_yLoLim)
        if Cur_yUpLim is not None:
            plt.ylim(ymax=Cur_yUpLim)
        if xLoLim is not None:
            plt.xlim(xmin=xLoLim)
        if xUpLim is not None:
            plt.xlim(xmax=xUpLim)

        if showplot is True:
            plt.show()

        plt.savefig("ViolinKTJointPlot.png", dpi=600)
        plt.close()

        # Plot individual plots

        # plt.close("All")
        # TORSION
        ax = sns.violinplot(x=self.data['res_n'], y=self.data['res_torsion'],
                            hue=None, data=None, split=True, inner=inner,
                            bw=bw, scale=scale, linewidth=linewidth)
        plt.xlabel(res_label)
        plt.ylabel(tor_label)
        plt.xticks(range(1, self.nres, xticks),
                   range(self.data['res_n'][0], self.data['res_n'][-1], xticks))
        plt.savefig("ViolinTOR.png", dpi=600)

        if showplot is True:
            plt.show()
        plt.close()

        # CURVATURE
        ax = sns.violinplot(x=self.data['res_n'], y=self.data['res_curvature'],
                            hue=None, data=None, split=True, inner=inner,
                            bw=bw, scale=scale, linewidth=linewidth)
        plt.xlabel(res_label)
        plt.ylabel(curv_label)
        plt.xticks(range(1, self.nres, xticks),
                   range(self.data['res_n'][0], self.data['res_n'][-1], xticks))
        plt.savefig("ViolinCURV.png", dpi=600)
        if showplot is True:
            plt.show()
        plt.close()

        # ARC LENGTH
        ax = sns.violinplot(x=self.data['res_n'], y=self.data['res_arc_len'],
                            hue=None, data=None, split=True, inner=inner,
                            bw=bw, scale=scale, linewidth=linewidth)
        plt.xlabel(res_label)
        plt.ylabel("Arc Length")
        plt.ylim(6.0)
        plt.xticks(range(1, self.nres, xticks),
                   range(self.data['res_n'][0], self.data['res_n'][-1], xticks)
                   )
        plt.savefig("ViolinARC.png", dpi=600)

        if showplot is True:
            plt.show()
        plt.close()

        # WRITHING
        ax = sns.violinplot(x=self.data['res_n'], y=self.data['res_writhing'],
                            hue=None, data=None, split=True, inner=inner,
                            bw=bw, scale=scale, linewidth=linewidth)
        plt.xlabel(res_label)
        plt.ylabel("Writhing")
        plt.xticks(range(1, self.nres, xticks),
                   range(self.data['res_n'][0], self.data['res_n'][-1], xticks)
                   )
        plt.savefig("ViolinWRI.png", dpi=600)
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

    # -----| methods to compare conformations |--------------------------------|
    def calcD2ConfperRes(self, cref_idx=0, ext_ref=None, showplot=False, xsticks=10, ysticks=10,
                     lang="EN", bw=.1, linewidth=.1, scale="count", inner=None,
                     out_prefix='LocalStructStab'):
        '''
        Calcualte KT Euclidean distance of all conformations to an arbitrary
        conformation on the ensemble.
        The default behaviour is calculate for the first conformation, but the
        user can specicify the index of the reference conformation specified by
        cref_idx (int).
        '''

        # Get language labels for plots
        if lang not in ("EN", "PT"):
            print(" WARNING: You need to provide a supported language for plot")
            print("          labels. Current valid options are 'EN' and 'PT'")
            print("          Setting to 'EN'")
            lang = "EN"

        if lang is "PT":
            reslabel = "Resíduo"
            conflabel = "Índice da conformação"
            d2label = "Distância Euclidiana"
        if lang is "EN":
            reslabel = "Residue"
            conflabel = "Conformation index"
            d2label = "Euclidean distance"

        # 0 - get the first conf data and set it as ref values
        if ext_ref is None:
            try:
                conf_ref = self.confsObjs[cref_idx][1]

            except(TypeError):
                print('ERROR: cref_idx provided is not a valid integer')
                exit(0)
            # TODO : add the error for cref_idx not in the range of the ensemble

        if ext_ref is not None:
            conf_ref = ext_ref.confsObjs[0][1]

        SerialMtx = np.zeros((self.nres, self.nconfs+1))

        # fill res columns
        for n, res_n in enumerate(self.reslist):
            SerialMtx[n][0] = res_n

        # calculate d2_ij for all residues and store results on Mtx
        for j_idx, conf_j in enumerate(self.confslist):
            # print j_idx, " - ", self.confsObjs[j_idx]
            conf_j = self.confsObjs[j_idx][1]
            d2_refjMTX = conf_ref.calcD2KTperResOf(conf_j)
            for each in d2_refjMTX:
                if each[1] < 0:
                    print(d2_refjMTX)
                    print("WARN: Negative value")
            for n, d2ij in enumerate(d2_refjMTX):
                SerialMtx[n][j_idx+1] = d2_refjMTX[n][1]

        self.d2perResSerialMTX = SerialMtx

        # Format data to plot
        data4plot = []
        for res_line in SerialMtx:
            res_d2_data = []
            res = res_line[0]
            isFirst = True

            for i, conf_i in enumerate(self.confslist):
                # ## DEBUG ###
                # print res, ", ", res_line[i+1],",", conf_i
                # ############

                d2_resn_ij = res_line[i+1]

                if isFirst is True:
                    res_d2_data = np.array([res, d2_resn_ij, conf_i])
                    isFirst = False
                    # print res_d2_data
                else:
                    data_point = np.array([res, d2_resn_ij, conf_i])
                    res_d2_data = np.vstack((res_d2_data, data_point))
                    # res_d2_data.append(data_point)

            data4plot.append(res_d2_data)

        # Generate 3D plot
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')

        # Get one color for each residue line
        colors = matplotlib.cm.rainbow(np.linspace(0, 1, len(data4plot)))

        # Generate individual plots
        for i, res_d in enumerate(data4plot):
            x = res_d[..., 0]  # res
            z = res_d[..., 1]  # d2
            y = res_d[..., 2]  # conf

            ax.plot(xs=x, ys=y, zs=z, c=colors[i], linewidth=0.5,
                    fillstyle='full')
        # set general plot labels and save .png
        ax.set_xlabel(reslabel)
        ax.set_ylabel(conflabel)
        ax.set_zlabel(d2label)

        plt.savefig("test.png", dpi=600)
        if showplot is True:
            plt.show()
        plt.close()

        # write csv out file
        out_f = open("LocalDistFromRef.csv", 'w')
        out_f.write('#res,')
        for idx, col in enumerate(SerialMtx[0]):
            if idx == 0:
                continue
            else:
                out_f.write('d2_n1'+str(idx)+', ')
            if idx == len(SerialMtx[0])-1:
                out_f.write('d2_n1'+str(idx)+'\n')

        for line in SerialMtx:
            for idx, col in enumerate(line):
                out_f.write(str(col)+',')
                if idx == len(line)-1:
                    out_f.write(str(col)+'\n')
        out_f.close()

        # PLOT MATRIX
        #print(SerialMtx)
        sns.heatmap(SerialMtx, vmax=.8, square=False,
                    xticklabels=xsticks,
                    yticklabels=ysticks,
                    cmap='Reds')
        plt.xlabel(conflabel)
        plt.ylabel(reslabel+" index")
        plt.savefig("LocalStrucStab_MTX.png", dpi=600)
        #if showplot is True:
        #    plt.show()
        plt.close()

        # ## WARNING ###
        # # GAMBIARRA ALERT ###########################################
        # # Para plotar o violin tive que reestruturar o SerialMTX   ##
        # # Future Antonio,
        # # I know you hate me by now, but try to rework the SerialMTX
        # # representation to avoid such gambiarras
        # ##############################################################
        # ViolinMTX = np.zeros((self.nres * self.nconfs, 2))

        ViolinArray = np.array([0, 0])

        for n, line in enumerate(SerialMtx):
            for i, col in enumerate(line):
                if i == 0:
                    continue
                # print line[0], line[i]

                dtapoint = np.array([line[0], line[i]])
                if line[i] < 0:
                    print(dtapoint)
                ViolinArray = np.vstack((ViolinArray, dtapoint))

        # Plot Violin d2
        #print(ViolinArray)
        x = ViolinArray[1:-1, 0]
        y = ViolinArray[1:-1, 1]

        ax = sns.violinplot(x=x, y=y, hue=None, data=None,
                            split=False,
                            scale=scale,
                            inner=inner,
                            bw=bw,
                            linewidth=linewidth)

        plt.xticks(range(1, self.nres, xsticks),
                   range(self.data['res_n'][0], self.data['res_n'][-1], xsticks))
        plt.xlabel("Residues")
        plt.ylabel(r"$\kappa-\tau$ Euclidean Distance")
        # plt.ylim(0,max(y))
        # plt.title("Violin plot of KT Euclidean Distance from a reference structure*")
        plt.savefig("ViolinKTDistFromRef.png", dpi=600)
        plt.show()
        plt.close()
