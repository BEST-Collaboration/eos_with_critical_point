#/usr/bin/env python3
# Copyright Chun Shen @ 2018

from numpy import *
import sys
from os import path
from scipy import interpolate
import matplotlib.pyplot as plt
import multiprocessing as mp


try:
    EOS_table_folder = str(sys.argv[1])
except:
    print("Usage: python3 mapEOS_TmuB_to_enB.py EOS_table_folder")
    exit(1)

ACCURACY = 1e-6
hbarC    = 0.19733

filename_string = EOS_table_folder.split("Files_")[1].split("/")[0]

ed_table = loadtxt(path.join(EOS_table_folder,
                             "EnerDens_Final_%s_3D.dat" % filename_string))
nB_table = loadtxt(path.join(EOS_table_folder,
                             "BarDens_Final_%s_3D.dat" % filename_string))
P_table  = loadtxt(path.join(EOS_table_folder,
                             "Press_Final_%s_3D.dat" % filename_string))

n_muB = 451
n_T   = 771
muB = linspace(0.000, 0.450, n_muB)
T   = linspace(0.030, 0.800, n_T)

muB_table, T_table = meshgrid(muB, T)
ed_table = ed_table[:, 2].reshape(n_T, n_muB)*(T_table**4.)/(hbarC**3.)
nB_table = nB_table[:, 2].reshape(n_T, n_muB)*(T_table**3.)/(hbarC**3.)
P_table  = P_table[:, 2].reshape(n_T, n_muB)*(T_table**4.)/(hbarC**3.)

print("e_min = %.5e GeV/fm^3, e_max = %.5e GeV/fm^3"
        % (matrix(ed_table).min(), matrix(ed_table).max()))
print("nB_min = %.5e 1/fm^3, nB_max = %.5e 1/fm^3"
        % (matrix(nB_table).min(), matrix(nB_table).max()))

f_p  = interpolate.interp2d(muB, T, P_table,  kind='cubic')
f_e  = interpolate.interp2d(muB, T, ed_table, kind='cubic')
f_nB = interpolate.interp2d(muB, T, nB_table, kind='cubic')

def binary_search_1d(ed_local, muB_local):
    T_min = 0.03; T_max = 0.80
    e_low = f_e(muB_local, T_min)
    e_up  = f_e(muB_local, T_max)
    if (ed_local < e_low):
        return(T_min)
    elif (ed_local > e_up):
        return(T_max)
    else:
        T_mid = (T_max + T_min)/2.
        e_mid = f_e(muB_local, T_mid)
        abs_err = abs(e_mid - ed_local)
        rel_err = abs_err/abs(e_mid + ed_local + 1e-15)
        while (rel_err > ACCURACY and abs_err > ACCURACY*1e-2):
            if (ed_local < e_mid):
                T_max = T_mid
            else:
                T_min = T_mid
            T_mid = (T_max + T_min)/2.
            e_mid = f_e(muB_local, T_mid)
            abs_err = abs(e_mid - ed_local)
            rel_err = abs_err/abs(e_mid + ed_local + 1e-15)
        return(T_mid)

def binary_search_2d(ed_local, nB_local):
    muB_min = 0.0; muB_max = 0.45
    T_max = binary_search_1d(ed_local, muB_min)
    nB_min = f_nB(muB_min, T_max)
    T_min = binary_search_1d(ed_local, muB_max)
    nB_max = f_nB(muB_max, T_min)
    if (nB_local < nB_min):
        return(T_max, muB_min)
    elif (nB_local > nB_max):
        return(T_min, muB_max)
    else:
        muB_mid = (muB_min + muB_max)/2.
        T_mid = binary_search_1d(ed_local, muB_mid)
        nB_mid = f_nB(muB_mid, T_mid)
        abs_err = abs(nB_mid - nB_local)
        rel_err = abs_err/abs(nB_mid + nB_local + 1e-15)
        while (rel_err > ACCURACY and abs_err > ACCURACY*1e-2):
            if (nB_local < nB_mid):
                muB_max = muB_mid
            else:
                muB_min = muB_mid
            muB_mid = (muB_max + muB_min)/2.
            T_mid = binary_search_1d(ed_local, muB_mid)
            nB_mid = f_nB(muB_mid, T_mid)
            abs_err = abs(nB_mid - nB_local)
            rel_err = abs_err/abs(nB_mid + nB_local)
        return(T_mid, muB_min)

#T_local, muB_local = binary_search_2d(1.0, 0.02)
#print(T_local, muB_local, f_e(muB_local, T_local), f_nB(muB_local, T_local))

ed_bounds = [0.0, 0.0036, 0.015, 0.045, 0.455, 20.355, 650.]
ne_list   = [13, 20, 31, 42, 200, 400]
nB_bounds = [0.005, 0.015, 0.045, 0.5, 3.5, 12.0]
nnB_list  = [501, 301, 181, 251, 351, 251]

# generate tables
for itable in range(len(ne_list)):
    print("Generating table %d (%d x %d)... "
            % (itable, ne_list[itable], nnB_list[itable]))

    # generate the grid arrays
    ed_list  = linspace(ed_bounds[itable], ed_bounds[itable+1], ne_list[itable])
    de       = ed_list[1] - ed_list[0]
    nB_list  = linspace(0.0, nB_bounds[itable], nnB_list[itable])
    drhob    = nB_list[1] - nB_list[0]
    p_list   = zeros([len(ed_list), len(nB_list)])
    T_list   = zeros([len(ed_list), len(nB_list)])
    muB_list = zeros([len(ed_list), len(nB_list)])
   
    def invert_EOS_tables(idx):
        ie       = int(idx/len(nB_list))
        inB      = idx % len(nB_list)
        ed_local = ed_list[ie]
        nB_local = nB_list[inB]
        T_local, muB_local = binary_search_2d(ed_local, nB_local)
        P_local            = f_p(muB_local, T_local)[0]
        return(ie, inB, P_local, T_local, muB_local)

    # calculate every points in parallel
    inputs  = range(len(ed_list)*len(nB_list))
    pool    = mp.Pool(processes=4)
    results = pool.map(invert_EOS_tables, inputs)

    # fill into the tables
    for idx in range(len(results)):
        ie                = results[idx][0]
        inB               = results[idx][1]
        p_list[ie, inB]   = results[idx][2]
        T_list[ie, inB]   = results[idx][3]
        muB_list[ie, inB] = results[idx][4]

    # save to files
    savetxt("BEST_eos_p_%d.dat" % itable, p_list, fmt='%.8e', delimiter="  ",
            header="e0= %.8e de= %.8e Ne= %d rhob0= %.8e drhob0= %.8e Nrhob= %d"
                    % (ed_bounds[itable], de, ne_list[itable],
                       0.0, drhob, nnB_list[itable]))
    savetxt("BEST_eos_T_%d.dat" % itable, T_list, fmt='%.8e', delimiter="  ")
    savetxt("BEST_eos_muB_%d.dat" % itable, muB_list, fmt='%.8e', delimiter="  ")
