#!/usr/bin/env python3
# Copyright Chun Shen @ 2018

from numpy import *
import sys
from os import path
import matplotlib.pyplot as plt


def inverse_mapping_2D(x, y, xt, yt, z1, z2, z3):
    """
    This function perform a 2D inverse mapping
    """
    nx, ny = xt.shape
    idx_x_array = zeros(ny)
    y_cut  = zeros(ny)
    z1_cut = zeros(ny)
    z2_cut = zeros(ny)
    z3_cut = zeros(ny)
    for iy in range(ny):
        idx_x = searchsorted(xt[:, iy], x, side='left')
        idx_x = max(0, min(nx - 2, idx_x - 1))
        x_1 = xt[idx_x, iy]
        x_2 = xt[idx_x + 1, iy]
        y_1 = yt[idx_x, iy]
        y_2 = yt[idx_x + 1, iy]
        z1_1 = z1[idx_x, iy]
        z1_2 = z1[idx_x + 1, iy]
        z2_1 = z2[idx_x, iy]
        z2_2 = z2[idx_x + 1, iy]
        z3_1 = z3[idx_x, iy]
        z3_2 = z3[idx_x + 1, iy]
        y_cut[iy]  = y_1 + (y_2 - y_1)/(x_2 - x_1)*(x - x_1)
        z1_cut[iy] = z1_1 + (z1_2 - z1_1)/(x_2 - x_1)*(x - x_1)
        z2_cut[iy] = z2_1 + (z2_2 - z2_1)/(x_2 - x_1)*(x - x_1)
        z3_cut[iy] = z3_1 + (z3_2 - z3_1)/(x_2 - x_1)*(x - x_1)
    if y > y_cut[-1]:
        z1_interp = z1_cut[-1]
        z2_interp = z2_cut[-1]
        z3_interp = z3_cut[-1]
    else:
        idx_y = min(ny - 2, max(0, searchsorted(y_cut, y, side='left') - 1))
        z1_interp = z1_cut[idx_y] + (z1_cut[idx_y + 1] - z1_cut[idx_y])/(y_cut[idx_y + 1] - y_cut[idx_y])*(y - y_cut[idx_y])
        z2_interp = z2_cut[idx_y] + (z2_cut[idx_y + 1] - z2_cut[idx_y])/(y_cut[idx_y + 1] - y_cut[idx_y])*(y - y_cut[idx_y])
        z3_interp = z3_cut[idx_y] + (z3_cut[idx_y + 1] - z3_cut[idx_y])/(y_cut[idx_y + 1] - y_cut[idx_y])*(y - y_cut[idx_y])
    z1_interp = max(0., z1_interp)
    z2_interp = max(0., z2_interp)
    z3_interp = max(0., z3_interp)
    return(z1_interp, z2_interp, z3_interp)


try:
    EOS_table_folder = str(sys.argv[1])
except:
    print("Usage: python3 mapEOS_TmuB_to_enB.py EOS_table_folder")
    exit(1)


hbarC = 0.19733

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

print("e_min = %.5e GeV/fm^3, e_max = %.5e GeV/fm^3" % (matrix(ed_table).min(), matrix(ed_table).max()))
print("nB_min = %.5e 1/fm^3, nB_max = %.5e 1/fm^3" % (matrix(nB_table).min(), matrix(nB_table).max()))

#ed_local = 1.0
#nB_local = 5.0
#P_local, T_local, muB_local = inverse_mapping_2D(
#    ed_local, nB_local, ed_table, nB_table, P_table, T_table, muB_table)
#s_local = (ed_local + P_local - muB_local*nB_local)/T_local
#print(ed_local, nB_local, P_local, T_local, muB_local, s_local, (ed_local + P_local)/T_local)


ed_bounds = [0.0, 0.0036, 0.015, 0.045, 0.455, 20.355, 650.]
ne_list   = [13, 20, 31, 42, 200, 400]
nB_bounds = [0.005, 0.015, 0.045, 0.5, 3.5, 12.0]
nnB_list  = [501, 301, 181, 251, 351, 251]

# generate tables
for itable in range(len(ne_list)):
    print("Generating table %d ... " % itable)
    ed_list = linspace(ed_bounds[itable], ed_bounds[itable+1], ne_list[itable])
    nB_list = linspace(0.0, nB_bounds[itable], nnB_list[itable])
    p_list = zeros([len(ed_list), len(nB_list)])
    T_list = zeros([len(ed_list), len(nB_list)])
    muB_list = zeros([len(ed_list), len(nB_list)])
    
    for ie in range(len(ed_list)):
        ed_local = ed_list[ie]
        for inB in range(len(nB_list)):
            nB_local = nB_list[inB]
            P_local, T_local, muB_local = inverse_mapping_2D(
                ed_local, nB_local, ed_table, nB_table, P_table, T_table, muB_table)
            p_list[ie, inB] = P_local
            T_list[ie, inB] = T_local
            muB_list[ie, inB] = muB_local
    savetxt("BEST_eos_p_%d.dat" % itable, p_list, fmt='%.8e', delimiter="  ")
    savetxt("BEST_eos_T_%d.dat" % itable, T_list, fmt='%.8e', delimiter="  ")
    savetxt("BEST_eos_muB_%d.dat" % itable, muB_list, fmt='%.8e', delimiter="  ")
