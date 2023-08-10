from os.path import join
import pandas as pd
from matplotlib import pyplot as plt
import numpy as np

def get_r():
    ct = pd.read_csv(join('..','data', 'mono_ss_report', 'ss_17_3_ct.csv'))[['B3']]
    ct['Time'] = np.arange(0, 72.1, 1/6)
    ct['species'] = 'Ct'
    ct.columns = ['OD','time','species']
    ct_rates = []
    for j in range(6):
        od = ct['OD'][j::6]
        od_r = np.gradient(od) /od
        ct_rates.append(max(od_r.dropna()))

    oa = pd.read_csv(join('..','data', 'mono_ss_report', 'oa_03_19_curves.csv'))[['Time','Citric acid']]
    oa['species'] = 'Oa'
    oa.columns = ['time','OD','species']
    oa_rates = []
    oa['time'] = np.arange(0, 44.6, 1/6)
    for j in range(6):
        od = oa['OD'][j::6]
        od_r = np.gradient(od) /od
        oa_rates.append(max(od_r.dropna()))
    
    return max(ct_rates),max(oa_rates)




def get_Y_cfus():
    ct_cfus = 2700000000.0
    oa_cfus = 13333333330.333333
    s_ct = 137.4976471744328 * 50.0 
    f_ct = 196.88529428319006 * 50.0 
    s_oa = 140.54647070640806 * 50.0 
    f_oa = 202.25294134652674 * 50.0 
    y_ct = ct_cfus / (f_ct - s_ct)
    y_oa = oa_cfus / (f_oa - s_oa)
    return y_ct, y_oa

def get_K():
    u_c,u_o = get_r()
    K_c = 137.4976471744328 * 50.0 * (u_c - 0.05) / 0.05
    K_o = 140.54647070640806 * 50.0 * (u_o - 0.05) / 0.05
    return K_c, K_o

"""r_oa,r_ct = get_r()
print(r_ct,r_oa)
y_ct,y_oa = get_Y_cfus()
print(y_ct,y_oa)
K_c,K_o = get_K()
print(K_c,K_o)"""

q = 0.075
sr = 10
x =0.2

s = (x-sr*q)/-q
k = s * (0.3-0.05)/0.05