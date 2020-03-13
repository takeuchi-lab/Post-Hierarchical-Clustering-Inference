import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import scipy.stats as stats
import sys
import subprocess
import os
from tqdm import tqdm
from mpmath import *

def calc_p(stat, df, final_interval):
    """
    近似計算なしの切断カイ二乗分布でのp値計算
    """
    l = final_interval[:, 0]
    u = final_interval[:, 1]
    deno = np.sum(stats.chi2.cdf(u,df) - stats.chi2.cdf(l,df))
    p = np.argmax( (l < stat) * (stat < u) )
    add = np.sum( stats.chi2.cdf(u[p+1:], df) - stats.chi2.cdf(l[p+1:], df) ) if p + 1 < len(l) else 0.0
    nume = stats.chi2.cdf(u[p], df) - stats.chi2.cdf(stat, df) + add
    return nume / deno

def calc_multi_p(tau, mean, std, L, U):
    """
    区間は昇順にソートされていることを想定
    """
    deno = np.sum(stats.norm.cdf(-(L - mean) / std) - stats.norm.cdf(-(U - mean) / std))
    p = np.argmax( (L < tau) * (tau < U) )
    add = np.sum( stats.norm.cdf(-(L[p+1:] - mean) / std) - stats.norm.cdf(-(U[p+1:] - mean) / std) ) if p + 1 < len(L) else 0.0
    nume = stats.norm.cdf(-(tau - mean) / std) - stats.norm.cdf(-(U[p] - mean) / std) + add 
    return nume / deno

def calc_p_approx(stat, df, final_interval):
    """
    切断カイ二乗分布でのp値計算 正規近似var
    """
    L = final_interval[:, 0]
    U = final_interval[:, 1]
    # 近似
    Lx = (stat / df)**(1. / 6) - (stat / df)**(1. / 3) / 2 + (stat / df)**(1. / 2) / 3
    mean = 5. / 6 - 1. / (9 * df) - 7. / (648 * df**2) + 25. / (2187 * df**3)
    var = 1. / (18 * df) + 1. / (162 * df**2) - 37. / (11664 * df**3)
    
    LL = (L / df)**(1. / 6) - (L / df)**(1. / 3) / 2 + (L / df)**(1. / 2) / 3
    LU = (U / df)**(1. / 6) - (U / df)**(1. / 3) / 2 + (U / df)**(1. / 2) / 3
    
    selective_p = calc_multi_p(Lx, mean, np.sqrt(var), LL, LU)
    return selective_p

def calc_norm_mp(stat, final_interval, digit=1000):
    """
    高精度な切断正規分布でのp値計算
    """
    mp.dps = digit
    l = final_interval[:, 0]
    u = final_interval[:, 1]
    deno = 0
    nume = 0
    for i in range(len(l)):
        deno += mp.ncdf(u[i]) - mp.ncdf(l[i])
        if l[i] < stat < u[i]:
            nume += mp.ncdf(u[i]) - mp.ncdf(stat)
        elif u[i] > stat :
            nume += mp.ncdf(u[i]) - mp.ncdf(l[i])
    return float(nume / deno)

def chi2cdf(stat, df, digit=1000):
    """
    カイ二乗分布の累積分布関数
    """
    mp.dps = digit
    if stat < 0:
        return 0.0
    else:
        gamma = mp.gamma(float(df) / 2)
        imcomp = mp.gammainc(float(df) / 2, 0, float(stat) / 2)
        return imcomp / gamma
    
def calc_chi2_mp(stat, df, final_interval, digit=1000):
    """
    高精度な切断カイ二乗分布でのp値計算
    """
    l = final_interval[:, 0]
    u = final_interval[:, 1]
    deno = 0
    nume = 0
    for i in range(len(l)):
        deno += chi2cdf(u[i], df, digit) - chi2cdf(l[i], df, digit)
        if l[i] < stat < u[i]:
            nume += chi2cdf(u[i], df, digit) - chi2cdf(stat, df, digit)
        elif u[i] > stat:
            nume += chi2cdf(u[i], df, digit) - chi2cdf(l[i], df, digit)
    return float(nume / deno)

if __name__ == '__main__':
    
    args = sys.argv
    datafile = args[1]
    sigmafile = args[2]
    xifile = args[3]
    threads = int(args[4])
    
    data = pd.read_csv(datafile, header=None)
    all_step = data.shape[0] - 1
    
    
    xi = pd.read_csv(xifile, header=None).values.reshape(-1, )[0]
    
    if not os.path.isdir("cluster_result"):
        os.mkdir("cluster_result")
    
    if not os.path.isdir("stat"):
        os.mkdir("stat")

    if not os.path.isdir("interval"):
        os.mkdir("interval")

    naive_p_all = []
    selective_p_all = []
    for i in tqdm(range(all_step)):
        statfile = "stat/test_stat" + "_step" + str(i) + ".csv"
        intervalfile = "interval/final_interval" + "_step" + str(i) + ".csv"
        if os.path.exists(statfile):
            os.remove(statfile)
        if os.path.exists(intervalfile):
            os.remove(intervalfile)
            
        if threads > 1:
            subprocess.run(['./pci_cluster_ex_parallel.exe', datafile, sigmafile, str(xi), str(i), str(threads)])
        else:
            subprocess.run(['./pci_cluster_ex.exe', datafile, sigmafile, str(xi), str(i), str(1)])

        # 結果 読み込み
        chi2, d = pd.read_csv(statfile, header=None).values.reshape(2, )
        final_interval = pd.read_csv(intervalfile, header=None).values.reshape(-1, 2)

        naive_p = 1 - stats.chi2.cdf(chi2, d)
        selective_p = calc_chi2_mp(chi2, d, final_interval)
        
        naive_p_all.append(naive_p)
        selective_p_all.append(selective_p)

    naive_p_all = pd.DataFrame(naive_p_all)
    selective_p_all = pd.DataFrame(selective_p_all)
    
    if not os.path.isdir("result"):
        os.mkdir("result")
    naive_p_all.to_csv("result/naive_p.csv", header=None, index=False)
    selective_p_all.to_csv("result/selective_p.csv", header=None, index=False)