import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.stats as stats
import sys
import subprocess
import os
from tqdm import tqdm


if __name__ == '__main__':
    
    args = sys.argv
    
    epo = args[1]
    n = args[2]
    d = args[3]
    step = args[4]
    threads = args[5]
    mu = args[6]

    
    if not os.path.isdir("result"):
        os.mkdir("result")

    lastflag = False
    if step == "last":
        lastflag = True
        step = int(n) - 2

    
    if int(step) > int(n) - 2:
        print("error!! please input step 0 ~ n - 2")
    else:
        if float(mu) > 0.0:
            loop_list = list(np.arange(0.5, float(mu) + 0.5, 0.5))
            step = str(int(n) - 2)
        else:
            loop_list = list(np.arange(10, int(n) + 10, 10))
        naive_p_all = []
        selective_p_all = []
        for i in tqdm(loop_list):
            if float(mu) > 0.0:
                sepfile = "./result/selective_p_n" + n + "_d" + d + "_step" + step + "_mu" + str(i) + ".csv"
                napfile = "./result/naive_p_n" + n + "_d" + d + "_step" + step + "_mu" + str(i) + ".csv"
                if os.path.exists(sepfile):
                    os.remove(sepfile)
                if os.path.exists(napfile):
                    os.remove(napfile)

                if int(threads) > 1:
                    subprocess.run(['./pci_cluster_dim_synthetic_parallel.exe', epo, str(n), d, step, threads, str(i)])
                else: 
                    subprocess.run(['./pci_cluster_dim_synthetic.exe', epo, str(n), d, step, str(1), str(i)])
            else:
                if lastflag:
                    step = str(i - 2) 
                sepfile = "./result/selective_p_n" + str(i) + "_d" + d + "_step" + step + ".csv"
                napfile = "./result/naive_p_n" + str(i) + "_d" + d + "_step" + step + ".csv"
                if os.path.exists(sepfile):
                    os.remove(sepfile)
                if os.path.exists(napfile):
                    os.remove(napfile)
                if int(threads) > 1:
                    subprocess.run(['./pci_cluster_dim_synthetic_parallel.exe', epo, str(i), d, step, threads, mu])
                else: 
                    subprocess.run(['./pci_cluster_dim_synthetic.exe', epo, str(i), d, step, str(1), mu])
            
                
            # 結果 読み込み
            selective_p = pd.read_csv(sepfile, header=None).values
            naive_p = pd.read_csv(napfile, header=None).values

            selective_p_all.append(selective_p)
            naive_p_all.append(naive_p)

        selective_p_all = np.array(selective_p_all)
        naive_p_all = np.array(naive_p_all)

        if not os.path.isdir("result"):
            os.mkdir("result")
        # if float(mu) > 0.0:
        #     pd.DataFrame(naive_p_all).to_csv("result/naive_p_TPR_epoch" + epo + "_step" + step + "_n" + n + "_d" + d + ".csv", header=None, index=False)
        #     pd.DataFrame(selective_p_all).to_csv("result/selective_p_TPR_epoch" + epo + "_step" + step + "_n" + n + "_d" + d + ".csv", header=None, index=False)    
        # else:
        #     pd.DataFrame(naive_p_all).to_csv("result/naive_p_FPR_epoch" + epo + "_step" + step + "_d" + d + ".csv", header=None, index=False)
        #     pd.DataFrame(selective_p_all).to_csv("result/selective_p_FPR_epoch" + epo + "_step" + step + "_d" + d + ".csv", header=None, index=False)
        

        pr_sep = np.sum(selective_p_all < 0.05, axis=1) / int(epo)
        pr_nap = np.sum(naive_p_all < 0.05, axis=1) / int(epo)

        if not os.path.isdir("fig"):
            os.mkdir("fig")  
        for dim in range(pr_sep.shape[1]):
            plt.figure()
            plt.ylim([-0.05, 1.05])
            plt.plot(loop_list, pr_nap[:, dim], '-o',label="naive")
            plt.plot(loop_list, pr_sep[:, dim], '-o', label="selective")

            if float(mu) > 0.0:
                plt.title(str(dim) + "dim")
                plt.xlabel("$\mu$")
                plt.xticks(loop_list)
                plt.ylabel("TPR")
                plt.legend()
                plt.savefig("fig/TPR_demo_dim" + str(dim) + ".pdf", pad_inches=0, bbox_inches='tight')
                plt.show()
            else:
                plt.title(str(dim) + "dim")
                plt.xlabel("$n$")
                plt.ylabel("FPR")
                plt.legend()
                plt.savefig("fig/FPR_demo_dim" + str(dim) + ".pdf", pad_inches=0, bbox_inches='tight')
                plt.show()