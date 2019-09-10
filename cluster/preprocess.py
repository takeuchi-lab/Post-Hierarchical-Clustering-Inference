import numpy as np 
import pandas as pd
import sys
import os

if __name__ == "__main__":
    
    if not os.path.isdir("stat"):
        os.mkdir("stat")

    if not os.path.isdir("interval"):
        os.mkdir("interval")

    args = sys.argv
    datafile = args[1]
    data = pd.read_csv(datafile, header=None)
    data_preprocess = ((data.T - data.T.mean()) / data.T.std()).T
    data_preprocess = (data_preprocess - data_preprocess.mean()) / data_preprocess.std()
    
    # 推定用のデータサイズが実験用のデータサイズと異なる場合
#     sigma = np.mean(estimate_set.iloc[:, :4].T.var())
#     Sigma = sigma * np.identity(data_preprocess.iloc[:, :4].shape[0])
    # 推定用と実験用のデータサイズが同じ場合
    Sigma = data_preprocess.T.var()
    
    Sigma = pd.DataFrame(Sigma)
    Sigma.to_csv("data/Sigma.csv", header=None, index=False)
    
    xi = np.mean(data_preprocess.var())
    xi = np.array([xi])
    xi = pd.DataFrame(xi)
    xi.to_csv("data/xi.csv", header=None, index=False)
