import numpy as np 
import pandas as pd
import sys

if __name__ == "__main__":
    args = sys.argv
    datafile = args[1]
    data = pd.read_csv(datafile, header=None)
    data_preprocess = ((data.T - data.T.mean()) / data.T.std()).T
    data_preprocess = (data_preprocess - data_preprocess.mean()) / data_preprocess.std()
    
    # 推定用のデータサイズが実験用のデータサイズと異なる場合
#     sigma = np.mean(estimate_set.iloc[:, :4].T.var())
#     Sigma = sigma * np.identity(data_preprocess.iloc[:, :4].shape[0])
#     xi = np.mean(estimate_set.iloc[:, :4].var())
#     Xi = xi * np.identity(data_preprocess.iloc[:, :4].shape[1])
    
    # 推定用と実験用のデータサイズが同じ場合
    Sigma = data_preprocess.T.cov()    
    Sigma = pd.DataFrame(Sigma)
    Xi = data_preprocess.cov()
    Xi = pd.DataFrame(Xi)

    Sigma.to_csv("data/Sigma.csv", header=None, index=False)
    Xi.to_csv("data/Xi.csv", header=None, index=False)
