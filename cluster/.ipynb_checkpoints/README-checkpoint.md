# Post Hierarchical Clustering

## Environmental Requirement

## Usage

### Preprocessing for real data

1. `interval`と`stat`というディレクトリを作成する 
2. 実データを分散推定用と$p$値計算用に分ける (目安 分散推定: $p$値計算 = 2 : 8)
3. `preprocess.py`を分散推定用のデータに使い, `sigma.csv`, `xi.csv`を出力する <br>
   例 `preprocess.py data_estimate.csv`


### Post_Hierarchical_Clustering
1. `pci_cluster_ex.exe`で`interval`と`stat`の直下にそれぞれファイルが作成される <br>
   例 `$ ./pci_cluster_ex.exe data.csv sigma.csv 0.926 4`
2. `calc_p.py`で最終的な$p$値が得られる