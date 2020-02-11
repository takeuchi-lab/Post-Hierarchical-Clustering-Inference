# Post Hierarchical Clustering Inference
画像(仮)
<div align="center">

![](figs/dendrogram_p_black.svg)

</div>

## Abstract
DNA マイクロアレイにより計測された遺伝子発現量データの解析が遺伝学分野で盛んにおこなわれている. 例えば, マイクロアレイデータに対してクラスタリングを適用し, 得られたクラスタを特徴づける遺伝子の統計的有意性を検証する. これにより, 疾患のサブタイプを同定するという試みである. しかし, このようなクラスタリング後の推論において古典的な検定手法を用いてしまうと, クラスタリングの影響を考慮できず, 検定におけるtype I error を制御できない. そこで, 本研究では, 階層型クラスタリングのウォード法に対して, Lee et al. [1] によって提案されたSelective Inference の枠組みを適用し, 妥当な検定を行う手法を提案する.
## Environmental Requirement
- c++ソースをコンパイルする場合はeigenをインストールする
- Python version 3.7.0
- Please install required packages when the python "ImportError" occurs
  

## Usage
### クラスタ中心間の検定
`cluster`ディレクトリ以下にあるものが相当する
#### コンパイル (実行ファイル作成)
次のコマンドで実行ファイルが`/cluster`直下に作成される
    
    $ cd /cluster/cpp_source
    $ make
`pci_cluster_ex`, および並列実行用`pci_cluster_ex_parallel`    

#### Preprocessing
1. `data`というディレクトリを作成する 
2. `data`ディレクトリに使いたいデータを入れ, 分散($\Sigma$, $\xi^2$)推定用とp値計算用に分ける (目安 分散推定: p値計算 = 2 : 8)
3. `preprocess.py`によって`data`ディレクトリ以下に `sigma.csv`, `xi.csv`が出力され, `stat`と`interval`のディレクトリが作成されていることを確認する <br>

使用例 <br>
    
    $ preprocess.py data_estimate.csv


#### クラスタリングと検定の実行

##### 一つの階層でのp値を計算したい場合
`calc_p.py`を実行する
以下の引数が必要である
- データファイルパス
- $\Sigma$のファイルパス
- $\xi^2$のファイルパス
- どの階層か

使用例: <br>

    $ pyton calc_p.py data.csv data/sigma.csv data/xi.csv 0 

#### 全ての階層でのp値を計算したい場合
`calc_p_all.py`を実行する
必要な引数
- データファイルパス
- $\Sigma$のファイルパス
- $\xi^2$のファイルパス
- 全部で何階層か? (サンプルサイズ-1が全階層となる)
使用例 $n=10, d=5$のデータの場合 (全ステップ数$n-1$を入力する必要があります) <br>

使用例
    
    $ pyton calc_p_all.py data.csv data/sigma.csv data/xi.csv 9



`calc_p.py`, `calc_p_all.py`両方とも, 次のcsvファイルが`result`ディレクトリ以下に出力される
- naive_p
- selective_p
  
また, `cluster_result`ディレクトリ以下に次のものが出力される
- output (クラスタリング結果 scikit-learnの[linkage](https://docs.scipy.org/doc/scipy-0.14.0/reference/generated/scipy.cluster.hierarchy.linkage.html)の出力`Z`と同じフォーマット)

#### デモ (FPR)
`demo_FPR`ディレクトリにFPRの計算プログラムがある <br>
`calc_FPR.py`を実行することで, 計算される, <font color="red"> 注意: $n$をあまり大きくすると, 膨大な計算時間がかかります (経験的に$\sim 100$) </font> <br>


必要な引数
- 繰り返し回数 (ここでは1000回)
- $n$の最大サイズ: ($10, 20, \ldots, n$)
- 次元 $d$
- 並列オプション : 2以上であれば, 並列計算を指定したコア数


例1 ($n = 10, \ldots, 50$, $d = 5$, $1000$回の計算 並列化しない)

    $ python calc_FPR.py 1000 50 5 1

例2 (4コアで並列)

    $ python calc_FPR.py 1000 50 5 4


### クラスタ中心の各次元での検定
`each_dim`ディレクトリ以下にあるものが相当する
#### 1. Preprocessing
1. `data`というディレクトリを作成する 
2. `data`ディレクトリに使いたいデータを入れ, 分散($\Sigma$, $\xi^2$)推定用とp値計算用に分ける (目安 分散推定: p値計算 = 2 : 8)
3. `preprocess.py`によって`data`ディレクトリ以下に `Sigma.csv`, `Xi.csv`が出力されることを確認する <br>

使用例 `preprocess.py data_estimate.csv`


#### 2. クラスタリングと検定の実行

##### 一つの階層でのp値を計算したい場合
`execute.py`を実行する
以下の引数が必要である
- データファイルパス
- $\Sigma$のファイルパス
- $\xi^2$のファイルパス
- どの階層か

使用例: `$ pyton calc_p.py data.csv data/sigma.csv data/xi.csv 0` 

#### 全ての階層でのp値を計算したい場合
`execute_allstep.py`を実行する
必要な引数
- データファイルパス
- $\Sigma$のファイルパス
- $\xi^2$のファイルパス
- 全部で何階層か? (サンプルサイズ-1が全階層となる)
使用例 $n=10, d=5$のデータの場合 (全ステップ数$n-1$を入力する必要があります) <br>
`$ pyton calc_p.py data.csv data/sigma.csv data/xi.csv 9`

`calc_p.py`, `calc_p_all.py`両方とも, 次のcsvファイルが`result`ディレクトリ以下に出力される
- output (クラスタリング結果 scikit-learnのlinkageの出力`Z`と同じフォーマット)
p-valueは各csvファイルに次元数分のp値が記載される
- naive_p 
- selective_p

##### toyデータ
 `cluster`, `each_dim`ともに, `data`ディレクトリにtoyデータを入れておいたので, 使い方の確認に使ってください

### デンドログラムにp値を付与した図の表示
`dendrogram_p.py`の中にある`pv_dendrogram`の関数をimportして使ってください

## ディレクトリ

```
root/
    |- cluster/
          |- cluster_result/
          |- cpp_source/
          |- data/
          |- interval/
          |- stat/
          |- calc_p_all.py
          |- calc_p.py
          |- ...
    |- each_dim/
          |- cpp_source/
          |- data/
          |- result/...
          |- execute_allstep.py
          |- execute.py
          |- ...
    |- figs
    |- README.md
```
### 注意事項
#### dataについて
- data.csvは値のみのformatとしてください. そうでない場合, errorとなります.

