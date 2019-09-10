#include <bits/stdc++.h>
#include <Eigen/Core>

using namespace Eigen;
using ind = Eigen::MatrixXd::Index;

// \bm{1}_k^{(t)}をつくる
std::vector<double> makeIndicator(std::vector<ind> &cluster_head, std::vector<ind> &cluster_next, int k, int n);
// C_k^{(t)}
std::vector<int> makeClusterindSet(VectorXd &indecator);
// \bar{x}_k^{(t)}
VectorXd calcCentroid(MatrixXd &data, std::vector<int> &c_ind_set);
// クラスタ内分散 ||x_i - \bar{x}_k ||^2
double IntraCVariance(MatrixXd &data, std::vector<int> &c_ind_set, VectorXd centroid);
// Ward法の距離計算 d(c_k, c_k')
double distWard(MatrixXd &data, std::vector<int> &c_k1_ind_set, std::vector<int> &c_k2_ind_set, VectorXd &k1_centroid, VectorXd &k2_centroid);
// vardelta 計算 関数化する意味ないかもしれない
// tステップ目のクラスタ数だけ要素数をもつSTLコンテナ
std::vector<double> calcVarDelta(MatrixXd &Sigma, VectorXd &delta, std::vector<std::vector<ind>> &c_ind_set_vec);
//nC2を計算するために実装
int combination(int n, int r);
//VectorXd型の条件にあうものだけを返す
VectorXd slice(VectorXd &vec, Array<bool, Dynamic, 1> &cond);
// 標準正規分布のcdf
double norm_cdf(double x);
// Eigen::Matrixのi列目の値でソート
MatrixXd matrix_icol_sort(MatrixXd &mat, int num);
// 区間の共通部分を計算する
std::vector<VectorXd> intersection(double L, double U, MatrixXd &new_interval);
//複数区間の場合のp値計算
// intersection()が区間を昇順に出力することを利用する
double calc_multisec_p(double tau, double sig, std::vector<VectorXd> &final_interval);
// PCI cluster
std::pair<double, std::vector<VectorXd>> PCI_cluster_ward_step(MatrixXd &data, std::vector<std::vector<ind>> &cluster_head_vec, std::vector<std::vector<ind>> &cluster_next_vec, std::vector<std::vector<ind>> &cluster_tail_vec, std::vector<std::pair<int, int>> &selected_c, MatrixXd &Sigma, double xi, int step);