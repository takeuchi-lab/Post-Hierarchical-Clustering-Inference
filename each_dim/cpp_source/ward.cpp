#include <bits/stdc++.h>
#include <Eigen/Dense>
#include "ward.h"

using namespace Eigen;
using ind = Eigen::MatrixXd::Index;

std::tuple<std::vector<std::vector<double>>, std::vector<std::vector<ind>>, std::vector<std::vector<ind>>, std::vector<std::vector<ind>>, std::vector<std::pair<int, int>>>ward(Eigen::MatrixXd &data){

    int n = data.rows();
	int d = data.cols();

	// PCIに必要なものセット こんなに用意するとメモリパンクしそう
	std::vector<std::vector<ind>> cluster_head_vec;
	std::vector<std::vector<ind>> cluster_tail_vec;
	std::vector<std::vector<ind>> cluster_next_vec;
	// クラスタリング自体の出力 dendrogram, fclusterに必要
    std::vector< std::vector<double> > result;

	std::vector<ind> cluster_head(n);
	std::vector<ind> cluster_tail(n);
	for(int i = 0; i < n; ++i){
		cluster_head[i] = i;
		cluster_tail[i] = i;
	}
	std::vector<ind> cluster_next(n, -1);
	std::vector<ind> cluster_elems(n, 1); // クラスタ内の要素数
    
    std::vector<ind> cluster_inds(n);
    for (int i = 0; i < cluster_inds.size(); i++){
        cluster_inds.at(i) = i;
    }

	MatrixXd cluster_centroid = data;

	// クラスタ内の、クラスタ中心からクラスタ内要素までの距離の二乗和。
	// 初期状態では0
	std::vector<double> cluster_internalDist(n, 0);

	// ward_mat(i, j) で、「i番目とj番目のクラスタを統合したと仮定したときのWard距離」を表す
	MatrixXd ward_mat = MatrixXd::Zero(n, n);

	int merged_i, merged_j;

	// tステップ目で選ばれたクラスタのid格納
	std::vector<std::pair<int, int>> selected_c;

	for (int t = 0; t < n - 1; t++) {
		// このとき、クラスターとして有効なのは0番目からn-t-1番目
		int last_cluster_id = n - t - 1;
		
		// 距離計算
		double minWard = -1;
		int minWard_i, minWard_j;
		double minWardCluster_internalDist;
		for (int i = 0; i < last_cluster_id; i++) {
			for (int j = i + 1; j <= last_cluster_id; j++) {
				if(t == 0){
					ward_mat(i, j) = (data.row(i) - data.row(j)).squaredNorm() / 2;
				}
				double new_ward = ward_mat(i, j);

				if(minWard < 0 || new_ward < minWard){
					minWard = new_ward;
					minWardCluster_internalDist = new_ward + cluster_internalDist[i] + cluster_internalDist[j];
					minWard_i = i;
					minWard_j = j;
				}
			}
		}
		std::pair<int, int> ij;
		ij.first = minWard_i;
		ij.second = minWard_j;
		selected_c.push_back(ij);

		cluster_head_vec.push_back(cluster_head);
		cluster_next_vec.push_back(cluster_next);
		cluster_tail_vec.push_back(cluster_tail);

		// ----------
		// クラスタ統合
		// 新たなクラスタの情報を、minWard_i 側に上書きする形で作成
		// 
		// また、事前に被統合クラスタを一番後ろに移動しておく
		// ※tの値によって、末尾の要素は順次無視されていくため
		// ※このほうが効率がよい
		// ※統合先クラスタのIDがlast_cluster_idよりは小さいことを前提としたコードを書いています
		// ----------

		// 移動
		// ※cluster_nextは動かさないことに注意 → なんでだろう
		// minward_jのヘッドを最後に持っていく
		// minward_jのテイルを最後に持っていく
		// 要素も交換
		// クラスタ内距離交換
		// クラスタの中心も移動
		std::swap(cluster_head[minWard_j], cluster_head[last_cluster_id]);
		std::swap(cluster_tail[minWard_j], cluster_tail[last_cluster_id]);
		std::swap(cluster_elems[minWard_j], cluster_elems[last_cluster_id]);
		std::swap(cluster_internalDist[minWard_j], cluster_internalDist[last_cluster_id]);
        std::swap(cluster_inds.at(minWard_j), cluster_inds.at(last_cluster_id));
		MatrixXd tmp = cluster_centroid.row(minWard_j);
		cluster_centroid.row(minWard_j) = cluster_centroid.row(last_cluster_id);
		cluster_centroid.row(last_cluster_id) = tmp;

		// ここのうまいやり方募集
		for(int k = 0; k < minWard_j; ++k){
			std::swap(ward_mat(k, minWard_j), ward_mat(k, last_cluster_id));
		}
		for(int k = minWard_j+1; k < last_cluster_id; ++k){
			std::swap(ward_mat(minWard_j, k), ward_mat(k, last_cluster_id));
		}

		minWard_j = last_cluster_id;

		// 統合されたクラスタの中心
		MatrixXd minWard_centroid = (cluster_elems[minWard_i] * cluster_centroid.row(minWard_i) + cluster_elems[minWard_j] * cluster_centroid.row(minWard_j)) / (cluster_elems[minWard_i] + cluster_elems[minWard_j]);

		// クラスタを統合した結果を書き込む
		cluster_next[cluster_tail[minWard_i]] = cluster_head[minWard_j];
		cluster_tail[minWard_i] = cluster_tail[minWard_j];
		cluster_elems[minWard_i] += cluster_elems[minWard_j];
		cluster_internalDist[minWard_i] = minWardCluster_internalDist;
		cluster_centroid.row(minWard_i) = minWard_centroid;


        std::vector<double> out;
		if (cluster_inds[minWard_i] > cluster_inds[minWard_j]){
			out.push_back(cluster_inds[minWard_j]);
			out.push_back(cluster_inds[minWard_i]);	
		}
		else {
			out.push_back(cluster_inds[minWard_i]);
			out.push_back(cluster_inds[minWard_j]);
		}
		out.push_back(minWard);
        out.push_back(cluster_elems[minWard_i]);

        cluster_inds.at(minWard_i) = n + t;

		// Ward距離の行列の更新
		for(int j_tmp = 0; j_tmp < last_cluster_id; ++j_tmp){
			// ここは "< last_cluster_id" でよいことに注意
			// （last_cluster_id番目は被統合クラスタなので）

			if(minWard_i == j_tmp) continue;
			int i = minWard_i, j = j_tmp;
			if(i > j) std::swap(i, j); // assure i < j

			double new_ward = - cluster_internalDist[i] - cluster_internalDist[j];
			MatrixXd tmp_centroid = (cluster_elems[minWard_i] * minWard_centroid + cluster_elems[j_tmp] * cluster_centroid.row(j_tmp)) / (cluster_elems[minWard_i] + cluster_elems[j_tmp]);

			for(int ci = cluster_head[minWard_i]; ci != -1; ci = cluster_next[ci]){
				new_ward += (data.row(ci) - tmp_centroid).squaredNorm();
			}
			for(int cj = cluster_head[j_tmp]; cj != -1; cj = cluster_next[cj]){
				new_ward += (data.row(cj) - tmp_centroid).squaredNorm();
			}

			ward_mat(i, j) = new_ward;
		}
		result.push_back(out);
	}
    return {result, cluster_head_vec, cluster_next_vec, cluster_tail_vec, selected_c};
}
