#include <bits/stdc++.h>
#include <Eigen/Core>
#include <fstream>

#include "ward.h"
#include "SI_step.h"

using namespace Eigen;

int main(){
    // データ作成
    std::vector<double> time_vec;
    for (int epoch = 0; epoch < 1; epoch++){
        std::cout << "------epoch: " << epoch << "------" << std::endl; 
        std::mt19937 seed_gen(epoch); 
        std::default_random_engine engine(seed_gen());
        std::normal_distribution<> dist(0.0, 1.0);
        int n = 155;
        int d = 30;
        std::cout << "n: " << n << std::endl;
        std::cout << "d: " << d << std::endl;
        MatrixXd data = MatrixXd::Random(n, d);
        // クラスタリング    
        // std::cout << "------------------Clustering------------------" << std::endl;
        std::chrono::system_clock::time_point  start, end; // 型は auto で可
        start = std::chrono::system_clock::now();
        std::vector< std::vector<double> > result;
	    std::vector< std::vector<ind> > cluster_head_vec, cluster_next_vec, cluster_tail_vec;
        std::vector< std::pair<int, int> > selected_c;
        std::tie(result, cluster_head_vec, cluster_next_vec, cluster_tail_vec, selected_c) = ward(data);
        // output ファイル出力
        // std::string fname = "output";
        // std::ofstream output(fname + "_n" + std::to_string(n) + "_d" + std::to_string(d) + ".csv");
        // for (int i = 0; i < n - 1; i++){
        //     for (int j = 0; j < 3; j++){
        //         output << result.at(i).at(j) << ",";
        //     }
        //     output << std::endl;
        // }
        // output.close();
        // result 出力
        // std::cout << "Clustering Result" << std::endl;
        // for (int i = 0; i < result.size(); i++) {
        //     for(int j = 0; j < result[i].size(); ++j){
        //         std::cout << result[i][j] << " ";
        //         }
        //     std::cout << std::endl;
        // }
    // std::cout << "-----------------Post Clustering Inference-----------------" << std::endl;
        std::vector<double> naive_p, selective_p;
        MatrixXd Sigma = MatrixXd::Identity(n, n);
        MatrixXd Xi = MatrixXd::Identity(d, d);
        int step = n - 2;
        std::tie(naive_p, selective_p) = PCI_ward_Degs_step(data, cluster_head_vec, cluster_next_vec, cluster_tail_vec, selected_c, Sigma, Xi, step);
        end = std::chrono::system_clock::now();
        double elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end-start).count();
        std::cout << elapsed << "[milisec]" << std::endl;
        time_vec.push_back(elapsed);
        for (int i = 0; i < selective_p.size(); i++){
            std::cout << selective_p.at(i) << std::endl;
        }
    }
    // std::string fn = "time_step";
    // std::ofstream timefile(fn + "_d" + ".csv");
    // for (int i = 0; i < time_vec.size(); i++){
    //     timefile << time_vec.at(i) << ", ";
    // }
    // timefile.close();
}