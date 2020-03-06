#include <bits/stdc++.h>
#include <Eigen/Core>
#include <fstream>

#include "ward.h"
#include "SI_step.h"

using namespace Eigen;

// 小数点数を指定してstring型に変える
template <typename ... Args>
std::string format(const std::string& fmt, Args ... args )
{
    size_t len = std::snprintf( nullptr, 0, fmt.c_str(), args ... );
    std::vector<char> buf(len + 1);
    std::snprintf(&buf[0], len + 1, fmt.c_str(), args ... );
    return std::string(&buf[0], &buf[0] + len);
}


int main(int argc, char ** argv){
    int epo = std::stod(argv[1]);
    int n = std::stod(argv[2]);
    int d = std::stod(argv[3]);
    int step = std::stod(argv[4]);
    int threads = std::stod(argv[5]);
    double mu = std::stod(argv[6]);

    std::vector<double> time_vec;
    int half_n = n / 2;
    for (int epoch = 0; epoch < epo; epoch++){

        // データ作成
        std::mt19937 seed_gen(epoch); 
        std::default_random_engine engine(seed_gen());
        std::normal_distribution<> dist_mu(mu, 1.0);
        MatrixXd data_mu(half_n, d);
        for (int i = 0; i < half_n; i++){
            for (int j = 0; j < d; j++){
                data_mu(i, j) = dist_mu(engine);
            }
        }
        std::normal_distribution<> dist(0, 1.0);
        MatrixXd data_zero(half_n, d);
        for (int i = 0; i < half_n; i++){
            for (int j = 0; j < d; j++){
                data_zero(i, j) = dist(engine);
            }
        }
        MatrixXd data(n, d);
        for (int i = 0; i < half_n; i++){
            data.row(i) = data_mu.row(i);
        }
        for (int i = half_n; i < n; i++){
            data.row(i) = data_zero.row(i - half_n);
        }

        // クラスタリング    
        std::chrono::system_clock::time_point  start, end; // 型は auto で可
        start = std::chrono::system_clock::now();
        std::vector< std::vector<double> > result;
	    std::vector< std::vector<ind> > cluster_head_vec, cluster_next_vec, cluster_tail_vec;
        std::vector< std::pair<int, int> > selected_c;
        std::tie(result, cluster_head_vec, cluster_next_vec, cluster_tail_vec, selected_c) = ward(data);
    
        // 検定
        std::vector<double> naive_p, selective_p;
        MatrixXd Sigma = MatrixXd::Identity(n, n);
        MatrixXd Xi = MatrixXd::Identity(d, d);
        std::tie(naive_p, selective_p) = PCI_ward_Degs_step(data, cluster_head_vec, cluster_next_vec, cluster_tail_vec, selected_c, Sigma, Xi, step, threads);
        end = std::chrono::system_clock::now();
        double elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end-start).count();
        // std::cout << elapsed << "[milisec]" << std::endl;
        

        // output ファイル出力
        std::string fname1 = "output";
        std::ofstream output("./result/" + fname1 + ".csv");
        for (int i = 0; i < n - 1; i++){
            for (int j = 0; j < result.at(i).size() - 1; j++){
                output << result.at(i).at(j) << ",";
            }
            output << result.at(i).at(result.at(i).size() - 1) << std::endl;
        }
        output.close();

        // selective_p ファイル出力
        std::string fname2;
        if (mu > 0.0) fname2 = "./result/selective_p_n" + std::to_string(n) + "_d" + std::to_string(d) + "_step" + std::to_string(step) + "_mu" + format("%.1f", mu) + ".csv";
        else fname2 = "./result/selective_p_n" + std::to_string(n) + "_d" + std::to_string(d) + "_step" + std::to_string(step) + ".csv";
        std::ofstream Sep(fname2, std::ios::app);
        for (int j = 0; j < d-1; j++){
            Sep << selective_p.at(j) << ",";
        }
        Sep << selective_p.at(d - 1) << std::endl;
        Sep.close();
        
        // naive_p ファイル出力
        std::string fname3;
        if (mu > 0.0) fname3 = "./result/naive_p_n" + std::to_string(n) + "_d" + std::to_string(d) + "_step" + std::to_string(step) + "_mu" + format("%.1f", mu) + ".csv";
        else fname3 = "./result/naive_p_n" + std::to_string(n) + "_d" + std::to_string(d) + "_step" + std::to_string(step) + ".csv";
        std::ofstream Nap(fname3, std::ios::app);
        for (int j = 0; j < d - 1; j++){
            Nap << naive_p.at(j) << ",";
        }
        Nap << naive_p.at(d - 1) << std::endl;    
        Nap.close();
    }
}