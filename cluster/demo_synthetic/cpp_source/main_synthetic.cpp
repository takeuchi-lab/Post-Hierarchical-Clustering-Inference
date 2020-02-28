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

int main(int argc, char* argv[]){
    if (argc >= 6){
        int epo = std::stod(argv[1]);
        int n = std::stod(argv[2]);
        int d = std::stod(argv[3]);
        int step = std::stod(argv[4]);
	    int threads = std::stod(argv[5]);
        double mu = std::stod(argv[6]);

        // データ作成
    	// std::chrono::system_clock::time_point  start, end; // 型は auto で可
        // start = std::chrono::system_clock::now();
        int half_n = n / 2;
        for (int epoch = 0; epoch < epo; epoch++){
            // std::cout << "------epoch: " << epoch << "------" << std::endl; 
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
            // std::cout << "------------------Clustering------------------" << std::endl;
            std::vector< std::vector<double> > result;
            std::vector< std::vector<ind> > cluster_head_vec, cluster_next_vec, cluster_tail_vec;
            std::vector< std::pair<int, int> > selected_c;
            std::tie(result, cluster_head_vec, cluster_next_vec, cluster_tail_vec, selected_c) = ward(data);
            // std::cout << "-----------------Post Clustering Inference-----------------" << std::endl;
            std::vector<VectorXd> final_interval;
            double chi2;
            MatrixXd Sigma = MatrixXd::Identity(n, n);
            double xi = 1;
            std::tie(chi2, final_interval) = PCI_cluster_ward_step(data, cluster_head_vec, cluster_next_vec, cluster_tail_vec, selected_c, Sigma, xi, step, threads);
            
            // 統計量 ファイル出力
            std::string fname1;
            if (mu > 0.0) fname1 = "stat/test_stat_d_epoch" + std::to_string(epo) + "_step" + std::to_string(step) + "_n" +  std::to_string(n) + "_d" + std::to_string(d) + "_mu" + format("%.1f" ,mu) + ".csv";
            else fname1 = "stat/test_stat_d_epoch" + std::to_string(epo) + "_step" + std::to_string(step) + "_n" +  std::to_string(n) + "_d" + std::to_string(d) + ".csv";
            
            std::ofstream Test_Stat;
            Test_Stat.open(fname1, std::ios::app);
            Test_Stat << chi2 << "," << d << std::endl;
            Test_Stat.close();

            // 区間 ファイル出力
            std::string fname2;
            if (mu > 0.0) fname2 = "interval/final_interval_epoch"+ std::to_string(epo) + "_step" + std::to_string(step) + "_n" +  std::to_string(n) + "_d" + std::to_string(d) + "_mu" + format("%.1f" ,mu) + ".csv";
            else fname2 = "interval/final_interval_epoch"+ std::to_string(epo) + "_step" + std::to_string(step) + "_n" +  std::to_string(n) + "_d" + std::to_string(d) + ".csv";
            std::ofstream Interval(fname2, std::ios::app);
            for (int i = 0; i < final_interval.size(); i++){
                if (i < final_interval.size() - 1){
                    Interval << pow(final_interval.at(i).array(), 2)(0) << "," << pow(final_interval.at(i).array(), 2)(1) << ",";
                }
                else if (i >= final_interval.size() - 1){
                    Interval << pow(final_interval.at(i).array(), 2)(0) << "," << pow(final_interval.at(i).array(), 2)(1);
                }
            }
            Interval << std::endl;
            Interval.close();
        }
	    // end = std::chrono::system_clock::now();
        // double elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end-start).count();
        // std::cout << elapsed << "[milisec]" << std::endl;

        // std::string fn = "time_step";
        // std::ofstream timefile(fn + "_d" + ".csv");
        // for (int i = 0; i < time_vec.size(); i++){
        //     timefile << time_vec.at(i) << ", ";
        // }
        // timefile.close();

    }
    else {
        std::cout << "Not corrected argument" << std::endl;
        return 0;
    }
}
