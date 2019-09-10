#include <bits/stdc++.h>
#include <Eigen/Core>
#include <fstream>

#include "ward.h"
#include "SI_step.h"

using namespace Eigen;

std::vector<double> split(std::string& input, char delimiter){
    std::istringstream stream(input);
    std::string field;
    std::vector<double> result;
    while (getline(stream, field, delimiter)) {
        result.push_back(std::stod(field));
    }
    return result;
}

MatrixXd load_csv(std::string fn){
    std::ifstream ifs(fn);
    
    std::string line;
    std::vector<std::vector<double>> data;
    while (getline(ifs, line)) {
        std::vector<double> numvec = split(line, ',');
        data.push_back(numvec);
    }
    MatrixXd A(data.size(), data[0].size());
    for (int i = 0; i < data.size(); i++){
        A.row(i) = Map<VectorXd> (data.at(i).data(), data.at(i).size());
    }
    return A;
}

int main(int argc, char ** argv){
    if (argc < 4){
        std::cout << "plese input below" << std::endl;
        std::cout << "./pci_cluster_ex datafile sigmafile xivalue step"
        std::cout << "example: ./pci_cluster_ex data/data.csv data/sigma.csv 0.9262638032548505 4" << std::endl;
        return 0;
    }
    // データ作成
    std::vector<double> time_vec;
    std::string datafile = argv[1];
    MatrixXd data = load_csv(datafile);
    int n = data.rows();
    int d = data.cols();
    std::cout << "n: " << n << std::endl;
    std::cout << "d: " << d << std::endl;

    // クラスタリング    
    std::chrono::system_clock::time_point  start, end; // 型は auto で可
    start = std::chrono::system_clock::now();
    std::vector< std::vector<double> > result;
    std::vector< std::vector<ind> > cluster_head_vec, cluster_next_vec, cluster_tail_vec;
    std::vector< std::pair<int, int> > selected_c;
    std::tie(result, cluster_head_vec, cluster_next_vec, cluster_tail_vec, selected_c) = ward(data);
    
    //検定
    std::string sigmafile = argv[2];
    // std::string xifile = argv[3];
    MatrixXd Sigma = load_csv(sigmafile);
    Sigma = MatrixXd::Identity(n, n);
    double xi = std::stod(argv[3]);
    // double xi = 0.9262638032548505;
    int step = std::stoi(argv[4]);
    std::vector<VectorXd> final_interval;
    double chi2;
    std::tie(chi2, final_interval) = PCI_cluster_ward_step(data, cluster_head_vec, cluster_next_vec, cluster_tail_vec, selected_c, Sigma, xi, step);
    end = std::chrono::system_clock::now();
    double elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end-start).count();
    std::cout << elapsed << "[milisec]" << std::endl;
    time_vec.push_back(elapsed);

    //  クラスタリング結果出力
    std::string fname1 = "output";
    std::ofstream output("./cluster_result/" + fname1 + ".csv");
    for (int i = 0; i < n - 1; i++){
        for (int j = 0; j < result.at(i).size(); j++){
            output << result.at(i).at(j) << ",";
        }
        output << std::endl;
    }
    output.close();

    // 統計量 ファイル出力
    std::string fname2 = "stat/test_stat_step" + std::to_string(step) + ".csv";
    std::ofstream Test_Stat;
    Test_Stat.open(fname2, std::ios::app);
    Test_Stat << chi2 << "," << d << std::endl;
    Test_Stat.close();

    // 区間 ファイル出力
    std::string fname3 = "interval/final_interval_step" + std::to_string(step) + ".csv";
    std::ofstream Interval(fname3, std::ios::app);
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