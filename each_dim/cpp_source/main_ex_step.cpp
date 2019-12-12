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
        std::cout << "./pci_cluster_step datafile Sigmafile Xifile step";
        std::cout << "example: ./pci_cluster_step_ex data/data.csv data/Sigma.csv data/Xi.csv 4" << std::endl;
        return 0;
    }
    // データ作成
    std::string datafile = argv[1];
    MatrixXd data = load_csv(datafile);
    int n = data.rows();
    int d = data.cols();

    // クラスタリング    
    std::chrono::system_clock::time_point  start, end; // 型は auto で可
    start = std::chrono::system_clock::now();
    std::vector< std::vector<double> > result;
    std::vector< std::vector<ind> > cluster_head_vec, cluster_next_vec, cluster_tail_vec;
    std::vector< std::pair<int, int> > selected_c;
    std::tie(result, cluster_head_vec, cluster_next_vec, cluster_tail_vec, selected_c) = ward(data);
    
    std::vector<double> naive_p, selective_p;
    std::string sigmafile = argv[2];
    std::string xifile = argv[3];
    MatrixXd Sigma = load_csv(sigmafile);
    MatrixXd Xi = load_csv(xifile);
    int step = std::stoi(argv[4]);
    std::tie(naive_p, selective_p) = PCI_ward_Degs_step(data, cluster_head_vec, cluster_next_vec, cluster_tail_vec, selected_c, Sigma, Xi, step);
    end = std::chrono::system_clock::now();
    double elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end-start).count();
    std::cout << elapsed << "[milisec]" << std::endl;
    
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
    std::string fname2 = "selective_p";
    std::ofstream csvfile("./result/" + fname2 + "_step" + std::to_string(step) + ".csv");
    for (int j = 0; j < d-1; j++){
        csvfile << selective_p.at(j) << ",";
    }
    csvfile << selective_p.at(d - 1);
    // csvfile << std::endl;
    csvfile.close();

        // naive_p ファイル出力
    std::string fname3 = "naive_p";
    std::ofstream naivefile("./result/" + fname3 + "_step" + std::to_string(step) + ".csv");
    for (int j = 0; j < d - 1; j++){
        naivefile << naive_p.at(j) << ",";
    }
    naivefile << naive_p.at(d - 1);    
    naivefile.close();
}