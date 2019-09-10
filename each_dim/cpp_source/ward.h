#include<bits/stdc++.h>
#include<Eigen/Core>

using namespace Eigen;
using ind = Eigen::MatrixXd::Index;

std::tuple<std::vector<std::vector<double>>, std::vector<std::vector<ind>>, std::vector<std::vector<ind>>, std::vector<std::vector<ind>>, std::vector<std::pair<int, int>>>ward(Eigen::MatrixXd &data);