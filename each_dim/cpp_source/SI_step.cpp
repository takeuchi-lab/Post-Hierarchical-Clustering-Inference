#include <bits/stdc++.h>
#include <Eigen/Core>
#include <fstream>

#include "SI_step.h"

using namespace Eigen;
using ind = Eigen::MatrixXd::Index;

// \bm{1}_k^{(t)}をつくる
std::vector<double> makeIndicator(std::vector<ind> &cluster_head, std::vector<ind> &cluster_next, int k, int n){
	std::vector<double> indicator(n, 0); 
	for(int ci = cluster_head.at(k); ci != -1; ci = cluster_next.at(ci)){
        indicator.at(ci) = 1;
        }
	return indicator;
}

// C_k^{(t)}
std::vector<int> makeClusterindSet(VectorXd &indecator){
    std::vector<int> c_ind_set;
    for (int i = 0; i < indecator.size(); i++){
        if (indecator(i) == 1) c_ind_set.push_back(i);
    }
    return c_ind_set;
}

// \bar{x}_k^{(t)}
VectorXd calcCentroid(MatrixXd &data, std::vector<int> &c_ind_set){
    int d = data.cols();
    VectorXd sumX = VectorXd::Zero(d);
    for (const int &i : c_ind_set){
        sumX += data.row(i);
    }

    VectorXd centroid(d);
    centroid = sumX / c_ind_set.size();
    return centroid;
}

// クラスタ内分散 ||x_i - \bar{x}_k ||^2
double IntraCVariance(MatrixXd &data, std::vector<int> &c_ind_set, VectorXd centroid){
    double energy = 0;
    for (const int &i : c_ind_set){
        energy += (data.row(i).transpose() - centroid).squaredNorm();
    }
    return energy;
}

// Ward法の距離計算 d(c_k, c_k')
double distWard(MatrixXd &data, std::vector<int> &c_k1_ind_set, std::vector<int> &c_k2_ind_set, VectorXd &k1_centroid, VectorXd &k2_centroid){
    double distance = 0.0;
    std::vector<int> c_k12_ind_set;
    for (const int &i : c_k1_ind_set){
        c_k12_ind_set.push_back(i);
    }
    for (const int &i : c_k2_ind_set){
        c_k12_ind_set.push_back(i);
    }
    VectorXd k12_centroid = (c_k1_ind_set.size() * k1_centroid + c_k2_ind_set.size() * k2_centroid) / (c_k1_ind_set.size() + c_k2_ind_set.size());
    distance = IntraCVariance(data, c_k12_ind_set, k12_centroid) - IntraCVariance(data, c_k1_ind_set, k1_centroid) - IntraCVariance(data, c_k2_ind_set, k2_centroid);
    return distance;
}

// vardelta 計算 関数化する意味ないかもしれない
// tステップ目のクラスタ数だけ要素数をもつSTLコンテナ
std::vector<double> calcVarDelta(MatrixXd &Sigma, VectorXd &delta, std::vector<std::vector<int>> &c_ind_set_vec){
    std::vector<double> vardelta_vec;
    VectorXd sigmadelta = Sigma * delta;
    for (int i = 0; i < c_ind_set_vec.size(); i++){
        double delta_k = 0.0;
        for (int j : c_ind_set_vec.at(i)){
            delta_k += sigmadelta(j);
        }
        delta_k = delta_k / c_ind_set_vec.at(i).size();
        vardelta_vec.push_back(delta_k);
    }
    return vardelta_vec;
}

//nC2を計算するために実装
int combination(int n, int r){
    if (r == 1) return n;
    if (n == r) return 1;

    return combination(n - 1, r - 1) + combination(n - 1, r);
}

//VectorXd型の条件にあうものだけを返す
VectorXd slice(VectorXd &vec, Array<bool, Dynamic, 1> &cond){
    std::vector<double> slicing;
    VectorXd slice_vec;
    for (double i = 0; i < cond.size(); i++){
        if (cond(i) != 0){
            slicing.push_back(vec.coeff(i));
        }
    }
    return Map<VectorXd>(slicing.data(), slicing.size());
}

// 標準正規分布のcdf
double norm_cdf(double x){
    double cdf = (1.0 + std::erf(x / sqrt(2.0))) / 2.0;
    return cdf;
}

// Eigen::Matrixのi列目の値でソート
MatrixXd matrix_icol_sort(MatrixXd &mat, int num){
    MatrixXd sorted(mat.rows(), mat.cols());
    VectorXd cols = mat.col(num);
    // argsort
    std::vector<size_t> indices(cols.size());
    std::iota(indices.begin(), indices.end(), 0);
    std::sort(indices.begin(), indices.end(), [&cols](size_t i1, size_t i2) {return cols(i1) < cols(i2);});
    for (int i = 0; i < indices.size(); i++){
        sorted.row(i) = mat.row(indices.at(i));
    }
    return sorted;
}

// 区間の共通部分を計算する
std::vector<VectorXd> intersection(double L, double U, MatrixXd &new_interval){
    std::vector<VectorXd> intersect;
    MatrixXd using_interval = new_interval;
    double lower = L;
    VectorXd interval(2);
    int step_num = 0;
    while (using_interval.size() > 0){
        VectorXd low = using_interval.col(0);
        VectorXd up = using_interval.col(1);
        VectorXd base = matrix_icol_sort(using_interval, 0).row(0);
        interval(0) = lower;
        interval(1) = base(0);
        intersect.push_back(interval);
        Array<bool, Dynamic, 1> inner_fg = low.array() < base(1);
        while (1){
            lower = slice(up, inner_fg).maxCoeff();
            if ((inner_fg.array() == (low.array() <= lower)).all()) break;
            else inner_fg = low.array() <= lower;
        }        
        Array<bool, Dynamic, 1> outer_fg = !inner_fg;
        double interval_size = outer_fg.cast<double>().matrix().transpose() * outer_fg.cast<double>().matrix();
        using_interval.resize(interval_size, 2);
        using_interval.col(0) = slice(low, outer_fg);
        using_interval.col(1) = slice(up, outer_fg);    
    }
    interval(0) = lower;
    interval(1) = U;
    intersect.push_back(interval);
    return intersect;
}

//複数区間の場合のp値計算
// intersection()が区間を昇順に出力することを利用する
double calc_multisec_p(double tau, double sig, std::vector<VectorXd> &final_interval){
    double denominator = 0.0;
    double numerator = 0.0;
    double tau_j = tau / sig;
    for (int i = final_interval.size() - 1; i > -1; i--){
        double l = final_interval.at(i).coeff(0) / sig;
        double u = final_interval.at(i).coeff(1) / sig;
        denominator += norm_cdf(-l) - norm_cdf(-u);
        // 区間がどの位置にあるか
        bool flag = (l <= tau_j) && (tau_j <= u);
        bool belong = false;
        if(flag){
            numerator += norm_cdf(-tau_j) - norm_cdf(-u);
            belong = true;
        }
        else if(!belong){
            numerator += norm_cdf(-l) - norm_cdf(-u);
        }
    }
    if (numerator / denominator < 0) std::cout << "In calc multi p-value under 0!!!!";
    if (numerator / denominator > 1) std::cout << "In calc multi p-value over 1!!!!";
    return numerator / denominator;
}

//Importance Sampling
double ImpSamp(double sig, double tau, double L, double U){
    int num_samp = pow(10, 6);
    double nume = 0.0;
    double deno = 0.0;
    std::random_device seed_gen;
    for (int i = 0; i < num_samp; i++){
        std::default_random_engine engine(seed_gen());
        std::normal_distribution<> dist(tau, sig);
        double x = dist(engine);
        if( (tau < x) && (x < U)){
            nume += exp(-tau * x / sig + pow(tau, 2)/ sig);
        }
        if((L < x) && (x < U)){
            deno += exp(-tau * x / sig + pow(tau, 2) / sig);
        } 
    }
    return nume / deno;
}

// PCI 各次元Ver 切断点保持
std::pair<std::vector<double>, std::vector<double>> PCI_ward_Degs_step(MatrixXd &data, std::vector<std::vector<ind>> &cluster_head_vec, std::vector<std::vector<ind>> &cluster_next_vec, std::vector<std::vector<ind>> &cluster_tail_vec, std::vector<std::pair<int, int>> &selected_c, MatrixXd &Sigma, MatrixXd &Xi, int step){
    int n = data.rows();
    int d = data.cols();
    assert(step >= 0 && step < n - 1);
    // 全ステップにおける区間 線形と下に凸だけ
    MatrixXd interval(d, 2);
    for (int i = 0; i < d; i++){
        Vector2d h(2);
        h << -std::numeric_limits<double>::infinity(), std::numeric_limits<double>::infinity();
        interval.row(i) = h;
    }
    // 上に凸だけ集める あとで共通部分をとる d * ? * 2
    std::vector<MatrixXd> section_mat(d);
    // 最終的な区間の計算に必要なもの
    std::vector<ind> cluster_head = cluster_head_vec.at(step);
    std::vector<ind> cluster_next = cluster_next_vec.at(step);
    std::vector<ind> cluster_tail = cluster_tail_vec.at(step);
    // Selection Eventの計算に必要なもの
    std::vector<VectorXd> indecator_vec;
    std::vector<std::vector<int>> c_ind_set_vec;
    for (int i = 0; i < n - step; i++){
        // 1_k
        VectorXd temp_indecator = Map<VectorXd> (makeIndicator(cluster_head,cluster_next, i, n).data(), makeIndicator(cluster_head, cluster_next, i, n).size());
        indecator_vec.push_back(temp_indecator);
        // C_k
        c_ind_set_vec.push_back(makeClusterindSet(temp_indecator));
    }
    // delta_ab
    int a = selected_c.at(step).first;
    int b = selected_c.at(step).second;
    VectorXd delta_ab = indecator_vec.at(a) / c_ind_set_vec.at(a).size() - indecator_vec.at(b) / c_ind_set_vec.at(b).size();

    // 検定統計量 \tau_ab
    VectorXd a_centroid = calcCentroid(data, c_ind_set_vec.at(a));
    VectorXd b_centroid = calcCentroid(data, c_ind_set_vec.at(b));
    VectorXd tau_ab = abs((a_centroid - b_centroid).array());
    VectorXd si = sign((a_centroid - b_centroid).array());
    VectorXd tilde_sigma2 = Xi.diagonal() * delta_ab.transpose() * Sigma * delta_ab;
    // メインループ
    for (int s = 0; s < step + 1; s++){
        std::cout << "-----------step: " << s << "---------------" << std::endl;
        std::vector<ind> cluster_head = cluster_head_vec.at(s);
        std::vector<ind> cluster_next = cluster_next_vec.at(s);
        std::vector<ind> cluster_tail = cluster_tail_vec.at(s);

        // Selection Eventの計算に必要なもの
        std::vector<VectorXd> indecator_vec;
        std::vector<std::vector<int>> c_ind_set_vec;
        
        for (int i = 0; i < n - s; i++){
            // 1_k
            VectorXd temp_indecator = Map<VectorXd> (makeIndicator(cluster_head,cluster_next, i, n).data(), makeIndicator(cluster_head, cluster_next, i, n).size());
            indecator_vec.push_back(temp_indecator);
            // C_k
            std::vector<int> hoge = makeClusterindSet(temp_indecator);
            c_ind_set_vec.push_back(makeClusterindSet(temp_indecator));
        }
        int a = selected_c.at(s).first;
        int b = selected_c.at(s).second;
        if (a > b){
            int temp = a;
            a = b;
            b = temp;
        }
        VectorXd a_centroid = calcCentroid(data, c_ind_set_vec.at(a));
        VectorXd b_centroid = calcCentroid(data, c_ind_set_vec.at(b));
        
        double d_ab = distWard(data, c_ind_set_vec.at(a), c_ind_set_vec.at(b), a_centroid, b_centroid);
        double w_ab = double(c_ind_set_vec.at(a).size() * c_ind_set_vec.at(b).size()) / double(c_ind_set_vec.at(a).size() + c_ind_set_vec.at(b).size());
        std::vector<double> varDelta = calcVarDelta(Sigma, delta_ab, c_ind_set_vec);
        std::vector<double> vBv;
        Matrix<double, Dynamic, Dynamic> cBv(combination(n - s, 2), d);
        Matrix<double, Dynamic, Dynamic> cBc(combination(n - s, 2), d);
        bool flag = false;
        int counter = 0;
        for (int k1 = 0; k1 < n - s - 1; k1++){
            for (int k2 = k1 + 1; k2 < n - s; k2++){
                VectorXd k1_centroid = calcCentroid(data, c_ind_set_vec.at(k1));
                VectorXd k2_centroid = calcCentroid(data, c_ind_set_vec.at(k2));
                if (d_ab - distWard(data, c_ind_set_vec.at(k1), c_ind_set_vec.at(k2), k1_centroid, k2_centroid) > 0){
                    flag = true;
                    break;
                }
                //vBv
                vBv.push_back(d_ab - distWard(data, c_ind_set_vec.at(k1), c_ind_set_vec.at(k2), k1_centroid, k2_centroid));
                //cBv
                double w_k12 = double(c_ind_set_vec.at(k1).size() * c_ind_set_vec.at(k2).size()) / double(c_ind_set_vec.at(k1).size() + c_ind_set_vec.at(k2).size());
                // \frac{s}{\tilde{\sigma}^2}
                VectorXd s_div_sigma = si.array() / tilde_sigma2.array();
                VectorXd d_vec = w_ab * (varDelta.at(a) - varDelta.at(b)) * (a_centroid - b_centroid) - (w_k12 * (varDelta.at(k1) - varDelta.at(k2))*(k1_centroid - k2_centroid));
                cBv.row(counter) = s_div_sigma.array() * (d_vec.transpose() * Xi).transpose().array();
                //cBc
                double scalar = w_ab * pow(varDelta.at(a) - varDelta.at(b), 2.0)  - w_k12 * pow(varDelta.at(k1) - varDelta.at(k2), 2.0);
                // 数値誤差処理
                if (std::abs(scalar) < pow(10, -15)) scalar = 0;
                // めちゃ値が小さくなってるときがある (10^-17 くらい)
                cBc.row(counter) = pow(s_div_sigma.array(), 2.0).array() * scalar * ((Xi.transpose() * Xi).diagonal().array());
                counter++;
            }
            if(flag) break;
        }
        if(flag){
            std::cout << "error occured!!!" << std::endl;
            break;
        }
        VectorXd vBv_vec = Map<VectorXd> (vBv.data(), vBv.size());
        // 区間計算
        VectorXd discriminant(combination(n - s, 2));
            
        for (int j = 0; j < d; j++){
            std::vector<double> forL;
            std::vector<double> forU;
            //コピーが発生しちゃう... 
            // VectorXd &h = cBc.col(j)ができない...
            VectorXd cBc_j = cBc.col(j);
            VectorXd cBv_j = cBv.col(j);
    
            discriminant = pow(cBv_j.array(), 2.0).array() - cBc_j.array() * vBv_vec.array();
            Array<bool, Dynamic, 1> cond1 = cBc_j.array() == 0;
            Array<bool, Dynamic, 1> cond1_l = (cBc_j.array() == 0) && (cBv_j.array() < 0.0);
            Array<bool, Dynamic, 1> cond1_u = (cBc_j.array() == 0) && (cBv_j.array() > 0.0);
            // 線形event
            if (cond1.any() == 1){
                if (cond1_l.any() == 1){
                    forL.push_back((-slice(vBv_vec, cond1_l).array() / (2 * slice(cBv_j, cond1_l)).array()).maxCoeff());
                }
                if (cond1_u.any() == 1){
                    forU.push_back((-slice(vBv_vec, cond1_u).array() / (2 * slice(cBv_j, cond1_u)).array()).minCoeff());
                }
            }
            // 下に凸
            Array<bool, Dynamic, 1> cond2 = cBc_j.array() > 0.0;
            if (cond2.any() == 1){
                forL.push_back(((-slice(cBv_j, cond2).array() - sqrt(slice(discriminant, cond2).array())).array() / slice(cBc_j, cond2).array()).maxCoeff());
                forU.push_back(((-slice(cBv_j, cond2).array() + sqrt(slice(discriminant, cond2).array())).array() / slice(cBc_j, cond2).array()).minCoeff());
            }
            // 線形と下に凸で共通部分をとる
            double maxforL, minforU;
            if (forL.size() > 0) maxforL = *std::max_element(forL.begin(), forL.end());
            else maxforL = -std::numeric_limits<double>::infinity();
        
            if (forU.size() > 0) minforU = *std::min_element(forU.begin(), forU.end());
            else minforU = std::numeric_limits<double>::infinity();

            maxforL = std::max(maxforL, -tau_ab(j));
            if (interval(j, 0) < maxforL + tau_ab(j)) {
                    interval(j, 0) = maxforL + tau_ab(j);
            }
            if (interval(j, 1) > minforU + tau_ab(j)) {
                    interval(j, 1) = minforU + tau_ab(j);
            }
            // 上に凸
            Array<bool, Dynamic, 1> cond3 = (cBc_j.array() < 0.0) && (discriminant.array() > 0.0);
            if (cond3.any() == 1){
                VectorXd leftside = (-slice(cBv_j, cond3).array() + sqrt(slice(discriminant, cond3).array())).array() / slice(cBc_j, cond3).array();
                VectorXd rightside = (-slice(cBv_j, cond3).array() - sqrt(slice(discriminant, cond3).array())).array() / slice(cBc_j, cond3).array();
                MatrixXd new_interval(leftside.size(), 2);
                new_interval.col(0) = leftside.array() + tau_ab(j);
                new_interval.col(1) = rightside.array() + tau_ab(j);
                if (section_mat.at(j).size() > 0){
                    MatrixXd temp_mat = section_mat.at(j);
                    section_mat.at(j).resize(temp_mat.rows() + new_interval.rows(), 2);
                    section_mat.at(j).block(0, 0, temp_mat.rows(), 2) = temp_mat;
                    section_mat.at(j).block(temp_mat.rows(), 0, new_interval.rows(), 2) = new_interval;
                }
                else {
                    section_mat.at(j) = new_interval;
                }
            }
        }
    }
    // 最終的な区間作成 下に凸で区間にすっぽり入ってしまうものと共通部分をとる
    std::vector<double> naive_p(d, 0);
    std::vector<double> selective_p(d, 0);
    for (int j = 0; j < d; j++){
        double L = interval(j, 0);
        double U = interval(j, 1);
        double sig = sqrt(tilde_sigma2.coeff(j));
        if (section_mat.at(j).size() > 0){
            VectorXd x_small = section_mat.at(j).col(0);
            VectorXd x_large = section_mat.at(j).col(1);
            // x_L < L < x_U
            Array<bool, Dynamic, 1> flag1 = (x_small.array() < L) && (L < x_large.array());
            while (flag1.any() == 1){
                    L = slice(x_large, flag1).maxCoeff();
                    flag1 = (x_small.array() < L) && (L < x_large.array());
            }   
            // x_L < U < x_U
            Array<bool, Dynamic, 1> flag2 = (x_small.array() < U) && (U < x_large.array());
            while (flag2.any() == 1){
                U = slice(x_small, flag2).minCoeff();
                flag2 = (x_small.array() < U) && (U < x_large.array());
            }      
            // L < x_L < x_U < U
            Array<bool, Dynamic, 1> flag3 = (L < x_small.array()) && (x_large.array() < U);
            if (flag3.any() == 1){
                MatrixXd new_interval(slice(x_small, flag3).size(), 2);
                new_interval.col(0) = slice(x_small, flag3);
                new_interval.col(1) = slice(x_large, flag3);
                std::vector<VectorXd> final_interval = intersection(L, U, new_interval);
                
                selective_p.at(j) = calc_multisec_p(tau_ab(j), sig, final_interval);
            }
            else {
                if ((norm_cdf(- tau_ab(j) / sig) - norm_cdf(-U / sig)) / (norm_cdf(-L / sig) - norm_cdf(-U / sig)) > 1) std::cout << "One interval p-value over 1!!!" << std::endl;
                if ((norm_cdf(- tau_ab(j) / sig) - norm_cdf(-U / sig)) / (norm_cdf(-L / sig) - norm_cdf(-U / sig)) < 0) std::cout << "One interval p-value under 0!!!" << std::endl;
                if ((norm_cdf(-L / sig) - norm_cdf(-U / sig)) == 0.000000000000000){
                    selective_p.at(j) = ImpSamp(sig, tau_ab(j), L, U);
                }
                else {
                selective_p.at(j) = (norm_cdf(- tau_ab(j) / sig) - norm_cdf(-U / sig)) / (norm_cdf(-L / sig) - norm_cdf(-U / sig));
                }
            }
        } 
        else {
            if ((norm_cdf(- tau_ab(j) / sig) - norm_cdf(-U / sig)) / (norm_cdf(-L / sig) - norm_cdf(-U / sig)) > 1) std::cout << "One interval p-value over 1!!!" << std::endl;
            if ((norm_cdf(- tau_ab(j) / sig) - norm_cdf(-U / sig)) / (norm_cdf(-L / sig) - norm_cdf(-U / sig)) < 0) std::cout << "One interval p-value under 0!!!" << std::endl;
            if ((norm_cdf(-L / sig) - norm_cdf(-U / sig)) == 0.000000000000000){
                    selective_p.at(j) = ImpSamp(sig, tau_ab(j), L, U);
                }
            else {
                selective_p.at(j) = (norm_cdf(- tau_ab(j) / sig) - norm_cdf(-U / sig)) / (norm_cdf(-L / sig) - norm_cdf(-U / sig));
            }
        }
        double sub1 = norm_cdf(tau_ab(j) / sig);
        double sub2 = norm_cdf(-tau_ab(j) / sig); 
        naive_p.at(j) = 2 * std::min(sub1, sub2); 
    }
    return {naive_p, selective_p};    
}