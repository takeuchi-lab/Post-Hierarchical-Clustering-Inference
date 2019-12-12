#include <bits/stdc++.h>
#include <Eigen/Core>

// #include <fstream>

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

// PCI クラスタ自体の検定
std::pair<double, std::vector<VectorXd>> PCI_cluster_ward_step(MatrixXd &data, std::vector<std::vector<ind>> &cluster_head_vec, std::vector<std::vector<ind>> &cluster_next_vec, std::vector<std::vector<ind>> &cluster_tail_vec, std::vector<std::pair<int, int>> &selected_c, MatrixXd &Sigma, double xi, int step){
    int n = data.rows();
    int d = data.cols();
    assert(step >= 0 && step < n - 1);
    // 全ステップにおける区間 線形と下に凸だけ
    VectorXd interval(2);
    interval << -std::numeric_limits<double>::infinity(), std::numeric_limits<double>::infinity();
    // 上に凸だけ集める あとで共通部分をとる ? * 2
    MatrixXd section_mat;
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
    // \tau_ab
    VectorXd a_centroid = calcCentroid(data, c_ind_set_vec.at(a));
    VectorXd b_centroid = calcCentroid(data, c_ind_set_vec.at(b));
    VectorXd tau_ab = a_centroid - b_centroid;

    // 検定統計量 \chi
    double delta2 = delta_ab.transpose() * delta_ab;
    double delsig = (delta_ab.transpose() * Sigma) * delta_ab;
    double check_sig = pow(xi, 2) * delsig / delta2;
    VectorXd Z = tau_ab.array() / (xi * sqrt(delsig));
    double chi2 = Z.transpose() * Z;
    double tau2 = tau_ab.transpose() * tau_ab;
    double P_vec_2 = sqrt(tau2 / delta2);
    double chi = sqrt(chi2);

    // メインループ
    for (int s = 0; s < step + 1; s++){
        // std::cout << "-----------step: " << s << "---------------" << std::endl;
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

        // 1_ab, \bar{x}_{ab}
        VectorXd indecator_ab = indecator_vec.at(a) + indecator_vec.at(b);
        VectorXd ab_centroid = (c_ind_set_vec.at(a).size() * a_centroid + c_ind_set_vec.at(b).size() * b_centroid) / (c_ind_set_vec.at(a).size() + c_ind_set_vec.at(b).size());
        
        double d_ab = distWard(data, c_ind_set_vec.at(a), c_ind_set_vec.at(b), a_centroid, b_centroid);
        // std::cout << "d_ab: " << d_ab << std::endl;
        std::vector<double> vBv;
        std::vector<double> wBv;
        std::vector<double> wBw;
        bool flag = false;
        int counter = 0;
        for (int k1 = 0; k1 < n - s - 1; k1++){
            for (int k2 = k1 + 1; k2 < n - s; k2++){
                if ((k1 == a) && (k2 == b)){
                    continue;
                }
                VectorXd k1_centroid = calcCentroid(data, c_ind_set_vec.at(k1));
                VectorXd k2_centroid = calcCentroid(data, c_ind_set_vec.at(k2));
                if (d_ab - distWard(data, c_ind_set_vec.at(k1), c_ind_set_vec.at(k2), k1_centroid, k2_centroid) > 0){
                    flag = true;
                    std::cout << "!!!! d_ab > d_kk' !!!!" << std::endl;
                    break;
                }
                // 1_kk', \bar{x}_kk'
                VectorXd indecator_k12 = indecator_vec.at(k1) + indecator_vec.at(k2);
                VectorXd k12_centroid = (c_ind_set_vec.at(k1).size() * k1_centroid + c_ind_set_vec.at(k2).size() * k2_centroid) / (c_ind_set_vec.at(k1).size() + c_ind_set_vec.at(k2).size());

                //vBv
                vBv.push_back(d_ab - distWard(data, c_ind_set_vec.at(k1), c_ind_set_vec.at(k2), k1_centroid, k2_centroid));
                //wBv
                MatrixXd mat = -indecator_ab * ab_centroid.transpose() + indecator_vec.at(a) * a_centroid.transpose() + indecator_vec.at(b) * b_centroid.transpose()
                             + indecator_k12 * k12_centroid.transpose() - indecator_vec.at(k1) * k1_centroid.transpose() - indecator_vec.at(k2) * k2_centroid.transpose();
                double scalar = delta_ab.transpose() * mat * tau_ab;
                wBv.push_back(scalar / (P_vec_2 * delta2));

                double nume_abab = delta_ab.transpose() * indecator_ab * indecator_ab.transpose() * delta_ab;
                double abab = nume_abab / (c_ind_set_vec.at(a).size() + c_ind_set_vec.at(b).size());
                double nume_aa = delta_ab.transpose() * indecator_vec.at(a) * indecator_vec.at(a).transpose() * delta_ab;
                double aa = nume_aa / c_ind_set_vec.at(a).size();
                double nume_bb = delta_ab.transpose() * indecator_vec.at(b) * indecator_vec.at(b).transpose() * delta_ab; 
                double bb = nume_bb / c_ind_set_vec.at(b).size();
                double nume_k1212 = delta_ab.transpose() * indecator_k12 * indecator_k12.transpose() * delta_ab;
                double k1212 = nume_k1212 / (c_ind_set_vec.at(k1).size() + c_ind_set_vec.at(k2).size());
                double nume_k11 = delta_ab.transpose() * indecator_vec.at(k1) * indecator_vec.at(k1).transpose() * delta_ab;
                double k11 = nume_k11 / c_ind_set_vec.at(k1).size();
                double nume_k22 = delta_ab.transpose() * indecator_vec.at(k2) * indecator_vec.at(k2).transpose() * delta_ab; 
                double k22 = nume_k22 / c_ind_set_vec.at(k2).size();
                double kappa = -abab + aa + bb + k1212 - k11 - k22;
                wBw.push_back(kappa / delta2);
            }
            if(flag) break;
        }
        if(flag){
            std::cout << "error occured!!!" << std::endl;
            break;
        }
        VectorXd vBv_vec = Map<VectorXd> (vBv.data(), vBv.size());
        VectorXd wBv_vec = Map<VectorXd> (wBv.data(), wBv.size());
        VectorXd wBw_vec = Map<VectorXd> (wBw.data(), wBw.size());

        // 区間計算
        VectorXd discriminant(combination(n - s, 2));
            
        std::vector<double> forL;
        std::vector<double> forU;

        VectorXd alpha = pow(check_sig, 2) * wBw_vec.array();
        VectorXd beta = 2 * check_sig * (wBv_vec.array() - check_sig * chi * wBw_vec.array());
        VectorXd gamma = vBv_vec.array() - 2 * chi * check_sig * wBv_vec.array() + pow(check_sig, 2) * chi2 * wBw_vec.array();

        discriminant = pow(beta.array() , 2) - 4 * alpha.array() * gamma.array();
        
        Array<bool, Dynamic, 1> cond1 = alpha.array() == 0;
        Array<bool, Dynamic, 1> cond1_l = (alpha.array() == 0) && (beta.array() < 0.0);
        Array<bool, Dynamic, 1> cond1_u = (alpha.array() == 0) && (beta.array() > 0.0);
        // 線形event
        if (cond1.any() == 1){
            if (cond1_l.any() == 1){
                // 0とmaxとる必要があるか
                double forL_1 = ((-slice(vBv_vec,cond1_l).array() / (2 * check_sig * slice(wBv_vec, cond1_l).array())).array() + chi).maxCoeff();
                forL.push_back(std::max(forL_1, 0.0));
            }
            if (cond1_u.any() == 1){
                forU.push_back(((-slice(vBv_vec,cond1_u).array() / (2 * check_sig * slice(wBv_vec, cond1_u).array())).array() + chi).minCoeff());
                }
        }
        // 下に凸
        Array<bool, Dynamic, 1> cond2 = wBw_vec.array() > 0.0;
        if (cond2.any() == 1){
            double forL_2 = ((-slice(beta, cond2).array() - sqrt(slice(discriminant, cond2).array())).array() / (2 * slice(alpha, cond2).array())).maxCoeff();
            forL.push_back(std::max(forL_2, 0.0));
            forU.push_back(((-slice(beta, cond2).array() + sqrt(slice(discriminant, cond2).array())).array() / (2 * slice(alpha, cond2).array())).minCoeff());
            }
        // 線形と下に凸で共通部分をとる
        double maxforL, minforU;
        if (forL.size() > 0) maxforL = *std::max_element(forL.begin(), forL.end());
        else maxforL = 0.0;
    
        if (forU.size() > 0) minforU = *std::min_element(forU.begin(), forU.end());
        else minforU = std::numeric_limits<double>::infinity();

        if (interval(0) < maxforL) {
                interval(0) = maxforL;
        }
        if (interval(1) > minforU) {
                interval(1) = minforU;
        }
        // 上に凸
        Array<bool, Dynamic, 1> cond3 = (wBw_vec.array() < 0.0) && (discriminant.array() > 0.0);
        if (cond3.any() == 1){
            VectorXd leftside = (-slice(beta, cond3).array() + sqrt(slice(discriminant, cond3).array())).array() / (2 * slice(alpha, cond3).array());
            VectorXd rightside = (-slice(beta, cond3).array() - sqrt(slice(discriminant, cond3).array())).array() / (2 * slice(alpha, cond3).array());
            MatrixXd new_interval(leftside.size(), 2);

            new_interval.col(0) = leftside;
            new_interval.col(1) = rightside;
            if (section_mat.size() > 0){
                MatrixXd temp_mat = section_mat;
                section_mat.resize(temp_mat.rows() + new_interval.rows(), 2);
                section_mat.block(0, 0, temp_mat.rows(), 2) = temp_mat;
                section_mat.block(temp_mat.rows(), 0, new_interval.rows(), 2) = new_interval;
            }
            else {
                section_mat = new_interval;
            }
        }    
    }
    // 最終的な区間作成 下に凸で区間にすっぽり入ってしまうものと共通部分をとる
    double L = interval(0);
    double U = interval(1);
    double sig = sqrt(check_sig);
    std::vector<VectorXd> final_interval;
    if (section_mat.size() > 0){
        VectorXd x_small = section_mat.col(0);
        VectorXd x_large = section_mat.col(1);
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
            final_interval = intersection(L, U, new_interval);
            return {chi2, final_interval};
        }
        else {
            VectorXd new_interval(2);
            new_interval(0) = L; 
            new_interval(1) = U;
            final_interval.push_back(new_interval);
            return {chi2, final_interval};
        }
    } 
    else {
        VectorXd new_interval(2);
        new_interval(0) = L;
        new_interval(1) = U;
        final_interval.push_back(new_interval);
        return {chi2, final_interval};
    }
}