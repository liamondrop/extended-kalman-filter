#include "./tools.h"
#include <iostream>
#include <vector>

using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
    VectorXd RMSE(4);
    RMSE << 0, 0, 0, 0;

    const std::size_t N = estimations.size();
    if (N == 0 || N != ground_truth.size()) {
        return RMSE;
    }

    // accumulate squared residuals
    for (std::size_t t=0; t < N; ++t) {
        VectorXd res = estimations[t] - ground_truth[t];
        res = res.array().square();
        RMSE += res;
    }

    RMSE = (RMSE / N).array().sqrt();
    return RMSE;
}

VectorXd Tools::CartesianToPolar(const VectorXd &x) {
    float px = x(0);
    float py = x(1);
    float vx = x(2);
    float vy = x(3);

    float rho = sqrt(px*px + py*py);
    float phi = (px == 0) ? 0.0 : atan2(py, px);
    float rho_dot = (rho == 0) ? 0.0 : (px*vx + py*vy) / rho;

    VectorXd out(3);
    out << rho, phi, rho_dot;
    return out;
}

MatrixXd Tools::CalculateJacobian(const VectorXd &x) {
    MatrixXd Hj(3, 4);

    // recover state parameters
    float px = x(0);
    float py = x(1);
    float vx = x(2);
    float vy = x(3);
    float norm_sq = px*px + py*py;

    // check division by vanishingly small numbers
    if (std::abs(norm_sq) < 1e-5) {
        std::cout << "CalculateJacobian() - Divide by Zero Error" << std::endl;
        return Hj;
    }

    float norm = sqrt(norm_sq);
    float norm_sq15 = norm * norm_sq;
    float vxpy = vx * py;
    float vypx = vy * px;

    // compute the Jacobian matrix
    float h00 = px/norm;
    float h01 = py/norm;
    float h10 = -py/norm_sq;
    float h11 = px/norm_sq;
    float h20 = py*(vxpy-vypx)/norm_sq15;
    float h21 = px*(vypx-vxpy)/norm_sq15;

    Hj << h00, h01, 0.0, 0.0,
          h10, h11, 0.0, 0.0,
          h20, h21, h00, h01;

    return Hj;
}
