#include "./kalman_filter.h"

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::init(VectorXd initial_x) {
    x = initial_x;

    // initialize the P and Identity Matrices
    size_t x_size = x.size();
    P = MatrixXd::Identity(x_size, x_size);
    P *= 1000.0;  // set the initial uncertainty high
    I_ = MatrixXd::Identity(x_size, x_size);
}

void KalmanFilter::predict(const MatrixXd &F, const MatrixXd &Q) {
    x = F * x;
    P = F * P * F.transpose() + Q;
}

void KalmanFilter::update(const VectorXd &z,
                          const MatrixXd &R,
                          const MatrixXd &H) {
    MatrixXd PH_t = P * H.transpose();
    MatrixXd S = H * PH_t + R;
    MatrixXd K = PH_t * S.inverse();
    VectorXd y = z - H * x;

    // new estimate
    x += K * y;
    P = (I_ - K * H) * P;
}

void KalmanFilter::update(const VectorXd &z, const MatrixXd &R) {
    MatrixXd Hj = tools.CalculateJacobian(x);
    MatrixXd PHj_t = P * Hj.transpose();
    MatrixXd S = Hj * PHj_t + R;
    MatrixXd K = PHj_t * S.inverse();
    VectorXd y = z - tools.CartesianToPolar(x);

    // normalize the measured heading to between -pi and pi
    float phi = y(1);
    if (phi < -M_PI || phi > M_PI) {
        float TWO_PI = 2 * M_PI;
        y(1) -= TWO_PI * floor((phi + M_PI) / TWO_PI);
    }

    // new estimate
    x += K * y;
    P = (I_ - K * Hj) * P;
}
