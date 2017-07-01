#ifndef SRC_KALMAN_FILTER_H_
#define SRC_KALMAN_FILTER_H_
#include "Eigen/Dense"
#include "./measurement_package.h"
#include "./tools.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;

class KalmanFilter {
 public:
    VectorXd x;  // state vector
    MatrixXd P;  // state covariance matrix

    KalmanFilter();

    virtual ~KalmanFilter();

    /**
     * Inititalize the KF
     */
    void init(VectorXd initial_x);

    /**
    * Predicts the state and the state covariance
    * using the process model
    * @param F The state transition matrix
    * @param Q The process covariance matrix
    */
    void predict(const MatrixXd &F, const MatrixXd &Q);

    /**
    * Updates the state by using standard linear Kalman Filter equations
    * @param z The measurement at k+1
    * @param R The measurement covariance matrix
    * @param H The matrix mapping the state to the measurements
    */
    void update(const VectorXd &z, const MatrixXd &R, const MatrixXd &H);

    /**
    * Updates the state by using Extended Kalman Filter equations
    * @param z The measurement at k+1
    * @param R The measurement covariance matrix
    */
    void update(const VectorXd &z, const MatrixXd &R);

 private:
    int sensor_type_;
    Tools tools;

    // identity matrix
    MatrixXd I_;
};

#endif  // SRC_KALMAN_FILTER_H_
