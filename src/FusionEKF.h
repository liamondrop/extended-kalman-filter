#ifndef SRC_FUSIONEKF_H_
#define SRC_FUSIONEKF_H_

#include "Eigen/Dense"
#include "./measurement_package.h"
#include "./kalman_filter.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

class FusionEKF {
 public:
    /**
    * Constructor.
    */
    FusionEKF();

    /**
    * Destructor.
    */
    virtual ~FusionEKF();

    /**
    * Run the whole flow of the Kalman Filter from here.
    */
    void ProcessMeasurement(const MeasurementPackage &measurement_pack);

    /**
    * Kalman Filter update and prediction math lives in here.
    */
    KalmanFilter ekf_;

 private:
    // check whether the tracking toolbox was initialized or not
    bool is_initialized_;
    int64_t previous_timestamp_;

    // acceleration noise components
    float noise_ax;
    float noise_ay;

    // tool object used to compute Jacobian and RMSE
    Tools tools;
    MatrixXd R_laser_;
    MatrixXd R_radar_;
    MatrixXd H_laser_;
    MatrixXd F_;

    VectorXd getInitialX_(const MeasurementPackage &measurement_pack);
    MatrixXd getQ_(const float &dt);
};

#endif  // SRC_FUSIONEKF_H_
