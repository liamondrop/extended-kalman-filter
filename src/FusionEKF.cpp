#include <iostream>
#include "./FusionEKF.h"

/**
 * Constructor.
 */
FusionEKF::FusionEKF() {
    is_initialized_ = false;

    previous_timestamp_ = 0;

    // initializing matrices
    R_laser_ = MatrixXd(2, 2);
    R_radar_ = MatrixXd(3, 3);

    // measurement covariance matrix - laser
    R_laser_ << 0.0225, 0,
                0, 0.0225;

    // measurement covariance matrix - radar
    R_radar_ << 0.09,  0,  0,
                0, 0.0009, 0,
                0,  0,  0.09;

    // measurement matrix
    H_laser_ = MatrixXd(2, 4);
    H_laser_ << 1, 0, 0, 0,
                0, 1, 0, 0;

    // the initial transition matrix F_
    F_ = MatrixXd(4, 4);
    F_ << 1, 0, 1, 0,
          0, 1, 0, 1,
          0, 0, 1, 0,
          0, 0, 0, 1;

    noise_ax = 9;
    noise_ay = 9;
}

/**
 * Destructor.
 */
FusionEKF::~FusionEKF() {}

VectorXd FusionEKF::getInitialX_(const MeasurementPackage &measurement_pack) {
    VectorXd initial_x(4);
    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
        float rho = measurement_pack.raw_measurements_[0];
        float phi = measurement_pack.raw_measurements_[1];
        float rho_dot = measurement_pack.raw_measurements_[2];
        initial_x << rho * sin(phi),
                     rho * cos(phi),
                     rho_dot * cos(phi),
                     rho_dot * sin(phi);
    } else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
        initial_x << measurement_pack.raw_measurements_[0],
                     measurement_pack.raw_measurements_[1],
                     0.0,
                     0.0;
    }
    return initial_x;
}

MatrixXd FusionEKF::getQ_(const float &dt) {
    float dt2 = (dt*dt < 1e-5) ? 0 : dt * dt;  // dt^2
    float dt3 = (dt2 * dt) / 2.0;              // dt^3/2
    float dt4 = (dt2 * dt2) / 4.0;             // dt^4/4
    MatrixXd Q = MatrixXd(4, 4);
    Q << dt4*noise_ax, 0.0, dt3*noise_ax, 0.0,
         0.0, dt4*noise_ay, 0.0, dt3*noise_ay,
         dt3*noise_ax, 0.0, dt2*noise_ax, 0.0,
         0.0, dt3*noise_ay, 0.0, dt2*noise_ay;
    return Q;
}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {
    // compute the time elapsed between measurements
    float dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1e+6;
    previous_timestamp_ = measurement_pack.timestamp_;

    // check if is initialized and delta_T is a sane value
    if (!is_initialized_ || dt <= 0 || dt > 100) {
        std::cout << "INIT" << std::endl;
        std::cout << "++++" << std::endl;
        VectorXd initial_x = getInitialX_(measurement_pack);
        ekf_.init(initial_x);

        // done initializing, no need to predict or update
        is_initialized_ = true;
    } else {
        // Modify the F matrix so that the time is integrated
        F_(0, 2) = dt;
        F_(1, 3) = dt;

        // Calculate the process covariance matrix Q
        MatrixXd Q = getQ_(dt);

        // Prediction step
        ekf_.predict(F_, Q);

        // Update step
        if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
            ekf_.update(measurement_pack.raw_measurements_, R_radar_);
        } else {
            ekf_.update(measurement_pack.raw_measurements_, R_laser_, H_laser_);
        }
    }

    std::cout << "x =\n" << ekf_.x << std::endl;
    std::cout << "P =\n" << ekf_.P << std::endl;
    std::cout << "=========================" << std::endl;
}
