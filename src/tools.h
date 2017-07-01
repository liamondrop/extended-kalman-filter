#ifndef SRC_TOOLS_H_
#define SRC_TOOLS_H_
#include <vector>
#include "Eigen/Dense"

using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

class Tools {
 public:
    /**
     * Constructor.
     */
    Tools();

    /**
     * Destructor.
     */
    virtual ~Tools();

    /**
     * A helper method to calculate RMSE.
     */
    VectorXd CalculateRMSE(const vector<VectorXd> &estimations,
                           const vector<VectorXd> &ground_truth);

    /**
     * Map cartesian coordinates px, py, vx, vy
     * to polar coords rho, phi, rho_dot
     */
    VectorXd CartesianToPolar(const VectorXd &x);

    /**
     * A helper method to calculate Jacobian Matrix for Radar Measurement
     */
    MatrixXd CalculateJacobian(const VectorXd &x);
};

#endif  // SRC_TOOLS_H_
