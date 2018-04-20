#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  /**
  TODO:
    * Calculate the RMSE here.
  */
  VectorXd rmse(4);
  rmse.fill(0);

  if (estimations.size() != ground_truth.size() ||
      estimations.size() == 0) {
      std::cerr << "Invalid estimation or ground truth data" << std::endl;
      return rmse;
  }

  for (unsigned int i=0; i<estimations.size(); i++) {
      VectorXd diff = estimations[i] - ground_truth[i];
      diff = diff.array()*diff.array();
      rmse += diff;
  }

  rmse /= estimations.size();
  rmse = rmse.array().sqrt();

  return rmse;
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
  /**
  TODO:
    * Calculate a Jacobian here.
  */
  MatrixXd Hj(3, 4);
  Hj.fill(0.0);

  float px = x_state(0);
  float py = x_state(1);
  float vx = x_state(2);
  float vy = x_state(3);

  float d1 = px*px + py*py;
  float d2 = sqrt(d1);
  float d3 = d1*d2;
  float n = vx*py - vy*px;

  if (fabs(d1) < 0.0001) {
      std::cerr << "CalculateJacobian () - Error - Division by Zero" << std::endl;
      return Hj;
  }

  Hj << px/d2,      py/d2,      0,      0,
        -py/d1,     px/d1,      0,      0,
        py*n/d3,    -px*n/d3,   px/d2,  py/d2;
  return Hj;
}
