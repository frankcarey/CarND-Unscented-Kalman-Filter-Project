#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  VectorXd rmse(4);
  rmse << 0,0,0,0;

  // check the validity of the following inputs:
  //  * the estimation vector size should not be zero
  //  * the estimation vector size should equal ground truth vector size
  if (estimations.empty()) {
    cout << "Err: No estimations found!";
    return rmse;
  }
  if (estimations.size() != ground_truth.size()) {
    cout << "Err: Number of estimations and ground_truth are not the same!";
    return rmse;
  }

  //accumulate squared residuals
  for(int i=0; i < estimations.size(); ++i){
    // ... your code here
    VectorXd diff;
    diff = estimations[i] - ground_truth[i];
    rmse << rmse.array() + diff.array() * diff.array();
  }

  // calculate the mean
  rmse /= estimations.size();

  //calculate the squared root
  rmse = rmse.array().sqrt();

  //return the result
  return rmse;
}

void Tools::NormalizeRadians(double &r) {
  // yaw angle normalization for angles that may be outside of +/- PI.
  double dif = fmod(r + M_PI, 2*M_PI);
  if (dif < 0) {
    dif += 2 * M_PI;
  }
  r = dif - M_PI;
}