#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>
#include "tools.h"

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() {
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd(5);
  x_ << 0, 0, 0, 0, 0;

  // initial covariance matrix
  P_ = MatrixXd(5, 5);
  P_ = MatrixXd::Identity(5, 5);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 1.6;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 0.9;

  // Laser measurement noise standard deviation position1 in m
  std_laspx_ = 0.15;

  // Laser measurement noise standard deviation position2 in m
  std_laspy_ = 0.15;

  // Radar measurement noise standard deviation radius in m
  std_radr_ = 0.3;

  // Radar measurement noise standard deviation angle in rad
  std_radphi_ = 0.03;

  // Radar measurement noise standard deviation radius change in m/s
  std_radrd_ = 0.3;

  /**
  TODO:

  Complete the initialization. See ukf.h for other member properties.

  Hint: one or more values initialized above might be wildly off...
  */

  time_us_ = 0;

  // State dimension
  n_x_ = 5;

  // Augmented state dimension
  n_aug_ = 7;

  // Define spreading parameter
  lambda_ = 3 - n_x_;

  weights_ = VectorXd(2*n_aug_+1);
  //set weights
  for(int i=0; i < (2 * n_aug_ + 1); i++) {
    if (i == 0) {
      weights_(i) = lambda_ / (lambda_ + n_aug_);
    } else {
      weights_(i) = 1 / (2 * (lambda_ + n_aug_));
    }
  }

  // Initialize X sigma predictions.
  Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ + 1);

}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Make sure you switch between lidar and radar
  measurements.
  */

  cout << "MEASUREMENT: ";
  if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
    cout << "RADAR\n\n";
  } else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
    cout << "LASER\n\n";
  }
  if (!is_initialized_) {
    /**
    TODO:
      * Initialize the state ekf_.x_ with the first measurement.
      * Create the covariance matrix.
      * Remember: you'll need to convert radar from polar to cartesian coordinates.
    */

    if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
      // Convert radar from polar to cartesian coordinates and initialize state.
      double rho = meas_package.raw_measurements_(0);
      double phi = meas_package.raw_measurements_(1);
      // Can we use ro-hat to set the predicted velocity? Per "tips" we don't get enough data for that,
      // so set x and y, but leave vx and vy at zero.
      x_(0) = rho * cos(phi);
      x_(1) = rho * sin(phi);
      std::cout << "RADAR_INIT\n\n";
    }
    else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
      // Initialize state but laser only gives position data, so leave vx and vy at zero.
      x_(0) = meas_package.raw_measurements_(0);
      x_(1) = meas_package.raw_measurements_(1);
      std::cout << "LIDAR_INIT\n\n";
    }
    else {

    }

    time_us_ = meas_package.timestamp_;

    // done initializing, no need to predict or update
    is_initialized_ = true;
    return;
  }
  // Compute the time elapsed between the current and previous measurements.
  double delta_t = (meas_package.timestamp_ - time_us_) / 1000000.0;	//dt - expressed in seconds

  if (meas_package.sensor_type_ == MeasurementPackage::RADAR && use_radar_) {
    std::cout << "PREDICTION\n\n";
    UKF::Prediction(delta_t);
    std::cout << "RADAR_UPDATE\n\n";
    UKF::UpdateRadar(meas_package);

  }
  else if (meas_package.sensor_type_ == MeasurementPackage::LASER && use_laser_) {
    std::cout << "PREDICTION\n\n";
    UKF::Prediction(delta_t);
    std::cout << "LASER_UPDATE\n\n";
    UKF::UpdateLidar(meas_package);

  }
  else if (meas_package.sensor_type_ != MeasurementPackage::RADAR && meas_package.sensor_type_ != MeasurementPackage::LASER ) {
    cout << "Unknown sensor type: " << meas_package.sensor_type_;
    throw std::exception();
  }

  time_us_ = meas_package.timestamp_;


}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
  /**
  TODO:

  Complete this function! Estimate the object's location. Modify the state
  vector, x_. Predict sigma points, the state, and the state covariance matrix.
  */

  //create augmented covariance matrix
  //create augmented state covariance
  MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);
  P_aug.fill(0.0);
  P_aug.topLeftCorner(5,5) = P_;
  P_aug(5,5) = std_a_*std_a_;
  P_aug(6,6) = std_yawdd_*std_yawdd_;
  cout << "P_aug: " << P_aug << "\n\n";


  //create square root matrix
  MatrixXd L = P_aug.llt().matrixL();
  cout << "L: " << L << "\n\n";
  //calculate sigma points ...

  //create augmented mean vector
  VectorXd x_aug = VectorXd(n_aug_);
  x_aug.fill(0);
  x_aug.head(5) << x_;
  cout << "x_aug: " << x_aug << "\n\n";

  //create sigma point matrix
  MatrixXd Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);
  Xsig_aug.fill(0.0);
  Xsig_aug.col(0)  = x_aug;
  for (int i = 0; i< n_aug_; i++)
  {
    Xsig_aug.col(i+1)        = x_aug + sqrt(lambda_+n_aug_) * L.col(i);
    Xsig_aug.col(i+1+n_aug_) = x_aug - sqrt(lambda_+n_aug_) * L.col(i);
  }

  for (int i = 0; i< 2*n_aug_+1; i++) {
    //predict sigma points
    double p_x = Xsig_aug(0, i);
    double p_y = Xsig_aug(1, i);
    double v = Xsig_aug(2, i);
    double yaw_angle = Xsig_aug(3, i);
    double yaw_rate = Xsig_aug(4, i);
    double accel_noise = Xsig_aug(5, i);
    double yaw_rate_noise = Xsig_aug(6, i);


    VectorXd pt = VectorXd(n_x_);
    VectorXd a = VectorXd(n_x_);
    VectorXd b = VectorXd(n_x_);
    pt << p_x, p_y, v, yaw_angle, yaw_rate;
    cout << "pt: " << i << "\n" << pt << "\n\n";

    //predicted state values
    double px_p, py_p;

    //avoid division by zero
    if (fabs(yaw_rate) > 0.001) {
      px_p = p_x + (v / yaw_rate) * (sin(yaw_angle+(yaw_rate * delta_t)) - sin(yaw_angle));
      py_p = p_y + (v / yaw_rate) * (-cos(yaw_angle+(yaw_rate * delta_t)) + cos(yaw_angle));
    }
    else {
      px_p = p_x + v * cos(yaw_angle) * delta_t;
      py_p = p_y + v * sin(yaw_angle) * delta_t;
    }

    double v_p = v;
    double yaw_angle_p = yaw_angle + yaw_rate * delta_t;
    double yaw_rate_p = yaw_rate;

    //add noise
    px_p = px_p + 0.5 * accel_noise * delta_t * delta_t * cos(yaw_angle);
    py_p = py_p + 0.5 * accel_noise * delta_t * delta_t * sin(yaw_angle);
    v_p = v_p + accel_noise*delta_t;

    yaw_angle_p = yaw_angle_p + 0.5*yaw_rate_noise*delta_t*delta_t;
    yaw_rate_p = yaw_rate_p + yaw_rate_noise*delta_t;

    //write predicted sigma point into right column
    Xsig_pred_(0, i) = px_p;
    Xsig_pred_(1, i) = py_p;
    Xsig_pred_(2, i) = v_p;
    Xsig_pred_(3, i) = yaw_angle_p;
    Xsig_pred_(4, i) = yaw_rate_p;
  }

  // Predicted state mean
  // The summation of predicted columns times the weights (which reverses the spread).
  //create vector for predicted state
  //VectorXd x_mean = VectorXd(n_x_);
  cout << "before x_: " << x_ << "\n\n";
  cout << "weights_: " << weights_ << "\n\n";
  cout << "Xsig_pred_: " << Xsig_pred_ << "\n\n";

  x_.fill(0.0);
  for(int i=0; i < (2 * n_aug_ + 1); i++) {
    x_ = x_ + weights_(i) * Xsig_pred_.col(i);
  }
  cout << "after x_: " << x_ << "\n\n";

  //predicted state covariance matrix
  //create covariance matrix for prediction
  P_.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {

    // The difference between the predicted sigma point and the mean (average) point.
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    cout << "x_diff: " << x_diff << "\n\n";

    Tools::NormalizeRadians(x_diff(3));

    // Covariance Matrix recovered from the predicted sigma points.
    P_ = P_ + weights_(i) * x_diff * x_diff.transpose() ;

  }
  cout << "P_: " << P_ << "\n\n";
  cout << "\n";


}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use lidar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the lidar NIS.
  */

  //set measurement dimension, laser can measure py and px;
  int n_z = 2;

  //create matrix for sigma points in measurement space
  MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);
  //Zsig.fill(0.0);

  //mean predicted measurement
  VectorXd z_pred = VectorXd(n_z);
  z_pred.fill(0.0);

  //measurement covariance matrix S
  MatrixXd S = MatrixXd(n_z,n_z);
  S.fill(0.0);

  cout << "LIDAR_A\n\n";
  //transform sigma points into measurement space
  for (int i=0; i<2*n_aug_+1; i++) {
    Zsig(0, i) = Xsig_pred_(0, i); //px
    Zsig(1, i) = Xsig_pred_(1, i); //py
  }
  cout << "Zsig: " << Zsig << "\n\n";
  cout << "LIDAR_B\n\n";
  //calculate mean predicted measurement
  for(int i=0; i < (2 * n_aug_ + 1); i++) {
    z_pred = z_pred + weights_(i) * Zsig.col(i);
  }
  cout << "z_pred: " << z_pred << "\n\n";
  cout << "LIDAR_C\n\n";
  //calculate measurement covariance matrix S
  //predicted state covariance matrix
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {

    // The difference between the predicted sigma point and the mean (average) point.
    VectorXd z_diff = Zsig.col(i) - z_pred;

    // Covariance Matrix recovered from the predicted sigma points.
    S = S + weights_(i) * z_diff * z_diff.transpose();
  }
  cout << "S: " << S << "\n\n";

  //add measurement noise covariance matrix
  MatrixXd R = MatrixXd(n_z,n_z);
  R <<    std_laspx_*std_laspx_, 0,
      0, std_laspy_*std_laspy_;

  S = S + R;
  cout << "S+R: " << S << "\n\n";


  //create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x_, n_z);
  Tc.fill(0.0);

  cout << "LIDAR_D\n\n";
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {

    VectorXd z_diff = Zsig.col(i) - z_pred;
    VectorXd x_diff = Xsig_pred_.col(i) - x_;

    // yaw angle normalization for angles that may be outside of +/- PI.
    Tools::NormalizeRadians(x_diff(3));
    Tools::NormalizeRadians(z_diff(1));

    // Covariance Matrix recovered from the predicted sigma points.
    Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
  }
  cout << "Tc: " << Tc << "\n\n";

  //calculate Kalman gain K;
  MatrixXd K = MatrixXd(n_x_, n_z);
  K = Tc * S.inverse();
  cout << "K: " << K << "\n\n";

  //update state mean and covariance matrix
  x_ = x_ + K * (meas_package.raw_measurements_ - z_pred);
  P_ = P_ - K * S * K.transpose();

  cout << "x_: " << x_ << "\n\n";
  cout << "P_: " << P_ << "\n\n";
  cout << "LIDAR_DONE\n\n";

}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use radar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the radar NIS.
  */

  //set measurement dimension, radar can measure r, phi, and r_dot
  int n_z = 3;

  //create matrix for sigma points in measurement space
  MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);
  Zsig.fill(0);
  cout << "Zsig: " << Zsig << "\n\n";

  //mean predicted measurement
  VectorXd z_pred = VectorXd(n_z);
  z_pred.fill(0);
  cout << "z_pred: " << z_pred << "\n\n";

  //measurement covariance matrix S
  MatrixXd S = MatrixXd(n_z,n_z);
  S.fill(0);
  cout << "S: " << S << "\n\n";

  //transform sigma points into measurement space
  for (int i=0; i<2*n_aug_+1; i++) {
    double px = Xsig_pred_(0, i);
    double py = Xsig_pred_(1, i);
    double v = Xsig_pred_(2, i);
    double yaw = Xsig_pred_(3, i);

    Zsig(0, i) = sqrt(px*px+py*py);
    Zsig(1, i) = atan2(py,px);
    if (Zsig(0, i) == 0) {
      continue;
    }
    Zsig(2, i) = (px * cos(yaw) * v + py * sin(yaw) * v) / Zsig(0, i);

  }
  cout << "ZSig: " << Zsig << "\n\n";
  cout << "RADAR_A\n\n";
  //calculate mean predicted measurement
  for(int i=0; i < (2 * n_aug_ + 1); i++) {
    z_pred = z_pred + weights_(i) * Zsig.col(i);
  }
  //z_pred  = Zsig * weights_;
  cout << "z_pred: " << z_pred << "\n\n";
  cout << "RADAR_B\n\n";
  //calculate measurement covariance matrix S
  //predicted state covariance matrix
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {

    // The difference between the predicted sigma point and the mean (average) point.
    VectorXd z_diff = Zsig.col(i) - z_pred;

    // yaw angle normalization for angles that may be outside of +/- PI.
    Tools::NormalizeRadians(z_diff(1));

    // Covariance Matrix recovered from the predicted sigma points.
    S = S + weights_(i) * z_diff * z_diff.transpose();
  }
  cout << "S: " << S << "\n\n";
  cout << "RADAR_C\n\n";
  //add measurement noise covariance matrix
  MatrixXd R = MatrixXd(n_z,n_z);
  R <<    std_radr_*std_radr_, 0, 0,
      0, std_radphi_*std_radphi_, 0,
      0, 0,std_radrd_*std_radrd_;

  S = S + R;
  cout << "S+R: " << S << "\n\n";

  //create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x_, n_z);
  Tc.fill(0.0);

  for (int i = 0; i < 2 * n_aug_ + 1; i++) {

    VectorXd z_diff = Zsig.col(i) - z_pred;
    VectorXd x_diff = Xsig_pred_.col(i) - x_;

    // yaw angle normalization for angles that may be outside of +/- PI.
    Tools::NormalizeRadians(z_diff(1));
    Tools::NormalizeRadians(x_diff(3));

    // Covariance Matrix recovered from the predicted sigma points.
    Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
  }

  VectorXd z_diff_full = meas_package.raw_measurements_ - z_pred;
  Tools::NormalizeRadians(z_diff_full(1));

  cout << "RADAR_D\n\n";
  //calculate Kalman gain K;
  MatrixXd K = MatrixXd(n_x_, n_z);
  K = Tc * S.inverse();

  //update state mean and covariance matrix
  x_ = x_ + K * z_diff_full;
  P_ = P_ - K * S * K.transpose();

  cout << "RADAR_DONE\n\n";
  cout << x_ << "\n\n";
}
