#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>

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

  // initial covariance matrix
  P_ = MatrixXd(5, 5);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 30;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 30;

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
  std_a_ = 6;

  std_yawdd_ = 6;

  n_x_ = 5;

  n_aug_ = n_x_ + 2;

  P_ << 1, 0, 0, 0, 0,
          0, 1, 0, 0, 0,
          0, 0, 1, 0, 0,
          0, 0, 0, 1, 0,
          0, 0, 0, 0, 1;

  lambda_ = 3 - n_x_;

  Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ + 1);
}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage measurement_pack) {
  /**
  TODO:

  Complete this function! Make sure you switch between lidar and radar
  measurements.
  */

  if (!is_initialized_) {
    x_ = VectorXd(4);
    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      /**
      * Convert radar from polar to cartesian coordinates and initialize state.
      * Remember: you'll need to convert radar from polar to cartesian coordinates.
      */
      double rho = measurement_pack.raw_measurements_[0];
      double phi = measurement_pack.raw_measurements_[1];

      double px = rho * cos(phi);
      double py = rho * sin(phi);
      x_ << px, py, 0, 0, 0;
    } else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      /**
      Initialize state.
      */
      double px = measurement_pack.raw_measurements_[0];
      double py = measurement_pack.raw_measurements_[1];
      x_ << px, py, 0, 0, 0;
    }
    time_us_ = measurement_pack.timestamp_;
    // done initializing, no need to predict or update
    is_initialized_ = true;
    return;
  }

  double dt = (measurement_pack.timestamp_ - time_us_) / 1000000.0;    //dt - expressed in seconds
  double dt2 = dt * dt;
  double dt3 = dt2 * dt;
  double dt4 = dt3 * dt;

  Prediction(dt);

  if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {

  } else if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {

  }
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
  // Using CTRV (Constant Turn Rate and Velocity Magnitude) model

  //calculate sigma points ...
  //set sigma points as columns of matrix Xsig

  //create augmented mean state
  VectorXd x_aug = VectorXd(7);
  x_aug.head(n_x_) = x_;

  //create augmented covariance matrix
  MatrixXd P_aug = MatrixXd(7, 7);
  P_aug.topLeftCorner(n_x_, n_x_) = P_;
  P_aug.bottomRightCorner(2, 2) << std_a_ * std_a_, 0,
          0, std_yawdd_ * std_yawdd_;

  //create square root matrix
  double coef = sqrt(lambda_ + n_aug_);

  //calculate square root of P
  MatrixXd A = P_aug.llt().matrixL();

  //create augmented sigma points
  MatrixXd Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);
  Xsig_aug.col(0) = x_aug;
  MatrixXd Aplus = A * coef;
  for (int i = 0; i < n_aug_; i++) {
    Xsig_aug.col(i + 1) = x_aug + Aplus.col(i);
    Xsig_aug.col(i + 1 + n_aug_) = x_aug - Aplus.col(i);
  }

  //predict sigma points
  //avoid division by zero
  //write predicted sigma points into right column
  for (int i = 0; i < Xsig_aug.cols(); i++) {
    VectorXd xk = Xsig_aug.col(i).head(n_x_);

    double v = Xsig_aug.col(i)(2);
    double yaw = Xsig_aug.col(i)(3);
    double yaw_rate = Xsig_aug.col(i)(4);
    double nu_acc = Xsig_aug.col(i)(n_x_);
    double nu_yaw = Xsig_aug.col(i)(n_x_ + 1);

    VectorXd vec1 = VectorXd(n_x_);
    if (fabs(yaw_rate) > 0.001) {
      vec1 << v / yaw_rate * (sin(yaw + yaw_rate * delta_t) - sin(yaw)),
              v / yaw_rate * (-cos(yaw + yaw_rate * delta_t) + cos(yaw)),
              0,
              yaw_rate * delta_t,
              0;
    } else {
      vec1 << v * cos(yaw) * delta_t,
              v * sin(yaw) * delta_t,
              0,
              yaw_rate * delta_t,
              0;
    }
    VectorXd vec2 = VectorXd(n_x_);
    vec2 << (1. / 2) * (delta_t * delta_t) * cos(yaw) * nu_acc,
            (1. / 2) * (delta_t * delta_t) * sin(yaw) * nu_acc,
            delta_t * nu_acc,
            (1. / 2) * (delta_t * delta_t) * nu_yaw,
            delta_t * nu_yaw;
    Xsig_pred_.col(i) = xk + vec1 + vec2;
  }
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
}
