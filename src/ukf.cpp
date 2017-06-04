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
  is_initialized_ = false;

  std_a_ = 1;

  std_yawdd_ = 1;

  n_x_ = 5;

  n_aug_ = n_x_ + 2;

  P_ << 1, 0, 0, 0, 0,
          0, 1, 0, 0, 0,
          0, 0, 1, 0, 0,
          0, 0, 0, 1, 0,
          0, 0, 0, 0, 1;

  lambda_ = 3 - n_x_;

  Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ + 1);

  //set vector for weights
  weights_ = VectorXd(2 * n_aug_ + 1);
  double weight_0 = lambda_ / (lambda_ + n_aug_);
  weights_(0) = weight_0;
  for (int i = 1; i < 2 * n_aug_ + 1; i++) {
    double weight = 0.5 / (n_aug_ + lambda_);
    weights_(i) = weight;
  }

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
    x_ = VectorXd(5);
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

  Prediction(dt);

  if (use_laser_ && measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
    UpdateLidar(measurement_pack);
  }
  if (use_radar_ && measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
//    cout << "RADAR = " << endl << measurement_pack.raw_measurements_ << endl;
    UpdateRadar(measurement_pack);
  }
  time_us_ = measurement_pack.timestamp_;
}


/**
 * Ensure angular values are within [-PI, PI]
 */
void FitRadian(double& rad) {
  while (rad < -M_PI) { rad += 2 * M_PI; }
  while (rad > M_PI) { rad -= 2 * M_PI; }
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
  /**
  Estimate the object's location. Modify the state
  vector, x_. Predict sigma points, the state, and the state covariance matrix.
  */
  // This prediction uses CTRV (Constant Turn Rate and Velocity Magnitude) model

  //calculate sigma points ...
  //set sigma points as columns of matrix Xsig

  //create augmented mean state
  VectorXd x_aug = VectorXd(n_aug_);
  x_aug.fill(0.);
  x_aug.head(n_x_) = x_;

  //create augmented covariance matrix
  MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);
  P_aug.fill(0.);
  P_aug.topLeftCorner(n_x_, n_x_) = P_;
  P_aug.bottomRightCorner(2, 2) << std_a_ * std_a_, 0,
          0, std_yawdd_ * std_yawdd_;

  //create square root matrix
  double coef = sqrt(lambda_ + n_aug_);

  //calculate square root of augmented covariance matrix
  MatrixXd L = P_aug.llt().matrixL();

  //create augmented sigma points
  MatrixXd Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);

  Xsig_aug.col(0) = x_aug;
  for (int i = 0; i < n_aug_; i++) {
    Xsig_aug.col(i + 1) = x_aug + L.col(i) * coef;
    Xsig_aug.col(i + 1 + n_aug_) = x_aug - L.col(i) * coef;
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

  //predict state mean
  VectorXd x = VectorXd(n_x_);
  x.fill(0.);

  for (int i = 0; i < Xsig_pred_.cols(); i++) {
    x += weights_(i) * Xsig_pred_.col(i);
    FitRadian(x(3));
  }
  //predict state covariance matrix
  MatrixXd P = MatrixXd(n_x_, n_x_);
  P.fill(0.);

  for (int i = 0; i < Xsig_pred_.cols(); i++) {
    VectorXd offset = (Xsig_pred_.col(i) - x);
    FitRadian(offset(3));
    P += weights_(i) * offset * offset.transpose();
  }

  //update predicted state and covariance Matrix
  x_ = x;
  P_ = P;
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

  Note that Lidar measurement, containing only target position, is inherently linear.
  We can use the original Kalman filter for linear estimation.

  You'll also need to calculate the lidar NIS.
  */
//  cout << "x_ = " << endl << x_ << endl;

  VectorXd z = meas_package.raw_measurements_;

  MatrixXd I = MatrixXd::Identity(n_x_, n_x_);
  MatrixXd H = MatrixXd(2, n_x_);
  H << 1, 0, 0, 0, 0,
          0, 1, 0, 0, 0;

  MatrixXd R = MatrixXd(2, 2);
  R << std_laspx_ * std_laspx_, 0,
          0, std_laspy_ * std_laspy_;
  VectorXd y = z - H * x_;
  MatrixXd S = H * P_ * H.transpose() + R;
  MatrixXd K = P_ * H.transpose() * S.inverse();

  x_ = x_ + (K * y);
  P_ = (I - K * H) * P_;
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
  long n_z = meas_package.raw_measurements_.size();

  //create matrix for sigma points in measurement space
  MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);

  //mean predicted measurement
  VectorXd z_pred = VectorXd(n_z);
  z_pred.fill(0.);

  //measurement covariance matrix S
  MatrixXd S = MatrixXd(n_z, n_z);
  S.fill(0.);

  MatrixXd R = MatrixXd(3, 3);
  R << std_radr_ * std_radr_, 0, 0,
          0, std_radphi_ * std_radphi_, 0,
          0, 0, std_radrd_ * std_radrd_;

  for (int i = 0; i < Xsig_pred_.cols(); i++) {
    double px = Xsig_pred_.col(i)[0];
    double py = Xsig_pred_.col(i)[1];
    double v = Xsig_pred_.col(i)[2];
    double yaw = Xsig_pred_.col(i)[3];
    double yaw_rate = Xsig_pred_.col(i)[4];

    MatrixXd H_x = VectorXd(3);
    double px2 = px * px;
    double py2 = py * py;
    H_x << sqrt(px2 + py2),
            atan2(py, px),
            (px * cos(yaw) * v + py * sin(yaw) * v) / sqrt(px2 + py2);
    Zsig.col(i) = H_x;
  }

  for (int i = 0; i < Zsig.cols(); i++) {
    z_pred += weights_(i) * Zsig.col(i);
  }

  for (int i = 0; i < Zsig.cols(); i++) {
    VectorXd offset = (Zsig.col(i) - z_pred);
    FitRadian(offset(1));
    S += weights_(i) * offset * offset.transpose();
  }
  S += R;

  //create matrix for cross correlation Tc
  VectorXd z = meas_package.raw_measurements_;

  //create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x_, n_z);
  Tc.fill(0.);

  //calculate cross correlation matrix
  //calculate Kalman gain K;
  //update state mean and covariance matrix
  for (int i = 0; i < Xsig_pred_.cols(); i++) {

    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    FitRadian(x_diff(3));
    VectorXd z_diff = Zsig.col(i) - z_pred;
    FitRadian(z_diff(1));
    Tc += x_diff * z_diff.transpose() * weights_(i);
  }

  MatrixXd K = Tc * S.inverse();
  VectorXd z_diff = z - z_pred;

  FitRadian(z_diff(1));
  x_ += K * z_diff;
  FitRadian(x_(3));
  P_ -= K * S * K.transpose();
}
