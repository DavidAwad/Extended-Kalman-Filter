#include <iostream>
#include "FusionEKF.h"
#include "tools.h"
#include "Eigen/Dense"

#define MIN 0.0001 // some min threshold for fabs()

using namespace std;
using std::vector;
using Eigen::MatrixXd;
using Eigen::VectorXd;

/*
 * Constructor.
 */
FusionEKF::FusionEKF() {
  is_initialized_ = false;
  previous_timestamp_ = 0;

  // initializing matrices
  H_laser_ = MatrixXd(2, 4);
  Hj_ = MatrixXd(3, 4);
  R_laser_ = MatrixXd(2, 2);
  R_laser_ << 0.0225, 0,
        0, 0.0225;

  //measurement covariance matrix - radar
  R_radar_ = MatrixXd(3, 3);
  R_radar_ << 0.09, 0,      0,
              0,    0.0009, 0,
              0,    0,      0.09;

  // Initialise P (object covariance) as zero matrix
  ekf_.P_ = MatrixXd(4, 4);
  ekf_.P_ <<  1, 0, 0,    0,
              0, 1, 0,    0,
              0, 0, 1000, 0,
              0, 0, 0,    1000;
  // Initialise F with dt = 0
  ekf_.F_ = MatrixXd(4, 4);
  ekf_.F_ << 1, 0, 0, 0,
             0, 1, 0, 0,
             0, 0, 1, 0,
             0, 0, 0, 1;

  // Initialise H_laser
  H_laser_ << 1, 0, 0, 0,
              0, 1, 0, 0;

  // Initialise H
  ekf_.H_ = MatrixXd(4, 4);
  ekf_.H_ << 1, 0, 0, 0,
             0, 1, 0, 0,
             0, 0, 1, 0,
             0, 0, 0, 1;
}

/**
* Destructor.
*/
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {

  /*****************************************************************************
   *  Initialization
   ****************************************************************************/
  if (!is_initialized_) {
    // first measurement
    ekf_.x_ = VectorXd(4);
    ekf_.x_ << 1, 1, 1, 1;

    // Create process covariance matrix
    ekf_.Q_ = Eigen::MatrixXd(4,4);

    float px;
    float py;

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      /**
      Convert radar from polar to cartesian coordinates and initialize state.
      */
      float rho = measurement_pack.raw_measurements_[0];
      float phi = measurement_pack.raw_measurements_[1];

      px = rho*cos(phi);
      py = rho*sin(phi);

    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      px = measurement_pack.raw_measurements_[0];
      py = measurement_pack.raw_measurements_[1];
    }

    // Handle small px, py
    if(fabs(px) < MIN){
        px = 0.1;
    }

    if(fabs(py) < MIN){
        py = 0.1;
    }

    ekf_.x_ << px, py, 0, 0;

    previous_timestamp_ = measurement_pack.timestamp_;
    // done initializing, no need to predict or update
    is_initialized_ = true;
    return;
  }

  /*****************************************************************************
   *  Prediction
   ****************************************************************************/

    float dt = (measurement_pack.timestamp_ - previous_timestamp_);

    dt = dt / 1000000.0;
    previous_timestamp_ = measurement_pack.timestamp_;

    // avoid repeating any work
    float dt_2 = dt * dt;
    float dt_3 = dt_2 * dt;
    float dt_4 = dt_3 * dt;

    // additional constants to avoid extra work
    float dt_3_2 = dt_3/2;
    float dt_4_4 = dt_4/4;

    // integrate with respect to time
    ekf_.F_(0, 2) = dt;
    ekf_.F_(1, 3) = dt;

    // using given noise values
    float noise_ax = 9;
    float noise_ay = 9;

    // update process covariance
    ekf_.Q_ <<  dt_4_4 * noise_ax, 0, dt_3_2 * noise_ax, 0,
                0, dt_4_4 * noise_ay, 0, dt_3_2 * noise_ay,
                dt_3_2 * noise_ax, 0, dt_2 * noise_ax, 0,
                0, dt_3_2 * noise_ay, 0, dt_2 * noise_ay;
  ekf_.Predict();

  /*****************************************************************************
   *  Update
   ****************************************************************************/

  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    // Radar updates
    ekf_.hx_ = VectorXd(3);

    float px = ekf_.x_[0];
    float py = ekf_.x_[1];
    float vx = ekf_.x_[2];
    float vy = ekf_.x_[3];

    float rho;
    float phi;
    float rho_dot;

    if(fabs(px) < MIN or fabs(py) < MIN){
      if(fabs(px) < MIN){
        px = MIN;
      }
      if(fabs(py) < MIN){
        py = MIN;
      }

      rho = sqrt(px*px + py*py);
      phi = 0;
      rho_dot = 0;

    } else {
      rho = sqrt(px*px + py*py);
      phi = atan2(py,px);
      rho_dot = (px*vx + py*vy) /rho;
    }

    ekf_.hx_ << rho, phi, rho_dot;

    // we're using a radar meaasurement so we need to calculate jacobian
    Hj_ = tools.CalculateJacobian(ekf_.x_);

    // don't update measurement if we can't compute the Jacobian
    if (Hj_.isZero(0)){ return; }

    ekf_.H_ = Hj_;
    ekf_.R_ = R_radar_;
    ekf_.UpdateEKF(measurement_pack.raw_measurements_);
  } else {
    // Laser updates
    ekf_.H_ = H_laser_;
    ekf_.R_ = R_laser_;
    ekf_.Update(measurement_pack.raw_measurements_);
  }
}
