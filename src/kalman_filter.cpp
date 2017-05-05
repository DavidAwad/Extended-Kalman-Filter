#include "kalman_filter.h"
#include <iostream>

#define PI  3.14159
#define TAU 6.28318

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
                        MatrixXd &H_in, MatrixXd &R_in, MatrixXd &Q_in)
{
  x_ = x_in;
  P_ = P_in; // object covariance matrix
  F_ = F_in;
  H_ = H_in;
  R_ = R_in;
  Q_ = Q_in;

}

void KalmanFilter::Predict() {
    x_ = F_ * x_;
    MatrixXd Ft = F_.transpose();
    P_ = F_ * P_ * Ft + Q_;
}

/*
 * Update
 *
 *
 * We take a new vector each update, so we perform a similar process to UpdateEKF
 * update our state matrices
 *
 * y = z - H_{x}'
 * S = HP'H^{T} + R
 * K = P'H^{T}S^{-1}
 * x = x' + Ky
 * P = (I-KH)P'
 */
void KalmanFilter::Update(const VectorXd &z) {

    VectorXd y = z - (H_ * x_);

    MatrixXd Ht = H_.transpose();

    MatrixXd S = H_ * P_ * Ht + R_;

    MatrixXd K = (P_ * Ht) * S.inverse();

    //compute our new estimate
    x_ = x_ + (K * y);
    long x_size = x_.size();
    MatrixXd I = MatrixXd::Identity(x_size, x_size);
    P_ = (I - K * H_) * P_;
}



/*
 * UpdateEKF
 *
 * We're going to make an update to our measurements
 * This involves the following equations: (see 13-17 in the included KF reference)
 *
 *
 * y = z - H_{x}'
 * S = HP'H^{T} + R
 * K = P'H^{T}S^{-1}
 * x = x' + Ky
 * P = (I-KH)P'
 */
void KalmanFilter::UpdateEKF(const VectorXd &z) {

    VectorXd y = z - hx_; // (13)

    // make sure phi is in range
    while (true) {
      if (y(1) > PI) {
        y(1) = y(1) - TAU;
      }
      else if (y(1) < - PI) {
        y(1) = y(1) + TAU;
      }
      else { break; }
    }

    // avoid dupicate work by setting h
    MatrixXd H_t = H_.transpose();

    MatrixXd S = H_ * P_ * H_t + R_; // (14)

    MatrixXd K = P_ * H_t * S.inverse(); // (15)

    //compute our new estimate!
    x_ = x_ + (K * y); // (16)
    long x_size = x_.size();

    MatrixXd I = MatrixXd::Identity(x_size, x_size);
    P_ = (I - K * H_) * P_; // (17)
}
