#include "kalman_filter.h"
#include <iostream>

using Eigen::MatrixXd;
using Eigen::VectorXd;

// Please note that the Eigen library does not initialize 
// VectorXd or MatrixXd objects with zeros upon creation.

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
                        MatrixXd &H_in, MatrixXd &R_in, MatrixXd &Q_in) {
  x_ = x_in;
  P_ = P_in;
  F_ = F_in;
  H_ = H_in;
  R_ = R_in;
  Q_ = Q_in;
}

void KalmanFilter::Predict() {
  /**
  TODO:
    * predict the state
  */
  x_ = F_ * x_;
  P_ = ((F_ * P_) * F_.transpose()) + Q_;
}

void KalmanFilter::Update(const VectorXd &z) {
  //Linear Kalman update step.
  //Use the linear update for the laser measurements.
  /**
  TODO:
    * update the state by using Kalman Filter equations
  */
  VectorXd y = z - (H_ * x_);
  MatrixXd Ht = H_.transpose();
  MatrixXd S = ((H_ * P_) * Ht) + R_;
  MatrixXd Si = S.inverse();
  MatrixXd PHt = P_ * Ht;
  MatrixXd K = PHt * Si;

  //new estimate
  x_ = x_ + (K * y);
  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = (I - (K * H_)) * P_;
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  //Nonlinear Kalman update step.
  //Use the nonlinear update for the radar measurments.
  /**
  TODO:
    * update the state by using Extended Kalman Filter equations
  */
  
  // Convert the cartesian state vector into the
  // polar coordinates of the radar measurements (Hx).
  float px = x_(0);
  float py = x_(1);
  float vx = x_(2);
  float vy = x_(3);
  float sqrt_pxpy = sqrt(px*px + py*py);
  VectorXd Hx = VectorXd(z.size());
  Hx << sqrt_pxpy,
        atan2(py, px),
        (px*vx + py*vy)/sqrt_pxpy;

  //calculate the prediction error y
  VectorXd y = z - Hx;
   
  //check if error is greater than Pi and if yes shift
  if(abs(y(1)) > M_PI)
      y(1) = z(1) + Hx(1);
    
  // Update filter equations.  Note these are the same equestions
  // as the Kalman Filter except the Jacobian Hj is used.
  MatrixXd Ht = H_.transpose();
  MatrixXd S = ((H_ * P_) * Ht) + R_;
  MatrixXd Si = S.inverse();
  MatrixXd PHt = P_ * Ht;
  MatrixXd K = PHt * Si;

  //new estimate
  x_ = x_ + (K * y);
  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = (I - (K * H_)) * P_;
}

/* this is a test */
