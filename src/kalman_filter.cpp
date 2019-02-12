#include "kalman_filter.h"
#include <iostream>
#include <math.h>

using Eigen::MatrixXd;
using Eigen::VectorXd;
using namespace std;

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

  // Kalman Filter Equations - Predictions
  // Lesson 5: Part 8 - Kalman Filter Equations in C++ 
  // Same for both Laser and Radar Measurements
  
  // Kalman Prediction Step
  x_ = F_ * x_;
  MatrixXd Ft_ = F_.transpose();
  P_ = F_ * P_ * Ft_ + Q_;
  
}

void KalmanFilter::Update(const VectorXd &z) {

  // Kalman Filter Equations Update
  // Lesson 5: Part 8 - Kalman Filter Equations in C++
  // Used only for Laser data (already linear)
  
  // Kalman Update Step
  VectorXd y_ = z - H_ * x_;
  MatrixXd Ht_ = H_.transpose();
  MatrixXd S_ = H_ * P_ * Ht_ + R_;
  MatrixXd Si_ = S_.inverse();
  MatrixXd K_ =  P_ * Ht_ * Si_;

  // Define Identity Matrix
  int x_size = x_.size();
  MatrixXd I_ = MatrixXd::Identity(x_size, x_size);
  
  // Update State and Covariance Matrix
  x_ = x_ + (K_ * y_);
  P_ = (I_ - K_ * H_) * P_;
  
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  
  // Extend Kalman Filter Equations Update
  // Lesson 5: Part 21 - EKF Algorithm Generalisation 
  // Used only for Radar measurements (non-linear)
  
  // Extract state values
  float x = x_(0);
  float y = x_(1);    
  float vx = x_(2);
  float vy = x_(3);
  
  // Convert to polar coordinates
  float rho = sqrt(x * x + y * y);
  float phi = atan2(y, x);
  float rho_dot;
  
  // Check division by zero
  if(fabs(rho) < 0.0001){ 
      rho_dot = 0;
  } else {
      rho_dot = (vx * x + vy * y) / rho;
  }
  
  VectorXd z_pred = VectorXd(3);
  z_pred << rho, phi, rho_dot;
  
  // Kalman Update Step
  VectorXd y_ = z - z_pred;
  
  // Normalise y-Phi between -PI, PI; 
  // Code referenced from #s-t2-p-extended-kalman slack channel
  if (y_(1) < -M_PI){
     while(y_(1) < -M_PI){
       y_(1) = y_(1) + 2 * M_PI;
     }
   }
   else if (y_(1) > M_PI){
     while(y_(1) > M_PI){
       y_(1) = y_(1) - 2 * M_PI;
     }
   }
  
  // Extended Kalman Update
  MatrixXd Ht_ = H_.transpose();
  MatrixXd S_ = H_ * P_ * Ht_ + R_;
  MatrixXd Si_ = S_.inverse();
  MatrixXd K_ =  P_ * Ht_ * Si_;

  // Define Identity Matrix
  int x_size = x_.size();
  MatrixXd I_ = MatrixXd::Identity(x_size, x_size);
  
  
  // Update State and Covariance Matrix
  x_ = x_ + (K_ * y_);
  P_ = (I_ - K_ * H_) * P_;  
  
}
