#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 * This is scaffolding, do not modify
 */

static void adjust_angle(double &angle){
  while(angle>M_PI) angle-=2*M_PI;
  while(angle<-M_PI) angle += 2*M_PI;
}

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
  std_a_ = 1;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 0.29;
  
  //DO NOT MODIFY measurement noise values below these are provided by the sensor manufacturer.
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
  //DO NOT MODIFY measurement noise values above these are provided by the sensor manufacturer.
  
  /**
  TODO:

  Complete the initialization. See ukf.h for other member properties.

  Hint: one or more values initialized above might be wildly off...
  */
  is_initialized_ = false;
  n_x_ = 5;
  n_aug_ = 7;
  lambda_ = 3 - n_aug_;
  Xsig_pred_ = MatrixXd(n_x_, 2*n_aug_+1);
  delta_Xsig_pred_ = MatrixXd(n_x_, 2*n_aug_+1);
  W_ = MatrixXd(2*n_aug_ + 1, 2*n_aug_+1);
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
  if (is_initialized_){
    double delta_t = meas_package.timestamp_ - time_us_;
    
    Prediction(delta_t);
    
    
    if(meas_package.sensor_type_ == MeasurementPackage::LASER){
      
      UpdateLidar(meas_package);
      
      

    }else if (meas_package.sensor_type_ == MeasurementPackage::RADAR){
      
      UpdateRadar(meas_package);
      

    }else{
      cout<<"The measurement type is not defined"<<endl;
    }
    time_us_ = meas_package.timestamp_;
  }else{
    //Initializing the state vector.
    x_.fill(0.0);
    P_.fill(0.0);
    P_.diagonal()<<1,1,1,1,1;
    if(meas_package.sensor_type_ == MeasurementPackage::LASER){
      x_.head(2) = meas_package.raw_measurements_;
      //P_(0,0) = pow(std_laspx_,2);
      //P_(1,1) = pow(std_laspy_,2);
    }else if (meas_package.sensor_type_ == MeasurementPackage::RADAR){
      double rho = meas_package.raw_measurements_(0);
      double phi = meas_package.raw_measurements_(1);
      //P_(0,0) = pow(std_radr_,2);
      //P_(1,1) = pow(std_radr_,2);
      x_.head(2)<<rho*cos(phi),rho*sin(phi);
    }else{
      cout<<"The measurement type is not defined"<<endl;
    }
    time_us_ = meas_package.timestamp_;
    is_initialized_ = true;
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
  //This finds the sigma points
  delta_t/=1000000;
  MatrixXd x_sig (n_aug_,2*n_aug_+1);
  MatrixXd P_aug (n_aug_, n_aug_);
  P_aug.fill(0.0);
  P_aug.topLeftCorner(n_x_,n_x_)=P_;
  P_aug.bottomRightCorner(2,2)<<pow(std_a_,2),0,
                                0,pow(std_yawdd_,2);
  MatrixXd root_P = P_aug.llt().matrixL();
  cout<<"P_aug"<<endl;
  cout<<P_aug<<endl;
  cout<<"root_P"<<endl;
  cout<<root_P<<endl;
  cout<<"root_P squared"<<endl;
  cout<<root_P.transpose()*root_P<<endl;
  VectorXd x_aug (n_aug_);
  x_aug.fill(0.0);
  x_aug.head(n_x_) = x_;
  x_sig.col(0)=x_aug;
  for (int i = 0; i<n_aug_; i++){
    x_sig.col(i+1) = x_aug + root_P.col(i)*sqrt(lambda_ + n_aug_);
    //adjust_angle(x_sig(3,i+1));
    x_sig.col(i+n_aug_+1) = x_aug - root_P.col(i)*sqrt(lambda_ + n_aug_);
    //adjust_angle(x_sig(3,i+n_aug_+1));
  }
  //This predicts the sigma points
  VectorXd temp(n_aug_);
  double px = 0, py = 0, v = 0, epsi = 0, epsi_dot = 0, neu_a = 0, neu_epsi = 0;
  for (int i = 0; i<x_sig.cols(); i++){
      temp = x_sig.col(i);
      px = temp(0);
      py = temp(1);
      v = temp(2);
      epsi = temp(3);
      epsi_dot = temp(4);
      neu_a = temp(5);
      neu_epsi = temp(6);
      Xsig_pred_.col(i) = temp.head(n_x_);
      if (fabs(epsi_dot) < EPSILON){
        Xsig_pred_.col(i) += (VectorXd(n_x_)<<v*cos(epsi)*delta_t,
                                            v*sin(epsi)*delta_t,
                                            0,
                                            epsi_dot*delta_t,
                                            0).finished();
        Xsig_pred_.col(i) += (VectorXd(n_x_)<<0.5*delta_t*delta_t*cos(epsi)*neu_a,
                                            0.5*delta_t*delta_t*sin(epsi)*neu_a,
                                            delta_t * neu_a,
                                            0.5 * delta_t * delta_t * neu_epsi,
                                            delta_t * neu_epsi).finished();
      }else{
        Xsig_pred_.col(i) += (VectorXd(n_x_)<<(v/epsi_dot)*(sin(epsi + epsi_dot * delta_t) - sin(epsi)),
                                            (v/epsi_dot)*(-cos(epsi + epsi_dot * delta_t) + cos(epsi)),
                                            0,
                                            epsi_dot*delta_t,
                                            0).finished();
        Xsig_pred_.col(i) += (VectorXd(n_x_)<<0.5 * delta_t * delta_t * cos(epsi) * neu_a,
                                            0.5 * delta_t * delta_t * sin(epsi) * neu_a,
                                            delta_t * neu_a,
                                            0.5 * delta_t * delta_t * neu_epsi,
                                            delta_t * neu_epsi).finished();
      }
      adjust_angle(Xsig_pred_(3,i));
  }
  W_.fill(0.0);
  W_(0,0) = lambda_/(lambda_ + n_aug_);
  for (int i = 1;i<2*n_aug_ + 1; i++) W_(i,i) = 0.5/(lambda_ + n_aug_); 
  x_ = Xsig_pred_ * W_.diagonal();
  adjust_angle(x_(3));
  for (int i = 0; i<Xsig_pred_.cols(); i++){
    delta_Xsig_pred_.col(i) = Xsig_pred_.col(i)-x_;
    adjust_angle(delta_Xsig_pred_(3,i));
  }
  cout<<"delta_Xsig_pred_"<<endl;
  cout<<delta_Xsig_pred_<<endl;
  P_ = delta_Xsig_pred_ * W_ * delta_Xsig_pred_.transpose(); 
  cout<<"Predicted P"<<endl;
  cout<<P_<<endl;
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
  
  double px=0,py=0;
  MatrixXd Zsig (2,Xsig_pred_.cols());
  
  
  MatrixXd H (2,n_x_);
  H.fill(0.0);
  H(0,0)=1;
  H(1,1)=1;
  Zsig = H*Xsig_pred_;
  
  
  VectorXd z_pred = Zsig * W_.diagonal();
  for (int i = 0; i<Zsig.cols(); i++) Zsig.col(i)-=z_pred;
  MatrixXd R (2,2);
  R.fill(0.0);
  R.diagonal()<<pow(std_laspx_,2), pow(std_laspy_,2);
  MatrixXd S = Zsig * W_ * Zsig.transpose() + R;
  MatrixXd T = delta_Xsig_pred_ * W_ * Zsig.transpose();
  MatrixXd K = T * S.inverse();
  x_ += K*(meas_package.raw_measurements_-z_pred);
  adjust_angle(x_(3));
  P_ -= K*S*K.transpose();
  cout<<"Updated P_ lidar:"<<endl;
  cout<<P_<<endl;
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
  
  double px=0,py=0,v=0,epsi=0,epsi_dot=0;
  MatrixXd Zsig (3,Xsig_pred_.cols());
  
  
  for (int i = 0; i<Xsig_pred_.cols(); i++){
    px=Xsig_pred_.col(i)(0);
    py=Xsig_pred_.col(i)(1);
    v=Xsig_pred_.col(i)(2);
    epsi=Xsig_pred_.col(i)(3);
    epsi_dot=Xsig_pred_.col(i)(4);
    Zsig.col(i)<<sqrt(px*px+py*py),atan2(py,px),(px*cos(epsi)*v + py*sin(epsi)*v)/(sqrt(px*px+py*py));
  }
  
  VectorXd z_pred = Zsig * W_.diagonal();
  adjust_angle(z_pred(1));
  for (int i = 0; i<Zsig.cols(); i++){
    Zsig.col(i)-=z_pred;
    adjust_angle(Zsig(1,i));
  }  
  MatrixXd R (3,3);
  R.fill(0.0);
  R.diagonal()<<pow(std_radr_,2), pow(std_radphi_,2), pow(std_radrd_,2);
  MatrixXd S = Zsig * W_ * Zsig.transpose() + R;
  

  MatrixXd T = delta_Xsig_pred_ * W_ * Zsig.transpose();
  
  MatrixXd K = T * S.inverse();
  
  x_ += K*(meas_package.raw_measurements_-z_pred);
  adjust_angle(x_(3));
  P_ -= K*S*K.transpose();
  cout<<"Updated P_ Radar"<<endl;
  cout<<P_<<endl;
}
