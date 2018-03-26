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
UKF::UKF() {
	
	is_initialized_ = false;
	// if this is false, laser measurements will be ignored (except during init)
	use_laser_ = true;

	// if this is false, radar measurements will be ignored (except during init)
	use_radar_ = true;
	
	//set state dimension
	n_x_ = 5;

	//set augmented dimension
	n_aug_ = 7;

	// initial state vector
	x_ = VectorXd(n_x_);

	// initial covariance matrix
	P_ = MatrixXd(n_x_, n_x_);

	// Process noise standard deviation longitudinal acceleration in m/s^2
	std_a_ = 1;

	// Process noise standard deviation yaw acceleration in rad/s^2
	std_yawdd_ = 0.5;
	
	// Noise matrix 
	Q_ = MatrixXd(n_aug_-n_x_,n_aug_-n_x_);
	
	Q_ << std_a_*std_a_, 0,
			0,      std_yawdd_*std_yawdd_;
			
	//define spreading parameter (design parameter)
	lambda_ = 3 - n_aug_;

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
	
	// time when the state is true, in us
	time_us_ = 0.0;

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
/*****************************************************************************
   *  Initialization
   ****************************************************************************/
  if (!is_initialized_) {
	 std::cout << "Initializing the UKF filter" << std::endl; 
	// Assuming initial states and Covariance matrix
      x_ << 1, 1, 2, 1, 0.1;

      // init covariance matrix
      P_ << 0.15,    0, 0, 0, 0,
               0, 0.15, 0, 0, 0,
               0,    0, 1, 0, 0,
               0,    0, 0, 1, 0,
			   0, 0, 0, 0, 1;
			   
	//create vector for weights
	weights_ = VectorXd(2*n_aug_+1);

	// set weights
	double weight_0 = lambda_/(lambda_+n_aug_);
	weights_(0) = weight_0;
	
	for (int i=1; i<2*n_aug_+1; i++) {  //2n+1 weights
		double weight = 0.5/(n_aug_+lambda_);
		weights_(i) = weight;
	}
	
	// initial time mark
	time_us_ = meas_package.timestamp_;

    if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
        /**
        Convert radar from polar to cartesian coordinates and initialize state.
        */
        float ro = meas_package.raw_measurements_(0);
        float phi = meas_package.raw_measurements_(1);
        float ro_dot = meas_package.raw_measurements_(2);
		
        x_(0) = ro * cos(phi);
		x_(1) = ro * sin(phi);
    }
    else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
		x_(0) = meas_package.raw_measurements_(0);
		x_(1) = meas_package.raw_measurements_(1);
	}

    // done initializing, no need to predict or update
    is_initialized_ = true;
	std::cout << "Initializion complete!" << std::endl;
	std::cout << "UKF algorithm running..." << std::endl;
	
	return;
  }
  
      //compute the time elapsed between the current and previous measurements
    delta_t_ = (meas_package.timestamp_ - time_us_) / 1000000.0;	//dt - expressed in seconds
	time_us_ = meas_package.timestamp_;

/**
 * Prediction and Update of the UKF filter
 */

	// Generate sigma points
	AugmentedSigmaPoints();
	// Predict sigma points for the next time step
	SigmaPointPrediction();
	// Compute the mean of the state and its Covariance matrix
	PredictMeanAndCovariance();
	
	// Predict and update sensor measurements
	if (meas_package.sensor_type_ == MeasurementPackage::RADAR && use_radar_== true) {
		PredictRadarMeasurement();
		VectorXd z = meas_package.raw_measurements_;
		UpdateState(z, 3); // 3 is the number of parameters measured from the radar
		
	}
	
	else if (meas_package.sensor_type_ == MeasurementPackage::LASER && use_laser_==true) {
		PredictLidarMeasurement();  
		VectorXd z = meas_package.raw_measurements_;
		UpdateState(z, 2);  // 2: p_x and p_y
		
	}
}

// FUNCTIONS FOR THE UKF


void UKF::AugmentedSigmaPoints() {

	//create augmented mean vector
	VectorXd x_aug = VectorXd(n_aug_);
	x_aug.fill(0.0);

	//create augmented state covariance
	MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);
	P_aug.fill(0.0);

	//create sigma point matrix
	Xsig_aug_ = MatrixXd(n_aug_, 2 * n_aug_ + 1);
	Xsig_aug_.fill(0.0);

	//create augmented mean state. The mean of the noises values is zero.
	x_aug.head(n_x_) = x_;
	x_aug(n_x_) = 0;
	x_aug(n_x_+1) = 0;


	//create augmented covariance matrix
	P_aug.topLeftCorner(n_x_,n_x_) = P_;

	P_aug.bottomRightCorner(n_aug_-n_x_,n_aug_-n_x_) = Q_;
	// std::cout << "Q" << Q_ << '\n';
	// std::cout << "P_aug:" << P_aug << '\n';

	//  create square root matrix
	MatrixXd L = P_aug.llt().matrixL();
	//std::cout << "square root of the P_aug matrix" << L << '\n';

	//create augmented sigma points
	// First column of the sigma points
	Xsig_aug_.col(0) = x_aug;
	double spreading_module = sqrt(lambda_ + n_aug_);

	//set remaining columns of sigma points (2nd to last column)

	for (int i = 0; i < n_aug_; i++)
	{
	Xsig_aug_.col(i+1)     = x_aug + spreading_module * L.col(i);
	Xsig_aug_.col(i+1+n_aug_) = x_aug - spreading_module * L.col(i);
	}
	
	//std::cout << P_ << std::endl;
	// *Xsig_out = Xsig_aug_;

}

void UKF::SigmaPointPrediction() {

   //create matrix with predicted sigma points as columns
	Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ + 1);
	Xsig_pred_.fill(0.0);

	// Predict the next states using the CTRV vehicle model
  for (int i = 0; i< 2*n_aug_+1; i++)
  {
    //extract values for better readability
    double p_x = Xsig_aug_(0,i);
    double p_y = Xsig_aug_(1,i);
    double v = Xsig_aug_(2,i);
    double yaw = Xsig_aug_(3,i);
    double yawd = Xsig_aug_(4,i);
    double nu_a = Xsig_aug_(5,i);
    double nu_yawdd = Xsig_aug_(6,i);

    //predicted state values
    double px_p, py_p;

    //avoid division by zero
    if (fabs(yawd) > 0.001) {
        px_p = p_x + v/yawd * ( sin (yaw + yawd*delta_t_) - sin(yaw));
        py_p = p_y + v/yawd * ( cos(yaw) - cos(yaw+yawd*delta_t_) );
    }
    else {
        px_p = p_x + v*delta_t_*cos(yaw);
        py_p = p_y + v*delta_t_*sin(yaw);
    }

    double v_p = v;
    double yaw_p = yaw + yawd*delta_t_;
    double yawd_p = yawd;
	
	long long delta_t2_ = 0.5*delta_t_*delta_t_;

    //add noise
    px_p = px_p + delta_t2_ * nu_a * cos(yaw);
    py_p = py_p + delta_t2_ * nu_a * sin(yaw);
    v_p = v_p + nu_a*delta_t_;

    yaw_p = yaw_p + delta_t2_ * nu_yawdd;
    yawd_p = yawd_p + nu_yawdd*delta_t_;

    //write predicted sigma point into right column
    Xsig_pred_(0,i) = px_p;
    Xsig_pred_(1,i) = py_p;
    Xsig_pred_(2,i) = v_p;
    Xsig_pred_(3,i) = yaw_p;
    Xsig_pred_(4,i) = yawd_p;
	
  }


  //write result
  //*Xsig_out = Xsig_pred_;

}

void UKF::PredictMeanAndCovariance() {
	
	VectorXd x = VectorXd(n_x_);
	MatrixXd P = MatrixXd(n_x_,n_x_);

	//predicted state mean
	x.fill(0.0);
	for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //iterate over sigma points
	x = x+ weights_(i) * Xsig_pred_.col(i);
	
	}

	//predicted state covariance matrix
	P.fill(0.0);
	for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //iterate over sigma points

    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x;
    //angle normalization
    while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
    while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;

    P = P + weights_(i) * x_diff * x_diff.transpose() ;
	
	x_ = x;
	P_ = P;
  }

	//write result
	x_ = x;
	P_ = P;
	
}

void UKF::PredictRadarMeasurement() {

	int n_z = 3; // rho, phi, rho_dot
	
	//create matrix for sigma points in measurement space
	Zsig_ = MatrixXd(n_z, 2 * n_aug_ + 1);

	//transform sigma points into measurement space
	for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points

		// extract values for better readibility
		double p_x = Xsig_pred_(0,i);
		double p_y = Xsig_pred_(1,i);
		double v  = Xsig_pred_(2,i);
		double yaw = Xsig_pred_(3,i);

		double v1 = cos(yaw)*v;
		double v2 = sin(yaw)*v;

		// measurement model
		Zsig_(0,i) = sqrt(p_x*p_x + p_y*p_y);                        //r
		Zsig_(1,i) = atan2(p_y,p_x);                                 //phi
		Zsig_(2,i) = (p_x*v1 + p_y*v2 ) / sqrt(p_x*p_x + p_y*p_y);   //r_dot
	}

	//mean predicted measurement
	z_pred_ = VectorXd(n_z);
	z_pred_.fill(0.0);
	
	for (int i=0; i < 2*n_aug_+1; i++) {
		z_pred_ = z_pred_ + weights_(i) * Zsig_.col(i);
	}

	//innovation covariance matrix S_
	S_ = MatrixXd(n_z,n_z);
	S_.fill(0.0);
	
	for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 sigma points
		//residual
		VectorXd z_diff = Zsig_.col(i) - z_pred_;

		//angle normalization
		while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
		while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

		S_ = S_ + weights_(i) * z_diff * z_diff.transpose();
	}

	//add measurement noise covariance matrix
	MatrixXd R = MatrixXd(n_z,n_z);
	
	R <<    std_radr_*std_radr_, 0, 0,
			  0, std_radphi_*std_radphi_, 0,
			  0, 0,std_radrd_*std_radrd_;

	//std::cout << "R: "  << std::endl << R << std::endl;
	
	S_ = S_ + R;

	//write result
	//*z_out = z_pred_;
	//*S_out = S_;
}  

void UKF::PredictLidarMeasurement() {

	int n_z = 2; // p_x, p_y
	
	//create matrix for sigma points in measurement space
	Zsig_ = MatrixXd(n_z, 2 * n_aug_ + 1);

	//transform sigma points into measurement space
	for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points

		// extract values for better readibility
		double p_x = Xsig_pred_(0,i);
		double p_y = Xsig_pred_(1,i);

		// measurement model
		Zsig_(0,i) = p_x;                        
		Zsig_(1,i) = p_y;                               
	}

	//mean predicted measurement
	z_pred_ = VectorXd(n_z);
	z_pred_.fill(0.0);
	
	for (int i=0; i < 2*n_aug_+1; i++) {
		z_pred_ = z_pred_ + weights_(i) * Zsig_.col(i);
	}

	//innovation covariance matrix S_
	S_ = MatrixXd(n_z,n_z);
	S_.fill(0.0);
	
	for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 sigma points
		//residual
		VectorXd z_diff = Zsig_.col(i) - z_pred_;

		S_ = S_ + weights_(i) * z_diff * z_diff.transpose();
	}

	//add measurement noise covariance matrix
	MatrixXd R = MatrixXd(n_z,n_z);
	
	R <<    std_laspx_*std_laspx_, 0,
			  0, std_laspy_*std_laspy_;

	//std::cout << "R: "  << std::endl << R << std::endl;
	
	S_ = S_ + R;

	//write result
	//*z_out = z_pred_;
	//*S_out = S_;
}  



void UKF::UpdateState(VectorXd z, int n_z){
	 
	//create matrix for cross correlation Tc
	MatrixXd Tc = MatrixXd(n_x_, n_z);

	//calculate cross correlation matrix
	Tc.fill(0.0);
	for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points

		//residual
		VectorXd z_diff = Zsig_.col(i) - z_pred_;
		//angle normalization
		while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
		while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

		// state difference
		VectorXd x_diff = Xsig_pred_.col(i) - x_;
		//angle normalization
		while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
		while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;

		Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
	}

	//Kalman gain K;
	MatrixXd K = Tc * S_.inverse();

	//residual
	VectorXd z_diff = z - z_pred_;

	//angle normalization
	while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
	while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

	//update state mean and covariance matrix
	x_ = x_ + K * z_diff;
	P_ = P_ - K*S_*K.transpose();
	
  
}