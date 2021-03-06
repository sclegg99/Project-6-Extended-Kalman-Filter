#include "FusionEKF.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/*
 * Constructor.
 */
FusionEKF::FusionEKF() {
  is_initialized_ = false;

  previous_timestamp_ = 0;

  // initializing matrices
  R_laser_ = MatrixXd(2, 2);
  R_radar_ = MatrixXd(3, 3);
  H_laser_ = MatrixXd(2, 4);
  Hj_ = MatrixXd(3, 4);

  //measurement covariance matrix - laser
  R_laser_ << 0.0225, 0,
        0, 0.0225;

  //measurement covariance matrix - radar
  R_radar_ << 0.09, 0, 0,
        0, 0.0009, 0,
        0, 0, 0.09;

  //initialize H_laser_
  H_laser_ << 1., 0., 0., 0.,
              0., 1., 0., 0.;

  /**
  TODO:
    * Finish initializing the FusionEKF.
    * Set the process and measurement noises
  */
  //initialize transition matrix F_
  ekf_.F_ = MatrixXd(4,4);
  ekf_.F_ << 1., 0., 0., 0.,
             0., 1., 0., 0.,
             0., 0., 1., 0.,
             0., 0., 0., 1.;

  //initialize process noise Q_
  ekf_.Q_ = MatrixXd(4,4);
  ekf_.Q_ << 0., 0., 0., 0.,
             0., 0., 0., 0.,
             0., 0., 0., 0.,
             0., 0., 0., 0;

  //initialize state covariance matrix P_
  ekf_.P_ = MatrixXd(4,4);
  ekf_.P_ << 1., 0., 0., 0.,
             0., 1., 0., 0.,
             0., 0., 1000., 0.,
             0., 0., 0., 1000.;
   
  //
  // Define values for process noise
  noise_ax = 9;
  noise_ay = 9;
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
    /**
    TODO:
      * Initialize the state ekf_.x_ with the first measurement.
      * Create the covariance matrix.
      * Remember: you'll need to convert radar from polar to cartesian coordinates.
    */
    // first measurement
    cout << "Initialize FusionEKF: " << endl;
    ekf_.x_ = VectorXd(4);
 
    // extract first of the raw measurement data m1, and m2
    float m1 = measurement_pack.raw_measurements_(0);
    float m2 = measurement_pack.raw_measurements_(1);
    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      /**
      Convert radar from polar to cartesian coordinates and initialize state.
      */
      //initial velocity vx, vy are assumed 0.
      ekf_.x_ << m1*cos(m2), m1*sin(m2), 0., 0.;
    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      /**
      Initialize state.
      */
      ekf_.x_ << m1, m2, 0., 0.;
    }

    // initialize the time stamp
    previous_timestamp_ = measurement_pack.timestamp_;

    cout << "x_ = " << ekf_.x_ << endl;
    // done initializing, no need to predict or update
    is_initialized_ = true;
    return;
  }
/*
  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      cout << "Radar measurement...skip" << endl;
      return;
  }
*/
  /*****************************************************************************
   *  Prediction
   ****************************************************************************/

  /**
   TODO:
     * Update the state transition matrix F according to the new elapsed time.
      - Time is measured in seconds.
     * Update the process noise covariance matrix.
     * Use noise_ax = 9 and noise_ay = 9 for your Q matrix.
   */
  // Calculate the time increment from the previous measurement to the
  // current measurement.  Then store the current time as the previous time
  float dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;	//dt - expressed in seconds
  previous_timestamp_ = measurement_pack.timestamp_;

  // set up the F and Q matrix for the time increment dt
  ekf_.F_(0,2) = dt;
  ekf_.F_(1,3) = dt;

  float dt_2 = dt*dt;
  float dt_3 = 0.5*dt*dt_2;
  float dt_4 = 0.5*dt*dt_3;
  ekf_.Q_(0,0) = dt_4*noise_ax;
  ekf_.Q_(0,2) = dt_3*noise_ax;
  ekf_.Q_(1,1) = dt_4*noise_ay;
  ekf_.Q_(1,3) = dt_3*noise_ay;
  ekf_.Q_(2,0) = dt_3*noise_ax;
  ekf_.Q_(2,2) = dt_2*noise_ax;
  ekf_.Q_(3,1) = dt_3*noise_ay;
  ekf_.Q_(3,3) = dt_2*noise_ay;

  // predict new state x_ and covariance P_
  ekf_.Predict();

  /*****************************************************************************
   *  Update
   ****************************************************************************/

  /**
   TODO:
     * Use the sensor type to perform the update step.
     * Update the state and covariance matrices.
   */

  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    // Radar updates
    // calcluate Jacobian for current state
    Hj_ = tools.CalculateJacobian(ekf_.x_);

    // set H_ to Jacobian and R_ to radar value and call update
    ekf_.H_ = Hj_;
    ekf_.R_ = R_radar_;
    ekf_.UpdateEKF(measurement_pack.raw_measurements_);
  } else {
    // Laser updates
    // set H_ and R_ to laser values then call update
    ekf_.H_ = H_laser_;
    ekf_.R_ = R_laser_;
    ekf_.Update(measurement_pack.raw_measurements_);
  }

  // print the output
  cout << "x_ = " << ekf_.x_ << endl;
  cout << "P_ = " << ekf_.P_ << endl;
}
