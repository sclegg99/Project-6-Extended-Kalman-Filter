#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  /**
  TODO:
    * Calculate the RMSE here.
  */
    VectorXd rmse(4);
    rmse << 0,0,0,0;

    // check the validity of the following inputs:
    //  * the estimation vector size should not be zero
    //  * the estimation vector size should equal ground truth vector size
    if(estimations.size() != ground_truth.size()
       || estimations.size() == 0){
        cout << "Invalid estimation or ground_truth data" << endl;
        return rmse;
    }

    //accumulate squared residuals
    for(unsigned int i=0; i < estimations.size(); ++i){

        VectorXd residual = estimations[i] - ground_truth[i];

        //coefficient-wise multiplication
        residual = residual.array()*residual.array();
        rmse += residual;
    }

    //calculate the mean
    rmse = rmse/estimations.size();

    //calculate the squared root
    rmse = rmse.array().sqrt();

    //return the result
    return rmse;
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
  /**
  TODO:
    * Calculate a Jacobian here.
  */

    // Initialize Hj
    MatrixXd Hj(3,4);

    // Store state in local variables p and v
    float px = x_state(0);
    float py = x_state(1);
    float vx = x_state(2);
    float vy = x_state(3);

    // Compute px**2 + py**2 and check if is small
    float pxpy = px*px+py*py;
    if(pxpy < 1.e-10) {
        cout << "Divide by zero in Jacobian" << endl;
        Hj << 0., 0., 0., 0.,
              0., 0., 0., 0.,
              0., 0., 0., 0;
        return Hj;
    }

    // Compute a inverse of pxpy and square root of the inverse
    pxpy = 1./pxpy;
    float sqrt_pxpy = sqrt(pxpy);

    // Now assign values to Hj
    Hj << px*sqrt_pxpy, py*sqrt_pxpy, 0., 0.,
         -py*pxpy,      px*pxpy,      0., 0.,
         (py*(vx*py-vy*px))*pxpy*sqrt_pxpy, (px*(vy*px-vx*py))*pxpy*sqrt_pxpy, px*sqrt_pxpy, py*sqrt_pxpy;

    return Hj;
}
