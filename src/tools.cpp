#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  VectorXd RMSE;
  if (estimations.size()!=ground_truth.size()){
    std::cout<<"The length of estimations is not the same as that of the ground truth"<< std::endl;
    return RMSE;
  }
  if (estimations.size()==0){
    std::cout<<"The length of estimations and ground truth is zero"<<std::endl;
    return RMSE;
  }
  RMSE = VectorXd(estimations[0].size());
  RMSE.fill(0.0);
  VectorXd diff(estimations[0].size());
  for (int i =0; i<estimations.size(); i++){
    diff = estimations[i] - ground_truth[i];
    RMSE += diff.array().square().matrix();
  }
  RMSE /= estimations.size();
  RMSE = RMSE.array().sqrt().matrix();
  return RMSE;
  /**
  TODO:
    * Calculate the RMSE here.
  */

}