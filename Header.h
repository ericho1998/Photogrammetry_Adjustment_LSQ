#ifndef HEADER_H
#define HEADER_H
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <math.h>
#include <vector>
#include<iomanip>
#define EIGEN_NO_DEBUG


#include <Eigen/Dense> //for Eigen library

using namespace std;
using namespace Eigen;

using Eigen::MatrixXd;
using Eigen::VectorXd;

void read(const string& f, MatrixXd& x); //read function
void write(const string& f, MatrixXd& x); // matrix output function
void writeicf(const string& f, MatrixXd& x); //write an icf file
void writegcp(const string& f, MatrixXd& x2); // write a gcp file
void writeori(const string& f, MatrixXd& eop); // create ori file
void writecam(const string& f);
MatrixXd interpolation(MatrixXd& gcps, MatrixXd& obs, MatrixXd& x, double& xp, double& yp);// interpolation to remove radial distortion
MatrixXd Amat(MatrixXd& gcp);
MatrixXd lmat(MatrixXd& meas);
double z(MatrixXd& meas);
#endif