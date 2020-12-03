#include "../include/angle_defect.h"
#include <igl/squared_edge_lengths.h>
#include "internal_angles.h"
#include <math.h>
void angle_defect(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  Eigen::VectorXd & D)
{ 
  // Set all values initially to 2*pi we will subtract afterwards
  D = Eigen::VectorXd::Ones(V.rows());
  D *= (2.0* M_PI);

  Eigen::MatrixXd l_sqr;
	igl::squared_edge_lengths(V, F, l_sqr);
	Eigen::MatrixXd A;
	internal_angles(l_sqr, A);

  for(int row = 0; row < F.rows(); row++){
    D(F(row,0)) -= A(row, 0);
    D(F(row,1)) -= A(row, 1);
    D(F(row,2)) -= A(row, 2);
  }
}
