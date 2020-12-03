#include "../include/internal_angles.h"
#include <math.h>
void internal_angles(
  const Eigen::MatrixXd & l_sqr,
  Eigen::MatrixXd & A)
{
  // Add with your code
  int rows = l_sqr.rows();
  A.resize(rows, 3);

  for (int row = 0; row < rows; row++){
    for (int col = 0; col < 3; col++){
      double a = l_sqr(row, col % 3);
      double b = l_sqr(row, (col + 1) % 3);
      double c = l_sqr(row, (col + 2) % 3);
      double cosTheta = (b+c-a)/(2.0*sqrt(b)*sqrt(c));
      A(row, col) = acos(cosTheta);
    }
  }
  
}
