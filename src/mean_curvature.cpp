#include "../include/mean_curvature.h"
#include <igl/cotmatrix.h>
#include <igl/massmatrix.h>
#include <igl/invert_diag.h>
#include <igl/per_vertex_normals.h>

void mean_curvature(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  Eigen::VectorXd & H)
{
  // Replace with your code
  Eigen::SparseMatrix<double> M, L, Minv;
  igl::massmatrix(V, F, igl::MASSMATRIX_TYPE_DEFAULT, M);
  igl::cotmatrix(V, F, L);
  igl::invert_diag(M, Minv);
  
  Eigen::MatrixXd curveNormals = Minv * L * V;

  Eigen::MatrixXd vertexNormals;
  igl::per_vertex_normals(V, F, vertexNormals);

  // Stripping the magnitude off the rows of the resulting matrix would give the unsigned mean curvature.
  // To make sure that the sign is preserved we can check 
  // whether each row in H agrees or disagrees with consistently oriented per-vertex normals in N
  H.resize(V.rows());
  for (int i = 0; i < V.rows(); i++){
    double test = curveNormals.row(i).dot(vertexNormals.row(i));
    H(i) = curveNormals.row(i).norm();
    if (test < 0){
      H(i) = -H(i);
    }
  }

}
