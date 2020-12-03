#include "../include/principal_curvatures.h"
#include <igl/adjacency_matrix.h>
#include <igl/per_vertex_normals.h>
#include <igl/pinv.h>
#include <math.h>


void principal_curvatures(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  Eigen::MatrixXd & D1,
  Eigen::MatrixXd & D2,
  Eigen::VectorXd & K1,
  Eigen::VectorXd & K2)
{
  // Replace with your code
  int rows = V.rows();
  K1 = Eigen::VectorXd(rows);
  K2 = Eigen::VectorXd(rows);
  D1 = Eigen::MatrixXd(rows, 3);
  D2 = Eigen::MatrixXd(rows, 3);

  // Normals for each vertex
	Eigen::MatrixXd N;
	igl::per_vertex_normals(V, F, N);

  // This produces the one ring relative to the vertex
  Eigen::SparseMatrix<int> AdjExactlyOne;
  igl::adjacency_matrix(F, AdjExactlyOne);
  // This produces the connection that are exactly two away
  Eigen::SparseMatrix<int> AdjExactlyTwo;
  AdjExactlyTwo = AdjExactlyOne * AdjExactlyOne;
  // If we add both, it will produce the desired "two ring" of v (disregarding diagonal)
  Eigen::SparseMatrix<int> A = AdjExactlyOne + AdjExactlyTwo;
  // This will have non-zero elements on the diagonal, but we don't want that
  for (int row = 0; row<A.rows(); row++){
      A.coeffRef(row,row) = 0;
  }

  for (int vertex = 0 ; vertex< V.rows(); vertex++){
    // Gather the positions of the k points relative to v (i.e., vi-v) into a matrix P
    int k = A.innerVector(vertex).nonZeros();
    // Fill P
    Eigen::MatrixXd P(k, 3);
    int row = 0;
		for (Eigen::SparseMatrix<int>::InnerIterator it(A, vertex); it; ++it) {
			P.row(row) = V.row(it.row()) - V.row(vertex);
			row++;
		}

    // PCA (eigen decomp on P.transpose()*P)
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> PCA(P.transpose() * P);
    Eigen::MatrixXd eigenV = PCA.eigenvectors();
    Eigen::Vector3d u = eigenV.col(2);
    Eigen::Vector3d v = eigenV.col(1);
    Eigen::Vector3d w = eigenV.col(0);

    // UV to find coefficients
    Eigen::MatrixXd UV(P.rows(), 5);
    UV.col(0) = P * u;
    UV.col(1) = P * v;
    UV.col(2) = UV.col(0).cwiseProduct(UV.col(0));
	  UV.col(3) = UV.col(0).cwiseProduct(UV.col(1));
	  UV.col(4) = UV.col(1).cwiseProduct(UV.col(1));

    // solve for coefficients
    Eigen::MatrixXd UV_inv;
		igl::pinv(UV, UV_inv);
		Eigen::VectorXd coeff = UV_inv * P * w;

    // Shape Operator 
    double E = 1 + coeff(0) * coeff(0);
    double F = coeff(1) * coeff(0);
    double G = 1 + coeff(1) * coeff(1);
    double denominator = sqrt(coeff(0) * coeff(0) + 1 + coeff(1) * coeff(1));
    double e = 2.0 * coeff(2) / denominator;
    double f = coeff(3) / denominator;
    double g = 2.0 * coeff(4) / denominator;

    Eigen::MatrixXd leftMatrixS(2, 2);
    leftMatrixS << e, f,
                  f, g;
    Eigen::MatrixXd rightMatrixS(2, 2);
    rightMatrixS << E, F, 
                   F, G;
    Eigen::MatrixXd S = -leftMatrixS * rightMatrixS.inverse();

    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solved_shape(S);
    Eigen::MatrixXd shape_eigenVector = solved_shape.eigenvectors();
    Eigen::VectorXd shape_eigenValues = solved_shape.eigenvalues();
    // K1 stores max, K2 stores minimum
    // I verified that shape_eigenValues is in increasing order ie: shape_eigenValues(0) < shape_eigenValues(1)
    int k1 = 1;
    int k2 = 0;
    K1(vertex) = shape_eigenValues(k1);
    K2(vertex) = shape_eigenValues(k2);

    D1.row(vertex) = shape_eigenVector(1,1)*u + shape_eigenVector(1,0)*v;
    D2.row(vertex) = shape_eigenVector(0,1)*u + shape_eigenVector(0,0)*v;
  }


}
