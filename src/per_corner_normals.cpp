#include "per_corner_normals.h"
#include "triangle_area_normal.h"
#include "vertex_triangle_adjacency.h"
#include <iostream>

void per_corner_normals(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  const double corner_threshold,
  Eigen::MatrixXd & N)
{
  N = Eigen::MatrixXd::Zero(F.rows()*3,3);
  ////////////////////////////////////////////////////////////////////////////
  std::vector<std::vector<int> > VF;
  vertex_triangle_adjacency(F, V.rows(), VF);

  const int F_rows = F.rows();
  // Namely: the final normal at N.row(i * 3 + j); the normal of a face; the normal of a face that is adjacent to normal1's face
  Eigen::RowVector3d normal, normal1, normal2;
  std::vector<int> faces;
  // The absolute deviation between normal1 and normal2 in degree
  double angle;
  // The absolute value of threshold
  const double threshold = abs(corner_threshold);
  
  for (int i = 0; i < F_rows; i++) {
    normal1 = triangle_area_normal(V.row(F(i, 0)), V.row(F(i, 1)), V.row(F(i, 2)));

    for (int j = 0; j < 3; j++) {

      normal = Eigen::RowVector3d(0, 0, 0);
      faces = VF[F(i, j)];

      for (int face_index : faces) {

        normal2 = triangle_area_normal(V.row(F(face_index, 0)), V.row(F(face_index, 1)), V.row(F(face_index, 2)));
        angle = abs(acos(normal1.dot(normal2) / normal1.norm() / normal2.norm()));

        if (angle > threshold) {
          normal += normal2;
        }
      }
      N.row(i * 3 + j) = normal.normalized();
    }
  }
  ////////////////////////////////////////////////////////////////////////////
}
