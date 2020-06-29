#include "sphere.h"
#include <iostream>

void sphere(
  const int num_faces_u,
  const int num_faces_v,
  Eigen::MatrixXd & V,
  Eigen::MatrixXi & F,
  Eigen::MatrixXd & UV,
  Eigen::MatrixXi & UF,
  Eigen::MatrixXd & NV,
  Eigen::MatrixXi & NF)
{
  ////////////////////////////////////////////////////////////////////////////

  // Reference: https://stackoverflow.com/questions/12732590/how-map-2d-grid-points-x-y-onto-sphere-as-3d-points-x-y-z

  const int num_vertex = (num_faces_u + 1) * (num_faces_v + 1);
  const int num_face = num_faces_u * num_faces_v;

  V.resize(num_vertex, 3);
  F.resize(num_face, 4);
  UV.resize(num_vertex, 2);
  UF.resize(num_face, 4);
  NV.resize(num_vertex, 3);
  NF.resize(num_face, 4);

  double x, y, z, longitude, latitude, p1, p2, p3, p4;
  int index;

  for (int i = 0; i <= num_faces_u; i++) {
    for (int j = 0; j <= num_faces_v; j++) {
      index = (num_faces_v + 1) * i + j;

      latitude = M_PI * (double) j / num_faces_v;
      longitude = 2 * M_PI * (double) i / num_faces_u;
      z = - sin(latitude) * cos(longitude);
      x = - sin(latitude) * sin(longitude);
      y = - cos(latitude);

      // update V, UV, NV
      V.row(index) = Eigen::RowVector3d(x, y, z);
      UV.row(index) = Eigen::RowVector2d((double) i / (num_faces_u + 1), (double) j / (num_faces_v + 1));
      NV.row(index) = Eigen::RowVector3d(x, y, z).normalized();
      
      // update F, UF, NF
      if (i != num_faces_u && j != num_faces_v) {
        index = num_faces_v * i + j;

        p1 = (num_faces_v + 1) * i + j;
        p2 = (num_faces_v + 1) * (i + 1) + j;
        p3 = (num_faces_v + 1) * (i + 1) + j + 1;
        p4 = (num_faces_v + 1) * i + j + 1;

        F.row(index) = Eigen::RowVector4i(p1, p2, p3, p4);
        UF.row(index) = Eigen::RowVector4i(p1, p2, p3, p4);
        NF.row(index) = Eigen::RowVector4i(p1, p2, p3, p4); 
      }
    }
  }

  ////////////////////////////////////////////////////////////////////////////
}
