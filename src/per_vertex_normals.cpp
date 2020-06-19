#include "per_vertex_normals.h"
#include "triangle_area_normal.h"

void per_vertex_normals(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  Eigen::MatrixXd & N)
{
  N = Eigen::MatrixXd::Zero(V.rows(),3);
  ////////////////////////////////////////////////////////////////////////////
  const int V_rows = V.rows();
  const int F_rows = F.rows();

  Eigen::RowVector3d n;
  for (int i = 0; i < V_rows; i++) {
    n = Eigen::RowVector3d(0, 0, 0);
    for (int j = 0; j < F_rows; j++) {
      // if the triangle at F.row(j) has V.row(i) as one of its vertices, add the area-weighted normal of the triangle to n
      if (F(j, 0) == i || F(j, 1) == i || F(j, 2) == i) {
        n += triangle_area_normal(V.row(F(i, 0)), V.row(F(i, 1)), V.row(F(i, 2)));
      }
    }
    N.row(i) = n.normalized();
  }
  ////////////////////////////////////////////////////////////////////////////
}
