#include "vertex_triangle_adjacency.h"

void vertex_triangle_adjacency(
  const Eigen::MatrixXi & F,
  const int num_vertices,
  std::vector<std::vector<int> > & VF)
{
  VF.resize(num_vertices);
  ////////////////////////////////////////////////////////////////////////////
  const int rows = F.rows();
  std::vector<int> VF_entry;

  for (int i = 0; i < num_vertices; i++) {
    VF_entry.clear();
    for (int j = 0; j < rows; j++) {
      if (F(j, 0) == i || F(j, 1) == i || F(j, 2) == i) {
        VF_entry.push_back(j);
        break;
      }
    }
    VF[i] = VF_entry;
  }
  ////////////////////////////////////////////////////////////////////////////
}

