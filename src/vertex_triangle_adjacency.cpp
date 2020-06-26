#include "vertex_triangle_adjacency.h"

void vertex_triangle_adjacency(
  const Eigen::MatrixXi & F,
  const int num_vertices,
  std::vector<std::vector<int> > & VF)
{
  VF.resize(num_vertices);
  ////////////////////////////////////////////////////////////////////////////
  const int rows = F.rows();
  for (int i = 0; i < num_vertices; i++) {
    std::vector<int> VF_entry;
    VF.at(i) = VF_entry;
    for (int j = 0; j < rows; j++) {
      if (F(j, 0) == i || F(j, 1) == i || F(j, 2) == i) {
        VF.at(i).push_back(j);
      }
    }
  }
  ////////////////////////////////////////////////////////////////////////////
}

