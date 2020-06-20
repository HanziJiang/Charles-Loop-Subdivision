#include "catmull_clark.h"
#include "vertex_triangle_adjacency.h"
#include <unordered_map>
#include <utility>
#include <functional>
#include <vector>
#include <iostream>

void catmull_clark(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  const int num_iters,
  Eigen::MatrixXd & SV,
  Eigen::MatrixXi & SF)
{
  ////////////////////////////////////////////////////////////////////////////
  // Reference: https://en.wikipedia.org/wiki/Catmull–Clark_subdivision_surface
  std::cout << "begin\n";
  if (num_iters <= 0) return;

  const int V_rows = V.rows();
  const int F_rows = F.rows();

  int SV_size = V_rows + F_rows;
  int SF_size = F_rows * 4;

  Eigen::Matrix3d face_points = Eigen::MatrixXd::Zero(F_rows, 3);
  for (int i = 0; i < F_rows; i++) {
    // For each face, add face point
    face_points.row(i) = (V.row(F(i, 0)) + V.row(F(i, 1)) + V.row(F(i, 2)) + V.row(F(i, 3))) / 4.0;
  }

  // For each face point, add edges from face point to adjacent edge point.
  std::unordered_map<int, std::vector<int>> new_edges;  // vertex index & indices of points sharing edge with vertex
  int vertex_index;
  // vertex index & edge points of the edges vertex-edges[vertex index]
  std::unordered_map<int, std::vector<Eigen::RowVector3d>> edge_points;
  int target_index_1, target_index_2;

  for (int i = 0; i < F_rows; i++) {
    for (int j = 0; j < 4; j++) {
      vertex_index = F(i, j);

      // if key not exist
      if (new_edges.find(vertex_index) == new_edges.end()) {
        std::vector<int> edge_entry;
        new_edges[vertex_index] = edge_entry;
        std::vector<Eigen::RowVector3d> edge_point_entry;
        edge_points[vertex_index] = edge_point_entry;
      }

      target_index_1 = F(i, (j + 1) % 4);
      target_index_2 = F(i, (j - 1 + 4) % 4);

      // if vertex_index-target_index_1 is not already stored
      if (std::find(new_edges[vertex_index].begin(), new_edges[vertex_index].end(), target_index_1) == new_edges[vertex_index].end() && (new_edges.find(target_index_1) == new_edges.end() || std::find(new_edges[target_index_1].begin(), new_edges[target_index_1].end(), vertex_index) == new_edges[target_index_1].end())) {
        new_edges[vertex_index].push_back(target_index_1);
        edge_points[vertex_index].push_back(V.row(target_index_1) + V.row(vertex_index) / 2.0);
        SV_size ++;
      }
      
      // if vertex_index-target_index_2 is not already stored
      if (std::find(new_edges[vertex_index].begin(), new_edges[vertex_index].end(), target_index_2) == new_edges[vertex_index].end() && (new_edges.find(target_index_2) == new_edges.end() || std::find(new_edges[target_index_2].begin(), new_edges[target_index_2].end(), vertex_index) == new_edges[target_index_2].end())) {
        new_edges[vertex_index].push_back(target_index_2);
        edge_points[vertex_index].push_back(V.row(target_index_2) + V.row(vertex_index) / 2.0);
        SV_size ++;
      }
      
    }
  }

  SV = Eigen::MatrixXd::Zero(SV_size, 3);
  SF = Eigen::MatrixXi::Zero(SF_size, 4);

  std::vector<std::vector<int> > VF;
  std::vector<int> face_indces;
  vertex_triangle_adjacency(F, V_rows, VF);

  Eigen::RowVector3d f(0, 0, 0), r(0, 0, 0);
  int num_edges;
  Eigen::RowVector3d barycenter(0, 0, 0);
  int num_points;
  for (int i = 0; i < V_rows; i++) {
    // calculate f, the average of all face points of faces adjacent to vertex at V.row(i)
    face_indces = VF[i];
    for (int face_index : face_indces) {
      f += face_points.row(face_index);
    }
    num_points = face_indces.size();

    // calculate r, the average of all edge points of edges touching vertex at V.row(i)
    num_edges = new_edges[i].size();
    for (int j = 0; j < num_edges; j++) {
      r += edge_points[i].at(j);
    }
    for (int j = 0; j < VF[i].size(); j++) {
      for (int k = 0; j < 4; k++) {
        for (int w = 0; w < new_edges[F(j, k)].size(); w++) {
          if (new_edges[F(j, k)].at(w) == i) {
            r += edge_points[F(j, k)].at(w);
            num_edges ++;
          }
        }
      }
    }
    r /= num_edges;
    
    barycenter = (f + 2 * r + (num_edges - 3) * V.row(i)) / num_edges;

    SV.row(i) = barycenter;
  }

  for (int i = 0; i < F_rows; i++) {
    SV.row(V_rows + i) = face_points.row(i);
  }

  int size;
  int index = V_rows + F_rows;
  for (int i = 0; i < V_rows; i++) {
    size = edge_points[i].size();
    for (int j = 0; j < size; j++) {
      SV.row(index) = edge_points[i].at(j);
      index++;
    }
  }

  int point_index_1, point_index_2, point_index_3, point_index_4, another_index_1, another_index_2; 
  int count_1, count_2;
  for (int i = 0; i < F_rows; i++) {
    for (int j = 0; j < 4; j++) {
      point_index_1 = F(i, j);
      point_index_3 = V_rows + i;

      another_index_1 = F(i, (j + 1) % 4);
      another_index_1 = F(i, (j - 1 + 4) % 4);
      count_1, count_2 = -1;
      for (int a = 0; a < new_edges.size() && (count_1 < 0 || count_2 < 0); i++) {
        for (int b = 0; b < new_edges[a].size() && (count_1 < 0 || count_2 < 0); b++) {
          if ((a == point_index_1 && new_edges[a].at(b) == another_index_1) || (a == another_index_1 && new_edges[a].at(b) == point_index_1)) {
            count_1 = a * new_edges[a].size() + b;
          }

          if ((a == point_index_2 && new_edges[a].at(b) == another_index_2) || (a == another_index_2 && new_edges[a].at(b) == point_index_2)) {
            count_2 = a * new_edges[a].size() + b;
          }
        }
      }

      point_index_2 = F_rows + V_rows + count_1;
      point_index_4 = F_rows + V_rows + count_2;
      
      SF.row(i) = Eigen::RowVector4i(point_index_1, point_index_2, point_index_3, point_index_4);
    }
  }

  catmull_clark(Eigen::MatrixXd(SV), Eigen::MatrixXi(SF), num_iters - 1, SV, SF);
  std::cout << "SV size:" << SV.size() << "\n";
  std::cout << "SF size:" << SF.size() << "\n";
  ////////////////////////////////////////////////////////////////////////////
}
