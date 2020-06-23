#include "catmull_clark.h"
#include "vertex_triangle_adjacency.h"
#include <unordered_map>
#include <utility>
#include <functional>
#include <vector>
#include <iostream>

struct hash_pair { 
    template <class T1, class T2> 
    size_t operator()(const std::pair<T1, T2>& p) const
    { 
        auto hash1 = std::hash<T1>{}(p.first); 
        auto hash2 = std::hash<T2>{}(p.second); 
        return hash1 ^ hash2; 
    } 
}; 

void catmull_clark(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  const int num_iters,
  Eigen::MatrixXd & SV,
  Eigen::MatrixXi & SF)
{
  ////////////////////////////////////////////////////////////////////////////
  // Reference: https://en.wikipedia.org/wiki/Catmull–Clark_subdivision_surface
  // https://www.geeksforgeeks.org/how-to-create-an-unordered_map-of-pairs-in-c/
  
  if (num_iters <= 0) return;

  const int V_rows = V.rows();
  const int F_rows = F.rows();

  int SV_size = V_rows + F_rows;
  int SF_size = F_rows * 4;

  Eigen::MatrixXd face_points = Eigen::MatrixXd::Zero(F_rows, 3);

  for (int i = 0; i < F_rows; i++) {
    // For each face, add face point
    face_points.row(i) = (V.row(F(i, 0)) + V.row(F(i, 1)) + V.row(F(i, 2)) + V.row(F(i, 3))) / 4.0;
  }

  std::unordered_map<int, std::vector<int>> new_edges;  // vertex index & indices of points sharing edge with vertex
  std::unordered_map<std::pair<int, int>, std::vector<int>, hash_pair> edge_to_faces;
  int current_vertex_index, another_vertex_index;
  // vertex index & edge points of the edges vertex-edges[vertex index]
  std::unordered_map<int, std::vector<Eigen::RowVector3d>> edge_points;
  std::pair<int, int> pair, reversed_pair;

  for (int i = 0; i < F_rows; i++) {
    for (int j = 0; j < 4; j++) {
      current_vertex_index = F(i, j);
      another_vertex_index = F(i, (j + 1) % 4);
      pair = std::make_pair(current_vertex_index, another_vertex_index);
      reversed_pair = std::make_pair(another_vertex_index, current_vertex_index);
      
      // if the edge has not been stored
      if (edge_to_faces.find(pair) == edge_to_faces.end() && edge_to_faces.find(reversed_pair) == edge_to_faces.end()) {
        std::vector<int> entry;
        edge_to_faces[pair] = entry;
        edge_to_faces[pair].push_back(i);
      } else if (edge_to_faces.find(pair) != edge_to_faces.end()) {
        edge_to_faces[pair].push_back(i);
      } else {
        edge_to_faces[reversed_pair].push_back(i);
      }
    }
  }

  int num_vertex;
  Eigen::RowVector3d edge_point_value(0, 0, 0);
  for (int i = 0; i < F_rows; i++) {
    for (int j = 0; j < 4; j++) {
      current_vertex_index = F(i, j);
      another_vertex_index = F(i, (j + 1) % 4);
      pair = std::make_pair(current_vertex_index, another_vertex_index);
      reversed_pair = std::make_pair(another_vertex_index, current_vertex_index);
    
      // if key not exist
      if (new_edges.find(current_vertex_index) == new_edges.end()) {
        std::vector<int> entry;
        new_edges[current_vertex_index] = entry;
        std::vector<Eigen::RowVector3d> edge_point_entry;
        edge_points[current_vertex_index] = edge_point_entry;
      }
      // if vertex_index-target_index_1 is not already stored
      if (std::find(new_edges[current_vertex_index].begin(), new_edges[current_vertex_index].end(), another_vertex_index) == new_edges[current_vertex_index].end() && (new_edges.find(another_vertex_index) == new_edges.end() || std::find(new_edges[another_vertex_index].begin(), new_edges[another_vertex_index].end(), current_vertex_index) == new_edges[another_vertex_index].end())) {
        new_edges[current_vertex_index].push_back(another_vertex_index);
        edge_point_value.setZero();
        if (edge_to_faces.find(pair) != edge_to_faces.end()) {
          num_vertex = edge_to_faces[pair].size();
          for (int m = 0; m < num_vertex; m++) {
            edge_point_value += face_points.row(edge_to_faces[pair].at(m));
          }
        } else {
          num_vertex = edge_to_faces[reversed_pair].size();
          for (int m = 0; m < num_vertex; m++) {
            edge_point_value += face_points.row(edge_to_faces[reversed_pair].at(m));
          }
        }
        num_vertex += 2;
        edge_points[current_vertex_index].push_back((V.row(another_vertex_index) + V.row(current_vertex_index) + edge_point_value) / (double) num_vertex);
        SV_size ++;
      }
    }
  }

  SV = Eigen::MatrixXd::Zero(SV_size, 3);
  SF = Eigen::MatrixXi::Zero(SF_size, 4);
  int num_edges;
  std::vector<std::vector<int> > VF;
  std::vector<int> face_indces;
  vertex_triangle_adjacency(F, V_rows, VF);

  Eigen::RowVector3d f(0, 0, 0), r(0, 0, 0);
  Eigen::RowVector3d barycenter(0, 0, 0);
  int num_points;
  for (int i = 0; i < V_rows; i++) {
    // calculate f, the average of all face points of faces adjacent to vertex at V.row(i)
    face_indces = VF[i];
    f.setZero();
    for (const int face_index : face_indces) {
      f += face_points.row(face_index);
    }
    f /= face_indces.size();
    
    num_edges = 0;
    // calculate r, the average of all edge points of edges touching vertex at V.row(i)
    r.setZero();
    for (auto entry : edge_to_faces ) {
      if (entry.first.first == i) {
        r += (V.row(entry.first.second) + V.row(i)) / 2.0;
        num_edges++;
      } else if (entry.first.second == i) {
        r += (V.row(entry.first.first) + V.row(i)) / 2.0;
        num_edges++;
      }
    }
    if (num_edges != 0) r /= num_edges;
    
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
      another_index_2 = F(i, (j - 1 + 4) % 4);
      count_1 = -1;
      count_2 = -1;
      index = 0;
      for (int a = 0; a < new_edges.size() && (count_1 < 0 || count_2 < 0); a++) {
        for (int b = 0; b < new_edges[a].size() && (count_1 < 0 || count_2 < 0); b++) {
          if ((a == point_index_1 && new_edges[a].at(b) == another_index_1) || (a == another_index_1 && new_edges[a].at(b) == point_index_1)) {
            count_1 = index;
          }

          if ((a == point_index_1 && new_edges[a].at(b) == another_index_2) || (a == another_index_2 && new_edges[a].at(b) == point_index_1)) {
            count_2 = index;
          }

          index++;
        }
      }

      point_index_2 = F_rows + V_rows + count_1;
      point_index_4 = F_rows + V_rows + count_2;
      
      SF.row(i * 4 + j) = Eigen::RowVector4i(point_index_1, point_index_2, point_index_3, point_index_4);
    }
  }

  catmull_clark(Eigen::MatrixXd(SV), Eigen::MatrixXi(SF), num_iters - 1, SV, SF);

  ////////////////////////////////////////////////////////////////////////////
}
