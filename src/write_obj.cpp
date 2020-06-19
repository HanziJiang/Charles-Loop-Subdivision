#include "write_obj.h"
#include <fstream>
#include <cassert>
#include <iostream>

bool write_obj(
  const std::string & filename,
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  const Eigen::MatrixXd & UV,
  const Eigen::MatrixXi & UF,
  const Eigen::MatrixXd & NV,
  const Eigen::MatrixXi & NF)
{
  assert((F.size() == 0 || F.cols() == 3 || F.cols() == 4) && "F must have 3 or 4 columns");
  ////////////////////////////////////////////////////////////////////////////
  
  // Open file and check if successful
  std::ofstream image_file;
  image_file.open(filename);
  if (image_file.fail() || !image_file.is_open()) return false;
  
  try {
    const int V_rows = V.rows();
    const int F_rows = F.rows();
    const int UV_rows = UV.rows();
    const int UF_rows = UF.rows();
    const int NV_rows = NV.rows();
    const int NF_rows = NF.rows();
    const int poly = F.cols();
    const int max_F = std::max(std::max(F_rows, UF_rows), NF_rows);

    int i;
    // Write geometric vertices (position information)
    for (i = 0; i < V_rows; i++) {
      image_file << "v " << (double) V(i, 0) << " " << (double) V(i, 1) << " " << (double) V(i, 2) << "\n" ;
    }

    // Write texture coordinates (parameterization information)
    for (i = 0; i < UV_rows; i++) {
      image_file << "vt " << (double) UV(i, 0) << " " << (double) UV(i, 1) << "\n" ;
    }

    // Write vertex normals (normal vector information)
    for (i = 0; i < NV_rows; i++) {
      image_file << "vn " << (double) NV(i, 0) << " " << (double) NV(i, 1) << " " << (double) NV(i, 2) << "\n" ;
    }

    // Write polygonal face elements
    for (i = 0; i < max_F; i++) {
      image_file << "f";
      for (int j = 0; j < poly; j++) {
        image_file << " ";
        if (i < F_rows) image_file << (int) F(i, j) + 1;
        image_file << "/";
        if (i < UF_rows) image_file << (int) UF(i, j) + 1;
        image_file << "/";
        if (i < NF_rows) image_file << (int) NF(i, j) + 1;
        image_file << "/";
      }
      image_file << "\n";
    }

    image_file.close();
    return true;
  } catch (const std::exception&) {
    if (image_file.is_open()) image_file.close();
  } catch (...) {
    if (image_file.is_open()) image_file.close();
  }
  ////////////////////////////////////////////////////////////////////////////
  return false;
}
