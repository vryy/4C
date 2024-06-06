/*----------------------------------------------------------------------*/
/*! \file

\brief Run time for locsys
\level 2

*----------------------------------------------------------------------*/

#include "4C_discretization_fem_general_largerotations.hpp"
#include "4C_linalg_fixedsizematrix.hpp"

#include <Teuchos_RCP.hpp>

#include <fstream>

/*======================================================================*/
/*======================================================================*/
int main(int argc, char** argv)
{
  using namespace FourC;

  // base vectors
  Core::LinAlg::Matrix<3, 1> vector1;
  Core::LinAlg::Matrix<3, 1> vector2;
  Core::LinAlg::Matrix<3, 1> vector3;
  vector1.Clear();
  vector2.Clear();
  vector3.Clear();

  for (int i = 0; i < 3; i++)
  {
    std::cout << "Enter component " << i + 1 << " of base vector 1: ";
    std::cin >> vector1(i);
  }

  for (int i = 0; i < 3; i++)
  {
    std::cout << "Enter component " << i + 1 << " of base vector 2: ";
    std::cin >> vector2(i);
  }

  // Check, if vectors are perpendicular to each other
  double scalarproduct = 0.0;
  for (int i = 0; i < 3; i++)
  {
    scalarproduct += vector1(i) * vector2(i);
  }
  if (scalarproduct > 1.0e-10)
  {
    std::cout << "The two base vectors are not perpendicular!" << std::endl;
    return (1);
  }

  // Scale vectors to unity
  double normvec1 = vector1.Norm2();
  vector1.Scale(1.0 / normvec1);
  double normvec2 = vector2.Norm2();
  vector2.Scale(1.0 / normvec2);

  // Compute third base vector
  Core::LinAlg::Matrix<3, 3> S_vector1;
  Core::LargeRotations::computespin(S_vector1, vector1);
  vector3.Multiply(S_vector1, vector2);

  // Compute rotation matrix
  Core::LinAlg::Matrix<3, 3> rotmatrix;
  for (int i = 0; i < 3; i++)
  {
    rotmatrix(i, 0) = vector1(i);
    rotmatrix(i, 1) = vector2(i);
    rotmatrix(i, 2) = vector3(i);
  }

  // Compute rotation angle via quaterion
  Core::LinAlg::Matrix<4, 1> quaterion;
  Core::LargeRotations::triadtoquaternion(rotmatrix, quaterion);
  Core::LinAlg::Matrix<3, 1> rotangle;
  Core::LargeRotations::quaterniontoangle(quaterion, rotangle);

  std::cout << std::endl << std::setprecision(10) << "Rotation vector: " << rotangle << std::endl;

  //  //Check via inverse mapping
  //  Core::LinAlg::Matrix<3,3> rotmatrix_test;
  //  Core::LargeRotations::angletotriad(rotangle, rotmatrix_test);
  //
  //  std::cout << endl << std::setprecision(10)<<  "rotmatrix_test: " << rotmatrix_test <<  endl;

  return 0;
}
