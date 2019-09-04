/*----------------------------------------------------------------------*/
/*! \file

\brief Run time for locsys
\level 2
\maintainer Martin Kronbichler

*----------------------------------------------------------------------*/

#include <string>
#include <fstream>

#include <Teuchos_CommandLineProcessor.hpp>
#include <Teuchos_RCP.hpp>

#include <cstdlib>

#ifdef PARALLEL
#include <Epetra_MpiComm.h>
#else
#include <Epetra_SerialComm.h>
#endif

#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_inputreader.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_lib/drt_parobjectregister.H"
#include "../linalg/linalg_fixedsizematrix.H"
#include "../drt_fem_general/largerotations.H"

/*======================================================================*/
/*======================================================================*/
int main(int argc, char** argv)
{
  // base vectors
  LINALG::Matrix<3, 1> vector1;
  LINALG::Matrix<3, 1> vector2;
  LINALG::Matrix<3, 1> vector3;
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
  LINALG::Matrix<3, 3> S_vector1;
  LARGEROTATIONS::computespin(S_vector1, vector1);
  vector3.Multiply(S_vector1, vector2);

  // Compute rotation matrix
  LINALG::Matrix<3, 3> rotmatrix;
  for (int i = 0; i < 3; i++)
  {
    rotmatrix(i, 0) = vector1(i);
    rotmatrix(i, 1) = vector2(i);
    rotmatrix(i, 2) = vector3(i);
  }

  // Compute rotation angle via quaterion
  LINALG::Matrix<4, 1> quaterion;
  LARGEROTATIONS::triadtoquaternion(rotmatrix, quaterion);
  LINALG::Matrix<3, 1> rotangle;
  LARGEROTATIONS::quaterniontoangle(quaterion, rotangle);

  std::cout << std::endl << std::setprecision(10) << "Rotation vector: " << rotangle << std::endl;

  //  //Check via inverse mapping
  //  LINALG::Matrix<3,3> rotmatrix_test;
  //  LARGEROTATIONS::angletotriad(rotangle, rotmatrix_test);
  //
  //  std::cout << endl << std::setprecision(10)<<  "rotmatrix_test: " << rotmatrix_test <<  endl;

  return 0;
}
