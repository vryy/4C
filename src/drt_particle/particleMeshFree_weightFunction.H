/*----------------------------------------------------------------------*/
/*!
\file particleMeshFree_weightFunction.H

\brief weight functions for MeshFree methods

\level 3

\maintainer Alessandro Cattabiani
*/

/*----------------------------------------------------------------------*/
#ifndef SRC_DRT_PARTICLE_PARTICLEMESHFREE_WEIGHTFUNCTION_H_
#define SRC_DRT_PARTICLE_PARTICLEMESHFREE_WEIGHTFUNCTION_H_

/*----------------------------------------------------------------------*/
/* headers */
#include "../linalg/linalg_fixedsizematrix.H"

/*----------------------------------------------------------------------*/
/* belongs to particle namespace */
namespace PARTICLE
{

  class WeightFunction_CubicBspline
  {
  public :

     //! @name Actions
     //@{

     //!
     static double Weight(const double &distRel, const double &radius);

     //!
     static double DerivativeWeight(const double &distRel, const double &radius);

     //!
     static LINALG::Matrix<3,1> GradientWeight(LINALG::Matrix<3,1> &WFGrad, const double &radius);

     //@}

  private :

     static double Resizer3D() {return 1.5 * M_1_PI; }

  };

}  // namespace PARTICLE

/*----------------------------------------------------------------------*/
#endif /* SRC_DRT_PARTICLE_PARTICLEMESHFREE_WEIGHTFUNCTION_H_ */
