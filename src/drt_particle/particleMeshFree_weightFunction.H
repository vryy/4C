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

   //! compute the cubicBspline weight function
   static double Weight(const double &disRel, const double &radius);

   //! compute the weight function derivative of the cubicBspline weight function
   static double WeightDerivative(const double &disRel, const double &radius);

   //@}

private :

   static double RszDisRel(const double& disRel, const double& radius) {return 2 * disRel / radius; }

   static double Rsz3D(const double& radius) {return 8 * M_1_PI / std::pow(radius,3); }

};

class WeightFunction_SqrtHyperbola
{
public :

   //! @name Actions
   //@{

   //! compute the cubicBspline weight function
   static double Weight(const double &disRel, const double &radius);

   //! compute the weight function derivative of the cubicBspline weight function
   static double WeightDerivative(const double &disRel, const double &radius);

   //@}

private :

   static double Rsz3D(const double& radius) {return 15.0 / (4.0 * M_PI * std::pow(radius,3)); }

};

}  // namespace PARTICLE

/*----------------------------------------------------------------------*/
#endif /* SRC_DRT_PARTICLE_PARTICLEMESHFREE_WEIGHTFUNCTION_H_ */
