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
#include "../drt_inpar/inpar_particle.H"

//TODO: remove this include as soon as the define flags in particle_algorithm.H have been set as input parameters!
#include "particle_algorithm.H"

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

   //Pre-factor of smoothing kernel depends on space dimension (see Monaghan2005)
   static double Rsz3D(const double& radius)
   {
     if(PARTICLE_DIM==3)
       return 8.0 * M_1_PI / std::pow(radius,3);
     else if(PARTICLE_DIM==2)
       return 40.0 * M_1_PI / (7.0*std::pow(radius,2));
     else if(PARTICLE_DIM==1)
       return 4.0 / (3.0*radius);
     else
       dserror("Only the problem dimensions 1, 2 and 3 are possible!");
   }

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


class WeightFunction_HyperbolaNoRsz
{
public :

   //! @name Actions
   //@{

   //! compute the cubicBspline weight function
   static double Weight(const double &disRel, const double &radius);

   //! compute the weight function derivative of the cubicBspline weight function
   static double WeightDerivative(const double &disRel, const double &radius);

   //@}

};


}  // namespace PARTICLE

/*----------------------------------------------------------------------*/
#endif /* SRC_DRT_PARTICLE_PARTICLEMESHFREE_WEIGHTFUNCTION_H_ */
