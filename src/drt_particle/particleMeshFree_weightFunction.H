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

class WeightFunction_Base
{
public :

  // destructor
  virtual ~WeightFunction_Base() { ; };

  // weight function
  virtual double Weight(const double &disRel, const double &radius) = 0;

  // weight function derivative
  virtual double WeightDerivative(const double &disRel, const double &radius) = 0;

};

class WeightFunction_CubicBspline : public WeightFunction_Base
{
public :

   //! @name Actions
   //@{

   //! compute the cubicBspline weight function
   virtual double Weight(const double &disRel, const double &radius);

   //! compute the weight function derivative of the cubicBspline weight function
   virtual double WeightDerivative(const double &disRel, const double &radius);

   //@}

private :

   virtual double RszDisRel(const double& disRel, const double& radius) {return 2 * disRel / radius; }

   virtual double Rsz3D(const double& radius) {return 8 * M_1_PI / std::pow(radius,3); }

};

class WeightFunction_SqrtHyperbola : public WeightFunction_Base
{
public :

   //! @name Actions
   //@{

   //! compute the cubicBspline weight function
   virtual double Weight(const double &disRel, const double &radius);

   //! compute the weight function derivative of the cubicBspline weight function
   virtual double WeightDerivative(const double &disRel, const double &radius);

   //@}

private :

   virtual double Rsz3D(const double& radius) {return 15.0 / (4.0 * M_PI * std::pow(radius,3)); }

};


class WeightFunction_HyperbolaNoRsz : public WeightFunction_Base
{
public :

   //! @name Actions
   //@{

   //! compute the cubicBspline weight function
   virtual double Weight(const double &disRel, const double &radius);

   //! compute the weight function derivative of the cubicBspline weight function
   virtual double WeightDerivative(const double &disRel, const double &radius);

   //@}

};


}  // namespace PARTICLE

/*----------------------------------------------------------------------*/
#endif /* SRC_DRT_PARTICLE_PARTICLEMESHFREE_WEIGHTFUNCTION_H_ */
