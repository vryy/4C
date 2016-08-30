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
/* belongs to particle namespace */
namespace PARTICLE
{

  class WeightFunctionBase
  {
  public :

    //! @name Life
     //@{

     //! Constructor
     WeightFunctionBase(){ ; }

     //! Destructor
     virtual ~WeightFunctionBase(){ ; }

     //@}

     //! @name Actions
     //@{

     //!
     virtual double ComputeWeight(const double dist_rel, const double radius) = 0;

     //!
     virtual double ComputeGradientWeight(const double dist_rel, const double radius) = 0;

     //@}

  };

  class WeightFunction_CubicBspline : public PARTICLE::WeightFunctionBase
  {
  public :
    //! @name Life
     //@{

     //! Constructor
     WeightFunction_CubicBspline(){ ; }

     //! Destructor
     virtual ~WeightFunction_CubicBspline(){ ; }

     //@}

     //! @name Actions
     //@{

     //!
     virtual double ComputeWeight(const double dist_rel, const double radius);

     //!
     virtual double ComputeGradientWeight(const double dist_rel, const double radius);

     //@}

  };

}  // namespace PARTICLE

/*----------------------------------------------------------------------*/
#endif /* SRC_DRT_PARTICLE_PARTICLEMESHFREE_WEIGHTFUNCTION_H_ */
