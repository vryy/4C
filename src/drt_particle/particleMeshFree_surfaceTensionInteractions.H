/*----------------------------------------------------------------------*/
/*!
\file particleMeshFree_surfaceTensionInteractions.H

\brief surface tension - related functions

\level 3

\maintainer  Christoph Meier
             meier@lnm.mw.tum.de
             http://www.lnm.mw.tum.de

// theory: http://doi.acm.org/10.1145/2508363.2508395\nhttp://dl.acm.org/ft_gateway.cfm?id=2508395&type=pdf R1

*-----------------------------------------------------------------------*/
#ifndef SRC_DRT_PARTICLE_PARTICLEMESHFREE_SURFACETENSIONINTERACTIONS_H_
#define SRC_DRT_PARTICLE_PARTICLEMESHFREE_SURFACETENSIONINTERACTIONS_H_

/*----------------------------------------------------------------------*/
/* headers */
#include "../linalg/linalg_fixedsizematrix.H"

/*----------------------------------------------------------------------*/
/* belongs to particle namespace */
namespace PARTICLE
{

class SurfaceTensionInteractions
{
public :

   //! @name Actions
   //@{

   //! cohesion interaction (E2)
   static double Cohesion(const double &disRel, const double &radius);

   //! adhesion interaction (E7)
   static double Adhesion(const double &disRel, const double &radius);

   //@}

};


}  // namespace PARTICLE

/*----------------------------------------------------------------------*/
#endif /* SRC_DRT_PARTICLE_PARTICLEMESHFREE_SURFACETENSIONINTERACTIONS_H_ */
