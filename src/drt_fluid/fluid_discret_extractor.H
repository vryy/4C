/*!-----------------------------------------------------------------------------------------------*
\file fluid_turb_inflow_gen.H

  \brief creates a second discretization as part of the complete discretization for inflow generation

<pre>
Maintainer: Ursula Rasthofer
            rasthofer@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15236
</pre>
 *------------------------------------------------------------------------------------------------*/

#ifdef CCADISCRET

#ifndef FLUID_DISCRET_EXTRACTOR_H
#define FLUID_DISCRET_EXTRACTOR_H

#include "../drt_lib/drt_discret.H"

namespace FLD {

class FluidDiscretExtractor
{
public:
   /*!
  \brief Constructor

  */
  FluidDiscretExtractor(
    RCP<DRT::Discretization> actdis,         //! parent discretization
    const string&            condition,      //! condition for separation of domain
    bool                     yescondition);  //! (unused) bool to distinguish between all nodes having the condition and all nodes not having it

   /*!
  \brief Destructor

  */
  virtual ~FluidDiscretExtractor();

  //! get child discretization
  const RCP<DRT::Discretization> GetChildDiscretization(){return childdiscret_;}
  //! get node to node coupling in case of periodic boundary conditions
  RCP<map<int,vector<int> > > GetCoupledColNodesChildDiscretization(){return col_pbcmapmastertoslave_;}
private:
  //! the parent discretization
  RCP<DRT::Discretization>    parentdiscret_;
  //! the child discretization
  RCP<DRT::Discretization>    childdiscret_;
  //! periodic boundary condition: node to node coupling
  RCP<map<int,vector<int> > > col_pbcmapmastertoslave_;
};

}

#endif

#endif  // #ifdef CCADISCRET
