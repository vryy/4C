/*!-----------------------------------------------------------------------------------------------*
\file xfem_edgestab.H

\brief provides the xfem fluid and ghost penalty stabilization based on EOS/CIP (edge-oriented, continuous interior penalty) scheme

<pre>
Maintainer: Benedikt Schott
            schott@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15241
</pre>
*------------------------------------------------------------------------------------------------*/


#ifndef XFEM_EDGESTAB_H_
#define XFEM_EDGESTAB_H_


#include <vector>
#include <Teuchos_RCP.hpp>


namespace DRT
{
  class AssembleStrategy;
  class Discretization;
  class Element;

  namespace ELEMENTS
  {
    class Fluid3;
  }
}

namespace GEO
{
  namespace CUT
  {
    class Side;
  }
}


namespace XFEM
{
  class FluidWizard;

/*!
\brief this class is the basic TIMEINT class for the projection, adaption or
       something else in XFEM-problems between consecutive time steps
 */
class XFEM_EdgeStab
{
public:

  //! constructor
  explicit XFEM_EdgeStab(Teuchos::RCP<XFEM::FluidWizard>              wizard);

  //! destructor
  ~XFEM_EdgeStab()
  {
    return;
  }

  void findSurface(std::vector<int>                nodeids,
                   DRT::ELEMENTS::Fluid3 *         actele,
                   std::vector<Teuchos::RCP<DRT::Element> > surfaces,
                   int & surfid);


  void findSide   ( std::vector<int>                nodeids, // node ids of surface
                    DRT::ELEMENTS::Fluid3 *         actele,
                    const std::vector<GEO::CUT::Side*> & sides,
                    int & side_index);

  void findNeighborElement(DRT::ELEMENTS::Fluid3 * actele,
                           Teuchos::RCP<DRT::Element> surface,
                           bool & neighbor_found,
                           int & neighbor_id);


  void stabilizeStandardElement     ( DRT::AssembleStrategy& strategy,
                                      DRT::ELEMENTS::Fluid3 *actele,
                                      DRT::Discretization & discret);

  void EvaluateEdgeStabandGhostPenalty( DRT::Discretization &  discret,
                                       DRT::AssembleStrategy&  strategy,
                                       DRT::ELEMENTS::Fluid3 * ele);

  void callEdgeStabandGhostPenalty(                           bool & edge_based_stab,
                                                              bool & ghost_penalty,
                                                              Teuchos::RCP<DRT::Element>       surface,
                                                              DRT::Element * fele_1,
                                                              DRT::Element * fele_2,
                                                              std::vector<int> & nds_1,
                                                              std::vector<int> & nds_2,
                                                              DRT::Discretization &   discret,
                                                              DRT::AssembleStrategy&  strategy);

  private:

  Teuchos::RCP<XFEM::FluidWizard> wizard_;

};

} // namespace XFEM

#endif /*XFEM_EDGESTAB_H_*/
