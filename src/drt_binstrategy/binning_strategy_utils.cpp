/*----------------------------------------------------------------------*/
/*! \file

\brief utils class for use of binning strategy

\level 2

\maintainer Jonas Eichinger
*----------------------------------------------------------------------*/


#include "binning_strategy_utils.H"

#include "../drt_scatra_ele/scatra_ele.H"
#include "../drt_fluid_ele/fluid_ele.H"
#include "../drt_bele3/bele3.H"
#include "../drt_beam3/beam3_base.H"
#include "../drt_rigidsphere/rigidsphere.H"
#include "../drt_so3/so_base.H"

#include "../drt_lib/drt_utils_parallel.H"

namespace BINSTRATEGY
{
  namespace UTILS
  {
    /*-----------------------------------------------------------------------------*
     *-----------------------------------------------------------------------------*/
    void ExtendDiscretizationGhosting(Teuchos::RCP<DRT::Discretization> discret,
        Teuchos::RCP<Epetra_Map> const& extendedelecolmap, bool assigndegreesoffreedom,
        bool initelements, bool doboundaryconditions)
    {
      // make sure that all procs are either filled or unfilled
      // oldmap in ExportColumnElements must be Reset() on every proc or nowhere
      discret->CheckFilledGlobally();

      // adapt layout to extended ghosting in discret
      // first export the elements according to the processor local element column maps
      discret->ExportColumnElements(*extendedelecolmap);

      // get the node ids of the elements that are to be ghosted
      // and create a proper node column map for their export
      std::set<int> nodes;
      for (int lid = 0; lid < extendedelecolmap->NumMyElements(); ++lid)
      {
        DRT::Element* ele = discret->gElement(extendedelecolmap->GID(lid));
        const int* nodeids = ele->NodeIds();
        for (int inode = 0; inode < ele->NumNode(); ++inode) nodes.insert(nodeids[inode]);
      }

      std::vector<int> colnodes(nodes.begin(), nodes.end());
      Teuchos::RCP<Epetra_Map> nodecolmap =
          Teuchos::rcp(new Epetra_Map(-1, (int)colnodes.size(), &colnodes[0], 0, discret->Comm()));

      // now ghost the nodes
      discret->ExportColumnNodes(*nodecolmap);

      // fillcomplete discret with extended ghosting
      discret->FillComplete(assigndegreesoffreedom, initelements, doboundaryconditions);

#ifdef DEBUG
      // print distribution after standard ghosting
      if (discret->Comm().MyPID() == 0)
        std::cout << "parallel distribution with extended ghosting" << std::endl;
      DRT::UTILS::PrintParallelDistribution(*discret);
#endif

      return;
    }

    /*-----------------------------------------------------------------------------*
     *-----------------------------------------------------------------------------*/
    BINSTRATEGY::UTILS::BinContentType ConvertElementToBinContentType(
        DRT::Element const* const eleptr)
    {
      // (Todo make this nicer and cheaper)

      if (dynamic_cast<DRT::ELEMENTS::Transport const*>(eleptr) != NULL)
      {
        return BINSTRATEGY::UTILS::Scatra;
      }
      else if (dynamic_cast<DRT::ELEMENTS::Fluid const*>(eleptr) != NULL)
      {
        return BINSTRATEGY::UTILS::Fluid;
      }
      else if (dynamic_cast<DRT::ELEMENTS::Bele3 const*>(eleptr) != NULL)
      {
        return BINSTRATEGY::UTILS::BELE3;
      }
      else if (dynamic_cast<DRT::ELEMENTS::Beam3Base const*>(eleptr) != NULL)
      {
        return BINSTRATEGY::UTILS::Beam;
      }
      else if (dynamic_cast<DRT::ELEMENTS::Rigidsphere const*>(eleptr) != NULL)
      {
        return BINSTRATEGY::UTILS::RigidSphere;
      }
      else if (dynamic_cast<DRT::ELEMENTS::So_base const*>(eleptr) != NULL)
      {
        return BINSTRATEGY::UTILS::Solid;
      }
      else
      {
        dserror(
            " Element you are about to assign to a bin could not be converted"
            " to a valid bin content type. ");
        exit(EXIT_FAILURE);
      }
    }



  }  // namespace UTILS
}  // namespace BINSTRATEGY
