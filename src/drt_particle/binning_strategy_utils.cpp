/*----------------------------------------------------------------------*/
/*!
\file binning_strategy_utils.cpp

\brief utils class for use of binning strategy

\level 2

\maintainer Georg Hammerl
*----------------------------------------------------------------------*/


#include "../drt_scatra_ele/scatra_ele.H"
#include "../drt_fluid_ele/fluid_ele.H"
#include "../drt_bele3/bele3.H"
#include "../drt_beam3/beam3_base.H"
#include "../drt_rigidsphere/rigidsphere.H"


#include "binning_strategy_utils.H"

namespace BINSTRATEGY
{
namespace UTILS
{

/*-----------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------*/
BINSTRATEGY::UTILS::BinContentType ConvertElementToBinContentType(
    DRT::Element const * const eleptr )
{
  if( dynamic_cast<const DRT::ELEMENTS::Transport*>(eleptr) != NULL )
  {
    return BINSTRATEGY::UTILS::Scatra;
  }
  else if( dynamic_cast<const DRT::ELEMENTS::Fluid*>(eleptr) != NULL )
  {
    return BINSTRATEGY::UTILS::Fluid;
  }
  else if( dynamic_cast<const DRT::ELEMENTS::Bele3*>(eleptr) != NULL )
  {
    return BINSTRATEGY::UTILS::BELE3;
  }
  else if( dynamic_cast<const DRT::ELEMENTS::Beam3Base*>(eleptr) != NULL )
  {
    return BINSTRATEGY::UTILS::Beam;
  }
  else if ( dynamic_cast<const DRT::ELEMENTS::Rigidsphere*>(eleptr) != NULL)
  {
    return BINSTRATEGY::UTILS::RigidSphere;
  }
  else
  {
    dserror(" Element you are about to assign to a bin could not be converted"
            " to a valid bin content type. ");
    // to make compiler happy
    return BINSTRATEGY::UTILS::enumsize;
  }
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void GetCurrentNodePos(
    Teuchos::RCP<DRT::Discretization> discret,
    const DRT::Node* node,
    Teuchos::RCP<Epetra_Vector> disnp,
    double* currpos )
{
  // Todo make this nicer
  // the problem is that we might have nodes without position DoFs
  // (e.g. for beam elements with 'interior' nodes that are only used for
  // triad interpolation)
  // instead of the node position itself, we return the position of the
  // first node of the  element here (for the sake of binning)

  // standard case
  DRT::Node const* node_with_position_Dofs = node;

  const DRT::Element* element = node->Elements()[0];
  const DRT::ELEMENTS::Beam3Base* beamelement =
      dynamic_cast<const DRT::ELEMENTS::Beam3Base*>(element);

  // if the node does not have position DoFs, we return the position of the first
  // node of the corresponding element
  if (beamelement != NULL and not beamelement->IsCenterlineNode(*node) )
  {
    node_with_position_Dofs = beamelement->Nodes()[0];
  }

  if(disnp!=Teuchos::null)
  {
    const int gid = discret->Dof(node_with_position_Dofs, 0);
    const int lid = disnp->Map().LID(gid);
    if( lid < 0 )
      dserror("Your displacement is incomplete (need to be based on a column map"
              " as this function is also called from a loop over elements and "
              "each proc does (usually) not own all nodes of his row elements ");
    for(int dim=0; dim<3; ++dim)
    {
      currpos[dim] = node_with_position_Dofs->X()[dim] + (*disnp)[lid+dim];
    }
  }
  else
  {
    for(int dim=0; dim<3; ++dim)
      currpos[dim] = node_with_position_Dofs->X()[dim];
  }

  return;
}

}
}
