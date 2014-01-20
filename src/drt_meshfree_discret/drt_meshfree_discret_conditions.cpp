/*!---------------------------------------------------------------------------*\
 * \file drt_meshfree_discret.H
 *
 * \brief discretisation with additional reference point vector
 *        (meshfree analysis)
 *
 * <pre>
 * Maintainer: Keijo Nissen (nis)
 *             nissen@lnm.mw.tum.de
 *             http://www.lnm.mw.tum.de
 *             089 - 289-15253
 * </pre>
 *
\*----------------------------------------------------------------------------*/

#include "drt_meshfree_discret.H"

/*----------------------------------------------------------------------------*
 |  Build meshfree line geometry in a condition (protected)         nis Jan14 |
 *----------------------------------------------------------------------------*/
void DRT::MESHFREE::MeshfreeDiscretization::BuildLinesinCondition(
  const std::string            name,
  Teuchos::RCP<DRT::Condition> cond)
{
  //----------------------------------------------------------------------------
  // call base class BuildLinesinCondition
  //----------------------------------------------------------------------------
  DRT::Discretization::BuildLinesinCondition(name, cond);

  //----------------------------------------------------------------------------
  // assign nodes to new elements
  //----------------------------------------------------------------------------
  AssignNodesToCells(cond->Geometry(), *(cond->Nodes()));
  return;
}

/*----------------------------------------------------------------------------*
 |  Build meshfree surface geometry in a condition (protected)      nis Jan14 |
 *----------------------------------------------------------------------------*/
void DRT::MESHFREE::MeshfreeDiscretization::BuildSurfacesinCondition(
  const std::string            name,
  Teuchos::RCP<DRT::Condition> cond)
{
  //----------------------------------------------------------------------------
  // call base class BuildSurfacesinCondition
  //----------------------------------------------------------------------------
  DRT::Discretization::BuildSurfacesinCondition(name, cond);

  //----------------------------------------------------------------------------
  // assign nodes to new elements
  //----------------------------------------------------------------------------
  AssignNodesToCells(cond->Geometry(), *(cond->Nodes()));

  return;
}
