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
 |  Overloading base class function to prevent storing condition node
 |  sets twice (protected)         nis Jan14 |
 *----------------------------------------------------------------------------*/
void DRT::MESHFREE::MeshfreeDiscretization::SetCondition(
  const std::string& name,
  Teuchos::RCP<Condition> cond
  )
{
  // get face of correct dimension
  std::vector<Face> faceset = domaintopology_[cond->GType()];

  // search to which face condition is enforced
  for(unsigned i=0; i<faceset.size(); ++i)
  {
    if (*(faceset[i].Nodes())== *(cond->Nodes()))
    {
      // let condition data point to same set of node ids as in face
      if (faceset[i].Nodes().get()==cond->Nodes())
        dserror("This move is obsolete since data is already the same.");
      else
      {
        // bend face nodeid pointer to condition nodeids to store node ids only once
        cond->Add("Node Ids",faceset[i].Nodes());
        // call base class method to finally add condition
        DRT::Discretization::SetCondition(name, cond);
      }

      // where done here
      return;
    }
  }

  // if condition could not be linked to a face
  dserror("No face could be identified on which condition in enforced.");

  return;
}

/*----------------------------------------------------------------------------*
 |  Build meshfree line geometry in a condition (protected)         nis Jan14 |
 *----------------------------------------------------------------------------*/
bool DRT::MESHFREE::MeshfreeDiscretization::BuildLinesinCondition(
  const std::string            name,
  Teuchos::RCP<DRT::Condition> cond)
{
  // get face of correct dimension - here 1 for lines
  std::vector<Face> faceset = domaintopology_[1];

  // search to which face the condition is enforced
  for(unsigned i=0; i<faceset.size(); ++i)
  {
    if (*(faceset[i].Nodes())== *(cond->Nodes()))
    {
      bool havenewelements = false;

      // if no geometry exists on this face...
      if (faceset[i].Geometry()==Teuchos::null)
      {
        // ... create new geometry for condition via base class call, ...
        havenewelements = DRT::Discretization::BuildLinesinCondition(name, cond);
        // ... make this geometry known to this face,
        faceset[i].AddGeometry(Teuchos::rcpFromRef(cond->Geometry()));
        // ... and assign nodes to cells
        faceset[i].AssignNodesToCells();
      }
      // ... else, the face already has a geometry
      else
      {
        // ... and we don't need to build a new geometry ...
        havenewelements = false;
        // ... but make this geometry known to current condition
        cond->AddGeometry(faceset[i].Geometry());
      }

      // where done here
      return havenewelements;
    }
  }

  // if condition could not be linked to a face
  dserror("No face could be identified on which condition in enforced.");

  return false;
}

/*----------------------------------------------------------------------------*
 |  Build meshfree surface geometry in a condition (protected)      nis Jan14 |
 *----------------------------------------------------------------------------*/
bool DRT::MESHFREE::MeshfreeDiscretization::BuildSurfacesinCondition(
  const std::string            name,
  Teuchos::RCP<DRT::Condition> cond)
{
  // get face of correct dimension - here 2 for surfaces
  std::vector<Face> faceset = domaintopology_[2];

  // search to which face the condition is enforced
  for(unsigned i=0; i<faceset.size(); ++i)
  {
    if (*(faceset[i].Nodes())== *(cond->Nodes()))
    {
      bool havenewelements = false;

      // if no geometry set in face...
      if (faceset[i].Geometry()==Teuchos::null)
      {
        // ... create new geometry for condition via base class call ...
        havenewelements = DRT::Discretization::BuildSurfacesinCondition(name, cond);
        // ... and make this geometry known to this face
        faceset[i].AddGeometry(Teuchos::rcpFromRef(cond->Geometry()));
        // ... and assign nodes to cells
        faceset[i].AssignNodesToCells();
      }
      // ... else the face already has a geometry
      else
      {
        // ... make this geometry known to current condition
        cond->AddGeometry(faceset[i].Geometry());
      }

      // where done here
      return havenewelements;
    }
  }

  // if condition could not be linked to a face
  dserror("No face could be identified on which condition in enforced.");

  return false;
}
