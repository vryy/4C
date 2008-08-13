/*!----------------------------------------------------------------------
\file multipointconstraint.cpp

\brief Basic constraint class, dealing with multi point constraints
<pre>
Maintainer: Thomas Kloeppel
            kloeppel@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/Members/kloeppel
            089 - 289-15257
</pre>

*----------------------------------------------------------------------*/
#ifdef CCADISCRET


#include "multipointconstraint.H"
#include "mpcdofset.H"
#include "constraint_element.H"

#include "../drt_lib/drt_discret.H"
#include "../drt_lib/linalg_utils.H"
#include "../drt_lib/linalg_sparsematrix.H"
#include "iostream"
#include "../drt_lib/drt_condition_utils.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_lib/drt_timecurve.H"



/*----------------------------------------------------------------------*
 |  ctor (public)                                               tk 07/08|
 *----------------------------------------------------------------------*/
UTILS::MPConstraint::MPConstraint(RCP<DRT::Discretization> discr,
        const string& conditionname,
        int& minID,
        int& maxID)
: UTILS::Constraint
  (
    discr,
    conditionname,
    minID,
    maxID
  )
{  
  return;
}

/*----------------------------------------------------------------------*
 |  ctor (public)                                               tk 07/08|
 *----------------------------------------------------------------------*/
UTILS::MPConstraint::MPConstraint(RCP<DRT::Discretization> discr,
        const string& conditionname)
: UTILS::Constraint
  (
    discr,
    conditionname
  )
{  
  return;
}


/*----------------------------------------------------------------------*
 |(private)                                                   tk 05/08  |
 |replace numdofs in elements of constraint discretization              |
 *----------------------------------------------------------------------*/
void UTILS::MPConstraint::ReplaceNumDof
(
  const RCP<DRT::Discretization> sourcedis,
  const RCP<DRT::Discretization> constraintdis
) const
{
  // find typical numdof of basis discretization (may not work for XFEM)
  const DRT::Element* actele = sourcedis->lColElement(0);
  const DRT::Node*const* nodes = actele->Nodes();
  const int mpc_numdof = sourcedis->NumDof(nodes[0]);

  // change numdof for all constraint elements
  const int numcolele = constraintdis->NumMyColElements();
  for (int i=0; i<numcolele; ++i)
  {
    DRT::ELEMENTS::ConstraintElement* mpcele = dynamic_cast<DRT::ELEMENTS::ConstraintElement*>(constraintdis->lColElement(i));
    mpcele->SetNumDofPerNode(mpc_numdof);
  }
  constraintdis->FillComplete();
  return;
 }


/*----------------------------------------------------------------------*
 |(private)                                                   tk 04/08  |
 |recompute nodecolmap of standard discretization to include constrained|
 |nodes as ghosted nodes                                                |
 *----------------------------------------------------------------------*/
RCP<Epetra_Map> UTILS::MPConstraint::ComputeNodeColMap(
        const RCP<DRT::Discretization> sourcedis,
        const RCP<DRT::Discretization> constraintdis
        ) const
{
    const Epetra_Map* oldcolnodemap = sourcedis->NodeColMap();

    vector<int> mycolnodes(oldcolnodemap->NumMyElements());
    oldcolnodemap->MyGlobalElements (&mycolnodes[0]);
    for (int inode = 0; inode != constraintdis->NumMyColNodes(); ++inode)
    {
        const DRT::Node* newnode = constraintdis->lColNode(inode);
        const int gid = newnode->Id();
        if (!(sourcedis->HaveGlobalNode(gid)))
        {
            mycolnodes.push_back(gid);
        }
    }

    // now reconstruct the extended colmap
    RCP<Epetra_Map> newcolnodemap = rcp(new Epetra_Map(-1,
                                       mycolnodes.size(),
                                       &mycolnodes[0],
                                       0,
                                       sourcedis->Comm()));
    return newcolnodemap;
}

#endif
