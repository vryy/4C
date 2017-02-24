/*----------------------------------------------------------------------------*/
/**
\file drt_utils_discret.cpp

\brief Utils methods concerning the discretization

\maintainer Martin Kronbichler

\level 1

*/
/*----------------------------------------------------------------------------*/
#include "drt_utils_discret.H"
#include "drt_discret_interface.H"
#include "drt_node.H"
#include "drt_globalproblem.H"

#include <Epetra_Map.h>


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void DRT::UTILS::EvaluateInitialField(
    const DRT::DiscretizationInterface & discret,
    const std::string&                  fieldstring,
    Teuchos::RCP<Epetra_Vector>         fieldvector,
    const std::vector<int>              locids)
{
  // check for valid input
  bool invalid = false;
  if (fieldstring=="Velocity" && (int)locids.size()!=3) invalid = true;
  if (fieldstring=="Pressure" && (int)locids.size()!=1) invalid = true;
  if (fieldstring=="Temperature" && (int)locids.size()!=1) invalid = true;
  if (fieldstring=="Porosity" && (int)locids.size()!=1) invalid = true;
  if (invalid) dserror("ERROR: Invalid input to EvaluateInitialField().");

  // get initial field conditions
  std::vector<DRT::Condition*> initfieldconditions(0);
  discret.GetCondition("Initfield",initfieldconditions);

  //--------------------------------------------------------
  // loop through Initfield conditions and evaluate them
  //--------------------------------------------------------
  // Note that this method does not sum up but 'sets' values in fieldvector.
  // For this reason, Initfield BCs are evaluated hierarchical meaning
  // in this order (just like Dirichlet BCs):
  //                VolumeInitfield
  //                SurfaceInitfield
  //                LineInitfield
  //                PointInitfield
  // This way, lower entities override higher ones. Whether
  // this is really useful for Initfield BCs, I don't know... (popp 06/11)

  // Do VolumeInitfield first
  for (int i=0; i<(int)initfieldconditions.size(); ++i)
  {

    if (initfieldconditions[i]->Type() != DRT::Condition::VolumeInitfield) continue;
    const std::string* condstring = initfieldconditions[i]->Get<std::string>("Field");
    if (*condstring != fieldstring) continue;
    DoInitialField(discret,*initfieldconditions[i],fieldvector,locids);
  }

  // Do SurfaceInitfield
  for (int i=0; i<(int)initfieldconditions.size(); ++i)
  {

    if (initfieldconditions[i]->Type() != DRT::Condition::SurfaceInitfield) continue;
    const std::string* condstring = initfieldconditions[i]->Get<std::string>("Field");
    if (*condstring != fieldstring) continue;
    DoInitialField(discret,*initfieldconditions[i],fieldvector,locids);
  }

  // Do LineInitfield
  for (int i=0; i<(int)initfieldconditions.size(); ++i)
  {

    if (initfieldconditions[i]->Type() != DRT::Condition::LineInitfield) continue;
    const std::string* condstring = initfieldconditions[i]->Get<std::string>("Field");
    if (*condstring != fieldstring) continue;
    DoInitialField(discret,*initfieldconditions[i],fieldvector,locids);
  }

  // Do PointInitfield
  for (int i=0; i<(int)initfieldconditions.size(); ++i)
  {

    if (initfieldconditions[i]->Type() != DRT::Condition::PointInitfield) continue;
    const std::string* condstring = initfieldconditions[i]->Get<std::string>("Field");
    if (*condstring != fieldstring) continue;
    DoInitialField(discret,*initfieldconditions[i],fieldvector,locids);
  }

  return;
}


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void DRT::UTILS::DoInitialField(
     const DRT::DiscretizationInterface & discret,
     DRT::Condition&                      cond,
     Teuchos::RCP<Epetra_Vector>          fieldvector,
     const std::vector<int>               locids)
{
  const std::vector<int>* nodeids = cond.Nodes();
  if (!nodeids) dserror("Initfield condition does not have nodal cloud");
  const int nnode = (*nodeids).size();
  const std::vector<int>* funct  = cond.Get<std::vector<int> >("funct");

  // check fieldvector
  if (fieldvector==Teuchos::null)
    dserror("ERROR: Fieldvector must not be Teuchos::null");

  // loop nodes to identify and evaluate spatial distributions
  // of Initfield boundary conditions
  for (int i=0; i<nnode; ++i)
  {
    // do only nodes in my row map
    int nlid = discret.NodeRowMap()->LID((*nodeids)[i]);
    if (nlid < 0) continue;
    DRT::Node* actnode = discret.lRowNode(nlid);

    // call explicitly the main dofset, i.e. the first column
    std::vector<int> dofs = discret.Dof(0,actnode);
    const unsigned total_numdf = dofs.size();

    // Get native number of dofs at this node. There might be multiple dofsets
    // (in xfem cases), thus the size of the dofs vector might be a multiple
    // of this.
    const int numele = actnode->NumElement();
    const DRT::Element * const * myele = actnode->Elements();
    int numdf = 0;
    for (int j=0; j<numele; ++j)
      numdf = std::max(numdf,myele[j]->NumDofPerNode(*actnode));

    if ( ( total_numdf % numdf ) != 0 )
      dserror( "illegal dof set number" );

    // now loop over all relevant DOFs
    for (unsigned j=0; j<total_numdf; ++j)
    {
      // check if something needs to be done for this DOF
      bool dosomething = false;
      int localdof = j % numdf;

      // something needs to be done if local DOF id exists
      // in the given locids vector
      for (int k=0;k<(int)(locids.size());++k)
        if (localdof == locids[k])
          dosomething = true;

      // evaluate function
      if (dosomething)
      {
        double time = 0.0; // dummy time here
        double functfac = 0.0;
        int funct_num = -1;
        if (funct)
        {
          funct_num = (*funct)[0];
          if (funct_num > 0)
            functfac = DRT::Problem::Instance()->Funct(funct_num-1).Evaluate(localdof,actnode->X(),time,&discret);
        }

        // assign value
        const int gid = dofs[j];
        const int lid = (*fieldvector).Map().LID(gid);
        if (lid<0) dserror("Global id %d not on this proc in system vector",gid);
        (*fieldvector)[lid] = functfac;

      } // if dosomething
    }  // loop over nodal DOFs
  }  // loop over nodes

  return;
}
