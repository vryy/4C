/*---------------------------------------------------------------------*/
/*! \file

\brief Utils methods concerning the discretization


\level 1

*/
/*---------------------------------------------------------------------*/

#include "baci_lib_utils_discret.H"

#include "baci_lib_discret.H"
#include "baci_lib_globalproblem.H"
#include "baci_lib_node.H"

#include <Epetra_Map.h>


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void DRT::UTILS::EvaluateInitialField(const DRT::Discretization& discret,
    const std::string& fieldstring, Teuchos::RCP<Epetra_Vector> fieldvector,
    const std::vector<int>& locids)
{
  // get initial field conditions
  std::vector<DRT::Condition*> initfieldconditions;
  discret.GetCondition("Initfield", initfieldconditions);

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
  // This way, lower entities override higher ones.
  const std::vector<Condition::ConditionType> evaluation_type_order = {
      DRT::Condition::VolumeInitfield, DRT::Condition::SurfaceInitfield,
      DRT::Condition::LineInitfield, DRT::Condition::PointInitfield};

  for (const auto& type : evaluation_type_order)
  {
    for (const auto& initfieldcondition : initfieldconditions)
    {
      if (initfieldcondition->Type() != type) continue;
      const std::string condstring = *initfieldcondition->Get<std::string>("Field");
      if (condstring != fieldstring) continue;
      DoInitialField(discret, *initfieldcondition, *fieldvector, locids);
    }
  }
}


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void DRT::UTILS::DoInitialField(const DRT::Discretization& discret, DRT::Condition& cond,
    Epetra_Vector& fieldvector, const std::vector<int>& locids)
{
  const std::vector<int> cond_nodeids = *cond.Nodes();
  if (cond_nodeids.empty()) dserror("Initfield condition does not have nodal cloud.");

  // loop nodes to identify and evaluate spatial distributions
  // of Initfield boundary conditions
  const auto* funct = cond.Get<std::vector<int>>("funct");
  if (funct->empty()) dserror("Cannot get function.");
  if (funct->size() != 1) dserror("Only one function expected function.");

  for (const int cond_nodeid : cond_nodeids)
  {
    // do only nodes in my row map
    int cond_node_lid = discret.NodeRowMap()->LID(cond_nodeid);
    if (cond_node_lid < 0) continue;
    DRT::Node* node = discret.lRowNode(cond_node_lid);

    // call explicitly the main dofset, i.e. the first column
    std::vector<int> node_dofs = discret.Dof(0, node);
    const int total_numdof = static_cast<int>(node_dofs.size());

    // Get native number of dofs at this node. There might be multiple dofsets
    // (in xfem cases), thus the size of the dofs vector might be a multiple
    // of this.
    auto* const myeles = node->Elements();
    auto* ele_with_max_dof = std::max_element(myeles, myeles + node->NumElement(),
        [&](Element* a, Element* b) { return a->NumDofPerNode(*node) < b->NumDofPerNode(*node); });
    const int numdof = (*ele_with_max_dof)->NumDofPerNode(*node);

    if ((total_numdof % numdof) != 0) dserror("illegal dof set number");

    // now loop over all relevant DOFs
    for (int j = 0; j < total_numdof; ++j)
    {
      int localdof = j % numdof;

      // evaluate function if local DOF id exists
      // in the given locids vector
      for (const int locid : locids)
      {
        if (localdof == locid)
        {
          const double time = 0.0;  // dummy time here
          const int funct_num = (*funct)[0];
          const double functfac =
              funct_num > 0 ? DRT::Problem::Instance()
                                  ->FunctionById<DRT::UTILS::FunctionOfSpaceTime>(funct_num - 1)
                                  .Evaluate(node->X(), time, localdof)
                            : 0.0;

          // assign value
          const int gid = node_dofs[j];
          const int lid = fieldvector.Map().LID(gid);
          if (lid < 0) dserror("Global id %d not on this proc in system vector", gid);
          fieldvector[lid] = functfac;
        }
      }
    }
  }
}
