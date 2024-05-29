/*---------------------------------------------------------------------*/
/*! \file

\brief Utils methods concerning the discretization


\level 1

*/
/*---------------------------------------------------------------------*/

#include "4C_lib_utils_discret.hpp"

#include "4C_global_data.hpp"
#include "4C_lib_discret.hpp"
#include "4C_lib_node.hpp"
#include "4C_utils_function.hpp"

#include <Epetra_Map.h>

FOUR_C_NAMESPACE_OPEN


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void DRT::UTILS::evaluate_initial_field(const DRT::Discretization& discret,
    const std::string& fieldstring, Teuchos::RCP<Epetra_Vector> fieldvector,
    const std::vector<int>& locids)
{
  // get initial field conditions
  std::vector<CORE::Conditions::Condition*> initfieldconditions;
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
  const std::vector<CORE::Conditions::ConditionType> evaluation_type_order = {
      CORE::Conditions::VolumeInitfield, CORE::Conditions::SurfaceInitfield,
      CORE::Conditions::LineInitfield, CORE::Conditions::PointInitfield};

  for (const auto& type : evaluation_type_order)
  {
    for (const auto& initfieldcondition : initfieldconditions)
    {
      if (initfieldcondition->Type() != type) continue;
      const std::string condstring = initfieldcondition->parameters().Get<std::string>("Field");
      if (condstring != fieldstring) continue;
      DoInitialField(discret, *initfieldcondition, *fieldvector, locids);
    }
  }
}


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void DRT::UTILS::DoInitialField(const DRT::Discretization& discret,
    CORE::Conditions::Condition& cond, Epetra_Vector& fieldvector, const std::vector<int>& locids)
{
  const std::vector<int> cond_nodeids = *cond.GetNodes();
  if (cond_nodeids.empty()) FOUR_C_THROW("Initfield condition does not have nodal cloud.");

  // loop nodes to identify and evaluate spatial distributions
  // of Initfield boundary conditions
  const auto funct_num = cond.parameters().Get<int>("funct");

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
        [&](CORE::Elements::Element* a, CORE::Elements::Element* b)
        { return a->NumDofPerNode(*node) < b->NumDofPerNode(*node); });
    const int numdof = (*ele_with_max_dof)->NumDofPerNode(*node);

    if ((total_numdof % numdof) != 0) FOUR_C_THROW("illegal dof set number");

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
          const double functfac =
              funct_num > 0 ? GLOBAL::Problem::Instance()
                                  ->FunctionById<CORE::UTILS::FunctionOfSpaceTime>(funct_num - 1)
                                  .Evaluate(node->X().data(), time, localdof)
                            : 0.0;

          // assign value
          const int gid = node_dofs[j];
          const int lid = fieldvector.Map().LID(gid);
          if (lid < 0) FOUR_C_THROW("Global id %d not on this proc in system vector", gid);
          fieldvector[lid] = functfac;
        }
      }
    }
  }
}

FOUR_C_NAMESPACE_CLOSE
