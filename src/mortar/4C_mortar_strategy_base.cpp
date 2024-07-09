/*----------------------------------------------------------------------*/
/*! \file
\brief Generic class for all mortar solution strategies


\level 2
*/
/*----------------------------------------------------------------------*/

#include "4C_mortar_strategy_base.hpp"

#include "4C_fem_discretization.hpp"
#include "4C_inpar_mortar.hpp"
#include "4C_inpar_structure.hpp"
#include "4C_linalg_sparsematrix.hpp"
#include "4C_linalg_utils_sparse_algebra_math.hpp"
#include "4C_mortar_defines.hpp"

#include <Teuchos_RCP.hpp>
#include <Teuchos_StandardParameterEntryValidators.hpp>

FOUR_C_NAMESPACE_OPEN


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Mortar::StratDataContainer::StratDataContainer()
    : probdofs_(Teuchos::null),
      probnodes_(Teuchos::null),
      comm_(Teuchos::null),
      scontact_(),
      dim_(0),
      alphaf_(0.0),
      parredist_(false),
      maxdof_(0),
      systype_(Inpar::CONTACT::system_none),
      dyntype_(Inpar::Solid::dyna_statics),
      dynparam_n_(0.0)
{
}

/*----------------------------------------------------------------------*
 | ctor (public)                                             popp 01/10 |
 *----------------------------------------------------------------------*/
Mortar::StrategyBase::StrategyBase(const Teuchos::RCP<Mortar::StratDataContainer>& data_ptr,
    const Epetra_Map* dof_row_map, const Epetra_Map* NodeRowMap,
    const Teuchos::ParameterList& params, const int spatialDim,
    const Teuchos::RCP<const Epetra_Comm>& comm, const double alphaf, const int maxdof)
    : probdofs_(data_ptr->prob_dofs_ptr()),
      probnodes_(data_ptr->prob_nodes_ptr()),
      comm_(data_ptr->comm_ptr()),
      scontact_(data_ptr->s_contact()),
      dim_(data_ptr->n_dim()),
      alphaf_(data_ptr->alpha_f()),
      parredist_(data_ptr->is_par_redist()),
      maxdof_(data_ptr->max_dof()),
      systype_(data_ptr->sys_type()),
      data_ptr_(data_ptr)
{
  // *** set data container variables
  data().prob_dofs_ptr() = Teuchos::rcp(new Epetra_Map(*(dof_row_map)));
  data().prob_nodes_ptr() = Teuchos::rcp(new Epetra_Map(*(NodeRowMap)));
  data().comm_ptr() = comm;
  data().s_contact() = params;
  data().n_dim() = spatialDim;
  data().alpha_f() = alphaf;
  data().max_dof() = maxdof;
  data().sys_type() = Core::UTILS::IntegralValue<Inpar::CONTACT::SystemType>(scontact_, "SYSTEM");
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void Mortar::StrategyBase::set_time_integration_info(
    const double time_fac, const Inpar::Solid::DynamicType dyntype)
{
  // Get weight for contribution from last time step

  data().set_dyn_type(dyntype);
  switch (dyntype)
  {
    case Inpar::Solid::dyna_statics:
      data().set_dyn_parameter_n(0.0);
      break;
    case Inpar::Solid::dyna_genalpha:
    case Inpar::Solid::dyna_onesteptheta:
      data().set_dyn_parameter_n(time_fac);
      break;
    default:
      FOUR_C_THROW(
          "Unsupported time integration detected! [\"%s\"]", DynamicTypeString(dyntype).c_str());
      exit(EXIT_FAILURE);
  }

  // Check if we only want to compute the contact force at the time endpoint
  if (Core::UTILS::IntegralValue<int>(data().s_contact(), "CONTACTFORCE_ENDTIME"))
    alphaf_ = 0.0;
  else
  {
    alphaf_ = data().get_dyn_parameter_n();
  }
}

FOUR_C_NAMESPACE_CLOSE
