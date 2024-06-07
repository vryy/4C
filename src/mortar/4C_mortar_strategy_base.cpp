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
      dyntype_(Inpar::STR::dyna_statics),
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
    : probdofs_(data_ptr->ProbDofsPtr()),
      probnodes_(data_ptr->ProbNodesPtr()),
      comm_(data_ptr->CommPtr()),
      scontact_(data_ptr->SContact()),
      dim_(data_ptr->Dim()),
      alphaf_(data_ptr->AlphaF()),
      parredist_(data_ptr->IsParRedist()),
      maxdof_(data_ptr->MaxDof()),
      systype_(data_ptr->SysType()),
      data_ptr_(data_ptr)
{
  // *** set data container variables
  data().ProbDofsPtr() = Teuchos::rcp(new Epetra_Map(*(dof_row_map)));
  data().ProbNodesPtr() = Teuchos::rcp(new Epetra_Map(*(NodeRowMap)));
  data().CommPtr() = comm;
  data().SContact() = params;
  data().Dim() = spatialDim;
  data().AlphaF() = alphaf;
  data().MaxDof() = maxdof;
  data().SysType() = Core::UTILS::IntegralValue<Inpar::CONTACT::SystemType>(scontact_, "SYSTEM");
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void Mortar::StrategyBase::set_time_integration_info(
    const double time_fac, const Inpar::STR::DynamicType dyntype)
{
  // Get weight for contribution from last time step

  data().SetDynType(dyntype);
  switch (dyntype)
  {
    case Inpar::STR::dyna_statics:
      data().SetDynParameterN(0.0);
      break;
    case Inpar::STR::dyna_genalpha:
    case Inpar::STR::dyna_onesteptheta:
      data().SetDynParameterN(time_fac);
      break;
    default:
      FOUR_C_THROW(
          "Unsupported time integration detected! [\"%s\"]", DynamicTypeString(dyntype).c_str());
      exit(EXIT_FAILURE);
  }

  // Check if we only want to compute the contact force at the time endpoint
  if (Core::UTILS::IntegralValue<int>(data().SContact(), "CONTACTFORCE_ENDTIME"))
    alphaf_ = 0.0;
  else
  {
    alphaf_ = data().GetDynParameterN();
  }
}

FOUR_C_NAMESPACE_CLOSE
