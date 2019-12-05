/*----------------------------------------------------------------------*/
/*! \file
\brief Generic class for all mortar solution strategies

\maintainer Matthias Mayr

\level 2
*/
/*----------------------------------------------------------------------*/

#include <Teuchos_StandardParameterEntryValidators.hpp>
#include "Teuchos_RCP.hpp"
#include "Epetra_SerialComm.h"
#include "mortar_strategy_base.H"
#include "mortar_defines.H"
#include "../linalg/linalg_utils_sparse_algebra_math.H"
#include "../linalg/linalg_sparsematrix.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_inpar/inpar_mortar.H"
#include "../drt_inpar/inpar_structure.H"


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
MORTAR::StratDataContainer::StratDataContainer()
    : probdofs_(Teuchos::null),
      probnodes_(Teuchos::null),
      comm_(Teuchos::null),
      scontact_(),
      dim_(0),
      alphaf_(0.0),
      parredist_(false),
      maxdof_(0),
      systype_(INPAR::CONTACT::system_none),
      dyntype_(INPAR::STR::dyna_statics),
      dynparamN_(0.0)
{
}

/*----------------------------------------------------------------------*
 | ctor (public)                                             popp 01/10 |
 *----------------------------------------------------------------------*/
MORTAR::StrategyBase::StrategyBase(const Teuchos::RCP<MORTAR::StratDataContainer>& data_ptr,
    const Epetra_Map* DofRowMap, const Epetra_Map* NodeRowMap, const Teuchos::ParameterList& params,
    const int spatialDim, const Teuchos::RCP<const Epetra_Comm>& comm, const double alphaf,
    const int maxdof)
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
  Data().ProbDofsPtr() = Teuchos::rcp(new Epetra_Map(*(DofRowMap)));
  Data().ProbNodesPtr() = Teuchos::rcp(new Epetra_Map(*(NodeRowMap)));
  Data().CommPtr() = comm;
  Data().SContact() = params;
  Data().Dim() = spatialDim;
  Data().AlphaF() = alphaf;
  Data().MaxDof() = maxdof;
  Data().SysType() = DRT::INPUT::IntegralValue<INPAR::CONTACT::SystemType>(scontact_, "SYSTEM");
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void MORTAR::StrategyBase::SetTimeIntegrationInfo(
    const double time_fac, const INPAR::STR::DynamicType dyntype)
{
  Data().SetDynType(dyntype);
  if (Data().AlphaF() != 0.0 and time_fac != Data().AlphaF())
    dserror(
        "The time integration factor has been changed already. This error is"
        " supposed to prohibit unintended changes.");
  switch (dyntype)
  {
    case INPAR::STR::dyna_statics:
      Data().SetDynParameterN(0.0);
      break;
    case INPAR::STR::dyna_genalpha:
    case INPAR::STR::dyna_onesteptheta:
      Data().SetDynParameterN(time_fac);
      break;
    default:
      dserror(
          "Unsupported time integration detected! [\"%s\"]", DynamicTypeString(dyntype).c_str());
      exit(EXIT_FAILURE);
  }
}
