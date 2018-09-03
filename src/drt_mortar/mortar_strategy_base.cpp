/*!---------------------------------------------------------------------*/
/*!
\file mortar_strategy_base.cpp

\brief Generic class for all mortar solution strategies

\maintainer Alexander Seitz

\level 2

*/
/*----------------------------------------------------------------------*/

#include <Teuchos_StandardParameterEntryValidators.hpp>
#include "Teuchos_RCP.hpp"
#include "Epetra_SerialComm.h"
#include "mortar_strategy_base.H"
#include "mortar_defines.H"
#include "../linalg/linalg_utils.H"
#include "../linalg/linalg_sparsematrix.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_inpar/inpar_mortar.H"


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
MORTAR::StratDataContainer::StratDataContainer()
    : probdofs_(Teuchos::null),
      probnodes_(Teuchos::null),
      comm_(Teuchos::null),
      scontact_(),
      dim_(0),
      alphaf_(0),
      parredist_(false),
      maxdof_(0),
      systype_(INPAR::CONTACT::system_none)
{
}

/*----------------------------------------------------------------------*
 | ctor (public)                                             popp 01/10 |
 *----------------------------------------------------------------------*/
MORTAR::StrategyBase::StrategyBase(const Teuchos::RCP<MORTAR::StratDataContainer>& data_ptr,
    const Epetra_Map* DofRowMap, const Epetra_Map* NodeRowMap, const Teuchos::ParameterList& params,
    int dim, const Teuchos::RCP<const Epetra_Comm>& comm, double alphaf, int maxdof)
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
  Data().Dim() = dim;
  Data().AlphaF() = alphaf;
  Data().MaxDof() = maxdof;
  Data().SysType() = DRT::INPUT::IntegralValue<INPAR::CONTACT::SystemType>(scontact_, "SYSTEM");
}
