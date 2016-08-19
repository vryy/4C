/*!------------------------------------------------------------------------------------------------*
 \file ssi_partitioned.cpp

 \brief base class for partitioned scalar structure interaction

 \level 1

 \maintainer Anh-Tu Vuong
             vuong@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289-15264

 *------------------------------------------------------------------------------------------------*/

#include "../linalg/linalg_utils.H"

#include "../drt_adapter/ad_str_wrapper.H"
#include "../drt_adapter/adapter_scatra_base_algorithm.H"

#include "../drt_scatra/scatra_timint_implicit.H"

#include "ssi_partitioned.H"

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
SSI::SSI_Part::SSI_Part(const Epetra_Comm& comm,
    const Teuchos::ParameterList& globaltimeparams)
  : SSI_Base(comm, globaltimeparams)
{
  // Keep this constructor empty!
  // First do everything on the more basic objects like the discretizations, like e.g. redistribution of elements.
  // Only then call the setup to this class. This will call the setup to all classes in the inheritance hierarchy.
  // This way, this class may also override a method that is called during Setup() in a base class.
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void SSI::SSI_Part::SetupSystem()
{
  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void SSI::SSI_Part::Setup(const Epetra_Comm& comm,
    const Teuchos::ParameterList& globaltimeparams,
    const Teuchos::ParameterList& scatraparams,
    const Teuchos::ParameterList& structparams,
    const std::string struct_disname,
    const std::string scatra_disname)
{
  // call setup of base class
  SSI::SSI_Base::Setup(comm, globaltimeparams,scatraparams,structparams,struct_disname, scatra_disname);

  return;
}
