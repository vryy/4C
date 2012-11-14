/*----------------------------------------------------------------------*/
/*!
\file adapter_structure_bio.cpp

<pre>
Maintainer: Mirella Coroneo
            coroneo@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15236
</pre>
*/
/*----------------------------------------------------------------------*/


#include "ad_str_bio.H"
#include "adapter_coupling.H"
#include "../drt_lib/drt_utils_createdis.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_fsi/fsi_utils.H"
#include "../drt_structure/stru_aux.H"
#include "../drt_adapter/ad_str_fsiwrapper.H"
#include "../drt_ale/ale_utils_mapextractor.H"

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
ADAPTER::StructureBio::StructureBio(
    const Epetra_Comm& comm,
    const Teuchos::ParameterList& prbdyn,
    bool isale,
    const std::string disname,
    const std::string condname
    )
:  AlgorithmBase(comm,prbdyn),
   params_(prbdyn)
{
  // access the structural discretization
  Teuchos::RCP<DRT::Discretization> structdis = DRT::Problem::Instance()->GetDis("structure");
  Teuchos::RCP<ADAPTER::StructureBaseAlgorithm> structurebase = Teuchos::rcp(new ADAPTER::StructureBaseAlgorithm(prbdyn, structdis));
  structure_ = Teuchos::rcp_dynamic_cast<FSIStructureWrapper>(structurebase->StructureFieldrcp());

  if(structure_ == Teuchos::null)
    dserror("cast from ADAPTER::Structure to ADAPTER::FSIStructureWrapper failed");

	Teuchos::RCP<ALE::AleBaseAlgorithm> ale = Teuchos::rcp(new ALE::AleBaseAlgorithm(prbdyn));

	 // set up ale-structure couplings
	const int ndim = DRT::Problem::Instance()->NDim();
        icoupsa_ = Teuchos::rcp(new Coupling());
	 icoupsa_->SetupConditionCoupling(*structure_->Discretization(),
									 structure_->Interface()->FSICondMap(),
									 *ale->AleField().Discretization(),
									 ale->AleField().Interface()->FSICondMap(),
									 condname,
		               ndim);

	  /// do we need this ???
	  //StructureField()->SetMeshMap(coupsa_.MasterDofMap());

	   // the ale matrix might be build just once ???
	  ale->AleField().BuildSystemMatrix();

  return;

}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
ADAPTER::StructureBio::~StructureBio()
{
}
