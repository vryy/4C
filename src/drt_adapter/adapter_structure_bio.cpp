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

#ifdef CCADISCRET

#include "adapter_structure_bio.H"
#include "../drt_lib/drt_utils_createdis.H"
#include "../drt_stru_ale/stru_ale_utils.H"
#include "../drt_lib/drt_utils.H"

#include "../drt_fsi/fsi_utils.H"

//#include "../drt_lib/drt_node.H"
//#include "../drt_adapter/adapter_structure.H"
//#include "../drt_ale/ale.H"


extern struct _GENPROB     genprob;


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
ADAPTER::StructureBio::StructureBio(
    Epetra_Comm& comm,
    const Teuchos::ParameterList& prbdyn,
    bool isale,
    const int disnum,
    const std::string condname
    )
:  AlgorithmBase(comm,prbdyn),
   //StructureBaseAlgorithm(prbdyn),
   //ScaTraBaseAlgorithm(prbdyn,isale,disnum),
   //AleBaseAlgorithm(prbdyn),
   params_(prbdyn)
{
	// create ale elements for the structure
	RefCountPtr<DRT::Discretization> aledis = null;

	RCP<DRT::Discretization> structuredis = DRT::Problem::Instance()->Dis(genprob.numsf,0);

	Teuchos::RCP<DRT::UTILS::DiscretizationCreator<STRU_ALE::UTILS::AleStructureCloneStrategy> > alecreator =
   	Teuchos::rcp(new DRT::UTILS::DiscretizationCreator<STRU_ALE::UTILS::AleStructureCloneStrategy>() );

	//alecreator->CreateMatchingDiscretization(structuredis,aledis,-1); //not working

	if (comm.MyPID()==0)
	{
   	  cout << "\n\nCreating ALE copy of the structure ....\n\n";
	}

    Teuchos::RCP<ADAPTER::StructureBaseAlgorithm> structure = Teuchos::rcp(new ADAPTER::StructureBaseAlgorithm(prbdyn));

	Teuchos::RCP<ALE::AleBaseAlgorithm> ale = Teuchos::rcp(new ALE::AleBaseAlgorithm(prbdyn));


	 // set up ale-structure couplings
	 icoupsa_.SetupConditionCoupling(*structure->StructureField().Discretization(),
									 structure->StructureField().Interface().FSICondMap(),
									 *ale->AleField().Discretization(),
									 ale->AleField().Interface().FSICondMap(),
									 condname,
		                             genprob.ndim);

	  // the structure-ale coupling always matches
	  const Epetra_Map* structurenodemap = structure->StructureField().Discretization()->NodeRowMap();
	  const Epetra_Map* alenodemap   = ale->AleField().Discretization()->NodeRowMap();

	  coupsa_.SetupCoupling(*structure->StructureField().Discretization(),
							*ale->AleField().Discretization(),
							*structurenodemap,
							*alenodemap,
	                        genprob.ndim);

	  /// do we need this ???
	  //StructureField().SetMeshMap(coupsa_.MasterDofMap());

	   // the ale matrix might be build just once ???
	  ale->AleField().BuildSystemMatrix();

  return;

}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
ADAPTER::StructureBio::~StructureBio()
{
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::StructureBio::StructAleSolve(
    Teuchos::RCP<Epetra_Vector> idisp)
{

  if (idisp!=Teuchos::null)
  {
    // if we have values at the interface we need to apply them
	  ale->AleField().ApplyInterfaceDisplacements(StructToAle(idisp));
  }

  ale->AleField().Solve();
  Teuchos::RCP<Epetra_Vector> structdisp = AleToStructField(ale->AleField().ExtractDisplacement());

  //change nodes reference position
  ChangeConfig(structdisp);

  // computation of structure solution
  structure->StructureField().Solve();

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::StructureBio::ApplyInterfaceDisplacements(Teuchos::RCP<Epetra_Vector> idisp)
{
	ALE::UTILS::MapExtractor interface_;
	interface_.Setup(*structuredis);
	interface_.InsertFSICondVector(idisp,dispnp_);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::StructureBio::ChangeConfig(Teuchos::RCP<Epetra_Vector> structdisp)
{
	RCP<DRT::Discretization> structuredis = structure->StructureField().Discretization();
	const int numnode = (structuredis->NodeRowMap())->NumMyElements();

	const Epetra_Vector& gvector =*structdisp;

	// loop over all nodes
	for (int i = 0; i < numnode; ++i)
	{
		// get current node
		int gid = (structuredis->NodeRowMap())->GID(i);
		DRT::Node* mynode = structuredis->gNode(gid);

		vector<int> globaldofs = structuredis->Dof(mynode);

		// extract local values from the global vector
		vector<double> nvector(globaldofs.size());
		DRT::UTILS::ExtractMyValues(gvector, nvector, globaldofs);

		mynode->ChangePos(nvector);

	}

	return;

}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::StructureBio::ReadRestart(int step)
{
return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> ADAPTER::StructureBio::AleToStructField(Teuchos::RCP<Epetra_Vector> iv) const
{
  return coupsa_.SlaveToMaster(iv);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> ADAPTER::StructureBio::AleToStructField(Teuchos::RCP<const Epetra_Vector> iv) const
{
  return coupsa_.SlaveToMaster(iv);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> ADAPTER::StructureBio::StructToAle(Teuchos::RCP<Epetra_Vector> iv) const
{
  return icoupsa_.MasterToSlave(iv);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> ADAPTER::StructureBio::StructToAle(Teuchos::RCP<const Epetra_Vector> iv) const
{
  return icoupsa_.MasterToSlave(iv);
}

#endif // CCADISCRET
