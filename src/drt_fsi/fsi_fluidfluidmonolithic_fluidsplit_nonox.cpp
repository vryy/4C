/*!----------------------------------------------------------------------
\file fsi_fluidfluidmonolithic_fluidsplit_nonox.cpp
\brief Class for monolithic fluid-fluid-FSI using XFEM
<pre>
Maintainer:  Shadan Shahmiri
             shahmiri@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289-15265
</pre>
*----------------------------------------------------------------------*/

#include "fsi_fluidfluidmonolithic_fluidsplit_nonox.H"

#include "../drt_adapter/adapter_coupling.H"
#include "../drt_adapter/ad_str_fsiwrapper.H"

#include "fsi_matrixtransform.H"

#include "fsi_debugwriter.H"
#include "fsi_statustest.H"
#include "fsi_monolithic_linearsystem.H"
#include "../linalg/linalg_solver.H"
#include "../linalg/linalg_utils.H"

#include "../drt_lib/drt_colors.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_inpar/inpar_fsi.H"
#include "../drt_inpar/inpar_xfem.H"
#include "../drt_inpar/inpar_ale.H"

#include "../drt_fluid/fluid_utils_mapextractor.H"
#include "../drt_structure/stru_aux.H"
#include "../drt_ale/ale_utils_mapextractor.H"
#include "../drt_ale/ale.H"

#include "../drt_io/io_control.H"
#include "../drt_io/io.H"
#include "../drt_constraint/constraint_manager.H"
#include "../drt_io/io_pstream.H"

/*----------------------------------------------------------------------
 * Constructor:
 *
 * - Remove fluid DBC-DOFs from FSI-interface
 * - Initialize members
 * - Extract and evaluate the parameter list for
 * 	 XFLUID_DYNAMIC and ALE_DYNAMIC parameters
 ----------------------------------------------------------------------*/
FSI::FluidFluidMonolithicFluidSplitNoNOX::FluidFluidMonolithicFluidSplitNoNOX(const Epetra_Comm& comm,
                                                                              const Teuchos::ParameterList& timeparams)
  : MonolithicNoNOX(comm,timeparams)
{

	//Map vector to store both fluid DBC- and fluid FSI-DOFs
	//These maps are supposed to be intersected later on!
	std::vector<Teuchos::RCP<const Epetra_Map> > intersectionmaps;

	//Get map of DBC-DOFs of the fluid field
	intersectionmaps.push_back(FluidField().GetDBCMapExtractor()->CondMap());
	//Get map of fluid FSI-DOFs
	intersectionmaps.push_back(FluidField().Interface()->FSICondMap());


	//Result: The complete set of Fluid-DBCs located on the FSI interface!
	Teuchos::RCP<Epetra_Map> intersectionmap=LINALG::MultiMapExtractor::IntersectMaps(intersectionmaps);

	//Remove the fluid DBCs on the FSI, if there are any!
	if (intersectionmap->NumGlobalElements() != 0)
	{
		// A method for removal of DBC DOFs was not implemented in class XFluidFluid, it had to be added!
		FluidField().RemoveDirichCond(intersectionmap);

		if(comm.MyPID() == 0)
		{
		     IO::cout << "  +---------------------------------------------------------------------------------------------+" << IO::endl;
		     IO::cout << "  |                                        PLEASE NOTE:                                         |" << IO::endl;
		     IO::cout << "  +---------------------------------------------------------------------------------------------+" << IO::endl;
		     IO::cout << "  | You run a monolithic fluid split scheme. Hence, there are no fluid interface DOFs.          |" << IO::endl;
		     IO::cout << "  | Fluid sided DBCs on the FSI interface will be neglected.         			                  |" << IO::endl;
		     IO::cout << "  | Check whether you have prescribed appropriate DBCs on fluid interface DOFs.                 |" << IO::endl;
		     IO::cout << "  +---------------------------------------------------------------------------------------------+" << IO::endl;
		 }

	}


	//Initialization of row/column transformation objects
	//These are needed for the system matrix setup,
	//as matrices from 3 different fields (S,F,A) are set together.

	fggtransform_=Teuchos::rcp(new UTILS::MatrixRowColTransform);
	//fmggtransform_=Teuchos::rcp(new UTILS::MatrixRowTransform);
	fgitransform_=Teuchos::rcp(new UTILS::MatrixRowTransform);
	figtransform_=Teuchos::rcp(new UTILS::MatrixColTransform);
	fmiitransform_=Teuchos::rcp(new UTILS::MatrixColTransform);
	fmgitransform_=Teuchos::rcp(new UTILS::MatrixRowColTransform);
	aigtransform_=Teuchos::rcp(new UTILS::MatrixColTransform);


	//Extract parameter list XFLUID_DYNAMIC
	const Teuchos::ParameterList& xfluiddyn  = DRT::Problem::Instance()->XFluidDynamicParams();


	monolithic_approach_ = DRT::INPUT::IntegralValue<INPAR::XFEM::Monolithic_xffsi_Approach>
							 (xfluiddyn.sublist("GENERAL"),"MONOLITHIC_XFFSI_APPROACH");

	//Initialize the current time step variable
	currentstep_ = 0;

	//Should ALE-relaxation take place? Get flag from parameter list!
	relaxing_ale_ = (bool)DRT::INPUT::IntegralValue<int>(xfluiddyn.sublist("GENERAL"),"RELAXING_ALE");
	//Get no. of steps, after which ALE field should be relaxed
	relaxing_ale_every_ = xfluiddyn.sublist("GENERAL").get<int>("RELAXING_ALE_EVERY");

	//Extract parameter list ALE_DYNAMIC
	const Teuchos::ParameterList& adyn     = DRT::Problem::Instance()->AleDynamicParams();
	//Get the ALE-type index, defining the underlying ALE-algorithm
	int aletype = DRT::INPUT::IntegralValue<int>(adyn,"ALE_TYPE");

	//Just ALE algorithm type 'incremental linear' and 'springs' supports ALE relaxation
	if ( (aletype!=INPAR::ALE::incr_lin and aletype!=INPAR::ALE::springs)
			and monolithic_approach_!=INPAR::XFEM::XFFSI_Full_Newton)
		dserror("Relaxing ALE approach is just possible with Ale-incr-lin!");


	//lambda_ is the Lagrange Multiplier (LM),
	//representing the FSI-Forces
	//It is added to the RHS vector and calculated at every time update.
	//lambda_ is filled in RecoverLangrangeMultiplier().
	//The underlying DOF-map equals to the fluid FSI-DOF-map.
	//The entries of the LM vector are initialized with 0:
	lambda_ = Teuchos::rcp(new Epetra_Vector (*FluidField().Interface()->FSICondMap(),true));

	//Storage for matrices from previous time steps
	fggcurr_  = Teuchos::null;
	fgicurr_  = Teuchos::null;
	fmggcurr_ = Teuchos::null;
	fmgicurr_ = Teuchos::null;

	//Structural predictor step, initially filled with zeros
	ddgpred_ = Teuchos::rcp(new Epetra_Vector(*StructureField()->Interface()->FSICondMap(),true));

	return;
}


/*----------------------------------------------------------------------
 * SetupSystem:
 *   - Extract and evaluate the parameter list for
 * 	   FSI_DYNAMIC parameters
 * 	 - Setup field coupling
 * 	 - System matrix initialization
 * 	 - Build of combined DOF map, including all DOFs
 * 	   except from the fluid FSI DOFs
 *----------------------------------------------------------------------*/
void FSI::FluidFluidMonolithicFluidSplitNoNOX::SetupSystem()
{
	//Extract parameter list FSI_DYNAMIC
	const Teuchos::ParameterList& fsidyn = DRT::Problem::Instance()->FSIDynamicParams();

	//Extract information about the linear block solver (FSIAMG / PreconditionedKrylov)
	linearsolverstrategy_ = DRT::INPUT::IntegralValue<INPAR::FSI::LinearBlockSolver>(fsidyn,"LINEARBLOCKSOLVER");

	//The field coupling methods from namespace ADAPTER
	ADAPTER::Coupling& coupsf  = StructureFluidCoupling();
	ADAPTER::Coupling& coupsa  = StructureAleCoupling();
	ADAPTER::Coupling& coupfa  = FluidAleCoupling();
	ADAPTER::Coupling& icoupfa = InterfaceFluidAleCoupling();

	//Dimensionality of problem
	const int ndim = DRT::Problem::Instance()->NDim();

	/*----------------------------------------------------------------------
	Setup of FSI coupling conditions
	----------------------------------------------------------------------*/

	//Structure to Fluid
	coupsf.SetupConditionCoupling(*StructureField()->Discretization(),
								  StructureField()->Interface()->FSICondMap(),
								  *FluidField().Discretization(),
								  FluidField().Interface()->FSICondMap(),
								  "FSICoupling",
								  ndim);
	//Structure to ALE
	coupsa.SetupConditionCoupling(*StructureField()->Discretization(),
								  StructureField()->Interface()->FSICondMap(),
								  *AleField().Discretization(),
								  AleField().Interface()->FSICondMap(),
								  "FSICoupling",
								  ndim);
	//Fluid to ALE @ interface
	icoupfa.SetupConditionCoupling(*FluidField().Discretization(),
								  FluidField().Interface()->FSICondMap(),
								  *AleField().Discretization(),
								  AleField().Interface()->FSICondMap(),
								  "FSICoupling",
								  ndim);

	//Fluid-ALE coupling always matches --> the maps of S-F and S-A coupling shouldn't differ!
	if ( not coupsf.MasterDofMap()->SameAs(*coupsa.MasterDofMap()) )
		dserror("Structure interface maps do not match");

	if ( coupsf.MasterDofMap()->NumGlobalElements() == 0 )
		dserror("Empty FSI coupling condition???");

	//Before coupling setup, get the node row maps of
	//embedded fluid field and ALE field!
	const Epetra_Map* embfluidnodemap = FluidField().Discretization()->NodeRowMap();
	const Epetra_Map* alenodemap = AleField().Discretization()->NodeRowMap();

	//Fluid to ALE
	coupfa.SetupCoupling(*FluidField().Discretization(),
						*AleField().Discretization(),
						*embfluidnodemap,
						*alenodemap,
						ndim);

	//The fluid field mesh map contains all velocity DOFs covered by an ALE node
	FluidField().SetMeshMap(coupfa.MasterDofMap());

	/*----------------------------------------------------------------------
	Create a combined map for Structure/Fluid/ALE-DOFs all in one!
	----------------------------------------------------------------------*/
	std::vector<Teuchos::RCP<const Epetra_Map> > vecSpaces;
	//Append the structural DOF map
	vecSpaces.push_back(StructureField()->DofRowMap());

	//Create the background fluid maps and the embedded fluid maps merged, without
	//fluid FSI-DOFs
	std::vector<Teuchos::RCP<const Epetra_Map> > vecSpacesFluid;
	//Background fluid DOFs
	vecSpacesFluid.push_back(FluidField().XFluidFluidMapExtractor()->XFluidMap());
	//Embedded fluid DOFs without FSI DOFs
	vecSpacesFluid.push_back(FluidField().Interface()->OtherMap());
	//Merged FSI-free fluid maps
	Teuchos::RCP<const Epetra_Map> fluidmaps=Teuchos::rcp(new Epetra_Map(*LINALG::MultiMapExtractor::MergeMaps(vecSpacesFluid)));
	//Append the final fluid DOF map
	vecSpaces.push_back(fluidmaps);

	//Append ALE DOF map
	vecSpaces.push_back(AleField().Interface()->OtherMap());

	//Non-FSI fluid maps empty??
	if (vecSpaces[1]->NumGlobalElements() == 0)
		dserror("No inner fluid equations. Can't split!");

	//The vector is complete, now fill the system's global block row map
	//with the maps previously set together!
	SetDofRowMaps(vecSpaces);

    // Use normal matrix for fluid equations but build (splitted) mesh movement
	// linearization (shapederivatives) if requested in the input file
	FluidField().UseBlockMatrix(true);

	//Build the ALE-matrix in already splitted system
	AleField().BuildSystemMatrix(false);

	//To avoid incremental ALE errors
	aleresidual_=Teuchos::rcp(new Epetra_Vector(*AleField().Interface()->OtherMap()));

	//Initialize the global system matrix!
	systemmatrix_=
			Teuchos::rcp(
					new LINALG::BlockSparseMatrix<LINALG::DefaultBlockMatrixStrategy>(
							Extractor(),
							Extractor(),
							81,
							false,
							true
							)
						);
}


/*----------------------------------------------------------------------
 * SetupRHS:
 *   - Build the RHS-vector for the Newton-loop!
 *   - The RHS-vector is made up of 2 parts:
 *     the single-field RHS-contributions and special
 *     terms resulting from condensation of fluid-DOFs from the FSI interface
 *     and predictor steps. These terms are added at the first Newton
 *     step only!
 *----------------------------------------------------------------------*/
void FSI::FluidFluidMonolithicFluidSplitNoNOX::SetupRHS(Epetra_Vector& f, bool firstcall)
{

	TEUCHOS_FUNC_TIME_MONITOR("FSI::MonolithicFluidSplit::SetupRHS");
	//f: RHS-vector

	#ifdef DEBUG
		if (FluidField().RHS() == Teuchos::null)
			dserror("empty fluid residual");
	#endif

	//Get the individual field RHS-vectors and adapt them!
	SetupVector(f,
		 	    StructureField()->RHS(),
		 	    FluidField().RHS(),
		    	AleField().RHS(),
		    	FluidField().ResidualScaling());

	//Add the additional ALE residual
	Extractor().AddVector(*aleresidual_,2,f);

	/*----------------------------------------------------------------------
	The following terms are added only at the first Newton iteration!
	----------------------------------------------------------------------*/

	if (firstcall)
	{
		//Store fluid interface velocity:
		Teuchos::RCP<const Epetra_Vector> fveln=FluidField().ExtractInterfaceVeln();

		/*----------------------------------------------------------------------*/
		//Time integration parameters
		/*----------------------------------------------------------------------*/
		//Structure:
		//a*x_n+(1-a)*x_n+1
		//Fluid:
		//b*y_n+(1-b)*y_n+1
		//a: stimintparam
		//b: ftimintparam
		/*----------------------------------------------------------------------*/

		const double stimintparam=StructureField()->TimIntParam();
		const double ftimintparam=FluidField().TimIntParam();

		//Fluid time scaling parameter
		//fluidtimescale: \tau
		// \tau = 1/\Delta t  for Backward-Euler;
		// \tau = 2/\Delta t  for Trapezoidal Rule
		const double fluidtimescale=FluidField().TimeScaling();
		const double fluidresidualscale=FluidField().ResidualScaling();

		//Get the fluid block system matrix "blockf"
		//shape derivative matrix (linearization
		//of Navier-Stokes with respect to mesh movement: moving mesh matrix "mmm")
		//ALE block matrix "blocka"
		//F_...

		//Fluid block system matrix
		const Teuchos::RCP<LINALG::BlockSparseMatrixBase> blockf=FluidField().BlockSystemMatrix();
		if (blockf == Teuchos::null)
			dserror("Expected fluid block matrix...");

		//F^G_...
		//Fluid shape derivative matrix
		const Teuchos::RCP<LINALG::BlockSparseMatrixBase> mmm=FluidField().ShapeDerivatives();

		//A_...
		//ALE block system matrix
		const Teuchos::RCP<LINALG::BlockSparseMatrixBase> blocka=AleField().BlockSystemMatrix();
		if (blocka == Teuchos::null)
			dserror("Expected ALE block matrix...");

		//Extracting submatrices for fluid & ALE from block field matrices:
		//F_{\Gamma\Gamma}, F_{I\Gamma} and A_{I\Gamma}
		const LINALG::SparseMatrix& fgg=blockf->Matrix(1,1);
		const LINALG::SparseMatrix& fig=blockf->Matrix(0,1);
		const LINALG::SparseMatrix& aig=blocka->Matrix(0,1);

	    //Vector for storage of the temporary result for the RHS vector
	    Teuchos::RCP<Epetra_Vector> rhs=Teuchos::null;


	    /*----------------------------------------------------------------------*/
	    //Starting the setup!
	    /*----------------------------------------------------------------------*/

	    //Step 1: Taking care of the  DOFs @ structural side
	    /*
			 * (1)  + (1-stintparam)/(1-flintparam)* \Delta t * timescale * F_{\Gamma\Gamma} * u^{n}_{\Gamma}
		     *
		     * (2)  - (1-stintparam)/(1-flintparam) * F^{G}_{\Gamma\Gamma} * \Delta d_{\Gamma,p}
		     *
		     * (3)  - (1-stintparam)/(1-flintparam) * timescale * F_{\Gamma\Gamma} * \Delta d_{\Gamma,p}
		     *
		     *
		     */

	    //------------------
	    //Build term (1):
	    //------------------

	    //Create zero-filled vector copy based on the row map of F_{\Gamma\Gamma}
	    rhs = Teuchos::rcp(new Epetra_Vector(fgg.RangeMap(),true));

	    //Compute F_{\Gamma\Gamma}*u^n_\Gamma
	    //Write into rhs
	    fgg.Apply(*fveln,*rhs);

	    //Apply scaling
	    rhs->Scale(fluidresidualscale * (1.0-stimintparam)/(1.0-ftimintparam) * Dt() *fluidtimescale);

	    //Insert into structural side of the interface
	    rhs=FluidToStruct(rhs);
	    rhs=StructureField()->Interface()->InsertFSICondVector(rhs);

	    //Add to the structure block (index 0) of the final RHS vector f
	    Extractor().AddVector(*rhs,0,f);

	    if (mmm != Teuchos::null)
	    {
		    //------------------
		    //Build term (2):
	    	//------------------

	    	//Extract the matrix F^G_\Gamma_\Gamma
	    	const LINALG::SparseMatrix& fmgg=mmm->Matrix(1,1);

	    	//Re-initialize rhs
	    	rhs=Teuchos::rcp(new Epetra_Vector(fmgg.RangeMap(),true));


		    //Compute F^{G}_{\Gamma\Gamma} * \Delta d_{\Gamma,p}
		    //Write into rhs
		    fmgg.Apply(*StructToFluid(ddgpred_),*rhs);

		    //Apply scaling
		    rhs->Scale(-1.0*(1.0-stimintparam)/(1.0-ftimintparam));

		    //Insert into structure side of the interface
		   	rhs=FluidToStruct(rhs);
		    rhs=StructureField()->Interface()->InsertFSICondVector(rhs);

		    Extractor().AddVector(*rhs,0,f);
	    }

	    //------------------
	    //Build term (3):
    	//------------------

	    //Re-initialize rhs
	    rhs=Teuchos::rcp(new Epetra_Vector(fgg.RangeMap(),true));

	    //Compute F_{\Gamma\Gamma} * \Delta d_{\Gamma,p}
	    //Write into rhs
	    fgg.Apply(*StructToFluid(ddgpred_),*rhs);

	    //Apply scaling
	    rhs->Scale(-1.0*fluidresidualscale*(1.0-stimintparam)/(1.0-ftimintparam)*fluidtimescale);

	    //Insert into structure side of the interface
	    rhs=FluidToStruct(rhs);
	    rhs=StructureField()->Interface()->InsertFSICondVector(rhs);

	    Extractor().AddVector(*rhs,0,f);


	    //Step 2: Fluid sided inner DOFs
	    /*
	     * (1) F_{I\Gamma} * \Delta t* timescale * u^{n}_{\Gamma}
	     *
	     * (2)  - timescale * F_{I\Gamma} * \Delta d_{\Gamma,p}
	     *
	     * (3) - F^{G}_{I\Gamma} * \Delta d_{\Gamma,p}
	     *
		 */

		//------------------
		//Build term (1):
		//------------------

	    //Re-initialize rhs
		rhs=Teuchos::rcp(new Epetra_Vector(fig.RangeMap(),true));

	    //Compute term F_{I\Gamma} *u^{n}_{\Gamma}
	    //Write into rhs
		fig.Apply(*fveln,*rhs);

		//Apply scaling
		rhs->Scale(fluidtimescale * Dt());

		//Add to the final vector f, to the fluid block (index 1)
		Extractor().AddVector(*rhs,1,f);


		//------------------
		//Build term (2):
		//------------------

		//Re-initialize rhs
		rhs=Teuchos::rcp(new Epetra_Vector(fig.RangeMap(),true));

	    //Compute term F_{I\Gamma} * \Delta d_{\Gamma,p}
	    //Write into rhs
		fig.Apply(*StructToFluid(ddgpred_),*rhs);

		//Apply scaling
		rhs->Scale(-fluidtimescale);

		Extractor().AddVector(*rhs,1,f);

		if (mmm != Teuchos::null)
		{
			//------------------
			//Build term (3):
			//------------------

			//Extract the matrix F^G_I\Gamma
			const LINALG::SparseMatrix& fmig=mmm->Matrix(0,1);

			//Re-initialize rhs
			rhs=Teuchos::rcp(new Epetra_Vector(fmig.RangeMap(),true));

		    //Compute F^{G}_{I\Gamma} * \Delta d_{\Gamma,p}
		    //Write into rhs
			fmig.Apply(*StructToFluid(ddgpred_),*rhs);

			rhs=FluidField().Interface()->InsertOtherVector(rhs);
			rhs=xfluidfluidsplitter_->InsertFluidVector(rhs);
			//Remove the FSI DOFs, to match the FSI-DOF-free map of the fluid block from f
			Teuchos::RCP<Epetra_Map> fsimap=Teuchos::rcp(new Epetra_Map(*FluidField().Interface()->FSICondMap()));
			Teuchos::RCP<LINALG::MapExtractor> ffsextractor=Teuchos::rcp(new LINALG::MapExtractor(*(FluidField().DofRowMap()),fsimap,true));
			rhs=ffsextractor->ExtractOtherVector(rhs);

			//Apply scaling
			rhs->Scale(-1.0);

			Extractor().AddVector(*rhs,1,f);
		}

		//Step 3: Inner ALE DOFs
		//
		//Adding terms to ALE-Part of RHS-vector
		//-A_{I\Gamma} * \Delta d_{\Gamma,p}

		//Re-initialize rhs
		rhs=Teuchos::rcp(new Epetra_Vector(aig.RangeMap(),true));

		//Compute term A_{I\Gamma} * \Delta d_{\Gamma,p}
		//Write into rhs
		aig.Apply(*StructToAle(ddgpred_),*rhs);

		//Apply scaling
		rhs->Scale(-1.0);

		Extractor().AddVector(*rhs,2,f);


		// Finally, apply the DBCs!
		// structure
		rhs=Extractor().ExtractVector(f,0);
		Teuchos::RCP<Epetra_Vector> zeros=Teuchos::rcp(new Epetra_Vector(rhs->Map(),true));
		LINALG::ApplyDirichlettoSystem(rhs,zeros,*(StructureField()->GetDBCMapExtractor()->CondMap()));
		Extractor().InsertVector(*rhs,0,f);

		// fluid
		rhs=Extractor().ExtractVector(f,1);
		zeros=Teuchos::rcp(new Epetra_Vector(rhs->Map(),true));
		LINALG::ApplyDirichlettoSystem(rhs,zeros,*(FluidField().GetDBCMapExtractor()->CondMap()));
		Extractor().InsertVector(*rhs,1,f);

		// ale
		rhs=Extractor().ExtractVector(f,2);
		zeros=Teuchos::rcp(new Epetra_Vector(rhs->Map(),true));
		LINALG::ApplyDirichlettoSystem(rhs,zeros,*(AleField().GetDBCMapExtractor()->CondMap()));
		Extractor().InsertVector(*rhs,2,f);
	}

}

/*----------------------------------------------------------------------
 * SetupSystemMatrix:
 *
 *   - Build the final system block matrix extracting, scaling
 *     & transforming the single field submatrices
 ----------------------------------------------------------------------*/
void FSI::FluidFluidMonolithicFluidSplitNoNOX::SetupSystemMatrix()
{
	TEUCHOS_FUNC_TIME_MONITOR("FSI::MonolithicFluidSplit::SetupSystemMatrix");

	//Short explanation of domain- & rangemap:
	//Ax=y; x-->domain, y-->range

	/*----------------------------------------------------------------------*/
	//Time integration parameters
	/*----------------------------------------------------------------------*/
	//Structure:
	//a*x_n+(1-a)*x_n+1
	//Fluid:
	//b*y_n+(1-b)*y_n+1
	//a: stimintparam
	//b: ftimintparam
	/*----------------------------------------------------------------------*/

	const double stimintparam=StructureField()->TimIntParam();
	const double ftimintparam=FluidField().TimIntParam();

	//Fluid time scaling parameter
	//fluidtimescale: \tau
	// \tau = 1/\Delta t  for Backward-Euler;
	// \tau = 2/\Delta t  for Trapezoidal Rule
	const double fluidtimescale=FluidField().TimeScaling();

	const double fluidresidualscale=FluidField().ResidualScaling();

	const ADAPTER::Coupling& coupsf=StructureFluidCoupling();
	const ADAPTER::Coupling& coupsa=StructureAleCoupling();
	const ADAPTER::Coupling& coupfa=FluidAleCoupling();

	/*----------------------------------------------------------------------*/
	//Extract Jacobians for all fields!
	/*----------------------------------------------------------------------*/
	//Structural matrix is SparseMatrix, as it is not splitted in this
	//algorithm!
	Teuchos::RCP<LINALG::SparseMatrix> s=StructureField()->SystemMatrix();

	//Extract fluid & ALE system matrices:
	//there has to be a reordering of blocks in the fluid split
	//scheme, therefore matrices f & a need to be of block type!
	const Teuchos::RCP<LINALG::BlockSparseMatrixBase> f=FluidField().BlockSystemMatrix();
	const Teuchos::RCP<LINALG::BlockSparseMatrixBase> a=AleField().BlockSystemMatrix();

	//Mesh motion matrix (known with subscript G)
	//Results from linearization with respect to mesh movement!
	const Teuchos::RCP<LINALG::BlockSparseMatrixBase> mmm=FluidField().ShapeDerivatives();

	#ifdef DEBUG
		if (s == Teuchos::null)
			dserror("expected structure matrix!");
		if (f == Teuchos::null)
			dserror("expected fluid matrix!");
		if (a == Teuchos::null)
			dserror("expected ale matrix!");
	#endif

	//Fluid submatrices
	//F_{II}
	LINALG::SparseMatrix& fii=f->Matrix(0,0);
	//F_{I\Gamma}
	LINALG::SparseMatrix& fig=f->Matrix(0,1);
	//F_{\GammaI}
	LINALG::SparseMatrix& fgi=f->Matrix(1,0);
	//F_{\Gamma\Gamma}
	LINALG::SparseMatrix& fgg=f->Matrix(1,1);

	//ALE-submatrices
	//A_{II}
	LINALG::SparseMatrix& aii=a->Matrix(0,0);
	//A_{I\Gamma}
	LINALG::SparseMatrix& aig=a->Matrix(0,1);

	/*----------------------------------------------------------------------*/
	//Build the final block system matrix
	/*----------------------------------------------------------------------*/

	s->UnComplete();

	//F_{\Gamma\Gamma} : Scale & transform!
	(*fggtransform_)(fgg,
					(1.0-stimintparam)/(1.0-ftimintparam)*fluidtimescale*fluidresidualscale,
					ADAPTER::CouplingSlaveConverter(coupsf),
					ADAPTER::CouplingSlaveConverter(coupsf),
					*s,
					true,
					true);

	//Creating the matrix F_{\GammaI}, the DBCs have to be inserted!
	//Therefore, an auxiliary matrix lfgi is created
	//based on a clone of s (inserted in the same rows as sgi)!
	Teuchos::RCP<LINALG::SparseMatrix> lfgi=Teuchos::rcp(new LINALG::SparseMatrix(s->RowMap(),81,false));
//---------------------------------------------
	//F_{\GammaI} : Scale & transform!
	(*fgitransform_)(fgi,
			 	 	 (1.0-stimintparam)/(1.0-ftimintparam)*fluidresidualscale,
			 	 	 ADAPTER::CouplingSlaveConverter(coupsf),
			 	 	 *lfgi);

	//Complete() does the FillComplete-call for the submatrices
	lfgi->Complete(fgi.DomainMap(),s->RangeMap());

	//Apply the (structural) DBCs before assembling into system matrix.
	lfgi->ApplyDirichlet(*(StructureField()->GetDBCMapExtractor()->CondMap()),false);

	//Insert into final matrix
	systemmatrix_->Assign(0,1,View,*lfgi);


	//F_{I\Gamma} : Scale, transform & directly insert!
    (*figtransform_)(f->FullRowMap(),
                     f->FullColMap(),
                     fig,
                     fluidtimescale,
                     ADAPTER::CouplingSlaveConverter(coupsf),
                     systemmatrix_->Matrix(1,0));

	//F_{II} : Insert!
	systemmatrix_->Assign(1,1,View,fii);

	//A_{I\Gamma} : Scale, transform & directly insert!
	(*aigtransform_)(a->FullRowMap(),
					 a->FullColMap(),
					 aig,
					 1.0,
					 ADAPTER::CouplingSlaveConverter(coupsa),
					 systemmatrix_->Matrix(2,0));

	//A_{II} : Insert!
	systemmatrix_->Assign(2,2,View,aii);

	/*----------------------------------------------------------------------*/
	// add optional fluid linearization with respect to mesh motion block
	if	(mmm != Teuchos::null)
	{
		//Extracting submatrices
		//F^G_{II}
		LINALG::SparseMatrix& fmii=mmm->Matrix(0,0);
		//F^G_{I\Gamma}
		LINALG::SparseMatrix& fmig=mmm->Matrix(0,1);
		//F^G_{\GammaI}
		LINALG::SparseMatrix& fmgi=mmm->Matrix(1,0);
		//F^G_{\Gamma\Gamma}
		LINALG::SparseMatrix& fmgg=mmm->Matrix(1,1);

		//F^G_{I\Gamma} : Scale, transform & directly insert!
		//no exact match, addition true!
		(*figtransform_) (f->FullRowMap(),
						  f->FullColMap(),
						  fmig,
						  1.0,
						  ADAPTER::CouplingSlaveConverter(coupsf),
						  systemmatrix_->Matrix(1,0),
						  false,
						  true);

		//F^G_{\Gamma\Gamma} : Scale & transform!
/*		(*fmggtransform_)(fmgg,
			   			(1.0-stimintparam)/(1.0-ftimintparam)*fluidresidualscale,
						ADAPTER::CouplingSlaveConverter(coupsf),
						*s,
						true);*/

		(*fggtransform_)(fmgg,
						(1.0-stimintparam)/(1.0-ftimintparam)*fluidresidualscale,
						ADAPTER::CouplingSlaveConverter(coupsf),
						ADAPTER::CouplingSlaveConverter(coupsf),
						*s,
						false,
						true);

		//F^G_{II} : Transform & directly insert!
		(*fmiitransform_)(mmm->FullRowMap(),
						  mmm->FullColMap(),
						  fmii,
						  1.0,
						  ADAPTER::CouplingMasterConverter(coupfa),
						  systemmatrix_->Matrix(1,2),
						  false);


		Teuchos::RCP<LINALG::SparseMatrix> lfmgi=Teuchos::rcp(new LINALG::SparseMatrix(s->RowMap(),81,false));
		//F^G_{\GammaI} : Scale & transform!
		(*fmgitransform_)(fmgi,
						 (1.0-stimintparam)/(1.0-ftimintparam)*fluidresidualscale,
						 ADAPTER::CouplingSlaveConverter(coupsf),
						 ADAPTER::CouplingMasterConverter(coupfa),
						 *lfmgi,
						 false,
						 false);

		lfmgi->Complete(aii.DomainMap(),s->RangeMap());
		lfmgi->ApplyDirichlet(*(StructureField()->GetDBCMapExtractor()->CondMap()),false);

		systemmatrix_->Assign(0,2,View,*lfmgi);
	}

	/*----------------------------------------------------------------------
	Prepare structure matrix & finally apply DBCs
    ----------------------------------------------------------------------*/

	s->Complete();
	s->ApplyDirichlet(*(StructureField()->GetDBCMapExtractor()->CondMap()),true);
	s->UnComplete();

	//And finally add the structural (unsplitted) matrix
	systemmatrix_->Assign(0,0,View,*s);

	//Done!
	systemmatrix_->Complete();

	//Store some submatrices required for RHS-Setup
	//to know them for LM-recovery!
	fgicurr_=Teuchos::rcp(new LINALG::SparseMatrix(f->Matrix(1,0)));
	fggcurr_=Teuchos::rcp(new LINALG::SparseMatrix(f->Matrix(1,1)));


	if (mmm != Teuchos::null)
	{
		fmgicurr_=Teuchos::rcp(new LINALG::SparseMatrix(mmm->Matrix(1,0)));
		fmggcurr_=Teuchos::rcp(new LINALG::SparseMatrix(mmm->Matrix(1,1)));
	}
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/

void FSI::FluidFluidMonolithicFluidSplitNoNOX::InitialGuess(Teuchos::RCP<Epetra_Vector> ig)
{
	TEUCHOS_FUNC_TIME_MONITOR("FSI::MonolithicFluidSplit::InitialGuess");

	SetupVector(*ig,
				StructureField()->InitialGuess(),
				FluidField().InitialGuess(),
				AleField().InitialGuess(),
				0.0);
}

/*----------------------------------------------------------------------
 * CombinedDBCMap:
 *
 *   - Creates map containing the DOFs with Dirichlet-BCs from all fields!
 *----------------------------------------------------------------------*/

Teuchos::RCP<Epetra_Map> FSI::FluidFluidMonolithicFluidSplitNoNOX::CombinedDBCMap()
{
	//As this is a fluid split scheme, the structural DBC map remains
	//unchanged. If there are any DBCs to be prescribed on the FSI,
	//this has to be done for the structure field.
	//Possible fluid DBCs on the FSI will be removed from the final map.

	/*----------------------------------------------------------------------*/
	//Get the DBC-maps for each field
	/*----------------------------------------------------------------------*/

	Teuchos::RCP<const Epetra_Map> scondmap=StructureField()->GetDBCMapExtractor()->CondMap();
	Teuchos::RCP<const Epetra_Map> ffcondmap=FluidField().FluidDirichMaps();
	Teuchos::RCP<const Epetra_Map> acondmap=AleField().GetDBCMapExtractor()->CondMap();

	//Create a combined map vector with the 3 field DBC maps
	std::vector<Teuchos::RCP<const Epetra_Map> > allfsidbcmapvector;
	allfsidbcmapvector.push_back(scondmap);
	allfsidbcmapvector.push_back(ffcondmap);
	allfsidbcmapvector.push_back(acondmap);

	//Merge the maps together
	Teuchos::RCP<Epetra_Map> alldbcmaps=LINALG::MultiMapExtractor::MergeMaps(allfsidbcmapvector);


	/*----------------------------------------------------------------------*/
	//Info: fullmap_  is built from the merged field DOF maps in SetDofRowMaps(),
	//which is called in Setup[New]System().
	//If the local ID (LID) in fullmap_ for a DOF characterized by a certain global ID (GID)
	//is -1 it means, that the DOF is not located on any processor: it's a fluid FSI-DOF then,
	//as we run a fluid split scheme!
	//This check is performed below.

	/*----------------------------------------------------------------------*/

	std::vector<int> otherdbcmapvector;
	const int allmapslen=alldbcmaps->NumMyElements(); //on each prossesor (lids)
	const int* allmapsgids=alldbcmaps->MyGlobalElements();

	for ( int i=0 ; i<allmapslen ; i++ )
	{
		//Run through all global IDs in alldbcmaps successively
		int gid=allmapsgids[i];
		//Get the local ID for every such GID
		int fullmaplid=fullmap_->LID(gid);
		//LID passed the check? Then it is NOT an FSI-DBC-DOF among the DBC-DOFs
		//and can be appended.
		if(fullmaplid>=0)
			otherdbcmapvector.push_back(gid);
	}

	Teuchos::RCP<Epetra_Map> otherdbcmap=Teuchos::rcp(new Epetra_Map(-1,otherdbcmapvector.size(),&otherdbcmapvector[0],0,Comm()));

	return otherdbcmap;

}

/*----------------------------------------------------------------------
 * SetupVector:
 *    - Called in every Newton iteration via SetupRHS()
 *    - Extract RHS-contribution from the fields,
 *      scale them and insert them.
 *    - Add the scaled Lagrange multiplier!
----------------------------------------------------------------------*/
void FSI::FluidFluidMonolithicFluidSplitNoNOX::SetupVector(Epetra_Vector &f,
                                                          Teuchos::RCP<const Epetra_Vector> sv,
                                                          Teuchos::RCP<const Epetra_Vector> fv,
                                                          Teuchos::RCP<const Epetra_Vector> av,
                                                          double fluidscale)
{
	TEUCHOS_FUNC_TIME_MONITOR("FSI::FluidFluidMonolithicFluidSplitNoNOX::SetupVector");


	//Writes the following entries into f :

	//[
	//f^S_I
	//-----
	//f^S_{\Gamma}+
	//(1-stintparam_)/(1-flintparam)*fluidscale*f^F_{\Gamma}
	//lambda*(stintparam_-flintparam*(1-stintparam_)/(1-flintparam)
	//-----
	//f^ F_I
	//-----
	//0
	//]

	/*----------------------------------------------------------------------*/
	//Time integration parameters
	/*----------------------------------------------------------------------*/
	//Structure:
	//a*x_n+(1-a)*x_n+1
	//Fluid:
	//b*y_n+(1-b)*y_n+1
	//a: stimintparam
	//b: ftimintparam
	/*----------------------------------------------------------------------*/

	const double stimintparam=StructureField()->TimIntParam();
	const double ftimintparam=FluidField().TimIntParam();


	//Extract inner DOFs for ALE-field
	Teuchos::RCP<Epetra_Vector> aov=AleField().Interface()->ExtractOtherVector(av);

	//Separate vectors into inner fluid and background fluid contribution!
	Teuchos::RCP<Epetra_Vector> ffluidv=xfluidfluidsplitter_->ExtractFluidVector(fv);
	Teuchos::RCP<Epetra_Vector> fxfluidv=xfluidfluidsplitter_->ExtractXFluidVector(fv);

	if (fabs(fluidscale) >= EPS15)
	{
		//Get the FSI-interface RHS-vector for the fluid side!
		Teuchos::RCP<Epetra_Vector> fcv=FluidField().Interface()->ExtractFSICondVector(ffluidv);

		//Convert previously extracted vector to structure !
		Teuchos::RCP<Epetra_Vector> modsv=StructureField()->Interface()->InsertFSICondVector(FluidToStruct(fcv));

		//Add the converted interface RHS-contributions (scaled) to the global structural RHS!
		int err=modsv->Update(1.0,*sv,(1.0-stimintparam)/(1.0-ftimintparam)*fluidscale);
		if (err)
			dserror("Update of structural residual vector failed! Error code %i",err);

		//Add the previous Lagrange Multiplier
		if (lambda_ != Teuchos::null)
		{
			Teuchos::RCP<Epetra_Vector> lambdaglob=StructureField()->Interface()->InsertFSICondVector(FluidToStruct(lambda_));
			err=modsv->Update(stimintparam-ftimintparam*(1.0-stimintparam)/(1.0-ftimintparam),*lambdaglob,1.0);
			if (err)
				dserror("Update of structural residual vector failed! Error code %i",err);
		}

		//Apply DBCs to system
		//Create zero-filled Epetra_Vector from modsv's map (bool 'zeros'= true!)
		Teuchos::RCP<const Epetra_Vector> zeros=Teuchos::rcp(new const Epetra_Vector(modsv->Map(),true));
		LINALG::ApplyDirichlettoSystem(modsv,zeros,*(StructureField()->GetDBCMapExtractor()->CondMap()));

		//Insert structural contribution
		Extractor().InsertVector(*modsv,0,f);
	}
	else
	{
		Extractor().InsertVector(*sv,0,f);
	}

	//Create a global temporary fluid vector to assemble all entries into
	Teuchos::RCP<Epetra_Vector> fglobalv=LINALG::CreateVector(*FluidField().DofRowMap(),true);

	//Insert fluid contribution
	xfluidfluidsplitter_->InsertXFluidVector(fxfluidv,fglobalv);
	xfluidfluidsplitter_->InsertFluidVector(ffluidv,fglobalv);

	//Performing a fluid split, the fluid FSI DOFs are removed from the system.
	//The fluid field RHS still contains these DOFs.
	//A MapExtractor() object is created in order to split
	//the fluid RHS into FSI- and non-FSI-DOFs!
	Teuchos::RCP<Epetra_Map> fsimap=Teuchos::rcp(new Epetra_Map(*FluidField().Interface()->FSICondMap()));
	Teuchos::RCP<LINALG::MapExtractor> ffsextractor=Teuchos::rcp(new LINALG::MapExtractor(*(FluidField().DofRowMap()),fsimap,true));
	//Remove the FSI DOFs from the global fluid vector
	fglobalv=ffsextractor->ExtractOtherVector(fglobalv);

	//Insert fluid contribution
	Extractor().InsertVector(*fglobalv,1,f);

	//Insert ALE contribution
	Extractor().InsertVector(*aov,2,f);

	return;
}
/*----------------------------------------------------------------------
 *
 * ExtractFieldVectors:
 *    - Called from Evaluate() method in Newton-loop with x=x_sum_
 *      (increment sum)
 *    - Field contributions sx,fx,ax are recovered from x
 *    - The embedded fluid FSI-DOFs are recovered from
 *    	structural interface DOFs
 *
 ----------------------------------------------------------------------*/

void FSI::FluidFluidMonolithicFluidSplitNoNOX::ExtractFieldVectors(Teuchos::RCP<const Epetra_Vector>  x,
																   Teuchos::RCP<const Epetra_Vector>& sx,
                                                                   Teuchos::RCP<const Epetra_Vector>& fx,
                                                                   Teuchos::RCP<const Epetra_Vector>& ax)
{
	TEUCHOS_FUNC_TIME_MONITOR("FSI::MonolithicFluidSplit::ExtractFieldVectors");

	/*----------------------------------------------------------------------*/
	//Process structure unknowns
	/*----------------------------------------------------------------------*/
	//Extract whole structure field vector
	sx=Extractor().ExtractVector(x,0);

	//Structural part of FSI interface
	Teuchos::RCP<Epetra_Vector> scx=StructureField()->Interface()->ExtractFSICondVector(sx);

	/*----------------------------------------------------------------------*/
	//Process ALE unknowns
	/*----------------------------------------------------------------------*/
	Teuchos::RCP<const Epetra_Vector> aox=Extractor().ExtractVector(x,2);
	//Update interface part of structure vector with predictor increment
	scx->Update(1.0,*ddgpred_,1.0);
	Teuchos::RCP<Epetra_Vector> acx=StructToAle(scx);

	Teuchos::RCP<Epetra_Vector> a=AleField().Interface()->InsertOtherVector(aox);
	//Insert the FSI-DOF vector into full vector a
	AleField().Interface()->InsertFSICondVector(acx,a);
	//Write a into passed argument ax
	ax=a;

	/*----------------------------------------------------------------------*/
	//Process fluid unknowns
	/*----------------------------------------------------------------------*/
	//Extract vector of fluid unknowns from x
	Teuchos::RCP <const Epetra_Vector> fox=Extractor().ExtractVector(x,1);

	//Conversion ALE displacement to fluid field:
	Teuchos::RCP<Epetra_Vector> fcx=AleToFluidInterface(acx);
	//fcx contains displacements now -
	//they have to be converted to velocities.
	FluidField().DisplacementToVelocity(fcx);


	//The previously computed fluid interface values have to be inserted into the fluid field vector
	//Again a FSI/non-FSI- MapExtractor() is necessary for this task (-->ffsextractor).
	Teuchos::RCP<Epetra_Map> fsimap=Teuchos::rcp(new Epetra_Map(*FluidField().Interface()->FSICondMap()));
	Teuchos::RCP<LINALG::MapExtractor> ffsextractor=Teuchos::rcp(new LINALG::MapExtractor(*(FluidField().DofRowMap()),fsimap,true));
	Teuchos::RCP<Epetra_Vector> f=ffsextractor->InsertCondVector(fcx);
	ffsextractor->InsertOtherVector(fox,f);

	fx=f;

	/*----------------------------------------------------------------------*/
	//Compute the increments needed for recovery of Lagrange Multiplier!
	/*----------------------------------------------------------------------*/
	//Update inner fluid velocity increment
	/*if (solipre_ != Teuchos::null)
		duiinc_->Update(1.0,*fox,-1.0,*solipre_,0.0);
	else
		solipre_=fox;

	//Update inner ale displacement increment
	if (solialepre_ != Teuchos::null)
		ddialeinc_->Update(1.0,*aox,-1.0,*solialepre_,0.0);
	else
		solialepre_=aox;

	//Update structural displacement increment at the interface
	if (solgpre_ != Teuchos::null)
		ddginc_->Update(1.0,*scx,-1.0,*solgpre_,0.0);
	else
		solgpre_=scx;*/
}

/*----------------------------------------------------------------------
  PrepareTimeStep:
 * 	- Prepare the Newton loop!
 * 	- Increment time step and value
 * 	- Call the internal PrepareTimeStep() method for each of the 3 fields
 *	- If monolithic type of approach is not XFFSI-FullNewton --> setup new system!
 *	- (options for the monolithic approach are: full newton,
 *	  fixed ale interpolation and fixed ale partitioned)
----------------------------------------------------------------------*/
void FSI::FluidFluidMonolithicFluidSplitNoNOX::PrepareTimeStep()
{
	TEUCHOS_FUNC_TIME_MONITOR("FSI::MonolithicFluidSplit::PrepareTimeStep");

	IncrementTimeAndStep();
	PrintHeader();

	StructureField()->PrepareTimeStep();
	FluidField().PrepareTimeStep();
	AleField().PrepareTimeStep();

	if ( monolithic_approach_ != INPAR::XFEM::XFFSI_Full_Newton )
		SetupNewSystem();

	xfluidfluidsplitter_=FluidField().XFluidFluidMapExtractor();
}

/*----------------------------------------------------------------------
 * Update:
 *    - Initiates LM-recovery
 *    - Calls Update() for each of the 3 fields
 *    - Initiates ALE-relaxation
----------------------------------------------------------------------*/
void FSI::FluidFluidMonolithicFluidSplitNoNOX::Update()
{
	TEUCHOS_FUNC_TIME_MONITOR("FSI::MonolithicFluidSplit::Update");

	//ALE relaxation flag
	bool aleupdate=false;

	//In case of ale relaxation: Check, if the current step is an ale relaxation step!
	if (relaxing_ale_)
		if (currentstep_%relaxing_ale_every_ == 0) aleupdate=true;

	//In case of ALE relaxation:
	if ( monolithic_approach_ != INPAR::XFEM::XFFSI_Full_Newton and aleupdate )
	{
		//Set the old state of ALE displacement before relaxation
		FluidField().ApplyEmbFixedMeshDisplacement(AleToFluid(AleField().WriteAccessDispnp()));
	}

	RecoverLagrangeMultiplier();

	currentstep_++;

	//In case of ALE relaxation
	if ( monolithic_approach_ != INPAR::XFEM::XFFSI_Full_Newton and aleupdate )
	{
		if (Comm().MyPID() == 0)
			IO::cout << "Relaxing ALE!" << IO::endl;
		//Set the ALE FSI-DOFs to Dirichlet and solve ALE system again
		//to obtain the true ALE displacement
		AleField().SolveAleXFluidFluidFSI();
		//Now apply the ALE-displacement to the (embedded!) fluid and update the
		//grid velocity!
		FluidField().ApplyMeshDisplacement(AleToFluid(AleField().WriteAccessDispnp()));
	}

	//Call Update() method for all fields
	StructureField()->Update();
	FluidField().Update();
	AleField().Update();

	if ( monolithic_approach_ != INPAR::XFEM::XFFSI_Full_Newton and aleupdate )
	{
		//Build the ALE-matrix after the update
		AleField().BuildSystemMatrix(false);
		//Create new ALE interface residual
		aleresidual_=Teuchos::rcp(new Epetra_Vector(*AleField().Interface()->OtherMap()));
	}
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::FluidFluidMonolithicFluidSplitNoNOX::ReadRestart(int step)
{
  //Read Lagrange Multiplier
  {
    Teuchos::RCP<Epetra_Vector> lambdafull = Teuchos::rcp(new Epetra_Vector(*(FluidField().XFluidFluidMapExtractor()->FluidMap()),true));
    //FluidField().Discretization() is embedded discretization!
    IO::DiscretizationReader reader = IO::DiscretizationReader(FluidField().Discretization(),step);
    reader.ReadVector(lambdafull, "fsilambda");
    lambda_=FluidField().Interface()->ExtractFSICondVector(lambdafull);
  }

  StructureField()->ReadRestart(step);
  FluidField().ReadRestart(step);
  AleField().ReadRestart(step);

  SetTimeStep(FluidField().Time(),FluidField().Step());
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::FluidFluidMonolithicFluidSplitNoNOX::Output()
{
	  StructureField()->Output();
	  FluidField().Output();

	  // Output of Lagrange Multiplier
	  {
	    /* 'lambda_' is only defined on the interface. So, insert 'lambda_' into
	     * 'lambdafull' that is defined on the entire (INNER!!) fluid field. Then, write
	     * output or restart data.
	     */
		Teuchos::RCP<Epetra_Vector> lambdafull = FluidField().Interface()->InsertFSICondVector(lambda_);
	    const Teuchos::ParameterList& fsidyn   = DRT::Problem::Instance()->FSIDynamicParams();
	    const int uprestart = fsidyn.get<int>("RESTARTEVRY");
	    const int upres = fsidyn.get<int>("UPRES");
	    if ((uprestart != 0 && FluidField().Step() % uprestart == 0) || FluidField().Step() % upres == 0)
	      FluidField().DiscWriter()->WriteVector("fsilambda", lambdafull);
	  }


	  AleField().Output();
	  FluidField().LiftDrag();

	  if (StructureField()->GetConstraintManager()->HaveMonitor())
	  {
	    StructureField()->GetConstraintManager()->ComputeMonitorValues(StructureField()->Dispnp());
	    if (Comm().MyPID() == 0)
	      StructureField()->GetConstraintManager()->PrintMonitorValues();
	  }

	return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::FluidFluidMonolithicFluidSplitNoNOX::SetupNewSystem()
{
	TEUCHOS_FUNC_TIME_MONITOR("FSI::MonolithicFluidSplit::SetupNewSystem()");

	/*----------------------------------------------------------------------
	Create a combined map for Structure/Fluid/ALE-DOFs all in one!
	----------------------------------------------------------------------*/
	//RC-Pointers to the single maps are stored in vecSpaces
	std::vector<Teuchos::RCP<const Epetra_Map> > vecSpaces;

	//Append the structural DOF map!
	vecSpaces.push_back(StructureField()->DofRowMap());

	//Create the background fluid maps and the embedded fluid maps merged, without
	//fluid FSI-DOFs!
	std::vector<Teuchos::RCP<const Epetra_Map> > vecSpacesFluid;

	//Therefore:
	/*----------------------------------------------------------------------*/
	//Background fluid DOFs
	vecSpacesFluid.push_back(FluidField().XFluidFluidMapExtractor()->XFluidMap());
	//Embedded fluid DOFs without FSI DOFs
	vecSpacesFluid.push_back(FluidField().Interface()->OtherMap());
	//Merged FSI-free fluid maps
	Teuchos::RCP<const Epetra_Map> fluidmaps=Teuchos::rcp(new Epetra_Map(*LINALG::MultiMapExtractor::MergeMaps(vecSpacesFluid)));
	//Append the final fluid DOF map!
	vecSpaces.push_back(fluidmaps);
	/*----------------------------------------------------------------------*/

	//Append ALE DOF map
	vecSpaces.push_back(AleField().Interface()->OtherMap());

	//If the non-FSI fluid maps are empty
	if (vecSpaces[1]->NumGlobalElements() == 0)
		dserror("No inner fluid equations. Can't split!");

	//The vector is complete, fill the system's global BlockRowMap
	//with the maps previously set together!
	SetDofRowMaps(vecSpaces);

	//Initializing the system matrix equivalent to SetupSystem!
	systemmatrix_= Teuchos::rcp(
					        new LINALG::BlockSparseMatrix<LINALG::DefaultBlockMatrixStrategy>(
							Extractor(),
							Extractor(),
							81,
							false,
							true
							)
				        );

}
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::FluidFluidMonolithicFluidSplitNoNOX::Newton()
{
	TEUCHOS_FUNC_TIME_MONITOR("FSI::MonolithicFluidSplit::Newton()");

	/*----------------------------------------------------------------------
	Initialization
	----------------------------------------------------------------------*/
	//Iteration counter
	iter_ = 1;

	//Create a new Epetra_Vector with a given row map
	//and initialize it!
	//DofRowMap() is a method inherited from MonolithicNoNOX base class
	//it returns blockrowdofmap_, a MapExtractor object for
	//a DOF map, splitted into blocks.

	//Increment sum vector
	x_sum_ = LINALG::CreateVector(*DofRowMap(), true);

	//Solution vector
	iterinc_ = LINALG::CreateVector(*DofRowMap(), true);

	//Residual vector
	rhs_ = LINALG::CreateVector(*DofRowMap(), true);

	zeros_ = LINALG::CreateVector(*DofRowMap(), true);

	//Flag for special treatment of RHS-setup for i==0
	firstcall_ = true;

	/*----------------------------------------------------------------------
	Extract predictor increments
	----------------------------------------------------------------------*/
	//Increment of structural interface displacement --> structural predictor!!
	ddgpred_=Teuchos::rcp(new Epetra_Vector(*StructureField()->ExtractInterfaceDispnp()));
	ddgpred_->Update(-1.0,*StructureField()->ExtractInterfaceDispn(),1.0);

	/*----------------------------------------------------------------------*/
	//Initialize the increment vectors, they are updated in Evaluate(...)->ExtractFieldVectors(...)
	//at every Newton iteration!

	//Initialization for 1st Newton call
	//structural interface predictor
	ddginc_    = Teuchos::rcp(new Epetra_Vector(*ddgpred_));
	ddialeinc_ = Teuchos::rcp(new Epetra_Vector(*AleField().Interface()->OtherMap()),true);
	duiinc_	   = Teuchos::rcp(new Epetra_Vector(*Extractor().ExtractVector(iterinc_,1)),true);


	//We want to make sure, that the loop is entered at least once!
	//We exit the loop if either the convergence criteria are met (checked in
	//Converged() method) OR the maximum number of iterations is exceeded!
	while ( (iter_ == 1) or ( (not Converged()) and (iter_ <= itermax_) ) )
	{
		//Evaluate()- call
		//This function calls evaluate methods for all fields,
		//assembles and transfers current ale mesh positions to the fluid field.
		Evaluate(iterinc_);

		//Did the fluid map change after increment?
		//Problem here: The fluid FSI DOFs are removed from the system
		//they are not in FSI algorithm's global DOF-map anymore -
		//the increment vector map is free of FSI DOFs.
		//But as they are still available in the internal fluid field DOF map,
		//a comparison of the 2 maps would always return 'false'.
		//Solution:
		//We need a splitter object to remove fluid FSI DOFs from the fluid DOF map
		//before the maps are compared!
		//We want to control, if the fluid map itself has really changed!

		Teuchos::RCP<Epetra_Map> fsimap=Teuchos::rcp(new Epetra_Map(*FluidField().Interface()->FSICondMap()));
		Teuchos::RCP<LINALG::MapExtractor> ffsextractor=Teuchos::rcp(new LINALG::MapExtractor(*(FluidField().DofRowMap()),fsimap,true));
		Teuchos::RCP<const Epetra_Map> innerfluidmap=ffsextractor->OtherMap();

		bool isnewfluidmap=(innerfluidmap->SameAs(Extractor().ExtractVector(iterinc_,1)->Map()));

		if ( not isnewfluidmap )
		{
			IO::cout << "New Map!" << IO::endl;

			//Save old sum of increments
			Teuchos::RCP<Epetra_Vector> x_sum_n=LINALG::CreateVector(*DofRowMap(), true);
			*x_sum_n = *x_sum_;
			//Extract structural increment sum
			Teuchos::RCP<const Epetra_Vector> sx_n;
			sx_n = Extractor().ExtractVector(x_sum_n,0);
			//Extract ALE increment sum
			Teuchos::RCP<const Epetra_Vector> ax_n;
			ax_n = Extractor().ExtractVector(x_sum_n,2);

			SetupNewSystem();
			//Get the new XFluidFluid-MapExtractor object
			xfluidfluidsplitter_ = FluidField().XFluidFluidMapExtractor();
			iterinc_ = LINALG::CreateVector(*DofRowMap(),true);
			rhs_ = LINALG::CreateVector(*DofRowMap(),true);
			zeros_ = LINALG::CreateVector(*DofRowMap(),true);
			x_sum_ = LINALG::CreateVector(*DofRowMap(),true);

			/*----------------------------------------------------------------------
			Set the new increment sum x_sum_ together
			----------------------------------------------------------------------*/

			Extractor().InsertVector(sx_n,0,x_sum_);

			//In FluidField().Stepinc() there are still FSI DOFs that need to be removed...
			Teuchos::RCP<Epetra_Map> fsimap = Teuchos::rcp(new Epetra_Map(*FluidField().Interface()->FSICondMap()));
			Teuchos::RCP<LINALG::MapExtractor> ffsextractor = Teuchos::rcp(new LINALG::MapExtractor(*(FluidField().DofRowMap()),fsimap,true));
			Teuchos::RCP<Epetra_Vector> ff_stepinc = ffsextractor->ExtractOtherVector(FluidField().Stepinc());
			Extractor().InsertVector(ff_stepinc,1,x_sum_);

			Extractor().InsertVector(ax_n,2,x_sum_);

			//The fluid length may have changed:
			nf_ = FluidField().RHS()->GlobalLength();
		}

		/*----------------------------------------------------------------------*/
		//Build the linear system
		//J^{n+1}(x_i) \Delta^{n+1}_{x_i+1}=-r^{n+1}(x_i)
		//i: Newton iteration counter
		//J: Jacobian
		//r: RHS-vector
		/*----------------------------------------------------------------------*/
		SetupSystemMatrix();

		if ( not systemmatrix_->Filled() )
			dserror("Unfilled system matrix! Fatal error!");

		//Create the RHS consisting of the field residuals, the Lagrange
		//multiplier and predictor-related terms added at the first Newton
		//iteration only.
		SetupRHS(*rhs_,firstcall_);

		//Solver call
		LinearSolve();

		//Adapt solver tolerance (important if Aztec solver is picked)
		solver_->ResetTolerance();

		//Build residual and incremental norms,
		//count the DOFs!
		BuildCovergenceNorms();

		//Give some output
		PrintNewtonIter();

		//Increment loop index
		iter_+=1;

		firstcall_ = false;
	}//End of Newton loop!

	//After the loop exit, the iteration counter is 1 higher than the true no. of
	//iterations! Correct that:
	iter_-=1;

	/*----------------------------------------------------------------------*/
	//Compute the increments needed for recovery of Lagrange Multiplier!
	//After the last Newton iteration, the increments are not updated.
	//We need the last increment for the recovery of lambda.
	/*----------------------------------------------------------------------*/
	//Fluid
	duiinc_->Update(1.0,*Extractor().ExtractVector(iterinc_,1),0.0);
	//Structure
	Teuchos::RCP<Epetra_Vector> ddinc=Extractor().ExtractVector(iterinc_,0);
	ddginc_->Update(1.0,*StructureField()->Interface()->ExtractFSICondVector(ddinc),0.0);
	//ALE
	ddialeinc_->Update(1.0,*Extractor().ExtractVector(iterinc_,2),0.0);

	if ( Converged() and Comm().MyPID() == 0 )
	{
	    IO::cout << "-------------------------------------------------------------------------------------------"
					"-------------------------------------------------------------------------"<< IO::endl;
	    IO::cout << "  Newton Converged! " << IO::endl;
	}
	else if ( iter_ >= itermax_ )
	{
		IO::cout << IO::endl;
	    IO::cout << " Newton unconverged in "<< iter_ << " iterations " <<  IO::endl;
	}
}

/*----------------------------------------------------------------------
 * BuildConvergenceNorms:
 *     - Calculate the residual and incremental norms required for
 *     	 the convergence test in Newton-loop
 *     - Implemented:
 *     		The (Euclidian) L2-Norm
 ----------------------------------------------------------------------*/
void FSI::FluidFluidMonolithicFluidSplitNoNOX::BuildCovergenceNorms()
{
	TEUCHOS_FUNC_TIME_MONITOR("FSI::MonolithicFluidSplit::BuildConvergenceNorms");

	/*----------------------------------------------------------------------
	Residual norms - L2 of:
	- global rhs
	- inner structural rhs
	- ale rhs
	- inner fluid-fluid velocity rhs
	- inner fluid-fluid pressure rhs
	- complete interface residual
	----------------------------------------------------------------------*/
	//Norm of global RHS vector
 	rhs_->Norm2(&normrhs_);

	//Inner structural RHS and interface RHS

	//RHS-vector from FluidField() without FSI-DOFs
	Teuchos::RCP<const Epetra_Vector> innerfluidfluidrhs=Extractor().ExtractVector(rhs_,1);
	//(Inner) ALE RHS
	Teuchos::RCP<const Epetra_Vector> alerhs=Extractor().ExtractVector(rhs_,2);

	//Norm of inner structural residual forces
	Teuchos::RCP<const Epetra_Vector> structrhs=Extractor().ExtractVector(rhs_,0);
	StructureField()->Interface()->ExtractOtherVector(structrhs)->Norm2(&normstrrhsL2_);
	StructureField()->Interface()->ExtractOtherVector(structrhs)->NormInf(&normstrrhsInf_);

	//Norm of ALE residual forces
	alerhs->Norm2(&normalerhsL2_);

	//Norm of fluid velocity residual
	//This requires an Epetra_Map of the inner fluid velocity DOFs first!
	Teuchos::RCP<const Epetra_Map>  innerfluidvel=FluidField().InnerVelocityRowMap();
	//The background fluid maps and the embedded fluid maps without
	//fluid FSI DOFs
	std::vector<Teuchos::RCP<const Epetra_Map> > vecSpacesFluid;
	//Background fluid DOFs
	vecSpacesFluid.push_back(FluidField().XFluidFluidMapExtractor()->XFluidMap());
	//Embedded fluid DOFs without FSI DOFs
	vecSpacesFluid.push_back(FluidField().Interface()->OtherMap());
	//Merged FSI-free fluid maps
	Teuchos::RCP<const Epetra_Map> fluidmaps=Teuchos::rcp(new Epetra_Map(*LINALG::MultiMapExtractor::MergeMaps(vecSpacesFluid)));
	//Create a MapExtractor to access the velocity DOFs from the FSI-free fluid map
	Teuchos::RCP<LINALG::MapExtractor> fluidvelextract=Teuchos::rcp(new LINALG::MapExtractor(*fluidmaps,innerfluidvel,true));

	//Finally, compute the fluid velocity RHS-norm
	fluidvelextract->ExtractCondVector(innerfluidfluidrhs)->Norm2(&normflvelrhsL2_);
	fluidvelextract->ExtractCondVector(innerfluidfluidrhs)->NormInf(&normflvelrhsInf_);

	//Norm of fluid pressure residual
	//This requires an Epetra_Map of the fluid pressure DOFs
	if ( FluidField().PressureRowMap() == Teuchos::null )
		dserror("Empty pressure row map!");

	//Finally, compute the fluid pressure RHS-norm
	fluidvelextract->ExtractOtherVector(innerfluidfluidrhs)->Norm2(&normflpresrhsL2_);
	fluidvelextract->ExtractOtherVector(innerfluidfluidrhs)->NormInf(&normflpresrhsInf_);

	//The true RHS for the FSI interface equation block consists of
	//more than just the structure residual, namely the scaled fluid interface residual and the
	//previous Lagrange multiplier. The first idea is, to test this whole term, which can be easily
	//extracted from rhs_. For a more strict testing, the L_inf-norm should be employed!
	Teuchos::RCP<Epetra_Vector> interfaceresidual=StructureField()->Interface()->ExtractFSICondVector(*structrhs);
	interfaceresidual->Norm2(&norminterfacerhsL2_);
	interfaceresidual->NormInf(&norminterfacerhsInf_);

	/*----------------------------------------------------------------------
	Incremental norms - L2
	- global inc
	- inner structural inc
    - complete interface inc
	- inner fluid-fluid velocity rhs
	- inner fluid-fluid pressure rhs
	----------------------------------------------------------------------*/
	//Norm of global increment vector
	iterinc_->Norm2(&norminc_);

	Teuchos::RCP<const Epetra_Vector> structinc=Extractor().ExtractVector(iterinc_,0);
	Teuchos::RCP<const Epetra_Vector> fluidinc=Extractor().ExtractVector(iterinc_,1);

	//Norm of inner structural increment vector
	StructureField()->Interface()->ExtractOtherVector(structinc)->Norm2(&normstrincL2_);
	StructureField()->Interface()->ExtractOtherVector(structinc)->NormInf(&normstrincInf_);

	//Norm of interface increment vector
	StructureField()->Interface()->ExtractFSICondVector(structinc)->Norm2(&norminterfaceincL2_);
	StructureField()->Interface()->ExtractFSICondVector(structinc)->NormInf(&norminterfaceincInf_);

	//Norm of fluid velocity increment
	fluidvelextract->ExtractCondVector(fluidinc)->Norm2(&normflvelincL2_);
	fluidvelextract->ExtractCondVector(fluidinc)->NormInf(&normflvelincInf_);

	//Norm of fluid pressure increment
	fluidvelextract->ExtractOtherVector(fluidinc)->Norm2(&normflpresincL2_);
	fluidvelextract->ExtractOtherVector(fluidinc)->NormInf(&normflpresincInf_);

	//Norm of ALE increment vector
	Extractor().ExtractVector(iterinc_,2)->Norm2(&normaleincL2_);

	//get length of the structural, fluid and ale vector
	ni_=(StructureField()->Interface()->ExtractFSICondVector(*structrhs))->GlobalLength();
	ns_=(StructureField()->Interface()->ExtractOtherVector(structrhs))->GlobalLength();
	nf_=innerfluidfluidrhs->GlobalLength();
	nfv_=fluidvelextract->ExtractCondVector(fluidinc)->GlobalLength();
	nfp_=fluidvelextract->ExtractOtherVector(fluidinc)->GlobalLength();
	na_=alerhs->GlobalLength();
	nall_=rhs_->GlobalLength();
}
/*----------------------------------------------------------------------
 * RecoverLagrangeMultiplier:
 *     - Compute the Lagrange multiplier at the FSI-interface
 ----------------------------------------------------------------------*/
 void FSI::FluidFluidMonolithicFluidSplitNoNOX::RecoverLagrangeMultiplier()
{
 	/*----------------------------------------------------------------------*/
	//Time integration parameter
	/*----------------------------------------------------------------------*/
	//Fluid
	//b*y_n+(1-b)*y_n+1
	//b:flintparam
	/*----------------------------------------------------------------------*/
	 const double ftimintparam=FluidField().TimIntParam();

	//Fluid time scaling parameter
	//fluidtimescale: \tau
	// \tau = 1/\Delta t  for Backward-Euler;
	// \tau = 2/\Delta t  for Trapezoidal Rule
	 const double fluidtimescale=FluidField().TimeScaling();

	 //Scaling factor for different fluid/structural units
	 const double fluidresidualscale=FluidField().ResidualScaling();

	/*----------------------------------------------------------------------
	  Lagrange Multiplier Setup
	------------------------------------------------------------------------
	 The Langrange Multiplier is updated as follows:

	 lambda_^{n+1}=

	 		1/(1-flintparam)*(

	 (1)	-flintparam*lambda_^n

	 (2)	-r_{\Gamma}^{F,n+1}

	 (3)	- 1/Tau *(F_{\Gamma\Gamma)})* \Delta d_{\Gamma}{S,n+1}

	 (4)	-F^{G}_{\Gamma\Gamma} * \Delta d_{\Gamma}{S,n+1}

	 (5)	-F_{\Gamma I} * \Delta u_{I}^{F,n+1}

	 (6)	-F^{G}_{\Gamma I} * \Delta d_{I,n+1}^{G,n+1}

			Only at first Newton iteration:
	 (7)	+dt / Tau * F_{\Gamma\Gamma} * u_{\Gamma}^n}
	        )
	  + tau: time scaling factor for interface time integration (tau = 1/FluidField().TimeScaling())
	  + the terms 2 to 7 are first saved in a tmpvec which will be added to lambda
	 ----------------------------------------------------------------------*/

	 //creating & initializing the storage vectors for the last four terms
	 Teuchos::RCP<Epetra_Vector> fggddg=Teuchos::null;
	 Teuchos::RCP<Epetra_Vector> fmggddg=Teuchos::null;
	 Teuchos::RCP<Epetra_Vector> fgidui=Teuchos::null;
	 Teuchos::RCP<Epetra_Vector> fmgiddia=Teuchos::null;

	 // stores intermediate result of terms (3)-(7)
	 Teuchos::RCP<Epetra_Vector> tmpvec = Teuchos::null;

	 // ---------Addressing term (2)
	 //store f^F_{\Gamma}! As Recover-LM is called after the Newton loop, the RHS will have changed!
	 Teuchos::RCP<Epetra_Vector> femb=xfluidfluidsplitter_->ExtractFluidVector(FluidField().RHS());
	 Teuchos::RCP<Epetra_Vector> fluidresidual = FluidField().Interface()->ExtractFSICondVector(femb);

	 // ---------Addressing term (1)
	 lambda_->Update(ftimintparam,*lambda_,0.0);

	 // ---------Addressing term (2)
	 tmpvec = Teuchos::rcp(new Epetra_Vector(*fluidresidual));
	 tmpvec->Scale(-1.0);


     // ---------Addressing term (3)
	 if (fggcurr_ != Teuchos::null)
	 {
		 fggddg=Teuchos::rcp(new Epetra_Vector(fggcurr_->RangeMap(),true));
		 fggcurr_->Apply(*StructToFluid(ddginc_),*fggddg);
		 tmpvec->Update(fluidtimescale,*fggddg,1.0);
	 }

	 //(4)
	 if (fmggcurr_ != Teuchos::null)
	 {
		 Teuchos::RCP<Epetra_Vector> fmggddg=Teuchos::rcp(new Epetra_Vector(fmggcurr_->RangeMap(),true));
		 fmggcurr_->Apply(*StructToFluid(ddginc_),*fmggddg);
		 tmpvec->Update(1.0,*fmggddg,1.0);
	 }

	 //(5)
	 if (fgicurr_ != Teuchos::null)
	 {
		 Teuchos::RCP<Epetra_Vector> fgidui=Teuchos::rcp(new Epetra_Vector(fgicurr_->RangeMap(),true));
		 fgicurr_->Apply(*duiinc_,*fgidui);
		 tmpvec->Update(1.0,*fgidui,1.0);
	 }

	 //(6)
	 if (fmgicurr_ != Teuchos::null)
	 {

		 //The DomainMap() of matrix fmgipre_ contains inner velocity and pressure DOFs!
		 //AleToFluid converts the inner ALE displacement increments to inner Fluid velocity DOFs.
		 //The underlying map of this vector has to match the DomainMap()!
		 //Hence, the missing pressure DOFs have to be appended.
		 Teuchos::RCP<Epetra_Vector> fmgiddia=Teuchos::rcp(new Epetra_Vector(fmgicurr_->RangeMap(),true));

		 //std::vector<Teuchos::RCP<const Epetra_Map> > fluidpresmaps;
		 //Global fluid pressure DOF map
		 //fluidpresmaps.push_back(FluidField().PressureRowMap());
		 //Inner FSI free pressure DOF map
		 //fluidpresmaps.push_back(FluidField().Interface()->OtherMap());

		 //Merged FSI-free fluid pressure DOF map
		 //Teuchos::RCP<const Epetra_Map> innerfluidpresmap=Teuchos::rcp(new Epetra_Map(*LINALG::MultiMapExtractor::IntersectMaps(fluidpresmaps)));

		 Teuchos::RCP<LINALG::MapExtractor> fluidpresextractor=Teuchos::rcp(new LINALG::MapExtractor(*xfluidfluidsplitter_->FluidMap(),FluidField().PressureRowMap(),false));
		 Teuchos::RCP<Epetra_Vector> tmp=Teuchos::rcp(new Epetra_Vector(fmgicurr_->DomainMap(),true));
		 Teuchos::RCP<Epetra_Vector> aux=AleToFluid(AleField().Interface()->InsertOtherVector(ddialeinc_));
		 // vector with pressure dofs
		 Teuchos::RCP<Epetra_Vector> auxaux=fluidpresextractor->InsertCondVector(aux);
		 tmp=FluidField().Interface()->ExtractOtherVector(auxaux);
		 fmgicurr_->Apply(*tmp,*fmgiddia);
		 tmpvec->Update(1.0,*fmgiddia,1.0);
	 }

	 //(7)
	 if (firstcall_)
	 {
		 if (fggcurr_ != Teuchos::null)
		 {
			 Teuchos::RCP<Epetra_Vector> tmp=Teuchos::rcp(new Epetra_Vector(fggcurr_->RangeMap(),true));
			 Teuchos::RCP<Epetra_Vector> fveln=FluidField().ExtractInterfaceVeln();
			 fggcurr_->Apply(*fveln,*tmp);
			 tmpvec->Update(Dt()*fluidtimescale,*tmp,1.0);
		 }
	 }


	// ---------Adding tmpvec to lambda_
	lambda_->Update(fluidresidualscale,*tmpvec,1.0); // scale with ResidualScaling() to get [N/m^2]

	//Scaling everything with -1/(1-flintparam_)
    lambda_->Scale(-1.0/(1.0-ftimintparam));

    return;


}

