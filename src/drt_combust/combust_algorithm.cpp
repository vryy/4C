/*!-----------------------------------------------------------------------------------------------*
\file combust_algorithm.cpp

\brief base combustion algorithm

    detailed description in header file combust_algorithm.H

\level 2

<pre>
\maintainer Benedikt Schott
            schott@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15241
</pre>

 *------------------------------------------------------------------------------------------------*/

#include "combust_algorithm.H"
#include "combust_defines.H"
#include "combust_flamefront.H"
#include "two_phase_defines.H"
#include "combust_fluidimplicitintegration.H"
#include "../drt_fluid_turbulence/turbulence_statistic_manager.H"
#include "../drt_fluid_turbulence/turbulence_statistics_mean_general.H"
#include "../drt_lib/drt_periodicbc.H"
#include "../drt_lib/drt_dofset.H"
#include "../drt_mat/newtonianfluid.H"
#include "../drt_mat/matlist.H"
#include "../drt_io/io_gmsh.H"
#include "../drt_io/io_pstream.H"
#include "../drt_io/io_control.H"
#include "../linalg/linalg_utils.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_levelset/levelset_algorithm.H"
#include "../drt_scatra/scatra_timint_ost.H"


/*------------------------------------------------------------------------------------------------*
 | constructor                                                                        henke 06/08 |
 *------------------------------------------------------------------------------------------------*/
COMBUST::Algorithm::Algorithm(const Epetra_Comm& comm, const Teuchos::ParameterList& combustdyn, const Teuchos::ParameterList& solverparams)
: ScaTraFluidCouplingAlgorithm(comm, combustdyn,false,"scatra",solverparams),
// initialize member variables
  fgiter_(0),
  fgitermax_(combustdyn.get<int>("ITEMAX")),
  convtol_(combustdyn.get<double>("CONVTOL")),
  combusttype_(DRT::INPUT::IntegralValue<INPAR::COMBUST::CombustionType>(combustdyn.sublist("COMBUSTION FLUID"),"COMBUSTTYPE")),
  evaltimeratio_(1.0),
  gmshoutput_(combustdyn.get<bool>("GMSH_OUTPUT")),
  combustdyn_(combustdyn),
  flamefront_(Teuchos::null),
  restart_(false),
  gen_flow_(DRT::INPUT::IntegralValue<int>(combustdyn,"GENERATE_FLOW_FIELD"))
{
  // Keep constructor empty
  return;
}

/*------------------------------------------------------------------------------------------------*
 | destructor                                                                         henke 06/08 |
 *------------------------------------------------------------------------------------------------*/
COMBUST::Algorithm::~Algorithm()
{
}

/*------------------------------------------------------------------------------------------------*
 | init                                                                               rauch 09/16 |
 *------------------------------------------------------------------------------------------------*/
void COMBUST::Algorithm::Init(
    const Teuchos::ParameterList&   prbdyn,         ///< parameter list for global problem
    const Teuchos::ParameterList&   scatradyn,      ///< parameter list for scalar transport subproblem
    const Teuchos::ParameterList&   solverparams,   ///< parameter list for scalar transport solver
    const std::string&              disname,        ///< name of scalar transport discretization
    const bool                      isale           ///< ALE flag
)
{
  // call Init() in base class
  ADAPTER::ScaTraFluidCouplingAlgorithm::Init(
      prbdyn,
      scatradyn,
      solverparams,
      disname,
      isale);

  if (Comm().MyPID()==0)
  {
    switch(combusttype_)
    {
    case INPAR::COMBUST::combusttype_premixedcombustion:
      IO::cout << "COMBUST::Algorithm: this is a premixed combustion problem" << IO::endl;
      break;
    case INPAR::COMBUST::combusttype_twophaseflow:
      IO::cout << "COMBUST::Algorithm: this is a two-phase flow problem" << IO::endl;
      break;
    case INPAR::COMBUST::combusttype_twophaseflow_surf:
      IO::cout << "COMBUST::Algorithm: this is a two-phase flow problem with kinks in vel and jumps in pres" << IO::endl;
      break;
    case INPAR::COMBUST::combusttype_twophaseflowjump:
      IO::cout << "COMBUST::Algorithm: this is a two-phase flow problem with jumps in vel and pres" << IO::endl;
      break;
    default:
      dserror("unknown type of combustion problem");
      break;
    }
  }

  return;
}

/*------------------------------------------------------------------------------------------------*
 | setup                                                                              rauch 09/16 |
 *------------------------------------------------------------------------------------------------*/
void COMBUST::Algorithm::Setup()
{
  // call Setup() in base class
  ADAPTER::ScaTraFluidCouplingAlgorithm::Setup();

  // get pointers to the discretizations from the time integration scheme of each field
  // remark: fluiddis cannot be of type "const Teuchos::RCP<const DRT::Dis...>", because parent
  // class. InterfaceHandle only accepts "const Teuchos::RCP<DRT::Dis...>"              henke 01/09
  const Teuchos::RCP<DRT::Discretization> fluiddis = FluidField()->Discretization();
  const Teuchos::RCP<DRT::Discretization> gfuncdis = ScaTraField()->Discretization();

  // FGI vectors are initialized
  velnpi_ = Teuchos::rcp(new Epetra_Vector(FluidField()->StdVelnp()->Map()),true);//*fluiddis->DofRowMap()),true);
  velnpi_->Update(1.0,*FluidField()->StdVelnp(),0.0);
  phinpi_ = Teuchos::rcp(new Epetra_Vector(*gfuncdis->DofRowMap()),true);
  phinpi_->Update(1.0,*ScaTraField()->Phinp(),0.0);

  /*----------------------------------------------------------------------------------------------*
   * initialize all data structures needed for the combustion algorithm
   *
   * - capture the flame front and create interface geometry (triangulation)
   * - determine initial enrichment (DofManager wird bereits mit dem Element d.h. Diskretisierung angelegt)
   * - ...
   *----------------------------------------------------------------------------------------------*/
  // construct initial flame front
  flamefront_ = Teuchos::rcp(new COMBUST::FlameFront(fluiddis,gfuncdis));
  flamefront_->UpdateFlameFront(combustdyn_,ScaTraField()->Phin(), ScaTraField()->Phinp());
  flamefront_->UpdateOldInterfaceHandle();

  //---------------------------------------------------
  // set initial fluid field (based on standard dofset)
  //---------------------------------------------------
  // update fluid interface with flamefront
  Teuchos::rcp_dynamic_cast<FLD::CombustFluidImplicitTimeInt>(FluidField())->ImportFlameFront(flamefront_,false);
  // read parameters for initial field
  const INPAR::COMBUST::InitialField initfield = DRT::INPUT::IntegralValue<INPAR::COMBUST::InitialField>(
      combustdyn_.sublist("COMBUSTION FLUID"),"INITIALFIELD");
  const int initfuncno = combustdyn_.sublist("COMBUSTION FLUID").get<int>("INITFUNCNO");
  // set initial flow field based on standard dofs only
  Teuchos::rcp_dynamic_cast<FLD::CombustFluidImplicitTimeInt>(FluidField())->SetInitialFlowField(initfield, initfuncno);
  // clear fluid's memory to flamefront
  Teuchos::rcp_dynamic_cast<FLD::CombustFluidImplicitTimeInt>(FluidField())->ImportFlameFront(Teuchos::null,false);

  //--------------------------------------------------------
  // incorporate the XFEM dofs into the fluid
  // remark: this includes setting the initial enriched dofs
  //--------------------------------------------------------
  // update fluid interface with flamefront
  Teuchos::rcp_dynamic_cast<FLD::CombustFluidImplicitTimeInt>(FluidField())->ImportFlameFront(flamefront_,true);
  // clear fluid's memory to flamefront
  Teuchos::rcp_dynamic_cast<FLD::CombustFluidImplicitTimeInt>(FluidField())->ImportFlameFront(Teuchos::null,true);

  // access the scatra discretization
  Teuchos::RCP<DRT::Discretization> scatradis = DRT::Problem::Instance()->GetDis("scatra");

  // add proxy of fluid transport degrees of freedom to scatra discretization
  if(scatradis->AddDofSet(Teuchos::rcp_const_cast<DRT::DofSet>(FluidField()->DofSet())) != 1)
    dserror("Scatra discretization has illegal number of dofsets!");

  // transfer the initial convective velocity from initial fluid field to scalar transport field
  // subgrid scales not transferred since they are zero at time t=0.0
  // this step has already been done in the ScaTraFluidCouplingAlgorithm(), however, for problems
  // with particles we have to set the old convective velocity
  // moreover, for combustion problems, we have to redo this step anyway, since we have to consider
  // the flame velocity here, which can not be supported by the base class
  // bool indicates initalization call
  SetVelocityLevelSet(true);

  return;
}


void COMBUST::Algorithm::DoAlgorithmSpecificInit()
{
  // check whether fluid and scatra discret still have the same maps
  // they may change due a modified ghosting required, i.e., for particle level-set methods
  if(Teuchos::rcp_dynamic_cast<SCATRA::LevelSetAlgorithm>(ScaTraField()) != Teuchos::null)
  {
    const Epetra_Map* scatraelecolmap = ScaTraField()->Discretization()->ElementColMap();
    if (not scatraelecolmap->PointSameAs(*FluidField()->Discretization()->ElementColMap()))
    {
      if (Comm().MyPID()==0)
        std::cout << "----- adaption of fluid ghosting to scatra ghosting ------" << std::endl;

      // adapt fluid ghosting to scatra ghosting
      FluidField()->Discretization()->ExtendedGhosting(*scatraelecolmap,false,false,false,false);
    }
  }
}


/*------------------------------------------------------------------------------------------------*
 | public: algorithm for a dynamic combustion problem                                 henke 06/08 |
 *------------------------------------------------------------------------------------------------*/
void COMBUST::Algorithm::TimeLoop()
{
  // safety check
  CheckIsInit();
  CheckIsSetup();

  // output of initial fields
  OutputInitialField();

  // get initial field by solving stationary problem first
  // however, calculate it if and only if the problem has not been restarted
  if(DRT::INPUT::IntegralValue<int>(combustdyn_.sublist("COMBUSTION FLUID"),"INITSTATSOL") == true and restart_==false)
    SolveInitialStationaryProblem();

  // time loop
  while (NotFinished())
  {
    // redistribute
    Redistribute();

    // prepare next time step; update field vectors
    PrepareTimeStep();

    // Fluid-G-function-Interaction loop
    // remark: - the G-function is solved first, then the fluid is solved
    //         - this guarantees that the interface geometry and the fluid solution match
    //         - this is a prerequisite for the semi-Lagrangean time integration method
    //         - changing this order would affect the restart
    //         -> do not change this order!
    while (NotConvergedFGI())
    {
      // prepare Fluid-G-function iteration
      PrepareFGIteration();

      // solve linear G-function equation
      if (not gen_flow_)
        DoGfuncField();
      else
      {
        // do not solve the level-field if flow generation for turbulent problems
        // is performed
        if (Comm().MyPID()==0)
          IO::cout << "Level set deactivated!" << IO::endl;
      }

      // update interface geometry
      // note: do not call this function at any other place
      UpdateInterface();

      // solve nonlinear Navier-Stokes system
      DoFluidField();

    } // Fluid-G-function-Interaction loop

    // solve for level-set field once more
    if (not gen_flow_)
      DoGfuncField();

    // update all field solvers
    UpdateTimeStep();
    //Remark (important for restart): the time level of phi (n+1, n or n-1) used to reconstruct the interface
    //                                conforming to the restart state of the fluid depends on the order
    //                                of Output() and UpdateTimeStep()

    // write output to screen and files
    Output();

  } // time loop

  return;
}

/*------------------------------------------------------------------------------------------------*
 | public: algorithm for a stationary combustion problem                              henke 10/09 |
 *------------------------------------------------------------------------------------------------*/
void COMBUST::Algorithm::SolveStationaryProblem()
{
  if (Comm().MyPID()==0)
  {
    printf("--------Stationary-Combustion-------  time step ----------------------------------------\n");
  }

  // check if ScaTraField()->initialvelset == true
  /* remark: initial velocity field has been transferred to scalar transport field in constructor of
   * ScaTraFluidCouplingAlgorithm (initialvelset_ == true). Time integration schemes, such as
   * the one-step-theta scheme, are thus initialized correctly.
   */

  // check time integration schemes of single fields
  // remark: this was already done in ScaTraFluidCouplingAlgorithm() before
  if (FluidField()->TimIntScheme() != INPAR::FLUID::timeint_stationary)
    dserror("Fluid time integration scheme is not stationary");
  if (ScaTraField()->MethodName() != INPAR::SCATRA::timeint_stationary)
    dserror("Scatra time integration scheme is not stationary");

  // solve nonlinear Navier-Stokes system
  DoFluidField();

  // solve (non)linear G-function equation
  if (Comm().MyPID()==0)
    IO::cout << "/!\\ warning === G-function field not solved for stationary problems" << IO::endl;
  //DoGfuncField();

  // write output to screen and files (and Gmsh)
  Output();

  return;
}


/*------------------------------------------------------------------------------------------------*
 | protected: output of initial field                                             rasthofer 09/13 |
 *------------------------------------------------------------------------------------------------*/
void COMBUST::Algorithm::OutputInitialField()
{
  if (Step() == 0)
  {
    // update fluid interface with flamefront
    Teuchos::rcp_dynamic_cast<FLD::CombustFluidImplicitTimeInt>(FluidField())->ImportFlameFront(flamefront_,false);
    // output fluid initial state
    if (FluidField()->TimIntScheme() != INPAR::FLUID::timeint_stationary
        and DRT::INPUT::IntegralValue<int>(combustdyn_.sublist("COMBUSTION FLUID"),"INITSTATSOL") == false)
      FluidField()->Output();
    // clear fluid's memory to flamefront
    Teuchos::rcp_dynamic_cast<FLD::CombustFluidImplicitTimeInt>(FluidField())->ImportFlameFront(Teuchos::null,false);

    // output G-function initial state
    if (ScaTraField()->MethodName() != INPAR::SCATRA::timeint_stationary and
        DRT::INPUT::IntegralValue<int>(combustdyn_.sublist("COMBUSTION FLUID"),"INITSTATSOL") == false )
      ScaTraField()->Output();

    // write center of mass
    if (DRT::INPUT::IntegralValue<bool>(combustdyn_,"WRITE_CENTER_OF_MASS"))
      CenterOfMass();
  }

  return;
}


/*------------------------------------------------------------------------------------------------*
 | protected: compute flame velocity                                                  henke 07/09 |
 *------------------------------------------------------------------------------------------------*/
const Teuchos::RCP<Epetra_Vector> COMBUST::Algorithm::ComputeFlameVel(
    const Teuchos::RCP<const Epetra_Vector>& convel,
    const Teuchos::RCP<const DRT::DofSet>& dofset
    //const Teuchos::RCP<const Epetra_Map >& dbcmap
    )
{
  if (Comm().MyPID()==0)
  {
    if((DRT::INPUT::IntegralValue<int>(combustdyn_.sublist("COMBUSTION FLUID"),"INITSTATSOL") == false) and
       (DRT::INPUT::IntegralValue<INPAR::COMBUST::InitialField>(combustdyn_.sublist("COMBUSTION FLUID"),"INITIALFIELD") == INPAR::COMBUST::initfield_zero_field))
      IO::cout << "/!\\ Compute an initial stationary fluid solution to avoid a non-zero initial flame velocity" << IO::endl;
  }

  // temporary vector for convective velocity (based on dofrowmap of standard (non-XFEM) dofset)
  // remark: operations must not be performed on 'convel', because the vector is accessed by both
  //         master and slave nodes, if periodic bounday conditions are present
  Teuchos::RCP<Epetra_Vector> conveltmp = LINALG::CreateVector(*(dofset->DofRowMap()),true);

  // get a pointer to the fluid discretization
  const Teuchos::RCP<const DRT::Discretization> fluiddis = FluidField()->Discretization();

  // get G-function value vector on fluid NodeColMap
  const Teuchos::RCP<Epetra_Vector> phinp = flamefront_->Phinp();
  const Teuchos::RCP<Epetra_MultiVector> gradphinp = flamefront_->GradPhi();
  const Teuchos::RCP<Epetra_Vector> curvature = flamefront_->Curvature();

#ifdef DEBUG
  // get map of this vector
  const Epetra_BlockMap& phimap = phinp->Map();
  // check, whether this map is still identical with the current node map in the discretization
  if (not phimap.SameAs(*fluiddis->NodeColMap())) dserror("node column map has changed!");
  // get map of this vector
  const Epetra_BlockMap& gradphimap = gradphinp->Map();
  // check, whether this map is still identical with the current node map in the discretization
  if (not gradphimap.SameAs(*fluiddis->NodeColMap())) dserror("node column map has changed!");
#endif

  // laminar flame speed
  const double sl = combustdyn_.sublist("COMBUSTION FLUID").get<double>("LAMINAR_FLAMESPEED");
  const double marksteinlength  = combustdyn_.sublist("COMBUSTION FLUID").get<double>("MARKSTEIN_LENGTH");

#ifdef COMBUST_GMSH_NORMALFIELD
  const std::string filestr = "normal_field";
  const std::string name_in_gmsh = "Normal field";

  const std::string filename = IO::GMSH::GetNewFileNameAndDeleteOldFiles(filestr, Step(), 500, true, fluiddis->Comm().MyPID());
  std::ofstream gmshfilecontent(filename.c_str());
  {
    gmshfilecontent << "View \" " << name_in_gmsh << "\" {\n";
#endif

    // loop over nodes on this processor
    for(int lnodeid=0; lnodeid < fluiddis->NumMyRowNodes(); ++lnodeid)
    {
      //-------------------------------------------------------------
      // get (smoothed) gradient of the G-function field at this node
      //-------------------------------------------------------------
      LINALG::Matrix<3,1> mygradphi(true);

      // get the processor local node
      DRT::Node*  lnode = fluiddis->lRowNode(lnodeid);
      // get the global id for current node
      const int gid = lnode->Id();

#ifdef ORACLES
      // skip the part of the domain left of the expansion
      if (lnode->X()[0] < 0.0)
        continue;
#endif
      // get local processor id according to global node id
      const int nodelid = (*gradphinp).Map().LID(gid);
      if (nodelid<0) dserror("Proc %d: Cannot find gid=%d in Epetra_Vector",(*gradphinp).Comm().MyPID(),gid);

      const int numcol = (*gradphinp).NumVectors();
      if( numcol != 3) dserror("number of columns in Epetra_MultiVector is not 3");

      // loop over dimensions (= number of columns in multivector)
      for(int col=0; col< numcol; col++)
      {
        // get columns vector of multivector
        double* globalcolumn = (*gradphinp)[col];
        // set smoothed gradient entry of phi into column of global multivector
        mygradphi(col) = globalcolumn[nodelid];
      }

      //------------------------------------
      // compute smoothed unit normal vector
      // n = - grad phi / |grad phi|
      //------------------------------------
      // smoothed normal vector at this node
      LINALG::Matrix<3,1> nvec = mygradphi;
      // compute norm of smoothed normal vector
      const double ngradphi = mygradphi.Norm2();

      // divide vector by its norm to get unit normal vector
      if (fabs(ngradphi < 1.0E-12))// 'ngradnorm' == 0.0
      {
        // length of smoothed normal is zero at this node -> node must be on a singularity of the
        // level set function (e.g. "regular level set cone"); all normals add up to zero normal vector
        // -> The fluid convective velocity 'fluidvel' alone constitutes the flame velocity, since the
        //    relative flame velocity 'flvelrel' turns out to be zero due to the zero average normal vector.
        IO::cout << "\n/!\\ phi gradient too small at node " << gid << " -> flame velocity is only the convective velocity" << IO::endl;
        nvec.PutScalar(0.0);
      }
      else
      {
        nvec.Scale(-1.0/ngradphi);
      }

#ifdef COLLAPSE_FLAME
      nvec(0) = lnode->X()[0];
      nvec(1) = lnode->X()[1];
      nvec(2) = lnode->X()[2];
#ifdef COMBUST_2D
      nvec(2) = 0.0;
#endif
      const double norm = nvec.Norm2(); // sqrt(normal(0)*normal(0) + normal(1)*normal(1) + normal(2)*normal(2))
      if (norm == 0.0)
      {
        nvec.PutScalar(0.0);
        //dserror("norm of normal vector is zero!");
      }
      else
      {
        nvec.Scale(1.0/norm);
      }
#endif

#ifdef COMBUST_GMSH_NORMALFIELD
      if (gmshoutput_)
      {
        LINALG::SerialDenseMatrix xyz(3,1);
        xyz(0,0) = lnode->X()[0];
        xyz(1,0) = lnode->X()[1];
        xyz(2,0) = lnode->X()[2];

        IO::GMSH::cellWithVectorFieldToStream(DRT::Element::point1, nvec, xyz, gmshfilecontent);
      }
#endif

      //------------------------
      // get material parameters
      //------------------------
      // get list of adjacent elements of this node
      DRT::Element** elelist = lnode->Elements();

      // get material from first (arbitrary!) element adjacent to this node
      const Teuchos::RCP<MAT::Material> matlistptr = elelist[0]->Material();
      dsassert(matlistptr->MaterialType() == INPAR::MAT::m_matlist, "material is not of type m_matlist");
      const MAT::MatList* matlist = static_cast<const MAT::MatList*>(matlistptr.get());

      // density burnt domain
      Teuchos::RCP<const MAT::Material> matptrplus = matlist->MaterialById(3);
      dsassert(matptrplus->MaterialType() == INPAR::MAT::m_fluid, "material is not of type m_fluid");
      const MAT::NewtonianFluid* matplus = static_cast<const MAT::NewtonianFluid*>(matptrplus.get());
      const double rhoplus = matplus->Density();

      // density unburnt domain
      Teuchos::RCP<const MAT::Material> matptrminus = matlist->MaterialById(4);
      dsassert(matptrminus->MaterialType() == INPAR::MAT::m_fluid, "material is not of type m_fluid");
      const MAT::NewtonianFluid* matminus = static_cast<const MAT::NewtonianFluid*>(matptrminus.get());
      const double rhominus = matminus->Density();

      //---------------------------------------------
      // compute relative flame velocity at this node
      //---------------------------------------------
      // get phi value for this node
      const double gfuncval = (*phinp)[nodelid];
      const double curv = (*curvature)[nodelid];

      double speedfac = 0.0;
      double wallfac = 1.0;
#ifdef ORACLES
      //---------------------------------------------------------
      // blend the flame speed close to walls for ORACLES problem
      //---------------------------------------------------------
      //    wall
      // 1.0 |     ______
      //     |    /
      // 0.0 |___/ ,
      //     |     H/6
      const double wallzone = 0.0299/6.0;
      if (lnode->X()[0] > 0.0) // inside combustion chamber
      {
        if ( (0.0653-abs(lnode->X()[1])) < wallzone or // close to top or bottom wall
                         lnode->X()[0]   < wallzone )  // close to step
        {
          // wall factor is 0 at the wall and 1 at H/6 or further away from the wall
          wallfac = 6.0/0.0299 * std::min(0.0653-abs(lnode->X()[1]),lnode->X()[0]);
        }
      }
#endif

      if (XFEM::plusDomain(gfuncval) == true){ // interface or burnt domain -> burnt material
        // flame speed factor = laminar flame speed * rho_unburnt / rho_burnt
        speedfac = sl*(1.0-marksteinlength*curv)* rhominus/rhoplus;
      }
      else{ // unburnt domain -> unburnt material
        // flame speed factor = laminar flame speed
        speedfac = sl*(1.0-marksteinlength*curv);
      }

      LINALG::Matrix<3,1> flvelrel(true);
      for (int icomp=0; icomp<3; ++icomp)
        flvelrel(icomp) = wallfac * speedfac * nvec(icomp);

      //-----------------------------------------------
      // compute (absolute) flame velocity at this node
      //-----------------------------------------------
      LINALG::Matrix<3,1> fluidvel(true);
      // get the set of dof IDs for this node (3 x vel + 1 x pressure) from standard FEM dofset
      const std::vector<int> dofids = (*dofset).Dof(lnode);
      std::vector<int> lids(3);

      // extract velocity values (no pressure!) from global velocity vector
      for (int icomp=0; icomp<3; ++icomp)
      {
        lids[icomp] = convel->Map().LID(dofids[icomp]);
        fluidvel(icomp) = (*convel)[lids[icomp]];
      }

#ifdef ORACLES
      // set relative convective velocity to zero at the walls
      // remark: - the wall factor does the same thing, this is just to be sure
      //         - this is a hack for ORACLES
      //         - this guarantees a zero transport velocity for no-slip Dirichlet walls
      //         - this could be done in clean way by using a vector holding the fluid Dirichlet conditions
      if ( (abs(lnode->X()[1])-0.0653) < 1.0E-9 or // top and bottom wall combustion chamber
          ( abs(lnode->X()[0]) < 1.0E-9 and abs(lnode->X()[1]) > (0.005-1.0E-9) ) ) // walls above and below step
        flvelrel.Clear();

      // this is an attempt to do this in a clean general way; it failed
      // set Dirichlet BC
      //for (int icomp=0; icomp<3; ++icomp)
      //{
      //  const int gid = dofids[icomp];
      //  IO::cout << gid << IO::endl;
      //  if(dbcmap->MyGID(gid))
      //  {
      //    IO::cout << "DBC overwritten" << IO::endl;
      //    IO::cout << nodecoord(0,0) << IO::endl;
      //    IO::cout << nodecoord(1,0) << IO::endl;
      //    IO::cout << nodecoord(2,0) << IO::endl;
      //    flvelrel(icomp) = 0.0;
      //  }
      //}
#endif

      LINALG::Matrix<3,1> flvelabs(true);
      // add fluid velocity (Navier Stokes solution) and relative flame velocity
      for (int icomp=0; icomp<3; ++icomp)
      {
        flvelabs(icomp) = fluidvel(icomp) + flvelrel(icomp);
        const int err = conveltmp->ReplaceMyValues(1,&flvelabs(icomp),&lids[icomp]);
        if (err) dserror("could not replace values for convective velocity");
      }
    }
#ifdef COMBUST_GMSH_NORMALFIELD
    gmshfilecontent << "};\n";
  }
  // Gmsh output curvature
  {
    gmshfilecontent << "View \" " << "Curvature \" {\n";
    // loop over nodes on this processor
    for(int lnodeid=0; lnodeid < fluiddis->NumMyRowNodes(); ++lnodeid)
    {
      // get the processor local node
      DRT::Node*  lnode = fluiddis->lRowNode(lnodeid);
      const int gid = lnode->Id();
      const int nodelid = (*curvature).Map().LID(gid);
      if (nodelid<0) dserror("Proc %d: Cannot find gid=%d in Epetra_Vector",(*gradphinp).Comm().MyPID(),gid);
      const double curv = (*curvature)[nodelid];
      LINALG::SerialDenseMatrix xyz(3,1);
      xyz(0,0) = lnode->X()[0];
      xyz(1,0) = lnode->X()[1];
      xyz(2,0) = lnode->X()[2];

      IO::GMSH::cellWithScalarToStream(DRT::Element::point1, curv, xyz, gmshfilecontent);
    }
    gmshfilecontent << "};\n";
  }
  gmshfilecontent.close();
  if (true) IO::cout << " done" << IO::endl;
#endif

  return conveltmp;
}


/*------------------------------------------------------------------------------------------------*
 | protected: FGI iteration converged?                                                henke 06/08 |
 *------------------------------------------------------------------------------------------------*/
bool COMBUST::Algorithm::NotConvergedFGI()
{
  bool notconverged = true;

  if (fgitermax_ == 0)
    dserror("A maximum of 0 FGI is not sensible!!!");

  if (fgiter_ > 0) // nothing to do for FGI = 0 (in particular, do not overwrite phinpi_, since it corresponds to the fluid solution of the previous time step)
  {
    const Teuchos::RCP<const Epetra_Vector> velnpip = FluidField()->StdVelnp();
    const Teuchos::RCP<const Epetra_Vector> phinpip = ScaTraField()->Phinp();

    double velnormL2 = 1.0;
    double gfuncnormL2 = 1.0;

    velnpip->Norm2(&velnormL2);
    phinpip->Norm2(&gfuncnormL2);

    if (velnormL2 < 1e-5) velnormL2 = 1.0;
    if (gfuncnormL2 < 1e-5) gfuncnormL2 = 1.0;

    double fgvelnormL2 = 1.0;
    double fggfuncnormL2 = 1.0;

    // compute increment and L2-norm of increment
    Teuchos::RCP<Epetra_Vector> incvel = Teuchos::rcp(new Epetra_Vector(velnpip->Map()),true);
    incvel->Update(1.0,*velnpip,-1.0,*velnpi_,0.0);
    incvel->Norm2(&fgvelnormL2);

    Teuchos::RCP<Epetra_Vector> incgfunc = Teuchos::rcp(new Epetra_Vector(phinpip->Map(),true));//*ScaTraField()->Discretization()->DofRowMap()),true);
    incgfunc->Update(1.0,*phinpip,-1.0,*phinpi_,0.0);
    incgfunc->Norm2(&fggfuncnormL2);

    if (fgitermax_ > 1)
    {
      if (Comm().MyPID()==0)
      {
        printf("\n|+------------------------ FGI ------------------------+|");
        printf("\n|iter/itermax|----tol-[Norm]--|-fluid inc--|-g-func inc-|");
        printf("\n|   %2d/%2d    | %10.3E[L2] | %10.3E | %10.3E |",
            fgiter_,fgitermax_,convtol_,fgvelnormL2/velnormL2,fggfuncnormL2/gfuncnormL2);
        printf("\n|+-----------------------------------------------------+|\n");
      } // end if processor 0 for output

      if ((fgvelnormL2/velnormL2 <= convtol_) and (fggfuncnormL2/gfuncnormL2 <= convtol_))
        notconverged = false;

      if (Comm().MyPID()==0)
      {
        if ((fgiter_ == fgitermax_) and (notconverged == true))
        {
          printf("|+---------------- not converged ----------------------+|");
          printf("\n|+-----------------------------------------------------+|\n");
        }
      }
    } // end if fgitermax > 1

    if (fgiter_ >= fgitermax_)
      notconverged = false;

    if (!notconverged)
      Teuchos::rcp_dynamic_cast<FLD::CombustFluidImplicitTimeInt>(FluidField())->ClearTimeInt(); // clear XFEM time integration

    // update fgi-vectors
    // phi vector used for final fluid solution
    velnpi_->Update(1.0,*velnpip,0.0);
    phinpi_->Update(1.0,*phinpip,0.0);
  }

  return notconverged;
}


/*------------------------------------------------------------------------------------------------*
 | protected: do a stationary first time step for combustion algorithm               schott 08/10 |
 *------------------------------------------------------------------------------------------------*/
void COMBUST::Algorithm::SolveInitialStationaryProblem()
{
  if (Comm().MyPID()==0)
  {
    printf("==============================================================================================\n");
    printf("----------------Stationary timestep prepares instationary algorithm---------------------------\n");
    printf("==============================================================================================\n");
  }
  //-----------------------------
  // prepare stationary algorithm
  //-----------------------------
  fgiter_ = 0;

  // check if ScaTraField()->initialvelset == true
  /* remark: initial velocity field has been transferred to scalar transport field in constructor of
   * ScaTraFluidCouplingAlgorithm (initialvelset_ == true). Time integration schemes, such as
   * the one-step-theta scheme, are thus initialized correctly.
   */

  // modify time and timestep for stationary timestep
  SetTimeStep(0.0,0); // algorithm timestep

  if (Comm().MyPID()==0)
  {
    IO::cout << "----------------------Combustion-------  time step "
             << std::setw(2) << Step() << " ----------------------------------------\nTIME: "
             << std::setw(11) << std::setprecision(4) << std::scientific << Time() << "/"
             << std::setw(11) << std::setprecision(4) << std::scientific << MaxTime() << "  DT = "
             << std::setw(11) << std::setprecision(4) << std::scientific << Dt() << " STEP = "
             << std::setw(4) << Step() << "/" << std::setw(4) << NStep() << IO::endl;
  }

  FluidField()->PrepareTimeStep();

  //-------------------------------------
  // solve nonlinear Navier-Stokes system
  //-------------------------------------
  DoFluidField();

  // update field vectors
  UpdateInterface();

  // assign the fluid velocity field to the G-function as convective velocity field
  SetVelocityLevelSet(true);

  // write output to screen and files (and Gmsh)
  Output();

  return;
}


/*------------------------------------------------------------------------------------------------*
 | protected: prepare time step for combustion algorithm                              henke 06/08 |
 *------------------------------------------------------------------------------------------------*/
void COMBUST::Algorithm::PrepareTimeStep()
{
  IncrementTimeAndStep();
  fgiter_ = 0;

  if (Comm().MyPID()==0)
  {
    IO::cout << "----------------------Combustion-------  time step "
             << std::setw(2) << Step() << " ----------------------------------------\nTIME: "
             << std::setw(11) << std::setprecision(4) << std::scientific << Time() << "/"
             << std::setw(11) << std::setprecision(4) << std::scientific << MaxTime() << "  DT = "
             << std::setw(11) << std::setprecision(4) << std::scientific << Dt() << " STEP = "
             << std::setw(4) << Step() << "/" << std::setw(4) << NStep() << IO::endl;
  }

  FluidField()->PrepareTimeStep();

  // prepare time step
  // remark: initial velocity field has been transferred to scalar transport field in constructor of
  //         ScaTraFluidCouplingAlgorithm (initialvelset_ == true). Time integration schemes, such
  //         as the one-step-theta scheme, are thus initialized correctly.
  ScaTraField()->PrepareTimeStep();

  // synchronicity check between combust algorithm and base algorithms
  if (FluidField()->Time() != Time())
    dserror("Time in Fluid time integration differs from time in combustion algorithm");
  if (ScaTraField()->Time() != Time())
    dserror("Time in ScaTra time integration  differs from time in combustion algorithm");

  return;
}

/*------------------------------------------------------------------------------------------------*
 | protected: prepare time step for combustion algorithm                              henke 06/08 |
 *------------------------------------------------------------------------------------------------*/
void COMBUST::Algorithm::PrepareFGIteration()
{
  fgiter_ += 1;
  if (Comm().MyPID()==0)
  {
    printf("\n---------------------------------------  FGI loop: iteration number: %2d ----------------------\n",fgiter_);
  }
}

/*------------------------------------------------------------------------------------------------*
 | protected: perform a fluid time integration step                                   henke 06/08 |
 *------------------------------------------------------------------------------------------------*/
void COMBUST::Algorithm::DoFluidField()
{
  if (Comm().MyPID()==0)
  {
    IO::cout<<"\n---------------------------------------  FLUID SOLVER  ---------------------------------------" << IO::endl;
  }

  // update fluid interface with flamefront
  Teuchos::rcp_dynamic_cast<FLD::CombustFluidImplicitTimeInt>(FluidField())->ImportFlameFront(flamefront_,true);

  // solve nonlinear Navier-Stokes equations
  FluidField()->Solve();

  // clear fluid's memory to flamefront
  Teuchos::rcp_dynamic_cast<FLD::CombustFluidImplicitTimeInt>(FluidField())->ImportFlameFront(Teuchos::null,false);


  return;
}

/*------------------------------------------------------------------------------------------------*
 | protected: perform a G-function time integration step                              henke 06/08 |
 *------------------------------------------------------------------------------------------------*/
void COMBUST::Algorithm::DoGfuncField()
{
  if (Comm().MyPID()==0)
  {
    IO::cout<<"\n---------------------------------------  G-FUNCTION SOLVER  ----------------------------------\n";
  }

  // assign the fluid velocity field to the G-function as convective velocity field
  SetVelocityLevelSet();

  //solve convection-diffusion equation
  ScaTraField()->Solve();

  return;
}

/*------------------------------------------------------------------------------------------------*
 | protected: update                                                                  henke 06/08 |
 *------------------------------------------------------------------------------------------------*/
void COMBUST::Algorithm::UpdateInterface()
{
  //overwrite old interfacehandle before updating flamefront in first FGI
  if (fgiter_<=1)
    flamefront_->UpdateOldInterfaceHandle();

  // update flame front according to evolved G-function field
  // remark: for only one FGI iteration, 'phinpip_' == ScaTraField()->Phin()
  // TODO: hier scheitern die Gen-Alpha-Rechnungen, wenn Scatra nicht genalpha verwendet
  if (FluidField()->TimIntScheme() == INPAR::FLUID::timeint_afgenalpha)
    flamefront_->UpdateFlameFront(combustdyn_, phinpi_, ScaTraField()->Phiaf());
  else
    flamefront_->UpdateFlameFront(combustdyn_, phinpi_, ScaTraField()->Phinp());

  // remark: Except for initialization (contructor), restart and redistribution, flamefront_->UpdateFlameFront()
  //         must not be called at any other places of the algorithm, that is within the time loop! Otherwise
  //         the history of the interface is lost and semi-Lagrange time-integration will fail.

  return;
}

/*------------------------------------------------------------------------------------------------*
 | protected: update                                                                  henke 06/08 |
 *------------------------------------------------------------------------------------------------*/
void COMBUST::Algorithm::UpdateTimeStep()
{
  FluidField()->Update();

  ScaTraField()->Update();

  return;
}

/*------------------------------------------------------------------------------------------------*
 | protected: output                                                                  henke 06/08 |
 *------------------------------------------------------------------------------------------------*/
void COMBUST::Algorithm::Output()
{
  // Note: The order is important here! In here control file entries are
  // written. And these entries define the order in which the filters handle
  // the Discretizations, which in turn defines the dof number ordering of the
  // Discretizations.
  //------------------------------------------------------------------------------------------------
  // this hack is necessary for the visualization of disconituities in Gmsh             henke 10/09
  //------------------------------------------------------------------------------------------------
  // show flame front to fluid time integration scheme
  Teuchos::rcp_dynamic_cast<FLD::CombustFluidImplicitTimeInt>(FluidField())->ImportFlameFront(flamefront_,false);
  // write fluid output
  FluidField()->Output();
  // clear fluid's memory to flamefront
  Teuchos::rcp_dynamic_cast<FLD::CombustFluidImplicitTimeInt>(FluidField())->ImportFlameFront(Teuchos::null,false);


  // causes error in DEBUG mode (trueresidual_ is null)
  //FluidField()->LiftDrag();
  ScaTraField()->Output();

  // we have to call the output of averaged fields for scatra separately
  if (FluidField()->TurbulenceStatisticManager() != Teuchos::null)
    FluidField()->TurbulenceStatisticManager()->DoOutputForScaTra(*ScaTraField()->DiscWriter(),ScaTraField()->Step());

  // write position of center of mass to file
  if (DRT::INPUT::IntegralValue<bool>(combustdyn_,"WRITE_CENTER_OF_MASS"))
    CenterOfMass();

  return;
}


/*------------------------------------------------------------------------------------------------*
 | set velocity field for level-set algorithm                                     rasthofer 03/14 |
 *------------------------------------------------------------------------------------------------*/
void COMBUST::Algorithm::SetVelocityLevelSet(bool init)
{
  // get the convel at the correct time
  const Teuchos::RCP<const Epetra_Vector> convel = (ScaTraField()->MethodName() == INPAR::SCATRA::timeint_gen_alpha)
                                                   ?(FluidField()->StdVelaf())
                                                   :(FluidField()->StdVelnp());

  // assign the fluid velocity field to the G-function as convective velocity field
  switch(combusttype_)
  {
    case INPAR::COMBUST::combusttype_twophaseflow:
    case INPAR::COMBUST::combusttype_twophaseflow_surf:
    case INPAR::COMBUST::combusttype_twophaseflowjump:
    {
      // for two-phase flow, the fluid velocity field is continuous; it can be directly transferred to
      // the scalar transport field

      if (ScaTraField()->MethodName() != INPAR::SCATRA::timeint_gen_alpha)
        Teuchos::rcp_dynamic_cast<SCATRA::LevelSetAlgorithm>(ScaTraField())->SetVelocityField(convel,
                                                                                              Teuchos::null,
                                                                                              Teuchos::null,
                                                                                              Teuchos::null,
                                                                                              false,
                                                                                              init);
      else // temporary solution, since level-set algorithm does not yet support gen-alpha
      {
        if (Comm().MyPID()==0)
          std::cout << "CORRECT THESE LINES WHEN LEVEL-SET ALOGITHM SUPPORTS GEN-ALPHA" << std::endl;

        ScaTraField()->SetVelocityField(convel,
                                        Teuchos::null,
                                        Teuchos::null,
                                        Teuchos::null,
                                        1);
      }

      // transfer history vector only for subgrid-velocity: remark: complete RBVMM with cross- and Reynolds-stress term
      // in level-set field does not work since the missing enrichment

      break;
    }
    case INPAR::COMBUST::combusttype_premixedcombustion:
    {
      // for combustion, the velocity field is discontinuous; the relative flame velocity is added
      if (ScaTraField()->MethodName() != INPAR::SCATRA::timeint_gen_alpha)
        Teuchos::rcp_dynamic_cast<SCATRA::LevelSetAlgorithm>(ScaTraField())->SetVelocityField(ComputeFlameVel(convel,FluidField()->DofSet()),
                                                                                              Teuchos::null,
                                                                                              Teuchos::null,
                                                                                              Teuchos::null,
                                                                                              false,
                                                                                              init);
      else // temporary solution, since level-set algorithm does not yet support gen-alpha
      {
        if (Comm().MyPID()==0)
          std::cout << "CORRECT THESE LINES WHEN LEVEL-SET ALOGITHM SUPPORTS GEN-ALPHA" << std::endl;

        ScaTraField()->SetVelocityField(ComputeFlameVel(convel,FluidField()->DofSet()),
                                        Teuchos::null,
                                        Teuchos::null,
                                        Teuchos::null,
                                        1);
      }

      break;
    }
    default:
      dserror("unknown type of combustion problem");
      break;
  }

  return;
}


/*------------------------------------------------------------------------------------------------*
 | compute center of mass                                                         rasthofer 12/13 |
 *------------------------------------------------------------------------------------------------*/
void COMBUST::Algorithm::CenterOfMass()
{
  // get interface handle
  const Teuchos::RCP<InterfaceHandleCombust> ih = flamefront_->InterfaceHandle();

  // initialize vector for local sum
  std::vector<double > my_part(3);
  for (int idim=0; idim<3; idim++)
    my_part[idim] = 0.0;
  // initialize local volume
  double my_vol = 0.0;

  const Epetra_Map* elerowmap = FluidField()->Discretization()->ElementRowMap();

  // loop all row elements
  for (int iele=0; iele<FluidField()->Discretization()->NumMyRowElements(); iele++)
  {
    const int ele_gid = elerowmap->GID(iele);
    // get integration cells
    GEO::DomainIntCells domainIntCells = ih->ElementDomainIntCells(ele_gid);
    for (GEO::DomainIntCells::const_iterator cell = domainIntCells.begin(); cell != domainIntCells.end(); ++cell)
    {
      // include all contributions to minus domain
      if (cell->getDomainPlus() == false)
      {
        // get center of cell
        const LINALG::Matrix<3,1> center_cell = cell->GetPhysicalCenterPosition();
        // get volume of cell
        const double vol_cell = cell->VolumeInPhysicalDomain();

        // sum up
        for (int idim=0; idim<3; idim++)
          my_part[idim] += center_cell(idim)*vol_cell;

        my_vol += vol_cell;
      }
    }
  }

  // initalize vector for global sum
  std::vector<double > center(3);
  for (int idim=0; idim<3; idim++)
    center[idim] = 0.0;
  // initialize global volume
  double vol = 0.0;

  // compute sum over all procs
  FluidField()->Discretization()->Comm().SumAll(&(my_part[0]),&(center[0]),3);
  FluidField()->Discretization()->Comm().SumAll(&my_vol,&vol,1);

  // finally compute center position
  for (int idim=0; idim<3; idim++)
    center[idim] /= vol;

  if (Comm().MyPID()==0)
  {

    // write to file
    const std::string simulation = DRT::Problem::Instance()->OutputControlFile()->FileName();
    const std::string fname = simulation+"_center_of_mass.txt";

    if (Step()==0)
    {
      std::ofstream f;
      f.open(fname.c_str());
      f << "#| Step | Time |       x       |       y       |       z       |\n";
      f << "  "<< std::setw(2) << std::setprecision(10) << Step() << "    " << std::setw(3)<< std::setprecision(5)
        << Time() << std::setw(4) << std::setprecision(8) << "  "
        << center[0] << "    " << std::setprecision(8) << center[1] << "    " << std::setprecision(8) << center[2] << " " <<"\n";

      f.flush();
      f.close();
    }
    else
    {
      std::ofstream f;
      f.open(fname.c_str(),std::fstream::ate | std::fstream::app);
      f << "  "<< std::setw(2) << std::setprecision(10) << Step() << "    " << std::setw(3)<< std::setprecision(5)
      << Time() << std::setw(4) << std::setprecision(8) << "  "
      << center[0] << "    " << std::setprecision(8) << center[1] << "    " << std::setprecision(8) << center[2] << " " <<"\n";

      f.flush();
      f.close();
    }

  }

  return;
}


/* -------------------------------------------------------------------------------*
 | Restart a combustion problem                                          rasthofer|
 | remark: G-function is solved before fluid                                      |
 |         switching the order would affect the restart                           |
 * -------------------------------------------------------------------------------*/
void COMBUST::Algorithm::Restart(int step, const bool restartscatrainput, const bool restartfromfluid)
{
  if (Comm().MyPID()==0)
  {
    IO::cout << "---------------------------------------------" << IO::endl;
    IO::cout << "| restart of combustion problem             |" << IO::endl;
    if (restartscatrainput)
      IO::cout << "| restart with scalar field from input file |" << IO::endl;
    if (restartfromfluid)
      IO::cout << "| restart from standard fluid problem       |" << IO::endl;
    IO::cout << "---------------------------------------------" << IO::endl;
  }
  if (restartfromfluid and !restartscatrainput)
    dserror("Scalar field has to be read from input file for restart from standard fluid");

  // read level-set field from input file instead of restart file
  Teuchos::RCP<Epetra_Vector> oldphidtn  = Teuchos::null;
  Teuchos::RCP<Epetra_Vector> oldphin    = Teuchos::null;
  Teuchos::RCP<Epetra_Vector> oldphinp   = Teuchos::null;
  if (restartscatrainput)
  {
    oldphin  = Teuchos::rcp(new Epetra_Vector(*(ScaTraField()->Phin())));
    oldphinp = Teuchos::rcp(new Epetra_Vector(*(ScaTraField()->Phinp())));
    if (ScaTraField()->MethodName() == INPAR::SCATRA::timeint_gen_alpha)
    {
      dserror("Check restart from scatra input for gen alpha!");
      // do we have to keep it?
      oldphidtn= Teuchos::rcp(new Epetra_Vector(*(ScaTraField()->Phidtn())));
    }
  }

  // restart of scalar transport (G-function) field
  if (!restartfromfluid) // not if restart is done from standard fluid field; there is no scalar field
  {
    // for particle level-set method, we cannot call the level-set restart, since the are no particles yet
    // (restart from flow generation via combustion fluid)
    if (!restartscatrainput)
      ScaTraField()->ReadRestart(step);
    else
    {
      // call restart function of base time integrator only
      if (ScaTraField()->MethodName() == INPAR::SCATRA::timeint_one_step_theta)
        Teuchos::rcp_dynamic_cast<SCATRA::TimIntOneStepTheta>(ScaTraField())->SCATRA::TimIntOneStepTheta::ReadRestart(step);
      else // particles are not yet supported by other time integration schemes
        ScaTraField()->ReadRestart(step);
    }
  }

  // get pointers to the discretizations from the time integration scheme of each field
  const Teuchos::RCP<DRT::Discretization> gfuncdis = ScaTraField()->Discretization();

  //-------------------------------------------------------------
  // create (old) flamefront conforming to restart state of fluid
  //-------------------------------------------------------------
  flamefront_->UpdateFlameFront(combustdyn_, ScaTraField()->Phin(), ScaTraField()->Phinp(),
               Teuchos::rcp_dynamic_cast<FLD::CombustFluidImplicitTimeInt>(FluidField())->ReadPhinp(step));

  // show flame front to fluid time integration scheme
  Teuchos::rcp_dynamic_cast<FLD::CombustFluidImplicitTimeInt>(FluidField())->ImportFlameFront(flamefront_,true);
  // delete fluid's memory of flame front; it should never have seen it in the first place!
  Teuchos::rcp_dynamic_cast<FLD::CombustFluidImplicitTimeInt>(FluidField())->ImportFlameFront(Teuchos::null,false);

  // restart of fluid field
  FluidField()->ReadRestart(step);

  // read level-set field from input file instead of restart file
  if (restartscatrainput)
  {
    if (Comm().MyPID()==0)
      IO::cout << "---  overwriting scalar field with field from input file... " << IO::endl;
    // now overwrite restart phis w/ the old phis
    ScaTraField()->Phinp()->Update(1.0, *(oldphinp), 0.0);
    ScaTraField()->Phin() ->Update(1.0, *(oldphin),  0.0);
    if (ScaTraField()->MethodName() == INPAR::SCATRA::timeint_gen_alpha)
    {
      ScaTraField()->Phidtn() ->Update(1.0, *(oldphidtn),  0.0);
      ScaTraField()->ComputeIntermediateValues();
    }
    else if (ScaTraField()->MethodName() == INPAR::SCATRA::timeint_one_step_theta)
    {
      ScaTraField()->Phidtnp()->PutScalar(0.0);
      // calls CalcInitialTimeDerivative()
      // ApplyDirichletBC() called within this function doesn't pose any problem since we
      // don't have any DBCs for level-set problems
      ScaTraField()->PrepareFirstTimeStep();
      // CalcInitialTimeDerivative() copies phidtnp to phidtn
    }
    else
      dserror("Unknown time-integration scheme!");
    if (Comm().MyPID()==0)
      IO::cout << "done" << IO::endl;

    //-------------------------------------------------------------
    // fill flamefront conforming to restart state
    //-------------------------------------------------------------
    flamefront_->UpdateFlameFront(combustdyn_, ScaTraField()->Phin(), ScaTraField()->Phinp(),
            Teuchos::rcp_dynamic_cast<FLD::CombustFluidImplicitTimeInt>(FluidField())->ReadPhinp(step));
    flamefront_->UpdateOldInterfaceHandle();

    // show flame front to fluid time integration scheme
    Teuchos::rcp_dynamic_cast<FLD::CombustFluidImplicitTimeInt>(FluidField())->ImportFlameFront(flamefront_,true);
    // delete fluid's memory of flame front; it should never have seen it in the first place!
    Teuchos::rcp_dynamic_cast<FLD::CombustFluidImplicitTimeInt>(FluidField())->ImportFlameFront(Teuchos::null,false);

    if (gmshoutput_)
    {
      //--------------------------
      // write output to Gmsh file
      //--------------------------
      const std::string filename = IO::GMSH::GetNewFileNameAndDeleteOldFiles("field_scalar_restart_before_reinit", Step(), 701, true, gfuncdis->Comm().MyPID());
      std::ofstream gmshfilecontent(filename.c_str());
      {
        // add 'View' to Gmsh postprocessing file
        gmshfilecontent << "View \" " << "Phinp \" {" << std::endl;
        // draw scalar field 'Phinp' for every element
        IO::GMSH::ScalarFieldToGmsh(gfuncdis,ScaTraField()->Phinp(),gmshfilecontent);
        gmshfilecontent << "};" << std::endl;
      }
      {
        // add 'View' to Gmsh postprocessing file
        gmshfilecontent << "View \" " << "Phin \" {" << std::endl;
        // draw scalar field 'Phin' for every element
        IO::GMSH::ScalarFieldToGmsh(gfuncdis,ScaTraField()->Phin(),gmshfilecontent);
        gmshfilecontent << "};" << std::endl;
      }
      if (ScaTraField()->MethodName() == INPAR::SCATRA::timeint_gen_alpha)
      {
        // add 'View' to Gmsh postprocessing file
        gmshfilecontent << "View \" " << "Phidtn \" {" << std::endl;
        // draw scalar field 'Phidtn' for every element
        IO::GMSH::ScalarFieldToGmsh(gfuncdis,ScaTraField()->Phidtn(),gmshfilecontent);
        gmshfilecontent << "};" << std::endl;
      }
      {
        //    // add 'View' to Gmsh postprocessing file
        //    gmshfilecontent << "View \" " << "Convective Velocity \" {" << std::endl;
        //    // draw vector field 'Convective Velocity' for every element
        //    IO::GMSH::VectorFieldNodeBasedToGmsh(gfuncdis,ScaTraField()->ConVel(),gmshfilecontent);
        //    gmshfilecontent << "};" << std::endl;
      }
      gmshfilecontent.close();
      IO::cout << " done" << IO::endl;
    }
  }

  phinpi_->Update(1.0,*(Teuchos::rcp_dynamic_cast<FLD::CombustFluidImplicitTimeInt>(FluidField())->ReadPhinp(step)),0.0);

  if (gmshoutput_)
  {
    //--------------------------
    // write output to Gmsh file
    //--------------------------
    const std::string filename = IO::GMSH::GetNewFileNameAndDeleteOldFiles("field_scalar_restart_after_reinit", Step(), 701, true, gfuncdis->Comm().MyPID());
    std::ofstream gmshfilecontent(filename.c_str());
    {
      // add 'View' to Gmsh postprocessing file
      gmshfilecontent << "View \" " << "Phinp \" {" << std::endl;
      // draw scalar field 'Phinp' for every element
      IO::GMSH::ScalarFieldToGmsh(gfuncdis,ScaTraField()->Phinp(),gmshfilecontent);
      gmshfilecontent << "};" << std::endl;
    }
    {
      // add 'View' to Gmsh postprocessing file
      gmshfilecontent << "View \" " << "Phin \" {" << std::endl;
      // draw scalar field 'Phin' for every element
      IO::GMSH::ScalarFieldToGmsh(gfuncdis,ScaTraField()->Phin(),gmshfilecontent);
      gmshfilecontent << "};" << std::endl;
    }
    if (ScaTraField()->MethodName() == INPAR::SCATRA::timeint_gen_alpha)
    {
      // add 'View' to Gmsh postprocessing file
      gmshfilecontent << "View \" " << "Phidtn \" {" << std::endl;
      // draw scalar field 'Phidtn' for every element
      IO::GMSH::ScalarFieldToGmsh(gfuncdis,ScaTraField()->Phidtn(),gmshfilecontent);
      gmshfilecontent << "};" << std::endl;
    }
    {
      //    // add 'View' to Gmsh postprocessing file
      //    gmshfilecontent << "View \" " << "Convective Velocity \" {" << std::endl;
      //    // draw vector field 'Convective Velocity' for every element
      //    IO::GMSH::VectorFieldNodeBasedToGmsh(gfuncdis,ScaTraField()->ConVel(),gmshfilecontent);
      //    gmshfilecontent << "};" << std::endl;
    }
    gmshfilecontent.close();
    IO::cout << " done" << IO::endl;
  }

  //FluidField()->Output();
  //ScaTraField()->Output();

  // set time in scalar transport time integration scheme
  // this has also to be done for restartscatrainput, since a particle field may now be included
  if(restartfromfluid or restartscatrainput)
    ScaTraField()->SetTimeStep(FluidField()->Time(),step);

  SetTimeStep(FluidField()->Time(),step);

  // assign the fluid velocity field to the G-function as convective velocity field
  // bool true allows for setting old convective velocity required for particle coupling
  if (not gen_flow_)
    SetVelocityLevelSet(true);

  //UpdateTimeStep();

  restart_ = true;

  return;
}

/*------------------------------------------------------------------------------------------------*
 | protected: redistribute scatra and fluid discretization                        rasthofer 07/11 |
 |                                                                                    DA wichmann |
 *------------------------------------------------------------------------------------------------*/
void COMBUST::Algorithm::Redistribute()
{
  const double ratiolimitfac = combustdyn_.get<double>("PARALLEL_REDIST_RATIO_FAC");

  if (Comm().NumProc() > 1 and ratiolimitfac > 1.0)
  {
    // for ease of changing some constants are set here
    double normalnodeweight   = 10.0;
    double enrichednodeweight = 2000.0;
    double slavenodeweight    = 1.0;

    // determine whether a redistribution is necessary
    // by comparing min and max evaltime of all procs
    double myprocevaltime  = FluidField()->EvalTime();
    double minprocevaltime = 0.0;
    double maxprocevaltime = 0.0;

    Comm().MinAll(&myprocevaltime, &minprocevaltime, 1);
    Comm().MaxAll(&myprocevaltime, &maxprocevaltime, 1);

    // Calculate the evaluation time ratio by which a further redistribution is determined.
    // For this the minimum ratio since the last redistribution is used
    if (minprocevaltime > 0.0 and maxprocevaltime / minprocevaltime < evaltimeratio_) // the order is important, as it avoids division by zero
      evaltimeratio_ = maxprocevaltime / minprocevaltime;

    if (minprocevaltime <= 0.0)
    {
      if (Comm().MyPID() == 0)
        IO::cout << "The max / min ratio could not be determined --> Not redistributing" << IO::endl;
    }
    else if (maxprocevaltime / minprocevaltime < ratiolimitfac * evaltimeratio_)
    {
      if (Comm().MyPID() == 0)
        IO::cout << "The max / min ratio is " << maxprocevaltime / minprocevaltime << " < " << ratiolimitfac * evaltimeratio_ << " --> Not redistributing" << IO::endl;
    }
    else
    {
      if (Comm().MyPID() == 0)
      {
        IO::cout << "-------------------------------Redistributing-------------------------------" << IO::endl;
        IO::cout << "The max / min ratio is " << maxprocevaltime / minprocevaltime << " > " << ratiolimitfac * evaltimeratio_ << " --> Redistributing" << IO::endl;
      }
      //--------------------------------------------------------------------------------------
      // Building graph for later use by parmetis
      //--------------------------------------------------------------------------------------
      const Teuchos::RCP<DRT::Discretization> fluiddis = FluidField()->Discretization();
      const Teuchos::RCP<DRT::Discretization> gfuncdis = ScaTraField()->Discretization();
      const Epetra_Map* noderowmap = fluiddis->NodeRowMap();
      Teuchos::RCP<std::map<int,std::vector<int> > > allcoupledcolnodes = gfuncdis->GetAllPBCCoupledColNodes();

      // weights for graph partition
      Epetra_Vector weights(*noderowmap,false);
      weights.PutScalar(normalnodeweight);

      // loop elements and raise the weight of nodes of cut elements to 1000
      {
        for (int iele = 0; iele < fluiddis->NumMyRowElements(); ++iele)
        {
          DRT::Element* ele = fluiddis->lRowElement(iele);
          if (flamefront_->InterfaceHandle()->ElementCutStatus(ele->Id()) != COMBUST::InterfaceHandleCombust::uncut)
          {
            const int* nodes = ele->NodeIds();
            for (int i = 0; i < ele->NumNode(); ++i)
            {
              int gid = nodes[i];
              weights.ReplaceGlobalValues(1,&enrichednodeweight,&gid);
            }
          }
        }
      }

      // ----------------------------------------
      // loop masternodes to adjust weights of slavenodes
      // they need a small weight since they do not contribute any dofs
      // to the linear system
      {
        std::map<int,std::vector<int> >::iterator masterslavepair;

        for(masterslavepair =allcoupledcolnodes->begin();
            masterslavepair!=allcoupledcolnodes->end()  ;
            ++masterslavepair)
        {
          // get masternode
          DRT::Node*  master = fluiddis->gNode(masterslavepair->first);

          if (master->Owner() != Comm().MyPID())
            continue;

          // loop slavenodes associated with master
          for(std::vector<int>::iterator iter=masterslavepair->second.begin();
              iter!=masterslavepair->second.end();++iter)
          {
            int gid       =*iter;
            weights.ReplaceGlobalValues(1,&slavenodeweight,&gid);
          }
        }
      }


      // allocate graph
      Teuchos::RCP<Epetra_CrsGraph> nodegraph = Teuchos::rcp(new Epetra_CrsGraph(Copy,*noderowmap,108,false));

      // -------------------------------------------------------------
      // iterate all elements on this proc including ghosted ones
      // compute connectivity

      // standard part without master<->slave coupling
      // Note:
      // if a proc stores the appropiate ghosted elements, the resulting
      // graph will be the correct and complete graph of the distributed
      // discretization even if nodes are not ghosted.

      for (int nele = 0; nele < fluiddis->NumMyColElements(); ++nele)
      {
        // get the element
        DRT::Element* ele = fluiddis->lColElement(nele);

        // get its nodes and nodeids
        const int  nnode   = ele->NumNode();
        const int* nodeids = ele->NodeIds();

        for (int row = 0; row < nnode; ++row)
        {
          const int rownode = nodeids[row];

          // insert into line of graph only when this proc owns the node
          if (!noderowmap->MyGID(rownode))
            continue;

          // insert all neighbors from element in the graph
          for (int col = 0; col < nnode; ++col)
          {
            int colnode = nodeids[col];
            int err = nodegraph->InsertGlobalIndices(rownode,1,&colnode);
            if (err<0)
              dserror("nodegraph->InsertGlobalIndices returned err=%d",err);
          }
        }
      }

      // -------------------------------------------------------------
      // additional coupling between master and slave
      // we do not only connect master and slave nodes but if a master/slave
      // is connected to a master/slave, we connect the corresponding slaves/master
      // as well

      for (int nele = 0; nele < fluiddis->NumMyColElements();++nele)
      {
        // get the element
        DRT::Element* ele = fluiddis->lColElement(nele);

        // get its nodes and nodeids
        const int  nnode   = ele->NumNode();
        const int* nodeids = ele->NodeIds();

        for (int row = 0; row < nnode; ++row)
        {
          const int rownode = nodeids[row];

          // insert into line of graph only when this proc owns the node
          if (!noderowmap->MyGID(rownode))
            continue;

          std::map<int,std::vector<int> >::iterator masterslavepair = allcoupledcolnodes->find(rownode);
          if(masterslavepair != allcoupledcolnodes->end())
          {
            // get all masternodes of this element
            for (int col = 0; col < nnode; ++col)
            {
              int colnode = nodeids[col];

              std::map<int,std::vector<int> >::iterator othermasterslavepair = allcoupledcolnodes->find(colnode);
              if(othermasterslavepair != allcoupledcolnodes->end())
              {
                // add connection to all slaves

                for(std::vector<int>::iterator iter = othermasterslavepair->second.begin();
                    iter != othermasterslavepair->second.end(); ++iter)
                {
                  int othermastersslaveindex = *iter;
                  int masterindex            = rownode;
                  int err = nodegraph->InsertGlobalIndices(rownode,1,&othermastersslaveindex);
                  if (err<0)
                    dserror("nodegraph->InsertGlobalIndices returned err=%d",err);

                  if (noderowmap->MyGID(*iter))
                  {
                    err = nodegraph->InsertGlobalIndices(*iter,1,&masterindex);
                    if (err<0)
                      dserror("nodegraph->InsertGlobalIndices returned err=%d",err);
                  }
                }
              }
            }
          }
        }
      }

      // finalize construction of initial graph
      int err = nodegraph->FillComplete();
      if (err)
        dserror("graph->FillComplete returned %d",err);


      //--------------------------------------------------------------------------------------
      // prepare and call METIS
      //--------------------------------------------------------------------------------------

      const int myrank   = nodegraph->Comm().MyPID();
      const int numproc  = nodegraph->Comm().NumProc();

      // proc that will do the serial partitioning
      // the graph is collapsed to this proc
      // Normally this would be proc 0 but 0 always has so much to do.... ;-)
      int workrank = 1;

      // get rowmap of the graph
      const Epetra_BlockMap& tmp = nodegraph->RowMap();
      Epetra_Map rowmap(tmp.NumGlobalElements(),tmp.NumMyElements(),
                        tmp.MyGlobalElements(),0,nodegraph->Comm());

      // -------------------------------------------------------------
      // build a target map that stores everything on proc workrank
      // We have arbirtary gids here and we do not tell metis about
      // them. So we have to keep rowrecv until the redistributed map is
      // build.

      // rowrecv is a fully redundant vector (size of number of nodes)
      std::vector<int> rowrecv(rowmap.NumGlobalElements());

      // after AllreduceEMap rowrecv contains
      //
      // *-+-+-    -+-+-*-+-+-    -+-+-*-           -*-+-+-    -+-+-*
      // * | | .... | | * | | .... | | * ..........  * | | .... | | *
      // *-+-+-    -+-+-*-+-+-    -+-+-*-           -*-+-+-    -+-+-*
      //   gids stored     gids stored                  gids stored
      //  on first proc  on second proc                 on last proc
      //
      // the ordering of the gids on the procs is arbitrary (as copied
      // from the map)
      LINALG::AllreduceEMap(rowrecv, rowmap);

      // construct an epetra map from the list of gids
      Epetra_Map tmap(rowmap.NumGlobalElements(),
                      // if ..........    then ............... else
                      (myrank == workrank) ? (int)rowrecv.size() : 0,
                      &rowrecv[0],
                      0,
                      rowmap.Comm());

      // export the graph to tmap
      Epetra_CrsGraph tgraph(Copy,tmap,108,false);
      Epetra_Export exporter(rowmap,tmap);
      {
        int err = tgraph.Export(*nodegraph,exporter,Add);
        if (err < 0)
          dserror("Graph export returned err=%d",err);
      }
      tgraph.FillComplete();
      tgraph.OptimizeStorage();

      // export the weights to tmap
      Epetra_Vector tweights(tmap,false);
      err = tweights.Export(weights,exporter,Insert);
      if (err < 0)
        dserror("Vector export returned err=%d",err);

      // metis requests indexes. So we need a reverse lookup from gids
      // to indexes.
      std::map<int,int> idxmap;
      // xadj points from index i to the index of the
      // first adjacent node
      std::vector<int> xadj  (rowmap.NumGlobalElements()+1);
      // a list of adjacent nodes, adressed using xadj
      std::vector<int> adjncy(tgraph.NumGlobalNonzeros()); // the size is an upper bound

      // This is a vector of size n that upon successful completion stores the partition vector of the graph
      std::vector<int> part(tmap.NumMyElements());

      // construct reverse lookup for all procs
      for (size_t i = 0; i < rowrecv.size(); ++i)
      {
        idxmap[rowrecv[i]] = i;
      }

      if (myrank == workrank)
      {
        // ----------------------------------------

        // rowrecv(i)       rowrecv(i+1)                      node gids
        //     ^                 ^
        //     |                 |
        //     | idxmap          | idxmap
        //     |                 |
        //     v                 v
        //     i                i+1                       equivalent indices
        //     |                 |
        //     | xadj            | xadj
        //     |                 |
        //     v                 v
        //    +-+-+-+-+-+-+-+-+-+-+                -+-+-+
        //    | | | | | | | | | | | ............... | | |      adjncy
        //    +-+-+-+-+-+-+-+-+-+-+                -+-+-+
        //
        //    |       i's       |    (i+1)'s
        //    |    neighbours   |   neighbours           (numbered by equivalent indices)
        //

        int count=0;
        xadj[0] = 0;
        for (int row=0; row<tgraph.NumMyRows(); ++row)
        {
          int grid = tgraph.RowMap().GID(row);
          int numindices;
          int* lindices;
          int err = tgraph.ExtractMyRowView(row,numindices,lindices);
          if (err)
            dserror("Epetra_CrsGraph::ExtractMyRowView returned err=%d",err);

          for (int col = 0; col < numindices; ++col)
          {
            int gcid = tgraph.ColMap().GID(lindices[col]);
            if (gcid == grid)
              continue;
            adjncy[count] = idxmap[gcid];
            ++count;
          }
          xadj[row+1] = count;
        }
      } // if (myrank == workrank)

      // broadcast xadj
      tmap.Comm().Broadcast(&xadj[0],xadj.size(),workrank);

      // broadcast adjacence (required for edge weights)
      int adjncysize = (int)adjncy.size();
      tmap.Comm().Broadcast(&adjncysize,1,workrank);
      adjncy.resize(adjncysize);
      tmap.Comm().Broadcast(&adjncy[0],adjncysize,workrank);

      // -------------------------------------------------------------
      // set a fully redundant vector of weights for edges
      std::vector<int> ladjwgt(adjncy.size(),0);
      std::vector<int>  adjwgt(adjncy.size(),0);

      for(std::vector<int>::iterator iter =ladjwgt.begin();
          iter!=ladjwgt.end();
          ++iter)
      {
        *iter=0;
      }

      // loop all master nodes on this proc
      {
        std::map<int,std::vector<int> >::iterator masterslavepair;

        for(masterslavepair =allcoupledcolnodes->begin();
            masterslavepair!=allcoupledcolnodes->end()  ;
            ++masterslavepair)
        {
          // get masternode
          DRT::Node*  master = fluiddis->gNode(masterslavepair->first);

          if(master->Owner()!=myrank)
          {
            continue;
          }

          std::map<int,int>::iterator paul = idxmap.find(master->Id());
          if (paul == idxmap.end())
            dserror("master not in reverse lookup");

          // inverse lookup
          int masterindex = idxmap[master->Id()];

          // loop slavenodes
          for(std::vector<int>::iterator iter = masterslavepair->second.begin();
              iter != masterslavepair->second.end(); ++iter)
          {
            DRT::Node*  slave = fluiddis->gNode(*iter);

            if(slave->Owner() != myrank)
              dserror("own master but not slave\n");

            int slaveindex = idxmap[slave->Id()];

            std::map<int,int>::iterator foo = idxmap.find(slave->Id());
            if (foo == idxmap.end())
              dserror("slave not in reverse lookup");

            // -------------------------------------------------------------
            // connections between master and slavenodes are very strong
            // we do not want to partition between master and slave nodes
            for(int j = xadj[masterindex]; j < xadj[masterindex+1]; ++j)
            {
              if(adjncy[j] == slaveindex)
              {
                ladjwgt[j] = 100;
              }
            }

            for(int j = xadj[slaveindex]; j < xadj[slaveindex+1]; ++j)
            {
              if(adjncy[j] == masterindex)
              {
                ladjwgt[j] = 100;
              }
            }
          }
        }
      }

      // do communication to aquire edge weight information from all procs
      tmap.Comm().SumAll(&ladjwgt[0], &adjwgt[0], adjwgt.size());

      // the standard edge weight is one
      for(std::vector<int>::iterator iter =adjwgt.begin();
          iter!=adjwgt.end();
          ++iter)
      {
        if(*iter==0)
        *iter=1;
      }

      // the reverse lookup is not required anymore
      idxmap.clear();

      // -------------------------------------------------------------
      // do partitioning using metis on workrank
      if (myrank == workrank)
      {
        // the vertex weights
        std::vector<int> vwgt(tweights.MyLength());
        for (int i=0; i<tweights.MyLength(); ++i) vwgt[i] = (int)tweights[i];

        // 0 No weights (vwgts and adjwgt are NULL)
        // 1 Weights on the edges only (vwgts = NULL)
        // 2 Weights on the vertices only (adjwgt = NULL)
        // 3 Weights both on vertices and edges.
        int wgtflag=3;
        // 0 C-style numbering is assumed that starts from 0
        // 1 Fortran-style numbering is assumed that starts from 1
        int numflag=0;
        // The number of parts to partition the graph.
        int npart=numproc;
        // This is an array of 5 integers that is used to pass parameters for the various phases of the algorithm.
        // If options[0]=0 then default values are used. If options[0]=1, then the remaining four elements of
        // options are interpreted as follows:
        // options[1]    Determines matching type. Possible values are:
        //               1 Random Matching (RM)
        //               2 Heavy-Edge Matching (HEM)
        //               3 Sorted Heavy-Edge Matching (SHEM) (Default)
        //               Experiments has shown that both HEM and SHEM perform quite well.
        // options[2]    Determines the algorithm used during initial partitioning. Possible values are:
        //               1 Region Growing (Default)
        // options[3]    Determines the algorithm used for re%Gï¬%@nement. Possible values are:
        //               1 Early-Exit Boundary FM re%Gï¬%@nement (Default)
        // options[4]    Used for debugging purposes. Always set it to 0 (Default).
        int options[5] = { 0,3,1,1,0 };
        // Upon successful completion, this variable stores the number of edges that are cut by the partition.
        int edgecut=0;
        // The number of vertices in the graph.
        int nummyele = tmap.NumMyElements();

        IO::cout << "proc " <<  myrank << " repartition graph using metis" << IO::endl;

        if (numproc<8) // better for smaller no. of partitions
        {
          METIS_PartGraphRecursive(&nummyele,
                                   &xadj[0],
                                   &adjncy[0],
                                   &vwgt[0],
                                   &adjwgt[0],
                                   &wgtflag,
                                   &numflag,
                                   &npart,
                                   options,
                                   &edgecut,
                                   &part[0]);

          IO::cout << "METIS_PartGraphRecursive produced edgecut of " << edgecut << IO::endl;
        }
        else
        {
          METIS_PartGraphKway(&nummyele,
                              &xadj[0],
                              &adjncy[0],
                              &vwgt[0],
                              &adjwgt[0],
                              &wgtflag,
                              &numflag,
                              &npart,
                              options,
                              &edgecut,
                              &part[0]);
          IO::cout << "METIS_PartGraphKway produced edgecut of " << edgecut << IO::endl;
        }
      } // if (myrank==workrank)

      // broadcast partitioning result
      int size = tmap.NumMyElements();
      tmap.Comm().Broadcast(&size,1,workrank);
      part.resize(size);
      tmap.Comm().Broadcast(&part[0],size,workrank);

      // loop part and count no. of nodes belonging to me
      // (we reuse part to save on memory)
      int count = 0;
      for (int i = 0; i < size; ++i)
      {
        if (part[i]==myrank)
        {
          part[count] = rowrecv[i];
          ++count;
        }
      }

      // rowrecv is done
      rowrecv.clear();

      // create map with new layout
      Epetra_Map newmap(size,count,&part[0],0,nodegraph->Comm());

      // create the new graph and export to it
      Teuchos::RCP<Epetra_CrsGraph> newnodegraph;

      newnodegraph = Teuchos::rcp(new Epetra_CrsGraph(Copy,newmap,108,false));
      Epetra_Export exporter2(nodegraph->RowMap(),newmap);
      err = newnodegraph->Export(*nodegraph,exporter2,Add);
      if (err<0)
        dserror("Graph export returned err=%d",err);
      newnodegraph->FillComplete();
      newnodegraph->OptimizeStorage();


      //--------------------------------------------------------------------------------------
      // Ensure the new distribution is valid
      //--------------------------------------------------------------------------------------
      {
        // the rowmap will become the new distribution of nodes
        const Epetra_BlockMap rntmp = newnodegraph->RowMap();
        Epetra_Map newnoderowmap(-1,rntmp.NumMyElements(),rntmp.MyGlobalElements(),0,Comm());

        // the column map will become the new ghosted distribution of nodes
        const Epetra_BlockMap Mcntmp = newnodegraph->ColMap();
        Epetra_Map newnodecolmap(-1,Mcntmp.NumMyElements(),Mcntmp.MyGlobalElements(),0,Comm());

        Teuchos::RCP<Epetra_Map> elerowmap;
        Teuchos::RCP<Epetra_Map> elecolmap;

        ScaTraField()->Discretization()->BuildElementRowColumn(newnoderowmap,newnodecolmap,elerowmap,elecolmap);

        int myrowele = elerowmap->NumMyElements();
        int minrowele = 1;
        Comm().MinAll(&myrowele, &minrowele, 1);
        if (minrowele <= 0)
        {
          if (Comm().MyPID() == 0)
            IO::cout << "A processor would be left without row elements --> Redistribution aborted\n"
                      << "----------------------------------------------------------------------------" << IO::endl;
          return;
        }
      }


      //--------------------------------------------------------------------------------------
      // call ScaTra and Fluid field and pass the new graph, so they can manage all
      // redistribution themselves
      //--------------------------------------------------------------------------------------
      if(Comm().MyPID()==0)
        IO::cout << "Redistributing ScaTra Discretization                                ... " << IO::endl;

      Teuchos::rcp_dynamic_cast<SCATRA::LevelSetAlgorithm>(ScaTraField())->Redistribute(newnodegraph);

      if(Comm().MyPID()==0)
        IO::cout << "Redistributing Fluid Discretization                                 ... " << IO::endl;

      FluidField()->Redistribute(newnodegraph);

      if(Comm().MyPID()==0)
        IO::cout << "done\nUpdating interface                                                  ... " << IO::endl;

      // update flame front according to evolved G-function field
      // remark: for only one FGI iteration, 'phinpip_' == ScaTraField().Phin()
      // TODO: passt das zu Martins XFEMTimInt?
      flamefront_->UpdateFlameFront(combustdyn_, ScaTraField()->Phin(), ScaTraField()->Phinp());

      if(Comm().MyPID()==0)
        IO::cout << "Transfering state vectors to new distribution                       ... " << IO::endl;

      FluidField()->TransferVectorsToNewDistribution(flamefront_);

      if(Comm().MyPID()==0)
        IO::cout << "done" << IO::endl;

      //--------------------------------------------------------------------------------------
      // with the scatra- and fluid-field updated it is time to do the same with the algorithm
      //--------------------------------------------------------------------------------------
      Teuchos::RCP<Epetra_Vector> old;

      if (velnpi_ != Teuchos::null)
      {
        old = velnpi_;
        velnpi_ = Teuchos::rcp(new Epetra_Vector(FluidField()->StdVelnp()->Map()),true);
        LINALG::Export(*old, *velnpi_);
      }

      if (phinpi_ != Teuchos::null)
      {
        old = phinpi_;
        phinpi_ = Teuchos::rcp(new Epetra_Vector(*gfuncdis->DofRowMap()),true);
        LINALG::Export(*old, *phinpi_);
      }

      //-----------------------------------------------------------------------------
      // redistibute vectors holding the means in the general mean statistics manager
      //-----------------------------------------------------------------------------
      // if applicable, provide scatra data to the turbulence statistics
      if (FluidField()->TurbulenceStatisticManager() != Teuchos::null)
      {
//        if(Comm().MyPID()==0)
//          IO::cout << "Updating pointer to ScaTra vector                                               ... " << IO::endl;

        // update the statistics manager to the new ScaTra discretization
//        FluidField()->TurbulenceStatisticManager()
//            ->AddScaTraResults(ScaTraField()->Discretization(),ScaTraField()->Phinp());

        if(Comm().MyPID()==0)
          IO::cout << "Redistributing General Mean Statistics Manager                      ... " << IO::endl;

//        // redistribute with redistributed standard fluid dofset
        if ( FluidField()->TurbulenceStatisticManager()->GetTurbulenceStatisticsGeneralMean()!=Teuchos::null)
          FluidField()->TurbulenceStatisticManager()->GetTurbulenceStatisticsGeneralMean()
              ->Redistribute(FluidField()->DofSet(), FluidField()->Discretization());

        if(Comm().MyPID()==0)
          IO::cout << "done" << IO::endl;
      }

      if (Comm().MyPID() == 0)
        IO::cout << "----------------------------------------------------------------------------" << IO::endl;

      // make sure the redistribution ratio will be reset
      evaltimeratio_ = 1e12;
    }
  }
  return;
}


/*----------------------------------------------------------------------*
 | perform result test                                  rasthofer 01/14 |
 *----------------------------------------------------------------------*/
void COMBUST::Algorithm::TestResults()
{
  DRT::Problem::Instance()->AddFieldTest(FluidField()->CreateFieldTest());

  if (ScaTraField()->MethodName() != INPAR::SCATRA::timeint_gen_alpha)
  {
    // DRT::Problem::Instance()->TestAll() is called in level-set field after adding particles
    Teuchos::rcp_dynamic_cast<SCATRA::LevelSetAlgorithm>(ScaTraField())->TestResults();
  }
  else
  {
    DRT::Problem::Instance()->AddFieldTest(CreateScaTraFieldTest());
    DRT::Problem::Instance()->TestAll(Comm());
  }

  return;
}
