#ifdef STKADAPTIVE

#include <limits>

#include <Teuchos_StandardParameterEntryValidators.hpp>
#include <Teuchos_Time.hpp>

#include <Epetra_CrsGraph.h>

#include "stk_algebra.H"
#include "stk_assemblestrategy.H"
#include "stk_errorestimate.H"
//#include "stk_fixedsparsematrix.H"
#include "stk_fluid.H"
#include "stk_fluidresulttest.H"
#include "stk_gmsh.H"
#include "../stk_lib/stk_gnuplot.H"
#include "../stk_lib/stk_iterator.H"
#include "../stk_lib/stk_mesh.H"
#include "../stk_lib/stk_types.H"
#include "../stk_lib/stk_utils.H"

#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_dserror.H"
#include "../drt_lib/drt_dirichletextractor.H"

#include "../drt_fluid/fluid_utils_mapextractor.H"

#include "../linalg/linalg_fixedsparsematrix.H"
#include "../linalg/linalg_solver.H"
#include "../linalg/linalg_sparsematrix.H"
#include "../linalg/linalg_utils.H"

#include "../drt_io/io.H"
#include "../drt_io/io_control.H"

#include "../drt_fsi/fsi_debugwriter.H"


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
STK::FluidDRT::FluidDRT( int counter, STK::Fluid & fluid )
  : counter_( counter )
{
  DRT::Discretization & dis = fluid.Discretization();
  STK::Mesh & mesh = fluid.GetMesh();

  // object holds maps/subsets for DOFs subjected to Dirichlet BCs and otherwise
  dbcmaps_ = Teuchos::rcp(new DRT::DirichletExtractor());
  dbcmaps_->Setup( dis );

  sysmat_ = Teuchos::rcp( new LINALG::FixedSparseMatrix( dbcmaps_->DirichletMap() ) );

  residual_ = Teuchos::rcp( new Epetra_Vector( *dis.DofRowMap() ) );

  assemblestrategy_ = Teuchos::rcp( new STK::AssembleStrategy( dis, mesh, dbcmaps_->DirichletMap(),
                                                               sysmat_, Teuchos::null,
                                                               residual_, Teuchos::null, Teuchos::null ) );

  sysmat_->SetMatrix( assemblestrategy_->MatrixGraph( dis, dbcmaps_->DirichletMap() ) );

  velpressplitter_ = Teuchos::rcp( new FLD::UTILS::VelPressExtractor() );
  velpressplitter_->Setup( dis );

  output_ = Teuchos::rcp( new IO::DiscretizationWriter( Teuchos::rcp( &dis, false ),
                                                        DRT::Problem::Instance()->OutputControlFile() ) );
  output_->WriteMesh( fluid.step_, fluid.time_ );
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void STK::FluidDRT::Output( STK::Fluid & fluid )
{
  if ( m_exo == Teuchos::null )
  {
    std::vector< const stk::mesh::FieldBase * > out_fields ;
    out_fields.push_back( fluid.velnp_ );
    out_fields.push_back( fluid.pressure_ );
    out_fields.push_back( fluid.constrained_ );
    out_fields.push_back( fluid.hanging_ );
    out_fields.push_back( fluid.error_ );

    std::string filename = DRT::Problem::Instance()->OutputControlFile()->FileName();
    std::stringstream str;
    str << filename << "-fluid-" << counter_ << ".exo";
    m_exo = fluid.GetMesh().OutputContext( str.str(), "Navier-Stokes Problem", out_fields );
  }

  m_exo->write( fluid.step_, fluid.time_ );

  output_->NewStep( fluid.step_, fluid.time_ );

  std::vector<stk::mesh::FieldBase*> v;
  v.push_back( fluid.velnp_ );
  v.push_back( fluid.pressure_ );

  Teuchos::RCP<Epetra_Vector> velnp = fluid.GatherFieldData( v );

  // velocity/pressure vector
  output_->WriteVector("velnp",velnp);

  // (hydrodynamic) pressure
  Teuchos::RCP<Epetra_Vector> pressure = velpressplitter_->ExtractPressureVector(velnp);
  //pressure->Scale(density_);
  output_->WriteVector("pressure", pressure);
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
STK::Fluid::Fluid( DRT::Discretization & dis,
                   Teuchos::RCP<LINALG::Solver> solver )
  : Adaptive( dis ),
    time_( 0.0 ),
    step_( 0 )
{
  const Teuchos::ParameterList& fdyn     = DRT::Problem::Instance()->FluidDynamicParams();

  physicaltype_ = Teuchos::getIntegralValue<INPAR::FLUID::PhysicalType>(fdyn,"PHYSICAL_TYPE");
  timealgo_     = Teuchos::getIntegralValue<FLUID_TIMEINTTYPE>(fdyn,"TIMEINTEGR");
  dyntype_      = Teuchos::getIntegralValue<int>(fdyn,"DYNAMICTYP");
  stepmax_      = fdyn.get<int>("NUMSTEP");
  maxtime_      = fdyn.get<double>("MAXTIME");
  dta_          = fdyn.get<double>("TIMESTEP");
  dtp_          = dta_;
  theta_        = fdyn.get<double>("THETA");
  alefluid_     = false;
  newton_       = fdyn.get<string>("NONLINITER");
  convform_     = fdyn.get<string>("CONVFORM");

  refinestep_   = 1;
  maxrefine_    = 5;

  max_error_    = 0.01;
  min_error_    = 0.001;
  min_volume_   = 1e-6;

  solver_ = solver;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void STK::Fluid::SetupDRTMesh()
{
  int counter = 0;
  if ( drt_ != Teuchos::null )
  {
    counter = drt_->Counter() + 1;
    drt_ = Teuchos::null;
  }
  drt_ = Teuchos::rcp( new FluidDRT( counter, *this ) );
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void STK::Fluid::DefineFields( stk::mesh::MetaData & meta_data )
{
  velnp_    = & declare_vector_field_on_all_nodes( meta_data, "velocity", genprob.ndim );
  veln_     = & declare_vector_field_on_all_nodes( meta_data, "veln",     genprob.ndim );
  velnm_    = & declare_vector_field_on_all_nodes( meta_data, "velnm",    genprob.ndim );
  accnp_    = & declare_vector_field_on_all_nodes( meta_data, "accnp",    genprob.ndim );
  accn_     = & declare_vector_field_on_all_nodes( meta_data, "accn",     genprob.ndim );
  hist_     = & declare_vector_field_on_all_nodes( meta_data, "hist",     genprob.ndim );
  pressure_ = & declare_scalar_field_on_all_nodes( meta_data, "pressure" );

  constrained_ = & declare_scalar_field_on_all_nodes( meta_data, "constrained" );
  hanging_     = & declare_scalar_field_on_all_nodes( meta_data, "hanging" );

  stk::mesh::Part & active = *meta_data.get_part( "active" );

  //error_    = & declare_scalar_field_on_all_elements( meta_data, "error" );
  error_    = & stk::mesh::put_field( meta_data.declare_field<stk::mesh::ScalarField>( "error" ),
                                      stk::mesh::Element,
                                      active );
  volume_   = & stk::mesh::put_field( meta_data.declare_field<stk::mesh::ScalarField>( "volume" ),
                                      stk::mesh::Element,
                                      active );
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void STK::Fluid::Integrate()
{
  stk::mesh::BulkData & bulk_data = GetMesh().BulkData();
  int p_rank = bulk_data.parallel_rank();

  while (step_<stepmax_ and time_<maxtime_)
  {
    PrepareTimeStep();
    // -------------------------------------------------------------------
    //                       output to screen
    // -------------------------------------------------------------------
    if (p_rank==0)
    {
      switch (timealgo_)
      {
      case timeint_one_step_theta:
        printf("TIME: %11.4E/%11.4E  DT = %11.4E   One-Step-Theta    STEP = %4d/%4d \n",time_,maxtime_,dta_,step_,stepmax_);
        break;
      case timeint_afgenalpha:
        printf("TIME: %11.4E/%11.4E  DT = %11.4E  Generalized-Alpha  STEP = %4d/%4d \n",time_,maxtime_,dta_,step_,stepmax_);
        break;
      case timeint_bdf2:
        printf("TIME: %11.4E/%11.4E  DT = %11.4E       BDF2          STEP = %4d/%4d \n",time_,maxtime_,dta_,step_,stepmax_);
        break;
      default:
        throw std::runtime_error("parameter out of range: IOP");
      }
    }

    switch (dyntype_)
    {
    case 0:
      // -----------------------------------------------------------------
      //                     solve nonlinear equation
      // -----------------------------------------------------------------
      //NonlinearSolve();
      AdaptiveNonlinearSolve();
      break;
    case 1:
      // -----------------------------------------------------------------
      //                     solve linearised equation
      // -----------------------------------------------------------------
      LinearSolve();
      break;
    default:
      dserror("type of dynamics unknown");
    }

    // -------------------------------------------------------------------
    //                         update solution
    //        current solution becomes old solution of next timestep
    // -------------------------------------------------------------------
    TimeUpdate();

    // -------------------------------------------------------------------
    //  lift'n'drag forces, statistics time sample and output of solution
    //  and statistics
    // -------------------------------------------------------------------
    StatisticsAndOutput();

    // -------------------------------------------------------------------
    //                       update time step sizes
    // -------------------------------------------------------------------
    dtp_ = dta_;
  }
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<DRT::ResultTest> STK::Fluid::CreateFieldTest()
{
  return Teuchos::rcp( new FluidResultTest( *this ) );
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void STK::Fluid::PrepareTimeStep()
{
  stk::mesh::BulkData & bulk_data = GetMesh().BulkData();

  // -------------------------------------------------------------------
  //              set time dependent parameters
  // -------------------------------------------------------------------
  step_ += 1;
  time_ += dta_;

  // for BDF2, theta is set by the time-step sizes, 2/3 for const. dt
  if (timealgo_==timeint_bdf2) theta_ = (dta_+dtp_)/(2.0*dta_ + dtp_);

  // -------------------------------------------------------------------
  // set part(s) of the rhs vector(s) belonging to the old timestep
  // (only meaningful for momentum part)
  //
  // stationary/af-generalized-alpha: hist_ = 0.0
  //
  // one-step-Theta:                  hist_ = veln_  + dt*(1-Theta)*accn_
  //
  // BDF2: for constant time step:    hist_ = 4/3 veln_  - 1/3 velnm_
  //
  // -------------------------------------------------------------------

  switch ( timealgo_ )
  {
  case timeint_stationary: /* Stationary algorithm */
  case timeint_afgenalpha: /* Af-generalized-alpha time integration */
    algebra::PutScalar( bulk_data, *hist_, 0.0 );
    break;

  case timeint_one_step_theta: /* One step Theta time integration */
    algebra::Update( bulk_data, *hist_, 0.0, *veln_, 1.0, *accn_, dta_*(1.0-theta_) );
    break;

  case timeint_bdf2:    /* 2nd order backward differencing BDF2 */
    algebra::Update( bulk_data, *hist_, 0.0, *veln_, 4./3., *velnm_, -1./3. );
    break;

  default:
    throw std::runtime_error("Time integration scheme unknown!");
  }


  // -------------------------------------------------------------------
  //                     do explicit predictor step
  //
  //                      +-                                      -+
  //                      | /     dta \          dta  veln_-velnm_ |
  // velnp_ = veln_ + dta | | 1 + --- | accn_ - ----- ------------ |
  //                      | \     dtp /          dtp     dtp       |
  //                      +-                                      -+
  //
  // -------------------------------------------------------------------
  //
  if ( step_>1 )
  {
    switch ( timealgo_ )
    {
    case timeint_stationary: /* Stationary algorithm */
      // do nothing
      break;
    case timeint_afgenalpha: /* Generalized-alpha time integration */
    {
      // do nothing for the time being, that is, steady-state predictor
      break;
    }
    case timeint_one_step_theta: /* One step Theta time integration */
    case timeint_bdf2:    /* 2nd order backward differencing BDF2 */
    {
      const double fact1 = dta_*(1.0+dta_/dtp_);
      const double fact2 = (dta_/dtp_)*(dta_/dtp_);

      algebra::Update( bulk_data, *velnp_, 1.0, *accn_, fact1, *veln_, -fact2, *velnm_, fact2 );
      break;
    }
    default:
      throw std::runtime_error("Time integration scheme unknown!");
    }
  }

  // -------------------------------------------------------------------
  //  evaluate Dirichlet and Neumann boundary conditions
  // -------------------------------------------------------------------
  {
    Teuchos::ParameterList eleparams;
    DRT::Discretization & discret = Discretization();

    std::vector<stk::mesh::FieldBase*> v;
    v.push_back( velnp_ );
    v.push_back( pressure_ );

    Teuchos::RCP<Epetra_Vector> velnp = GatherFieldData( v );

    // total time required for Dirichlet conditions
    eleparams.set("total time",time_);

    // set vector values needed by elements
    discret.ClearState();
    discret.SetState("velnp",velnp);

    // predicted dirichlet values
    // velnp then also holds prescribed new dirichlet values
    discret.EvaluateDirichlet(eleparams,velnp,null,null,null);
    discret.ClearState();

    ScatterFieldData( velnp, v );

#if 0
    // set all parameters and states required for Neumann conditions
    eleparams.set("Physical Type",physicaltype_);
    eleparams.set("thermpress at n+1",thermpressaf_);
    if (timealgo_==timeint_afgenalpha)
    {
      eleparams.set("total time",time_-(1-alphaF_)*dta_);
      eleparams.set("thsl",1.0);
    }
    else
    {
      eleparams.set("total time",time_);
      eleparams.set("thsl",theta_*dta_);
    }

    // evaluate Neumann conditions
    neumann_loads_->PutScalar(0.0);
    discret_->SetState("scanp",scaaf_);
    discret_->EvaluateNeumann(eleparams,*neumann_loads_);
    discret_->ClearState();
#endif
  }

  // -------------------------------------------------------------------
  //  For af-generalized-alpha time-integration scheme:
  //  set "pseudo-theta", calculate initial accelerations according to
  //  prescribed Dirichlet values for generalized-alpha time
  //  integration and values at intermediate time steps
  // -------------------------------------------------------------------
  if (timealgo_==timeint_afgenalpha)
  {
    throw std::runtime_error( "not supported" );
  }
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void STK::Fluid::LinearSolve()
{
  stk::mesh::BulkData & bulk_data = GetMesh().BulkData();

  // -------------------------------------------------------------------
  // call elements to calculate system matrix
  // -------------------------------------------------------------------

  // get cpu time
  const double tcpuele = Teuchos::Time::wallTime();

  // create the parameters for the discretization
  Teuchos::ParameterList eleparams;
  DRT::Discretization & discret = Discretization();

  std::vector<stk::mesh::FieldBase*> v;
  v.push_back( velnp_ );
  v.push_back( pressure_ );

  Teuchos::RCP<Epetra_Vector> velnp = GatherFieldData( v );

  v[0] = accnp_;
  Teuchos::RCP<Epetra_Vector> accnp = GatherFieldData( v );

  v[0] = hist_;
  Teuchos::RCP<Epetra_Vector> hist  = GatherFieldData( v );

  Teuchos::RCP<Epetra_Vector> rhs   = Teuchos::rcp( new Epetra_Vector( *discret.DofRowMap() ) );

  // action for elements
  eleparams.set("action","calc_linear_fluid");

  // other parameters that might be needed by the elements
  eleparams.set("total time",time_);
  eleparams.set("thsl",theta_*dta_);
  eleparams.set("Physical Type",physicaltype_);

//   eleparams.set("thermpress at n+alpha_F/n+1",thermpressaf_);
//   eleparams.set("thermpress at n+alpha_M/n",thermpressam_);
//   eleparams.set("thermpressderiv at n+alpha_M/n+1",thermpressdtam_);

  // set vector values needed by elements
  discret.ClearState();
  discret.SetState("velaf",velnp);
  //discret.SetState("scaaf",scaaf);
  discret.SetState("accam",accnp);
  discret.SetState("hist" ,hist );

  // call standard loop over linear elements
  discret.Evaluate( eleparams, drt_->sysmat_, Teuchos::null, rhs, Teuchos::null, Teuchos::null );
  discret.ClearState();

  // finalize the complete matrix
  drt_->sysmat_->Complete();

  // end time measurement for element
  const double dtele = Teuchos::Time::wallTime() - tcpuele;

  //-------solve for total new velocities and pressures
  // get cpu time
  const double tcpusolve = Teuchos::Time::wallTime();

  LINALG::ApplyDirichlettoSystem( drt_->sysmat_, velnp, rhs, velnp,
                                  *drt_->dbcmaps_->DirichletMap() );

  solver_->Solve( drt_->sysmat_->EpetraOperator(), velnp, rhs, true, true );

  v[0] = velnp_;
  ScatterFieldData( velnp, v );

  // end time measurement for solver
  const double dtsolve = Teuchos::Time::wallTime() - tcpusolve;

  if ( bulk_data.parallel_rank() == 0 )
    std::cout << "te=" << dtele << ", ts=" << dtsolve << "\n\n" ;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void STK::Fluid::NonlinearSolve()
{
  stk::mesh::BulkData & bulk_data = GetMesh().BulkData();
  int p_rank = bulk_data.parallel_rank();

  DRT::Discretization & discret = Discretization();

  const Teuchos::ParameterList& fdyn = DRT::Problem::Instance()->FluidDynamicParams();

  const double  ittol  = fdyn.get<double>("CONVTOL");
  const int     itemax = fdyn.get<int>("ITEMAX");

  int  itnum = 0;
  bool stopnonliniter = false;

  double dtsolve = 0;
  double dtele = 0;

  if (p_rank == 0)
  {
    printf("+------------+-------------------+--------------+--------------+--------------+--------------+\n");
    printf("|- step/max -|- tol      [norm] -|-- vel-res ---|-- pre-res ---|-- vel-inc ---|-- pre-inc ---|\n");
  }

  std::vector<stk::mesh::FieldBase*> v;
  v.push_back( velnp_ );
  v.push_back( pressure_ );

  Teuchos::RCP<Epetra_Vector> velnp = GatherFieldData( v );

  v[0] = accnp_;
  Teuchos::RCP<Epetra_Vector> accnp = GatherFieldData( v );

  v[0] = hist_;
  Teuchos::RCP<Epetra_Vector> hist  = GatherFieldData( v );

  //Teuchos::RCP<Epetra_Vector> residual = Teuchos::rcp( new Epetra_Vector( *discret.DofRowMap() ) );

  Teuchos::RCP<Epetra_Vector> incvel = Teuchos::rcp( new Epetra_Vector( *discret.DofRowMap() ) );
  incvel->PutScalar( 0.0 );

  Teuchos::RCP<Epetra_Vector> zeros = Teuchos::rcp( new Epetra_Vector( *discret.DofRowMap() ) );
  zeros->PutScalar( 0.0 );

// #ifdef DEBUG

//   std::vector< const stk::mesh::FieldBase * > out_fields ;
//   out_fields.push_back( velnp_ );
//   out_fields.push_back( pressure_ );

//   std::stringstream s;
//   s << "fluid_" << step_ << ".exo";

//   Teuchos::RCP<phdmesh::exodus::FileOutput> exo = GetMesh().OutputContext( s.str(), "Navier-Stokes Problem", out_fields );

//   FSI::UTILS::SimpleDebugWriter sdw( Teuchos::rcp( &Discretization(), false ), "fluid" );
//   sdw.NewLinearSystem( step_, "fluid" );

// #endif

  while (stopnonliniter==false)
  {
    itnum++;

    // get cpu time
    const double tcpu=Teuchos::Time::wallTime();

    drt_->sysmat_->Zero();

    // create the parameters for the discretization
    Teuchos::ParameterList eleparams;

    drt_->residual_->PutScalar( 0.0 );

    // add Neumann loads
    //residual_->Update(1.0,*neumann_loads_,0.0);

    // set general element parameters
    eleparams.set("dt",dta_);
    eleparams.set("theta",theta_);
//     eleparams.set("omtheta",omtheta_);
    eleparams.set("form of convective term",convform_);
//     eleparams.set("fs subgrid viscosity",fssgv_);
    eleparams.set("Linearisation", newton_);
    eleparams.set("Physical Type", physicaltype_);
//     eleparams.set("mixed_formulation", params_.get<bool>("mixed_formulation", false));

    // parameters for stabilization
    eleparams.sublist("STABILIZATION") = fdyn.sublist("STABILIZATION");

    // parameters for stabilization
    eleparams.sublist("TURBULENCE MODEL") = fdyn.sublist("TURBULENCE MODEL");

    eleparams.set("thermpress at n+alpha_F/n+1",0.0);
    eleparams.set("thermpress at n+alpha_M/n",0.0);
    eleparams.set("thermpressderiv at n+alpha_M/n+1",0.0);

    // set general vector values needed by elements
    discret.ClearState();
    discret.SetState("hist" ,hist );
    discret.SetState("accam",accnp);
    discret.SetState("scaaf",zeros);
    discret.SetState("scaam",zeros);
    if (alefluid_)
    {
//       discret_->SetState("dispnp", dispnp_);
//       discret_->SetState("gridv", gridv_);
    }

    // set scheme-specific element parameters and vector values
    if (timealgo_==timeint_stationary)
    {
      eleparams.set("action","calc_fluid_stationary_systemmat_and_residual");
      eleparams.set("using generalized-alpha time integration",false);
      eleparams.set("total time",time_);
      eleparams.set("is stationary", true);

      discret.SetState("velaf",velnp);
    }
    else if (timealgo_==timeint_afgenalpha)
    {
//       eleparams.set("action","calc_fluid_afgenalpha_systemmat_and_residual");
//       eleparams.set("using generalized-alpha time integration",true);
//       eleparams.set("total time",time_-(1-alphaF_)*dta_);
//       eleparams.set("is stationary", false);
//       eleparams.set("alphaF",alphaF_);
//       eleparams.set("alphaM",alphaM_);
//       eleparams.set("gamma",gamma_);

      //discret_->SetState("velaf",velaf_);
    }
    else
    {
      eleparams.set("action","calc_fluid_systemmat_and_residual");
      eleparams.set("using generalized-alpha time integration",false);
      eleparams.set("total time",time_);
      eleparams.set("is stationary", false);

      discret.SetState("velaf",velnp);
    }

    // call standard loop over elements
    discret.Evaluate(eleparams,*drt_->assemblestrategy_);

    discret.ClearState();

    // finalize the complete matrix
    drt_->sysmat_->Complete();

    // end time measurement for element
    dtele = Teuchos::Time::wallTime()-tcpu;

    drt_->dbcmaps_->ZeroDirichlets( drt_->residual_ );

    double incvelnorm_L2 = drt_->velpressplitter_->VelocityNorm2( *incvel );
    double incprenorm_L2 = drt_->velpressplitter_->PressureNorm2( *incvel );

    double velnorm_L2 = algebra::Norm2( bulk_data, *velnp_ );
    double prenorm_L2 = algebra::Norm2( bulk_data, *pressure_ );

    double vresnorm = drt_->velpressplitter_->VelocityNorm2( *drt_->residual_ );
    double presnorm = drt_->velpressplitter_->PressureNorm2( *drt_->residual_ );

    // care for the case that nothing really happens in the velocity
    // or pressure field
    if (velnorm_L2 < 1e-5) velnorm_L2 = 1.0;
    if (prenorm_L2 < 1e-5) prenorm_L2 = 1.0;

    //-------------------------------------------------- output to screen
    /* special case of very first iteration step:
        - solution increment is not yet available
        - convergence check is not required (we solve at least once!)    */
    if (itnum == 1)
    {
      if (p_rank == 0)
      {
        printf("|  %3d/%3d   | %10.3E[L_2 ]  | %10.3E   | %10.3E   |      --      |      --      |",
               itnum,itemax,ittol,vresnorm,presnorm);
        printf(" (      --     ,te=%10.3E",dtele);
        printf(")\n");
      }
    }
    /* ordinary case later iteration steps:
        - solution increment can be printed
        - convergence check should be done*/
    else
    {
      // this is the convergence check
      // We always require at least one solve. Otherwise the
      // perturbation at the FSI interface might get by unnoticed.
      if (vresnorm <= ittol and presnorm <= ittol and
          incvelnorm_L2/velnorm_L2 <= ittol and incprenorm_L2/prenorm_L2 <= ittol)
      {
        stopnonliniter=true;
        if (p_rank == 0)
        {
          printf("|  %3d/%3d   | %10.3E[L_2 ]  | %10.3E   | %10.3E   | %10.3E   | %10.3E   |",
                 itnum,itemax,ittol,vresnorm,presnorm,
                 incvelnorm_L2/velnorm_L2,incprenorm_L2/prenorm_L2);
          printf(" (ts=%10.3E,te=%10.3E",dtsolve,dtele);
          printf(")\n");
          printf("+------------+-------------------+--------------+--------------+--------------+--------------+\n");
        }
        break;
      }
      else
      {
        // if not yet converged
        if (p_rank == 0)
        {
          printf("|  %3d/%3d   | %10.3E[L_2 ]  | %10.3E   | %10.3E   | %10.3E   | %10.3E   |",
                 itnum,itemax,ittol,vresnorm,presnorm,
                 incvelnorm_L2/velnorm_L2,incprenorm_L2/prenorm_L2);
          printf(" (ts=%10.3E,te=%10.3E",dtsolve,dtele);
          printf(")\n");
        }
      }
    }

    // warn if itemax is reached without convergence, but proceed to
    // next timestep...
    if ((itnum == itemax) and (vresnorm > ittol or presnorm > ittol or
                               incvelnorm_L2/velnorm_L2 > ittol or
                               incprenorm_L2/prenorm_L2 > ittol))
    {
      stopnonliniter=true;
      if (p_rank == 0)
      {
        printf("+---------------------------------------------------------------+\n");
        printf("|            >>>>>> not converged in itemax steps!              |\n");
        printf("+---------------------------------------------------------------+\n");
      }
      break;
    }

    //--------- Apply Dirichlet boundary conditions to system of equations
    //          residual displacements are supposed to be zero at
    //          boundary conditions
    if ( itnum==1 )
    {
      incvel->Update( 1.0, *velnp, 0.0 );
      LINALG::ApplyDirichlettoSystem( drt_->sysmat_, incvel, drt_->residual_, incvel,
                                      *drt_->dbcmaps_->DirichletMap() );
    }
    else
    {
      incvel->PutScalar( 0.0 );
      LINALG::ApplyDirichlettoSystem( drt_->sysmat_, incvel, drt_->residual_, zeros,
                                      *drt_->dbcmaps_->DirichletMap() );
    }

    //sysmat_->EpetraMatrix()->Print( std::cout );
    //residual_->Print( std::cout );

    const double tcpusolve=Teuchos::Time::wallTime();

    solver_->Solve( drt_->sysmat_->EpetraOperator(), incvel, drt_->residual_, true, itnum==1 );
    //solver_->ResetTolerance();

    // end time measurement for solver
    dtsolve = Teuchos::Time::wallTime()-tcpusolve;

// #ifdef DEBUG
//     sdw.NewIteration();
//     sdw.WriteVector( "velnp", *velnp );
//     sdw.WriteVector( "residual", *residual_ );
//     sdw.WriteVector( "incvel", *incvel );
// #endif

    // -------------------------------------------------------------------
    // update velocity and pressure values by increments
    // -------------------------------------------------------------------
    if ( itnum==1 )
    {
      incvel->Update( -1.0, *velnp, 1.0 );
      velnp->Update( 1.0, *incvel, 1.0 );
    }
    else
    {
      velnp->Update( 1.0, *incvel, 1.0 );
    }

    v[0] = velnp_;
    ScatterFieldData( velnp, v );

// #ifdef DEBUG
//     exo->write( itnum, itnum );
// #endif
  }

// #ifdef DEBUG
//   std::stringstream gs;
//   gs << "fluid_" << step_ << ".gp";

//   GnuplotDump( GetMesh(), gs.str(), *velnp_ );

//   DumpGmsh( Discretization(), "fluid", *velnp, step_ );
// #endif
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void STK::Fluid::TimeUpdate()
{
  stk::mesh::BulkData & bulk_data = GetMesh().BulkData();

  if ( step_ == 1 )
  {
    switch ( timealgo_ )
    {
      case timeint_stationary: /* no accelerations for stationary problems*/
      {
        algebra::PutScalar( bulk_data, *accnp_, 0.0 );
        break;
      }
      case timeint_one_step_theta: /* One step Theta time integration */
      case timeint_bdf2:    /* 2nd order backward differencing BDF2 */
      {
        // do just a linear interpolation within the first timestep
        algebra::Update( bulk_data, *accnp_, 0.0, *velnp_, 1.0/dta_, *veln_, -1.0/dta_ );
        break;
      }
      case timeint_afgenalpha: /* Af-generalized-alpha time integration */
      {
        // startup is handled separately
        break;
      }
      default:
        throw std::runtime_error("Time integration scheme unknown!");
    }
  }
  else
  {
    /*

    Following formulations are for n+1; acceleration values, however, are
    directly stored in vectors at time n (velocity has not yet been updated).

    One-step-Theta:

     acc(n+1) = (vel(n+1)-vel(n)) / (Theta * dt(n)) - (1/Theta -1) * acc(n)


    BDF2:

                   2*dt(n)+dt(n-1)                  dt(n)+dt(n-1)
     acc(n+1) = --------------------- vel(n+1) - --------------- vel(n)
                 dt(n)*[dt(n)+dt(n-1)]              dt(n)*dt(n-1)

                         dt(n)
               + ----------------------- vel(n-1)
                 dt(n-1)*[dt(n)+dt(n-1)]

    */

    switch ( timealgo_ )
    {
      case timeint_stationary: /* no accelerations for stationary problems*/
      {
        algebra::PutScalar( bulk_data, *accnp_, 0.0 );
        break;
      }
      case timeint_one_step_theta: /* One-step-theta time integration */
      {
        const double fact1 = 1.0/(theta_*dta_);
        const double fact2 =-1.0/theta_ +1.0;   /* = -1/Theta + 1 */

        algebra::Update( bulk_data, *accnp_, 0.0, *velnp_, fact1, *veln_, -fact1, *accn_, fact2 );
        break;
      }
      case timeint_bdf2:    /* 2nd order backward differencing BDF2 */
      {
        if (dta_*dtp_ < 1e-15) throw std::runtime_error("zero time step size");
        const double sum = dta_ + dtp_;

        algebra::Update( bulk_data,
                         *accnp_, 0.0,
                         *velnp_, (2.0*dta_+dtp_)/(dta_*sum),
                         *veln_, -sum/(dta_*dtp_),
                         *velnm_, dta_/(dtp_*sum) );
        break;
      }
      case timeint_afgenalpha: /* Af-generalized-alpha time integration */
      {
        // do nothing: new acceleration is calculated at beginning of next time step
        break;
      }
      default:
        throw std::runtime_error("Time integration scheme unknown!");
    }
  }

  // update old acceleration
  algebra::Update( bulk_data, *accn_, 0.0, *accnp_, 1.0 );

  // velocities/pressures of this step become most recent
  // velocities/pressures of the last step
  algebra::Update( bulk_data, *velnm_, 0.0, *veln_ , 1.0 );
  algebra::Update( bulk_data, *veln_ , 0.0, *velnp_, 1.0 );
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void STK::Fluid::StatisticsAndOutput()
{
  stk::mesh::BulkData & bulk = GetMesh().BulkData();

  // mark current constraints

  algebra::PutScalar( bulk, *constrained_, 0.0 );
  algebra::PutScalar( bulk, *hanging_, 0.0 );

  const std::vector<stk::mesh::Bucket*> & constraints = bulk.buckets( stk::mesh::Constraint );
  for ( std::vector<stk::mesh::Bucket*>::const_iterator i=constraints.begin();
        i!=constraints.end();
        ++i )
  {
    stk::mesh::Bucket & bucket = **i;

    for ( stk::mesh::Bucket::iterator j=bucket.begin();
          j!=bucket.end();
          ++j )
    {
      stk::mesh::Entity & c = *j;

      stk::mesh::PairIterRelation rel = c.relations( stk::mesh::Node );

      double & hnv = * stk::mesh::field_data( *hanging_, *rel->entity() );
      hnv += 1;

      for ( ; not rel.empty(); ++rel )
      {
        stk::mesh::Entity & n = *rel->entity();
        double & v = * stk::mesh::field_data( *constrained_, n );
        v += 1;
      }
    }
  }

  drt_->Output( *this );
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void STK::Fluid::AdaptiveNonlinearSolve()
{
  NonlinearSolve();

  if ( step_ % refinestep_ == 0 )
  {
    std::vector<stk::mesh::EntityKey> refine;
    std::vector<stk::mesh::EntityKey> unrefine;

    ErrorEstimate( refine, unrefine );

    int r = 0;
    for ( ; r < maxrefine_; ++r )
    {
      //if ( refine.size()==0 and unrefine.size()==0 )
      if ( refine.size()==0 )
        break;

      Refine( refine );
      Unrefine( unrefine );

      SetupDRTMesh();

      NonlinearSolve();

      refine.clear();
      unrefine.clear();
      ErrorEstimate( refine, unrefine );
    }
    if ( refine.size()!=0 or unrefine.size()!=0 )
    {
      if ( r==maxrefine_ )
      {
        stk::mesh::BulkData & bulk = GetMesh().BulkData();
        if ( bulk.parallel_rank()==0 )
        {
          std::cout << "Warning: failed to reduce error\n";
        }
      }
    }
  }
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void STK::Fluid::ErrorEstimate( std::vector<stk::mesh::EntityKey> & refine,
                                std::vector<stk::mesh::EntityKey> & unrefine )
{
  stk::mesh::BulkData & bulk = GetMesh().BulkData();
  stk::mesh::Part & active = GetMesh().ActivePart();

  //PrintMesh( std::cout, bulk, active );

  algebra::PutScalar( bulk, *error_ , 0.0 );
  algebra::PutScalar( bulk, *volume_, 0.0 );

  FluidJumpIntegrator integrator( Discretization(), bulk, active );
  integrator.Integrate( GetMesh().Coordinates(), *velnp_, *pressure_, *error_, *volume_, drt_->residual_ );

  double min_err = std::numeric_limits<double>::max();
  double max_err = std::numeric_limits<double>::min();

  for ( ElementIterator i( bulk, active );
        not i.done();
        ++i )
  {
    stk::mesh::Entity & e = i.element();
    double & err = *reinterpret_cast<double*>( stk::mesh::field_data( *error_ , e ) );

    min_err = std::min( err, min_err );
    max_err = std::max( err, max_err );
  }

  stk::all_reduce( bulk.parallel(), stk::ReduceMin<1>( &min_err ) );
  stk::all_reduce( bulk.parallel(), stk::ReduceMax<1>( &max_err ) );

//   double omega = 0.5;
//   double med_err = min_err + omega*( max_err-min_err );

  if ( bulk.parallel_rank()==0 )
  {
    std::cout << "    max error = " << max_err
              << "    min error = " << min_err
//               << "    med error = " << med_err
              << "\n";
  }

  unsigned passed = 0;

  for ( ElementIterator i( bulk, active );
        not i.done();
        ++i )
  {
    stk::mesh::Entity & e = i.element();
    double & err = *reinterpret_cast<double*>( stk::mesh::field_data( *error_ , e ) );

    if ( err > max_error_ )
    {
      double & vol = *reinterpret_cast<double*>( stk::mesh::field_data( *volume_ , e ) );
      if ( vol > min_volume_ )
      {
        refine.push_back( e.key() );
      }
      else
      {
        passed += 1;
      }
    }
    else if ( err < min_error_ )
    {
      unrefine.push_back( e.key() );
    }
  }

  if ( bulk.parallel_rank()==0 )
  {
    std::cout << "    #refine elements: " << refine.size()
              << "    #unrefine elements: " << unrefine.size()
              << "    #passed elements: " << passed
              << "\n";
  }
}

#endif
