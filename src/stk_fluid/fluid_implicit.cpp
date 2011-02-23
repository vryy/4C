
#ifdef STKADAPTIVE

#include <Teuchos_StandardParameterEntryValidators.hpp>
#include <Teuchos_Time.hpp>

#include "fluid_implicit.H"
#include "fluid_errorestimate.H"
#include "fluid_resulttest.H"

#include "../stk_lib/stk_discret.H"
#include "../stk_refine/stk_mesh.H"
#include "../stk_lib/stk_iterator.H"

#include "../stk_lib/stk_algebra.H"

#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/standardtypes_cpp.H"
#include "../drt_lib/drt_element.H"

#include "../linalg/linalg_fixedsparsematrix.H"
#include "../linalg/linalg_solver.H"
#include "../linalg/linalg_sparsematrix.H"
#include "../linalg/linalg_utils.H"

#include "../drt_io/io_control.H"

#include "../drt_f3_impl/fluid3_impl.H"

#include "../drt_mat/matpar_parameter.H"
#include "../drt_mat/material.H"


  /*----------------------------------------------------------------------*
   *----------------------------------------------------------------------*/
  template <int nsd, int nen>
  void ExtractValues( stk::mesh::Entity & e,
                      const stk::mesh::VectorField * v,
                      typename LINALG::Matrix<nsd,nen> & ev )
  {
    stk::mesh::PairIterRelation nodes = e.relations( stk::mesh::Node );
    unsigned size = nodes.size();
    for ( unsigned i=0; i<size; ++i )
    {
      stk::mesh::Entity & n = * nodes[i].entity();
      double * data = field_data( *v, n );
      for ( int j=0; j<nsd; ++j )
      {
        ev( j, i ) = data[j];
      }
    }
  }


  /*----------------------------------------------------------------------*
   *----------------------------------------------------------------------*/
  template <int nen>
  void ExtractValues( stk::mesh::Entity & e,
                      const stk::mesh::ScalarField * s,
                      typename LINALG::Matrix<nen, 1> & es )
  {
    stk::mesh::PairIterRelation nodes = e.relations( stk::mesh::Node );
    unsigned size = nodes.size();
    for ( unsigned i=0; i<size; ++i )
    {
      stk::mesh::Entity & n = * nodes[i].entity();
      double * data = field_data( *s, n );
      es( i, 0 ) = *data;
    }
  }


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
STK::FLD::FluidState::FluidState( int counter, Fluid & fluid )
  : counter_( counter ),
    system_( fluid.Discretization() )
{
//   STK::Discretization & dis = fluid.Discretization();
//   const DRT::DirichletExtractor & extractor = dis.DirichletExtractor();
//   const Teuchos::RCP<const Epetra_Map> & dirichletmap = extractor.DirichletMap();

//   sysmat_ = Teuchos::rcp( new LINALG::FixedSparseMatrix( dirichletmap ) );

//   residual_ = Teuchos::rcp( new Epetra_Vector( dis.DofRowMap() ) );

//   assemblestrategy_ = Teuchos::rcp( new STK::FEI::AssembleStrategy( dis,
//                                                                     sysmat_, Teuchos::null,
//                                                                     residual_, Teuchos::null, Teuchos::null ) );

//   sysmat_->SetMatrix( assemblestrategy_->MatrixGraph( dis, dirichletmap ) );
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void STK::FLD::FluidState::Output( STK::FLD::Fluid & fluid )
{
  if ( m_exo==Teuchos::null )
  {
    std::vector< const stk::mesh::FieldBase * > out_fields ;
    out_fields.push_back( fluid.velnp_ );
    out_fields.push_back( fluid.pressure_ );
    out_fields.push_back( fluid.constrained_ );
    out_fields.push_back( fluid.hanging_ );
    out_fields.push_back( fluid.error_ );

    STK::Discretization & dis = fluid.Discretization();

    std::string filename = DRT::Problem::Instance()->OutputControlFile()->FileName();
    std::stringstream str;
    str << filename << "-fluid-" << counter_ << ".exo";
    m_exo = dis.GetMesh().OutputContext( str.str(), "Navier-Stokes Problem", out_fields );
  }
  m_exo->write( fluid.step_, fluid.time_ );
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
STK::FLD::Fluid::Fluid( STK::Discretization & dis, Teuchos::RCP<LINALG::Solver> solver )
  : dis_( dis ),
    time_( 0.0 ),
    step_( 0 )
{
  const Teuchos::ParameterList& fdyn     = DRT::Problem::Instance()->FluidDynamicParams();

  physicaltype_ = DRT::INPUT::IntegralValue<INPAR::FLUID::PhysicalType>(fdyn,"PHYSICAL_TYPE");
  timealgo_     = DRT::INPUT::IntegralValue<INPAR::FLUID::TimeIntegrationScheme>(fdyn,"TIMEINTEGR");
  stepmax_      = fdyn.get<int>("NUMSTEP");
  maxtime_      = fdyn.get<double>("MAXTIME");
  dta_          = fdyn.get<double>("TIMESTEP");
  dtp_          = dta_;
  theta_        = fdyn.get<double>("THETA");
  alefluid_     = false;
  newton_       = DRT::INPUT::IntegralValue<INPAR::FLUID::LinearisationAction>(fdyn,"NONLINITER");
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
Teuchos::RCP<DRT::ResultTest> STK::FLD::Fluid::CreateFieldTest()
{
  return Teuchos::rcp( new STK::FLD::FluidResultTest( *this ) );
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void STK::FLD::Fluid::declare_fields( stk::mesh::MetaData & meta )
{
  incvel_   = & declare_vector_field_on_all_nodes( meta, "dv",       genprob.ndim );
  resvel_   = & declare_vector_field_on_all_nodes( meta, "rv",       genprob.ndim );
  velnp_    = & declare_vector_field_on_all_nodes( meta, "velocity", genprob.ndim );
  veln_     = & declare_vector_field_on_all_nodes( meta, "veln",     genprob.ndim );
  velnm_    = & declare_vector_field_on_all_nodes( meta, "velnm",    genprob.ndim );
  accnp_    = & declare_vector_field_on_all_nodes( meta, "accnp",    genprob.ndim );
  accn_     = & declare_vector_field_on_all_nodes( meta, "accn",     genprob.ndim );
  hist_     = & declare_vector_field_on_all_nodes( meta, "hist",     genprob.ndim );

  incpres_  = & declare_scalar_field_on_all_nodes( meta, "dp" );
  respres_  = & declare_scalar_field_on_all_nodes( meta, "rp" );
  pressure_ = & declare_scalar_field_on_all_nodes( meta, "pressure" );

  constrained_ = & declare_scalar_field_on_all_nodes( meta, "constrained" );
  hanging_     = & declare_scalar_field_on_all_nodes( meta, "hanging" );

  stk::mesh::Part & active = *meta.get_part( "active" );

  //error_    = & declare_scalar_field_on_all_elements( meta, "error" );
  error_    = & stk::mesh::put_field( meta.declare_field<stk::mesh::ScalarField>( "error" ),
                                      stk::mesh::Element,
                                      active );
  volume_   = & stk::mesh::put_field( meta.declare_field<stk::mesh::ScalarField>( "volume" ),
                                      stk::mesh::Element,
                                      active );
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void STK::FLD::Fluid::collect_unknowns( std::vector<stk::mesh::FieldBase*> & fields )
{
  fields.resize( 2 );
  fields[0] = incvel_;
  fields[1] = incpres_;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void STK::FLD::Fluid::notify_state_changed()
{
  int counter = 0;
  if ( state_ != Teuchos::null )
  {
    counter = state_->Counter() + 1;
    state_ = Teuchos::null;
  }
  state_ = Teuchos::rcp( new FluidState( counter, *this ) );

  if ( Discretization().GetMesh().parallel_rank()==0 )
  {
    std::cout << "mesh number " << counter << "\n";
  }
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void STK::FLD::Fluid::Integrate()
{
  stk::mesh::BulkData & bulk = dis_.GetMesh().BulkData();
  int p_rank = bulk.parallel_rank();

  StatisticsAndOutput();

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
      case INPAR::FLUID::timeint_one_step_theta:
        printf("TIME: %11.4E/%11.4E  DT = %11.4E   One-Step-Theta    STEP = %4d/%4d \n",time_,maxtime_,dta_,step_,stepmax_);
        break;
      case INPAR::FLUID::timeint_afgenalpha:
        printf("TIME: %11.4E/%11.4E  DT = %11.4E  Generalized-Alpha  STEP = %4d/%4d \n",time_,maxtime_,dta_,step_,stepmax_);
        break;
      case INPAR::FLUID::timeint_bdf2:
        printf("TIME: %11.4E/%11.4E  DT = %11.4E       BDF2          STEP = %4d/%4d \n",time_,maxtime_,dta_,step_,stepmax_);
        break;
      default:
        throw std::runtime_error("parameter out of range: IOP");
      }
    }

    // -----------------------------------------------------------------
    //                     solve nonlinear equation
    // -----------------------------------------------------------------
    //NonlinearSolve();
    AdaptiveNonlinearSolve();

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
void STK::FLD::Fluid::PrepareTimeStep()
{
  stk::mesh::BulkData & bulk = dis_.GetMesh().BulkData();

  // -------------------------------------------------------------------
  //              set time dependent parameters
  // -------------------------------------------------------------------
  step_ += 1;
  time_ += dta_;

  // for BDF2, theta is set by the time-step sizes, 2/3 for const. dt
  if (timealgo_==INPAR::FLUID::timeint_bdf2) theta_ = (dta_+dtp_)/(2.0*dta_ + dtp_);

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
  case INPAR::FLUID::timeint_stationary: /* Stationary algorithm */
  case INPAR::FLUID::timeint_afgenalpha: /* Af-generalized-alpha time integration */
    algebra::PutScalar( bulk, *hist_, 0.0 );
    break;

  case INPAR::FLUID::timeint_one_step_theta: /* One step Theta time integration */
    algebra::Update( bulk, *hist_, 0.0, *veln_, 1.0, *accn_, dta_*(1.0-theta_) );
    break;

  case INPAR::FLUID::timeint_bdf2:    /* 2nd order backward differencing BDF2 */
    algebra::Update( bulk, *hist_, 0.0, *veln_, 4./3., *velnm_, -1./3. );
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
    case INPAR::FLUID::timeint_stationary: /* Stationary algorithm */
      // do nothing
      break;
    case INPAR::FLUID::timeint_afgenalpha: /* Generalized-alpha time integration */
    {
      // do nothing for the time being, that is, steady-state predictor
      break;
    }
    case INPAR::FLUID::timeint_one_step_theta: /* One step Theta time integration */
    case INPAR::FLUID::timeint_bdf2:    /* 2nd order backward differencing BDF2 */
    {
      const double fact1 = dta_*(1.0+dta_/dtp_);
      const double fact2 = (dta_/dtp_)*(dta_/dtp_);

      algebra::Update( bulk, *velnp_, 1.0, *accn_, fact1, *veln_, -fact2, *velnm_, fact2 );
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
    std::vector<stk::mesh::FieldBase*> v;
    v.push_back( velnp_ );
    v.push_back( pressure_ );

    dis_.EvaluateDirichlet( time_, &v );

#if 0
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
    eleparams.set<int>("Physical Type",physicaltype_);
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
#endif
  }

  // -------------------------------------------------------------------
  //  For af-generalized-alpha time-integration scheme:
  //  set "pseudo-theta", calculate initial accelerations according to
  //  prescribed Dirichlet values for generalized-alpha time
  //  integration and values at intermediate time steps
  // -------------------------------------------------------------------
  if (timealgo_==INPAR::FLUID::timeint_afgenalpha)
  {
    throw std::runtime_error( "not supported" );
  }
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void STK::FLD::Fluid::LinearSolve()
{
  dserror( "not implemented" );
}


namespace STK
{
  namespace FLD
  {

    /*----------------------------------------------------------------------*
     *----------------------------------------------------------------------*/
    template<DRT::Element::DiscretizationType distype>
    void EvaluateFluidElements( STK::Discretization & dis,
                                stk::mesh::Bucket & bucket,
                                STK::LinearSystem & system,
                                Teuchos::RCP<MAT::Material> mat,
                                Teuchos::ParameterList & eleparams,
                                stk::mesh::VectorField * velnp,
                                stk::mesh::ScalarField * pressure,
                                stk::mesh::VectorField * accnp,
                                stk::mesh::VectorField * hist,
                                bool alefluid )
    {
      typedef typename DRT::ELEMENTS::Fluid3Impl<distype> FluidElementType;

      const int nen = FluidElementType::nen_;
      const int nsd = FluidElementType::nsd_;
      const int dim = (nsd+1)*nen;

      Epetra_SerialDenseMatrix epetra_elemat1( dim, dim );
      Epetra_SerialDenseMatrix epetra_elemat2( dim, dim );
      Epetra_SerialDenseVector epetra_elevec1( dim );

      LINALG::Matrix<dim,dim> elemat1( epetra_elemat1, true );
      LINALG::Matrix<dim,dim> elemat2( epetra_elemat2, true );
      LINALG::Matrix<dim,  1> elevec1( epetra_elevec1, true );

      LINALG::Matrix<nsd,nen> edeadaf( true ); // body force

      LINALG::Matrix<nsd,nen> evelaf( true );
      LINALG::Matrix<nen,1>   epreaf( true );
      LINALG::Matrix<nen,1>   escaaf( true );
      LINALG::Matrix<nsd,nen> emhist( true );
      LINALG::Matrix<nsd,nen> eaccam( true );
      LINALG::Matrix<nen,1>   escadtam( true );
      LINALG::Matrix<nsd,nen> eveln( true );
      LINALG::Matrix<nen,1>   escaam( true );
      LINALG::Matrix<nsd,nen> edispnp( true );
      LINALG::Matrix<nsd,nen> egridv( true );
      LINALG::Matrix<nsd,nen> fsevelaf( true );

      std::vector<int> lm( dim );
      std::vector<int> lmowner( dim );

      FluidElementType * f3 = FluidElementType::Instance();

      for ( stk::mesh::Bucket::iterator j=bucket.begin();
            j!=bucket.end();
            ++j )
      {
        stk::mesh::Entity & e = *j;

        dis.LocationVector( e, lm, lmowner );

        elemat1.PutScalar( 0. );
        //elemat2.PutScalar( 0. );
        elevec1.PutScalar( 0. );

        // We need just a few element matrices for normal fluid

        ExtractValues<nsd,nen>( e, & dis.GetMesh().Coordinates(), f3->xyz() );

        ExtractValues<nsd,nen>( e, velnp, evelaf );
        ExtractValues<nen>( e, pressure, epreaf );

        ExtractValues<nsd,nen>( e, hist, emhist );

        ExtractValues<nsd,nen>( e, accnp, eaccam );
        //ExtractValues<nen>( e, _, escadtam );

        if (alefluid)
        {
          //ExtractValues<nsd,nen>( e, _, edispnp );
          //ExtractValues<nsd,nen>( e, _, egridv );
        }

        int eid = static_cast<int>( e.key().id() );

        int result = f3->Evaluate(
          eid,
          eleparams,
          edeadaf,
          elemat1,
          elemat2,
          elevec1,
          evelaf,
          epreaf,
          escaaf,
          emhist,
          eaccam,
          escadtam,
          eveln,
          escaam,
          edispnp,
          egridv,
          fsevelaf,
          mat,
          alefluid,
          e.owner_rank()==dis.GetMesh().parallel_rank(),
          0.0,
          NULL,
          NULL,
          NULL);

        // Assemble

        system.Assemble(eid,epetra_elemat1,lm,lmowner,lm);
        system.Assemble(eid,epetra_elevec1,lm);
      }
    }

  }
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void STK::FLD::Fluid::NonlinearSolve()
{
  STK::Mesh & mesh = dis_.GetMesh();
  stk::mesh::BulkData & bulk = mesh.BulkData();
  int p_rank = bulk.parallel_rank();

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

  Teuchos::RCP<Epetra_Vector> incvel = Teuchos::rcp( new Epetra_Vector( dis_.DofRowMap() ) );
  Teuchos::RCP<Epetra_Vector> tmpvec = Teuchos::rcp( new Epetra_Vector( dis_.DofRowMap() ) );

//   algebra::Update( bulk, *incvel_,  0.0, *velnp_,    1. );
//   algebra::Update( bulk, *incpres_, 0.0, *pressure_, 1. );

#ifdef WRITEDEBUG

  std::vector< const stk::mesh::FieldBase * > out_fields ;
  out_fields.push_back( velnp_ );
  out_fields.push_back( pressure_ );
  out_fields.push_back( incvel_ );
  out_fields.push_back( incpres_ );
  out_fields.push_back( resvel_ );
  out_fields.push_back( respres_ );

  std::stringstream s;
  s << "fluid_" << step_ << ".exo";

  Teuchos::RCP<phdmesh::exodus::FileOutput> exo = dis_.GetMesh().OutputContext( s.str(), "Navier-Stokes Problem", out_fields );

#endif

  LinearSystem & system = state_->system_;
  stk::mesh::Selector dirichlet_select = dis_.DirichletSelector();

  while (stopnonliniter==false)
  {
    itnum++;

    // get cpu time
    const double tcpu=Teuchos::Time::wallTime();

    system.Zero();

    // add Neumann loads
    //residual_->Update(1.0,*neumann_loads_,0.0);

    // element loop

    // set general element parameters
    Teuchos::ParameterList eleparams;
    eleparams.set("dt",dta_);
    eleparams.set("theta",theta_);
//     eleparams.set("omtheta",omtheta_);
    eleparams.set<int>("form of convective term",convform_);
//     eleparams.set("fs subgrid viscosity",fssgv_);
    eleparams.set<int>("Linearisation", newton_);
    eleparams.set<int>("Physical Type", physicaltype_);
//     eleparams.set("mixed_formulation", params_.get<bool>("mixed_formulation", false));

    // parameters for stabilization
    eleparams.sublist("STABILIZATION") = fdyn.sublist("STABILIZATION");

    // parameters for stabilization
    eleparams.sublist("TURBULENCE MODEL") = fdyn.sublist("TURBULENCE MODEL");

    eleparams.set("thermpress at n+alpha_F/n+1",0.0);
    eleparams.set("thermpress at n+alpha_M/n",0.0);
    eleparams.set("thermpressderiv at n+alpha_M/n+1",0.0);

    // set scheme-specific element parameters and vector values
    if (timealgo_==INPAR::FLUID::timeint_stationary)
    {
      eleparams.set("action","calc_fluid_stationary_systemmat_and_residual");
      eleparams.set("using generalized-alpha time integration",false);
      eleparams.set("total time",time_);
      eleparams.set("is stationary", true);
    }
    else if (timealgo_==INPAR::FLUID::timeint_afgenalpha)
    {
//       eleparams.set("action","calc_fluid_afgenalpha_systemmat_and_residual");
//       eleparams.set("using generalized-alpha time integration",true);
//       eleparams.set("total time",time_-(1-alphaF_)*dta_);
//       eleparams.set("is stationary", false);
//       eleparams.set("alphaF",alphaF_);
//       eleparams.set("alphaM",alphaM_);
//       eleparams.set("gamma",gamma_);
    }
    else
    {
      eleparams.set("action","calc_fluid_systemmat_and_residual");
      eleparams.set("using generalized-alpha time integration",false);
      eleparams.set("total time",time_);
      eleparams.set("is stationary", false);
    }

    // no need to call preevaluate here
    DRT::ELEMENTS::Fluid3ImplParameter* f3Parameter = DRT::ELEMENTS::Fluid3ImplParameter::Instance();
    f3Parameter->SetParameter( eleparams );

    // We need to evaluate all active elements including ghosted ones. There
    // might even be elements that do not have any node owned by the local
    // processor...

    //stk::mesh::Selector activeselector = mesh.ActivePart();
    stk::mesh::Selector ownedselector = mesh.OwnedPart() & mesh.ActivePart();

    const std::vector<stk::mesh::Bucket*> & elements = bulk.buckets( stk::mesh::Element );
    for ( std::vector<stk::mesh::Bucket*>::const_iterator j=elements.begin();
          j!=elements.end();
          ++j )
    {
      stk::mesh::Bucket & bucket = **j;

      if ( ownedselector( bucket ) )
      {
        // create material on-demand
        MAT::PAR::Parameter * matpat = dis_.MaterialParameter( bucket );
        if ( matpat==NULL )
          dserror( "no material parameters" );
        Teuchos::RCP<MAT::Material> mat = matpat->CreateMaterial();

        if ( stk::mesh::has_superset( bucket, mesh.Quad4() ) )
        {
          STK::FLD::EvaluateFluidElements<DRT::Element::quad4>( dis_, bucket,
                                                                system,
                                                                mat, eleparams,
                                                                velnp_, pressure_, accnp_, hist_,
                                                                alefluid_ );
        }
        else
        {
          dserror( "unknown bucket type" );
        }
      }
    }

    // finalize the complete matrix
    system.Complete();

    // end time measurement for element
    dtele = Teuchos::Time::wallTime()-tcpu;

    ScatterFieldData( *system.RHS(), resvel_, respres_ );

    // zero dirichlet values
    algebra::PutScalar( bulk, dirichlet_select, *resvel_, 0.0 );
    algebra::PutScalar( bulk, dirichlet_select, *respres_, 0.0 );

    double incvelnorm_L2 = algebra::Norm2( bulk, *incvel_ );
    double incprenorm_L2 = algebra::Norm2( bulk, *incpres_ );

    double velnorm_L2 = algebra::Norm2( bulk, *velnp_ );
    double prenorm_L2 = algebra::Norm2( bulk, *pressure_ );

    double vresnorm = algebra::Norm2( bulk, *resvel_ );
    double presnorm = algebra::Norm2( bulk, *respres_ );

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
      GatherFieldData( velnp_, pressure_, *incvel );
    }
    else
    {
      incvel->PutScalar( 0.0 );
    }

    system.ApplyDirichlet( *incvel );

    const double tcpusolve=Teuchos::Time::wallTime();

    solver_->Solve( system.Matrix(), incvel, system.RHS(), true, itnum==1 );
    //solver_->ResetTolerance();

    ScatterFieldData( *incvel, incvel_, incpres_ );

#if 0
    {
      state_->sysmat_->Multiply( false, *incvel, *tmpvec );
      tmpvec->Update( -1, *state_->residual_, 1 );
      double norm;
      tmpvec->Norm2( &norm );
      std::cout << "solve test: " << norm << "\n";
    }
#endif

    // end time measurement for solver
    dtsolve = Teuchos::Time::wallTime()-tcpusolve;

    // -------------------------------------------------------------------
    // update velocity and pressure values by increments
    // -------------------------------------------------------------------
    if ( itnum==1 )
    {
      //incvel->Update( -1.0, *velnp, 1.0 );
      algebra::Update( bulk, *incvel_,  1., *velnp_,    -1. );
      algebra::Update( bulk, *incpres_, 1., *pressure_, -1. );

      //velnp->Update( 1.0, *incvel, 1.0 );
      algebra::Update( bulk, *velnp_,    1., *incvel_,  1. );
      algebra::Update( bulk, *pressure_, 1., *incpres_, 1. );
    }
    else
    {
      //velnp->Update( 1.0, *incvel, 1.0 );
      algebra::Update( bulk, *velnp_,    1., *incvel_,  1. );
      algebra::Update( bulk, *pressure_, 1., *incpres_, 1. );
    }

#ifdef WRITEDEBUG
    exo->write( itnum, itnum );
#endif
  }
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void STK::FLD::Fluid::TimeUpdate()
{
  stk::mesh::BulkData & bulk = dis_.GetMesh().BulkData();

  if ( step_ == 1 )
  {
    switch ( timealgo_ )
    {
      case INPAR::FLUID::timeint_stationary: /* no accelerations for stationary problems*/
      {
        algebra::PutScalar( bulk, *accnp_, 0.0 );
        break;
      }
      case INPAR::FLUID::timeint_one_step_theta: /* One step Theta time integration */
      case INPAR::FLUID::timeint_bdf2:    /* 2nd order backward differencing BDF2 */
      {
        // do just a linear interpolation within the first timestep
        algebra::Update( bulk, *accnp_, 0.0, *velnp_, 1.0/dta_, *veln_, -1.0/dta_ );
        break;
      }
      case INPAR::FLUID::timeint_afgenalpha: /* Af-generalized-alpha time integration */
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
      case INPAR::FLUID::timeint_stationary: /* no accelerations for stationary problems*/
      {
        algebra::PutScalar( bulk, *accnp_, 0.0 );
        break;
      }
      case INPAR::FLUID::timeint_one_step_theta: /* One-step-theta time integration */
      {
        const double fact1 = 1.0/(theta_*dta_);
        const double fact2 =-1.0/theta_ +1.0;   /* = -1/Theta + 1 */

        algebra::Update( bulk, *accnp_, 0.0, *velnp_, fact1, *veln_, -fact1, *accn_, fact2 );
        break;
      }
      case INPAR::FLUID::timeint_bdf2:    /* 2nd order backward differencing BDF2 */
      {
        if (dta_*dtp_ < 1e-15) throw std::runtime_error("zero time step size");
        const double sum = dta_ + dtp_;

        algebra::Update( bulk,
                         *accnp_, 0.0,
                         *velnp_, (2.0*dta_+dtp_)/(dta_*sum),
                         *veln_, -sum/(dta_*dtp_),
                         *velnm_, dta_/(dtp_*sum) );
        break;
      }
      case INPAR::FLUID::timeint_afgenalpha: /* Af-generalized-alpha time integration */
      {
        // do nothing: new acceleration is calculated at beginning of next time step
        break;
      }
      default:
        throw std::runtime_error("Time integration scheme unknown!");
    }
  }

  // update old acceleration
  algebra::Update( bulk, *accn_, 0.0, *accnp_, 1.0 );

  // velocities/pressures of this step become most recent
  // velocities/pressures of the last step
  algebra::Update( bulk, *velnm_, 0.0, *veln_ , 1.0 );
  algebra::Update( bulk, *veln_ , 0.0, *velnp_, 1.0 );
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void STK::FLD::Fluid::StatisticsAndOutput()
{
  stk::mesh::BulkData & bulk = dis_.GetMesh().BulkData();

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

  state_->Output( *this );
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void STK::FLD::Fluid::AdaptiveNonlinearSolve()
{
  NonlinearSolve();

  if ( step_ % refinestep_ == 0 )
  {
    std::vector<stk::mesh::EntityKey> refine;
    std::vector<stk::mesh::EntityKey> unrefine;

    ErrorEstimate( refine, unrefine );

    unsigned refine_size   = refine.size();
    unsigned unrefine_size = unrefine.size();

    int r = 0;
    for ( ; r < maxrefine_; ++r )
    {
      stk::all_reduce( Discretization().GetMesh().parallel(), stk::ReduceSum<1>( &refine_size ) );
      //stk::all_reduce( Discretization().GetMesh().parallel(), stk::ReduceSum<1>( &unrefine_size ) );

      if ( refine_size==0 )
        break;

      dis_.AdaptMesh( refine, unrefine );
      dis_.GetMesh().Statistics();

      NonlinearSolve();

      refine.clear();
      unrefine.clear();
      ErrorEstimate( refine, unrefine );

      refine_size   = refine.size();
      unrefine_size = unrefine.size();
    }

    stk::all_reduce( Discretization().GetMesh().parallel(), stk::ReduceSum<1>( &refine_size ) );
    stk::all_reduce( Discretization().GetMesh().parallel(), stk::ReduceSum<1>( &unrefine_size ) );

    if ( refine_size!=0 or unrefine_size!=0 )
    {
      if ( r==maxrefine_ )
      {
        stk::mesh::BulkData & bulk = Discretization().GetMesh().BulkData();
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
void STK::FLD::Fluid::ErrorEstimate( std::vector<stk::mesh::EntityKey> & refine,
                                     std::vector<stk::mesh::EntityKey> & unrefine )
{
  stk::mesh::BulkData & bulk = Discretization().GetMesh().BulkData();
  stk::mesh::Part & active = Discretization().GetMesh().ActivePart();
  stk::mesh::Part & owned  = Discretization().GetMesh().OwnedPart();

  //PrintMesh( std::cout, bulk, active );

  algebra::PutScalar( bulk, *error_ , 0.0 );
  algebra::PutScalar( bulk, *volume_, 0.0 );

  //Discretization().GetMesh().Dump( "mesh-debug" );

  STK::FLD::FluidJumpIntegrator integrator( Discretization() );
  integrator.Integrate( Discretization().GetMesh().Coordinates(),
                        *velnp_, *pressure_,
                        *error_, *volume_,
                        *resvel_, *respres_ );

  double min_err = std::numeric_limits<double>::max();
  double max_err = std::numeric_limits<double>::min();

  for ( STK::ElementIterator i( bulk, stk::mesh::Selector( owned & active ) );
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

  for ( STK::ElementIterator i( bulk, stk::mesh::Selector( active ) );
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

#if 0
  if ( bulk.parallel_rank()==0 )
  {
    std::cout << "    #refine elements: " << refine.size()
              << "    #unrefine elements: " << unrefine.size()
              << "    #passed elements: " << passed
              << "\n";
  }
#endif
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void STK::FLD::Fluid::ScatterFieldData( const Epetra_Vector & v,
                                        stk::mesh::VectorField * vel,
                                        stk::mesh::ScalarField * pres )
{
  std::vector<stk::mesh::FieldBase*> fields( 2 );
  fields[0] = vel;
  fields[1] = pres;
  dis_.ScatterFieldData( v, fields );
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void STK::FLD::Fluid::GatherFieldData( stk::mesh::VectorField * vel,
                                       stk::mesh::ScalarField * pres,
                                       Epetra_Vector & v )
{
  std::vector<stk::mesh::FieldBase*> fields( 2 );
  fields[0] = vel;
  fields[1] = pres;
  dis_.GatherFieldData( fields, v );
}


#endif
