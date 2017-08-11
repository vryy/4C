/*-----------------------------------------------------------*/
/*!
\file str_timint_basedataio.cpp

\brief Input/output data container for the structural (time)
       integration

\maintainer Michael Hiermeier

\date Jan 11, 2016

\level 3

*/
/*-----------------------------------------------------------*/


#include "str_timint_basedataio.H"
#include "str_timint_basedataio_runtime_vtk_output.H"
#include "str_timint_basedataio_runtime_vtp_output.H"

#include "../drt_io/every_iteration_writer.H"
#include "../solver_nonlin_nox/nox_nln_aux.H"

#include <NOX_Solver_Generic.H>

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
STR::TIMINT::BaseDataIO::BaseDataIO()
    : isinit_(false),
      issetup_(false),
      output_(Teuchos::null),
      writer_every_iter_(Teuchos::null),
      params_runtime_vtk_output_(Teuchos::null),
      params_runtime_vtp_output_(Teuchos::null),
      energyfile_(Teuchos::null),
      errfile_(NULL),
      gmsh_out_(false),
      printlogo_(false),
      printerrfile_(false),
      printiter_(false),
      outputeveryiter_(false),
      writesurfactant_(false),
      writestate_(false),
      writevelacc_(false),
      writejac2matlab_(false),
      firstoutputofrun_(false),
      printscreen_(-1),
      outputcounter_(-1),
      writerestartevery_(-1),
      writereducedrestart_(-1),
      writeresultsevery_(-1),
      writeenergyevery_(-1),
      kinergy_(-1.0),
      intergy_(-1.0),
      extergy_(-1.0),
      writestress_(INPAR::STR::stress_none),
      writecouplstress_(INPAR::STR::stress_none),
      writestrain_(INPAR::STR::strain_none),
      writeplstrain_(INPAR::STR::strain_none)
{
  // empty constructor
}


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::TIMINT::BaseDataIO::Init(const Teuchos::ParameterList& ioparams,
    const Teuchos::ParameterList& sdynparams,
    const Teuchos::ParameterList& xparams,
    Teuchos::RCP<IO::DiscretizationWriter> output
    )
{
  // We have to call Setup() after Init()
  issetup_ = false;

  // ---------------------------------------------------------------------------
  // initialize the printing and output parameters
  // ---------------------------------------------------------------------------
  {
    output_ = output;
    printscreen_ = ioparams.get<int>("STDOUTEVRY");
    printlogo_ = (printscreen_>0 ? true : false);
    errfile_ = xparams.get<FILE*>("err file");
    gmsh_out_ = (bool) DRT::INPUT::IntegralValue<int>(ioparams,"OUTPUT_GMSH");
    printerrfile_ = (true and errfile_);
    printiter_ = true;
    p_io_every_iteration_ = Teuchos::rcp( new Teuchos::ParameterList(
        ioparams.sublist( "EVERY ITERATION" ) ) );
    outputeveryiter_ = DRT::INPUT::IntegralValue<bool>( *p_io_every_iteration_,
        "OUTPUT_EVERY_ITER" );
    writerestartevery_ = sdynparams.get<int>("RESTARTEVRY");
    writereducedrestart_ = xparams.get<int>("REDUCED_OUTPUT");
    writestate_ = (bool) DRT::INPUT::IntegralValue<int>(ioparams,"STRUCT_DISP");
    writevelacc_ = (bool) DRT::INPUT::IntegralValue<int>(ioparams,"STRUCT_VEL_ACC");
    writejac2matlab_ = (bool) DRT::INPUT::IntegralValue<int>(ioparams,"STRUCT_JACOBIAN_MATLAB");
    firstoutputofrun_ = true;
    writeresultsevery_ = sdynparams.get<int>("RESULTSEVRY");
    writestress_ = DRT::INPUT::IntegralValue<INPAR::STR::StressType>(ioparams,"STRUCT_STRESS");
    writecouplstress_ = DRT::INPUT::IntegralValue<INPAR::STR::StressType>(ioparams,"STRUCT_COUPLING_STRESS");
    writestrain_ = DRT::INPUT::IntegralValue<INPAR::STR::StrainType>(ioparams,"STRUCT_STRAIN");
    writeplstrain_ = DRT::INPUT::IntegralValue<INPAR::STR::StrainType>(ioparams,"STRUCT_PLASTIC_STRAIN");
    writeenergyevery_ = sdynparams.get<int>("RESEVRYERGY");
    writesurfactant_ = (bool) DRT::INPUT::IntegralValue<int>(ioparams,"STRUCT_SURFACTANT");

    // check whether VTK output at runtime is desired
    if ( ioparams.sublist("RUNTIME VTK OUTPUT").get<int>("INTERVAL_STEPS") != -1 )
    {
      params_runtime_vtk_output_ = Teuchos::rcp( new ParamsRuntimeVtkOutput() );

      params_runtime_vtk_output_->Init( ioparams.sublist("RUNTIME VTK OUTPUT") );
      params_runtime_vtk_output_->Setup();
    }

    // check whether VTP output at runtime is desired
    if ( ioparams.sublist("RUNTIME VTP OUTPUT STRUCTURE").get<int>("INTERVAL_STEPS") != -1 )
    {
      params_runtime_vtp_output_ = Teuchos::rcp( new ParamsRuntimeVtpOutput() );

      params_runtime_vtp_output_->Init( ioparams.sublist("RUNTIME VTP OUTPUT STRUCTURE") );
      params_runtime_vtp_output_->Setup();
    }
  }

  isinit_ = true;

  return;
}


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::TIMINT::BaseDataIO::Setup()
{
  // safety check
  if (!IsInit())
    dserror("Init() has not been called, yet!");

  if ( outputeveryiter_ )
    writer_every_iter_ = Teuchos::rcp( new IO::EveryIterationWriter() );

  issetup_ = true;

  // Good bye
  return;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::TIMINT::BaseDataIO::InitSetupEveryIterationWriter(
    IO::EveryIterationWriterInterface* interface,
    Teuchos::ParameterList& p_nox )
{
  if ( not outputeveryiter_ )
    return;

  writer_every_iter_->Init( output_.get(), interface, *p_io_every_iteration_ );
  writer_every_iter_->Setup();

  Teuchos::ParameterList& p_sol_opt = p_nox.sublist("Solver Options");

  Teuchos::RCP<NOX::Abstract::PrePostOperator> prepost_solver_ptr =
      Teuchos::rcp( new NOX::NLN::Solver::PrePostOp::TIMINT::WriteOutputEveryIteration(
          *writer_every_iter_ ) );

  NOX::NLN::AUX::AddToPrePostOpVector( p_sol_opt, prepost_solver_ptr );
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
NOX::NLN::Solver::PrePostOp::TIMINT::WriteOutputEveryIteration::WriteOutputEveryIteration(
    IO::EveryIterationWriter& every_iter_writer )
    : every_iter_writer_( every_iter_writer )
{
  /* empty */
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void NOX::NLN::Solver::PrePostOp::TIMINT::WriteOutputEveryIteration::runPreSolve(
    const NOX::Solver::Generic& solver )
{
  every_iter_writer_.InitNewtonIteration();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void NOX::NLN::Solver::PrePostOp::TIMINT::WriteOutputEveryIteration::runPostIterate(
    const NOX::Solver::Generic& solver )
{
  const int newton_iteration = solver.getNumIterations();
  every_iter_writer_.AddNewtonIteration( newton_iteration );
}
