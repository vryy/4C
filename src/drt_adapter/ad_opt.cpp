/*!------------------------------------------------------------------------------------------------*
\file ad_opt.cpp

\brief topology optimization adapter

<pre>
Maintainer: Martin Winklmaier
            winklmaier@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15241
</pre>
 *------------------------------------------------------------------------------------------------*/


#include "ad_opt.H"
#include "../drt_io/io.H"
#include "../drt_io/io_control.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_opti/topopt_optimizer.H"
#include "../drt_opti/topopt_utils.H"
#include "../drt_inpar/inpar_parameterlist_utils.H"


/// constructor
ADAPTER::TopOptBaseAlgorithm::TopOptBaseAlgorithm(
    const Teuchos::ParameterList& prbdyn,  ///< problem-dependent parameters
    const std::string disname  ///< optimization field discretization name(default: "opti")
)
{
  DRT::Problem* problem = DRT::Problem::Instance();

  // -------------------------------------------------------------------
  // access the fluid and the optimization discretization
  // -------------------------------------------------------------------
  Teuchos::RCP<DRT::Discretization> optidis = Teuchos::null;
  optidis = problem->GetDis(disname);
  Teuchos::RCP<DRT::Discretization> fluiddis = Teuchos::null;
  fluiddis = problem->GetDis("fluid");

  // -------------------------------------------------------------------
  // check degrees of freedom in the discretization
  // -------------------------------------------------------------------
  if (!optidis->Filled()) dserror("optimization discretization should be filled before");
  if (!fluiddis->Filled()) dserror("fluid discretization should be filled before");

  // -------------------------------------------------------------------
  // context for output and restart
  // -------------------------------------------------------------------
  // output control for optimization field
  // equal to output for fluid equations except for the filename
  // and the - not necessary - input file name
  Teuchos::RCP<IO::OutputControl> optioutput =
      Teuchos::rcp(new IO::OutputControl(optidis->Comm(), problem->ProblemName(),
          problem->SpatialApproximation(), problem->OutputControlFile()->InputFileName(),
          TOPOPT::modifyFilename(problem->OutputControlFile()->FileName(), "xxx_opti_",
              (bool)DRT::Problem::Instance()->Restart(), true),
          problem->NDim(), problem->Restart(), problem->OutputControlFile()->FileSteps(),
          DRT::INPUT::IntegralValue<int>(problem->IOParams(), "OUTPUT_BIN")));

  Teuchos::RCP<IO::DiscretizationWriter> output = optidis->Writer();
  output->SetOutput(optioutput);
  output->WriteMesh(0, 0.0);

  // -------------------------------------------------------------------
  // create instance of the optimization class (call the constructor)
  // -------------------------------------------------------------------
  optimizer_ = Teuchos::rcp(new TOPOPT::Optimizer(optidis, fluiddis, prbdyn, output));

  return;
}
