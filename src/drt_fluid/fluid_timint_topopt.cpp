/*----------------------------------------------------------------------*/
/*!
\file fluid_timint_topopt.cpp

\brief TimIntTopOpt

\level 2
<pre>
\maintainer Benjamin Krank & Martin Kronbichler
             {krank,kronbichler}@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289-15252/-235
</pre>
*/
/*----------------------------------------------------------------------*/

#include "fluid_timint_topopt.H"
#include "../drt_opti/topopt_optimizer.H"
#include "../drt_fluid_ele/fluid_ele_action.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_mat/matpar_bundle.H"
#include "../drt_mat/optimization_density.H"
#include "../drt_io/io.H"
#include "../linalg/linalg_utils.H"


/*----------------------------------------------------------------------*
 |  Constructor (public)                                       bk 11/13 |
 *----------------------------------------------------------------------*/
FLD::TimIntTopOpt::TimIntTopOpt(
        const Teuchos::RCP<DRT::Discretization>&      actdis,
        const Teuchos::RCP<LINALG::Solver>&           solver,
        const Teuchos::RCP<Teuchos::ParameterList>&   params,
        const Teuchos::RCP<IO::DiscretizationWriter>& output,
        bool                                          alefluid /*= false*/)
    : FluidImplicitTimeInt(actdis,solver,params,output,alefluid),
      optimizer_(Teuchos::null)
{
  return;
}


/*----------------------------------------------------------------------*
 |  initialize algorithm                                rasthofer 04/14 |
 *----------------------------------------------------------------------*/
void FLD::TimIntTopOpt::Init()
{
  //set some Topopt-specific parameters
  SetElementCustomParameter();
  return;
}


/*----------------------------------------------------------------------*
| Destructor dtor (public)                                    bk 11/13 |
*----------------------------------------------------------------------*/
FLD::TimIntTopOpt::~TimIntTopOpt()
{
  return;
}

/*----------------------------------------------------------------------*
| Import fluid data to optimizer at end of Solve()             bk 12/13 |
*----------------------------------------------------------------------*/
void FLD::TimIntTopOpt::Solve()
{

  FluidImplicitTimeInt::Solve();

  optimizer_->ImportFluidData(velnp_,step_);

  // initial solution (=u_0) is old solution at time step 1
  if (step_==1 and timealgo_!=INPAR::FLUID::timeint_stationary)
    optimizer_->ImportFluidData(veln_,0);

  return;
}

/*----------------------------------------------------------------------*
| TopOptSetDensity                                             bk 12/13 |
*----------------------------------------------------------------------*/
void FLD::TimIntTopOpt::SetCustomEleParamsAssembleMatAndRHS(Teuchos::ParameterList& eleparams)
{
  eleparams.set("topopt_density",Teuchos::rcp_const_cast<const Epetra_Vector>(density_scaling_));
  eleparams.set("dens_type",DRT::INPUT::IntegralValue<INPAR::TOPOPT::DensityField>(optimizer_->OptiParams(),"DENS_TYPE"));
  return;
}

/*----------------------------------------------------------------------*
 | sent density field for topology optimization         winklmaier 12/11|
 *----------------------------------------------------------------------*/
void FLD::TimIntTopOpt::SetTopOptData(
    Teuchos::RCP<Epetra_Vector> density,
    Teuchos::RCP<TOPOPT::Optimizer>& optimizer
)
{
  // currently the maps have to fit and, thus, this works
  // in the future this simple procedure may have to be altered,
  // see setiterlomafields or settimelomafields for examples
  density_scaling_ = density;
  optimizer_=optimizer;
}

// -------------------------------------------------------------------
// set topopt parameters                           winklmaier  07/2013
// -------------------------------------------------------------------
void FLD::TimIntTopOpt::SetElementCustomParameter()
{
  Teuchos::ParameterList eleparams;

  eleparams.set<int>("action",FLD::set_topopt_parameter);

  // get the material
  const int nummat = DRT::Problem::Instance()->Materials()->Num();
  for (int id = 1; id-1 < nummat; ++id)
  {
    Teuchos::RCP<const MAT::PAR::Material> imat = DRT::Problem::Instance()->Materials()->ById(id);

    if (imat == Teuchos::null)
      dserror("Could not find material Id %d", id);
    else
    {
      switch (imat->Type())
      {
      case INPAR::MAT::m_fluid:
        break;
      case INPAR::MAT::m_opti_dens:
      {
        const MAT::PAR::Parameter* matparam = imat->Parameter();
        const MAT::PAR::TopOptDens* mat = static_cast<const MAT::PAR::TopOptDens* >(matparam);

        eleparams.set("MIN_PORO",mat->PoroBdDown());
        eleparams.set("MAX_PORO",mat->PoroBdUp());


        const INPAR::TOPOPT::OptiCase testcase = (INPAR::TOPOPT::OptiCase)(params_->get<int>("opti testcase"));
        switch (testcase)
        {
        case INPAR::TOPOPT::optitest_channel:
        case INPAR::TOPOPT::optitest_channel_with_step:
        case INPAR::TOPOPT::optitest_cornerflow:
        case INPAR::TOPOPT::optitest_lin_poro:
        case INPAR::TOPOPT::optitest_quad_poro:
        case INPAR::TOPOPT::optitest_cub_poro:
        {
          eleparams.set("SMEAR_FAC",(double)(-(int)testcase));
          break;
        }
        default:
        {
          eleparams.set("SMEAR_FAC",mat->SmearFac());
          break;
        }
        }

        break;
      }
      default:
      {
        dserror("unknown material %s",imat->Name().c_str());
        break;
      }
      }
    }
  }

  // call standard loop over elements
  discret_->Evaluate(eleparams,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null);
  return;
}


