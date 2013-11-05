/*----------------------------------------------------------------------*/
/*!
 \file fluid_ele_poro_evaluate.cpp

 \brief

 <pre>
   Maintainer: Anh-Tu Vuong
               vuong@lnm.mw.tum.de
               http://www.lnm.mw.tum.de
               089 - 289-15264
 </pre>
 *----------------------------------------------------------------------*/


#include "fluid_ele_poro.H"
#include "fluid_ele_action.H"
#include "fluid_ele_parameter.H"

#include "fluid_ele_factory.H"
#include "fluid_ele_interface.H"

#include "../drt_inpar/inpar_fpsi.H"

#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_discret.H"


/*---------------------------------------------------------------------*
|  Call the element to set all basic parameter                         |
*----------------------------------------------------------------------*/
void DRT::ELEMENTS::FluidPoroEleType::PreEvaluate(DRT::Discretization&            dis,
                                            Teuchos::ParameterList&               p,
                                            Teuchos::RCP<LINALG::SparseOperator>  systemmatrix1,
                                            Teuchos::RCP<LINALG::SparseOperator>  systemmatrix2,
                                            Teuchos::RCP<Epetra_Vector>           systemvector1,
                                            Teuchos::RCP<Epetra_Vector>           systemvector2,
                                            Teuchos::RCP<Epetra_Vector>           systemvector3)
{
  const FLD::Action action = DRT::INPUT::get<FLD::Action>(p,"action");

  int num=0;
  if(DRT::Problem::Instance()->ProblemType()==prb_fpsi)
  {
    if(dis.Name()=="fluid")
    {
      num=INPAR::FPSI::fluid;
    }
    else if(dis.Name()=="porofluid")
    {
      num=INPAR::FPSI::porofluid;
    }
   }
  p.set<int>("numfield",num);

  if (action == FLD::set_general_fluid_parameter)
  {//this should be unnecessary when paramter lists are separated
    Teuchos::RCP<DRT::ELEMENTS::FluidEleParameter> fldpara = DRT::ELEMENTS::FluidEleParameter::Instance(num);
    fldpara->SetElementGeneralFluidParameter(p,dis.Comm().MyPID());
  }
  else if (action == FLD::set_time_parameter)
  {//this should be unnecessary when paramter lists are separated
    Teuchos::RCP<DRT::ELEMENTS::FluidEleParameter> fldpara = DRT::ELEMENTS::FluidEleParameter::Instance(num);
    fldpara->SetElementTimeParameter(p);
  }
  //poro specific actions
  else if (action == FLD::set_poro_parameter)
  {
    Teuchos::RCP<DRT::ELEMENTS::FluidEleParameter> fldpara = DRT::ELEMENTS::FluidEleParameter::Instance(num);
    fldpara->SetElementPoroParameter(p);
  }
  else
    //call standard fluid type
    FluidType::PreEvaluate(dis, p, systemmatrix1, systemmatrix2, systemvector1, systemvector2, systemvector3);

  return;
}

/*----------------------------------------------------------------------*
 |  evaluate the element (public)                            g.bau 03/07|
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::FluidPoro::Evaluate( Teuchos::ParameterList&   params,
                                        DRT::Discretization&      discretization,
                                        std::vector<int>&         lm,
                                        Epetra_SerialDenseMatrix& elemat1,
                                        Epetra_SerialDenseMatrix& elemat2,
                                        Epetra_SerialDenseVector& elevec1,
                                        Epetra_SerialDenseVector& elevec2,
                                        Epetra_SerialDenseVector& elevec3)
{
  // get the action required
  const FLD::Action act = DRT::INPUT::get<FLD::Action>(params,"action");

  // get material
  RCP<MAT::Material> mat = Material();

  // switch between different physical types as used below
  std::string impltype = "std";
  switch(params.get<int>("physical type",INPAR::FLUID::incompressible))
  {
  //incompressible case included for setting general fluid parameters
  case INPAR::FLUID::incompressible:    impltype = "std";                 break;
  case INPAR::FLUID::poro:              impltype = "poro";                break;
  case INPAR::FLUID::poro_p1:           impltype = "poro_p1";             break;
  case INPAR::FLUID::poro_p2:           impltype = "poro_p2";             break;
  default: dserror("invalid physical type '%s' for porous fluid!",impltype.c_str() );  break;
  }

  switch(act)
  {
    //-----------------------------------------------------------------------
    // standard implementation enabling time-integration schemes such as
    // one-step-theta, BDF2, and generalized-alpha (n+alpha_F and n+1)
    //-----------------------------------------------------------------------
    case FLD::calc_fluid_systemmat_and_residual:
    {
      return DRT::ELEMENTS::FluidFactory::ProvideImpl(Shape(), impltype)->Evaluate(
              this,
              discretization,
              lm,
              params,
              mat,
              elemat1,
              elemat2,
              elevec1,
              elevec2,
              elevec3 );
    }
    break;
    //-----------------------------------------------------------------------
    // standard implementation enabling time-integration schemes such as
    // one-step-theta, BDF2, and generalized-alpha (n+alpha_F and n+1)
    // for the particular case of porous flow
    //-----------------------------------------------------------------------
    /***********************************************/
    case FLD::calc_porousflow_fluid_coupling:
    {
      return DRT::ELEMENTS::FluidFactory::ProvideImpl(Shape(), impltype)->Evaluate(
            this,
            discretization,
            lm,
            params,
            mat,
            elemat1,
            elemat2,
            elevec1,
            elevec2,
            elevec3,
            true);
    }
    break;
    //-----------------------------------------------------------------------
    // standard implementation enabling time-integration schemes such as
    // one-step-theta, BDF2, and generalized-alpha (n+alpha_F and n+1)
    // for evaluation of off-diagonal matrix block for monolithic
    // porous flow with scalar transport
    //-----------------------------------------------------------------------
    case FLD::calc_poroscatra_mono_odblock:
    {
      switch(params.get<int>("physical type",INPAR::FLUID::incompressible))
      {
      case INPAR::FLUID::poro_p2:
      {
      return DRT::ELEMENTS::FluidFactory::ProvideImpl(Shape(), "poro_p2")->Evaluate(
          this,
          discretization,
          lm,
          params,
          mat,
          elemat1,
          elemat2,
          elevec1,
          elevec2,
          elevec3,
          true);
      }
      break;
      default:
        dserror("Unknown physical type for poroelasticity\n");
      break;
      }
    }
    break;
    case FLD::calc_volume:
    {
      return DRT::ELEMENTS::FluidFactory::ProvideImpl(Shape(), impltype)->EvaluateService(this,
                                                                       params,
                                                                       mat,
                                                                       discretization,
                                                                       lm,
                                                                       elemat1,
                                                                       elemat2,
                                                                       elevec1,
                                                                       elevec2,
                                                                       elevec3);
      break;
    }
    case FLD::set_poro_parameter:
      break;
    default:
      //call evaluate of standard fluid
      return Fluid::Evaluate(params,
                             discretization,
                             lm,
                             elemat1,
                             elemat2,
                             elevec1,
                             elevec2,
                             elevec3);
      break;
  } // end of switch(act)

  return 0;
} // end of DRT::ELEMENTS::Fluid::Evaluate
