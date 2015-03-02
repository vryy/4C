/*!---------------------------------------------------------------------------

\file meshfree_fluid_cell_evaluate.cpp

\brief

<pre>
Maintainer: Keijo Nissen
            nissen@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15253
</pre>
*---------------------------------------------------------------------------*/

#include "meshfree_fluid_cell.H"
#include "meshfree_fluid_cell_interface.H"

#include "../drt_fluid_ele/fluid_ele_factory.H"
#include "../drt_fluid_ele/fluid_ele_parameter.H"
#include "../drt_fluid_ele/fluid_ele_parameter_std.H"
#include "../drt_fluid_ele/fluid_ele_parameter_timint.H"
#include "../drt_fluid_ele/fluid_ele_action.H"

#include "../drt_lib/drt_discret.H"

/*
  Depending on the type of action and the element type (tet, hex etc.),
  the elements allocate common static arrays.

  */

/*--------------------------------------------------------------------------*
 |  Call the element to set all basic parameter          (public) nis Jan13 |
 *--------------------------------------------------------------------------*/
void DRT::ELEMENTS::MeshfreeFluidType::PreEvaluate(DRT::Discretization&           dis,
                                            Teuchos::ParameterList&               p,
                                            Teuchos::RCP<LINALG::SparseOperator>  systemmatrix1,
                                            Teuchos::RCP<LINALG::SparseOperator>  systemmatrix2,
                                            Teuchos::RCP<Epetra_Vector>           systemvector1,
                                            Teuchos::RCP<Epetra_Vector>           systemvector2,
                                            Teuchos::RCP<Epetra_Vector>           systemvector3)
{
  const FLD::Action action = DRT::INPUT::get<FLD::Action>(p,"action");

  if (action==FLD::set_general_fluid_parameter)
  {
    DRT::ELEMENTS::FluidEleParameterStd* fldpara = DRT::ELEMENTS::FluidEleParameterStd::Instance();
    fldpara->SetElementGeneralFluidParameter(p,dis.Comm().MyPID());
  }
  else if (action==FLD::set_time_parameter)
  {
    DRT::ELEMENTS::FluidEleParameterTimInt* fldpara = DRT::ELEMENTS::FluidEleParameterTimInt::Instance();
    fldpara->SetElementTimeParameter(p);
  }

  return;
}


 /*----------------------------------------------------------------------*
 |  evaluate the element (public)                            g.bau 03/07|
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::MeshfreeFluid::Evaluate(Teuchos::ParameterList&            params,
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
  Teuchos::RCP<MAT::Material> mat = Material();

  switch(act)
  {
    //-----------------------------------------------------------------------
    // standard implementation enabling time-integration schemes such as
    // one-step-theta, BDF2, and generalized-alpha (n+alpha_F and n+1)
    //-----------------------------------------------------------------------
    case FLD::calc_fluid_systemmat_and_residual:
    {
      switch(params.get<int>("Physical Type",INPAR::FLUID::incompressible))
      {
      default:
        return DRT::ELEMENTS::FluidFactory::ProvideImplMeshfree(Shape(), "std_meshfree")->Evaluate(
            this,
            discretization,
            lm,
            params,
            mat,
            elemat1,
            elemat2,
            elevec1,
            elevec2,
            elevec3);
      }
      break;
    }
    case FLD::integrate_shape:
    {
      return DRT::ELEMENTS::FluidFactory::ProvideImplMeshfree(Shape(), "std_meshfree")->EvaluateService(this,
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
    case FLD::calc_fluid_error:
    case FLD::set_general_fluid_parameter:
    case FLD::set_time_parameter:
    case FLD::set_turbulence_parameter:
    case FLD::set_loma_parameter:
    case FLD::set_topopt_parameter:
    case FLD::set_general_adjoint_parameter:
    case FLD::set_adjoint_time_parameter:
      break;
    default:
    {
      std::cout << "action type: " << act << std::endl;
      dserror("Unknown type of action for MeshfreeFluid");
      break;
    }
  } // end of switch(act)

  return 0;
} // end of DRT::ELEMENTS::MeshfreeFluid::Evaluate


/*----------------------------------------------------------------------*
 |  do nothing (public)                                      gammi 04/07|
 |                                                                      |
 |  The function is just a dummy. For fluid elements, the integration   |
 |  integration of volume Neumann conditions (body forces) takes place  |
 |  in the element. We need it there for the stabilisation terms!       |
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::MeshfreeFluid::EvaluateNeumann(
  Teuchos::ParameterList&   params,
  DRT::Discretization&      discretization,
  DRT::Condition&           condition,
  std::vector<int>&         lm,
  Epetra_SerialDenseVector& elevec1,
  Epetra_SerialDenseMatrix* elemat1)
{
  return 0;
}


