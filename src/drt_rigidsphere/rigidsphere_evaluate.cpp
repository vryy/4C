/*----------------------------------------------------------------------------*/
/*!

\brief spherical particle element for brownian dynamics

\level 3

\maintainer Maximilian Grill
*/
/*----------------------------------------------------------------------------*/

#include "rigidsphere.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_exporter.H"
#include "../drt_lib/drt_dserror.H"
#include "../drt_lib/drt_utils.H"
#include "../linalg/linalg_utils.H"
#include "../drt_fem_general/drt_utils_fem_shapefunctions.H"
#include "../drt_mat/stvenantkirchhoff.H"
#include "../linalg/linalg_fixedsizematrix.H"
#include "../drt_fem_general/largerotations.H"
#include "../drt_fem_general/drt_utils_integration.H"
#include "../drt_inpar/inpar_structure.H"
#include <Epetra_CrsMatrix.h>

#include "../drt_inpar/inpar_browniandyn.H"
#include "../drt_lib/standardtypes_cpp.H"
#include "../drt_structure_new/str_elements_paramsinterface.H"

/*-----------------------------------------------------------------------------------------------------------*
 |  evaluate the element (public) meier 02/14|
 *----------------------------------------------------------------------------------------------------------*/
int DRT::ELEMENTS::Rigidsphere::Evaluate(Teuchos::ParameterList& params,
    DRT::Discretization& discretization, std::vector<int>& lm, Epetra_SerialDenseMatrix& elemat1,
    Epetra_SerialDenseMatrix& elemat2, Epetra_SerialDenseVector& elevec1,
    Epetra_SerialDenseVector& elevec2, Epetra_SerialDenseVector& elevec3)
{
  SetParamsInterfacePtr(params);

  // start with "none"
  ELEMENTS::ActionType act = ParamsInterface().GetActionType();

  switch (act)
  {
    case ELEMENTS::struct_calc_linstiff:
    case ELEMENTS::struct_calc_nlnstiff:
    case ELEMENTS::struct_calc_internalforce:
    case ELEMENTS::struct_calc_linstiffmass:
    case ELEMENTS::struct_calc_nlnstiffmass:
    case ELEMENTS::struct_calc_nlnstifflmass:
    case ELEMENTS::struct_calc_internalinertiaforce:
    {
      // need current global displacement and residual forces and get them from discretization
      // making use of the local-to-global map lm one can extract current displacement and residual
      // values for each degree of freedom

      // get element displacements
      Teuchos::RCP<const Epetra_Vector> disp = discretization.GetState("displacement");
      if (disp == Teuchos::null) dserror("Cannot get state vectors 'displacement'");
      std::vector<double> mydisp(lm.size());
      DRT::UTILS::ExtractMyValues(*disp, mydisp, lm);

      Teuchos::RCP<const Epetra_Vector> vel;
      std::vector<double> myvel(lm.size());
      myvel.clear();

      // get element acceleration
      std::vector<double> myacc(lm.size());
      myacc.clear();

      if (act == ELEMENTS::struct_calc_nlnstiffmass or act == ELEMENTS::struct_calc_nlnstifflmass or
          act == ELEMENTS::struct_calc_linstiffmass)
      {
        nlnstiffmass(params, myacc, myvel, mydisp, &elemat1, &elemat2, &elevec1, &elevec2);
      }
      else if (act == ELEMENTS::struct_calc_linstiff or act == ELEMENTS::struct_calc_nlnstiff)
      {
        nlnstiffmass(params, myacc, myvel, mydisp, &elemat1, NULL, &elevec1, NULL);
      }
      else if (act == ELEMENTS::struct_calc_internalforce)
      {
        nlnstiffmass(params, myacc, myvel, mydisp, NULL, NULL, &elevec1, NULL);
      }
      else if (act == ELEMENTS::struct_calc_internalinertiaforce)
      {
        nlnstiffmass(params, myacc, myvel, mydisp, NULL, NULL, &elevec1, &elevec2);
      }
    }
    break;

    case ELEMENTS::struct_calc_brownianforce:
    case ELEMENTS::struct_calc_brownianstiff:
    {
      // get element displacements
      Teuchos::RCP<const Epetra_Vector> disp = discretization.GetState("displacement");
      if (disp == Teuchos::null) dserror("Cannot get state vectors 'displacement'");
      std::vector<double> mydisp(lm.size());
      DRT::UTILS::ExtractMyValues(*disp, mydisp, lm);

      // get element velocity
      Teuchos::RCP<const Epetra_Vector> vel = discretization.GetState("velocity");
      if (vel == Teuchos::null) dserror("Cannot get state vectors 'velocity'");
      std::vector<double> myvel(lm.size());
      DRT::UTILS::ExtractMyValues(*vel, myvel, lm);

      if (act == ELEMENTS::struct_calc_brownianforce)
        CalcBrownianForcesAndStiff(params, myvel, mydisp, NULL, &elevec1);
      else if (act == ELEMENTS::struct_calc_brownianstiff)
        CalcBrownianForcesAndStiff(params, myvel, mydisp, &elemat1, &elevec1);
      else
        dserror("You shouldn't be here.");

      break;
    }

    case ELEMENTS::struct_calc_stress:
    {
      dserror("No stress output implemented for beam3 elements");
      break;
    }
    case ELEMENTS::struct_calc_update_istep:
    case ELEMENTS::struct_calc_reset_istep:
    case ELEMENTS::struct_calc_recover:
    {
      // not necessary since no class variables are modified in predicting steps
      break;
    }

    case ELEMENTS::struct_calc_predict:
    {
      // do nothing here
      break;
    }

    // element based PTC scaling
    case ELEMENTS::struct_calc_addjacPTC:
    {
      // nothing to do here
      break;
    }

    case ELEMENTS::struct_calc_energy:
    {
      // no contribution of rigid sphere to system energy
      break;
    }

    default:
    {
      dserror("Unknown type of action for Rigidsphere %d", act);
      break;
    }
  }

  return (0);
}

/*------------------------------------------------------------------------------------------------------------*
 | nonlinear stiffness and mass matrix (private) meier 05/12|
 *-----------------------------------------------------------------------------------------------------------*/
void DRT::ELEMENTS::Rigidsphere::nlnstiffmass(Teuchos::ParameterList& params,
    std::vector<double>& acc, std::vector<double>& vel, std::vector<double>& disp,
    Epetra_SerialDenseMatrix* stiffmatrix, Epetra_SerialDenseMatrix* massmatrix,
    Epetra_SerialDenseVector* force, Epetra_SerialDenseVector* inertia_force)
{
  // assemble internal force vector if requested
  if (force != NULL)
  {
    for (int i = 0; i < 3; ++i) (*force)(i) = 0.0;
  }

  // assemble stiffmatrix if requested
  if (stiffmatrix != NULL)
  {
    for (int i = 0; i < 3; ++i)
      for (int j = 0; j < 3; ++j) (*stiffmatrix)(i, j) = 0.0;
  }

  // assemble massmatrix if requested
  if (massmatrix != NULL)
  {
    double m = rho_ * 4.0 / 3.0 * PI * radius_ * radius_ * radius_;
    for (int i = 0; i < 3; ++i) (*massmatrix)(i, i) = m;
  }

  //    //assemble inertia force vector if requested
  //    if ( inertia_force != NULL and massmatrix != NULL )
  //    {
  //      for ( int i = 0; i < 3; ++i )
  //        (*inertia_force)(i) = acc[i] * (*massmatrix)(i,i);
  //    }

  return;
}

/*------------------------------------------------------------------------------------------------------------*
 | calculation of thermal (i.e. stochastic) and damping forces according to Brownian dynamics grill
 03/14|
 *------------------------------------------------------------------------------------------------------------*/
void DRT::ELEMENTS::Rigidsphere::CalcBrownianForcesAndStiff(Teuchos::ParameterList& params,
    std::vector<double>& vel, std::vector<double>& disp, Epetra_SerialDenseMatrix* stiffmatrix,
    Epetra_SerialDenseVector* force)
{
  CalcDragForce(params, vel, disp, stiffmatrix, force);
  CalcStochasticForce(params, vel, disp, stiffmatrix, force);
}

/*------------------------------------------------------------------------------------------------------------*
 | compute drag forces and contribution to stiffness matrix  (private) grill 03/14|
 *-----------------------------------------------------------------------------------------------------------*/
void DRT::ELEMENTS::Rigidsphere::CalcDragForce(Teuchos::ParameterList& params,
    const std::vector<double>& vel,         //!< element velocity vector
    const std::vector<double>& disp,        //!< element displacement vector
    Epetra_SerialDenseMatrix* stiffmatrix,  //!< element stiffness matrix
    Epetra_SerialDenseVector* force)        //!< element internal force vector
{
  double gamma = MyDampingConstant();

  // get time step size
  double dt = ParamsInterface().GetDeltaTime();

  // velocity and gradient of background velocity field
  LINALG::Matrix<3, 1> velbackground;
  LINALG::Matrix<3, 3> velbackgroundgrad;  // is a dummy so far

  // Compute background velocity
  GetBackgroundVelocity(params, velbackground, velbackgroundgrad);

  // Drag force contribution
  if (force != NULL)
    for (int i = 0; i < 3; ++i) (*force)(i) += gamma * (vel[i] - velbackground(i));


  // contribution to stiffness matrix
  // depends on TIME INTEGRATION SCHEME (so far, damping is allowed for StatMech only => Backward
  // Euler) GenAlpha would require scaling with gamma_genalpha/beta_genalpha
  if (stiffmatrix != NULL)
  {
    // StatMech: Backward Euler
    for (int l = 0; l < 3; l++)
    {
      (*stiffmatrix)(l, l) += gamma / dt;
    }
  }

  return;
}

/*-----------------------------------------------------------------------------------------------------------*
 |computes velocity of background fluid and gradient of that velocity at a certain evaluation point
 in       | |the physical space                                                         (public)
 grill   03/14|
 *----------------------------------------------------------------------------------------------------------*/
void DRT::ELEMENTS::Rigidsphere::GetBackgroundVelocity(
    Teuchos::ParameterList& params,           //!< parameter list
    LINALG::Matrix<3, 1>& velbackground,      //!< velocity of background fluid
    LINALG::Matrix<3, 3>& velbackgroundgrad)  //!< gradient of velocity of background fluid
{
  // only constant background velocity implemented yet. for case of shear flow, see beam3r


  // default values for background velocity and its gradient
  velbackground.PutScalar(0);
  velbackgroundgrad.PutScalar(0);

  //  double time = params.get<double>("total time",0.0);
  //  double starttime = params.get<double>("STARTTIMEACT",0.0);
  //  double dt = params.get<double>("delta time");
  //
  //  Teuchos::RCP<std::vector<double> > defvalues = Teuchos::rcp(new std::vector<double>(3,0.0));
  //  Teuchos::RCP<std::vector<double> > periodlength = params.get("PERIODLENGTH", defvalues);
  //
  //  // check and throw error if shear flow is applied
  //  INPAR::STATMECH::DBCType dbctype = params.get<INPAR::STATMECH::DBCType>("DBCTYPE",
  //  INPAR::STATMECH::dbctype_std); bool shearflow = false;
  //  if(dbctype==INPAR::STATMECH::dbctype_shearfixed ||
  //     dbctype==INPAR::STATMECH::dbctype_shearfixeddel ||
  //     dbctype==INPAR::STATMECH::dbctype_sheartrans ||
  //     dbctype==INPAR::STATMECH::dbctype_affineshear||
  //     dbctype==INPAR::STATMECH::dbctype_affinesheardel)
  //  {
  //    shearflow = true;
  //    dserror("Shear flow not implemented yet for rigid spherical particles!");
  //  }
  //
  //  // constant background velocity specified in input file?
  //  Teuchos::RCP<std::vector<double> > constbackgroundvel = params.get("CONSTBACKGROUNDVEL",
  //  defvalues);
  //
  //  if (constbackgroundvel->size() != 3) dserror("\nSpecified vector for constant background
  //  velocity has wrong dimension! Check input file!"); bool constflow = false; for (int i=0; i<3;
  //  ++i)
  //  {
  //    if (constbackgroundvel->at(i)!=0.0) constflow=true;
  //  }
  //
  //  if(periodlength->at(0) > 0.0)
  //  {
  //    if(constflow && time>starttime && fabs(time-starttime)>dt/1e4)
  //    {
  //      for (int i=0; i<3; ++i) velbackground(i) = constbackgroundvel->at(i);
  //
  //      // shear flow AND constant background flow not implemented
  //      if(shearflow) dserror("Conflict in input parameters: shearflow AND constant background
  //      velocity specified. Not implemented!\n");
  //    }
  //  }
}

/*-----------------------------------------------------------------------------------------------------------*
 | computes damping coefficient                                             (private) grill   03/14|
 *----------------------------------------------------------------------------------------------------------*/
double DRT::ELEMENTS::Rigidsphere::MyDampingConstant()
{
  // (dynamic) viscosity of background fluid
  double eta = ParamsInterface().GetBrownianDynParamInterface()->GetViscosity();

  // damping/friction coefficient of a rigid sphere (Stokes' law for very small Reynolds numbers)
  return 6 * PI * eta * radius_;
}

/*-----------------------------------------------------------------------------------------------------------*
 |computes the number of different random numbers required in each time step for generation of
 stochastic    | |forces; (public)           grill   03/14|
 *----------------------------------------------------------------------------------------------------------*/
int DRT::ELEMENTS::Rigidsphere::HowManyRandomNumbersINeed()
{
  /*three randomly excited (translational) DOFs for Rigidsphere element*/
  return 3;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void DRT::ELEMENTS::Rigidsphere::GetGeneralizedInterpolationMatrixVariationsAtXi(
    LINALG::SerialDenseMatrix& Ivar, const double& dummy1, const std::vector<double>& dummy2) const
{
  LINALG::TMatrix<double, 6, 3> Ivar_fixedsize(&Ivar(0, 0), true);
  for (unsigned int i = 0; i < 3; ++i) Ivar_fixedsize(i, i) = 1.0;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void DRT::ELEMENTS::Rigidsphere::GetGeneralizedInterpolationMatrixIncrementsAtXi(
    LINALG::SerialDenseMatrix& Iinc, const double& dummy1, const std::vector<double>& dummy2) const
{
  LINALG::TMatrix<double, 6, 3> Iinc_fixedsize(&Iinc(0, 0), true);
  for (unsigned int i = 0; i < 3; ++i) Iinc_fixedsize(i, i) = 1.0;
}

/*-----------------------------------------------------------------------------------------------------------*
 | computes stochastic forces and resulting stiffness (public) grill   03/14|
 *----------------------------------------------------------------------------------------------------------*/
void DRT::ELEMENTS::Rigidsphere::CalcStochasticForce(
    Teuchos::ParameterList& params,         //!< parameter list
    const std::vector<double>& vel,         //!< element velocity vector
    const std::vector<double>& disp,        //!< element disp vector
    Epetra_SerialDenseMatrix* stiffmatrix,  //!< element stiffness matrix
    Epetra_SerialDenseVector* force)        //!< element internal force vector
{
  // damping coefficient
  double gamma = MyDampingConstant();

  /*get pointer at Epetra multivector in parameter list linking to random numbers for stochastic
   * forces with zero mean and standard deviation (2*kT / dt)^0.5*/
  Teuchos::RCP<Epetra_MultiVector> randomnumbers =
      ParamsInterface().GetBrownianDynParamInterface()->GetRandomForces();

  if (force != NULL)
  {
    for (unsigned int k = 0; k < 3; ++k)
    {
      (*force)(k) -= sqrt(gamma) * (*randomnumbers)[k][LID()];
    }
  }

  // no contribution to stiffmatrix

  return;
}
