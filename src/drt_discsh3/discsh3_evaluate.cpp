/*-----------------------------------------------------------*/
/*! \file
\brief evaluate routines for discsh3 element

\level 3

\maintainer Christoph Meier
*/
/*-----------------------------------------------------------*/

#include "discsh3.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_lib/drt_discret.H"
//#include "../drt_lib/drt_element.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_dserror.H"
#include "../linalg/linalg_utils.H"
#include "../drt_mat/stvenantkirchhoff.H"
#include "../linalg/linalg_serialdensevector.H"
#include "Epetra_SerialDenseSolver.h"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_fem_general/drt_utils_integration.H"
#include "../drt_fem_general/drt_utils_fem_shapefunctions.H"

#include <Teuchos_StandardParameterEntryValidators.hpp>

#include "Sacado.hpp"

#include "../drt_inpar/inpar_browniandyn.H"
typedef Sacado::Fad::DFad<double> FAD;


/*----------------------------------------------------------------------*
 |  evaluate the element (public)                        mukherjee 04/15|
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::DiscSh3::Evaluate(Teuchos::ParameterList& params,
    DRT::Discretization& discretization, std::vector<int>& lm, Epetra_SerialDenseMatrix& elemat1,
    Epetra_SerialDenseMatrix& elemat2, Epetra_SerialDenseVector& elevec1,
    Epetra_SerialDenseVector& elevec2, Epetra_SerialDenseVector& elevec3)
{
  // start with "none"
  DRT::ELEMENTS::DiscSh3::ActionType act = DiscSh3::none;

  // get the required action
  std::string action = params.get<std::string>("action", "none");
  if (action == "none")
    dserror("No action supplied");
  else if (action == "calc_struct_linstiff")
    act = DiscSh3::calc_struct_linstiff;
  else if (action == "calc_struct_nlnstiff")
    act = DiscSh3::calc_struct_nlnstiff;
  else if (action == "calc_struct_internalforce")
    act = DiscSh3::calc_struct_internalforce;
  else if (action == "calc_struct_linstiffmass")
    act = DiscSh3::calc_struct_linstiffmass;
  else if (action == "calc_struct_nlnstiffmass")
    act = DiscSh3::calc_struct_nlnstiffmass;
  else if (action == "calc_struct_nlnstifflmass")
    act = DiscSh3::calc_struct_nlnstifflmass;
  else if (action == "calc_struct_ptcstiff")
    act = DiscSh3::calc_struct_ptcstiff;
  else if (action == "calc_struct_stress")
    act = DiscSh3::calc_struct_stress;
  else if (action == "calc_struct_eleload")
    act = DiscSh3::calc_struct_eleload;
  else if (action == "calc_struct_fsiload")
    act = DiscSh3::calc_struct_fsiload;
  else if (action == "calc_struct_refvol")
    act = DiscSh3::calc_struct_refvol;
  else if (action == "calc_struct_currvol")
    act = DiscSh3::calc_struct_currvol;
  else if (action == "calc_struct_refCG")
    act = DiscSh3::calc_struct_refCG;
  else if (action == "calc_struct_currCG")
    act = DiscSh3::calc_struct_currCG;
  else if (action == "calc_struct_refarea")
    act = DiscSh3::calc_struct_refarea;
  else if (action == "calc_struct_currarea")
    act = DiscSh3::calc_struct_currarea;
  else if (action == "calc_struct_update_istep")
    act = DiscSh3::calc_struct_update_istep;
  else if (action == "calc_struct_reset_istep")
    act = DiscSh3::calc_struct_reset_istep;
  else
    dserror("Unknown type of action for DiscSh3: %s", action.c_str());


  // what should the element do
  switch (act)
  {
    case DiscSh3::calc_struct_ptcstiff:
    {
      EvaluatePTC(params, elemat1);
    }
    break;
    //==================================================================================
    // nonlinear stiffness and internal force vector
    case calc_struct_nlnstiff:
    case calc_struct_linstiff:
    case calc_struct_nlnstiffmass:
    case calc_struct_nlnstifflmass:
    case calc_struct_linstiffmass:
    {
      // need current displacement and residual forces
      Teuchos::RCP<const Epetra_Vector> disp = discretization.GetState("displacement");
      Teuchos::RCP<const Epetra_Vector> res = discretization.GetState("residual displacement");
      if (disp == Teuchos::null || res == Teuchos::null)
        dserror("Cannot get state vectors 'displacement' and/or residual");
      std::vector<double> mydisp(lm.size());
      DRT::UTILS::ExtractMyValues(*disp, mydisp, lm);
      std::vector<double> myres(lm.size());
      DRT::UTILS::ExtractMyValues(*res, myres, lm);
      Teuchos::RCP<const Epetra_Vector> vel;
      std::vector<double> myvel(lm.size(), 0.0);
      myvel.clear();
      const Teuchos::ParameterList& sdyn = DRT::Problem::Instance()->StructuralDynamicParams();

      if (DRT::INPUT::IntegralValue<INPAR::STR::DynamicType>(sdyn, "DYNAMICTYP") !=
          INPAR::STR::dyna_statics)
      {
        vel = discretization.GetState("velocity");
        if (vel == Teuchos::null) dserror("Cannot get state vectors 'velocity'");
        DRT::UTILS::ExtractMyValues(*vel, myvel, lm);
      }

      sh3_nlnstiffmass(params, discretization, lm, myvel, mydisp, &elemat1, &elemat2, &elevec1);

      break;
    }


    //==================================================================================
    // internal force vector only
    case calc_struct_internalforce:
    {
      // need current displacement and residual forces
      break;
    }

      //==================================================================================
      // nonlinear stiffness, internal force vector, and consistent mass matrix

      {
        // need current displacement and residual forces
        Teuchos::RCP<const Epetra_Vector> disp = discretization.GetState("displacement");
        Teuchos::RCP<const Epetra_Vector> res = discretization.GetState("residual displacement");
        // need current velocities and accelerations (for non constant mass matrix)
        Teuchos::RCP<const Epetra_Vector> vel = discretization.GetState("velocity");
        Teuchos::RCP<const Epetra_Vector> acc = discretization.GetState("acceleration");
        if (disp == Teuchos::null || res == Teuchos::null)
          dserror("Cannot get state vectors 'displacement' and/or residual");
        if (vel == Teuchos::null) dserror("Cannot get state vectors 'velocity'");
        if (acc == Teuchos::null) dserror("Cannot get state vectors 'acceleration'");

        std::vector<double> mydisp(lm.size());
        DRT::UTILS::ExtractMyValues(*disp, mydisp, lm);
        std::vector<double> myvel(lm.size());
        DRT::UTILS::ExtractMyValues(*vel, myvel, lm);
        std::vector<double> myacc(lm.size());
        DRT::UTILS::ExtractMyValues(*acc, myacc, lm);
        std::vector<double> myres(lm.size());
        DRT::UTILS::ExtractMyValues(*res, myres, lm);

        std::vector<double> mydispmat(lm.size());
        break;
      }


    //==================================================================================
    // evaluate stresses and strains at gauss points
    case calc_struct_stress:
      dserror("Case not yet implemented");
      break;



    //==================================================================================
    case calc_struct_eleload:
      dserror("this method is not supposed to evaluate a load, use EvaluateNeumann(...)");
      break;

    //==================================================================================
    case calc_struct_fsiload:
      dserror("Case not yet implemented");
      break;

    //==================================================================================
    case calc_struct_refvol:
    {
      LINALG::Matrix<1, NUMDOF_DISCSH3> x(true);
      x = MaterialConfiguration();

      // vectors for vertices of pyramid
      LINALG::Matrix<1, 3> v1(true);
      LINALG::Matrix<1, 3> v2(true);
      LINALG::Matrix<1, 3> v3(true);

      for (int j = 0; j < 3; j++)
      {
        v1(j) = x(j);
        v2(j) = x(j + 3);
        v3(j) = x(j + 6);
      }

      double vol_ref =
          (-(v3(0) * v2(1) * v1(2)) + (v2(0) * v3(1) * v1(2)) + (v3(0) * v1(1) * v2(2)) -
              (v1(0) * v3(1) * v2(2)) - (v2(0) * v1(1) * v3(2)) + (v1(0) * v2(1) * v3(2))) /
          6;

      elevec1[0] += vol_ref;
    }
    break;

    //==================================================================================
    case calc_struct_currvol:
    {
      // need current displacement
      Teuchos::RCP<const Epetra_Vector> disp = discretization.GetState("displacement");

      if (disp == Teuchos::null) dserror("Cannot get state vectors 'displacement'");

      std::vector<double> mydisp(lm.size());
      DRT::UTILS::ExtractMyValues(*disp, mydisp, lm);

      LINALG::Matrix<1, NUMDOF_DISCSH3> x(true);
      x = SpatialConfiguration(mydisp);

      // vectors for vertices of pyramid
      LINALG::Matrix<1, 3> v1(true);
      LINALG::Matrix<1, 3> v2(true);
      LINALG::Matrix<1, 3> v3(true);


      for (int j = 0; j < 3; j++)
      {
        v1(j) = x(j);
        v2(j) = x(j + 3);
        v3(j) = x(j + 6);
      }

      double vol_curr =
          (-(v3(0) * v2(1) * v1(2)) + (v2(0) * v3(1) * v1(2)) + (v3(0) * v1(1) * v2(2)) -
              (v1(0) * v3(1) * v2(2)) - (v2(0) * v1(1) * v3(2)) + (v1(0) * v2(1) * v3(2))) /
          6;

      elevec1[0] += vol_curr;
    }
    break;
    //==================================================================================
    case calc_struct_refCG:
    {
      LINALG::Matrix<1, NUMDOF_DISCSH3> x_ref(true);
      x_ref = MaterialConfiguration();

      LINALG::Matrix<1, 3> bary_ref(true);

      for (int i = 0; i < 3; i++) bary_ref(i) = (x_ref(i) + x_ref(i + 3) + x_ref(i + 6)) / 3;


      elevec1[0] += bary_ref(0);
      elevec1[1] += bary_ref(1);
      elevec1[2] += bary_ref(2);
    }
    break;

    //==================================================================================
    case calc_struct_currCG:
    {
      // need current displacement
      Teuchos::RCP<const Epetra_Vector> disp = discretization.GetState("displacement");

      if (disp == Teuchos::null) dserror("Cannot get state vectors 'displacement'");

      std::vector<double> mydisp(lm.size());
      DRT::UTILS::ExtractMyValues(*disp, mydisp, lm);

      // Get Spatial positions
      LINALG::Matrix<1, NUMDOF_DISCSH3> x_curr(true);
      x_curr = SpatialConfiguration(mydisp);

      LINALG::Matrix<1, 3> bary_curr(true);

      for (int i = 0; i < 3; i++) bary_curr(i) = (x_curr(i) + x_curr(i + 3) + x_curr(i + 6)) / 3;


      elevec1[0] += bary_curr(0);
      elevec1[1] += bary_curr(1);
      elevec1[2] += bary_curr(2);
    }
    break;

    //==================================================================================
    case calc_struct_refarea:
    {
      LINALG::Matrix<1, NUMDOF_DISCSH3> x(true);
      x = MaterialConfiguration();

      // vectors for vertices of pyramid
      LINALG::Matrix<1, 3> v1(true);
      LINALG::Matrix<1, 3> v2(true);
      LINALG::Matrix<1, 3> v3(true);

      for (int j = 0; j < 3; j++)
      {
        v1(j) = x(j);
        v2(j) = x(j + 3);
        v3(j) = x(j + 6);
      }

      double area_ref = CalcSurfaceArea(v1, v2, v3);

      elevec1[0] += area_ref;
    }
    break;

    //==================================================================================
    case calc_struct_currarea:
    {
      // need current displacement
      Teuchos::RCP<const Epetra_Vector> disp = discretization.GetState("displacement");

      if (disp == Teuchos::null) dserror("Cannot get state vectors 'displacement'");

      std::vector<double> mydisp(lm.size());
      DRT::UTILS::ExtractMyValues(*disp, mydisp, lm);

      LINALG::Matrix<1, NUMDOF_DISCSH3> x(true);
      x = SpatialConfiguration(mydisp);

      // vectors for vertices of pyramid
      LINALG::Matrix<1, 3> v1(true);
      LINALG::Matrix<1, 3> v2(true);
      LINALG::Matrix<1, 3> v3(true);

      for (int j = 0; j < 3; j++)
      {
        v1(j) = x(j);
        v2(j) = x(j + 3);
        v3(j) = x(j + 6);
      }

      double area_curr = CalcSurfaceArea(v1, v2, v3);

      elevec1[0] += area_curr;
    }
    break;

    //==================================================================================
    case calc_struct_update_istep:
    {
      for (int i = 0; i < NUMDOF_DISCSH3; i++) x_n_1_(i) = x_n_(i);
    }

    break;

    //==================================================================================
    case calc_struct_reset_istep:
      dserror("Case not yet implemented");
      break;


    //==================================================================================
    default:
      dserror("Unknown type of action for DiscSh3");
      break;
  }
  return 0;
}

/*-----------------------------------------------------------------------------------------------------------*
 | Evaluate PTC damping (public) mukherjee 12/15|
 *----------------------------------------------------------------------------------------------------------*/
void DRT::ELEMENTS::DiscSh3::EvaluatePTC(
    Teuchos::ParameterList& params, Epetra_SerialDenseMatrix& elemat1)
{
  // PTC for translational degrees of freedom; the Lobatto integration weight is 0.5 for 3-noded
  // elements
  for (int i = 0; i < 9; i++)
  {
    elemat1(i, i) += params.get<double>("ctransptc", 0.0);
  }

  return;
}  // DRT::ELEMENTS::Beam3r::EvaluatePTC


/*----------------------------------------------------------------------------*
 |  Integrate a Volume Neumann boundary condition (public)     mukherjee 05/15|
 *----------------------------------------------------------------------------*/
int DRT::ELEMENTS::DiscSh3::EvaluateNeumann(Teuchos::ParameterList& params,
    DRT::Discretization& discretization, DRT::Condition& condition, std::vector<int>& lm,
    Epetra_SerialDenseVector& elevec1, Epetra_SerialDenseMatrix* elemat1)
{
  dserror("Method not configured yet!");
  // get values and switches from the condition
  const std::vector<int>* onoff = condition.Get<std::vector<int>>("onoff");
  //   const std::vector<double>* val   = condition.Get<std::vector<double> >("val"  );

  // ensure that at least as many curves/functs as dofs are available
  if (int(onoff->size()) < NUMDIM_DISCSH3)
    dserror("Fewer functions or curves defined than the element has dofs.");

  for (int checkdof = NUMDIM_DISCSH3; checkdof < int(onoff->size()); ++checkdof)
  {
    if ((*onoff)[checkdof] != 0)
      dserror("Number of Dimensions in Neumann_Evalutaion is 3. Further DoFs are not considered.");
  }

  // (SPATIAL) FUNCTION BUSINESS
  const std::vector<int>* funct = condition.Get<std::vector<int>>("funct");
  LINALG::Matrix<NUMDIM_DISCSH3, 1> xrefegp(false);
  bool havefunct = false;
  if (funct)
    for (int dim = 0; dim < NUMDIM_DISCSH3; dim++)
      if ((*funct)[dim] > 0) havefunct = havefunct or true;

  return 0;
}  // DRT::ELEMENTS::DiscSh3::EvaluateNeumann


/*----------------------------------------------------------------------*
 |  evaluate the element (private)                       mukherjee 05/15|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::DiscSh3::sh3_nlnstiffmass(Teuchos::ParameterList& params,
    DRT::Discretization& discretization, std::vector<int>& lm, const std::vector<double>& vel,
    const std::vector<double>& disp, Epetra_SerialDenseMatrix* stiffmatrix,
    Epetra_SerialDenseMatrix* massmatrix, Epetra_SerialDenseVector* force)
{
  x_n_ = SpatialConfiguration(disp);
  //  Teuchos::ParameterList StatMechParams =
  //  DRT::Problem::Instance()->StatisticalMechanicsParams();

  // Calculate the stiffness and viscous contribution arising from area constraint
  const int NumGElements = discretization.NumGlobalElements();

  double time = params.get<double>("total time", 0.0);
  if (time == 0.0)  // Only at first time step
    CheckIfOutwardsNormal(params, NumGElements);

  // Local way of computing constraint on Center of Gravity
  //  CGConstrtStiffmass(params,disp,vel,stiffmatrix,force,NumGElements);

  // Global way of computing constraint on Center of Gravity
  //   CGGlobalConstrtStiffmass(params,disp,vel,stiffmatrix,force,NumGElements);

  // Global way of computing constraint on Barycenter
  //   BaryConstrtStiffmass(params,disp,vel,stiffmatrix,force,NumGElements);

  //  INPAR::STATMECH::AreaPenaltyType area_pen_type =
  //  DRT::INPUT::IntegralValue<INPAR::STATMECH::AreaPenaltyType>(StatMechParams,"AREA_PENALTY_TYPE");
  //
  //  if(area_pen_type==INPAR::STATMECH::areapenalty_local)
  //  {
  //    // Local way of computing constraint on element area
  //    AreaConstrtStiffmass(params,disp,vel,stiffmatrix,force);
  //  }
  //  else if(area_pen_type==INPAR::STATMECH::areapenalty_global)
  //  {
  //    // Global way of computing constraint on element area
  //    AreaConstrtGlobalStiff(params,discretization,disp,stiffmatrix,force);
  //  }

  // Local way of computing constraint on element area (Quadrature based)
  //  AreaConstrtQuadStiffmass(params,disp,vel,stiffmatrix,force);

  // Local way of Calculating the force & stiffness contribution
  // arising from volume constraint (Only applicable in case of enclosed volume)
  //   VolConstrtStiffmass(params,disp,stiffmatrix,force);

  // Global way of Calculating the force & stiffness contribution
  // arising from volume constraint
  VolConstrtGlobalStiff(params, discretization, disp, stiffmatrix, force);

  // Calculate the mass matrix
  sh3_lumpmass(disp, massmatrix);

  // in statistical mechanics simulations, a deletion influenced by the values of the internal force
  // vector might occur
  //  if(params.get<std::string>("internalforces","no")=="yes" && force != NULL)
  //  internalforces_ = *force;
  /*the following function call applied statistical forces and damping matrix according to the
   * fluctuation dissipation theorem; it is dedicated to the application of beam3 elements in the
   * frame of statistical mechanics problems; for these problems a special vector has to be passed
   * to the element packed in the params parameter list; in case that the control routine calling
   * the element does not attach this special vector to params the following method is just doing
   * nothing, which means that for any ordinary problem of structural mechanics it may be ignored*/

  //  // Get if normal dynamics problem or statmech problem
  //  if(DRT::INPUT::IntegralValue<int>(DRT::Problem::Instance()->StatisticalMechanicsParams(),
  //  "BROWNDYNPROB"))
  //  {
  //      #ifdef INEXTENSIBLE
  //        dserror("INEXTENSIBLE formulation not possible for statmech so far. Adapt vector vel ->
  //        myvel like above!");
  //      #endif
  //        CalcBrownian(params,vel,disp,stiffmatrix,force);
  //
  //        // Internal damping forces
  ////          MyDampingForces(params,vel,disp,stiffmatrix,force);
  //    }

  return;
}  // DRT::ELEMENTS::DiscSh3::sh3_nlnstiffmass

/*-----------------------------------------------------------------------------------------------------------*
 | Assemble stochastic and viscous forces and respective stiffness according to fluctuation
 dissipation      | | theorem (private) Mukherjee 08/15|
 *----------------------------------------------------------------------------------------------------------*/
inline void DRT::ELEMENTS::DiscSh3::CalcBrownian(Teuchos::ParameterList& params,
    const std::vector<double>& vel,         //!< element velocity vector
    const std::vector<double>& disp,        //!< element displacement vector
    Epetra_SerialDenseMatrix* stiffmatrix,  //!< element stiffness matrix
    Epetra_SerialDenseVector* force)        //!< element internal force vector
{
  // if no random numbers for generation of stochastic forces are passed to the element no Brownian
  // dynamics calculations are conducted
  if (params.get<Teuchos::RCP<Epetra_MultiVector>>("RandomNumbers", Teuchos::null) == Teuchos::null)
    return;

  // Evaluation of force vectors and stiffness matrices

  // add stiffness and forces due to translational damping effects
  // Approach proposed by Baraff1998
  //  MyDampingForces(params,vel,disp,stiffmatrix,force);

  // Modified approach (Modify damping with disk in fluid)
  MyTranslationalDamping(params, vel, disp, stiffmatrix, force);

  // add stochastic forces and (if required) resulting stiffness
  //  MyStochasticForces(params,vel,disp,stiffmatrix,force);

  return;

}  // DRT::ELEMENTS::DiscSh3::CalcBrownian(.)


/*-----------------------------------------------------------------------------------------------------------*
 | computes translational damping forces and stiffness (private) Mukherjee   08/15|
 *----------------------------------------------------------------------------------------------------------*/
inline void DRT::ELEMENTS::DiscSh3::MyTranslationalDamping(
    Teuchos::ParameterList& params,   //!< parameter list
    const std::vector<double>& vel,   //!< vector containing first order time derivative of nodal
                                      //!< positions and nodal tangents of an element
    const std::vector<double>& disp,  //!< vector containing change in nodal positions and nodal
                                      //!< tangents of an element w.r.t. inital config.
    Epetra_SerialDenseMatrix* stiffmatrix,  //!< element stiffness matrix
    Epetra_SerialDenseVector* force)        //!< element internal force vector
{
  // get time step size
  double dt = params.get<double>("delta time", 0.0);

  // damping coefficients for translational and rotational degrees of freedom
  LINALG::Matrix<NUMDIM_DISCSH3, 1> gamma(true);
  MyDampingConstants(params, gamma);

  // compute normal vector t_{\par} at current Gauss point
  std::vector<FAD> x_FAD(NUMDOF_DISCSH3, 0.0);
  // Get Spatial positions
  LINALG::Matrix<1, NUMDOF_DISCSH3> x = SpatialConfiguration(disp);
  for (int i = 0; i < NUMDOF_DISCSH3; i++)
  {
    x_FAD[i] = x(i);
    x_FAD[i].diff(i, NUMDOF_DISCSH3);
  }

  LINALG::TMatrix<FAD, 1, 3> normal = CalcSurfaceNormal(x_FAD);

  /** Calc damping matrix **/
  LINALG::TMatrix<FAD, 3, 3> Damping_mat_small(true);
  LINALG::TMatrix<FAD, NUMDOF_DISCSH3, NUMDOF_DISCSH3> Damping_mat(true);
  // Compute outer product of normals
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
    {
      Damping_mat_small(i, j) = (i == j) * gamma(1) + (gamma(0) - gamma(1)) * normal(i) * normal(j);
      Damping_mat(i, j) = Damping_mat_small(i, j);
      Damping_mat(i + 3, j + 3) = Damping_mat_small(i, j);
      Damping_mat(i + 6, j + 6) = Damping_mat_small(i, j);
    }

  LINALG::TMatrix<FAD, 1, NUMDOF_DISCSH3> FAD_force(true);
  LINALG::TMatrix<FAD, NUMDOF_DISCSH3, NUMDOF_DISCSH3> FAD_StiffMat(true);

  // loop over columns of matrix t_{\par} \otimes t_{\par}

  LINALG::TMatrix<LINALG::TMatrix<FAD, NUMDOF_DISCSH3, NUMDOF_DISCSH3>, 1, NUMDOF_DISCSH3>
      Grad_DampingMat(true);
  for (int i = 0; i < NUMDOF_DISCSH3; i++)
    for (int j = 0; j < NUMDOF_DISCSH3; j++)
      for (int k = 0; k < NUMDOF_DISCSH3; k++) Grad_DampingMat(k)(i, j) = Damping_mat(i, j).dx(k);

  for (int i = 0; i < NUMDOF_DISCSH3; i++)
    for (int j = 0; j < NUMDOF_DISCSH3; j++)
    {
      if (force != NULL)
        FAD_force(i) += Damping_mat(i, j) * vel[j] / 3;  // 1/3 because of the nodal weight
      if (stiffmatrix != NULL)
        for (int k = 0; k < NUMDOF_DISCSH3; k++)
        {
          FAD_StiffMat(i, j) +=
              Grad_DampingMat(k)(i, j) * vel[k] / 3 + Damping_mat(i, j) / (3 * dt);
        }
    }

  for (int i = 0; i < NUMDOF_DISCSH3; i++)
  {
    if (force != NULL) (*force)(i) += FAD_force(i).val();
    for (int j = 0; j < NUMDOF_DISCSH3; j++)
    {
      if (stiffmatrix != NULL) (*stiffmatrix)(i, j) += FAD_StiffMat(i, j).val();
    }
  }

  return;
}  // DRT::ELEMENTS::DiscSh3::MyTranslationalDamping(.)


/*-----------------------------------------------------------------------------------------------------------*
 | computes stochastic forces and resulting stiffness (private) mukherjee   10/15|
 *----------------------------------------------------------------------------------------------------------*/
inline void DRT::ELEMENTS::DiscSh3::MyStochasticForces(
    Teuchos::ParameterList& params,         //!< parameter list
    const std::vector<double>& vel,         //!< element velocity vector
    const std::vector<double>& disp,        //!< element disp vector
    Epetra_SerialDenseMatrix* stiffmatrix,  //!< element stiffness matrix
    Epetra_SerialDenseVector* force)        //!< element internal force vector
{
  // damping coefficients for three translational and one rotatinal degree of freedom
  LINALG::Matrix<3, 1> gamma(true);

  MyDampingConstants(params, gamma);


  /*get pointer at Epetra multivector in parameter list linking to random numbers for stochastic
   * forces with zero mean and standard deviation (2*kT / dt)^0.5; note carefully: a space between
   * the two subsequal ">" signs is mandatory for the C++ parser in order to avoid confusion with
   * ">>" for streams*/
  Teuchos::RCP<Epetra_MultiVector> randomnumbers =
      params.get<Teuchos::RCP<Epetra_MultiVector>>("RandomNumbers", Teuchos::null);

  // compute normal vector t_{\par} at current Gauss point
  std::vector<FAD> x_FAD(NUMDOF_DISCSH3, 0.0);
  // Get Spatial positions
  LINALG::Matrix<1, NUMDOF_DISCSH3> x(true);
  LINALG::Matrix<1, NUMDOF_DISCSH3> x_mat(true);
  x = SpatialConfiguration(disp);
  for (int i = 0; i < NUMDOF_DISCSH3; i++)
  {
    x_FAD[i] = x(i);
    x_FAD[i].diff(i, NUMDOF_DISCSH3);
  }

  LINALG::TMatrix<FAD, 1, 3> normal = CalcSurfaceNormal(x_FAD);

  /** Calc damping matrix **/
  LINALG::TMatrix<FAD, 3, 3> Stoch_mat_small(true);
  LINALG::TMatrix<FAD, NUMDOF_DISCSH3, NUMDOF_DISCSH3> Stoch_mat(true);
  // Compute outer product of normals
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
    {
      Stoch_mat_small(i, j) = (i == j) * std::sqrt(gamma(1)) +
                              (std::sqrt(gamma(0)) - std::sqrt(gamma(1))) * normal(i) * normal(j);
      Stoch_mat(i, j) = Stoch_mat_small(i, j);
      Stoch_mat(i + 3, j + 3) = Stoch_mat_small(i, j);
      Stoch_mat(i + 6, j + 6) = Stoch_mat_small(i, j);
    }

  LINALG::TMatrix<FAD, 1, NUMDOF_DISCSH3> FAD_force(true);
  LINALG::TMatrix<FAD, NUMDOF_DISCSH3, NUMDOF_DISCSH3> FAD_StiffMat(true);

  LINALG::TMatrix<LINALG::TMatrix<FAD, NUMDOF_DISCSH3, NUMDOF_DISCSH3>, 1, NUMDOF_DISCSH3>
      Grad_StochMat(true);
  for (int i = 0; i < NUMDOF_DISCSH3; i++)
    for (int j = 0; j < NUMDOF_DISCSH3; j++)
      for (int k = 0; k < NUMDOF_DISCSH3; k++) Grad_StochMat(k)(i, j) = Stoch_mat(i, j).dx(k);

  // Sort Random numbers
  LINALG::Matrix<1, NUMDOF_DISCSH3> RandomNumber(true);
  for (int i = 0; i < NUMNOD_DISCSH3; i++)
    for (int j = 0; j < NUMDIM_DISCSH3; j++)
    {
      RandomNumber(3 * i + j) = (*randomnumbers)[j][this->Nodes()[i]->LID()];
    }


  for (int i = 0; i < NUMDOF_DISCSH3; i++)
    for (int j = 0; j < NUMDOF_DISCSH3; j++)
    {
      if (force != NULL)
        FAD_force(i) += Stoch_mat(i, j) * RandomNumber(j) / 3;  // weight 1/3 per node

      if (stiffmatrix != NULL)  //    loop over all column nodes
        for (int k = 0; k < NUMNOD_DISCSH3; k++)
        {
          FAD_StiffMat(i, j) += Grad_StochMat(k)(i, j) * RandomNumber(j) / 3;
        }
    }

  // loop dimensions with respect to columns
  for (int i = 0; i < NUMDIM_DISCSH3; i++)
    if (force != NULL)
    {
      (*force)(i) -= FAD_force(i).val();
    }


  for (int i = 0; i < NUMDIM_DISCSH3; i++)
    for (int j = 0; j < NUMDIM_DISCSH3; j++)
    {
      if (stiffmatrix != NULL) (*stiffmatrix)(i, j) -= FAD_StiffMat(i, j).val();
    }


  return;
}  // DRT::ELEMENTS::DiscSh3::MyStochasticForces(.)


/*-----------------------------------------------------------------------------------------------------------*
 |computes the number of different random numbers required in each time step for generation of
 stochastic    | |forces; (public)       Mukherjee   10/15|
 *----------------------------------------------------------------------------------------------------------*/
int DRT::ELEMENTS::DiscSh3::HowManyRandomNumbersINeed()
{
  /*at each node one needs as many random numbers as randomly excited degrees of freedom, i.e. three
   *random numbers for the translational degrees of freedom and one random number for the rotation
   *around the element axis*/
  return (3 * NumNode());
}
/*-----------------------------------------------------------------------------------------------------------*
 |computes the number of different random numbers required in each time step for generation of
 stochastic    | |forces; (public)       Mukherjee   10/15|
 *----------------------------------------------------------------------------------------------------------*/
int DRT::ELEMENTS::DiscSh3::HowManyRandomNumbersPerNode()
{
  /*at each node one needs as many random numbers as randomly excited degrees of freedom, i.e. three
   *random numbers for the translational degrees of freedom and one random number for the rotation
   *around the element axis*/
  return (3);
}


/*-----------------------------------------------------------------------------------------------------------*
 | computes damping coefficients per lengthand stores them in a matrix in the following order:
 damping of    | | translation parallel to filament axis, damping of translation orthogonal to
 filament axis, damping of     | | rotation around filament axis (private)       Mukherjee  10/15|
 *----------------------------------------------------------------------------------------------------------*/
inline void DRT::ELEMENTS::DiscSh3::MyDampingConstants(
    Teuchos::ParameterList& params, LINALG::Matrix<3, 1>& gamma)
{
  // translational damping coefficients according to Ahmadi
  gamma(0) = 16 * params.get<double>("ETA", 0.0);      // Parallel to normal direction
  gamma(1) = 32 / 3 * params.get<double>("ETA", 0.0);  // Perpendicular to normal direction

  gamma(0) = gamma(1);
}

/*-----------------------------------------------------------------------------------------------------------*
 | computes translational damping forces and stiffness (private) Mukherjee   10/15|
 *----------------------------------------------------------------------------------------------------------*/
inline void DRT::ELEMENTS::DiscSh3::MyDampingForces(
    Teuchos::ParameterList& params,   //!< parameter list
    const std::vector<double>& vel,   //!< vector containing first order time derivative of nodal
                                      //!< positions and nodal tangents of an element
    const std::vector<double>& disp,  //!< vector containing change in nodal positions and nodal
                                      //!< tangents of an element w.r.t. inital config.
    Epetra_SerialDenseMatrix* stiffmatrix,  //!< element stiffness matrix
    Epetra_SerialDenseVector* force)        //!< element internal force vector
{
  // Calculate damping force contribution from curvature constraint
  CalcDampingForces_ = true;

  // Calculate damping force contribution from area constraint

  MyDampingForcesArea(params, vel, disp, stiffmatrix, force);

  // Calculate damping force contribution from volume constraint

  //  MyDampingForcesVol(params,vel,disp,stiffmatrix,force);


  return;
}  // DRT::ELEMENTS::DiscSh3::MyDampingForces()

/*-----------------------------------------------------------------------------------------------------------*
 | computes damping forces arising from strectch (area) (private) Mukherjee   10/15|
 *----------------------------------------------------------------------------------------------------------*/
inline void DRT::ELEMENTS::DiscSh3::MyDampingForcesArea(
    Teuchos::ParameterList& params,   //!< parameter list
    const std::vector<double>& vel,   //!< vector containing first order time derivative of nodal
                                      //!< positions and nodal tangents of an element
    const std::vector<double>& disp,  //!< vector containing change in nodal positions and nodal
                                      //!< tangents of an element w.r.t. inital config.
    Epetra_SerialDenseMatrix* stiffmatrix,  //!< element stiffness matrix
    Epetra_SerialDenseVector* force)        //!< element internal force vector
{
  // get time step size
  FAD dt = params.get<double>("delta time", 0.0);

  // damping coefficients for translational and rotational degrees of freedom
  FAD gamma = params.get<double>("ETA", 0.0);


  /****** Calculate damping force contribution from area constraint ******/

  FAD BehaviorFunctArea;

  CalcBehaviorFunctArea(params, disp, BehaviorFunctArea);

  LINALG::TMatrix<FAD, 1, 9> FAD_force_viscous(true);
  for (int i = 0; i < NUMDOF_DISCSH3; i++)
  {
    for (int j = 0; j < NUMDOF_DISCSH3; j++)
    {
      FAD_force_viscous(i) += gamma * BehaviorFunctArea.dx(i) * BehaviorFunctArea.dx(j) * vel[j];
    }
    (*force)(i) += FAD_force_viscous(i).val();
  }


  /* Analytical Calculations */

  std::vector<FAD> x_FAD(NUMDOF_DISCSH3, 0.0);
  // Get Spatial positions
  LINALG::Matrix<1, NUMDOF_DISCSH3> x(true);
  x = SpatialConfiguration(disp);
  LINALG::TMatrix<FAD, 1, 9> crossprod2(true);
  for (int i = 0; i < NUMDOF_DISCSH3; i++)
  {
    x_FAD[i] = x(i);
    x_FAD[i].diff(i, NUMDOF_DISCSH3);
  }

  // Calculate analytical expressions for derivative of behavior function
  LINALG::TMatrix<FAD, 1, 9> BehaviorFunctArea_dx(true);

  for (int j = 0; j < NUMNOD_DISCSH3; j++)
  {
    // For calculation of cross-product
    LINALG::TMatrix<FAD, 1, 3> side1(true);
    LINALG::TMatrix<FAD, 1, 3> side2(true);
    // Auxilary vector
    LINALG::TMatrix<FAD, 1, 3> side3(true);

    if (j == 0)  // 1st node
    {
      for (int i = 0; i < NODDOF_DISCSH3; i++)
      {
        side1(i) = x_FAD[3 + i] - x_FAD[i];
        side2(i) = x_FAD[6 + i] - x_FAD[3 + i];
        side3(i) = -x_FAD[6 + i] + x_FAD[3 + i];
      }
    }
    else if (j == 1)  // 2nd node
    {
      for (int i = 0; i < NODDOF_DISCSH3; i++)
      {
        side1(i) = x_FAD[6 + i] - x_FAD[3 + i];
        side2(i) = x_FAD[i] - x_FAD[6 + i];
        side3(i) = -x_FAD[i] + x_FAD[6 + i];
      }
    }
    else  // 3rd node
    {
      for (int i = 0; i < NODDOF_DISCSH3; i++)
      {
        side1(i) = x_FAD[i] - x_FAD[6 + i];
        side2(i) = x_FAD[3 + i] - x_FAD[i];
        side3(i) = -x_FAD[3 + i] + x_FAD[i];
      }
    }

    LINALG::TMatrix<FAD, 1, 3> crossprod1(true);

    // Cross Product side1xside2
    crossprod1 = CalcCrossProduct(side1, side2);

    FAD norm_crossprod = 0.0;
    for (int i = 0; i < 3; i++)
    {
      norm_crossprod += pow(crossprod1(i), 2);
    }
    norm_crossprod = pow(norm_crossprod, 0.5);

    LINALG::TMatrix<FAD, 1, 3> crossprod_aux(true);

    // Cross Product side3xcrossprod1
    crossprod_aux = CalcCrossProduct(side3, crossprod1);

    for (int i = 0; i < 3; i++) crossprod2(3 * j + i) = crossprod_aux(i);
  }

  // Calculate surface area at reference config FAD
  FAD area_ref = CalcSurfaceArea(disp, true);

  // Calculate surface area at spatial config FAD
  FAD area_curr = CalcSurfaceArea(disp, false);

  // Calculate the analytical expression of force
  for (int i = 0; i < NUMDOF_DISCSH3; i++)
    BehaviorFunctArea_dx(i) = 0.25 * crossprod2(i) / (area_curr);


  LINALG::TMatrix<FAD, 9, 9> stiffmat_viscous(true);
  for (int i = 0; i < NUMDOF_DISCSH3; i++)
    for (int j = 0; j < NUMDOF_DISCSH3; j++)
    {
      for (int k = 0; k < NUMDOF_DISCSH3; k++)
      {
        stiffmat_viscous(i, j) +=
            gamma * BehaviorFunctArea_dx(i).dx(j) * BehaviorFunctArea_dx(k) * vel[k] +
            gamma * BehaviorFunctArea_dx(i) * BehaviorFunctArea_dx(k).dx(j) * vel[k] +
            gamma * BehaviorFunctArea_dx(i) * BehaviorFunctArea_dx(k) * (j == k) / dt;
      }

      (*stiffmatrix)(i, j) += stiffmat_viscous(i, j).val();
    }

  return;
}  // DRT::ELEMENTS::DiscSh3::MyDampingForcesArea()

/*-----------------------------------------------------------------------------------------------------------*
 | computes damping forces arising from strectch (area) (private                           Mukherjee
 10/15|
 *----------------------------------------------------------------------------------------------------------*/
inline void DRT::ELEMENTS::DiscSh3::MyDampingForcesVol(
    Teuchos::ParameterList& params,   //!< parameter list
    const std::vector<double>& vel,   //!< vector containing first order time derivative of nodal
                                      //!< positions and nodal tangents of an element
    const std::vector<double>& disp,  //!< vector containing change in nodal positions and nodal
                                      //!< tangents of an element w.r.t. inital config.
    Epetra_SerialDenseMatrix* stiffmatrix,  //!< element stiffness matrix
    Epetra_SerialDenseVector* force)        //!< element internal force vector
{
  // get time step size
  FAD dt = params.get<double>("delta time", 0.0);

  // damping coefficients for translational and rotational degrees of freedom
  FAD gamma = params.get<double>("ETA", 0.0);

  LINALG::TMatrix<FAD, 1, 9> FAD_force_viscous(true);
  LINALG::TMatrix<FAD, 9, 9> stiffmat_viscous(true);


  /****** Calculate damping force contribution from volume constraint ******/

  FAD BehaviorFunctVol;

  CalcBehaviorFunctVol(params, disp, BehaviorFunctVol);

  for (int i = 0; i < NUMDOF_DISCSH3; i++)
  {
    for (int j = 0; j < NUMDOF_DISCSH3; j++)
    {
      FAD_force_viscous(i) += gamma * BehaviorFunctVol.dx(i) * BehaviorFunctVol.dx(j) * vel[j];
    }
    (*force)(i) += FAD_force_viscous(i).val();
  }

  std::vector<FAD> x_FAD(NUMDOF_DISCSH3, 0.0);
  // Get Spatial positions
  LINALG::Matrix<1, NUMDOF_DISCSH3> x(true);
  x = SpatialConfiguration(disp);
  LINALG::TMatrix<FAD, 1, 9> crossprod(true);
  LINALG::TMatrix<FAD, 1, 9> sign(true);
  for (int i = 0; i < NUMDOF_DISCSH3; i++)
  {
    x_FAD[i] = x(i);
    x_FAD[i].diff(i, NUMDOF_DISCSH3);
  }

  // Calculate analytical expressions for derivative of behavior function
  LINALG::TMatrix<FAD, 1, 9> BehaviorFunctVol_dx(true);
  LINALG::Matrix<1, 9> tol(true);

  for (int j = 0; j < NUMNOD_DISCSH3; j++)
  {
    // vectors for vertices of pyramid
    LINALG::TMatrix<FAD, 1, 3> vertex1(true);
    LINALG::TMatrix<FAD, 1, 3> vertex2(true);
    LINALG::TMatrix<FAD, 1, 3> vertex3(true);

    if (j == 0)  // 1st node
    {
      for (int i = 0; i < NODDOF_DISCSH3; i++)
      {
        vertex1(i) = x_FAD[3 + i];
        vertex2(i) = x_FAD[6 + i];
        vertex3(i) = x_FAD[i];
      }
    }
    else if (j == 1)  // 2nd node
    {
      for (int i = 0; i < NODDOF_DISCSH3; i++)
      {
        vertex1(i) = x_FAD[6 + i];
        vertex2(i) = x_FAD[i];
        vertex3(i) = x_FAD[3 + i];
      }
    }
    else  // 3rd node
    {
      for (int i = 0; i < NODDOF_DISCSH3; i++)
      {
        vertex1(i) = x_FAD[i];
        vertex2(i) = x_FAD[3 + i];
        vertex3(i) = x_FAD[6 + i];
      }
    }
    LINALG::TMatrix<FAD, 1, 3> crossprod_aux(true);

    // Cross Product side3xcrossprod1
    crossprod_aux = CalcCrossProduct(vertex1, vertex2);

    for (int i = 0; i < 3; i++) crossprod(3 * j + i) = crossprod_aux(i);

    FAD dotprod = 0;
    // Dot Product
    for (int i = 0; i < NUMNOD_DISCSH3; i++)
    {
      dotprod += crossprod(3 * j + i) * vertex3(i);
    }
    if (dotprod != 0)
    {
      sign(3 * j + 0) = dotprod / abs(dotprod);
      sign(3 * j + 1) = dotprod / abs(dotprod);
      sign(3 * j + 2) = dotprod / abs(dotprod);
    }
  }

  FAD vol_ref;
  // Calculate volume at reference config FAD
  CalcVolume(vol_ref, disp, true);

  FAD vol_curr;
  // Calculate surface area at spatial config FAD
  CalcVolume(vol_curr, disp, false);

  // Calculate the analytical expression of force
  for (int i = 0; i < NUMDOF_DISCSH3; i++)
  {
    BehaviorFunctVol_dx(i) = crossprod(i) * sign(i) / 6;
    tol(i) = BehaviorFunctVol_dx(i).val() - FAD_force_viscous(i).val();
  }


  for (int i = 0; i < NUMDOF_DISCSH3; i++)
    for (int j = 0; j < NUMDOF_DISCSH3; j++)
    {
      for (int k = 0; k < NUMDOF_DISCSH3; k++)
      {
        stiffmat_viscous(i, j) +=
            gamma * BehaviorFunctVol_dx(i).dx(j) * BehaviorFunctVol_dx(k) * vel[k] +
            gamma * BehaviorFunctVol_dx(i) * BehaviorFunctVol_dx(k).dx(j) * vel[k] +
            gamma * BehaviorFunctVol_dx(i) * BehaviorFunctVol_dx(k) * (j == k) / dt;
      }

      (*stiffmatrix)(i, j) += stiffmat_viscous(i, j).val();
    }


  return;

}  // DRT::ELEMENTS::DiscSh3::MyDampingForcesVol()

/*----------------------------------------------------------------------*
 |  evaluate the element (private)                       mukherjee 06/15|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::DiscSh3::CGConstrtStiffmass(Teuchos::ParameterList& params,
    const std::vector<double>& disp, const std::vector<double>& vel,
    Epetra_SerialDenseMatrix* stiffmatrix, Epetra_SerialDenseVector* force, const int NumGlobalEles)
{
  FAD BehaviorFunctArea;

  //  Teuchos::ParameterList StatMechParams =
  //  DRT::Problem::Instance()->StatisticalMechanicsParams(); FAD
  //  cg_penalty=StatMechParams.get<double>("CG_PENALTY",0.0); if (cg_penalty==0)
  //    dserror("Please enter a non-zero CG_PENALTY");

  FAD cg_penalty = 0.0;

  std::vector<FAD> x_FAD(NUMDOF_DISCSH3, 0.0);
  // Get Spatial positions
  LINALG::Matrix<1, NUMDOF_DISCSH3> x_curr(true);
  LINALG::Matrix<1, NUMDOF_DISCSH3> x_ref(true);

  LINALG::TMatrix<FAD, 1, 3> bary_curr(true);
  LINALG::TMatrix<FAD, 1, 3> bary_ref(true);

  x_ref = MaterialConfiguration();
  x_curr = SpatialConfiguration(disp);

  LINALG::TMatrix<FAD, 1, 9> crossprod2(true);
  for (int i = 0; i < NUMDOF_DISCSH3; i++)
  {
    x_FAD[i] = x_curr(i);
    x_FAD[i].diff(i, NUMDOF_DISCSH3);
  }

  for (int i = 0; i < 3; i++)
  {
    bary_curr(i) = (x_FAD[i] + x_FAD[i + 3] + x_FAD[i + 6]) / 3;
    bary_ref(i) = (x_ref(i) + x_ref(i + 3) + x_ref(i + 6)) / 3;
  }

  FAD EnergyFuctional = 0;


  LINALG::TMatrix<FAD, 1, 9> aux_vect(true);
  for (int i = 0; i < 3; i++)
  {
    EnergyFuctional += cg_penalty * (bary_curr(i) - bary_ref(i)) * (bary_curr(i) - bary_ref(i)) /
                       (2 * NumGlobalEles);
    aux_vect(i) = (bary_curr(i) - bary_ref(i)) / 3;
    aux_vect(i + 3) = (bary_curr(i) - bary_ref(i)) / 3;
    aux_vect(i + 6) = (bary_curr(i) - bary_ref(i)) / 3;
  }


  // Calculate analytical expressions for FAD force
  LINALG::Matrix<1, 9> tol(true);
  LINALG::TMatrix<FAD, 1, 9> FAD_force(true);
  LINALG::TMatrix<FAD, 1, 9> ana_force(true);

  // Calculate the FAD expression of force
  for (int i = 0; i < NUMDOF_DISCSH3; i++)
  {
    ana_force(i) = cg_penalty * aux_vect(i) / (NumGlobalEles);
    FAD_force(i) = EnergyFuctional.dx(i);
    //     tol(i)= ana_force(i).val()-FAD_force(i).val();
  }


  LINALG::Matrix<9, 9> FAD_stiff_val(true);
  for (int i = 0; i < NUMDOF_DISCSH3; i++)
    for (int j = 0; j < NUMDOF_DISCSH3; j++)
    {
      FAD_stiff_val(i, j) = ana_force(i).dx(j);
    }

  return;
}  // DRT::ELEMENTS::DiscSh3::AreaConstrtStiffmass


/*----------------------------------------------------------------------*
 |  evaluate the element (private)                       mukherjee 07/15|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::DiscSh3::CGGlobalConstrtStiffmass(Teuchos::ParameterList& params,
    const std::vector<double>& disp, const std::vector<double>& vel,
    Epetra_SerialDenseMatrix* stiffmatrix, Epetra_SerialDenseVector* force, const int NumGlobalEles)
{
  FAD BehaviorFunctArea;

  //  Teuchos::ParameterList StatMechParams =
  //  DRT::Problem::Instance()->StatisticalMechanicsParams(); double
  //  cg_penalty=StatMechParams.get<double>("CG_PENALTY",0.0); if (cg_penalty==0)
  //    dserror("Please enter a non-zero CG_PENALTY");
  double cg_penalty = 0.0;

  Teuchos::RCP<Epetra_SerialDenseVector> CG_ref_rcp =
      params.get<Teuchos::RCP<Epetra_SerialDenseVector>>("reference CG", Teuchos::null);
  Teuchos::RCP<Epetra_SerialDenseVector> CG_curr_rcp =
      params.get<Teuchos::RCP<Epetra_SerialDenseVector>>("current CG", Teuchos::null);

  FAD EnergyFuctional = 0;

  LINALG::Matrix<1, 3> CG_curr(true);
  LINALG::Matrix<1, 3> CG_ref(true);
  for (int i = 0; i < 3; i++)
  {
    CG_ref(i) = (*CG_ref_rcp)(i) / NumGlobalEles;
    CG_curr(i) = (*CG_curr_rcp)(i) / NumGlobalEles;
  }


  LINALG::Matrix<1, 9> aux_vect(true);
  for (int i = 0; i < 3; i++)
  {
    EnergyFuctional += 0.5 * cg_penalty * (CG_curr(i) - CG_ref(i)) * (CG_curr(i) - CG_ref(i));
    aux_vect(i) = (CG_curr(i) - CG_ref(i)) / 3;
    aux_vect(i + 3) = (CG_curr(i) - CG_ref(i)) / 3;
    aux_vect(i + 6) = (CG_curr(i) - CG_ref(i)) / 3;
  }

  // Calculate analytical expressions for FAD force
  LINALG::Matrix<1, 9> tol(true);
  LINALG::TMatrix<FAD, 1, 9> FAD_force(true);
  LINALG::Matrix<1, 9> ana_force(true);
  LINALG::Matrix<9, 9> stiff_aux(true);

  // Calculate the FAD expression of force
  for (int i = 0; i < NUMDOF_DISCSH3; i++)
  {
    for (int j = 0; j < 9; j++) ana_force(i) = cg_penalty * aux_vect(i);
  }

  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
    {
      stiff_aux(i, j) = cg_penalty / 9;
      stiff_aux(i + 3, j + 3) = cg_penalty / 9;
      stiff_aux(i + 6, j + 6) = cg_penalty / 9;
    }

  LINALG::Matrix<9, 9> FAD_stiff_val(true);
  for (int i = 0; i < NUMDOF_DISCSH3; i++)
  {
    (*force)(i) += ana_force(i);
    for (int j = 0; j < NUMDOF_DISCSH3; j++)
    {
      (*stiffmatrix)(i, j) += stiff_aux(i, j);
    }
  }

  return;
}  // DRT::ELEMENTS::DiscSh3::CGGlobalConstrtStiffmass

/*----------------------------------------------------------------------*
 |  evaluate the element (private)                       mukherjee 09/15|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::DiscSh3::BaryConstrtStiffmass(Teuchos::ParameterList& params,
    const std::vector<double>& disp, const std::vector<double>& vel,
    Epetra_SerialDenseMatrix* stiffmatrix, Epetra_SerialDenseVector* force, const int NumGlobalEles)
{
  //  Teuchos::ParameterList StatMechParams =
  //  DRT::Problem::Instance()->StatisticalMechanicsParams(); FAD
  //  bary_penalty=StatMechParams.get<double>("BARY_PENALTY",0.0); if (bary_penalty==0)
  //    dserror("Please enter a non-zero BARY_PENALTY");

  FAD bary_penalty = 0.0;

  std::vector<FAD> x_FAD(NUMDOF_DISCSH3, 0.0);
  // Get Spatial positions
  LINALG::Matrix<1, NUMDOF_DISCSH3> x_curr(true);
  LINALG::Matrix<1, NUMDOF_DISCSH3> x_ref(true);

  LINALG::TMatrix<FAD, 1, 3> bary_curr(true);
  LINALG::TMatrix<FAD, 1, 3> bary_ref(true);


  // Get coordinates at previous time-step
  x_ref = this->x_n_1_;
  x_curr = SpatialConfiguration(disp);

  LINALG::TMatrix<FAD, 1, 9> crossprod2(true);
  for (int i = 0; i < NUMDOF_DISCSH3; i++)
  {
    x_FAD[i] = x_curr(i);
    x_FAD[i].diff(i, NUMDOF_DISCSH3);
  }

  for (int i = 0; i < 3; i++)
  {
    bary_curr(i) = (x_FAD[i] + x_FAD[i + 3] + x_FAD[i + 6]) / 3;
    bary_ref(i) = (x_ref(i) + x_ref(i + 3) + x_ref(i + 6)) / 3;
  }

  FAD EnergyFuctional = 0;
  LINALG::TMatrix<FAD, 1, 9> aux_vect(true);
  for (int i = 0; i < 3; i++)
  {
    EnergyFuctional +=
        bary_penalty * (bary_curr(i) - bary_ref(i)) * (bary_curr(i) - bary_ref(i)) / 2;
    aux_vect(i) = (bary_curr(i) - bary_ref(i)) / 3;
    aux_vect(i + 3) = (bary_curr(i) - bary_ref(i)) / 3;
    aux_vect(i + 6) = (bary_curr(i) - bary_ref(i)) / 3;
  }

  // Calculate analytical expressions for FAD force
  LINALG::Matrix<1, 9> tol(true);

  LINALG::TMatrix<FAD, 1, 9> FAD_force(true);
  LINALG::TMatrix<FAD, 1, 9> ana_force(true);
  // Calculate the FAD expression of force
  for (int i = 0; i < NUMDOF_DISCSH3; i++)
  {
    ana_force(i) = bary_penalty * aux_vect(i);
    FAD_force(i) = EnergyFuctional.dx(i);
    tol(i) = ana_force(i).val() - FAD_force(i).val();
  }

  LINALG::Matrix<9, 9> FAD_stiff_val(true);
  for (int i = 0; i < NUMDOF_DISCSH3; i++)
  {
    (*force)(i) += ana_force(i).val();
    for (int j = 0; j < NUMDOF_DISCSH3; j++)
    {
      (*stiffmatrix)(i, j) += ana_force(i).dx(j);
    }
  }

  return;
}  // DRT::ELEMENTS::DiscSh3::AreaConstrtStiffmass

/*----------------------------------------------------------------------*
 |  evaluate the element (private)                       mukherjee 06/15|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::DiscSh3::AreaConstrtStiffmass(Teuchos::ParameterList& params,
    const std::vector<double>& disp, const std::vector<double>& vel,
    Epetra_SerialDenseMatrix* stiffmatrix, Epetra_SerialDenseVector* force)
{
  FAD BehaviorFunctArea;

  //  Teuchos::ParameterList StatMechParams =
  //  DRT::Problem::Instance()->StatisticalMechanicsParams(); FAD
  //  area_penalty=StatMechParams.get<double>("AREA_PENALTY",0.0); if (area_penalty==0)
  //    dserror("Please enter a non-zero AREA_PENALTY");

  FAD area_penalty = 0.0;

  CalcBehaviorFunctArea(params, disp, BehaviorFunctArea);

  LINALG::TMatrix<FAD, 1, 9> FAD_force(true);
  //  LINALG::TMatrix<FAD,1,9> FAD_force_viscous(true);
  // Calculate FAD force
  for (int i = 0; i < NUMDOF_DISCSH3; i++)
  {
    FAD_force(i) = area_penalty * BehaviorFunctArea * BehaviorFunctArea.dx(i);
    //    FAD_force(i) = BehaviorFunctArea.dx(i);
    //    if(CalcDampingForces)
    //    {
    //      FAD v_FAD=vel[i];
    //      FAD_force_viscous(i) = gamma*BehaviorFunctArea.dx(i)*BehaviorFunctArea*v_FAD;
    //    }
  }


  std::vector<FAD> x_FAD(NUMDOF_DISCSH3, 0.0);
  // Get Spatial positions
  LINALG::Matrix<1, NUMDOF_DISCSH3> x(true);
  x = SpatialConfiguration(disp);
  LINALG::TMatrix<FAD, 1, 9> crossprod2(true);
  for (int i = 0; i < NUMDOF_DISCSH3; i++)
  {
    x_FAD[i] = x(i);
    x_FAD[i].diff(i, NUMDOF_DISCSH3);
  }

  // Calculate analytical expressions for FAD force
  LINALG::TMatrix<FAD, 1, 9> ana_force(true);
  LINALG::Matrix<1, 9> tol(true);

  for (int j = 0; j < NUMNOD_DISCSH3; j++)
  {
    // For calculation of cross-product
    LINALG::TMatrix<FAD, 1, 3> side1(true);
    LINALG::TMatrix<FAD, 1, 3> side2(true);
    // Auxilary vector
    LINALG::TMatrix<FAD, 1, 3> side3(true);

    if (j == 0)  // 1st node
    {
      for (int i = 0; i < NODDOF_DISCSH3; i++)
      {
        side1(i) = x_FAD[3 + i] - x_FAD[i];
        side2(i) = x_FAD[6 + i] - x_FAD[3 + i];
        side3(i) = -x_FAD[6 + i] + x_FAD[3 + i];
      }
    }
    else if (j == 1)  // 2nd node
    {
      for (int i = 0; i < NODDOF_DISCSH3; i++)
      {
        side1(i) = x_FAD[6 + i] - x_FAD[3 + i];
        side2(i) = x_FAD[i] - x_FAD[6 + i];
        side3(i) = -x_FAD[i] + x_FAD[6 + i];
      }
    }
    else  // 3rd node
    {
      for (int i = 0; i < NODDOF_DISCSH3; i++)
      {
        side1(i) = x_FAD[i] - x_FAD[6 + i];
        side2(i) = x_FAD[3 + i] - x_FAD[i];
        side3(i) = -x_FAD[3 + i] + x_FAD[i];
      }
    }
    LINALG::TMatrix<FAD, 1, 3> crossprod1(true);

    // Cross Product side1xside2
    crossprod1 = CalcCrossProduct(side1, side2);

    LINALG::TMatrix<FAD, 1, 3> crossprod_aux(true);

    // Cross Product side3xcrossprod1
    crossprod_aux = CalcCrossProduct(side3, crossprod1);

    for (int i = 0; i < 3; i++) crossprod2(3 * j + i) = crossprod_aux(i);
  }

  // Calculate surface area at reference config FAD
  FAD area_ref = CalcSurfaceArea(disp, true);

  // Calculate surface area at spatial config FAD
  FAD area_curr = CalcSurfaceArea(disp, false);

  // Calculate the analytical expression of force
  for (int i = 0; i < NUMDOF_DISCSH3; i++)
  {
    ana_force(i) = -area_penalty * BehaviorFunctArea * crossprod2(i) /
                   (4 * std::pow(area_ref, 0.5) * area_curr);
    //     tol(i)= ana_force(i).val()-FAD_force(i).val();
  }

  // FAD force vector for computation of forces
  for (int i = 0; i < NUMDOF_DISCSH3; i++)
  {
    (*force)(i) += ana_force(i).val();
    //      if(CalcDampingForces)
    //      {
    //        (*force)(i)+=FAD_force_viscous(i).val();
    //      }
  }
  LINALG::Matrix<9, 9> stiff_aux(true);
  for (int i = 0; i < NUMDOF_DISCSH3; i++)
    for (int j = 0; j < NUMDOF_DISCSH3; j++)
    {
      (*stiffmatrix)(i, j) += ana_force(i).dx(j);
    }

  return;
}  // DRT::ELEMENTS::DiscSh3::AreaConstrtStiffmass

/*----------------------------------------------------------------------*
 |  evaluate the element (private)                             maf 04/07|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::DiscSh3::AreaConstrtQuadStiffmass(Teuchos::ParameterList& params,
    const std::vector<double>& disp, const std::vector<double>& vel,
    Epetra_SerialDenseMatrix* stiffmatrix, Epetra_SerialDenseVector* force)
{
  FAD BehaviorFunctArea;

  //  Teuchos::ParameterList StatMechParams =
  //  DRT::Problem::Instance()->StatisticalMechanicsParams(); double
  //  area_penalty=StatMechParams.get<double>("AREA_PENALTY",0.0); if (area_penalty==0)
  //    dserror("Please enter a non-zero AREA_PENALTY");
  //  // element geometry update

  double area_penalty = 0.0;

  const int numnode = NumNode();
  LINALG::SerialDenseMatrix x(numnode, 3);
  LINALG::SerialDenseMatrix x0(numnode, 3);

  SpatialConfiguration(x, disp);
  MaterialConfiguration(x0);

  // set up matrices and parameters needed for the evaluation of current
  // interfacial area and its derivatives w.r.t. the displacements

  int ndof = 3 * numnode;  // overall number of surface dofs
  double A = 0;            // interfacial area
  double A0 = 0;           // interfacial area Ref
  // first partial derivatives
  Teuchos::RCP<Epetra_SerialDenseVector> Adiff = Teuchos::rcp(new Epetra_SerialDenseVector);
  // second partial derivatives
  Teuchos::RCP<Epetra_SerialDenseMatrix> Adiff2 = Teuchos::rcp(new Epetra_SerialDenseMatrix);

  ComputeAreaRef(x0, numnode, ndof, A0);
  ComputeAreaDeriv(x, numnode, ndof, A, Adiff, Adiff2);


  LINALG::Matrix<1, 9> dummy_force(true);
  LINALG::Matrix<9, 9> dummy_stiff(true);
  for (int i = 0; i < ndof; ++i)
  {
    (*force)(i) += -area_penalty * (1 - A / A0) * (*Adiff)[i];
    for (int j = 0; j < ndof; ++j)
    {
      (*stiffmatrix)(i, j) += (area_penalty * (*Adiff)[i] * (*Adiff)[j] / A0 -
                               area_penalty * (1 - A / A0) * (*Adiff2)(i, j));
    }
  }

  return;
}  // DRT::ELEMENTS::DiscSh3::AreaConstrtQuadStiffmass


/*----------------------------------------------------------------------*
 |  evaluate the element (private)                       mukherjee 09/15|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::DiscSh3::VolConstrtStiffmass(Teuchos::ParameterList& params,
    const std::vector<double>& disp, Epetra_SerialDenseMatrix* stiffmatrix,
    Epetra_SerialDenseVector* force)
{
  dserror("stop");
  FAD BehaviorFunctVol;


  //  Teuchos::ParameterList StatMechParams =
  //  DRT::Problem::Instance()->StatisticalMechanicsParams(); FAD
  //  vol_penalty=StatMechParams.get<double>("VOL_PENALTY",0.0); if (vol_penalty==0)
  //    dserror("Please enter a non-zero VOL_PENALTY");

  FAD vol_penalty = 0.0;

  CalcBehaviorFunctVol(params, disp, BehaviorFunctVol);

  LINALG::TMatrix<FAD, 1, 9> FAD_force(true);
  LINALG::TMatrix<FAD, 1, 9> FAD_force_new(true);

  // Calculate FAD force
  for (int i = 0; i < NUMDOF_DISCSH3; i++)
  {
    FAD_force(i) = vol_penalty * BehaviorFunctVol * BehaviorFunctVol.dx(i);
  }

  std::vector<FAD> x_FAD(NUMDOF_DISCSH3, 0.0);
  // Get Spatial positions
  LINALG::Matrix<1, NUMDOF_DISCSH3> x(true);
  x = SpatialConfiguration(disp);
  LINALG::TMatrix<FAD, 1, 9> crossprod(true);
  LINALG::TMatrix<FAD, 1, 9> sign(true);
  for (int i = 0; i < NUMDOF_DISCSH3; i++)
  {
    x_FAD[i] = x(i);
    x_FAD[i].diff(i, NUMDOF_DISCSH3);
  }


  // Calculate analytical expressions for FAD force
  LINALG::TMatrix<FAD, 1, 9> ana_force(true);
  LINALG::Matrix<1, 9> tol(true);

  for (int j = 0; j < NUMNOD_DISCSH3; j++)
  {
    // vectors for vertices of pyramid
    LINALG::TMatrix<FAD, 1, 3> vertex1(true);
    LINALG::TMatrix<FAD, 1, 3> vertex2(true);
    LINALG::TMatrix<FAD, 1, 3> vertex3(true);

    if (j == 0)  // 1st node
    {
      for (int i = 0; i < NODDOF_DISCSH3; i++)
      {
        vertex1(i) = x_FAD[3 + i];
        vertex2(i) = x_FAD[6 + i];
        vertex3(i) = x_FAD[i];
      }
    }
    else if (j == 1)  // 2nd node
    {
      for (int i = 0; i < NODDOF_DISCSH3; i++)
      {
        vertex1(i) = x_FAD[6 + i];
        vertex2(i) = x_FAD[i];
        vertex3(i) = x_FAD[3 + i];
      }
    }
    else  // 3rd node
    {
      for (int i = 0; i < NODDOF_DISCSH3; i++)
      {
        vertex1(i) = x_FAD[i];
        vertex2(i) = x_FAD[3 + i];
        vertex3(i) = x_FAD[6 + i];
      }
    }
    LINALG::TMatrix<FAD, 1, 3> crossprod_aux(true);

    // Cross Product side3xcrossprod1
    crossprod_aux = CalcCrossProduct(vertex1, vertex2);

    for (int i = 0; i < 3; i++) crossprod(3 * j + i) = crossprod_aux(i);

    FAD dotprod = 0;
    // Dot Product
    for (int i = 0; i < NUMNOD_DISCSH3; i++)
    {
      dotprod += crossprod(3 * j + i) * vertex3(i);
    }
    if (dotprod != 0)
    {
      sign(3 * j + 0) = dotprod / abs(dotprod);
      sign(3 * j + 1) = dotprod / abs(dotprod);
      sign(3 * j + 2) = dotprod / abs(dotprod);
    }
  }

  FAD vol_ref;
  // Calculate surface area at reference config FAD
  CalcVolume(vol_ref, disp, true);

  //  INPAR::STATMECH::StatOutput StatOut =
  //  params.get<INPAR::STATMECH::StatOutput>("SPECIALOUTPUT",INPAR::STATMECH::statout_none); if
  //  (StatOut == INPAR::STATMECH::statout_vesceqshapes)
  //  {
  //    double time = params.get<double>("total time",0.0);
  //
  //
  //    if(time>=1.0)
  //      dserror("Stopping simulation! Volume can't be zero");
  //
  //    vol_ref=(1-time)*vol_ref;
  //  }


  FAD vol_curr;
  // Calculate surface area at spatial config FAD
  CalcVolume(vol_curr, disp, false);

  // vectors for vertices of pyramid
  LINALG::TMatrix<FAD, 1, 3> v1(true);
  LINALG::TMatrix<FAD, 1, 3> v2(true);
  LINALG::TMatrix<FAD, 1, 3> v3(true);
  for (int i = 0; i < NODDOF_DISCSH3; i++)
  {
    v1(i) = x_FAD[i];
    v2(i) = x_FAD[3 + i];
    v3(i) = x_FAD[6 + i];
  }

  // Calculate the analytical expression of force
  for (int i = 0; i < NUMDOF_DISCSH3; i++)
  {
    ana_force(i) = vol_penalty * 1 / std::pow(abs(vol_ref), 0.5) * BehaviorFunctVol * crossprod(i) *
                   sign(i) / 6;
    //    tol(i)=ana_force(i).val()-FAD_force(i).val();
  }

  // FAD force vector for computation of forces
  for (int i = 0; i < NUMDOF_DISCSH3; i++)
  {
    (*force)(i) += FAD_force(i).val();
  }

  for (int i = 0; i < NUMDOF_DISCSH3; i++)
    for (int j = 0; j < NUMDOF_DISCSH3; j++) (*stiffmatrix)(i, j) += ana_force(i).dx(j);

  return;
}  // DRT::ELEMENTS::DiscSh3::VolConstrtStiffmass

/*----------------------------------------------------------------------*
 |  evaluate the element (private)                       mukherjee 10/15|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::DiscSh3::VolConstrtGlobalStiff(Teuchos::ParameterList& params,
    DRT::Discretization& discretization, const std::vector<double>& disp,
    Epetra_SerialDenseMatrix* stiffmatrix, Epetra_SerialDenseVector* force)
{
  //  Teuchos::ParameterList StatMechParams =
  //  DRT::Problem::Instance()->StatisticalMechanicsParams(); FAD
  //  vol_penalty=StatMechParams.get<double>("VOL_PENALTY",0.0); if (vol_penalty==0)
  //    dserror("Please enter a non-zero VOL_PENALTY");

  FAD vol_penalty = 0.0;

  double vol_ref = params.get<double>("reference volume", 0.0);
  double vol_curr = params.get<double>("current volume", 0.0);

  //  INPAR::STATMECH::StatOutput StatOut =
  //  params.get<INPAR::STATMECH::StatOutput>("SPECIALOUTPUT",INPAR::STATMECH::statout_none); if
  //  (StatOut == INPAR::STATMECH::statout_vesceqshapes)
  //  {
  //    double time = params.get<double>("total time",0.0);
  //
  //    if(time>=1.0)
  //      dserror("Stopping simulation! Volume can't be zero");
  //
  //    vol_ref=(1-time)*vol_ref;
  //  }

  // Calculate behavior function
  FAD BehaviorFunctVol = std::pow(vol_ref, 0.5) * (1 - vol_curr / vol_ref);

  LINALG::TMatrix<FAD, 1, 9> FAD_force(true);

  // Calculate FAD force
  /* Important:  Use only for single processor
   * for verification purposes    */
  FAD FAD_vol_curr = 0;
  CalcVolume(FAD_vol_curr, disp, false);
  FAD BehaviorFunctVolFAD = std::pow(vol_ref, 0.5) * (1 - FAD_vol_curr / vol_ref);

  for (int i = 0; i < NUMDOF_DISCSH3; i++)
  {
    FAD_force(i) = vol_penalty * BehaviorFunctVol * BehaviorFunctVolFAD.dx(i);
  }

  std::vector<FAD> x_FAD(NUMDOF_DISCSH3, 0.0);
  // Get Spatial positions
  LINALG::Matrix<1, NUMDOF_DISCSH3> x(true);
  x = SpatialConfiguration(disp);
  LINALG::TMatrix<FAD, 1, 9> crossprod(true);
  for (int i = 0; i < NUMDOF_DISCSH3; i++)
  {
    x_FAD[i] = x(i);
    x_FAD[i].diff(i, NUMDOF_DISCSH3);
  }


  // Calculate analytical expressions for FAD force
  LINALG::TMatrix<FAD, 1, NUMDOF_DISCSH3> ana_force(true);
  LINALG::Matrix<1, NUMDOF_DISCSH3> tol(true);

  for (int j = 0; j < NUMNOD_DISCSH3; j++)
  {
    // vectors for vertices of pyramid
    LINALG::TMatrix<FAD, 1, 3> vertex1(true);
    LINALG::TMatrix<FAD, 1, 3> vertex2(true);

    if (j == 0)  // 1st node
    {
      for (int i = 0; i < NODDOF_DISCSH3; i++)
      {
        vertex1(i) = x_FAD[3 + i];
        vertex2(i) = x_FAD[6 + i];
      }
    }
    else if (j == 1)  // 2nd node
    {
      for (int i = 0; i < NODDOF_DISCSH3; i++)
      {
        vertex1(i) = x_FAD[6 + i];
        vertex2(i) = x_FAD[i];
      }
    }
    else  // 3rd node
    {
      for (int i = 0; i < NODDOF_DISCSH3; i++)
      {
        vertex1(i) = x_FAD[i];
        vertex2(i) = x_FAD[3 + i];
      }
    }
    LINALG::TMatrix<FAD, 1, 3> crossprod_aux(true);


    // Cross Product side3xcrossprod1
    crossprod_aux = CalcCrossProduct(vertex1, vertex2);

    for (int i = 0; i < 3; i++) crossprod(3 * j + i) = crossprod_aux(i);
  }


  LINALG::TMatrix<FAD, 1, 9> grad_vol(true);
  // Calculate the analytical expression of force
  for (int i = 0; i < NUMDOF_DISCSH3; i++)
  {
    grad_vol(i) = crossprod(i) / 6;
    ana_force(i) = -vol_penalty * (1 - vol_curr / vol_ref) * grad_vol(i);
    tol(i) = ana_force(i).val() - FAD_force(i).val();
  }

  // FAD force vector for computation of forces
  for (int i = 0; i < NUMDOF_DISCSH3; i++)
  {
    (*force)(i) += ana_force(i).val();
  }

  LINALG::TMatrix<FAD, 9, 9> grad2_vol(true);
  LINALG::TMatrix<FAD, 9, 9> stiff_vol(true);
  LINALG::Matrix<9, 9> stiff_aux(true);
  // Corrected linearization
  for (int i = 0; i < NUMDOF_DISCSH3; i++)
    for (int j = 0; j < NUMDOF_DISCSH3; j++)
    {
      grad2_vol(i, j) = grad_vol(i).dx(j);
      stiff_vol(i, j) = vol_penalty * grad_vol(i) * grad_vol(j) / vol_ref -
                        vol_penalty * (1 - vol_curr / vol_ref) * grad2_vol(i, j);
      (*stiffmatrix)(i, j) += stiff_vol(i, j).val();
      stiff_aux(i, j) = stiff_vol(i, j).val();
    }

  return;
}  // DRT::ELEMENTS::DiscSh3::VolConstrtStiffmass

/*----------------------------------------------------------------------*
 |  evaluate the element (private)                             maf 04/07|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::DiscSh3::AreaConstrtGlobalStiff(Teuchos::ParameterList& params,
    DRT::Discretization& discretization, const std::vector<double>& disp,
    Epetra_SerialDenseMatrix* stiffmatrix, Epetra_SerialDenseVector* force)
{
  //  Teuchos::ParameterList StatMechParams =
  //  DRT::Problem::Instance()->StatisticalMechanicsParams(); FAD
  //  area_penalty=StatMechParams.get<double>("AREA_PENALTY",0.0);
  //
  //  if (area_penalty==0)
  //    dserror("Please enter a non-zero AREA_PENALTY");

  FAD area_penalty = 0.0;

  double area_ref = params.get<double>("reference area", 0.0);
  double area_curr = params.get<double>("current area", 0.0);

  // Calculate behavior function
  FAD BehaviorFunctArea = std::pow(area_ref, 0.5) * (1 - area_curr / area_ref);

  LINALG::TMatrix<FAD, 1, 9> FAD_force(true);

  // Calculate FAD force
  /* Important:  Use only for single processor
   * for verification purposes    */
  FAD FAD_area_curr = CalcSurfaceArea(disp, false);
  FAD BehaviorFunctAreaFAD = std::pow(area_ref, 0.5) * (1 - FAD_area_curr / area_ref);

  for (int i = 0; i < NUMDOF_DISCSH3; i++)
  {
    FAD_force(i) = area_penalty * BehaviorFunctArea * BehaviorFunctAreaFAD.dx(i);
  }

  std::vector<FAD> x_FAD(NUMDOF_DISCSH3, 0.0);
  // Get Spatial positions
  LINALG::Matrix<1, NUMDOF_DISCSH3> x(true);
  x = SpatialConfiguration(disp);
  LINALG::TMatrix<FAD, 1, 9> crossprod(true);
  for (int i = 0; i < NUMDOF_DISCSH3; i++)
  {
    x_FAD[i] = x(i);
    x_FAD[i].diff(i, NUMDOF_DISCSH3);
  }


  // Calculate analytical expressions for FAD force
  LINALG::TMatrix<FAD, 1, 9> ana_force(true);
  LINALG::Matrix<1, 9> tol(true);

  // Calculate analytical expressions for FAD force
  LINALG::TMatrix<FAD, 1, 9> crossprod2(true);
  for (int j = 0; j < NUMNOD_DISCSH3; j++)
  {
    // For calculation of cross-product
    LINALG::TMatrix<FAD, 1, 3> side1(true);
    LINALG::TMatrix<FAD, 1, 3> side2(true);
    // Auxilary vector
    LINALG::TMatrix<FAD, 1, 3> side3(true);

    if (j == 0)  // 1st node
    {
      for (int i = 0; i < NODDOF_DISCSH3; i++)
      {
        side1(i) = x_FAD[3 + i] - x_FAD[i];
        side2(i) = x_FAD[6 + i] - x_FAD[3 + i];
        side3(i) = -x_FAD[6 + i] + x_FAD[3 + i];
      }
    }
    else if (j == 1)  // 2nd node
    {
      for (int i = 0; i < NODDOF_DISCSH3; i++)
      {
        side1(i) = x_FAD[6 + i] - x_FAD[3 + i];
        side2(i) = x_FAD[i] - x_FAD[6 + i];
        side3(i) = -x_FAD[i] + x_FAD[6 + i];
      }
    }
    else  // 3rd node
    {
      for (int i = 0; i < NODDOF_DISCSH3; i++)
      {
        side1(i) = x_FAD[i] - x_FAD[6 + i];
        side2(i) = x_FAD[3 + i] - x_FAD[i];
        side3(i) = -x_FAD[3 + i] + x_FAD[i];
      }
    }
    LINALG::TMatrix<FAD, 1, 3> crossprod1(true);

    // Cross Product side1xside2
    crossprod1 = CalcCrossProduct(side1, side2);

    LINALG::TMatrix<FAD, 1, 3> crossprod_aux(true);

    // Cross Product side3xcrossprod1
    crossprod_aux = CalcCrossProduct(side3, crossprod1);

    for (int i = 0; i < 3; i++) crossprod2(3 * j + i) = crossprod_aux(i);
  }


  LINALG::TMatrix<FAD, 1, 9> grad_area(true);
  // Calculate the analytical expression of force
  for (int i = 0; i < NUMDOF_DISCSH3; i++)
  {
    grad_area(i) = crossprod2(i) / (4 * FAD_area_curr);
    ana_force(i) = -area_penalty * (1 - area_curr / area_ref) * grad_area(i);
    tol(i) = ana_force(i).val() - FAD_force(i).val();
  }


  // FAD force vector for computation of forces
  for (int i = 0; i < NUMDOF_DISCSH3; i++)
  {
    (*force)(i) += ana_force(i).val();
  }

  LINALG::TMatrix<FAD, 9, 9> grad2_area(true);
  LINALG::TMatrix<FAD, 9, 9> stiff_area(true);

  // Corrected linearization
  for (int i = 0; i < NUMDOF_DISCSH3; i++)
    for (int j = 0; j < NUMDOF_DISCSH3; j++)
    {
      grad2_area(i, j) = grad_area(i).dx(j);
      //      grad2_area(i,j)= crossprod2(i).dx(j)/(4*area_curr)-
      //      crossprod2(i)*crossprod2(j)/(16*std::pow(area_curr,3));
      stiff_area(i, j) = area_penalty * grad_area(i) * grad_area(j) / area_ref -
                         area_penalty * (1 - area_curr / area_ref) * grad2_area(i, j);
      (*stiffmatrix)(i, j) += stiff_area(i, j).val();
    }


  return;
}  // DRT::ELEMENTS::DiscSh3::VolConstrtStiffmass

/*--------------------------------------------------------------------*
 |Calculate energy functional arising from the area constraint        |
 |of the element (private)                             mukherjee 07/15|
 *--------------------------------------------------------------------*/
void DRT::ELEMENTS::DiscSh3::CalcBehaviorFunctArea(
    Teuchos::ParameterList& params, const std::vector<double>& disp, FAD& BehaviorFunctArea)
{
  //  Teuchos::ParameterList StatMechParams =
  //  DRT::Problem::Instance()->StatisticalMechanicsParams(); FAD
  //  area_penalty=StatMechParams.get<double>("AREA_PENALTY",0.0); if (area_penalty==0)
  //    dserror("Please enter a non-zero AREA_PENALTY");

  FAD area_penalty = 0.0;

  // Calculate surface area at reference config FAD
  FAD area_ref = CalcSurfaceArea(disp, true);

  // Calculate surface area at spatial config FAD
  FAD area_curr = CalcSurfaceArea(disp, false);
  // Normalised
  BehaviorFunctArea = std::pow(area_ref, 0.5) * (1 - (area_curr) / (area_ref));

  return;
}  // DRT::ELEMENTS::DiscSh3::CalcEnergyAreaConst


/*----------------------------------------------------------------------*
 |Calculate energy functional arising from the vol constraint           |
 |of the element (private)                               mukherjee 07/15|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::DiscSh3::CalcBehaviorFunctVol(
    Teuchos::ParameterList& params, const std::vector<double>& disp, FAD& BehaviorFunctVol)
{
  // Volume at reference configuration
  FAD vol_ref;

  // Calculate volume at reference config FAD
  CalcVolume(vol_ref, disp, true);

  //  INPAR::STATMECH::StatOutput StatOut =
  //  params.get<INPAR::STATMECH::StatOutput>("SPECIALOUTPUT",INPAR::STATMECH::statout_none); if
  //  (StatOut == INPAR::STATMECH::statout_vesceqshapes)
  //  {
  //    double time = params.get<double>("total time",0.0);
  //    if(time>=1.0)
  //      dserror("Stopping simulation! Volume can't be zero");
  //    vol_ref=(1-time)*vol_ref;
  //  }

  // Volume at spatial configuration
  FAD vol_curr;

  // Calculate voluem at spatial config
  CalcVolume(vol_curr, disp, false);

  if (vol_curr.val() == 0.0) dserror("Volume can't be zero!");
  // Calculate energy function arising from volume constraint
  BehaviorFunctVol = std::pow(abs(vol_curr), 0.5) * (1 - vol_ref / vol_curr);

  return;
}  // DRT::ELEMENTS::DiscSh3::CalcEnergyVolConst


/*------------------------------------------------------------------------*
 |  Calculate surface area of the element (private)        mukherjee 04/15|
 *------------------------------------------------------------------------*/
FAD DRT::ELEMENTS::DiscSh3::CalcSurfaceArea(const std::vector<double>& disp, bool refconfig)
{
  FAD area;
  std::vector<FAD> x_FAD(NUMDOF_DISCSH3, 0.0);

  // 3 nodes, 3 dimensions
  LINALG::Matrix<1, NUMDOF_DISCSH3> x(true);
  if (refconfig)  // reference config
  {
    x = MaterialConfiguration();
    for (int i = 0; i < NUMDOF_DISCSH3; i++)
    {
      x_FAD[i] = x(i);
    }
  }
  else
  {
    x = SpatialConfiguration(disp);
    for (int i = 0; i < NUMDOF_DISCSH3; i++)
    {
      x_FAD[i] = x(i);
      x_FAD[i].diff(i, NUMDOF_DISCSH3);
    }
  }

  LINALG::TMatrix<FAD, 1, 3> side1(true);
  LINALG::TMatrix<FAD, 1, 3> side2(true);
  for (int j = 0; j < 3; j++)
  {
    side1(j) = x_FAD[j + 3] - x_FAD[j];
    side2(j) = x_FAD[j + 6] - x_FAD[j + 3];
  }

  LINALG::TMatrix<FAD, 1, 3> crossprod(true);

  // Cross Product
  crossprod = CalcCrossProduct(side1, side2);

  FAD norm_crossprod = 0.0;
  for (int i = 0; i < 3; i++)
  {
    norm_crossprod += pow(crossprod(i), 2);
  }
  norm_crossprod = pow(norm_crossprod, 0.5);

  area = 0.5 * norm_crossprod;  // Always positive

  return area;
}

/*------------------------------------------------------------------------*
 |  Calculate surface area of the element (private)        mukherjee 04/15|
 *------------------------------------------------------------------------*/
FAD DRT::ELEMENTS::DiscSh3::CalcSurfaceAreaFAD(const std::vector<double>& disp, bool refconfig)
{
  FAD area;
  LINALG::TMatrix<FAD, 1, 3> area_vector(true);
  LINALG::TMatrix<FAD, 1, NUMDOF_DISCSH3> FAD_force_vector(true);
  std::vector<FAD> x_FAD(NUMDOF_DISCSH3, 0.0);

  // 3 nodes, 3 dimensions
  LINALG::Matrix<1, NUMDOF_DISCSH3> x(true);
  if (refconfig)  // reference config
  {
    x = MaterialConfiguration();
    for (int i = 0; i < NUMDOF_DISCSH3; i++)
    {
      x_FAD[i] = x(i);
    }
  }
  else
  {
    x = SpatialConfiguration(disp);
    for (int i = 0; i < NUMDOF_DISCSH3; i++)
    {
      x_FAD[i] = x(i);
      x_FAD[i].diff(i, NUMDOF_DISCSH3);
    }
  }

  LINALG::TMatrix<FAD, 1, 3> side1(true);
  LINALG::TMatrix<FAD, 1, 3> side2(true);
  for (int j = 0; j < NUMNOD_DISCSH3; j++)
  {
    // For calculation of cross-product
    LINALG::TMatrix<FAD, 1, 3> side1(true);
    LINALG::TMatrix<FAD, 1, 3> side2(true);
    // Auxilary vector
    LINALG::TMatrix<FAD, 1, 3> side3(true);

    if (j == 0)  // 1st node
    {
      for (int i = 0; i < NODDOF_DISCSH3; i++)
      {
        side1(i) = x_FAD[3 + i] - x_FAD[i];
        side2(i) = x_FAD[6 + i] - x_FAD[3 + i];
        side3(i) = -x_FAD[6 + i] + x_FAD[3 + i];
      }
    }
    else if (j == 1)  // 2nd node
    {
      for (int i = 0; i < NODDOF_DISCSH3; i++)
      {
        side1(i) = x_FAD[6 + i] - x_FAD[3 + i];
        side2(i) = x_FAD[i] - x_FAD[6 + i];
        side3(i) = -x_FAD[i] + x_FAD[6 + i];
      }
    }
    else  // 3rd node
    {
      for (int i = 0; i < NODDOF_DISCSH3; i++)
      {
        side1(i) = x_FAD[i] - x_FAD[6 + i];
        side2(i) = x_FAD[3 + i] - x_FAD[i];
        side3(i) = -x_FAD[3 + i] + x_FAD[i];
      }
    }
    LINALG::TMatrix<FAD, 1, 3> crossprod(true);

    // Cross Product side1xside2
    crossprod = CalcCrossProduct(side1, side2);

    FAD norm_crossprod = 0.0;
    for (int i = 0; i < 3; i++)
    {
      norm_crossprod += pow(crossprod(i), 2);
    }
    norm_crossprod = pow(norm_crossprod, 0.5);

    area = 0.5 * norm_crossprod;  // Always positive

    area_vector(j) = area;

    for (int k = 0; k < 3; k++) FAD_force_vector(k + 3 * j) = area.dx(k + 3 * j);
  }

  return area;
}

/*------------------------------------------------------------------------*
 |  Calculate surface area of the element (private)        mukherjee 04/07|
 *------------------------------------------------------------------------*/
FAD DRT::ELEMENTS::DiscSh3::CalcSurfaceArea(DRT::Discretization& discretization, bool refconfig)
{
  FAD area;
  std::vector<FAD> x_FAD(NUMDOF_DISCSH3, 0.0);

  // 3 nodes, 3 dimensions
  LINALG::Matrix<1, NUMDOF_DISCSH3> x(true);
  if (refconfig)  // reference config
  {
    x = MaterialConfiguration();
    for (int i = 0; i < NUMDOF_DISCSH3; i++)
    {
      x_FAD[i] = x(i);
    }
  }
  else
  {
    x = SpatialConfiguration(discretization);
    for (int i = 0; i < NUMDOF_DISCSH3; i++)
    {
      x_FAD[i] = x(i);
      x_FAD[i].diff(i, NUMDOF_DISCSH3);
    }
  }

  LINALG::TMatrix<FAD, 1, 3> side1(true);
  LINALG::TMatrix<FAD, 1, 3> side2(true);
  for (int j = 0; j < 3; j++)
  {
    side1(j) = x_FAD[j + 3] - x_FAD[j];
    side2(j) = x_FAD[j + 6] - x_FAD[j + 3];
  }

  LINALG::TMatrix<FAD, 1, 3> crossprod(true);

  // Cross Product
  crossprod = CalcCrossProduct(side1, side2);

  FAD norm_crossprod = 0.0;
  for (int i = 0; i < 3; i++)
  {
    norm_crossprod += pow(crossprod(i), 2);
  }
  norm_crossprod = pow(norm_crossprod, 0.5);

  area = 0.5 * norm_crossprod;  // Always positive

  return area;
}

/*------------------------------------------------------------------------*
 |  Calculate surface area of the element (public)         mukherjee 04/15|
 *------------------------------------------------------------------------*/
double DRT::ELEMENTS::DiscSh3::CalcSurfArea(
    DRT::Discretization& dis, const Epetra_Vector& discol, bool refconfig)
{
  FAD area;
  std::vector<FAD> x_FAD(NUMDOF_DISCSH3, 0.0);

  // 3 nodes, 3 dimensions
  LINALG::Matrix<1, NUMDOF_DISCSH3> x(true);
  if (refconfig)  // reference config
  {
    x = MaterialConfiguration();
    for (int i = 0; i < NUMDOF_DISCSH3; i++) x_FAD[i] = x(i);
  }
  else
  {
    x = SpatialConfiguration(dis, discol);
    for (int i = 0; i < NUMDOF_DISCSH3; i++)
    {
      x_FAD[i] = x(i);
      x_FAD[i].diff(i, NUMDOF_DISCSH3);
    }
  }

  LINALG::TMatrix<FAD, 1, 3> side1(true);
  LINALG::TMatrix<FAD, 1, 3> side2(true);
  for (int j = 0; j < 3; j++)
  {
    side1(j) = x_FAD[j + 3] - x_FAD[j];
    side2(j) = x_FAD[j + 6] - x_FAD[j + 3];
  }

  LINALG::TMatrix<FAD, 1, 3> crossprod(true);

  // Cross Product
  crossprod = CalcCrossProduct(side1, side2);

  FAD norm_crossprod = 0.0;
  for (int i = 0; i < 3; i++)
  {
    norm_crossprod += pow(crossprod(i), 2);
  }
  norm_crossprod = pow(norm_crossprod, 0.5);

  area = 0.5 * norm_crossprod;  // Always positive

  return area.val();
}


/*------------------------------------------------------------------------*
 |  Calculate surface area of a triangle (private)        mukherjee 04/15|
 *------------------------------------------------------------------------*/
FAD DRT::ELEMENTS::DiscSh3::CalcSurfaceArea(LINALG::TMatrix<FAD, 1, 3>& vertex1,
    LINALG::TMatrix<FAD, 1, 3>& vertex2, LINALG::TMatrix<FAD, 1, 3>& vertex3)
{
  FAD area;

  LINALG::TMatrix<FAD, 1, 3> side1(true);
  LINALG::TMatrix<FAD, 1, 3> side2(true);
  for (int j = 0; j < 3; j++)
  {
    side1(j) = vertex1(j) - vertex2(j);
    side2(j) = vertex3(j) - vertex2(j);
  }

  LINALG::TMatrix<FAD, 1, 3> crossprod(true);

  // Cross Product
  crossprod = CalcCrossProduct(side1, side2);

  FAD norm_crossprod = 0.0;
  for (int i = 0; i < 3; i++)
  {
    norm_crossprod += pow(crossprod(i), 2);
  }
  norm_crossprod = pow(norm_crossprod, 0.5);

  area = 0.5 * norm_crossprod;  // Always positive

  return area;
}
/*------------------------------------------------------------------------*
 |  Calculate surface area of a triangle (private)        mukherjee 04/15|
 *------------------------------------------------------------------------*/
double DRT::ELEMENTS::DiscSh3::CalcSurfaceArea(
    LINALG::Matrix<1, 3>& vertex1, LINALG::Matrix<1, 3>& vertex2, LINALG::Matrix<1, 3>& vertex3)
{
  FAD area;

  LINALG::TMatrix<FAD, 1, 3> side1(true);
  LINALG::TMatrix<FAD, 1, 3> side2(true);
  for (int j = 0; j < 3; j++)
  {
    side1(j) = vertex1(j) - vertex2(j);
    side2(j) = vertex3(j) - vertex2(j);
  }

  LINALG::TMatrix<FAD, 1, 3> crossprod(true);

  // Cross Product
  crossprod = CalcCrossProduct(side1, side2);

  FAD norm_crossprod = 0.0;
  for (int i = 0; i < 3; i++)
  {
    norm_crossprod += pow(crossprod(i), 2);
  }
  norm_crossprod = pow(norm_crossprod, 0.5);

  area = 0.5 * norm_crossprod;  // Always positive

  return area.val();
}

/*----------------------------------------------------------------------------*
 |Calculate vol encompassed by surface of the el. (private)    mukherjee 04/07|
 *----------------------------------------------------------------------------*/
void DRT::ELEMENTS::DiscSh3::CalcVolume(
    FAD& volume, const std::vector<double>& disp, bool refconfig)
{
  std::vector<FAD> x_FAD(NUMDOF_DISCSH3, 0.0);

  // 3 nodes, 3 dimensions
  LINALG::Matrix<1, NUMDOF_DISCSH3> x(true);
  if (refconfig)
  {
    x = MaterialConfiguration();
    for (int i = 0; i < NUMDOF_DISCSH3; i++)
    {
      x_FAD[i] = x(i);
    }
  }
  else
  {
    x = SpatialConfiguration(disp);
    for (int i = 0; i < NUMDOF_DISCSH3; i++)
    {
      x_FAD[i] = x(i);
      x_FAD[i].diff(i, NUMDOF_DISCSH3);
    }
  }
  /* calculations according to Nystroem2002 */
  // Calculate unit normal
  LINALG::TMatrix<FAD, 1, 3> normal = CalcSurfaceNormal(x_FAD);

  // vectors for vertices of pyramid
  LINALG::TMatrix<FAD, 1, 3> v1(true);
  LINALG::TMatrix<FAD, 1, 3> v2(true);
  LINALG::TMatrix<FAD, 1, 3> v3(true);

  for (int j = 0; j < 3; j++)
  {
    v1(j) = x_FAD[j];
    v2(j) = x_FAD[j + 3];
    v3(j) = x_FAD[j + 6];
  }

  /**********Old Way ***********

//  FAD area = CalcSurfaceArea(vertex1, vertex2, vertex3);

// Projection of triangular area along x axis onto y-z plane

//  FAD area_projected= 0.5* normal(0);
//
  LINALG::TMatrix  <FAD,1,3> crossprod(true);
//
//
//  //Cross Product
  crossprod= CalcCrossProduct(vertex1, vertex2);

//  FAD norm_crossprod= 0.0;
  for (int i=0; i<3; i++)
  {
    volume += crossprod(i)*vertex3(i)/6;
//    volume += normal(i)*(vertex1(i)+vertex2(i)+vertex3(i))/3;
  }

  ***********End Old Way*************/
  //  volume = 1/6*-(v3(0));

  volume = (-(v3(0) * v2(1) * v1(2)) + (v2(0) * v3(1) * v1(2)) + (v3(0) * v1(1) * v2(2)) -
               (v1(0) * v3(1) * v2(2)) - (v2(0) * v1(1) * v3(2)) + (v1(0) * v2(1) * v3(2))) /
           6;

  //  std::cout<<"v1="<<v1<<std::endl;
  //  std::cout<<"v2="<<v2<<std::endl;
  //  std::cout<<"v3="<<v3<<std::endl;
  //
  //  std::cout<<"volume="<<volume<<std::endl;
  //  dserror("stop");
  //  volume= abs(volume);

  return;
}


/*----------------------------------------------------------------------------*
 |Calculate vol encompassed by surface of the el. (public)     mukherjee 06/15|
 *----------------------------------------------------------------------------*/
double DRT::ELEMENTS::DiscSh3::CalcVolume(
    DRT::Discretization& dis, const Epetra_Vector& discol, bool refconfig)
{
  FAD volume = 0;
  std::vector<FAD> x_FAD(NUMDOF_DISCSH3, 0.0);

  // 3 nodes, 3 dimensions
  LINALG::Matrix<1, NUMDOF_DISCSH3> x(true);
  if (refconfig)
  {
    x = MaterialConfiguration();
    for (int i = 0; i < NUMDOF_DISCSH3; i++)
    {
      x_FAD[i] = x(i);
    }
  }
  else
  {
    x = SpatialConfiguration(dis, discol);
    //    x=SpatialConfiguration(discretization);
    for (int i = 0; i < NUMDOF_DISCSH3; i++)
    {
      x_FAD[i] = x(i);
      x_FAD[i].diff(i, NUMDOF_DISCSH3);
    }
  }

  /* calculations according to Nystroem2002 */
  // Calculate unit normal
  LINALG::TMatrix<FAD, 1, 3> normal = CalcSurfaceNormal(x_FAD);

  // vectors for vertices of pyramid
  LINALG::TMatrix<FAD, 1, 3> v1(true);
  LINALG::TMatrix<FAD, 1, 3> v2(true);
  LINALG::TMatrix<FAD, 1, 3> v3(true);

  for (int j = 0; j < 3; j++)
  {
    v1(j) = x_FAD[j];
    v2(j) = x_FAD[j + 3];
    v3(j) = x_FAD[j + 6];
  }

  /**********Old Way ***********

//  FAD area = CalcSurfaceArea(vertex1, vertex2, vertex3);

// Projection of triangular area along x axis onto y-z plane

//  FAD area_projected= 0.5* normal(0);
//
  LINALG::TMatrix  <FAD,1,3> crossprod(true);
//
//
//  //Cross Product
  crossprod= CalcCrossProduct(vertex1, vertex2);

//  FAD norm_crossprod= 0.0;
  for (int i=0; i<3; i++)
  {
    volume += crossprod(i)*vertex3(i)/6;
//    volume += normal(i)*(vertex1(i)+vertex2(i)+vertex3(i))/3;
  }
//  volume= abs(volume);

  ***********End Old Way*************/

  volume = (-(v3(0) * v2(1) * v1(2)) + (v2(0) * v3(1) * v1(2)) + (v3(0) * v1(1) * v2(2)) -
               (v1(0) * v3(1) * v2(2)) - (v2(0) * v1(1) * v3(2)) + (v1(0) * v2(1) * v3(2))) /
           6;

  return volume.val();
}


/*----------------------------------------------------------------------*
 |  Calculate surface normal of the element (private)    mukherjee 04/07|
 *----------------------------------------------------------------------*/
LINALG::TMatrix<FAD, 1, 3> DRT::ELEMENTS::DiscSh3::CalcSurfaceNormal(std::vector<FAD>& x_FAD)
{
  LINALG::TMatrix<FAD, 1, 3> normal(true);
  LINALG::Matrix<1, 3> normal_val(true);

  LINALG::TMatrix<FAD, 1, 3> side1(true);
  LINALG::TMatrix<FAD, 1, 3> side2(true);
  for (int j = 0; j < 3; j++)
  {
    side1(j) = x_FAD[j + 3] - x_FAD[j];
    side2(j) = x_FAD[j + 6] - x_FAD[j];
  }

  LINALG::TMatrix<FAD, 1, 3> crossprod(true);

  // Cross Product
  crossprod = CalcCrossProduct(side1, side2);

  FAD norm_crossprod = pow(crossprod.Dot(crossprod), 0.5);

  for (int i = 0; i < 3; i++)
  {
    normal(i) = crossprod(i) / norm_crossprod;
    normal_val(i) = normal(i).val();
  }


  return normal;
}



/*----------------------------------------------------------------------*
 |  lump mass matrix (private)                           mukherjee 07/15|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::DiscSh3::sh3_lumpmass(
    const std::vector<double>& disp, Epetra_SerialDenseMatrix* massmatrix)
{
  // lump mass matrix (In this case, dummy. There is no massterm. )
  if (massmatrix != NULL)
  {
    // we assume #elemat2 is a square matrix
    for (int c = 0; c < (*massmatrix).N(); c++)  // parse columns
    {
      for (int r = 0; r < (*massmatrix).M(); r++)  // parse rows
      {
        if (r == c)
        {
          (*massmatrix)(r, c) = 1;
        }
        else
          (*massmatrix)(r, c) = 0.0;
      }
    }
  }
}

/*----------------------------------------------------------------------*
 |  init the element (public)                            mukherjee 04/15|
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::DiscSh3Type::Initialize(DRT::Discretization& dis)
{
  for (int i = 0; i < dis.NumMyColElements(); ++i)
  {
    if (dis.lColElement(i)->ElementType() != *this) continue;
    DRT::ELEMENTS::DiscSh3* actele = dynamic_cast<DRT::ELEMENTS::DiscSh3*>(dis.lColElement(i));
    if (!actele) dserror("cast to DiscSh3* failed");
    LINALG::Matrix<1, 9> X = actele->MaterialConfiguration();
    for (int i = 0; i < 9; i++)
    {
      actele->x_n_1_(i) = X(i);
    }

    actele->CalcDampingForces_ = false;
  }

  return 0;
}
