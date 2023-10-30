/*----------------------------------------------------------------------*/
/*! \file
\brief averaged nodal volume tet4
\level 3
*----------------------------------------------------------------------*/
#include "baci_contact_analytical.H"
#include "baci_discretization_fem_general_utils_fem_shapefunctions.H"
#include "baci_discretization_fem_general_utils_integration.H"
#include "baci_lib_discret.H"
#include "baci_lib_globalproblem.H"
#include "baci_lib_utils.H"
#include "baci_linalg_serialdensematrix.H"
#include "baci_linalg_serialdensevector.H"
#include "baci_linalg_utils_sparse_algebra_math.H"
#include "baci_mat_so3_material.H"
#include "baci_so3_tet4av.H"
#include "baci_utils_exceptions.H"
#include "baci_utils_function.H"

#include <Teuchos_SerialDenseSolver.hpp>
#include <Teuchos_StandardParameterEntryValidators.hpp>



/*----------------------------------------------------------------------*
 |  evaluate the element (public)                           seitz 03/16 |
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::So_tet4av::Evaluate(Teuchos::ParameterList& params,
    DRT::Discretization& discretization, std::vector<int>& lm,
    CORE::LINALG::SerialDenseMatrix& elemat1_epetra,
    CORE::LINALG::SerialDenseMatrix& elemat2_epetra,
    CORE::LINALG::SerialDenseVector& elevec1_epetra,
    CORE::LINALG::SerialDenseVector& elevec2_epetra,
    CORE::LINALG::SerialDenseVector& elevec3_epetra)
{
  // Check whether the solid material PostSetup() routine has already been called and call it if not
  EnsureMaterialPostSetup(params);

  CORE::LINALG::Matrix<NUMDOF_SOTET4av, NUMDOF_SOTET4av> elemat1(elemat1_epetra.values(), true);
  CORE::LINALG::Matrix<NUMDOF_SOTET4av, NUMDOF_SOTET4av> elemat2(elemat2_epetra.values(), true);
  CORE::LINALG::Matrix<NUMDOF_SOTET4av, 1> elevec1(elevec1_epetra.values(), true);
  CORE::LINALG::Matrix<NUMDOF_SOTET4av, 1> elevec2(elevec2_epetra.values(), true);

  // start with "none"
  DRT::ELEMENTS::So_tet4av::ActionType act = So_tet4av::none;

  // get the required action
  std::string action = params.get<std::string>("action", "none");
  if (action == "none")
    dserror("No action supplied");
  else if (action == "calc_struct_nlnstiff")
    act = So_tet4av::calc_struct_nlnstiff;
  else if (action == "calc_struct_internalforce")
    act = So_tet4av::calc_struct_internalforce;
  else if (action == "calc_struct_nlnstiffmass")
    act = So_tet4av::calc_struct_nlnstiffmass;
  else if (action == "calc_struct_stress")
    act = So_tet4av::calc_struct_stress;
  else if (action == "calc_struct_update_istep")
    act = So_tet4av::calc_struct_update_istep;
  else if (action == "calc_struct_reset_istep")
    act = So_tet4av::calc_struct_reset_istep;
  else if (action == "calc_struct_reset_all")
    act = So_tet4av::calc_struct_reset_all;
  else if (action == "calc_struct_recover")
    return 0;
  else if (action == "calc_struct_predict")
    return 0;
  else
    dserror("Unknown type of action for So_tet4av");


  // what should the element do
  switch (act)
  {
    //==================================================================================
    // nonlinear stiffness and internal force vector
    case calc_struct_nlnstiff:
    {
      // need current displacement and residual forces
      Teuchos::RCP<const Epetra_Vector> disp = discretization.GetState("displacement");
      if (disp == Teuchos::null) dserror("Cannot get state vectors 'displacement' and/or residual");
      std::vector<double> mydisp(lm.size());
      DRT::UTILS::ExtractMyValues(*disp, mydisp, lm);

      nlnstiffmass(lm, mydisp, &elemat1, nullptr, &elevec1, nullptr, nullptr, params,
          INPAR::STR::stress_none, INPAR::STR::strain_none);
    }
    break;

    //==================================================================================
    // internal force vector only
    case calc_struct_internalforce:
    {
      // need current displacement and residual forces
      Teuchos::RCP<const Epetra_Vector> disp = discretization.GetState("displacement");
      if (disp == Teuchos::null) dserror("Cannot get state vectors 'displacement' and/or residual");
      std::vector<double> mydisp(lm.size());
      DRT::UTILS::ExtractMyValues(*disp, mydisp, lm);
      nlnstiffmass(lm, mydisp, nullptr, nullptr, &elevec1, nullptr, nullptr, params,
          INPAR::STR::stress_none, INPAR::STR::strain_none);
    }
    break;

    //==================================================================================
    // nonlinear stiffness, internal force vector, and consistent mass matrix
    case calc_struct_nlnstiffmass:
    {
      // need current displacement and residual forces
      Teuchos::RCP<const Epetra_Vector> disp = discretization.GetState("displacement");
      if (disp == Teuchos::null) dserror("Cannot get state vectors 'displacement' and/or residual");
      std::vector<double> mydisp(lm.size());
      DRT::UTILS::ExtractMyValues(*disp, mydisp, lm);

      nlnstiffmass(lm, mydisp, &elemat1, &elemat2, &elevec1, nullptr, nullptr, params,
          INPAR::STR::stress_none, INPAR::STR::strain_none);
    }
    break;

    //==================================================================================
    // evaluate stresses and strains at gauss points
    case calc_struct_stress:
    {
      // nothing to do for ghost elements
      if (discretization.Comm().MyPID() == Owner())
      {
        Teuchos::RCP<const Epetra_Vector> disp = discretization.GetState("displacement");
        Teuchos::RCP<std::vector<char>> stressdata =
            params.get<Teuchos::RCP<std::vector<char>>>("stress", Teuchos::null);
        Teuchos::RCP<std::vector<char>> straindata =
            params.get<Teuchos::RCP<std::vector<char>>>("strain", Teuchos::null);
        if (disp == Teuchos::null) dserror("Cannot get state vectors 'displacement'");
        if (stressdata == Teuchos::null) dserror("Cannot get 'stress' data");
        if (straindata == Teuchos::null) dserror("Cannot get 'strain' data");
        std::vector<double> mydisp(lm.size());
        DRT::UTILS::ExtractMyValues(*disp, mydisp, lm);
        CORE::LINALG::Matrix<NUMGPT_SOTET4av, MAT::NUM_STRESS_3D> stress(true);  // set to zero
        CORE::LINALG::Matrix<NUMGPT_SOTET4av, MAT::NUM_STRESS_3D> strain(true);
        auto iostress =
            DRT::INPUT::get<INPAR::STR::StressType>(params, "iostress", INPAR::STR::stress_none);
        auto iostrain =
            DRT::INPUT::get<INPAR::STR::StrainType>(params, "iostrain", INPAR::STR::strain_none);

        nlnstiffmass(
            lm, mydisp, nullptr, nullptr, nullptr, &stress, &strain, params, iostress, iostrain);

        {
          DRT::PackBuffer data;
          AddtoPack(data, stress);
          data.StartPacking();
          AddtoPack(data, stress);
          std::copy(data().begin(), data().end(), std::back_inserter(*stressdata));
        }
        {
          DRT::PackBuffer data;
          AddtoPack(data, strain);
          data.StartPacking();
          AddtoPack(data, strain);
          std::copy(data().begin(), data().end(), std::back_inserter(*straindata));
        }
      }
    }
    break;

    //==================================================================================
    case calc_struct_update_istep:
    {
      // Update of history for materials
      SolidMaterial()->Update();
    }
    break;

    //==================================================================================
    case calc_struct_reset_istep:
    {
      // Reset of history (if needed)
      SolidMaterial()->ResetStep();
    }
    break;

    default:
      dserror("Unknown type of action for so_tet4");
      break;
  }

  return 0;
}


/*----------------------------------------------------------------------*
 |  Integrate a Volume Neumann boundary condition (public)  seitz 03/16 |
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::So_tet4av::EvaluateNeumann(Teuchos::ParameterList& params,
    DRT::Discretization& discretization, DRT::Condition& condition, std::vector<int>& lm,
    CORE::LINALG::SerialDenseVector& elevec1, CORE::LINALG::SerialDenseMatrix* elemat1)
{
  // get values and switches from the condition
  const auto* onoff = condition.Get<std::vector<int>>("onoff");
  const auto* val = condition.Get<std::vector<double>>("val");

  /*
  **    TIME CURVE BUSINESS
  */
  // find out whether we will use a time curve
  const double time = params.get("total time", -1.0);

  // ensure that at least as many curves/functs as dofs are available
  if (int(onoff->size()) < NUMDIM_SOTET4av)
    dserror("Fewer functions or curves defined than the element has dofs.");

  for (int checkdof = NUMDIM_SOTET4av; checkdof < int(onoff->size()); ++checkdof)
  {
    if ((*onoff)[checkdof] != 0)
      dserror("Number of Dimensions in Neumann_Evalutaion is 3. Further DoFs are not considered.");
  }

  // (SPATIAL) FUNCTION BUSINESS
  static_assert(NUMGPT_SOTET4av == 1);
  const auto* funct = condition.Get<std::vector<int>>("funct");
  CORE::LINALG::Matrix<NUMDIM_SOTET4av, 1> xrefegp(false);
  bool havefunct = false;
  if (funct)
    for (int dim = 0; dim < NUMDIM_SOTET4av; dim++)
      if ((*funct)[dim] > 0)
      {
        havefunct = true;
        break;
      }
  /* ============================================================================*/

  // update element geometry
  CORE::LINALG::Matrix<NUMNOD_SOTET4av, NUMDIM_SOTET4av> xrefe;
  DRT::Node** nodes = Nodes();
  for (int i = 0; i < NUMNOD_SOTET4av; ++i)
  {
    const double* x = nodes[i]->X();
    xrefe(i, 0) = x[0];
    xrefe(i, 1) = x[1];
    xrefe(i, 2) = x[2];
  }

  CORE::LINALG::Matrix<NUMDIM_SOTET4av, NUMDIM_SOTET4av> jac;
  CORE::LINALG::Matrix<NUMNOD_SOTET4av, 1> shapefct;
  CORE::LINALG::Matrix<NUMDIM_SOTET4av, NUMNOD_SOTET4av> deriv;

  /* ================================================= Loop over Gauss Points */
  for (int gp = 0; gp < NUMGPT_SOTET4av; gp++)
  {
    CORE::DRT::UTILS::shape_function<CORE::FE::CellType::tet4>(xsi_[gp], shapefct);
    CORE::DRT::UTILS::shape_function_deriv1<CORE::FE::CellType::tet4>(xsi_[gp], deriv);
    jac.Multiply(deriv, xrefe);

    // material/reference co-ordinates of Gauss point
    if (havefunct)
    {
      for (int dim = 0; dim < NUMDIM_SOTET4av; dim++)
      {
        xrefegp(dim) = 0.0;
        for (int nodid = 0; nodid < NUMNOD_SOTET4av; ++nodid)
          xrefegp(dim) += shapefct(nodid) * xrefe(nodid, dim);
      }
    }

    // integration factor
    double fac = wgt_[gp] * detJ_[gp];
    // distribute/add over element load vector
    for (int dim = 0; dim < NUMDIM_SOTET4av; dim++)
    {
      if ((*onoff)[dim])
      {
        // function evaluation
        const int functnum = (funct) ? (*funct)[dim] : -1;
        const double functfac =
            (functnum > 0) ? DRT::Problem::Instance()
                                 ->FunctionById<CORE::UTILS::FunctionOfSpaceTime>(functnum - 1)
                                 .Evaluate(xrefegp.A(), time, dim)
                           : 1.0;
        const double dim_fac = (*val)[dim] * fac * functfac;
        for (int nodid = 0; nodid < NUMNOD_SOTET4av; ++nodid)
        {
          elevec1[nodid * NUMDIM_SOTET4av + NODDOF_SOTET4av] += shapefct(nodid) * dim_fac;
        }
      }
    }


  } /* ==================================================== end of Loop over GP */

  return 0;
}  // DRT::ELEMENTS::So_tet4av::EvaluateNeumann


/*----------------------------------------------------------------------*
 |  init the element jacobian mapping and integration       seitz 03/16 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_tet4av::InitJacobianMapping()
{
  CORE::LINALG::Matrix<NUMNOD_SOTET4av, NUMDIM_SOTET4av> xrefe;
  DRT::Node** nodes = Nodes();
  for (int i = 0; i < NUMNOD_SOTET4av; ++i)
  {
    const double* x = nodes[i]->X();
    xrefe(i, 0) = x[0];
    xrefe(i, 1) = x[1];
    xrefe(i, 2) = x[2];
  }
  CORE::LINALG::Matrix<NUMDIM_SOTET4av, NUMNOD_SOTET4av> deriv;

  CORE::DRT::UTILS::IntPointsAndWeights<3> intpoints(CORE::DRT::UTILS::GaussRule3D::tet_1point);
  numgpt_ = intpoints.IP().nquad;
  xsi_.resize(numgpt_);
  wgt_.resize(numgpt_);
  invJ_.resize(numgpt_);
  detJ_.resize(numgpt_);
  for (int gp = 0; gp < numgpt_; ++gp)
  {
    wgt_[gp] = (intpoints.IP().qwgt)[gp];
    const double* gpcoord = (intpoints.IP().qxg)[gp];
    for (int idim = 0; idim < 3; idim++) xsi_[gp](idim) = gpcoord[idim];

    CORE::DRT::UTILS::shape_function_deriv1<CORE::FE::CellType::tet4>(xsi_[gp], deriv);

    invJ_[gp].Multiply(deriv, xrefe);

    // xij_ = ds/dx
    detJ_[gp] = invJ_[gp].Invert();
    if (detJ_[gp] < 1.0E-16) dserror("ZERO OR NEGATIVE JACOBIAN DETERMINANT: %f", detJ_[gp]);
  }

  return;
}


/*----------------------------------------------------------------------*
 |  evaluate the element (private)                          seitz 03/16 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_tet4av::nlnstiffmass(std::vector<int>& lm,  // location matrix
    std::vector<double>& disp,                                     // current displacements
    CORE::LINALG::Matrix<NUMDOF_SOTET4av, NUMDOF_SOTET4av>*
        stiffmatrix,                                                     // element stiffness matrix
    CORE::LINALG::Matrix<NUMDOF_SOTET4av, NUMDOF_SOTET4av>* massmatrix,  // element mass matrix
    CORE::LINALG::Matrix<NUMDOF_SOTET4av, 1>* force,  // element internal force vector
    CORE::LINALG::Matrix<NUMGPT_SOTET4av, MAT::NUM_STRESS_3D>* elestress,  // stresses at GP
    CORE::LINALG::Matrix<NUMGPT_SOTET4av, MAT::NUM_STRESS_3D>* elestrain,  // strains at GP
    Teuchos::ParameterList& params,         // algorithmic parameters e.g. time
    const INPAR::STR::StressType iostress,  // stress output option
    const INPAR::STR::StrainType iostrain   // strain output option
)
{
  // current  displacements of element
  CORE::LINALG::Matrix<NUMNOD_SOTET4av, NUMDIM_SOTET4av> xrefe;
  CORE::LINALG::Matrix<NUMNOD_SOTET4av, NUMDIM_SOTET4av> xcurr;
  DRT::Node** nodes = Nodes();
  CORE::LINALG::Matrix<NUMNOD_SOTET4av, 1> nodalVol;
  for (int i = 0; i < NUMNOD_SOTET4av; ++i)
  {
    const double* x = nodes[i]->X();
    xrefe(i, 0) = x[0];
    xrefe(i, 1) = x[1];
    xrefe(i, 2) = x[2];

    xcurr(i, 0) = x[0] + disp[i * NODDOF_SOTET4av + 0];
    xcurr(i, 1) = x[1] + disp[i * NODDOF_SOTET4av + 1];
    xcurr(i, 2) = x[2] + disp[i * NODDOF_SOTET4av + 2];
    nodalVol(i) = 1. + disp[i * NODDOF_SOTET4av + 3];
  }

  CORE::LINALG::Matrix<NUMNOD_SOTET4av, 1> shapefct;
  CORE::LINALG::Matrix<NUMDIM_SOTET4av, NUMNOD_SOTET4av> deriv;
  CORE::LINALG::Matrix<NUMDIM_SOTET4av, NUMNOD_SOTET4av> N_XYZ;
  CORE::LINALG::Matrix<NUMDIM_SOTET4av, NUMDIM_SOTET4av> defgrd;
  CORE::LINALG::Matrix<NUMDIM_SOTET4av, NUMDIM_SOTET4av> defgrd_bar;
  CORE::LINALG::Matrix<NUMDIM_SOTET4av, NUMDIM_SOTET4av> rcg_bar;
  CORE::LINALG::Matrix<MAT::NUM_STRESS_3D, 1> gl_bar;
  CORE::LINALG::Matrix<MAT::NUM_STRESS_3D, NUMDOF_SOTET4av> bop;

  /* =========================================================================*/
  /* ============================================== Loop over Gauss Points ===*/
  /* =========================================================================*/
  for (int gp = 0; gp < numgpt_; gp++)
  {
    // shape functions (shapefunct) and their first derivatives (deriv)
    CORE::DRT::UTILS::shape_function<CORE::FE::CellType::tet4>(xsi_[gp], shapefct);
    CORE::DRT::UTILS::shape_function_deriv1<CORE::FE::CellType::tet4>(xsi_[gp], deriv);


    /* get the inverse of the Jacobian matrix which looks like:
    **            [ x_,r  y_,r  z_,r ]^-1
    **     J^-1 = [ x_,s  y_,s  z_,s ]
    **            [ x_,t  y_,t  z_,t ]
    */
    // compute derivatives N_XYZ at gp w.r.t. material coordinates
    // by N_XYZ = J^-1 * N_rst
    N_XYZ.Multiply(invJ_[gp], deriv);
    const double detJ = detJ_[gp];
    defgrd.MultiplyTT(xcurr, N_XYZ);
    CORE::LINALG::Matrix<NUMDIM_SOTET4av, NUMDIM_SOTET4av> invdefgrd;
    const double detF = invdefgrd.Invert(defgrd);
    const double intNodalVol = shapefct.Dot(nodalVol);
    if (intNodalVol < 0.) dserror("intNodalVol < 0");
    if (detF < 0.) dserror("detF < 0");
    const double fbar_fac = std::pow(intNodalVol / detF, 1. / 3.);
    defgrd_bar.Update(fbar_fac, defgrd, 0.);
    rcg_bar.MultiplyTN(defgrd_bar, defgrd_bar);
    for (int i = 0; i < NUMDIM_SOTET4av; ++i) gl_bar(i) = .5 * (rcg_bar(i, i) - 1.);
    gl_bar(3) = rcg_bar(0, 1);
    gl_bar(4) = rcg_bar(1, 2);
    gl_bar(5) = rcg_bar(0, 2);

    CORE::LINALG::Matrix<6, 1> pk2;
    CORE::LINALG::Matrix<6, 6> cmat;
    SolidMaterial()->Evaluate(&defgrd_bar, &gl_bar, params, &pk2, &cmat, gp, Id());

    // return gp stresses
    switch (iostress)
    {
      case INPAR::STR::stress_2pk:
      {
        if (elestress == nullptr) dserror("stress data not available");
        for (int i = 0; i < MAT::NUM_STRESS_3D; ++i) (*elestress)(gp, i) = pk2(i);
      }
      break;
      case INPAR::STR::stress_cauchy:
      {
        if (elestress == nullptr) dserror("stress data not available");
        const double detF_bar = defgrd_bar.Determinant();

        CORE::LINALG::Matrix<3, 3> pkstress_bar;
        pkstress_bar(0, 0) = pk2(0);
        pkstress_bar(0, 1) = pk2(3);
        pkstress_bar(0, 2) = pk2(5);
        pkstress_bar(1, 0) = pkstress_bar(0, 1);
        pkstress_bar(1, 1) = pk2(1);
        pkstress_bar(1, 2) = pk2(4);
        pkstress_bar(2, 0) = pkstress_bar(0, 2);
        pkstress_bar(2, 1) = pkstress_bar(1, 2);
        pkstress_bar(2, 2) = pk2(2);

        CORE::LINALG::Matrix<3, 3> temp;
        CORE::LINALG::Matrix<3, 3> cauchystress_bar;
        temp.Multiply(1.0 / detF_bar, defgrd_bar, pkstress_bar);
        cauchystress_bar.MultiplyNT(temp, defgrd_bar);

        (*elestress)(gp, 0) = cauchystress_bar(0, 0);
        (*elestress)(gp, 1) = cauchystress_bar(1, 1);
        (*elestress)(gp, 2) = cauchystress_bar(2, 2);
        (*elestress)(gp, 3) = cauchystress_bar(0, 1);
        (*elestress)(gp, 4) = cauchystress_bar(1, 2);
        (*elestress)(gp, 5) = cauchystress_bar(0, 2);
      }
      break;
      case INPAR::STR::stress_none:
        break;
      default:
        dserror("requested stress type not available");
        break;
    }

    for (int i = 0; i < NUMNOD_SOTET4av; ++i)
    {
      for (int a = 0; a < NUMDIM_SOTET4av; ++a)
        for (int b = 0; b < NUMDIM_SOTET4av; ++b)
          bop(a, NODDOF_SOTET4av * i + b) = defgrd(b, a) * N_XYZ(a, i);
      /* ~~~ */
      bop(3, NODDOF_SOTET4av * i + 0) = defgrd(0, 0) * N_XYZ(1, i) + defgrd(0, 1) * N_XYZ(0, i);
      bop(3, NODDOF_SOTET4av * i + 1) = defgrd(1, 0) * N_XYZ(1, i) + defgrd(1, 1) * N_XYZ(0, i);
      bop(3, NODDOF_SOTET4av * i + 2) = defgrd(2, 0) * N_XYZ(1, i) + defgrd(2, 1) * N_XYZ(0, i);
      bop(4, NODDOF_SOTET4av * i + 0) = defgrd(0, 1) * N_XYZ(2, i) + defgrd(0, 2) * N_XYZ(1, i);
      bop(4, NODDOF_SOTET4av * i + 1) = defgrd(1, 1) * N_XYZ(2, i) + defgrd(1, 2) * N_XYZ(1, i);
      bop(4, NODDOF_SOTET4av * i + 2) = defgrd(2, 1) * N_XYZ(2, i) + defgrd(2, 2) * N_XYZ(1, i);
      bop(5, NODDOF_SOTET4av * i + 0) = defgrd(0, 2) * N_XYZ(0, i) + defgrd(0, 0) * N_XYZ(2, i);
      bop(5, NODDOF_SOTET4av * i + 1) = defgrd(1, 2) * N_XYZ(0, i) + defgrd(1, 0) * N_XYZ(2, i);
      bop(5, NODDOF_SOTET4av * i + 2) = defgrd(2, 2) * N_XYZ(0, i) + defgrd(2, 0) * N_XYZ(2, i);
    }

    const double detJ_w = detJ * wgt_[gp];

    if (force != nullptr)
    {
      force->MultiplyTN(detJ_w / fbar_fac, bop, pk2, 1.0);

      if (gp == 0)
        for (int i = 0; i < NUMNOD_SOTET4av; ++i)
          (*force)(i * NODDOF_SOTET4av + 3) += nodalVol(i) - detF;
    }

    if (stiffmatrix != nullptr)
    {
      // integrate `elastic' and `initial-displacement' stiffness matrix
      // keu = keu + (B^T . C . B) * detJ * w(gp)
      CORE::LINALG::Matrix<MAT::NUM_STRESS_3D, NUMDOF_SOTET4av> cb;
      cb.Multiply(cmat, bop);
      stiffmatrix->MultiplyTN(detJ_w * fbar_fac, bop, cb, 1.0);

      // integrate `geometric' stiffness matrix and add to keu *****************
      CORE::LINALG::Matrix<6, 1> sfac(pk2);  // auxiliary integrated stress
      sfac.Scale(detJ_w / fbar_fac);         // detJ*w(gp)*[S11,S22,S33,S12=S21,S23=S32,S13=S31]
      std::vector<double> SmB_L(3);          // intermediate Sm.B_L
      // kgeo += (B_L^T . sigma . B_L) * detJ * w(gp)  with B_L = Ni,Xj see NiliFEM-Skript
      for (int inod = 0; inod < NUMNOD_SOTET4av; ++inod)
      {
        SmB_L[0] = sfac(0) * N_XYZ(0, inod) + sfac(3) * N_XYZ(1, inod) + sfac(5) * N_XYZ(2, inod);
        SmB_L[1] = sfac(3) * N_XYZ(0, inod) + sfac(1) * N_XYZ(1, inod) + sfac(4) * N_XYZ(2, inod);
        SmB_L[2] = sfac(5) * N_XYZ(0, inod) + sfac(4) * N_XYZ(1, inod) + sfac(2) * N_XYZ(2, inod);
        for (int jnod = 0; jnod < NUMNOD_SOTET4av; ++jnod)
        {
          double bopstrbop = 0.0;  // intermediate value
          for (int idim = 0; idim < NUMDIM_SOTET4av; ++idim)
            bopstrbop += N_XYZ(idim, jnod) * SmB_L[idim];
          (*stiffmatrix)(4 * inod + 0, 4 * jnod + 0) += bopstrbop;
          (*stiffmatrix)(4 * inod + 1, 4 * jnod + 1) += bopstrbop;
          (*stiffmatrix)(4 * inod + 2, 4 * jnod + 2) += bopstrbop;
        }
      }  // end of integrate `geometric' stiffness******************************

      // integrate additional fbar matrix
      CORE::LINALG::Matrix<NUMDIM_SOTET4av, NUMDIM_SOTET4av> cauchygreen;
      cauchygreen.MultiplyTN(defgrd, defgrd);
      CORE::LINALG::Matrix<MAT::NUM_STRESS_3D, 1> cauchygreenvector;
      cauchygreenvector(0) = cauchygreen(0, 0);
      cauchygreenvector(1) = cauchygreen(1, 1);
      cauchygreenvector(2) = cauchygreen(2, 2);
      cauchygreenvector(3) = 2 * cauchygreen(0, 1);
      cauchygreenvector(4) = 2 * cauchygreen(1, 2);
      cauchygreenvector(5) = 2 * cauchygreen(2, 0);

      CORE::LINALG::Matrix<4 * 4, 1> htensor;
      for (int n = 0; n < 4 * 3; n++)
        for (int i = 0; i < 3; i++)
          htensor(n + n / 3) -= invdefgrd(i, n % 3) * N_XYZ(i, n / 3) / detF * intNodalVol;
      for (int i = 0; i < 4; ++i) htensor(i * 4 + 3) += shapefct(i) / detF;

      CORE::LINALG::Matrix<4 * 4, 1> bops;
      bops.MultiplyTN(bop, pk2);
      stiffmatrix->MultiplyNT(-1. / 3. * std::pow(fbar_fac, -4.) * detJ_w, bops, htensor, 1.);

      CORE::LINALG::Matrix<6, 1> ccg;
      ccg.Multiply(cmat, cauchygreenvector);
      CORE::LINALG::Matrix<4 * 4, 1> bopccg;
      bopccg.MultiplyTN(bop, ccg);
      stiffmatrix->MultiplyNT(detJ_w * std::pow(fbar_fac, -2.) / 3., bopccg, htensor, 1.);

      if (gp == 0)
        for (int inod = 0; inod < 4; ++inod)
        {
          (*stiffmatrix)(inod * 4 + 3, inod * 4 + 3) += 1.;
          for (int n = 0; n < 4 * 3; n++)
            for (int i = 0; i < 3; i++)
              (*stiffmatrix)(inod * 4 + 3, n + n / 3) -=
                  detF * invdefgrd(i, n % 3) * N_XYZ(i, n / 3);
        }
    }
  }  // end gp loop

  if (massmatrix != nullptr)  // evaluate mass matrix +++++++++++++++++++++++++
  {
    double density = Material()->Density(0);  // density at the only Gauss point the material has!
    // integrate consistent mass matrix
    // jacobian is constatnt
    double ifactor, massfactor;
    // needs more than one gauss point
    CORE::DRT::UTILS::IntPointsAndWeights<3> intpoints(CORE::DRT::UTILS::GaussRule3D::tet_4point);
    CORE::LINALG::Matrix<3, 1> xsi;

    for (int gp = 0; gp < intpoints.IP().nquad; gp++)
    {
      for (int i = 0; i < NUMDIM_SOTET4av; ++i) xsi(i) = (intpoints.IP().qxg)[gp][i];
      CORE::DRT::UTILS::shape_function<CORE::FE::CellType::tet4>(xsi_[gp], shapefct);
      const double factor = detJ_[0] * density * (intpoints.IP().qwgt)[gp];
      for (int inod = 0; inod < 4; ++inod)
      {
        ifactor = shapefct(inod) * factor;
        for (int jnod = 0; jnod < 4; ++jnod)
        {
          massfactor = shapefct(jnod) * ifactor;  // intermediate factor
          (*massmatrix)(4 * inod + 0, 4 * jnod + 0) += massfactor;
          (*massmatrix)(4 * inod + 1, 4 * jnod + 1) += massfactor;
          (*massmatrix)(4 * inod + 2, 4 * jnod + 2) += massfactor;
        }
      }
    }
  }  // end of mass matrix +++++++++++++++++++++++++++++++++++++++++++++++++++

  return;
}  // DRT::ELEMENTS::So_tet4av::nlnstiffmass


/*----------------------------------------------------------------------*
 |  init the element (public)                               seitz 03/16 |
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::So_tet4avType::Initialize(DRT::Discretization& dis)
{
  for (int i = 0; i < dis.NumMyColElements(); ++i)
  {
    if (dis.lColElement(i)->ElementType() != *this) continue;
    auto* actele = dynamic_cast<DRT::ELEMENTS::So_tet4av*>(dis.lColElement(i));
    if (!actele) dserror("cast to So_tet4av* failed");
    actele->InitJacobianMapping();
  }
  return 0;
}
