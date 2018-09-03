/*----------------------------------------------------------------------*/
/*!
 \file wall1_poro_p1_evaluate.cpp

 \brief evaluate methods for 2D wall element for structure part of porous medium
        using p1 approach (mixed approach)

 \level 2

 \maintainer Andreas Rauch
             rauch@lnm.mw.tum.de
             http://www.lnm.mw.tum.de

 *----------------------------------------------------------------------*/

#include "wall1_poro_p1.H"
#include "wall1_poro_p1_eletypes.H"

#include "../drt_lib/drt_utils.H"
#include "../drt_lib/drt_discret.H"

#include "../drt_mat/fluidporo.H"
#include "../drt_mat/structporo.H"

#include "../drt_lib/drt_globalproblem.H"

#include "../drt_fem_general/drt_utils_gder2.H"

#include "../drt_structure_new/str_elements_paramsinterface.H"

/*----------------------------------------------------------------------*
 *                                                            vuong 07/13|
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::Wall1_PoroP1<distype>::ComputePorosityAndLinearization(
    Teuchos::ParameterList& params, const double& press, const double& J, const int& gp,
    const LINALG::Matrix<my::numnod_, 1>& shapfct, const LINALG::Matrix<my::numnod_, 1>* myporosity,
    const LINALG::Matrix<1, my::numdof_>& dJ_dus, double& porosity,
    LINALG::Matrix<1, my::numdof_>& dphi_dus)
{
  if (myporosity == NULL)
    dserror("no porosity values given!");
  else
    porosity = shapfct.Dot(*myporosity);

  dphi_dus.PutScalar(0.0);

  return;
}

/*----------------------------------------------------------------------*
 *                                                            vuong 07/13|
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::Wall1_PoroP1<distype>::ComputePorosityAndLinearizationOD(
    Teuchos::ParameterList& params, const double& press, const double& J, const int& gp,
    const LINALG::Matrix<my::numnod_, 1>& shapfct, const LINALG::Matrix<my::numnod_, 1>* myporosity,
    double& porosity, double& dphi_dp)
{
  dphi_dp = 0.0;

  if (myporosity == NULL)
    dserror("no porosity values given!");
  else
    porosity = shapfct.Dot(*myporosity);

  return;
}

/*----------------------------------------------------------------------*
 |  evaluate the element (public)                                       |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
int DRT::ELEMENTS::Wall1_PoroP1<distype>::Evaluate(Teuchos::ParameterList& params,
    DRT::Discretization& discretization, DRT::Element::LocationArray& la,
    Epetra_SerialDenseMatrix& elemat1_epetra, Epetra_SerialDenseMatrix& elemat2_epetra,
    Epetra_SerialDenseVector& elevec1_epetra, Epetra_SerialDenseVector& elevec2_epetra,
    Epetra_SerialDenseVector& elevec3_epetra)
{
  this->SetParamsInterfacePtr(params);
  ELEMENTS::ActionType act = ELEMENTS::none;

  if (this->IsParamsInterface())
  {
    act = this->ParamsInterface().GetActionType();
  }
  else
  {
    // get the required action
    std::string action = params.get<std::string>("action", "none");
    if (action == "none")
      dserror("No action supplied");
    else if (action == "struct_poro_calc_fluidcoupling")
      act = ELEMENTS::struct_poro_calc_fluidcoupling;
    else if (action == "calc_struct_energy")
      act = ELEMENTS::struct_calc_energy;
  }

  // what should the element do
  switch (act)
  {
    //==================================================================================
    // off diagonal terms in stiffness matrix for monolithic coupling
    case ELEMENTS::struct_poro_calc_fluidcoupling:
    {
      // in some cases we need to write/change some data before evaluating
      my::PreEvaluate(params, discretization, la);

      MyEvaluate(params, discretization, la, elemat1_epetra, elemat2_epetra, elevec1_epetra,
          elevec2_epetra, elevec3_epetra);
    }
    break;
    //==================================================================================
    case ELEMENTS::struct_calc_energy:
    {
      // in some cases we need to write/change some data before evaluating
      my::PreEvaluate(params, discretization, la);

      // evaluate parent solid element
      Wall1::Evaluate(params, discretization, la[0].lm_, elemat1_epetra, elemat2_epetra,
          elevec1_epetra, elevec2_epetra, elevec3_epetra);
    }
    break;
    //==================================================================================
    default:
    {
      // in some cases we need to write/change some data before evaluating
      my::PreEvaluate(params, discretization, la);

      Epetra_SerialDenseMatrix elemat1_sub;
      Epetra_SerialDenseMatrix elemat2_sub;
      Epetra_SerialDenseVector elevec1_sub;
      Epetra_SerialDenseVector elevec2_sub;
      Epetra_SerialDenseVector elevec3_sub;

      if (elemat1_epetra.A()) elemat1_sub.Shape(my::numdof_, my::numdof_);
      if (elemat2_epetra.A()) elemat2_sub.Shape(my::numdof_, my::numdof_);
      if (elevec1_epetra.A()) elevec1_sub.Resize(my::numdof_);
      if (elevec2_epetra.A()) elevec2_sub.Resize(my::numdof_);
      if (elevec3_epetra.A()) elevec3_sub.Resize(my::numdof_);

      std::vector<int> lm_sub;
      for (int i = 0; i < my::numnod_; i++)
        for (int j = 0; j < my::numdim_; j++) lm_sub.push_back(la[0].lm_[i * noddof_ + j]);

      // evaluate parent solid element
      Wall1::Evaluate(params, discretization, lm_sub, elemat1_sub, elemat2_sub, elevec1_sub,
          elevec2_sub, elevec3_sub);


      if (elemat1_epetra.A())
        for (int i = 0; i < my::numnod_; i++)
          for (int j = 0; j < my::numdim_; j++)
            for (int k = 0; k < my::numnod_; k++)
              for (int l = 0; l < my::numdim_; l++)
                elemat1_epetra(i * noddof_ + j, k * noddof_ + l) =
                    elemat1_sub(i * my::noddof_ + j, k * my::noddof_ + l);

      if (elemat2_epetra.A())
        for (int i = 0; i < my::numnod_; i++)
          for (int j = 0; j < my::numdim_; j++)
            for (int k = 0; k < my::numnod_; k++)
              for (int l = 0; l < my::numdim_; l++)
                elemat2_epetra(i * noddof_ + j, k * noddof_ + l) =
                    elemat2_sub(i * my::noddof_ + j, k * my::noddof_ + l);

      if (elevec1_epetra.A())
        for (int i = 0; i < my::numnod_; i++)
          for (int j = 0; j < my::numdim_; j++)
            elevec1_epetra(i * noddof_ + j) = elevec1_sub(i * my::noddof_ + j);

      if (elevec2_epetra.A())
        for (int i = 0; i < my::numnod_; i++)
          for (int j = 0; j < my::numdim_; j++)
            elevec2_epetra(i * noddof_ + j) = elevec2_sub(i * my::noddof_ + j);

      if (elevec3_epetra.A())
        for (int i = 0; i < my::numnod_; i++)
          for (int j = 0; j < my::numdim_; j++)
            elevec3_epetra(i * noddof_ + j) = elevec3_sub(i * my::noddof_ + j);

      // add volume coupling specific terms
      MyEvaluate(params, discretization, la, elemat1_epetra, elemat2_epetra, elevec1_epetra,
          elevec2_epetra, elevec3_epetra);
    }
    break;
  }  // action

  return 0;
}

/*----------------------------------------------------------------------*
 |  evaluate the element (protected)                                       |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
int DRT::ELEMENTS::Wall1_PoroP1<distype>::MyEvaluate(Teuchos::ParameterList& params,
    DRT::Discretization& discretization, DRT::Element::LocationArray& la,
    Epetra_SerialDenseMatrix& elemat1_epetra, Epetra_SerialDenseMatrix& elemat2_epetra,
    Epetra_SerialDenseVector& elevec1_epetra, Epetra_SerialDenseVector& elevec2_epetra,
    Epetra_SerialDenseVector& elevec3_epetra)
{
  this->SetParamsInterfacePtr(params);
  ELEMENTS::ActionType act = ELEMENTS::none;

  if (this->IsParamsInterface())
  {
    act = this->ParamsInterface().GetActionType();
  }
  else
  {
    // get the required action
    std::string action = params.get<std::string>("action", "none");
    if (action == "none")
      dserror("No action supplied");
    else if (action == "calc_struct_internalforce")
      act = ELEMENTS::struct_calc_internalforce;
    else if (action == "calc_struct_nlnstiff")
      act = ELEMENTS::struct_calc_nlnstiff;
    else if (action == "calc_struct_nlnstiffmass")
      act = ELEMENTS::struct_calc_nlnstiffmass;
    else if (action == "struct_poro_calc_fluidcoupling")
      act = ELEMENTS::struct_poro_calc_fluidcoupling;
  }

  // what should the element do
  switch (act)
  {
    //==================================================================================
    // nonlinear stiffness, damping and internal force vector for poroelasticity
    case ELEMENTS::struct_calc_nlnstiff:
    case ELEMENTS::struct_calc_nlnstiffmass:
    {
      if (la.Size() > 1)
      {
        // stiffness
        LINALG::Matrix<numdof_, numdof_> elemat1(elemat1_epetra.A(), true);
        // damping
        LINALG::Matrix<numdof_, numdof_> elemat2(elemat2_epetra.A(), true);
        // internal force vector
        LINALG::Matrix<numdof_, 1> elevec1(elevec1_epetra.A(), true);
        // LINALG::Matrix<numdof_,1> elevec2(elevec2_epetra.A(),true);
        // elemat2,elevec2+3 are not used anyway

        // build the location vector only for the structure field
        std::vector<int> lm = la[0].lm_;

        LINALG::Matrix<my::numdim_, my::numnod_> mydisp(true);
        LINALG::Matrix<my::numnod_, 1> myporosity(true);
        my::ExtractValuesFromGlobalVector(
            discretization, 0, la[0].lm_, &mydisp, &myporosity, "displacement");

        LINALG::Matrix<numdof_, numdof_>* matptr = NULL;
        if (elemat1.IsInitialized()) matptr = &elemat1;

        enum INPAR::STR::DampKind damping =
            params.get<enum INPAR::STR::DampKind>("damping", INPAR::STR::damp_none);
        LINALG::Matrix<numdof_, numdof_>* matptr2 = NULL;
        if (elemat2.IsInitialized() and (damping == INPAR::STR::damp_material)) matptr2 = &elemat2;

        // need current fluid state,
        // call the fluid discretization: fluid equates 2nd dofset
        // disassemble velocities and pressures

        LINALG::Matrix<my::numdim_, my::numnod_> myvel(true);

        LINALG::Matrix<my::numdim_, my::numnod_> myfluidvel(true);
        LINALG::Matrix<my::numnod_, 1> myepreaf(true);

        if (discretization.HasState(0, "velocity"))
          my::ExtractValuesFromGlobalVector(discretization, 0, la[0].lm_, &myvel, NULL, "velocity");

        if (discretization.HasState(1, "fluidvel"))
          // extract local values of the global vectors
          my::ExtractValuesFromGlobalVector(
              discretization, 1, la[1].lm_, &myfluidvel, &myepreaf, "fluidvel");

        // calculate tangent stiffness matrix
        nlnstiff_poroelast(lm, mydisp, myvel, &myporosity, myfluidvel, myepreaf, matptr, matptr2,
            &elevec1, params);
      }
    }
    break;

    //==================================================================================
    // coupling terms in force-vector and stiffness matrix for poroelasticity
    case ELEMENTS::struct_poro_calc_fluidcoupling:
    {
      // stiffness
      LINALG::Matrix<numdof_, (my::numdim_ + 1) * my::numnod_> elemat1(elemat1_epetra.A(), true);
      // LINALG::Matrix<numdof_,(numdim_+1)*numnod_> elemat2(elemat2_epetra.A(),true);

      // internal force vector
      // LINALG::Matrix<numdof_,1> elevec1(elevec1_epetra.A(),true);
      // LINALG::Matrix<numdof_,1> elevec2(elevec2_epetra.A(),true);

      // elemat2,elevec2+3 are not used anyway

      // build the location vector only for the structure field
      std::vector<int> lm = la[0].lm_;

      LINALG::Matrix<numdof_, (my::numdim_ + 1)* my::numnod_>* matptr = NULL;
      if (elemat1.IsInitialized()) matptr = &elemat1;

      // need current fluid state,
      // call the fluid discretization: fluid equates 2nd dofset
      // disassemble velocities and pressures
      if (discretization.HasState(1, "fluidvel"))
      {
        LINALG::Matrix<my::numdim_, my::numnod_> myvel(true);
        LINALG::Matrix<my::numdim_, my::numnod_> myfluidvel(true);
        LINALG::Matrix<my::numnod_, 1> myepreaf(true);

        LINALG::Matrix<my::numdim_, my::numnod_> mydisp(true);
        LINALG::Matrix<my::numnod_, 1> myporosity(true);
        my::ExtractValuesFromGlobalVector(
            discretization, 0, la[0].lm_, &mydisp, &myporosity, "displacement");

        if (discretization.HasState(0, "velocity"))
          my::ExtractValuesFromGlobalVector(discretization, 0, la[0].lm_, &myvel, NULL, "velocity");

        if (discretization.HasState(1, "fluidvel"))
          // extract local values of the global vectors
          my::ExtractValuesFromGlobalVector(
              discretization, 1, la[1].lm_, &myfluidvel, &myepreaf, "fluidvel");

        coupling_poroelast(lm, mydisp, myvel, &myporosity, myfluidvel, myepreaf, matptr,  // NULL,
            NULL, NULL, params);
      }
    }
    break;

    //==================================================================================
    // nonlinear stiffness and internal force vector for poroelasticity
    case ELEMENTS::struct_calc_internalforce:
    {
      // internal force vector
      LINALG::Matrix<numdof_, 1> elevec1(elevec1_epetra.A(), true);
      // elemat2,elevec2+3 are not used anyway

      // build the location vector only for the structure field
      std::vector<int> lm = la[0].lm_;

      LINALG::Matrix<my::numdim_, my::numnod_> mydisp(true);
      LINALG::Matrix<my::numnod_, 1> myporosity(true);
      my::ExtractValuesFromGlobalVector(
          discretization, 0, la[0].lm_, &mydisp, &myporosity, "displacement");

      LINALG::Matrix<my::numdim_, my::numnod_> myvel(true);

      LINALG::Matrix<my::numdim_, my::numnod_> myfluidvel(true);
      LINALG::Matrix<my::numnod_, 1> myepreaf(true);

      // need current fluid state,
      // call the fluid discretization: fluid equates 2nd dofset
      // disassemble velocities and pressures
      if (discretization.HasState(1, "fluidvel"))
      {
        // extract local values of the global vectors
        my::ExtractValuesFromGlobalVector(
            discretization, 1, la[1].lm_, &myfluidvel, &myepreaf, "fluidvel");

        my::ExtractValuesFromGlobalVector(discretization, 0, la[0].lm_, &myvel, NULL, "velocity");

        nlnstiff_poroelast(
            lm, mydisp, myvel, &myporosity, myfluidvel, myepreaf, NULL, NULL, &elevec1, params);
      }
    }
    break;
    //==================================================================================
    default:
      // do nothing (no error because there are some actions the poro element is supposed to ignore)
      break;
  }  // action
  return 0;
}

/*----------------------------------------------------------------------*
 *                                                            vuong 07/13|
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::Wall1_PoroP1<distype>::InitElement()
{
  // initialize base element
  my::InitElement();
  return;
}

/*----------------------------------------------------------------------*
 *                                                            vuong 07/13|
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::Wall1_PoroP1<distype>::nlnstiff_poroelast(
    std::vector<int>& lm,                            ///< location matrix
    LINALG::Matrix<my::numdim_, my::numnod_>& disp,  ///< current displacements
    LINALG::Matrix<my::numdim_, my::numnod_>& vel,   ///< current velocities
    LINALG::Matrix<my::numnod_, 1>* porosity_dof,
    LINALG::Matrix<my::numdim_, my::numnod_>& evelnp,  ///< fluid velocity of element
    LINALG::Matrix<my::numnod_, 1>& epreaf,            ///< fluid pressure of element
    LINALG::Matrix<numdof_, numdof_>* stiffmatrix,     ///< element stiffness matrix
    LINALG::Matrix<numdof_, numdof_>* reamatrix,       // element reactive matrix
    LINALG::Matrix<numdof_, 1>* force,                 ///< element internal force vector
    Teuchos::ParameterList& params                     ///< algorithmic parameters e.g. time
)
{
  my::GetMaterials();

  // update element geometry
  LINALG::Matrix<my::numdim_, my::numnod_> xrefe;  // material coord. of element
  LINALG::Matrix<my::numdim_, my::numnod_> xcurr;  // current  coord. of element

  DRT::Node** nodes = my::Nodes();
  for (int i = 0; i < my::numnod_; ++i)
  {
    const double* x = nodes[i]->X();
    for (int j = 0; j < my::numdim_; j++)
    {
      xrefe(j, i) = x[j];
      xcurr(j, i) = xrefe(j, i) + disp(j, i);
    }
  }

  // initialize element matrizes and vectors
  LINALG::Matrix<my::numdof_, my::numdof_> erea_v(true);
  LINALG::Matrix<my::numdof_, my::numdof_> sub_stiff(true);
  LINALG::Matrix<my::numdof_, 1> sub_force(true);

  LINALG::Matrix<my::numdof_, my::numnod_> ecoupl_p1(true);
  LINALG::Matrix<my::numnod_, numdof_> estiff_p1(true);
  LINALG::Matrix<my::numnod_, 1> ecoupl_force_p1(true);

  /* =========================================================================*/
  /* ================================================= Loop over Gauss Points */
  /* =========================================================================*/
  GaussPointLoopP1(params, xrefe, xcurr, disp, vel, evelnp, epreaf, porosity_dof, erea_v,
      &sub_stiff, &sub_force, ecoupl_p1, estiff_p1, ecoupl_force_p1);

  // update stiffness matrix
  if (stiffmatrix != NULL)
  {
    if (reamatrix != NULL)
    {
      /* additional "reactive darcy-term"
       detJ * w(gp) * ( J * reacoeff * phi^2  ) * D(v_s)
       */
      for (int k = 0; k < my::numnod_; k++)
        for (int l = 0; l < my::numdim_; l++)
          for (int i = 0; i < my::numnod_; i++)
            for (int j = 0; j < my::numdim_; j++)
              (*reamatrix)(i * noddof_ + j, k * noddof_ + l) +=
                  erea_v(i * my::numdim_ + j, k * my::numdim_ + l);
    }

    for (int k = 0; k < my::numnod_; k++)
    {
      for (int l = 0; l < my::numdim_; l++)
        for (int i = 0; i < my::numnod_; i++)
        {
          for (int j = 0; j < my::numdim_; j++)
            (*stiffmatrix)(i * noddof_ + j, k * noddof_ + l) +=
                sub_stiff(i * my::numdim_ + j, k * my::numdim_ + l);
        }
      for (int i = 0; i < my::numnod_; i++)
        for (int j = 0; j < my::numdim_; j++)
          (*stiffmatrix)(i * noddof_ + j, k * noddof_ + my::numdim_) +=
              ecoupl_p1(i * my::noddof_ + j, k);
    }

    for (int i = 0; i < my::numnod_; i++)
      for (int j = 0; j < my::numnod_; j++)
        for (int k = 0; k < noddof_; k++)
          (*stiffmatrix)(i * noddof_ + my::numdim_, j * noddof_ + k) +=
              estiff_p1(i, j * noddof_ + k);
  }

  // update internal force vector
  if (force != NULL)
  {
    for (int i = 0; i < my::numnod_; i++)
    {
      for (int j = 0; j < my::numdim_; j++)
        (*force)(i * noddof_ + j) += sub_force(i * my::numdim_ + j);

      (*force)(i * noddof_ + my::numdim_) += ecoupl_force_p1(i);
    }
  }  // if (force != NULL )

  return;
}  // nlnstiff_poroelast()

/*----------------------------------------------------------------------*
 |  evaluate only the poroelasticity fraction for the element (protected) vuong 07/13 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::Wall1_PoroP1<distype>::GaussPointLoopP1(Teuchos::ParameterList& params,
    const LINALG::Matrix<my::numdim_, my::numnod_>& xrefe,
    const LINALG::Matrix<my::numdim_, my::numnod_>& xcurr,
    const LINALG::Matrix<my::numdim_, my::numnod_>& nodaldisp,
    const LINALG::Matrix<my::numdim_, my::numnod_>& nodalvel,
    const LINALG::Matrix<my::numdim_, my::numnod_>& evelnp,
    const LINALG::Matrix<my::numnod_, 1>& epreaf,
    const LINALG::Matrix<my::numnod_, 1>* porosity_dof,
    LINALG::Matrix<my::numdof_, my::numdof_>& erea_v,
    LINALG::Matrix<my::numdof_, my::numdof_>* sub_stiff, LINALG::Matrix<my::numdof_, 1>* sub_force,
    LINALG::Matrix<my::numdof_, my::numnod_>& ecoupl_p1,
    LINALG::Matrix<my::numnod_, numdof_>& estiff_p1,
    LINALG::Matrix<my::numnod_, 1>& ecoupl_force_p1)
{
  LINALG::Matrix<my::numdim_, my::numnod_> N_XYZ;
  // build deformation gradient wrt to material configuration
  // in case of prestressing, build defgrd wrt to last stored configuration
  // CAUTION: defgrd(true): filled with zeros!
  LINALG::Matrix<my::numdim_, my::numdim_> defgrd(true);
  LINALG::Matrix<my::numnod_, 1> shapefct;
  LINALG::Matrix<my::numdim_, my::numnod_> deriv;

  LINALG::Matrix<my::numstr_, 1> fstress(true);

  for (int gp = 0; gp < my::numgpt_; ++gp)
  {
    // evaluate shape functions and derivatives at integration point
    my::ComputeShapeFunctionsAndDerivatives(gp, shapefct, deriv, N_XYZ);

    // compute deformation gradient
    my::ComputeDefGradient(defgrd, N_XYZ, xcurr);

    // inverse deformation gradient F^-1
    LINALG::Matrix<my::numdim_, my::numdim_> defgrd_inv(false);
    defgrd_inv.Invert(defgrd);

    // jacobian determinant of transformation between spatial and material space "|dx/dX|"
    double J = 0.0;
    //------linearization of jacobi determinant detF=J w.r.t. structure displacement   dJ/d(us) =
    //dJ/dF : dF/dus = J * F^-T * N,X
    static LINALG::Matrix<1, my::numdof_> dJ_dus;
    // volume change (used for porosity law). Same as J in nonlinear theory.
    double volchange = 0.0;
    //------linearization of volume change w.r.t. structure displacement
    static LINALG::Matrix<1, my::numdof_> dvolchange_dus;

    // compute J, the volume change and the respctive linearizations w.r.t. structure displacement
    my::ComputeJacobianDeterminantVolumeChangeAndLinearizations(
        J, volchange, dJ_dus, dvolchange_dus, defgrd, defgrd_inv, N_XYZ, nodaldisp);

    //----------------------------------------------------
    // pressure at integration point
    double press = shapefct.Dot(epreaf);

    // structure velocity at integration point
    LINALG::Matrix<my::numdim_, 1> velint(true);

    for (int i = 0; i < my::numnod_; i++)
      for (int j = 0; j < my::numdim_; j++) velint(j) += nodalvel(j, i) * shapefct(i);

    // fluid velocity at integration point
    LINALG::Matrix<my::numdim_, 1> fvelint;
    fvelint.Multiply(evelnp, shapefct);

    // material fluid velocity gradient at integration point
    LINALG::Matrix<my::numdim_, my::numdim_> fvelder;
    fvelder.MultiplyNT(evelnp, N_XYZ);

    // pressure gradient at integration point
    LINALG::Matrix<my::numdim_, 1> Gradp;
    Gradp.Multiply(N_XYZ, epreaf);

    // non-linear B-operator
    LINALG::Matrix<my::numstr_, my::numdof_> bop;
    my::ComputeBOperator(bop, defgrd, N_XYZ);

    // Right Cauchy-Green tensor = F^T * F
    LINALG::Matrix<my::numdim_, my::numdim_> cauchygreen;
    cauchygreen.MultiplyTN(defgrd, defgrd);

    // inverse Right Cauchy-Green tensor
    LINALG::Matrix<my::numdim_, my::numdim_> C_inv(false);
    C_inv.Invert(cauchygreen);

    //------linearization of material gradient of jacobi determinant GradJ  w.r.t. strucuture
    //displacement d(GradJ)/d(us)
    //---------------------d(GradJ)/dus =  dJ/dus * F^-T . : dF/dX + J * dF^-T/dus : dF/dX + J *
    //F^-T : N_X_X

    // compute some auxiliary matrixes for computation of linearization
    // dF^-T/dus
    LINALG::Matrix<my::numdim_ * my::numdim_, my::numdof_> dFinvTdus(true);
    // F^-T * Grad p
    LINALG::Matrix<my::numdim_, 1> Finvgradp;
    // dF^-T/dus * Grad p
    LINALG::Matrix<my::numdim_, my::numdof_> dFinvdus_gradp(true);
    // dC^-1/dus * Grad p
    LINALG::Matrix<my::numstr_, my::numdof_> dCinv_dus(true);

    my::ComputeAuxiliaryValues(
        N_XYZ, defgrd_inv, C_inv, Gradp, dFinvTdus, Finvgradp, dFinvdus_gradp, dCinv_dus);

    //--------------------------------------------------------------------

    // linearization of porosity w.r.t structure displacement d\phi/d(us) = d\phi/dJ*dJ/d(us)
    LINALG::Matrix<1, my::numdof_> dphi_dus;
    double porosity = 0.0;

    ComputePorosityAndLinearization(
        params, press, volchange, gp, shapefct, porosity_dof, dvolchange_dus, porosity, dphi_dus);

    double dW_dphi = 0.0;
    double dW_dJ = 0.0;
    double dW_dp = 0.0;
    double W = 0.0;
    my::structmat_->ConstitutiveDerivatives(params, press, volchange, porosity,
        &dW_dp,  // dW_dp not needed
        &dW_dphi, &dW_dJ, NULL, &W);

    //--------------------------------------------------------

    // **********************evaluate stiffness matrix and force vector+++++++++++++++++++++++++
    if (my::fluidmat_->Type() == MAT::PAR::darcy_brinkman)
    {
      my::FillMatrixAndVectorsBrinkman(gp, J, porosity, fvelder, defgrd_inv, bop, C_inv, dphi_dus,
          dJ_dus, dCinv_dus, dFinvTdus, sub_stiff, sub_force, fstress);
    }

    my::FillMatrixAndVectors(gp, shapefct, N_XYZ, J, press, porosity, velint, fvelint, fvelder,
        defgrd_inv, bop, C_inv, Finvgradp, dphi_dus, dJ_dus, dCinv_dus, dFinvdus_gradp, dFinvTdus,
        erea_v, sub_stiff, sub_force, fstress);

    // **********************evaluate stiffness matrix and force vector+++++++++++++++++++++++++
    double detJ_w = my::detJ_[gp] * my::intpoints_.Weight(gp);  // gpweights[gp];

    const double reacoeff = my::fluidmat_->ComputeReactionCoeff();
    // if (force != NULL or stiffmatrix != NULL or reamatrix != NULL )
    {
      for (int k = 0; k < my::numnod_; k++)
      {
        const double fac = detJ_w * shapefct(k);

        ecoupl_force_p1(k) += fac * W;

        for (int i = 0; i < my::numnod_; i++)
        {
          for (int j = 0; j < my::numdim_; j++)
          {
            estiff_p1(k, i * noddof_ + j) += fac * dW_dJ * dJ_dus(i * my::numdim_ + j);

            ecoupl_p1(i * my::numdim_ + j, k) +=
                fac * (2 * J * reacoeff * porosity * (velint(j) - fvelint(j)) + J * Finvgradp(j)) *
                shapefct(i);
          }
          estiff_p1(k, i * noddof_ + my::numdim_) += fac * dW_dphi * shapefct(i);
        }
      }
    }

    if (my::fluidmat_->Type() == MAT::PAR::darcy_brinkman)
    {
      double visc = my::fluidmat_->Viscosity();
      LINALG::Matrix<my::numdim_, my::numdim_> CinvFvel(true);
      LINALG::Matrix<my::numdim_, my::numdim_> visctress1(true);
      CinvFvel.Multiply(C_inv, fvelder);
      visctress1.MultiplyNT(CinvFvel, defgrd_inv);
      LINALG::Matrix<my::numdim_, my::numdim_> visctress2(visctress1);
      visctress1.UpdateT(1.0, visctress2, 1.0);

      fstress(0) = visctress1(0, 0);
      fstress(1) = visctress1(1, 1);
      fstress(2) = visctress1(0, 1);

      fstress.Scale(detJ_w * visc * J);

      // B^T . C^-1
      LINALG::Matrix<my::numdof_, 1> fstressb(true);
      fstressb.MultiplyTN(bop, fstress);

      for (int k = 0; k < my::numnod_; k++)
      {
        const double fac = detJ_w * shapefct(k);
        for (int i = 0; i < my::numnod_; i++)
          for (int j = 0; j < my::numdim_; j++)
            ecoupl_p1(i * my::numdim_ + j, k) += fac * fstressb(i * my::numdim_ + j) * shapefct(i);
      }
    }
    /* =========================================================================*/
  } /* ==================================================== end of Loop over GP */
  /* =========================================================================*/
}

/*----------------------------------------------------------------------*
 |  evaluate only the poroelasticity fraction for the element (protected) |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::Wall1_PoroP1<distype>::coupling_poroelast(
    std::vector<int>& lm,                            // location matrix
    LINALG::Matrix<my::numdim_, my::numnod_>& disp,  // current displacements
    LINALG::Matrix<my::numdim_, my::numnod_>& vel,   // current velocities
    LINALG::Matrix<my::numnod_, 1>* porosity,
    LINALG::Matrix<my::numdim_, my::numnod_>& evelnp,  // current fluid velocity
    LINALG::Matrix<my::numnod_, 1>& epreaf,            // current fluid pressure
    LINALG::Matrix<numdof_, (my::numdim_ + 1) * my::numnod_>*
        stiffmatrix,  // element stiffness matrix
    LINALG::Matrix<numdof_, (my::numdim_ + 1) * my::numnod_>* reamatrix,  // element reactive matrix
    LINALG::Matrix<numdof_, 1>* force,  // element internal force vector
    Teuchos::ParameterList& params)     // algorithmic parameters e.g. time
{
  my::GetMaterials();

  //=======================================================================

  // update element geometry
  LINALG::Matrix<my::numdim_, my::numnod_> xrefe;  // material coord. of element
  LINALG::Matrix<my::numdim_, my::numnod_> xcurr;  // current  coord. of element

  DRT::Node** nodes = my::Nodes();
  for (int i = 0; i < my::numnod_; ++i)
  {
    const double* x = nodes[i]->X();
    for (int j = 0; j < my::numdim_; j++)
    {
      xrefe(j, i) = x[j];
      xcurr(j, i) = xrefe(j, i) + disp(j, i);
    }
  }
  // initialize element matrizes
  LINALG::Matrix<my::numdof_, (my::numdim_ + 1) * my::numnod_> ecoupl(true);

  LINALG::Matrix<my::numnod_, my::numnod_> ecoupl_p1_p(true);

  /* =========================================================================*/
  /* ================================================= Loop over Gauss Points */
  /* =========================================================================*/
  GaussPointLoopP1OD(
      params, xrefe, xcurr, disp, vel, evelnp, epreaf, porosity, ecoupl_p1_p, ecoupl);

  // if (force != NULL )
  //{
  // all rhs terms are added in nlnstiff_poroelast
  //}

  if (stiffmatrix != NULL)
  {
    for (int k = 0; k < my::numnod_; k++)
      for (int l = 0; l < (my::numdim_ + 1); l++)
        for (int i = 0; i < my::numnod_; i++)
          for (int j = 0; j < my::numdim_; j++)
            (*stiffmatrix)(i * noddof_ + j, k * (my::numdim_ + 1) + l) +=
                ecoupl(i * my::numdim_ + j, k * (my::numdim_ + 1) + l);

    for (int ui = 0; ui < my::numnod_; ++ui)
      for (int ni = 0; ni < my::numnod_; ++ni)
        (*stiffmatrix)(noddof_ * ui + my::numdim_, (my::numdim_ + 1) * ni + my::numdim_) +=
            ecoupl_p1_p(ui, ni);
  }

  return;

}  // coupling_poroelast()

/*----------------------------------------------------------------------*
 |  evaluate only the poroelasticity fraction for the element (protected) vuong 07/13|
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::Wall1_PoroP1<distype>::GaussPointLoopP1OD(Teuchos::ParameterList& params,
    const LINALG::Matrix<my::numdim_, my::numnod_>& xrefe,
    const LINALG::Matrix<my::numdim_, my::numnod_>& xcurr,
    const LINALG::Matrix<my::numdim_, my::numnod_>& nodaldisp,
    const LINALG::Matrix<my::numdim_, my::numnod_>& nodalvel,
    const LINALG::Matrix<my::numdim_, my::numnod_>& evelnp,
    const LINALG::Matrix<my::numnod_, 1>& epreaf,
    const LINALG::Matrix<my::numnod_, 1>* porosity_dof,
    LINALG::Matrix<my::numnod_, my::numnod_>& ecoupl_p1,
    LINALG::Matrix<my::numdof_, (my::numdim_ + 1) * my::numnod_>& sub_stiff)
{
  LINALG::Matrix<my::numdim_, my::numnod_> N_XYZ;  //  first derivatives at gausspoint w.r.t. X, Y,Z
  // build deformation gradient wrt to material configuration
  // in case of prestressing, build defgrd wrt to last stored configuration
  // CAUTION: defgrd(true): filled with zeros!
  LINALG::Matrix<my::numdim_, my::numdim_> defgrd(
      true);                                //  deformation gradiant evaluated at gauss point
  LINALG::Matrix<my::numnod_, 1> shapefct;  //  shape functions evalulated at gauss point
  LINALG::Matrix<my::numdim_, my::numnod_> deriv(
      true);  //  first derivatives at gausspoint w.r.t. r,s,t

  for (int gp = 0; gp < my::numgpt_; ++gp)
  {
    // evaluate shape functions and derivatives at integration point
    my::ComputeShapeFunctionsAndDerivatives(gp, shapefct, deriv, N_XYZ);
    // evaluate second derivatives of shape functions at integration point
    // ComputeSecondDerivativesOfShapeFunctions(gp,xrefe,deriv,deriv2,N_XYZ,N_XYZ2);

    // (material) deformation gradient F = d xcurr / d xrefe = xcurr * N_XYZ^T
    my::ComputeDefGradient(defgrd, N_XYZ, xcurr);

    // inverse deformation gradient F^-1
    LINALG::Matrix<my::numdim_, my::numdim_> defgrd_inv(false);
    defgrd_inv.Invert(defgrd);

    // jacobian determinant of transformation between spatial and material space "|dx/dX|"
    double J = 0.0;
    // volume change (used for porosity law). Same as J in nonlinear theory.
    double volchange = 0.0;

    // compute J, the volume change and the respctive linearizations w.r.t. structure displacement
    my::ComputeJacobianDeterminantVolumeChange(J, volchange, defgrd, N_XYZ, nodaldisp);

    // non-linear B-operator
    LINALG::Matrix<my::numstr_, my::numdof_> bop;
    my::ComputeBOperator(bop, defgrd, N_XYZ);

    // -----------------Right Cauchy-Green tensor = F^T * F
    LINALG::Matrix<my::numdim_, my::numdim_> cauchygreen;
    cauchygreen.MultiplyTN(defgrd, defgrd);

    //------------------ inverse Right Cauchy-Green tensor
    LINALG::Matrix<my::numdim_, my::numdim_> C_inv(false);
    C_inv.Invert(cauchygreen);

    //---------------- get pressure at integration point
    double press = shapefct.Dot(epreaf);

    //------------------ get material pressure gradient at integration point
    LINALG::Matrix<my::numdim_, 1> Gradp;
    Gradp.Multiply(N_XYZ, epreaf);

    //--------------------- get fluid velocity at integration point
    LINALG::Matrix<my::numdim_, 1> fvelint;
    fvelint.Multiply(evelnp, shapefct);

    // material fluid velocity gradient at integration point
    LINALG::Matrix<my::numdim_, my::numdim_> fvelder;
    fvelder.MultiplyNT(evelnp, N_XYZ);

    //! ----------------structure velocity at integration point
    LINALG::Matrix<my::numdim_, 1> velint(true);
    for (int i = 0; i < my::numnod_; i++)
      for (int j = 0; j < my::numdim_; j++) velint(j) += nodalvel(j, i) * shapefct(i);

    //**************************************************+auxilary variables for computing the
    //porosity and linearization
    double dphi_dp = 0.0;
    double porosity = 0.0;

    ComputePorosityAndLinearizationOD(
        params, press, volchange, gp, shapefct, porosity_dof, porosity, dphi_dp);

    // **********************evaluate stiffness matrix and force vector+++++++++++++++++++++++++

    my::FillMatrixAndVectorsOD(gp, shapefct, N_XYZ, J, porosity, dphi_dp, velint, fvelint,
        defgrd_inv, Gradp, bop, C_inv, sub_stiff);

    if (my::fluidmat_->Type() == MAT::PAR::darcy_brinkman)
    {
      my::FillMatrixAndVectorsBrinkmanOD(
          gp, shapefct, N_XYZ, J, porosity, dphi_dp, fvelder, defgrd_inv, bop, C_inv, sub_stiff);
    }  // darcy-brinkman

    double dW_dp = 0.0;
    my::structmat_->ConstitutiveDerivatives(params, press, J, porosity, &dW_dp,
        NULL,  // not needed
        NULL,  // not needed
        NULL,  // not needed
        NULL   // not needed
    );
    //--------------------------------------------------------

    // **********************evaluate stiffness matrix and force vector+++++++++++++++++++++++++
    double detJ_w = my::detJ_[gp] * my::intpoints_.Weight(gp);  // gpweights[gp];

    for (int k = 0; k < my::numnod_; k++)
    {
      const double fac = detJ_w * shapefct(k);

      for (int i = 0; i < my::numnod_; i++) ecoupl_p1(k, i) += fac * dW_dp * shapefct(i);
    }

    /* =========================================================================*/
  } /* ==================================================== end of Loop over GP */
  /* =========================================================================*/
}

/*----------------------------------------------------------------------*
 |  Integrate a Surface Neumann boundary condition (public)  |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
int DRT::ELEMENTS::Wall1_PoroP1<distype>::EvaluateNeumann(Teuchos::ParameterList& params,
    DRT::Discretization& discretization, DRT::Condition& condition, std::vector<int>& lm,
    Epetra_SerialDenseVector& elevec1, Epetra_SerialDenseMatrix* elemat1)
{
  LINALG::Matrix<my::numdim_, my::numnod_> disp(true);
  LINALG::Matrix<my::numnod_, 1> myporosity(true);
  my::ExtractValuesFromGlobalVector(discretization, 0, lm, &disp, &myporosity, "displacement");

  // find out whether we will use a time curve
  const double time = params.get("total time", -1.0);

  /*----------------------------------------------------- geometry update */
  // update element geometry
  LINALG::Matrix<my::numdim_, my::numnod_> xrefe;  // material coord. of element
  LINALG::Matrix<my::numdim_, my::numnod_> xcurr;  // current  coord. of element

  DRT::Node** nodes = my::Nodes();
  for (int i = 0; i < my::numnod_; ++i)
  {
    const double* x = nodes[i]->X();
    for (int j = 0; j < my::numdim_; j++)
    {
      xrefe(j, i) = x[j];
      xcurr(j, i) = xrefe(j, i) + disp(j, i);
    }
  }


  // get values and switches from the condition
  const std::vector<int>* onoff = condition.Get<std::vector<int>>("onoff");
  const std::vector<double>* val = condition.Get<std::vector<double>>("val");
  const std::vector<int>* funct = condition.Get<std::vector<int>>("funct");


  LINALG::Matrix<my::numdim_, my::numnod_> N_XYZ;
  // build deformation gradient wrt to material configuration
  // in case of prestressing, build defgrd wrt to last stored configuration
  // CAUTION: defgrd(true): filled with zeros!
  LINALG::Matrix<my::numdim_, my::numdim_> defgrd(true);
  LINALG::Matrix<my::numnod_, 1> shapefcts;
  LINALG::Matrix<my::numdim_, my::numnod_> deriv;

  LINALG::Matrix<my::numstr_, 1> fstress(true);

  for (int gp = 0; gp < my::numgpt_; ++gp)
  {
    // evaluate shape functions and derivatives at integration point
    my::ComputeShapeFunctionsAndDerivatives(gp, shapefcts, deriv, N_XYZ);

    // jacobian determinant of transformation between spatial and material space "|dx/dX|"
    my::ComputeJacobianDeterminant(gp, xcurr, deriv);

    /*------------------------------------ integration factor  -------*/
    double fac = my::detJ_[gp] * my::intpoints_.Weight(gp);

    // load vector ar
    double ar[my::numdim_] = {0.0, 0.0};
    // loop the dofs of a node
    // ar[i] = ar[i] * facr * ds * onoff[i] * val[i]
    for (int i = 0; i < my::numdim_; ++i)
    {
      if ((*onoff)[i])
      {
        // factor given by spatial function
        const int functnum = (funct) ? (*funct)[i] : -1;
        double functfac = 1.0;
        if (functnum > 0)
        {
          // calculate reference position of GP
          LINALG::Matrix<1, my::numdim_> gp_coord;
          gp_coord.MultiplyTT(shapefcts, xrefe);

          // write coordinates in another datatype
          double gp_coord2[3];  // the position vector has to be given in 3D!!!
          for (int k = 0; k < my::numdim_; k++) gp_coord2[k] = gp_coord(0, k);
          for (int k = my::numdim_; k < 3;
               k++)  // set a zero value for the remaining spatial directions
            gp_coord2[k] = 0.0;
          const double* coordgpref = &gp_coord2[0];  // needed for function evaluation

          // evaluate function at current gauss point
          functfac = DRT::Problem::Instance()->Funct(functnum - 1).Evaluate(i, coordgpref, time);
        }

        ar[i] = fac * (*val)[i] * functfac;
      }
    }

    // add load components
    for (int node = 0; node < my::numnod_; ++node)
      for (int dim = 0; dim < my::numdim_; ++dim)
        elevec1[node * noddof_ + dim] += shapefcts(node) * ar[dim];

  }  // for (int ip=0; ip<totngp; ++ip)

  // finished
  return 0;
}

/*----------------------------------------------------------------------*
 *                                                            vuong 07/13|
 *----------------------------------------------------------------------*/
template class DRT::ELEMENTS::Wall1_PoroP1<DRT::Element::tri3>;
template class DRT::ELEMENTS::Wall1_PoroP1<DRT::Element::quad4>;
template class DRT::ELEMENTS::Wall1_PoroP1<DRT::Element::quad9>;
