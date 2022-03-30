/*----------------------------------------------------------------------------*/
/*! \file
\brief Evaluate methods for 2D wall element for structure part of porous medium
       using p1 approach (mixed approach).

\level 2


*/
/*---------------------------------------------------------------------------*/

#include "wall1_poro_p1.H"
#include "wall1_poro_p1_eletypes.H"

#include "../drt_lib/drt_discret.H"

#include "../drt_mat/fluidporo.H"
#include "../drt_mat/structporo.H"

#include "../drt_structure_new/str_elements_paramsinterface.H"

template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::Wall1_PoroP1<distype>::ComputePorosityAndLinearization(
    Teuchos::ParameterList& params, const double& press, const double& J, const int& gp,
    const LINALG::Matrix<Base::numnod_, 1>& shapfct,
    const LINALG::Matrix<Base::numnod_, 1>* myporosity,
    const LINALG::Matrix<1, Base::numdof_>& dJ_dus, double& porosity,
    LINALG::Matrix<1, Base::numdof_>& dphi_dus)
{
  if (myporosity == nullptr)
    dserror("no porosity values given!");
  else
    porosity = shapfct.Dot(*myporosity);

  dphi_dus.PutScalar(0.0);
}

template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::Wall1_PoroP1<distype>::ComputePorosityAndLinearizationOD(
    Teuchos::ParameterList& params, const double& press, const double& J, const int& gp,
    const LINALG::Matrix<Base::numnod_, 1>& shapfct,
    const LINALG::Matrix<Base::numnod_, 1>* myporosity, double& porosity, double& dphi_dp)
{
  dphi_dp = 0.0;

  if (myporosity == nullptr)
    dserror("no porosity values given!");
  else
    porosity = shapfct.Dot(*myporosity);
}

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
      Base::PreEvaluate(params, discretization, la);

      MyEvaluate(params, discretization, la, elemat1_epetra, elemat2_epetra, elevec1_epetra,
          elevec2_epetra, elevec3_epetra);
    }
    break;
    //==================================================================================
    case ELEMENTS::struct_calc_energy:
    {
      // in some cases we need to write/change some data before evaluating
      Base::PreEvaluate(params, discretization, la);

      // evaluate parent solid element
      Wall1::Evaluate(params, discretization, la[0].lm_, elemat1_epetra, elemat2_epetra,
          elevec1_epetra, elevec2_epetra, elevec3_epetra);
    }
    break;
    //==================================================================================
    default:
    {
      // in some cases we need to write/change some data before evaluating
      Base::PreEvaluate(params, discretization, la);

      Epetra_SerialDenseMatrix elemat1_sub;
      Epetra_SerialDenseMatrix elemat2_sub;
      Epetra_SerialDenseVector elevec1_sub;
      Epetra_SerialDenseVector elevec2_sub;
      Epetra_SerialDenseVector elevec3_sub;

      if (elemat1_epetra.A()) elemat1_sub.Shape(Base::numdof_, Base::numdof_);
      if (elemat2_epetra.A()) elemat2_sub.Shape(Base::numdof_, Base::numdof_);
      if (elevec1_epetra.A()) elevec1_sub.Resize(Base::numdof_);
      if (elevec2_epetra.A()) elevec2_sub.Resize(Base::numdof_);
      if (elevec3_epetra.A()) elevec3_sub.Resize(Base::numdof_);

      std::vector<int> lm_sub;
      for (int i = 0; i < Base::numnod_; i++)
        for (int j = 0; j < Base::numdim_; j++) lm_sub.push_back(la[0].lm_[i * noddof_ + j]);

      // evaluate parent solid element
      Wall1::Evaluate(params, discretization, lm_sub, elemat1_sub, elemat2_sub, elevec1_sub,
          elevec2_sub, elevec3_sub);


      if (elemat1_epetra.A())
      {
        for (int i = 0; i < Base::numnod_; i++)
        {
          for (int j = 0; j < Base::numdim_; j++)
          {
            for (int k = 0; k < Base::numnod_; k++)
            {
              for (int l = 0; l < Base::numdim_; l++)
                elemat1_epetra(i * noddof_ + j, k * noddof_ + l) =
                    elemat1_sub(i * Base::noddof_ + j, k * Base::noddof_ + l);
            }
          }
        }
      }

      if (elemat2_epetra.A())
      {
        for (int i = 0; i < Base::numnod_; i++)
        {
          for (int j = 0; j < Base::numdim_; j++)
          {
            for (int k = 0; k < Base::numnod_; k++)
            {
              for (int l = 0; l < Base::numdim_; l++)
                elemat2_epetra(i * noddof_ + j, k * noddof_ + l) =
                    elemat2_sub(i * Base::noddof_ + j, k * Base::noddof_ + l);
            }
          }
        }
      }

      if (elevec1_epetra.A())
      {
        for (int i = 0; i < Base::numnod_; i++)
          for (int j = 0; j < Base::numdim_; j++)
            elevec1_epetra(i * noddof_ + j) = elevec1_sub(i * Base::noddof_ + j);
      }

      if (elevec2_epetra.A())
      {
        for (int i = 0; i < Base::numnod_; i++)
          for (int j = 0; j < Base::numdim_; j++)
            elevec2_epetra(i * noddof_ + j) = elevec2_sub(i * Base::noddof_ + j);
      }

      if (elevec3_epetra.A())
      {
        for (int i = 0; i < Base::numnod_; i++)
          for (int j = 0; j < Base::numdim_; j++)
            elevec3_epetra(i * noddof_ + j) = elevec3_sub(i * Base::noddof_ + j);
      }

      // add volume coupling specific terms
      MyEvaluate(params, discretization, la, elemat1_epetra, elemat2_epetra, elevec1_epetra,
          elevec2_epetra, elevec3_epetra);
    }
    break;
  }

  return 0;
}

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
        // elevec2+3 are not used anyway

        // build the location vector only for the structure field
        std::vector<int> lm = la[0].lm_;

        LINALG::Matrix<Base::numdim_, Base::numnod_> mydisp(true);
        LINALG::Matrix<Base::numnod_, 1> myporosity(true);
        Base::ExtractValuesFromGlobalVector(
            discretization, 0, la[0].lm_, &mydisp, &myporosity, "displacement");

        LINALG::Matrix<numdof_, numdof_>* matptr = nullptr;
        if (elemat1.IsInitialized()) matptr = &elemat1;

        enum INPAR::STR::DampKind damping =
            params.get<enum INPAR::STR::DampKind>("damping", INPAR::STR::damp_none);
        LINALG::Matrix<numdof_, numdof_>* matptr2 = nullptr;
        if (elemat2.IsInitialized() and (damping == INPAR::STR::damp_material)) matptr2 = &elemat2;

        // need current fluid state,
        // call the fluid discretization: fluid equates 2nd dofset
        // disassemble velocities and pressures

        LINALG::Matrix<Base::numdim_, Base::numnod_> myvel(true);

        LINALG::Matrix<Base::numdim_, Base::numnod_> myfluidvel(true);
        LINALG::Matrix<Base::numnod_, 1> myepreaf(true);

        if (discretization.HasState(0, "velocity"))
          Base::ExtractValuesFromGlobalVector(
              discretization, 0, la[0].lm_, &myvel, nullptr, "velocity");

        if (discretization.HasState(1, "fluidvel"))
        {
          // extract local values of the global vectors
          Base::ExtractValuesFromGlobalVector(
              discretization, 1, la[1].lm_, &myfluidvel, &myepreaf, "fluidvel");
        }

        // calculate tangent stiffness matrix
        NonlinearStiffnessPoroelast(lm, mydisp, myvel, &myporosity, myfluidvel, myepreaf, matptr,
            matptr2, &elevec1, params);
      }
    }
    break;

    //==================================================================================
    // coupling terms in force-vector and stiffness matrix for poroelasticity
    case ELEMENTS::struct_poro_calc_fluidcoupling:
    {
      // stiffness
      LINALG::Matrix<numdof_, (Base::numdim_ + 1) * Base::numnod_> elemat1(
          elemat1_epetra.A(), true);

      // elemat2,elevec1-3 are not used anyway

      // build the location vector only for the structure field
      std::vector<int> lm = la[0].lm_;

      LINALG::Matrix<numdof_, (Base::numdim_ + 1)* Base::numnod_>* matptr = nullptr;
      if (elemat1.IsInitialized()) matptr = &elemat1;

      // need current fluid state,
      // call the fluid discretization: fluid equates 2nd dofset
      // disassemble velocities and pressures
      if (discretization.HasState(1, "fluidvel"))
      {
        LINALG::Matrix<Base::numdim_, Base::numnod_> myvel(true);
        LINALG::Matrix<Base::numdim_, Base::numnod_> myfluidvel(true);
        LINALG::Matrix<Base::numnod_, 1> myepreaf(true);

        LINALG::Matrix<Base::numdim_, Base::numnod_> mydisp(true);
        LINALG::Matrix<Base::numnod_, 1> myporosity(true);
        Base::ExtractValuesFromGlobalVector(
            discretization, 0, la[0].lm_, &mydisp, &myporosity, "displacement");

        if (discretization.HasState(0, "velocity"))
          Base::ExtractValuesFromGlobalVector(
              discretization, 0, la[0].lm_, &myvel, nullptr, "velocity");

        if (discretization.HasState(1, "fluidvel"))
        {
          // extract local values of the global vectors
          Base::ExtractValuesFromGlobalVector(
              discretization, 1, la[1].lm_, &myfluidvel, &myepreaf, "fluidvel");
        }

        CouplingPoroelast(lm, mydisp, myvel, &myporosity, myfluidvel, myepreaf,
            matptr,  // nullptr,
            nullptr, nullptr, params);
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

      LINALG::Matrix<Base::numdim_, Base::numnod_> mydisp(true);
      LINALG::Matrix<Base::numnod_, 1> myporosity(true);
      Base::ExtractValuesFromGlobalVector(
          discretization, 0, la[0].lm_, &mydisp, &myporosity, "displacement");

      LINALG::Matrix<Base::numdim_, Base::numnod_> myvel(true);

      LINALG::Matrix<Base::numdim_, Base::numnod_> myfluidvel(true);
      LINALG::Matrix<Base::numnod_, 1> myepreaf(true);

      // need current fluid state,
      // call the fluid discretization: fluid equates 2nd dofset
      // disassemble velocities and pressures
      if (discretization.HasState(1, "fluidvel"))
      {
        // extract local values of the global vectors
        Base::ExtractValuesFromGlobalVector(
            discretization, 1, la[1].lm_, &myfluidvel, &myepreaf, "fluidvel");

        Base::ExtractValuesFromGlobalVector(
            discretization, 0, la[0].lm_, &myvel, nullptr, "velocity");

        NonlinearStiffnessPoroelast(lm, mydisp, myvel, &myporosity, myfluidvel, myepreaf, nullptr,
            nullptr, &elevec1, params);
      }
    }
    break;
    //==================================================================================
    default:
      // do nothing (no error because there are some actions the poro element is supposed to ignore)
      break;
  }
  return 0;
}

template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::Wall1_PoroP1<distype>::InitElement()
{
  // initialize base element
  Base::InitElement();
}

template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::Wall1_PoroP1<distype>::NonlinearStiffnessPoroelast(std::vector<int>& lm,
    LINALG::Matrix<Base::numdim_, Base::numnod_>& disp,
    LINALG::Matrix<Base::numdim_, Base::numnod_>& vel,
    LINALG::Matrix<Base::numnod_, 1>* porosity_dof,
    LINALG::Matrix<Base::numdim_, Base::numnod_>& evelnp, LINALG::Matrix<Base::numnod_, 1>& epreaf,
    LINALG::Matrix<numdof_, numdof_>* stiffmatrix, LINALG::Matrix<numdof_, numdof_>* reamatrix,
    LINALG::Matrix<numdof_, 1>* force, Teuchos::ParameterList& params)
{
  Base::GetMaterials();

  // update element geometry
  LINALG::Matrix<Base::numdim_, Base::numnod_> xrefe;  // material coord. of element
  LINALG::Matrix<Base::numdim_, Base::numnod_> xcurr;  // current  coord. of element

  DRT::Node** nodes = Base::Nodes();
  for (int i = 0; i < Base::numnod_; ++i)
  {
    const double* x = nodes[i]->X();
    for (int j = 0; j < Base::numdim_; j++)
    {
      xrefe(j, i) = x[j];
      xcurr(j, i) = xrefe(j, i) + disp(j, i);
    }
  }

  // initialize element matrizes and vectors
  LINALG::Matrix<Base::numdof_, Base::numdof_> erea_v(true);
  LINALG::Matrix<Base::numdof_, Base::numdof_> sub_stiff(true);
  LINALG::Matrix<Base::numdof_, 1> sub_force(true);

  LINALG::Matrix<Base::numdof_, Base::numnod_> ecoupl_p1(true);
  LINALG::Matrix<Base::numnod_, numdof_> estiff_p1(true);
  LINALG::Matrix<Base::numnod_, 1> ecoupl_force_p1(true);

  /* =========================================================================*/
  /* ================================================= Loop over Gauss Points */
  /* =========================================================================*/
  GaussPointLoopP1(params, xrefe, xcurr, disp, vel, evelnp, epreaf, porosity_dof, erea_v,
      &sub_stiff, &sub_force, ecoupl_p1, estiff_p1, ecoupl_force_p1);

  // update stiffness matrix
  if (stiffmatrix != nullptr)
  {
    if (reamatrix != nullptr)
    {
      /* additional "reactive darcy-term"
       detJ * w(gp) * ( J * reacoeff * phi^2  ) * D(v_s)
       */
      for (int k = 0; k < Base::numnod_; k++)
      {
        for (int l = 0; l < Base::numdim_; l++)
        {
          for (int i = 0; i < Base::numnod_; i++)
          {
            for (int j = 0; j < Base::numdim_; j++)
              (*reamatrix)(i * noddof_ + j, k * noddof_ + l) +=
                  erea_v(i * Base::numdim_ + j, k * Base::numdim_ + l);
          }
        }
      }
    }

    for (int k = 0; k < Base::numnod_; k++)
    {
      for (int l = 0; l < Base::numdim_; l++)
      {
        for (int i = 0; i < Base::numnod_; i++)
        {
          for (int j = 0; j < Base::numdim_; j++)
            (*stiffmatrix)(i * noddof_ + j, k * noddof_ + l) +=
                sub_stiff(i * Base::numdim_ + j, k * Base::numdim_ + l);
        }
      }
      for (int i = 0; i < Base::numnod_; i++)
      {
        for (int j = 0; j < Base::numdim_; j++)
          (*stiffmatrix)(i * noddof_ + j, k * noddof_ + Base::numdim_) +=
              ecoupl_p1(i * Base::noddof_ + j, k);
      }
    }

    for (int i = 0; i < Base::numnod_; i++)
    {
      for (int j = 0; j < Base::numnod_; j++)
      {
        for (int k = 0; k < noddof_; k++)
          (*stiffmatrix)(i * noddof_ + Base::numdim_, j * noddof_ + k) +=
              estiff_p1(i, j * noddof_ + k);
      }
    }
  }

  // update internal force vector
  if (force != nullptr)
  {
    for (int i = 0; i < Base::numnod_; i++)
    {
      for (int j = 0; j < Base::numdim_; j++)
        (*force)(i * noddof_ + j) += sub_force(i * Base::numdim_ + j);

      (*force)(i * noddof_ + Base::numdim_) += ecoupl_force_p1(i);
    }
  }
}

template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::Wall1_PoroP1<distype>::GaussPointLoopP1(Teuchos::ParameterList& params,
    const LINALG::Matrix<Base::numdim_, Base::numnod_>& xrefe,
    const LINALG::Matrix<Base::numdim_, Base::numnod_>& xcurr,
    const LINALG::Matrix<Base::numdim_, Base::numnod_>& nodaldisp,
    const LINALG::Matrix<Base::numdim_, Base::numnod_>& nodalvel,
    const LINALG::Matrix<Base::numdim_, Base::numnod_>& evelnp,
    const LINALG::Matrix<Base::numnod_, 1>& epreaf,
    const LINALG::Matrix<Base::numnod_, 1>* porosity_dof,
    LINALG::Matrix<Base::numdof_, Base::numdof_>& erea_v,
    LINALG::Matrix<Base::numdof_, Base::numdof_>* sub_stiff,
    LINALG::Matrix<Base::numdof_, 1>* sub_force,
    LINALG::Matrix<Base::numdof_, Base::numnod_>& ecoupl_p1,
    LINALG::Matrix<Base::numnod_, numdof_>& estiff_p1,
    LINALG::Matrix<Base::numnod_, 1>& ecoupl_force_p1)
{
  LINALG::Matrix<Base::numdim_, Base::numnod_> N_XYZ;
  // build deformation gradient wrt to material configuration
  // in case of prestressing, build defgrd wrt to last stored configuration
  // CAUTION: defgrd(true): filled with zeros!
  LINALG::Matrix<Base::numdim_, Base::numdim_> defgrd(true);
  LINALG::Matrix<Base::numnod_, 1> shapefct;
  LINALG::Matrix<Base::numdim_, Base::numnod_> deriv;

  LINALG::Matrix<Base::numstr_, 1> fstress(true);

  for (int gp = 0; gp < Base::numgpt_; ++gp)
  {
    // evaluate shape functions and derivatives at integration point
    Base::ComputeShapeFunctionsAndDerivatives(gp, shapefct, deriv, N_XYZ);

    // compute deformation gradient
    Base::ComputeDefGradient(defgrd, N_XYZ, xcurr);

    // inverse deformation gradient F^-1
    LINALG::Matrix<Base::numdim_, Base::numdim_> defgrd_inv(false);
    defgrd_inv.Invert(defgrd);

    // jacobian determinant of transformation between spatial and material space "|dx/dX|"
    double J = 0.0;
    //------linearization of jacobi determinant detF=J w.r.t. structure displacement   dJ/d(us) =
    // dJ/dF : dF/dus = J * F^-T * N,X
    static LINALG::Matrix<1, Base::numdof_> dJ_dus;
    // volume change (used for porosity law). Same as J in nonlinear theory.
    double volchange = 0.0;
    //------linearization of volume change w.r.t. structure displacement
    static LINALG::Matrix<1, Base::numdof_> dvolchange_dus;

    // compute J, the volume change and the respctive linearizations w.r.t. structure displacement
    Base::ComputeJacobianDeterminantVolumeChangeAndLinearizations(
        J, volchange, dJ_dus, dvolchange_dus, defgrd, defgrd_inv, N_XYZ, nodaldisp);

    //----------------------------------------------------
    // pressure at integration point
    double press = shapefct.Dot(epreaf);

    // structure velocity at integration point
    LINALG::Matrix<Base::numdim_, 1> velint(true);

    for (int i = 0; i < Base::numnod_; i++)
      for (int j = 0; j < Base::numdim_; j++) velint(j) += nodalvel(j, i) * shapefct(i);

    // fluid velocity at integration point
    LINALG::Matrix<Base::numdim_, 1> fvelint;
    fvelint.Multiply(evelnp, shapefct);

    // material fluid velocity gradient at integration point
    LINALG::Matrix<Base::numdim_, Base::numdim_> fvelder;
    fvelder.MultiplyNT(evelnp, N_XYZ);

    // pressure gradient at integration point
    LINALG::Matrix<Base::numdim_, 1> Gradp;
    Gradp.Multiply(N_XYZ, epreaf);

    // non-linear B-operator
    LINALG::Matrix<Base::numstr_, Base::numdof_> bop;
    Base::ComputeBOperator(bop, defgrd, N_XYZ);

    // Right Cauchy-Green tensor = F^T * F
    LINALG::Matrix<Base::numdim_, Base::numdim_> cauchygreen;
    cauchygreen.MultiplyTN(defgrd, defgrd);

    // inverse Right Cauchy-Green tensor
    LINALG::Matrix<Base::numdim_, Base::numdim_> C_inv(false);
    C_inv.Invert(cauchygreen);

    //------linearization of material gradient of jacobi determinant GradJ  w.r.t. strucuture
    // displacement d(GradJ)/d(us)
    //---------------------d(GradJ)/dus =  dJ/dus * F^-T . : dF/dX + J * dF^-T/dus : dF/dX + J *
    // F^-T : N_X_X

    // compute some auxiliary matrixes for computation of linearization
    // dF^-T/dus
    LINALG::Matrix<Base::numdim_ * Base::numdim_, Base::numdof_> dFinvTdus(true);
    // F^-T * Grad p
    LINALG::Matrix<Base::numdim_, 1> Finvgradp;
    // dF^-T/dus * Grad p
    LINALG::Matrix<Base::numdim_, Base::numdof_> dFinvdus_gradp(true);
    // dC^-1/dus * Grad p
    LINALG::Matrix<Base::numstr_, Base::numdof_> dCinv_dus(true);

    Base::ComputeAuxiliaryValues(
        N_XYZ, defgrd_inv, C_inv, Gradp, dFinvTdus, Finvgradp, dFinvdus_gradp, dCinv_dus);

    //--------------------------------------------------------------------

    // linearization of porosity w.r.t structure displacement d\phi/d(us) = d\phi/dJ*dJ/d(us)
    LINALG::Matrix<1, Base::numdof_> dphi_dus;
    double porosity = 0.0;

    ComputePorosityAndLinearization(
        params, press, volchange, gp, shapefct, porosity_dof, dvolchange_dus, porosity, dphi_dus);

    double dW_dphi = 0.0;
    double dW_dJ = 0.0;
    double dW_dp = 0.0;
    double W = 0.0;
    Base::struct_mat_->ConstitutiveDerivatives(params, press, volchange, porosity,
        &dW_dp,  // dW_dp not needed
        &dW_dphi, &dW_dJ, nullptr, &W);

    //--------------------------------------------------------

    // **********************evaluate stiffness matrix and force vector+++++++++++++++++++++++++
    if (Base::fluid_mat_->Type() == MAT::PAR::darcy_brinkman)
    {
      Base::FillMatrixAndVectorsBrinkman(gp, J, porosity, fvelder, defgrd_inv, bop, C_inv, dphi_dus,
          dJ_dus, dCinv_dus, dFinvTdus, sub_stiff, sub_force, fstress);
    }

    Base::FillMatrixAndVectors(gp, shapefct, N_XYZ, J, press, porosity, velint, fvelint, fvelder,
        defgrd_inv, bop, C_inv, Finvgradp, dphi_dus, dJ_dus, dCinv_dus, dFinvdus_gradp, dFinvTdus,
        erea_v, sub_stiff, sub_force, fstress);

    // **********************evaluate stiffness matrix and force vector+++++++++++++++++++++++++
    double detJ_w = Base::detJ_[gp] * Base::intpoints_.Weight(gp);  // gpweights[gp];

    const double reacoeff = Base::fluid_mat_->ComputeReactionCoeff();
    {
      for (int k = 0; k < Base::numnod_; k++)
      {
        const double fac = detJ_w * shapefct(k);

        ecoupl_force_p1(k) += fac * W;

        for (int i = 0; i < Base::numnod_; i++)
        {
          for (int j = 0; j < Base::numdim_; j++)
          {
            estiff_p1(k, i * noddof_ + j) += fac * dW_dJ * dJ_dus(i * Base::numdim_ + j);

            ecoupl_p1(i * Base::numdim_ + j, k) +=
                fac * (2 * J * reacoeff * porosity * (velint(j) - fvelint(j)) + J * Finvgradp(j)) *
                shapefct(i);
          }
          estiff_p1(k, i * noddof_ + Base::numdim_) += fac * dW_dphi * shapefct(i);
        }
      }
    }

    if (Base::fluid_mat_->Type() == MAT::PAR::darcy_brinkman)
    {
      double visc = Base::fluid_mat_->Viscosity();
      LINALG::Matrix<Base::numdim_, Base::numdim_> CinvFvel(true);
      LINALG::Matrix<Base::numdim_, Base::numdim_> visctress1(true);
      CinvFvel.Multiply(C_inv, fvelder);
      visctress1.MultiplyNT(CinvFvel, defgrd_inv);
      LINALG::Matrix<Base::numdim_, Base::numdim_> visctress2(visctress1);
      visctress1.UpdateT(1.0, visctress2, 1.0);

      fstress(0) = visctress1(0, 0);
      fstress(1) = visctress1(1, 1);
      fstress(2) = visctress1(0, 1);

      fstress.Scale(detJ_w * visc * J);

      // B^T . C^-1
      LINALG::Matrix<Base::numdof_, 1> fstressb(true);
      fstressb.MultiplyTN(bop, fstress);

      for (int k = 0; k < Base::numnod_; k++)
      {
        const double fac = detJ_w * shapefct(k);
        for (int i = 0; i < Base::numnod_; i++)
          for (int j = 0; j < Base::numdim_; j++)
            ecoupl_p1(i * Base::numdim_ + j, k) +=
                fac * fstressb(i * Base::numdim_ + j) * shapefct(i);
      }
    }
    /* =========================================================================*/
  } /* ==================================================== end of Loop over GP */
  /* =========================================================================*/
}

template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::Wall1_PoroP1<distype>::CouplingPoroelast(
    std::vector<int>& lm,                                // location matrix
    LINALG::Matrix<Base::numdim_, Base::numnod_>& disp,  // current displacements
    LINALG::Matrix<Base::numdim_, Base::numnod_>& vel,   // current velocities
    LINALG::Matrix<Base::numnod_, 1>* porosity,
    LINALG::Matrix<Base::numdim_, Base::numnod_>& evelnp,  // current fluid velocity
    LINALG::Matrix<Base::numnod_, 1>& epreaf,              // current fluid pressure
    LINALG::Matrix<numdof_, (Base::numdim_ + 1) * Base::numnod_>*
        stiffmatrix,  // element stiffness matrix
    LINALG::Matrix<numdof_, (Base::numdim_ + 1) * Base::numnod_>*
        reamatrix,                      // element reactive matrix
    LINALG::Matrix<numdof_, 1>* force,  // element internal force vector
    Teuchos::ParameterList& params)     // algorithmic parameters e.g. time
{
  Base::GetMaterials();

  //=======================================================================

  // update element geometry
  LINALG::Matrix<Base::numdim_, Base::numnod_> xrefe;  // material coord. of element
  LINALG::Matrix<Base::numdim_, Base::numnod_> xcurr;  // current  coord. of element

  DRT::Node** nodes = Base::Nodes();
  for (int i = 0; i < Base::numnod_; ++i)
  {
    const double* x = nodes[i]->X();
    for (int j = 0; j < Base::numdim_; j++)
    {
      xrefe(j, i) = x[j];
      xcurr(j, i) = xrefe(j, i) + disp(j, i);
    }
  }
  // initialize element matrizes
  LINALG::Matrix<Base::numdof_, (Base::numdim_ + 1) * Base::numnod_> ecoupl(true);

  LINALG::Matrix<Base::numnod_, Base::numnod_> ecoupl_p1_p(true);

  /* =========================================================================*/
  /* ================================================= Loop over Gauss Points */
  /* =========================================================================*/
  GaussPointLoopP1OD(
      params, xrefe, xcurr, disp, vel, evelnp, epreaf, porosity, ecoupl_p1_p, ecoupl);

  if (stiffmatrix != nullptr)
  {
    for (int k = 0; k < Base::numnod_; k++)
    {
      for (int l = 0; l < (Base::numdim_ + 1); l++)
      {
        for (int i = 0; i < Base::numnod_; i++)
        {
          for (int j = 0; j < Base::numdim_; j++)
            (*stiffmatrix)(i * noddof_ + j, k * (Base::numdim_ + 1) + l) +=
                ecoupl(i * Base::numdim_ + j, k * (Base::numdim_ + 1) + l);
        }
      }
    }

    for (int ui = 0; ui < Base::numnod_; ++ui)
    {
      for (int ni = 0; ni < Base::numnod_; ++ni)
        (*stiffmatrix)(noddof_ * ui + Base::numdim_, (Base::numdim_ + 1) * ni + Base::numdim_) +=
            ecoupl_p1_p(ui, ni);
    }
  }
}

template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::Wall1_PoroP1<distype>::GaussPointLoopP1OD(Teuchos::ParameterList& params,
    const LINALG::Matrix<Base::numdim_, Base::numnod_>& xrefe,
    const LINALG::Matrix<Base::numdim_, Base::numnod_>& xcurr,
    const LINALG::Matrix<Base::numdim_, Base::numnod_>& nodaldisp,
    const LINALG::Matrix<Base::numdim_, Base::numnod_>& nodalvel,
    const LINALG::Matrix<Base::numdim_, Base::numnod_>& evelnp,
    const LINALG::Matrix<Base::numnod_, 1>& epreaf,
    const LINALG::Matrix<Base::numnod_, 1>* porosity_dof,
    LINALG::Matrix<Base::numnod_, Base::numnod_>& ecoupl_p1,
    LINALG::Matrix<Base::numdof_, (Base::numdim_ + 1) * Base::numnod_>& sub_stiff)
{
  LINALG::Matrix<Base::numdim_, Base::numnod_>
      N_XYZ;  //  first derivatives at gausspoint w.r.t. X,Y,Z
  // build deformation gradient wrt to material configuration
  // in case of prestressing, build defgrd wrt to last stored configuration
  // CAUTION: defgrd(true): filled with zeros!
  LINALG::Matrix<Base::numdim_, Base::numdim_> defgrd(
      true);                                  //  deformation gradiant evaluated at gauss point
  LINALG::Matrix<Base::numnod_, 1> shapefct;  //  shape functions evalulated at gauss point
  LINALG::Matrix<Base::numdim_, Base::numnod_> deriv(
      true);  //  first derivatives at gausspoint w.r.t. r,s,t

  for (int gp = 0; gp < Base::numgpt_; ++gp)
  {
    // evaluate shape functions and derivatives at integration point
    Base::ComputeShapeFunctionsAndDerivatives(gp, shapefct, deriv, N_XYZ);
    // evaluate second derivatives of shape functions at integration point
    // ComputeSecondDerivativesOfShapeFunctions(gp,xrefe,deriv,deriv2,N_XYZ,N_XYZ2);

    // (material) deformation gradient F = d xcurr / d xrefe = xcurr * N_XYZ^T
    Base::ComputeDefGradient(defgrd, N_XYZ, xcurr);

    // inverse deformation gradient F^-1
    LINALG::Matrix<Base::numdim_, Base::numdim_> defgrd_inv(false);
    defgrd_inv.Invert(defgrd);

    // jacobian determinant of transformation between spatial and material space "|dx/dX|"
    double J = 0.0;
    // volume change (used for porosity law). Same as J in nonlinear theory.
    double volchange = 0.0;

    // compute J, the volume change and the respctive linearizations w.r.t. structure displacement
    Base::ComputeJacobianDeterminantVolumeChange(J, volchange, defgrd, N_XYZ, nodaldisp);

    // non-linear B-operator
    LINALG::Matrix<Base::numstr_, Base::numdof_> bop;
    Base::ComputeBOperator(bop, defgrd, N_XYZ);

    // -----------------Right Cauchy-Green tensor = F^T * F
    LINALG::Matrix<Base::numdim_, Base::numdim_> cauchygreen;
    cauchygreen.MultiplyTN(defgrd, defgrd);

    //------------------ inverse Right Cauchy-Green tensor
    LINALG::Matrix<Base::numdim_, Base::numdim_> C_inv(false);
    C_inv.Invert(cauchygreen);

    //---------------- get pressure at integration point
    double press = shapefct.Dot(epreaf);

    //------------------ get material pressure gradient at integration point
    LINALG::Matrix<Base::numdim_, 1> Gradp;
    Gradp.Multiply(N_XYZ, epreaf);

    //--------------------- get fluid velocity at integration point
    LINALG::Matrix<Base::numdim_, 1> fvelint;
    fvelint.Multiply(evelnp, shapefct);

    //---------------- material fluid velocity gradient at integration point
    LINALG::Matrix<Base::numdim_, Base::numdim_> fvelder;
    fvelder.MultiplyNT(evelnp, N_XYZ);

    //---------------- structure velocity at integration point
    LINALG::Matrix<Base::numdim_, 1> velint(true);
    for (int i = 0; i < Base::numnod_; i++)
      for (int j = 0; j < Base::numdim_; j++) velint(j) += nodalvel(j, i) * shapefct(i);

    //**************************************************+auxilary variables for computing the
    // porosity and linearization
    double dphi_dp = 0.0;
    double porosity = 0.0;

    ComputePorosityAndLinearizationOD(
        params, press, volchange, gp, shapefct, porosity_dof, porosity, dphi_dp);

    // **********************evaluate stiffness matrix and force vector+++++++++++++++++++++++++

    Base::FillMatrixAndVectorsOD(gp, shapefct, N_XYZ, J, porosity, dphi_dp, velint, fvelint,
        defgrd_inv, Gradp, bop, C_inv, sub_stiff);

    if (Base::fluid_mat_->Type() == MAT::PAR::darcy_brinkman)
    {
      Base::FillMatrixAndVectorsBrinkmanOD(
          gp, shapefct, N_XYZ, J, porosity, dphi_dp, fvelder, defgrd_inv, bop, C_inv, sub_stiff);
    }  // darcy-brinkman

    double dW_dp = 0.0;
    Base::struct_mat_->ConstitutiveDerivatives(params, press, J, porosity, &dW_dp,
        nullptr,  // not needed
        nullptr,  // not needed
        nullptr,  // not needed
        nullptr   // not needed
    );
    //--------------------------------------------------------

    // **********************evaluate stiffness matrix and force vector+++++++++++++++++++++++++
    double detJ_w = Base::detJ_[gp] * Base::intpoints_.Weight(gp);

    for (int k = 0; k < Base::numnod_; k++)
    {
      const double fac = detJ_w * shapefct(k);

      for (int i = 0; i < Base::numnod_; i++) ecoupl_p1(k, i) += fac * dW_dp * shapefct(i);
    }

    /* =========================================================================*/
  } /* ==================================================== end of Loop over GP */
  /* =========================================================================*/
}

template <DRT::Element::DiscretizationType distype>
int DRT::ELEMENTS::Wall1_PoroP1<distype>::EvaluateNeumann(Teuchos::ParameterList& params,
    DRT::Discretization& discretization, DRT::Condition& condition, std::vector<int>& lm,
    Epetra_SerialDenseVector& elevec1, Epetra_SerialDenseMatrix* elemat1)
{
  LINALG::Matrix<Base::numdim_, Base::numnod_> disp(true);
  LINALG::Matrix<Base::numnod_, 1> myporosity(true);
  Base::ExtractValuesFromGlobalVector(discretization, 0, lm, &disp, &myporosity, "displacement");

  // find out whether we will use a time curve
  const double time = params.get("total time", -1.0);

  /*----------------------------------------------------- geometry update */
  // update element geometry
  LINALG::Matrix<Base::numdim_, Base::numnod_> xrefe;  // material coord. of element
  LINALG::Matrix<Base::numdim_, Base::numnod_> xcurr;  // current  coord. of element

  DRT::Node** nodes = Base::Nodes();
  for (int i = 0; i < Base::numnod_; ++i)
  {
    const double* x = nodes[i]->X();
    for (int j = 0; j < Base::numdim_; j++)
    {
      xrefe(j, i) = x[j];
      xcurr(j, i) = xrefe(j, i) + disp(j, i);
    }
  }


  // get values and switches from the condition
  const auto* onoff = condition.Get<std::vector<int>>("onoff");
  const auto* val = condition.Get<std::vector<double>>("val");
  const auto* funct = condition.Get<std::vector<int>>("funct");


  LINALG::Matrix<Base::numdim_, Base::numnod_> N_XYZ;
  // build deformation gradient wrt to material configuration
  // in case of prestressing, build defgrd wrt to last stored configuration
  // CAUTION: defgrd(true): filled with zeros!
  LINALG::Matrix<Base::numdim_, Base::numdim_> defgrd(true);
  LINALG::Matrix<Base::numnod_, 1> shapefcts;
  LINALG::Matrix<Base::numdim_, Base::numnod_> deriv;

  LINALG::Matrix<Base::numstr_, 1> fstress(true);

  for (int gp = 0; gp < Base::numgpt_; ++gp)
  {
    // evaluate shape functions and derivatives at integration point
    Base::ComputeShapeFunctionsAndDerivatives(gp, shapefcts, deriv, N_XYZ);

    // jacobian determinant of transformation between spatial and material space "|dx/dX|"
    Base::ComputeJacobianDeterminant(gp, xcurr, deriv);

    /*------------------------------------ integration factor  -------*/
    double fac = Base::detJ_[gp] * Base::intpoints_.Weight(gp);

    // load vector ar
    std::array<double, Base::numdim_> ar = {0.0, 0.0};
    // loop the dofs of a node
    for (int i = 0; i < Base::numdim_; ++i)
    {
      if ((*onoff)[i])
      {
        // factor given by spatial function
        const int functnum = (funct) ? (*funct)[i] : -1;
        double functfac = 1.0;
        if (functnum > 0)
        {
          // calculate reference position of GP
          LINALG::Matrix<1, Base::numdim_> gp_coord;
          gp_coord.MultiplyTT(shapefcts, xrefe);

          // write coordinates in another datatype
          std::array<double, 3> gp_coord2;  // the position vector has to be given in 3D!!!
          for (int k = 0; k < Base::numdim_; k++) gp_coord2[k] = gp_coord(0, k);
          for (int k = Base::numdim_; k < 3;
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
    for (int node = 0; node < Base::numnod_; ++node)
      for (int dim = 0; dim < Base::numdim_; ++dim)
        elevec1[node * noddof_ + dim] += shapefcts(node) * ar[dim];
  }

  return 0;
}

template class DRT::ELEMENTS::Wall1_PoroP1<DRT::Element::tri3>;
template class DRT::ELEMENTS::Wall1_PoroP1<DRT::Element::quad4>;
template class DRT::ELEMENTS::Wall1_PoroP1<DRT::Element::quad9>;
