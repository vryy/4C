/*----------------------------------------------------------------------*/
/*! \file

 \brief evaluate methods for porous media using the p1 approach (mixed approach)

 \level 2

 *----------------------------------------------------------------------*/

#include "baci_lib_globalproblem.H"
#include "baci_mat_fluidporo.H"
#include "baci_mat_structporo.H"
#include "baci_so3_poro_p1.H"
#include "baci_so3_poro_p1_eletypes.H"


template <class so3_ele, DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::So3_Poro_P1<so3_ele, distype>::ComputePorosityAndLinearization(
    Teuchos::ParameterList& params, const double& press, const double& J, const int& gp,
    const CORE::LINALG::Matrix<Base::numnod_, 1>& shapfct,
    const CORE::LINALG::Matrix<Base::numnod_, 1>* myporosity,
    const CORE::LINALG::Matrix<1, Base::numdof_>& dJ_dus, double& porosity,
    CORE::LINALG::Matrix<1, Base::numdof_>& dphi_dus)
{
  if (myporosity == nullptr)
    dserror("no porosity values given!");
  else
    porosity = shapfct.Dot(*myporosity);

  dphi_dus.PutScalar(0.0);
}

template <class so3_ele, DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::So3_Poro_P1<so3_ele, distype>::ComputePorosityAndLinearizationOD(
    Teuchos::ParameterList& params, const double& press, const double& J, const int& gp,
    const CORE::LINALG::Matrix<Base::numnod_, 1>& shapfct,
    const CORE::LINALG::Matrix<Base::numnod_, 1>* myporosity, double& porosity, double& dphi_dp)
{
  dphi_dp = 0.0;

  if (myporosity == nullptr)
    dserror("no porosity values given!");
  else
    porosity = shapfct.Dot(*myporosity);
}

template <class so3_ele, DRT::Element::DiscretizationType distype>
int DRT::ELEMENTS::So3_Poro_P1<so3_ele, distype>::Evaluate(Teuchos::ParameterList& params,
    DRT::Discretization& discretization, DRT::Element::LocationArray& la,
    CORE::LINALG::SerialDenseMatrix& elemat1_epetra,
    CORE::LINALG::SerialDenseMatrix& elemat2_epetra,
    CORE::LINALG::SerialDenseVector& elevec1_epetra,
    CORE::LINALG::SerialDenseVector& elevec2_epetra,
    CORE::LINALG::SerialDenseVector& elevec3_epetra)
{
  // start with "none"
  typename Base::ActionType act = Base::none;

  // get the required action
  std::string action = params.get<std::string>("action", "none");
  if (action == "none")
    dserror("No action supplied");
  else if (action == "struct_poro_calc_fluidcoupling")
    act = Base::calc_struct_multidofsetcoupling;
  else if (action == "interpolate_porosity_to_given_point")
    act = Base::interpolate_porosity_to_given_point;

  // what should the element do
  switch (act)
  {
    //==================================================================================
    // off diagonal terms in stiffness matrix for monolithic coupling
    case Base::calc_struct_multidofsetcoupling:
    {
      // in some cases we need to write/change some data before evaluating
      PreEvaluate(params, discretization, la);

      MyEvaluate(params, discretization, la, elemat1_epetra, elemat2_epetra, elevec1_epetra,
          elevec2_epetra, elevec3_epetra);
    }
    break;
    case Base::interpolate_porosity_to_given_point:
    {
      // given point
      std::vector<CORE::LINALG::Matrix<Base::numdim_, 1>> xsi(1);
      (xsi[0])(0, 0) = elevec2_epetra(0);
      (xsi[0])(1, 0) = elevec2_epetra(1);
      (xsi[0])(2, 0) = elevec2_epetra(2);

      //  evalulate shape functions at given point
      CORE::LINALG::Matrix<Base::numnod_, 1> shapefct;
      CORE::DRT::UTILS::shape_function<distype>(xsi[0], shapefct);

      CORE::LINALG::Matrix<Base::numdim_, Base::numnod_> mydisp(true);
      CORE::LINALG::Matrix<Base::numnod_, 1> myporosity(true);
      Base::ExtractValuesFromGlobalVector(
          discretization, 0, la[0].lm_, nullptr, &myporosity, "displacement");

      elevec1_epetra(0) = shapefct.Dot(myporosity);
    }
    break;
    //==================================================================================
    default:
    {
      // in some cases we need to write/change some data before evaluating
      PreEvaluate(params, discretization, la);

      CORE::LINALG::SerialDenseMatrix elemat1_sub;
      CORE::LINALG::SerialDenseMatrix elemat2_sub;
      CORE::LINALG::SerialDenseVector elevec1_sub;
      CORE::LINALG::SerialDenseVector elevec2_sub;
      CORE::LINALG::SerialDenseVector elevec3_sub;

      if (elemat1_epetra.values()) elemat1_sub.shape(Base::numdof_, Base::numdof_);
      if (elemat2_epetra.values()) elemat2_sub.shape(Base::numdof_, Base::numdof_);
      if (elevec1_epetra.values()) elevec1_sub.resize(Base::numdof_);
      if (elevec2_epetra.values()) elevec2_sub.resize(Base::numdof_);
      if (elevec3_epetra.values()) elevec3_sub.resize(Base::numdof_);

      std::vector<int> lm_sub;
      for (int i = 0; i < Base::numnod_; i++)
        for (int j = 0; j < Base::numdim_; j++) lm_sub.push_back(la[0].lm_[i * noddof_ + j]);

      // evaluate parent solid element
      so3_ele::Evaluate(params, discretization, lm_sub, elemat1_sub, elemat2_sub, elevec1_sub,
          elevec2_sub, elevec3_sub);

      if (elemat1_epetra.values())
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

      if (elemat2_epetra.values())
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

      if (elevec1_epetra.values())
      {
        for (int i = 0; i < Base::numnod_; i++)
          for (int j = 0; j < Base::numdim_; j++)
            elevec1_epetra(i * noddof_ + j) = elevec1_sub(i * Base::noddof_ + j);
      }

      if (elevec2_epetra.values())
      {
        for (int i = 0; i < Base::numnod_; i++)
          for (int j = 0; j < Base::numdim_; j++)
            elevec2_epetra(i * noddof_ + j) = elevec2_sub(i * Base::noddof_ + j);
      }

      if (elevec3_epetra.values())
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

template <class so3_ele, DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::So3_Poro_P1<so3_ele, distype>::PreEvaluate(Teuchos::ParameterList& params,
    DRT::Discretization& discretization, DRT::Element::LocationArray& la)
{
  Base::PreEvaluate(params, discretization, la);

  if (discretization.HasState(0, "displacement") and (not is_init_porosity_))
  {
    init_porosity_ = Teuchos::rcp(new CORE::LINALG::Matrix<Base::numnod_, 1>(true));
    Base::ExtractValuesFromGlobalVector(
        discretization, 0, la[0].lm_, nullptr, &(*init_porosity_), "displacement");
    is_init_porosity_ = true;
  }
}

template <class so3_ele, DRT::Element::DiscretizationType distype>
int DRT::ELEMENTS::So3_Poro_P1<so3_ele, distype>::MyEvaluate(Teuchos::ParameterList& params,
    DRT::Discretization& discretization, DRT::Element::LocationArray& la,
    CORE::LINALG::SerialDenseMatrix& elemat1_epetra,
    CORE::LINALG::SerialDenseMatrix& elemat2_epetra,
    CORE::LINALG::SerialDenseVector& elevec1_epetra,
    CORE::LINALG::SerialDenseVector& elevec2_epetra,
    CORE::LINALG::SerialDenseVector& elevec3_epetra)
{
  // start with "none"
  typename Base::ActionType act = Base::none;

  // get the required action
  std::string action = params.get<std::string>("action", "none");
  if (action == "none")
    dserror("No action supplied");
  else if (action == "calc_struct_internalforce")
    act = Base::calc_struct_internalforce;
  else if (action == "calc_struct_nlnstiff")
    act = Base::calc_struct_nlnstiff;
  else if (action == "calc_struct_nlnstiffmass")
    act = Base::calc_struct_nlnstiffmass;
  else if (action == "struct_poro_calc_fluidcoupling")
    act = Base::calc_struct_multidofsetcoupling;

  // what should the element do
  switch (act)
  {
    //==================================================================================
    // nonlinear stiffness, damping and internal force vector for poroelasticity
    case Base::calc_struct_nlnstiff:
    case Base::calc_struct_nlnstiffmass:
    {
      if (la.Size() > 1)
      {
        // stiffness
        CORE::LINALG::Matrix<numdof_, numdof_> elemat1(elemat1_epetra.values(), true);
        // damping
        CORE::LINALG::Matrix<numdof_, numdof_> elemat2(elemat2_epetra.values(), true);
        // internal force vector
        CORE::LINALG::Matrix<numdof_, 1> elevec1(elevec1_epetra.values(), true);
        // elevec2+3 are not used anyway

        // build the location vector only for the structure field
        std::vector<int> lm = la[0].lm_;

        CORE::LINALG::Matrix<Base::numdim_, Base::numnod_> mydisp(true);
        CORE::LINALG::Matrix<Base::numnod_, 1> myporosity(true);
        Base::ExtractValuesFromGlobalVector(
            discretization, 0, la[0].lm_, &mydisp, &myporosity, "displacement");

        CORE::LINALG::Matrix<numdof_, numdof_>* matptr = nullptr;
        if (elemat1.IsInitialized()) matptr = &elemat1;

        enum INPAR::STR::DampKind damping =
            params.get<enum INPAR::STR::DampKind>("damping", INPAR::STR::damp_none);
        CORE::LINALG::Matrix<numdof_, numdof_>* matptr2 = nullptr;
        if (elemat2.IsInitialized() and (damping == INPAR::STR::damp_material)) matptr2 = &elemat2;

        // need current fluid state,
        // call the fluid discretization: fluid equates 2nd dofset
        // disassemble velocities and pressures

        CORE::LINALG::Matrix<Base::numdim_, Base::numnod_> myvel(true);

        CORE::LINALG::Matrix<Base::numdim_, Base::numnod_> myfluidvel(true);
        CORE::LINALG::Matrix<Base::numnod_, 1> myepreaf(true);

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
    case Base::calc_struct_multidofsetcoupling:
    {
      // stiffness
      CORE::LINALG::Matrix<numdof_, (Base::numdim_ + 1) * Base::numnod_> elemat1(
          elemat1_epetra.values(), true);
      // elemat2,elevec1-3 are not used anyway

      // build the location vector only for the structure field
      std::vector<int> lm = la[0].lm_;

      CORE::LINALG::Matrix<numdof_, (Base::numdim_ + 1)* Base::numnod_>* matptr = nullptr;
      if (elemat1.IsInitialized()) matptr = &elemat1;

      // need current fluid state,
      // call the fluid discretization: fluid equates 2nd dofset
      // disassemble velocities and pressures
      if (discretization.HasState(1, "fluidvel"))
      {
        CORE::LINALG::Matrix<Base::numdim_, Base::numnod_> myvel(true);
        CORE::LINALG::Matrix<Base::numdim_, Base::numnod_> myfluidvel(true);
        CORE::LINALG::Matrix<Base::numnod_, 1> myepreaf(true);

        CORE::LINALG::Matrix<Base::numdim_, Base::numnod_> mydisp(true);
        CORE::LINALG::Matrix<Base::numnod_, 1> myporosity(true);
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

        CouplingPoroelast(
            lm, mydisp, myvel, &myporosity, myfluidvel, myepreaf, matptr, nullptr, nullptr, params);
      }
    }
    break;

    //==================================================================================
    // nonlinear stiffness and internal force vector for poroelasticity
    case Base::calc_struct_internalforce:
    {
      // internal force vector
      CORE::LINALG::Matrix<numdof_, 1> elevec1(elevec1_epetra.values(), true);
      // elemat1+2,elevec2+3 are not used anyway

      // build the location vector only for the structure field
      std::vector<int> lm = la[0].lm_;

      CORE::LINALG::Matrix<Base::numdim_, Base::numnod_> mydisp(true);
      CORE::LINALG::Matrix<Base::numnod_, 1> myporosity(true);
      Base::ExtractValuesFromGlobalVector(
          discretization, 0, la[0].lm_, &mydisp, &myporosity, "displacement");

      CORE::LINALG::Matrix<Base::numdim_, Base::numnod_> myvel(true);

      CORE::LINALG::Matrix<Base::numdim_, Base::numnod_> myfluidvel(true);
      CORE::LINALG::Matrix<Base::numnod_, 1> myepreaf(true);

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

template <class so3_ele, DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::So3_Poro_P1<so3_ele, distype>::InitElement()
{
  // initialize base element
  Base::InitElement();
}

template <class so3_ele, DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::So3_Poro_P1<so3_ele, distype>::NonlinearStiffnessPoroelast(std::vector<int>& lm,
    CORE::LINALG::Matrix<Base::numdim_, Base::numnod_>& disp,
    CORE::LINALG::Matrix<Base::numdim_, Base::numnod_>& vel,
    CORE::LINALG::Matrix<Base::numnod_, 1>* porosity_dof,
    CORE::LINALG::Matrix<Base::numdim_, Base::numnod_>& evelnp,
    CORE::LINALG::Matrix<Base::numnod_, 1>& epreaf,
    CORE::LINALG::Matrix<numdof_, numdof_>* stiffmatrix,
    CORE::LINALG::Matrix<numdof_, numdof_>* reamatrix, CORE::LINALG::Matrix<numdof_, 1>* force,
    Teuchos::ParameterList& params)
{
  Base::GetMaterials();

  // update element geometry
  CORE::LINALG::Matrix<Base::numdim_, Base::numnod_> xrefe;  // material coord. of element
  CORE::LINALG::Matrix<Base::numdim_, Base::numnod_> xcurr;  // current  coord. of element

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
  CORE::LINALG::Matrix<Base::numdof_, Base::numdof_> erea_v(true);
  CORE::LINALG::Matrix<Base::numdof_, Base::numdof_> sub_stiff(true);
  CORE::LINALG::Matrix<Base::numdof_, 1> sub_force(true);

  CORE::LINALG::Matrix<Base::numdof_, Base::numnod_> ecoupl_p1(true);
  CORE::LINALG::Matrix<Base::numnod_, numdof_> estiff_p1(true);
  CORE::LINALG::Matrix<Base::numnod_, 1> ecoupl_force_p1(true);

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
    else
    {
      const double dt = params.get<double>("delta time");
      // if the reaction part is not supposed to be computed separately, we add it to the stiffness
      //(this is not the best way to do it, but it only happens once during initialization)
      for (int k = 0; k < Base::numnod_; k++)
      {
        for (int l = 0; l < Base::numdim_; l++)
        {
          for (int i = 0; i < Base::numnod_; i++)
          {
            for (int j = 0; j < Base::numdim_; j++)
              (*stiffmatrix)(i * noddof_ + j, k * noddof_ + l) +=
                  erea_v(i * Base::numdim_ + j, k * Base::numdim_ + l) / dt;
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

template <class so3_ele, DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::So3_Poro_P1<so3_ele, distype>::GaussPointLoopP1(Teuchos::ParameterList& params,
    const CORE::LINALG::Matrix<Base::numdim_, Base::numnod_>& xrefe,
    const CORE::LINALG::Matrix<Base::numdim_, Base::numnod_>& xcurr,
    const CORE::LINALG::Matrix<Base::numdim_, Base::numnod_>& nodaldisp,
    const CORE::LINALG::Matrix<Base::numdim_, Base::numnod_>& nodalvel,
    const CORE::LINALG::Matrix<Base::numdim_, Base::numnod_>& evelnp,
    const CORE::LINALG::Matrix<Base::numnod_, 1>& epreaf,
    const CORE::LINALG::Matrix<Base::numnod_, 1>* porosity_dof,
    CORE::LINALG::Matrix<Base::numdof_, Base::numdof_>& erea_v,
    CORE::LINALG::Matrix<Base::numdof_, Base::numdof_>* sub_stiff,
    CORE::LINALG::Matrix<Base::numdof_, 1>* sub_force,
    CORE::LINALG::Matrix<Base::numdof_, Base::numnod_>& ecoupl_p1,
    CORE::LINALG::Matrix<Base::numnod_, numdof_>& estiff_p1,
    CORE::LINALG::Matrix<Base::numnod_, 1>& ecoupl_force_p1)
{
  CORE::LINALG::Matrix<Base::numdim_, Base::numnod_> N_XYZ;
  // build deformation gradient wrt to material configuration
  // in case of prestressing, build defgrd wrt to last stored configuration
  // CAUTION: defgrd(true): filled with zeros!
  CORE::LINALG::Matrix<Base::numdim_, Base::numdim_> defgrd(true);
  CORE::LINALG::Matrix<Base::numnod_, 1> shapefct;
  CORE::LINALG::Matrix<Base::numdim_, Base::numnod_> deriv;

  CORE::LINALG::Matrix<Base::numstr_, 1> fstress(true);

  for (int gp = 0; gp < Base::numgpt_; ++gp)
  {
    // evaluate shape functions and derivatives at integration point
    Base::ComputeShapeFunctionsAndDerivatives(gp, shapefct, deriv, N_XYZ);

    // (material) deformation gradient F = d xcurr / d xrefe = xcurr * N_XYZ^T
    Base::ComputeDefGradient(defgrd, N_XYZ, xcurr);

    // inverse deformation gradient F^-1
    static CORE::LINALG::Matrix<Base::numdim_, Base::numdim_> defgrd_inv(false);
    defgrd_inv.Invert(defgrd);

    // jacobian determinant of transformation between spatial and material space "|dx/dX|"
    double J = 0.0;
    //------linearization of jacobi determinant detF=J w.r.t. structure displacement   dJ/d(us) =
    // dJ/dF : dF/dus = J * F^-T * N,X
    static CORE::LINALG::Matrix<1, Base::numdof_> dJ_dus;
    // volume change (used for porosity law). Same as J in nonlinear theory.
    double volchange = 0.0;
    //------linearization of volume change w.r.t. structure displacement
    static CORE::LINALG::Matrix<1, Base::numdof_> dvolchange_dus;

    // compute J, the volume change and the respctive linearizations w.r.t. structure displacement
    Base::ComputeJacobianDeterminantVolumeChangeAndLinearizations(
        J, volchange, dJ_dus, dvolchange_dus, defgrd, defgrd_inv, N_XYZ, nodaldisp);

    // pressure at integration point
    double press = shapefct.Dot(epreaf);

    // structure velocity at integration point
    CORE::LINALG::Matrix<Base::numdim_, 1> velint(true);

    for (int i = 0; i < Base::numnod_; i++)
      for (int j = 0; j < Base::numdim_; j++) velint(j) += nodalvel(j, i) * shapefct(i);

    // fluid velocity at integration point
    CORE::LINALG::Matrix<Base::numdim_, 1> fvelint;
    fvelint.Multiply(evelnp, shapefct);

    // material fluid velocity gradient at integration point
    CORE::LINALG::Matrix<Base::numdim_, Base::numdim_> fvelder;
    fvelder.MultiplyNT(evelnp, N_XYZ);

    // pressure gradient at integration point
    CORE::LINALG::Matrix<Base::numdim_, 1> Gradp;
    Gradp.Multiply(N_XYZ, epreaf);

    // non-linear B-operator
    CORE::LINALG::Matrix<Base::numstr_, Base::numdof_> bop;
    Base::ComputeBOperator(bop, defgrd, N_XYZ);

    // Right Cauchy-Green tensor = F^T * F
    CORE::LINALG::Matrix<Base::numdim_, Base::numdim_> cauchygreen;
    cauchygreen.MultiplyTN(defgrd, defgrd);

    // inverse Right Cauchy-Green tensor
    CORE::LINALG::Matrix<Base::numdim_, Base::numdim_> C_inv(false);
    C_inv.Invert(cauchygreen);

    //------linearization of material gradient of jacobi determinant GradJ  w.r.t. strucuture
    // displacement d(GradJ)/d(us)
    //---------------------d(GradJ)/dus =  dJ/dus * F^-T . : dF/dX + J * dF^-T/dus : dF/dX + J *
    // F^-T : N_X_X

    // compute some auxiliary matrixes for computation of linearization
    // dF^-T/dus
    CORE::LINALG::Matrix<Base::numdim_ * Base::numdim_, Base::numdof_> dFinvTdus(true);
    // F^-T * Grad p
    CORE::LINALG::Matrix<Base::numdim_, 1> Finvgradp;
    // dF^-T/dus * Grad p
    CORE::LINALG::Matrix<Base::numdim_, Base::numdof_> dFinvdus_gradp(true);
    // dC^-1/dus * Grad p
    CORE::LINALG::Matrix<Base::numstr_, Base::numdof_> dCinv_dus(true);

    Base::ComputeAuxiliaryValues(
        N_XYZ, defgrd_inv, C_inv, Gradp, dFinvTdus, Finvgradp, dFinvdus_gradp, dCinv_dus);

    //--------------------------------------------------------------------

    // linearization of porosity w.r.t structure displacement d\phi/d(us) = d\phi/dJ*dJ/d(us)
    CORE::LINALG::Matrix<1, Base::numdof_> dphi_dus;
    double porosity = 0.0;

    ComputePorosityAndLinearization(
        params, press, volchange, gp, shapefct, porosity_dof, dvolchange_dus, porosity, dphi_dus);

    double dW_dphi = 0.0;
    double dW_dJ = 0.0;
    double dW_dp = 0.0;
    double W = 0.0;

    if (init_porosity_ == Teuchos::null)
      dserror("Failed to create vector of nodal intial porosity");
    double init_porosity = shapefct.Dot(*init_porosity_);
    Base::struct_mat_->ConstitutiveDerivatives(params, press, volchange, porosity, init_porosity,
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
    // if (force != nullptr or stiffmatrix != nullptr or reamatrix != nullptr )
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
      CORE::LINALG::Matrix<Base::numdim_, Base::numdim_> CinvFvel;
      CORE::LINALG::Matrix<Base::numdim_, Base::numdim_> visctress1;
      CinvFvel.Multiply(C_inv, fvelder);
      visctress1.MultiplyNT(CinvFvel, defgrd_inv);
      CORE::LINALG::Matrix<Base::numdim_, Base::numdim_> visctress2(visctress1);
      visctress1.UpdateT(1.0, visctress2, 1.0);

      fstress(0) = visctress1(0, 0);
      fstress(1) = visctress1(1, 1);
      fstress(2) = visctress1(2, 2);
      fstress(3) = visctress1(0, 1);
      fstress(4) = visctress1(1, 2);
      fstress(5) = visctress1(2, 0);

      fstress.Scale(detJ_w * visc * J);

      // B^T . C^-1
      CORE::LINALG::Matrix<Base::numdof_, 1> fstressb(true);
      fstressb.MultiplyTN(bop, fstress);

      for (int k = 0; k < Base::numnod_; k++)
      {
        const double fac = detJ_w * shapefct(k);
        for (int i = 0; i < Base::numnod_; i++)
        {
          for (int j = 0; j < Base::numdim_; j++)
            ecoupl_p1(i * Base::numdim_ + j, k) +=
                fac * fstressb(i * Base::numdim_ + j) * shapefct(i);
        }
      }
    }
  }
}

template <class so3_ele, DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::So3_Poro_P1<so3_ele, distype>::CouplingPoroelast(std::vector<int>& lm,
    CORE::LINALG::Matrix<Base::numdim_, Base::numnod_>& disp,
    CORE::LINALG::Matrix<Base::numdim_, Base::numnod_>& vel,
    CORE::LINALG::Matrix<Base::numnod_, 1>* porosity,
    CORE::LINALG::Matrix<Base::numdim_, Base::numnod_>& evelnp,
    CORE::LINALG::Matrix<Base::numnod_, 1>& epreaf,
    CORE::LINALG::Matrix<numdof_, (Base::numdim_ + 1) * Base::numnod_>* stiffmatrix,
    [[maybe_unused]] CORE::LINALG::Matrix<numdof_, (Base::numdim_ + 1) * Base::numnod_>* reamatrix,
    [[maybe_unused]] CORE::LINALG::Matrix<numdof_, 1>* force, Teuchos::ParameterList& params)
{
  //=============================get parameters

  Base::GetMaterials();

  //=======================================================================

  // update element geometry
  CORE::LINALG::Matrix<Base::numdim_, Base::numnod_> xrefe;  // material coord. of element
  CORE::LINALG::Matrix<Base::numdim_, Base::numnod_> xcurr;  // current  coord. of element

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
  CORE::LINALG::Matrix<Base::numdof_, (Base::numdim_ + 1) * Base::numnod_> ecoupl(true);

  CORE::LINALG::Matrix<Base::numnod_, Base::numnod_> ecoupl_p1_p(true);

  GaussPointLoopP1OD(
      params, xrefe, xcurr, disp, vel, evelnp, epreaf, porosity, ecoupl_p1_p, &ecoupl);

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

template <class so3_ele, DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::So3_Poro_P1<so3_ele, distype>::GaussPointLoopP1OD(
    Teuchos::ParameterList& params, const CORE::LINALG::Matrix<Base::numdim_, Base::numnod_>& xrefe,
    const CORE::LINALG::Matrix<Base::numdim_, Base::numnod_>& xcurr,
    const CORE::LINALG::Matrix<Base::numdim_, Base::numnod_>& nodaldisp,
    const CORE::LINALG::Matrix<Base::numdim_, Base::numnod_>& nodalvel,
    const CORE::LINALG::Matrix<Base::numdim_, Base::numnod_>& evelnp,
    const CORE::LINALG::Matrix<Base::numnod_, 1>& epreaf,
    const CORE::LINALG::Matrix<Base::numnod_, 1>* porosity_dof,
    CORE::LINALG::Matrix<Base::numnod_, Base::numnod_>& ecoupl_p1,
    CORE::LINALG::Matrix<Base::numdof_, (Base::numdim_ + 1) * Base::numnod_>* sub_stiff)
{
  CORE::LINALG::Matrix<Base::numdim_, Base::numnod_>
      N_XYZ;  //  first derivatives at gausspoint w.r.t. X, Y,Z
  // build deformation gradient wrt to material configuration
  // in case of prestressing, build defgrd wrt to last stored configuration
  // CAUTION: defgrd(true): filled with zeros!
  CORE::LINALG::Matrix<Base::numdim_, Base::numdim_> defgrd(
      true);  //  deformation gradiant evaluated at gauss point
  CORE::LINALG::Matrix<Base::numnod_, 1> shapefct;  //  shape functions evalulated at gauss point
  CORE::LINALG::Matrix<Base::numdim_, Base::numnod_> deriv(
      true);  //  first derivatives at gausspoint w.r.t. r,s,t

  for (int gp = 0; gp < Base::numgpt_; ++gp)
  {
    // evaluate shape functions and derivatives at integration point
    Base::ComputeShapeFunctionsAndDerivatives(gp, shapefct, deriv, N_XYZ);

    // (material) deformation gradient F = d xcurr / d xrefe = xcurr * N_XYZ^T
    defgrd.MultiplyNT(xcurr, N_XYZ);  //  (6.17)

    // inverse deformation gradient F^-1
    static CORE::LINALG::Matrix<Base::numdim_, Base::numdim_> defgrd_inv(false);
    defgrd_inv.Invert(defgrd);

    // jacobian determinant of transformation between spatial and material space "|dx/dX|"
    double J = 0.0;
    //------linearization of jacobi determinant detF=J w.r.t. structure displacement   dJ/d(us) =
    // dJ/dF : dF/dus = J * F^-T * N,X
    static CORE::LINALG::Matrix<1, numdof_> dJ_dus;
    // volume change (used for porosity law). Same as J in nonlinear theory.
    double volchange = 0.0;
    //------linearization of volume change w.r.t. structure displacement
    static CORE::LINALG::Matrix<1, numdof_> dvolchange_dus;

    // compute J, the volume change and the respctive linearizations w.r.t. structure displacement
    Base::ComputeJacobianDeterminantVolumeChange(J, volchange, defgrd, N_XYZ, nodaldisp);

    // non-linear B-operator
    CORE::LINALG::Matrix<Base::numstr_, Base::numdof_> bop;
    Base::ComputeBOperator(bop, defgrd, N_XYZ);

    // -----------------Right Cauchy-Green tensor = F^T * F
    CORE::LINALG::Matrix<Base::numdim_, Base::numdim_> cauchygreen;
    cauchygreen.MultiplyTN(defgrd, defgrd);

    //------------------ inverse Right Cauchy-Green tensor
    CORE::LINALG::Matrix<Base::numdim_, Base::numdim_> C_inv(false);
    C_inv.Invert(cauchygreen);

    //---------------- get pressure at integration point
    double press = shapefct.Dot(epreaf);

    //------------------ get material pressure gradient at integration point
    CORE::LINALG::Matrix<Base::numdim_, 1> Gradp;
    Gradp.Multiply(N_XYZ, epreaf);

    //--------------------- get fluid velocity at integration point
    CORE::LINALG::Matrix<Base::numdim_, 1> fvelint;
    fvelint.Multiply(evelnp, shapefct);

    //----------------------- material fluid velocity gradient at integration point
    CORE::LINALG::Matrix<Base::numdim_, Base::numdim_> fvelder;
    fvelder.MultiplyNT(evelnp, N_XYZ);

    //---------------- structure velocity at integration point
    CORE::LINALG::Matrix<Base::numdim_, 1> velint;
    velint.Multiply(nodalvel, shapefct);

    //**************************************************+auxilary variables for computing the
    // porosity and linearization
    double dphi_dp = 0.0;
    double porosity = 0.0;

    ComputePorosityAndLinearizationOD(
        params, press, volchange, gp, shapefct, porosity_dof, porosity, dphi_dp);

    // **********************evaluate stiffness matrix and force vector**********************

    Base::FillMatrixAndVectorsOD(gp, shapefct, N_XYZ, J, porosity, dphi_dp, velint, fvelint,
        defgrd_inv, Gradp, bop, C_inv, sub_stiff);

    if (Base::fluid_mat_->Type() == MAT::PAR::darcy_brinkman)
    {
      Base::FillMatrixAndVectorsBrinkmanOD(
          gp, shapefct, N_XYZ, J, porosity, dphi_dp, fvelder, defgrd_inv, bop, C_inv, sub_stiff);
    }

    double dW_dp = 0.0;
    double init_porosity = shapefct.Dot(*init_porosity_);
    Base::struct_mat_->ConstitutiveDerivatives(params, press, J, porosity, init_porosity, &dW_dp,
        nullptr,  // not needed
        nullptr,  // not needed
        nullptr,  // not needed
        nullptr   // not needed
    );
    //--------------------------------------------------------

    // **********************evaluate stiffness matrix and force vector**********************
    double detJ_w = Base::detJ_[gp] * Base::intpoints_.Weight(gp);

    for (int k = 0; k < Base::numnod_; k++)
    {
      const double fac = detJ_w * shapefct(k);

      for (int i = 0; i < Base::numnod_; i++) ecoupl_p1(k, i) += fac * dW_dp * shapefct(i);
    }
  }
}

#include "baci_so3_poro_p1_fwd.hpp"
