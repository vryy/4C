/*----------------------------------------------------------------------*/
/*! \file

 \brief evaluate methods for porous media using the p1 approach (mixed approach)

 \level 2

 *----------------------------------------------------------------------*/

#include "4C_global_data.hpp"
#include "4C_mat_fluidporo.hpp"
#include "4C_mat_structporo.hpp"
#include "4C_so3_poro_p1.hpp"
#include "4C_so3_poro_p1_eletypes.hpp"

FOUR_C_NAMESPACE_OPEN

template <class so3_ele, Core::FE::CellType distype>
void Discret::ELEMENTS::So3PoroP1<so3_ele, distype>::compute_porosity_and_linearization(
    Teuchos::ParameterList& params, const double& press, const double& J, const int& gp,
    const Core::LinAlg::Matrix<Base::numnod_, 1>& shapfct,
    const Core::LinAlg::Matrix<Base::numnod_, 1>* myporosity,
    const Core::LinAlg::Matrix<1, Base::numdof_>& dJ_dus, double& porosity,
    Core::LinAlg::Matrix<1, Base::numdof_>& dphi_dus)
{
  if (myporosity == nullptr)
    FOUR_C_THROW("no porosity values given!");
  else
    porosity = shapfct.Dot(*myporosity);

  dphi_dus.PutScalar(0.0);
}

template <class so3_ele, Core::FE::CellType distype>
void Discret::ELEMENTS::So3PoroP1<so3_ele, distype>::compute_porosity_and_linearization_od(
    Teuchos::ParameterList& params, const double& press, const double& J, const int& gp,
    const Core::LinAlg::Matrix<Base::numnod_, 1>& shapfct,
    const Core::LinAlg::Matrix<Base::numnod_, 1>* myporosity, double& porosity, double& dphi_dp)
{
  dphi_dp = 0.0;

  if (myporosity == nullptr)
    FOUR_C_THROW("no porosity values given!");
  else
    porosity = shapfct.Dot(*myporosity);
}

template <class so3_ele, Core::FE::CellType distype>
int Discret::ELEMENTS::So3PoroP1<so3_ele, distype>::Evaluate(Teuchos::ParameterList& params,
    Core::FE::Discretization& discretization, Core::Elements::Element::LocationArray& la,
    Core::LinAlg::SerialDenseMatrix& elemat1_epetra,
    Core::LinAlg::SerialDenseMatrix& elemat2_epetra,
    Core::LinAlg::SerialDenseVector& elevec1_epetra,
    Core::LinAlg::SerialDenseVector& elevec2_epetra,
    Core::LinAlg::SerialDenseVector& elevec3_epetra)
{
  // start with "none"
  typename Base::ActionType act = Base::none;

  // get the required action
  std::string action = params.get<std::string>("action", "none");
  if (action == "none")
    FOUR_C_THROW("No action supplied");
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
      pre_evaluate(params, discretization, la);

      my_evaluate(params, discretization, la, elemat1_epetra, elemat2_epetra, elevec1_epetra,
          elevec2_epetra, elevec3_epetra);
    }
    break;
    case Base::interpolate_porosity_to_given_point:
    {
      // given point
      std::vector<Core::LinAlg::Matrix<Base::numdim_, 1>> xsi(1);
      (xsi[0])(0, 0) = elevec2_epetra(0);
      (xsi[0])(1, 0) = elevec2_epetra(1);
      (xsi[0])(2, 0) = elevec2_epetra(2);

      //  evalulate shape functions at given point
      Core::LinAlg::Matrix<Base::numnod_, 1> shapefct;
      Core::FE::shape_function<distype>(xsi[0], shapefct);

      Core::LinAlg::Matrix<Base::numdim_, Base::numnod_> mydisp(true);
      Core::LinAlg::Matrix<Base::numnod_, 1> myporosity(true);
      Base::extract_values_from_global_vector(
          discretization, 0, la[0].lm_, nullptr, &myporosity, "displacement");

      elevec1_epetra(0) = shapefct.Dot(myporosity);
    }
    break;
    //==================================================================================
    default:
    {
      // in some cases we need to write/change some data before evaluating
      pre_evaluate(params, discretization, la);

      Core::LinAlg::SerialDenseMatrix elemat1_sub;
      Core::LinAlg::SerialDenseMatrix elemat2_sub;
      Core::LinAlg::SerialDenseVector elevec1_sub;
      Core::LinAlg::SerialDenseVector elevec2_sub;
      Core::LinAlg::SerialDenseVector elevec3_sub;

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
      my_evaluate(params, discretization, la, elemat1_epetra, elemat2_epetra, elevec1_epetra,
          elevec2_epetra, elevec3_epetra);
    }
    break;
  }

  return 0;
}

template <class so3_ele, Core::FE::CellType distype>
void Discret::ELEMENTS::So3PoroP1<so3_ele, distype>::pre_evaluate(Teuchos::ParameterList& params,
    Core::FE::Discretization& discretization, Core::Elements::Element::LocationArray& la)
{
  Base::pre_evaluate(params, discretization, la);

  if (discretization.HasState(0, "displacement") and (not is_init_porosity_))
  {
    init_porosity_ = Teuchos::rcp(new Core::LinAlg::Matrix<Base::numnod_, 1>(true));
    Base::extract_values_from_global_vector(
        discretization, 0, la[0].lm_, nullptr, &(*init_porosity_), "displacement");
    is_init_porosity_ = true;
  }
}

template <class so3_ele, Core::FE::CellType distype>
int Discret::ELEMENTS::So3PoroP1<so3_ele, distype>::my_evaluate(Teuchos::ParameterList& params,
    Core::FE::Discretization& discretization, Core::Elements::Element::LocationArray& la,
    Core::LinAlg::SerialDenseMatrix& elemat1_epetra,
    Core::LinAlg::SerialDenseMatrix& elemat2_epetra,
    Core::LinAlg::SerialDenseVector& elevec1_epetra,
    Core::LinAlg::SerialDenseVector& elevec2_epetra,
    Core::LinAlg::SerialDenseVector& elevec3_epetra)
{
  // start with "none"
  typename Base::ActionType act = Base::none;

  // get the required action
  std::string action = params.get<std::string>("action", "none");
  if (action == "none")
    FOUR_C_THROW("No action supplied");
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
        Core::LinAlg::Matrix<numdof_, numdof_> elemat1(elemat1_epetra.values(), true);
        // damping
        Core::LinAlg::Matrix<numdof_, numdof_> elemat2(elemat2_epetra.values(), true);
        // internal force vector
        Core::LinAlg::Matrix<numdof_, 1> elevec1(elevec1_epetra.values(), true);
        // elevec2+3 are not used anyway

        // build the location vector only for the structure field
        std::vector<int> lm = la[0].lm_;

        Core::LinAlg::Matrix<Base::numdim_, Base::numnod_> mydisp(true);
        Core::LinAlg::Matrix<Base::numnod_, 1> myporosity(true);
        Base::extract_values_from_global_vector(
            discretization, 0, la[0].lm_, &mydisp, &myporosity, "displacement");

        Core::LinAlg::Matrix<numdof_, numdof_>* matptr = nullptr;
        if (elemat1.IsInitialized()) matptr = &elemat1;

        enum Inpar::STR::DampKind damping =
            params.get<enum Inpar::STR::DampKind>("damping", Inpar::STR::damp_none);
        Core::LinAlg::Matrix<numdof_, numdof_>* matptr2 = nullptr;
        if (elemat2.IsInitialized() and (damping == Inpar::STR::damp_material)) matptr2 = &elemat2;

        // need current fluid state,
        // call the fluid discretization: fluid equates 2nd dofset
        // disassemble velocities and pressures

        Core::LinAlg::Matrix<Base::numdim_, Base::numnod_> myvel(true);

        Core::LinAlg::Matrix<Base::numdim_, Base::numnod_> myfluidvel(true);
        Core::LinAlg::Matrix<Base::numnod_, 1> myepreaf(true);

        if (discretization.HasState(0, "velocity"))
          Base::extract_values_from_global_vector(
              discretization, 0, la[0].lm_, &myvel, nullptr, "velocity");

        if (discretization.HasState(1, "fluidvel"))
        {
          // extract local values of the global vectors
          Base::extract_values_from_global_vector(
              discretization, 1, la[1].lm_, &myfluidvel, &myepreaf, "fluidvel");
        }

        // calculate tangent stiffness matrix
        nonlinear_stiffness_poroelast(lm, mydisp, myvel, &myporosity, myfluidvel, myepreaf, matptr,
            matptr2, &elevec1, params);
      }
    }
    break;

    //==================================================================================
    // coupling terms in force-vector and stiffness matrix for poroelasticity
    case Base::calc_struct_multidofsetcoupling:
    {
      // stiffness
      Core::LinAlg::Matrix<numdof_, (Base::numdim_ + 1) * Base::numnod_> elemat1(
          elemat1_epetra.values(), true);
      // elemat2,elevec1-3 are not used anyway

      // build the location vector only for the structure field
      std::vector<int> lm = la[0].lm_;

      Core::LinAlg::Matrix<numdof_, (Base::numdim_ + 1)* Base::numnod_>* matptr = nullptr;
      if (elemat1.IsInitialized()) matptr = &elemat1;

      // need current fluid state,
      // call the fluid discretization: fluid equates 2nd dofset
      // disassemble velocities and pressures
      if (discretization.HasState(1, "fluidvel"))
      {
        Core::LinAlg::Matrix<Base::numdim_, Base::numnod_> myvel(true);
        Core::LinAlg::Matrix<Base::numdim_, Base::numnod_> myfluidvel(true);
        Core::LinAlg::Matrix<Base::numnod_, 1> myepreaf(true);

        Core::LinAlg::Matrix<Base::numdim_, Base::numnod_> mydisp(true);
        Core::LinAlg::Matrix<Base::numnod_, 1> myporosity(true);
        Base::extract_values_from_global_vector(
            discretization, 0, la[0].lm_, &mydisp, &myporosity, "displacement");

        if (discretization.HasState(0, "velocity"))
          Base::extract_values_from_global_vector(
              discretization, 0, la[0].lm_, &myvel, nullptr, "velocity");

        if (discretization.HasState(1, "fluidvel"))
        {
          // extract local values of the global vectors
          Base::extract_values_from_global_vector(
              discretization, 1, la[1].lm_, &myfluidvel, &myepreaf, "fluidvel");
        }

        coupling_poroelast(
            lm, mydisp, myvel, &myporosity, myfluidvel, myepreaf, matptr, nullptr, nullptr, params);
      }
    }
    break;

    //==================================================================================
    // nonlinear stiffness and internal force vector for poroelasticity
    case Base::calc_struct_internalforce:
    {
      // internal force vector
      Core::LinAlg::Matrix<numdof_, 1> elevec1(elevec1_epetra.values(), true);
      // elemat1+2,elevec2+3 are not used anyway

      // build the location vector only for the structure field
      std::vector<int> lm = la[0].lm_;

      Core::LinAlg::Matrix<Base::numdim_, Base::numnod_> mydisp(true);
      Core::LinAlg::Matrix<Base::numnod_, 1> myporosity(true);
      Base::extract_values_from_global_vector(
          discretization, 0, la[0].lm_, &mydisp, &myporosity, "displacement");

      Core::LinAlg::Matrix<Base::numdim_, Base::numnod_> myvel(true);

      Core::LinAlg::Matrix<Base::numdim_, Base::numnod_> myfluidvel(true);
      Core::LinAlg::Matrix<Base::numnod_, 1> myepreaf(true);

      // need current fluid state,
      // call the fluid discretization: fluid equates 2nd dofset
      // disassemble velocities and pressures
      if (discretization.HasState(1, "fluidvel"))
      {
        // extract local values of the global vectors
        Base::extract_values_from_global_vector(
            discretization, 1, la[1].lm_, &myfluidvel, &myepreaf, "fluidvel");

        Base::extract_values_from_global_vector(
            discretization, 0, la[0].lm_, &myvel, nullptr, "velocity");

        nonlinear_stiffness_poroelast(lm, mydisp, myvel, &myporosity, myfluidvel, myepreaf, nullptr,
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

template <class so3_ele, Core::FE::CellType distype>
void Discret::ELEMENTS::So3PoroP1<so3_ele, distype>::InitElement()
{
  // initialize base element
  Base::InitElement();
}

template <class so3_ele, Core::FE::CellType distype>
void Discret::ELEMENTS::So3PoroP1<so3_ele, distype>::nonlinear_stiffness_poroelast(
    std::vector<int>& lm, Core::LinAlg::Matrix<Base::numdim_, Base::numnod_>& disp,
    Core::LinAlg::Matrix<Base::numdim_, Base::numnod_>& vel,
    Core::LinAlg::Matrix<Base::numnod_, 1>* porosity_dof,
    Core::LinAlg::Matrix<Base::numdim_, Base::numnod_>& evelnp,
    Core::LinAlg::Matrix<Base::numnod_, 1>& epreaf,
    Core::LinAlg::Matrix<numdof_, numdof_>* stiffmatrix,
    Core::LinAlg::Matrix<numdof_, numdof_>* reamatrix, Core::LinAlg::Matrix<numdof_, 1>* force,
    Teuchos::ParameterList& params)
{
  Base::get_materials();

  // update element geometry
  Core::LinAlg::Matrix<Base::numdim_, Base::numnod_> xrefe;  // material coord. of element
  Core::LinAlg::Matrix<Base::numdim_, Base::numnod_> xcurr;  // current  coord. of element

  Core::Nodes::Node** nodes = Base::Nodes();
  for (int i = 0; i < Base::numnod_; ++i)
  {
    const auto& x = nodes[i]->X();
    for (int j = 0; j < Base::numdim_; j++)
    {
      xrefe(j, i) = x[j];
      xcurr(j, i) = xrefe(j, i) + disp(j, i);
    }
  }

  // initialize element matrizes and vectors
  Core::LinAlg::Matrix<Base::numdof_, Base::numdof_> erea_v(true);
  Core::LinAlg::Matrix<Base::numdof_, Base::numdof_> sub_stiff(true);
  Core::LinAlg::Matrix<Base::numdof_, 1> sub_force(true);

  Core::LinAlg::Matrix<Base::numdof_, Base::numnod_> ecoupl_p1(true);
  Core::LinAlg::Matrix<Base::numnod_, numdof_> estiff_p1(true);
  Core::LinAlg::Matrix<Base::numnod_, 1> ecoupl_force_p1(true);

  gauss_point_loop_p1(params, xrefe, xcurr, disp, vel, evelnp, epreaf, porosity_dof, erea_v,
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

template <class so3_ele, Core::FE::CellType distype>
void Discret::ELEMENTS::So3PoroP1<so3_ele, distype>::gauss_point_loop_p1(
    Teuchos::ParameterList& params, const Core::LinAlg::Matrix<Base::numdim_, Base::numnod_>& xrefe,
    const Core::LinAlg::Matrix<Base::numdim_, Base::numnod_>& xcurr,
    const Core::LinAlg::Matrix<Base::numdim_, Base::numnod_>& nodaldisp,
    const Core::LinAlg::Matrix<Base::numdim_, Base::numnod_>& nodalvel,
    const Core::LinAlg::Matrix<Base::numdim_, Base::numnod_>& evelnp,
    const Core::LinAlg::Matrix<Base::numnod_, 1>& epreaf,
    const Core::LinAlg::Matrix<Base::numnod_, 1>* porosity_dof,
    Core::LinAlg::Matrix<Base::numdof_, Base::numdof_>& erea_v,
    Core::LinAlg::Matrix<Base::numdof_, Base::numdof_>* sub_stiff,
    Core::LinAlg::Matrix<Base::numdof_, 1>* sub_force,
    Core::LinAlg::Matrix<Base::numdof_, Base::numnod_>& ecoupl_p1,
    Core::LinAlg::Matrix<Base::numnod_, numdof_>& estiff_p1,
    Core::LinAlg::Matrix<Base::numnod_, 1>& ecoupl_force_p1)
{
  Core::LinAlg::Matrix<Base::numdim_, Base::numnod_> N_XYZ;
  // build deformation gradient wrt to material configuration
  // in case of prestressing, build defgrd wrt to last stored configuration
  // CAUTION: defgrd(true): filled with zeros!
  Core::LinAlg::Matrix<Base::numdim_, Base::numdim_> defgrd(true);
  Core::LinAlg::Matrix<Base::numnod_, 1> shapefct;
  Core::LinAlg::Matrix<Base::numdim_, Base::numnod_> deriv;

  Core::LinAlg::Matrix<Base::numstr_, 1> fstress(true);

  for (int gp = 0; gp < Base::numgpt_; ++gp)
  {
    // evaluate shape functions and derivatives at integration point
    Base::compute_shape_functions_and_derivatives(gp, shapefct, deriv, N_XYZ);

    // (material) deformation gradient F = d xcurr / d xrefe = xcurr * N_XYZ^T
    Base::compute_def_gradient(defgrd, N_XYZ, xcurr);

    // inverse deformation gradient F^-1
    static Core::LinAlg::Matrix<Base::numdim_, Base::numdim_> defgrd_inv(false);
    defgrd_inv.Invert(defgrd);

    // jacobian determinant of transformation between spatial and material space "|dx/dX|"
    double J = 0.0;
    //------linearization of jacobi determinant detF=J w.r.t. structure displacement   dJ/d(us) =
    // dJ/dF : dF/dus = J * F^-T * N,X
    static Core::LinAlg::Matrix<1, Base::numdof_> dJ_dus;
    // volume change (used for porosity law). Same as J in nonlinear theory.
    double volchange = 0.0;
    //------linearization of volume change w.r.t. structure displacement
    static Core::LinAlg::Matrix<1, Base::numdof_> dvolchange_dus;

    // compute J, the volume change and the respctive linearizations w.r.t. structure displacement
    Base::compute_jacobian_determinant_volume_change_and_linearizations(
        J, volchange, dJ_dus, dvolchange_dus, defgrd, defgrd_inv, N_XYZ, nodaldisp);

    // pressure at integration point
    double press = shapefct.Dot(epreaf);

    // structure velocity at integration point
    Core::LinAlg::Matrix<Base::numdim_, 1> velint(true);

    for (int i = 0; i < Base::numnod_; i++)
      for (int j = 0; j < Base::numdim_; j++) velint(j) += nodalvel(j, i) * shapefct(i);

    // fluid velocity at integration point
    Core::LinAlg::Matrix<Base::numdim_, 1> fvelint;
    fvelint.Multiply(evelnp, shapefct);

    // material fluid velocity gradient at integration point
    Core::LinAlg::Matrix<Base::numdim_, Base::numdim_> fvelder;
    fvelder.MultiplyNT(evelnp, N_XYZ);

    // pressure gradient at integration point
    Core::LinAlg::Matrix<Base::numdim_, 1> Gradp;
    Gradp.Multiply(N_XYZ, epreaf);

    // non-linear B-operator
    Core::LinAlg::Matrix<Base::numstr_, Base::numdof_> bop;
    Base::compute_b_operator(bop, defgrd, N_XYZ);

    // Right Cauchy-Green tensor = F^T * F
    Core::LinAlg::Matrix<Base::numdim_, Base::numdim_> cauchygreen;
    cauchygreen.MultiplyTN(defgrd, defgrd);

    // inverse Right Cauchy-Green tensor
    Core::LinAlg::Matrix<Base::numdim_, Base::numdim_> C_inv(false);
    C_inv.Invert(cauchygreen);

    //------linearization of material gradient of jacobi determinant GradJ  w.r.t. strucuture
    // displacement d(GradJ)/d(us)
    //---------------------d(GradJ)/dus =  dJ/dus * F^-T . : dF/dX + J * dF^-T/dus : dF/dX + J *
    // F^-T : N_X_X

    // compute some auxiliary matrixes for computation of linearization
    // dF^-T/dus
    Core::LinAlg::Matrix<Base::numdim_ * Base::numdim_, Base::numdof_> dFinvTdus(true);
    // F^-T * Grad p
    Core::LinAlg::Matrix<Base::numdim_, 1> Finvgradp;
    // dF^-T/dus * Grad p
    Core::LinAlg::Matrix<Base::numdim_, Base::numdof_> dFinvdus_gradp(true);
    // dC^-1/dus * Grad p
    Core::LinAlg::Matrix<Base::numstr_, Base::numdof_> dCinv_dus(true);

    Base::compute_auxiliary_values(
        N_XYZ, defgrd_inv, C_inv, Gradp, dFinvTdus, Finvgradp, dFinvdus_gradp, dCinv_dus);

    //--------------------------------------------------------------------

    // linearization of porosity w.r.t structure displacement d\phi/d(us) = d\phi/dJ*dJ/d(us)
    Core::LinAlg::Matrix<1, Base::numdof_> dphi_dus;
    double porosity = 0.0;

    compute_porosity_and_linearization(
        params, press, volchange, gp, shapefct, porosity_dof, dvolchange_dus, porosity, dphi_dus);

    double dW_dphi = 0.0;
    double dW_dJ = 0.0;
    double dW_dp = 0.0;
    double W = 0.0;

    if (init_porosity_ == Teuchos::null)
      FOUR_C_THROW("Failed to create vector of nodal intial porosity");
    double init_porosity = shapefct.Dot(*init_porosity_);
    Base::struct_mat_->constitutive_derivatives(params, press, volchange, porosity, init_porosity,
        &dW_dp,  // dW_dp not needed
        &dW_dphi, &dW_dJ, nullptr, &W);
    //--------------------------------------------------------

    // **********************evaluate stiffness matrix and force vector+++++++++++++++++++++++++
    if (Base::fluid_mat_->Type() == Mat::PAR::darcy_brinkman)
    {
      Base::fill_matrix_and_vectors_brinkman(gp, J, porosity, fvelder, defgrd_inv, bop, C_inv,
          dphi_dus, dJ_dus, dCinv_dus, dFinvTdus, sub_stiff, sub_force, fstress);
    }

    Base::fill_matrix_and_vectors(gp, shapefct, N_XYZ, J, press, porosity, velint, fvelint, fvelder,
        defgrd_inv, bop, C_inv, Finvgradp, dphi_dus, dJ_dus, dCinv_dus, dFinvdus_gradp, dFinvTdus,
        erea_v, sub_stiff, sub_force, fstress);

    // **********************evaluate stiffness matrix and force vector+++++++++++++++++++++++++
    double detJ_w = Base::detJ_[gp] * Base::intpoints_.Weight(gp);  // gpweights[gp];

    const double reacoeff = Base::fluid_mat_->compute_reaction_coeff();
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

    if (Base::fluid_mat_->Type() == Mat::PAR::darcy_brinkman)
    {
      double visc = Base::fluid_mat_->Viscosity();
      Core::LinAlg::Matrix<Base::numdim_, Base::numdim_> CinvFvel;
      Core::LinAlg::Matrix<Base::numdim_, Base::numdim_> visctress1;
      CinvFvel.Multiply(C_inv, fvelder);
      visctress1.MultiplyNT(CinvFvel, defgrd_inv);
      Core::LinAlg::Matrix<Base::numdim_, Base::numdim_> visctress2(visctress1);
      visctress1.UpdateT(1.0, visctress2, 1.0);

      fstress(0) = visctress1(0, 0);
      fstress(1) = visctress1(1, 1);
      fstress(2) = visctress1(2, 2);
      fstress(3) = visctress1(0, 1);
      fstress(4) = visctress1(1, 2);
      fstress(5) = visctress1(2, 0);

      fstress.Scale(detJ_w * visc * J);

      // B^T . C^-1
      Core::LinAlg::Matrix<Base::numdof_, 1> fstressb(true);
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

template <class so3_ele, Core::FE::CellType distype>
void Discret::ELEMENTS::So3PoroP1<so3_ele, distype>::coupling_poroelast(std::vector<int>& lm,
    Core::LinAlg::Matrix<Base::numdim_, Base::numnod_>& disp,
    Core::LinAlg::Matrix<Base::numdim_, Base::numnod_>& vel,
    Core::LinAlg::Matrix<Base::numnod_, 1>* porosity,
    Core::LinAlg::Matrix<Base::numdim_, Base::numnod_>& evelnp,
    Core::LinAlg::Matrix<Base::numnod_, 1>& epreaf,
    Core::LinAlg::Matrix<numdof_, (Base::numdim_ + 1) * Base::numnod_>* stiffmatrix,
    Core::LinAlg::Matrix<numdof_, (Base::numdim_ + 1) * Base::numnod_>* reamatrix,
    Core::LinAlg::Matrix<numdof_, 1>* force, Teuchos::ParameterList& params)
{
  //=============================get parameters

  Base::get_materials();

  //=======================================================================

  // update element geometry
  Core::LinAlg::Matrix<Base::numdim_, Base::numnod_> xrefe;  // material coord. of element
  Core::LinAlg::Matrix<Base::numdim_, Base::numnod_> xcurr;  // current  coord. of element

  Core::Nodes::Node** nodes = Base::Nodes();
  for (int i = 0; i < Base::numnod_; ++i)
  {
    const auto& x = nodes[i]->X();
    for (int j = 0; j < Base::numdim_; j++)
    {
      xrefe(j, i) = x[j];
      xcurr(j, i) = xrefe(j, i) + disp(j, i);
    }
  }
  // initialize element matrizes
  Core::LinAlg::Matrix<Base::numdof_, (Base::numdim_ + 1) * Base::numnod_> ecoupl(true);

  Core::LinAlg::Matrix<Base::numnod_, Base::numnod_> ecoupl_p1_p(true);

  gauss_point_loop_p1_od(
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

template <class so3_ele, Core::FE::CellType distype>
void Discret::ELEMENTS::So3PoroP1<so3_ele, distype>::gauss_point_loop_p1_od(
    Teuchos::ParameterList& params, const Core::LinAlg::Matrix<Base::numdim_, Base::numnod_>& xrefe,
    const Core::LinAlg::Matrix<Base::numdim_, Base::numnod_>& xcurr,
    const Core::LinAlg::Matrix<Base::numdim_, Base::numnod_>& nodaldisp,
    const Core::LinAlg::Matrix<Base::numdim_, Base::numnod_>& nodalvel,
    const Core::LinAlg::Matrix<Base::numdim_, Base::numnod_>& evelnp,
    const Core::LinAlg::Matrix<Base::numnod_, 1>& epreaf,
    const Core::LinAlg::Matrix<Base::numnod_, 1>* porosity_dof,
    Core::LinAlg::Matrix<Base::numnod_, Base::numnod_>& ecoupl_p1,
    Core::LinAlg::Matrix<Base::numdof_, (Base::numdim_ + 1) * Base::numnod_>* sub_stiff)
{
  Core::LinAlg::Matrix<Base::numdim_, Base::numnod_>
      N_XYZ;  //  first derivatives at gausspoint w.r.t. X, Y,Z
  // build deformation gradient wrt to material configuration
  // in case of prestressing, build defgrd wrt to last stored configuration
  // CAUTION: defgrd(true): filled with zeros!
  Core::LinAlg::Matrix<Base::numdim_, Base::numdim_> defgrd(
      true);  //  deformation gradiant evaluated at gauss point
  Core::LinAlg::Matrix<Base::numnod_, 1> shapefct;  //  shape functions evalulated at gauss point
  Core::LinAlg::Matrix<Base::numdim_, Base::numnod_> deriv(
      true);  //  first derivatives at gausspoint w.r.t. r,s,t

  for (int gp = 0; gp < Base::numgpt_; ++gp)
  {
    // evaluate shape functions and derivatives at integration point
    Base::compute_shape_functions_and_derivatives(gp, shapefct, deriv, N_XYZ);

    // (material) deformation gradient F = d xcurr / d xrefe = xcurr * N_XYZ^T
    defgrd.MultiplyNT(xcurr, N_XYZ);  //  (6.17)

    // inverse deformation gradient F^-1
    static Core::LinAlg::Matrix<Base::numdim_, Base::numdim_> defgrd_inv(false);
    defgrd_inv.Invert(defgrd);

    // jacobian determinant of transformation between spatial and material space "|dx/dX|"
    double J = 0.0;
    //------linearization of jacobi determinant detF=J w.r.t. structure displacement   dJ/d(us) =
    // dJ/dF : dF/dus = J * F^-T * N,X
    static Core::LinAlg::Matrix<1, numdof_> dJ_dus;
    // volume change (used for porosity law). Same as J in nonlinear theory.
    double volchange = 0.0;
    //------linearization of volume change w.r.t. structure displacement
    static Core::LinAlg::Matrix<1, numdof_> dvolchange_dus;

    // compute J, the volume change and the respctive linearizations w.r.t. structure displacement
    Base::compute_jacobian_determinant_volume_change(J, volchange, defgrd, N_XYZ, nodaldisp);

    // non-linear B-operator
    Core::LinAlg::Matrix<Base::numstr_, Base::numdof_> bop;
    Base::compute_b_operator(bop, defgrd, N_XYZ);

    // -----------------Right Cauchy-Green tensor = F^T * F
    Core::LinAlg::Matrix<Base::numdim_, Base::numdim_> cauchygreen;
    cauchygreen.MultiplyTN(defgrd, defgrd);

    //------------------ inverse Right Cauchy-Green tensor
    Core::LinAlg::Matrix<Base::numdim_, Base::numdim_> C_inv(false);
    C_inv.Invert(cauchygreen);

    //---------------- get pressure at integration point
    double press = shapefct.Dot(epreaf);

    //------------------ get material pressure gradient at integration point
    Core::LinAlg::Matrix<Base::numdim_, 1> Gradp;
    Gradp.Multiply(N_XYZ, epreaf);

    //--------------------- get fluid velocity at integration point
    Core::LinAlg::Matrix<Base::numdim_, 1> fvelint;
    fvelint.Multiply(evelnp, shapefct);

    //----------------------- material fluid velocity gradient at integration point
    Core::LinAlg::Matrix<Base::numdim_, Base::numdim_> fvelder;
    fvelder.MultiplyNT(evelnp, N_XYZ);

    //---------------- structure velocity at integration point
    Core::LinAlg::Matrix<Base::numdim_, 1> velint;
    velint.Multiply(nodalvel, shapefct);

    //**************************************************+auxilary variables for computing the
    // porosity and linearization
    double dphi_dp = 0.0;
    double porosity = 0.0;

    compute_porosity_and_linearization_od(
        params, press, volchange, gp, shapefct, porosity_dof, porosity, dphi_dp);

    // **********************evaluate stiffness matrix and force vector**********************

    Base::fill_matrix_and_vectors_od(gp, shapefct, N_XYZ, J, porosity, dphi_dp, velint, fvelint,
        defgrd_inv, Gradp, bop, C_inv, sub_stiff);

    if (Base::fluid_mat_->Type() == Mat::PAR::darcy_brinkman)
    {
      Base::fill_matrix_and_vectors_brinkman_od(
          gp, shapefct, N_XYZ, J, porosity, dphi_dp, fvelder, defgrd_inv, bop, C_inv, sub_stiff);
    }

    double dW_dp = 0.0;
    double init_porosity = shapefct.Dot(*init_porosity_);
    Base::struct_mat_->constitutive_derivatives(params, press, J, porosity, init_porosity, &dW_dp,
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

FOUR_C_NAMESPACE_CLOSE

#include "4C_so3_poro_p1_fwd.hpp"
