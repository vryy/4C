/*----------------------------------------------------------------------*/
/*! \file

\brief Evaluation methods for the 3D structural poro element


\level 2

*----------------------------------------------------------------------*/

#include "baci_discretization_fem_general_utils_nurbs_shapefunctions.H"
#include "baci_lib_globalproblem.H"
#include "baci_lib_utils.H"
#include "baci_linalg_serialdensevector.H"
#include "baci_linalg_utils_densematrix_multiply.H"
#include "baci_linalg_utils_sparse_algebra_math.H"
#include "baci_mat_fluidporo.H"
#include "baci_mat_fluidporo_multiphase.H"
#include "baci_mat_list.H"
#include "baci_mat_structporo.H"
#include "baci_nurbs_discret_nurbs_utils.H"
#include "baci_so3_poro.H"
#include "baci_so3_poro_eletypes.H"
#include "baci_structure_new_elements_paramsinterface.H"
#include "baci_utils_function.H"

#include <Teuchos_SerialDenseSolver.hpp>

#include <iterator>

BACI_NAMESPACE_OPEN

template <class so3_ele, CORE::FE::CellType distype>
void DRT::ELEMENTS::So3_Poro<so3_ele, distype>::PreEvaluate(Teuchos::ParameterList& params,
    DRT::Discretization& discretization, DRT::Element::LocationArray& la)
{
  if (scatra_coupling_)
  {
    if (la.Size() > 2)
    {
      if (discretization.HasState(2, "scalar"))
      {
        // check if you can get the scalar state
        Teuchos::RCP<const Epetra_Vector> scalarnp = discretization.GetState(2, "scalar");

        // extract local values of the global vectors
        std::vector<double> myscalar(la[2].lm_.size());
        DRT::UTILS::ExtractMyValues(*scalarnp, myscalar, la[2].lm_);

        if (so3_ele::NumMaterial() < 2)
          dserror("no second material defined for Wall poro element!");
        Teuchos::RCP<MAT::Material> scatramat = so3_ele::Material(2);

        int numscal = 1;
        if (scatramat->MaterialType() == INPAR::MAT::m_matlist or
            scatramat->MaterialType() == INPAR::MAT::m_matlist_reactions)
        {
          Teuchos::RCP<MAT::MatList> matlist = Teuchos::rcp_dynamic_cast<MAT::MatList>(scatramat);
          numscal = matlist->NumMat();
        }

        Teuchos::RCP<std::vector<double>> scalar =
            Teuchos::rcp(new std::vector<double>(numscal, 0.0));
        if ((int)myscalar.size() != numscal * numnod_) dserror("sizes do not match!");

        for (int i = 0; i < numnod_; i++)
          for (int j = 0; j < numscal; j++) scalar->at(j) += myscalar[numscal * i + j] / numnod_;

        params.set("scalar", scalar);
      }
    }
    else
    {
      const double time = params.get("total time", 0.0);
      // find out whether we will use a time curve and get the factor
      int num = 0;  // TO BE READ FROM INPUTFILE AT EACH ELEMENT!!!
      std::vector<double> xrefe;
      xrefe.resize(3);
      DRT::Node** nodes = Nodes();
      // get displacements of this element
      //  DRT::UTILS::ExtractMyValues(*disp,mydisp,lm);
      for (int i = 0; i < numnod_; ++i)
      {
        const auto& x = nodes[i]->X();
        xrefe[0] += x[0] / numnod_;
        xrefe[1] += x[1] / numnod_;
        xrefe[2] += x[2] / numnod_;
      }
      const double* coordgpref = xrefe.data();
      double functfac =
          DRT::Problem::Instance()->FunctionById<CORE::UTILS::FunctionOfSpaceTime>(num).Evaluate(
              coordgpref, time, 0);
      params.set<double>("scalar", functfac);
    }
  }
}

template <class so3_ele, CORE::FE::CellType distype>
int DRT::ELEMENTS::So3_Poro<so3_ele, distype>::Evaluate(Teuchos::ParameterList& params,
    DRT::Discretization& discretization, DRT::Element::LocationArray& la,
    CORE::LINALG::SerialDenseMatrix& elemat1_epetra,
    CORE::LINALG::SerialDenseMatrix& elemat2_epetra,
    CORE::LINALG::SerialDenseVector& elevec1_epetra,
    CORE::LINALG::SerialDenseVector& elevec2_epetra,
    CORE::LINALG::SerialDenseVector& elevec3_epetra)
{
  if (not init_) dserror("internal element data not initialized!");

  // set the pointer to the parameter list in element
  so3_ele::SetParamsInterfacePtr(params);

  // start with "none"
  ELEMENTS::ActionType act = ELEMENTS::none;

  if (so3_ele::IsParamsInterface())
  {
    act = so3_ele::ParamsInterface().GetActionType();
  }
  else
  {
    // get the required action
    std::string action = params.get<std::string>("action", "none");
    if (action == "none")
      dserror("No action supplied");
    else if (action == "struct_poro_calc_fluidcoupling")
      act = ELEMENTS::struct_poro_calc_fluidcoupling;
    else if (action == "struct_poro_calc_scatracoupling")
      act = ELEMENTS::struct_poro_calc_scatracoupling;
  }

  // what should the element do
  switch (act)
  {
    //==================================================================================
    // off diagonal terms in stiffness matrix for monolithic coupling
    case ELEMENTS::struct_poro_calc_fluidcoupling:
    {
      MyEvaluate(params, discretization, la, elemat1_epetra, elemat2_epetra, elevec1_epetra,
          elevec2_epetra, elevec3_epetra);
    }
    break;
    case ELEMENTS::struct_poro_calc_scatracoupling:
      // no coupling-> return
      break;
    //==================================================================================
    default:
    {
      // in some cases we need to write/change some data before evaluating
      PreEvaluate(params, discretization, la);

      // evaluate parent solid element
      so3_ele::Evaluate(params, discretization, la[0].lm_, elemat1_epetra, elemat2_epetra,
          elevec1_epetra, elevec2_epetra, elevec3_epetra);

      // add volume coupling specific terms
      MyEvaluate(params, discretization, la, elemat1_epetra, elemat2_epetra, elevec1_epetra,
          elevec2_epetra, elevec3_epetra);
    }
    break;
  }

  return 0;
}

template <class so3_ele, CORE::FE::CellType distype>
int DRT::ELEMENTS::So3_Poro<so3_ele, distype>::MyEvaluate(Teuchos::ParameterList& params,
    DRT::Discretization& discretization, DRT::Element::LocationArray& la,
    CORE::LINALG::SerialDenseMatrix& elemat1_epetra,
    CORE::LINALG::SerialDenseMatrix& elemat2_epetra,
    CORE::LINALG::SerialDenseVector& elevec1_epetra,
    CORE::LINALG::SerialDenseVector& elevec2_epetra,
    CORE::LINALG::SerialDenseVector& elevec3_epetra)
{
  // start with "none"
  // ActionType act = none;
  ELEMENTS::ActionType act = ELEMENTS::none;

  if (so3_ele::IsParamsInterface())
  {
    act = so3_ele::ParamsInterface().GetActionType();
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
    else if (action == "calc_struct_stress")
      act = ELEMENTS::struct_calc_stress;
  }

  // what should the element do
  switch (act)
  {
    //==================================================================================
    // nonlinear stiffness, damping and internal force vector for poroelasticity
    case ELEMENTS::struct_calc_nlnstiff:
    {
      // stiffness
      CORE::LINALG::Matrix<numdof_, numdof_> elemat1(elemat1_epetra.values(), true);
      // damping
      CORE::LINALG::Matrix<numdof_, numdof_> elemat2(elemat2_epetra.values(), true);
      // internal force vector
      CORE::LINALG::Matrix<numdof_, 1> elevec1(elevec1_epetra.values(), true);
      // elevec2+3 are not used anyway

      std::vector<int> lm = la[0].lm_;

      CORE::LINALG::Matrix<numdim_, numnod_> mydisp(true);
      ExtractValuesFromGlobalVector(discretization, 0, lm, &mydisp, nullptr, "displacement");

      CORE::LINALG::Matrix<numdof_, numdof_>* matptr = nullptr;
      if (elemat1.IsInitialized()) matptr = &elemat1;

      enum INPAR::STR::DampKind damping =
          params.get<enum INPAR::STR::DampKind>("damping", INPAR::STR::damp_none);
      CORE::LINALG::Matrix<numdof_, numdof_>* matptr2 = nullptr;
      if (elemat2.IsInitialized() and (damping == INPAR::STR::damp_material)) matptr2 = &elemat2;

      if (la.Size() > 1)
      {
        if (discretization.HasState(1, "fluidvel"))
        {
          // need current fluid state,
          // call the fluid discretization: fluid equates 2nd dofset
          // disassemble velocities and pressures
          CORE::LINALG::Matrix<numdim_, numnod_> myvel(true);
          CORE::LINALG::Matrix<numdim_, numnod_> myfluidvel(true);
          CORE::LINALG::Matrix<numnod_, 1> myepreaf(true);

          if (discretization.HasState(0, "velocity"))
            ExtractValuesFromGlobalVector(
                discretization, 0, la[0].lm_, &myvel, nullptr, "velocity");

          if (discretization.HasState(1, "fluidvel"))
          {
            // extract local values of the global vectors
            ExtractValuesFromGlobalVector(
                discretization, 1, la[1].lm_, &myfluidvel, &myepreaf, "fluidvel");
          }

          // calculate tangent stiffness matrix
          NonlinearStiffnessPoroelast(
              lm, mydisp, myvel, myfluidvel, myepreaf, matptr, matptr2, &elevec1, params);
        }
        else if (la.Size() > 2)
        {
          if (discretization.HasState(1, "porofluid"))
          {
            // get primary variables of multiphase porous medium flow
            std::vector<double> myephi(la[1].Size());
            Teuchos::RCP<const Epetra_Vector> matrix_state =
                discretization.GetState(1, "porofluid");
            DRT::UTILS::ExtractMyValues(*matrix_state, myephi, la[1].lm_);

            // calculate tangent stiffness matrix
            NonlinearStiffnessPoroelastPressureBased(lm, mydisp, myephi, matptr, &elevec1, params);
          }
        }
      }
    }
    break;

    //==================================================================================
    // nonlinear stiffness, mass matrix and internal force vector for poroelasticity
    case ELEMENTS::struct_calc_nlnstiffmass:
    {
      // stiffness
      CORE::LINALG::Matrix<numdof_, numdof_> elemat1(elemat1_epetra.values(), true);
      // internal force vector
      CORE::LINALG::Matrix<numdof_, 1> elevec1(elevec1_epetra.values(), true);
      // elemat2,elevec2+3 are not used anyway

      // build the location vector only for the structure field
      std::vector<int> lm = la[0].lm_;

      CORE::LINALG::Matrix<numdim_, numnod_> mydisp(true);
      ExtractValuesFromGlobalVector(discretization, 0, la[0].lm_, &mydisp, nullptr, "displacement");

      CORE::LINALG::Matrix<numdof_, numdof_>* matptr = nullptr;
      if (elemat1.IsInitialized()) matptr = &elemat1;

      if (isNurbs_)
      {
        // access knots and weights for this element
        bool zero_size =
            DRT::NURBS::GetMyNurbsKnotsAndWeights(discretization, this, myknots_, weights_);

        // if we have a zero sized element due to a interpolated point -> exit here
        if (zero_size) return 0;
      }

      // we skip this evaluation if the coupling is not setup yet, i.e.
      // if the secondary dofset or the secondary material was not set
      // this can happen during setup of the time integrator or restart
      // there might be a better way. For instance do not evaluate
      // before the setup of the multiphysics problem is completed.
      if (la.Size() > 1 and so3_ele::NumMaterial() > 1)
      {
        // need current fluid state,
        // call the fluid discretization: fluid equates 2nd dofset
        // disassemble velocities and pressures
        CORE::LINALG::Matrix<numdim_, numnod_> myvel(true);
        CORE::LINALG::Matrix<numdim_, numnod_> myfluidvel(true);
        CORE::LINALG::Matrix<numnod_, 1> myepreaf(true);

        if (discretization.HasState(0, "velocity"))
          ExtractValuesFromGlobalVector(discretization, 0, la[0].lm_, &myvel, nullptr, "velocity");

        // this is kind of a hack. Find a better way! (e.g. move the pressure based variant
        // into own element)
        if (discretization.HasState(1, "fluidvel"))
        {
          // extract local values of the global vectors
          ExtractValuesFromGlobalVector(
              discretization, 1, la[1].lm_, &myfluidvel, &myepreaf, "fluidvel");

          // calculate tangent stiffness matrix
          NonlinearStiffnessPoroelast(
              lm, mydisp, myvel, myfluidvel, myepreaf, matptr, nullptr, &elevec1, params);
        }
        else if (la.Size() > 2)
        {
          if (discretization.HasState(1, "porofluid"))
          {
            // get primary variables of multiphase porous medium flow
            std::vector<double> myephi(la[1].Size());
            Teuchos::RCP<const Epetra_Vector> matrix_state =
                discretization.GetState(1, "porofluid");
            DRT::UTILS::ExtractMyValues(*matrix_state, myephi, la[1].lm_);

            // calculate tangent stiffness matrix
            NonlinearStiffnessPoroelastPressureBased(lm, mydisp, myephi, matptr, &elevec1, params);
          }
        }
      }
    }
    break;

    //==================================================================================
    // coupling terms in force-vector and stiffness matrix for poroelasticity
    case ELEMENTS::struct_poro_calc_fluidcoupling:
    {
      // stiffness
      CORE::LINALG::Matrix<numdof_, (numdim_ + 1) * numnod_> elemat1(elemat1_epetra.values(), true);

      // build the location vector only for the structure field
      std::vector<int> lm = la[0].lm_;

      CORE::LINALG::Matrix<numdof_, (numdim_ + 1)* numnod_>* matptr = nullptr;
      if (elemat1.IsInitialized()) matptr = &elemat1;

      if (isNurbs_)
      {
        // access knots and weights for this element
        bool zero_size =
            DRT::NURBS::GetMyNurbsKnotsAndWeights(discretization, this, myknots_, weights_);

        // if we have a zero sized element due to a interpolated point -> exit here
        if (zero_size) return 0;
      }

      // need current fluid state,
      // call the fluid discretization: fluid equates 2nd dofset
      // disassemble velocities and pressures
      if (discretization.HasState(1, "fluidvel"))
      {
        CORE::LINALG::Matrix<numdim_, numnod_> myvel(true);
        CORE::LINALG::Matrix<numdim_, numnod_> myfluidvel(true);
        CORE::LINALG::Matrix<numnod_, 1> myepreaf(true);

        CORE::LINALG::Matrix<numdim_, numnod_> mydisp(true);
        ExtractValuesFromGlobalVector(
            discretization, 0, la[0].lm_, &mydisp, nullptr, "displacement");

        if (discretization.HasState(0, "velocity"))
          ExtractValuesFromGlobalVector(discretization, 0, la[0].lm_, &myvel, nullptr, "velocity");

        if (discretization.HasState(1, "fluidvel"))
        {
          // extract local values of the global vectors
          ExtractValuesFromGlobalVector(
              discretization, 1, la[1].lm_, &myfluidvel, &myepreaf, "fluidvel");
        }

        CouplingPoroelast(
            lm, mydisp, myvel, myfluidvel, myepreaf, matptr, nullptr, nullptr, params);
      }
      else if (la.Size() > 2)
      {
        if (discretization.HasState(1, "porofluid"))
        {
          // get primary variables of multiphase porous medium flow
          std::vector<double> myephi(la[1].Size());
          Teuchos::RCP<const Epetra_Vector> matrix_state = discretization.GetState(1, "porofluid");
          DRT::UTILS::ExtractMyValues(*matrix_state, myephi, la[1].lm_);

          CORE::LINALG::Matrix<numdim_, numnod_> mydisp(true);
          ExtractValuesFromGlobalVector(
              discretization, 0, la[0].lm_, &mydisp, nullptr, "displacement");

          // calculate OD-Matrix
          CouplingPoroelastPressureBased(lm, mydisp, myephi, elemat1_epetra, params);
        }
        else
          dserror("cannot find global states displacement or solidpressure");
      }
    }
    break;

    //==================================================================================
    // nonlinear stiffness and internal force vector for poroelasticity
    case ELEMENTS::struct_calc_internalforce:
    {
      // internal force vector
      CORE::LINALG::Matrix<numdof_, 1> elevec1(elevec1_epetra.values(), true);
      // elemat1+2,elevec2+3 are not used anyway

      // build the location vector only for the structure field
      std::vector<int> lm = la[0].lm_;

      CORE::LINALG::Matrix<numdim_, numnod_> mydisp(true);
      ExtractValuesFromGlobalVector(discretization, 0, lm, &mydisp, nullptr, "displacement");

      // need current fluid state,
      // call the fluid discretization: fluid equates 2nd dofset
      // disassemble velocities and pressures
      if (discretization.HasState(1, "fluidvel"))
      {
        // extract local values of the global vectors
        CORE::LINALG::Matrix<numdim_, numnod_> myfluidvel(true);
        CORE::LINALG::Matrix<numnod_, 1> myepreaf(true);
        ExtractValuesFromGlobalVector(
            discretization, 1, la[1].lm_, &myfluidvel, &myepreaf, "fluidvel");

        CORE::LINALG::Matrix<numdim_, numnod_> myvel(true);
        ExtractValuesFromGlobalVector(discretization, 0, la[0].lm_, &myvel, nullptr, "velocity");

        // calculate tangent stiffness matrix
        NonlinearStiffnessPoroelast(
            lm, mydisp, myvel, myfluidvel, myepreaf, nullptr, nullptr, &elevec1, params);
      }
      else if (la.Size() > 2)
      {
        if (discretization.HasState(1, "porofluid"))
        {
          // get primary variables of multiphase porous medium flow
          std::vector<double> myephi(la[1].Size());
          Teuchos::RCP<const Epetra_Vector> matrix_state = discretization.GetState(1, "porofluid");
          DRT::UTILS::ExtractMyValues(*matrix_state, myephi, la[1].lm_);

          // calculate tangent stiffness matrix
          NonlinearStiffnessPoroelastPressureBased(lm, mydisp, myephi, nullptr, &elevec1, params);
        }
      }
    }
    break;

    //==================================================================================
    // evaluate stresses and strains at gauss points
    case ELEMENTS::struct_calc_stress:
    {
      // nothing to do for ghost elements
      if (discretization.Comm().MyPID() == so3_ele::Owner())
      {
        // get the location vector only for the structure field
        std::vector<int> lm = la[0].lm_;

        CORE::LINALG::Matrix<numdim_, numnod_> mydisp(true);
        ExtractValuesFromGlobalVector(discretization, 0, lm, &mydisp, nullptr, "displacement");

        Teuchos::RCP<std::vector<char>> couplingstressdata = Teuchos::null;
        INPAR::STR::StressType iocouplingstress = INPAR::STR::stress_none;
        if (this->IsParamsInterface())
        {
          couplingstressdata = this->StrParamsInterface().CouplingStressDataPtr();
          iocouplingstress = this->StrParamsInterface().GetCouplingStressOutputType();
        }
        else
        {
          iocouplingstress = DRT::INPUT::get<INPAR::STR::StressType>(
              params, "iocouplstress", INPAR::STR::stress_none);

          couplingstressdata =
              params.get<Teuchos::RCP<std::vector<char>>>("couplstress", Teuchos::null);

          if (couplingstressdata == Teuchos::null) dserror("Cannot get 'couplstress' data");
        }

        // initialize the coupling stress
        CORE::LINALG::SerialDenseMatrix couplstress(numgpt_, numstr_);
        // need current fluid state,
        // call the fluid discretization: fluid equates 2nd dofset
        // disassemble velocities and pressures
        if (iocouplingstress != INPAR::STR::stress_none)
        {
          if (discretization.HasState(1, "fluidvel"))
          {
            // extract local values of the global vectors
            CORE::LINALG::Matrix<numdim_, numnod_> myfluidvel(true);
            CORE::LINALG::Matrix<numnod_, 1> myepreaf(true);
            ExtractValuesFromGlobalVector(
                discretization, 1, la[1].lm_, &myfluidvel, &myepreaf, "fluidvel");

            CouplingStressPoroelast(
                mydisp, myfluidvel, myepreaf, &couplstress, nullptr, params, iocouplingstress);
          }
          else if (la.Size() > 2)
          {
            if (discretization.HasState(1, "porofluid"))
            {
              dserror("coupl stress poroelast not yet implemented for pressure-based variant");
            }
          }
        }

        // pack the data for postprocessing
        {
          CORE::COMM::PackBuffer data;
          // get the size of stress
          so3_ele::AddtoPack(data, couplstress);
          data.StartPacking();
          // pack the stresses
          so3_ele::AddtoPack(data, couplstress);
          std::copy(data().begin(), data().end(), std::back_inserter(*couplingstressdata));
        }
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

template <class so3_ele, CORE::FE::CellType distype>
void DRT::ELEMENTS::So3_Poro<so3_ele, distype>::NonlinearStiffnessPoroelast(
    std::vector<int>& lm,                                 // location matrix
    CORE::LINALG::Matrix<numdim_, numnod_>& disp,         // current displacements
    CORE::LINALG::Matrix<numdim_, numnod_>& vel,          // current velocities
    CORE::LINALG::Matrix<numdim_, numnod_>& evelnp,       // current fluid velocities
    CORE::LINALG::Matrix<numnod_, 1>& epreaf,             // current fluid pressure
    CORE::LINALG::Matrix<numdof_, numdof_>* stiffmatrix,  // element stiffness matrix
    CORE::LINALG::Matrix<numdof_, numdof_>* reamatrix,    // element reactive matrix
    CORE::LINALG::Matrix<numdof_, 1>* force,              // element internal force vector
    // CORE::LINALG::Matrix<numgptpar_, numstr_>* elestress, // stresses at GP
    // CORE::LINALG::Matrix<numgptpar_, numstr_>* elestrain, // strains at GP
    Teuchos::ParameterList& params  // algorithmic parameters e.g. time
    //   const INPAR::STR::StressType       iostress     // stress output option
)
{
  GetMaterials();

  // update element geometry
  CORE::LINALG::Matrix<numdim_, numnod_> xrefe;  // material coord. of element
  CORE::LINALG::Matrix<numdim_, numnod_> xcurr;  // current  coord. of element

  DRT::Node** nodes = Nodes();
  for (int i = 0; i < numnod_; ++i)
  {
    const auto& x = nodes[i]->X();
    for (int j = 0; j < numdim_; j++)
    {
      xrefe(j, i) = x[j];
      xcurr(j, i) = xrefe(j, i) + disp(j, i);
    }
  }

  // initialize element matrizes and vectors
  CORE::LINALG::Matrix<numdof_, numdof_> erea_v(true);

  /* =========================================================================*/
  /* ================================================= Loop over Gauss Points */
  /* =========================================================================*/
  GaussPointLoop(
      params, xrefe, xcurr, disp, vel, evelnp, epreaf, nullptr, erea_v, stiffmatrix, force);

  // update stiffness matrix
  if (stiffmatrix != nullptr)
  {
    if (reamatrix != nullptr)
    {
      /* additional "reactive darcy-term"
       detJ * w(gp) * ( J * reacoeff * phi^2  ) * D(v_s)
       */
      reamatrix->Update(1.0, erea_v, 1.0);
    }
  }
}

template <class so3_ele, CORE::FE::CellType distype>
void DRT::ELEMENTS::So3_Poro<so3_ele, distype>::NonlinearStiffnessPoroelastPressureBased(
    std::vector<int>& lm,                          // location matrix
    CORE::LINALG::Matrix<numdim_, numnod_>& disp,  // current displacements
    const std::vector<double>& ephi,               // primary variable for poro-multiphase flow
    CORE::LINALG::Matrix<numdof_, numdof_>* stiffmatrix,  // element stiffness matrix
    CORE::LINALG::Matrix<numdof_, 1>* force,              // element internal force vector
    Teuchos::ParameterList& params                        // algorithmic parameters e.g. time
)
{
  GetMaterialsPressureBased();

  // update element geometry
  CORE::LINALG::Matrix<numdim_, numnod_> xrefe;  // material coord. of element
  CORE::LINALG::Matrix<numdim_, numnod_> xcurr;  // current  coord. of element

  DRT::Node** nodes = Nodes();
  for (int i = 0; i < numnod_; ++i)
  {
    const auto& x = nodes[i]->X();
    for (int j = 0; j < numdim_; j++)
    {
      xrefe(j, i) = x[j];
      xcurr(j, i) = xrefe(j, i) + disp(j, i);
    }
  }

  /* =========================================================================*/
  /* ================================================= Loop over Gauss Points */
  /* =========================================================================*/
  GaussPointLoopPressureBased(params, xrefe, xcurr, disp, ephi, stiffmatrix, force);
}

template <class so3_ele, CORE::FE::CellType distype>
void DRT::ELEMENTS::So3_Poro<so3_ele, distype>::GaussPointLoop(Teuchos::ParameterList& params,
    const CORE::LINALG::Matrix<numdim_, numnod_>& xrefe,
    const CORE::LINALG::Matrix<numdim_, numnod_>& xcurr,
    const CORE::LINALG::Matrix<numdim_, numnod_>& nodaldisp,
    const CORE::LINALG::Matrix<numdim_, numnod_>& nodalvel,
    const CORE::LINALG::Matrix<numdim_, numnod_>& evelnp,
    const CORE::LINALG::Matrix<numnod_, 1>& epreaf,
    const CORE::LINALG::Matrix<numnod_, 1>* porosity_dof,
    CORE::LINALG::Matrix<numdof_, numdof_>& erea_v,
    CORE::LINALG::Matrix<numdof_, numdof_>* stiffmatrix, CORE::LINALG::Matrix<numdof_, 1>* force)
{
  static CORE::LINALG::Matrix<numdim_, numnod_> N_XYZ;
  // build deformation gradient wrt to material configuration
  // in case of prestressing, build defgrd wrt to last stored configuration
  // CAUTION: defgrd(true): filled with zeros!
  static CORE::LINALG::Matrix<numdim_, numdim_> defgrd(true);
  static CORE::LINALG::Matrix<numnod_, 1> shapefct;
  static CORE::LINALG::Matrix<numdim_, numnod_> deriv;

  static CORE::LINALG::Matrix<numstr_, 1> fstress(true);

  for (int gp = 0; gp < numgpt_; ++gp)
  {
    // evaluate shape functions and derivatives at integration point
    ComputeShapeFunctionsAndDerivatives(gp, shapefct, deriv, N_XYZ);

    // (material) deformation gradient F = d xcurr / d xrefe = xcurr * N_XYZ^T
    ComputeDefGradient(defgrd, N_XYZ, xcurr);

    // inverse deformation gradient F^-1
    static CORE::LINALG::Matrix<numdim_, numdim_> defgrd_inv(false);
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
    ComputeJacobianDeterminantVolumeChangeAndLinearizations(
        J, volchange, dJ_dus, dvolchange_dus, defgrd, defgrd_inv, N_XYZ, nodaldisp);

    // pressure at integration point
    double press = shapefct.Dot(epreaf);

    // structure displacement and velocity at integration point
    static CORE::LINALG::Matrix<numdim_, 1> velint;
    velint.Multiply(nodalvel, shapefct);

    // fluid velocity at integration point
    static CORE::LINALG::Matrix<numdim_, 1> fvelint;
    fvelint.Multiply(evelnp, shapefct);

    // material fluid velocity gradient at integration point
    static CORE::LINALG::Matrix<numdim_, numdim_> fvelder;
    fvelder.MultiplyNT(evelnp, N_XYZ);

    // pressure gradient at integration point
    static CORE::LINALG::Matrix<numdim_, 1> Gradp;
    Gradp.Multiply(N_XYZ, epreaf);

    // non-linear B-operator
    static CORE::LINALG::Matrix<numstr_, numdof_> bop;
    bop.Clear();
    ComputeBOperator(bop, defgrd, N_XYZ);

    // Right Cauchy-Green tensor = F^T * F
    static CORE::LINALG::Matrix<numdim_, numdim_> cauchygreen;
    cauchygreen.MultiplyTN(defgrd, defgrd);

    // inverse Right Cauchy-Green tensor
    static CORE::LINALG::Matrix<numdim_, numdim_> C_inv(false);
    C_inv.Invert(cauchygreen);

    // compute some auxiliary matrixes for computation of linearization
    // dF^-T/dus
    static CORE::LINALG::Matrix<numdim_ * numdim_, numdof_> dFinvTdus(true);
    // F^-T * Grad p
    static CORE::LINALG::Matrix<numdim_, 1> Finvgradp;
    // dF^-T/dus * Grad p
    static CORE::LINALG::Matrix<numdim_, numdof_> dFinvdus_gradp(true);
    // dC^-1/dus * Grad p
    static CORE::LINALG::Matrix<numstr_, numdof_> dCinv_dus(true);

    ComputeAuxiliaryValues(
        N_XYZ, defgrd_inv, C_inv, Gradp, dFinvTdus, Finvgradp, dFinvdus_gradp, dCinv_dus);

    // linearization of porosity w.r.t structure displacement d\phi/d(us) = d\phi/dJ*dJ/d(us)
    static CORE::LINALG::Matrix<1, numdof_> dphi_dus;
    double porosity = 0.0;

    ComputePorosityAndLinearization(
        params, press, volchange, gp, shapefct, porosity_dof, dvolchange_dus, porosity, dphi_dus);

    // **********************fill stiffness matrix and force vector+++++++++++++++++++++++++
    if (fluid_mat_->Type() == MAT::PAR::darcy_brinkman)
    {
      FillMatrixAndVectorsBrinkman(gp, J, porosity, fvelder, defgrd_inv, bop, C_inv, dphi_dus,
          dJ_dus, dCinv_dus, dFinvTdus, stiffmatrix, force, fstress);
    }

    FillMatrixAndVectors(gp, shapefct, N_XYZ, J, press, porosity, velint, fvelint, fvelder,
        defgrd_inv, bop, C_inv, Finvgradp, dphi_dus, dJ_dus, dCinv_dus, dFinvdus_gradp, dFinvTdus,
        erea_v, stiffmatrix, force, fstress);
  }
}

template <class so3_ele, CORE::FE::CellType distype>
void DRT::ELEMENTS::So3_Poro<so3_ele, distype>::GaussPointLoopPressureBased(
    Teuchos::ParameterList& params, const CORE::LINALG::Matrix<numdim_, numnod_>& xrefe,
    const CORE::LINALG::Matrix<numdim_, numnod_>& xcurr,
    const CORE::LINALG::Matrix<numdim_, numnod_>& nodaldisp, const std::vector<double>& ephi,
    CORE::LINALG::Matrix<numdof_, numdof_>* stiffmatrix, CORE::LINALG::Matrix<numdof_, 1>* force)
{
  /*--------------------------------- get node weights for nurbs elements */
  if (distype == CORE::FE::CellType::nurbs4 || distype == CORE::FE::CellType::nurbs9)
  {
    for (int inode = 0; inode < numnod_; ++inode)
    {
      auto* cp = dynamic_cast<DRT::NURBS::ControlPoint*>(Nodes()[inode]);

      weights_(inode) = cp->W();
    }
  }

  /* =========================================================================*/
  /* ================================================= Loop over Gauss Points */
  /* =========================================================================*/
  // first derivatives N_XYZ at gp w.r.t. material coordinates
  CORE::LINALG::Matrix<numdim_, numnod_> N_XYZ;
  // build deformation gradient wrt to material configuration
  // in case of prestressing, build defgrd wrt to last stored configuration
  // CAUTION: defgrd(true): filled with zeros!
  CORE::LINALG::Matrix<numdim_, numdim_> defgrd(true);
  // shape function at gp w.r.t. reference coordinates
  CORE::LINALG::Matrix<numnod_, 1> shapefct;
  // first derivatives at gp w.r.t. reference coordinates
  CORE::LINALG::Matrix<numdim_, numnod_> deriv;

  // Initialize
  const int totalnumdofpernode = fluidmulti_mat_->NumMat();
  const int numfluidphases = fluidmulti_mat_->NumFluidPhases();
  const int numvolfrac = fluidmulti_mat_->NumVolFrac();
  const bool hasvolfracs = (totalnumdofpernode > numfluidphases);
  std::vector<double> phiAtGP(totalnumdofpernode);

  for (int gp = 0; gp < numgpt_; ++gp)
  {
    // evaluate shape functions and derivatives at integration point
    ComputeShapeFunctionsAndDerivatives(gp, shapefct, deriv, N_XYZ);

    // compute deformation gradient
    ComputeDefGradient(defgrd, N_XYZ, xcurr);

    // inverse deformation gradient F^-1
    CORE::LINALG::Matrix<numdim_, numdim_> defgrd_inv(false);
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
    ComputeJacobianDeterminantVolumeChangeAndLinearizations(
        J, volchange, dJ_dus, dvolchange_dus, defgrd, defgrd_inv, N_XYZ, nodaldisp);

    // non-linear B-operator
    CORE::LINALG::Matrix<numstr_, numdof_> bop;
    bop.Clear();
    ComputeBOperator(bop, defgrd, N_XYZ);

    // derivative of press w.r.t. displacements (only in case of vol fracs)
    CORE::LINALG::Matrix<1, numdof_> dps_dus(true);

    //----------------------------------------------------
    // pressure at integration point
    ComputePrimaryVariableAtGP(ephi, totalnumdofpernode, shapefct, phiAtGP);
    double press = ComputeSolPressureAtGP(totalnumdofpernode, numfluidphases, phiAtGP);
    // recalculate for the case of volume fractions
    if (hasvolfracs)
    {
      CORE::LINALG::Matrix<1, numdof_> dphi_dus;
      double porosity = 0.0;

      ComputePorosityAndLinearization(
          params, press, volchange, gp, shapefct, nullptr, dvolchange_dus, porosity, dphi_dus);
      // save the pressure coming from the fluid S_i*p_i
      const double fluidpress = press;
      press = RecalculateSolPressureAtGP(
          fluidpress, porosity, totalnumdofpernode, numfluidphases, numvolfrac, phiAtGP);
      ComputeLinearizationOfSolPressWrtDisp(fluidpress, porosity, totalnumdofpernode,
          numfluidphases, numvolfrac, phiAtGP, dphi_dus, dps_dus);
    }

    // Right Cauchy-Green tensor = F^T * F
    CORE::LINALG::Matrix<numdim_, numdim_> cauchygreen;
    cauchygreen.MultiplyTN(defgrd, defgrd);

    // inverse Right Cauchy-Green tensor
    CORE::LINALG::Matrix<numdim_, numdim_> C_inv(false);
    C_inv.Invert(cauchygreen);

    // compute some auxiliary matrixes for computation of linearization
    // dC^-1/dus
    CORE::LINALG::Matrix<numstr_, numdof_> dCinv_dus(true);
    for (int n = 0; n < numnod_; ++n)
    {
      for (int k = 0; k < numdim_; ++k)
      {
        const int gid = n * numdim_ + k;
        for (int i = 0; i < numdim_; ++i)
        {
          dCinv_dus(0, gid) += -2 * C_inv(0, i) * N_XYZ(i, n) * defgrd_inv(0, k);
          dCinv_dus(1, gid) += -2 * C_inv(1, i) * N_XYZ(i, n) * defgrd_inv(1, k);
          dCinv_dus(2, gid) += -2 * C_inv(2, i) * N_XYZ(i, n) * defgrd_inv(2, k);
          /* ~~~ */
          dCinv_dus(3, gid) += -C_inv(0, i) * N_XYZ(i, n) * defgrd_inv(1, k) -
                               defgrd_inv(0, k) * N_XYZ(i, n) * C_inv(1, i);
          dCinv_dus(4, gid) += -C_inv(1, i) * N_XYZ(i, n) * defgrd_inv(2, k) -
                               defgrd_inv(1, k) * N_XYZ(i, n) * C_inv(2, i);
          dCinv_dus(5, gid) += -C_inv(2, i) * N_XYZ(i, n) * defgrd_inv(0, k) -
                               defgrd_inv(2, k) * N_XYZ(i, n) * C_inv(0, i);
        }
      }
    }

    // **********************evaluate stiffness matrix and force vector+++++++++++++++++++++++++
    FillMatrixAndVectorsPressureBased(
        gp, shapefct, N_XYZ, J, press, bop, C_inv, dJ_dus, dCinv_dus, dps_dus, stiffmatrix, force);
  }
}

template <class so3_ele, CORE::FE::CellType distype>
void DRT::ELEMENTS::So3_Poro<so3_ele, distype>::CouplingPoroelast(
    std::vector<int>& lm,                            // location matrix
    CORE::LINALG::Matrix<numdim_, numnod_>& disp,    // current displacements
    CORE::LINALG::Matrix<numdim_, numnod_>& vel,     // current velocities
    CORE::LINALG::Matrix<numdim_, numnod_>& evelnp,  // current fluid velocity
    CORE::LINALG::Matrix<numnod_, 1>& epreaf,        // current fluid pressure
    CORE::LINALG::Matrix<numdof_, (numdim_ + 1) * numnod_>*
        stiffmatrix,                                                    // element stiffness matrix
    CORE::LINALG::Matrix<numdof_, (numdim_ + 1) * numnod_>* reamatrix,  // element reactive matrix
    CORE::LINALG::Matrix<numdof_, 1>* force,  // element internal force vector
    Teuchos::ParameterList& params)           // algorithmic parameters e.g. time
{
  //=============================get parameters

  GetMaterials();

  //=======================================================================

  // update element geometry
  static CORE::LINALG::Matrix<numdim_, numnod_> xrefe;  // material coord. of element
  static CORE::LINALG::Matrix<numdim_, numnod_> xcurr;  // current  coord. of element

  DRT::Node** nodes = Nodes();
  for (int i = 0; i < numnod_; ++i)
  {
    const auto& x = nodes[i]->X();
    for (int j = 0; j < numdim_; j++)
    {
      xrefe(j, i) = x[j];
      xcurr(j, i) = xrefe(j, i) + disp(j, i);
    }
  }

  /* =========================================================================*/
  /* ================================================= Loop over Gauss Points */
  /* =========================================================================*/
  if (stiffmatrix != nullptr)
    GaussPointLoopOD(params, xrefe, xcurr, disp, vel, evelnp, epreaf, stiffmatrix);
}

template <class so3_ele, CORE::FE::CellType distype>
void DRT::ELEMENTS::So3_Poro<so3_ele, distype>::CouplingPoroelastPressureBased(
    std::vector<int>& lm,                          // location matrix
    CORE::LINALG::Matrix<numdim_, numnod_>& disp,  // current displacements
    const std::vector<double>& ephi,            // current primary variable for poro-multiphase flow
    CORE::LINALG::SerialDenseMatrix& couplmat,  // element stiffness matrix
    Teuchos::ParameterList& params)             // algorithmic parameters e.g. time
{
  GetMaterialsPressureBased();

  //=======================================================================

  // update element geometry
  CORE::LINALG::Matrix<numdim_, numnod_> xrefe;  // material coord. of element
  CORE::LINALG::Matrix<numdim_, numnod_> xcurr;  // current  coord. of element

  DRT::Node** nodes = Nodes();
  for (int i = 0; i < numnod_; ++i)
  {
    const auto& x = nodes[i]->X();
    for (int j = 0; j < numdim_; j++)
    {
      xrefe(j, i) = x[j];
      xcurr(j, i) = xrefe(j, i) + disp(j, i);
    }
  }

  /* =========================================================================*/
  /* ================================================= Loop over Gauss Points */
  /* =========================================================================*/

  GaussPointLoopODPressureBased(params, xrefe, xcurr, disp, ephi, couplmat);
}

template <class so3_ele, CORE::FE::CellType distype>
void DRT::ELEMENTS::So3_Poro<so3_ele, distype>::GaussPointLoopOD(Teuchos::ParameterList& params,
    const CORE::LINALG::Matrix<numdim_, numnod_>& xrefe,
    const CORE::LINALG::Matrix<numdim_, numnod_>& xcurr,
    const CORE::LINALG::Matrix<numdim_, numnod_>& nodaldisp,
    const CORE::LINALG::Matrix<numdim_, numnod_>& nodalvel,
    const CORE::LINALG::Matrix<numdim_, numnod_>& evelnp,
    const CORE::LINALG::Matrix<numnod_, 1>& epreaf,
    CORE::LINALG::Matrix<numdof_, (numdim_ + 1) * numnod_>* stiffmatrix)
{
  static CORE::LINALG::Matrix<numdim_, numnod_>
      N_XYZ;  //  first derivatives at gausspoint w.r.t. X, Y,Z
  // build deformation gradient wrt to material configuration
  // in case of prestressing, build defgrd wrt to last stored configuration
  // CAUTION: defgrd(true): filled with zeros!
  static CORE::LINALG::Matrix<numdim_, numdim_> defgrd(
      true);  //  deformation gradiant evaluated at gauss point
  static CORE::LINALG::Matrix<numnod_, 1> shapefct;  //  shape functions evalulated at gauss point
  static CORE::LINALG::Matrix<numdim_, numnod_> deriv(
      true);  //  first derivatives at gausspoint w.r.t. r,s,t

  for (int gp = 0; gp < numgpt_; ++gp)
  {
    // evaluate shape functions and derivatives at integration point
    ComputeShapeFunctionsAndDerivatives(gp, shapefct, deriv, N_XYZ);
    // evaluate second derivatives of shape functions at integration point
    // ComputeSecondDerivativesOfShapeFunctions(gp,xrefe,deriv,deriv2,N_XYZ,N_XYZ2);

    // (material) deformation gradient F = d xcurr / d xrefe = xcurr * N_XYZ^T
    ComputeDefGradient(defgrd, N_XYZ, xcurr);

    // inverse deformation gradient F^-1
    static CORE::LINALG::Matrix<numdim_, numdim_> defgrd_inv(false);
    defgrd_inv.Invert(defgrd);

    // jacobian determinant of transformation between spatial and material space "|dx/dX|"
    double J = 0.0;
    // volume change (used for porosity law). Same as J in nonlinear theory.
    double volchange = 0.0;

    // compute J and the volume change
    ComputeJacobianDeterminantVolumeChange(J, volchange, defgrd, N_XYZ, nodaldisp);

    // non-linear B-operator
    static CORE::LINALG::Matrix<numstr_, numdof_> bop;
    ComputeBOperator(bop, defgrd, N_XYZ);

    // -----------------Right Cauchy-Green tensor = F^T * F
    static CORE::LINALG::Matrix<numdim_, numdim_> cauchygreen;
    cauchygreen.MultiplyTN(defgrd, defgrd);

    //------------------ inverse Right Cauchy-Green tensor
    static CORE::LINALG::Matrix<numdim_, numdim_> C_inv(false);
    C_inv.Invert(cauchygreen);

    //---------------- get pressure at integration point
    double press = shapefct.Dot(epreaf);

    //------------------ get material pressure gradient at integration point
    static CORE::LINALG::Matrix<numdim_, 1> Gradp;
    Gradp.Multiply(N_XYZ, epreaf);

    //--------------------- get fluid velocity at integration point
    static CORE::LINALG::Matrix<numdim_, 1> fvelint;
    fvelint.Multiply(evelnp, shapefct);

    // material fluid velocity gradient at integration point
    static CORE::LINALG::Matrix<numdim_, numdim_> fvelder;
    fvelder.MultiplyNT(evelnp, N_XYZ);

    //----------------structure velocity at integration point
    static CORE::LINALG::Matrix<numdim_, 1> velint;
    velint.Multiply(nodalvel, shapefct);

    //**************************************************+auxilary variables for computing the
    // porosity and linearization
    double dphi_dp = 0.0;
    double porosity = 0.0;

    ComputePorosityAndLinearizationOD(
        params, press, volchange, gp, shapefct, nullptr, porosity, dphi_dp);

    // **********************evaluate stiffness matrix and force vector+++++++++++++++++++++++++

    FillMatrixAndVectorsOD(gp, shapefct, N_XYZ, J, porosity, dphi_dp, velint, fvelint, defgrd_inv,
        Gradp, bop, C_inv, stiffmatrix);

    if (fluid_mat_->Type() == MAT::PAR::darcy_brinkman)
    {
      FillMatrixAndVectorsBrinkmanOD(
          gp, shapefct, N_XYZ, J, porosity, dphi_dp, fvelder, defgrd_inv, bop, C_inv, stiffmatrix);
    }
  }
}

template <class so3_ele, CORE::FE::CellType distype>
void DRT::ELEMENTS::So3_Poro<so3_ele, distype>::GaussPointLoopODPressureBased(
    Teuchos::ParameterList& params, const CORE::LINALG::Matrix<numdim_, numnod_>& xrefe,
    const CORE::LINALG::Matrix<numdim_, numnod_>& xcurr,
    const CORE::LINALG::Matrix<numdim_, numnod_>& nodaldisp, const std::vector<double>& ephi,
    CORE::LINALG::SerialDenseMatrix& couplmat)
{
  /* =========================================================================*/
  /* ================================================= Loop over Gauss Points */
  /* =========================================================================*/
  CORE::LINALG::Matrix<numdim_, numnod_> N_XYZ;  //  first derivatives at gausspoint w.r.t. X, Y,Z
  // build deformation gradient wrt to material configuration
  // in case of prestressing, build defgrd wrt to last stored configuration
  // CAUTION: defgrd(true): filled with zeros!
  CORE::LINALG::Matrix<numdim_, numdim_> defgrd(
      true);                                  //  deformation gradiant evaluated at gauss point
  CORE::LINALG::Matrix<numnod_, 1> shapefct;  //  shape functions evalulated at gauss point
  CORE::LINALG::Matrix<numdim_, numnod_> deriv(
      true);  //  first derivatives at gausspoint w.r.t. r,s,t

  // Initialize
  const int numfluidphases = fluidmulti_mat_->NumFluidPhases();
  const int totalnumdofpernode = fluidmulti_mat_->NumMat();
  const int numvolfrac = fluidmulti_mat_->NumVolFrac();
  const bool hasvolfracs = (totalnumdofpernode - numfluidphases);
  std::vector<double> phiAtGP(totalnumdofpernode);
  std::vector<double> solpressderiv(totalnumdofpernode);

  for (int gp = 0; gp < numgpt_; ++gp)
  {
    // evaluate shape functions and derivatives at integration point
    ComputeShapeFunctionsAndDerivatives(gp, shapefct, deriv, N_XYZ);

    // (material) deformation gradient F = d xcurr / d xrefe = xcurr * N_XYZ^T
    ComputeDefGradient(defgrd, N_XYZ, xcurr);

    // jacobian determinant of transformation between spatial and material space "|dx/dX|"
    double J = 0.0;
    // volume change (used for porosity law). Same as J in nonlinear theory.
    double volchange = 0.0;

    // compute J, the volume change and the respctive linearizations w.r.t. structure displacement
    ComputeJacobianDeterminantVolumeChange(J, volchange, defgrd, N_XYZ, nodaldisp);

    // non-linear B-operator (may so be called, meaning
    CORE::LINALG::Matrix<numstr_, numdof_> bop;
    ComputeBOperator(bop, defgrd, N_XYZ);

    // Right Cauchy-Green tensor = F^T * F
    CORE::LINALG::Matrix<numdim_, numdim_> cauchygreen;
    cauchygreen.MultiplyTN(defgrd, defgrd);

    // inverse Right Cauchy-Green tensor
    CORE::LINALG::Matrix<numdim_, numdim_> C_inv(false);
    C_inv.Invert(cauchygreen);

    // compute derivative of solid pressure w.r.t primary variable phi at node
    ComputePrimaryVariableAtGP(ephi, totalnumdofpernode, shapefct, phiAtGP);
    ComputeSolPressureDeriv(phiAtGP, numfluidphases, solpressderiv);
    // in case of volume fractions --> recalculate
    if (hasvolfracs)
    {
      double dphi_dp = 0.0;
      double porosity = 0.0;

      double press = ComputeSolPressureAtGP(totalnumdofpernode, numfluidphases, phiAtGP);

      ComputePorosityAndLinearizationOD(
          params, press, volchange, gp, shapefct, nullptr, porosity, dphi_dp);

      RecalculateSolPressureDeriv(
          phiAtGP, totalnumdofpernode, numfluidphases, numvolfrac, press, porosity, solpressderiv);
    }

    // **********************evaluate stiffness matrix and force vector+++++++++++++++++++++++++
    FillMatrixAndVectorsODPressureBased(
        gp, shapefct, N_XYZ, J, bop, C_inv, solpressderiv, couplmat);
  }
}

template <class so3_ele, CORE::FE::CellType distype>
void DRT::ELEMENTS::So3_Poro<so3_ele, distype>::CouplingStressPoroelast(
    CORE::LINALG::Matrix<numdim_, numnod_>& disp,    // current displacements
    CORE::LINALG::Matrix<numdim_, numnod_>& evelnp,  // current fluid velocities
    CORE::LINALG::Matrix<numnod_, 1>& epreaf,        // current fluid pressure
    CORE::LINALG::SerialDenseMatrix* elestress,      // stresses at GP
    CORE::LINALG::SerialDenseMatrix* elestrain,      // strains at GP
    Teuchos::ParameterList& params,                  // algorithmic parameters e.g. time
    const INPAR::STR::StressType iostress            // stress output option
)
{
  // update element geometry
  CORE::LINALG::Matrix<numdim_, numnod_> xrefe;  // material coord. of element
  CORE::LINALG::Matrix<numdim_, numnod_> xcurr;  // current  coord. of element

  DRT::Node** nodes = Nodes();
  for (int i = 0; i < numnod_; ++i)
  {
    const auto& x = nodes[i]->X();
    for (int j = 0; j < numdim_; j++)
    {
      xrefe(j, i) = x[j];
      xcurr(j, i) = xrefe(j, i) + disp(j, i);
    }
  }

  // get structure material
  Teuchos::RCP<MAT::StructPoro> structmat = Teuchos::rcp_dynamic_cast<MAT::StructPoro>(Material());
  if (structmat->MaterialType() != INPAR::MAT::m_structporo)
    dserror("invalid structure material for poroelasticity");

  CORE::LINALG::Matrix<numnod_, 1> shapefct;
  CORE::LINALG::Matrix<numdim_, numdim_> defgrd(true);
  CORE::LINALG::Matrix<numdim_, numnod_> N_XYZ;
  CORE::LINALG::Matrix<numdim_, numnod_> deriv;

  for (int gp = 0; gp < numgpt_; ++gp)
  {
    // evaluate shape functions and derivatives at integration point
    ComputeShapeFunctionsAndDerivatives(gp, shapefct, deriv, N_XYZ);

    // (material) deformation gradient F = d xcurr / d xrefe = xcurr * N_XYZ^T
    ComputeDefGradient(defgrd, N_XYZ, xcurr);

    //----------------------------------------------------
    // pressure at integration point
    double press = shapefct.Dot(epreaf);

    CORE::LINALG::Matrix<numstr_, 1> couplstress(true);

    structmat->CouplStress(defgrd, press, couplstress);

    // return gp stresses
    switch (iostress)
    {
      case INPAR::STR::stress_2pk:
      {
        if (elestress == nullptr) dserror("stress data not available");
        for (int i = 0; i < numstr_; ++i) (*elestress)(gp, i) = couplstress(i);
      }
      break;
      case INPAR::STR::stress_cauchy:
      {
        if (elestress == nullptr) dserror("stress data not available");

        // push forward of material stress to the spatial configuration
        CORE::LINALG::Matrix<numdim_, numdim_> cauchycouplstress;
        PK2toCauchy(couplstress, defgrd, cauchycouplstress);

        (*elestress)(gp, 0) = cauchycouplstress(0, 0);
        (*elestress)(gp, 1) = cauchycouplstress(1, 1);
        (*elestress)(gp, 2) = cauchycouplstress(2, 2);
        (*elestress)(gp, 3) = cauchycouplstress(0, 1);
        (*elestress)(gp, 4) = cauchycouplstress(1, 2);
        (*elestress)(gp, 5) = cauchycouplstress(0, 2);
      }
      break;
      case INPAR::STR::stress_none:
        break;

      default:
        dserror("requested stress type not available");
        break;
    }
  }
}

template <class so3_ele, CORE::FE::CellType distype>
void DRT::ELEMENTS::So3_Poro<so3_ele, distype>::InitElement()
{
  CORE::LINALG::Matrix<numdim_, numnod_> deriv;
  CORE::LINALG::Matrix<numnod_, numdim_> xrefe;
  for (int i = 0; i < numnod_; ++i)
  {
    Node** nodes = Nodes();
    if (!nodes) dserror("Nodes() returned null pointer");
    xrefe(i, 0) = Nodes()[i]->X()[0];
    xrefe(i, 1) = Nodes()[i]->X()[1];
    xrefe(i, 2) = Nodes()[i]->X()[2];
  }

  if (distype == CORE::FE::CellType::nurbs27) isNurbs_ = true;

  invJ_.resize(numgpt_);
  detJ_.resize(numgpt_);
  xsi_.resize(numgpt_);

  for (int gp = 0; gp < numgpt_; ++gp)
  {
    const double* gpcoord = intpoints_.Point(gp);
    for (int idim = 0; idim < numdim_; idim++)
    {
      xsi_[gp](idim) = gpcoord[idim];
    }

    if (not isNurbs_)
    {
      CORE::DRT::UTILS::shape_function_deriv1<distype>(xsi_[gp], deriv);

      invJ_[gp].Multiply(deriv, xrefe);
      detJ_[gp] = invJ_[gp].Invert();
      if (detJ_[gp] <= 0.0) dserror("Element Jacobian mapping %10.5e <= 0.0", detJ_[gp]);
    }
  }

  init_ = true;

  scatra_coupling_ = false;

  ProblemType probtype = DRT::Problem::Instance()->GetProblemType();
  if (probtype == ProblemType::poroscatra) scatra_coupling_ = true;
}

template <class so3_ele, CORE::FE::CellType distype>
void DRT::ELEMENTS::So3_Poro<so3_ele, distype>::PK2toCauchy(
    CORE::LINALG::Matrix<numstr_, 1>& stress, CORE::LINALG::Matrix<numdim_, numdim_>& defgrd,
    CORE::LINALG::Matrix<numdim_, numdim_>& cauchystress)
{
  // calculate the Jacobi-determinant
  const double detF = (defgrd).Determinant();

  // sigma = 1/J . F . S . F^T
  CORE::LINALG::Matrix<numdim_, numdim_> pkstress;
  pkstress(0, 0) = (stress)(0);
  pkstress(0, 1) = (stress)(3);
  pkstress(0, 2) = (stress)(5);
  pkstress(1, 0) = pkstress(0, 1);
  pkstress(1, 1) = (stress)(1);
  pkstress(1, 2) = (stress)(4);
  pkstress(2, 0) = pkstress(0, 2);
  pkstress(2, 1) = pkstress(1, 2);
  pkstress(2, 2) = (stress)(2);

  CORE::LINALG::Matrix<numdim_, numdim_> temp;
  temp.Multiply((1.0 / detF), (defgrd), pkstress);
  (cauchystress).MultiplyNT(temp, (defgrd));
}

template <class so3_ele, CORE::FE::CellType distype>
void DRT::ELEMENTS::So3_Poro<so3_ele, distype>::ComputePorosityAndLinearization(
    Teuchos::ParameterList& params, const double& press, const double& J, const int& gp,
    const CORE::LINALG::Matrix<numnod_, 1>& shapfct,
    const CORE::LINALG::Matrix<numnod_, 1>* myporosity,
    const CORE::LINALG::Matrix<1, numdof_>& dJ_dus, double& porosity,
    CORE::LINALG::Matrix<1, numdof_>& dphi_dus)
{
  double dphi_dJ = 0.0;

  struct_mat_->ComputePorosity(params, press, J, gp, porosity,
      nullptr,  // dphi_dp not needed
      &dphi_dJ,
      nullptr,  // dphi_dJdp not needed
      nullptr,  // dphi_dJJ not needed
      nullptr   // dphi_dpp not needed
  );

  dphi_dus.Update(dphi_dJ, dJ_dus);
}

template <class so3_ele, CORE::FE::CellType distype>
void DRT::ELEMENTS::So3_Poro<so3_ele, distype>::ComputePorosityAndLinearizationOD(
    Teuchos::ParameterList& params, const double& press, const double& J, const int& gp,
    const CORE::LINALG::Matrix<numnod_, 1>& shapfct,
    const CORE::LINALG::Matrix<numnod_, 1>* myporosity, double& porosity, double& dphi_dp)
{
  struct_mat_->ComputePorosity(params, press, J, gp, porosity, &dphi_dp,
      nullptr,  // dphi_dJ not needed
      nullptr,  // dphi_dJdp not needed
      nullptr,  // dphi_dJJ not needed
      nullptr   // dphi_dpp not needed
  );
}

template <class so3_ele, CORE::FE::CellType distype>
void DRT::ELEMENTS::So3_Poro<so3_ele, distype>::ExtractValuesFromGlobalVector(
    const DRT::Discretization& discretization, const int& dofset, const std::vector<int>& lm,
    CORE::LINALG::Matrix<numdim_, numnod_>* matrixtofill,
    CORE::LINALG::Matrix<numnod_, 1>* vectortofill, const std::string& state)
{
  // get state of the global vector
  Teuchos::RCP<const Epetra_Vector> matrix_state = discretization.GetState(dofset, state);
  if (matrix_state == Teuchos::null) dserror("Cannot get state vector %s", state.c_str());

  // ask for the number of dofs of dofset
  const int numdofpernode = discretization.NumDof(dofset, Nodes()[0]);

  // extract local values of the global vectors
  std::vector<double> mymatrix(lm.size());
  DRT::UTILS::ExtractMyValues(*matrix_state, mymatrix, lm);

  if (numdofpernode == numdim_ + 1)
  {
    for (int inode = 0; inode < numnod_; ++inode)  // number of nodes
    {
      // fill a vector field via a pointer
      if (matrixtofill != nullptr)
      {
        for (int idim = 0; idim < numdim_; ++idim)  // number of dimensions
        {
          (*matrixtofill)(idim, inode) = mymatrix[idim + (inode * numdofpernode)];
        }
      }
      // fill a scalar field via a pointer
      if (vectortofill != nullptr)
        (*vectortofill)(inode, 0) = mymatrix[numdim_ + (inode * numdofpernode)];
    }
  }
  else if (numdofpernode == numdim_)
  {
    for (int inode = 0; inode < numnod_; ++inode)  // number of nodes
    {
      // fill a vector field via a pointer
      if (matrixtofill != nullptr)
      {
        for (int idim = 0; idim < numdim_; ++idim)  // number of dimensions
        {
          (*matrixtofill)(idim, inode) = mymatrix[idim + (inode * numdofpernode)];
        }
      }
    }
  }
  else if (numdofpernode == 1)
  {
    for (int inode = 0; inode < numnod_; ++inode)  // number of nodes
    {
      if (vectortofill != nullptr) (*vectortofill)(inode, 0) = mymatrix[inode * numdofpernode];
    }
  }
  else
  {
    for (int inode = 0; inode < numnod_; ++inode)  // number of nodes
    {
      if (vectortofill != nullptr) (*vectortofill)(inode, 0) = mymatrix[inode * numdofpernode];
    }
  }
}

template <class so3_ele, CORE::FE::CellType distype>
void DRT::ELEMENTS::So3_Poro<so3_ele, distype>::ComputeSolPressureDeriv(
    const std::vector<double>& phiAtGP, const int numfluidphases,
    std::vector<double>& solidpressderiv)
{
  // zero out everything
  std::fill(solidpressderiv.begin(), solidpressderiv.end(), 0.0);

  // initialize auxiliary variables
  std::vector<double> genpress(numfluidphases);
  std::vector<double> press(numfluidphases);
  std::vector<double> sat(numfluidphases);
  CORE::LINALG::SerialDenseMatrix helpderiv(numfluidphases, numfluidphases, true);
  CORE::LINALG::SerialDenseMatrix satderiv(numfluidphases, numfluidphases, true);
  CORE::LINALG::SerialDenseMatrix pressderiv(numfluidphases, numfluidphases, true);
  std::vector<double> fluidphi(phiAtGP.data(), phiAtGP.data() + numfluidphases);

  // evaluate the pressures
  fluidmulti_mat_->EvaluateGenPressure(genpress, fluidphi);

  // transform generalized pressures to true pressure values
  fluidmulti_mat_->TransformGenPresToTruePres(genpress, press);

  // explicit evaluation of saturation
  fluidmulti_mat_->EvaluateSaturation(sat, fluidphi, press);

  // calculate the derivative of the pressure (actually first its inverse)
  fluidmulti_mat_->EvaluateDerivOfDofWrtPressure(pressderiv, fluidphi);

  // now invert the derivatives of the dofs w.r.t. pressure to get the derivatives
  // of the pressure w.r.t. the dofs
  {
    using ordinalType = CORE::LINALG::SerialDenseMatrix::ordinalType;
    using scalarType = CORE::LINALG::SerialDenseMatrix::scalarType;
    Teuchos::SerialDenseSolver<ordinalType, scalarType> inverse;
    inverse.setMatrix(Teuchos::rcpFromRef(pressderiv));
    int err = inverse.invert();
    if (err != 0)
      dserror("Inversion of matrix for pressure derivative failed with error code %d.", err);
  }

  // calculate derivatives of saturation w.r.t. pressure
  fluidmulti_mat_->EvaluateDerivOfSaturationWrtPressure(helpderiv, press);

  // chain rule: the derivative of saturation w.r.t. dof =
  // (derivative of saturation w.r.t. pressure) * (derivative of pressure w.r.t. dof)
  CORE::LINALG::multiply(satderiv, helpderiv, pressderiv);

  // compute derivative of solid pressure w.r.t. dofs with product rule
  // standard derivative: no volume fractions present
  for (int iphase = 0; iphase < numfluidphases; iphase++)
  {
    for (int jphase = 0; jphase < numfluidphases; jphase++)
      solidpressderiv[iphase] +=
          pressderiv(jphase, iphase) * sat[jphase] + satderiv(jphase, iphase) * press[jphase];
  }
}

template <class so3_ele, CORE::FE::CellType distype>
void DRT::ELEMENTS::So3_Poro<so3_ele, distype>::ComputeLinearizationOfSolPressWrtDisp(
    const double fluidpress, const double porosity, const int totalnumdofpernode,
    const int numfluidphases, const int numvolfrac, const std::vector<double>& phiAtGP,
    const CORE::LINALG::Matrix<1, numdof_>& dphi_dus, CORE::LINALG::Matrix<1, numdof_>& dps_dus)
{
  // get volume fraction primary variables
  std::vector<double> volfracphi(
      phiAtGP.data() + numfluidphases, phiAtGP.data() + numfluidphases + numvolfrac);
  double sumaddvolfrac = 0.0;
  for (int ivolfrac = 0; ivolfrac < numvolfrac; ivolfrac++) sumaddvolfrac += volfracphi[ivolfrac];

  // get volume fraction pressure at [numfluidphases+numvolfrac...totalnumdofpernode-1]
  std::vector<double> volfracpressure(
      phiAtGP.data() + numfluidphases + numvolfrac, phiAtGP.data() + totalnumdofpernode);

  // p_s = (porosity - sumaddvolfrac)/porosity * fluidpress
  //       + 1.0 / porosity sum_i=1^numvolfrac (volfrac_i*pressure_i)
  // d (p_s) / d porosity = + sumaddvolfrac/porosity/porosity * fluidpress
  double dps_dphi = sumaddvolfrac / (porosity * porosity) * fluidpress;

  // ... + 1.0 / porosity / porosity sum_i=1^numvolfrac (volfrac_i*pressure_i)
  for (int ivolfrac = 0; ivolfrac < numvolfrac; ivolfrac++)
    dps_dphi -= volfracphi[ivolfrac] * volfracpressure[ivolfrac] / (porosity * porosity);

  // d (p_s) / d u_s = d (p_s) / d porosity * d porosity / d u_s
  dps_dus.Update(dps_dphi, dphi_dus);
}

template <class so3_ele, CORE::FE::CellType distype>
void DRT::ELEMENTS::So3_Poro<so3_ele, distype>::RecalculateSolPressureDeriv(
    const std::vector<double>& phiAtGP, const int totalnumdofpernode, const int numfluidphases,
    const int numvolfrac, const double press, const double porosity,
    std::vector<double>& solidpressderiv)
{
  // get volume fraction primary variables
  std::vector<double> volfracphi(
      phiAtGP.data() + numfluidphases, phiAtGP.data() + numfluidphases + numvolfrac);
  double sumaddvolfrac = 0.0;
  for (int ivolfrac = 0; ivolfrac < numvolfrac; ivolfrac++) sumaddvolfrac += volfracphi[ivolfrac];

  // p_s = (porosity - sumaddvolfrac)/porosity * fluidpress
  //      + 1.0 / porosity sum_i=1^numvolfrac (volfrac_i*pressure_i)
  const double scale = (porosity - sumaddvolfrac) / porosity;

  // scale original fluid press deriv with (porosity - sumaddvolfrac)/porosity
  for (int iphase = 0; iphase < numfluidphases; iphase++) solidpressderiv[iphase] *= scale;

  // get volfrac pressures at [numfluidphases+numvolfrac...totalnumdofpernode-1]
  std::vector<double> volfracpressure(
      phiAtGP.data() + numfluidphases + numvolfrac, phiAtGP.data() + totalnumdofpernode);

  for (int ivolfrac = 0; ivolfrac < numvolfrac; ivolfrac++)
  {
    // d p_s / d volfrac = - fluidpress/porosity + volfracpressure/porosity
    solidpressderiv[ivolfrac + numfluidphases] =
        -1.0 / porosity * press + 1.0 / porosity * volfracpressure[ivolfrac];
    // d p_s / d volfracpress = + volfracphi/porosity
    solidpressderiv[ivolfrac + numfluidphases + numvolfrac] = volfracphi[ivolfrac] / porosity;
  }
}

template <class so3_ele, CORE::FE::CellType distype>
double DRT::ELEMENTS::So3_Poro<so3_ele, distype>::ComputeSolPressureAtGP(
    const int totalnumdofpernode, const int numfluidphases, const std::vector<double>& phiAtGP)
{
  // initialize auxiliary variables
  std::vector<double> genpress(numfluidphases, 0.0);
  std::vector<double> sat(numfluidphases, 0.0);
  std::vector<double> press(numfluidphases, 0.0);
  std::vector<double> fluidphi(phiAtGP.data(), phiAtGP.data() + numfluidphases);

  // evaluate the pressures
  fluidmulti_mat_->EvaluateGenPressure(genpress, fluidphi);

  //! transform generalized pressures to true pressure values
  fluidmulti_mat_->TransformGenPresToTruePres(genpress, press);

  // explicit evaluation of saturation
  fluidmulti_mat_->EvaluateSaturation(sat, fluidphi, press);

  // solid pressure = sum (S_i*p_i)
  const double solidpressure = std::inner_product(sat.begin(), sat.end(), press.begin(), 0.0);

  return solidpressure;
}

template <class so3_ele, CORE::FE::CellType distype>
double DRT::ELEMENTS::So3_Poro<so3_ele, distype>::RecalculateSolPressureAtGP(double press,
    const double porosity, const int totalnumdofpernode, const int numfluidphases,
    const int numvolfrac, const std::vector<double>& phiAtGP)
{
  // get volume fraction primary variables at [numfluidphases-1...numfluidphase-1+numvolfrac]
  std::vector<double> volfracphi(
      phiAtGP.data() + numfluidphases, phiAtGP.data() + numfluidphases + numvolfrac);
  double sumaddvolfrac = 0.0;
  for (int ivolfrac = 0; ivolfrac < numvolfrac; ivolfrac++) sumaddvolfrac += volfracphi[ivolfrac];

  // p_s = (porosity - sumaddvolfrac)/porosity * fluidpress
  //      + 1.0 / porosity * sum_i=1^numvolfrac (volfrac_i*pressure_i)
  // first part
  press *= (porosity - sumaddvolfrac) / porosity;

  // get volfrac pressures at [numfluidphases+numvolfrac...totalnumdofpernode-1]
  std::vector<double> volfracpressure(
      phiAtGP.data() + numfluidphases + numvolfrac, phiAtGP.data() + totalnumdofpernode);

  // second part
  for (int ivolfrac = 0; ivolfrac < numvolfrac; ivolfrac++)
    press += volfracphi[ivolfrac] / porosity * volfracpressure[ivolfrac];

  // note: in RecalculateSolidPressure in porofluid_phasemanager calculation is performed a bit
  //       differently since we already pass porosity = porosity - sumaddvolfrac, but result is
  //       equivalent

  return press;
}

template <class so3_ele, CORE::FE::CellType distype>
void DRT::ELEMENTS::So3_Poro<so3_ele, distype>::ComputePrimaryVariableAtGP(
    const std::vector<double>& ephi, const int totalnumdofpernode,
    const CORE::LINALG::Matrix<numnod_, 1>& shapefct, std::vector<double>& phiAtGP)
{
  // zero out everything
  std::fill(phiAtGP.begin(), phiAtGP.end(), 0.0);
  // compute phi at GP = phi * shapefunction
  for (int i = 0; i < numnod_; i++)
  {
    for (int j = 0; j < totalnumdofpernode; j++)
    {
      phiAtGP[j] += shapefct(i) * ephi[i * totalnumdofpernode + j];
    }
  }
}

template <class so3_ele, CORE::FE::CellType distype>
void DRT::ELEMENTS::So3_Poro<so3_ele, distype>::GetMaterials()
{
  // get structure material
  if (struct_mat_ == Teuchos::null)
  {
    struct_mat_ = Teuchos::rcp_dynamic_cast<MAT::StructPoro>(Material());
    if (struct_mat_->MaterialType() != INPAR::MAT::m_structporo and
        struct_mat_->MaterialType() != INPAR::MAT::m_structpororeaction and
        struct_mat_->MaterialType() != INPAR::MAT::m_structpororeactionECM)
      dserror("invalid structure material for poroelasticity");
  }

  // get fluid material
  if (fluid_mat_ == Teuchos::null)
  {
    // access second material in structure element
    if (so3_ele::NumMaterial() > 1)
    {
      fluid_mat_ = Teuchos::rcp_dynamic_cast<MAT::FluidPoro>(so3_ele::Material(1));
      if (fluid_mat_->MaterialType() != INPAR::MAT::m_fluidporo)
        dserror("invalid fluid material for poroelasticity");
    }
    else
      dserror("no second material defined for element %i", Id());
  }
}

template <class so3_ele, CORE::FE::CellType distype>
void DRT::ELEMENTS::So3_Poro<so3_ele, distype>::GetMaterialsPressureBased()
{
  // get structure material
  if (struct_mat_ == Teuchos::null)
  {
    struct_mat_ = Teuchos::rcp_dynamic_cast<MAT::StructPoro>(Material());
    if (struct_mat_ == Teuchos::null) dserror("cast to poro material failed");

    if (struct_mat_->MaterialType() != INPAR::MAT::m_structporo and
        struct_mat_->MaterialType() != INPAR::MAT::m_structpororeaction and
        struct_mat_->MaterialType() != INPAR::MAT::m_structpororeactionECM)
      dserror("invalid structure material for poroelasticity");
  }

  // Get Fluid-multiphase-Material
  if (fluidmulti_mat_ == Teuchos::null)
  {
    // access second material in structure element
    if (so3_ele::NumMaterial() > 1)
    {
      fluidmulti_mat_ = Teuchos::rcp_dynamic_cast<MAT::FluidPoroMultiPhase>(so3_ele::Material(1));
      if (fluidmulti_mat_ == Teuchos::null)
        dserror("cast to multiphase fluid poro material failed");
      if (fluidmulti_mat_->MaterialType() != INPAR::MAT::m_fluidporo_multiphase and
          fluidmulti_mat_->MaterialType() != INPAR::MAT::m_fluidporo_multiphase_reactions)
        dserror("invalid fluid material for poro-multiphase-elasticity");
      if (fluidmulti_mat_->NumFluidPhases() == 0)
      {
        dserror(
            "NUMFLUIDPHASES_IN_MULTIPHASEPORESPACE = 0 currently not supported since this requires "
            "an adaption of the definition of the solid pressure");
      }
    }
    else
      dserror("no second material defined for element %i", Id());
  }
}

template <class so3_ele, CORE::FE::CellType distype>
void DRT::ELEMENTS::So3_Poro<so3_ele, distype>::ComputePorosity(Teuchos::ParameterList& params,
    double press, double J, int gp, double& porosity, double* dphi_dp, double* dphi_dJ,
    double* dphi_dJdp, double* dphi_dJJ, double* dphi_dpp, bool save)
{
  struct_mat_->ComputePorosity(
      params, press, J, gp, porosity, dphi_dp, dphi_dJ, dphi_dJdp, dphi_dJJ, dphi_dpp, save);
}

template <class so3_ele, CORE::FE::CellType distype>
void DRT::ELEMENTS::So3_Poro<so3_ele, distype>::ComputeSurfPorosity(Teuchos::ParameterList& params,
    double press, double J, int surfnum, int gp, double& porosity, double* dphi_dp, double* dphi_dJ,
    double* dphi_dJdp, double* dphi_dJJ, double* dphi_dpp, bool save)
{
  struct_mat_->ComputeSurfPorosity(params, press, J, surfnum, gp, porosity, dphi_dp, dphi_dJ,
      dphi_dJdp, dphi_dJJ, dphi_dpp, save);
}

template <class so3_ele, CORE::FE::CellType distype>
double DRT::ELEMENTS::So3_Poro<so3_ele, distype>::RefPorosityTimeDeriv()
{
  return struct_mat_->RefPorosityTimeDeriv();
}

template <class so3_ele, CORE::FE::CellType distype>
void DRT::ELEMENTS::So3_Poro<so3_ele, distype>::ComputeShapeFunctionsAndDerivatives(const int& gp,
    CORE::LINALG::Matrix<numnod_, 1>& shapefct, CORE::LINALG::Matrix<numdim_, numnod_>& deriv,
    CORE::LINALG::Matrix<numdim_, numnod_>& N_XYZ)
{
  if (!isNurbs_)
  {
    CORE::DRT::UTILS::shape_function<distype>(xsi_[gp], shapefct);
    CORE::DRT::UTILS::shape_function_deriv1<distype>(xsi_[gp], deriv);
  }
  else
  {
    CORE::DRT::NURBS::UTILS::nurbs_get_funct_deriv(
        shapefct, deriv, xsi_[gp], myknots_, weights_, distype);

    CORE::LINALG::Matrix<numnod_, numdim_> xrefe;
    for (int i = 0; i < numnod_; ++i)
    {
      Node** nodes = Nodes();
      if (!nodes) dserror("Nodes() returned null pointer");
      xrefe(i, 0) = Nodes()[i]->X()[0];
      xrefe(i, 1) = Nodes()[i]->X()[1];
      xrefe(i, 2) = Nodes()[i]->X()[2];
    }

    invJ_[gp].Multiply(deriv, xrefe);
    detJ_[gp] = invJ_[gp].Invert();
    if (detJ_[gp] <= 0.0) dserror("Element Jacobian mapping %10.5e <= 0.0", detJ_[gp]);
  }

  /* get the inverse of the Jacobian matrix which looks like:
   **            [ X_,r  Y_,r  Z_,r ]^-1
   **     J^-1 = [ X_,s  Y_,s  Z_,s ]
   **            [ X_,t  Y_,t  Z_,t ]
   */

  // compute derivatives N_XYZ at gp w.r.t. material coordinates
  // by N_XYZ = J^-1 * N_rst
  N_XYZ.Multiply(invJ_[gp], deriv);  // (6.21)
}

template <class so3_ele, CORE::FE::CellType distype>
void DRT::ELEMENTS::So3_Poro<so3_ele, distype>::ComputeJacobianDeterminantVolumeChange(double& J,
    double& volchange, const CORE::LINALG::Matrix<numdim_, numdim_>& defgrd,
    const CORE::LINALG::Matrix<numdim_, numnod_>& N_XYZ,
    const CORE::LINALG::Matrix<numdim_, numnod_>& nodaldisp)
{
  // compute J
  J = defgrd.Determinant();

  if (so3_ele::kintype_ == INPAR::STR::kinem_nonlinearTotLag)  // total lagrange (nonlinear)
  {
    // for nonlinear kinematics the Jacobian of the deformation gradient is the volume change
    volchange = J;
  }
  else if (so3_ele::kintype_ == INPAR::STR::kinem_linear)  // linear kinematics
  {
    // for linear kinematics the volume change is the trace of the linearized strains

    // gradient of displacements
    static CORE::LINALG::Matrix<numdim_, numdim_> dispgrad;
    dispgrad.Clear();
    // gradient of displacements
    dispgrad.MultiplyNT(nodaldisp, N_XYZ);

    volchange = 1.0;
    // volchange = 1 + trace of the linearized strains (= trace of displacement gradient)
    for (int i = 0; i < numdim_; ++i) volchange += dispgrad(i, i);
  }
  else
    dserror("invalid kinematic type!");
}

template <class so3_ele, CORE::FE::CellType distype>
void DRT::ELEMENTS::So3_Poro<so3_ele,
    distype>::ComputeJacobianDeterminantVolumeChangeAndLinearizations(double& J, double& volchange,
    CORE::LINALG::Matrix<1, numdof_>& dJ_dus, CORE::LINALG::Matrix<1, numdof_>& dvolchange_dus,
    const CORE::LINALG::Matrix<numdim_, numdim_>& defgrd,
    const CORE::LINALG::Matrix<numdim_, numdim_>& defgrd_inv,
    const CORE::LINALG::Matrix<numdim_, numnod_>& N_XYZ,
    const CORE::LINALG::Matrix<numdim_, numnod_>& nodaldisp)
{
  // compute J
  J = defgrd.Determinant();
  // compute linearization of J
  ComputeLinearizationOfJacobian(dJ_dus, J, N_XYZ, defgrd_inv);

  if (so3_ele::kintype_ == INPAR::STR::kinem_nonlinearTotLag)  // total lagrange (nonlinear)
  {
    // for nonlinear kinematics the Jacobian of the deformation gradient is the volume change
    volchange = J;
    dvolchange_dus = dJ_dus;
  }
  else if (so3_ele::kintype_ == INPAR::STR::kinem_linear)  // linear kinematics
  {
    // for linear kinematics the volume change is the trace of the linearized strains

    // gradient of displacements
    static CORE::LINALG::Matrix<numdim_, numdim_> dispgrad;
    dispgrad.Clear();
    // gradient of displacements
    dispgrad.MultiplyNT(nodaldisp, N_XYZ);

    volchange = 1.0;
    // volchange = 1 + trace of the linearized strains (= trace of displacement gradient)
    for (int i = 0; i < numdim_; ++i) volchange += dispgrad(i, i);

    for (int i = 0; i < numdim_; ++i)
      for (int j = 0; j < numnod_; ++j) dvolchange_dus(numdim_ * j + i) = N_XYZ(i, j);
  }
  else
    dserror("invalid kinematic type!");
}

template <class so3_ele, CORE::FE::CellType distype>
void DRT::ELEMENTS::So3_Poro<so3_ele, distype>::ComputeAuxiliaryValues(
    const CORE::LINALG::Matrix<numdim_, numnod_>& N_XYZ,
    const CORE::LINALG::Matrix<numdim_, numdim_>& defgrd_inv,
    const CORE::LINALG::Matrix<numdim_, numdim_>& C_inv,
    const CORE::LINALG::Matrix<numdim_, 1>& Gradp,
    CORE::LINALG::Matrix<numdim_ * numdim_, numdof_>& dFinvTdus,
    CORE::LINALG::Matrix<numdim_, 1>& Finvgradp,
    CORE::LINALG::Matrix<numdim_, numdof_>& dFinvdus_gradp,
    CORE::LINALG::Matrix<numstr_, numdof_>& dCinv_dus)
{
  // F^-T * Grad p
  Finvgradp.MultiplyTN(defgrd_inv, Gradp);

  if (so3_ele::kintype_ != INPAR::STR::kinem_linear)
  {
    // dF^-T/dus
    dFinvTdus.Clear();
    for (int i = 0; i < numdim_; i++)
    {
      for (int n = 0; n < numnod_; n++)
      {
        for (int j = 0; j < numdim_; j++)
        {
          const int gid = numdim_ * n + j;
          for (int k = 0; k < numdim_; k++)
            for (int l = 0; l < numdim_; l++)
              dFinvTdus(i * numdim_ + l, gid) += -defgrd_inv(l, j) * N_XYZ(k, n) * defgrd_inv(k, i);
        }
      }
    }

    // dF^-T/dus * Grad p
    dFinvdus_gradp.Clear();
    for (int i = 0; i < numdim_; i++)
    {
      for (int n = 0; n < numnod_; n++)
      {
        for (int j = 0; j < numdim_; j++)
        {
          const int gid = numdim_ * n + j;
          for (int l = 0; l < numdim_; l++)
            dFinvdus_gradp(i, gid) += dFinvTdus(i * numdim_ + l, gid) * Gradp(l);
        }
      }
    }
  }

  dCinv_dus.Clear();
  for (int n = 0; n < numnod_; ++n)
  {
    for (int k = 0; k < numdim_; ++k)
    {
      const int gid = n * numdim_ + k;
      for (int i = 0; i < numdim_; ++i)
      {
        dCinv_dus(0, gid) += -2 * C_inv(0, i) * N_XYZ(i, n) * defgrd_inv(0, k);
        dCinv_dus(1, gid) += -2 * C_inv(1, i) * N_XYZ(i, n) * defgrd_inv(1, k);
        dCinv_dus(2, gid) += -2 * C_inv(2, i) * N_XYZ(i, n) * defgrd_inv(2, k);
        /* ~~~ */
        dCinv_dus(3, gid) += -C_inv(0, i) * N_XYZ(i, n) * defgrd_inv(1, k) -
                             defgrd_inv(0, k) * N_XYZ(i, n) * C_inv(1, i);
        dCinv_dus(4, gid) += -C_inv(1, i) * N_XYZ(i, n) * defgrd_inv(2, k) -
                             defgrd_inv(1, k) * N_XYZ(i, n) * C_inv(2, i);
        dCinv_dus(5, gid) += -C_inv(2, i) * N_XYZ(i, n) * defgrd_inv(0, k) -
                             defgrd_inv(2, k) * N_XYZ(i, n) * C_inv(0, i);
      }
    }
  }
}

template <class so3_ele, CORE::FE::CellType distype>
inline void DRT::ELEMENTS::So3_Poro<so3_ele, distype>::ComputeBOperator(
    CORE::LINALG::Matrix<numstr_, numdof_>& bop,
    const CORE::LINALG::Matrix<numdim_, numdim_>& defgrd,
    const CORE::LINALG::Matrix<numdim_, numnod_>& N_XYZ)
{
  /* non-linear B-operator (may so be called, meaning
   ** of B-operator is not so sharp in the non-linear realm) *
   ** B = F . Bl *
   **
   **      [ ... | F_11*N_{,1}^k  F_21*N_{,1}^k  F_31*N_{,1}^k | ... ]
   **      [ ... | F_12*N_{,2}^k  F_22*N_{,2}^k  F_32*N_{,2}^k | ... ]
   **      [ ... | F_13*N_{,3}^k  F_23*N_{,3}^k  F_33*N_{,3}^k | ... ]
   ** B =  [ ~~~   ~~~~~~~~~~~~~  ~~~~~~~~~~~~~  ~~~~~~~~~~~~~   ~~~ ]
   **      [       F_11*N_{,2}^k+F_12*N_{,1}^k                       ]
   **      [ ... |          F_21*N_{,2}^k+F_22*N_{,1}^k        | ... ]
   **      [                       F_31*N_{,2}^k+F_32*N_{,1}^k       ]
   **      [                                                         ]
   **      [       F_12*N_{,3}^k+F_13*N_{,2}^k                       ]
   **      [ ... |          F_22*N_{,3}^k+F_23*N_{,2}^k        | ... ]
   **      [                       F_32*N_{,3}^k+F_33*N_{,2}^k       ]
   **      [                                                         ]
   **      [       F_13*N_{,1}^k+F_11*N_{,3}^k                       ]
   **      [ ... |          F_23*N_{,1}^k+F_21*N_{,3}^k        | ... ]
   **      [                       F_33*N_{,1}^k+F_31*N_{,3}^k       ]
   */
  for (int i = 0; i < numnod_; ++i)
  {
    bop(0, noddof_ * i + 0) = defgrd(0, 0) * N_XYZ(0, i);
    bop(0, noddof_ * i + 1) = defgrd(1, 0) * N_XYZ(0, i);
    bop(0, noddof_ * i + 2) = defgrd(2, 0) * N_XYZ(0, i);
    bop(1, noddof_ * i + 0) = defgrd(0, 1) * N_XYZ(1, i);
    bop(1, noddof_ * i + 1) = defgrd(1, 1) * N_XYZ(1, i);
    bop(1, noddof_ * i + 2) = defgrd(2, 1) * N_XYZ(1, i);
    bop(2, noddof_ * i + 0) = defgrd(0, 2) * N_XYZ(2, i);
    bop(2, noddof_ * i + 1) = defgrd(1, 2) * N_XYZ(2, i);
    bop(2, noddof_ * i + 2) = defgrd(2, 2) * N_XYZ(2, i);
    /* ~~~ */
    bop(3, noddof_ * i + 0) = defgrd(0, 0) * N_XYZ(1, i) + defgrd(0, 1) * N_XYZ(0, i);
    bop(3, noddof_ * i + 1) = defgrd(1, 0) * N_XYZ(1, i) + defgrd(1, 1) * N_XYZ(0, i);
    bop(3, noddof_ * i + 2) = defgrd(2, 0) * N_XYZ(1, i) + defgrd(2, 1) * N_XYZ(0, i);
    bop(4, noddof_ * i + 0) = defgrd(0, 1) * N_XYZ(2, i) + defgrd(0, 2) * N_XYZ(1, i);
    bop(4, noddof_ * i + 1) = defgrd(1, 1) * N_XYZ(2, i) + defgrd(1, 2) * N_XYZ(1, i);
    bop(4, noddof_ * i + 2) = defgrd(2, 1) * N_XYZ(2, i) + defgrd(2, 2) * N_XYZ(1, i);
    bop(5, noddof_ * i + 0) = defgrd(0, 2) * N_XYZ(0, i) + defgrd(0, 0) * N_XYZ(2, i);
    bop(5, noddof_ * i + 1) = defgrd(1, 2) * N_XYZ(0, i) + defgrd(1, 0) * N_XYZ(2, i);
    bop(5, noddof_ * i + 2) = defgrd(2, 2) * N_XYZ(0, i) + defgrd(2, 0) * N_XYZ(2, i);
  }
}

template <class so3_ele, CORE::FE::CellType distype>
inline void DRT::ELEMENTS::So3_Poro<so3_ele, distype>::ComputeLinearizationOfJacobian(
    CORE::LINALG::Matrix<1, numdof_>& dJ_dus, const double& J,
    const CORE::LINALG::Matrix<numdim_, numnod_>& N_XYZ,
    const CORE::LINALG::Matrix<numdim_, numdim_>& defgrd_inv)
{
  if (so3_ele::kintype_ == INPAR::STR::kinem_nonlinearTotLag)  // total lagrange (nonlinear)
  {
    //------------------------------------ build F^-1 as vector 9x1
    CORE::LINALG::Matrix<numdim_ * numdim_, 1> defgrd_inv_vec;
    defgrd_inv_vec(0) = defgrd_inv(0, 0);
    defgrd_inv_vec(1) = defgrd_inv(0, 1);
    defgrd_inv_vec(2) = defgrd_inv(0, 2);
    defgrd_inv_vec(3) = defgrd_inv(1, 0);
    defgrd_inv_vec(4) = defgrd_inv(1, 1);
    defgrd_inv_vec(5) = defgrd_inv(1, 2);
    defgrd_inv_vec(6) = defgrd_inv(2, 0);
    defgrd_inv_vec(7) = defgrd_inv(2, 1);
    defgrd_inv_vec(8) = defgrd_inv(2, 2);

    //--------------------------- build N_X operator (wrt material config)
    CORE::LINALG::Matrix<9, numdof_> N_X(true);  // set to zero
    for (int i = 0; i < numnod_; ++i)
    {
      N_X(0, 3 * i + 0) = N_XYZ(0, i);
      N_X(1, 3 * i + 1) = N_XYZ(0, i);
      N_X(2, 3 * i + 2) = N_XYZ(0, i);

      N_X(3, 3 * i + 0) = N_XYZ(1, i);
      N_X(4, 3 * i + 1) = N_XYZ(1, i);
      N_X(5, 3 * i + 2) = N_XYZ(1, i);

      N_X(6, 3 * i + 0) = N_XYZ(2, i);
      N_X(7, 3 * i + 1) = N_XYZ(2, i);
      N_X(8, 3 * i + 2) = N_XYZ(2, i);
    }

    //------linearization of jacobi determinant detF=J w.r.t. strucuture displacement   dJ/d(us) =
    // dJ/dF : dF/dus = J * F^-T * N,X
    dJ_dus.MultiplyTN(J, defgrd_inv_vec, N_X);
  }
  else if (so3_ele::kintype_ == INPAR::STR::kinem_linear)  // linear kinematics
  {
    // J=1 -> no linearization
    dJ_dus.Clear();
  }
  else
    dserror("invalid kinematic type!");
}

template <class so3_ele, CORE::FE::CellType distype>
void DRT::ELEMENTS::So3_Poro<so3_ele, distype>::FillMatrixAndVectors(const int& gp,
    const CORE::LINALG::Matrix<numnod_, 1>& shapefct,
    const CORE::LINALG::Matrix<numdim_, numnod_>& N_XYZ, const double& J, const double& press,
    const double& porosity, const CORE::LINALG::Matrix<numdim_, 1>& velint,
    const CORE::LINALG::Matrix<numdim_, 1>& fvelint,
    const CORE::LINALG::Matrix<numdim_, numdim_>& fvelder,
    const CORE::LINALG::Matrix<numdim_, numdim_>& defgrd_inv,
    const CORE::LINALG::Matrix<numstr_, numdof_>& bop,
    const CORE::LINALG::Matrix<numdim_, numdim_>& C_inv,
    const CORE::LINALG::Matrix<numdim_, 1>& Finvgradp,
    const CORE::LINALG::Matrix<1, numdof_>& dphi_dus,
    const CORE::LINALG::Matrix<1, numdof_>& dJ_dus,
    const CORE::LINALG::Matrix<numstr_, numdof_>& dCinv_dus,
    const CORE::LINALG::Matrix<numdim_, numdof_>& dFinvdus_gradp,
    const CORE::LINALG::Matrix<numdim_ * numdim_, numdof_>& dFinvTdus,
    CORE::LINALG::Matrix<numdof_, numdof_>& erea_v,
    CORE::LINALG::Matrix<numdof_, numdof_>* stiffmatrix, CORE::LINALG::Matrix<numdof_, 1>* force,
    CORE::LINALG::Matrix<numstr_, 1>& fstress)
{
  const double detJ_w = detJ_[gp] * intpoints_.Weight(gp);

  {
    static CORE::LINALG::Matrix<numdim_, numdim_> matreatensor(true);
    static CORE::LINALG::Matrix<numdim_, numdim_> reatensor(true);
    static CORE::LINALG::Matrix<numdim_, numdim_> linreac_dphi(true);
    static CORE::LINALG::Matrix<numdim_, numdim_> linreac_dJ(true);
    static CORE::LINALG::Matrix<numdim_, 1> reafvel(true);
    static CORE::LINALG::Matrix<numdim_, 1> reavel(true);
    {
      static CORE::LINALG::Matrix<numdim_, numdim_> temp(true);
      std::vector<double> anisotropic_permeability_coeffs =
          ComputeAnisotropicPermeabilityCoeffsAtGP(shapefct);
      fluid_mat_->ComputeReactionTensor(matreatensor, J, porosity,
          anisotropic_permeability_directions_, anisotropic_permeability_coeffs);
      fluid_mat_->ComputeLinMatReactionTensor(linreac_dphi, linreac_dJ, J, porosity);
      temp.Multiply(1.0, matreatensor, defgrd_inv);
      reatensor.MultiplyTN(defgrd_inv, temp);
      reavel.Multiply(reatensor, velint);
      reafvel.Multiply(reatensor, fvelint);
    }

    for (int idim = 0; idim < numdim_; idim++)
    {
      const double reafvel_idim = reafvel(idim);
      const double reac_vel_idim = reavel(idim);
      const double Finvgradp_idim = Finvgradp(idim);

      for (int inode = 0; inode < numnod_; inode++)
      {
        const double fac = detJ_w * shapefct(inode);
        const double v = fac * porosity * porosity * J * J;
        const int fk = numdim_ * inode;

        /*-------structure- fluid velocity coupling:  RHS
         "darcy-terms"
         - reacoeff * J^2 *  phi^2 *  v^f
         */
        (*force)(fk + idim) += -v * reafvel_idim;

        /* "reactive darcy-terms"
         reacoeff * J^2 *  phi^2 *  v^s
         */
        (*force)(fk + idim) += v * reac_vel_idim;

        /*-------structure- fluid pressure coupling: RHS
         *                        "pressure gradient terms"
         - J *  F^-T * Grad(p) * phi
         */
        (*force)(fk + idim) += fac * J * Finvgradp_idim * (-porosity);
      }
    }

    for (int idim = 0; idim < numdim_; idim++)
    {
      for (int jdim = 0; jdim < numdim_; jdim++)
      {
        const double reatensor_i_j = reatensor(idim, jdim);

        for (int inode = 0; inode < numnod_; inode++)
        {
          const int fk = numdim_ * inode;
          const double v = detJ_w * shapefct(inode) * porosity * porosity * J * J;

          for (int jnode = 0; jnode < numnod_; jnode++)
          {
            const int fi = numdim_ * jnode;

            /* additional "reactive darcy-term"
             detJ * w(gp) * ( J^2 * reacoeff * phi^2  ) * D(v_s)
             */
            erea_v(fk + idim, fi + jdim) += v * reatensor_i_j * shapefct(jnode);
          }
        }
      }
    }

    for (int idim = 0; idim < numdim_; idim++)
    {
      const double Finvgradp_j = Finvgradp(idim);

      for (int jdim = 0; jdim < numdim_; jdim++)
      {
        for (int jnode = 0; jnode < numnod_; jnode++)
        {
          const int fi = numdim_ * jnode;

          const double val = detJ_w * (-porosity * dJ_dus(fi + jdim) * Finvgradp_j -
                                          porosity * J * dFinvdus_gradp(idim, fi + jdim) -
                                          dphi_dus(fi + jdim) * J * Finvgradp_j);

          for (int inode = 0; inode < numnod_; inode++)
          {
            /* additional "pressure gradient term"
             -  detJ * w(gp) * phi *  ( dJ/d(us) * F^-T * Grad(p) - J * d(F^-T)/d(us) *Grad(p) ) *
             D(us)
             - detJ * w(gp) * d(phi)/d(us) * J * F^-T * Grad(p) * D(us)
             */
            (*stiffmatrix)(numdim_ * inode + idim, fi + jdim) += shapefct(inode) * val;
          }
        }
      }
    }

    for (int idim = 0; idim < numdim_; idim++)
    {
      const double reac_vel_j = reavel(idim);
      const double reafvel_j = reafvel(idim);

      for (int jdim = 0; jdim < numdim_; jdim++)
      {
        for (int jnode = 0; jnode < numnod_; jnode++)
        {
          const int fi = numdim_ * jnode;
          const double val = detJ_w * J * porosity * 2 * (reac_vel_j - reafvel_j) *
                             (porosity * dJ_dus(fi + jdim) + J * dphi_dus(fi + jdim));

          for (int inode = 0; inode < numnod_; inode++)
          {
            /* additional "reactive darcy-term"
               detJ * w(gp) * 2 * ( dJ/d(us) * vs * reacoeff * phi^2 + J * reacoeff * phi *
             d(phi)/d(us) * vs ) * D(us)
             - detJ * w(gp) *  2 * ( J * dJ/d(us) * v^f * reacoeff * phi^2 + J * reacoeff * phi *
             d(phi)/d(us) * v^f ) * D(us)
             */
            (*stiffmatrix)(numdim_ * inode + idim, fi + jdim) += shapefct(inode) * val;
          }
        }
      }
    }

    // check if derivatives of reaction tensor are zero --> significant speed up
    if (fluid_mat_->PermeabilityFunction() == MAT::PAR::constant)
    {
      const double fac = detJ_w * porosity * porosity * J * J;
      for (int idim = 0; idim < numdim_; idim++)
      {
        for (int jdim = 0; jdim < numdim_; jdim++)
        {
          for (int jnode = 0; jnode < numnod_; jnode++)
          {
            const int fi = numdim_ * jnode;

            for (int inode = 0; inode < numnod_; inode++)
            {
              double val = 0.0;
              for (int p = 0; p < numdim_; ++p)
              {
                const double velint_p = velint(p);
                const double fvelint_p = fvelint(p);
                for (int n = 0; n < numdim_; ++n)
                {
                  const double defgrd_inv_n_p = defgrd_inv(n, p);
                  const double dFinvTdus_n_p = dFinvTdus(p * numdim_ + n, fi + jdim);
                  for (int m = 0; m < numdim_; ++m)
                  {
                    val += fac * (velint_p - fvelint_p) *
                           (dFinvTdus(idim * numdim_ + m, fi + jdim) * matreatensor(m, n) *
                                   defgrd_inv_n_p +
                               defgrd_inv(m, idim) * matreatensor(m, n) * dFinvTdus_n_p);
                  }
                }
              }

              (*stiffmatrix)(numdim_ * inode + idim, fi + jdim) += shapefct(inode) * val;
            }
          }
        }
      }
    }
    else
    {
      const double fac = detJ_w * porosity * porosity * J * J;
      for (int idim = 0; idim < numdim_; idim++)
      {
        for (int jdim = 0; jdim < numdim_; jdim++)
        {
          for (int jnode = 0; jnode < numnod_; jnode++)
          {
            const int fi = numdim_ * jnode;
            const double dphi_dus_fi_l = dphi_dus(fi + jdim);
            const double dJ_dus_fi_l = dJ_dus(fi + jdim);

            for (int inode = 0; inode < numnod_; inode++)
            {
              double val = 0.0;
              for (int m = 0; m < numdim_; ++m)
              {
                const double dFinvTdus_idim_m_fi_jdim = dFinvTdus(idim * numdim_ + m, fi + jdim);
                const double defgrd_inv_m_idim = defgrd_inv(m, idim);
                for (int n = 0; n < numdim_; ++n)
                {
                  const double matreatensor_m_n = matreatensor(m, n);
                  const double linreac_dphi_m_n = linreac_dphi(m, n);
                  const double linreac_dJ_m_n = linreac_dJ(m, n);

                  for (int p = 0; p < numdim_; ++p)
                  {
                    val +=
                        fac * (velint(p) - fvelint(p)) *
                        (dFinvTdus_idim_m_fi_jdim * matreatensor_m_n * defgrd_inv(n, p) +
                            defgrd_inv_m_idim * matreatensor_m_n *
                                dFinvTdus(p * numdim_ + n, fi + jdim) +
                            defgrd_inv_m_idim *
                                (linreac_dphi_m_n * dphi_dus_fi_l + linreac_dJ_m_n * dJ_dus_fi_l) *
                                defgrd_inv(n, p));
                  }
                }
              }
              (*stiffmatrix)(numdim_ * inode + idim, fi + jdim) += val * shapefct(inode);
            }
          }
        }
      }
    }

    // inverse Right Cauchy-Green tensor as vector
    static CORE::LINALG::Matrix<numstr_, 1> C_inv_vec;
    for (int i = 0, k = 0; i < numdim_; i++)
      for (int j = 0; j < numdim_ - i; j++, k++) C_inv_vec(k) = C_inv(i + j, j);

    // B^T . C^-1
    static CORE::LINALG::Matrix<numdof_, 1> cinvb(true);
    cinvb.MultiplyTN(bop, C_inv_vec);

    const double fac1 = -detJ_w * press;
    const double fac2 = fac1 * J;

    // additional fluid stress term -(B^T . C^-1 * J * p^f * detJ * w(gp))
    force->Update(fac2, cinvb, 1.0);

    static CORE::LINALG::Matrix<numdof_, numdof_> tmp1;
    static CORE::LINALG::Matrix<numdof_, numdof_> tmp2;

    tmp1.Multiply(fac1, cinvb, dJ_dus);
    tmp2.MultiplyTN(fac2, bop, dCinv_dus);

    // additional fluid stress- stiffness term -(B^T . C^-1 . dJ/d(us) * p^f * detJ * w(gp))
    stiffmatrix->Update(1.0, tmp1, 1.0);

    // additional fluid stress- stiffness term -(B^T .  dC^-1/d(us) * J * p^f * detJ * w(gp))
    stiffmatrix->Update(1.0, tmp2, 1.0);

    // integrate `geometric' stiffness matrix and add to keu *****************
    CORE::LINALG::Matrix<numstr_, 1> sfac(C_inv_vec);  // auxiliary integrated stress

    // scale and add viscous stress
    sfac.Update(detJ_w, fstress, fac2);  // detJ*w(gp)*[S11,S22,S33,S12=S21,S23=S32,S13=S31]

    std::vector<double> SmB_L(3);  // intermediate Sm.B_L
    // kgeo += (B_L^T . sigma . B_L) * detJ * w(gp)  with B_L = Ni,Xj see NiliFEM-Skript
    for (int inod = 0; inod < numnod_; ++inod)
    {
      SmB_L[0] = sfac(0) * N_XYZ(0, inod) + sfac(3) * N_XYZ(1, inod) + sfac(5) * N_XYZ(2, inod);
      SmB_L[1] = sfac(3) * N_XYZ(0, inod) + sfac(1) * N_XYZ(1, inod) + sfac(4) * N_XYZ(2, inod);
      SmB_L[2] = sfac(5) * N_XYZ(0, inod) + sfac(4) * N_XYZ(1, inod) + sfac(2) * N_XYZ(2, inod);
      for (int jnod = 0; jnod < numnod_; ++jnod)
      {
        double bopstrbop = 0.0;  // intermediate value
        for (int idim = 0; idim < numdim_; ++idim) bopstrbop += N_XYZ(idim, jnod) * SmB_L[idim];
        (*stiffmatrix)(numdim_ * inod + 0, numdim_ * jnod + 0) += bopstrbop;
        (*stiffmatrix)(numdim_ * inod + 1, numdim_ * jnod + 1) += bopstrbop;
        (*stiffmatrix)(numdim_ * inod + 2, numdim_ * jnod + 2) += bopstrbop;
      }
    }
  }
}

template <class so3_ele, CORE::FE::CellType distype>
void DRT::ELEMENTS::So3_Poro<so3_ele, distype>::FillMatrixAndVectorsPressureBased(const int& gp,
    const CORE::LINALG::Matrix<numnod_, 1>& shapefct,
    const CORE::LINALG::Matrix<numdim_, numnod_>& N_XYZ, const double& J, const double& press,
    const CORE::LINALG::Matrix<numstr_, numdof_>& bop,
    const CORE::LINALG::Matrix<numdim_, numdim_>& C_inv,
    const CORE::LINALG::Matrix<1, numdof_>& dJ_dus,
    const CORE::LINALG::Matrix<numstr_, numdof_>& dCinv_dus,
    const CORE::LINALG::Matrix<1, numdof_>& dps_dus,
    CORE::LINALG::Matrix<numdof_, numdof_>* stiffmatrix, CORE::LINALG::Matrix<numdof_, 1>* force)
{
  const double detJ_w = detJ_[gp] * intpoints_.Weight(gp);

  //-----------inverse Right Cauchy-Green tensor as vector in voigt notation
  static CORE::LINALG::Matrix<numstr_, 1> C_inv_vec(true);
  for (int i = 0, k = 0; i < numdim_; i++)
    for (int j = 0; j < numdim_ - i; j++, k++) C_inv_vec(k) = C_inv(i + j, j);

  // B^T . C^-1
  static CORE::LINALG::Matrix<numdof_, 1> cinvb(true);
  cinvb.MultiplyTN(bop, C_inv_vec);

  const double fac1 = -detJ_w * press;
  const double fac2 = fac1 * J;

  // update internal force vector
  if (force != nullptr)
  {
    // additional fluid stress- stiffness term RHS -(B^T .  C^-1  * J * p^f * detJ * w(gp))
    force->Update(fac2, cinvb, 1.0);
  }

  // update stiffness matrix
  if (stiffmatrix != nullptr)
  {
    static CORE::LINALG::Matrix<numdof_, numdof_> tmp;

    // additional fluid stress- stiffness term -(B^T . C^-1 . dJ/d(us) * p^f * detJ * w(gp))
    tmp.Multiply(fac1, cinvb, dJ_dus);
    stiffmatrix->Update(1.0, tmp, 1.0);

    // additional fluid stress- stiffness term -(B^T .  dC^-1/d(us) * J * p^f * detJ * w(gp))
    tmp.MultiplyTN(fac2, bop, dCinv_dus);
    stiffmatrix->Update(1.0, tmp, 1.0);

    // additional fluid stress- stiffness term -(B^T .  dC^-1 * J * dp^s/d(us) * detJ * w(gp))
    tmp.Multiply(-detJ_w * J, cinvb, dps_dus);
    stiffmatrix->Update(1.0, tmp, 1.0);

    // integrate `geometric' stiffness matrix and add to keu *****************
    CORE::LINALG::Matrix<numstr_, 1> sfac(C_inv_vec);  // auxiliary integrated stress

    // scale
    sfac.Scale(fac2);

    std::vector<double> SmB_L(3);  // intermediate Sm.B_L
    // kgeo += (B_L^T . sigma . B_L) * detJ * w(gp)  with B_L = Ni,Xj see NiliFEM-Skript
    for (int inod = 0; inod < numnod_; ++inod)
    {
      SmB_L[0] = sfac(0) * N_XYZ(0, inod) + sfac(3) * N_XYZ(1, inod) + sfac(5) * N_XYZ(2, inod);
      SmB_L[1] = sfac(3) * N_XYZ(0, inod) + sfac(1) * N_XYZ(1, inod) + sfac(4) * N_XYZ(2, inod);
      SmB_L[2] = sfac(5) * N_XYZ(0, inod) + sfac(4) * N_XYZ(1, inod) + sfac(2) * N_XYZ(2, inod);
      for (int jnod = 0; jnod < numnod_; ++jnod)
      {
        double bopstrbop = 0.0;  // intermediate value
        for (int idim = 0; idim < numdim_; ++idim) bopstrbop += N_XYZ(idim, jnod) * SmB_L[idim];
        (*stiffmatrix)(numdim_ * inod + 0, numdim_ * jnod + 0) += bopstrbop;
        (*stiffmatrix)(numdim_ * inod + 1, numdim_ * jnod + 1) += bopstrbop;
        (*stiffmatrix)(numdim_ * inod + 2, numdim_ * jnod + 2) += bopstrbop;
      }
    }
  }
}

template <class so3_ele, CORE::FE::CellType distype>
void DRT::ELEMENTS::So3_Poro<so3_ele, distype>::FillMatrixAndVectorsBrinkman(const int& gp,
    const double& J, const double& porosity, const CORE::LINALG::Matrix<numdim_, numdim_>& fvelder,
    const CORE::LINALG::Matrix<numdim_, numdim_>& defgrd_inv,
    const CORE::LINALG::Matrix<numstr_, numdof_>& bop,
    const CORE::LINALG::Matrix<numdim_, numdim_>& C_inv,
    const CORE::LINALG::Matrix<1, numdof_>& dphi_dus,
    const CORE::LINALG::Matrix<1, numdof_>& dJ_dus,
    const CORE::LINALG::Matrix<numstr_, numdof_>& dCinv_dus,
    const CORE::LINALG::Matrix<numdim_ * numdim_, numdof_>& dFinvTdus,
    CORE::LINALG::Matrix<numdof_, numdof_>* stiffmatrix, CORE::LINALG::Matrix<numdof_, 1>* force,
    CORE::LINALG::Matrix<numstr_, 1>& fstress)
{
  double detJ_w = detJ_[gp] * intpoints_.Weight(gp);

  double visc = fluid_mat_->Viscosity();
  CORE::LINALG::Matrix<numdim_, numdim_> CinvFvel;
  CORE::LINALG::Matrix<numdim_, numdim_> visctress1;
  CinvFvel.Multiply(C_inv, fvelder);
  visctress1.MultiplyNT(CinvFvel, defgrd_inv);
  CORE::LINALG::Matrix<numdim_, numdim_> visctress2(visctress1);
  visctress1.UpdateT(1.0, visctress2, 1.0);

  fstress(0) = visctress1(0, 0);
  fstress(1) = visctress1(1, 1);
  fstress(2) = visctress1(2, 2);
  fstress(3) = visctress1(0, 1);
  fstress(4) = visctress1(1, 2);
  fstress(5) = visctress1(2, 0);

  fstress.Scale(detJ_w * visc * J * porosity);

  // B^T . C^-1
  static CORE::LINALG::Matrix<numdof_, 1> fstressb(true);
  fstressb.MultiplyTN(bop, fstress);

  force->Update(1.0, fstressb, 1.0);

  // evaluate viscous terms (for darcy-brinkman flow only)
  {
    static CORE::LINALG::Matrix<numdim_, numdim_> tmp;
    tmp.MultiplyNT(fvelder, defgrd_inv);

    double fac = detJ_w * visc;

    CORE::LINALG::Matrix<numstr_, numdof_> fstress_dus(true);
    {
      const double tmp_0_0 = tmp(0, 0);
      const double tmp_0_1 = tmp(0, 1);
      const double tmp_0_2 = tmp(0, 2);
      const double tmp_1_0 = tmp(1, 0);
      const double tmp_1_1 = tmp(1, 1);
      const double tmp_1_2 = tmp(1, 2);
      const double tmp_2_0 = tmp(2, 0);
      const double tmp_2_1 = tmp(2, 1);
      const double tmp_2_2 = tmp(2, 2);

      const double CinvFvel_0_0 = CinvFvel(0, 0);
      const double CinvFvel_0_1 = CinvFvel(0, 1);
      const double CinvFvel_0_2 = CinvFvel(0, 2);
      const double CinvFvel_1_0 = CinvFvel(1, 0);
      const double CinvFvel_1_1 = CinvFvel(1, 1);
      const double CinvFvel_1_2 = CinvFvel(1, 2);
      const double CinvFvel_2_0 = CinvFvel(2, 0);
      const double CinvFvel_2_1 = CinvFvel(2, 1);
      const double CinvFvel_2_2 = CinvFvel(2, 2);

      for (int n = 0; n < numnod_; ++n)
      {
        for (int k = 0; k < numdim_; ++k)
        {
          const int gid = n * numdim_ + k;

          fstress_dus(0, gid) += 2 * (dCinv_dus(0, gid) * tmp_0_0 + dCinv_dus(3, gid) * tmp_1_0 +
                                         dCinv_dus(5, gid) * tmp_2_0);
          fstress_dus(1, gid) += 2 * (dCinv_dus(3, gid) * tmp_0_1 + dCinv_dus(1, gid) * tmp_1_1 +
                                         dCinv_dus(4, gid) * tmp_2_1);
          fstress_dus(2, gid) += 2 * (dCinv_dus(5, gid) * tmp_0_2 + dCinv_dus(4, gid) * tmp_1_2 +
                                         dCinv_dus(2, gid) * tmp_2_2);
          /* ~~~ */
          fstress_dus(3, gid) += +dCinv_dus(0, gid) * tmp_0_1 + dCinv_dus(3, gid) * tmp_1_1 +
                                 dCinv_dus(5, gid) * tmp_2_1 + dCinv_dus(3, gid) * tmp_0_0 +
                                 dCinv_dus(1, gid) * tmp_1_0 + dCinv_dus(4, gid) * tmp_2_0;
          fstress_dus(4, gid) += +dCinv_dus(3, gid) * tmp_0_2 + dCinv_dus(1, gid) * tmp_1_2 +
                                 dCinv_dus(4, gid) * tmp_2_2 + dCinv_dus(5, gid) * tmp_0_1 +
                                 dCinv_dus(4, gid) * tmp_1_1 + dCinv_dus(2, gid) * tmp_2_1;
          fstress_dus(5, gid) += +dCinv_dus(5, gid) * tmp_0_0 + dCinv_dus(4, gid) * tmp_1_0 +
                                 dCinv_dus(2, gid) * tmp_2_0 + dCinv_dus(0, gid) * tmp_0_2 +
                                 dCinv_dus(3, gid) * tmp_1_2 + dCinv_dus(5, gid) * tmp_2_2;

          fstress_dus(0, gid) += 2 * CinvFvel_0_0 * dFinvTdus(0 * numdim_, gid) +
                                 2 * CinvFvel_0_1 * dFinvTdus(1 * numdim_, gid) +
                                 2 * CinvFvel_0_2 * dFinvTdus(2 * numdim_, gid);
          fstress_dus(1, gid) += 2 * CinvFvel_1_0 * dFinvTdus(0 * numdim_ + 1, gid) +
                                 2 * CinvFvel_1_1 * dFinvTdus(1 * numdim_ + 1, gid) +
                                 2 * CinvFvel_1_2 * dFinvTdus(2 * numdim_ + 1, gid);
          fstress_dus(2, gid) += 2 * CinvFvel_2_0 * dFinvTdus(0 * numdim_ + 2, gid) +
                                 2 * CinvFvel_2_1 * dFinvTdus(1 * numdim_ + 2, gid) +
                                 2 * CinvFvel_2_2 * dFinvTdus(2 * numdim_ + 2, gid);
          /* ~~~ */
          fstress_dus(3, gid) += +CinvFvel_0_0 * dFinvTdus(0 * numdim_ + 1, gid) +
                                 CinvFvel_1_0 * dFinvTdus(0 * numdim_, gid) +
                                 CinvFvel_0_1 * dFinvTdus(1 * numdim_ + 1, gid) +
                                 CinvFvel_1_1 * dFinvTdus(1 * numdim_, gid) +
                                 CinvFvel_0_2 * dFinvTdus(2 * numdim_ + 1, gid) +
                                 CinvFvel_1_2 * dFinvTdus(2 * numdim_, gid);
          fstress_dus(4, gid) += +CinvFvel_1_0 * dFinvTdus(0 * numdim_ + 2, gid) +
                                 CinvFvel_2_0 * dFinvTdus(0 * numdim_ + 1, gid) +
                                 CinvFvel_1_1 * dFinvTdus(1 * numdim_ + 2, gid) +
                                 CinvFvel_2_1 * dFinvTdus(1 * numdim_ + 1, gid) +
                                 CinvFvel_1_2 * dFinvTdus(2 * numdim_ + 2, gid) +
                                 CinvFvel_2_2 * dFinvTdus(2 * numdim_ + 1, gid);
          fstress_dus(5, gid) += +CinvFvel_2_0 * dFinvTdus(0 * numdim_, gid) +
                                 CinvFvel_0_0 * dFinvTdus(0 * numdim_ + 2, gid) +
                                 CinvFvel_2_1 * dFinvTdus(1 * numdim_, gid) +
                                 CinvFvel_0_1 * dFinvTdus(1 * numdim_ + 2, gid) +
                                 CinvFvel_2_2 * dFinvTdus(2 * numdim_, gid) +
                                 CinvFvel_0_2 * dFinvTdus(2 * numdim_ + 2, gid);
        }
      }
    }

    static CORE::LINALG::Matrix<numdof_, numdof_> fluidstress_part;

    // additional viscous fluid stress- stiffness term (B^T . fstress . dJ/d(us) * porosity * detJ *
    // w(gp))
    fluidstress_part.Multiply(fac * porosity, fstressb, dJ_dus);
    stiffmatrix->Update(1.0, fluidstress_part, 1.0);

    // additional fluid stress- stiffness term (B^T .  d\phi/d(us) . fstress  * J * w(gp))
    fluidstress_part.Multiply(fac * J, fstressb, dphi_dus);
    stiffmatrix->Update(1.0, fluidstress_part, 1.0);

    // additional fluid stress- stiffness term (B^T .  phi . dfstress/d(us)  * J * w(gp))
    fluidstress_part.MultiplyTN(detJ_w * visc * J * porosity, bop, fstress_dus);
    stiffmatrix->Update(1.0, fluidstress_part, 1.0);
  }
}

template <class so3_ele, CORE::FE::CellType distype>
void DRT::ELEMENTS::So3_Poro<so3_ele, distype>::FillMatrixAndVectorsOD(const int& gp,
    const CORE::LINALG::Matrix<numnod_, 1>& shapefct,
    const CORE::LINALG::Matrix<numdim_, numnod_>& N_XYZ, const double& J, const double& porosity,
    const double& dphi_dp, const CORE::LINALG::Matrix<numdim_, 1>& velint,
    const CORE::LINALG::Matrix<numdim_, 1>& fvelint,
    const CORE::LINALG::Matrix<numdim_, numdim_>& defgrd_inv,
    const CORE::LINALG::Matrix<numdim_, 1>& Gradp,
    const CORE::LINALG::Matrix<numstr_, numdof_>& bop,
    const CORE::LINALG::Matrix<numdim_, numdim_>& C_inv,
    CORE::LINALG::Matrix<numdof_, (numdim_ + 1) * numnod_>* stiffmatrix)
{
  double detJ_w = detJ_[gp] * intpoints_.Weight(gp);

  static CORE::LINALG::Matrix<numdim_, numdim_> matreatensor(true);
  static CORE::LINALG::Matrix<numdim_, numdim_> reatensor(true);
  static CORE::LINALG::Matrix<numdim_, numdim_> linreac_dphi(true);
  static CORE::LINALG::Matrix<numdim_, numdim_> linreac_dJ(true);
  static CORE::LINALG::Matrix<numdim_, 1> reafvel(true);
  static CORE::LINALG::Matrix<numdim_, 1> reavel(true);
  {
    CORE::LINALG::Matrix<numdim_, numdim_> temp(true);
    std::vector<double> anisotropic_permeability_coeffs =
        ComputeAnisotropicPermeabilityCoeffsAtGP(shapefct);
    fluid_mat_->ComputeReactionTensor(matreatensor, J, porosity,
        anisotropic_permeability_directions_, anisotropic_permeability_coeffs);
    fluid_mat_->ComputeLinMatReactionTensor(linreac_dphi, linreac_dJ, J, porosity);
    temp.Multiply(1.0, matreatensor, defgrd_inv);
    reatensor.MultiplyTN(defgrd_inv, temp);
    reavel.Multiply(reatensor, velint);
    reafvel.Multiply(reatensor, fvelint);
  }

  //-----------inverse Right Cauchy-Green tensor as vector in voigt notation
  static CORE::LINALG::Matrix<numstr_, 1> C_inv_vec(true);
  for (int i = 0, k = 0; i < numdim_; i++)
    for (int j = 0; j < numdim_ - i; j++, k++) C_inv_vec(k) = C_inv(i + j, j);

  // B^T . C^-1
  static CORE::LINALG::Matrix<numdof_, 1> cinvb(true);
  cinvb.MultiplyTN(bop, C_inv_vec);

  // F^-T * grad p
  static CORE::LINALG::Matrix<numdim_, 1> Finvgradp;
  Finvgradp.MultiplyTN(defgrd_inv, Gradp);

  // F^-T * N_XYZ
  static CORE::LINALG::Matrix<numdim_, numnod_> FinvNXYZ;
  FinvNXYZ.MultiplyTN(defgrd_inv, N_XYZ);

  {
    const double fac = detJ_w * J * J * 2 * porosity * dphi_dp;
    for (int idim = 0; idim < numdim_; idim++)
    {
      const double reafvel_idim = reafvel(idim);
      const double reac_vel_idim = reavel(idim);

      for (int jnode = 0; jnode < numnod_; jnode++)
      {
        const int fkp1 = (numdim_ + 1) * jnode;

        const double val = fac * shapefct(jnode) * (reac_vel_idim - reafvel_idim);
        for (int inode = 0; inode < numnod_; inode++)
        {
          /*-------structure- fluid pressure coupling:  "dracy-terms" + "reactive darcy-terms"
           - 2 * reacoeff * J * v^f * phi * d(phi)/dp  Dp
           + 2 * reacoeff * J * v^s * phi * d(phi)/dp  Dp
           */
          (*stiffmatrix)(numdim_ * inode + idim, fkp1 + numdim_) += shapefct(inode) * val;
        }
      }
    }
  }

  {
    for (int idim = 0; idim < numdim_; idim++)
    {
      const double Finvgradp_idim = Finvgradp(idim);
      for (int jnode = 0; jnode < numnod_; jnode++)
      {
        const int fkp1 = (numdim_ + 1) * jnode;

        const double val1 = detJ_w * (-1.0) * J * shapefct(jnode);
        const double val2 =
            -1.0 * detJ_w * J *
            (Finvgradp_idim * dphi_dp * shapefct(jnode) + porosity * FinvNXYZ(idim, jnode));

        for (int inode = 0; inode < numnod_; inode++)
        {
          /*-------structure- fluid pressure coupling: "stress terms" + "pressure gradient terms"
           -B^T . ( -1*J*C^-1 ) * Dp
           - J * F^-T * dphi/dp * Dp - J * F^-T * d(Grad((p))/(dp) * phi * Dp
           */
          (*stiffmatrix)(numdim_ * inode + idim, fkp1 + numdim_) +=
              val1 * cinvb(numdim_ * inode + idim) + val2 * shapefct(inode);
        }
      }
    }
  }

  // check if derivatives of reaction tensor are zero --> significant speed up
  if (fluid_mat_->PermeabilityFunction() != MAT::PAR::constant)
  {
    const double fac = detJ_w * J * J * porosity * porosity * dphi_dp;
    for (int idim = 0; idim < numdim_; idim++)
    {
      for (int jnode = 0; jnode < numnod_; jnode++)
      {
        const int fkp1 = (numdim_ + 1) * jnode;
        const double shapefct_jnode = shapefct(jnode);

        for (int inode = 0; inode < numnod_; inode++)
        {
          double val = 0.0;
          for (int p = 0; p < numdim_; ++p)
          {
            const double velint_fvelint_p = velint(p) - fvelint(p);
            for (int n = 0; n < numdim_; ++n)
            {
              const double defgrd_inv_n_p = defgrd_inv(n, p);
              for (int m = 0; m < numdim_; ++m)
              {
                val += fac * defgrd_inv(m, idim) * linreac_dphi(m, n) * defgrd_inv_n_p *
                       velint_fvelint_p;
              }
            }
          }
          val *= shapefct_jnode;

          /*-------structure- fluid pressure coupling:   "reactive darcy-terms"
           + J * J * phi * phi * defgrd_^-T * d(mat_reacoeff)/d(phi) * defgrd_^-1 * (v^s-v^f) *
           d(phi)/dp Dp
           */
          (*stiffmatrix)(numdim_ * inode + idim, fkp1 + numdim_) += shapefct(inode) * val;
        }
      }
    }
  }

  {
    const double fac = detJ_w * J * J * porosity * porosity;
    for (int idim = 0; idim < numdim_; idim++)
    {
      for (int jdim = 0; jdim < numdim_; jdim++)
      {
        const double reatensor_idim_jdim = reatensor(idim, jdim);
        for (int jnode = 0; jnode < numnod_; jnode++)
        {
          const double val = -1.0 * fac * shapefct(jnode) * reatensor_idim_jdim;

          /*-------structure- fluid velocity coupling:  "darcy-terms"
           -reacoeff * J * J *  phi^2 *  Dv^f
           */
          for (int inode = 0; inode < numnod_; inode++)
            (*stiffmatrix)(numdim_ * inode + idim, (numdim_ + 1) * jnode + jdim) +=
                val * shapefct(inode);
        }
      }
    }
  }
}

template <class so3_ele, CORE::FE::CellType distype>
void DRT::ELEMENTS::So3_Poro<so3_ele, distype>::FillMatrixAndVectorsODPressureBased(const int& gp,
    const CORE::LINALG::Matrix<numnod_, 1>& shapefct,
    const CORE::LINALG::Matrix<numdim_, numnod_>& N_XYZ, const double& J,
    const CORE::LINALG::Matrix<numstr_, numdof_>& bop,
    const CORE::LINALG::Matrix<numdim_, numdim_>& C_inv, const std::vector<double>& solpressderiv,
    CORE::LINALG::SerialDenseMatrix& couplmat)
{
  const double detJ_w = detJ_[gp] * intpoints_.Weight(gp);

  // inverse Right Cauchy-Green tensor as vector
  static CORE::LINALG::Matrix<numstr_, 1> C_inv_vec;
  for (int i = 0, k = 0; i < numdim_; i++)
    for (int j = 0; j < numdim_ - i; j++, k++) C_inv_vec(k) = C_inv(i + j, j);

  // B^T . C^-1
  static CORE::LINALG::Matrix<numdof_, 1> cinvb(true);
  cinvb.MultiplyTN(bop, C_inv_vec);

  const int totalnumdofpernode = fluidmulti_mat_->NumMat();

  {
    for (int i = 0; i < numnod_; i++)
    {
      const int fi = numdim_ * i;

      for (int j = 0; j < numdim_; j++)
      {
        for (int k = 0; k < numnod_; k++)
        {
          for (int iphase = 0; iphase < totalnumdofpernode; iphase++)
          {
            int fk_press = k * totalnumdofpernode + iphase;

            /*-------structure- fluid pressure coupling: "stress term"
             -B^T . ( -1*J*C^-1 ) * Dp
             */
            couplmat(fi + j, fk_press) +=
                detJ_w * cinvb(fi + j) * (-1.0) * J * shapefct(k) * solpressderiv[iphase];
          }
        }
      }
    }
  }
}

template <class so3_ele, CORE::FE::CellType distype>
void DRT::ELEMENTS::So3_Poro<so3_ele, distype>::FillMatrixAndVectorsBrinkmanOD(const int& gp,
    const CORE::LINALG::Matrix<numnod_, 1>& shapefct,
    const CORE::LINALG::Matrix<numdim_, numnod_>& N_XYZ, const double& J, const double& porosity,
    const double& dphi_dp, const CORE::LINALG::Matrix<numdim_, numdim_>& fvelder,
    const CORE::LINALG::Matrix<numdim_, numdim_>& defgrd_inv,
    const CORE::LINALG::Matrix<numstr_, numdof_>& bop,
    const CORE::LINALG::Matrix<numdim_, numdim_>& C_inv,
    CORE::LINALG::Matrix<numdof_, (numdim_ + 1) * numnod_>* stiffmatrix)
{
  double detJ_w = detJ_[gp] * intpoints_.Weight(gp);  // gpweights[gp];

  static CORE::LINALG::Matrix<numstr_, 1> fstress;

  double visc = fluid_mat_->Viscosity();
  static CORE::LINALG::Matrix<numdim_, numdim_> CinvFvel;
  static CORE::LINALG::Matrix<numdim_, numdim_> tmp;
  CinvFvel.Multiply(C_inv, fvelder);
  tmp.MultiplyNT(CinvFvel, defgrd_inv);
  static CORE::LINALG::Matrix<numdim_, numdim_> tmp2(tmp);
  tmp.UpdateT(1.0, tmp2, 1.0);

  fstress(0) = tmp(0, 0);
  fstress(1) = tmp(1, 1);
  fstress(2) = tmp(2, 2);
  fstress(3) = tmp(0, 1);
  fstress(4) = tmp(1, 2);
  fstress(5) = tmp(2, 0);

  // B^T . \sigma
  static CORE::LINALG::Matrix<numdof_, 1> fstressb;
  fstressb.MultiplyTN(bop, fstress);
  static CORE::LINALG::Matrix<numdim_, numnod_> N_XYZ_Finv;
  N_XYZ_Finv.Multiply(defgrd_inv, N_XYZ);

  // dfstress/dv^f
  CORE::LINALG::Matrix<numstr_, numdof_> dfstressb_dv;
  for (int j = 0; j < numdim_; j++)
  {
    const double C_inv_0_j = C_inv(0, j);
    const double C_inv_1_j = C_inv(0, j);
    const double C_inv_2_j = C_inv(0, j);

    for (int i = 0; i < numnod_; i++)
    {
      const int k = numdim_ * i + j;
      const double N_XYZ_Finv_0_i = N_XYZ_Finv(0, i);
      const double N_XYZ_Finv_1_i = N_XYZ_Finv(0, i);
      const double N_XYZ_Finv_2_i = N_XYZ_Finv(0, i);

      dfstressb_dv(0, k) = 2 * N_XYZ_Finv_0_i * C_inv_0_j;
      dfstressb_dv(1, k) = 2 * N_XYZ_Finv_1_i * C_inv_1_j;
      dfstressb_dv(2, k) = 2 * N_XYZ_Finv_2_i * C_inv_2_j;
      //**********************************
      dfstressb_dv(3, k) = N_XYZ_Finv_0_i * C_inv_1_j + N_XYZ_Finv_1_i * C_inv_0_j;
      dfstressb_dv(4, k) = N_XYZ_Finv_1_i * C_inv_2_j + N_XYZ_Finv_2_i * C_inv_1_j;
      dfstressb_dv(5, k) = N_XYZ_Finv_2_i * C_inv_0_j + N_XYZ_Finv_0_i * C_inv_2_j;
    }
  }

  // B^T . dfstress/dv^f
  CORE::LINALG::Matrix<numdof_, numdof_> dfstressb_dv_bop(true);
  dfstressb_dv_bop.MultiplyTN(bop, dfstressb_dv);

  for (int i = 0; i < numnod_; i++)
  {
    const int fi = noddof_ * i;

    for (int j = 0; j < numdim_; j++)
    {
      const double fstressb_i_j = fstressb(fi + j);

      for (int k = 0; k < numnod_; k++)
      {
        const int fk = noddof_ * k;
        const int fkp1 = (numdim_ + 1) * k;

        /*-------structure- fluid pressure coupling: "darcy-brinkman stress terms"
         B^T . ( \mu*J - d(phi)/(dp) * fstress ) * Dp
         */
        (*stiffmatrix)(fi + j, fkp1 + numdim_) +=
            detJ_w * fstressb_i_j * dphi_dp * visc * J * shapefct(k);
        for (int l = 0; l < noddof_; l++)
        {
          /*-------structure- fluid velocity coupling: "darcy-brinkman stress terms"
           B^T . ( \mu*J - phi * dfstress/dv^f ) * Dp
           */
          (*stiffmatrix)(fi + j, fkp1 + l) +=
              detJ_w * visc * J * porosity * dfstressb_dv_bop(fi + j, fk + l);
        }
      }
    }
  }
}

template <class so3_ele, CORE::FE::CellType distype>
void DRT::ELEMENTS::So3_Poro<so3_ele, distype>::ComputeDefGradient(
    CORE::LINALG::Matrix<numdim_, numdim_>& defgrd,
    const CORE::LINALG::Matrix<numdim_, numnod_>& N_XYZ,
    const CORE::LINALG::Matrix<numdim_, numnod_>& xcurr)
{
  if (so3_ele::kintype_ == INPAR::STR::kinem_nonlinearTotLag)  // total lagrange (nonlinear)
  {
    // (material) deformation gradient F = d xcurr / d xrefe = xcurr * N_XYZ^T
    defgrd.MultiplyNT(xcurr, N_XYZ);  //  (6.17)
  }
  else if (so3_ele::kintype_ == INPAR::STR::kinem_linear)  // linear kinematics
  {
    defgrd.Clear();
    for (int i = 0; i < numdim_; i++) defgrd(i, i) = 1.0;
  }
  else
    dserror("invalid kinematic type!");
}

template <class so3_ele, CORE::FE::CellType distype>
void DRT::ELEMENTS::So3_Poro<so3_ele, distype>::GetCauchyNDirAndDerivativesAtXi(
    const CORE::LINALG::Matrix<3, 1>& xi, const std::vector<double>& disp,
    const std::vector<double>& pres, const CORE::LINALG::Matrix<3, 1>& n,
    const CORE::LINALG::Matrix<3, 1>& dir, double& cauchy_n_dir,
    CORE::LINALG::SerialDenseMatrix* d_cauchyndir_dd,
    CORE::LINALG::SerialDenseMatrix* d_cauchyndir_dp, CORE::LINALG::Matrix<3, 1>* d_cauchyndir_dn,
    CORE::LINALG::Matrix<3, 1>* d_cauchyndir_ddir, CORE::LINALG::Matrix<3, 1>* d_cauchyndir_dxi)
{
  if (fluid_mat_->Type() != MAT::PAR::darcy)
    dserror("GetCauchyAtXi just implemented for pure Darcy flow!");

  if (distype != CORE::FE::CellType::hex8)
    dserror("GetCauchyAtXi for Poro just implemented for hex8!");

  so3_ele::GetCauchyNDirAndDerivativesAtXi(xi, disp, n, dir, cauchy_n_dir, d_cauchyndir_dd, nullptr,
      nullptr, nullptr, nullptr, d_cauchyndir_dn, d_cauchyndir_ddir, d_cauchyndir_dxi, nullptr,
      nullptr, nullptr, nullptr, nullptr);

  // Add pressure to sigma_nt
  const double dot = n(0, 0) * dir(0, 0) + n(1, 0) * dir(1, 0) + n(2, 0) * dir(2, 0);
  if (fabs(dot) > 1e-30)
  {
    CORE::LINALG::Matrix<NUMNOD_SOH8, 1> shapefcts;
    CORE::DRT::UTILS::shape_function<CORE::FE::CellType::hex8>(xi, shapefcts);

    for (unsigned nlid = 0; nlid < NUMNOD_SOH8; ++nlid)
      cauchy_n_dir -= pres[nlid] * shapefcts(nlid, 0) * dot;

    if (d_cauchyndir_dp || d_cauchyndir_dn || d_cauchyndir_ddir || d_cauchyndir_dxi)
    {
      CORE::LINALG::Matrix<NUMDIM_SOH8, NUMNOD_SOH8> deriv;
      CORE::DRT::UTILS::shape_function_deriv1<CORE::FE::CellType::hex8>(xi, deriv);

      d_cauchyndir_dp->reshape(NUMNOD_SOH8, 1);
      CORE::LINALG::Matrix<NUMNOD_SOH8, 1> dsntdp_m(d_cauchyndir_dp->values(), true);

      for (unsigned nlid = 0; nlid < NUMNOD_SOH8; ++nlid)
      {
        dsntdp_m(nlid, 0) = -dot * shapefcts(nlid, 0);
        for (unsigned dim = 0; dim < 3; ++dim)
        {
          (*d_cauchyndir_dn)(dim, 0) -= pres[nlid] * shapefcts(nlid, 0) * dir(dim, 0);
          (*d_cauchyndir_ddir)(dim, 0) -= pres[nlid] * shapefcts(nlid, 0) * n(dim, 0);
          (*d_cauchyndir_dxi)(dim, 0) -= pres[nlid] * deriv(dim, nlid) * dot;
        }
      }
    }
  }
}

template <class so3_ele, CORE::FE::CellType distype>
std::vector<double>
DRT::ELEMENTS::So3_Poro<so3_ele, distype>::ComputeAnisotropicPermeabilityCoeffsAtGP(
    const CORE::LINALG::Matrix<numnod_, 1>& shapefct) const
{
  std::vector<double> anisotropic_permeability_coeffs(numdim_, 0.0);

  for (int node = 0; node < numnod_; ++node)
  {
    const double shape_val = shapefct(node);
    for (int dim = 0; dim < numdim_; ++dim)
    {
      anisotropic_permeability_coeffs[dim] +=
          shape_val * anisotropic_permeability_nodal_coeffs_[dim][node];
    }
  }

  return anisotropic_permeability_coeffs;
}

BACI_NAMESPACE_CLOSE

#include "baci_so3_poro_fwd.hpp"
