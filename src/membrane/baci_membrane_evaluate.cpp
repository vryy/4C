/*----------------------------------------------------------------------*/
/*! \file
\brief

\level 3


\brief Nonlinear Membrane Finite Element evaluation

*----------------------------------------------------------------------*/
#include "baci_discretization_fem_general_utils_fem_shapefunctions.H"
#include "baci_global_data.H"
#include "baci_lib_discret.H"
#include "baci_lib_utils.H"
#include "baci_linalg_fixedsizematrix.H"
#include "baci_linalg_fixedsizematrix_tensor_transformation.H"
#include "baci_linalg_utils_densematrix_eigen.H"
#include "baci_mat_material.H"
#include "baci_mat_membrane_elasthyper.H"
#include "baci_mat_membrane_material_interfaces.H"
#include "baci_membrane.H"
#include "baci_membrane_service.H"
#include "baci_structure_new_elements_paramsinterface.H"
#include "baci_utils_function_of_time.H"

#include <Teuchos_RCP.hpp>

BACI_NAMESPACE_OPEN


/*----------------------------------------------------------------------*
 |  evaluate the element (public)                          fbraeu 06/16 |
 *----------------------------------------------------------------------*/
template <CORE::FE::CellType distype>
int DRT::ELEMENTS::Membrane<distype>::Evaluate(Teuchos::ParameterList& params,
    DRT::Discretization& discretization, std::vector<int>& lm,
    CORE::LINALG::SerialDenseMatrix& elemat1_epetra,
    CORE::LINALG::SerialDenseMatrix& elemat2_epetra,
    CORE::LINALG::SerialDenseVector& elevec1_epetra,
    CORE::LINALG::SerialDenseVector& elevec2_epetra,
    CORE::LINALG::SerialDenseVector& elevec3_epetra)
{
  // determine size of each element matrix
  CORE::LINALG::Matrix<numdof_, numdof_> elemat1(elemat1_epetra.values(), true);
  CORE::LINALG::Matrix<numdof_, numdof_> elemat2(elemat2_epetra.values(), true);
  CORE::LINALG::Matrix<numdof_, 1> elevec1(elevec1_epetra.values(), true);
  CORE::LINALG::Matrix<numdof_, 1> elevec2(elevec2_epetra.values(), true);
  CORE::LINALG::Matrix<numdof_, 1> elevec3(elevec3_epetra.values(), true);

  // set params interface pointer
  SetParamsInterfacePtr(params);

  // start with ActionType none
  ELEMENTS::ActionType act = ELEMENTS::none;

  if (IsParamsInterface())  // new structural time integration
  {
    act = ParamsInterface().GetActionType();
  }
  else  // old structural time integration
  {
    // get the action required
    std::string action = params.get<std::string>("action", "none");
    if (action == "none")
      dserror("No action supplied");
    else if (action == "calc_struct_nlnstiff")
      act = ELEMENTS::struct_calc_nlnstiff;
    else if (action == "calc_struct_nlnstiffmass")
      act = ELEMENTS::struct_calc_nlnstiffmass;
    else if (action == "calc_struct_update_istep")
      act = ELEMENTS::struct_calc_update_istep;
    else if (action == "calc_struct_reset_istep")
      act = ELEMENTS::struct_calc_reset_istep;
    else if (action == "calc_struct_stress")
      act = ELEMENTS::struct_calc_stress;
    else if (action == "calc_struct_thickness")
      act = ELEMENTS::struct_calc_thickness;
    else if (action == "calc_struct_energy")
      act = ELEMENTS::struct_calc_energy;
    else if (action == "postprocess_thickness")
      act = ELEMENTS::struct_postprocess_thickness;
    else
    {
      dserror("Unknown type of action for Membrane: %s", action.c_str());
    }
  }

  switch (act)
  {
    /*===============================================================================*
     | struct_calc_nlnstiff                                                          |
     *===============================================================================*/
    case ELEMENTS::struct_calc_nlnstiff:
    {
      // need current displacement
      Teuchos::RCP<const Epetra_Vector> disp = discretization.GetState("displacement");
      if (disp == Teuchos::null) dserror("Cannot get state vector 'displacement'");
      std::vector<double> mydisp(lm.size());
      DRT::UTILS::ExtractMyValues(*disp, mydisp, lm);
      CORE::LINALG::Matrix<numdof_, numdof_>* matptr = nullptr;
      if (elemat1.IsInitialized()) matptr = &elemat1;

      mem_nlnstiffmass(lm, mydisp, matptr, nullptr, &elevec1, nullptr, nullptr, params,
          INPAR::STR::stress_none, INPAR::STR::strain_none);
    }
    break;

    /*===============================================================================*
     | struct_calc_nlnstiffmass                                                      |
     *===============================================================================*/
    case ELEMENTS::struct_calc_nlnstiffmass:  // do mass, stiffness and internal forces
    {
      // need current displacement
      Teuchos::RCP<const Epetra_Vector> disp = discretization.GetState("displacement");
      if (disp == Teuchos::null) dserror("Cannot get state vector 'displacement'");
      std::vector<double> mydisp(lm.size());
      DRT::UTILS::ExtractMyValues(*disp, mydisp, lm);
      CORE::LINALG::Matrix<numdof_, numdof_>* matptr = nullptr;
      if (elemat1.IsInitialized()) matptr = &elemat1;

      mem_nlnstiffmass(lm, mydisp, matptr, &elemat2, &elevec1, nullptr, nullptr, params,
          INPAR::STR::stress_none, INPAR::STR::strain_none);
    }
    break;

    /*===============================================================================*
     | struct_calc_internalforce                                                     |
     *===============================================================================*/
    case ELEMENTS::struct_calc_internalforce:
    {
      // need current displacement
      Teuchos::RCP<const Epetra_Vector> disp = discretization.GetState("displacement");
      if (disp == Teuchos::null) dserror("Cannot get state vector 'displacement'");
      std::vector<double> mydisp(lm.size());
      DRT::UTILS::ExtractMyValues(*disp, mydisp, lm);

      mem_nlnstiffmass(lm, mydisp, nullptr, nullptr, &elevec1, nullptr, nullptr, params,
          INPAR::STR::stress_none, INPAR::STR::strain_none);
    }
    break;

    /*===============================================================================*
     | struct_calc_update_istep                                                      |
     *===============================================================================*/
    case ELEMENTS::struct_calc_update_istep:
    {
      // Update materials
      Teuchos::RCP<const Epetra_Vector> disp = discretization.GetState("displacement");
      if (disp == Teuchos::null) dserror("Cannot get state vectors 'displacement'");
      std::vector<double> mydisp(lm.size());
      DRT::UTILS::ExtractMyValues(*disp, mydisp, lm);
      Update_element(mydisp, params, Material());
    }
    break;

    /*===============================================================================*
     | struct_calc_reset_istep                                                       |
     *===============================================================================*/
    case ELEMENTS::struct_calc_reset_istep:
    {
      // Reset of history (if needed)
      SolidMaterial()->ResetStep();
    }
    break;

    /*===============================================================================*
     | struct_calc_stress                                                            |
     *===============================================================================*/
    case ELEMENTS::struct_calc_stress:
    {
      // nothing to do for ghost elements
      if (discretization.Comm().MyPID() == Owner())
      {
        // need current displacement
        Teuchos::RCP<const Epetra_Vector> disp = discretization.GetState("displacement");
        if (disp == Teuchos::null) dserror("Cannot get state vectors 'displacement'");
        std::vector<double> mydisp(lm.size());
        DRT::UTILS::ExtractMyValues(*disp, mydisp, lm);

        Teuchos::RCP<std::vector<char>> stressdata = Teuchos::null;
        Teuchos::RCP<std::vector<char>> straindata = Teuchos::null;

        INPAR::STR::StressType iostress = INPAR::STR::stress_none;
        INPAR::STR::StrainType iostrain = INPAR::STR::strain_none;

        if (IsParamsInterface())  // new structural time integration
        {
          stressdata = StrParamsInterface().StressDataPtr();
          straindata = StrParamsInterface().StrainDataPtr();

          iostress = StrParamsInterface().GetStressOutputType();
          iostrain = StrParamsInterface().GetStrainOutputType();
        }
        else  // old structural time integration
        {
          stressdata = params.get<Teuchos::RCP<std::vector<char>>>("stress", Teuchos::null);
          straindata = params.get<Teuchos::RCP<std::vector<char>>>("strain", Teuchos::null);

          iostress =
              INPUT::get<INPAR::STR::StressType>(params, "iostress", INPAR::STR::stress_none);
          iostrain =
              INPUT::get<INPAR::STR::StrainType>(params, "iostrain", INPAR::STR::strain_none);
        }

        if (stressdata == Teuchos::null) dserror("Cannot get 'stress' data");
        if (straindata == Teuchos::null) dserror("Cannot get 'strain' data");

        CORE::LINALG::Matrix<numgpt_post_, 6> stress;
        CORE::LINALG::Matrix<numgpt_post_, 6> strain;

        // determine strains and/or stresses
        mem_nlnstiffmass(
            lm, mydisp, nullptr, nullptr, nullptr, &stress, &strain, params, iostress, iostrain);

        // add data to pack
        {
          CORE::COMM::PackBuffer data;
          AddtoPack(data, stress);
          data.StartPacking();
          AddtoPack(data, stress);
          std::copy(data().begin(), data().end(), std::back_inserter(*stressdata));
        }

        {
          CORE::COMM::PackBuffer data;
          AddtoPack(data, strain);
          data.StartPacking();
          AddtoPack(data, strain);
          std::copy(data().begin(), data().end(), std::back_inserter(*straindata));
        }
      }
    }
    break;

    /*===============================================================================*
     | struct_calc_thickness                                                         |
     *===============================================================================*/
    case ELEMENTS::struct_calc_thickness:
    {
      // nothing to do for ghost elements
      if (discretization.Comm().MyPID() == Owner())
      {
        Teuchos::RCP<std::vector<char>> thickdata = Teuchos::null;

        if (IsParamsInterface())  // new structural time integration
          thickdata = StrParamsInterface().OptQuantityDataPtr();
        else  // old structural time integration
          thickdata = params.get<Teuchos::RCP<std::vector<char>>>("optquantity", Teuchos::null);

        if (thickdata == Teuchos::null) dserror("Cannot get 'thickness' data");

        CORE::LINALG::Matrix<numgpt_post_, 1> thickness;
        for (int i = 0; i < numgpt_post_; ++i) thickness(i) = cur_thickness_[i];

        // add data to pack
        {
          CORE::COMM::PackBuffer data;
          AddtoPack(data, thickness);
          data.StartPacking();
          AddtoPack(data, thickness);
          std::copy(data().begin(), data().end(), std::back_inserter(*thickdata));
        }
      }
    }
    break;

    /*===============================================================================*
     | struct_calc_energy                                                            |
     *===============================================================================*/
    case ELEMENTS::struct_calc_energy:
    {
      // initialization of internal energy
      double intenergy = 0.0;

      // need current displacement
      Teuchos::RCP<const Epetra_Vector> disp = discretization.GetState("displacement");
      if (disp == Teuchos::null) dserror("Cannot get state vectors 'displacement'");
      std::vector<double> mydisp(lm.size());
      DRT::UTILS::ExtractMyValues(*disp, mydisp, lm);

      // get reference configuration and determine current configuration
      CORE::LINALG::Matrix<numnod_, noddof_> xrefe(true);
      CORE::LINALG::Matrix<numnod_, noddof_> xcurr(true);

      mem_configuration(mydisp, xrefe, xcurr);

      /*===============================================================================*
       | loop over the gauss points                                                    |
       *===============================================================================*/

      // allocate matrix for shape function derivatives at gp
      CORE::LINALG::Matrix<numdim_, numnod_> derivs(true);

      for (int gp = 0; gp < intpoints_.nquad; ++gp)
      {
        // get gauss points from integration rule
        double xi_gp = intpoints_.qxg[gp][0];
        double eta_gp = intpoints_.qxg[gp][1];

        // get gauss weight at current gp
        double gpweight = intpoints_.qwgt[gp];

        // get shape function derivatives in the plane of the element
        CORE::FE::shape_function_2D_deriv1(derivs, xi_gp, eta_gp, Shape());

        /*===============================================================================*
         | orthonormal base (t1,t2,tn) in the undeformed configuration at current GP     |
         *===============================================================================*/

        CORE::LINALG::Matrix<numdim_, numnod_> derivs_ortho(true);
        double G1G2_cn;
        CORE::LINALG::Matrix<noddof_, 1> dXds1(true);
        CORE::LINALG::Matrix<noddof_, 1> dXds2(true);
        CORE::LINALG::Matrix<noddof_, 1> dxds1(true);
        CORE::LINALG::Matrix<noddof_, 1> dxds2(true);
        CORE::LINALG::Matrix<noddof_, noddof_> Q_localToGlobal(true);

        mem_orthonormalbase(xrefe, xcurr, derivs, derivs_ortho, G1G2_cn, dXds1, dXds2, dxds1, dxds2,
            Q_localToGlobal);

        /*===============================================================================*
         | surface deformation gradient                                                  |
         *===============================================================================*/

        // surface deformation gradient in 3 dimensions in global coordinates
        CORE::LINALG::Matrix<noddof_, noddof_> defgrd_glob(true);

        // surface deformation gradient in 3 dimensions in local coordinates
        CORE::LINALG::Matrix<noddof_, noddof_> defgrd_loc(true);

        // principle stretch in thickness direction
        double lambda3 = 1.0;

        // standard evaluation (incompressible, plane stress)
        if (Material()->MaterialType() == INPAR::MAT::m_membrane_elasthyper)
        {
          // incompressibility condition to get principle stretch in thickness direction
          lambda3 = std::sqrt(
              1.0 / (dxds1.Dot(dxds1) * dxds2.Dot(dxds2) - std::pow(dxds1.Dot(dxds2), 2.0)));
        }
        else
          dserror(
              "Type of material not implemented for evaluation of strain energy for membranes!");

        // surface deformation gradient in 3 dimensions in global coordinates
        mem_defgrd_global(dXds1, dXds2, dxds1, dxds2, lambda3, defgrd_glob);

        // surface deformation gradient in 3 dimensions in local coordinates
        CORE::LINALG::TENSOR::InverseTensorRotation<3>(Q_localToGlobal, defgrd_glob, defgrd_loc);

        /*===============================================================================*
         | right cauchygreen tensor in local coordinates                                 |
         *===============================================================================*/

        // calculate three dimensional right cauchy-green strain tensor in orthonormal base
        CORE::LINALG::Matrix<noddof_, noddof_> cauchygreen_loc(true);
        cauchygreen_loc.MultiplyTN(1.0, defgrd_loc, defgrd_loc, 0.0);

        /*===============================================================================*
         | call material law for evaluation of strain energy                            |
         *===============================================================================*/

        double psi = 0.0;

        // standard evaluation (incompressible, plane stress)
        if (Material()->MaterialType() == INPAR::MAT::m_membrane_elasthyper)
        {
          Teuchos::rcp_dynamic_cast<MAT::Membrane_ElastHyper>(DRT::Element::Material(), true)
              ->StrainEnergy(cauchygreen_loc, psi, gp, Id());
        }
        else
          dserror(
              "Type of material not implemented for evaluation of strain energy for membranes!");

        // add gauss point contribution to internal energy
        double fac = gpweight * thickness_ * G1G2_cn;
        intenergy += fac * psi;
      }

      if (IsParamsInterface())  // new structural time integration
      {
        // only add contributions from row elements to avoid counting them on more than one proc
        if (discretization.Comm().MyPID() == Owner())
          StrParamsInterface().AddContributionToEnergyType(intenergy, STR::internal_energy);
      }
      else  // old structural time integration
      {
        // check length of elevec1
        if (elevec1_epetra.length() < 1) dserror("The given result vector is too short.");

        elevec1_epetra(0) = intenergy;
      }
    }
    break;

    /*===============================================================================*
     | struct_postprocess_thickness                                                  |
     *===============================================================================*/
    case ELEMENTS::struct_postprocess_thickness:
    {
      const Teuchos::RCP<std::map<int, Teuchos::RCP<CORE::LINALG::SerialDenseMatrix>>> gpthickmap =
          params.get<Teuchos::RCP<std::map<int, Teuchos::RCP<CORE::LINALG::SerialDenseMatrix>>>>(
              "gpthickmap", Teuchos::null);
      if (gpthickmap == Teuchos::null) dserror("no gp thickness map available for postprocessing");

      std::string optquantitytype = params.get<std::string>("optquantitytype", "ndxyz");

      int gid = Id();
      CORE::LINALG::Matrix<numgpt_post_, 1> gpthick(((*gpthickmap)[gid])->values(), true);

      Teuchos::RCP<Epetra_MultiVector> postthick =
          params.get<Teuchos::RCP<Epetra_MultiVector>>("postthick", Teuchos::null);
      if (postthick == Teuchos::null) dserror("No element thickness vector available");

      if (optquantitytype == "ndxyz")
      {
        // extrapolation matrix: static because equal for all elements of the same discretization
        // type
        static CORE::LINALG::Matrix<numnod_, numgpt_post_> extrapol(mem_extrapolmat());

        // extrapolate the nodal thickness for current element
        CORE::LINALG::Matrix<numnod_, 1> nodalthickness;
        nodalthickness.Multiply(1.0, extrapol, gpthick, 0.0);

        // "assembly" of extrapolated nodal thickness
        for (int i = 0; i < numnod_; ++i)
        {
          int gid = NodeIds()[i];
          if (postthick->Map().MyGID(NodeIds()[i]))  // rownode
          {
            int lid = postthick->Map().LID(gid);
            int myadjele = Nodes()[i]->NumElement();
            (*((*postthick)(0)))[lid] += nodalthickness(i) / myadjele;
          }
        }
      }
      else
        dserror("unknown type of thickness output on element level");
    }
    break;

    /*===============================================================================*
     | struct_calc_recover                                                           |
     *===============================================================================*/
    case ELEMENTS::struct_calc_recover:
    {
      // do nothing here
    }
    break;

    case ELEMENTS::struct_calc_predict:
    {
      // do nothing here
      break;
    }

    /*===============================================================================*
     | default                                                                       |
     *===============================================================================*/
    default:
      dserror("Unknown type of action for Membrane: %s", ActionType2String(act).c_str());
      break;
  }

  return 0;
}


/*-----------------------------------------------------------------------*
 |  Integrate a Surface Neumann boundary condition (public) fbraeu 06/16 |
 *-----------------------------------------------------------------------*/
template <CORE::FE::CellType distype>
int DRT::ELEMENTS::Membrane<distype>::EvaluateNeumann(Teuchos::ParameterList& params,
    DRT::Discretization& discretization, DRT::Condition& condition, std::vector<int>& lm,
    CORE::LINALG::SerialDenseVector& elevec1_epetra,
    CORE::LINALG::SerialDenseMatrix* elemat1_epetra)
{
  // set params interface pointer
  SetParamsInterfacePtr(params);

  // get values and switches from the condition
  const std::vector<int>* onoff = condition.Get<std::vector<int>>("onoff");
  const std::vector<double>* val = condition.Get<std::vector<double>>("val");

  // find out whether we will use a time curve
  double time = -1.0;

  if (IsParamsInterface())  // new structural time integration
    time = ParamsInterface().GetTotalTime();
  else  // old structural time integration
    time = params.get("total time", -1.0);

  // ensure that at least as many curves/functs as dofs are available
  if (int(onoff->size()) < noddof_)
    dserror("Fewer functions or curves defined than the element has dofs.");

  // check membrane pressure input
  for (int checkdof = 1; checkdof < int(onoff->size()); ++checkdof)
    if ((*onoff)[checkdof] != 0) dserror("membrane pressure on 1st dof only!");

  // find out whether we will use time curves and get the factors
  const std::vector<int>* tmp_funct = condition.Get<std::vector<int>>("funct");
  std::vector<double> functfacs(noddof_, 1.0);
  for (int i = 0; i < noddof_; ++i)
  {
    const int functnum = (tmp_funct) ? (*tmp_funct)[i] : -1;
    if (functnum > 0)
      functfacs[i] = DRT::Problem::Instance()
                         ->FunctionById<CORE::UTILS::FunctionOfTime>(functnum - 1)
                         .Evaluate(time);
  }

  // determine current pressure
  double pressure;
  if ((*onoff)[0])
    pressure = (*val)[0] * functfacs[0];
  else
    pressure = 0.0;

  // need displacement new
  Teuchos::RCP<const Epetra_Vector> disp = discretization.GetState("displacement new");
  if (disp == Teuchos::null) dserror("Cannot get state vector 'displacement new'");
  std::vector<double> mydisp(lm.size());
  DRT::UTILS::ExtractMyValues(*disp, mydisp, lm);

  // get reference configuration and determine current configuration
  CORE::LINALG::Matrix<numnod_, noddof_> xrefe(true);
  CORE::LINALG::Matrix<numnod_, noddof_> xcurr(true);

  mem_configuration(mydisp, xrefe, xcurr);

  /*===============================================================================*
   | loop over the gauss points                                                    |
   *===============================================================================*/

  // allocate vector for shape functions and matrix for derivatives at gp
  CORE::LINALG::Matrix<numnod_, 1> shapefcts(true);
  CORE::LINALG::Matrix<numdim_, numnod_> derivs(true);

  for (int gp = 0; gp < intpoints_.nquad; ++gp)
  {
    // get gauss points from integration rule
    double xi_gp = intpoints_.qxg[gp][0];
    double eta_gp = intpoints_.qxg[gp][1];

    // get gauss weight at current gp
    double gpweight = intpoints_.qwgt[gp];

    // get shape functions and derivatives in the plane of the element
    CORE::FE::shape_function_2D(shapefcts, xi_gp, eta_gp, Shape());
    CORE::FE::shape_function_2D_deriv1(derivs, xi_gp, eta_gp, Shape());

    /*===============================================================================*
     | orthonormal base (t1,t2,tn) in the undeformed configuration at current GP     |
     *===============================================================================*/

    CORE::LINALG::Matrix<numdim_, numnod_> derivs_ortho(true);
    double G1G2_cn;
    CORE::LINALG::Matrix<noddof_, 1> dXds1(true);
    CORE::LINALG::Matrix<noddof_, 1> dXds2(true);
    CORE::LINALG::Matrix<noddof_, 1> dxds1(true);
    CORE::LINALG::Matrix<noddof_, 1> dxds2(true);
    CORE::LINALG::Matrix<noddof_, noddof_> Q_localToGlobal(true);

    mem_orthonormalbase(
        xrefe, xcurr, derivs, derivs_ortho, G1G2_cn, dXds1, dXds2, dxds1, dxds2, Q_localToGlobal);

    // determine cross product x,1 x x,2
    CORE::LINALG::Matrix<noddof_, 1> xcurr_cross(true);
    xcurr_cross(0) = dxds1(1) * dxds2(2) - dxds1(2) * dxds2(1);
    xcurr_cross(1) = dxds1(2) * dxds2(0) - dxds1(0) * dxds2(2);
    xcurr_cross(2) = dxds1(0) * dxds2(1) - dxds1(1) * dxds2(0);

    // determine cross product X,1 x X,2
    CORE::LINALG::Matrix<noddof_, 1> xrefe_cross(true);
    xrefe_cross(0) = dXds1(1) * dXds2(2) - dXds1(2) * dXds2(1);
    xrefe_cross(1) = dXds1(2) * dXds2(0) - dXds1(0) * dXds2(2);
    xrefe_cross(2) = dXds1(0) * dXds2(1) - dXds1(1) * dXds2(0);

    // euclidian norm of xref_cross
    double xrefe_cn = xrefe_cross.Norm2();

    // integration factor
    double fac = (pressure * G1G2_cn * gpweight) / xrefe_cn;

    // loop over all 4 nodes
    for (int i = 0; i < numnod_; ++i)
    {
      // assemble external force vector
      elevec1_epetra[noddof_ * i + 0] += fac * xcurr_cross(0) * (shapefcts)(i);
      elevec1_epetra[noddof_ * i + 1] += fac * xcurr_cross(1) * (shapefcts)(i);
      elevec1_epetra[noddof_ * i + 2] += fac * xcurr_cross(2) * (shapefcts)(i);

      // evaluate external stiffness matrix if needed
      if (elemat1_epetra != nullptr)
      {
        // determine P matrix for all 4 nodes, Gruttmann92 equation (41) and directly fill up
        // elemat1_epetra
        for (int j = 0; j < numnod_; ++j)
        {
          double p1_ij =
              (dxds1(0) * derivs_ortho(1, i) - dxds2(0) * derivs_ortho(0, i)) * (shapefcts)(j);
          double p2_ij =
              (dxds1(1) * derivs_ortho(1, i) - dxds2(1) * derivs_ortho(0, i)) * (shapefcts)(j);
          double p3_ij =
              (dxds1(2) * derivs_ortho(1, i) - dxds2(2) * derivs_ortho(0, i)) * (shapefcts)(j);

          // entries of P matrix are in round brackets
          (*elemat1_epetra)(noddof_ * i + 0, noddof_ * j + 1) += fac * -p3_ij;
          (*elemat1_epetra)(noddof_ * i + 0, noddof_ * j + 2) += fac * +p2_ij;
          (*elemat1_epetra)(noddof_ * i + 1, noddof_ * j + 0) += fac * +p3_ij;
          (*elemat1_epetra)(noddof_ * i + 1, noddof_ * j + 2) += fac * -p1_ij;
          (*elemat1_epetra)(noddof_ * i + 2, noddof_ * j + 0) += fac * -p2_ij;
          (*elemat1_epetra)(noddof_ * i + 2, noddof_ * j + 1) += fac * +p1_ij;
        }
      }
    }
  }
  return 0;
}

/*----------------------------------------------------------------------*
 |  evaluate the element (private)                         fbraeu 06/16 |
 *----------------------------------------------------------------------*/
template <CORE::FE::CellType distype>
void DRT::ELEMENTS::Membrane<distype>::mem_nlnstiffmass(std::vector<int>& lm,  // location matrix
    std::vector<double>& disp,                            // current displacements
    CORE::LINALG::Matrix<numdof_, numdof_>* stiffmatrix,  // element stiffness matrix
    CORE::LINALG::Matrix<numdof_, numdof_>* massmatrix,   // element mass matrix
    CORE::LINALG::Matrix<numdof_, 1>* force,              // element internal force vector
    CORE::LINALG::Matrix<numgpt_post_, 6>* elestress,     // stresses at GP
    CORE::LINALG::Matrix<numgpt_post_, 6>* elestrain,     // strains at GP
    Teuchos::ParameterList& params,                       // algorithmic parameters e.g. time
    const INPAR::STR::StressType iostress,                // stress output option
    const INPAR::STR::StrainType iostrain)                // strain output option
{
  // get reference configuration and determine current configuration
  CORE::LINALG::Matrix<numnod_, noddof_> xrefe(true);
  CORE::LINALG::Matrix<numnod_, noddof_> xcurr(true);

  auto material_local_coordinates =
      Teuchos::rcp_dynamic_cast<MAT::MembraneMaterialLocalCoordinates>(DRT::Element::Material());
  auto material_global_coordinates =
      Teuchos::rcp_dynamic_cast<MAT::MembraneMaterialGlobalCoordinates>(DRT::Element::Material());
  auto material_inelastic_thickness =
      Teuchos::rcp_dynamic_cast<MAT::MembraneMaterialInelasticThickness>(DRT::Element::Material());

  mem_configuration(disp, xrefe, xcurr);

  /*===============================================================================*
   | loop over the gauss points                                                    |
   *===============================================================================*/

  // allocate vector for shape functions and matrix for derivatives at gp
  CORE::LINALG::Matrix<numnod_, 1> shapefcts(true);
  CORE::LINALG::Matrix<numdim_, numnod_> derivs(true);

  for (int gp = 0; gp < intpoints_.nquad; ++gp)
  {
    // get gauss points from integration rule
    double xi_gp = intpoints_.qxg[gp][0];
    double eta_gp = intpoints_.qxg[gp][1];

    // get gauss weight at current gp
    double gpweight = intpoints_.qwgt[gp];

    // get shape functions and derivatives in the plane of the element
    CORE::FE::shape_function_2D(shapefcts, xi_gp, eta_gp, Shape());
    CORE::FE::shape_function_2D_deriv1(derivs, xi_gp, eta_gp, Shape());

    /*===============================================================================*
     | orthonormal base (t1,t2,tn) in the undeformed configuration at current GP     |
     *===============================================================================*/

    CORE::LINALG::Matrix<numdim_, numnod_> derivs_ortho(true);
    double G1G2_cn;
    CORE::LINALG::Matrix<noddof_, 1> dXds1(true);
    CORE::LINALG::Matrix<noddof_, 1> dXds2(true);
    CORE::LINALG::Matrix<noddof_, 1> dxds1(true);
    CORE::LINALG::Matrix<noddof_, 1> dxds2(true);
    CORE::LINALG::Matrix<noddof_, noddof_> Q_localToGlobal(true);

    mem_orthonormalbase(
        xrefe, xcurr, derivs, derivs_ortho, G1G2_cn, dXds1, dXds2, dxds1, dxds2, Q_localToGlobal);

    /*===============================================================================*
     | surface deformation gradient                                                  |
     *===============================================================================*/

    // surface deformation gradient in 3 dimensions in global coordinates
    CORE::LINALG::Matrix<noddof_, noddof_> defgrd_glob(true);

    // surface deformation gradient in 3 dimensions in local coordinates
    CORE::LINALG::Matrix<noddof_, noddof_> defgrd_loc(true);

    // principle stretch in thickness direction
    double lambda3 = 1.0;

    if (material_inelastic_thickness != Teuchos::null)
    {
      // incompressibility is just valid for the elastic quantities, therefore
      // use thickness from previous iteration step to get principle stretch in thickness
      // direction.
      // Stretch in thickness direction is evaluated later by the material
      lambda3 = cur_thickness_[gp] / thickness_;
    }
    else
    {
      // standard evaluation (incompressible, plane stress)
      // incompressibility condition to get principle stretch in thickness direction
      lambda3 =
          std::sqrt(1.0 / (dxds1.Dot(dxds1) * dxds2.Dot(dxds2) - std::pow(dxds1.Dot(dxds2), 2.0)));
    }

    // surface deformation gradient in 3 dimensions in global coordinates
    mem_defgrd_global(dXds1, dXds2, dxds1, dxds2, lambda3, defgrd_glob);

    // surface deformation gradient in 3 dimensions in local coordinates
    CORE::LINALG::TENSOR::InverseTensorRotation<3>(Q_localToGlobal, defgrd_glob, defgrd_loc);

    /*===============================================================================*
     | right cauchygreen tensor in local coordinates                                 |
     *===============================================================================*/

    // calculate three dimensional right cauchy-green strain tensor in orthonormal base
    CORE::LINALG::Matrix<noddof_, noddof_> cauchygreen_loc(true);
    cauchygreen_loc.MultiplyTN(1.0, defgrd_loc, defgrd_loc, 0.0);

    /*===============================================================================*
     | call material law                                                             |
     *===============================================================================*/

    // 2nd piola kirchhoff stress vector under plane stress assumption
    CORE::LINALG::Matrix<3, 1> pk2red_loc(true);

    // material tangent matrix for plane stress
    CORE::LINALG::Matrix<3, 3> cmatred_loc(true);

    // The growth remodel elast hyper material needs some special quantities for its evaluation
    if (Material()->MaterialType() == INPAR::MAT::m_growthremodel_elasthyper)
    {
      // Gauss-point coordinates in reference configuration
      CORE::LINALG::Matrix<1, noddof_> gprefecoord(true);
      gprefecoord.MultiplyTN(shapefcts, xrefe);
      params.set("gprefecoord", gprefecoord);

      // center of element in reference configuration
      CORE::LINALG::Matrix<numnod_, 1> funct_center;
      CORE::FE::shape_function_2D(funct_center, 0.0, 0.0, distype);
      CORE::LINALG::Matrix<1, noddof_> midpoint;
      midpoint.MultiplyTN(funct_center, xrefe);
      params.set("elecenter", midpoint);
    }

    if (material_inelastic_thickness != Teuchos::null)
    {
      // Let material decide the total stretch in thickness direction
      lambda3 = material_inelastic_thickness->EvaluateMembraneThicknessStretch(
          defgrd_glob, params, gp, Id());

      // update surface deformation gradient in 3 dimensions in global coordinates
      mem_defgrd_global(dXds1, dXds2, dxds1, dxds2, lambda3, defgrd_glob);

      // update surface deformation gradient in 3 dimensions in local coordinates
      CORE::LINALG::TENSOR::InverseTensorRotation<3>(Q_localToGlobal, defgrd_glob, defgrd_loc);

      // update three dimensional right cauchy-green strain tensor in orthonormal base
      cauchygreen_loc.MultiplyTN(1.0, defgrd_loc, defgrd_loc, 0.0);
    }

    // standard evaluation (incompressible, plane stress)
    if (material_local_coordinates != Teuchos::null)
    {
      material_local_coordinates->EvaluateMembrane(
          defgrd_loc, cauchygreen_loc, params, Q_localToGlobal, pk2red_loc, cmatred_loc, gp, Id());
    }
    else if (material_global_coordinates != Teuchos::null)
    {
      CORE::LINALG::Matrix<3, 3> pk2M_glob(true);
      CORE::LINALG::Matrix<6, 6> cmat_glob(true);

      // Evaluate material with quantities in the global coordinate system
      material_global_coordinates->EvaluateMembrane(
          defgrd_glob, params, pk2M_glob, cmat_glob, gp, Id());

      // Transform stress and elasticity into the local membrane coordinate system
      CORE::LINALG::Matrix<3, 3> pk2M_loc(true);
      CORE::LINALG::TENSOR::InverseTensorRotation<3>(Q_localToGlobal, pk2M_glob, pk2M_loc);
      MEMBRANE::LocalPlaneStressToStressLikeVoigt(pk2M_loc, pk2red_loc);

      CORE::LINALG::Matrix<6, 6> cmat_loc(true);
      CORE::LINALG::TENSOR::InverseFourthTensorRotation(Q_localToGlobal, cmat_glob, cmat_loc);
      MEMBRANE::LocalFourthTensorPlaneStressToStressLikeVoigt(cmat_loc, cmatred_loc);
    }
    else
    {
      dserror("The material does not support the evaluation of membranes");
    }

    /*===============================================================================*
     | update current thickness at gp                                                |
     *===============================================================================*/
    cur_thickness_[gp] = lambda3 * thickness_;

    /*===============================================================================*
     | calculate force, stiffness matrix and mass matrix                             |
     *===============================================================================*/
    // evaluate just force vector (stiffness matrix not needed)
    if (stiffmatrix == nullptr && force != nullptr)
    {
      // determine B matrix for all 4 nodes, Gruttmann1992 equation (36)
      CORE::LINALG::Matrix<noddof_, numdof_> B_matrix(true);

      for (int i = 0; i < numnod_; ++i)
      {
        B_matrix(0, noddof_ * i + 0) = derivs_ortho(0, i) * dxds1(0);
        B_matrix(1, noddof_ * i + 0) = derivs_ortho(1, i) * dxds2(0);
        B_matrix(2, noddof_ * i + 0) =
            derivs_ortho(0, i) * dxds2(0) + derivs_ortho(1, i) * dxds1(0);

        B_matrix(0, noddof_ * i + 1) = derivs_ortho(0, i) * dxds1(1);
        B_matrix(1, noddof_ * i + 1) = derivs_ortho(1, i) * dxds2(1);
        B_matrix(2, noddof_ * i + 1) =
            derivs_ortho(0, i) * dxds2(1) + derivs_ortho(1, i) * dxds1(1);

        B_matrix(0, noddof_ * i + 2) = derivs_ortho(0, i) * dxds1(2);
        B_matrix(1, noddof_ * i + 2) = derivs_ortho(1, i) * dxds2(2);
        B_matrix(2, noddof_ * i + 2) =
            derivs_ortho(0, i) * dxds2(2) + derivs_ortho(1, i) * dxds1(2);
      }

      double fac = gpweight * thickness_ * G1G2_cn;

      // determine force and stiffness matrix, Gruttmann1992 equation (37) and (39)
      force->MultiplyTN(fac, B_matrix, pk2red_loc, 1.0);
    }

    // evaluate stiffness matrix and force vector if needed
    if (stiffmatrix != nullptr && force != nullptr)
    {
      // determine B matrix and G matrix for all 4 nodes, Gruttmann1992 equation (36) and (40)
      CORE::LINALG::Matrix<noddof_, numdof_> B_matrix(true);
      CORE::LINALG::Matrix<numdof_, numdof_> G_matrix(true);
      double g_ij;

      for (int i = 0; i < numnod_; ++i)
      {
        B_matrix(0, noddof_ * i + 0) = derivs_ortho(0, i) * dxds1(0);
        B_matrix(1, noddof_ * i + 0) = derivs_ortho(1, i) * dxds2(0);
        B_matrix(2, noddof_ * i + 0) =
            derivs_ortho(0, i) * dxds2(0) + derivs_ortho(1, i) * dxds1(0);

        B_matrix(0, noddof_ * i + 1) = derivs_ortho(0, i) * dxds1(1);
        B_matrix(1, noddof_ * i + 1) = derivs_ortho(1, i) * dxds2(1);
        B_matrix(2, noddof_ * i + 1) =
            derivs_ortho(0, i) * dxds2(1) + derivs_ortho(1, i) * dxds1(1);

        B_matrix(0, noddof_ * i + 2) = derivs_ortho(0, i) * dxds1(2);
        B_matrix(1, noddof_ * i + 2) = derivs_ortho(1, i) * dxds2(2);
        B_matrix(2, noddof_ * i + 2) =
            derivs_ortho(0, i) * dxds2(2) + derivs_ortho(1, i) * dxds1(2);

        for (int j = 0; j < numnod_; ++j)
        {
          g_ij = pk2red_loc(0) * derivs_ortho(0, i) * derivs_ortho(0, j) +
                 pk2red_loc(1) * derivs_ortho(1, i) * derivs_ortho(1, j) +
                 pk2red_loc(2) * (derivs_ortho(0, i) * derivs_ortho(1, j) +
                                     derivs_ortho(1, i) * derivs_ortho(0, j));

          G_matrix(noddof_ * i + 0, noddof_ * j + 0) = g_ij;
          G_matrix(noddof_ * i + 1, noddof_ * j + 1) = g_ij;
          G_matrix(noddof_ * i + 2, noddof_ * j + 2) = g_ij;
        }
      }

      double fac = gpweight * thickness_ * G1G2_cn;

      // determine force and stiffness matrix, Gruttmann1992 equation (37) and (39)
      force->MultiplyTN(fac, B_matrix, pk2red_loc, 1.0);

      CORE::LINALG::Matrix<numdof_, noddof_> temp(true);
      temp.MultiplyTN(1.0, B_matrix, cmatred_loc, 0.0);
      CORE::LINALG::Matrix<numdof_, numdof_> temp2(true);
      temp2.Multiply(1.0, temp, B_matrix, 0.0);
      temp2.Update(1.0, G_matrix, 1.0);

      stiffmatrix->Update(fac, temp2, 1.0);
    }

    // evaluate massmatrix if needed, just valid for a constant density
    if (massmatrix != nullptr)
    {
      // get density
      double density = SolidMaterial()->Density();

      // integrate consistent mass matrix
      const double factor = gpweight * thickness_ * G1G2_cn * density;
      double ifactor = 0.0;
      double massfactor = 0.0;

      for (int i = 0; i < numnod_; ++i)
      {
        ifactor = shapefcts(i) * factor;

        for (int j = 0; j < numnod_; ++j)
        {
          massfactor = shapefcts(j) * ifactor;  // intermediate factor

          (*massmatrix)(noddof_ * i + 0, noddof_ * j + 0) += massfactor;
          (*massmatrix)(noddof_ * i + 1, noddof_ * j + 1) += massfactor;
          (*massmatrix)(noddof_ * i + 2, noddof_ * j + 2) += massfactor;
        }
      }

      // check for non constant mass matrix
      if (SolidMaterial()->VaryingDensity())
      {
        dserror("Varying Density not supported for Membrane");
      }
    }

    /*===============================================================================*
     | return gp strains (only in case of stress/strain output)                      |
     *===============================================================================*/
    switch (iostrain)
    {
      // Green-Lagrange strains
      case INPAR::STR::strain_gl:
      {
        if (elestrain == nullptr) dserror("strain data not available");

        // transform local cauchygreen to global coordinates
        CORE::LINALG::Matrix<noddof_, noddof_> cauchygreen_glob(true);
        CORE::LINALG::TENSOR::TensorRotation<3>(Q_localToGlobal, cauchygreen_loc, cauchygreen_glob);

        // green-lagrange strain tensor in global coordinates
        CORE::LINALG::Matrix<noddof_, noddof_> glstrain_glob(true);
        glstrain_glob(0, 0) = 0.5 * (cauchygreen_glob(0, 0) - 1.0);
        glstrain_glob(1, 1) = 0.5 * (cauchygreen_glob(1, 1) - 1.0);
        glstrain_glob(2, 2) = 0.5 * (cauchygreen_glob(2, 2) - 1.0);
        glstrain_glob(0, 1) = 0.5 * cauchygreen_glob(0, 1);
        glstrain_glob(0, 2) = 0.5 * cauchygreen_glob(0, 2);
        glstrain_glob(1, 2) = 0.5 * cauchygreen_glob(1, 2);
        glstrain_glob(1, 0) = glstrain_glob(0, 1);
        glstrain_glob(2, 0) = glstrain_glob(0, 2);
        glstrain_glob(2, 1) = glstrain_glob(1, 2);

        (*elestrain)(gp, 0) = glstrain_glob(0, 0);
        (*elestrain)(gp, 1) = glstrain_glob(1, 1);
        (*elestrain)(gp, 2) = glstrain_glob(2, 2);
        (*elestrain)(gp, 3) = glstrain_glob(0, 1);
        (*elestrain)(gp, 4) = glstrain_glob(1, 2);
        (*elestrain)(gp, 5) = glstrain_glob(0, 2);
      }
      break;
      // Euler-Almansi strains
      case INPAR::STR::strain_ea:
      {
        if (elestrain == nullptr) dserror("strain data not available");

        // transform local cauchygreen to global coordinates
        CORE::LINALG::Matrix<noddof_, noddof_> cauchygreen_glob(true);
        CORE::LINALG::TENSOR::TensorRotation<3>(Q_localToGlobal, cauchygreen_loc, cauchygreen_glob);

        // green-lagrange strain tensor in global coordinates
        CORE::LINALG::Matrix<noddof_, noddof_> glstrain_glob(true);
        glstrain_glob(0, 0) = 0.5 * (cauchygreen_glob(0, 0) - 1);
        glstrain_glob(1, 1) = 0.5 * (cauchygreen_glob(1, 1) - 1);
        glstrain_glob(2, 2) = 0.5 * (cauchygreen_glob(2, 2) - 1);
        glstrain_glob(0, 1) = 0.5 * cauchygreen_glob(0, 1);
        glstrain_glob(0, 2) = 0.5 * cauchygreen_glob(0, 2);
        glstrain_glob(1, 2) = 0.5 * cauchygreen_glob(1, 2);
        glstrain_glob(1, 0) = glstrain_glob(0, 1);
        glstrain_glob(2, 0) = glstrain_glob(0, 2);
        glstrain_glob(2, 1) = glstrain_glob(1, 2);

        // pushforward of gl strains to ea strains
        CORE::LINALG::Matrix<noddof_, noddof_> euler_almansi(true);
        mem_GLtoEA(glstrain_glob, defgrd_glob, euler_almansi);

        (*elestrain)(gp, 0) = euler_almansi(0, 0);
        (*elestrain)(gp, 1) = euler_almansi(1, 1);
        (*elestrain)(gp, 2) = euler_almansi(2, 2);
        (*elestrain)(gp, 3) = euler_almansi(0, 1);
        (*elestrain)(gp, 4) = euler_almansi(1, 2);
        (*elestrain)(gp, 5) = euler_almansi(0, 2);
      }
      break;
      // Logarithmic strains
      case INPAR::STR::strain_log:
      {
        if (elestrain == nullptr) dserror("strain data not available");

        // the Eularian logarithmic strain is defined as the natural logarithm of the left stretch
        // tensor [1,2]: e_{log} = e_{hencky} = ln (\mathbf{V}) = \sum_{i=1}^3 (ln \lambda_i)
        // \mathbf{n}_i \otimes \mathbf{n}_i References: [1] H. Xiao, Beijing, China, O. T. Bruhns
        // and A. Meyers (1997) Logarithmic strain, logarithmic spin and logarithmic rate, Eq. 5 [2]
        // Caminero et al. (2011) Modeling large strain anisotropic elasto-plasticity with
        // logarithmic strain and stress measures, Eq. 70

        // transform local cauchygreen to global coordinates
        CORE::LINALG::Matrix<noddof_, noddof_> cauchygreen_glob(true);
        CORE::LINALG::TENSOR::TensorRotation<3>(Q_localToGlobal, cauchygreen_loc, cauchygreen_glob);

        // eigenvalue decomposition (from elasthyper.cpp)
        CORE::LINALG::Matrix<noddof_, noddof_> prstr2(true);  // squared principal stretches
        CORE::LINALG::Matrix<noddof_, 1> prstr(true);         // principal stretch
        CORE::LINALG::Matrix<noddof_, noddof_> prdir(true);   // principal directions
        CORE::LINALG::SYEV(cauchygreen_glob, prstr2, prdir);

        // THE principal stretches
        for (int al = 0; al < 3; ++al) prstr(al) = std::sqrt(prstr2(al, al));

        // populating the logarithmic strain matrix
        CORE::LINALG::Matrix<noddof_, noddof_> lnv(true);

        // checking if cauchy green is correctly determined to ensure eigenvectors in correct
        // direction i.e. a flipped eigenvector is also a valid solution C = \sum_{i=1}^3
        // (\lambda_i^2) \mathbf{n}_i \otimes \mathbf{n}_i
        CORE::LINALG::Matrix<noddof_, noddof_> tempCG(true);

        for (int k = 0; k < 3; ++k)
        {
          double n_00, n_01, n_02, n_11, n_12, n_22 = 0.0;

          n_00 = prdir(0, k) * prdir(0, k);
          n_01 = prdir(0, k) * prdir(1, k);
          n_02 = prdir(0, k) * prdir(2, k);
          n_11 = prdir(1, k) * prdir(1, k);
          n_12 = prdir(1, k) * prdir(2, k);
          n_22 = prdir(2, k) * prdir(2, k);

          // only compute the symmetric components from a single eigenvector,
          // because eigenvalue directions are not consistent (it can be flipped)
          tempCG(0, 0) += (prstr(k)) * (prstr(k)) * n_00;
          tempCG(0, 1) += (prstr(k)) * (prstr(k)) * n_01;
          tempCG(0, 2) += (prstr(k)) * (prstr(k)) * n_02;
          tempCG(1, 0) += (prstr(k)) * (prstr(k)) * n_01;  // symmetry
          tempCG(1, 1) += (prstr(k)) * (prstr(k)) * n_11;
          tempCG(1, 2) += (prstr(k)) * (prstr(k)) * n_12;
          tempCG(2, 0) += (prstr(k)) * (prstr(k)) * n_02;  // symmetry
          tempCG(2, 1) += (prstr(k)) * (prstr(k)) * n_12;  // symmetry
          tempCG(2, 2) += (prstr(k)) * (prstr(k)) * n_22;

          // Computation of the Logarithmic strain tensor

          lnv(0, 0) += (std::log(prstr(k))) * n_00;
          lnv(0, 1) += (std::log(prstr(k))) * n_01;
          lnv(0, 2) += (std::log(prstr(k))) * n_02;
          lnv(1, 0) += (std::log(prstr(k))) * n_01;  // symmetry
          lnv(1, 1) += (std::log(prstr(k))) * n_11;
          lnv(1, 2) += (std::log(prstr(k))) * n_12;
          lnv(2, 0) += (std::log(prstr(k))) * n_02;  // symmetry
          lnv(2, 1) += (std::log(prstr(k))) * n_12;  // symmetry
          lnv(2, 2) += (std::log(prstr(k))) * n_22;
        }

        // compare CG computed with deformation gradient with CG computed
        // with eigenvalues and -vectors to determine/ensure the correct
        // orientation of the eigen vectors
        CORE::LINALG::Matrix<noddof_, noddof_> diffCG(true);

        for (int i = 0; i < 3; ++i)
        {
          for (int j = 0; j < 3; ++j)
          {
            diffCG(i, j) = cauchygreen_glob(i, j) - tempCG(i, j);
            // the solution to this problem is to evaluate the cauchygreen tensor with
            // tempCG computed with every combination of eigenvector orientations -- up to nine
            // comparisons
            if (diffCG(i, j) > 1e-10)
              dserror(
                  "eigenvector orientation error with the diffCG giving problems: %10.5e \n BUILD "
                  "SOLUTION TO FIX IT",
                  diffCG(i, j));
          }
        }

        (*elestrain)(gp, 0) = lnv(0, 0);
        (*elestrain)(gp, 1) = lnv(1, 1);
        (*elestrain)(gp, 2) = lnv(2, 2);
        (*elestrain)(gp, 3) = lnv(0, 1);
        (*elestrain)(gp, 4) = lnv(1, 2);
        (*elestrain)(gp, 5) = lnv(0, 2);
      }
      break;
      // no strain output
      case INPAR::STR::strain_none:
        break;
      default:
        dserror("requested strain type not available");
        break;
    }

    /*===============================================================================*
     | return gp stresses (only in case of stress/strain output)                     |
     *===============================================================================*/
    switch (iostress)
    {
      // 2nd Piola-Kirchhoff stresses
      case INPAR::STR::stress_2pk:
      {
        if (elestress == nullptr) dserror("stress data not available");

        // 2nd Piola-Kirchhoff stress in tensor notation, plane stress meaning entries in 2i and i2
        // are zero for i=0,1,2
        CORE::LINALG::Matrix<noddof_, noddof_> pkstressM_local(true);
        pkstressM_local(0, 0) = pk2red_loc(0);
        pkstressM_local(1, 1) = pk2red_loc(1);
        pkstressM_local(0, 1) = pk2red_loc(2);
        pkstressM_local(1, 0) = pk2red_loc(2);

        // determine 2nd Piola-Kirchhoff stresses in global coordinates
        CORE::LINALG::Matrix<noddof_, noddof_> pkstress_glob(true);
        CORE::LINALG::TENSOR::TensorRotation<3>(Q_localToGlobal, pkstressM_local, pkstress_glob);

        (*elestress)(gp, 0) = pkstress_glob(0, 0);
        (*elestress)(gp, 1) = pkstress_glob(1, 1);
        (*elestress)(gp, 2) = pkstress_glob(2, 2);
        (*elestress)(gp, 3) = pkstress_glob(0, 1);
        (*elestress)(gp, 4) = pkstress_glob(1, 2);
        (*elestress)(gp, 5) = pkstress_glob(0, 2);
      }
      break;
      // Cauchy stresses
      case INPAR::STR::stress_cauchy:
      {
        if (elestress == nullptr) dserror("stress data not available");

        // 2nd Piola-Kirchhoff stress in tensor notation, plane stress meaning entries in 2i and i2
        // are zero for i=0,1,2
        CORE::LINALG::Matrix<noddof_, noddof_> pkstressM_loc(true);
        pkstressM_loc(0, 0) = pk2red_loc(0);
        pkstressM_loc(1, 1) = pk2red_loc(1);
        pkstressM_loc(0, 1) = pk2red_loc(2);
        pkstressM_loc(1, 0) = pk2red_loc(2);

        // determine 2nd Piola-Kirchhoff stresses in global coordinates
        CORE::LINALG::Matrix<noddof_, noddof_> pkstress_glob(true);
        CORE::LINALG::TENSOR::TensorRotation<3>(Q_localToGlobal, pkstressM_loc, pkstress_glob);

        CORE::LINALG::Matrix<noddof_, noddof_> cauchy_glob(true);
        mem_PK2toCauchy(pkstress_glob, defgrd_glob, cauchy_glob);

        (*elestress)(gp, 0) = cauchy_glob(0, 0);
        (*elestress)(gp, 1) = cauchy_glob(1, 1);
        (*elestress)(gp, 2) = cauchy_glob(2, 2);
        (*elestress)(gp, 3) = cauchy_glob(0, 1);
        (*elestress)(gp, 4) = cauchy_glob(1, 2);
        (*elestress)(gp, 5) = cauchy_glob(0, 2);
      }
      break;
      // no stress output
      case INPAR::STR::stress_none:
        break;
      default:
        dserror("requested stress type not available");
        break;
    }
  }
  return;

}  // DRT::ELEMENTS::Membrane::membrane_nlnstiffmass

/*----------------------------------------------------------------------*
 |  Return names of visualization data (public)                fb 09/15 |
 *----------------------------------------------------------------------*/
template <CORE::FE::CellType distype>
void DRT::ELEMENTS::Membrane<distype>::VisNames(std::map<std::string, int>& names)
{
  std::string result_thickness = "thickness";

  names[result_thickness] = 1;


  SolidMaterial()->VisNames(names);

  return;

}  // DRT::ELEMENTS::Membrane::VisNames

/*----------------------------------------------------------------------*
 |  Return visualization data (public)                     fbraeu 06/16 |
 *----------------------------------------------------------------------*/
template <CORE::FE::CellType distype>
bool DRT::ELEMENTS::Membrane<distype>::VisData(const std::string& name, std::vector<double>& data)
{
  // Put the owner of this element into the file (use base class method for this)
  if (DRT::Element::VisData(name, data)) return true;

  if (name == "thickness")
  {
    if (data.size() != 1) dserror("size mismatch");
    for (int gp = 0; gp < intpoints_.nquad; gp++)
    {
      data[0] += cur_thickness_[gp];
    }
    data[0] = data[0] / intpoints_.nquad;

    return true;
  }

  return SolidMaterial()->VisData(name, data, intpoints_.nquad, this->Id());

}  // DRT::ELEMENTS::Membrane::VisData

/*----------------------------------------------------------------------*
 |  get reference and current configuration                fbraeu 06/16 |
 *----------------------------------------------------------------------*/
template <CORE::FE::CellType distype>
void DRT::ELEMENTS::Membrane<distype>::mem_configuration(const std::vector<double>& disp,
    CORE::LINALG::Matrix<numnod_, noddof_>& xrefe, CORE::LINALG::Matrix<numnod_, noddof_>& xcurr)
{
  // get reference configuration and determine current configuration
  DRT::Node** nodes = Nodes();
  if (!nodes) dserror("Nodes() returned null pointer");

  for (int i = 0; i < numnod_; ++i)
  {
    const auto& x = nodes[i]->X();
    xrefe(i, 0) = x[0];
    xrefe(i, 1) = x[1];
    xrefe(i, 2) = x[2];

    xcurr(i, 0) = xrefe(i, 0) + disp[i * noddof_ + 0];
    xcurr(i, 1) = xrefe(i, 1) + disp[i * noddof_ + 1];
    xcurr(i, 2) = xrefe(i, 2) + disp[i * noddof_ + 2];
  }

  return;

}  // DRT::ELEMENTS::Membrane::mem_configuration

/*------------------------------------------------------------------------------------------------------*
 |  introduce an orthonormal base in the undeformed configuration at current Gauss point   fbraeu
 06/16 |
 *------------------------------------------------------------------------------------------------------*/
template <CORE::FE::CellType distype>
void DRT::ELEMENTS::Membrane<distype>::mem_orthonormalbase(
    const CORE::LINALG::Matrix<numnod_, noddof_>& xrefe,
    const CORE::LINALG::Matrix<numnod_, noddof_>& xcurr,
    const CORE::LINALG::Matrix<numdim_, numnod_>& derivs,
    CORE::LINALG::Matrix<numdim_, numnod_>& derivs_ortho, double& G1G2_cn,
    CORE::LINALG::Matrix<noddof_, 1>& dXds1, CORE::LINALG::Matrix<noddof_, 1>& dXds2,
    CORE::LINALG::Matrix<noddof_, 1>& dxds1, CORE::LINALG::Matrix<noddof_, 1>& dxds2,
    CORE::LINALG::Matrix<noddof_, noddof_>& Q_localToGlobal) const
{
  /*===============================================================================*
   | introduce an orthonormal base in the undeformed configuration as proposed in: |
   | Gruttmann, "Theory and finite element formulation of rubberlike membrane      |
   | shells using principal stretches", 1992                                       |
   *===============================================================================*/

  CORE::LINALG::Matrix<noddof_, numdim_> G12(true);
  G12.MultiplyTT(1.0, xrefe, derivs, 0.0);

  // G1 and G2 Gruttmann1992 equation (43)
  CORE::LINALG::Matrix<noddof_, 1> G1(true);
  G1(0) = G12(0, 0);
  G1(1) = G12(1, 0);
  G1(2) = G12(2, 0);

  CORE::LINALG::Matrix<noddof_, 1> G2(true);
  G2(0) = G12(0, 1);
  G2(1) = G12(1, 1);
  G2(2) = G12(2, 1);

  // cross product G1xG2
  CORE::LINALG::Matrix<noddof_, 1> G1G2_cross(true);
  G1G2_cross(0) = G1(1) * G2(2) - G1(2) * G2(1);
  G1G2_cross(1) = G1(2) * G2(0) - G1(0) * G2(2);
  G1G2_cross(2) = G1(0) * G2(1) - G1(1) * G2(0);

  // 2 norm of vectors
  G1G2_cn = G1G2_cross.Norm2();
  double G1_n = G1.Norm2();

  // Gruttmann1992 equation (44), orthonormal base vectors
  CORE::LINALG::Matrix<noddof_, 1> tn(true);
  tn(0) = G1G2_cross(0) / G1G2_cn;
  tn(1) = G1G2_cross(1) / G1G2_cn;
  tn(2) = G1G2_cross(2) / G1G2_cn;

  CORE::LINALG::Matrix<noddof_, 1> t1(true);
  t1(0) = G1(0) / G1_n;
  t1(1) = G1(1) / G1_n;
  t1(2) = G1(2) / G1_n;

  CORE::LINALG::Matrix<noddof_, 1> t2(true);
  t2(0) = tn(1) * t1(2) - tn(2) * t1(1);
  t2(1) = tn(2) * t1(0) - tn(0) * t1(2);
  t2(2) = tn(0) * t1(1) - tn(1) * t1(0);

  CORE::LINALG::Matrix<noddof_, numdim_> t12(true);
  t12(0, 0) = t1(0);
  t12(1, 0) = t1(1);
  t12(2, 0) = t1(2);
  t12(0, 1) = t2(0);
  t12(1, 1) = t2(1);
  t12(2, 1) = t2(2);

  // Jacobian transformation matrix and its inverse, Gruttmann1992 equation (44b)
  // for the Trafo from local membrane orthonormal coordinates to global coordinates
  // It is not the Jacobian for the Trafo from the parameter space xi, eta to the global coords!
  CORE::LINALG::Matrix<numdim_, numdim_> J(true);
  J.MultiplyTN(1.0, G12, t12, 0.0);

  CORE::LINALG::Matrix<numdim_, numdim_> Jinv(true);
  Jinv.Invert(J);

  // calclate derivatives of shape functions in orthonormal base, Gruttmann1992 equation (42)
  derivs_ortho.Multiply(1.0, Jinv, derivs, 0.0);

  // derivative of the reference position wrt the orthonormal base
  CORE::LINALG::Matrix<noddof_, numdim_> dXds(true);
  dXds.MultiplyTT(1.0, xrefe, derivs_ortho, 0.0);

  dXds1(0) = dXds(0, 0);
  dXds1(1) = dXds(1, 0);
  dXds1(2) = dXds(2, 0);

  dXds2(0) = dXds(0, 1);
  dXds2(1) = dXds(1, 1);
  dXds2(2) = dXds(2, 1);

  // derivative of the current position wrt the orthonormal base
  CORE::LINALG::Matrix<noddof_, numdim_> dxds(true);
  dxds.MultiplyTT(1.0, xcurr, derivs_ortho, 0.0);

  dxds1(0) = dxds(0, 0);
  dxds1(1) = dxds(1, 0);
  dxds1(2) = dxds(2, 0);

  dxds2(0) = dxds(0, 1);
  dxds2(1) = dxds(1, 1);
  dxds2(2) = dxds(2, 1);

  // determine Trafo from local membrane orthonormal coordinates to global coordinates
  Q_localToGlobal(0, 0) = t1(0);
  Q_localToGlobal(1, 0) = t1(1);
  Q_localToGlobal(2, 0) = t1(2);
  Q_localToGlobal(0, 1) = t2(0);
  Q_localToGlobal(1, 1) = t2(1);
  Q_localToGlobal(2, 1) = t2(2);
  Q_localToGlobal(0, 2) = tn(0);
  Q_localToGlobal(1, 2) = tn(1);
  Q_localToGlobal(2, 2) = tn(2);

  return;

}  // DRT::ELEMENTS::Membrane::mem_orthonormalbase

/*-------------------------------------------------------------------------------------------------*
 |  pushforward of 2nd PK stresses to Cauchy stresses at gp                           fbraeu 06/16 |
 *-------------------------------------------------------------------------------------------------*/
template <CORE::FE::CellType distype>
void DRT::ELEMENTS::Membrane<distype>::mem_PK2toCauchy(
    const CORE::LINALG::Matrix<noddof_, noddof_>& pkstress_global,
    const CORE::LINALG::Matrix<noddof_, noddof_>& defgrd,
    CORE::LINALG::Matrix<noddof_, noddof_>& cauchy) const
{
  // calculate the Jacobi-deterinant
  const double detF = defgrd.Determinant();

  // check determinant of deformation gradient
  if (detF == 0) dserror("Zero Determinant of Deformation Gradient.");

  // determine the cauchy stresses
  CORE::LINALG::Matrix<noddof_, noddof_> temp;
  temp.Multiply((1.0 / detF), defgrd, pkstress_global, 0.0);
  cauchy.MultiplyNT(1.0, temp, defgrd, 1.0);

  return;

}  // DRT::ELEMENTS::Membrane::mem_PK2toCauchy

/*-------------------------------------------------------------------------------------------------*
 |  pushforward of Green-Lagrange to Euler-Almansi strains at gp                      fbraeu 06/16 |
 *-------------------------------------------------------------------------------------------------*/
template <CORE::FE::CellType distype>
void DRT::ELEMENTS::Membrane<distype>::mem_GLtoEA(
    const CORE::LINALG::Matrix<noddof_, noddof_>& glstrain_global,
    const CORE::LINALG::Matrix<noddof_, noddof_>& defgrd,
    CORE::LINALG::Matrix<noddof_, noddof_>& euler_almansi) const
{
  // check determinant of deformation gradient
  if (defgrd.Determinant() == 0)
    dserror("Inverse of Deformation Gradient can not be calcualated due to a zero Determinant.");

  // inverse of deformation gradient
  CORE::LINALG::Matrix<noddof_, noddof_> invdefgrd(true);
  invdefgrd.Invert(defgrd);

  // determine the euler-almansi strains
  CORE::LINALG::Matrix<noddof_, noddof_> temp;
  temp.Multiply(1.0, glstrain_global, invdefgrd, 0.0);
  euler_almansi.MultiplyTN(1.0, invdefgrd, temp, 1.0);

  return;

}  // DRT::ELEMENTS::Membrane::mem_GLtoEA

/*-------------------------------------------------------------------------------------------------*
 |  determine deformation gradient in global coordinates                              fbraeu 06/16 |
 *-------------------------------------------------------------------------------------------------*/
template <CORE::FE::CellType distype>
void DRT::ELEMENTS::Membrane<distype>::mem_defgrd_global(
    const CORE::LINALG::Matrix<noddof_, 1>& dXds1, const CORE::LINALG::Matrix<noddof_, 1>& dXds2,
    const CORE::LINALG::Matrix<noddof_, 1>& dxds1, const CORE::LINALG::Matrix<noddof_, 1>& dxds2,
    const double& lambda3, CORE::LINALG::Matrix<noddof_, noddof_>& defgrd_glob) const
{
  // clear
  defgrd_glob.Clear();

  // determine cross product x,1 x x,2
  CORE::LINALG::Matrix<noddof_, 1> xcurr_cross(true);
  xcurr_cross(0) = dxds1(1) * dxds2(2) - dxds1(2) * dxds2(1);
  xcurr_cross(1) = dxds1(2) * dxds2(0) - dxds1(0) * dxds2(2);
  xcurr_cross(2) = dxds1(0) * dxds2(1) - dxds1(1) * dxds2(0);

  // normalize the cross product for the current configuration
  xcurr_cross.Scale(1.0 / xcurr_cross.Norm2());

  // determine cross product X,1 x X,2, has unit length due to orthonormal basis
  CORE::LINALG::Matrix<noddof_, 1> xrefe_cross(true);
  xrefe_cross(0) = dXds1(1) * dXds2(2) - dXds1(2) * dXds2(1);
  xrefe_cross(1) = dXds1(2) * dXds2(0) - dXds1(0) * dXds2(2);
  xrefe_cross(2) = dXds1(0) * dXds2(1) - dXds1(1) * dXds2(0);

  defgrd_glob.MultiplyNT(1.0, dxds1, dXds1, 0.0);
  defgrd_glob.MultiplyNT(1.0, dxds2, dXds2, 1.0);
  // scale third dimension by sqrt(rcg33), that equals the principle stretch lambda_3
  defgrd_glob.MultiplyNT(lambda3, xcurr_cross, xrefe_cross, 1.0);

  return;

}  // DRT::ELEMENTS::Membrane::mem_defgrd_global

/*-------------------------------------------------------------------------------------------------*
 |  determine extrapolation matrix                                                    sfuchs 02/18 |
 *-------------------------------------------------------------------------------------------------*/
template <CORE::FE::CellType distype>
CORE::LINALG::Matrix<CORE::FE::num_nodes<distype>, THR::DisTypeToNumGaussPoints<distype>::nquad>
DRT::ELEMENTS::Membrane<distype>::mem_extrapolmat() const
{
  // extrapolation matrix
  // note: equal for all elements of the same discretization type
  CORE::LINALG::Matrix<numnod_, numgpt_post_> extrapol;

  // check for correct gaussrule
  if (intpoints_.nquad != numgpt_post_)
    dserror(
        "number of gauss points of gaussrule_ does not match numgpt_post_ used for postprocessing");

  // allocate vector for shape functions and matrix for derivatives at gp
  CORE::LINALG::Matrix<numnod_, 1> shapefcts(true);

  // loop over the nodes and gauss points
  // interpolation matrix, inverted later to be the extrapolation matrix
  for (int nd = 0; nd < numnod_; ++nd)
  {
    // gaussian coordinates
    const double e1 = intpoints_.qxg[nd][0];
    const double e2 = intpoints_.qxg[nd][1];

    // shape functions for the extrapolated coordinates
    CORE::LINALG::Matrix<numgpt_post_, 1> funct;
    CORE::FE::shape_function_2D(funct, e1, e2, Shape());

    for (int i = 0; i < numgpt_post_; ++i) extrapol(nd, i) = funct(i);
  }

  // fixedsizesolver for inverting extrapol
  CORE::LINALG::FixedSizeSerialDenseSolver<numnod_, numgpt_post_, 1> solver;
  solver.SetMatrix(extrapol);
  int err = solver.Invert();
  if (err != 0.) dserror("Matrix extrapol is not invertible");

  return extrapol;

}  // DRT::ELEMENTS::Membrane::mem_extrapolmat

/*---------------------------------------------------------------------------------------------*
 |  Update history variables (e.g. remodeling of fiber directions) (protected)      braeu 07/16|
 *---------------------------------------------------------------------------------------------*/
template <CORE::FE::CellType distype>
void DRT::ELEMENTS::Membrane<distype>::Update_element(
    std::vector<double>& disp, Teuchos::ParameterList& params, Teuchos::RCP<MAT::Material> mat)
{
  // Calculate current deformation gradient
  if (SolidMaterial()->UsesExtendedUpdate())
  {
    // get reference configuration and determine current configuration
    CORE::LINALG::Matrix<numnod_, noddof_> xrefe(true);
    CORE::LINALG::Matrix<numnod_, noddof_> xcurr(true);

    mem_configuration(disp, xrefe, xcurr);

    /*===============================================================================*
     | loop over the gauss points                                                    |
     *===============================================================================*/

    // allocate vector for shape functions and matrix for derivatives at gp
    CORE::LINALG::Matrix<numdim_, numnod_> derivs(true);

    for (int gp = 0; gp < intpoints_.nquad; ++gp)
    {
      // get gauss points from integration rule
      double xi_gp = intpoints_.qxg[gp][0];
      double eta_gp = intpoints_.qxg[gp][1];

      // get derivatives in the plane of the element
      CORE::FE::shape_function_2D_deriv1(derivs, xi_gp, eta_gp, Shape());

      /*===============================================================================*
       | orthonormal base (t1,t2,tn) in the undeformed configuration at current GP     |
       *===============================================================================*/

      CORE::LINALG::Matrix<numdim_, numnod_> derivs_ortho(true);
      double G1G2_cn;
      CORE::LINALG::Matrix<noddof_, 1> dXds1(true);
      CORE::LINALG::Matrix<noddof_, 1> dXds2(true);
      CORE::LINALG::Matrix<noddof_, 1> dxds1(true);
      CORE::LINALG::Matrix<noddof_, 1> dxds2(true);
      CORE::LINALG::Matrix<noddof_, noddof_> Q_localToGlobal(true);

      mem_orthonormalbase(
          xrefe, xcurr, derivs, derivs_ortho, G1G2_cn, dXds1, dXds2, dxds1, dxds2, Q_localToGlobal);

      /*===============================================================================*
       | surface deformation gradient                                                  |
       *===============================================================================*/

      // surface deformation gradient in 3 dimensions in global coordinates
      CORE::LINALG::Matrix<noddof_, noddof_> defgrd_glob(true);
      CORE::LINALG::Matrix<noddof_, noddof_> defgrd_loc(true);

      // principle stretch in thickness direction
      double lambda3 = cur_thickness_[gp] / thickness_;

      // surface deformation gradient in 3 dimensions in global coordinates
      mem_defgrd_global(dXds1, dXds2, dxds1, dxds2, lambda3, defgrd_glob);

      CORE::LINALG::TENSOR::InverseTensorRotation<3>(Q_localToGlobal, defgrd_glob, defgrd_loc);

      auto material_local_coordinates =
          Teuchos::rcp_dynamic_cast<MAT::MembraneMaterialLocalCoordinates>(
              DRT::Element::Material());
      auto material_global_coordinates =
          Teuchos::rcp_dynamic_cast<MAT::MembraneMaterialGlobalCoordinates>(
              DRT::Element::Material());
      if (material_local_coordinates != Teuchos::null)
      {
        material_local_coordinates->UpdateMembrane(defgrd_loc, params, Q_localToGlobal, gp, Id());
      }
      else if (material_global_coordinates != Teuchos::null)
      {
        SolidMaterial()->Update(defgrd_glob, gp, params, Id());
      }
    }
  }

  SolidMaterial()->Update();

  return;
}

template class DRT::ELEMENTS::Membrane<CORE::FE::CellType::tri3>;
template class DRT::ELEMENTS::Membrane<CORE::FE::CellType::tri6>;
template class DRT::ELEMENTS::Membrane<CORE::FE::CellType::quad4>;
template class DRT::ELEMENTS::Membrane<CORE::FE::CellType::quad9>;

BACI_NAMESPACE_CLOSE
