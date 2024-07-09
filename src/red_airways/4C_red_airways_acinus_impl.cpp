/*---------------------------------------------------------------------*/
/*! \file

\brief Internal implementation of RedAcinus element. Methods implemented here
       are called by acinus_evaluate.cpp by Discret::ELEMENTS::RedAcinus::evaluate()
       with the corresponding action.


\level 3

*/
/*---------------------------------------------------------------------*/



#include "4C_red_airways_acinus_impl.hpp"

#include "4C_fem_discretization.hpp"
#include "4C_fem_general_extract_values.hpp"
#include "4C_fem_general_utils_fem_shapefunctions.hpp"
#include "4C_fem_general_utils_gder2.hpp"
#include "4C_global_data.hpp"
#include "4C_mat_maxwell_0d_acinus.hpp"
#include "4C_mat_par_bundle.hpp"
#include "4C_red_airways_elem_params.hpp"
#include "4C_red_airways_evaluation_data.hpp"
#include "4C_utils_function.hpp"
#include "4C_utils_function_of_time.hpp"

#include <fstream>
#include <iomanip>

FOUR_C_NAMESPACE_OPEN


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Discret::ELEMENTS::RedAcinusImplInterface* Discret::ELEMENTS::RedAcinusImplInterface::impl(
    Discret::ELEMENTS::RedAcinus* red_acinus)
{
  switch (red_acinus->shape())
  {
    case Core::FE::CellType::line2:
    {
      static AcinusImpl<Core::FE::CellType::line2>* acinus;
      if (acinus == nullptr)
      {
        acinus = new AcinusImpl<Core::FE::CellType::line2>;
      }
      return acinus;
    }
    default:
      FOUR_C_THROW(
          "shape %d (%d nodes) not supported", red_acinus->shape(), red_acinus->num_node());
      break;
  }
  return nullptr;
}


/*----------------------------------------------------------------------*
  | constructor (public)                                    ismail 01/10 |
 *----------------------------------------------------------------------*/
template <Core::FE::CellType distype>
Discret::ELEMENTS::AcinusImpl<distype>::AcinusImpl()
{
}


/*----------------------------------------------------------------------*
 |  calculate element matrix and right hand side (private)  ismail 01/10|
 |                                                                      |
 *----------------------------------------------------------------------*/
/*!
        \brief calculate element matrix and rhs

\param ele              (i) the element those matrix is calculated
        \param eqnp             (i) nodal volumetric flow rate at n+1
        \param evelnp           (i) nodal velocity at n+1
        \param eareanp          (i) nodal cross-sectional area at n+1
        \param eprenp           (i) nodal pressure at n+1
        \param estif            (o) element matrix to calculate
        \param eforce           (o) element rhs to calculate
        \param material         (i) acinus material/dimesion
        \param time             (i) current simulation time
        \param dt               (i) timestep
        */
template <Core::FE::CellType distype>
void Sysmat(Discret::ELEMENTS::RedAcinus* ele, Core::LinAlg::SerialDenseVector& epnp,
    Core::LinAlg::SerialDenseVector& epn, Core::LinAlg::SerialDenseVector& epnm,
    Core::LinAlg::SerialDenseMatrix& sysmat, Core::LinAlg::SerialDenseVector& rhs,
    Teuchos::RCP<const Core::Mat::Material> material, Discret::ReducedLung::ElemParams& params,
    double time, double dt)
{
  const auto acinus_params = ele->get_acinus_params();

  // Decide which acinus material should be used
  if ((material->material_type() == Core::Materials::m_0d_maxwell_acinus_neohookean) ||
      (material->material_type() == Core::Materials::m_0d_maxwell_acinus_exponential) ||
      (material->material_type() == Core::Materials::m_0d_maxwell_acinus_doubleexponential) ||
      (material->material_type() == Core::Materials::m_0d_maxwell_acinus_ogden))
  {
    const double VolAcinus = acinus_params.volume_relaxed;
    const double volAlvDuct = acinus_params.alveolar_duct_volume;
    const auto NumOfAcini = double(floor(VolAcinus / volAlvDuct));

    const Teuchos::RCP<Mat::Maxwell0dAcinus> acinus_mat =
        Teuchos::rcp_dynamic_cast<Mat::Maxwell0dAcinus>(ele->material());

    // Evaluate material law for acinus
    acinus_mat->evaluate(epnp, epn, epnm, sysmat, rhs, params, NumOfAcini, volAlvDuct, time, dt);
  }
  else
  {
    FOUR_C_THROW("Material law is not a valid reduced dimensional lung acinus material.");
    exit(1);
  }
}


/*----------------------------------------------------------------------*
 | evaluate (public)                                       ismail 01/10 |
 *----------------------------------------------------------------------*/
template <Core::FE::CellType distype>
int Discret::ELEMENTS::AcinusImpl<distype>::evaluate(RedAcinus* ele, Teuchos::ParameterList& params,
    Core::FE::Discretization& discretization, std::vector<int>& lm,
    Core::LinAlg::SerialDenseMatrix& elemat1_epetra,
    Core::LinAlg::SerialDenseMatrix& elemat2_epetra,
    Core::LinAlg::SerialDenseVector& elevec1_epetra,
    Core::LinAlg::SerialDenseVector& elevec2_epetra,
    Core::LinAlg::SerialDenseVector& elevec3_epetra, Teuchos::RCP<Core::Mat::Material> mat)
{
  const int elemVecdim = elevec1_epetra.length();

  Discret::ReducedLung::EvaluationData& evaluation_data =
      Discret::ReducedLung::EvaluationData::get();

  // Get control parameters for time integration
  // get time-step size
  const double dt = evaluation_data.dt;
  // get time
  const double time = evaluation_data.time;

  // Get control parameters for stabilization and higher-order elements (currently unused)
  // flag for higher order elements
  // bool higher_order_ele = ele->is_higher_order_element(distype);

  // Get all general state vectors: flow, pressure,
  Teuchos::RCP<const Epetra_Vector> pnp = discretization.get_state("pnp");
  Teuchos::RCP<const Epetra_Vector> pn = discretization.get_state("pn");
  Teuchos::RCP<const Epetra_Vector> pnm = discretization.get_state("pnm");

  Teuchos::RCP<const Epetra_Vector> ial = discretization.get_state("intr_ac_link");

  if (pnp == Teuchos::null || pn == Teuchos::null || pnm == Teuchos::null)
    FOUR_C_THROW("Cannot get state vectors 'pnp', 'pn', and/or 'pnm''");

  // Extract local values from the global vectors
  std::vector<double> mypnp(lm.size());
  Core::FE::ExtractMyValues(*pnp, mypnp, lm);

  // Extract local values from the global vectors
  std::vector<double> mypn(lm.size());
  Core::FE::ExtractMyValues(*pn, mypn, lm);

  // Extract local values from the global vectors
  std::vector<double> mypnm(lm.size());
  Core::FE::ExtractMyValues(*pnm, mypnm, lm);

  // Extract local values from the global vectors
  std::vector<double> myial(lm.size());
  Core::FE::ExtractMyValues(*ial, myial, lm);

  // Create objects for element arrays
  Core::LinAlg::SerialDenseVector epnp(elemVecdim);
  Core::LinAlg::SerialDenseVector epn(elemVecdim);
  Core::LinAlg::SerialDenseVector epnm(elemVecdim);
  for (int i = 0; i < elemVecdim; ++i)
  {
    // Split area and volumetric flow rate, insert into element arrays
    epnp(i) = mypnp[i];
    epn(i) = mypn[i];
    epnm(i) = mypnm[i];
  }

  double e_acin_e_vnp;
  double e_acin_e_vn;

  // Split area and volumetric flow rate, insert into element arrays
  e_acin_e_vnp = (*evaluation_data.acinar_vnp)[ele->lid()];
  e_acin_e_vn = (*evaluation_data.acinar_vn)[ele->lid()];

  // Get the volumetric flow rate from the previous time step
  Discret::ReducedLung::ElemParams elem_params;
  elem_params.qout_np = (*evaluation_data.qout_np)[ele->lid()];
  elem_params.qout_n = (*evaluation_data.qout_n)[ele->lid()];
  elem_params.qout_nm = (*evaluation_data.qout_nm)[ele->lid()];
  elem_params.qin_np = (*evaluation_data.qin_np)[ele->lid()];
  elem_params.qin_n = (*evaluation_data.qin_n)[ele->lid()];
  elem_params.qin_nm = (*evaluation_data.qin_nm)[ele->lid()];

  elem_params.acin_vnp = e_acin_e_vnp;
  elem_params.acin_vn = e_acin_e_vn;

  elem_params.lungVolume_np = evaluation_data.lungVolume_np;
  elem_params.lungVolume_n = evaluation_data.lungVolume_n;
  elem_params.lungVolume_nm = evaluation_data.lungVolume_nm;

  // Call routine for calculating element matrix and right hand side
  Sysmat<distype>(ele, epnp, epn, epnm, elemat1_epetra, elevec1_epetra, mat, elem_params, time, dt);

  // Put zeros on second line of matrix and rhs in case of interacinar linker
  if (myial[1] > 0.0)
  {
    elemat1_epetra(1, 0) = 0.0;
    elemat1_epetra(1, 1) = 0.0;
    elevec1_epetra(1) = 0.0;
  }

  return 0;
}


/*----------------------------------------------------------------------*
 |  calculate element matrix and right hand side (private)  ismail 01/10|
 *----------------------------------------------------------------------*/
template <Core::FE::CellType distype>
void Discret::ELEMENTS::AcinusImpl<distype>::initial(RedAcinus* ele, Teuchos::ParameterList& params,
    Core::FE::Discretization& discretization, std::vector<int>& lm,
    Teuchos::RCP<const Core::Mat::Material> material)
{
  const int myrank = discretization.get_comm().MyPID();

  Discret::ReducedLung::EvaluationData& evaluation_data =
      Discret::ReducedLung::EvaluationData::get();
  const auto acinus_params = ele->get_acinus_params();

  std::vector<int> lmstride;
  Teuchos::RCP<std::vector<int>> lmowner = Teuchos::rcp(new std::vector<int>);
  ele->location_vector(discretization, lm, *lmowner, lmstride);

  // Initialize the pressure vectors
  if (myrank == (*lmowner)[0])
  {
    int gid = lm[0];
    double val = 0.0;
    evaluation_data.p0np->ReplaceGlobalValues(1, &val, &gid);
    evaluation_data.p0n->ReplaceGlobalValues(1, &val, &gid);
    evaluation_data.p0nm->ReplaceGlobalValues(1, &val, &gid);
  }

  // Find the volume of an acinus element
  {
    int gid2 = ele->id();
    double acin_vol = acinus_params.volume_relaxed;
    evaluation_data.acini_e_volume->ReplaceGlobalValues(1, &acin_vol, &gid2);
  }

  // Get the generation numbers
  for (int i = 0; i < 2; i++)
  {
    if (ele->nodes()[i]->get_condition("RedAirwayEvalLungVolCond"))
    {
      // find the acinus condition
      int gid = ele->id();
      double val = 1.0;
      evaluation_data.acini_bc->ReplaceGlobalValues(1, &val, &gid);
    }
  }
  {
    int gid = ele->id();
    int generation = -1;
    double val = double(generation);
    evaluation_data.generations->ReplaceGlobalValues(1, &val, &gid);
  }

}  // AcinusImpl::Initial


/*----------------------------------------------------------------------*
 |  Evaluate the values of the degrees of freedom           ismail 01/10|
 |  at terminal nodes.                                                  |
 *----------------------------------------------------------------------*/
template <Core::FE::CellType distype>
void Discret::ELEMENTS::AcinusImpl<distype>::evaluate_terminal_bc(RedAcinus* ele,
    Teuchos::ParameterList& params, Core::FE::Discretization& discretization, std::vector<int>& lm,
    Core::LinAlg::SerialDenseVector& rhs, Teuchos::RCP<Core::Mat::Material> material)
{
  const int myrank = discretization.get_comm().MyPID();

  Discret::ReducedLung::EvaluationData& evaluation_data =
      Discret::ReducedLung::EvaluationData::get();

  // Get total time
  const double time = evaluation_data.time;

  // The number of nodes
  const int numnode = lm.size();

  Teuchos::RCP<const Epetra_Vector> pnp = discretization.get_state("pnp");

  if (pnp == Teuchos::null) FOUR_C_THROW("Cannot get state vectors 'pnp'");

  // Extract local values from the global vectors
  std::vector<double> mypnp(lm.size());
  Core::FE::ExtractMyValues(*pnp, mypnp, lm);

  // Create objects for element arrays
  Core::LinAlg::SerialDenseVector epnp(numnode);

  // Get all values at the last computed time step
  for (int i = 0; i < numnode; ++i)
  {
    // Split area and volumetric flow rate, insert into element arrays
    epnp(i) = mypnp[i];
  }

  /**
   * Resolve the BCs
   **/
  for (int i = 0; i < ele->num_node(); i++)
  {
    if (ele->nodes()[i]->owner() == myrank)
    {
      if (ele->nodes()[i]->get_condition("RedAirwayPrescribedCond") ||
          ele->nodes()[i]->get_condition("Art_redD_3D_CouplingCond") ||
          ele->nodes()[i]->get_condition("RedAcinusVentilatorCond"))
      {
        std::string Bc;
        double BCin = 0.0;
        if (ele->nodes()[i]->get_condition("RedAirwayPrescribedCond"))
        {
          Core::Conditions::Condition* condition =
              ele->nodes()[i]->get_condition("RedAirwayPrescribedCond");
          // Get the type of prescribed bc
          Bc = (condition->parameters().get<std::string>("boundarycond"));

          const auto* vals = &condition->parameters().get<std::vector<double>>("val");
          const auto* curve = &condition->parameters().get<std::vector<int>>("curve");
          const auto* functions = &condition->parameters().get<std::vector<int>>("funct");

          // Read in the value of the applied BC
          // Get factor of first CURVE
          double curvefac = 1.0;
          if ((*curve)[0] >= 0)
          {
            curvefac = Global::Problem::instance()
                           ->function_by_id<Core::UTILS::FunctionOfTime>((*curve)[0])
                           .evaluate(time);
            BCin = (*vals)[0] * curvefac;
          }
          else
          {
            FOUR_C_THROW("no boundary condition defined!");
            exit(1);
          }

          // Get factor of FUNCT
          int functnum = -1;
          if (functions)
            functnum = (*functions)[0];
          else
            functnum = -1;

          double functionfac = 0.0;
          if (functnum > 0)
          {
            functionfac = Global::Problem::instance()
                              ->function_by_id<Core::UTILS::FunctionOfSpaceTime>(functnum - 1)
                              .evaluate((ele->nodes()[i])->x().data(), time, 0);
          }

          // Get factor of second CURVE
          int curve2num = -1;
          double curve2fac = 1.0;
          if (curve) curve2num = (*curve)[1];
          if (curve2num >= 0)
            curve2fac = Global::Problem::instance()
                            ->function_by_id<Core::UTILS::FunctionOfTime>(curve2num)
                            .evaluate(time);

          // Add first_CURVE + FUNCTION * second_CURVE
          BCin += functionfac * curve2fac;

          // Get the local id of the node to whom the bc is prescribed
          int local_id = discretization.node_row_map()->LID(ele->nodes()[i]->id());
          if (local_id < 0)
          {
            FOUR_C_THROW("node (%d) doesn't exist on proc(%d)", ele->nodes()[i]->id(),
                discretization.get_comm().MyPID());
            exit(1);
          }
        }
        /**
         * For Art_redD_3D_CouplingCond bc
         **/
        else if (ele->nodes()[i]->get_condition("Art_redD_3D_CouplingCond"))
        {
          const Core::Conditions::Condition* condition =
              ele->nodes()[i]->get_condition("Art_redD_3D_CouplingCond");

          Teuchos::RCP<Teuchos::ParameterList> CoupledTo3DParams =
              params.get<Teuchos::RCP<Teuchos::ParameterList>>("coupling with 3D fluid params");
          // -----------------------------------------------------------------
          // If the parameter list is empty, then something is wrong!
          // -----------------------------------------------------------------
          if (CoupledTo3DParams.get() == nullptr)
          {
            FOUR_C_THROW(
                "Cannot prescribe a boundary condition from 3D to reduced D, if the parameters "
                "passed don't exist");
            exit(1);
          }

          // -----------------------------------------------------------------
          // Read in Condition type
          // -----------------------------------------------------------------
          //        Type = (condition->parameters().get<std::string>("CouplingType"));
          // -----------------------------------------------------------------
          // Read in coupling variable rescribed by the 3D simulation
          //
          //     In this case a map called map3D has the following form:
          //     +-----------------------------------------------------------+
          //     |           std::map< string               ,  double        >    |
          //     |     +------------------------------------------------+    |
          //     |     |  ID  | coupling variable name | variable value |    |
          //     |     +------------------------------------------------+    |
          //     |     |  1   |   flow1                |     0.12116    |    |
          //     |     +------+------------------------+----------------+    |
          //     |     |  2   |   pressure2            |    10.23400    |    |
          //     |     +------+------------------------+----------------+    |
          //     |     .  .   .   ....                 .     .......    .    |
          //     |     +------+------------------------+----------------+    |
          //     |     |  N   |   variableN            |    value(N)    |    |
          //     |     +------+------------------------+----------------+    |
          //     +-----------------------------------------------------------+
          // -----------------------------------------------------------------

          int ID = condition->parameters().get<int>("ConditionID");
          Teuchos::RCP<std::map<std::string, double>> map3D;
          map3D = CoupledTo3DParams->get<Teuchos::RCP<std::map<std::string, double>>>(
              "3D map of values");

          // find the applied boundary variable
          std::stringstream stringID;
          stringID << "_" << ID;
          for (std::map<std::string, double>::iterator itr = map3D->begin(); itr != map3D->end();
               itr++)
          {
            std::string VariableWithId = itr->first;
            size_t found;
            found = VariableWithId.rfind(stringID.str());
            if (found != std::string::npos)
            {
              Bc = std::string(VariableWithId, 0, found);
              BCin = itr->second;
              break;
            }
          }
        }
        /**
         * For RedAcinusVentilatorCond bc
         **/
        else if (ele->nodes()[i]->get_condition("RedAcinusVentilatorCond"))
        {
          Core::Conditions::Condition* condition =
              ele->nodes()[i]->get_condition("RedAcinusVentilatorCond");
          // Get the type of prescribed bc
          Bc = (condition->parameters().get<std::string>("phase1"));

          double period = condition->parameters().get<double>("period");
          double period1 = condition->parameters().get<double>("phase1_period");

          unsigned int phase_number = 0;

          if (fmod(time, period) > period1)
          {
            phase_number = 1;
            Bc = (condition->parameters().get<std::string>("phase2"));
          }

          const auto* curve = &condition->parameters().get<std::vector<int>>("curve");
          double curvefac = 1.0;
          const auto* vals = &condition->parameters().get<std::vector<double>>("val");

          // Read in the value of the applied BC
          if ((*curve)[phase_number] >= 0)
          {
            curvefac = Global::Problem::instance()
                           ->function_by_id<Core::UTILS::FunctionOfTime>((*curve)[phase_number])
                           .evaluate(time);
            BCin = (*vals)[phase_number] * curvefac;
          }
          else
          {
            FOUR_C_THROW("no boundary condition defined!");
            exit(1);
          }

          // Get the local id of the node to whom the bc is prescribed
          int local_id = discretization.node_row_map()->LID(ele->nodes()[i]->id());
          if (local_id < 0)
          {
            FOUR_C_THROW("node (%d) doesn't exist on proc(%d)", ele->nodes()[i]->id(),
                discretization.get_comm().MyPID());
            exit(1);
          }
        }
        else
        {
        }

        /**
         * For pressure or VolumeDependentPleuralPressure bc
         **/
        if (Bc == "pressure" || Bc == "VolumeDependentPleuralPressure")
        {
          if (Bc == "VolumeDependentPleuralPressure")
          {
            Core::Conditions::Condition* pplCond =
                ele->nodes()[i]->get_condition("RedAirwayVolDependentPleuralPressureCond");
            double Pp_np = 0.0;
            if (pplCond)
            {
              const auto* curve = &pplCond->parameters().get<std::vector<int>>("curve");
              double curvefac = 1.0;
              const auto* vals = &pplCond->parameters().get<std::vector<double>>("val");

              // Read in the value of the applied BC
              if ((*curve)[0] >= 0)
              {
                curvefac = Global::Problem::instance()
                               ->function_by_id<Core::UTILS::FunctionOfTime>((*curve)[0])
                               .evaluate(time);
              }

              // Get parameters for VolumeDependentPleuralPressure condition
              std::string ppl_Type = (pplCond->parameters().get<std::string>("TYPE"));
              double ap = pplCond->parameters().get<double>("P_PLEURAL_0");
              double bp = pplCond->parameters().get<double>("P_PLEURAL_LIN");
              double cp = pplCond->parameters().get<double>("P_PLEURAL_NONLIN");
              double dp = pplCond->parameters().get<double>("TAU");
              double RV = pplCond->parameters().get<double>("RV");
              double TLC = pplCond->parameters().get<double>("TLC");

              Discret::ReducedLung::EvaluationData& evaluation_data =
                  Discret::ReducedLung::EvaluationData::get();

              // Safety check: in case of polynomial TLC is not used
              if (((ppl_Type == "Linear_Polynomial") or (ppl_Type == "Nonlinear_Polynomial")) and
                  (TLC != 0.0))
              {
                FOUR_C_THROW(
                    "TLC is not used for the following type of VolumeDependentPleuralPressure BC: "
                    "%s.\n Set TLC = 0.0",
                    ppl_Type.c_str());
              }

              // Safety check: in case of Ogden TLC, P_PLEURAL_0, and P_PLEURAL_LIN
              if ((ppl_Type == "Nonlinear_Ogden") and
                  ((TLC != 0.0) or (ap != 0.0) or (bp != 0.0) or (dp == 0.0)))
              {
                FOUR_C_THROW(
                    "Parameters are not set correctly for Nonlinear_Ogden. Only P_PLEURAL_NONLIN, "
                    "TAU and RV are used. Set all others to zero. TAU is not allowed to be zero.");
              }

              if (ppl_Type == "Linear_Polynomial")
              {
                const double lungVolumenp = evaluation_data.lungVolume_n;
                Pp_np = ap + bp * (lungVolumenp - RV) + cp * pow((lungVolumenp - RV), dp);
              }
              else if (ppl_Type == "Linear_Exponential")
              {
                const double lungVolumenp = evaluation_data.lungVolume_n;
                const double TLCnp = (lungVolumenp - RV) / (TLC - RV);
                Pp_np = ap + bp * TLCnp + cp * exp(dp * TLCnp);
              }
              else if (ppl_Type == "Linear_Ogden")
              {
                const double lungVolumenp = evaluation_data.lungVolume_n;
                Pp_np = RV / lungVolumenp * cp / dp * (1 - pow(RV / lungVolumenp, dp));
              }
              else if (ppl_Type == "Nonlinear_Polynomial")
              {
                const double lungVolumenp = evaluation_data.lungVolume_np;
                Pp_np = ap + bp * (lungVolumenp - RV) + cp * pow((lungVolumenp - RV), dp);
              }
              else if (ppl_Type == "Nonlinear_Exponential")
              {
                const double lungVolumenp = evaluation_data.lungVolume_np;
                const double TLCnp = (lungVolumenp - RV) / (TLC - RV);
                Pp_np = ap + bp * TLCnp + cp * exp(dp * TLCnp);
              }
              else if (ppl_Type == "Nonlinear_Ogden")
              {
                const double lungVolumenp = evaluation_data.lungVolume_np;
                Pp_np = RV / lungVolumenp * cp / dp * (1 - pow(RV / lungVolumenp, dp));
              }
              else
              {
                FOUR_C_THROW("Unknown volume pleural pressure type: %s", ppl_Type.c_str());
              }
              Pp_np *= curvefac * ((*vals)[0]);
            }
            else
            {
              FOUR_C_THROW("No volume dependent pleural pressure condition was defined");
            }
            BCin += Pp_np;
          }

          Discret::ReducedLung::EvaluationData& evaluation_data =
              Discret::ReducedLung::EvaluationData::get();

          // Set pressure at node i
          int gid;
          double val;

          gid = lm[i];
          val = BCin;
          evaluation_data.bcval->ReplaceGlobalValues(1, &val, &gid);

          gid = lm[i];
          val = 1;
          evaluation_data.dbctog->ReplaceGlobalValues(1, &val, &gid);
        }
        /**
         * For flow bc
         **/
        else if (Bc == "flow")
        {
          // ----------------------------------------------------------
          // Since a node might belong to multiple elements then the
          // flow might be added to the rhs multiple time.
          // To fix this the flow is divided by the number of elements
          // (which is the number of branches). Thus the sum of the
          // final added values is the actual prescribed flow.
          // ----------------------------------------------------------
          int numOfElems = (ele->nodes()[i])->num_element();
          BCin /= double(numOfElems);
          rhs(i) += -BCin + rhs(i);
        }
        else
        {
          FOUR_C_THROW("prescribed [%s] is not defined for reduced acinuss", Bc.c_str());
          exit(1);
        }
      }
      /**
       * If the node is a terminal node, but no b.c is prescribed to it
       * then a zero output pressure is assumed
       **/
      else
      {
        if (ele->nodes()[i]->num_element() == 1)
        {
          // Get the local id of the node to whome the bc is prescribed
          int local_id = discretization.node_row_map()->LID(ele->nodes()[i]->id());
          if (local_id < 0)
          {
            FOUR_C_THROW("node (%d) doesn't exist on proc(%d)", ele->nodes()[i],
                discretization.get_comm().MyPID());
            exit(1);
          }

          Discret::ReducedLung::EvaluationData& evaluation_data =
              Discret::ReducedLung::EvaluationData::get();

          // Set pressure=0.0 at node i
          int gid;
          double val;

          gid = lm[i];
          val = 0.0;
          evaluation_data.bcval->ReplaceGlobalValues(1, &val, &gid);

          gid = lm[i];
          val = 1;
          evaluation_data.dbctog->ReplaceGlobalValues(1, &val, &gid);
        }
      }  // END of if there is no BC but the node still is at the terminal

    }  // END of if node is available on this processor
  }    // End of node i has a condition
}


/*----------------------------------------------------------------------*
 |  Calculate flowrate at current iteration step and the    ismail 01/10|
 |  correponding acinus volume via dV = 0.5*(qnp+qn)*dt                 |
 *----------------------------------------------------------------------*/
template <Core::FE::CellType distype>
void Discret::ELEMENTS::AcinusImpl<distype>::calc_flow_rates(RedAcinus* ele,
    Teuchos::ParameterList& params, Core::FE::Discretization& discretization, std::vector<int>& lm,
    Teuchos::RCP<Core::Mat::Material> material)

{
  const int elemVecdim = lm.size();

  Discret::ReducedLung::EvaluationData& evaluation_data =
      Discret::ReducedLung::EvaluationData::get();
  const auto acinus_params = ele->get_acinus_params();

  // Get control parameters for time integration
  // Get time-step size
  const double dt = evaluation_data.dt;
  // Get time
  const double time = evaluation_data.time;

  // Get control parameters for stabilization and higher-order elements
  // flag for higher order elements
  //  bool higher_order_ele = ele->is_higher_order_element(distype);

  // Get all general state vectors: flow, pressure,
  Teuchos::RCP<const Epetra_Vector> pnp = discretization.get_state("pnp");
  Teuchos::RCP<const Epetra_Vector> pn = discretization.get_state("pn");
  Teuchos::RCP<const Epetra_Vector> pnm = discretization.get_state("pnm");

  if (pnp == Teuchos::null || pn == Teuchos::null || pnm == Teuchos::null)
    FOUR_C_THROW("Cannot get state vectors 'pnp', 'pn', and/or 'pnm''");

  // Extract local values from the global vectors
  std::vector<double> mypnp(lm.size());
  Core::FE::ExtractMyValues(*pnp, mypnp, lm);

  // Extract local values from the global vectors
  std::vector<double> mypn(lm.size());
  Core::FE::ExtractMyValues(*pn, mypn, lm);

  // Extract local values from the global vectors
  std::vector<double> mypnm(lm.size());
  Core::FE::ExtractMyValues(*pnm, mypnm, lm);

  // Create objects for element arrays
  Core::LinAlg::SerialDenseVector epnp(elemVecdim);
  Core::LinAlg::SerialDenseVector epn(elemVecdim);
  Core::LinAlg::SerialDenseVector epnm(elemVecdim);
  for (int i = 0; i < elemVecdim; ++i)
  {
    // Split area and volumetric flow rate, insert into element arrays
    epnp(i) = mypnp[i];
    epn(i) = mypn[i];
    epnm(i) = mypnm[i];
  }

  double e_acin_vnp = 0.0;
  double e_acin_vn = 0.0;

  for (int i = 0; i < elemVecdim; ++i)
  {
    // Split area and volumetric flow rate, insert into element arrays
    e_acin_vnp = (*evaluation_data.acinar_vnp)[ele->lid()];
    e_acin_vn = (*evaluation_data.acinar_vn)[ele->lid()];
  }

  // Get the volumetric flow rate from the previous time step
  Discret::ReducedLung::ElemParams elem_params;
  elem_params.qout_np = (*evaluation_data.qout_np)[ele->lid()];
  elem_params.qout_n = (*evaluation_data.qout_n)[ele->lid()];
  elem_params.qout_nm = (*evaluation_data.qout_nm)[ele->lid()];
  elem_params.qin_np = (*evaluation_data.qin_np)[ele->lid()];
  elem_params.qin_n = (*evaluation_data.qin_n)[ele->lid()];
  elem_params.qin_nm = (*evaluation_data.qin_nm)[ele->lid()];

  elem_params.acin_vnp = e_acin_vnp;
  elem_params.acin_vn = e_acin_vn;

  Core::LinAlg::SerialDenseMatrix sysmat(elemVecdim, elemVecdim, true);
  Core::LinAlg::SerialDenseVector rhs(elemVecdim);

  // Call routine for calculating element matrix and right hand side
  Sysmat<distype>(ele, epnp, epn, epnm, sysmat, rhs, material, elem_params, time, dt);

  double qn = (*evaluation_data.qin_n)[ele->lid()];
  double qnp = -1.0 * (sysmat(0, 0) * epnp(0) + sysmat(0, 1) * epnp(1) - rhs(0));

  int gid = ele->id();

  evaluation_data.qin_np->ReplaceGlobalValues(1, &qnp, &gid);
  evaluation_data.qout_np->ReplaceGlobalValues(1, &qnp, &gid);

  // Calculate the new volume of the acinus due to the incoming flow; 0.5*(qnp+qn)*dt
  {
    double acinus_volume = e_acin_vn;
    acinus_volume += 0.5 * (qnp + qn) * dt;
    evaluation_data.acinar_vnp->ReplaceGlobalValues(1, &acinus_volume, &gid);

    // Calculate correponding acinar strain
    const double vo = acinus_params.volume_relaxed;
    double avs_np = (acinus_volume - vo) / vo;
    evaluation_data.acinar_vnp_strain->ReplaceGlobalValues(1, &avs_np, &gid);
  }
}


/*----------------------------------------------------------------------*
 |  Calculate element volume                                ismail 07/13|
 |                                                                      |
 *----------------------------------------------------------------------*/
template <Core::FE::CellType distype>
void Discret::ELEMENTS::AcinusImpl<distype>::calc_elem_volume(RedAcinus* ele,
    Teuchos::ParameterList& params, Core::FE::Discretization& discretization, std::vector<int>& lm,
    Teuchos::RCP<Core::Mat::Material> material)

{
  Discret::ReducedLung::EvaluationData& evaluation_data =
      Discret::ReducedLung::EvaluationData::get();

  // Get acinus size
  double evolnp = (*evaluation_data.elemVolumenp)[ele->lid()];

  // Get element global ID
  int gid = ele->id();

  // Update elem
  evaluation_data.elemVolumenp->ReplaceGlobalValues(1, &evolnp, &gid);

  // calculate and update element radius
  double eRadiusnp = std::pow(evolnp * 0.75 * M_1_PI, 1.0 / 3.0);
  evaluation_data.elemRadiusnp->ReplaceGlobalValues(1, &eRadiusnp, &gid);
}

/*----------------------------------------------------------------------*
 |  Get the coupled the values on the coupling interface    ismail 07/10|
 |  of the 3D/reduced-D problem                                         |
 *----------------------------------------------------------------------*/
template <Core::FE::CellType distype>
void Discret::ELEMENTS::AcinusImpl<distype>::get_coupled_values(RedAcinus* ele,
    Teuchos::ParameterList& params, Core::FE::Discretization& discretization, std::vector<int>& lm,
    Teuchos::RCP<Core::Mat::Material> material)
{
  const int myrank = discretization.get_comm().MyPID();

  // The number of nodes
  const int numnode = lm.size();

  Teuchos::RCP<const Epetra_Vector> pnp = discretization.get_state("pnp");

  if (pnp == Teuchos::null) FOUR_C_THROW("Cannot get state vectors 'pnp'");

  // extract local values from the global vectors
  std::vector<double> mypnp(lm.size());
  Core::FE::ExtractMyValues(*pnp, mypnp, lm);

  // create objects for element arrays
  Core::LinAlg::SerialDenseVector epnp(numnode);

  // get all values at the last computed time step
  for (int i = 0; i < numnode; ++i)
  {
    // split area and volumetric flow rate, insert into element arrays
    epnp(i) = mypnp[i];
  }

  // ---------------------------------------------------------------------------------
  // Resolve the BCs
  // ---------------------------------------------------------------------------------
  for (int i = 0; i < ele->num_node(); i++)
  {
    if (ele->nodes()[i]->owner() == myrank)
    {
      if (ele->nodes()[i]->get_condition("Art_redD_3D_CouplingCond"))
      {
        const Core::Conditions::Condition* condition =
            ele->nodes()[i]->get_condition("Art_redD_3D_CouplingCond");
        Teuchos::RCP<Teuchos::ParameterList> CoupledTo3DParams =
            params.get<Teuchos::RCP<Teuchos::ParameterList>>("coupling with 3D fluid params");
        // -----------------------------------------------------------------
        // If the parameter list is empty, then something is wrong!
        // -----------------------------------------------------------------
        if (CoupledTo3DParams.get() == nullptr)
        {
          FOUR_C_THROW(
              "Cannot prescribe a boundary condition from 3D to reduced D, if the parameters "
              "passed don't exist");
          exit(1);
        }


        // -----------------------------------------------------------------
        // Compute the variable solved by the reduced D simulation to be
        // passed to the 3D simulation
        //
        //     In this case a map called map1D has the following form:
        //     +-----------------------------------------------------------+
        //     |              std::map< string            ,  double        > >  |
        //     |     +------------------------------------------------+    |
        //     |     |  ID  | coupling variable name | variable value |    |
        //     |     +------------------------------------------------+    |
        //     |     |  1   |   flow1                |     xxxxxxx    |    |
        //     |     +------+------------------------+----------------+    |
        //     |     |  2   |   pressure2            |     xxxxxxx    |    |
        //     |     +------+------------------------+----------------+    |
        //     |     .  .   .   ....                 .     .......    .    |
        //     |     +------+------------------------+----------------+    |
        //     |     |  N   |   variable(N)          | trash value(N) |    |
        //     |     +------+------------------------+----------------+    |
        //     +-----------------------------------------------------------+
        // -----------------------------------------------------------------

        int ID = condition->parameters().get<int>("ConditionID");
        Teuchos::RCP<std::map<std::string, double>> map1D;
        map1D = CoupledTo3DParams->get<Teuchos::RCP<std::map<std::string, double>>>(
            "reducedD map of values");

        std::string returnedBC = (condition->parameters().get<std::string>("ReturnedVariable"));

        double BC3d = 0.0;
        if (returnedBC == "flow")
        {
          // MUST BE DONE
        }
        else if (returnedBC == "pressure")
        {
          BC3d = epnp(i);
        }
        else
        {
          std::string str = (condition->parameters().get<std::string>("ReturnedVariable"));
          FOUR_C_THROW("%s, is an unimplimented type of coupling", str.c_str());
          exit(1);
        }
        std::stringstream returnedBCwithId;
        returnedBCwithId << returnedBC << "_" << ID;

        //        std::cout<<"Return ["<<returnedBC<<"] form 1D problem to 3D SURFACE of
        //        ID["<<ID<<"]: "<<BC3d<<std::endl;

        // -----------------------------------------------------------------
        // Check whether the coupling wrapper has already initialized this
        // map else wise we will have problems with parallelization, that's
        // because of the preassumption that the map is filled and sorted
        // Thus we can use parallel addition
        // -----------------------------------------------------------------

        std::map<std::string, double>::iterator itrMap1D;
        itrMap1D = map1D->find(returnedBCwithId.str());
        if (itrMap1D == map1D->end())
        {
          FOUR_C_THROW("The 3D map for (1D - 3D coupling) has no variable (%s) for ID [%d]",
              returnedBC.c_str(), ID);
          exit(1);
        }

        // update the 1D map
        (*map1D)[returnedBCwithId.str()] = BC3d;
      }
    }  // END of if node is available on this processor
  }    // End of node i has a condition
}


/*----------------------------------------------------------------------*
 |  calculate the amount of fluid mixing inside a           ismail 02/13|
 |  junction                                                            |
 *----------------------------------------------------------------------*/
template <Core::FE::CellType distype>
void Discret::ELEMENTS::AcinusImpl<distype>::get_junction_volume_mix(RedAcinus* ele,
    Teuchos::ParameterList& params, Core::FE::Discretization& discretization,
    Core::LinAlg::SerialDenseVector& volumeMix_np, std::vector<int>& lm,
    Teuchos::RCP<Core::Mat::Material> material)
{
  Discret::ReducedLung::EvaluationData& evaluation_data =
      Discret::ReducedLung::EvaluationData::get();
  const auto acinus_params = ele->get_acinus_params();

  // get the element qout
  double q_out = (*evaluation_data.qout_np)[ele->lid()];

  // if transport is flowing into the acinus
  if (q_out >= 0.0)
  {
    volumeMix_np(1) = acinus_params.area;
  }
  // els if transport is flowing out of the acinus
  else
  {
    volumeMix_np(0) = acinus_params.area;
    volumeMix_np(1) = acinus_params.area;
  }

  // extra treatment if an acinus is not connected to anything else
  for (int i = 0; i < iel; i++)
  {
    if (ele->nodes()[i]->num_element() == 1)
    {
      volumeMix_np(1) = acinus_params.area;
    }
  }
}


/*----------------------------------------------------------------------*
 |  calculate the scalar transport                          ismail 02/13|
 |                                                                      |
 *----------------------------------------------------------------------*/
template <Core::FE::CellType distype>
void Discret::ELEMENTS::AcinusImpl<distype>::update_scatra(RedAcinus* ele,
    Teuchos::ParameterList& params, Core::FE::Discretization& discretization, std::vector<int>& lm,
    Teuchos::RCP<Core::Mat::Material> material)
{
  const int myrank = discretization.get_comm().MyPID();

  Teuchos::RCP<const Epetra_Vector> dscatranp = discretization.get_state("dscatranp");
  Discret::ReducedLung::EvaluationData& evaluation_data =
      Discret::ReducedLung::EvaluationData::get();

  // get flowrate
  double qin = (*evaluation_data.qin_np)[ele->lid()];

  // extract local values from the global vectors
  std::vector<double> mydscatra(lm.size());
  Core::FE::ExtractMyValues(*dscatranp, mydscatra, lm);

  //--------------------------------------------------------------------
  // if vel>=0 then node(2) is analytically evaluated;
  //  ---> node(1) is either prescribed or comes from the junction
  // if vel< 0 then node(1) is analytically evaluated;
  //  ---> node(2) is either prescribed or comes from the junction
  //--------------------------------------------------------------------
  if (qin < 0.0)
  {
    int gid = lm[1];
    double val = mydscatra[1];
    if (myrank == ele->nodes()[1]->owner())
    {
      evaluation_data.dscatranp->ReplaceGlobalValues(1, &val, &gid);
    }
  }
}  // update_scatra


template <Core::FE::CellType distype>
void Discret::ELEMENTS::AcinusImpl<distype>::update_elem12_scatra(RedAcinus* ele,
    Teuchos::ParameterList& params, Core::FE::Discretization& discretization, std::vector<int>& lm,
    Teuchos::RCP<Core::Mat::Material> material)
{
  Teuchos::RCP<const Epetra_Vector> scatranp = discretization.get_state("scatranp");
  Teuchos::RCP<const Epetra_Vector> dscatranp = discretization.get_state("dscatranp");
  Teuchos::RCP<const Epetra_Vector> volumeMix = discretization.get_state("junctionVolumeInMix");

  Discret::ReducedLung::EvaluationData& evaluation_data =
      Discret::ReducedLung::EvaluationData::get();

  // extract local values from the global vectors
  std::vector<double> myscatranp(lm.size());
  Core::FE::ExtractMyValues(*scatranp, myscatranp, lm);

  // extract local values from the global vectors
  std::vector<double> mydscatranp(lm.size());
  Core::FE::ExtractMyValues(*dscatranp, mydscatranp, lm);

  // extract local values from the global vectors
  std::vector<double> myvolmix(lm.size());
  Core::FE::ExtractMyValues(*volumeMix, myvolmix, lm);

  // get flowrate
  double qin = (*evaluation_data.qin_np)[ele->lid()];
  // Get the average concentration

  // ---------------------------------------------------------------------
  // element scatra must be updated only at the capillary nodes.
  // ---------------------------------------------------------------------
  //  double e2s = (*e2scatranp)[ele->LID()] + mydscatranp[1]*myvolmix[1];
  double e2s = myscatranp[1];

  int gid = ele->id();
  evaluation_data.e2scatranp->ReplaceGlobalValues(1, &e2s, &gid);
  if (qin < 0.0)
  {
    evaluation_data.e1scatranp->ReplaceGlobalValues(1, &e2s, &gid);
  }
}


/*----------------------------------------------------------------------*
 |  calculate essential nodal values                        ismail 06/13|
 |                                                                      |
 *----------------------------------------------------------------------*/
template <Core::FE::CellType distype>
void Discret::ELEMENTS::AcinusImpl<distype>::eval_nodal_essential_values(RedAcinus* ele,
    Teuchos::ParameterList& params, Core::FE::Discretization& discretization,
    Core::LinAlg::SerialDenseVector& nodal_surface, Core::LinAlg::SerialDenseVector& nodal_volume,
    Core::LinAlg::SerialDenseVector& nodal_avg_scatra, std::vector<int>& lm,
    Teuchos::RCP<Core::Mat::Material> material)
{
  // Get all general state vectors: flow, pressure,
  Discret::ReducedLung::EvaluationData& evaluation_data =
      Discret::ReducedLung::EvaluationData::get();
  const auto acinus_params = ele->get_acinus_params();

  Teuchos::RCP<const Epetra_Vector> scatranp = discretization.get_state("scatranp");

  // Extract scatra values
  // Extract local values from the global vectors
  std::vector<double> myscatranp(lm.size());
  Core::FE::ExtractMyValues(*scatranp, myscatranp, lm);

  // Find the volume of an acinus
  // Get the current acinar volume
  double volAcinus = (*evaluation_data.acinar_v)[ele->lid()];
  // Set nodal volume
  nodal_volume[1] = volAcinus;

  // Find the average scalar transport concentration
  // Set nodal flowrate
  nodal_avg_scatra[0] = myscatranp[1];
  nodal_avg_scatra[1] = myscatranp[1];

  // Find the total gas exchange surface inside an acinus
  // get the initial volume of an acinus
  const double volAcinus0 = acinus_params.volume_relaxed;
  // get the initial volume of an alveolar duct
  const double volAlvDuct0 = acinus_params.alveolar_duct_volume;
  // find the  number of alveolar duct
  const double numOfAlvDucts = double(floor(volAcinus0 / volAlvDuct0));
  // define the number of alveoli per alveolar duct
  const double nAlveoliPerAlveolarDuct = 36.0;
  // define the number of alveoli per duct
  const double nAlveoliPerDuct = 4.0;
  // find the volume of one alveolus
  const double volAlveolus = volAcinus / numOfAlvDucts / nAlveoliPerAlveolarDuct;
  // find the surface of one alveolus
  const double surfAlveolus =
      (6.0 + 12.0 * sqrt(3.0)) * pow(volAlveolus / (8.0 * sqrt(2.0)), 2.0 / 3.0);
  // get the length of an edge of an alveolus
  const double al = pow(volAlveolus / (8.0 * sqrt(2.0)), 1.0 / 3.0) / 3.0;
  // find the surface of an alveolar duct
  const double surfAlveolarDuct =
      (nAlveoliPerAlveolarDuct - 2.0 * nAlveoliPerDuct) * surfAlveolus + 6.0 * (al * al);
  // find the surface of an acinus
  const double surfAcinus = surfAlveolarDuct * numOfAlvDucts;
  // set nodal surface area
  nodal_surface[1] = surfAcinus;
}

FOUR_C_NAMESPACE_CLOSE
