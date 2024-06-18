/*---------------------------------------------------------------------*/
/*! \file

\brief Internal implementation of RedInterAcinarDep element. Methods implemented here
       are called by inter_acinar_dep_evaluate.cpp by
Discret::ELEMENTS::RedInterAcinarDep::evaluate() with the corresponding action.


\level 3

*/
/*---------------------------------------------------------------------*/



#include "4C_red_airways_interacinardep_impl.hpp"

#include "4C_fem_discretization.hpp"
#include "4C_fem_general_extract_values.hpp"
#include "4C_fem_general_utils_fem_shapefunctions.hpp"
#include "4C_global_data.hpp"
#include "4C_mat_newtonianfluid.hpp"
#include "4C_red_airways_evaluation_data.hpp"
#include "4C_utils_function.hpp"
#include "4C_utils_function_of_time.hpp"

#include <fstream>

FOUR_C_NAMESPACE_OPEN


/*----------------------------------------------------------------------*
 |                                                         ismail 01/10 |
 *----------------------------------------------------------------------*/
Discret::ELEMENTS::RedInterAcinarDepImplInterface*
Discret::ELEMENTS::RedInterAcinarDepImplInterface::Impl(
    Discret::ELEMENTS::RedInterAcinarDep* red_acinus)
{
  switch (red_acinus->Shape())
  {
    case Core::FE::CellType::line2:
    {
      static InterAcinarDepImpl<Core::FE::CellType::line2>* acinus;
      if (acinus == nullptr)
      {
        acinus = new InterAcinarDepImpl<Core::FE::CellType::line2>;
      }
      return acinus;
    }
    default:
      FOUR_C_THROW(
          "shape %d (%d nodes) not supported", red_acinus->Shape(), red_acinus->num_node());
      break;
  }
  return nullptr;
}


/*----------------------------------------------------------------------*
 | Constructor (public)                                    ismail 01/10 |
 *----------------------------------------------------------------------*/
template <Core::FE::CellType distype>
Discret::ELEMENTS::InterAcinarDepImpl<distype>::InterAcinarDepImpl()
{
}


/*----------------------------------------------------------------------*
 | Evaluate (public)                                       ismail 01/10 |
 *----------------------------------------------------------------------*/
template <Core::FE::CellType distype>
int Discret::ELEMENTS::InterAcinarDepImpl<distype>::evaluate(RedInterAcinarDep* ele,
    Teuchos::ParameterList& params, Core::FE::Discretization& discretization, std::vector<int>& lm,
    Core::LinAlg::SerialDenseMatrix& elemat1_epetra,
    Core::LinAlg::SerialDenseMatrix& elemat2_epetra,
    Core::LinAlg::SerialDenseVector& elevec1_epetra,
    Core::LinAlg::SerialDenseVector& elevec2_epetra,
    Core::LinAlg::SerialDenseVector& elevec3_epetra, Teuchos::RCP<Core::Mat::Material> mat)
{
  // Get the vector with inter-acinar linkers
  Teuchos::RCP<const Epetra_Vector> ial = discretization.GetState("intr_ac_link");

  // Extract local values from the global vectors
  std::vector<double> myial(lm.size());
  Core::FE::ExtractMyValues(*ial, myial, lm);

  // Calculate the system matrix for inter-acinar linkers
  sysmat(myial, elemat1_epetra, elevec1_epetra);

  return 0;
}


/*----------------------------------------------------------------------*
 | Initial routine, sets generation number for inter-acinar linker      |
 | element to -2.0 and sets the number of linkers per node in this      |
 | element to 1.0. The final sum of linkers for each node is auto-      |
 | matically evaluated during the assembly process later.               |
 |                                              (private)  ismail 01/10 |
 *----------------------------------------------------------------------*/
template <Core::FE::CellType distype>
void Discret::ELEMENTS::InterAcinarDepImpl<distype>::Initial(RedInterAcinarDep* ele,
    Teuchos::ParameterList& params, Core::FE::Discretization& discretization, std::vector<int>& lm,
    Core::LinAlg::SerialDenseVector& n_intr_acn_l, Teuchos::RCP<const Core::Mat::Material> material)
{
  Discret::ReducedLung::EvaluationData& evaluation_data =
      Discret::ReducedLung::EvaluationData::get();

  // Set the generation number for the inter-acinar linker element to -2.0
  int gid = ele->Id();
  double val = -2.0;
  evaluation_data.generations->ReplaceGlobalValues(1, &val, &gid);

  // In this element, each node of an inter-acinar linker element has
  // one linker. The final sum of linkers for each node is automatically
  // evaluated during the assembly process.
  n_intr_acn_l(0) = 1.0;
  n_intr_acn_l(1) = 1.0;

}  // InterAcinarDepImpl::Initial


/*----------------------------------------------------------------------*
 | Calculate element matrix and right hand side (private). The system   |
 | matrix of an inter-acinar linker element is +/-1/(number of linkers  |
 | per node). The right hand side is zero.                              |
 |                                                         ismail 01/10 |
 *----------------------------------------------------------------------*/
template <Core::FE::CellType distype>
void Discret::ELEMENTS::InterAcinarDepImpl<distype>::sysmat(std::vector<double>& ial,
    Core::LinAlg::SerialDenseMatrix& sysmat, Core::LinAlg::SerialDenseVector& rhs)
{
  // Get the number of inter_acinar linkers on the 1st node (N0)
  double N0 = ial[0];
  // Get the number of inter_acinar linkers on the 2nd node (N1)
  double N1 = ial[1];
  if (N0 > 0)
  {
    sysmat(0, 0) = 1.0 / (N0);
    sysmat(0, 1) = -1.0 / (N0);
  }
  if (N1 > 0)
  {
    sysmat(1, 0) = -1.0 / (N1);
    sysmat(1, 1) = 1.0 / (N1);
  }
  rhs.putScalar(0.0);
}


/*----------------------------------------------------------------------*
 |  Evaluate the values of the degrees of freedom           ismail 04/13|
 |  at terminal nodes.                                                  |
 *----------------------------------------------------------------------*/
template <Core::FE::CellType distype>
void Discret::ELEMENTS::InterAcinarDepImpl<distype>::EvaluateTerminalBC(RedInterAcinarDep* ele,
    Teuchos::ParameterList& params, Core::FE::Discretization& discretization, std::vector<int>& lm,
    Core::LinAlg::SerialDenseVector& rhs, Teuchos::RCP<Core::Mat::Material> material)
{
  const int myrank = discretization.Comm().MyPID();

  Discret::ReducedLung::EvaluationData& evaluation_data =
      Discret::ReducedLung::EvaluationData::get();

  // Get total time
  const double time = evaluation_data.time;

  // Get the number of nodes
  const int numnode = lm.size();

  // Get state for pressure
  Teuchos::RCP<const Epetra_Vector> pnp = discretization.GetState("pnp");
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
    if (ele->Nodes()[i]->Owner() == myrank)
    {
      if (ele->Nodes()[i]->GetCondition("RedAirwayPrescribedCond"))
      {
        std::string Bc;
        double BCin = 0.0;
        if (ele->Nodes()[i]->GetCondition("RedAirwayPrescribedCond"))
        {
          Core::Conditions::Condition* condition =
              ele->Nodes()[i]->GetCondition("RedAirwayPrescribedCond");
          // Get the type of prescribed bc
          Bc = (condition->parameters().Get<std::string>("boundarycond"));

          const auto* curve = &condition->parameters().Get<std::vector<int>>("curve");
          double curvefac = 1.0;
          const auto* vals = &condition->parameters().Get<std::vector<double>>("val");
          const auto* functions = &condition->parameters().Get<std::vector<int>>("funct");

          // Read in the value of the applied BC
          // Get factor of first CURVE
          if ((*curve)[0] >= 0)
          {
            curvefac = Global::Problem::Instance()
                           ->FunctionById<Core::UTILS::FunctionOfTime>((*curve)[0])
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
            functionfac = Global::Problem::Instance()
                              ->FunctionById<Core::UTILS::FunctionOfSpaceTime>(functnum - 1)
                              .evaluate((ele->Nodes()[i])->X().data(), time, 0);
          }

          // Get factor of second CURVE
          int curve2num = -1;
          double curve2fac = 1.0;
          if (curve) curve2num = (*curve)[1];
          if (curve2num >= 0)
            curve2fac = Global::Problem::Instance()
                            ->FunctionById<Core::UTILS::FunctionOfTime>(curve2num)
                            .evaluate(time);

          // Add first_CURVE + FUNCTION * second_CURVE
          BCin += functionfac * curve2fac;

          // Get the local id of the node to whom the bc is prescribed
          int local_id = discretization.NodeRowMap()->LID(ele->Nodes()[i]->Id());
          if (local_id < 0)
          {
            FOUR_C_THROW("node (%d) doesn't exist on proc(%d)", ele->Nodes()[i]->Id(),
                discretization.Comm().MyPID());
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
                ele->Nodes()[i]->GetCondition("RedAirwayVolDependentPleuralPressureCond");
            double Pp_np = 0.0;
            if (pplCond)
            {
              const auto* curve = &pplCond->parameters().Get<std::vector<int>>("curve");
              double curvefac = 1.0;
              const auto* vals = &pplCond->parameters().Get<std::vector<double>>("val");

              // Read in the value of the applied BC
              if ((*curve)[0] >= 0)
              {
                curvefac = Global::Problem::Instance()
                               ->FunctionById<Core::UTILS::FunctionOfTime>((*curve)[0])
                               .evaluate(time);
              }

              // Get parameters for VolumeDependentPleuralPressure condition
              std::string ppl_Type = (pplCond->parameters().Get<std::string>("TYPE"));
              auto ap = pplCond->parameters().Get<double>("P_PLEURAL_0");
              auto bp = pplCond->parameters().Get<double>("P_PLEURAL_LIN");
              auto cp = pplCond->parameters().Get<double>("P_PLEURAL_NONLIN");
              auto dp = pplCond->parameters().Get<double>("TAU");
              auto RV = pplCond->parameters().Get<double>("RV");
              auto TLC = pplCond->parameters().Get<double>("TLC");

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

              Discret::ReducedLung::EvaluationData& evaluation_data =
                  Discret::ReducedLung::EvaluationData::get();

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
              std::cout << "Node " << ele->Nodes()[i]->Id() + 1 << "is not on corresponding DLINE "
                        << std::endl;
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
        else
        {
          FOUR_C_THROW(
              "Prescribed [%s] is not defined for reduced-inter-acinar linkers", Bc.c_str());
          exit(1);
        }
      }
      /**
       * If the node is a terminal node, but no b.c is prescribed to it
       * then a zero output pressure is assumed
       **/
      else
      {
        if (ele->Nodes()[i]->NumElement() == 1)
        {
          // Get the local id of the node to whom the bc is prescribed
          int local_id = discretization.NodeRowMap()->LID(ele->Nodes()[i]->Id());
          if (local_id < 0)
          {
            FOUR_C_THROW("node (%d) doesn't exist on proc(%d)", ele->Nodes()[i],
                discretization.Comm().MyPID());
            exit(1);
          }

          Discret::ReducedLung::EvaluationData& evaluation_data =
              Discret::ReducedLung::EvaluationData::get();

          // Set pressure at node i
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
    }    // END of if node is available on this processor
  }      // End of node i has a condition
}

FOUR_C_NAMESPACE_CLOSE
