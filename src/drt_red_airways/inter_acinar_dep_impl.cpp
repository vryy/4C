/*----------------------------------------------------------------------*/
/*!
\file inter_acinar_dep_impl.cpp

\brief Internal implementation of RedInterAcinarDep element. Routines
       from inter_acinar_dep_evaluate are re-implemented here.

\maintainer Lena Yoshihara

\level 3
*/
/*----------------------------------------------------------------------*/



#include "acinus_impl.H"
#include "inter_acinar_dep_impl.H"

#include "../drt_mat/newtonianfluid.H"
#include "../drt_mat/maxwell_0d_acinus.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_function.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/standardtypes_cpp.H"
#include "../drt_fem_general/drt_utils_fem_shapefunctions.H"
#include <fstream>
#include <iomanip>


/*----------------------------------------------------------------------*
 |                                                         ismail 01/10 |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::RedInterAcinarDepImplInterface* DRT::ELEMENTS::RedInterAcinarDepImplInterface::Impl(
    DRT::ELEMENTS::RedInterAcinarDep* red_acinus)
{
  switch (red_acinus->Shape())
  {
    case DRT::Element::line2:
    {
      static InterAcinarDepImpl<DRT::Element::line2>* acinus;
      if (acinus == NULL)
      {
        acinus = new InterAcinarDepImpl<DRT::Element::line2>;
      }
      return acinus;
    }
    default:
      dserror("shape %d (%d nodes) not supported", red_acinus->Shape(), red_acinus->NumNode());
      break;
  }
  return NULL;
}


/*----------------------------------------------------------------------*
 | Constructor (public)                                    ismail 01/10 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::InterAcinarDepImpl<distype>::InterAcinarDepImpl()
{
}


/*----------------------------------------------------------------------*
 | Evaluate (public)                                       ismail 01/10 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
int DRT::ELEMENTS::InterAcinarDepImpl<distype>::Evaluate(RedInterAcinarDep* ele,
    Teuchos::ParameterList& params, DRT::Discretization& discretization, std::vector<int>& lm,
    Epetra_SerialDenseMatrix& elemat1_epetra, Epetra_SerialDenseMatrix& elemat2_epetra,
    Epetra_SerialDenseVector& elevec1_epetra, Epetra_SerialDenseVector& elevec2_epetra,
    Epetra_SerialDenseVector& elevec3_epetra, Teuchos::RCP<MAT::Material> mat)
{
  // Get the vector with inter-acinar linkers
  Teuchos::RCP<const Epetra_Vector> ial = discretization.GetState("intr_ac_link");

  // Extract local values from the global vectors
  std::vector<double> myial(lm.size());
  DRT::UTILS::ExtractMyValues(*ial, myial, lm);

  // Calculate the system matrix for inter-acinar linkers
  Sysmat(myial, elemat1_epetra, elevec1_epetra);

  return 0;
}


/*----------------------------------------------------------------------*
 | Initial routine, sets generation number for inter-acinar linker      |
 | element to -2.0 and sets the number of linkers per node in this      |
 | element to 1.0. The final sum of linkers for each node is auto-      |
 | matically evaluated during the assembly process later.               |
 |                                              (private)  ismail 01/10 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::InterAcinarDepImpl<distype>::Initial(RedInterAcinarDep* ele,
    Teuchos::ParameterList& params, DRT::Discretization& discretization, std::vector<int>& lm,
    Epetra_SerialDenseVector& n_intr_acn_l, Teuchos::RCP<const MAT::Material> material)
{
  Teuchos::RCP<Epetra_Vector> generations = params.get<Teuchos::RCP<Epetra_Vector>>("generations");

  // Set the generation number for the inter-acinar linker element to -2.0
  int gid = ele->Id();
  double val = -2.0;
  generations->ReplaceGlobalValues(1, &val, &gid);

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
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::InterAcinarDepImpl<distype>::Sysmat(
    std::vector<double>& ial, Epetra_SerialDenseMatrix& sysmat, Epetra_SerialDenseVector& rhs)
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
  rhs.Scale(0.0);
}


/*----------------------------------------------------------------------*
 |  Evaluate the values of the degrees of freedom           ismail 04/13|
 |  at terminal nodes.                                                  |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::InterAcinarDepImpl<distype>::EvaluateTerminalBC(RedInterAcinarDep* ele,
    Teuchos::ParameterList& params, DRT::Discretization& discretization, std::vector<int>& lm,
    Epetra_SerialDenseVector& rhs, Teuchos::RCP<MAT::Material> material)
{
  const int myrank = discretization.Comm().MyPID();

  // Get total time
  const double time = params.get<double>("total time");

  // Get the number of nodes
  const int numnode = lm.size();
  std::vector<int>::iterator it_vcr;

  // Get state for pressure
  Teuchos::RCP<const Epetra_Vector> pnp = discretization.GetState("pnp");
  if (pnp == Teuchos::null) dserror("Cannot get state vectors 'pnp'");

  // Extract local values from the global vectors
  std::vector<double> mypnp(lm.size());
  DRT::UTILS::ExtractMyValues(*pnp, mypnp, lm);

  // Create objects for element arrays
  Epetra_SerialDenseVector epnp(numnode);

  // Get all values at the last computed time step
  for (int i = 0; i < numnode; ++i)
  {
    // Split area and volumetric flow rate, insert into element arrays
    epnp(i) = mypnp[i];
  }

  /**
   * Resolve the BCs
   **/
  for (int i = 0; i < ele->NumNode(); i++)
  {
    if (ele->Nodes()[i]->Owner() == myrank)
    {
      if (ele->Nodes()[i]->GetCondition("RedAirwayPrescribedCond"))
      {
        std::string Bc;
        double BCin = 0.0;
        if (ele->Nodes()[i]->GetCondition("RedAirwayPrescribedCond"))
        {
          DRT::Condition* condition = ele->Nodes()[i]->GetCondition("RedAirwayPrescribedCond");
          // Get the type of prescribed bc
          Bc = *(condition->Get<std::string>("boundarycond"));

          const std::vector<int>* curve = condition->Get<std::vector<int>>("curve");
          double curvefac = 1.0;
          const std::vector<double>* vals = condition->Get<std::vector<double>>("val");
          const std::vector<int>* functions = condition->Get<std::vector<int>>("funct");

          // Read in the value of the applied BC
          // Get factor of first CURVE
          if ((*curve)[0] >= 0)
          {
            curvefac = DRT::Problem::Instance()->Funct((*curve)[0]).EvaluateTime(time);
            BCin = (*vals)[0] * curvefac;
          }
          else
          {
            dserror("no boundary condition defined!");
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
            functionfac = DRT::Problem::Instance()
                              ->Funct(functnum - 1)
                              .Evaluate(0, (ele->Nodes()[i])->X(), time);
          }

          // Get factor of second CURVE
          int curve2num = -1;
          double curve2fac = 1.0;
          if (curve) curve2num = (*curve)[1];
          if (curve2num >= 0)
            curve2fac = DRT::Problem::Instance()->Funct(curve2num).EvaluateTime(time);

          // Add first_CURVE + FUNCTION * second_CURVE
          BCin += functionfac * curve2fac;

          // Get the local id of the node to whom the bc is prescribed
          int local_id = discretization.NodeRowMap()->LID(ele->Nodes()[i]->Id());
          if (local_id < 0)
          {
            dserror("node (%d) doesn't exist on proc(%d)", ele->Nodes()[i]->Id(),
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
          Teuchos::RCP<Epetra_Vector> bcval = params.get<Teuchos::RCP<Epetra_Vector>>("bcval");
          Teuchos::RCP<Epetra_Vector> dbctog = params.get<Teuchos::RCP<Epetra_Vector>>("dbctog");

          if (bcval == Teuchos::null || dbctog == Teuchos::null)
          {
            dserror("Cannot get state vectors 'bcval' and 'dbctog'");
            exit(1);
          }

          if (Bc == "VolumeDependentPleuralPressure")
          {
            DRT::Condition* pplCond =
                ele->Nodes()[i]->GetCondition("RedAirwayVolDependentPleuralPressureCond");
            double Pp_np = 0.0;
            if (pplCond)
            {
              const std::vector<int>* curve = pplCond->Get<std::vector<int>>("curve");
              double curvefac = 1.0;
              const std::vector<double>* vals = pplCond->Get<std::vector<double>>("val");

              // Read in the value of the applied BC
              if ((*curve)[0] >= 0)
              {
                curvefac = DRT::Problem::Instance()->Funct((*curve)[0]).EvaluateTime(time);
              }

              // Get parameters for VolumeDependentPleuralPressure condition
              std::string ppl_Type = *(pplCond->Get<std::string>("TYPE"));
              double ap = pplCond->GetDouble("P_PLEURAL_0");
              double bp = pplCond->GetDouble("P_PLEURAL_LIN");
              double cp = pplCond->GetDouble("P_PLEURAL_NONLIN");
              double dp = pplCond->GetDouble("TAU");
              double RV = pplCond->GetDouble("RV");
              double TLC = pplCond->GetDouble("TLC");

              // Safety check: in case of polynomial TLC is not used
              if (((ppl_Type == "Linear_Polynomial") or (ppl_Type == "Nonlinear_Polynomial")) and
                  (TLC != 0.0))
              {
                dserror(
                    "TLC is not used for the following type of VolumeDependentPleuralPressure BC: "
                    "%s.\n Set TLC = 0.0",
                    ppl_Type.c_str());
              }
              // Safety check: in case of Ogden TLC, P_PLEURAL_0, and P_PLEURAL_LIN
              if ((ppl_Type == "Nonlinear_Ogden") and
                  ((TLC != 0.0) or (ap != 0.0) or (bp != 0.0) or (dp == 0.0)))
              {
                dserror(
                    "Parameters are not set correctly for Nonlinear_Ogden. Only P_PLEURAL_NONLIN, "
                    "TAU and RV are used. Set all others to zero. TAU is not allowed to be zero.");
              }

              if (ppl_Type == "Linear_Polynomial")
              {
                const double lungVolumenp = params.get<double>("lungVolume_n");
                Pp_np = ap + bp * (lungVolumenp - RV) + cp * pow((lungVolumenp - RV), dp);
              }
              else if (ppl_Type == "Linear_Exponential")
              {
                const double lungVolumenp = params.get<double>("lungVolume_n");
                const double TLCnp = (lungVolumenp - RV) / (TLC - RV);
                Pp_np = ap + bp * TLCnp + cp * exp(dp * TLCnp);
              }
              else if (ppl_Type == "Linear_Ogden")
              {
                const double lungVolumenp = params.get<double>("lungVolume_n");
                Pp_np = RV / lungVolumenp * cp / dp * (1 - pow(RV / lungVolumenp, dp));
              }
              else if (ppl_Type == "Nonlinear_Polynomial")
              {
                const double lungVolumenp = params.get<double>("lungVolume_np");
                Pp_np = ap + bp * (lungVolumenp - RV) + cp * pow((lungVolumenp - RV), dp);
              }
              else if (ppl_Type == "Nonlinear_Exponential")
              {
                const double lungVolumenp = params.get<double>("lungVolume_np");
                const double TLCnp = (lungVolumenp - RV) / (TLC - RV);
                Pp_np = ap + bp * TLCnp + cp * exp(dp * TLCnp);
              }
              else if (ppl_Type == "Nonlinear_Ogden")
              {
                const double lungVolumenp = params.get<double>("lungVolume_np");
                Pp_np = RV / lungVolumenp * cp / dp * (1 - pow(RV / lungVolumenp, dp));
              }
              else
              {
                dserror("Unknown volume pleural pressure type: %s", ppl_Type.c_str());
              }
              Pp_np *= curvefac * ((*vals)[0]);
            }
            else
            {
              std::cout << "Node " << ele->Nodes()[i]->Id() + 1 << "is not on corresponding DLINE "
                        << std::endl;
              dserror("No volume dependent pleural pressure condition was defined");
            }

            BCin += Pp_np;
          }

          // Set pressure at node i
          int gid;
          double val;

          gid = lm[i];
          val = BCin;
          bcval->ReplaceGlobalValues(1, &val, &gid);

          gid = lm[i];
          val = 1;
          dbctog->ReplaceGlobalValues(1, &val, &gid);
        }
        else
        {
          dserror("Prescribed [%s] is not defined for reduced-inter-acinar linkers", Bc.c_str());
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
            dserror("node (%d) doesn't exist on proc(%d)", ele->Nodes()[i],
                discretization.Comm().MyPID());
            exit(1);
          }

          Teuchos::RCP<Epetra_Vector> bcval = params.get<Teuchos::RCP<Epetra_Vector>>("bcval");
          Teuchos::RCP<Epetra_Vector> dbctog = params.get<Teuchos::RCP<Epetra_Vector>>("dbctog");

          if (bcval == Teuchos::null || dbctog == Teuchos::null)
          {
            dserror("Cannot get state vectors 'bcval' and 'dbctog'");
            exit(1);
          }

          // Set pressure at node i
          int gid;
          double val;

          gid = lm[i];
          val = 0.0;
          bcval->ReplaceGlobalValues(1, &val, &gid);

          gid = lm[i];
          val = 1;
          dbctog->ReplaceGlobalValues(1, &val, &gid);
        }
      }  // END of if there is no BC but the node still is at the terminal
    }    // END of if node is available on this processor
  }      // End of node i has a condition
}


/*----------------------------------------------------------------------*
 |  CalcFlowRates                                                       |
 |                                                         ismail 01/10 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::InterAcinarDepImpl<distype>::CalcFlowRates(RedInterAcinarDep* ele,
    Teuchos::ParameterList& params, DRT::Discretization& discretization,
    Epetra_SerialDenseVector& elevec1,  // a_volumenp,
    Epetra_SerialDenseVector& elevec2,  // a_volume_strain_np,
    std::vector<int>& lm, Teuchos::RCP<MAT::Material> material)

{
}


/*----------------------------------------------------------------------*
 |  Get the coupled the values on the coupling interface    ismail 07/10|
 |  of the 3D/reduced-D problem                                         |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::InterAcinarDepImpl<distype>::GetCoupledValues(RedInterAcinarDep* ele,
    Teuchos::ParameterList& params, DRT::Discretization& discretization, std::vector<int>& lm,
    Teuchos::RCP<MAT::Material> material)
{
}
