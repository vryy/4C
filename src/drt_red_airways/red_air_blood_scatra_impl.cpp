/*---------------------------------------------------------------------*/
/*!

\brief Incomplete! - Purpose: Internal implementation of RedAirBloodScatra element

\maintainer Carolin Geitner

\level 3

*/
/*---------------------------------------------------------------------*/



#include "red_air_blood_scatra_impl.H"

#include "../drt_mat/newtonianfluid.H"
#include "../drt_mat/matpar_bundle.H"
#include "../drt_mat/matlist.H"
#include "../drt_mat/air_0d_O2_saturation.H"
#include "../drt_mat/hemoglobin_0d_O2_saturation.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_function.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/standardtypes_cpp.H"
#include "../drt_fem_general/drt_utils_fem_shapefunctions.H"
#include "../drt_fem_general/drt_utils_gder2.H"
#include <fstream>
#include <iomanip>


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::RedAirBloodScatraImplInterface* DRT::ELEMENTS::RedAirBloodScatraImplInterface::Impl(
    DRT::ELEMENTS::RedAirBloodScatra* red_acinus)
{
  switch (red_acinus->Shape())
  {
    case DRT::Element::line2:
    {
      static RedAirBloodScatraImpl<DRT::Element::line2>* acinus;
      if (acinus == NULL)
      {
        acinus = new RedAirBloodScatraImpl<DRT::Element::line2>;
      }
      return acinus;
    }
    default:
      dserror("shape %d (%d nodes) not supported", red_acinus->Shape(), red_acinus->NumNode());
  }
  return NULL;
}



/*----------------------------------------------------------------------*
  | constructor (public)                                    ismail 01/10 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::RedAirBloodScatraImpl<distype>::RedAirBloodScatraImpl()
{
}

/*----------------------------------------------------------------------*
 | evaluate (public)                                       ismail 01/10 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
int DRT::ELEMENTS::RedAirBloodScatraImpl<distype>::Evaluate(RedAirBloodScatra* ele,
    Teuchos::ParameterList& params, DRT::Discretization& discretization, std::vector<int>& lm,
    Epetra_SerialDenseMatrix& elemat1_epetra, Epetra_SerialDenseMatrix& elemat2_epetra,
    Epetra_SerialDenseVector& elevec1_epetra, Epetra_SerialDenseVector& elevec2_epetra,
    Epetra_SerialDenseVector& elevec3_epetra, Teuchos::RCP<MAT::Material> mat)
{
  return 0;
}


/*----------------------------------------------------------------------*
 |  calculate element matrix and right hand side (private)  ismail 01/10|
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::RedAirBloodScatraImpl<distype>::Initial(RedAirBloodScatra* ele,
    Teuchos::ParameterList& params, DRT::Discretization& discretization, std::vector<int>& lm,
    Teuchos::RCP<const MAT::Material> material)
{
  Teuchos::RCP<Epetra_Vector> generations = params.get<Teuchos::RCP<Epetra_Vector>>("generations");

  //--------------------------------------------------------------------
  // get the generation numbers
  //--------------------------------------------------------------------
  //  if(myrank == ele->Owner())
  {
    int gid = ele->Id();
    double val = -2.0;
    generations->ReplaceGlobalValues(1, &val, &gid);
  }

}  // RedAirBloodScatraImpl::Initial

/*----------------------------------------------------------------------*
 |  calculate element matrix and right hand side (private)  ismail 01/10|
 |                                                                      |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::RedAirBloodScatraImpl<distype>::Sysmat(RedAirBloodScatra* ele,
    Epetra_SerialDenseVector& epnp, Epetra_SerialDenseVector& epn, Epetra_SerialDenseVector& epnm,
    Epetra_SerialDenseMatrix& sysmat, Epetra_SerialDenseVector& rhs,
    Teuchos::RCP<const MAT::Material> material, Teuchos::ParameterList& params, double time,
    double dt)
{
}

/*----------------------------------------------------------------------*
 |  Evaluate the values of the degrees of freedom           ismail 01/10|
 |  at terminal nodes.                                                  |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::RedAirBloodScatraImpl<distype>::EvaluateTerminalBC(RedAirBloodScatra* ele,
    Teuchos::ParameterList& params, DRT::Discretization& discretization, std::vector<int>& lm,
    Epetra_SerialDenseVector& rhs, Teuchos::RCP<MAT::Material> material)
{
}


/*----------------------------------------------------------------------*
 |  Evaluate the values of the degrees of freedom           ismail 01/10|
 |  at terminal nodes.                                                  |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::RedAirBloodScatraImpl<distype>::CalcFlowRates(RedAirBloodScatra* ele,
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
void DRT::ELEMENTS::RedAirBloodScatraImpl<distype>::GetCoupledValues(RedAirBloodScatra* ele,
    Teuchos::ParameterList& params, DRT::Discretization& discretization, std::vector<int>& lm,
    Teuchos::RCP<MAT::Material> material)
{
}

/*----------------------------------------------------------------------*
 |  solve the transport of O2 from air to blood             ismail 06/13|
 |                                                                      |
 | Example of use (1):                                                  |
 |--------------------                                                  |
 |                                                                      |
 |        [RED_AIRWAY element]          [RED_ACINUS element]            |
 |             |                             |                          |
 |             |                             |                          |
 |  (node1)    V     (node2)        (node2)  V   (node1)                |
 |     o======>>>=======o              o============o                   |
 |     |       ^                       |                                |
 |     |       |                       |                                |
 |     |(flow direction)               |                                |
 |     |                               |                                |
 |     V                               V                                |
 |     o===============================o                                |
 |  (node1)     ^                   (node2)                             |
 |    or        |                     or                                |
 |  (node2)     |                   (node1)                             |
 |              |                                                       |
 |              |                                                       |
 |    [RED_AIR_BLOOD_SCATRA element]                                    |
 |                                                                      |
 | Example of use (2):                                                  |
 |--------------------                                                  |
 |                                                                      |
 |        [RED_AIRWAY element]          [RED_ACINUS element]            |
 |             |                             |                          |
 |             |                             |                          |
 |  (node1)    V     (node2)        (node2)  V   (node1)                |
 |     o======>>>=======o              o============o                   |
 |             ^        |              |                                |
 |             |        |              |                                |
 |      (flow direction)|              |                                |
 |                      |              |                                |
 |                      V              V                                |
 |                      o==============o                                |
 |                   (node1)    ^   (node2)                             |
 |                     or       |     or                                |
 |                   (node2)    |   (node1)                             |
 |                              |                                       |
 |                              |                                       |
 |                [RED_AIR_BLOOD_SCATRA element]                        |
 |                                                                      |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::RedAirBloodScatraImpl<distype>::SolveBloodAirTransport(RedAirBloodScatra* ele,
    Epetra_SerialDenseVector& dscatra, Epetra_SerialDenseVector& dvo2,
    Epetra_SerialDenseVector& scatra_acinus, Teuchos::ParameterList& params,
    DRT::Discretization& discretization, std::vector<int>& lm, Teuchos::RCP<MAT::Material> material)

{
  // const int   myrank  = discretization.Comm().MyPID();

  // get time-step size
  const double dt = params.get<double>("time step size");


  Teuchos::RCP<const Epetra_Vector> volnp = discretization.GetState("volumenp");
  Teuchos::RCP<const Epetra_Vector> areanp = discretization.GetState("areanp");
  Teuchos::RCP<const Epetra_Vector> scatranp = discretization.GetState("scatranp");

  // extract local values from the global vectors
  std::vector<double> myscatranp(lm.size());
  DRT::UTILS::ExtractMyValues(*scatranp, myscatranp, lm);

  // extract local values from the global vectors
  std::vector<double> myvolnp(lm.size());
  DRT::UTILS::ExtractMyValues(*volnp, myvolnp, lm);

  // extract local values from the global vectors
  std::vector<double> myareanp(lm.size());
  DRT::UTILS::ExtractMyValues(*areanp, myareanp, lm);

  //--------------------------------------------------------------------
  // define the nodes connected to an acinus
  std::vector<unsigned int> ai;
  std::vector<unsigned int> ci;

  for (unsigned int i = 0; i < lm.size(); i++)
  {
    if (ele->Nodes()[i]->GetCondition("RedAirwayScatraAirCond") != NULL)
    {
      ai.push_back(i);
    }
    else if (ele->Nodes()[i]->GetCondition("RedAirwayScatraHemoglobinCond") != NULL)
    {
      ci.push_back(i);
    }
    else
    {
    }
  }

  // define an empty vO2 vector
  std::vector<double> vO2(lm.size(), 0.0);

  // -------------------------------------------------------------------
  // Convert O2 concentration to PO2
  // -------------------------------------------------------------------

  // -------------------------------------------------------------------
  // Get O2 properties in air and blood
  // -------------------------------------------------------------------

  // --------------------------------
  // Get O2 properties in air
  // --------------------------------
  int aid =
      DRT::Problem::Instance()->Materials()->FirstIdByType(INPAR::MAT::m_0d_o2_air_saturation);
  // check if O2 properties material exists
  if (aid == -1)
  {
    dserror("A material defining O2 properties in air could not be found");
    exit(1);
  }
  const MAT::PAR::Parameter* amat = DRT::Problem::Instance()->Materials()->ParameterById(aid);
  const MAT::PAR::Air_0d_O2_saturation* aactmat =
      static_cast<const MAT::PAR::Air_0d_O2_saturation*>(amat);

  // get atmospheric pressure
  const double patm = aactmat->atmospheric_p_;
  // get number of O2 moles per unit volume of O2
  const double nO2perVO2a = aactmat->nO2_per_VO2_;

  // --------------------------------
  // Get O2 properties in Blood
  // --------------------------------
  int bid = DRT::Problem::Instance()->Materials()->FirstIdByType(
      INPAR::MAT::m_0d_o2_hemoglobin_saturation);

  // check if O2 properties material exists
  if (bid == -1)
  {
    dserror("A material defining O2 properties in blood could not be found");
    exit(1);
  }
  const MAT::PAR::Parameter* bmat = DRT::Problem::Instance()->Materials()->ParameterById(bid);
  const MAT::PAR::Hemoglobin_0d_O2_saturation* bactmat =
      static_cast<const MAT::PAR::Hemoglobin_0d_O2_saturation*>(bmat);

  // how much of blood satisfies this rule
  double per_volume_blood = bactmat->per_volume_blood_;
  double o2_sat_per_vol_blood = bactmat->o2_sat_per_vol_blood_;
  double ph = bactmat->p_half_;
  double power = bactmat->power_;
  double nO2perVO2b = bactmat->nO2_per_VO2_;

  // -------------------------------------------------------------------
  // Evaluate VO2 properties in air
  // -------------------------------------------------------------------
  // get acinar volume
  double vAir = myvolnp[ai[0]];
  {
    // calculate the VO2 at nodes
    vO2[ai[0]] = (myscatranp[ai[0]] * vAir) / nO2perVO2a;
  }

  // -------------------------------------------------------------------
  // Evaluate VO2 in blood
  // -------------------------------------------------------------------
  // get capillary volume
  double vBlood = 0.0;
  for (unsigned int i = 0; i < ci.size(); i++)
  {
    vBlood += myvolnp[ci[i]];
  }
  vBlood /= float(ci.size());

  // get the ratio of blood volume to the reference saturation volume
  double alpha = vBlood / per_volume_blood;
  double kv = o2_sat_per_vol_blood * alpha;
  {
    // calculate the VO2 at nodes
    for (unsigned int i = 0; i < ci.size(); i++)
    {
      vO2[ci[i]] = (myscatranp[ci[i]] * vBlood) / nO2perVO2b;
    }
  }

  // -------------------------------------------------------------------
  // Find the deltaVO2 that flowed from air to blood
  // -------------------------------------------------------------------
  const double tol = 1e-7;
  const unsigned int maxItr = 200;
  double error = 1e7;

  // get Diffussion coefficient
  double D = 0.0;
  ele->getParams("DiffusionCoefficient", D);
  // get wall thickness
  double th = 0.0;
  ele->getParams("WallThickness", th);
  // get percentage of diffusion area
  double percDiffArea = 0.0;
  ele->getParams("PercentageOfDiffusionArea", percDiffArea);
  {
    D *= percDiffArea * (myareanp[ai[0]] / th);
  }
  // define vO2 in air
  double vO2a = vO2[ai[0]];
  // define vO2 in capillary
  double vO2b = 0.0;
  for (unsigned int i = 0; i < ci.size(); i++)
  {
    vO2b += vO2[ci[i]];
  }
  vO2b /= float(ci.size());

  // define vO2ab which is the volume of blood flowing from air to blood
  double vO2ab = 0.0;  // D*(pO2[ai[0]] - tempPO2b)*dt;
  // loop till convergence
  for (unsigned int itr = 0; fabs(error) > tol; itr++)
  {
    // -----------------------------------------------------------------
    // Evaluate all terms associated with flow rate of O2
    // -----------------------------------------------------------------
    // flow rate of O2 from air to blood
    double qO2 = vO2ab / dt;
    // derivative of O2 flowrate w.r.t volume of O2 flowing into blood
    double dqO2_dvO2ab = 1.0 / dt;
    // -----------------------------------------------------------------
    // Evaluate all terms associated with PO2 in air
    // -----------------------------------------------------------------
    // new volume of O2 in air is old volume minus the one flowing into
    // blood
    double vO2ap = vO2a - vO2ab;
    // pO2 in air
    double pO2a = vO2ap * patm / vAir;
    // derivative of pO2a w.r.t volume of O2 flowing into blood
    double dpO2a_dvO2 = -patm / vAir;
    // -----------------------------------------------------------------
    // Evaluate all terms associated with PO2 in blood
    // -----------------------------------------------------------------
    // volumen of O2 in blood is old volume plus the one flowing inot
    // blood
    double vO2bp = vO2b + vO2ab;

    double pO2b = pow(vO2bp / kv, 1.0 / power) * pow(1.0 - vO2bp / kv, -1.0 / power) * ph;
    double dpO2b_dvO2ab = (1.0 / power) * pow(vO2bp / kv, 1.0 / power - 1.0) *
                              pow(1.0 - vO2bp / kv, -1.0 / power) * ph / kv -
                          (1.0 / power) * pow(vO2bp / kv, 1.0 / power) *
                              pow(1.0 - vO2bp / kv, -1.0 / power - 1.0) * ph / (-kv);


    // -----------------------------------------------------------------
    // perform the Newton-Raphson step
    // -----------------------------------------------------------------
    // calculate residual
    double f = qO2 - D * (pO2a - pO2b);
    // calculate df/dvO2ab
    double df_dvO2ab = dqO2_dvO2ab - D * (dpO2a_dvO2 - dpO2b_dvO2ab);
    // calculate the corrector
    double dvO2ab = f / (df_dvO2ab);
    // corrector vO2ab
    vO2ab -= dvO2ab;

    // evaluate a dimensionless error
    error = fabs(f * dt);

    // check append iteration step
    if (itr > maxItr)
    {
      dserror("[Warning in ELEMENT(%d)]: solving for VO2(air/blood) did not converge Error: %f\n",
          ele->Id(), fabs(error));
      exit(1);
    }
  }
  // -----------------------------------------------------------------
  // update the change of vO2 into dvO2 vector:
  // -----------------------------------------------------------------
  // update dscatra in air
  dscatra[ai[0]] = -vO2ab * nO2perVO2a / vAir;
  dvo2[ai[0]] = -vO2ab;
  // update dscatra in blood
  dscatra[ci[0]] = vO2ab * nO2perVO2b / vBlood;
  dvo2[ci[0]] = vO2ab;

  scatra_acinus[ai[0]] = (vO2[ai[0]] - vO2ab) * nO2perVO2a / vAir;
  scatra_acinus[ci[0]] = (vO2[ci[0]] + vO2ab) * nO2perVO2b / vBlood;
}
