/*----------------------------------------------------------------------*/
/*! \file
 \brief bond material

\level 3

\maintainer Christoph Schmidt

 *----------------------------------------------------------------------*/


#include <vector>
#include "scatra_reaction_mat.H"
#include "scatra_reaction_coupling.H"
#include "scatra_bondreac_mat.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_mat/matpar_bundle.H"
#include "../drt_comm/comm_utils.H"


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::PAR::ScatraBondReacMat::ScatraBondReacMat(Teuchos::RCP<MAT::PAR::Material> matdata)
    : ScatraReactionMat(matdata),
      // get thermal energy
      kBT_(DRT::Problem::Instance()->CellMigrationParams().get<double>("kBT")),
      penalty_(DRT::Problem::Instance()
                   ->CellMigrationParams()
                   .sublist("ADHESION MODULE")
                   .get<double>("PENALTY")),
      bondtype_(SetBondType(matdata)),
      slipcoeff_(matdata->GetDouble("SLIPCOEFF")),
      catchcoeff1_(matdata->GetDouble("CATCHCOEFF1")),
      catchcoeff2_(matdata->GetDouble("CATCHCOEFF2")),
      catchcoeff3_(matdata->GetDouble("CATCHCOEFF3")),
      catchcoeff4_(matdata->GetDouble("CATCHCOEFF4")),
      r_bind_(matdata->GetDouble("BINDING_RADIUS")),
      fiber_diameter_(DRT::Problem::Instance()
                          ->CellMigrationParams()
                          .sublist("ADHESION MODULE")
                          .get<double>("ECM_FIBER_DIAMETER"))
{
  // Some checks for more safety
  if (bondtype_ == bondtype_none)
    dserror(
        "The bond type '%s' is not a valid type. Valid bond types are 'no_bond', 'slip_bond' and "
        "'integrin_binding' and 'integrin_rupture'.",
        (matdata->Get<std::string>("TYPE"))->c_str());

  // check if slipcoeff is defined
  if (bondtype_ == bondtype_slip_bond and slipcoeff_ == -1.0)
    dserror("No binding length defined for the slip bond!");

  // check if catchcoeffs are defined
  if (bondtype_ == bondtype_catch_bond and (catchcoeff1_ == -123.0 or catchcoeff2_ == -123.0 or
                                               catchcoeff3_ == -123.0 or catchcoeff4_ == -123.0))
    dserror("No binding length defined for the slip bond!");

  // only one species is allowed to disassemble
  if (bondtype_ == bondtype_slip_bond or bondtype_ == bondtype_catch_bond)
  {
    int counter = 0;
    for (int ii = 0; ii < numscal_; ii++)
    {
      if (stoich_->at(ii) == -1) counter += 1;
    }
    if (counter != 1) dserror("More than one species is disassembled! Fix input file!");
  }

  // the reaction coefficient has to be set to 1.0 for bondtype_integrin_rupture
  // as the reaction coefficient is calculated using the catchbondcoeffs.
  if (bondtype_ == bondtype_catch_bond and reaccoeff_ != 1.0)
    dserror("The reaction coefficient has to be set to 1.0 for bond type bondtype_catch_bond!");

  // check if binding radius is defined
  if (bondtype_ == bondtype_integrin_binding and r_bind_ == -1.0)
    dserror("No binding radius defined for cell-ECM interaction!");

  // check if ECM fiber diameter is defined
  if (bondtype_ == bondtype_integrin_binding and fiber_diameter_ == -1.0)
    dserror("No fiber diameter defined for cell-ECM interaction!");

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<MAT::Material> MAT::PAR::ScatraBondReacMat::CreateMaterial()
{
  return Teuchos::rcp(new MAT::ScatraBondReacMat(this));
}


MAT::ScatraBondReacMatType MAT::ScatraBondReacMatType::instance_;


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::PAR::bond_type MAT::PAR::ScatraBondReacMat::SetBondType(
    Teuchos::RCP<MAT::PAR::Material> matdata)
{
  if (*(matdata->Get<std::string>("BONDTYPE")) == "no_bond")
  {
    return bondtype_no_bond;
  }
  else if (*(matdata->Get<std::string>("BONDTYPE")) == "slip_bond")
  {
    return bondtype_slip_bond;
  }
  else if (*(matdata->Get<std::string>("BONDTYPE")) == "catch_bond")
  {
    return bondtype_catch_bond;
  }
  else if (*(matdata->Get<std::string>("BONDTYPE")) == "integrin_binding")
  {
    return bondtype_integrin_binding;
  }
  else if (*(matdata->Get<std::string>("BONDTYPE")) == "no_bondtype")
  {
    return bondtype_none;
  }
  else
  {
    return bondtype_none;
  }
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
DRT::ParObject* MAT::ScatraBondReacMatType::Create(const std::vector<char>& data)
{
  MAT::ScatraBondReacMat* scatra_bond_mat = new MAT::ScatraBondReacMat();
  scatra_bond_mat->Unpack(data);
  return scatra_bond_mat;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::ScatraBondReacMat::ScatraBondReacMat() : ScatraReactionMat(), params_(NULL) {}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::ScatraBondReacMat::ScatraBondReacMat(MAT::PAR::ScatraBondReacMat* params)
    : ScatraReactionMat(params), params_(params)
{
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ScatraBondReacMat::Pack(DRT::PackBuffer& data) const
{
  DRT::PackBuffer::SizeMarker sm(data);
  sm.Insert();

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data, type);

  // matid
  int matid = -1;
  if (params_ != NULL) matid = params_->Id();  // in case we are in post-process mode
  AddtoPack(data, matid);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ScatraBondReacMat::Unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;
  // extract type
  int type = 0;
  ExtractfromPack(position, data, type);
  if (type != UniqueParObjectId())
    dserror(
        "wrong instance type data. type = %d, UniqueParObjectId()=%d", type, UniqueParObjectId());
  // matid and recover params_
  int matid;
  ExtractfromPack(position, data, matid);
  params_ = NULL;
  if (DRT::Problem::Instance()->Materials() != Teuchos::null)
    if (DRT::Problem::Instance()->Materials()->Num() != 0)
    {
      const int probinst = DRT::Problem::Instance()->Materials()->GetReadFromProblem();
      MAT::PAR::Parameter* mat =
          DRT::Problem::Instance(probinst)->Materials()->ParameterById(matid);
      if (mat->Type() == MaterialType())
        params_ = static_cast<MAT::PAR::ScatraBondReacMat*>(mat);
      else
        dserror("Type of parameter material %d does not fit to calling type %d", mat->Type(),
            MaterialType());
    }

  if (position != data.size()) dserror("Mismatch in size of data %d <-> %d", data.size(), position);
}


/*----------------------------------------------------------------------/
 | calculate advanced reaction terms                        Thon 08/16 |
/----------------------------------------------------------------------*/
double MAT::ScatraBondReacMat::CalcReaBodyForceTerm(const int k,  //!< current scalar id
    const std::vector<double>& phinp,                             //!< scalar values at t_(n+1)
    const std::vector<double>& phin,                              //!< scalar values at t_n
    const double violation,                                       //!< traction at curren gp
    const double porosity,  //!< porosity of background element
    const double
        scale_phi,         //!< scaling factor for scalar values (used for reference concentrations)
    const double* gpcoord  //!< Gauss-point coordinates
    ) const
{
  // add time and space coordinates
  std::vector<std::pair<std::string, double>> constants;
  constants.push_back(std::pair<std::string, double>("t", 0.0));
  constants.push_back(std::pair<std::string, double>("x", gpcoord[0]));
  constants.push_back(std::pair<std::string, double>("y", gpcoord[1]));
  constants.push_back(std::pair<std::string, double>("z", gpcoord[2]));

  if (Stoich()->at(k) != 0 and fabs(ReacCoeff(constants)) > 1.0e-14)
  {
    double ReacCoeffFactor = CalcReacCoeffFactor(violation, porosity, phin);

    return CalcReaBodyForceTerm(k, phinp, ReacCoeffFactor * ReacCoeff(constants) * Stoich()->at(k),
        scale_phi);  // scalar at integration point np
  }
  else
    return 0.0;
}


/*----------------------------------------------------------------------/
 | calculate advanced reaction term derivatives             Thon 08/16 |
/----------------------------------------------------------------------*/
void MAT::ScatraBondReacMat::CalcReaBodyForceDerivMatrix(const int k,  //!< current scalar id
    std::vector<double>& derivs,       //!< vector with derivatives (to be filled)
    const std::vector<double>& phinp,  //!< scalar values at t_(n+1)
    const std::vector<double>& phin,   //!< scalar values at t_n
    const double violation,            //!< penalty violation at current gp
    const double porosity,             //!< porosity of background element
    const double
        scale_phi,         //!< scaling factor for scalar values (used for reference concentrations)
    const double* gpcoord  //!< Gauss-point coordinates
    ) const
{
  // set time and space coordinates
  std::vector<std::pair<std::string, double>> constants;
  constants.push_back(std::pair<std::string, double>("t", 0.0));
  constants.push_back(std::pair<std::string, double>("x", gpcoord[0]));
  constants.push_back(std::pair<std::string, double>("y", gpcoord[1]));
  constants.push_back(std::pair<std::string, double>("z", gpcoord[2]));

  if (Stoich()->at(k) != 0 and fabs(ReacCoeff(constants)) > 1.0e-14)
  {
    // get reaction coefficient derivative force scale factor
    double ReacCoeffDerivFactor = CalcReacCoeffFactor(violation, porosity, phinp);

    CalcReaBodyForceDeriv(k, derivs, phinp, constants,
        ReacCoeffDerivFactor * ReacCoeff(constants) * Stoich()->at(k), scale_phi);
  }

  return;
}


/*----------------------------------------------------------------------/
 | calculate advanced reaction terms                        Thon 08/16 |
/----------------------------------------------------------------------*/
double MAT::ScatraBondReacMat::CalcReaBodyForceTerm(const int k,  //!< current scalar id
    const std::vector<double>& phinp,                             //!< scalar values at t_(n+1)
    const std::vector<double>& phin,                              //!< scalar values at t_n
    const std::vector<std::pair<std::string, double>>&
        constants,           //!< vector containing values which are independent of the scalars
    const double violation,  //!< traction at curren gp
    const double porosity,   //!< porosity of background element
    const double
        scale_phi,         //!< scaling factor for scalar values (used for reference concentrations)
    const double* gpcoord  //!< Gauss-point coordinates
    ) const
{
  // add time and space coordinates
  std::vector<std::pair<std::string, double>> constants_mod(constants);
  constants_mod.push_back(std::pair<std::string, double>("t", 0.0));
  constants_mod.push_back(std::pair<std::string, double>("x", gpcoord[0]));
  constants_mod.push_back(std::pair<std::string, double>("y", gpcoord[1]));
  constants_mod.push_back(std::pair<std::string, double>("z", gpcoord[2]));

  if (Stoich()->at(k) != 0 and fabs(ReacCoeff(constants)) > 1.0e-14)
  {
    double ReacCoeffFactor = CalcReacCoeffFactor(violation, porosity, phin);

    return CalcReaBodyForceTerm(k, phinp, constants_mod,
        ReacCoeffFactor * ReacCoeff(constants_mod) * Stoich()->at(k),
        scale_phi);  // scalar at integration point np
  }
  else
    return 0.0;
}


/*----------------------------------------------------------------------/
 | calculate advanced reaction term derivatives             Thon 08/16 |
/----------------------------------------------------------------------*/
void MAT::ScatraBondReacMat::CalcReaBodyForceDerivMatrix(const int k,  //!< current scalar id
    std::vector<double>& derivs,       //!< vector with derivatives (to be filled)
    const std::vector<double>& phinp,  //!< scalar values at t_(n+1)
    const std::vector<double>& phin,   //!< scalar values at t_n
    const std::vector<std::pair<std::string, double>>&
        constants,          //!< vector containing values which are independent of the scalars
    const double traction,  //!< traction at curren gp
    const double porosity,  //!< porosity of background element
    const double
        scale_phi,         //!< scaling factor for scalar values (used for reference concentrations)
    const double* gpcoord  //!< Gauss-point coordinates
    ) const
{
  std::vector<std::pair<std::string, double>> constants_mod(constants);
  constants_mod.push_back(std::pair<std::string, double>("t", 0.0));
  constants_mod.push_back(std::pair<std::string, double>("x", gpcoord[0]));
  constants_mod.push_back(std::pair<std::string, double>("y", gpcoord[1]));
  constants_mod.push_back(std::pair<std::string, double>("z", gpcoord[2]));

  if (Stoich()->at(k) != 0 and fabs(ReacCoeff(constants_mod)) > 1.0e-14)
  {
    double ReacCoeffFactor = CalcReacCoeffFactor(traction, porosity, phin);

    CalcReaBodyForceDeriv(k, derivs, phinp, constants_mod,
        ReacCoeffFactor * ReacCoeff(constants_mod) * Stoich()->at(k), scale_phi);
  }

  return;
}


/*----------------------------------------------------------------------*
 |  helper for calculating advanced reaction terms           thon 08/16 |
 *----------------------------------------------------------------------*/
double MAT::ScatraBondReacMat::CalcReaBodyForceTerm(int k,  //!< current scalar id
    const std::vector<double>& phinp,                       //!< scalar values at t_(n+1)
    double
        scale_reac,   //!< scaling factor for reaction term (= reaction coefficient * stoichometry)
    double scale_phi  //!< scaling factor for scalar values (used for reference concentrations)
    ) const
{
  // set time and space coordinates
  std::vector<std::pair<std::string, double>> constants;
  constants.push_back(std::pair<std::string, double>("t", 0.0));
  constants.push_back(std::pair<std::string, double>("x", 0.0));
  constants.push_back(std::pair<std::string, double>("y", 0.0));
  constants.push_back(std::pair<std::string, double>("z", 0.0));

  return params_->reaction_->CalcReaBodyForceTerm(
      k, NumScal(), phinp, constants, *Couprole(), scale_reac, scale_phi);
}


/*----------------------------------------------------------------------*
 |  helper for calculating advanced reaction terms           thon 08/16 |
 *----------------------------------------------------------------------*/
double MAT::ScatraBondReacMat::CalcReaBodyForceTerm(int k,  //!< current scalar id
    const std::vector<double>& phinp,                       //!< scalar values at t_(n+1)
    const std::vector<std::pair<std::string, double>>&
        constants,  //!< vector containing values which are independent of the scalars
    double
        scale_reac,   //!< scaling factor for reaction term (= reaction coefficient * stoichometry)
    double scale_phi  //!< scaling factor for scalar values (used for reference concentrations)
    ) const
{
  return params_->reaction_->CalcReaBodyForceTerm(
      k, NumScal(), phinp, constants, *Couprole(), scale_reac, scale_phi);
}


/*----------------------------------------------------------------------/
 | calculate the force/porosity dependency of the reaction coefficient  |
/----------------------------------------------------------------------*/
double MAT::ScatraBondReacMat::CalcReacCoeffFactor(
    const double violation,           //!< penalty violation at current gauss point
    const double porosity,            //!< porosity of background element
    const std::vector<double>& phinp  //!< scalar values at t_(n+1)
    ) const
{
  double fac = 1.0;  //!< reaction rate scaling factor

  // thermal energy
  const double kBT = params_->kBT_;

  // adhesion penalty parameter
  double penalty = params_->penalty_;

  switch (params_->bondtype_)
  {
    case MAT::PAR::bondtype_no_bond:
    {
      // --> no force-dependency --> do nothing here!
      break;
    }
    case MAT::PAR::bondtype_slip_bond:
    {
      // slip bond model according to Bell, 1978
      // k_0: reaction rate of unstressed bond
      // gamma: binding length of specific bond
      // f: traction of single bond (penalty * violation)
      // kBT: thermal energy
      // reaction rate k = k_0 * exp(gamma * f / kBT)

      // get slip bond coefficient
      double gamma = params_->slipcoeff_;

      // limit exponent to avoid errors in exponential function
      if (gamma * penalty * violation / kBT > 45) penalty = 45 * kBT / gamma / violation;

      // scaling factor according to Bell, 1978
      fac = exp(gamma * penalty * violation / kBT);

      break;
    }
    case MAT::PAR::bondtype_catch_bond:
    {
      // model according to Bell model (see slip bond), but with two factors:
      // one slip pathway and one catch pathway. Small forces therefore prolong
      // bond lifetime, whereas high forces reduce it.
      // k_c: reaction rate constant of catch pathway
      // x_c: binding length of catch pathway
      // k_s: reaction rate constant of slip pathway
      // x_s: binding length of slip pathway
      // f: single bond traction (f = F/conc = penalty * violation)
      // reaction rate k = k_c * exp(x_c * f / kBT) + k_s * exp(x_s * f / kBT)

      // get catch bond constants
      const double coeff_a = params_->catchcoeff1_;
      const double coeff_b = params_->catchcoeff2_;
      const double coeff_c = params_->catchcoeff3_;
      const double coeff_d = params_->catchcoeff4_;

      // limit exponent to avoid errors in exponential function
      if (penalty * violation > 400) penalty = 400 / violation;

      fac = coeff_a * exp(coeff_b * penalty * violation / kBT) +
            coeff_c * exp(coeff_d * penalty * violation / kBT);

      break;
    }
    case MAT::PAR::bondtype_integrin_binding:
    {
      // integrin binding can be seen as two step process: first, the ligand and receptor
      // have to get close to each other within the binding distance l_bind; then, the
      // actual reaction takes place, characterized by the reaction rate k.
      // Here, the first step is modeled via a probability rho of the ligand and receptor
      // to be close to each other, which mainly depends on the porosity.
      // The actual reaction rate is scaled with this probability.

      // define constants
      const double l_bind = params_->r_bind_;  //!< integrin binding radius in microns
      const double fiber_diameter =
          params_->fiber_diameter_;  //!< effective ECM fiber diameter in microns

      // the probability is modeled according to the mikado model of Metzner et al., 2011
      fac = 1.0 -
            exp(-4.0 * (l_bind * l_bind) / (fiber_diameter * fiber_diameter) * (1.0 - porosity));

      break;
    }
    default:
    {
      dserror("Force dependency for selected bond type not implemented!");
    }
  }  // switch(params_->bondtype_)

  return fac;
}
