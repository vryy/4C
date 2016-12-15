/*----------------------------------------------------------------------*/
/*!
 \file scatra_bondreac_mat.cpp

 \brief bond material

\level 3

 \maintainer  Andreas Rauch

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
MAT::PAR::ScatraBondReacMat::ScatraBondReacMat(
    Teuchos::RCP<MAT::PAR::Material> matdata
)
: ScatraReactionMat(matdata),
  bondtype_(SetBondType(matdata)),
  bindinglength_((matdata->GetDouble("GAMMA")))
{
  //Some checks for more safety
  if (bondtype_ == bondtype_none)
    dserror("The bond type '%s' is not a valid type. Valid bond types are 'no_bond', 'slip_bond' and 'integrin_binding' and 'integrin_rupture'.",(matdata->Get<std::string >("TYPE"))->c_str() );

  if (bondtype_ == bondtype_slip_bond  &&  bindinglength_==-1.0)
    dserror("No binding length defined for the slip bond!");

  // only one species is allowed to disassemble
  if ( bondtype_ != bondtype_no_bond and bondtype_ != bondtype_integrin_binding )
  {
    int counter = 0;
    for (int ii=0; ii < numscal_; ii++)
    {
      if (stoich_->at(ii) == -1)
        counter += 1;
    }
    if (counter != 1)
      dserror("reac_coup_power_multiplicative must contain at least one positive entry in the ROLE list");
  }

  // the reaction coefficient has to be set to 1.0 for bondtype_integrin_rupture,
  // as the real reaction coefficient is set in MAT::ScatraBondReacMat::AdjustReacCoeff.
  if ( bondtype_ == bondtype_integrin_rupture and  reaccoeff_!=1.0)
    dserror("The reaction coefficient has to be set to 1.0 for bond types integrin_binding and integrin_rupture!");

  return;
}


Teuchos::RCP<MAT::Material> MAT::PAR::ScatraBondReacMat::CreateMaterial()
{
  return Teuchos::rcp(new MAT::ScatraBondReacMat(this));
}


MAT::ScatraBondReacMatType MAT::ScatraBondReacMatType::instance_;


MAT::PAR::bond_type MAT::PAR::ScatraBondReacMat::SetBondType( Teuchos::RCP<MAT::PAR::Material> matdata )
{
  if ( *(matdata->Get<std::string >("BONDTYPE")) == "no_bond" )
  {
    return bondtype_no_bond;
  }
  else if ( *(matdata->Get<std::string >("BONDTYPE")) == "slip_bond" )
  {
    return bondtype_slip_bond;
  }
  else if ( *(matdata->Get<std::string >("BONDTYPE")) == "integrin_binding")
  {
    return bondtype_integrin_binding;
  }
  else if ( *(matdata->Get<std::string >("BONDTYPE")) == "integrin_rupture")
  {
    return bondtype_integrin_rupture;
  }
  else if ( *(matdata->Get<std::string >("BONDTYPE")) == "no_bondtype")
  {
    return bondtype_none;
  }
  else
  {
    return bondtype_none;
  }
}

DRT::ParObject* MAT::ScatraBondReacMatType::Create( const std::vector<char> & data )
{
  MAT::ScatraBondReacMat* scatra_bond_mat = new MAT::ScatraBondReacMat();
  scatra_bond_mat->Unpack(data);
  return scatra_bond_mat;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::ScatraBondReacMat::ScatraBondReacMat()
: ScatraReactionMat(),
  params_(NULL)
{
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::ScatraBondReacMat::ScatraBondReacMat(MAT::PAR::ScatraBondReacMat* params)
: ScatraReactionMat(params),
  params_(params)
{
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ScatraBondReacMat::Pack(DRT::PackBuffer& data) const
{
  DRT::PackBuffer::SizeMarker sm( data );
  sm.Insert();

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data,type);

  // matid
  int matid = -1;
  if (params_ != NULL) matid = params_->Id();  // in case we are in post-process mode
  AddtoPack(data,matid);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ScatraBondReacMat::Unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;
  // extract type
  int type = 0;
  ExtractfromPack(position,data,type);
  if (type != UniqueParObjectId()) dserror("wrong instance type data. type = %d, UniqueParObjectId()=%d",type,UniqueParObjectId());
  // matid and recover params_
  int matid;
  ExtractfromPack(position,data,matid);
  params_ = NULL;
  if (DRT::Problem::Instance()->Materials() != Teuchos::null)
    if (DRT::Problem::Instance()->Materials()->Num() != 0)
    {
      const int probinst = DRT::Problem::Instance()->Materials()->GetReadFromProblem();
      MAT::PAR::Parameter* mat = DRT::Problem::Instance(probinst)->Materials()->ParameterById(matid);
      if (mat->Type() == MaterialType())
        params_ = static_cast<MAT::PAR::ScatraBondReacMat*>(mat);
      else
        dserror("Type of parameter material %d does not fit to calling type %d", mat->Type(), MaterialType());
    }

  if (position != data.size())
    dserror("Mismatch in size of data %d <-> %d",data.size(),position);
}


/*----------------------------------------------------------------------/
 | calculate advanced reaction terms                        Thon 08/16 |
/----------------------------------------------------------------------*/
double MAT::ScatraBondReacMat::CalcReaBodyForceTerm(
    const int k,                         //!< current scalar id
    const std::vector<double>& phinp,    //!< scalar values at t_(n+1)
    const std::vector<double>& phin,     //!< scalar values at t_n
    const double traction,               //!< traction at curren gp
    const double porosity,               //!< porosity of background element
    const double scale_phi,              //!< scaling factor for scalar values (used for reference concentrations)
    const double* gpcoord                //!< Gauss-point coordinates
    ) const
{
  if ( Stoich()->at(k)!=0 and fabs(ReacCoeff(gpcoord))>1.0e-14)
  {
    double ReacCoeffFactor = AdjustReacCoeff(traction,porosity,phin,k);

    return CalcReaBodyForceTerm(k,phinp,ReacCoeffFactor*ReacCoeff(gpcoord)*Stoich()->at(k),scale_phi);// scalar at integration point np
  }
  else
    return 0.0;
}

/*----------------------------------------------------------------------/
 | calculate advanced reaction term derivatives             Thon 08/16 |
/----------------------------------------------------------------------*/
void MAT::ScatraBondReacMat::CalcReaBodyForceDerivMatrix(
    const int k,                         //!< current scalar id
    std::vector<double>& derivs,         //!< vector with derivatives (to be filled)
    const std::vector<double>& phinp,    //!< scalar values at t_(n+1)
    const std::vector<double>& phin,     //!< scalar values at t_n
    const double traction,               //!< traction at curren gp
    const double porosity,               //!< porosity of background element
    const double scale_phi,              //!< scaling factor for scalar values (used for reference concentrations)
    const double* gpcoord                //!< Gauss-point coordinates
    ) const
{
  if ( Stoich()->at(k)!=0 and fabs(ReacCoeff(gpcoord))>1.0e-14)
  {
    double ReacCoeffFactor = AdjustReacCoeff(traction,porosity,phin,k);

    CalcReaBodyForceDeriv(k,derivs,phinp,ReacCoeffFactor*ReacCoeff(gpcoord)*Stoich()->at(k),scale_phi);
  }

  return;
}


/*----------------------------------------------------------------------/
 | calculate advanced reaction terms                        Thon 08/16 |
/----------------------------------------------------------------------*/
double MAT::ScatraBondReacMat::CalcReaBodyForceTerm(
    const int k,                                                     //!< current scalar id
    const std::vector<double>& phinp,                                //!< scalar values at t_(n+1)
    const std::vector<double>& phin,                                 //!< scalar values at t_n
    const std::vector<std::pair<std::string,double> >& constants,    //!< vector containing values which are independent of the scalars
    const double traction,                                           //!< traction at curren gp
    const double porosity,                                           //!< porosity of background element
    const double scale_phi,                                          //!< scaling factor for scalar values (used for reference concentrations)
    const double* gpcoord                                            //!< Gauss-point coordinates
    ) const
{
  if ( Stoich()->at(k)!=0 and fabs(ReacCoeff(gpcoord))>1.0e-14)
  {
    double ReacCoeffFactor = AdjustReacCoeff(traction,porosity,phin,k);

    return  CalcReaBodyForceTerm(k,phinp,constants,ReacCoeffFactor*ReacCoeff(gpcoord)*Stoich()->at(k),scale_phi);// scalar at integration point np
  }
  else
    return 0.0;
}

/*----------------------------------------------------------------------/
 | calculate advanced reaction term derivatives             Thon 08/16 |
/----------------------------------------------------------------------*/
void MAT::ScatraBondReacMat::CalcReaBodyForceDerivMatrix(
    const int k,                                                     //!< current scalar id
    std::vector<double>& derivs,                                     //!< vector with derivatives (to be filled)
    const std::vector<double>& phinp,                                //!< scalar values at t_(n+1)
    const std::vector<double>& phin,                                 //!< scalar values at t_n
    const std::vector<std::pair<std::string,double> >& constants,    //!< vector containing values which are independent of the scalars
    const double traction,                                           //!< traction at curren gp
    const double porosity,                                           //!< porosity of background element
    const double scale_phi,                                          //!< scaling factor for scalar values (used for reference concentrations)
    const double* gpcoord                //!< Gauss-point coordinates
    ) const
{
  if ( Stoich()->at(k)!=0 and fabs(ReacCoeff(gpcoord))>1.0e-14)
  {
    double ReacCoeffFactor = AdjustReacCoeff(traction,porosity,phin,k);

    CalcReaBodyForceDeriv(k,derivs,phinp,constants,ReacCoeffFactor*ReacCoeff(gpcoord)*Stoich()->at(k),scale_phi);
  }

  return;
}






/*----------------------------------------------------------------------*
 |  helper for calculating advanced reaction terms           thon 08/16 |
 *----------------------------------------------------------------------*/
double MAT::ScatraBondReacMat::CalcReaBodyForceTerm(
    int k,                               //!< current scalar id
    const std::vector<double>& phinp,    //!< scalar values at t_(n+1)
    double scale_reac,                   //!< scaling factor for reaction term (= reaction coefficient * stoichometry)
    double scale_phi                     //!< scaling factor for scalar values (used for reference concentrations)
    ) const
{
  return params_->reaction_->CalcReaBodyForceTerm(k,NumScal(),phinp,*Couprole(),scale_reac,scale_phi);
}


/*----------------------------------------------------------------------*
 |  helper for calculating advanced reaction terms           thon 08/16 |
 *----------------------------------------------------------------------*/
double MAT::ScatraBondReacMat::CalcReaBodyForceTerm(
    int k,                                                           //!< current scalar id
    const std::vector<double>& phinp,                                //!< scalar values at t_(n+1)
    const std::vector<std::pair<std::string,double> >& constants,    //!< vector containing values which are independent of the scalars
    double scale_reac,                                               //!< scaling factor for reaction term (= reaction coefficient * stoichometry)
    double scale_phi                                                 //!< scaling factor for scalar values (used for reference concentrations)
    ) const
{
  return params_->reaction_->CalcReaBodyForceTerm(k,NumScal(),phinp,constants,*Couprole(),scale_reac,scale_phi);
}







/*----------------------------------------------------------------------/
 | calculate the force/porosity dependency of the reaction coefficient  |
/----------------------------------------------------------------------*/
double MAT::ScatraBondReacMat::AdjustReacCoeff(
    const double traction,              //!< traction at current gauss point
    const double porosity,              //!< porosity of background element
    const std::vector<double>& phin,    //!< scalar values at t_(n)
    const int k                         //!< current scalar id
) const
{

  double fac= 1.0;                      //!< reaction rate scaling factor

  switch( BondType() )
  {
  case MAT::PAR::bondtype_no_bond:
  {
    // in the case of bondtype_no_bond, the reaction coefficient is constant
    // --> do nothing here!
    break;
  }
  case MAT::PAR::bondtype_slip_bond:
  {
    // slip bond model according to Bell 1978
    // k_0: reaction rate of unstressed bond
    // gamma: binding length of specific bond
    // f: traction of single bond
    // kBT: thermal energy
    // reaction rate k = k_0 * exp(gamma * f / kBT)

    // define constants
    const double therm_nrg= 4.04530016e-3;      //!< thermal energy at 293K

    // get traction of one single bond
    double single_bond_trac = CalcSingleBondTrac(traction,phin);

    // restrict single_bond_trac to values <= 1.0e100/BindingLength/therm_nrg
    // (Otherwise, the calculated vector norm for the concentration is going to be infinity!)
    if (BindingLength()*single_bond_trac/therm_nrg>20)
    {
      single_bond_trac = 20.0*therm_nrg/BindingLength();
    }

    // scaling factor according to Bell 1978
    fac = exp(BindingLength()*single_bond_trac/therm_nrg);

    break;
  }
  case MAT::PAR::bondtype_integrin_binding:
  {
    // integrin binding can be seen as two step process: first, the ligand and receptor
    // have to get close to each other within the binding distance l_bind; then, the
    // actual reaction takes place, characterized by the reaction rate k.
    // Here, the first step is modeled via a probability rho of the ligand and receptor
    // to be close to each other, which mainly depends on the porosity.
    // The actual reaction rate is therefore scaled with this probability.



    // define constants
    const double l_bind = 23.0e-3;             //!< binding radius in micrometers
    const double fiber_diameter = 25.0e-3;     //!< effective ECM fiber diameter in micrometers
    const double coll_spec_vol = 1.89e9;      //!< effective specific volume of collagen in micrometers3/mg
    const double coll_mol_mass = 300.0;       //!< collagen molar mass in mg/micromoles
    const double avogadro = 6.022141e17;      //!< Avogadro's number in micromoles-1

    // the probability is modeled according to the mikado model of Metzner et al. 2011
    double rho = 1.0 - exp( -4.0 * pow(l_bind/fiber_diameter,2) * (1.0-porosity) );

    // calculate the collagen concentration based on the porosity
    // (multiplied by avogadro's number to get number of species from moles)
    double coll_conc = (1.0-porosity)/(coll_spec_vol * coll_mol_mass) * avogadro;

    // reaction rate has to be multiplied by both the probability of the reactants
    // to be close and the unbound ligand concentration
    fac = rho * coll_conc;

    break;
  }
  case MAT::PAR::bondtype_integrin_rupture:
  {
    // scaling factor curve fitted to data from Kong et al. 2007
    // model according to Bell model (see slip bond), but with two factors:
    // one slip pathway and one catch pathway. Small forces therefore prolong
    // bond lifetime, whereas high forces reduce it.
    // k_c: reaction rate constant of catch pathway
    // x_c: binding length of catch pathway
    // k_s: reaction rate constant of slip pathway
    // x_s: binding length of slip pathway
    // f: single bond traction
    // reaction rate k = k_c * exp(x_c * f / kBT) + k_s * exp(x_s * f / kBT)

    // define constants
    const double coeff_a =  1.082719e+01;
    const double coeff_b = -1.719530e-01;
    const double coeff_c =  1.825432e-03;
    const double coeff_d =  1.262438e-01;

    // get the traction of one single bond
    double single_bond_trac = CalcSingleBondTrac(traction,phin);

    // restrict single_bond_trac to values <= 1.0e100/x_s/therm_nrg
    // (Otherwise, the calculated vector norm for the concentration is going to be infinity in some cases!)
    if (single_bond_trac>200)
      single_bond_trac = 200;

    fac = coeff_a*exp(coeff_b*single_bond_trac) + coeff_c*exp(coeff_d*single_bond_trac);

    break;
  }
  default:
  {
    dserror("Reaction rate dependency for selected bond type not implemented!");
  }
  }

  return fac;
}


/*----------------------------------------------------------------------/
 | calculate the single bond traction                                   |
/----------------------------------------------------------------------*/
double MAT::ScatraBondReacMat::CalcSingleBondTrac(
    const double traction,              //!< traction at current gauss point
    const std::vector<double>& phin     //!< scalar values at t_(n)
) const
{
  // get scalar of species that is disassembled ( Stoich()==-1 )
  int bond_scalar=-1234;

  for (int ii=0; ii < NumScal(); ii++)
  {
    if (Stoich()->at(ii) == -1)
    {
      // check if there is only one species disassembled
      if (bond_scalar != -1234)
        dserror("More than one species is disassembled in integrin rupture reaction. Check definition!");

      bond_scalar = ii;
      break;
    }
  }
  if ( bond_scalar==-1234 )
    dserror("Couldn't get scalar of disassembled bond. Check your reaction definition!");

  // traction experienced by one single bond = traction/concentration
  // the concentration of t_(n) has to be used, otherwise the process is self-energizing
  // check whether phin > 0 (too high forces otherwise!)
  double single_bond_trac=0;

  if ( phin.at(bond_scalar) > 1e-10 )
  {
    single_bond_trac = traction/phin.at(bond_scalar);
  }

  return single_bond_trac;
}
