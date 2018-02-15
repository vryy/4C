/*!----------------------------------------------------------------------*/
/*!
\file electrode.cpp

\brief electrode material carrying concentration and electric potential as degrees of freedom

\level 2

<pre>
\maintainer Rui Fang
            fang@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/
            089-289-15251
</pre>
*/
/*----------------------------------------------------------------------*/
#include "electrode.H"

#include "../drt_lib/drt_globalproblem.H"

#include "../drt_mat/matpar_bundle.H"

/*----------------------------------------------------------------------*
 | constructor                                               fang 08/15 |
 *----------------------------------------------------------------------*/
MAT::PAR::Electrode::Electrode(
    Teuchos::RCP<MAT::PAR::Material> matdata
    ) :
  ElchSingleMat(matdata),
  cmax_(matdata->GetDouble("C_MAX")),
  ocpmodel_(StringToOCPModel(*matdata->Get<std::string>("OCP_MODEL"))),
  ocpparanum_(matdata->GetInt("OCP_PARA_NUM")),
  ocppara_(*matdata->Get<std::vector<double> >("OCP_PARA")),
  xmin_(matdata->GetDouble("X_MIN")),
  xmax_(matdata->GetDouble("X_MAX"))
{
  // safety checks
  if(cmax_ < 1.e-12)
    dserror("Saturation value c_max of intercalated Lithium concentration is too small!");
  if(ocpparanum_ < 1)
    dserror("No parameters found for electrode half cell open circuit potential!");
  if(ocpmodel_ == ocp_taralov and ocpparanum_ != 13)
    dserror("Electrode half cell open circuit potential according to Taralov, Taralova, Popov, Iliev, Latz, and Zausch (2012) needs to be specified by exactly 13 coefficients!");
  if((int) ocppara_.size() != ocpparanum_)
    dserror("Length of coefficient vector for electrode half cell open circuit potential doesn't match prescribed number of coefficients!");
  if((xmin_ > 1.0) or (xmax_ > 1.0))
    dserror("Lower bound (X_MIN) and upper bound (X_MAX) of range of validity for ocp calculation model cannot be larger than one since X "
        "is calculated as c/c_max! If you do not want to prescribe bounds, you have to set the two variables to negative values. "
        "If you set the bounds to realistic values (i.e. [0,1]) you will get a warning printed to the screen if bounds are violated throughout the simulation time!");
  if(xmin_ > xmax_)
    dserror("X_MIN cannot be larger than X_MAX!");

  return;
}


/*----------------------------------------------------------------------*
 | create instance of electrode material                     fang 02/15 |
 *----------------------------------------------------------------------*/
Teuchos::RCP<MAT::Material> MAT::PAR::Electrode::CreateMaterial()
{
  return Teuchos::rcp(new MAT::Electrode(this));
}


/*---------------------------------------------------------------------------*
 | convert string to model for half cell open circuit potential   fang 08/15 |
 *---------------------------------------------------------------------------*/
MAT::PAR::OCPModels MAT::PAR::Electrode::StringToOCPModel(const std::string& ocpmodelstring) const
{
  OCPModels ocpmodelenum(ocp_undefined);

  // Redlich-Kister expansion
  if(ocpmodelstring == "Redlich-Kister")
    ocpmodelenum = ocp_redlichkister;

  // empirical correlation given in Taralov, Taralova, Popov, Iliev, Latz, and Zausch (2012)
  else if(ocpmodelstring == "Taralov")
    ocpmodelenum = ocp_taralov;

  // polynomial
  else if(ocpmodelstring == "Polynomial")
    ocpmodelenum = ocp_polynomial;

  // unknown model
  else
    ocpmodelenum = ocp_undefined;

  return ocpmodelenum;
}


MAT::ElectrodeType MAT::ElectrodeType::instance_;


DRT::ParObject* MAT::ElectrodeType::Create(const std::vector<char>& data)
{
  MAT::Electrode* electrode = new MAT::Electrode();
  electrode->Unpack(data);
  return electrode;
}


/*----------------------------------------------------------------------*
 | construct empty electrode material                        fang 02/15 |
 *----------------------------------------------------------------------*/
MAT::Electrode::Electrode() :
  params_(NULL)
{
  return;
}


/*-----------------------------------------------------------------------------*
 | construct electrode material with specific material parameters   fang 02/15 |
 *-----------------------------------------------------------------------------*/
MAT::Electrode::Electrode(MAT::PAR::Electrode* params) :
  params_(params)
{
  return;
}


/*----------------------------------------------------------------------*
 | pack material for communication purposes                  fang 02/15 |
 *----------------------------------------------------------------------*/
void MAT::Electrode::Pack(DRT::PackBuffer& data) const
{
  DRT::PackBuffer::SizeMarker sm( data );
  sm.Insert();

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data,type);

  int matid = -1;
  if(params_ != NULL)
    matid = params_->Id();  // in case we are in post-process mode
  AddtoPack(data,matid);

  return;
}


/*----------------------------------------------------------------------*
 | unpack data from a char vector                            fang 02/15 |
 *----------------------------------------------------------------------*/
void MAT::Electrode::Unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;

  // extract type
  int type = 0;
  ExtractfromPack(position,data,type);
  if(type != UniqueParObjectId())
    dserror("Wrong instance type data!");

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
        params_ = static_cast<MAT::PAR::Electrode*>(mat);
      else
        dserror("Type of parameter material %d does not fit to calling type %d", mat->Type(), MaterialType());
    }

  if(position != data.size())
    dserror("Mismatch in size of data %d <-> %d",data.size(),position);

  return;
}


/*----------------------------------------------------------------------*
 | compute half cell open circuit potential                  fang 08/15 |
 *----------------------------------------------------------------------*/
double MAT::Electrode::ComputeOpenCircuitPotential(
    const double concentration,   //!< concentration
    const double faraday,         //!< Faraday constant
    const double frt              //!< factor F/RT
    ) const
{
  double ocp(0.);

  // intercalation fraction
  const double X = concentration/params_->cmax_;

  // print warning to screen if prescribed interval of validity for ocp calculation model is given but not satisfied
  if(((X < params_->xmin_) or (X > params_->xmax_)) and !(params_->xmax_ < 0.))
  {
    std::cout << "WARNING: intercalation fraction X = c/c_max is violating prescribed bounds of ocp calculation model. Calculated "
        "values might therefore not be reasonable!" << std::endl;
    std::cout << "X: " << X << " lower bound is: " << params_->xmin_ << "  upper bound is: " << params_->xmax_ << std::endl << std::endl;
  }

  // physically reasonable intercalation fraction
  if(X > 0. and X < 1.)
  {
    switch(params_->ocpmodel_)
    {
      // half cell open circuit potential according to Redlich-Kister expansion
      case MAT::PAR::ocp_redlichkister:
      {
        // cf. Colclasure and Kee, Electrochimica Acta 55 (2010) 8960:
        // ocppara_[0]         = DeltaG
        // ocppara_[1,2,3,...] = Redlich-Kister coefficients

        // terms not associated with any Redlich-Kister coefficient
        ocp = params_->ocppara_[0] + faraday/frt*log((1.-X)/X);

        // terms associated with first and second Redlich-Kister coefficients
        // these two terms are separated from the remaining sum and simplified thereafter to remove singularities in the expansion in case X == 0.5
        ocp += params_->ocppara_[1]*(2.*X-1.)+params_->ocppara_[2]*(6.*X*X-6.*X+1.);

        // terms associated with remaining Redlich-Kister coefficients
        for(int i=2; i<params_->ocpparanum_-1; ++i)
          ocp += params_->ocppara_[i+1]*(pow(2.*X-1.,i+1)-2.*i*X*(1.-X)*pow(2.*X-1.,i-1));

        // final scaling
        ocp /= faraday;

        break;
      }

      case MAT::PAR::ocp_taralov:
      {
        // cf. Taralov, Taralova, Popov, Iliev, Latz, and Zausch (2012)
        ocp = params_->ocppara_[0]+params_->ocppara_[1]*tanh(params_->ocppara_[2]*X+params_->ocppara_[3])
              +params_->ocppara_[4]*exp(params_->ocppara_[5]*pow(X,8.0))
              +params_->ocppara_[6]*(1/(pow((params_->ocppara_[7]-X),params_->ocppara_[8]))+params_->ocppara_[9])
              +params_->ocppara_[10]*exp(params_->ocppara_[11]*(X+params_->ocppara_[12]));

        break;
      }

      // polynomial ocp
      case MAT::PAR::ocp_polynomial:
      {
        // add constant
        ocp = params_->ocppara_[0];
        // add higher polynomial order terms
        for(int i=1; i<params_->ocpparanum_; ++i)
        {
          ocp += params_->ocppara_[i]*pow(X,i);
        }

        break;
      }

      default:
      {
        dserror("Model for half cell open circuit potential not recognized!");
        break;
      }
    }
  }

  // non-physical intercalation fraction
  else
    ocp = std::numeric_limits<double>::infinity();

  return ocp;
} // MAT::Electrode::ComputeOpenCircuitPotential


/*---------------------------------------------------------------------------------------------------------*
 | compute first derivative of half cell open circuit potential with respect to concentration   fang 08/15 |
 *---------------------------------------------------------------------------------------------------------*/
double MAT::Electrode::ComputeFirstDerivOpenCircuitPotential(
    const double concentration,   //!< concentration
    const double faraday,         //!< Faraday constant
    const double frt              //!< factor F/RT
    ) const
{
  double ocpderiv(0.);

  // intercalation fraction
  const double X = concentration/params_->cmax_;

  // physically reasonable intercalation fraction
  if(X > 0. and X < 1.)
  {
    switch(params_->ocpmodel_)
    {
      // derivative of half cell open circuit potential w.r.t. concentration according to Redlich-Kister expansion
      case MAT::PAR::ocp_redlichkister:
      {
        // cf. Colclasure and Kee, Electrochimica Acta 55 (2010) 8960:
        // ocppara_[0]         = DeltaG
        // ocppara_[1,2,3,...] = Redlich-Kister coefficients

        // term not associated with any Redlich-Kister coefficient
        ocpderiv = faraday/(2.*frt*X*(X-1.));

        // terms associated with first, second, and third Redlich-Kister coefficients
        // these three terms are separated from the remaining sum and simplified thereafter to remove singularities in the derivative of the expansion in case X == 0.5
        ocpderiv += params_->ocppara_[1]+params_->ocppara_[2]*(6.*X-3.)+params_->ocppara_[3]*(24.*X*X-24.*X+5.);

        // terms associated with remaining Redlich-Kister coefficients
        for(int i=3; i<params_->ocpparanum_-1; ++i)
          ocpderiv += params_->ocppara_[i+1]*((2.*i+1.)*pow(2.*X-1.,i)+2.*X*i*(X-1.)*(i-1.)*pow(2.*X-1.,i-2));

        // final scaling
        ocpderiv *= 2./(faraday*params_->cmax_);

        break;
      }

      // derivative of half cell open circuit potential w.r.t. concentration according to Taralov, Taralova, Popov, Iliev, Latz, and Zausch (2012)
      case MAT::PAR::ocp_taralov:
      {
        ocpderiv = params_->ocppara_[1]*params_->ocppara_[2]/pow(cosh(params_->ocppara_[2]*X+params_->ocppara_[3]),2)
                   +8.*params_->ocppara_[4]*params_->ocppara_[5]*exp(params_->ocppara_[5]*pow(X,8))*pow(X,7)
                   +params_->ocppara_[6]*params_->ocppara_[8]/pow(params_->ocppara_[7]-X,params_->ocppara_[8]+1.)
                   +params_->ocppara_[10]*params_->ocppara_[11]*exp(params_->ocppara_[11]*(X+params_->ocppara_[12]));
        ocpderiv /= params_->cmax_;

        break;
      }

      // derivative of polynomial half cell open circuit potential w.r.t. concentration
      case MAT::PAR::ocp_polynomial:
      {
        for(int i=1; i<params_->ocpparanum_; ++i)
          ocpderiv += i*params_->ocppara_[i]*pow(X,i-1);
        ocpderiv /= params_->cmax_;

        break;
      }

      default:
      {
        dserror("Model for half cell open circuit potential not recognized!");
        break;
      }
    }
  }

  // non-physical intercalation fraction
  else
    ocpderiv = std::numeric_limits<double>::infinity();

  return ocpderiv;
} // MAT::Electrode::ComputeFirstDerivOpenCircuitPotential


/*----------------------------------------------------------------------------------------------------------*
 | compute second derivative of half cell open circuit potential with respect to concentration   fang 08/15 |
 *----------------------------------------------------------------------------------------------------------*/
double MAT::Electrode::ComputeSecondDerivOpenCircuitPotential(
    const double concentration,   //!< concentration
    const double faraday,         //!< Faraday constant
    const double frt              //!< factor F/RT
    ) const
{
  double ocpderiv2(0.);

  // intercalation fraction
  const double X = concentration/params_->cmax_;

  // physically reasonable intercalation fraction
  if(X > 0. and X < 1.)
  {
    switch(params_->ocpmodel_)
    {
      // second derivative of half cell open circuit potential w.r.t. concentration according to Redlich-Kister expansion
      case MAT::PAR::ocp_redlichkister:
      {
        // cf. Colclasure and Kee, Electrochimica Acta 55 (2010) 8960:
        // ocppara_[0]         = DeltaG
        // ocppara_[1,2,3,...] = Redlich-Kister coefficients

        // term not associated with any Redlich-Kister coefficient
        ocpderiv2 = -faraday*(2.*X-1.)/(4.*frt*X*X*(X-1.)*(X-1.));

        // term associated with first Redlich-Kister coefficient vanishes

        // terms associated with second, third, and fourth Redlich-Kister coefficients
        // these three terms are separated from the remaining sum and simplified thereafter to remove singularities in the second derivative of the expansion in case X == 0.5
        ocpderiv2 += 3.*params_->ocppara_[2]+params_->ocppara_[3]*(24.*X-12.)+params_->ocppara_[4]*(120.*X*X-120.*X+27.);

        // terms associated with remaining Redlich-Kister coefficients
        for(int i=4; i<params_->ocpparanum_-1; ++i)
          ocpderiv2 += params_->ocppara_[i+1]*(3.*i*i*pow(2.*X-1.,i-1)+2.*i*(i-1.)*(i-2.)*X*(X-1.)*pow(2.*X-1.,i-3));

        // final scaling
        ocpderiv2 *= 4./(faraday*pow(params_->cmax_,2));

        break;
      }

      // second derivative of half cell open circuit potential w.r.t. concentration according to Taralov, Taralova, Popov, Iliev, Latz, and Zausch (2012)
      case MAT::PAR::ocp_taralov:
      {
        ocpderiv2 = -2.*params_->ocppara_[1]*pow(params_->ocppara_[2],2)/pow(cosh(params_->ocppara_[2]*X+params_->ocppara_[3]),2)*tanh(params_->ocppara_[2]*X+params_->ocppara_[3])
                    +8.*params_->ocppara_[4]*params_->ocppara_[5]*pow(X,6)*exp(params_->ocppara_[5]*pow(X,8))*(7.+8.*params_->ocppara_[5]*pow(X,8))
                    +params_->ocppara_[6]*params_->ocppara_[8]*(params_->ocppara_[8]+1.)/pow(params_->ocppara_[7]-X,params_->ocppara_[8]+2.)
                    +params_->ocppara_[10]*pow(params_->ocppara_[11],2)*exp(params_->ocppara_[11]*(X+params_->ocppara_[12]));
        ocpderiv2 /= pow(params_->cmax_,2);

        break;
      }

      // second derivative of polynomial half cell open circuit potential w.r.t. concentration
      case MAT::PAR::ocp_polynomial:
      {
        for(int i=2; i<params_->ocpparanum_; ++i)
          ocpderiv2 += i*(i-1)*params_->ocppara_[i]*pow(X,i-2);
        ocpderiv2 /= pow(params_->cmax_,2);

        break;
      }

      default:
      {
        dserror("Model for half cell open circuit potential not recognized!");
        break;
      }
    }
  }

  // non-physical intercalation fraction
  else
    ocpderiv2 = std::numeric_limits<double>::infinity();

  return ocpderiv2;
} // MAT::Electrode::ComputeSecondDerivOpenCircuitPotential
