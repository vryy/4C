/*----------------------------------------------------------------------*/
/*!
\file elast_isoyeoh.cpp
\brief


the input line should read
  MAT 1 ELAST_IsoYeoh C1 100 C2 0 C3 200 

<pre>
Maintainer: Sophie Rausch & Thomas Kloeppel
            {rausch,kloeppel}@lnm.mw.tum.de
            089/289 15257
</pre>
*/

/*----------------------------------------------------------------------*/
/* macros */
#ifdef CCADISCRET

/*----------------------------------------------------------------------*/
/* headers */
#include "elast_isoyeoh.H"

/*----------------------------------------------------------------------*
 |                                                                      |
 *----------------------------------------------------------------------*/
MAT::ELAST::PAR::IsoYeoh::IsoYeoh(
  Teuchos::RCP<MAT::PAR::Material> matdata
  )
: Parameter(matdata),
  c1_(matdata->GetDouble("C1")),
  c2_(matdata->GetDouble("C2")),
  c3_(matdata->GetDouble("C3"))
{
}


/*----------------------------------------------------------------------*
 |  Constructor                                   (public)  bborn 04/09 |
 *----------------------------------------------------------------------*/
MAT::ELAST::IsoYeoh::IsoYeoh()
  : Summand(),
    params_(NULL)
{
}


/*----------------------------------------------------------------------*
 |  Constructor                             (public)   bborn 04/09 |
 *----------------------------------------------------------------------*/
MAT::ELAST::IsoYeoh::IsoYeoh(MAT::ELAST::PAR::IsoYeoh* params)
  : params_(params)
{
}

/*----------------------------------------------------------------------*
 |  Pack                                          (public)  bborn 04/09 |
 *----------------------------------------------------------------------*/
/*
void MAT::ELAST::IsoYeoh::Pack(std::vector<char>& data) const
{
  data.resize(0);

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data,type);

  // matid
  int matid = -1;
  if (params_ != NULL) matid = params_->Id();  // in case we are in post-process mode
  AddtoPack(data,matid);
}
*/

/*----------------------------------------------------------------------*
 |  Unpack                                        (public)  bborn 04/09 |
 *----------------------------------------------------------------------*/
/*
void MAT::ELAST::IsoYeoh::Unpack(const std::vector<char>& data)
{
  int position = 0;
  // extract type
  int type = 0;
  ExtractfromPack(position,data,type);
  if (type != UniqueParObjectId()) dserror("wrong instance type data");

  // matid
  int matid;
  ExtractfromPack(position,data,matid);
  // in post-process mode we do not have any instance of DRT::Problem
  if (DRT::Problem::NumInstances() > 0)
  {
    const int probinst = DRT::Problem::Instance()->Materials()->GetReadFromProblem();
    MAT::PAR::Parameter* mat = DRT::Problem::Instance(probinst)->Materials()->ParameterById(matid);
    if (mat->Type() == MaterialType())
      params_ = static_cast<MAT::ELAST::PAR::IsoYeoh*>(mat);
    else
      dserror("Type of parameter material %d does not fit to calling type %d", mat->Type(), MaterialType());
  }
  else
  {
    params_ = NULL;
  }

  if (position != (int)data.size())
    dserror("Mismatch in size of data %d <-> %d",(int)data.size(),position);
}
*/

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ELAST::IsoYeoh::AddCoefficientsModified(
  bool& havecoefficients,
  LINALG::Matrix<3,1>& gamma,
  LINALG::Matrix<5,1>& delta,
  const LINALG::Matrix<3,1>& modinv
  )
{
  havecoefficients = havecoefficients or true;

  const double c1 = params_ -> c1_;
  const double c2 = params_ -> c2_;
  const double c3 = params_ -> c2_;
  
  gamma(0) += 2*c1+4*c2*(modinv(0)-3)+6*c3*(modinv(0)-3)*(modinv(0)-3);
  
  delta(0) += 8*(c2+3*c3*(modinv(0)-3));
  
  return;
}


/*----------------------------------------------------------------------*/
#endif // CCADISCRET
