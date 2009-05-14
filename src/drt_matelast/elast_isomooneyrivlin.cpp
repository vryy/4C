/*----------------------------------------------------------------------*/
/*!
\file elast_isomooneyrivlin.cpp
\brief
This file contains the routines required to calculate the isochoric contribution 
of a Mooney-Rivlin-type material

<pre>
Maintainer: Sophie Rausch & Thomas Kloeppel
            {rausch,kloeppel}@lnm.mw.tum.de
            089/289 15255
</pre>

*----------------------------------------------------------------------*/
/* macros */
#ifdef CCADISCRET

/*----------------------------------------------------------------------*/
/* headers */
#include "elast_isomooneyrivlin.H"

/*----------------------------------------------------------------------*
 |                                                                      |
 *----------------------------------------------------------------------*/
MAT::ELAST::PAR::IsoMooneyRivlin::IsoMooneyRivlin(
  Teuchos::RCP<MAT::PAR::Material> matdata
  )
: Parameter(matdata),
  c1_(matdata->GetDouble("C1")),
  c2_(matdata->GetDouble("C2"))
{
}


/*----------------------------------------------------------------------*
 |  Constructor                                   (public)  bborn 04/09 |
 *----------------------------------------------------------------------*/
MAT::ELAST::IsoMooneyRivlin::IsoMooneyRivlin()
  : Summand(),
    params_(NULL)
{
}


/*----------------------------------------------------------------------*
 |  Constructor                             (public)   bborn 04/09 |
 *----------------------------------------------------------------------*/
MAT::ELAST::IsoMooneyRivlin::IsoMooneyRivlin(MAT::ELAST::PAR::IsoMooneyRivlin* params)
  : params_(params)
{
}

/*----------------------------------------------------------------------*
 |  Pack                                          (public)  bborn 04/09 |
 *----------------------------------------------------------------------*/
/*
void MAT::ELAST::IsoMooneyRivlin::Pack(std::vector<char>& data) const
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
void MAT::ELAST::IsoMooneyRivlin::Unpack(const std::vector<char>& data)
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
      params_ = static_cast<MAT::ELAST::PAR::IsoMooneyRivlin*>(mat);
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
void MAT::ELAST::IsoMooneyRivlin::AddCoefficientsModified(
  bool& havecoefficients,
  LINALG::Matrix<3,1>& gamma,
  LINALG::Matrix<5,1>& delta,
  const LINALG::Matrix<3,1>& modinv
  )
{
  havecoefficients = havecoefficients or true;

  const double c1 = params_ -> c1_;
  const double c2 = params_ -> c2_;
  
  gamma(0) += 2*c1+2*modinv(0)*c2;
  gamma(1) += - 2*c2;
  
  delta(0) += 4*c2;
  delta(3) += - 4*c2;
  
  return;
}


/*----------------------------------------------------------------------*/
#endif // CCADISCRET
