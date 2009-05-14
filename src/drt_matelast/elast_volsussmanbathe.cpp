/*----------------------------------------------------------------------*/
/*!
\file elast_volsussmanbathe.cpp
\brief


the input line should read
  MAT 1 ELAST_VolSussmanBathe KAPPA 100 

<pre>
Maintainer: Sophie Rausch & Thomas Kloeppel
            {rausch,kloeppel}@lnm.mw.tum.de
            089/289 15255
</pre>
*/

/*----------------------------------------------------------------------*/
/* macros */
#ifdef CCADISCRET

/*----------------------------------------------------------------------*/
/* headers */
#include "elast_volsussmanbathe.H"

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
MAT::ELAST::PAR::VolSussmanBathe::VolSussmanBathe(
  Teuchos::RCP<MAT::PAR::Material> matdata
  )
: Parameter(matdata),
  kappa_(matdata->GetDouble("KAPPA"))
{
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
MAT::ELAST::VolSussmanBathe::VolSussmanBathe()
  : Summand(),
    params_(NULL)
{
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
MAT::ELAST::VolSussmanBathe::VolSussmanBathe(MAT::ELAST::PAR::VolSussmanBathe* params)
  : params_(params)
{
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
/*
void MAT::ELAST::VolSussmanBathe::Pack(std::vector<char>& data) const
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
 *----------------------------------------------------------------------*/
/*
void MAT::ELAST::VolSussmanBathe::Unpack(const std::vector<char>& data)
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
      params_ = static_cast<MAT::ELAST::PAR::VolSussmanBathe*>(mat);
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
void MAT::ELAST::VolSussmanBathe::AddCoefficientsModified(
  bool& havecoefficients,
  LINALG::Matrix<3,1>& gamma,
  LINALG::Matrix<5,1>& delta,
  const LINALG::Matrix<3,1>& modinv
  )
{
  havecoefficients = havecoefficients or true;

  const double kappa = params_ -> kappa_;
  
  gamma(2) += kappa*(modinv(2)-1);
  delta(4) += gamma(2) + modinv(2)*kappa;
  
  return;
}


/*----------------------------------------------------------------------*/
#endif // CCADISCRET
