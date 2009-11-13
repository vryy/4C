/*----------------------------------------------------------------------*/
/*!
\file fourieriso.cpp
\brief


the input line should read
  MAT 1 MAT_Struct_FourierIso YOUNG 1.044E7 NUE 0.3 DENS 1.0

<pre>
Maintainer: Burkhard Bornemann
            bornemann@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15237
</pre>
*/

/*----------------------------------------------------------------------*
 |  definitions                                              dano 09/09 |
 *----------------------------------------------------------------------*/
#ifdef CCADISCRET

/*----------------------------------------------------------------------*
 |  headers                                                  dano 09/09 |
 *----------------------------------------------------------------------*/
#include "fourieriso.H"

/*----------------------------------------------------------------------*
 |                                                                      |
 *----------------------------------------------------------------------*/
MAT::PAR::FourierIso::FourierIso(
  Teuchos::RCP<MAT::PAR::Material> matdata
  )
: Parameter(matdata),
  capa_(matdata->GetDouble("CAPA")),
  conduct_(matdata->GetDouble("CONDUCT"))
{
}

/*----------------------------------------------------------------------*
 |  Constructor                                   (public)  bborn 04/09 |
 *----------------------------------------------------------------------*/
MAT::FourierIso::FourierIso()
  : params_(NULL)
{
}

/*----------------------------------------------------------------------*
 |  Constructor                                  (public)   bborn 04/09 |
 *----------------------------------------------------------------------*/
MAT::FourierIso::FourierIso(MAT::PAR::FourierIso* params)
  : params_(params)
{
}

/*----------------------------------------------------------------------*
 |  Pack                                          (public)  bborn 04/09 |
 *----------------------------------------------------------------------*/
void MAT::FourierIso::Pack(std::vector<char>& data) const
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

/*----------------------------------------------------------------------*
 |  Unpack                                        (public)  bborn 04/09 |
 *----------------------------------------------------------------------*/
void MAT::FourierIso::Unpack(const std::vector<char>& data)
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
      params_ = static_cast<MAT::PAR::FourierIso*>(mat);
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

/*----------------------------------------------------------------------*
 |  calculate for 1D                                         dano 09/09 |
 *----------------------------------------------------------------------*/
void MAT::FourierIso::Evaluate(
  const LINALG::Matrix<1,1>& gradtemp,
  LINALG::Matrix<1,1>& cmat,
  LINALG::Matrix<1,1>& heatflux
  ) const
{
  // conductivity tensor
  cmat(0,0) = params_->conduct_;

  // heatflux
  heatflux.MultiplyNN(cmat,gradtemp);

  // done
  return;
}

/*----------------------------------------------------------------------*
 |  calculate for 2D                                         dano 09/09 |
 *----------------------------------------------------------------------*/
void MAT::FourierIso::Evaluate(
  const LINALG::Matrix<2,1>& gradtemp,
  LINALG::Matrix<2,2>& cmat,
  LINALG::Matrix<2,1>& heatflux
  ) const
{
  // conductivity tensor
  cmat.Clear();
  for (int i=0; i<2; ++i) cmat(i,i) = params_->conduct_;

  // heatflux
  heatflux.MultiplyNN(cmat,gradtemp);

  // done
  return;
}

/*----------------------------------------------------------------------*
 |  calculate for 3D                                         dano 09/09 |
 *----------------------------------------------------------------------*/
void MAT::FourierIso::Evaluate(
  const LINALG::Matrix<3,1>& gradtemp,
  LINALG::Matrix<3,3>& cmat,
  LINALG::Matrix<3,1>& heatflux
  ) const
{
  // conductivity tensor
  cmat.Clear();
  for (int i=0; i<3; ++i) cmat(i,i) = params_->conduct_;

  // heatflux
  heatflux.MultiplyNN(cmat,gradtemp);

  // done
  return;
}

/*----------------------------------------------------------------------*/
#endif // CCADISCRET
