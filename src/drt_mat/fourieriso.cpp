/*----------------------------------------------------------------------*/
/*! \file
\brief heat condution according to fourier's law


the input line should read
  MAT 1   THERM_FourierIso   CAPA 420 CONDUCT 52

\level 1
*/

/*----------------------------------------------------------------------*
 |  headers                                                  dano 09/09 |
 *----------------------------------------------------------------------*/
#include "fourieriso.H"
#include "matpar_bundle.H"

#include "../drt_lib/drt_globalproblem.H"


/*----------------------------------------------------------------------*
 |                                                                      |
 *----------------------------------------------------------------------*/
MAT::PAR::FourierIso::FourierIso(Teuchos::RCP<MAT::PAR::Material> matdata)
    : Parameter(matdata),
      // be careful: capa_ := rho * C_V, e.g contains the density
      capa_(matdata->GetDouble("CAPA")),
      conduct_(matdata->GetDouble("CONDUCT"))
{
}

Teuchos::RCP<MAT::Material> MAT::PAR::FourierIso::CreateMaterial()
{
  return Teuchos::rcp(new MAT::FourierIso(this));
}


MAT::FourierIsoType MAT::FourierIsoType::instance_;


DRT::ParObject* MAT::FourierIsoType::Create(const std::vector<char>& data)
{
  MAT::FourierIso* fourieriso = new MAT::FourierIso();
  fourieriso->Unpack(data);
  return fourieriso;
}


/*----------------------------------------------------------------------*
 |  Constructor                                   (public)  bborn 04/09 |
 *----------------------------------------------------------------------*/
MAT::FourierIso::FourierIso() : params_(NULL) {}

/*----------------------------------------------------------------------*
 |  Constructor                                  (public)   bborn 04/09 |
 *----------------------------------------------------------------------*/
MAT::FourierIso::FourierIso(MAT::PAR::FourierIso* params) : params_(params) {}

/*----------------------------------------------------------------------*
 |  Pack                                          (public)  bborn 04/09 |
 *----------------------------------------------------------------------*/
void MAT::FourierIso::Pack(DRT::PackBuffer& data) const
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

/*----------------------------------------------------------------------*
 |  Unpack                                        (public)  bborn 04/09 |
 *----------------------------------------------------------------------*/
void MAT::FourierIso::Unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;
  // extract type
  int type = 0;
  ExtractfromPack(position, data, type);
  if (type != UniqueParObjectId()) dserror("wrong instance type data");

  // matid
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
        params_ = static_cast<MAT::PAR::FourierIso*>(mat);
      else
        dserror("Type of parameter material %d does not fit to calling type %d", mat->Type(),
            MaterialType());
    }

  if (position != data.size()) dserror("Mismatch in size of data %d <-> %d", data.size(), position);
}

/*----------------------------------------------------------------------*
 |  calculate for 1D                                         dano 09/09 |
 *----------------------------------------------------------------------*/
void MAT::FourierIso::Evaluate(const LINALG::Matrix<1, 1>& gradtemp, LINALG::Matrix<1, 1>& cmat,
    LINALG::Matrix<1, 1>& heatflux) const
{
  // conductivity tensor
  cmat(0, 0) = params_->conduct_;

  // heatflux
  heatflux.MultiplyNN(cmat, gradtemp);
}

/*----------------------------------------------------------------------*
 |  calculate for 2D                                         dano 09/09 |
 *----------------------------------------------------------------------*/
void MAT::FourierIso::Evaluate(const LINALG::Matrix<2, 1>& gradtemp, LINALG::Matrix<2, 2>& cmat,
    LINALG::Matrix<2, 1>& heatflux) const
{
  // conductivity tensor
  cmat.Clear();
  for (int i = 0; i < 2; ++i) cmat(i, i) = params_->conduct_;

  // heatflux
  heatflux.MultiplyNN(cmat, gradtemp);
}

/*----------------------------------------------------------------------*
 |  calculate for 3D                                         dano 09/09 |
 *----------------------------------------------------------------------*/
void MAT::FourierIso::Evaluate(const LINALG::Matrix<3, 1>& gradtemp, LINALG::Matrix<3, 3>& cmat,
    LINALG::Matrix<3, 1>& heatflux) const
{
  // conductivity tensor
  cmat.Clear();
  for (int i = 0; i < 3; ++i) cmat(i, i) = params_->conduct_;

  // heatflux
  heatflux.MultiplyNN(cmat, gradtemp);
}
