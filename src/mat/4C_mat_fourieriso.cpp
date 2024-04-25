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
#include "4C_mat_fourieriso.hpp"

#include "4C_global_data.hpp"
#include "4C_mat_par_bundle.hpp"

FOUR_C_NAMESPACE_OPEN


/*----------------------------------------------------------------------*
 |                                                                      |
 *----------------------------------------------------------------------*/
MAT::PAR::FourierIso::FourierIso(Teuchos::RCP<MAT::PAR::Material> matdata)
    : Parameter(matdata),
      // be careful: capa_ := rho * C_V, e.g contains the density
      capa_(*matdata->Get<double>("CAPA")),
      conduct_(*matdata->Get<double>("CONDUCT"))
{
}

Teuchos::RCP<MAT::Material> MAT::PAR::FourierIso::CreateMaterial()
{
  return Teuchos::rcp(new MAT::FourierIso(this));
}


MAT::FourierIsoType MAT::FourierIsoType::instance_;


CORE::COMM::ParObject* MAT::FourierIsoType::Create(const std::vector<char>& data)
{
  MAT::FourierIso* fourieriso = new MAT::FourierIso();
  fourieriso->Unpack(data);
  return fourieriso;
}


/*----------------------------------------------------------------------*
 |  Constructor                                   (public)  bborn 04/09 |
 *----------------------------------------------------------------------*/
MAT::FourierIso::FourierIso() : params_(nullptr) {}

/*----------------------------------------------------------------------*
 |  Constructor                                  (public)   bborn 04/09 |
 *----------------------------------------------------------------------*/
MAT::FourierIso::FourierIso(MAT::PAR::FourierIso* params) : params_(params) {}

/*----------------------------------------------------------------------*
 |  Pack                                          (public)  bborn 04/09 |
 *----------------------------------------------------------------------*/
void MAT::FourierIso::Pack(CORE::COMM::PackBuffer& data) const
{
  CORE::COMM::PackBuffer::SizeMarker sm(data);
  sm.Insert();

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data, type);

  // matid
  int matid = -1;
  if (params_ != nullptr) matid = params_->Id();  // in case we are in post-process mode
  AddtoPack(data, matid);
}

/*----------------------------------------------------------------------*
 |  Unpack                                        (public)  bborn 04/09 |
 *----------------------------------------------------------------------*/
void MAT::FourierIso::Unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;

  CORE::COMM::ExtractAndAssertId(position, data, UniqueParObjectId());

  // matid
  int matid;
  ExtractfromPack(position, data, matid);
  params_ = nullptr;
  if (GLOBAL::Problem::Instance()->Materials() != Teuchos::null)
    if (GLOBAL::Problem::Instance()->Materials()->Num() != 0)
    {
      const int probinst = GLOBAL::Problem::Instance()->Materials()->GetReadFromProblem();
      MAT::PAR::Parameter* mat =
          GLOBAL::Problem::Instance(probinst)->Materials()->ParameterById(matid);
      if (mat->Type() == MaterialType())
        params_ = static_cast<MAT::PAR::FourierIso*>(mat);
      else
        FOUR_C_THROW("Type of parameter material %d does not fit to calling type %d", mat->Type(),
            MaterialType());
    }

  if (position != data.size())
    FOUR_C_THROW("Mismatch in size of data %d <-> %d", data.size(), position);
}

/*----------------------------------------------------------------------*
 |  calculate for 1D                                         dano 09/09 |
 *----------------------------------------------------------------------*/
void MAT::FourierIso::Evaluate(const CORE::LINALG::Matrix<1, 1>& gradtemp,
    CORE::LINALG::Matrix<1, 1>& cmat, CORE::LINALG::Matrix<1, 1>& heatflux) const
{
  // conductivity tensor
  cmat(0, 0) = params_->conduct_;

  // heatflux
  heatflux.MultiplyNN(cmat, gradtemp);
}

/*----------------------------------------------------------------------*
 |  calculate for 2D                                         dano 09/09 |
 *----------------------------------------------------------------------*/
void MAT::FourierIso::Evaluate(const CORE::LINALG::Matrix<2, 1>& gradtemp,
    CORE::LINALG::Matrix<2, 2>& cmat, CORE::LINALG::Matrix<2, 1>& heatflux) const
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
void MAT::FourierIso::Evaluate(const CORE::LINALG::Matrix<3, 1>& gradtemp,
    CORE::LINALG::Matrix<3, 3>& cmat, CORE::LINALG::Matrix<3, 1>& heatflux) const
{
  // conductivity tensor
  cmat.Clear();
  for (int i = 0; i < 3; ++i) cmat(i, i) = params_->conduct_;

  // heatflux
  heatflux.MultiplyNN(cmat, gradtemp);
}

FOUR_C_NAMESPACE_CLOSE
