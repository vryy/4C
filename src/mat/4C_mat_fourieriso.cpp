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
Mat::PAR::FourierIso::FourierIso(const Core::Mat::PAR::Parameter::Data& matdata)
    : Parameter(matdata),
      // be careful: capa_ := rho * C_V, e.g contains the density
      capa_(matdata.parameters.get<double>("CAPA")),
      conduct_(matdata.parameters.get<double>("CONDUCT"))
{
}

Teuchos::RCP<Core::Mat::Material> Mat::PAR::FourierIso::create_material()
{
  return Teuchos::rcp(new Mat::FourierIso(this));
}


Mat::FourierIsoType Mat::FourierIsoType::instance_;


Core::Communication::ParObject* Mat::FourierIsoType::create(const std::vector<char>& data)
{
  Mat::FourierIso* fourieriso = new Mat::FourierIso();
  fourieriso->unpack(data);
  return fourieriso;
}


/*----------------------------------------------------------------------*
 |  Constructor                                   (public)  bborn 04/09 |
 *----------------------------------------------------------------------*/
Mat::FourierIso::FourierIso() : params_(nullptr) {}

/*----------------------------------------------------------------------*
 |  Constructor                                  (public)   bborn 04/09 |
 *----------------------------------------------------------------------*/
Mat::FourierIso::FourierIso(Mat::PAR::FourierIso* params) : params_(params) {}

/*----------------------------------------------------------------------*
 |  Pack                                          (public)  bborn 04/09 |
 *----------------------------------------------------------------------*/
void Mat::FourierIso::pack(Core::Communication::PackBuffer& data) const
{
  Core::Communication::PackBuffer::SizeMarker sm(data);

  // pack type of this instance of ParObject
  int type = unique_par_object_id();
  add_to_pack(data, type);

  // matid
  int matid = -1;
  if (params_ != nullptr) matid = params_->id();  // in case we are in post-process mode
  add_to_pack(data, matid);
}

/*----------------------------------------------------------------------*
 |  Unpack                                        (public)  bborn 04/09 |
 *----------------------------------------------------------------------*/
void Mat::FourierIso::unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;

  Core::Communication::ExtractAndAssertId(position, data, unique_par_object_id());

  // matid
  int matid;
  extract_from_pack(position, data, matid);
  params_ = nullptr;
  if (Global::Problem::instance()->materials() != Teuchos::null)
    if (Global::Problem::instance()->materials()->num() != 0)
    {
      const int probinst = Global::Problem::instance()->materials()->get_read_from_problem();
      Core::Mat::PAR::Parameter* mat =
          Global::Problem::instance(probinst)->materials()->parameter_by_id(matid);
      if (mat->type() == material_type())
        params_ = static_cast<Mat::PAR::FourierIso*>(mat);
      else
        FOUR_C_THROW("Type of parameter material %d does not fit to calling type %d", mat->type(),
            material_type());
    }

  if (position != data.size())
    FOUR_C_THROW("Mismatch in size of data %d <-> %d", data.size(), position);
}

/*----------------------------------------------------------------------*
 |  calculate for 1D                                         dano 09/09 |
 *----------------------------------------------------------------------*/
void Mat::FourierIso::evaluate(const Core::LinAlg::Matrix<1, 1>& gradtemp,
    Core::LinAlg::Matrix<1, 1>& cmat, Core::LinAlg::Matrix<1, 1>& heatflux) const
{
  // conductivity tensor
  cmat(0, 0) = params_->conduct_;

  // heatflux
  heatflux.multiply_nn(cmat, gradtemp);
}

/*----------------------------------------------------------------------*
 |  calculate for 2D                                         dano 09/09 |
 *----------------------------------------------------------------------*/
void Mat::FourierIso::evaluate(const Core::LinAlg::Matrix<2, 1>& gradtemp,
    Core::LinAlg::Matrix<2, 2>& cmat, Core::LinAlg::Matrix<2, 1>& heatflux) const
{
  // conductivity tensor
  cmat.clear();
  for (int i = 0; i < 2; ++i) cmat(i, i) = params_->conduct_;

  // heatflux
  heatflux.multiply_nn(cmat, gradtemp);
}

/*----------------------------------------------------------------------*
 |  calculate for 3D                                         dano 09/09 |
 *----------------------------------------------------------------------*/
void Mat::FourierIso::evaluate(const Core::LinAlg::Matrix<3, 1>& gradtemp,
    Core::LinAlg::Matrix<3, 3>& cmat, Core::LinAlg::Matrix<3, 1>& heatflux) const
{
  // conductivity tensor
  cmat.clear();
  for (int i = 0; i < 3; ++i) cmat(i, i) = params_->conduct_;

  // heatflux
  heatflux.multiply_nn(cmat, gradtemp);
}

FOUR_C_NAMESPACE_CLOSE
