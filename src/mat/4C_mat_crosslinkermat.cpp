/*-----------------------------------------------------------*/
/*! \file
\brief A class for a crosslinker material


\level 3

*/
/*-----------------------------------------------------------*/

#include "4C_mat_crosslinkermat.hpp"

#include "4C_global_data.hpp"
#include "4C_mat_par_bundle.hpp"

FOUR_C_NAMESPACE_OPEN


/*----------------------------------------------------------------------*
 |                                                                      |
 *----------------------------------------------------------------------*/
Mat::PAR::CrosslinkerMat::CrosslinkerMat(const Core::Mat::PAR::Parameter::Data& matdata)
    : Parameter(matdata),
      link_element_matnum_(matdata.parameters.get<double>("MATNUM")),
      jointtype_(Inpar::BEAMINTERACTION::String2JointType(
          (matdata.parameters.get<std::string>("JOINTTYPE")))),
      linkinglength_(matdata.parameters.get<double>("LINKINGLENGTH")),
      linkinglengthtol_(matdata.parameters.get<double>("LINKINGLENGTHTOL")),
      linkingangle_(matdata.parameters.get<double>("LINKINGANGLE")),
      linkingangletol_(matdata.parameters.get<double>("LINKINGANGLETOL")),
      k_on_(matdata.parameters.get<double>("K_ON")),
      k_off_(matdata.parameters.get<double>("K_OFF")),
      deltabelleq_(matdata.parameters.get<double>("DELTABELLEQ")),
      nobonddistsphere(matdata.parameters.get<double>("NOBONDDISTSPHERE")),
      linkertype_(Inpar::BEAMINTERACTION::String2CrosslinkerType(
          (matdata.parameters.get<std::string>("TYPE"))))
{
  if (link_element_matnum_ < 0)
    FOUR_C_THROW(
        "Material number for underlying linker element of this crosslinker"
        "must be greater than zero");
  if (linkinglength_ < 1e-08)
    FOUR_C_THROW(
        "Linking length (distance of two binding spots of a linker) must be\n"
        "greater than zero (as you need to divide by it during crosslinker diffusion).");
  if (linkinglengthtol_ < 0.0 || linkinglengthtol_ > linkinglength_)
    FOUR_C_THROW(" Value for tolerance of linking does not make sense.");
  if ((linkinglength_ - linkinglengthtol_) < 1e-08)
    FOUR_C_THROW(
        "choose linkinglengthtol < linkinglength_, otherwise a linker with length 0.0 maybe be "
        "possible.");
}

Teuchos::RCP<Core::Mat::Material> Mat::PAR::CrosslinkerMat::create_material()
{
  return Teuchos::rcp(new Mat::CrosslinkerMat(this));
}

Mat::CrosslinkerMatType Mat::CrosslinkerMatType::instance_;


Core::Communication::ParObject* Mat::CrosslinkerMatType::create(const std::vector<char>& data)
{
  Mat::CrosslinkerMat* linkermat = new Mat::CrosslinkerMat();
  linkermat->unpack(data);
  return linkermat;
}


/*----------------------------------------------------------------------*
 |                                                                      |
 *----------------------------------------------------------------------*/
Mat::CrosslinkerMat::CrosslinkerMat() : params_(nullptr) {}


/*----------------------------------------------------------------------*
 |                                                                      |
 *----------------------------------------------------------------------*/
Mat::CrosslinkerMat::CrosslinkerMat(Mat::PAR::CrosslinkerMat* params) : params_(params) {}


/*----------------------------------------------------------------------*
 |                                                                      |
 *----------------------------------------------------------------------*/
void Mat::CrosslinkerMat::pack(Core::Communication::PackBuffer& data) const
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
 |                                                                      |
 *----------------------------------------------------------------------*/
void Mat::CrosslinkerMat::unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;

  Core::Communication::ExtractAndAssertId(position, data, unique_par_object_id());

  // matid and recover params_
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
        params_ = static_cast<Mat::PAR::CrosslinkerMat*>(mat);
      else
        FOUR_C_THROW("Type of parameter material %d does not fit to calling type %d", mat->type(),
            material_type());
    }

  if (position != data.size())
    FOUR_C_THROW("Mismatch in size of data %d <-> %d", data.size(), position);
}

FOUR_C_NAMESPACE_CLOSE
