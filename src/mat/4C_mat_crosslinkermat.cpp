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
MAT::PAR::CrosslinkerMat::CrosslinkerMat(Teuchos::RCP<MAT::PAR::Material> matdata)
    : Parameter(matdata),
      link_element_matnum_(matdata->Get<double>("MATNUM")),
      jointtype_(
          INPAR::BEAMINTERACTION::String2JointType((matdata->Get<std::string>("JOINTTYPE")))),
      linkinglength_(matdata->Get<double>("LINKINGLENGTH")),
      linkinglengthtol_(matdata->Get<double>("LINKINGLENGTHTOL")),
      linkingangle_(matdata->Get<double>("LINKINGANGLE")),
      linkingangletol_(matdata->Get<double>("LINKINGANGLETOL")),
      k_on_(matdata->Get<double>("K_ON")),
      k_off_(matdata->Get<double>("K_OFF")),
      deltabelleq_(matdata->Get<double>("DELTABELLEQ")),
      nobonddistsphere(matdata->Get<double>("NOBONDDISTSPHERE")),
      linkertype_(
          INPAR::BEAMINTERACTION::String2CrosslinkerType((matdata->Get<std::string>("TYPE"))))
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

Teuchos::RCP<MAT::Material> MAT::PAR::CrosslinkerMat::CreateMaterial()
{
  return Teuchos::rcp(new MAT::CrosslinkerMat(this));
}

MAT::CrosslinkerMatType MAT::CrosslinkerMatType::instance_;


CORE::COMM::ParObject* MAT::CrosslinkerMatType::Create(const std::vector<char>& data)
{
  MAT::CrosslinkerMat* linkermat = new MAT::CrosslinkerMat();
  linkermat->Unpack(data);
  return linkermat;
}


/*----------------------------------------------------------------------*
 |                                                                      |
 *----------------------------------------------------------------------*/
MAT::CrosslinkerMat::CrosslinkerMat() : params_(nullptr) {}


/*----------------------------------------------------------------------*
 |                                                                      |
 *----------------------------------------------------------------------*/
MAT::CrosslinkerMat::CrosslinkerMat(MAT::PAR::CrosslinkerMat* params) : params_(params) {}


/*----------------------------------------------------------------------*
 |                                                                      |
 *----------------------------------------------------------------------*/
void MAT::CrosslinkerMat::Pack(CORE::COMM::PackBuffer& data) const
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
 |                                                                      |
 *----------------------------------------------------------------------*/
void MAT::CrosslinkerMat::Unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;

  CORE::COMM::ExtractAndAssertId(position, data, UniqueParObjectId());

  // matid and recover params_
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
        params_ = static_cast<MAT::PAR::CrosslinkerMat*>(mat);
      else
        FOUR_C_THROW("Type of parameter material %d does not fit to calling type %d", mat->Type(),
            MaterialType());
    }

  if (position != data.size())
    FOUR_C_THROW("Mismatch in size of data %d <-> %d", data.size(), position);
}

FOUR_C_NAMESPACE_CLOSE
