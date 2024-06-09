/*----------------------------------------------------------------------*/
/*! \file
\brief material for heat transport due to Fourier-type thermal conduction and the Soret effect


\level 2
*/
/*----------------------------------------------------------------------*/
#include "4C_mat_soret.hpp"

#include "4C_global_data.hpp"
#include "4C_mat_par_bundle.hpp"

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*
 | constructor                                               fang 06/15 |
 *----------------------------------------------------------------------*/
Mat::PAR::Soret::Soret(Teuchos::RCP<Core::Mat::PAR::Material> matdata)
    : FourierIso(matdata), soretcoefficient_(matdata->Get<double>("SORET"))
{
  return;
}


/*------------------------------------------------------------------*
 | create instance of Soret material                     fang 06/15 |
 *------------------------------------------------------------------*/
Teuchos::RCP<Core::Mat::Material> Mat::PAR::Soret::create_material()
{
  return Teuchos::rcp(new Mat::Soret(this));
}

Mat::SoretType Mat::SoretType::instance_;

Core::Communication::ParObject* Mat::SoretType::Create(const std::vector<char>& data)
{
  Mat::Soret* soret = new Mat::Soret();
  soret->Unpack(data);
  return soret;
}


/*------------------------------------------------------------------*
 | construct empty Soret material                        fang 06/15 |
 *------------------------------------------------------------------*/
Mat::Soret::Soret() : params_(nullptr) { return; }


/*-------------------------------------------------------------------------*
 | construct Soret material with specific material parameters   fang 06/15 |
 *-------------------------------------------------------------------------*/
Mat::Soret::Soret(Mat::PAR::Soret* params) : FourierIso(params), params_(params) { return; }


/*----------------------------------------------------------------------*
 | pack material for communication purposes                  fang 06/15 |
 *----------------------------------------------------------------------*/
void Mat::Soret::Pack(Core::Communication::PackBuffer& data) const
{
  Core::Communication::PackBuffer::SizeMarker sm(data);

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  add_to_pack(data, type);

  int matid = -1;
  if (params_ != nullptr) matid = params_->Id();  // in case we are in post-process mode
  add_to_pack(data, matid);

  // pack base class material
  FourierIso::Pack(data);

  return;
}


/*----------------------------------------------------------------------*
 | unpack data from a char vector                            fang 06/15 |
 *----------------------------------------------------------------------*/
void Mat::Soret::Unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;

  Core::Communication::ExtractAndAssertId(position, data, UniqueParObjectId());

  // matid and recover params_
  int matid;
  extract_from_pack(position, data, matid);
  params_ = nullptr;
  if (Global::Problem::Instance()->Materials() != Teuchos::null)
    if (Global::Problem::Instance()->Materials()->Num() != 0)
    {
      const int probinst = Global::Problem::Instance()->Materials()->GetReadFromProblem();
      Core::Mat::PAR::Parameter* mat =
          Global::Problem::Instance(probinst)->Materials()->ParameterById(matid);
      if (mat->Type() == MaterialType())
        params_ = static_cast<Mat::PAR::Soret*>(mat);
      else
        FOUR_C_THROW("Type of parameter material %d does not match calling type %d!", mat->Type(),
            MaterialType());
    }

  // extract base class material
  std::vector<char> basedata(0);
  extract_from_pack(position, data, basedata);
  FourierIso::Unpack(basedata);

  // final safety check
  if (position != data.size())
    FOUR_C_THROW("Mismatch in size of data %d <-> %d!", data.size(), position);

  return;
}

FOUR_C_NAMESPACE_CLOSE
