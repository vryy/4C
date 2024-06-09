/*----------------------------------------------------------------------*/
/*! \file
\brief material for macro-scale elements in multi-scale simulations of scalar transport problems

\level 2

*/
/*----------------------------------------------------------------------*/
#include "4C_mat_scatra_multiscale.hpp"

#include "4C_global_data.hpp"
#include "4C_mat_par_bundle.hpp"

FOUR_C_NAMESPACE_OPEN

/*--------------------------------------------------------------------*
 *--------------------------------------------------------------------*/
Mat::PAR::ScatraMultiScale::ScatraMultiScale(Teuchos::RCP<Core::Mat::PAR::Material> matdata)
    : ScatraMat(matdata),
      ScatraMicroMacroCoupling(matdata),
      porosity_(matdata->Get<double>("POROSITY")),
      tortuosity_(matdata->Get<double>("TORTUOSITY"))
{
  return;
}


/*--------------------------------------------------------------------*
 *--------------------------------------------------------------------*/
Teuchos::RCP<Core::Mat::Material> Mat::PAR::ScatraMultiScale::create_material()
{
  return Teuchos::rcp(new Mat::ScatraMultiScale(this));
}


Mat::ScatraMultiScaleType Mat::ScatraMultiScaleType::instance_;


Core::Communication::ParObject* Mat::ScatraMultiScaleType::Create(const std::vector<char>& data)
{
  Mat::ScatraMultiScale* ScatraMatMultiScale = new Mat::ScatraMultiScale();
  ScatraMatMultiScale->Unpack(data);
  return ScatraMatMultiScale;
}


/*--------------------------------------------------------------------*
 *--------------------------------------------------------------------*/
Mat::ScatraMultiScale::ScatraMultiScale() : params_(nullptr) { return; }


/*--------------------------------------------------------------------*
 *--------------------------------------------------------------------*/
Mat::ScatraMultiScale::ScatraMultiScale(Mat::PAR::ScatraMultiScale* params)
    : ScatraMat(params), params_(params)
{
  return;
}


/*--------------------------------------------------------------------*
 *--------------------------------------------------------------------*/
void Mat::ScatraMultiScale::Pack(Core::Communication::PackBuffer& data) const
{
  Core::Communication::PackBuffer::SizeMarker sm(data);
  sm.Insert();

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  add_to_pack(data, type);

  int matid = -1;
  if (params_ != nullptr) matid = params_->Id();  // in case we are in post-process mode
  add_to_pack(data, matid);

  // pack base class material
  ScatraMat::Pack(data);

  return;
}


/*--------------------------------------------------------------------*
 *--------------------------------------------------------------------*/
void Mat::ScatraMultiScale::Unpack(const std::vector<char>& data)
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
        params_ = static_cast<Mat::PAR::ScatraMultiScale*>(mat);
      else
        FOUR_C_THROW("Type of parameter material %d does not match calling type %d!", mat->Type(),
            MaterialType());
    }

  // extract base class material
  std::vector<char> basedata(0);
  extract_from_pack(position, data, basedata);
  ScatraMat::Unpack(basedata);

  // final safety check
  if (position != data.size())
    FOUR_C_THROW("Mismatch in size of data %d <-> %d!", data.size(), position);

  return;
}

FOUR_C_NAMESPACE_CLOSE
