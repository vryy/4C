/*----------------------------------------------------------------------*/
/*! \file
\brief Material law for elastic spring (either translational or rotational spring)

\level 3

*----------------------------------------------------------------------*/


#include "4C_mat_spring.hpp"

#include "4C_global_data.hpp"
#include "4C_mat_par_bundle.hpp"

#include <vector>

FOUR_C_NAMESPACE_OPEN


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Mat::PAR::Spring::Spring(Teuchos::RCP<Core::Mat::PAR::Material> matdata)
    : Parameter(matdata),
      stiffness_(matdata->Get<double>("STIFFNESS")),
      density_(matdata->Get<double>("DENS"))
{
}


Teuchos::RCP<Core::Mat::Material> Mat::PAR::Spring::create_material()
{
  return Teuchos::rcp(new Mat::Spring(this));
}

Mat::SpringType Mat::SpringType::instance_;


Core::Communication::ParObject* Mat::SpringType::Create(const std::vector<char>& data)
{
  Mat::Spring* spring = new Mat::Spring();
  spring->Unpack(data);
  return spring;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Mat::Spring::Spring() : params_(nullptr) {}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Mat::Spring::Spring(Mat::PAR::Spring* params) : params_(params) {}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void Mat::Spring::Pack(Core::Communication::PackBuffer& data) const
{
  Core::Communication::PackBuffer::SizeMarker sm(data);

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  add_to_pack(data, type);

  // matid
  int matid = -1;
  if (params_ != nullptr) matid = params_->Id();  // in case we are in post-process mode
  add_to_pack(data, matid);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void Mat::Spring::Unpack(const std::vector<char>& data)
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
        params_ = static_cast<Mat::PAR::Spring*>(mat);
      else
        FOUR_C_THROW("Type of parameter material %d does not fit to calling type %d", mat->Type(),
            MaterialType());
    }

  if (position != data.size())
    FOUR_C_THROW("Mismatch in size of data %d <-> %d", data.size(), position);
}

FOUR_C_NAMESPACE_CLOSE
