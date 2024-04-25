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
MAT::PAR::ScatraMultiScale::ScatraMultiScale(Teuchos::RCP<MAT::PAR::Material> matdata)
    : ScatraMat(matdata),
      ScatraMicroMacroCoupling(matdata),
      porosity_(*matdata->Get<double>("POROSITY")),
      tortuosity_(*matdata->Get<double>("TORTUOSITY"))
{
  return;
}


/*--------------------------------------------------------------------*
 *--------------------------------------------------------------------*/
Teuchos::RCP<MAT::Material> MAT::PAR::ScatraMultiScale::CreateMaterial()
{
  return Teuchos::rcp(new MAT::ScatraMultiScale(this));
}


MAT::ScatraMultiScaleType MAT::ScatraMultiScaleType::instance_;


CORE::COMM::ParObject* MAT::ScatraMultiScaleType::Create(const std::vector<char>& data)
{
  MAT::ScatraMultiScale* ScatraMatMultiScale = new MAT::ScatraMultiScale();
  ScatraMatMultiScale->Unpack(data);
  return ScatraMatMultiScale;
}


/*--------------------------------------------------------------------*
 *--------------------------------------------------------------------*/
MAT::ScatraMultiScale::ScatraMultiScale() : params_(nullptr) { return; }


/*--------------------------------------------------------------------*
 *--------------------------------------------------------------------*/
MAT::ScatraMultiScale::ScatraMultiScale(MAT::PAR::ScatraMultiScale* params)
    : ScatraMat(params), params_(params)
{
  return;
}


/*--------------------------------------------------------------------*
 *--------------------------------------------------------------------*/
void MAT::ScatraMultiScale::Pack(CORE::COMM::PackBuffer& data) const
{
  CORE::COMM::PackBuffer::SizeMarker sm(data);
  sm.Insert();

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data, type);

  int matid = -1;
  if (params_ != nullptr) matid = params_->Id();  // in case we are in post-process mode
  AddtoPack(data, matid);

  // pack base class material
  ScatraMat::Pack(data);

  return;
}


/*--------------------------------------------------------------------*
 *--------------------------------------------------------------------*/
void MAT::ScatraMultiScale::Unpack(const std::vector<char>& data)
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
        params_ = static_cast<MAT::PAR::ScatraMultiScale*>(mat);
      else
        FOUR_C_THROW("Type of parameter material %d does not match calling type %d!", mat->Type(),
            MaterialType());
    }

  // extract base class material
  std::vector<char> basedata(0);
  ExtractfromPack(position, data, basedata);
  ScatraMat::Unpack(basedata);

  // final safety check
  if (position != data.size())
    FOUR_C_THROW("Mismatch in size of data %d <-> %d!", data.size(), position);

  return;
}

FOUR_C_NAMESPACE_CLOSE
