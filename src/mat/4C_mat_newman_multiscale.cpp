/*----------------------------------------------------------------------*/
/*! \file
\brief material for macro-scale elements in multi-scale simulations of electrochemistry problems

\level 2

*/
/*----------------------------------------------------------------------*/
#include "4C_mat_newman_multiscale.hpp"

#include "4C_global_data.hpp"
#include "4C_mat_par_bundle.hpp"
#include "4C_utils_function.hpp"

FOUR_C_NAMESPACE_OPEN

/*--------------------------------------------------------------------*
 | constructor                                             fang 07/17 |
 *--------------------------------------------------------------------*/
Mat::PAR::NewmanMultiScale::NewmanMultiScale(const Core::Mat::PAR::Parameter::Data& matdata)
    : Newman(matdata),
      ScatraMicroMacroCoupling(matdata),
      electronic_cond_(matdata.parameters.get<double>("ELECTRONIC_COND")),
      conc_dep_scale_func_num_(matdata.parameters.get<int>("ELECTRONIC_COND_CONC_SCALE_FUNC_NUM"))
{
}


/*--------------------------------------------------------------------*
 | create instance of Newman multi-scale material          fang 07/17 |
 *--------------------------------------------------------------------*/
Teuchos::RCP<Core::Mat::Material> Mat::PAR::NewmanMultiScale::create_material()
{
  return Teuchos::rcp(new Mat::NewmanMultiScale(this));
}


Mat::NewmanMultiScaleType Mat::NewmanMultiScaleType::instance_;


/*--------------------------------------------------------------------*
 | unpack instance of Newman multi-scale material          fang 07/17 |
 *--------------------------------------------------------------------*/
Core::Communication::ParObject* Mat::NewmanMultiScaleType::Create(const std::vector<char>& data)
{
  Mat::NewmanMultiScale* NewmanMultiScale = new Mat::NewmanMultiScale();
  NewmanMultiScale->unpack(data);
  return NewmanMultiScale;
}


/*--------------------------------------------------------------------*
 | construct empty Newman multi-scale material             fang 07/17 |
 *--------------------------------------------------------------------*/
Mat::NewmanMultiScale::NewmanMultiScale() : params_(nullptr) {}


/*--------------------------------------------------------------------------------------*
 | construct Newman multi-scale material with specific material parameters   fang 07/17 |
 *--------------------------------------------------------------------------------------*/
Mat::NewmanMultiScale::NewmanMultiScale(Mat::PAR::NewmanMultiScale* params)
    : Newman(params), params_(params)
{
}


/*--------------------------------------------------------------------*
 | pack material for communication purposes                fang 07/17 |
 *--------------------------------------------------------------------*/
void Mat::NewmanMultiScale::pack(Core::Communication::PackBuffer& data) const
{
  Core::Communication::PackBuffer::SizeMarker sm(data);

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  add_to_pack(data, type);

  int matid = -1;
  if (params_ != nullptr) matid = params_->Id();  // in case we are in post-process mode
  add_to_pack(data, matid);

  // pack base class material
  Newman::pack(data);
}


/*--------------------------------------------------------------------*
 | unpack data from a char vector                          fang 07/17 |
 *--------------------------------------------------------------------*/
void Mat::NewmanMultiScale::unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;

  Core::Communication::ExtractAndAssertId(position, data, UniqueParObjectId());

  // matid and recover params_
  int matid;
  extract_from_pack(position, data, matid);
  params_ = nullptr;
  if (Global::Problem::Instance()->Materials() != Teuchos::null)
  {
    if (Global::Problem::Instance()->Materials()->Num() != 0)
    {
      const int probinst = Global::Problem::Instance()->Materials()->GetReadFromProblem();
      Core::Mat::PAR::Parameter* mat =
          Global::Problem::Instance(probinst)->Materials()->ParameterById(matid);
      if (mat->Type() == MaterialType())
        params_ = static_cast<Mat::PAR::NewmanMultiScale*>(mat);
      else
        FOUR_C_THROW("Type of parameter material %d does not match calling type %d!", mat->Type(),
            MaterialType());
    }
  }

  // extract base class material
  std::vector<char> basedata(0);
  extract_from_pack(position, data, basedata);
  Newman::unpack(basedata);

  // final safety check
  if (position != data.size())
    FOUR_C_THROW("Mismatch in size of data %d <-> %d!", data.size(), position);
}

/*--------------------------------------------------------------------*
 *--------------------------------------------------------------------*/
double Mat::NewmanMultiScale::electronic_cond(const int gp) const
{
  const int func_num = params_->conc_dep_scale_func_num();
  if (func_num > 0)
  {
    return Global::Problem::Instance()
               ->FunctionById<Core::UTILS::FunctionOfAnything>(func_num - 1)
               .evaluate({{"c", evaluate_mean_concentration(gp)}}, {}, 0) *
           params_->electronic_cond();
  }
  else
  {
    return params_->electronic_cond();
  }
}

FOUR_C_NAMESPACE_CLOSE
