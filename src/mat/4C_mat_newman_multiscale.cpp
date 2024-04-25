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
MAT::PAR::NewmanMultiScale::NewmanMultiScale(Teuchos::RCP<MAT::PAR::Material> matdata)
    : Newman(matdata),
      ScatraMicroMacroCoupling(matdata),
      electronic_cond_(*matdata->Get<double>("ELECTRONIC_COND")),
      conc_dep_scale_func_num_(*matdata->Get<int>("ELECTRONIC_COND_CONC_SCALE_FUNC_NUM"))
{
}


/*--------------------------------------------------------------------*
 | create instance of Newman multi-scale material          fang 07/17 |
 *--------------------------------------------------------------------*/
Teuchos::RCP<MAT::Material> MAT::PAR::NewmanMultiScale::CreateMaterial()
{
  return Teuchos::rcp(new MAT::NewmanMultiScale(this));
}


MAT::NewmanMultiScaleType MAT::NewmanMultiScaleType::instance_;


/*--------------------------------------------------------------------*
 | unpack instance of Newman multi-scale material          fang 07/17 |
 *--------------------------------------------------------------------*/
CORE::COMM::ParObject* MAT::NewmanMultiScaleType::Create(const std::vector<char>& data)
{
  MAT::NewmanMultiScale* NewmanMultiScale = new MAT::NewmanMultiScale();
  NewmanMultiScale->Unpack(data);
  return NewmanMultiScale;
}


/*--------------------------------------------------------------------*
 | construct empty Newman multi-scale material             fang 07/17 |
 *--------------------------------------------------------------------*/
MAT::NewmanMultiScale::NewmanMultiScale() : params_(nullptr) {}


/*--------------------------------------------------------------------------------------*
 | construct Newman multi-scale material with specific material parameters   fang 07/17 |
 *--------------------------------------------------------------------------------------*/
MAT::NewmanMultiScale::NewmanMultiScale(MAT::PAR::NewmanMultiScale* params)
    : Newman(params), params_(params)
{
}


/*--------------------------------------------------------------------*
 | pack material for communication purposes                fang 07/17 |
 *--------------------------------------------------------------------*/
void MAT::NewmanMultiScale::Pack(CORE::COMM::PackBuffer& data) const
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
  Newman::Pack(data);
}


/*--------------------------------------------------------------------*
 | unpack data from a char vector                          fang 07/17 |
 *--------------------------------------------------------------------*/
void MAT::NewmanMultiScale::Unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;

  CORE::COMM::ExtractAndAssertId(position, data, UniqueParObjectId());

  // matid and recover params_
  int matid;
  ExtractfromPack(position, data, matid);
  params_ = nullptr;
  if (GLOBAL::Problem::Instance()->Materials() != Teuchos::null)
  {
    if (GLOBAL::Problem::Instance()->Materials()->Num() != 0)
    {
      const int probinst = GLOBAL::Problem::Instance()->Materials()->GetReadFromProblem();
      MAT::PAR::Parameter* mat =
          GLOBAL::Problem::Instance(probinst)->Materials()->ParameterById(matid);
      if (mat->Type() == MaterialType())
        params_ = static_cast<MAT::PAR::NewmanMultiScale*>(mat);
      else
        FOUR_C_THROW("Type of parameter material %d does not match calling type %d!", mat->Type(),
            MaterialType());
    }
  }

  // extract base class material
  std::vector<char> basedata(0);
  ExtractfromPack(position, data, basedata);
  Newman::Unpack(basedata);

  // final safety check
  if (position != data.size())
    FOUR_C_THROW("Mismatch in size of data %d <-> %d!", data.size(), position);
}

/*--------------------------------------------------------------------*
 *--------------------------------------------------------------------*/
double MAT::NewmanMultiScale::electronic_cond(const int gp) const
{
  const int func_num = params_->conc_dep_scale_func_num();
  if (func_num > 0)
  {
    return GLOBAL::Problem::Instance()
               ->FunctionById<CORE::UTILS::FunctionOfAnything>(func_num - 1)
               .Evaluate({{"c", EvaluateMeanConcentration(gp)}}, {}, 0) *
           params_->electronic_cond();
  }
  else
  {
    return params_->electronic_cond();
  }
}

FOUR_C_NAMESPACE_CLOSE
