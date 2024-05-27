/*----------------------------------------------------------------------------*/
/*! \file
\brief material stores parameters for ion species in electrolyte solution. The newman material is
derived for a binary electrolyte using the electroneutrality condition to condense the non-reacting
species

\level 2


*/
/*----------------------------------------------------------------------------*/

#include "4C_mat_newman.hpp"

#include "4C_global_data.hpp"
#include "4C_mat_par_bundle.hpp"
#include "4C_utils_function_of_time.hpp"

#include <vector>

FOUR_C_NAMESPACE_OPEN

// TODO: math.H was included automatically

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::PAR::Newman::Newman(Teuchos::RCP<CORE::MAT::PAR::Material> matdata)
    : ElchSingleMat(matdata),
      valence_(matdata->Get<double>("VALENCE")),
      transnrcurve_(matdata->Get<int>("TRANSNR")),
      thermfaccurve_(matdata->Get<int>("THERMFAC")),
      transnrparanum_(matdata->Get<int>("TRANS_PARA_NUM")),
      transnrpara_(matdata->Get<std::vector<double>>("TRANS_PARA")),
      thermfacparanum_(matdata->Get<int>("THERM_PARA_NUM")),
      thermfacpara_(matdata->Get<std::vector<double>>("THERM_PARA"))
{
  if (transnrparanum_ != (int)transnrpara_.size())
    FOUR_C_THROW("number of materials %d does not fit to size of material vector %d",
        transnrparanum_, transnrpara_.size());
  if (thermfacparanum_ != (int)thermfacpara_.size())
    FOUR_C_THROW("number of materials %d does not fit to size of material vector %d",
        thermfacparanum_, thermfacpara_.size());

  // check if number of provided parameter is valid for a the chosen predefined function
  check_provided_params(transnrcurve_, transnrpara_);
  check_provided_params(thermfaccurve_, thermfacpara_);
}


Teuchos::RCP<CORE::MAT::Material> MAT::PAR::Newman::create_material()
{
  return Teuchos::rcp(new MAT::Newman(this));
}

MAT::NewmanType MAT::NewmanType::instance_;


CORE::COMM::ParObject* MAT::NewmanType::Create(const std::vector<char>& data)
{
  MAT::Newman* newman = new MAT::Newman();
  newman->Unpack(data);
  return newman;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::Newman::Newman() : params_(nullptr) { return; }


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::Newman::Newman(MAT::PAR::Newman* params) : params_(params) { return; }


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::Newman::Pack(CORE::COMM::PackBuffer& data) const
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

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::Newman::Unpack(const std::vector<char>& data)
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
      CORE::MAT::PAR::Parameter* mat =
          GLOBAL::Problem::Instance(probinst)->Materials()->ParameterById(matid);
      if (mat->Type() == MaterialType())
        params_ = static_cast<MAT::PAR::Newman*>(mat);
      else
        FOUR_C_THROW("Type of parameter material %d does not fit to calling type %d", mat->Type(),
            MaterialType());
    }

  if (position != data.size())
    FOUR_C_THROW("Mismatch in size of data %d <-> %d", data.size(), position);

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
double MAT::Newman::compute_transference_number(const double cint) const
{
  double trans = 0.0;

  if (trans_nr_curve() < 0)
    trans = eval_pre_defined_funct(trans_nr_curve(), cint, trans_nr_params());
  else if (trans_nr_curve() == 0)
    trans = eval_pre_defined_funct(-1, cint, trans_nr_params());
  else
    trans = GLOBAL::Problem::Instance()
                ->FunctionById<CORE::UTILS::FunctionOfTime>(trans_nr_curve() - 1)
                .Evaluate(cint);

  return trans;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
double MAT::Newman::compute_first_deriv_trans(const double cint) const
{
  double firstderiv = 0.0;

  if (trans_nr_curve() < 0)
    firstderiv = eval_first_deriv_pre_defined_funct(trans_nr_curve(), cint, trans_nr_params());
  else if (trans_nr_curve() == 0)
    firstderiv = eval_first_deriv_pre_defined_funct(-1, cint, trans_nr_params());
  else
    firstderiv = GLOBAL::Problem::Instance()
                     ->FunctionById<CORE::UTILS::FunctionOfTime>(trans_nr_curve() - 1)
                     .EvaluateDerivative(cint);

  return firstderiv;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
double MAT::Newman::ComputeThermFac(const double cint) const
{
  double therm = 0.0;

  if (therm_fac_curve() < 0)
    therm = eval_pre_defined_funct(therm_fac_curve(), cint, therm_fac_params());
  else if (therm_fac_curve() == 0)
    // thermodynamic factor has to be one if not defined
    therm = 1.0;
  else
    therm = GLOBAL::Problem::Instance()
                ->FunctionById<CORE::UTILS::FunctionOfTime>(therm_fac_curve() - 1)
                .Evaluate(cint);

  return therm;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
double MAT::Newman::compute_first_deriv_therm_fac(const double cint) const
{
  double firstderiv = 0.0;

  if (therm_fac_curve() < 0)
    firstderiv = eval_first_deriv_pre_defined_funct(therm_fac_curve(), cint, therm_fac_params());
  else if (therm_fac_curve() == 0)
    // thermodynamic factor has to be one if not defined
    // -> first derivative = 0.0
    firstderiv = 0.0;
  else
    firstderiv = GLOBAL::Problem::Instance()
                     ->FunctionById<CORE::UTILS::FunctionOfTime>(therm_fac_curve() - 1)
                     .EvaluateDerivative(cint);

  return firstderiv;
}

FOUR_C_NAMESPACE_CLOSE
