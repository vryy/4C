/*----------------------------------------------------------------------*/
/*! \file

\brief Base of four-element Maxwell material model for reduced dimensional
acinus elements

Four-element Maxwell model consists of a parallel configuration of a spring (Stiffness1),
spring-dashpot (Stiffness2 and Viscosity1) and dashpot (Viscosity2) element
(derivation: see Ismail Mahmoud's dissertation, chapter 3.4)

Input line reads:
(material section)
MAT 3 MAT_0D_MAXWELL_ACINUS_OGDEN Stiffness1 1.0 Stiffness2 5249.1 Viscosity1 3221.86 Viscosity2
1000.0 // acinus properties;


\level 3
*/
/*----------------------------------------------------------------------*/


#include "baci_mat_maxwell_0d_acinus.hpp"

#include "baci_global_data.hpp"
#include "baci_mat_par_bundle.hpp"

#include <vector>

FOUR_C_NAMESPACE_OPEN



/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::PAR::Maxwell0dAcinus::Maxwell0dAcinus(Teuchos::RCP<MAT::PAR::Material> matdata)
    : Parameter(matdata),
      stiffness1_(*matdata->Get<double>("Stiffness1")),
      stiffness2_(*matdata->Get<double>("Stiffness2")),
      viscosity1_(*matdata->Get<double>("Viscosity1")),
      viscosity2_(*matdata->Get<double>("Viscosity2"))
{
}

Teuchos::RCP<MAT::Material> MAT::PAR::Maxwell0dAcinus::CreateMaterial()
{
  return Teuchos::rcp(new MAT::Maxwell0dAcinus(this));
}


MAT::Maxwell0dAcinusType MAT::Maxwell0dAcinusType::instance_;


CORE::COMM::ParObject* MAT::Maxwell0dAcinusType::Create(const std::vector<char>& data)
{
  MAT::Maxwell0dAcinus* mxwll_0d_acin = new MAT::Maxwell0dAcinus();
  mxwll_0d_acin->Unpack(data);
  return mxwll_0d_acin;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::Maxwell0dAcinus::Maxwell0dAcinus() : params_(nullptr) {}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::Maxwell0dAcinus::Maxwell0dAcinus(MAT::PAR::Maxwell0dAcinus* params) : params_(params) {}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::Maxwell0dAcinus::Pack(CORE::COMM::PackBuffer& data) const
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


/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void MAT::Maxwell0dAcinus::Unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;

  CORE::COMM::ExtractAndAssertId(position, data, UniqueParObjectId());

  // matid
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
        params_ = static_cast<MAT::PAR::Maxwell0dAcinus*>(mat);
      else
        FOUR_C_THROW("Type of parameter material %d does not fit to calling type %d", mat->Type(),
            MaterialType());
    }

  if (position != data.size())
    FOUR_C_THROW("Mismatch in size of data %d <-> %d", data.size(), position);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
double MAT::Maxwell0dAcinus::GetParams(std::string parametername)
{
  FOUR_C_THROW("GetParams not implemented yet for this material!");
  return 0;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::Maxwell0dAcinus::SetParams(std::string parametername, double new_value)
{
  FOUR_C_THROW("SetParams not implemented yet for this material!");
}

FOUR_C_NAMESPACE_CLOSE
