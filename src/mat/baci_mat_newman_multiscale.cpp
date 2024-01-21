/*----------------------------------------------------------------------*/
/*! \file
\brief material for macro-scale elements in multi-scale simulations of electrochemistry problems

\level 2

*/
/*----------------------------------------------------------------------*/
#include "baci_mat_newman_multiscale.H"

#include "baci_global_data.H"
#include "baci_mat_par_bundle.H"

BACI_NAMESPACE_OPEN

/*--------------------------------------------------------------------*
 | constructor                                             fang 07/17 |
 *--------------------------------------------------------------------*/
MAT::PAR::NewmanMultiScale::NewmanMultiScale(Teuchos::RCP<MAT::PAR::Material> matdata)
    : Newman(matdata), ScatraMultiScale(matdata), sigma_(matdata->GetDouble("SIGMA"))
{
  return;
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
MAT::NewmanMultiScale::NewmanMultiScale() : params_(nullptr) { return; }


/*--------------------------------------------------------------------------------------*
 | construct Newman multi-scale material with specific material parameters   fang 07/17 |
 *--------------------------------------------------------------------------------------*/
MAT::NewmanMultiScale::NewmanMultiScale(MAT::PAR::NewmanMultiScale* params)
    : Newman(params), params_(params)
{
  return;
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

  return;
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
  if (DRT::Problem::Instance()->Materials() != Teuchos::null)
    if (DRT::Problem::Instance()->Materials()->Num() != 0)
    {
      const int probinst = DRT::Problem::Instance()->Materials()->GetReadFromProblem();
      MAT::PAR::Parameter* mat =
          DRT::Problem::Instance(probinst)->Materials()->ParameterById(matid);
      if (mat->Type() == MaterialType())
        params_ = static_cast<MAT::PAR::NewmanMultiScale*>(mat);
      else
        dserror("Type of parameter material %d does not match calling type %d!", mat->Type(),
            MaterialType());
    }

  // extract base class material
  std::vector<char> basedata(0);
  ExtractfromPack(position, data, basedata);
  Newman::Unpack(basedata);

  // final safety check
  if (position != data.size())
    dserror("Mismatch in size of data %d <-> %d!", data.size(), position);

  return;
}

BACI_NAMESPACE_CLOSE
