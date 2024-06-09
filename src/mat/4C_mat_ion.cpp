/*----------------------------------------------------------------------------*/
/*! \file
\brief material that stores parameters of ion species in electrolyte solution. Former file of Georg
Bauer

\level 2


*/
/*----------------------------------------------------------------------------*/

#include "4C_mat_ion.hpp"

#include "4C_global_data.hpp"
#include "4C_mat_par_bundle.hpp"

#include <vector>

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Mat::PAR::Ion::Ion(Teuchos::RCP<Core::Mat::PAR::Material> matdata)
    : Parameter(matdata),
      valence_(matdata->Get<double>("VALENCE")),
      diffusivity_(matdata->Get<double>("DIFFUSIVITY")),
      densification_(matdata->Get<double>("DENSIFICATION")),
      elimvalence_(matdata->Get<double>("ELIM_VALENCE")),
      elimdiffusivity_(matdata->Get<double>("ELIM_DIFFUSIVITY"))
{
}


Teuchos::RCP<Core::Mat::Material> Mat::PAR::Ion::create_material()
{
  return Teuchos::rcp(new Mat::Ion(this));
}

Mat::IonType Mat::IonType::instance_;


Core::Communication::ParObject* Mat::IonType::Create(const std::vector<char>& data)
{
  Mat::Ion* ion = new Mat::Ion();
  ion->Unpack(data);
  return ion;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Mat::Ion::Ion() : params_(nullptr) {}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Mat::Ion::Ion(Mat::PAR::Ion* params) : params_(params) {}



/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void Mat::Ion::Pack(Core::Communication::PackBuffer& data) const
{
  Core::Communication::PackBuffer::SizeMarker sm(data);
  sm.Insert();

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  add_to_pack(data, type);

  // matid
  int matid = -1;
  if (params_ != nullptr) matid = params_->Id();  // in case we are in post-process mode
  add_to_pack(data, matid);

  /*
  for (unsigned i=0;i<data().size();i++)
  std::cout<<"Pack ION: pb["<<i<<"] = "<<(data())[i]<<std::endl;
*/
  /*
  // extract type
  std::vector<char>::size_type posit = 0;
  std::vector<char> pbtest;
  int typio = 0;
  extract_from_pack(posit,data(),typio);
  std::cout<<"ION Pack: Type will be "<<typio<<std::endl;
*/
  // std::cout<<"Ion Pack: "<<data().size()<<std::endl;

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void Mat::Ion::Unpack(const std::vector<char>& data)
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
        params_ = static_cast<Mat::PAR::Ion*>(mat);
      else
        FOUR_C_THROW("Type of parameter material %d does not fit to calling type %d", mat->Type(),
            MaterialType());
    }

  if (position != data.size())
    FOUR_C_THROW("Mismatch in size of data %d <-> %d", data.size(), position);

  return;
}

FOUR_C_NAMESPACE_CLOSE
