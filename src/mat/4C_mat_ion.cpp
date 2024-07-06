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
Mat::PAR::Ion::Ion(const Core::Mat::PAR::Parameter::Data& matdata)
    : Parameter(matdata),
      valence_(matdata.parameters.get<double>("VALENCE")),
      diffusivity_(matdata.parameters.get<double>("DIFFUSIVITY")),
      densification_(matdata.parameters.get<double>("DENSIFICATION")),
      elimvalence_(matdata.parameters.get<double>("ELIM_VALENCE")),
      elimdiffusivity_(matdata.parameters.get<double>("ELIM_DIFFUSIVITY"))
{
}


Teuchos::RCP<Core::Mat::Material> Mat::PAR::Ion::create_material()
{
  return Teuchos::rcp(new Mat::Ion(this));
}

Mat::IonType Mat::IonType::instance_;


Core::Communication::ParObject* Mat::IonType::create(const std::vector<char>& data)
{
  Mat::Ion* ion = new Mat::Ion();
  ion->unpack(data);
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
void Mat::Ion::pack(Core::Communication::PackBuffer& data) const
{
  Core::Communication::PackBuffer::SizeMarker sm(data);

  // pack type of this instance of ParObject
  int type = unique_par_object_id();
  add_to_pack(data, type);

  // matid
  int matid = -1;
  if (params_ != nullptr) matid = params_->id();  // in case we are in post-process mode
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
void Mat::Ion::unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;

  Core::Communication::ExtractAndAssertId(position, data, unique_par_object_id());

  // matid and recover params_
  int matid;
  extract_from_pack(position, data, matid);
  params_ = nullptr;
  if (Global::Problem::instance()->materials() != Teuchos::null)
    if (Global::Problem::instance()->materials()->num() != 0)
    {
      const int probinst = Global::Problem::instance()->materials()->get_read_from_problem();
      Core::Mat::PAR::Parameter* mat =
          Global::Problem::instance(probinst)->materials()->parameter_by_id(matid);
      if (mat->type() == material_type())
        params_ = static_cast<Mat::PAR::Ion*>(mat);
      else
        FOUR_C_THROW("Type of parameter material %d does not fit to calling type %d", mat->type(),
            material_type());
    }

  if (position != data.size())
    FOUR_C_THROW("Mismatch in size of data %d <-> %d", data.size(), position);

  return;
}

FOUR_C_NAMESPACE_CLOSE
