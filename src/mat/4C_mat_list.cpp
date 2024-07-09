/*----------------------------------------------------------------------------*/
/*! \file
\brief generic material that stores a list of materials, where each material itself defines the
properties of e.g. one species in a scalar transport problem, or one phase in a fluid problem

\level 1


*/
/*----------------------------------------------------------------------------*/

#include "4C_mat_list.hpp"

#include "4C_global_data.hpp"
#include "4C_mat_par_bundle.hpp"

#include <vector>

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Mat::PAR::MatList::MatList(const Core::Mat::PAR::Parameter::Data& matdata)
    : Parameter(matdata),
      nummat_(matdata.parameters.get<int>("NUMMAT")),
      matids_(matdata.parameters.get<std::vector<int>>("MATIDS")),
      local_((matdata.parameters.get<bool>("LOCAL")))
{
  // check if sizes fit
  if (nummat_ != (int)matids_.size())
    FOUR_C_THROW("number of materials %d does not fit to size of material vector %d", nummat_,
        matids_.size());

  if (not local_)
  {
    // make sure the referenced materials in material list have quick access parameters
    std::vector<int>::const_iterator m;
    for (m = matids_.begin(); m != matids_.end(); ++m)
    {
      const int matid = *m;
      Teuchos::RCP<Core::Mat::Material> mat = Mat::Factory(matid);
      mat_.insert(std::pair<int, Teuchos::RCP<Core::Mat::Material>>(matid, mat));
    }
  }
}

Teuchos::RCP<Core::Mat::Material> Mat::PAR::MatList::create_material()
{
  return Teuchos::rcp(new Mat::MatList(this));
}

Teuchos::RCP<Core::Mat::Material> Mat::PAR::MatList::material_by_id(const int id) const
{
  if (not local_)
  {
    std::map<int, Teuchos::RCP<Core::Mat::Material>>::const_iterator m = mat_.find(id);

    if (m == mat_.end())
    {
      FOUR_C_THROW("Material %d could not be found", id);
      return Teuchos::null;
    }
    else
      return m->second;
  }
  else
    FOUR_C_THROW("This is not allowed");

  return Teuchos::null;
}


Mat::MatListType Mat::MatListType::instance_;


Core::Communication::ParObject* Mat::MatListType::create(const std::vector<char>& data)
{
  Mat::MatList* matlist = new Mat::MatList();
  matlist->unpack(data);
  return matlist;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Mat::MatList::MatList() : params_(nullptr) {}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Mat::MatList::MatList(Mat::PAR::MatList* params) : params_(params)
{
  // setup of material map
  if (params_->local_)
  {
    setup_mat_map();
  }
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void Mat::MatList::setup_mat_map()
{
  // safety first
  mat_.clear();
  if (not mat_.empty()) FOUR_C_THROW("What's going wrong here?");

  // make sure the referenced materials in material list have quick access parameters

  // here's the recursive creation of materials
  std::vector<int>::const_iterator m;
  for (m = params_->mat_ids()->begin(); m != params_->mat_ids()->end(); ++m)
  {
    const int matid = *m;
    Teuchos::RCP<Core::Mat::Material> mat = Mat::Factory(matid);
    if (mat == Teuchos::null) FOUR_C_THROW("Failed to allocate this material");
    mat_.insert(std::pair<int, Teuchos::RCP<Core::Mat::Material>>(matid, mat));
  }
  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void Mat::MatList::clear()
{
  params_ = nullptr;
  mat_.clear();
  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void Mat::MatList::pack(Core::Communication::PackBuffer& data) const
{
  Core::Communication::PackBuffer::SizeMarker sm(data);

  // pack type of this instance of ParObject
  int type = unique_par_object_id();
  add_to_pack(data, type);
  // matid
  int matid = -1;
  if (params_ != nullptr) matid = params_->id();  // in case we are in post-process mode
  add_to_pack(data, matid);

  if (params_ != nullptr)
    if (params_->local_)
    {
      // loop map of associated local materials
      if (params_ != nullptr)
      {
        // std::map<int, Teuchos::RCP<Core::Mat::Material> >::const_iterator m;
        std::vector<int>::const_iterator m;
        for (m = params_->mat_ids()->begin(); m != params_->mat_ids()->end(); m++)
        {
          (mat_.find(*m))->second->pack(data);
        }
      }
    }
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void Mat::MatList::unpack(const std::vector<char>& data)
{
  // make sure we have a pristine material
  clear();

  std::vector<char>::size_type position = 0;

  Core::Communication::ExtractAndAssertId(position, data, unique_par_object_id());

  // matid and recover params_
  int matid(-1);
  extract_from_pack(position, data, matid);
  params_ = nullptr;
  if (Global::Problem::instance()->materials() != Teuchos::null)
    if (Global::Problem::instance()->materials()->num() != 0)
    {
      const int probinst = Global::Problem::instance()->materials()->get_read_from_problem();
      Core::Mat::PAR::Parameter* mat =
          Global::Problem::instance(probinst)->materials()->parameter_by_id(matid);
      if (mat->type() == material_type())
        params_ = static_cast<Mat::PAR::MatList*>(mat);
      else
        FOUR_C_THROW("Type of parameter material %d does not fit to calling type %d", mat->type(),
            material_type());
    }

  if (params_ != nullptr)  // params_ are not accessible in postprocessing mode
  {
    // make sure the referenced materials in material list have quick access parameters
    std::vector<int>::const_iterator m;
    for (m = params_->mat_ids()->begin(); m != params_->mat_ids()->end(); m++)
    {
      const int actmatid = *m;
      Teuchos::RCP<Core::Mat::Material> mat = Mat::Factory(actmatid);
      if (mat == Teuchos::null) FOUR_C_THROW("Failed to allocate this material");
      mat_.insert(std::pair<int, Teuchos::RCP<Core::Mat::Material>>(actmatid, mat));
    }

    if (params_->local_)
    {
      // loop map of associated local materials
      for (m = params_->mat_ids()->begin(); m != params_->mat_ids()->end(); m++)
      {
        std::vector<char> pbtest;
        extract_from_pack(position, data, pbtest);
        (mat_.find(*m))->second->unpack(pbtest);
      }
    }
    // in the postprocessing mode, we do not unpack everything we have packed
    // -> position check cannot be done in this case
    if (position != data.size())
      FOUR_C_THROW("Mismatch in size of data %d <-> %d", data.size(), position);
  }  // if (params_ != nullptr)
}

/*----------------------------------------------------------------------*
 | material ID by Index                                      thon 11/14 |
 *----------------------------------------------------------------------*/
int Mat::MatList::mat_id(const unsigned index) const
{
  if ((int)index < params_->nummat_)
    return params_->matids_.at(index);
  else
  {
    FOUR_C_THROW("Index too large");
    return -1;
  }
}

/*----------------------------------------------------------------------*
 | provide access to material by its ID                      thon 11/14 |
 *----------------------------------------------------------------------*/
///
Teuchos::RCP<Core::Mat::Material> Mat::MatList::material_by_id(const int id) const
{
  if (params_->local_)
  {
    std::map<int, Teuchos::RCP<Core::Mat::Material>>::const_iterator m =
        material_map_read()->find(id);
    if (m == mat_.end())
    {
      FOUR_C_THROW("Material %d could not be found", id);
      return Teuchos::null;
    }
    else
      return m->second;
  }
  else  // material is global (stored in material parameters)
    return params_->material_by_id(id);
}

FOUR_C_NAMESPACE_CLOSE
