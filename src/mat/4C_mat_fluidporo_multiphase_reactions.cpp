/*----------------------------------------------------------------------*/
/*! \file
 \brief a fluid material for porous multiphase flow with reactions (mass sources and sinks)

   \level 3

 *----------------------------------------------------------------------*/


#include "4C_mat_fluidporo_multiphase_reactions.hpp"

#include "4C_global_data.hpp"
#include "4C_mat_fluidporo_multiphase_singlereaction.hpp"
#include "4C_mat_par_bundle.hpp"

#include <vector>

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*
 | rstandard constructor                                     vuong 08/16 |
 *----------------------------------------------------------------------*/
Mat::PAR::FluidPoroMultiPhaseReactions::FluidPoroMultiPhaseReactions(
    const Core::Mat::PAR::Parameter::Data& matdata)
    : FluidPoroMultiPhase(matdata),
      numreac_((matdata.parameters.Get<int>("NUMREAC"))),
      reacids_((matdata.parameters.Get<std::vector<int>>("REACIDS")))
{
  // check if sizes fit
  if (numreac_ != (int)reacids_.size())
    FOUR_C_THROW("number of materials %d does not fit to size of material vector %d", nummat_,
        reacids_.size());

  if (numreac_ < 1)
    FOUR_C_THROW("if you don't have reactions, use MAT_matlist instead of MAT_matlist_reactions!");

  if (not local_)
  {
    // make sure the referenced materials in material list have quick access parameters
    std::vector<int>::const_iterator m;
    for (m = reacids_.begin(); m != reacids_.end(); ++m)
    {
      const int reacid = *m;
      Teuchos::RCP<Core::Mat::Material> mat = Mat::Factory(reacid);

      // safety check and cast
      if (mat->MaterialType() != Core::Materials::m_fluidporo_singlereaction)
        FOUR_C_THROW("only MAT_FluidPoroSingleReaction material valid");
      Mat::FluidPoroSingleReaction singlereacmat = static_cast<Mat::FluidPoroSingleReaction&>(*mat);
      if (singlereacmat.TotalNumDof() != this->nummat_)
        FOUR_C_THROW(
            "TOTALNUMDOF in MAT_FluidPoroSingleReaction does not correspond to NUMMAT in "
            "MAT_FluidPoroMultiPhaseReactions");

      material_map_write()->insert(std::pair<int, Teuchos::RCP<Core::Mat::Material>>(reacid, mat));
    }
  }
}

Teuchos::RCP<Core::Mat::Material> Mat::PAR::FluidPoroMultiPhaseReactions::create_material()
{
  return Teuchos::rcp(new Mat::FluidPoroMultiPhaseReactions(this));
}


Mat::FluidPoroMultiPhaseReactionsType Mat::FluidPoroMultiPhaseReactionsType::instance_;


Core::Communication::ParObject* Mat::FluidPoroMultiPhaseReactionsType::Create(
    const std::vector<char>& data)
{
  Mat::FluidPoroMultiPhaseReactions* FluidPoroMultiPhaseReactions =
      new Mat::FluidPoroMultiPhaseReactions();
  FluidPoroMultiPhaseReactions->Unpack(data);
  return FluidPoroMultiPhaseReactions;
}


/*----------------------------------------------------------------------*
 | construct empty material object                           vuong 08/16 |
 *----------------------------------------------------------------------*/
Mat::FluidPoroMultiPhaseReactions::FluidPoroMultiPhaseReactions()
    : FluidPoroMultiPhase(), paramsreac_(nullptr)
{
}

/*----------------------------------------------------------------------*
 | construct the material object given material paramete     vuong 08/16 |
 *----------------------------------------------------------------------*/
Mat::FluidPoroMultiPhaseReactions::FluidPoroMultiPhaseReactions(
    Mat::PAR::FluidPoroMultiPhaseReactions* params)
    : FluidPoroMultiPhase(params), paramsreac_(params)
{
  // setup of material map
  if (paramsreac_->local_)
  {
    setup_mat_map();
  }
}

/*----------------------------------------------------------------------*
 | setup of material map                                     vuong 08/16 |
 *----------------------------------------------------------------------*/
void Mat::FluidPoroMultiPhaseReactions::setup_mat_map()
{
  // We just have to add the reaction materials, since the rest is already done in
  // Mat::MatList::setup_mat_map() called from the MatList constructor

  // here's the recursive creation of materials
  std::vector<int>::const_iterator m;
  for (m = paramsreac_->ReacIds()->begin(); m != paramsreac_->ReacIds()->end(); ++m)
  {
    const int reacid = *m;
    Teuchos::RCP<Core::Mat::Material> mat = Mat::Factory(reacid);
    if (mat == Teuchos::null) FOUR_C_THROW("Failed to allocate this material");
    material_map_write()->insert(std::pair<int, Teuchos::RCP<Core::Mat::Material>>(reacid, mat));
  }
  return;
}

/*----------------------------------------------------------------------*
 | reset everything                                          vuong 08/16 |
 *----------------------------------------------------------------------*/
void Mat::FluidPoroMultiPhaseReactions::clear()
{
  paramsreac_ = nullptr;
  return;
}

/*----------------------------------------------------------------------*
 | Unpack data from a char vector into this class            vuong 08/16 |
 *----------------------------------------------------------------------*/
void Mat::FluidPoroMultiPhaseReactions::Pack(Core::Communication::PackBuffer& data) const
{
  Core::Communication::PackBuffer::SizeMarker sm(data);

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  add_to_pack(data, type);

  // matid
  int matid = -1;
  if (paramsreac_ != nullptr) matid = paramsreac_->Id();  // in case we are in post-process mode

  add_to_pack(data, matid);

  // Pack base class material
  Mat::FluidPoroMultiPhase::Pack(data);
}

/*----------------------------------------------------------------------*
 | Unpack data from a char vector into this class            vuong 08/16 |
 *----------------------------------------------------------------------*/
void Mat::FluidPoroMultiPhaseReactions::Unpack(const std::vector<char>& data)
{
  // make sure we have a pristine material
  clear();

  std::vector<char>::size_type position = 0;

  Core::Communication::ExtractAndAssertId(position, data, UniqueParObjectId());

  // matid and recover paramsreac_
  int matid(-1);
  extract_from_pack(position, data, matid);
  paramsreac_ = nullptr;
  if (Global::Problem::Instance()->Materials() != Teuchos::null)
    if (Global::Problem::Instance()->Materials()->Num() != 0)
    {
      const int probinst = Global::Problem::Instance()->Materials()->GetReadFromProblem();
      Core::Mat::PAR::Parameter* mat =
          Global::Problem::Instance(probinst)->Materials()->ParameterById(matid);
      if (mat->Type() == MaterialType())
      {
        // Note: We need to do a dynamic_cast here since Chemotaxis, Reaction, and Chemo-reaction
        // are in a diamond inheritance structure
        paramsreac_ = dynamic_cast<Mat::PAR::FluidPoroMultiPhaseReactions*>(mat);
      }
      else
        FOUR_C_THROW("Type of parameter material %d does not fit to calling type %d", mat->Type(),
            MaterialType());
    }

  // extract base class material
  std::vector<char> basedata(0);
  Mat::FluidPoroMultiPhase::extract_from_pack(position, data, basedata);
  Mat::FluidPoroMultiPhase::Unpack(basedata);

  // in the postprocessing mode, we do not unpack everything we have packed
  // -> position check cannot be done in this case
  if (position != data.size())
    FOUR_C_THROW("Mismatch in size of data %d <-> %d", data.size(), position);
}

/*----------------------------------------------------------------------*
 | reaction ID by Index                                      vuong 08/16 |
 *----------------------------------------------------------------------*/
int Mat::FluidPoroMultiPhaseReactions::ReacID(const unsigned index) const
{
  if ((int)index < paramsreac_->numreac_)
    return paramsreac_->reacids_.at(index);
  else
  {
    FOUR_C_THROW("Index too large");
    return -1;
  }
}

FOUR_C_NAMESPACE_CLOSE
