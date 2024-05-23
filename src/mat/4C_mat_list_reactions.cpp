/*----------------------------------------------------------------------*/
/*! \file
\brief This file contains the material for reactive scalars. It derives from MAT_matlist
and adds everything to supervise all the MAT_scatra_raction materials. The reactions
itself are defined inside the MAT_scatra_raction materials. So MAT_matlist_reactions
is just a "control instance".

\level 2

*----------------------------------------------------------------------*/


#include "4C_mat_list_reactions.hpp"

#include "4C_global_data.hpp"
#include "4C_mat_par_bundle.hpp"
#include "4C_mat_scatra_reaction.hpp"

#include <vector>

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*
 | standard constructor                                     thon 11/14 |
 *----------------------------------------------------------------------*/
MAT::PAR::MatListReactions::MatListReactions(Teuchos::RCP<CORE::MAT::PAR::Material> matdata)
    : MatList(matdata),
      numreac_((matdata->Get<int>("NUMREAC"))),
      reacids_((matdata->Get<std::vector<int>>("REACIDS")))
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
      Teuchos::RCP<CORE::MAT::Material> mat = MAT::Factory(reacid);
      MaterialMapWrite()->insert(std::pair<int, Teuchos::RCP<CORE::MAT::Material>>(reacid, mat));
    }
  }
}

Teuchos::RCP<CORE::MAT::Material> MAT::PAR::MatListReactions::CreateMaterial()
{
  return Teuchos::rcp(new MAT::MatListReactions(this));
}


MAT::MatListReactionsType MAT::MatListReactionsType::instance_;


CORE::COMM::ParObject* MAT::MatListReactionsType::Create(const std::vector<char>& data)
{
  MAT::MatListReactions* MatListReactions = new MAT::MatListReactions();
  MatListReactions->Unpack(data);
  return MatListReactions;
}


/*----------------------------------------------------------------------*
 | construct empty material object                           thon 11/14 |
 *----------------------------------------------------------------------*/
MAT::MatListReactions::MatListReactions() : MatList(), paramsreac_(nullptr) {}

/*----------------------------------------------------------------------*
 | construct the material object given material parameter    thon 11/14 |
 *----------------------------------------------------------------------*/
MAT::MatListReactions::MatListReactions(MAT::PAR::MatListReactions* params)
    : MatList(params), paramsreac_(params)
{
  // setup of material map
  if (paramsreac_->local_)
  {
    SetupMatMap();
  }
}

/*----------------------------------------------------------------------*
 | setup of material map                                     vuong 08/16 |
 *----------------------------------------------------------------------*/
void MAT::MatListReactions::Initialize()
{
  if (paramsreac_ != nullptr)
  {
    std::vector<int>::const_iterator m;
    for (m = paramsreac_->ReacIds()->begin(); m != paramsreac_->ReacIds()->end(); ++m)
    {
      const int reacid = *m;
      Teuchos::RCP<CORE::MAT::Material> mat = MaterialById(reacid);
      if (mat == Teuchos::null) FOUR_C_THROW("Failed to allocate this material");
      Teuchos::RCP<MAT::ScatraReactionMat> reacmat =
          Teuchos::rcp_dynamic_cast<MAT::ScatraReactionMat>(mat, true);
      reacmat->Initialize();
    }
  }
  return;
}

/*----------------------------------------------------------------------*
 | setup of material map                                     thon 11/14 |
 *----------------------------------------------------------------------*/
void MAT::MatListReactions::SetupMatMap()
{
  // We just have to add the reaction materials, since the rest is already done in
  // MAT::MatList::SetupMatMap() called from the MatList constructor

  // here's the recursive creation of materials
  std::vector<int>::const_iterator m;
  for (m = paramsreac_->ReacIds()->begin(); m != paramsreac_->ReacIds()->end(); ++m)
  {
    const int reacid = *m;
    Teuchos::RCP<CORE::MAT::Material> mat = MAT::Factory(reacid);
    if (mat == Teuchos::null) FOUR_C_THROW("Failed to allocate this material");
    MaterialMapWrite()->insert(std::pair<int, Teuchos::RCP<CORE::MAT::Material>>(reacid, mat));
  }
  return;
}

/*----------------------------------------------------------------------*
 | reset everything                                          thon 11/14 |
 *----------------------------------------------------------------------*/
void MAT::MatListReactions::Clear()
{
  paramsreac_ = nullptr;
  return;
}

/*----------------------------------------------------------------------*
 | Unpack data from a char vector into this class            thon 11/14 |
 *----------------------------------------------------------------------*/
void MAT::MatListReactions::Pack(CORE::COMM::PackBuffer& data) const
{
  CORE::COMM::PackBuffer::SizeMarker sm(data);
  sm.Insert();

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data, type);

  // matid
  int matid = -1;
  if (paramsreac_ != nullptr) matid = paramsreac_->Id();  // in case we are in post-process mode

  AddtoPack(data, matid);

  // Pack base class material
  MAT::MatList::Pack(data);

  if (paramsreac_ != nullptr)
  {
    if (paramsreac_->local_)
    {
      std::vector<int>::const_iterator m;
      for (m = paramsreac_->ReacIds()->begin(); m != paramsreac_->ReacIds()->end(); m++)
      {
        (MaterialMapRead()->find(*m))->second->Pack(data);
      }
    }
  }
}

/*----------------------------------------------------------------------*
 | Unpack data from a char vector into this class            thon 11/14 |
 *----------------------------------------------------------------------*/
void MAT::MatListReactions::Unpack(const std::vector<char>& data)
{
  // make sure we have a pristine material
  Clear();

  std::vector<char>::size_type position = 0;

  CORE::COMM::ExtractAndAssertId(position, data, UniqueParObjectId());

  // matid and recover paramsreac_
  int matid(-1);
  ExtractfromPack(position, data, matid);
  paramsreac_ = nullptr;
  if (GLOBAL::Problem::Instance()->Materials() != Teuchos::null)
    if (GLOBAL::Problem::Instance()->Materials()->Num() != 0)
    {
      const int probinst = GLOBAL::Problem::Instance()->Materials()->GetReadFromProblem();
      CORE::MAT::PAR::Parameter* mat =
          GLOBAL::Problem::Instance(probinst)->Materials()->ParameterById(matid);
      if (mat->Type() == MaterialType())
      {
        // Note: We need to do a dynamic_cast here since Chemotaxis, Reaction, and Chemo-reaction
        // are in a diamond inheritance structure
        paramsreac_ = dynamic_cast<MAT::PAR::MatListReactions*>(mat);
      }
      else
        FOUR_C_THROW("Type of parameter material %d does not fit to calling type %d", mat->Type(),
            MaterialType());
    }

  // extract base class material
  std::vector<char> basedata(0);
  MAT::MatList::ExtractfromPack(position, data, basedata);
  MAT::MatList::Unpack(basedata);

  if (paramsreac_ != nullptr)  // paramsreac_ are not accessible in postprocessing mode
  {
    // make sure the referenced materials in material list have quick access parameters
    std::vector<int>::const_iterator m;
    for (m = paramsreac_->ReacIds()->begin(); m != paramsreac_->ReacIds()->end(); m++)
    {
      const int actmatid = *m;
      Teuchos::RCP<CORE::MAT::Material> mat = MAT::Factory(actmatid);
      if (mat == Teuchos::null) FOUR_C_THROW("Failed to allocate this material");
      MaterialMapWrite()->insert(std::pair<int, Teuchos::RCP<CORE::MAT::Material>>(actmatid, mat));
    }

    if (paramsreac_->local_)
    {
      // loop map of associated local materials
      for (m = paramsreac_->ReacIds()->begin(); m != paramsreac_->ReacIds()->end(); m++)
      {
        std::vector<char> pbtest;
        ExtractfromPack(position, data, pbtest);
        (MaterialMapWrite()->find(*m))->second->Unpack(pbtest);
      }
    }
    // in the postprocessing mode, we do not unpack everything we have packed
    // -> position check cannot be done in this case
    if (position != data.size())
      FOUR_C_THROW("Mismatch in size of data %d <-> %d", data.size(), position);
  }
}

/*----------------------------------------------------------------------*
 | reaction ID by Index                                      thon 11/14 |
 *----------------------------------------------------------------------*/
int MAT::MatListReactions::ReacID(const unsigned index) const
{
  if ((int)index < paramsreac_->numreac_)
    return paramsreac_->reacids_.at(index);
  else
  {
    FOUR_C_THROW("Index too large");
    return -1;
  }
}

/*----------------------------------------------------------------------*
 | calculate advanced reaction terms                         thon 08/16 |
 *----------------------------------------------------------------------*/
double MAT::MatListReactions::calc_rea_body_force_term(
    const int k, const std::vector<double>& phinp, const double* gpcoord, const double scale) const
{
  // set time and space coordinates
  std::vector<std::pair<std::string, double>> constants;
  constants.push_back(std::pair<std::string, double>("t", 0.0));
  constants.push_back(std::pair<std::string, double>("x", gpcoord[0]));
  constants.push_back(std::pair<std::string, double>("y", gpcoord[1]));
  constants.push_back(std::pair<std::string, double>("z", gpcoord[2]));

  double bodyforcetermK = 0.0;

  for (int condnum = 0; condnum < NumReac(); condnum++)
  {
    const int reacid = ReacID(condnum);
    const Teuchos::RCP<const MAT::ScatraReactionMat> reacmat =
        Teuchos::rcp_static_cast<const MAT::ScatraReactionMat>(MaterialById(reacid));

    bodyforcetermK += reacmat->calc_rea_body_force_term(k, phinp, constants, scale);
  }

  return bodyforcetermK;
}

/*----------------------------------------------------------------------*
 | calculate advanced reaction term derivatives              thon 08/16 |
 *----------------------------------------------------------------------*/
void MAT::MatListReactions::calc_rea_body_force_deriv_matrix(const int k,
    std::vector<double>& derivs, const std::vector<double>& phinp, const double* gpcoord,
    const double scale) const
{
  // set time and space coordinates
  std::vector<std::pair<std::string, double>> constants;
  constants.push_back(std::pair<std::string, double>("t", 0.0));
  constants.push_back(std::pair<std::string, double>("x", gpcoord[0]));
  constants.push_back(std::pair<std::string, double>("y", gpcoord[1]));
  constants.push_back(std::pair<std::string, double>("z", gpcoord[2]));

  for (int condnum = 0; condnum < NumReac(); condnum++)
  {
    const int reacid = ReacID(condnum);
    const Teuchos::RCP<const MAT::ScatraReactionMat> reacmat =
        Teuchos::rcp_static_cast<const MAT::ScatraReactionMat>(MaterialById(reacid));

    reacmat->calc_rea_body_force_deriv_matrix(k, derivs, phinp, constants, scale);
  }
  // gpcoord_
  return;
}

/*----------------------------------------------------------------------*
 | calculate advanced reaction terms                         thon 08/16 |
 *----------------------------------------------------------------------*/
double MAT::MatListReactions::calc_rea_body_force_term(const int k,
    const std::vector<double>& phinp, const std::vector<std::pair<std::string, double>>& constants,
    const double* gpcoord, const double scale) const
{
  // add time and space coordinates
  std::vector<std::pair<std::string, double>> constants_mod(constants);
  constants_mod.push_back(std::pair<std::string, double>("t", 0.0));
  constants_mod.push_back(std::pair<std::string, double>("x", gpcoord[0]));
  constants_mod.push_back(std::pair<std::string, double>("y", gpcoord[1]));
  constants_mod.push_back(std::pair<std::string, double>("z", gpcoord[2]));

  double bodyforcetermK = 0.0;

  for (int condnum = 0; condnum < NumReac(); condnum++)
  {
    const int reacid = ReacID(condnum);
    const Teuchos::RCP<const MAT::ScatraReactionMat> reacmat =
        Teuchos::rcp_static_cast<const MAT::ScatraReactionMat>(MaterialById(reacid));

    bodyforcetermK += reacmat->calc_rea_body_force_term(k, phinp, constants_mod, scale);
  }

  return bodyforcetermK;
}

/*----------------------------------------------------------------------*
 | calculate advanced reaction term derivatives              thon 08/16 |
 *----------------------------------------------------------------------*/
void MAT::MatListReactions::calc_rea_body_force_deriv_matrix(const int k,
    std::vector<double>& derivs, const std::vector<double>& phinp,
    const std::vector<std::pair<std::string, double>>& constants, const double* gpcoord,
    const double scale) const
{
  // add time and space coordinates
  std::vector<std::pair<std::string, double>> constants_mod(constants);
  constants_mod.push_back(std::pair<std::string, double>("t", 0.0));
  constants_mod.push_back(std::pair<std::string, double>("x", gpcoord[0]));
  constants_mod.push_back(std::pair<std::string, double>("y", gpcoord[1]));
  constants_mod.push_back(std::pair<std::string, double>("z", gpcoord[2]));

  for (int condnum = 0; condnum < NumReac(); condnum++)
  {
    const int reacid = ReacID(condnum);
    const Teuchos::RCP<const MAT::ScatraReactionMat> reacmat =
        Teuchos::rcp_static_cast<const MAT::ScatraReactionMat>(MaterialById(reacid));

    reacmat->calc_rea_body_force_deriv_matrix(k, derivs, phinp, constants_mod, scale);
  }

  return;
}

/*----------------------------------------------------------------------*
 | add additional variables to reaction                kremheller 07/17 |
 *----------------------------------------------------------------------*/
void MAT::MatListReactions::add_additional_variables(
    const int k, const std::vector<std::pair<std::string, double>>& variables) const
{
  for (int condnum = 0; condnum < NumReac(); condnum++)
  {
    const int reacid = ReacID(condnum);
    const Teuchos::RCP<const MAT::ScatraReactionMat> reacmat =
        Teuchos::rcp_static_cast<const MAT::ScatraReactionMat>(MaterialById(reacid));

    reacmat->add_additional_variables(k, variables);
  }
}

/*--------------------------------------------------------------------------------*
 |  calculating advanced reaction term derivatives after additional variables     |
 |  (e.g. for monolithic coupling)                               kremheller 07/17 |
 *--------------------------------------------------------------------------------*/
void MAT::MatListReactions::calc_rea_body_force_deriv_matrix_add_variables(const int k,
    std::vector<double>& derivs, const std::vector<double>& phinp,
    const std::vector<std::pair<std::string, double>>& variables,
    const std::vector<std::pair<std::string, double>>& constants, const double* gpcoord,
    const double scale) const
{
  // add time and space coordinates
  std::vector<std::pair<std::string, double>> constants_mod(constants);

  // add scalar values as constants
  for (unsigned iscal = 0; iscal < phinp.size(); iscal++)
  {
    std::ostringstream temp;
    temp << iscal + 1;
    constants_mod.push_back(std::pair<std::string, double>("phi" + temp.str(), phinp[iscal]));
  }
  constants_mod.push_back(std::pair<std::string, double>("t", 0.0));
  constants_mod.push_back(std::pair<std::string, double>("x", gpcoord[0]));
  constants_mod.push_back(std::pair<std::string, double>("y", gpcoord[1]));
  constants_mod.push_back(std::pair<std::string, double>("z", gpcoord[2]));

  for (int condnum = 0; condnum < NumReac(); condnum++)
  {
    const int reacid = ReacID(condnum);
    const Teuchos::RCP<const MAT::ScatraReactionMat> reacmat =
        Teuchos::rcp_static_cast<const MAT::ScatraReactionMat>(MaterialById(reacid));

    reacmat->calc_rea_body_force_deriv_matrix_add_variables(
        k, derivs, variables, constants_mod, scale);
  }

  return;
}

FOUR_C_NAMESPACE_CLOSE
