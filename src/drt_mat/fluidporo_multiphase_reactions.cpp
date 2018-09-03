/*----------------------------------------------------------------------*/
/*!
 \file fluidporo_multiphase_reactions.cpp

 \brief a fluid material for porous multiphase flow with reactions (mass sources and sinks)

   \level 3

   \maintainer  Lena Yoshihara
                yoshihara@lnm.mw.tum.de
                http://www.lnm.mw.tum.de
 *----------------------------------------------------------------------*/


#include <vector>
#include "fluidporo_multiphase_reactions.H"
#include "fluidporo_multiphase_singlereaction.H"
#include "matpar_bundle.H"
#include "../drt_lib/drt_globalproblem.H"

/*----------------------------------------------------------------------*
 | rstandard constructor                                     vuong 08/16 |
 *----------------------------------------------------------------------*/
MAT::PAR::FluidPoroMultiPhaseReactions::FluidPoroMultiPhaseReactions(
    Teuchos::RCP<MAT::PAR::Material> matdata)
    : FluidPoroMultiPhase(matdata),
      numreac_((matdata->GetInt("NUMREAC"))),
      reacids_((matdata->Get<std::vector<int>>("REACIDS")))
{
  // check if sizes fit
  if (numreac_ != (int)reacids_->size())
    dserror("number of materials %d does not fit to size of material vector %d", nummat_,
        reacids_->size());

  if (numreac_ < 1)
    dserror("if you don't have reactions, use MAT_matlist instead of MAT_matlist_reactions!");

  if (not local_)
  {
    // make sure the referenced materials in material list have quick access parameters
    std::vector<int>::const_iterator m;
    for (m = reacids_->begin(); m != reacids_->end(); ++m)
    {
      const int reacid = *m;
      Teuchos::RCP<MAT::Material> mat = MAT::Material::Factory(reacid);

      // safety check and cast
      if (mat->MaterialType() != INPAR::MAT::m_fluidporo_singlereaction)
        dserror("only MAT_FluidPoroSingleReaction material valid");
      MAT::FluidPoroSingleReaction singlereacmat = static_cast<MAT::FluidPoroSingleReaction&>(*mat);
      if (singlereacmat.TotalNumDof() != this->nummat_)
        dserror(
            "TOTALNUMDOF in MAT_FluidPoroSingleReaction does not correspond to NUMMAT in "
            "MAT_FluidPoroMultiPhaseReactions");

      MaterialMapWrite()->insert(std::pair<int, Teuchos::RCP<MAT::Material>>(reacid, mat));
    }
  }
}

Teuchos::RCP<MAT::Material> MAT::PAR::FluidPoroMultiPhaseReactions::CreateMaterial()
{
  return Teuchos::rcp(new MAT::FluidPoroMultiPhaseReactions(this));
}


MAT::FluidPoroMultiPhaseReactionsType MAT::FluidPoroMultiPhaseReactionsType::instance_;


DRT::ParObject* MAT::FluidPoroMultiPhaseReactionsType::Create(const std::vector<char>& data)
{
  MAT::FluidPoroMultiPhaseReactions* FluidPoroMultiPhaseReactions =
      new MAT::FluidPoroMultiPhaseReactions();
  FluidPoroMultiPhaseReactions->Unpack(data);
  return FluidPoroMultiPhaseReactions;
}


/*----------------------------------------------------------------------*
 | construct empty material object                           vuong 08/16 |
 *----------------------------------------------------------------------*/
MAT::FluidPoroMultiPhaseReactions::FluidPoroMultiPhaseReactions()
    : FluidPoroMultiPhase(), paramsreac_(NULL)
{
}

/*----------------------------------------------------------------------*
 | construct the material object given material paramete     vuong 08/16 |
 *----------------------------------------------------------------------*/
MAT::FluidPoroMultiPhaseReactions::FluidPoroMultiPhaseReactions(
    MAT::PAR::FluidPoroMultiPhaseReactions* params)
    : FluidPoroMultiPhase(params), paramsreac_(params)
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
void MAT::FluidPoroMultiPhaseReactions::SetupMatMap()
{
  // We just have to add the reaction materials, since the rest is already done in
  // MAT::MatList::SetupMatMap() called from the MatList constructor

  // here's the recursive creation of materials
  std::vector<int>::const_iterator m;
  for (m = paramsreac_->ReacIds()->begin(); m != paramsreac_->ReacIds()->end(); ++m)
  {
    const int reacid = *m;
    Teuchos::RCP<MAT::Material> mat = MAT::Material::Factory(reacid);
    if (mat == Teuchos::null) dserror("Failed to allocate this material");
    MaterialMapWrite()->insert(std::pair<int, Teuchos::RCP<MAT::Material>>(reacid, mat));
  }
  return;
}

/*----------------------------------------------------------------------*
 | reset everything                                          vuong 08/16 |
 *----------------------------------------------------------------------*/
void MAT::FluidPoroMultiPhaseReactions::Clear()
{
  paramsreac_ = NULL;
  return;
}

/*----------------------------------------------------------------------*
 | Unpack data from a char vector into this class            vuong 08/16 |
 *----------------------------------------------------------------------*/
void MAT::FluidPoroMultiPhaseReactions::Pack(DRT::PackBuffer& data) const
{
  DRT::PackBuffer::SizeMarker sm(data);
  sm.Insert();

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data, type);

  // matid
  int matid = -1;
  if (paramsreac_ != NULL) matid = paramsreac_->Id();  // in case we are in post-process mode

  AddtoPack(data, matid);

  // Pack base class material
  MAT::FluidPoroMultiPhase::Pack(data);
}

/*----------------------------------------------------------------------*
 | Unpack data from a char vector into this class            vuong 08/16 |
 *----------------------------------------------------------------------*/
void MAT::FluidPoroMultiPhaseReactions::Unpack(const std::vector<char>& data)
{
  // make sure we have a pristine material
  Clear();

  std::vector<char>::size_type position = 0;
  // extract type
  int type = 0;
  ExtractfromPack(position, data, type);
  if (type != UniqueParObjectId()) dserror("wrong instance type data");

  // matid and recover paramsreac_
  int matid(-1);
  ExtractfromPack(position, data, matid);
  paramsreac_ = NULL;
  if (DRT::Problem::Instance()->Materials() != Teuchos::null)
    if (DRT::Problem::Instance()->Materials()->Num() != 0)
    {
      const int probinst = DRT::Problem::Instance()->Materials()->GetReadFromProblem();
      MAT::PAR::Parameter* mat =
          DRT::Problem::Instance(probinst)->Materials()->ParameterById(matid);
      if (mat->Type() == MaterialType())
      {
        // Note: We need to do a dynamic_cast here since Chemotaxis, Reaction, and Chemo-reaction
        // are in a diamond inheritance structure
        paramsreac_ = dynamic_cast<MAT::PAR::FluidPoroMultiPhaseReactions*>(mat);
      }
      else
        dserror("Type of parameter material %d does not fit to calling type %d", mat->Type(),
            MaterialType());
    }

  // extract base class material
  std::vector<char> basedata(0);
  MAT::FluidPoroMultiPhase::ExtractfromPack(position, data, basedata);
  MAT::FluidPoroMultiPhase::Unpack(basedata);

  // in the postprocessing mode, we do not unpack everything we have packed
  // -> position check cannot be done in this case
  if (position != data.size()) dserror("Mismatch in size of data %d <-> %d", data.size(), position);
}

/*----------------------------------------------------------------------*
 | reaction ID by Index                                      vuong 08/16 |
 *----------------------------------------------------------------------*/
int MAT::FluidPoroMultiPhaseReactions::ReacID(const unsigned index) const
{
  if ((int)index < paramsreac_->numreac_)
    return paramsreac_->reacids_->at(index);
  else
  {
    dserror("Index too large");
    return -1;
  }
}
