/*!----------------------------------------------------------------------
\file scatra_reaction_mat.cpp

 \brief

This file contains the base material for reactive scalars.

<pre>
Maintainer: Moritz Thon
            thon@mhpc.mw.tum.de
            http://www.lnm.mw.tum.de
            089-289-10364
</pre>
*----------------------------------------------------------------------*/


#include <vector>
#include "scatra_reaction_mat.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_mat/matpar_bundle.H"
#include "../drt_comm/comm_utils.H"


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::PAR::ScatraReactionMat::ScatraReactionMat(
  Teuchos::RCP<MAT::PAR::Material> matdata
  )
: Parameter(matdata),
  numscal_(matdata->GetInt("NUMSCAL")),
  stoich_(matdata->Get<std::vector<int> >("STOICH")),
  reaccoeff_(matdata->GetDouble("REACCOEFF")),
  coupling_(SetCouplingType(matdata)),
  couprole_(matdata->Get<std::vector<double> >("ROLE")),
  reacstart_((matdata->GetDouble("REACSTART")))
{
  //Some checks for more safety
  if (coupling_ == reac_coup_none)
    dserror("The coupling '%s' is not a valid reaction coupling. Valid couplings are 'simple_multiplicative', 'constand'and michaelis_menten.",(matdata->Get<std::string >("COUPLING"))->c_str() );

  if (numscal_ != (int)stoich_->size())
    dserror("number of scalars %d does not fit to size of the STOICH vector %d", numscal_, stoich_->size());

  if (numscal_ != (int)couprole_->size())
    dserror("number of scalars %d does not fit to size of the ROLE vector %d", numscal_, stoich_->size());

  switch (coupling_)
  {
    case MAT::PAR::reac_coup_simple_multiplicative: //reaction of type A*B*C:
    {
      bool allpositiv = true;
      for (int ii=0; ii < numscal_; ii++)
      {
        if (stoich_->at(ii)<0)
          allpositiv = false;
      }
      if (allpositiv)
        dserror("In the case of simple_multiplicative there must be at least one negative entry in each STOICH list!");
      break;
    }

    case MAT::PAR::reac_coup_power_multiplicative: //reaction of type A*B*C:
    {
      bool allpositiv = true;
      bool rolezero = false;
      for (int ii=0; ii < numscal_; ii++)
      {
        if (stoich_->at(ii)<0)
          allpositiv = false;
        if (stoich_->at(ii)!=0 and couprole_->at(ii) == 0)
          rolezero = true;
      }
      if (allpositiv)
        dserror("In the case of reac_coup_potential_multiplicative there must be at least one negative entry in each STOICH list!");
      if (rolezero)
        dserror("There is one reacting scalar with a zero exponent STOICH list. This does not make sense!");
      break;
    }

    case MAT::PAR::reac_coup_constant: //constant source term:
    {
      bool issomepositiv = false;
      for (int ii=0; ii < numscal_; ii++)
        {
          if (stoich_->at(ii)<0)
            dserror("reac_coup_constant must only contain positive entries in the STOICH list");
          if (stoich_->at(ii)>0)
            issomepositiv=true;
        }
      if (not issomepositiv)
        dserror("reac_coup_constant must contain at least one positive entry in the STOICH list");
      break;
    }

    case MAT::PAR::reac_coup_michaelis_menten: //reaction of type A*B/(B+4)
    {
      bool stoichallzero = true;
      bool roleallzero = true;
      for (int ii=0; ii < numscal_; ii++)
        {
          if (stoich_->at(ii) != 0)
            stoichallzero = false;
          if (couprole_->at(ii) != 0)
            roleallzero = false;
        }
      if (stoichallzero or roleallzero)
        dserror("reac_coup_michaelis_menten must contain at least one non-zero entry in the STOICH and ROLE list");
      break;
    }

    case MAT::PAR::reac_coup_none:
      dserror("reac_coup_none is not a valid coupling");
      break;

    default:
      dserror("The couplingtype %i is not a valid coupling type.", coupling_);
      break;
  }

  return;
}


Teuchos::RCP<MAT::Material> MAT::PAR::ScatraReactionMat::CreateMaterial()
{
  return Teuchos::rcp(new MAT::ScatraReactionMat(this));
}


MAT::ScatraReactionMatType MAT::ScatraReactionMatType::instance_;


void MAT::PAR::ScatraReactionMat::OptParams(std::map<std::string,int>* pnames)
{
  pnames->insert(std::pair<std::string,int>("REACSTART", reacstart_));
}

MAT::PAR::reaction_coupling MAT::PAR::ScatraReactionMat::SetCouplingType( Teuchos::RCP<MAT::PAR::Material> matdata )
{
  if ( *(matdata->Get<std::string >("COUPLING")) == "simple_multiplicative" )
  {
    return reac_coup_simple_multiplicative;
  }
  else if ( *(matdata->Get<std::string >("COUPLING")) == "power_multiplicative")
  {
    return reac_coup_power_multiplicative;
  }
  else if ( *(matdata->Get<std::string >("COUPLING")) == "constant")
  {
    return reac_coup_constant;
  }
  else if ( *(matdata->Get<std::string >("COUPLING")) == "michaelis_menten")
  {
    return reac_coup_michaelis_menten;
  }
  else if ( *(matdata->Get<std::string >("COUPLING")) == "no_coupling")
  {
    return reac_coup_none;
  }
  else
  {
    return reac_coup_none;
  }
}

DRT::ParObject* MAT::ScatraReactionMatType::Create( const std::vector<char> & data )
{
  MAT::ScatraReactionMat* scatra_reaction_mat = new MAT::ScatraReactionMat();
  scatra_reaction_mat->Unpack(data);
  return scatra_reaction_mat;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::ScatraReactionMat::ScatraReactionMat()
  : params_(NULL)
{
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::ScatraReactionMat::ScatraReactionMat(MAT::PAR::ScatraReactionMat* params)
  : params_(params)
{
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ScatraReactionMat::Pack(DRT::PackBuffer& data) const
{
  DRT::PackBuffer::SizeMarker sm( data );
  sm.Insert();

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data,type);

  // matid
  int matid = -1;
  if (params_ != NULL) matid = params_->Id();  // in case we are in post-process mode
  AddtoPack(data,matid);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ScatraReactionMat::Unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;
  // extract type
  int type = 0;
  ExtractfromPack(position,data,type);
  if (type != UniqueParObjectId()) dserror("wrong instance type data. type = %d, UniqueParObjectId()=%d",type,UniqueParObjectId());
  // matid and recover params_
  int matid;
  ExtractfromPack(position,data,matid);
  params_ = NULL;
  if (DRT::Problem::Instance()->Materials() != Teuchos::null)
    if (DRT::Problem::Instance()->Materials()->Num() != 0)
    {
      const int probinst = DRT::Problem::Instance()->Materials()->GetReadFromProblem();
      MAT::PAR::Parameter* mat = DRT::Problem::Instance(probinst)->Materials()->ParameterById(matid);
      if (mat->Type() == MaterialType())
        params_ = static_cast<MAT::PAR::ScatraReactionMat*>(mat);
      else
        dserror("Type of parameter material %d does not fit to calling type %d", mat->Type(), MaterialType());
    }

  if (position != data.size())
    dserror("Mismatch in size of data %d <-> %d",data.size(),position);
}
