/*!----------------------------------------------------------------------
\file scatra_reaction_mat.cpp

 \brief This file contains the base material for reactive scalars. This includes all
       calculations of the reactions terms and all its derivatives.

\level 2
<pre>
\maintainer Moritz Thon
            thon@mhpc.mw.tum.de
            http://www.lnm.mw.tum.de
            089-289-10364
</pre>
*----------------------------------------------------------------------*/


#include <vector>
#include "scatra_reaction_mat.H"
#include "scatra_reaction_coupling.H"
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
  reacstart_(matdata->Get<std::vector<double> >("REACSTART")),
  isreacstart_(false),
  isinit_(false)
{
  //Some checks for more safety
  if (coupling_ == reac_coup_none)
    dserror("The coupling '%s' is not a valid reaction coupling. Valid couplings are 'simple_multiplicative', 'constand'and michaelis_menten.",(matdata->Get<std::string >("COUPLING"))->c_str() );

  if (numscal_ != (int)stoich_->size())
    dserror("number of scalars %d does not fit to size of the STOICH vector %d", numscal_, stoich_->size());

  if (numscal_ != (int)couprole_->size())
    dserror("number of scalars %d does not fit to size of the ROLE vector %d", numscal_, stoich_->size());

  if (numscal_ != (int)reacstart_->size())
    dserror("number of scalars %d does not fit to size of the REACSTART vector %d", numscal_, reacstart_->size());

  for (int ii=0; ii < numscal_; ii++)
  {
    if (reacstart_->at(ii)<0)
    {
      dserror("In the REACSTART vector only non-negative values are allowed!");
    }
    else if (reacstart_->at(ii)>0)
    {
      isreacstart_=true;
    }
  }

  // do some more input checks depending on coupling type
  {
    switch (coupling_)
    {
      case MAT::PAR::reac_coup_simple_multiplicative: //reaction of type A*B*C:
      {
        bool stoichallzero = true;
        bool roleallzero = true;
        for (int ii=0; ii < numscal_; ii++)
          {
            if (stoich_->at(ii) != 0)
              stoichallzero = false;
            if (couprole_->at(ii) != 0.0)
              roleallzero = false;
          }
        if (roleallzero)
          dserror("reac_coup_simple_multiplicative must contain at least one non-zero entry in the ROLE list");
        if (stoichallzero)
          dserror("reac_coup_simple_multiplicative must contain at least one non-zero entry in the STOICH list");

        break;
      }

      case MAT::PAR::reac_coup_power_multiplicative: //reaction of type A^2*B^-1.5*C:
      {
        bool stoichallzero = true;
        bool roleallzero = true;
        for (int ii=0; ii < numscal_; ii++)
          {
            if (stoich_->at(ii) != 0)
              stoichallzero = false;
            if (couprole_->at(ii) != 0.0)
              roleallzero = false;
          }
        if (roleallzero)
          dserror("reac_coup_power_multiplicative must contain at least one positive entry in the ROLE list");
        if (stoichallzero)
          dserror("reac_coup_michaelis_menten must contain at least one non-zero entry in the STOICH list");

        break;
      }

      case MAT::PAR::reac_coup_constant: //constant source term:
      {
        bool issomepositiv = false;
        for (int ii=0; ii < numscal_; ii++)
          {
            if (stoich_->at(ii)<0)
              dserror("reac_coup_constant must only contain positive entries in the STOICH list");
            if (couprole_->at(ii) != 0.0)
              dserror("reac_coup_constant must only contain zero entries in the ROLE list");
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
        if (roleallzero)
          dserror("reac_coup_michaelis_menten must contain at least one non-zero entry in the ROLE list");
        if (stoichallzero)
          dserror("reac_coup_michaelis_menten must contain at least one non-zero entry in the STOICH list");

        break;
      }

      case MAT::PAR::reac_coup_byfunction: //reaction by function
      {
        int functID = -1;
        for (int ii=0; ii < numscal_; ii++)
          {
            if (stoich_->at(ii) != 0)
            {
              if (round(couprole_->at(ii)) < 1)
                dserror("reac_coup_byfunction: no function defined in the ROLE list for scalar with positive entry in the STOICH list");
              if(functID==-1)
                functID=round(couprole_->at(ii));
              else if(functID!=round(couprole_->at(ii)))
                dserror("The FUNC IDs defined in the ROLE list should all match");
            }
          }
        if(functID==-1)
          dserror("reac_coup_byfunction must contain at least one positive entry in the STOICH list");

        break;
      }

      case MAT::PAR::reac_coup_none:
        dserror("reac_coup_none is not a valid coupling");
        break;

      default:
        dserror("The couplingtype %i is not a valid coupling type.", coupling_);
        break;
    }
  }

  // if all checks are passed, we can build the reaction class
  reaction_ = MAT::PAR::REACTIONCOUPLING::ReactionInterface::CreateReaction(coupling_,isreacstart_,*reacstart_);

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void MAT::PAR::ScatraReactionMat::Initialize()
{
  reaction_->Initialize(numscal_,*couprole_);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<MAT::Material> MAT::PAR::ScatraReactionMat::CreateMaterial()
{
  return Teuchos::rcp(new MAT::ScatraReactionMat(this));
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
MAT::ScatraReactionMatType MAT::ScatraReactionMatType::instance_;

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
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
  else if ( *(matdata->Get<std::string >("COUPLING")) == "by_function")
  {
    return reac_coup_byfunction;
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

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
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

/*----------------------------------------------------------------------/
 | calculate advanced reaction terms                        Thon 08/16 |
/----------------------------------------------------------------------*/
double MAT::ScatraReactionMat::CalcReaBodyForceTerm(
    const int k,                         //!< current scalar id
    const std::vector<double>& phinp,    //!< scalar values at t_(n+1)
    const double scale_phi               //!< scaling factor for scalar values (used for reference concentrations)
    ) const
{
  if ( Stoich()->at(k)!=0 and fabs(ReacCoeff())>1.0e-14)
  {
    return CalcReaBodyForceTerm(k,phinp,ReacCoeff()*Stoich()->at(k),scale_phi);// scalar at integration point np
  }
  else
    return 0.0;
}

/*----------------------------------------------------------------------/
 | calculate advanced reaction term derivatives             Thon 08/16 |
/----------------------------------------------------------------------*/
void MAT::ScatraReactionMat::CalcReaBodyForceDerivMatrix(
    const int k,                         //!< current scalar id
    std::vector<double>& derivs,         //!< vector with derivatives (to be filled)
    const std::vector<double>& phinp,    //!< scalar values at t_(n+1)
    const double scale_phi               //!< scaling factor for scalar values (used for reference concentrations)
    ) const
{
  if ( Stoich()->at(k)!=0 and fabs(ReacCoeff())>1.0e-14)
  {
    CalcReaBodyForceDeriv(k,derivs,phinp,ReacCoeff()*Stoich()->at(k),scale_phi);
  }

  return;
}

/*----------------------------------------------------------------------*
 |  helper for calculating advanced reaction terms           thon 08/16 |
 *----------------------------------------------------------------------*/
double MAT::ScatraReactionMat::CalcReaBodyForceTerm(
    int k,                               //!< current scalar id
    const std::vector<double>& phinp,    //!< scalar values at t_(n+1)
    double scale_reac,                   //!< scaling factor for reaction term (= reaction coefficient * stoichometry)
    double scale_phi                     //!< scaling factor for scalar values (used for reference concentrations)
    ) const
{
  return params_->reaction_->CalcReaBodyForceTerm(k,NumScal(),phinp,*Couprole(),scale_reac,scale_phi);
}

/*--------------------------------------------------------------------------------*
 |  helper for calculating advanced reaction term derivatives          thon 08/16 |
 *--------------------------------------------------------------------------------*/
void MAT::ScatraReactionMat::CalcReaBodyForceDeriv(
    int k,                               //!< current scalar id
    std::vector<double>& derivs,         //!< vector with derivatives (to be filled)
    const std::vector<double>& phinp,    //!< scalar values at t_(n+1)
    double scale_reac,                   //!< scaling factor for reaction term (= reaction coefficient * stoichometry)
    double scale_phi                    //!< scaling factor for scalar values (used for reference concentrations)
    ) const
{
  params_->reaction_->CalcReaBodyForceDeriv(k,NumScal(),derivs,phinp,*Couprole(),scale_reac,scale_phi);

  return;
}

/*---------------------------------------------------------------------------------/
 | Calculate influence factor for scalar dependent membrane transport   Thon 08/16 |
/--------------------------------------------------------------------------------- */
double MAT::ScatraReactionMat::CalcPermInfluence(
    const int k,                         //!< current scalar id
    const std::vector<double>& phinp,    //!< scalar values at t_(n+1)
    const double scale                   //!< scaling factor for reference concentrations
    ) const
{
  if ( not (Stoich()->at(k) > 0) )
    dserror("You need to specify a positive STOICH entry for scalar %i",k);
  if ( fabs(ReacCoeff()) > 1.0e-14 )
    dserror("You need to set REACOEFF to 0.0!");

  return (CalcReaBodyForceTerm(k,phinp,Stoich()->at(k),scale));
}

/*---------------------------------------------------------------------------------/
 | Calculate influence factor for scalar dependent membrane transport   Thon 08/16 |
/--------------------------------------------------------------------------------- */
void MAT::ScatraReactionMat::CalcPermInfluenceDeriv(
    const int k,                         //!< current scalar id
    std::vector<double>& derivs,         //!< vector with derivatives (to be filled)
    const std::vector<double>& phinp,    //!< scalar values at t_(n+1)
    const double scale                   //!< scaling factor for reference concentrations
    ) const
{
  CalcReaBodyForceDeriv(k,derivs,phinp,Stoich()->at(k),scale);
}
