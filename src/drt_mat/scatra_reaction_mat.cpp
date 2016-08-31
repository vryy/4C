/*!----------------------------------------------------------------------
\file scatra_reaction_mat.cpp

 \brief This file contains the base material for reactive scalars.

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
  isreacstart_(false)
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

  switch (coupling_)
  {
    case MAT::PAR::reac_coup_simple_multiplicative: //reaction of type A*B*C:
    {
      bool stoichallzero = true;
      bool isonerolepositive = false;
      for (int ii=0; ii < numscal_; ii++)
        {
          if (stoich_->at(ii) != 0)
            stoichallzero = false;
          if (couprole_->at(ii) > 0.0)
            isonerolepositive = true;
          if (couprole_->at(ii) < 0.0)
            dserror("reac_coup_simple_multiplicative must only contain positive entries in the ROLE list");
        }
      if (not isonerolepositive)
        dserror("reac_coup_simple_multiplicative must contain at least one positive entry in the ROLE list");
      if (stoichallzero and fabs(reaccoeff_)>1.0e-14)
        dserror("reac_coup_michaelis_menten must contain at least one non-zero entry in the STOICH list");

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
      if (stoichallzero and fabs(reaccoeff_)>1.0e-14)
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
                      dserror("reac_coup_constant must only contain positive entries in the STOICH list");
          if (stoich_->at(ii)>0)
            issomepositiv=true;
        }
      if ( (not issomepositiv) and fabs(reaccoeff_)>1.0e-14)
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
      if (stoichallzero and fabs(reaccoeff_)>1.0e-14 )
        dserror("reac_coup_michaelis_menten must contain at least one non-zero entry in the STOICH list");
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

/*----------------------------------------------------------------------/
 | Calculate K(c)                                           Thon 08/16 |
/----------------------------------------------------------------------*/
double MAT::ScatraReactionMat::CalcReaCoeff(
    const int k,                         //!< current scalar id
    const std::vector<double>& phinp,    //!< scalar values at t_(n+1)
    const double scale                   //!< scaling factor for reference concentrations
    ) const
{
  double reactermK=0.0;

  if ( Stoich()->at(k)!=0 and Couprole()->at(k)!=0.0 and fabs(ReacCoeff())>1.0e-14)
  {
    const double rcfac= CalcReaCoeffFac(k,phinp,scale);

    reactermK -= ReacCoeff()*Stoich()->at(k)*rcfac; // scalar at integration point np
  }

  return reactermK;
}

/*----------------------------------------------------------------------/
 | calculate \frac{partial}{\partial c} K(c)                 Thon 08/16 |
/----------------------------------------------------------------------*/
double MAT::ScatraReactionMat::CalcReaCoeffDerivMatrix(
    const int k,                         //!< current scalar id
    const int toderive,                  //!< current id to derive to
    const std::vector<double>& phinp,    //!< scalar values at t_(n+1)
    const double scale                   //!< scaling factor for reference concentrations
    ) const
{
  double reacoeffderivmatrixKJ=0.0;

  if ( Stoich()->at(k)!=0 and Couprole()->at(k)!=0.0 and fabs(ReacCoeff())>1.0e-14)
  {
    const double rcdmfac = CalcReaCoeffDerivFac(k,toderive,phinp,scale);

    reacoeffderivmatrixKJ -= ReacCoeff()*Stoich()->at(k)*rcdmfac;
  }

  return reacoeffderivmatrixKJ;
}


/*----------------------------------------------------------------------/
 | calculate f(c)                                           Thon 08/16 |
/----------------------------------------------------------------------*/
double MAT::ScatraReactionMat::CalcReaBodyForceTerm(
    const int k,                         //!< current scalar id
    const std::vector<double>& phinp,    //!< scalar values at t_(n+1)
    const double scale                   //!< scaling factor for reference concentrations
    ) const
{
  double bodyforcetermK = 0.0;

  if ( Stoich()->at(k)!=0 and Couprole()->at(k)==0.0 and fabs(ReacCoeff())>1.0e-14)
  {
    const double bftfac = CalcReaBodyForceTermFac(k,phinp,scale);// scalar at integration point np

    bodyforcetermK += ReacCoeff()*Stoich()->at(k)*bftfac;
  }

  return bodyforcetermK;
}

/*----------------------------------------------------------------------/
 | calculate \frac{partial}{\partial c} f(c)                 Thon 08/16 |
/----------------------------------------------------------------------*/
double MAT::ScatraReactionMat::CalcReaBodyForceDerivMatrix(
    const int k,                         //!< current scalar id
    const int toderive,                  //!< current id to derive to
    const std::vector<double>& phinp,    //!< scalar values at t_(n+1)
    const double scale                   //!< scaling factor for reference concentrations
    ) const
{
  double reabodyforcederivmatrixKJ=0.0;

  if ( Stoich()->at(k)!=0 and Couprole()->at(k)==0.0 and fabs(ReacCoeff())>1.0e-14)
  {
    const double bfdmfac = CalcReaBodyForceDerivFac(k,toderive,phinp,scale);

    reabodyforcederivmatrixKJ += ReacCoeff()*Stoich()->at(k)*bfdmfac;
  }

  return reabodyforcederivmatrixKJ;
}

/*----------------------------------------------------------------------*
 |  Modify concentrations according to reacstart vector      thon 08/16 |
 *----------------------------------------------------------------------*/
void MAT::ScatraReactionMat::ApplyReacStartAndScaling(
        std::vector<double>& phinp,
        const std::vector<double>* reacstart,
        const double scale
        ) const
{

  for (unsigned int ii=0; ii < phinp.size(); ii++)
  {
    phinp.at(ii) -= reacstart->at(ii);
    if (phinp.at(ii) < 0.0)
      phinp.at(ii)=0.0;

    phinp.at(ii)*=scale;
  }

  return;
}

/*----------------------------------------------------------------------*
 |  helper for calculating K(c)                              thon 08/16 |
 *----------------------------------------------------------------------*/
double MAT::ScatraReactionMat::CalcReaCoeffFac(
    const int k,                         //!< current scalar id
    const std::vector<double>& phinp_org,//!< scalar values at t_(n+1)
    const double scale,                  //!< scaling factor for reference concentrations
    const bool skipreacstart
    ) const
{
  const std::vector<double> couprole = *Couprole();

  std::vector<double> phinp(phinp_org);
  if ( (IsReacStart() or scale != 1.0) and not skipreacstart )
    ApplyReacStartAndScaling( phinp,ReacStart(),scale );

  double rcfac=1.0;

  switch ( Coupling() )
  {
    case MAT::PAR::reac_coup_simple_multiplicative: //reaction of type A*B*C:
    {
      for (int ii=0; ii < NumScal(); ii++)
      {
        if (couprole[ii]!=0.0)
        {
          if (ii!=k)
            rcfac *=phinp.at(ii);
        }
      }
      break;
    }

    case MAT::PAR::reac_coup_power_multiplicative: //reaction of type A^2*B^-1.5*C:
    {
      for (int ii=0; ii < NumScal(); ii++)
      {
        if (couprole[ii]!=0.0)
        {
          if (ii!=k)
            rcfac *= std::pow(phinp.at(ii),couprole[ii]);
          else
            rcfac *= std::pow(phinp.at(ii),couprole[ii]-1.0);
        }
      }

      if ( IsReacStart() )
        dserror("The reacstart feature is only tested for reactions of type simple_multiplicative. It should work, but be careful!");
      break;
    }

    case MAT::PAR::reac_coup_constant: //constant source term:
    {
      rcfac = 0.0;

      if ( IsReacStart() )
        dserror("The reacstart feature is only tested for reactions of type simple_multiplicative. It should work, but be careful!");
      break;
    }

     case MAT::PAR::reac_coup_michaelis_menten: //reaction of type A*B/(B+4)
    {
      for (int ii=0; ii < NumScal(); ii++)
      {
        if (couprole[k]<0.0)
        {
          if ((couprole[ii] < 0.0) and (ii != k))
            rcfac *= phinp.at(ii);
          else if ((couprole[ii] > 0.0) and (ii!=k))
            rcfac *= (phinp.at(ii)/(phinp.at(ii) + couprole[ii]));
        }
        else if (couprole[k]>0.0)
        {
          if (ii == k)
            rcfac *= (1/(phinp.at(ii)+couprole[ii]));
          else if (couprole[ii] < 0.0)
            rcfac *= phinp.at(ii);
          else if (couprole[ii] > 0.0)
            rcfac *= (phinp.at(ii)/(phinp.at(ii) + couprole[ii]));
        }
      }
      break;
    }

    case MAT::PAR::reac_coup_none:
      dserror("reac_coup_none is not a valid coupling");
      break;

    default:
      dserror("The couplingtype %i is not a valid coupling type.", Coupling());
      break;
  }

  //reaction start feature
  if ( ReacStart()->at(k) > 0.0 and not skipreacstart )
  {
    if ( phinp.at(k) > 0.0 ) //! Calculate (K(c)-reacstart(c)./c)_{+}
      rcfac *= (1-ReacStart()->at(k)/(phinp_org.at(k)*scale) );
    else
      rcfac = 0.0;
  }

    return rcfac;
}

/*----------------------------------------------------------------------*
 |  helper for calculating \frac{partial}{\partial c} K(c)   thon 08/16 |
 *----------------------------------------------------------------------*/
double MAT::ScatraReactionMat::CalcReaCoeffDerivFac(
    const int k,                         //!< current scalar id
    const int toderive,                  //!< current id to derive to
    const std::vector<double>& phinp_org,//!< scalar values at t_(n+1)
    const double scale                   //!< scaling factor for reference concentrations
    ) const
{
  const std::vector<double> couprole = *Couprole();

  std::vector<double> phinp(phinp_org);
  if ( IsReacStart() or scale != 1.0 )
    ApplyReacStartAndScaling( phinp,ReacStart(),scale );

  double rcdmfac=1.0;

  switch (Coupling())
  {
    case MAT::PAR::reac_coup_simple_multiplicative: //reaction of type A*B*C:
    {
      if (couprole[toderive]!=0 and toderive!=k)
      {
        for (int ii=0; ii < NumScal(); ii++)
        {
          if (couprole[ii]!=0 and ii!=k and ii!= toderive)
            rcdmfac *= phinp.at(ii);
        }
      }
      else
        rcdmfac=0.0;
      break;
    }

    case MAT::PAR::reac_coup_power_multiplicative: //reaction of type A^2*B^-1.5*C:
    {
      if (couprole[toderive]!=0 and toderive!=k)
      {
        for (int ii=0; ii < NumScal(); ii++)
        {
          if (couprole[ii]!=0 and ii!=k and ii!= toderive)
            rcdmfac *= std::pow(phinp.at(ii),couprole[ii]);
          else if(couprole[ii]!=0 and ii!=k and ii== toderive)
            rcdmfac *= couprole[ii]*std::pow(phinp.at(ii),couprole[ii]-1.0);
          else if(couprole[ii]!=0 and ii==k and ii== toderive and couprole[ii]!=1.0)
            rcdmfac *= (couprole[ii]-1.0)*std::pow(phinp.at(ii),couprole[ii]-2.0);
        }
      }
      else
        rcdmfac=0.0;
      break;
    }

    case MAT::PAR::reac_coup_constant: //constant source term:
    {
      rcdmfac = 0.0;
      break;
    }

    case MAT::PAR::reac_coup_michaelis_menten: //reaction of type A*B/(B+4)
    {
      for (int ii=0; ii < NumScal(); ii++)
      {
        if (couprole[k] < 0.0)
        {
          if ((couprole[toderive]==0.0) or (k==toderive))
            rcdmfac = 0.0;
          else if ( (k!= toderive) and (couprole[toderive]<0.0) )
          {
            if (ii==k)
              rcdmfac *= 1.0;
            else if (ii!=toderive and couprole[ii]<0.0)
              rcdmfac *= phinp.at(ii);
            else if (ii!=toderive and couprole[ii]>0.0)
              rcdmfac *= phinp.at(ii)/(couprole[ii]+phinp.at(ii));
            else if (ii!=toderive and couprole[ii] == 0.0)
              rcdmfac *= 1.0;
            else if (ii == toderive)
              rcdmfac *= 1.0;
          }
          else if ( (k!= toderive) and (couprole[toderive]>00.) )
          {
            if (ii==k)
              rcdmfac *= 1;
            else if (ii!=toderive and couprole[ii]<0.0)
              rcdmfac *= phinp.at(ii);
            else if (ii!=toderive and couprole[ii]>0.0)
              rcdmfac *= phinp.at(ii)/(couprole[ii]+phinp.at(ii));
            else if (ii!=toderive and couprole[ii]==0.0)
              rcdmfac *= 1.0;
            else if (ii == toderive)
              rcdmfac *= couprole[ii]/std::pow((couprole[ii]+phinp.at(ii)),2);
          }
        }
        else if (couprole[k] > 0.0 )
        {
          if (couprole[toderive]==0)
            rcdmfac = 0.0;
          else if ( (k!=toderive) and couprole[toderive]<0.0)
          {
            if (ii==k)
              rcdmfac *= 1/(couprole[ii]+phinp.at(ii));
            else if ((ii !=toderive) and (couprole[ii]<0.0))
              rcdmfac *= phinp.at(ii);
            else if ((ii!=toderive) and (couprole[ii]>0.0))
              rcdmfac *= phinp.at(ii)/(couprole[ii]+phinp.at(ii));
            else if ((ii!=toderive) and (couprole[ii]==0.0))
              rcdmfac *= 1.0;
            else if (ii==toderive)
              rcdmfac *= 1.0;
          }
          else if ( (k!=toderive) and couprole[toderive]>0.0)
          {
            if (ii==k)
              rcdmfac *= 1/(couprole[ii]+phinp.at(ii));
            else if ((ii!=toderive) and couprole[ii]<0.0)
              rcdmfac *= phinp.at(ii);
            else if ((ii!=toderive) and couprole[ii]>0.0)
              rcdmfac *= phinp.at(ii)/(couprole[ii] + phinp.at(ii));
            else if ((ii!=toderive) and (couprole[ii]==0.0))
              rcdmfac *= 1.0;
            else if (ii==toderive)
              rcdmfac *= couprole[ii]/std::pow((couprole[ii]+phinp.at(ii)),2);
          }
          else if (k==toderive)
          {
            if (ii==k)
              rcdmfac *= -1/std::pow((couprole[ii]+phinp.at(ii)),2);
            else if ((ii!=toderive) and (couprole[ii]<0.0))
              rcdmfac *= phinp.at(ii);
            else if ((ii!=toderive) and (couprole[ii]>0.0))
              rcdmfac *= phinp.at(ii)/(couprole[ii] + phinp.at(ii));
            else if ((ii!=toderive) and (couprole[ii]==0.0))
              rcdmfac *= 1.0;
          }
        }
      }

      break;
    }

    case MAT::PAR::reac_coup_none:
      dserror("reac_coup_none is not a valid coupling");
      break;

    default:
      dserror("The couplingtype %i is not a valid coupling type.", Coupling());
      break;
  }

  //reaction start feature
  if ( ReacStart()->at(k) > 0.0 )
  {
    if ( phinp.at(k) > 0.0 ) //! Calculate (K(c)-reacstart(c)./c)_{+}
    {
      rcdmfac *= (1-ReacStart()->at(k)/(phinp_org.at(k)*scale) );

      if ( k==toderive )
      {
        const double rcfac= CalcReaCoeffFac(k,phinp,1.0,true);

        rcdmfac += rcfac*ReacStart()->at(k)/std::pow(phinp_org.at(k)*scale,2);
      }
    }
    else
      rcdmfac = 0.0;
  }

  return rcdmfac;
}

/*----------------------------------------------------------------------*
 |  helper for calculating f(c)                              thon 08/16 |
 *----------------------------------------------------------------------*/
double MAT::ScatraReactionMat::CalcReaBodyForceTermFac(
    const int k,                         //!< current scalar id
    const std::vector<double>& phinp_org,//!< scalar values at t_(n+1)
    const double scale                   //!< scaling factor for reference concentrations
    ) const
{
  const std::vector<double> couprole = *Couprole();

  std::vector<double> phinp(phinp_org);
  if ( IsReacStart() or scale != 1.0 )
    ApplyReacStartAndScaling( phinp,ReacStart(),scale );

  double bftfac=1.0;

  switch ( Coupling() )
  {
    case MAT::PAR::reac_coup_simple_multiplicative: //reaction of type A*B*C:
    {
      for (int ii=0; ii < NumScal(); ii++)
      {
        if (couprole[ii]!=0)
        {
          bftfac *=phinp[ii];
        }
      }
      break;
    }

    case MAT::PAR::reac_coup_power_multiplicative: //reaction of type A^2*B^-1.5*C:
    {
      for (int ii=0; ii < NumScal(); ii++)
      {
        if (couprole[ii]!=0)
        {
          bftfac *= std::pow(phinp[ii],couprole[ii]);
        }
      }
      break;
    }

    case MAT::PAR::reac_coup_constant: //constant source term:
    {
      bftfac = 1.0;
      break;
    }

    case MAT::PAR::reac_coup_michaelis_menten: //reaction of type A*B/(B+4)
    {
      for (int ii=0; ii < NumScal() ; ii++)
      {
        if (couprole[ii]>0.0) //and (ii!=k))
          bftfac *= phinp[ii]/(couprole[ii]+phinp[ii]);
        else if (couprole[ii]<0.0) //and (ii!=k))
          bftfac *= phinp[ii];
      }
      break;
    }

    case MAT::PAR::reac_coup_none:
      dserror("reac_coup_none is not a valid coupling");
      break;

    default:
      dserror("The couplingtype %i is not a valid coupling type.", Coupling() );
      break;
  }

  return bftfac;
}

/*--------------------------------------------------------------------------------*
 |  helper for calculating calculate \frac{partial}{\partial c} f(c)   thon 08/16 |
 *--------------------------------------------------------------------------------*/
double MAT::ScatraReactionMat::CalcReaBodyForceDerivFac(
    const int k,                         //!< current scalar id
    const int toderive,                  //!< current id to derive to
    const std::vector<double>& phinp_org,//!< scalar values at t_(n+1)
    const double scale                   //!< scaling factor for reference concentrations
    ) const
{
  const std::vector<double> couprole = *Couprole();

  std::vector<double> phinp(phinp_org);
  if ( IsReacStart() or scale != 1.0 )
    ApplyReacStartAndScaling( phinp,ReacStart(),scale );

  double bfdmfac=1.0;

  switch (Coupling())
  {
    case MAT::PAR::reac_coup_simple_multiplicative: //reaction of type A*B*C:
    {
      if (couprole[toderive]!=0)
        {
          for (int ii=0; ii < NumScal(); ii++)
          {
            if (couprole[ii]!=0 and ii!=toderive)
              bfdmfac *= phinp.at(ii);
          }
        }
      else
        bfdmfac=0.0;
      break;
    }

    case MAT::PAR::reac_coup_power_multiplicative: //reaction of type A^2*B^-1.5*C:
    {
      if (couprole[toderive]!=0)
        {
          for (int ii=0; ii < NumScal(); ii++)
          {
            if (couprole[ii]!=0 and ii!=toderive)
              bfdmfac *= std::pow(phinp.at(ii),couprole[ii]);
            else if(couprole[ii]!=0 and ii==toderive)
              bfdmfac *= std::pow(phinp.at(ii),couprole[ii]-1.0);
          }
        }
      else
        bfdmfac=0.0;
      break;
    }

    case MAT::PAR::reac_coup_constant: //constant source term:
    {
      bfdmfac = 0.0;
      break;
    }

    case MAT::PAR::reac_coup_michaelis_menten: //reaction of type A*B/(B+4)
    {
      for (int ii=0; ii < NumScal(); ii++)
      {
        if (ii != toderive)
        {
          if (couprole[ii] > 0.0)
            bfdmfac *= phinp.at(ii)/(couprole[ii] + phinp.at(ii));
          else if (couprole[ii] < 0.0)
            bfdmfac *= phinp.at(ii);
          else
            bfdmfac *= 1;
        }
        else
        {
          if (couprole[ii] > 0.0)
            bfdmfac *= couprole[ii]/(std::pow((phinp.at(ii)+couprole[ii]), 2));
          else if (couprole[ii] < 0.0)
            bfdmfac *= 1;
          else
            bfdmfac = 0;
        }
      }
      break;
    }

    case MAT::PAR::reac_coup_none:
      dserror("reac_coup_none is not a valid coupling");
      break;

    default:
      dserror("The couplingtype %i is not a valid coupling type.", Coupling());
      break;
  }

  if ( ReacStart()->at(toderive) > 0 and phinp.at(toderive)==0.0 )
    bfdmfac=0.0;

  return bfdmfac;
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

  return (Stoich()->at(k)*CalcReaBodyForceTermFac(k,phinp,scale));
}

/*---------------------------------------------------------------------------------/
 | Calculate influence factor for scalar dependent membrane transport   Thon 08/16 |
/--------------------------------------------------------------------------------- */
double MAT::ScatraReactionMat::CalcPermInfluenceDeriv(
    const int k,                         //!< current scalar id
    const int toderive,                  //!< current id to derive to
    const std::vector<double>& phinp,    //!< scalar values at t_(n+1)
    const double scale                   //!< scaling factor for reference concentrations
    ) const
{
  return (Stoich()->at(k)*CalcReaBodyForceDerivFac(k,toderive,phinp,scale));
}
