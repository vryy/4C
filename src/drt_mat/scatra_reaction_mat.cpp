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

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void MAT::PAR::ScatraReactionMat::Initialize()
{
  if(not isinit_)
  {
    switch (coupling_)
    {
      case MAT::PAR::reac_coup_byfunction: //reaction by function
      {
        for(int ii=0;ii<numscal_;ii++)
        {
          // we take the value in couprole list as function ID
          const int functID = round(couprole_->at(ii));
          if (functID!=0)
          {
            if(Function(functID-1).NumberComponents()!=1)
              dserror("expected only one component for the reaction evaluation");

            for(int k=0;k<numscal_;k++)
            {
              // construct the strings for scalar
              std::ostringstream temp;
              temp << k+1;
              std::string name = "phi"+temp.str();

              // add the variable name to the parser
              if(not Function(functID-1).IsVariable(0,name))
                Function(functID-1).AddVariable(0,name,0.0);
            }
          }
        }
        break;
      }
      default:
        break;
    }
  }

  isinit_=true;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
inline DRT::UTILS::VariableExprFunction& MAT::PAR::ScatraReactionMat::Function(int functnum) const
{
  // try to cast to variable expression function for phi evaluation
  // try-catch because of cast of references
  try
  {
    DRT::UTILS::VariableExprFunction& funct =
        dynamic_cast<DRT::UTILS::VariableExprFunction&>(DRT::Problem::Instance()->Funct(functnum));

    return funct;
  }
  catch(std::bad_cast & exp)
  {
    dserror("Cast to VarExp Function failed! For reaction law definition only 'VAREXPR' functions are allowed!\n"
        "Check your input file!");
    return dynamic_cast<DRT::UTILS::VariableExprFunction&>(DRT::Problem::Instance()->Funct(functnum));
  }
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
    const double scale                   //!< scaling factor for reference concentrations
    ) const
{
  double bodyforcetermK = 0.0;

  if ( Stoich()->at(k)!=0 and fabs(ReacCoeff())>1.0e-14)
  {
    const double bftfac = CalcReaBodyForceTermFac(k,phinp,scale);// scalar at integration point np

    bodyforcetermK += ReacCoeff()*Stoich()->at(k)*bftfac;
  }

  return bodyforcetermK;
}

/*----------------------------------------------------------------------/
 | calculate advanced reaction term derivatives             Thon 08/16 |
/----------------------------------------------------------------------*/
double MAT::ScatraReactionMat::CalcReaBodyForceDerivMatrix(
    const int k,                         //!< current scalar id
    const int toderive,                  //!< current id to derive to
    const std::vector<double>& phinp,    //!< scalar values at t_(n+1)
    const double scale                   //!< scaling factor for reference concentrations
    ) const
{
  double reabodyforcederivmatrixKJ=0.0;

  if ( Stoich()->at(k)!=0 and fabs(ReacCoeff())>1.0e-14)
  {
    const double bfdmfac = CalcReaBodyForceDerivFac(k,toderive,phinp,scale);

    reabodyforcederivmatrixKJ += ReacCoeff()*Stoich()->at(k)*bfdmfac;
  }

  return reabodyforcederivmatrixKJ;
}

/*----------------------------------------------------------------------------------*
 |  Modify concentrations according to reacstart vector and do scaling   thon 08/16 |
 *----------------------------------------------------------------------------------*/
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
 |  helper for calculating advanced reaction terms           thon 08/16 |
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

  // TODO: Split the reaction types into classes
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

    case MAT::PAR::reac_coup_byfunction: //reaction by function
    {
      // copy phi vector in different format to be read by the function
      std::vector<std::pair<std::string,double> > variables = BuildPhiVectorForFunction(phinp);
      // evaluate reaction term
      bftfac= Function(round(couprole[k])-1).Evaluate(0,variables);

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
 |  helper for calculating advanced reaction term derivatives          thon 08/16 |
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

  // TODO: Split the reaction types into classes
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

    case MAT::PAR::reac_coup_byfunction: //reaction by function
    {
      // copy phi vector in different format to be read by the function
      std::vector<std::pair<std::string,double> > variables = BuildPhiVectorForFunction(phinp);
      // evaluate the derivatives of the reaction term
      std::vector<std::vector<double> > deriv = Function(round(couprole[k])-1).FctDer(0,variables);

      // the derivative needed is the derivative of the first function (index 0) w.r.t. to
      // the index of the scalar to derive (index toderive)
      bfdmfac = deriv[0][toderive];

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

/*---------------------------------------------------------------------------------/
 | helper for evaluation by function                                     Vuong 08/16 |
/--------------------------------------------------------------------------------- */
std::vector<std::pair<std::string,double> >  MAT::ScatraReactionMat::BuildPhiVectorForFunction(
    const std::vector<double>& phinp    //!< scalar values at t_(n+1)
    ) const
{
  std::vector<std::pair<std::string,double> > variables;
  for (int ii=0; ii < NumScal(); ii++)
  {
    // construct the strings for scalar
    std::ostringstream temp;
    temp << ii+1;
    std::string name = "phi"+temp.str();

    // save the phi values with correct name
    variables.push_back(std::pair<std::string,double>(name,phinp[ii]));
  }
  return variables;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
inline DRT::UTILS::VariableExprFunction& MAT::ScatraReactionMat::Function(int functnum) const
{
  // try to cast to variable expression function for phi evaluation
  // try-catch because of cast of references
  try
  {
    DRT::UTILS::VariableExprFunction& funct =
        dynamic_cast<DRT::UTILS::VariableExprFunction&>(DRT::Problem::Instance()->Funct(functnum));

    return funct;
  }
  catch(std::bad_cast & exp)
  {
    dserror("Cast to VarExp Function failed! For phase law definition only 'VAREXPR' functions are allowed!\n"
        "Check your input file!");
    return dynamic_cast<DRT::UTILS::VariableExprFunction&>(DRT::Problem::Instance()->Funct(functnum));
  }
}
