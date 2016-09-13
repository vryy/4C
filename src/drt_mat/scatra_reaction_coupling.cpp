/*----------------------------------------------------------------------*/
/*!
 \file scatra_reaction_coupling.cpp

 \brief helper class encapsulating the reaction terms and its derivatives

   \level 3

   \maintainer  Anh-Tu Vuong
                vuong@lnm.mw.tum.de
                http://www.lnm.mw.tum.de
                089 - 289-15251
 *----------------------------------------------------------------------*/

#include "scatra_reaction_coupling.H"

#include "../drt_lib/drt_globalproblem.H"

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void MAT::PAR::REACTIONCOUPLING::SimpleMultiplicative::Initialize(
    int numscal,                         //!< number of scalars
    const std::vector<double>& couprole //!< coupling role vector
    )
{
  ReactionBase::Initialize(numscal,couprole);
}

/*----------------------------------------------------------------------*
 |  helper for calculating advanced reaction terms           thon 08/16 |
 *----------------------------------------------------------------------*/
double MAT::PAR::REACTIONCOUPLING::SimpleMultiplicative::CalcReaBodyForceTerm(
    const int k,                         //!< current scalar id
    int numscal,                         //!< number of scalars
    const std::vector<double>& phinp,    //!< scalar values at t_(n+1)
    const std::vector<double>& couprole, //!< coupling role vector
    const double scale_reac              //!< scaling factor for reaction term (= reaction coefficient * stoichometry)
    )
{

  double bftfac=1.0;

  for (int ii=0; ii < numscal; ii++)
  {
    if (couprole[ii]!=0)
    {
      bftfac *=phinp[ii];
    }
  }

  return scale_reac*bftfac;
}

/*--------------------------------------------------------------------------------*
 |  helper for calculating advanced reaction term derivatives          thon 08/16 |
 *--------------------------------------------------------------------------------*/
void MAT::PAR::REACTIONCOUPLING::SimpleMultiplicative::CalcReaBodyForceDeriv(
    const int k,                         //!< current scalar id
    int numscal,                         //!< number of scalars
    std::vector<double>& derivs,         //!< vector with derivatives (to be filled)
    const std::vector<double>& phinp,    //!< scalar values at t_(n+1)
    const std::vector<double>& couprole, //!< coupling role vector
    const std::vector<double>& reacstart,//!< reaction start vector
    const double scale_reac              //!< scaling factor for reaction term (= reaction coefficient * stoichometry)
    )
{

  for (int toderive=0; toderive<numscal ;toderive++)
  {
    double bfdmfac=1.0;
    if (couprole[toderive]!=0)
    {
      for (int ii=0; ii < numscal; ii++)
      {
        if (couprole[ii]!=0 and ii!=toderive)
          bfdmfac *= phinp[ii];
      }
    }
    else
      bfdmfac=0.0;

    if ( not(reacstart[toderive] > 0 and phinp[toderive]==0.0) )
      derivs[toderive]+=scale_reac*bfdmfac;
  }

  return;
}

/*----------------------------------------------------------------------*
 * **********************************************************************
 *----------------------------------------------------------------------*/


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void MAT::PAR::REACTIONCOUPLING::PowerMultiplicative::Initialize(
    int numscal,                         //!< number of scalars
    const std::vector<double>& couprole  //!< coupling role vector
    )
{
  ReactionBase::Initialize(numscal,couprole);
}

/*----------------------------------------------------------------------*
 |  helper for calculating advanced reaction terms           thon 08/16 |
 *----------------------------------------------------------------------*/
double MAT::PAR::REACTIONCOUPLING::PowerMultiplicative::CalcReaBodyForceTerm(
    const int k,                         //!< current scalar id
    int numscal,                         //!< number of scalars
    const std::vector<double>& phinp,    //!< scalar values at t_(n+1)
    const std::vector<double>& couprole, //!< coupling role vector
    const double scale_reac              //!< scaling factor for reaction term (= reaction coefficient * stoichometry)
    )
{

  double bftfac=1.0;

  for (int ii=0; ii < numscal; ii++)
  {
    if (couprole[ii]!=0)
    {
      bftfac *= std::pow(phinp[ii],couprole[ii]);
    }
  }

  return scale_reac*bftfac;
}

/*--------------------------------------------------------------------------------*
 |  helper for calculating advanced reaction term derivatives          thon 08/16 |
 *--------------------------------------------------------------------------------*/
void MAT::PAR::REACTIONCOUPLING::PowerMultiplicative::CalcReaBodyForceDeriv(
    const int k,                         //!< current scalar id
    int numscal,                         //!< number of scalars
    std::vector<double>& derivs,         //!< vector with derivatives (to be filled)
    const std::vector<double>& phinp,    //!< scalar values at t_(n+1)
    const std::vector<double>& couprole, //!< coupling role vector
    const std::vector<double>& reacstart,//!< reaction start vector
    const double scale_reac              //!< scaling factor for reaction term (= reaction coefficient * stoichometry)
    )
{

  for (int toderive=0; toderive<numscal ;toderive++)
  {
    double bfdmfac=1.0;
    if (couprole[toderive]!=0)
      {
        for (int ii=0; ii < numscal; ii++)
        {
          if (couprole[ii]!=0 and ii!=toderive)
            bfdmfac *= std::pow(phinp[ii],couprole[ii]);
          else if(couprole[ii]!=0 and ii==toderive)
            bfdmfac *= std::pow(phinp[ii],couprole[ii]-1.0);
        }
      }
    else
      bfdmfac=0.0;

    if ( not(reacstart[toderive] > 0 and phinp[toderive]==0.0) )
      derivs[toderive]+=scale_reac*bfdmfac;
  }

  return;
}

/*----------------------------------------------------------------------*
 * **********************************************************************
 *----------------------------------------------------------------------*/


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void MAT::PAR::REACTIONCOUPLING::Constant::Initialize(
    int numscal,                         //!< number of scalars
    const std::vector<double>& couprole  //!< coupling role vector
    )
{
  ReactionBase::Initialize(numscal,couprole);
}

/*----------------------------------------------------------------------*
 |  helper for calculating advanced reaction terms           thon 08/16 |
 *----------------------------------------------------------------------*/
double MAT::PAR::REACTIONCOUPLING::Constant::CalcReaBodyForceTerm(
    const int k,                         //!< current scalar id
    int numscal,                         //!< number of scalars
    const std::vector<double>& phinp,    //!< scalar values at t_(n+1)
    const std::vector<double>& couprole, //!< coupling role vector
    const double scale_reac              //!< scaling factor for reaction term (= reaction coefficient * stoichometry)
    )
{
  return scale_reac;
}

/*--------------------------------------------------------------------------------*
 |  helper for calculating advanced reaction term derivatives          thon 08/16 |
 *--------------------------------------------------------------------------------*/
void MAT::PAR::REACTIONCOUPLING::Constant::CalcReaBodyForceDeriv(
    const int k,                         //!< current scalar id
    int numscal,                         //!< number of scalars
    std::vector<double>& derivs,         //!< vector with derivatives (to be filled)
    const std::vector<double>& phinp,    //!< scalar values at t_(n+1)
    const std::vector<double>& couprole, //!< coupling role vector
    const std::vector<double>& reacstart,//!< reaction start vector
    const double scale_reac              //!< scaling factor for reaction term (= reaction coefficient * stoichometry)
    )
{
  // zero derivative -> do nothing
  return;
}

/*----------------------------------------------------------------------*
 * **********************************************************************
 *----------------------------------------------------------------------*/


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void MAT::PAR::REACTIONCOUPLING::MichaelisMenten::Initialize(
    int numscal,                         //!< number of scalars
    const std::vector<double>& couprole  //!< coupling role vector
    )
{
  ReactionBase::Initialize(numscal,couprole);
}

/*----------------------------------------------------------------------*
 |  helper for calculating advanced reaction terms           thon 08/16 |
 *----------------------------------------------------------------------*/
double MAT::PAR::REACTIONCOUPLING::MichaelisMenten::CalcReaBodyForceTerm(
    const int k,                         //!< current scalar id
    int numscal,                         //!< number of scalars
    const std::vector<double>& phinp,    //!< scalar values at t_(n+1)
    const std::vector<double>& couprole, //!< coupling role vector
    const double scale_reac              //!< scaling factor for reaction term (= reaction coefficient * stoichometry)
    )
{

  double bftfac=1.0;

  for (int ii=0; ii < numscal ; ii++)
  {
    if (couprole[ii]>0.0) //and (ii!=k))
      bftfac *= phinp[ii]/(couprole[ii]+phinp[ii]);
    else if (couprole[ii]<0.0) //and (ii!=k))
      bftfac *= phinp[ii];
  }

  return scale_reac*bftfac;
}

/*--------------------------------------------------------------------------------*
 |  helper for calculating advanced reaction term derivatives          thon 08/16 |
 *--------------------------------------------------------------------------------*/
void MAT::PAR::REACTIONCOUPLING::MichaelisMenten::CalcReaBodyForceDeriv(
    const int k,                         //!< current scalar id
    int numscal,                         //!< number of scalars
    std::vector<double>& derivs,         //!< vector with derivatives (to be filled)
    const std::vector<double>& phinp,    //!< scalar values at t_(n+1)
    const std::vector<double>& couprole, //!< coupling role vector
    const std::vector<double>& reacstart,//!< reaction start vector
    const double scale_reac              //!< scaling factor for reaction term (= reaction coefficient * stoichometry)
    )
{

  for (int toderive=0; toderive<numscal ;toderive++)
  {
    double bfdmfac=1.0;
    for (int ii=0; ii < numscal; ii++)
    {
      if (ii != toderive)
      {
        if (couprole[ii] > 0.0)
          bfdmfac *= phinp.at(ii)/(couprole[ii] + phinp[ii]);
        else if (couprole[ii] < 0.0)
          bfdmfac *= phinp.at(ii);
        else
          bfdmfac *= 1;
      }
      else
      {
        if (couprole[ii] > 0.0)
          bfdmfac *= couprole[ii]/(std::pow((phinp[ii]+couprole[ii]), 2));
        else if (couprole[ii] < 0.0)
          bfdmfac *= 1;
        else
          bfdmfac = 0;
      }
    }
    if ( not(reacstart[toderive] > 0 and phinp[toderive]==0.0) )
      derivs[toderive]+=scale_reac*bfdmfac;
  }

  return;
}


/*----------------------------------------------------------------------*
 * **********************************************************************
 *----------------------------------------------------------------------*/


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void MAT::PAR::REACTIONCOUPLING::ByFunction::Initialize(
    int numscal,                         //!< number of scalars
    const std::vector<double>& couprole  //!< coupling role vector
    )
{
  if(not IsInit())
  {
    variables_.reserve(numscal);
    for(int ii=0;ii<numscal;ii++)
    {
      // we take the value in couprole list as function ID
      const int functID = round(couprole[ii]);
      if (functID!=0)
      {
        if(Function(functID-1).NumberComponents()!=1)
          dserror("expected only one component for the reaction evaluation");

        for(int k=0;k<numscal;k++)
        {
          // construct the strings for scalar
          std::ostringstream temp;
          temp << k+1;
          std::string name = "phi"+temp.str();

          // add the variable name to the parser
          if(not Function(functID-1).IsVariable(0,name))
          {
            Function(functID-1).AddVariable(0,name,0.0);

            // save the phi values with correct name
            variables_.push_back(std::pair<std::string,double>(name,0.0));
          }
        }
      }
    }
  }

  ReactionBase::Initialize(numscal,couprole);
}

/*----------------------------------------------------------------------*
 |  helper for calculating advanced reaction terms           thon 08/16 |
 *----------------------------------------------------------------------*/
double MAT::PAR::REACTIONCOUPLING::ByFunction::CalcReaBodyForceTerm(
    const int k,                         //!< current scalar id
    int numscal,                         //!< number of scalars
    const std::vector<double>& phinp,    //!< scalar values at t_(n+1)
    const std::vector<double>& couprole, //!< coupling role vector
    const double scale_reac              //!< scaling factor for reaction term (= reaction coefficient * stoichometry)
    )
{

  // copy phi vector in different format to be read by the function
  BuildPhiVectorForFunction(phinp,numscal);
  // evaluate reaction term
  double bftfac= Function(round(couprole[k])-1).Evaluate(0,variables_);

  return scale_reac*bftfac;
}

/*--------------------------------------------------------------------------------*
 |  helper for calculating advanced reaction term derivatives          thon 08/16 |
 *--------------------------------------------------------------------------------*/
void MAT::PAR::REACTIONCOUPLING::ByFunction::CalcReaBodyForceDeriv(
    const int k,                         //!< current scalar id
    int numscal,                         //!< number of scalars
    std::vector<double>& derivs,         //!< vector with derivatives (to be filled)
    const std::vector<double>& phinp,    //!< scalar values at t_(n+1)
    const std::vector<double>& couprole, //!< coupling role vector
    const std::vector<double>& reacstart,//!< reaction start vector
    const double scale_reac              //!< scaling factor for reaction term (= reaction coefficient * stoichometry)
    )
{

  // copy phi vector in different format to be read by the function
  BuildPhiVectorForFunction(phinp,numscal);
  // evaluate the derivatives of the reaction term
  std::vector<double> myderiv = Function(round(couprole[k])-1).FctDer(0,variables_);

  for (int toderive=0; toderive<numscal ;toderive++)
  {
    // the derivative needed is the derivative of the first function (index 0) w.r.t. to
    // the index of the scalar to derive (index toderive)
    double bfdmfac = myderiv[toderive];

    if ( not(reacstart[toderive] > 0 and phinp[toderive]==0.0) )
      derivs[toderive]+=scale_reac*bfdmfac;
  }

  return;
}

/*---------------------------------------------------------------------------------/
 | helper for evaluation by function                                     vuong 08/16 |
/--------------------------------------------------------------------------------- */
void MAT::PAR::REACTIONCOUPLING::ByFunction::BuildPhiVectorForFunction(
    const std::vector<double>& phinp,    //!< scalar values at t_(n+1)
    int numscal                          //!< number of scalars
    )
{
  // note: we use the fact that the 'variables_' vector is ordered in the same way
  //       as the phi vector!
  for (int ii=0; ii < numscal; ii++)
  {
    variables_[ii].second = phinp[ii];
  }
  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
inline DRT::UTILS::VariableExprFunction& MAT::PAR::REACTIONCOUPLING::ByFunction::Function(int functnum) const
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
