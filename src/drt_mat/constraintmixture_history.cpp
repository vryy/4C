/*----------------------------------------------------------------------------*/
/*! \file
\brief This file contains the history class of the constraintmixture material

\level 2

\maintainer Amadeus Gebauer

*/
/*----------------------------------------------------------------------------*/

#include "constraintmixture_history.H"
#include "../drt_lib/drt_dserror.H"

MAT::ConstraintMixtureHistoryType MAT::ConstraintMixtureHistoryType::instance_;

DRT::ParObject* MAT::ConstraintMixtureHistoryType::Create(const std::vector<char>& data)
{
  MAT::ConstraintMixtureHistory* cmhis = new MAT::ConstraintMixtureHistory();
  cmhis->Unpack(data);
  return cmhis;
}

/*----------------------------------------------------------------------*
 |  History: Pack                                 (public)         03/11|
 *----------------------------------------------------------------------*/
void MAT::ConstraintMixtureHistory::Pack(DRT::PackBuffer& data) const
{
  DRT::PackBuffer::SizeMarker sm(data);
  sm.Insert();

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data, type);

  // Pack internal variables
  AddtoPack(data, depositiontime_);
  AddtoPack(data, dt_);
  AddtoPack(data, numgp_);
  AddtoPack(data, expvar_);
  for (int gp = 0; gp < numgp_; ++gp)
  {
    AddtoPack(data, collagenstretch1_->at(gp));
    AddtoPack(data, collagenstretch2_->at(gp));
    AddtoPack(data, collagenstretch3_->at(gp));
    AddtoPack(data, collagenstretch4_->at(gp));
    AddtoPack(data, massprod1_->at(gp));
    AddtoPack(data, massprod2_->at(gp));
    AddtoPack(data, massprod3_->at(gp));
    AddtoPack(data, massprod4_->at(gp));
    if (expvar_)
    {
      AddtoPack(data, vardegrad1_->at(gp));
      AddtoPack(data, vardegrad2_->at(gp));
      AddtoPack(data, vardegrad3_->at(gp));
      AddtoPack(data, vardegrad4_->at(gp));
    }
  }

  return;
}

/*----------------------------------------------------------------------*
 |  History: Unpack                               (public)         03/11|
 *----------------------------------------------------------------------*/
void MAT::ConstraintMixtureHistory::Unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;
  // extract type
  int type = 0;
  ExtractfromPack(position, data, type);
  if (type != UniqueParObjectId()) dserror("wrong instance type data");

  // unpack internal variables
  double a;
  ExtractfromPack(position, data, a);
  depositiontime_ = a;
  ExtractfromPack(position, data, a);
  dt_ = a;
  int b;
  ExtractfromPack(position, data, b);
  numgp_ = b;
  ExtractfromPack(position, data, b);
  expvar_ = b;

  collagenstretch1_ = Teuchos::rcp(new std::vector<double>(numgp_));
  collagenstretch2_ = Teuchos::rcp(new std::vector<double>(numgp_));
  collagenstretch3_ = Teuchos::rcp(new std::vector<double>(numgp_));
  collagenstretch4_ = Teuchos::rcp(new std::vector<double>(numgp_));
  massprod1_ = Teuchos::rcp(new std::vector<double>(numgp_));
  massprod2_ = Teuchos::rcp(new std::vector<double>(numgp_));
  massprod3_ = Teuchos::rcp(new std::vector<double>(numgp_));
  massprod4_ = Teuchos::rcp(new std::vector<double>(numgp_));
  if (expvar_)
  {
    vardegrad1_ = Teuchos::rcp(new std::vector<double>(numgp_));
    vardegrad2_ = Teuchos::rcp(new std::vector<double>(numgp_));
    vardegrad3_ = Teuchos::rcp(new std::vector<double>(numgp_));
    vardegrad4_ = Teuchos::rcp(new std::vector<double>(numgp_));
  }

  for (int gp = 0; gp < numgp_; ++gp)
  {
    ExtractfromPack(position, data, a);
    collagenstretch1_->at(gp) = a;
    ExtractfromPack(position, data, a);
    collagenstretch2_->at(gp) = a;
    ExtractfromPack(position, data, a);
    collagenstretch3_->at(gp) = a;
    ExtractfromPack(position, data, a);
    collagenstretch4_->at(gp) = a;
    ExtractfromPack(position, data, a);
    massprod1_->at(gp) = a;
    ExtractfromPack(position, data, a);
    massprod2_->at(gp) = a;
    ExtractfromPack(position, data, a);
    massprod3_->at(gp) = a;
    ExtractfromPack(position, data, a);
    massprod4_->at(gp) = a;
    if (expvar_)
    {
      ExtractfromPack(position, data, a);
      vardegrad1_->at(gp) = a;
      ExtractfromPack(position, data, a);
      vardegrad2_->at(gp) = a;
      ExtractfromPack(position, data, a);
      vardegrad3_->at(gp) = a;
      ExtractfromPack(position, data, a);
      vardegrad4_->at(gp) = a;
    }
  }

  if (position != data.size()) dserror("Mismatch in size of data %d <-> %d", data.size(), position);

  return;
}

/*----------------------------------------------------------------------*
 |  History: Setup                                (public)         03/11|
 *----------------------------------------------------------------------*/
void MAT::ConstraintMixtureHistory::Setup(const int ngp, const double massprodbasal, bool expvar)
{
  dt_ = 0.0;
  depositiontime_ = 0.0;

  numgp_ = ngp;
  expvar_ = expvar;
  // history variables
  collagenstretch1_ = Teuchos::rcp(new std::vector<double>(numgp_));
  collagenstretch2_ = Teuchos::rcp(new std::vector<double>(numgp_));
  collagenstretch3_ = Teuchos::rcp(new std::vector<double>(numgp_));
  collagenstretch4_ = Teuchos::rcp(new std::vector<double>(numgp_));
  massprod1_ = Teuchos::rcp(new std::vector<double>(numgp_));
  massprod2_ = Teuchos::rcp(new std::vector<double>(numgp_));
  massprod3_ = Teuchos::rcp(new std::vector<double>(numgp_));
  massprod4_ = Teuchos::rcp(new std::vector<double>(numgp_));
  if (expvar_)
  {
    vardegrad1_ = Teuchos::rcp(new std::vector<double>(numgp_));
    vardegrad2_ = Teuchos::rcp(new std::vector<double>(numgp_));
    vardegrad3_ = Teuchos::rcp(new std::vector<double>(numgp_));
    vardegrad4_ = Teuchos::rcp(new std::vector<double>(numgp_));
  }

  for (int gp = 0; gp < numgp_; gp++)
  {
    collagenstretch1_->at(gp) = 1.0;
    collagenstretch2_->at(gp) = 1.0;
    collagenstretch3_->at(gp) = 1.0;
    collagenstretch4_->at(gp) = 1.0;
    massprod1_->at(gp) = massprodbasal;
    massprod2_->at(gp) = massprodbasal;
    massprod3_->at(gp) = massprodbasal;  //*4.;
    massprod4_->at(gp) = massprodbasal;  //*4.;
    if (expvar_)
    {
      vardegrad1_->at(gp) = 1.0;
      vardegrad2_->at(gp) = 1.0;
      vardegrad3_->at(gp) = 1.0;
      vardegrad4_->at(gp) = 1.0;
    }
  }
}

/*----------------------------------------------------------------------*
 |  History: SetStretches                         (private)        03/11|
 *----------------------------------------------------------------------*/
void MAT::ConstraintMixtureHistory::SetStretches(int gp, LINALG::Matrix<4, 1> stretches)
{
  if (gp < numgp_)
  {
    collagenstretch1_->at(gp) = stretches(0);
    collagenstretch2_->at(gp) = stretches(1);
    collagenstretch3_->at(gp) = stretches(2);
    collagenstretch4_->at(gp) = stretches(3);
  }
  else
    dserror("gp out of range in SetStretches");
}

/*----------------------------------------------------------------------*
 |  History: GetStretches                         (private)        03/11|
 *----------------------------------------------------------------------*/
void MAT::ConstraintMixtureHistory::GetStretches(int gp, LINALG::Matrix<4, 1>* stretches)
{
  if (gp < numgp_)
  {
    (*stretches)(0) = collagenstretch1_->at(gp);
    (*stretches)(1) = collagenstretch2_->at(gp);
    (*stretches)(2) = collagenstretch3_->at(gp);
    (*stretches)(3) = collagenstretch4_->at(gp);
  }
  else
    dserror("gp out of range in GetStretches");
}

/*----------------------------------------------------------------------*
 |  History: SetMass                              (private)        03/11|
 *----------------------------------------------------------------------*/
void MAT::ConstraintMixtureHistory::SetMass(int gp, LINALG::Matrix<4, 1> massprod)
{
  if (gp < numgp_)
  {
    massprod1_->at(gp) = massprod(0);
    massprod2_->at(gp) = massprod(1);
    massprod3_->at(gp) = massprod(2);
    massprod4_->at(gp) = massprod(3);
  }
  else
    dserror("gp out of range in SetMass");
}

/*----------------------------------------------------------------------*
 |  History: SetMass                              (private)        04/13|
 *----------------------------------------------------------------------*/
void MAT::ConstraintMixtureHistory::SetMass(int gp, double massprod, int idfiber)
{
  if (gp < numgp_)
  {
    if (idfiber == 0)
    {
      massprod1_->at(gp) = massprod;
    }
    else if (idfiber == 1)
    {
      massprod2_->at(gp) = massprod;
    }
    else if (idfiber == 2)
    {
      massprod3_->at(gp) = massprod;
    }
    else if (idfiber == 3)
    {
      massprod4_->at(gp) = massprod;
    }
    else
      dserror("no valid fiber id: %d", idfiber);
  }
  else
    dserror("gp out of range in SetMass");
}

/*----------------------------------------------------------------------*
 |  History: GetMass                              (private)        03/11|
 *----------------------------------------------------------------------*/
void MAT::ConstraintMixtureHistory::GetMass(int gp, LINALG::Matrix<4, 1>* massprod)
{
  if (gp < numgp_)
  {
    (*massprod)(0) = massprod1_->at(gp);
    (*massprod)(1) = massprod2_->at(gp);
    (*massprod)(2) = massprod3_->at(gp);
    (*massprod)(3) = massprod4_->at(gp);
  }
  else
    dserror("gp out of range in GetMass");
}

/*----------------------------------------------------------------------*
 |  History: SetVarDegrad                         (private)        07/13|
 *----------------------------------------------------------------------*/
void MAT::ConstraintMixtureHistory::SetVarDegrad(int gp, int idfiber, double vardegrad)
{
  if (gp < numgp_)
  {
    if (idfiber == 0)
    {
      vardegrad1_->at(gp) = vardegrad;
    }
    else if (idfiber == 1)
    {
      vardegrad2_->at(gp) = vardegrad;
    }
    else if (idfiber == 2)
    {
      vardegrad3_->at(gp) = vardegrad;
    }
    else if (idfiber == 3)
    {
      vardegrad4_->at(gp) = vardegrad;
    }
    else
      dserror("no valid fiber id: %d", idfiber);
  }
  else
    dserror("gp out of range in SetVarDegrad");
}

/*----------------------------------------------------------------------*
 |  History: GetVarDegrad                         (private)        07/13|
 *----------------------------------------------------------------------*/
void MAT::ConstraintMixtureHistory::GetVarDegrad(int gp, int idfiber, double* vardegrad)
{
  if (gp < numgp_)
  {
    if (idfiber == 0)
    {
      *vardegrad = vardegrad1_->at(gp);
    }
    else if (idfiber == 1)
    {
      *vardegrad = vardegrad2_->at(gp);
    }
    else if (idfiber == 2)
    {
      *vardegrad = vardegrad3_->at(gp);
    }
    else if (idfiber == 3)
    {
      *vardegrad = vardegrad4_->at(gp);
    }
    else
      dserror("no valid fiber id: %d", idfiber);
  }
  else
    dserror("gp out of range in GetVarDegrad");
}
