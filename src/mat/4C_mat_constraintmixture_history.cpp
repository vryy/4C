/*----------------------------------------------------------------------------*/
/*! \file
\brief This file contains the history class of the constraintmixture material

\level 2


*/
/*----------------------------------------------------------------------------*/

#include "4C_mat_constraintmixture_history.hpp"

#include "4C_utils_exceptions.hpp"

FOUR_C_NAMESPACE_OPEN

Mat::ConstraintMixtureHistoryType Mat::ConstraintMixtureHistoryType::instance_;

Core::Communication::ParObject* Mat::ConstraintMixtureHistoryType::Create(
    const std::vector<char>& data)
{
  Mat::ConstraintMixtureHistory* cmhis = new Mat::ConstraintMixtureHistory();
  cmhis->Unpack(data);
  return cmhis;
}

/*----------------------------------------------------------------------*
 |  History: Pack                                 (public)         03/11|
 *----------------------------------------------------------------------*/
void Mat::ConstraintMixtureHistory::Pack(Core::Communication::PackBuffer& data) const
{
  Core::Communication::PackBuffer::SizeMarker sm(data);

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  add_to_pack(data, type);

  // Pack internal variables
  add_to_pack(data, depositiontime_);
  add_to_pack(data, dt_);
  add_to_pack(data, numgp_);
  add_to_pack(data, expvar_);
  for (int gp = 0; gp < numgp_; ++gp)
  {
    add_to_pack(data, collagenstretch1_->at(gp));
    add_to_pack(data, collagenstretch2_->at(gp));
    add_to_pack(data, collagenstretch3_->at(gp));
    add_to_pack(data, collagenstretch4_->at(gp));
    add_to_pack(data, massprod1_->at(gp));
    add_to_pack(data, massprod2_->at(gp));
    add_to_pack(data, massprod3_->at(gp));
    add_to_pack(data, massprod4_->at(gp));
    if (expvar_)
    {
      add_to_pack(data, vardegrad1_->at(gp));
      add_to_pack(data, vardegrad2_->at(gp));
      add_to_pack(data, vardegrad3_->at(gp));
      add_to_pack(data, vardegrad4_->at(gp));
    }
  }

  return;
}

/*----------------------------------------------------------------------*
 |  History: Unpack                               (public)         03/11|
 *----------------------------------------------------------------------*/
void Mat::ConstraintMixtureHistory::Unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;

  Core::Communication::ExtractAndAssertId(position, data, UniqueParObjectId());

  // unpack internal variables
  double a;
  extract_from_pack(position, data, a);
  depositiontime_ = a;
  extract_from_pack(position, data, a);
  dt_ = a;
  int b;
  extract_from_pack(position, data, b);
  numgp_ = b;
  extract_from_pack(position, data, b);
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
    extract_from_pack(position, data, a);
    collagenstretch1_->at(gp) = a;
    extract_from_pack(position, data, a);
    collagenstretch2_->at(gp) = a;
    extract_from_pack(position, data, a);
    collagenstretch3_->at(gp) = a;
    extract_from_pack(position, data, a);
    collagenstretch4_->at(gp) = a;
    extract_from_pack(position, data, a);
    massprod1_->at(gp) = a;
    extract_from_pack(position, data, a);
    massprod2_->at(gp) = a;
    extract_from_pack(position, data, a);
    massprod3_->at(gp) = a;
    extract_from_pack(position, data, a);
    massprod4_->at(gp) = a;
    if (expvar_)
    {
      extract_from_pack(position, data, a);
      vardegrad1_->at(gp) = a;
      extract_from_pack(position, data, a);
      vardegrad2_->at(gp) = a;
      extract_from_pack(position, data, a);
      vardegrad3_->at(gp) = a;
      extract_from_pack(position, data, a);
      vardegrad4_->at(gp) = a;
    }
  }

  if (position != data.size())
    FOUR_C_THROW("Mismatch in size of data %d <-> %d", data.size(), position);

  return;
}

/*----------------------------------------------------------------------*
 |  History: Setup                                (public)         03/11|
 *----------------------------------------------------------------------*/
void Mat::ConstraintMixtureHistory::setup(const int ngp, const double massprodbasal, bool expvar)
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
 |  History: set_stretches                         (private)        03/11|
 *----------------------------------------------------------------------*/
void Mat::ConstraintMixtureHistory::set_stretches(int gp, Core::LinAlg::Matrix<4, 1> stretches)
{
  if (gp < numgp_)
  {
    collagenstretch1_->at(gp) = stretches(0);
    collagenstretch2_->at(gp) = stretches(1);
    collagenstretch3_->at(gp) = stretches(2);
    collagenstretch4_->at(gp) = stretches(3);
  }
  else
    FOUR_C_THROW("gp out of range in set_stretches");
}

/*----------------------------------------------------------------------*
 |  History: get_stretches                         (private)        03/11|
 *----------------------------------------------------------------------*/
void Mat::ConstraintMixtureHistory::get_stretches(int gp, Core::LinAlg::Matrix<4, 1>* stretches)
{
  if (gp < numgp_)
  {
    (*stretches)(0) = collagenstretch1_->at(gp);
    (*stretches)(1) = collagenstretch2_->at(gp);
    (*stretches)(2) = collagenstretch3_->at(gp);
    (*stretches)(3) = collagenstretch4_->at(gp);
  }
  else
    FOUR_C_THROW("gp out of range in get_stretches");
}

/*----------------------------------------------------------------------*
 |  History: set_mass                              (private)        03/11|
 *----------------------------------------------------------------------*/
void Mat::ConstraintMixtureHistory::set_mass(int gp, Core::LinAlg::Matrix<4, 1> massprod)
{
  if (gp < numgp_)
  {
    massprod1_->at(gp) = massprod(0);
    massprod2_->at(gp) = massprod(1);
    massprod3_->at(gp) = massprod(2);
    massprod4_->at(gp) = massprod(3);
  }
  else
    FOUR_C_THROW("gp out of range in set_mass");
}

/*----------------------------------------------------------------------*
 |  History: set_mass                              (private)        04/13|
 *----------------------------------------------------------------------*/
void Mat::ConstraintMixtureHistory::set_mass(int gp, double massprod, int idfiber)
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
      FOUR_C_THROW("no valid fiber id: %d", idfiber);
  }
  else
    FOUR_C_THROW("gp out of range in set_mass");
}

/*----------------------------------------------------------------------*
 |  History: get_mass                              (private)        03/11|
 *----------------------------------------------------------------------*/
void Mat::ConstraintMixtureHistory::get_mass(int gp, Core::LinAlg::Matrix<4, 1>* massprod)
{
  if (gp < numgp_)
  {
    (*massprod)(0) = massprod1_->at(gp);
    (*massprod)(1) = massprod2_->at(gp);
    (*massprod)(2) = massprod3_->at(gp);
    (*massprod)(3) = massprod4_->at(gp);
  }
  else
    FOUR_C_THROW("gp out of range in get_mass");
}

/*----------------------------------------------------------------------*
 |  History: set_var_degrad                         (private)        07/13|
 *----------------------------------------------------------------------*/
void Mat::ConstraintMixtureHistory::set_var_degrad(int gp, int idfiber, double vardegrad)
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
      FOUR_C_THROW("no valid fiber id: %d", idfiber);
  }
  else
    FOUR_C_THROW("gp out of range in set_var_degrad");
}

/*----------------------------------------------------------------------*
 |  History: get_var_degrad                         (private)        07/13|
 *----------------------------------------------------------------------*/
void Mat::ConstraintMixtureHistory::get_var_degrad(int gp, int idfiber, double* vardegrad)
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
      FOUR_C_THROW("no valid fiber id: %d", idfiber);
  }
  else
    FOUR_C_THROW("gp out of range in get_var_degrad");
}

FOUR_C_NAMESPACE_CLOSE
