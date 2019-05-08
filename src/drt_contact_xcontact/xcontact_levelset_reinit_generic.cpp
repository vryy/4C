/*----------------------------------------------------------------------------*/
/**
\file xcontact_levelset_reinit_generic.cpp

\brief xcontact level-set generic reinitialization algorithm

\maintainer Matthias Mayr

\date Dec 1, 2016

\level 3

*/
/*----------------------------------------------------------------------------*/


#include "xcontact_levelset_reinit_generic.H"

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
XCONTACT::LEVELSET::REINIT::Generic::Generic() : isinit_(false), issetup_(false), algorithm_(NULL)
{
  /* left blank */
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void XCONTACT::LEVELSET::REINIT::Generic::Init(XCONTACT::LEVELSET::Algorithm* const algorithm)
{
  issetup_ = false;

  algorithm_ = algorithm;

  isinit_ = true;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void XCONTACT::LEVELSET::REINIT::Generic::Compute(const Epetra_Vector& phinp)
{
  CheckInitSetup();

  PreSolve();

  Solve(phinp);

  PostSolve();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
XCONTACT::LEVELSET::Algorithm& XCONTACT::LEVELSET::REINIT::Generic::Algorithm()
{
  return *algorithm_;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
const XCONTACT::LEVELSET::Algorithm& XCONTACT::LEVELSET::REINIT::Generic::Algorithm() const
{
  return *algorithm_;
}
