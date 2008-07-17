/*----------------------------------------------------------------------*/
/*!
\file strutimint_vector.cpp
\brief Dummy object holding vector utility functions

<pre>
Maintainer: Burkhard Bornemann
            bornemann@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15237
</pre>
*/

/*----------------------------------------------------------------------*/
/* definitions */
#ifdef CCADISCRET

/*----------------------------------------------------------------------*/
/* headers */
#include "strutimint_vector.H"

/*----------------------------------------------------------------------*/
/* Map vector norm identifaction string to enum term */
enum STR::StruTimIntVector::NormEnum STR::StruTimIntVector::MapNormStringToEnum
(
  const std::string name  //!< name identification string
)
{
  if (name == "Vague")
  {
    return norm_vague;
  }
  else if (name == "L1")
  {
    return norm_l1;
  }
  else if (name == "L2")
  {
    return norm_l2;
  }
  else if (name == "Rms")
  {
    return norm_rms;
  }
  else if (name == "Inf")
  {
    return norm_inf;
  }
  else
  {
    dserror("Unknown kind of predictor %s", name.c_str());
    return norm_vague;
  }
}

/*----------------------------------------------------------------------*/
/* map enum term to string */
std::string STR::StruTimIntVector::MapNormEnumToString
(
  const enum NormEnum norm  //!< input enum term
)
{
  // holds results
  std::string str;

  // select it
  switch (norm)
  {
  case norm_vague:
    str = "Vague";
    break;
  case norm_l1:
    str = "L1";
    break;
  case norm_l2:
    str = "L2";
    break;
  case norm_rms:
    str = "Rms";
    break;
  case norm_inf:
    str = "Inf";
    break;
  default:
    str = "norm is undefined";
    break;
  }

  // send to hell
  return str;
}

/*----------------------------------------------------------------------*/
/* Calculate vector norm */
double STR::StruTimIntVector::CalculateNorm
(
  const enum NormEnum norm,
  const RCP<Epetra_Vector> vect
)
{
  // average norm
  if (norm == norm_l1)
  {
    double vectnorm;
    vect->Norm1(&vectnorm);
    return vectnorm;
  }
  // quadratic norm
  else if (norm == norm_l2)
  {
    double vectnorm;
    vect->Norm2(&vectnorm);
    return vectnorm;
  }
  // RMS norm
  else if (norm == norm_rms)
  {
    double vectnorm;
    vect->Norm2(&vectnorm);
    return vectnorm/sqrt((double) vect->GlobalLength());
  }
  // infinity/maximum norm
  else if (norm == norm_inf)
  {
    double vectnorm;
    vect->NormInf(&vectnorm);
    return vectnorm;
  }
  else
  {
    dserror("Cannot handle vector norm");
    return NaN;
  }
}

/*----------------------------------------------------------------------*/
#endif  // #ifdef CCADISCRET
