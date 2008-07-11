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
#ifdef CCADISCRET

/*----------------------------------------------------------------------*/
/* headers */
#include "strutimint_vector.H"

/*----------------------------------------------------------------------*/
/* Map vector norm identifaction string to enum term */
enum StruTimIntVector::NormEnum StruTimIntVector::MapNormStringToEnum
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
/* Calculate vector norm */
double StruTimIntVector::CalculateNorm
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
    return -1;
  }
}

/*----------------------------------------------------------------------*/
#endif  // #ifdef CCADISCRET
