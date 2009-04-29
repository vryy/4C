/*----------------------------------------------------------------------*/
/*!
\file str_utils.cpp
\brief various auxiliar methods needed in structural analysis

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
#include "str_aux.H"

/*----------------------------------------------------------------------*/
/* Calculate vector norm */
double STR::AUX::CalculateVectorNorm
(
  const enum INPAR::STR::VectorNorm norm,
  const Teuchos::RCP<Epetra_Vector> vect
)
{
  // L1 norm
  if (norm == INPAR::STR::norm_l1)
  {
    double vectnorm;
    vect->Norm1(&vectnorm);
    return vectnorm;
  }
  // L2/Euclidian norm
  else if (norm == INPAR::STR::norm_l2)
  {
    double vectnorm;
    vect->Norm2(&vectnorm);
    return vectnorm;
  }
  // RMS norm
  else if (norm == INPAR::STR::norm_rms)
  {
    double vectnorm;
    vect->Norm2(&vectnorm);
    return vectnorm/sqrt((double) vect->GlobalLength());
  }
  // infinity/maximum norm
  else if (norm == INPAR::STR::norm_inf)
  {
    double vectnorm;
    vect->NormInf(&vectnorm);
    return vectnorm;
  }
  else
  {
    dserror("Cannot handle vector norm");
    return 0;
  }
}

/*----------------------------------------------------------------------*/
#endif  // #ifdef CCADISCRET
