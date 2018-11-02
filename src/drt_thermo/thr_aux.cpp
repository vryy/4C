/*----------------------------------------------------------------------*/
/*!
\file thr_aux.cpp
\brief various auxiliar methods needed in thermal analysis
\level 1
\maintainer Christoph Meier
*/


/*----------------------------------------------------------------------*
 | headers                                                  bborn 08/09 |
 *----------------------------------------------------------------------*/
#include "thr_aux.H"


/*----------------------------------------------------------------------*
 | Calculate vector norm                                    bborn 08/09 |
 *----------------------------------------------------------------------*/
double THR::AUX::CalculateVectorNorm(
    const enum INPAR::THR::VectorNorm norm, const Teuchos::RCP<Epetra_Vector> vect)
{
  // L1 norm
  if (norm == INPAR::THR::norm_l1)
  {
    double vectnorm;
    vect->Norm1(&vectnorm);
    return vectnorm;
  }
  // L2/Euclidian norm
  else if (norm == INPAR::THR::norm_l2)
  {
    double vectnorm;
    vect->Norm2(&vectnorm);
    return vectnorm;
  }
  // RMS norm
  else if (norm == INPAR::THR::norm_rms)
  {
    double vectnorm;
    vect->Norm2(&vectnorm);
    return vectnorm / std::sqrt((double)vect->GlobalLength());
  }
  // infinity/maximum norm
  else if (norm == INPAR::THR::norm_inf)
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
}  // CalculateVectorNorm()


/*----------------------------------------------------------------------*/
