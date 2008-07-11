/*----------------------------------------------------------------------*/
/*!
\file strutimint_vector.H

<pre>
Maintainer: Burkhard Bornemann
            bornemann@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15237
</pre>
*/

/*----------------------------------------------------------------------*/
#ifdef CCADISCRET
#ifndef STRUTIMINT_VECTOR_H
#define STRUTIMINT_VECTOR_H

/*----------------------------------------------------------------------*/
/* headers */
#include "strutimint.H"

/*======================================================================*/
/*!
\brief Dummy object holding vector utility functions
\author bborn
\date 07/08
*/
class StruTimIntVector
{

public:

  //! friendship with StruTimInt
  friend class StruTimInt;

  //! Kind of vector norms available for convergence check
  enum NormEnum
  {
    norm_vague,  //!< undetermined
    norm_l1,  //!< L1/linear/average norm
    norm_l2,  //!< L2/Euclidean norm
    norm_rms,  //!< RMS/weighted-L2 norm
    norm_inf  //!< Maximum/infinity norm
  };

  //! Map vector norm identifaction string to enum term
  static enum NormEnum MapNormStringToEnum
  (
    const std::string name  //!< name identification string
  );

  //! Determine norm of force residual
  static double CalculateNorm
  (
    const enum StruTimIntVector::NormEnum norm,  //!< norm to use, see #StruTimIntVector::NormEnum
    const RCP<Epetra_Vector> vect  //!< the vector
  );

private:

  //! Constructor
  StruTimIntVector() { ; }

  //! Destructor
  ~StruTimIntVector() { ; }

};

/*----------------------------------------------------------------------*/
#endif  // #ifndef STRUTIMINT_VECTOR_H
#endif  // #ifdef CCADISCRET
