/*----------------------------------------------------------------------*/
/*!
\file strtimint_vector.H

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
#ifndef STRTIMINT_VECTOR_H
#define STRTIMINT_VECTOR_H

/*----------------------------------------------------------------------*/
/* headers */
#include <string>

#include "strtimint.H"

/*----------------------------------------------------------------------*/
/* belongs to structural dynamics namespace */
namespace STR
{

  /*====================================================================*/
  /*!
   * \brief Dummy object holding vector utility functions
   * \author bborn
   * \date 07/08
   */
  class TimIntVector
  {

    //! friendship with TimInt
    friend class TimInt;

  public:

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

    //! Map enum term to string
    static std::string MapNormEnumToString
    (
      const enum NormEnum norm  //!< input enum term
    );

    //! Determine norm of force residual
    static double CalculateNorm
    (
      const enum TimIntVector::NormEnum norm,  //!< norm to use
      const Teuchos::RCP<Epetra_Vector> vect  //!< the vector
    );

  private:

    //! Constructor
    TimIntVector() { ; }

    //! Destructor
    ~TimIntVector() { ; }

  };

}  // namespace STR

/*----------------------------------------------------------------------*/
#endif  // #ifndef STRTIMINT_VECTOR_H
#endif  // #ifdef CCADISCRET
