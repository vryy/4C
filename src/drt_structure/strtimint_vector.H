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

#include "../drt_inpar/inpar_structure.H"

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

/*
    //! Map vector norm identifaction string to enum term
    static enum INPAR::STR::VectorNorm MapNormStringToEnum
    (
      const std::string name  //!< name identification string
    );
*/

/*
    //! Map enum term to string
    static std::string MapNormEnumToString
    (
      const enum INPAR::STR::VectorNorm norm  //!< input enum term
    );
*/

    //! Determine norm of force residual
    static double CalculateNorm
    (
      const enum INPAR::STR::VectorNorm norm,  //!< norm to use
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
