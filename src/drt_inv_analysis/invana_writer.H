/*----------------------------------------------------------------------------*/
/*!
\file invana_writer.H

\brief Output for the inverse analysis control routine

<pre>
\level 3
\maintainer Sebastian Brandstaeter
            brandstaeter@lnm.mw.tum.de
            089 - 289-15276
</pre>
*/

/*----------------------------------------------------------------------------*/

#ifndef INVANA_WRITER_H
#define INVANA_WRITER_H

/*----------------------------------------------------------------------------*/
/* headers */
#include <Teuchos_RCP.hpp>
#include "Epetra_MultiVector.h"

#include "../drt_io/io.H"

/*----------------------------------------------------------------------------*/
/* forward declarations */

/*----------------------------------------------------------------------------*/
namespace INVANA
{
  /*! \class InvanaWriter
   *  \brief Class to write solver output to restart from
   *
   *  This writer is to be passed to the INVANA::Optimizer to be used
   *  as output writer.
   *
   *  \author kehl \date 12/2015
   */
  class InvanaWriter
  {
   public:
    //! @name Construction/destruction
    //@{

    //! Constructor (empty)
    InvanaWriter();

    //! Destructor
    virtual ~InvanaWriter() {}

    //@}

    //! @name Initialization
    //@{

    /*! \brief Initialization
     *
     * TODO: Documentation (kehl)
     *
     *  \author kehl \date 12/2015
     */
    void Init(Teuchos::RCP<IO::DiscretizationWriter> writer);

    //@}

    //! Write methods
    //@{

    //! Write a new step
    virtual void WriteNewStep(const int step,  /// [in] new step to be written
        const double time                      /// [in] time for this step
    );

    //! Write named int
    virtual void WriteNamedInt(const std::string name,  /// [in] an name for the int
        const int aint                                  ///[int] integer to write
    );

    //! Write named int
    virtual void WriteNamedDouble(const std::string name,  /// [in] an name for the d
        const double adouble                               ///[int] double to write
    );

    //! Write named vector
    virtual void WriteNamedVector(const std::string name,  ///< [in] a name for the vector
        Teuchos::RCP<const Epetra_MultiVector> vector,     ///< [in] the vector to write
        IO::VectorType type = IO::elementvector            ///<[in] vector type
        ) const;

    //! Write named vectors (each vector in MV with a number)
    virtual void WriteNamedVectors(const std::string name,  ///< [in] a name for the vector
        Teuchos::RCP<const Epetra_MultiVector> vector,      ///< [in] the vector to write
        IO::VectorType type = IO::elementvector             ///<[in] vector type
        ) const;

    //! Write mesh
    virtual void WriteMesh(int step, double time);

   private:
    //! output writer
    Teuchos::RCP<IO::DiscretizationWriter> output_;

  };  // class InvanaWriter
}  // namespace INVANA

#endif /* INVANA_WRITER_H */
