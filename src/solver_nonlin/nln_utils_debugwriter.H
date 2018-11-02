/*----------------------------------------------------------------------------*/
/*!
\file nln_utils_debugwriter.H

\brief Write output for debugging purposes

\maintainer Matthias Mayr

\level 3
*/

/*----------------------------------------------------------------------------*/

#ifndef NLN_UTILS_DEBUGWRITER_H
#define NLN_UTILS_DEBUGWRITER_H

/*----------------------------------------------------------------------------*/
/* headers */

// std
#include <string.h>

// Teuchos
#include <Teuchos_RCP.hpp>
#include <Teuchos_VerboseObject.hpp>

// baci
#include "../drt_io/io.H"

/*----------------------------------------------------------------------------*/
/* forward declarations */

// DRT
namespace DRT
{
  class Discretization;
}  // namespace DRT


// Epetra
class Epetra_Comm;
class Epetra_MultiVector;

// IO
namespace IO
{
  class DiscretizationWriter;
  class OutputControl;
}  // namespace IO

// Teuchos
namespace Teuchos
{
  class ParameterList;
}  // namespace Teuchos

/*----------------------------------------------------------------------------*/
/* definition of methods */
namespace NLNSOL
{
  namespace UTILS
  {
    /*! \class DebugWriterBase
     *  \brief Base class of output functionality for nonlinear solvers
     *
     *  Base class for single-field and multi-field debug writers, that allows for
     *  writing debug output to a separate control file.
     *
     *  The names for output data are generated automatically by GenerateName().
     *  Since the output names need to be unique, an output name is created using
     *  the debug output name and a unique ID delivered by #counter_.
     *
     *  The debug writer creates an additional output file #log_, that encodes
     *  the auto-generated name and a user-provided description of the output
     *  quantity.
     *
     *  \author mayr.mt \date 05/2015
     */
    class DebugWriterBase : public Teuchos::VerboseObject<DebugWriterBase>
    {
     public:
      //! @name Construction / Destruction
      //@{

      //! destructor
      virtual ~DebugWriterBase(){};

      //@}

      //! @name Setup
      //@{

      //! Initialize class members
      virtual void Init(const Epetra_Comm& comm,        ///< communicator
          Teuchos::RCP<const DRT::Discretization> dis,  ///< discretization of single field
          const std::string& name,                      ///< name for output file
          const unsigned int step                       ///< time step number
          ) = 0;

      //! Setup the object
      virtual void Setup() = 0;

      //@}

      //! @name Access to setup status indicators
      //@{

      //! Has Init() already been called?
      virtual bool IsInit() const { return isinit_; }

      //! Has Setup() already been called?
      virtual bool IsSetup() const { return issetup_; }

      //@}

      //! Output of data
      //@{

      //! Write vector to output
      virtual void WriteVector(Teuchos::RCP<const Epetra_MultiVector>
                                   vec,  ///< vector with result data (to be outputted)
          const std::string&
              description,  ///< user-provided description of output quantity for log-file #log_
          const IO::VectorType vt = IO::dofvector  ///< vector type
          ) const = 0;

      //@}

      //! @name Utilities
      //@{

      //! Returns a character string describing the operator
      virtual const char* Label() const { return "NLNSOL::UTILS::DebugWriterBase"; }

      //! Returns a character string describing the operator briefly
      virtual const char* LabelShort() const { return "DebugWriterBase"; }

      //@}

     protected:
      //! @name Construction / Destruction
      //@{

      //! constructor
      DebugWriterBase();

      //@}

      //! @name Setup
      //@{

      //! Set communicator
      virtual void SetComm(const Epetra_Comm& comm);

      //! Set output filename
      virtual void SetOutputFileName(const std::string& name);

      //! Set time step number
      virtual void SetStep(const unsigned int step);

      /*! \brief Create the log file and write header into file
       *
       *  \note Only processor 0 is writing to the log file.
       */
      virtual void CreateLogFile();

      //@}

      //! Set setup flags
      //@{

      //! Init() has been called.
      virtual void SetIsInit()
      {
        isinit_ = true;
        return;
      }

      //! Setup() has been called.
      virtual void SetIsSetup()
      {
        issetup_ = true;
        return;
      }

      //@}


      //! @name Access to class members
      //@{

      //! Return name of debug output appendix
      const std::string& GetName() const { return name_; }

      //! Return time step number
      unsigned int GetStep() const { return step_; }

      //! communicator
      virtual Teuchos::RCP<const Epetra_Comm> Comm() const;

      //@}
      //! @name Utilities
      //@{

      /*! \brief Generate unique output name
       *
       *  Use either default prefix "nlndbg_" or a user-provided prefix.
       *
       *  Append counter #counter_ to name to make name unique. Increase counter
       *  #counter_ afterwards to ensure that numbering is always unique.
       */
      virtual const std::string GenerateName() const;

      /*! \brief Write entry to log file
       *
       *  Write the auto-generated output name \p name and the user-provided
       *  description \p description to the log file.
       *
       *  \note Only processor 0 is writing to the log file.
       */
      virtual void WriteLogFileEntry(
          const std::string& name,        ///< auto-generated name of output data
          const std::string& description  ///< user-provided description
          ) const;

      //@}

     private:
      //! Has Init() been called?
      bool isinit_;

      //! Has Setup() been called?
      bool issetup_;

      //! name of debug output file
      std::string name_;

      /*! \brief Counter for automated name generation
       *
       *  \note Needs to be mutable to enable increasing the counter every time
       *  the debug writer is used for writing output.
       */
      mutable unsigned int counter_;

      /*! \brief Output stream for log-file
       *
       *  Log file to list encodings of auto-generated names and user-provided
       *  descriptions.
       */
      Teuchos::RCP<std::ofstream> log_;

      //! Communicator
      Teuchos::RCP<const Epetra_Comm> comm_;

      //! time step number
      unsigned int step_;

    };  // class DebugWriterBase

    // ---------------------------------------------------------------------------
    // ---------------------------------------------------------------------------

    /*! \class DebugWriterSingleField
     *  \brief Single-field output functionality for nonlinear solvers
     *
     *  \author mayr.mt \date 05/2015
     */
    class DebugWriterSingleField : public DebugWriterBase
    {
     public:
      //! @name Construction / Destruction
      //@{

      //! constructor
      DebugWriterSingleField();

      //! destructor
      virtual ~DebugWriterSingleField(){};

      //@}

      //! @name Setup
      //@{

      //! Initialize class members
      virtual void Init(const Epetra_Comm& comm,        ///< communicator
          Teuchos::RCP<const DRT::Discretization> dis,  ///< discretization of single field
          const std::string& name,                      ///< name for output file
          const unsigned int step                       ///< time step number
      );

      //! Setup the output tool
      virtual void Setup();

      //@}

      //! @name Output of data
      //@{

      //! Write vector to output
      virtual void WriteVector(Teuchos::RCP<const Epetra_MultiVector>
                                   vec,  ///< vector with result data (to be outputted)
          const std::string&
              description,  ///< user-provided description of output quantity for log-file #log_
          const IO::VectorType vt = IO::dofvector  ///< vector type
          ) const;

      //@}

      //! @name Utilities
      //@{

      //! Returns a character string describing the operator
      virtual const char* Label() const { return "NLNSOL::UTILS::DebugWriterSingleField"; }

      //! Returns a character string describing the operator briefly
      virtual const char* LabelShort() const { return "DebugWriterSingleField"; }

      //@}

     protected:
     private:
      //! discretization
      Teuchos::RCP<const DRT::Discretization> dis_;

      //! current control file
      Teuchos::RCP<IO::OutputControl> control_;

      //! writer to control file
      Teuchos::RCP<IO::DiscretizationWriter> writer_;

    };  // class DebugWriterSingleField
  }     // namespace UTILS
}  // namespace NLNSOL

#endif /* NLN_UTILS_DEBUGWRITER_H */
