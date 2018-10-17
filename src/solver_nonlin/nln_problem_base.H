/*----------------------------------------------------------------------------*/
/*!
\file nln_problem_base.H

\brief Base class for interface of nonlinear solver to BACI

\maintainer Matthias Mayr

\level 3
*/

/*----------------------------------------------------------------------------*/

#ifndef NLN_PROBLEM_BASE_H
#define NLN_PROBLEM_BASE_H

/*----------------------------------------------------------------------------*/
/* headers */

// Teuchos
#include <Teuchos_RCP.hpp>
#include <Teuchos_VerboseObject.hpp>

// baci
#include "../drt_io/io.H"
#include "../drt_lib/drt_dserror.H"

/*----------------------------------------------------------------------------*/
/* forward declarations */

class Epetra_Comm;
class Epetra_Map;
class Epetra_MultiVector;
class Epetra_Operator;

namespace NLNSOL
{
  namespace UTILS
  {
    class NlnConfig;
    class DebugWriterBase;
  }  // namespace UTILS
}  // namespace NLNSOL

namespace NOX
{
  namespace Abstract
  {
    class Group;
  }  // namespace Abstract
}  // namespace NOX

namespace Teuchos
{
  class ParameterList;
}  // namespace Teuchos

/*----------------------------------------------------------------------------*/
/* definition of classes */
namespace NLNSOL
{
  /*! \class NlnProblemBase
   *  \brief Base class for interface of nonlinear solver to BACI
   *
   *  Interface class between the nonlinear solver and BACI. We utilize this set
   *  of classes to access the time integration within BACI for Evaluate() /
   *  computeF() / computeJacobian / ...
   *
   *  Whenever a nonlinear operator wants to call Evaluate() or
   *  ConvergenceCheck(), the NLNSOL::NlnProblemBase and its derived classes
   *  provide access to the time integration scheme in BACI.
   *
   *  NLNSOL::NlnProblem and derived classes act on the DofRowMap() of the
   *  fine level, only. For coarse level representations, we use
   *  NLNSOL::NlnProblemCoarseLevel().
   *
   *  \author mayr.mt \date 11/2015
   */
  class NlnProblemBase : public Teuchos::VerboseObject<NlnProblemBase>
  {
   public:
    //! @name Construction / Destruction
    //@{

    //! constructor (empty)
    NlnProblemBase();

    //! Destructor
    virtual ~NlnProblemBase(){};

    //@}

    //! @name Setup
    //@{

    //! Initialize member variables
    virtual void Init(const Epetra_Comm& comm,                ///< communicator
        Teuchos::RCP<const NLNSOL::UTILS::NlnConfig> config,  ///< configuration
        const std::string listname,                           ///< name of sublist for this operator
        const Teuchos::ParameterList&
            params,  ///< parameter list with Teuchos::RCP to BACI evaluator class
        Teuchos::RCP<const Epetra_Map> dofrowmap,  ///< dof row map
        Teuchos::RCP<NLNSOL::UTILS::DebugWriterBase> dbgwriter =
            Teuchos::null  ///< Debug writer for output of intermediate solutions
    );

    //! Setup of nonlinear problem
    virtual void Setup() = 0;

    /*! \brief Set actual model evaluator
     *
     *  Extract Teuchos::RCP to the actual BACI model evaluator form
     *  Teuchos::ParameterList and store it.
     *
     *  \return Bool: true if successful, false otherwise
     *
     *  \author mayr.mt \date 11/2015
     */
    virtual bool SetModelEvaluator() = 0;

    /*! \brief Set Jacobian operator
     *
     *  Extract Jacobian operator form Teuchos::ParameterList and store it.
     *
     *  \return Bool: true if successful, false otherwise
     *
     *  \author mayr.mt \date 11/2015
     */
    virtual bool SetJacobianOperator() = 0;

    /*! \brief Set coarse level residual modifications for FAS
     *
     *  \note Only needed for coarse levels. Implemented here only for purpose
     *  of inheritance.
     *
     *  \sa NLNSOL::NlnProblemCoarseLevel::SetupResidualModification
     */
    virtual void SetupResidualModification(
        Teuchos::RCP<const Epetra_MultiVector> xbar,  ///< \f$\bar{x}\f$
        Teuchos::RCP<const Epetra_MultiVector> fbar   ///< \f$\bar{F}\f$
        ) = 0;

    //@}

    //! @name Actions on the nonlinear problem at hand
    //@{

    /*! \brief Compute nonlinear residual \f$F(x)\f$
     *
     *  Set the most recent solution 'x' in the NOX group #noxgrp_. Afterwards,
     *  the marching time integrator is asked to evaluate the residual via the
     *  computeF()-call.
     *
     *  Updating the solution variables in the marching time integrator is done
     *  during computeF().
     *
     *  Return the most recent residual in 'f' as a descent direction.
     *
     *  \note Solvers implemented in the NLNSOL framework require a residual
     *  that points into descent direction.
     *
     *  \note One can define the 'mechanical' residual as
     *  \f$F(x) = f_{int} - f_{ext}\f$ or as \f$F(x) = f_{ext} - f_{int}\f$.
     *  This choice is arbitrary, since both formulations result in exactly the
     *  same solution \f$x\f$. In many Finite Element Codes as well as in baci,
     *  often the first version is utilized which produces an ascent direction.
     *  For a descent direction, one has to use the latter formula.
     *  To be consistent with mathematical literature, the residual needs
     *  to be a descent direction. Hence, we need to switch the sign of the
     *  'mechanical' baci-residual.
     *
     *  \author mayr.mt \date 10/2015
     */
    virtual void ComputeF(const Epetra_MultiVector& x,  ///< current solution
        Epetra_MultiVector& f                           ///< current residual (to be filled)
        ) const = 0;

    /*! \brief Compute nonlinear residual \f$F(x)\f$
     *
     *  \author mayr.mt \date 10/2015
     */
    virtual Teuchos::RCP<const Epetra_MultiVector> ComputeF(
        const Epetra_MultiVector& x  ///< current solution
        ) const = 0;

    /*! \brief Evaluate Jacobian (linearization of residual)
     *
     *  Updating the solution variables in the marching time integrator has
     *  already been done during computeF().
     *
     *  \author mayr.mt \date 06/2014
     */
    virtual void ComputeJacobian() const = 0;

    //! Check convergence
    virtual bool ConvergenceCheck(const Epetra_MultiVector& f  ///< current residual
        ) const;

    //! Check convergence
    virtual bool ConvergenceCheck(const Epetra_MultiVector& f,  ///< current residual
        double& fnorm2                                          ///< residual L2 norm (to be filled)
        ) const;

    //@}

    //! @name Get attributes
    //@{

    //! communicator
    virtual const Epetra_Comm& Comm() const;

    //! Access to global configuration (given in xml-file)
    virtual Teuchos::RCP<const NLNSOL::UTILS::NlnConfig> Configuration() const;

    //! Get name of this operator's parameter list in the global configuration
    virtual const std::string& MyListName() const { return listname_; }

    //! Access to runtime parameter list
    virtual Teuchos::RCP<const Teuchos::ParameterList> Params() const;

    /*! \brief Level ID of level where this nonlinear problem lives on
     *
     *  \author mayr.mt \date 11/2015
     */
    virtual int LevelID() const = 0;

    //@}

    //! Debug output of data
    //@{

    /*! \brief Access to debug writer
     *
     *  Return RCP-pointer to debug writer. We allow for returning
     *  Teuchos::null, since the debug writer will not be set for 'production'
     *  code and, thus, it will be mostly Teuchos::null. This case is taken care
     *  of during output.
     *
     *  \author mayr.mt \date 06/2015
     */
    Teuchos::RCP<NLNSOL::UTILS::DebugWriterBase> DebugWriter() const;

    /*! \brief Write Teuchos::RCP<Epetra_MultiVector> to output
     *
     *  \author mayr.mt \date 06/2015
     */
    virtual void WriteVector(
        Teuchos::RCP<const Epetra_MultiVector> vec,  ///< vector with result data (to be outputted)
        const std::string& description,              ///< description of output quantity
        const IO::VectorType vt = IO::dofvector      ///< vector type
        ) const;

    /*! \brief Wrapper to write plain Epetra_MultiVector to output
     *
     *  \author mayr.mt \date 06/2015
     */
    virtual void WriteVector(
        const Epetra_MultiVector& vec,           ///< vector with result data (to be outputted)
        const std::string& description,          ///< description of output quantity
        const IO::VectorType vt = IO::dofvector  ///< vector type
        ) const;

    //@}

    //! get Jacobian as Epetra_Operator
    virtual Teuchos::RCP<Epetra_Operator> GetJacobianOperator() const = 0;

    //! @name Utilities
    //@{

    //! Returns a character string describing the operator
    virtual const char* Label() const { return "NLNSOL::NlnProblemBase"; }

    //! Returns Teuchos::RCP to DOF row map of underlying problem
    virtual Teuchos::RCP<const Epetra_Map> DofRowMap() const;

    //@}

   protected:
    //! @name Attribute set methods
    //@{

    //! Init() has been called
    virtual void SetIsInit()
    {
      isinit_ = true;
      return;
    }

    //! Setup() has been called
    virtual void SetIsSetup()
    {
      issetup_ = true;
      return;
    }

    //@}

    //! @name Access to status indicators
    //@{

    //! Has Init() already been called?
    virtual bool IsInit() const { return isinit_; }

    //! Has Setup() already been called?
    virtual bool IsSetup() const { return issetup_; }

    /*! \brief Is there a valid debug writer?
     *
     *  The debug writer #dbgwriter_ is considered as valid and as ready to use
     *  if it's RCP_pointer != Teuchos::null and if it's setup routine has been
     *  called.
     *
     *  \author mayr.mt \date 06/2015
     */
    virtual bool HaveDebugWriter() const;

    //@}

   private:
    //! @name Status variables
    //@{

    //! Init() has been called
    bool isinit_;

    //! Setup() has been called
    bool issetup_;

    //@}

    //! communicator
    Teuchos::RCP<const Epetra_Comm> comm_;

    //! @name Configuration
    //@{

    //! Configuration of this nonlinear operator (given in xml-file)
    Teuchos::RCP<const NLNSOL::UTILS::NlnConfig> config_;

    //! Name of sublist in #config_ for this operator
    std::string listname_;

    //! Parameter list for parameters that are configured during runtime
    Teuchos::RCP<const Teuchos::ParameterList> params_;

    //@}

    //! Row map of nonlinear problem at hand
    Teuchos::RCP<const Epetra_Map> dofrowmap_;

    //! Debug writer for output of intermediate solutions
    Teuchos::RCP<NLNSOL::UTILS::DebugWriterBase> dbgwriter_;

    //! @name Convergence check
    //@{

    //! residual L2-norm tolerance
    double tolresl2_;

    //! Do we scale norms by \f$\sqrt{n}\f$?
    bool lengthscaling_;

    //@}

  };  // class NlnProblemBase
}  // namespace NLNSOL

#endif /* NLN_PROBLEM_BASE_H */
