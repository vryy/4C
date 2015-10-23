/*----------------------------------------------------------------------------*/
/*!
\file nln_problem.H

\brief Base class for interface of nonlinear solver to BACI

<pre>
Maintainer: Matthias Mayr
            mayr@mhpc.mw.tum.de
            089 - 289-10362
</pre>
*/

/*----------------------------------------------------------------------------*/

#ifndef NLN_PROBLEM_H
#define NLN_PROBLEM_H

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

namespace LINALG
{
  class SparseOperator;
} // namespace LINALG

namespace NLNSOL
{
  namespace UTILS
  {
    class DebugWriterBase;
  }  // namespace UTILS
} // namespace NLNSOL

namespace NOX
{
  namespace Abstract
  {
    class Group;
  } // namespace Abstract
} // namespace NOX

namespace Teuchos
{
  class ParameterList;
} // namespace Teuchos

/*----------------------------------------------------------------------------*/
/* definition of classes */
namespace NLNSOL
{
  /*! \class NlnProblem
   *  \brief Base class for interface of nonlinear solver to BACI
   *
   *  Interface class between the nonlinear solver and BACI. We utilize the
   *  NOX::Abstract::Group() to access the time integration within BACI for
   *  Evaluate() / computeF().
   *
   *  Whenever a nonlinear operator wants to call Evaluate() or
   *  ConvergenceCheck(), the NLNSOL::NlnProblem provides access to the time
   *  integration scheme in BACI.
   *
   *  NLNSOL::NlnProblem acts on the DofRowMap() of the (fine) level, only. For
   *  coarse level representations, we use NLNSOL::NlnProblemCoarseLevel().
   *
   *  \sa NlnProblemCoarseLevel
   *
   *  \author mayr.mt \date 06/2014
   */
  class NlnProblem : public Teuchos::VerboseObject<NlnProblem>
  {
  public:

    //! @name Construction / Destruction
    //@{

    //! constructor (empty)
    NlnProblem();

    //! Destructor
    virtual ~NlnProblem(){};

    //@}

    //! @name Setup
    //@{

    //! Initialize member variables
    virtual void Init(const Epetra_Comm& comm, ///< communicator
        const Teuchos::ParameterList& params, ///< parameter list
        NOX::Abstract::Group& noxgrp, ///< NOX group for call back routines to time integration
        Teuchos::RCP<LINALG::SparseOperator> jac, ///< Jacobian matrix of underlying problem
        Teuchos::RCP<NLNSOL::UTILS::DebugWriterBase> dbgwriter = Teuchos::null ///< Debug writer for output of intermediate solutions
        );

    //! Setup of nonlinear problem
    virtual void Setup();

    /*! \brief Set coarse level residual modifications for FAS
     *
     *  \note Only needed for coarse levels. Implemented here only for purpose
     *  of inheritance.
     *
     *  \sa NLNSOL::NlnProblemCoarseLevel::SetupResidualModification
     */
    virtual void SetupResidualModification(
        Teuchos::RCP<const Epetra_MultiVector> xbar, ///< \f$\bar{x}\f$
        Teuchos::RCP<const Epetra_MultiVector> fbar ///< \f$\bar{F}\f$
        )
    {
      dserror("Not implemented in base class. This routine is needed only in "
          "case of a coarse level problem within FAS-AMG.");
    }

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
     *  \author mayr.mt \date 09/2014
     */
    virtual void ComputeF(const Epetra_MultiVector& x, ///< current solution
        Epetra_MultiVector& f ///< current residual (to be filled)
        ) const;

    /*! \brief Compute nonlinear residual \f$F(x)\f$
     *
     *  \author mayr.mt \date 10/2015
     */
    virtual Teuchos::RCP<const Epetra_MultiVector> ComputeF(
        const Epetra_MultiVector& x ///< current solution
        ) const;

    /*! \brief Compute the plain residual on a coarse level
     *
     *  Since we are already on the fine level, we just call ComputeF().
     *
     *  \note Only needed for coarse levels. Implemented here only for purpose
     *  of inheritance.
     *
     *  \sa NLNSOL::NlnProblemCoarseLevel::ComputePlainF
     *
     *  \author mayr.mt \date 10/2015
     */
    virtual Teuchos::RCP<const Epetra_MultiVector> ComputePlainF(
        const Epetra_MultiVector& xc ///< current solution on coarse level
        ) const;

    /*! \brief Evaluate Jacobian (linearization of residual)
     *
     *  Updating the solution variables in the marching time integrator has
     *  already been done during computeF().
     *
     *  \author mayr.mt \date 06/2014
     */
    virtual void ComputeJacobian() const;

    //! Check convergence
    virtual bool ConvergenceCheck(const Epetra_MultiVector& f ///< current residual
        ) const;

    //! Check convergence
    virtual bool ConvergenceCheck(const Epetra_MultiVector& f, ///< current residual
        double& fnorm2 ///< residual L2 norm (to be filled)
        ) const;

    //@}

    //! @name Get attributes
    //@{

    //! communicator
    virtual const Epetra_Comm& Comm() const;

    //! parameter list
    virtual const Teuchos::ParameterList& Params() const;

    //! NOX group
    virtual NOX::Abstract::Group& NOXGroup() const;

    /*! \brief Level ID of level where this nonlinear problem lives on
     *
     *  Since this is the fine level problem, the level ID is always 0.
     *
     *  \author mayr.mt \date 03/2015
     */
    virtual const int LevelID() const { return 0; }

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
    virtual void WriteVector(Teuchos::RCP<const Epetra_MultiVector> vec, ///< vector with result data (to be outputted)
        const std::string& description, ///< description of output quantity
        const IO::DiscretizationWriter::VectorType vt = IO::DiscretizationWriter::dofvector ///< vector type
        ) const;

    /*! \brief Wrapper to write plain Epetra_MultiVector to output
     *
     *  \author mayr.mt \date 06/2015
     */
    virtual void WriteVector(const Epetra_MultiVector& vec, ///< vector with result data (to be outputted)
        const std::string& description, ///< description of output quantity
        const IO::DiscretizationWriter::VectorType vt = IO::DiscretizationWriter::dofvector ///< vector type
        ) const;

    //@}

    //! @name Access to underlying NOX group elements
    //@{

    //! get Jacobian as Epetra_Operator
    virtual Teuchos::RCP<Epetra_Operator> GetJacobianOperator() const;

    //@}

    //! @name Utilities
    //@{

    //! Returns a character string describing the operator
    virtual const char* Label() const { return "NLNSOL::NlnProblem"; }

    //! Returns reference to DOF row map of underlying problem
    virtual const Epetra_Map& DofRowMap() const;

    //@}

  protected:

    //! @name Attribute set methods
    //@{

    //! Init() has been called
    virtual void SetIsInit() { isinit_ = true; return; }

    //! Setup() has been called
    virtual void SetIsSetup() { issetup_ = true; return; }

    //@}

    //! @name Access to status indicators
    //@{

    //! Has Init() already been called?
    virtual const bool IsInit() const { return isinit_; }

    //! Has Setup() already been called?
    virtual const bool IsSetup() const { return issetup_; }

    /*! \brief Is there a valid debug writer?
     *
     *  The debug writer #dbgwriter_ is considered as valid and as ready to use
     *  if it's RCP_pointer != Teuchos::null and if it's setup routine has been
     *  called.
     *
     *  \author mayr.mt \date 06/2015
     */
    virtual const bool HaveDebugWriter() const;

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

    //! parameter list
    Teuchos::RCP<const Teuchos::ParameterList> params_;

    //! NOX group for access to time integrator via call-back routines
    Teuchos::RCP<NOX::Abstract::Group> noxgrp_;

    /*! \brief Pointer to Jacobian matrix of underlying problem
     *
     *  This needs to be a LINALG::SparseOperator since this is the base class
     *  for all matrix types in BACI. We can always extract an Epetra_Operator
     *  when we want to call a Trilinos package/solver.
     *
     *  \warning Access to time integration WIHTOUT using the NOX group
     *  #noxgrp_. Seems to be inconsistent since all other accesses to time
     *  integration are done via the NOX group #noxgrp_. However, usage of NOX
     *  group might go away at some point.
     */
    Teuchos::RCP<LINALG::SparseOperator> jac_;

    //! Debug writer for output of intermediate solutions
    Teuchos::RCP<NLNSOL::UTILS::DebugWriterBase> dbgwriter_;

    //! @name Convergence check
    //@{

    //! residual L2-norm tolerance
    double tolresl2_;

    //! Do we scale norms by \f$\sqrt{n}\f$?
    bool lengthscaling_;

    //@}

  }; // class NlnProblem
} // namespace NLNSOL

#endif /* NLN_PROBLEM_H */
