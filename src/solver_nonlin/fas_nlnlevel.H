/*----------------------------------------------------------------------------*/
/*!

\brief Nonlinear level containing the representation and operations for a single
nonlinear multigrid level

\level 3

\maintainer Matthias Mayr
*/

/*----------------------------------------------------------------------------*/

#ifndef FAS_NLNLEVEL_H
#define FAS_NLNLEVEL_H

/*----------------------------------------------------------------------------*/
/* headers */

// Teuchos
#include <Teuchos_RCP.hpp>
#include <Teuchos_VerboseObject.hpp>

/*----------------------------------------------------------------------------*/
/* forward declarations */

class Epetra_Comm;
class Epetra_CrsMatrix;
class Epetra_MultiVector;
class Epetra_Vector;

namespace NLNSOL
{
  class NlnOperatorBase;
  class NlnProblemBase;

  namespace UTILS
  {
    class NlnConfig;
  }  // namespace UTILS
}  // namespace NLNSOL

/*----------------------------------------------------------------------------*/
/* definition of classes */
namespace NLNSOL
{
  namespace FAS
  {
    /*! \class NlnLevel
     *  \brief Nonlinear level
     *
     *  Here, we store all data that belongs to a multigrid level like transfer
     *  operators #rop_ and #pop_, the coarse level matrix representation #A_,
     *  level smoothers (#presmoother_, #postsmoother_, #coarsesolver_), ...
     *
     *  We have access to the nonlinear problem (#nlnproblem_) either on the
     *  fine level (if #levelid_ == 0) or on the coarse level (if #levelid_ > 0).
     *
     *  Each level can exchange / transfer data to the next finer level, only.
     *  This transfer is managed and initiated by the NLNSOL::FAS::AMGHierarchy.
     *
     *  \sa NLNSOL::FAS::AMGHierarchy
     *
     * \author mayr.mt \date 06/2014
     */
    class NlnLevel : public Teuchos::VerboseObject<NlnLevel>
    {
     public:
      //! @name Construction / Destruction
      //@{

      //! Constructor (empty)
      NlnLevel();

      //! Destructor
      virtual ~NlnLevel() {}

      //@}

      //! @name Setup
      //@{

      //! Initialization
      virtual void Init(const int levelid,         ///< ID of current level
          const int numlevels,                     ///< total number of levels
          Teuchos::RCP<const Epetra_CrsMatrix> A,  ///< matrix
          Teuchos::RCP<const Epetra_CrsMatrix> R,  ///< restriction operator
          Teuchos::RCP<const Epetra_CrsMatrix> P,  ///< prolongation operator
          const Epetra_Comm& comm,                 ///< communicator
          Teuchos::RCP<const NLNSOL::UTILS::NlnConfig>
              config,                  ///< configuration of nonlinear solver
          const std::string listname,  ///< name of sublist with configuration of this operator
          Teuchos::RCP<NLNSOL::NlnProblemBase>
              nlnproblem,  ///< reference to nonlinear problem evaluator
          Teuchos::RCP<const Epetra_MultiVector> nullspace =
              Teuchos::null  ///< set of null space vectors
      );

      //! Setup
      virtual void Setup();

      /*! \brief Update matrix and rebuild preconditioners
       *
       *  Update the coarse level matrix by the recently computes RAP. Rebuild
       *  possible preconditioners of #presmoother_, #postsmoother_, and
       *  #coarsesolver_.
       *
       *  \author mayr.mt \date 03/2015
       */
      virtual void UpdateMatrix(Teuchos::RCP<const Epetra_CrsMatrix> A  ///< matrix
      );

      //@}

      //! @name VCycle-related stuff
      //@{

      /*! \brief Restriction vector to next coarser level
       *
       *  Restriction is performed on coarser level, i.e. from 'level-1' to
       *  'level', since only the coarse level knows the transfer operator.
       *
       *  \return coarse vector
       *
       *  \sa ProlongateToNextFinerLevel()
       *
       *  \author mayr.mt \date 10/2013
       */
      virtual Teuchos::RCP<Epetra_MultiVector> RestrictToNextCoarserLevel(
          Teuchos::RCP<const Epetra_MultiVector> vf  ///< fine vector to be restricted
          ) const;

      /*! \brief Prolongate vector to next finer level
       *
       *  Prolongation is performed on coarser level, i.e. from 'level' to
       *  'level-1', since only the coarse level knows the transfer operator.
       *
       *  \return fine vector
       *
       *  \sa RestrictToNextCoarserLevel()
       *
       *  \author mayr.mt \date 11/2013
       */
      virtual Teuchos::RCP<Epetra_MultiVector> ProlongateToNextFinerLevel(
          Teuchos::RCP<const Epetra_MultiVector> vc  ///< coarse vector to be prolongated
          ) const;

      //@}

      /*! \brief Check convergence
       *
       *  \return bool: true if converged, false otherwise
       */
      virtual bool ConvergenceCheck(const Epetra_MultiVector& f  ///< residual
          ) const;

      /*! \brief Check convergence
       *
       *  Provide L2-norm of residual as well.
       *
       *  \return bool: true if converged, false otherwise
       */
      virtual bool ConvergenceCheck(const Epetra_MultiVector& f,  ///< residual
          double& fnorm2  /// L2-norm of residual (to be filled)
          ) const;

      //! @name Utilities
      //@{

      //! Print method
      virtual void Print(std::ostream& os  ///< output stream used for printing
          ) const;

      //! Returns a character string describing the operator
      virtual const char* Label() const { return "NLNSOL::FAS::NlnLevel"; }

      //@}

      //! @name Access methods
      //@{

      //! Current level
      virtual int LevelID() const { return levelid_; }

      //! Max number of levels
      virtual int NumLevels() const { return numlevels_; }

      //! Is this level the coarsest level?
      virtual bool IsCoarsestLevel() const;

      //! Get matrix
      virtual Teuchos::RCP<const Epetra_CrsMatrix> GetMatrix() const;

      //! Get restriction operator
      virtual Teuchos::RCP<const Epetra_CrsMatrix> GetROp() const;

      //! Get prolongation operator
      virtual Teuchos::RCP<const Epetra_CrsMatrix> GetPOp() const;

      //! Get null space
      virtual Teuchos::RCP<const Epetra_MultiVector> GetNullSpace() const;

      //! Get presmoother
      virtual Teuchos::RCP<NLNSOL::NlnOperatorBase> GetPreSmoother() const;

      //! Get postsmoother
      virtual Teuchos::RCP<NLNSOL::NlnOperatorBase> GetPostSmoother() const;

      //! DofRowMap of this level
      const Epetra_BlockMap& DofRowMap() const;

      //@}

      //! @name Check status of member variables
      //@{

      //! Matrix set?
      virtual bool HaveMatrix() const { return (not A_.is_null()); }

      //! Restriction operator set?
      virtual bool HaveROp() const { return (not rop_.is_null()); }

      //! Prolongation operator set?
      virtual bool HavePOp() const { return (not pop_.is_null()); }

      //! Null space set?
      virtual bool HaveNullSpace() const { return (not nullspace_.is_null()); }

      //@}

      //! @name Smoothing methods
      //@{

      /*! \brief Do presmoothing sweeps on this level
       *
       *  \return Integer error code, set to 0 if successful.
       *
       *  \author mayr.mt \date 10/2013
       */
      virtual int DoPreSmoothing(const Epetra_MultiVector& f,  ///< residual
          Epetra_MultiVector& x                                ///< solution
          ) const;

      /*! \brief Do postsmoothing sweeps on this level
       *
       *  \return Integer error code, set to 0 if successful.
       *
       *  \author mayr.mt \date 10/2013
       */
      virtual int DoPostSmoothing(const Epetra_MultiVector& f,  ///< residual
          Epetra_MultiVector& x                                 ///< solution
          ) const;

      /*! \brief Do coarse level solve
       *
       *  \return Integer error code, set to 0 if successful.
       *
       *  \author mayr.mt \date 10/2013
       */
      virtual int DoCoarseLevelSolve(const Epetra_MultiVector& f,  ///< current residual
          Epetra_MultiVector& x                                    ///< solution
          ) const;

      //@}

      //! Access to problem where the smoother acts on
      virtual Teuchos::RCP<NLNSOL::NlnProblemBase> NlnProblem() const;

      /*! \brief Check for stagnation of level smoother
       *
       *  We check for stagnation in all pre- and poststmoothers and in the coarse
       *  level solver.
       *
       *  One might don't care for stagtnation in the presmoother since this might
       *  be cured by the coarser levels.
       *
       *  \return bool: true if stagnation is detected, false otherwise
       *
       *  \author mayr.mt \date 03/2015
       */
      virtual bool CheckStagnation() const;

     protected:
     private:
      //! Returns pointer to Epetra_Comm associated with this operator.
      virtual const Epetra_Comm& Comm() const;

      //! Access to global configuration
      virtual Teuchos::RCP<const NLNSOL::UTILS::NlnConfig> Configuration() const;

      //! Get name of this operator's parameter list in the global configuration
      virtual const std::string& MyListName() const { return listname_; }

      //! Has Init() already been called?
      virtual bool IsInit() const { return isinit_; }

      //! Has Setup() already been called?
      virtual bool IsSetup() const { return issetup_; }

      //! @name Initialization
      //@{

      //! Set level ID
      virtual void SetLevelID(const int levelid  ///< level ID (to be set)
      )
      {
        levelid_ = levelid;
        return;
      }

      //! Set total number of levels
      virtual void SetNumLevels(const int numlevels  ///< total number of levels (to be set)
      )
      {
        numlevels_ = numlevels;
        return;
      }

      //! Set matrix
      virtual void SetMatrix(Teuchos::RCP<const Epetra_CrsMatrix> A  ///< matrix for this level
      );

      //! Set restriction operator
      virtual void SetROp(Teuchos::RCP<const Epetra_CrsMatrix>
              R  ///< restriction operator from next finer level to this level
      );

      //! Set prolongation operator
      virtual void SetPOp(Teuchos::RCP<const Epetra_CrsMatrix>
              P  ///< prolongation operator from this level to next finer level
      );

      //! Set null space
      virtual void SetNullSpace(Teuchos::RCP<const Epetra_MultiVector>
              nullspace  ///< set of null space vectors on this level
      );

      //! Init() has been called?
      bool isinit_;

      //! Setup() has been called?
      bool issetup_;

      //@}

      //! @name AMG operators and hierarchy
      //@{

      //! current level
      int levelid_;

      //! total number of levels
      int numlevels_;

      //! matrix
      Teuchos::RCP<Epetra_CrsMatrix> A_;

      //! restriction operator
      Teuchos::RCP<Epetra_CrsMatrix> rop_;

      //! prolongation operator
      Teuchos::RCP<Epetra_CrsMatrix> pop_;

      //@}

      //! @name FAS quantities
      //@{

      //! restriction of fine level solution: xbar_ = rop_ * x_fine
      Teuchos::RCP<Epetra_Vector> xbar_;

      //! restriction of fine level residual: fbar_ = rop_ * f_fine
      Teuchos::RCP<Epetra_Vector> fbar_;

      //! variationally evaluated coarse grid residual: fhat_ = rop_ * F(P*R*x_fine)
      Teuchos::RCP<Epetra_Vector> fhat_;

      //!

      //@}

      //! Access to nonlinear problem
      Teuchos::RCP<NLNSOL::NlnProblemBase> nlnproblem_;

      //! communicator
      Teuchos::RCP<const Epetra_Comm> comm_;

      //! @name Configuration
      //@{

      //! Configuration of this nonlinear operator
      Teuchos::RCP<const NLNSOL::UTILS::NlnConfig> config_;

      //! Name of sublist in #config_ for this operator
      std::string listname_;

      //@}

      //! @name Smoothers for this level
      //@{

      //! Presmoother
      Teuchos::RCP<NLNSOL::NlnOperatorBase> presmoother_;

      //! Postsmoother
      Teuchos::RCP<NLNSOL::NlnOperatorBase> postsmoother_;

      //@}

      //! @name Coarse level solver
      //@{

      //! Coarse level solver
      Teuchos::RCP<NLNSOL::NlnOperatorBase> coarsesolver_;

      //@}

      //! @name Debug tools
      //@{

      /*! \brief Set of null space vectors
       *
       *  The null space is needed for setup of the hierarchy only. However, for
       *  debugging purposes we store it in the level. This can go away at some
       *  point in future.
       *
       *  \author mayr.mt \date 10/2015
       */
      Teuchos::RCP<Epetra_MultiVector> nullspace_;

      //@}
    };  // class NlnLevel
  }     // namespace FAS
}  // namespace NLNSOL


#endif /* FAS_NLNLEVEL_H */
