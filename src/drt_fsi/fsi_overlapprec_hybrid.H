/*----------------------------------------------------------------------------*/
/*!
\file fsi_overlapprec_hybrid.H

\brief Hybrid Additive/Multiplicative Schwarz Block Preconditioner for FSI

\level 1

\maintainer Matthias Mayr
*/
/*----------------------------------------------------------------------------*/

#ifndef FSI_OVERLAPPREC_HYBRID_H
#define FSI_OVERLAPPREC_HYBRID_H

/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/* forward declarations */
namespace LINALG
{
  namespace SOLVER
  {
    class IFPACKPreconditioner;
  }  // namespace SOLVER
}  // namespace LINALG


/*----------------------------------------------------------------------------*/
/* headers */
// Ifpack
#include <Ifpack.h>

// Teuchos
#include <Teuchos_RCP.hpp>

// baci
#include "fsi_overlapprec.H"

namespace FSI
{
  /*! \class OverlappingBlockMatrixHybridSchwarz
   *  \brief Block Matrix including a Hybrid Additive/Multiplicative Schwarz Block Preconditioner
   *
   *  Any type of existing AMG preconditioner (formally a multiplicative Schwarz
   *  preconditioner) is hybridized with a interface-local additive Schwarz
   *  preconditioner in order to reduce the error close to the interface more
   *  efficiently.
   *
   *  \note A first implementation has been done by Maximilian Noll during his
   *  Semesterarbeit "Domain Decomposition/Redistribution and Hybrid Additive/
   *  Multiplicative Schwarz Preconditioning for Monolithic Fluid-Structure
   *  Interaction Solvers" (2015).
   *
   *  \author mayr.mt \date 11/2015
   */
  class OverlappingBlockMatrixHybridSchwarz : public OverlappingBlockMatrix
  {
   public:
    //! Constructor
    OverlappingBlockMatrixHybridSchwarz(const LINALG::MultiMapExtractor& maps,
        ADAPTER::FSIStructureWrapper& structure,  ///< structure field
        ADAPTER::Fluid& fluid,                    ///< fluid field
        ADAPTER::AleFsiWrapper& ale,              ///< ale field
        bool structuresplit,                      ///< structure split or fluid split?
        int symmetric, std::vector<std::string>& blocksmoother, std::vector<double>& schuromega,
        std::vector<double>& omega, std::vector<int>& iterations, std::vector<double>& somega,
        std::vector<int>& siterations, std::vector<double>& fomega, std::vector<int>& fiterations,
        std::vector<double>& aomega, std::vector<int>& aiterations,
        int analyze,                             ///< Run FSIAMG anaylzer?
        INPAR::FSI::LinearBlockSolver strategy,  ///< type of linear preconditioner
        std::list<int> interfaceprocs,           ///< IDs of processors that own interface nodes
        INPAR::FSI::Verbosity verbosity,         ///< verbosity of FSI algorithm
        FILE* err = NULL                         ///< error file handle
    );

    //! @name Setup
    //@{

    //! setup of block preconditioners
    virtual void SetupPreconditioner();

    //@}

    //! @name Applying the preconditioner
    //@{

    /*! \brief Apply preconditioner
     *
     *  The hybrid preconditioner \f$M^{-1}_{H}\f$ is a chain of additive and
     *  multiplicative Schwarz preconditioners \f$M^{-1}_{AS}\f$ and
     *  \f$M^{-1}_{MS}\f$, respectively:
     *  \f[
     *    M^{-1}_{H} = M^{-1}_{AS,pre} \circ M^{-1}_{MS} \circ M^{-1}_{AS,post}
     *  \f]
     *
     *  It is applied with 3 stationary Richardson iterations (cf. eq. (4.12) in
     *  term paper by Maximilian Noll). Therefore, we embed ApplyInverse()-calls
     *  to the inner preconditioners #ifpackprec_ and #amgprec_ within a dampped
     *  stationary Richardson iteration, where each ApplyInverse() might also be
     *  a stationary Richardson iteration by itself.
     */
    virtual int ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const;

    //@}

    //! @name Attribute access functions
    //@{

    /// Returns a character string describing the operator.
    virtual const char* Label() const;

    //@}

   protected:
   private:
    //! type of preconditioner
    INPAR::FSI::LinearBlockSolver strategy_;

    //! IFPACK preconditioner (additive Schwarz)
    Teuchos::RCP<LINALG::SOLVER::IFPACKPreconditioner> ifpackprec_;

    //! IFPACK preconditioner (additive Schwarz)
    Teuchos::RCP<Ifpack_Preconditioner> directifpackprec_;

    //! AMG preconditioner (multiplicative Schwarz)
    Teuchos::RCP<FSI::OverlappingBlockMatrix> amgprec_;

    //! IDs of processors that own interface nodes
    std::list<int> interfaceprocs_;

    //! Apply #ifpackprec_ on all procs or on interface procs only
    bool additiveschwarzeverywhere_;
  };
}  // namespace FSI

#endif  // #ifndef FSI_OVERLAPPREC_HYBRID_H
