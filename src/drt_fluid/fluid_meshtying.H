/*!----------------------------------------------------------------------
\file fluidrmeshtying.H

\brief Methods to apply meshtying to fluid and scatra systems

<pre>
Maintainer: Andreas Ehrl
            ehrl@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15252
</pre>

 *----------------------------------------------------------------------*/
#ifndef FLUID_MESHTYING_H_
#define FLUID_MESHTYING_H_

#include "../drt_inpar/inpar_fluid.H"

#include <Teuchos_RCP.hpp>
#include <Epetra_Vector.h>


//forward declarations
namespace DRT
{
  class Discretization;

  namespace ELEMENTS
  {
    class FluidEleParameter;
  }
}

namespace ADAPTER
{
  class CouplingMortar;
}

namespace LINALG
{
  class Solver;
  class SparseMatrix;
  class MultiMapExtractor;
  class MapExtractor;
  class BlockSparseMatrixBase;
  class SparseOperator;
}


/*!
\brief  This class manages local transformations(rotation) of velocity fields

        o used for rotationally symmetric boundary conditions
        o used for rotationally symmetric boundary conditions

\author ehrl

 */

namespace FLD
{
namespace UTILS
{
  class MapExtractor;
  class InterfaceSplitStrategy;
}
class Meshtying
{
  friend class FluidEleParameter;

  public:

    //! Constructor
    Meshtying(Teuchos::RCP<DRT::Discretization>      dis,                      ///> actual discretisation
              LINALG::Solver&               solver,
              Teuchos::ParameterList&       params,                   ///> meshting parameter list
              const UTILS::MapExtractor*    surfacesplitter = NULL);  ///> surface splitter


    virtual ~Meshtying(){return;};

    //! Set up mesh-tying framework
    Teuchos::RCP<LINALG::SparseOperator> Setup();

    //! Prepare matrix and residual for meshtying
    void PrepareMeshtyingSystem(Teuchos::RCP<LINALG::SparseOperator>&  sysmat,     ///> sysmat established by the element routine
                      Teuchos::RCP<Epetra_Vector>&                     residual);  ///> residual established by the element routine

    //! Solve mesh-tying problem
    void SolveMeshtying(LINALG::Solver&                 solver,
                        Teuchos::RCP<LINALG::SparseOperator>     sysmat,
                        Teuchos::RCP<Epetra_Vector>&              incvel,
                        Teuchos::RCP<Epetra_Vector>              residual,
                        int                             itnum,
                        Teuchos::RCP<Epetra_MultiVector>         w,
                        Teuchos::RCP<Epetra_MultiVector>         c,
                        bool                            project
                        );

    //! Adjust null-space for Krylov projector (slave node are in-active)
    void KrylovProjection(Teuchos::RCP<Epetra_Vector>   c);

    //! Time update for Lagrange multiplier
    void UpdateLag();

    //! Output: maps & projection matrix
    void OutputSetUp();

    //! Output: split sparse matrix
    void OutputSparseMatrixSplit(Teuchos::RCP<LINALG::SparseOperator>  conmat);

    //! Output: single blocks of the block matrix
    void OutputBlockMatrix(Teuchos::RCP<LINALG::SparseOperator>      blockmatrix,
                           Teuchos::RCP<Epetra_Vector>               residual);

    //! Output: split vector
    void OutputVectorSplit(Teuchos::RCP<Epetra_Vector>  vector);

    //! Analyze system matrix
    void AnalyzeMatrix(
        Teuchos::RCP<LINALG::SparseMatrix>        sparsematrix); ///> sparse matrix to analyze

    //! Replace matrix entries
    /// Replace computed identity matrix by a real identity matrix
    void ReplaceMatrixEntries(
        Teuchos::RCP<LINALG::SparseMatrix>        sparsematrix); ///> sparse matrix to analyze

  private:

    //! Adjust the residual for a saddle point system (PrepareMeshtyingSystem)
    void ResidualSaddlePointSystem(Teuchos::RCP<Epetra_Vector>  residual);   ///> residual established by the element routine

    //! Prepare condensation for a block matrix (PrepareMeshtyingSystem)
    void CondensationBlockMatrix(Teuchos::RCP<LINALG::SparseOperator>&   sysmat,     ///> sysmat established by the element routine
                                 Teuchos::RCP<Epetra_Vector>&            residual);  ///> residual established by the element routine

    //! Prepare condensation for sparse matrix (PrepareMeshtyingSystem)
    void CondensationSparseMatrix(
        Teuchos::RCP<LINALG::SparseOperator>&   sysmat,      ///> sysmat established by the element routine
        Teuchos::RCP<Epetra_Vector>&            residual);   ///> residual established by the element routine

    //! Split sparse matrix and save blocks in a std::vector<Teuchos::RCP<Epetra_Vector> >
    void SplitMatrix(Teuchos::RCP<LINALG::SparseOperator>                matrix,         ///> sysmat established by the element routine
                     std::vector<Teuchos::RCP<LINALG::SparseMatrix> >&   splitmatrix);   ///> container for the split blocks

    //! Split vector and save parts in a std::vector<Teuchos::RCP<Epetra_Vector> >
    void SplitVector(Teuchos::RCP<Epetra_Vector>                   vector,           ///> vector to split
                     std::vector<Teuchos::RCP<Epetra_Vector> >&    splitvector);     ///> container for the split vector

    //! Split vector and save parts in a std::vector<Teuchos::RCP<Epetra_Vector> >
    void SplitVectorBasedOn3x3(Teuchos::RCP<Epetra_Vector> orgvector,        ///> original vector based on 3x3 blockmatrix
                               Teuchos::RCP<Epetra_Vector> vectorbasedon2x2);     ///> split vector based on 2x2 blockmatrix

    //! Condensation operation for a sparse matrix:
    /// the sysmat is manipulated via a second sparse matrix
    /// Assembling is slower, since the graph cannot be saved
    void CondensationOperationSparseMatrix(
        Teuchos::RCP<LINALG::SparseOperator>&              sysmat,         ///> sysmat established by the element routine
        Teuchos::RCP<Epetra_Vector>&                       residual,       ///> residual established by the element routine
        std::vector<Teuchos::RCP<LINALG::SparseMatrix> > & splitmatrix,    ///> container with split original sysmat
        std::vector<Teuchos::RCP<Epetra_Vector> >&         splitvector);   ///> container with split original residual

    //! Condensation operation for a block matrix
    /// the original blocks (nn, nm, mn, mm) are manipulated directly;
    /// the remaining blocks (ns, ms, ss, sn, sm) are not touched at all,
    /// since finally a 2x2 block matrix is solved
    void CondensationOperationBlockMatrix(
        Teuchos::RCP<LINALG::SparseOperator>&              sysmat,         ///> sysmat established by the element routine
        Teuchos::RCP<Epetra_Vector>&                       residual,       ///> residual established by the element routine
        std::vector<Teuchos::RCP<Epetra_Vector> >&         splitvector);   ///> container with split original residual

    //! Compute and update the increments of the slave node (Condensation)
    void UpdateSlaveDOF(Teuchos::RCP<Epetra_Vector>&   inc);

    //! Prepare the saddle point system in sparse matrix
    /// the system can only be solved by Umfpack
    void PrepareSaddlePointSystem(
        Teuchos::RCP<LINALG::SparseOperator>    sysmat,            ///> sysmat established by the element routine
        Teuchos::RCP<LINALG::SparseMatrix>      mergedsysmat,      ///> merged sysmat (including Lagrange multiplier)
        Teuchos::RCP<Epetra_Vector>             residual,          ///> residual established by the element routine
        Teuchos::RCP<Epetra_Vector>             mergedresidual);   ///> merged residual (including Lagrange multiplier)

    //! Prepare the saddle point system in block matrix
    /// the system can only be solved by SIMPLER algorithm
    void PrepareSaddlePointSystemPC(
        Teuchos::RCP<LINALG::SparseOperator>        sysmat,            ///> sysmat established by the element routine
        Teuchos::RCP<LINALG::BlockSparseMatrixBase> blocksysmat,       ///> merged sysmat (including Lagrange multiplier)
        Teuchos::RCP<Epetra_Vector>                 residual,          ///> residual established by the element routine
        Teuchos::RCP<Epetra_Vector>                 mergedresidual);   ///> merged residual (including Lagrange multiplier)

    //! Update increment and Lagrange multiplier after each Newton iteration
    void UpdateSaddlePointSystem(
        Teuchos::RCP<Epetra_Vector>      inc,            ///> original increment
        Teuchos::RCP<Epetra_Vector>      mergedinc);     ///> merged increment

  private:

    //! discretisation
    Teuchos::RCP<DRT::Discretization> discret_;

    LINALG::Solver& solver_;                        // standard solver object

    //! meshting options
    /// 0: no_meshtying     -> no mesh-tying
    /// 0: condensed_smat   -> condensation in a sparse matrix
    /// 1: condensed_bmat   -> condensation in a block matrix
    /// 2: sps_coupled      -> saddle point system in a sparse matrix
    /// 3: sps_pc           -> saddle point system in a block matrix
    int msht_;

    //! the processor ID from the communicator
    int                           myrank_;

    // interface splitter
    const UTILS::MapExtractor* surfacesplitter_;

    //! dof row map of the complete system
    const Epetra_Map* dofrowmap_;

    //! dof row map of the complete system
    Teuchos::RCP<Epetra_Map> problemrowmap_;

    //! dof rowmap of all nodes, which are not on the interface
    Teuchos::RCP<Epetra_Map> gndofrowmap_;

    //! slave & master dof rowmap
    Teuchos::RCP<Epetra_Map> gsmdofrowmap_;

    //! slave dof rowmap
    Teuchos::RCP<Epetra_Map> gsdofrowmap_;

    //! master dof rowmap
    Teuchos::RCP<Epetra_Map> gmdofrowmap_;

    //! dof row map of the complete system and the lagrange multiplier
    Teuchos::RCP<Epetra_Map> glmdofrowmap_;

    //! merged map for saddle point system and 2x2 block matrix
    Teuchos::RCP<Epetra_Map> mergedmap_;

    //! lagrange multiplier at time n+1
    Teuchos::RCP<Epetra_Vector> lag_;

    //! lagrange multiplier at time n
    Teuchos::RCP<Epetra_Vector> lagold_;

    //! time integration factor (only necessary for saddle point system)
    double theta_;

    //! adapter to mortar framework
    Teuchos::RCP<ADAPTER::CouplingMortar> adaptermeshtying_;

    //! 2x2 (3x3) block matrix for solving condensed system (3x3 block matrix)
    Teuchos::RCP<LINALG::SparseOperator> sysmatsolve_;

    //! flag defining pressure coupling
    bool pcoupled_;

};  // end  class Meshtying
} // end namespace FLD

#endif   /* FLUID_MESHTYING_H_ */
