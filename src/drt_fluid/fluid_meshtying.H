/*-----------------------------------------------------------*/
/*! \file

\brief Methods to apply meshtying to fluid and scatra systems

\maintainer Martin Kronbichler

\level 2

*/
/*-----------------------------------------------------------*/

#ifndef FLUID_MESHTYING_H_
#define FLUID_MESHTYING_H_

#include "../drt_inpar/inpar_fluid.H"

#include "../linalg/linalg_mapextractor.H"

#include <Teuchos_RCP.hpp>
#include <Epetra_Vector.h>


// forward declarations
namespace DRT
{
  class Discretization;

  namespace ELEMENTS
  {
    class FluidEleParameter;
  }
}  // namespace DRT

namespace MORTAR
{
  class MortarInterface;
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
  class KrylovProjector;
}  // namespace LINALG

namespace FLD
{
  namespace UTILS
  {
    class MapExtractor;
    class InterfaceSplitStrategy;
  }  // namespace UTILS
  class Meshtying
  {
    friend class FluidEleParameter;

   public:
    //! Constructor
    Meshtying(Teuchos::RCP<DRT::Discretization> dis,         ///> actual discretisation
        LINALG::Solver& solver,                              ///> solver
        int msht,                                            ///> meshting parameter list
        int nsd,                                             ///> number space dimensions
        const UTILS::MapExtractor* surfacesplitter = NULL);  ///> surface splitter

    const Epetra_Map* GetMergedMap();
    virtual ~Meshtying() { return; };

    //! Set up mesh-tying framework
    Teuchos::RCP<LINALG::SparseOperator> Setup(std::vector<int> coupleddof);

    //! Applied Dirichlet values are adapted on the slave side of the internal interface
    //! in order to avoid an over-constraint problem setup
    void CheckOverlappingBC(Teuchos::RCP<Epetra_Map> map  ///> map of boundary condition
    );

    //! Old routine handling Dirichlet conditions on the master side of the internal interface
    /// During PrepareTimeStep() DC are projected from the master to the slave
    void ProjectMasterToSlaveForOverlappingBC(
        Teuchos::RCP<Epetra_Vector>& velnp,   ///> solution vector n+1
        Teuchos::RCP<const Epetra_Map> bmaps  ///> map of boundary condition
    );

    //! Check whether Dirichlet BC are defined on the master
    void DirichletOnMaster(Teuchos::RCP<const Epetra_Map> bmaps  ///> map of boundary condition
    );

    //! Preparation for including Dirichlet conditions in the condensation process
    void IncludeDirichletInCondensation(
        const Teuchos::RCP<Epetra_Vector>& velnp,  ///> solution vector n+1
        const Teuchos::RCP<Epetra_Vector>& veln    ///> solution vector n
    );

    //! evaluation of matrix P with potential mesh relocation in ALE case
    void EvaluateWithMeshRelocation(
        Teuchos::RCP<Epetra_Vector>& dispnp);  ///> current ALE displacement vector

    //! Prepare matrix, shapederivatives and residual for meshtying
    void PrepareMeshtying(Teuchos::RCP<LINALG::SparseOperator>&
                              sysmat,  ///> sysmat established by the element routine
        const Teuchos::RCP<Epetra_Vector>&
            residual,                              ///> residual established by the element routine
        const Teuchos::RCP<Epetra_Vector>& velnp,  ///> current ALE displacement vector
        Teuchos::RCP<LINALG::BlockSparseMatrixBase>&
            shapederivatives);  ///> shapederivatives established by the element routine

    //! Prepare matrix and residual for meshtying
    void PrepareMeshtyingSystem(const Teuchos::RCP<LINALG::SparseOperator>&
                                    sysmat,  ///> sysmat established by the element routine
        const Teuchos::RCP<Epetra_Vector>&
            residual,                               ///> residual established by the element routine
        const Teuchos::RCP<Epetra_Vector>& velnp);  ///> current ALE displacement vector

    //! The residual has another length in case of bmat_merged --> residual has to be calculated in
    //! split form
    void ApplyPTToResidual(Teuchos::RCP<LINALG::SparseOperator> sysmat,
        Teuchos::RCP<Epetra_Vector> residual, Teuchos::RCP<LINALG::KrylovProjector> projector);

    //! Solve mesh-tying problem (including ALE case)
    void SolveMeshtying(LINALG::Solver& solver, const Teuchos::RCP<LINALG::SparseOperator>& sysmat,
        const Teuchos::RCP<Epetra_Vector>& incvel, const Teuchos::RCP<Epetra_Vector>& residual,
        const Teuchos::RCP<Epetra_Vector>& velnp, const int& itnum,
        const Teuchos::RCP<LINALG::KrylovProjector>& projector);

    //! Adjust null-space for Krylov projector (slave node are in-active)
    Teuchos::RCP<Epetra_Vector> AdaptKrylovProjector(Teuchos::RCP<Epetra_Vector> vec);

    //! Output: maps & projection matrix
    void OutputSetUp();

    //! Output: split sparse matrix
    void OutputSparseMatrixSplit(Teuchos::RCP<LINALG::SparseOperator> conmat);

    //! Output: single blocks of the block matrix
    void OutputBlockMatrix(
        Teuchos::RCP<LINALG::SparseOperator> blockmatrix, Teuchos::RCP<Epetra_Vector> residual);

    //! Output: split vector
    void OutputVectorSplit(Teuchos::RCP<Epetra_Vector> vector);

    //! Analyze system matrix
    void AnalyzeMatrix(
        Teuchos::RCP<LINALG::SparseMatrix> sparsematrix);  ///> sparse matrix to analyze

    //! Replace matrix entries
    /// Replace computed identity matrix by a real identity matrix
    void ReplaceMatrixEntries(
        Teuchos::RCP<LINALG::SparseMatrix> sparsematrix);  ///> sparse matrix to analyze

    //! Compute and update the increments of the slave node (including ALE case)
    void UpdateSlaveDOF(
        const Teuchos::RCP<Epetra_Vector>& inc, const Teuchos::RCP<Epetra_Vector>& velnp);

    //! Set the flag for multifield problems
    void IsMultifield(Teuchos::RCP<std::set<int>> condelements,  ///< conditioned elements of fluid
        const LINALG::MultiMapExtractor& domainmaps,     ///< domain maps for split of fluid matrix
        const LINALG::MultiMapExtractor& rangemaps,      ///< range maps for split of fluid matrix
        Teuchos::RCP<std::set<int>> condelements_shape,  ///< conditioned elements
        const LINALG::MultiMapExtractor&
            domainmaps_shape,  ///< domain maps for split of shape deriv. matrix
        const LINALG::MultiMapExtractor&
            rangemaps_shape,  ///< domain maps for split of shape deriv. matrix
        bool splitmatrix,     ///< flag for split of matrices
        bool ismultifield     ///< flag for multifield problems
    );

    //! Use the split of the fluid mesh tying for the sysmat
    void MshtSplit(Teuchos::RCP<LINALG::SparseOperator>& sysmat,
        Teuchos::RCP<LINALG::BlockSparseMatrixBase>& shapederivatives);

    //! Use the split of the fluid mesh tying for the shape derivatives
    void MshtSplitShape(Teuchos::RCP<LINALG::BlockSparseMatrixBase>& shapederivatives);

    //! Use the split of the multifield problem for the sysmat
    void MultifieldSplit(Teuchos::RCP<LINALG::SparseOperator>& sysmat);

    //! Use the split of the multifield problem for the shape derivatives
    void MultifieldSplitShape(Teuchos::RCP<LINALG::BlockSparseMatrixBase>& shapederivatives);

    //! Prepare condensation of the shape derivatives
    void CondensationOperationBlockMatrixShape(
        Teuchos::RCP<LINALG::BlockSparseMatrixBase>& shapederivatives);

   private:
    //! Prepare condensation for sparse matrix (including ALE case)
    void CondensationSparseMatrix(const Teuchos::RCP<LINALG::SparseOperator>&
                                      sysmat,  ///> sysmat established by the element routine
        const Teuchos::RCP<Epetra_Vector>&
            residual,  ///> residual established by the element routine
        const Teuchos::RCP<Epetra_Vector>& velnp);

    //! Prepare condensation for a block matrix (including ALE case)
    void CondensationBlockMatrix(const Teuchos::RCP<LINALG::SparseOperator>&
                                     sysmat,  ///> sysmat established by the element routine
        const Teuchos::RCP<Epetra_Vector>&
            residual,                               ///> residual established by the element routine
        const Teuchos::RCP<Epetra_Vector>& velnp);  ///> current velocity vector

    //! split sparse global system matrix into 3x3 block sparse matrix associated with interior,
    //! master, and slave dofs
    void SplitMatrix(Teuchos::RCP<LINALG::SparseOperator>
                         matrix,  //!< original sparse global system matrix before split
        Teuchos::RCP<LINALG::BlockSparseMatrixBase>&
            splitmatrix  //!< resulting block sparse matrix after split
    );

    //! Split vector and save parts in a std::vector<Teuchos::RCP<Epetra_Vector> >
    void SplitVector(Teuchos::RCP<Epetra_Vector> vector,         ///> vector to split
        std::vector<Teuchos::RCP<Epetra_Vector>>& splitvector);  ///> container for the split vector

    //! Split vector and save parts in a std::vector<Teuchos::RCP<Epetra_Vector> >
    void SplitVectorBasedOn3x3(
        Teuchos::RCP<Epetra_Vector> orgvector,  ///> original vector based on 3x3 blockmatrix
        Teuchos::RCP<Epetra_Vector> vectorbasedon2x2);  ///> split vector based on 2x2 blockmatrix

    //! Condensation operation for a sparse matrix (including ALE case):
    /// the sysmat is manipulated via a second sparse matrix
    /// Assembling is slower, since the graph cannot be saved
    void CondensationOperationSparseMatrix(
        const Teuchos::RCP<LINALG::SparseOperator>&
            sysmat,  ///> sysmat established by the element routine
        const Teuchos::RCP<Epetra_Vector>&
            residual,  ///> residual established by the element routine
        const Teuchos::RCP<LINALG::BlockSparseMatrixBase>&
            splitmatrix,  ///> container with split original sysmat
        const std::vector<Teuchos::RCP<Epetra_Vector>>&
            splitres,  ///> container with split original residual
        const std::vector<Teuchos::RCP<Epetra_Vector>>&
            splitvel  ///> container with split velocity vector
    );

    //! Condensation operation for a block matrix (including ALE case):
    /// the original blocks (nn, nm, mn, mm) are manipulated directly;
    /// the remaining blocks (ns, ms, ss, sn, sm) are not touched at all,
    /// since finally a 2x2 block matrix is solved
    void CondensationOperationBlockMatrix(
        const Teuchos::RCP<LINALG::SparseOperator>&
            sysmat,  ///> sysmat established by the element routine
        const Teuchos::RCP<Epetra_Vector>&
            residual,  ///> residual established by the element routine
        const std::vector<Teuchos::RCP<Epetra_Vector>>&
            splitres,  ///> container with split original residual
        const std::vector<Teuchos::RCP<Epetra_Vector>>&
            splitvel);  ///> container with split velocity vector

   private:
    //! discretisation
    Teuchos::RCP<DRT::Discretization> discret_;

    LINALG::Solver& solver_;  // standard solver object

    //! meshting options
    /// 0: no_meshtying     -> no mesh-tying
    /// 0: condensed_smat   -> condensation in a sparse matrix
    /// 1: condensed_bmat   -> condensation in a block matrix
    /// 2: condensed_bmat_merged   -> condensation in a block matrix

    /// deactivated:
    /// 3: sps_coupled      -> saddle point system in a sparse matrix
    /// 4: sps_pc           -> saddle point system in a block matrix
    int msht_;

    //! the processor ID from the communicator
    int myrank_;

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
    Teuchos::RCP<const Epetra_Map> gsdofrowmap_;

    //! master dof rowmap
    Teuchos::RCP<const Epetra_Map> gmdofrowmap_;

    //! merged map for saddle point system and 2x2 block matrix
    Teuchos::RCP<Epetra_Map> mergedmap_;

    //! vector containing time-depending values of the dirichlet condition
    /// valuesdc_ = (velnp after applying DC) - (veln)
    Teuchos::RCP<Epetra_Vector> valuesdc_;

    //! adapter to mortar framework
    Teuchos::RCP<ADAPTER::CouplingMortar> adaptermeshtying_;

    //! 2x2 (3x3) block matrix for solving condensed system (3x3 block matrix)
    Teuchos::RCP<LINALG::SparseOperator> sysmatsolve_;
    Teuchos::RCP<Epetra_Vector> residual_;
    //! flag defining pressure coupling
    bool pcoupled_;

    //! flag defining if Dirichlet  or Dirichlet-like boundary conditions are defined on the master
    //! side of the internal interface
    bool dconmaster_;

    //! flag for identifying first Newton iteration in each time step
    bool firstnonliniter_;

    //! number of space dimensions
    int nsd_;

    //! conditioned elements of fluid in multifield simulation
    Teuchos::RCP<std::set<int>> multifield_condelements_;

    //! domain maps for split of fluid matrix in multifield simulation
    LINALG::MultiMapExtractor multifield_domainmaps_;

    //! range maps for split of fluid matrix in multifield simulation
    LINALG::MultiMapExtractor multifield_rangemaps_;

    //! conditioned elements in multifield simulation
    Teuchos::RCP<std::set<int>> multifield_condelements_shape_;

    //! domain maps for split of shape deriv. matrix in multifield simulation
    LINALG::MultiMapExtractor multifield_domainmaps_shape_;

    //! domain maps for split of shape deriv. matrix in multifield simulation
    LINALG::MultiMapExtractor multifield_rangemaps_shape_;

    //! flag for split of matrices in multifield simulation
    bool multifield_splitmatrix_;

    //! flag for multifield problems in multifield simulation
    bool is_multifield_;

  };  // end  class Meshtying
}  // end namespace FLD

#endif /* FLUID_MESHTYING_H_ */
