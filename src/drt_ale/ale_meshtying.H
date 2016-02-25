/*--------------------------------------------------------------------------*/
/*!
\file ale_meshtying.H

\brief Mesh tying for ale problems

<pre>
Maintainer: Andy Wirtz
            wirtz@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089-289-15270
</pre>
*/
/*--------------------------------------------------------------------------*/
#ifndef SRC_DRT_ALE_ALE_MESHTYING_H_
#define SRC_DRT_ALE_ALE_MESHTYING_H_


#include "../drt_inpar/inpar_ale.H"

#include "../linalg/linalg_mapextractor.H"

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
}

namespace ALE
{
namespace UTILS
{
  class MapExtractor;
  class InterfaceSplitStrategy;
}
class Meshtying
{

  public:

    //! Constructor
    Meshtying(Teuchos::RCP<DRT::Discretization>      dis,             ///> actual discretisation
              LINALG::Solver&               solver,                   ///> solver
              int                           msht,                   ///> meshting parameter list
              int                           nsd,                      ///> number space dimensions
              const UTILS::MapExtractor*    surfacesplitter = NULL);  ///> surface splitter

    virtual ~Meshtying(){return;};

    //! Set up mesh-tying framework
    virtual Teuchos::RCP<LINALG::SparseOperator> Setup(std::vector<int> coupleddof, Teuchos::RCP<Epetra_Vector>& dispnp);

    //! Use the split of the ale mesh tying for the sysmat
    Teuchos::RCP<LINALG::SparseOperator> MshtSplit();

    //! Check weather Dirichlet BC are defined on the master
    void DirichletOnMaster(
        Teuchos::RCP<const Epetra_Map>  bmaps   ///> map of boundary condition
        );

    //! Prepare matrix and residual for meshtying
    void PrepareMeshtyingSystem(Teuchos::RCP<LINALG::SparseOperator>&  sysmat,     ///> sysmat established by the element routine
                      Teuchos::RCP<Epetra_Vector>&                     residual, ///> residual established by the element routine
                      Teuchos::RCP<Epetra_Vector>&                     dispnp); ///> current ALE displacement vector

    //! Set the flag for multifield problems
    void IsMultifield(const LINALG::MultiMapExtractor& interface,         ///< interface maps for split of ale matrix
                      bool ismultifield                                   ///< flag for multifield problems
                      );

    //! Use the split of the ale mesh tying for the sysmat
    void MshtSplit(Teuchos::RCP<LINALG::SparseOperator>& sysmat);

    //! Use the split of the multifield problem for the sysmat
    void MultifieldSplit(Teuchos::RCP<LINALG::SparseOperator>& sysmat);

    //! Call the constructor and the setup of the mortar coupling adapter
    virtual void AdapterMortar(std::vector<int> coupleddof);

    //! Compare the size of the slave and master dof row map
    virtual void CompareNumDof();

    //! Get function for the slave and master dof row map
    virtual void DofRowMaps();

    //! Get function for the P matrix
    virtual Teuchos::RCP<LINALG::SparseMatrix> GetMortarTrafo();

    //! Condensation operation for a block matrix (including ALE case):
    /// the original blocks (nn, nm, mn, mm) are manipulated directly;
    /// the remaining blocks (ns, ms, ss, sn, sm) are not touched at all,
    /// since finally a 2x2 block matrix is solved
    virtual void CondensationOperationBlockMatrix(
        Teuchos::RCP<LINALG::SparseOperator>&  sysmat,         ///> sysmat established by the element routine
        Teuchos::RCP<Epetra_Vector>&           residual,       ///> residual established by the element routine
        Teuchos::RCP<Epetra_Vector>&           dispnp);        ///> current displacement vector

    //! Compute and update the increments of the slave node
    virtual void UpdateSlaveDOF(Teuchos::RCP<Epetra_Vector>&   inc,
                        Teuchos::RCP<Epetra_Vector>&   dispnp);

    //! Recover method for Lagrange multipliers (do nothing in mesh tying case)
    virtual void Recover(Teuchos::RCP<Epetra_Vector>& inc){};

    //! Solve ALE mesh tying problem
    virtual int SolveMeshtying(
      LINALG::Solver&                       solver,
      Teuchos::RCP<LINALG::SparseOperator>  sysmat,
      Teuchos::RCP<Epetra_Vector>&          disi,
      Teuchos::RCP<Epetra_Vector>           residual,
      Teuchos::RCP<Epetra_Vector>&          dispnp);

    //! Split vector and save parts in a std::vector<Teuchos::RCP<Epetra_Vector> >
    void SplitVector(Teuchos::RCP<Epetra_Vector>                   vector,           ///> vector to split
                     std::vector<Teuchos::RCP<Epetra_Vector> >&    splitvector);     ///> container for the split vector

  protected:

    //! discretisation
    Teuchos::RCP<DRT::Discretization> discret_;

    LINALG::Solver& solver_;                        // standard solver object

    //! dof row map of the complete system
    const Epetra_Map* dofrowmap_;

    //! slave dof rowmap
    Teuchos::RCP<const Epetra_Map> gsdofrowmap_;

    //! master dof rowmap
    Teuchos::RCP<const Epetra_Map> gmdofrowmap_;

    //! merged map for saddle point system and 2x2 block matrix
    Teuchos::RCP<Epetra_Map> mergedmap_;

  private:

    //! Split vector and save parts in a std::vector<Teuchos::RCP<Epetra_Vector> >
    void SplitVectorBasedOn3x3(Teuchos::RCP<Epetra_Vector> orgvector,        ///> original vector based on 3x3 blockmatrix
                               Teuchos::RCP<Epetra_Vector> vectorbasedon2x2);     ///> split vector based on 2x2 blockmatrix

  private:

    //! meshting options
    /// 0: no_meshtying     -> no mesh-tying
    /// 1: yes_meshtying   -> condensation in a block matrix
    int msht_;

    //! the processor ID from the communicator
    int                           myrank_;

    // interface splitter
    const UTILS::MapExtractor* surfacesplitter_;

    //! dof row map of the complete system
    Teuchos::RCP<Epetra_Map> problemrowmap_;

    //! dof rowmap of all nodes, which are not on the interface
    Teuchos::RCP<Epetra_Map> gndofrowmap_;

    //! slave & master dof rowmap
    Teuchos::RCP<Epetra_Map> gsmdofrowmap_;

    //! vector containing time-depending values of the dirichlet condition
    /// valuesdc_ = (dispnp after applying DC) - (dispn)
    Teuchos::RCP<Epetra_Vector> valuesdc_;

    //! adapter to mortar framework
    Teuchos::RCP<ADAPTER::CouplingMortar> adaptermeshtying_;

    //! 2x2 (3x3) block matrix for solving condensed system (3x3 block matrix)
    Teuchos::RCP<LINALG::SparseOperator> sysmatsolve_;
    Teuchos::RCP<Epetra_Vector> residual_;

    //! flag defining if Dirichlet  or Dirichlet-like boundary conditions are defined on the master
    //! side of the internal interface
    bool dconmaster_;

    //! flag for identifying first Newton iteration in each time step
    bool firstnonliniter_;

    //! number of space dimensions
    int nsd_;

    //! interface maps for split of ale matrix in multifield simulation
    LINALG::MultiMapExtractor multifield_interface_;

    //! flag for multifield problems in multifield simulation
    bool is_multifield_;

};  // end  class Meshtying
} // end namespace ALE


#endif /* SRC_DRT_ALE_ALE_MESHTYING_H_ */
