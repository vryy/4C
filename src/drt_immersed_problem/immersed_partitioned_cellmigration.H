  /*!----------------------------------------------------------------------
\file immersed_partitioned_fsi_dirichletneumann.H

\brief partitioned immersed fsi algorithm for dirichlet-neumann coupling (volume dirichlet coupling)

<pre>
Maintainers: Andreas Rauch
             rauch@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289 -15240
</pre>
*----------------------------------------------------------------------*/
#ifndef IMMERSED_PARTITIONED_CELLMIGRATION_H
#define IMMERSED_PARTITIONED_CELLMIGRATION_H

#include "immersed_partitioned.H"

/*----------------------------------------------------------------------*
 | forward declarations                                                  |
 *----------------------------------------------------------------------*/
namespace POROELAST
{
  class Monolithic;
}

namespace ADAPTER
{
  class FSIStructureWrapper;
}

/*----------------------------------------------------------------------*/

namespace IMMERSED {

class ImmersedPartitionedCellMigration : public ImmersedPartitioned
{

public:

  /// setup partitioned ImmersedFSI algorithm
  explicit ImmersedPartitionedCellMigration(const Epetra_Comm& comm);

  /// initialize search tree for structure discretization and ghost structure redundantly
  void SetupImmersedDiscretization();

  /// initialize search tree for poro discretization
  void SetupBackgroundDiscretization();

protected:

  /// composed FSI operator
  void CouplingOp();

  /// background operator
  Teuchos::RCP<Epetra_Vector> BackgroundOp();

  /// immersed operator
  Teuchos::RCP<Epetra_Vector> ImmersedOp();
  /// initial guess
  virtual Teuchos::RCP<Epetra_Vector> InitialGuess();

  /// get immersed nodes and determine their dofs
  void BuildImmersedDirichMap();

  /// add immersed dirichlet values from immersed dis to systemvector of background dis
  void DoImmersedDirichletCond();

  /// determine elements cut by the boundary
  void PrepareImmersedOp();

  /// dummy read restart
  virtual void ReadRestart(int step){dserror("NOT PROPERLY IMPLEMENTED YET!!"); return;};

  virtual void PrepareOutput();

  //! @name Pointers to Discretizations
    //@{
  Teuchos::RCP<DRT::Discretization> backgroundfluiddis_;
  Teuchos::RCP<DRT::Discretization> backgroundstructuredis_;
  Teuchos::RCP<DRT::Discretization> immerseddis_;
  //@}

  /// cell structure of the immersed problem
  Teuchos::RCP< ::ADAPTER::FSIStructureWrapper>    cellstructure_;

  /// poroelast subproblem
  Teuchos::RCP<POROELAST::Monolithic>       poroelast_subproblem_;

  int myrank_;
  int numproc_;

  DRT::Problem* globalproblem_;

  bool displacementcoupling_;
  bool multicellmigration_;

  //! @name Search Trees
    //@{
    Teuchos::RCP<GEO::SearchTree> fluid_SearchTree_;  //!< search tree for fluid domain
    Teuchos::RCP<GEO::SearchTree> cell_SearchTree_;   //!< search tree for cell domain
    //@}

  //! @name current nodal positions of the discretizations
    //@{
    std::map<int,LINALG::Matrix<3,1> > currpositions_cell_; //!< map of vectors for search tree containing current structural positions
    std::map<int,LINALG::Matrix<3,1> > currpositions_ECM_;  //!< map of vectors for search tree containing current fluid positions
   //@}


};// class ImmersedFSI
} // namspace IMMERSED

#endif
