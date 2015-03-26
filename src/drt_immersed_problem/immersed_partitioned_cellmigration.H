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
  void CouplingOp(const Epetra_Vector &x, Epetra_Vector &F, const FillType fillFlag);

  /// background operator
  void BackgroundOp(Teuchos::RCP<Epetra_Vector> backgrd_dirichlet_values,
                    const FillType fillFlag);

  /// immersed operator
  Teuchos::RCP<Epetra_Vector> ImmersedOp(Teuchos::RCP<Epetra_Vector> bdry_traction,
                                         const FillType fillFlag);
  /// initial guess
  virtual Teuchos::RCP<Epetra_Vector> InitialGuess();

  /// get immersed nodes and determine their dofs
  void BuildImmersedDirichMap(Teuchos::RCP<DRT::Discretization> dis,
                              Teuchos::RCP<Epetra_Map>& dirichmap,
                              const Teuchos::RCP<const Epetra_Map>& dirichmap_original,
                              int dofsetnum);

  /// add immersed dirichlet values from immersed dis to systemvector of background dis
  void DoImmersedDirichletCond(Teuchos::RCP<Epetra_Vector> statevector, Teuchos::RCP<Epetra_Vector> dirichvals, Teuchos::RCP<Epetra_Map> dbcmap);

  /// prepare the immersed solve
  void PrepareImmersedOp();

  /// prepare the background solve
  void PrepareBackgroundOp();

  /// dummy read restart
  virtual void ReadRestart(int step){dserror("NOT PROPERLY IMPLEMENTED YET!!"); return;};

  /// prepare output of data to HD
  virtual void PrepareOutput();

  /// output of data to HD
  virtual void Output();

  /// update fields
  virtual void Update();

  /// initialize the relaxation parameters and print output
  void SetupRelaxation();

  /// setup poro algorithm (coupling, combined dofmaps, etc.)
  void SetupPoroAlgorithm();

  //! @name Various global vectors
  //@{
  Teuchos::RCP<Epetra_Vector> porofluid_artificial_velocity_;//!< background velocity interpolated from immersed dis
  Teuchos::RCP<Epetra_Vector> cell_bdry_traction_;           //!< boundary traction rhs on immersed cell
  Teuchos::RCP<Epetra_Vector> cell_bdry_traction_old_;       //!< old boundary traction rhs on immersed cell
  //@}

  //! @name Pointers to Discretizations
    //@{
  Teuchos::RCP<DRT::Discretization> backgroundfluiddis_;
  Teuchos::RCP<DRT::Discretization> backgroundstructuredis_;
  Teuchos::RCP<DRT::Discretization> immerseddis_;
  //@}

  //! @name parameters and variables for relaxation
    //@{
  Teuchos::RCP<NOX::LineSearch::UserDefinedFactory> aitkenfactory_;
  bool relaxforceglobally_;
  bool relaxvelglobally_;
  double forcerelax_;
  double velrelax_;
  int coupalgo_;
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

    //! @name current subsets of discretizations
      //@{
    std::map<int,std::set<int> > curr_subset_of_backgrounddis_;    //!< backgrd elements to evaluate the dirichlet interpolation
    //@}

    //! @name current subsets of discretizations
      //@{
    Teuchos::RCP<const Epetra_Map> dbcmap_ordinary_;       //!< dirichlet bc map before consideration of immersed values
    Teuchos::RCP<Epetra_Map> dbcmap_immersed_;             //!< dirichlet bc map of immersed values
    //@}

};// class ImmersedFSI
} // namspace IMMERSED

#endif
