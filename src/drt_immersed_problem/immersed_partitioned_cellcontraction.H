/*---------------------------------------------------------------------*/
/*! \file
\brief partitioned immersed contraction algorithm

\level 3

\maintainer Jonas Eichinger
*/
/*---------------------------------------------------------------------*/
#ifndef IMMERSED_PARTITIONED_CELLCONTRACTION_H
#define IMMERSED_PARTITIONED_CELLCONTRACTION_H

#include "immersed_partitioned.H"
#include "immersed_field_exchange_manager.H"

/*----------------------------------------------------------------------*
 | forward declarations                                                  |
 *----------------------------------------------------------------------*/
namespace POROELAST
{
  // class Monolithic;
  class PoroScatraBase;
}  // namespace POROELAST

namespace ADAPTER
{
  class FPSIStructureWrapper;
  class FSIStructureWrapperImmersed;
}  // namespace ADAPTER

namespace SSI
{
  class SSI_Part2WC_BIOCHEMOMECHANO;
}

/*----------------------------------------------------------------------*/
namespace IMMERSED
{
  class ImmersedPartitionedCellContraction : public ImmersedPartitioned
  {
   public:
    /// setup partitioned immersed cell migration algorithm
    explicit ImmersedPartitionedCellContraction(
        const Teuchos::ParameterList& params, const Epetra_Comm& comm);

    /*! \brief Initialize this object

    \warning none
    \return int
    \date 08/16
    \author rauch  */
    virtual int Init(const Teuchos::ParameterList& params);

    /*! \brief Setup all class internal objects and members

    \warning none
    \return void
    \date 08/16
    \author rauch  */
    virtual void Setup();

    /// read restart
    virtual void ReadRestart(int step);

   protected:
    /// composed Field-Field-Interaction operator
    void CouplingOp(const Epetra_Vector& x, Epetra_Vector& F, const FillType fillFlag);

    /// background operator
    void BackgroundOp(
        Teuchos::RCP<Epetra_Vector> backgrd_dirichlet_values, const FillType fillFlag);

    /// immersed operator
    Teuchos::RCP<Epetra_Vector> ImmersedOp(
        Teuchos::RCP<Epetra_Vector> bdry_traction, const FillType fillFlag);

    /// preparation of background field solve
    void PrepareBackgroundOp();

    /// preparation of immersed field solve
    void PrepareImmersedOp();

    /// prepare output of data to HD
    virtual void PrepareOutput();

    /// output of data to HD
    virtual void Output();

    /// update fields
    virtual void Update();

    /// prepare time step of single fields and print header
    virtual void PrepareTimeStep();

    /// initial guess
    virtual Teuchos::RCP<Epetra_Vector> InitialGuess();

   private:
    //! @name Data Management
    //@{
    DRT::ImmersedFieldExchangeManager*
        exchange_manager_;  //!< manager for data exchange between fields
    //@}

    //! @name Search Trees
    //@{
    Teuchos::RCP<GEO::SearchTree> fluid_SearchTree_;  //!< search tree for fluid domain
    Teuchos::RCP<GEO::SearchTree> cell_SearchTree_;   //!< search tree for cell domain
    //@}

    //! @name current nodal positions of the discretizations
    //@{
    std::map<int, LINALG::Matrix<3, 1>>*
        currpositions_cell_;  //!< pointer to map of vectors for search tree containing current
                              //!< structural positions
    std::map<int, LINALG::Matrix<3, 1>>*
        currpositions_ECM_;  //!< pointer to map of vectors for search tree containing current fluid
                             //!< positions
    //@}

    /// cell structure of the immersed problem
    Teuchos::RCP<::ADAPTER::FSIStructureWrapperImmersed> cellstructure_;

    /// poro structure of the immersed problem
    Teuchos::RCP<::ADAPTER::FPSIStructureWrapper> porostructure_;

    /// pointer to poroelast subproblem (poro-scatra interaction)
    Teuchos::RCP<POROELAST::PoroScatraBase> poroscatra_subproblem_;

    /// pointer to cellscatra subproblem (structure-scatra interaction)
    Teuchos::RCP<SSI::SSI_Part2WC_BIOCHEMOMECHANO> cellscatra_subproblem_;

    //! @name basic information for parallelism
    //@{
    int myrank_;   //!< processor id
    int numproc_;  //!< total number of processors
    //@}

    /// pointer to global problem
    DRT::Problem* globalproblem_;

    /// coupling scheme
    bool displacementcoupling_;

    //! @name Pointers to Discretizations
    //@{
    Teuchos::RCP<DRT::Discretization> backgroundstructuredis_;  //!< background ecm discretisation
    Teuchos::RCP<DRT::Discretization>
        backgroundfluiddis_;  //!< background interstitial flow discretisation
    Teuchos::RCP<DRT::Discretization> immerseddis_;  //!< immersed cell discretisation
    Teuchos::RCP<DRT::Discretization> scatradis_;    //!< cell scatra discretisation
    //@}

  };  // class ImmersedPartitionedProtrusionFormation
}  // namespace IMMERSED

#endif /* IMMERSED_PARTITIONED_CELLCONTRACTION_H */
