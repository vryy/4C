/*----------------------------------------------------------------------*/
/*!
\brief Wear interface implementation.

\level 2

\maintainer Matthias Mayr
*/
/*----------------------------------------------------------------------*/
#ifndef CONTACT_WEAR_INTERFACE_H_
#define CONTACT_WEAR_INTERFACE_H_

#include "../drt_inpar/inpar_contact.H"
#include "../drt_inpar/inpar_wear.H"
#include "contact_interface.H"
#include <Epetra_FEVector.h>

namespace WEAR
{
  class WearInterface : public CONTACT::CoInterface
  {
   public:
    /*!
    \brief Constructor

    */
    WearInterface(const Teuchos::RCP<MORTAR::IDataContainer>& idata_ptr, const int id,
        const Epetra_Comm& comm, const int dim, const Teuchos::ParameterList& icontact,
        bool selfcontact, INPAR::MORTAR::RedundantStorage redundant);

    /*!
    \brief Destructor

    */
    virtual ~WearInterface() {}

    /*!
    \brief Assemble second mortar D matrix for both-sided wear

    */
    virtual void AssembleD2(LINALG::SparseMatrix& dglobal);

    /*!
    \brief Assemble Mortar wear matrices T and E

    */
    virtual void AssembleTE(LINALG::SparseMatrix& tglobal, LINALG::SparseMatrix& eglobal);

    /*!
    \brief Assemble Mortar wear matrices T and E (maser side)

    */
    virtual void AssembleTE_Master(LINALG::SparseMatrix& tglobal, LINALG::SparseMatrix& eglobal);

    /*!
    \brief Assemble matrices LinT containing linearizations
           w.r.t. displacements
    */
    virtual void AssembleLinT_D(LINALG::SparseMatrix& lintglobal);

    /*!
    \brief Assemble matrices LinT containing linearizations
           w.r.t. displacements (for master side)
    */
    virtual void AssembleLinT_D_Master(LINALG::SparseMatrix& lintglobal);

    /*!
    \brief Assemble matrices LinT containing linearizations
           w.r.t. LM
    */
    virtual void AssembleLinT_LM(LINALG::SparseMatrix& lintglobal);

    /*!
    \brief Assemble matrices LinT containing linearizations
           w.r.t. LM
    */
    virtual void AssembleLinT_LM_Master(LINALG::SparseMatrix& lintglobal);

    /*!
    \brief Assemble matrices LinE containing linearizations
           w.r.t. displacements
    */
    virtual void AssembleLinE_D(LINALG::SparseMatrix& lineglobal);

    /*!
    \brief Assemble matrices LinE containing linearizations
           w.r.t. displacements (for master side)
    */
    virtual void AssembleLinE_D_Master(LINALG::SparseMatrix& lineglobal);

    /*!
    \brief Assemble matrix S containing linearizations

    This method builds an algebraic form of the FULL linearization
    of the normal contact condition g~ = 0. Concretely, this
    includes assembling the linearizations of the slave side
    nodal normals and of the Mortar matrices D  and M.

    */
    virtual void AssembleS(LINALG::SparseMatrix& sglobal);

    /*!
    \brief Assemble matrix S containing linearizations w

    */
    virtual void AssembleLinG_W(LINALG::SparseMatrix& sglobal);

    /*!
    \brief Assemble matrix LinStick containing linearizations

    This method builds an algebraic form of the FULL linearization
    of the tangential stick condition delta tg = 0. Concretely, this
    includes assembling the linearizations of the slave side
    nodal tangents and of the Mortar matrices D  and M.

    */
    virtual void AssembleLinStick(LINALG::SparseMatrix& linstickLMglobal,
        LINALG::SparseMatrix& linstickDISglobal, Epetra_Vector& linstickRHSglobal);
    /*!
    \brief Assemble matrix LinSlip containing linearizations

    This method builds an algebraic form of the FULL linearization
    of the tangential slip condition. Concretely, this
    includes assembling the linearizations of the slave side
    nodal tangents and of the Mortar matrices D  and M.

    */
    virtual void AssembleLinSlip(LINALG::SparseMatrix& linslipLMglobal,
        LINALG::SparseMatrix& linslipDISglobal, Epetra_Vector& linslipRHSglobal);

    /*!
    \brief Assemble matrix LinSlip containing w linearizations

    */
    virtual void AssembleLinSlip_W(LINALG::SparseMatrix& linslipWglobal);

    /*!
    \brief Assemble matrices W containing linearizations

    This method builds an algebraic form of the FULL linearization
    of the normal contact and slip contact condition for ~w.
    --> w.r.t. lagr. mult.

    */
    virtual void AssembleLinWLm(LINALG::SparseMatrix& sglobal);
    virtual void AssembleLinWLmSl(LINALG::SparseMatrix& sglobal);
    virtual void AssembleLinWLmSt(LINALG::SparseMatrix& sglobal);

    /*!
    \brief Assemble wear w

     This method assembles the weighted wear vector.
     */
    virtual void AssembleWear(Epetra_Vector& wglobal);

    /*!
    \brief Build active set (nodes / dofs) of this interface

    If the flag init==true, the active set is initialized (for t=0)
    according to the contact initialization defined in the input file.

    */
    virtual bool BuildActiveSet(bool init = false);

    /*!
    \brief Build corresponding active set for master side

    */
    virtual bool BuildActiveSetMaster();

    /*!
    \brief Check mortar wear T derivatives with finite differences

    */
    void FDCheckMortarTDeriv();

    /*!
    \brief Check mortar wear T derivatives with finite differences (Master)

    */
    void FDCheckMortarT_Master_Deriv();

    /*!
    \brief Check mortar wear E derivatives with finite differences

    */
    void FDCheckMortarEDeriv();

    /*!
    \brief Check mortar wear E derivatives with finite differences (for master)

    */
    void FDCheckMortarE_Master_Deriv();

    /*!
    \brief Check mortar wear T derivatives with finite differences
      --> for wear condition

    */
    void FDCheckDerivT_D(LINALG::SparseMatrix& lintdis);

    /*!
    \brief Check mortar wear T derivatives with finite differences
      --> for wear condition (Master)

    */
    void FDCheckDerivT_D_Master(LINALG::SparseMatrix& lintdis);

    /*!
    \brief Check mortar wear E derivatives with finite differences
      --> for wear condition

    */
    void FDCheckDerivE_D(LINALG::SparseMatrix& linedis);

    /*!
    \brief Check mortar wear E derivatives with finite differences
      --> for wear condition (Master)

    */
    void FDCheckDerivE_D_Master(LINALG::SparseMatrix& linedis);
    /*!
    \brief Check weighted gap g derivatives with finite differences

    */
    void FDCheckGapDeriv();

    /*!
    \brief Check weighted gap g derivatives with finite differences

    */
    void FDCheckGapDeriv_W();

    /*!
    \brief Check weighted wear ~w derivatives with finite differences
           derivation w.r.t. displ.

    */
    void FDCheckWearDeriv();

    /*!
    \brief Check weighted wear ~w derivatives with finite differences
           derivation w.r.t. lagr.-mult.

    */
    void FDCheckWearDerivLm();

    /*!
    \brief Check slip condition derivatives with finite differences

    */
    virtual void FDCheckSlipDeriv(LINALG::SparseMatrix& linslipLMglobal,
        LINALG::SparseMatrix& linslipDISglobal, LINALG::SparseMatrix& linslipWglobal);

    /*!
    \brief Assemble inactive rhs (incremental delta_w_)
    */
    virtual void AssembleInactiveWearRhs(Epetra_Vector& inactiverhs);

    /*!
    \brief Assemble inactive rhs (incremental delta_w_)
    */
    virtual void AssembleInactiveWearRhs_Master(Epetra_FEVector& inactiverhs);

    /*!
    \brief Assemble wear-cond. rhs
    */
    virtual void AssembleWearCondRhs(Epetra_Vector& rhs);

    /*!
    \brief Assemble wear-cond. rhs
    */
    virtual void AssembleWearCondRhs_Master(Epetra_FEVector& rhs);

    /*!
    \brief Initialize / reset interface for contact

    Derived version!

    */
    virtual void Initialize();


    /*!
    \brief Returning dofs for both-sided wear mapping

    */
    virtual const Teuchos::RCP<Epetra_Map> InvolvedDofs() const { return involveddofs_; }

    virtual const Teuchos::RCP<Epetra_Map> InvolvedNodes() const { return involvednodes_; }

    /*!
    \brief Set element areas

    Derived version!

    */
    void SplitSlaveDofs();
    void SplitMasterDofs();
    /*!
    \brief Set element areas

    Derived version!

    */
    void SetElementAreas();

    /*!
    \brief Evaluate nodal normals

    */
    virtual void EvaluateNodalNormals() const;


    /*!
    \brief Evaluate nodal normals

    */
    virtual void ExportNodalNormals() const;

    /*!
    \brief Update interface Wear variable sets

    This update is usually only done ONCE in the initialization phase
    and sets up the wear unknowns (only dofs) for the whole
    simulation.

    */
    virtual void UpdateWSets(int offset_if, int maxdofwear, bool bothdiscr);

    /*!
    \brief Get map of slave wear dofs (Filled()==true is prerequisite)

    */
    virtual const Teuchos::RCP<Epetra_Map> WDofs() const
    {
      if (Filled())
        return wdofmap_;
      else
        dserror("CONTACT::WearInterface::FillComplete was not called");
      exit(EXIT_FAILURE);  // calm down the compiler
    }

    /*!
    \brief Get map of master wear dofs (Filled()==true is prerequisite)

    */
    virtual const Teuchos::RCP<Epetra_Map> WMDofs() const
    {
      if (Filled())
        return wmdofmap_;
      else
        dserror("CONTACT::WearInterface::FillComplete was not called");
      exit(EXIT_FAILURE);  // calm down the compiler
    }

    /*!
    \brief Get map of Lagrange multiplier dofs (Filled()==true is prerequisite)

    */
    virtual const Teuchos::RCP<Epetra_Map> SNDofs() const
    {
      if (Filled())
        return sndofmap_;
      else
        dserror("CONTACT::WearInterface::FillComplete was not called");
      exit(EXIT_FAILURE);  // calm down the compiler
    }

    /*!
    \brief Get map of Lagrange multiplier dofs (Filled()==true is prerequisite)

    */
    virtual const Teuchos::RCP<Epetra_Map> MNDofs() const
    {
      if (Filled())
        return mndofmap_;
      else
        dserror("CONTACT::WearInterface::FillComplete was not called");
      exit(EXIT_FAILURE);  // calm down the compiler
    }

    /*!
    \brief Get row map of active nodes (Filled()==true is prerequisite)

    */
    virtual const Teuchos::RCP<Epetra_Map> ActiveMasterNodes() const
    {
      if (Filled())
        return activmasternodes_;
      else
        dserror("CONTACT::Interface::FillComplete was not called");
      exit(EXIT_FAILURE);  // calm down the compiler
    }

    /*!
    \brief Get row map of active nodes (Filled()==true is prerequisite)

    */
    virtual const Teuchos::RCP<Epetra_Map> SlipMasterNodes() const
    {
      if (Filled())
        return slipmasternodes_;
      else
        dserror("CONTACT::Interface::FillComplete was not called");
      exit(EXIT_FAILURE);  // calm down the compiler
    }

    /*!
    \brief Get row map of active nodes (Filled()==true is prerequisite)

    */
    virtual const Teuchos::RCP<Epetra_Map> SlipMasterNDofs() const
    {
      if (Filled())
        return slipmn_;
      else
        dserror("CONTACT::Interface::FillComplete was not called");
      exit(EXIT_FAILURE);  // calm down the compiler
    }


    /*!
    \brief Get type of wear shapefnct

    */
    INPAR::WEAR::WearShape WearShapeFcn()
    {
      return DRT::INPUT::IntegralValue<INPAR::WEAR::WearShape>(imortar_, "WEAR_SHAPEFCN");
    }

   private:
    /*!
    \brief initialize node and element data container

    Derived verision!

    */
    virtual void InitializeDataContainer();


    // both-sided wear specific stuff
    Teuchos::RCP<Epetra_Map> involvednodes_;  // row map of all involved master nodes
    Teuchos::RCP<Epetra_Map> involveddofs_;   // row map of all involved master dofs

    Teuchos::RCP<Epetra_Map> wdofmap_;   // row map of all slave wear dofs
    Teuchos::RCP<Epetra_Map> wmdofmap_;  // row map of all master wear dofs

    Teuchos::RCP<Epetra_Map> sndofmap_;  // row map of all slave dofs (first entries)
    Teuchos::RCP<Epetra_Map> mndofmap_;  // row map of all master dofs (first entries)

    Teuchos::RCP<Epetra_Map>
        activmasternodes_;  // row map of all active master nodes (first entries)
    Teuchos::RCP<Epetra_Map>
        slipmasternodes_;              // row map of all active master nodes (first entries)
    Teuchos::RCP<Epetra_Map> slipmn_;  // row map of all active master nodes (first entries)

    bool wear_;      // bool for wear
    bool wearimpl_;  // bool for implicit wear
    bool wearpv_;    // bool for wear with own discretization
    bool wearboth_;  // bool for wear on both sides
    bool sswear_;    // bool steady state wear

  };  // class


}  // namespace WEAR

#endif /* CONTACT_WEAR_INTERFACE_H_ */
