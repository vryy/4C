/*----------------------------------------------------------------------*/
/*! \file
\brief Wear interface implementation.

\level 2

*/
/*----------------------------------------------------------------------*/
#ifndef FOUR_C_CONTACT_WEAR_INTERFACE_HPP
#define FOUR_C_CONTACT_WEAR_INTERFACE_HPP

#include "4C_config.hpp"

#include "4C_contact_interface.hpp"
#include "4C_inpar_contact.hpp"
#include "4C_inpar_wear.hpp"

#include <Epetra_FEVector.h>

FOUR_C_NAMESPACE_OPEN

namespace Wear
{
  class WearInterface : public CONTACT::Interface
  {
   public:
    /*!
    \brief Constructor

    */
    WearInterface(const Teuchos::RCP<Mortar::InterfaceDataContainer>& interfaceData_ptr,
        const int id, const Epetra_Comm& comm, const int dim,
        const Teuchos::ParameterList& icontact, bool selfcontact);

    /*!
    \brief Assemble second mortar D matrix for both-sided wear

    */
    virtual void assemble_d2(Core::LinAlg::SparseMatrix& dglobal);

    /*!
    \brief Assemble Mortar wear matrices T and E

    */
    virtual void assemble_te(
        Core::LinAlg::SparseMatrix& tglobal, Core::LinAlg::SparseMatrix& eglobal);

    /*!
    \brief Assemble Mortar wear matrices T and E (maser side)

    */
    virtual void assemble_te_master(
        Core::LinAlg::SparseMatrix& tglobal, Core::LinAlg::SparseMatrix& eglobal);

    /*!
    \brief Assemble matrices LinT containing linearizations
           w.r.t. displacements
    */
    virtual void assemble_lin_t_d(Core::LinAlg::SparseMatrix& lintglobal);

    /*!
    \brief Assemble matrices LinT containing linearizations
           w.r.t. displacements (for master side)
    */
    virtual void assemble_lin_t_d_master(Core::LinAlg::SparseMatrix& lintglobal);

    /*!
    \brief Assemble matrices LinT containing linearizations
           w.r.t. LM
    */
    virtual void assemble_lin_t_lm(Core::LinAlg::SparseMatrix& lintglobal);

    /*!
    \brief Assemble matrices LinT containing linearizations
           w.r.t. LM
    */
    virtual void assemble_lin_t_lm_master(Core::LinAlg::SparseMatrix& lintglobal);

    /*!
    \brief Assemble matrices LinE containing linearizations
           w.r.t. displacements
    */
    virtual void assemble_lin_e_d(Core::LinAlg::SparseMatrix& lineglobal);

    /*!
    \brief Assemble matrices LinE containing linearizations
           w.r.t. displacements (for master side)
    */
    virtual void assemble_lin_e_d_master(Core::LinAlg::SparseMatrix& lineglobal);

    /*!
    \brief Assemble matrix S containing linearizations

    This method builds an algebraic form of the FULL linearization
    of the normal contact condition g~ = 0. Concretely, this
    includes assembling the linearizations of the slave side
    nodal normals and of the Mortar matrices D  and M.

    */
    void assemble_s(Core::LinAlg::SparseMatrix& sglobal) override;

    /*!
    \brief Assemble matrix S containing linearizations w

    */
    virtual void assemble_lin_g_w(Core::LinAlg::SparseMatrix& sglobal);

    /*!
    \brief Assemble matrix LinStick containing linearizations

    This method builds an algebraic form of the FULL linearization
    of the tangential stick condition delta tg = 0. Concretely, this
    includes assembling the linearizations of the slave side
    nodal tangents and of the Mortar matrices D  and M.

    */
    void assemble_lin_stick(Core::LinAlg::SparseMatrix& linstickLMglobal,
        Core::LinAlg::SparseMatrix& linstickDISglobal, Epetra_Vector& linstickRHSglobal) override;
    /*!
    \brief Assemble matrix LinSlip containing linearizations

    This method builds an algebraic form of the FULL linearization
    of the tangential slip condition. Concretely, this
    includes assembling the linearizations of the slave side
    nodal tangents and of the Mortar matrices D  and M.

    */
    void assemble_lin_slip(Core::LinAlg::SparseMatrix& linslipLMglobal,
        Core::LinAlg::SparseMatrix& linslipDISglobal, Epetra_Vector& linslipRHSglobal) override;

    /*!
    \brief Assemble matrix LinSlip containing w linearizations

    */
    virtual void assemble_lin_slip_w(Core::LinAlg::SparseMatrix& linslipWglobal);

    /*!
    \brief Assemble matrices W containing linearizations

    This method builds an algebraic form of the FULL linearization
    of the normal contact and slip contact condition for ~w.
    --> w.r.t. lagr. mult.

    */
    virtual void assemble_lin_w_lm(Core::LinAlg::SparseMatrix& sglobal);
    virtual void assemble_lin_w_lm_sl(Core::LinAlg::SparseMatrix& sglobal);
    virtual void assemble_lin_w_lm_st(Core::LinAlg::SparseMatrix& sglobal);

    /*!
    \brief Assemble wear w

     This method assembles the weighted wear vector.
     */
    virtual void assemble_wear(Epetra_Vector& wglobal);

    /*!
    \brief Build active set (nodes / dofs) of this interface

    If the flag init==true, the active set is initialized (for t=0)
    according to the contact initialization defined in the input file.

    */
    bool build_active_set(bool init = false) override;

    /*!
    \brief Build corresponding active set for master side

    */
    virtual bool build_active_set_master();

    /*!
    \brief Check mortar wear T derivatives with finite differences

    */
    void fd_check_mortar_t_deriv();

    /*!
    \brief Check mortar wear T derivatives with finite differences (Master)

    */
    void fd_check_mortar_t_master_deriv();

    /*!
    \brief Check mortar wear E derivatives with finite differences

    */
    void fd_check_mortar_e_deriv();

    /*!
    \brief Check mortar wear E derivatives with finite differences (for master)

    */
    void fd_check_mortar_e_master_deriv();

    /*!
    \brief Check mortar wear T derivatives with finite differences
      --> for wear condition

    */
    void fd_check_deriv_t_d(Core::LinAlg::SparseMatrix& lintdis);

    /*!
    \brief Check mortar wear T derivatives with finite differences
      --> for wear condition (Master)

    */
    void fd_check_deriv_t_d_master(Core::LinAlg::SparseMatrix& lintdis);

    /*!
    \brief Check mortar wear E derivatives with finite differences
      --> for wear condition

    */
    void fd_check_deriv_e_d(Core::LinAlg::SparseMatrix& linedis);

    /*!
    \brief Check mortar wear E derivatives with finite differences
      --> for wear condition (Master)

    */
    void fd_check_deriv_e_d_master(Core::LinAlg::SparseMatrix& linedis);
    /*!
    \brief Check weighted gap g derivatives with finite differences

    */
    void fd_check_gap_deriv();

    /*!
    \brief Check weighted gap g derivatives with finite differences

    */
    void fd_check_gap_deriv_w();

    /*!
    \brief Check weighted wear ~w derivatives with finite differences
           derivation w.r.t. displ.

    */
    void fd_check_wear_deriv();

    /*!
    \brief Check weighted wear ~w derivatives with finite differences
           derivation w.r.t. lagr.-mult.

    */
    void fd_check_wear_deriv_lm();

    /*!
    \brief Check slip condition derivatives with finite differences

    */
    virtual void fd_check_slip_deriv(Core::LinAlg::SparseMatrix& linslipLMglobal,
        Core::LinAlg::SparseMatrix& linslipDISglobal, Core::LinAlg::SparseMatrix& linslipWglobal);

    /*!
    \brief Assemble inactive rhs (incremental delta_w_)
    */
    virtual void assemble_inactive_wear_rhs(Epetra_Vector& inactiverhs);

    /*!
    \brief Assemble inactive rhs (incremental delta_w_)
    */
    virtual void assemble_inactive_wear_rhs_master(Epetra_FEVector& inactiverhs);

    /*!
    \brief Assemble wear-cond. rhs
    */
    virtual void assemble_wear_cond_rhs(Epetra_Vector& rhs);

    /*!
    \brief Assemble wear-cond. rhs
    */
    virtual void assemble_wear_cond_rhs_master(Epetra_FEVector& rhs);

    /*!
    \brief Initialize / reset interface for contact

    Derived version!

    */
    void initialize() final;


    /*!
    \brief Returning dofs for both-sided wear mapping

    */
    virtual Teuchos::RCP<Epetra_Map> involved_dofs() const { return involveddofs_; }

    virtual Teuchos::RCP<Epetra_Map> involved_nodes() const { return involvednodes_; }

    /*!
    \brief Set element areas

    Derived version!

    */
    void split_slave_dofs();
    void split_master_dofs();
    /*!
    \brief Set element areas

    Derived version!

    */
    void set_element_areas() final;

    /*!
    \brief Evaluate nodal normals

    */
    void evaluate_nodal_normals() const final;


    /*!
    \brief Evaluate nodal normals

    */
    void export_nodal_normals() const final;

    /*!
    \brief Update interface Wear variable sets

    This update is usually only done ONCE in the initialization phase
    and sets up the wear unknowns (only dofs) for the whole
    simulation.

    */
    virtual void update_w_sets(int offset_if, int maxdofwear, bool bothdiscr);

    /*!
    \brief Get map of slave wear dofs (Filled()==true is prerequisite)

    */
    virtual Teuchos::RCP<Epetra_Map> w_dofs() const
    {
      if (filled())
        return wdofmap_;
      else
        FOUR_C_THROW("CONTACT::WearInterface::fill_complete was not called");
      exit(EXIT_FAILURE);  // calm down the compiler
    }

    /*!
    \brief Get map of master wear dofs (Filled()==true is prerequisite)

    */
    virtual Teuchos::RCP<Epetra_Map> wm_dofs() const
    {
      if (filled())
        return wmdofmap_;
      else
        FOUR_C_THROW("CONTACT::WearInterface::fill_complete was not called");
      exit(EXIT_FAILURE);  // calm down the compiler
    }

    /*!
    \brief Get map of Lagrange multiplier dofs (Filled()==true is prerequisite)

    */
    virtual Teuchos::RCP<Epetra_Map> sn_dofs() const
    {
      if (filled())
        return sndofmap_;
      else
        FOUR_C_THROW("CONTACT::WearInterface::fill_complete was not called");
      exit(EXIT_FAILURE);  // calm down the compiler
    }

    /*!
    \brief Get map of Lagrange multiplier dofs (Filled()==true is prerequisite)

    */
    virtual Teuchos::RCP<Epetra_Map> mn_dofs() const
    {
      if (filled())
        return mndofmap_;
      else
        FOUR_C_THROW("CONTACT::WearInterface::fill_complete was not called");
      exit(EXIT_FAILURE);  // calm down the compiler
    }

    /*!
    \brief Get row map of active nodes (Filled()==true is prerequisite)

    */
    virtual Teuchos::RCP<Epetra_Map> active_master_nodes() const
    {
      if (filled())
        return activmasternodes_;
      else
        FOUR_C_THROW("CONTACT::Interface::fill_complete was not called");
      exit(EXIT_FAILURE);  // calm down the compiler
    }

    /*!
    \brief Get row map of active nodes (Filled()==true is prerequisite)

    */
    virtual Teuchos::RCP<Epetra_Map> slip_master_nodes() const
    {
      if (filled())
        return slipmasternodes_;
      else
        FOUR_C_THROW("CONTACT::Interface::fill_complete was not called");
      exit(EXIT_FAILURE);  // calm down the compiler
    }

    /*!
    \brief Get row map of active nodes (Filled()==true is prerequisite)

    */
    virtual Teuchos::RCP<Epetra_Map> slip_master_n_dofs() const
    {
      if (filled())
        return slipmn_;
      else
        FOUR_C_THROW("CONTACT::Interface::fill_complete was not called");
      exit(EXIT_FAILURE);  // calm down the compiler
    }


    /*!
    \brief Get type of wear shapefnct

    */
    Inpar::Wear::WearShape wear_shape_fcn()
    {
      return Core::UTILS::IntegralValue<Inpar::Wear::WearShape>(imortar_, "WEAR_SHAPEFCN");
    }

   private:
    /*!
    \brief initialize node and element data container

    Derived verision!

    */
    void initialize_data_container() final;


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


}  // namespace Wear

FOUR_C_NAMESPACE_CLOSE

#endif
