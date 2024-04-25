/*---------------------------------------------------------------------*/
/*! \file
\brief Some helpers for nitsche contact

\level 3


*/
/*---------------------------------------------------------------------*/
#ifndef FOUR_C_CONTACT_NITSCHE_UTILS_HPP
#define FOUR_C_CONTACT_NITSCHE_UTILS_HPP

#include "4C_config.hpp"

#include "4C_contact_utils.hpp"
#include "4C_discretization_fem_general_utils_local_connectivity_matrices.hpp"
#include "4C_linalg_fixedsizematrix.hpp"
#include "4C_linalg_serialdensevector.hpp"
#include "4C_linalg_sparsematrix.hpp"
#include "4C_mortar_element.hpp"

#include <Epetra_CrsMatrix.h>
#include <Epetra_FEVector.h>
#include <Teuchos_RCP.hpp>

#include <unordered_map>

FOUR_C_NAMESPACE_OPEN

namespace MORTAR
{
  class ElementNitscheContainer
  {
   public:
    /**
     * Virtual destructor.
     */
    virtual ~ElementNitscheContainer() = default;

    ElementNitscheContainer() = default;

    virtual void Clear() = 0;

    virtual void AssembleRHS(MORTAR::Element* mele, CONTACT::VecBlockType row,
        Teuchos::RCP<Epetra_FEVector> fc) const = 0;

    virtual void AssembleMatrix(MORTAR::Element* mele, CONTACT::MatBlockType block,
        Teuchos::RCP<CORE::LINALG::SparseMatrix> kc) const = 0;

    virtual double* Rhs(int dof) = 0;
    virtual double* Rhs() = 0;
    virtual double* K(int col) = 0;
    virtual double* K(int col, int dof) = 0;

    virtual double* RhsT(int dof) = 0;
    virtual double* RhsT() = 0;
    virtual double* Ktt(int col) = 0;
    virtual double* Ktd(int col) = 0;
    virtual double* Kdt(int col) = 0;

    virtual double* RhsP(int dof) = 0;
    virtual double* Kpp(int col) = 0;
    virtual double* Kpd(int col) = 0;
    virtual double* Kdp(int col) = 0;

    virtual double* RhsS(int dof) = 0;
    virtual double* Kss(int col) = 0;
    virtual double* Ksd(int col) = 0;
    virtual double* Kds(int col) = 0;

    virtual double* RhsE(int dof) = 0;
    virtual double* Kee(int col) = 0;
    virtual double* Ked(int col) = 0;
    virtual double* Ked(int col, int dof) = 0;
    virtual double* Kde(int col) = 0;
  };

  template <CORE::FE::CellType parent_distype>
  class ElementNitscheDataTsi
  {
   public:
    void Clear()
    {
      rhs_t_.Clear();
      k_tt_.clear();
      k_td_.clear();
      k_dt_.clear();
    }

    static constexpr int num_parent_disp_dof =
        CORE::FE::num_nodes<parent_distype> * CORE::FE::dim<parent_distype>;
    static constexpr int num_parent_thermo_dof = CORE::FE::num_nodes<parent_distype>;

    CORE::LINALG::Matrix<num_parent_thermo_dof, 1> rhs_t_;
    std::unordered_map<int, CORE::LINALG::Matrix<num_parent_thermo_dof, 1>> k_tt_;
    std::unordered_map<int, CORE::LINALG::Matrix<num_parent_thermo_dof, 1>> k_td_;
    std::unordered_map<int, CORE::LINALG::Matrix<num_parent_disp_dof, 1>> k_dt_;
  };

  template <CORE::FE::CellType parent_distype>
  class ElementNitscheDataPoro
  {
   public:
    void Clear()
    {
      rhs_p_.Clear();
      k_pp_.clear();
      k_pd_.clear();
      k_dp_.clear();
    }

    static constexpr int num_parent_disp_dof =
        CORE::FE::num_nodes<parent_distype> * CORE::FE::dim<parent_distype>;
    static constexpr int num_parent_pf_dof =
        CORE::FE::num_nodes<parent_distype> * (CORE::FE::dim<parent_distype> + 1);

    CORE::LINALG::Matrix<num_parent_pf_dof, 1> rhs_p_;
    std::unordered_map<int, CORE::LINALG::Matrix<num_parent_pf_dof, 1>> k_pp_;
    std::unordered_map<int, CORE::LINALG::Matrix<num_parent_pf_dof, 1>> k_pd_;
    std::unordered_map<int, CORE::LINALG::Matrix<num_parent_disp_dof, 1>> k_dp_;
  };

  template <CORE::FE::CellType parent_distype>
  class ElementNitscheDataSsi
  {
   public:
    void Clear()
    {
      rhs_s_.Clear();
      k_ss_.clear();
      k_sd_.clear();
      k_ds_.clear();
    }

    static constexpr int num_parent_disp_dof =
        CORE::FE::num_nodes<parent_distype> * CORE::FE::dim<parent_distype>;
    static constexpr int num_parent_scatra_dof = CORE::FE::num_nodes<parent_distype>;

    CORE::LINALG::Matrix<num_parent_scatra_dof, 1> rhs_s_;
    std::unordered_map<int, CORE::LINALG::Matrix<num_parent_scatra_dof, 1>> k_ss_;
    std::unordered_map<int, CORE::LINALG::Matrix<num_parent_scatra_dof, 1>> k_sd_;
    std::unordered_map<int, CORE::LINALG::Matrix<num_parent_disp_dof, 1>> k_ds_;
  };

  template <CORE::FE::CellType parent_distype>
  class ElementNitscheDataSsiElch
  {
   public:
    void Clear()
    {
      rhs_e_.Clear();
      k_ee_.clear();
      k_ed_.clear();
      k_de_.clear();
    }

    static constexpr int num_parent_disp_dof =
        CORE::FE::num_nodes<parent_distype> * CORE::FE::dim<parent_distype>;
    static constexpr int num_parent_elch_dof = CORE::FE::num_nodes<parent_distype> * 2;

    CORE::LINALG::Matrix<num_parent_elch_dof, 1> rhs_e_;
    std::unordered_map<int, CORE::LINALG::Matrix<num_parent_elch_dof, 1>> k_ee_;
    std::unordered_map<int, CORE::LINALG::Matrix<num_parent_elch_dof, 1>> k_ed_;
    std::unordered_map<int, CORE::LINALG::Matrix<num_parent_disp_dof, 1>> k_de_;
  };

  template <CORE::FE::CellType parent_distype>
  class ElementNitscheData : public ElementNitscheContainer
  {
    using VectorType =
        CORE::LINALG::Matrix<CORE::FE::num_nodes<parent_distype> * CORE::FE::dim<parent_distype>,
            1>;

   public:
    const VectorType& RhsVec() { return rhs_; }
    double* Rhs(int dof) override { return &rhs_(dof); }
    double* Rhs() override { return rhs_.A(); }
    double* K(int col) override { return k_[col].A(); }
    double* K(int col, int dof) override { return &k_[col](dof); }

    double* RhsT(int dof) override { return &tsi_data_.rhs_t_(dof); }
    double* RhsT() override { return tsi_data_.rhs_t_.A(); }
    double* Ktt(int col) override { return tsi_data_.k_tt_[col].A(); }
    double* Ktd(int col) override { return tsi_data_.k_td_[col].A(); }
    double* Kdt(int col) override { return tsi_data_.k_dt_[col].A(); }

    double* RhsP(int dof) override { return &poro_data_.rhs_p_(dof); }
    double* Kpp(int col) override { return poro_data_.k_pp_[col].A(); }
    double* Kpd(int col) override { return poro_data_.k_pd_[col].A(); }
    double* Kdp(int col) override { return poro_data_.k_dp_[col].A(); }

    double* RhsS(int dof) override { return &ssi_data_.rhs_s_(dof); }
    double* Kss(int col) override { return ssi_data_.k_ss_[col].A(); }
    double* Ksd(int col) override { return ssi_data_.k_sd_[col].A(); }
    double* Kds(int col) override { return ssi_data_.k_ds_[col].A(); }

    double* RhsE(int dof) override { return &ssi_elch_data_.rhs_e_(dof); }
    double* Kee(int col) override { return ssi_elch_data_.k_ee_[col].A(); }
    double* Ked(int col) override { return ssi_elch_data_.k_ed_[col].A(); }
    double* Ked(int col, int dof) override { return &ssi_elch_data_.k_ed_[col](dof); }
    double* Kde(int col) override { return ssi_elch_data_.k_de_[col].A(); }

    void AssembleRHS(MORTAR::Element* mele, CONTACT::VecBlockType row,
        Teuchos::RCP<Epetra_FEVector> fc) const override;

    void AssembleMatrix(MORTAR::Element* mele, CONTACT::MatBlockType block,
        Teuchos::RCP<CORE::LINALG::SparseMatrix> kc) const override;

    template <int num_dof_per_node>
    void AssembleRHS(MORTAR::Element* mele,
        const CORE::LINALG::Matrix<CORE::FE::num_nodes<parent_distype> * num_dof_per_node, 1>& rhs,
        std::vector<int>& dofs, Teuchos::RCP<Epetra_FEVector> fc) const;

    template <int num_dof_per_node>
    void AssembleMatrix(MORTAR::Element* mele,
        const std::unordered_map<int,
            CORE::LINALG::Matrix<CORE::FE::num_nodes<parent_distype> * num_dof_per_node, 1>>& k,
        std::vector<int>& dofs, Teuchos::RCP<CORE::LINALG::SparseMatrix> kc) const;


    void Clear() override
    {
      rhs_.Clear();
      k_.clear();
      tsi_data_.Clear();
      poro_data_.Clear();
      ssi_data_.Clear();
      ssi_elch_data_.Clear();
    }

   private:
    VectorType rhs_;
    std::unordered_map<int, VectorType> k_;
    MORTAR::ElementNitscheDataTsi<parent_distype> tsi_data_;
    MORTAR::ElementNitscheDataPoro<parent_distype> poro_data_;
    MORTAR::ElementNitscheDataSsi<parent_distype> ssi_data_;
    MORTAR::ElementNitscheDataSsiElch<parent_distype> ssi_elch_data_;
  };

}  // namespace MORTAR

FOUR_C_NAMESPACE_CLOSE

#endif
