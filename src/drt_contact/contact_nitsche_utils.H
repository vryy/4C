/*---------------------------------------------------------------------*/
/*! \file
\brief Some helpers for nitsche contact

\level 3

\maintainer Christoph Ager

*/
/*---------------------------------------------------------------------*/
#ifndef CONTACT_NITSCHE_UTILS_H
#define CONTACT_NITSCHE_UTILS_H

#include <Epetra_FEVector.h>
#include <Epetra_CrsMatrix.h>
#include <Teuchos_RCP.hpp>
#include "../linalg/linalg_serialdensevector.H"
#include "../linalg/linalg_sparsematrix.H"
#include "../linalg/linalg_fixedsizematrix.H"
#include "../drt_mortar/mortar_element.H"
#include "../drt_fem_general/drt_utils_local_connectivity_matrices.H"
#include "../drt_lib/drt_utils_matrix_vector_enums.H"
#include <unordered_map>

namespace MORTAR
{
  class MortarElementNitscheContainer
  {
   public:
    MortarElementNitscheContainer(){};

    virtual void Clear() = 0;

    virtual void AssembleRHS(MORTAR::MortarElement* mele, DRT::UTILS::VecBlockType row,
        Teuchos::RCP<Epetra_FEVector> fc) = 0;

    virtual void AssembleMatrix(MORTAR::MortarElement* mele, DRT::UTILS::MatBlockType block,
        Teuchos::RCP<LINALG::SparseMatrix> kc) = 0;

    virtual double* rhs(int dof) = 0;
    virtual double* rhs() = 0;
    virtual double* k(int col) = 0;
    virtual double* k(int col, int dof) = 0;

    virtual double* rhs_t(int dof) = 0;
    virtual double* rhs_t() = 0;
    virtual double* k_tt(int col) = 0;
    virtual double* k_tt(int col, int dof) = 0;
    virtual double* k_td(int col) = 0;
    virtual double* k_td(int col, int dof) = 0;
    virtual double* k_dt(int col) = 0;
    virtual double* k_dt(int col, int dof) = 0;

    virtual double* rhs_p(int dof) = 0;
    virtual double* rhs_p() = 0;
    virtual double* k_pp(int col) = 0;
    virtual double* k_pp(int col, int dof) = 0;
    virtual double* k_pd(int col) = 0;
    virtual double* k_pd(int col, int dof) = 0;
    virtual double* k_dp(int col) = 0;
    virtual double* k_dp(int col, int dof) = 0;
  };

  template <DRT::Element::DiscretizationType parent_distype>
  class MortarElementNitscheDataTsi
  {
   public:
    void Clear()
    {
      rhs_t_.Clear();
      k_tt_.clear();
      k_td_.clear();
      k_dt_.clear();
    }

    static const int num_parent_disp_dof =
        DRT::UTILS::DisTypeToNumNodePerEle<parent_distype>::numNodePerElement *
        DRT::UTILS::DisTypeToDim<parent_distype>::dim;
    static const int num_parent_thermo_dof =
        DRT::UTILS::DisTypeToNumNodePerEle<parent_distype>::numNodePerElement;

    LINALG::Matrix<num_parent_thermo_dof, 1> rhs_t_;
    std::unordered_map<int, LINALG::Matrix<num_parent_thermo_dof, 1>> k_tt_;
    std::unordered_map<int, LINALG::Matrix<num_parent_thermo_dof, 1>> k_td_;
    std::unordered_map<int, LINALG::Matrix<num_parent_disp_dof, 1>> k_dt_;
  };

  template <DRT::Element::DiscretizationType parent_distype>
  class MortarElementNitscheDataPoro
  {
   public:
    void Clear()
    {
      rhs_p_.Clear();
      k_pp_.clear();
      k_pd_.clear();
      k_dp_.clear();
    }

    static const int num_parent_disp_dof =
        DRT::UTILS::DisTypeToNumNodePerEle<parent_distype>::numNodePerElement *
        DRT::UTILS::DisTypeToDim<parent_distype>::dim;
    static const int num_parent_pf_dof =
        DRT::UTILS::DisTypeToNumNodePerEle<parent_distype>::numNodePerElement *
        (DRT::UTILS::DisTypeToDim<parent_distype>::dim + 1);

    LINALG::Matrix<num_parent_pf_dof, 1> rhs_p_;
    std::unordered_map<int, LINALG::Matrix<num_parent_pf_dof, 1>> k_pp_;
    std::unordered_map<int, LINALG::Matrix<num_parent_pf_dof, 1>> k_pd_;
    std::unordered_map<int, LINALG::Matrix<num_parent_disp_dof, 1>> k_dp_;
  };


  template <DRT::Element::DiscretizationType parent_distype>
  class MortarElementNitscheData : public MortarElementNitscheContainer
  {
   public:
    virtual double* rhs(int dof) { return &rhs_(dof); }
    virtual double* rhs() { return rhs_.A(); }
    virtual double* k(int col) { return k_[col].A(); }
    virtual double* k(int col, int dof) { return &k_[col](dof); }

    virtual double* rhs_t(int dof) { return &tsi_data_.rhs_t_(dof); }
    virtual double* rhs_t() { return tsi_data_.rhs_t_.A(); }
    virtual double* k_tt(int col) { return tsi_data_.k_tt_[col].A(); }
    virtual double* k_tt(int col, int dof) { return &tsi_data_.k_tt_[col](dof); }
    virtual double* k_td(int col) { return tsi_data_.k_td_[col].A(); }
    virtual double* k_td(int col, int dof) { return &tsi_data_.k_td_[col](dof); }
    virtual double* k_dt(int col) { return tsi_data_.k_dt_[col].A(); }
    virtual double* k_dt(int col, int dof) { return &tsi_data_.k_dt_[col](dof); }

    virtual double* rhs_p(int dof) { return &poro_data_.rhs_p_(dof); }
    virtual double* rhs_p() { return poro_data_.rhs_p_.A(); }
    virtual double* k_pp(int col) { return poro_data_.k_pp_[col].A(); }
    virtual double* k_pp(int col, int dof) { return &poro_data_.k_pp_[col](dof); }
    virtual double* k_pd(int col) { return poro_data_.k_pd_[col].A(); }
    virtual double* k_pd(int col, int dof) { return &poro_data_.k_pd_[col](dof); }
    virtual double* k_dp(int col) { return poro_data_.k_dp_[col].A(); }
    virtual double* k_dp(int col, int dof) { return &poro_data_.k_dp_[col](dof); }

    virtual void AssembleRHS(MORTAR::MortarElement* mele, DRT::UTILS::VecBlockType row,
        Teuchos::RCP<Epetra_FEVector> fc);

    virtual void AssembleMatrix(MORTAR::MortarElement* mele, DRT::UTILS::MatBlockType block,
        Teuchos::RCP<LINALG::SparseMatrix> kc);

    template <int num_dof_per_node>
    void AssembleRHS(MORTAR::MortarElement* mele,
        const LINALG::Matrix<DRT::UTILS::DisTypeToNumNodePerEle<parent_distype>::numNodePerElement *
                                 num_dof_per_node,
            1>& rhs,
        std::vector<int>& dofs, Teuchos::RCP<Epetra_FEVector> fc);

    template <int num_dof_per_node>
    void AssembleMatrix(MORTAR::MortarElement* mele,
        const std::unordered_map<int,
            LINALG::Matrix<DRT::UTILS::DisTypeToNumNodePerEle<parent_distype>::numNodePerElement *
                               num_dof_per_node,
                1>>& k,
        std::vector<int>& dofs, Teuchos::RCP<LINALG::SparseMatrix> kc);


    virtual void Clear()
    {
      rhs_.Clear();
      k_.clear();
      tsi_data_.Clear();
      poro_data_.Clear();
    }

   private:
    LINALG::Matrix<DRT::UTILS::DisTypeToNumNodePerEle<parent_distype>::numNodePerElement *
                       DRT::UTILS::DisTypeToDim<parent_distype>::dim,
        1>
        rhs_;
    std::unordered_map<int,
        LINALG::Matrix<DRT::UTILS::DisTypeToNumNodePerEle<parent_distype>::numNodePerElement *
                           DRT::UTILS::DisTypeToDim<parent_distype>::dim,
            1>>
        k_;
    MORTAR::MortarElementNitscheDataTsi<parent_distype> tsi_data_;
    MORTAR::MortarElementNitscheDataPoro<parent_distype> poro_data_;
  };

}  // namespace MORTAR

#endif  // CONTACT_NITSCHE_UTILS_H
