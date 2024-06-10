/*----------------------------------------------------------------------------*/
/*! \file
\brief (augmented) contact integration policies

\level 3

*/
/*----------------------------------------------------------------------------*/

#ifndef FOUR_C_CONTACT_AUG_INTEGRATOR_POLICY_HPP
#define FOUR_C_CONTACT_AUG_INTEGRATOR_POLICY_HPP

#include "4C_config.hpp"

#include "4C_contact_aug_contact_integrator_utils.hpp"
#include "4C_contact_aug_timemonitor.hpp"
#include "4C_fem_general_node.hpp"

FOUR_C_NAMESPACE_OPEN

namespace CONTACT
{
  namespace Aug
  {
    // forward declarations
    class NodeDataContainer;

    /*--------------------------------------------------------------------------*/
    template <unsigned probdim, Core::FE::CellType slavetype>
    class BaseSlaveIntPolicy
    {
     protected:
      static constexpr unsigned SLAVEDIM = Core::FE::dim<slavetype>;
      static constexpr unsigned SLAVENUMNODE = Core::FE::num_nodes<slavetype>;

     public:
      /// constructor
      BaseSlaveIntPolicy(){/* empty */};


      double unit_slave_element_normal(const Mortar::Element& sele,
          const Core::LinAlg::Matrix<3, 2>& tau,
          Core::LinAlg::Matrix<probdim, 1>& unit_normal) const;

      /// @name First derivatives
      /// @{

      void deriv1st_non_unit_slave_element_normal(const Mortar::Element& sele,
          const Core::LinAlg::Matrix<probdim, SLAVENUMNODE, int>& nodal_dofs,
          const Core::LinAlg::Matrix<SLAVEDIM, SLAVENUMNODE>& deriv,
          const Core::LinAlg::Matrix<3, 2>& tau, Deriv1stVecMap& d_non_unit_normal) const;

      void deriv1st_unit_slave_element_normal(const Core::LinAlg::Matrix<probdim, 1>& unit_normal,
          const double length_n_inv, const Deriv1stVecMap& d_non_unit_normal,
          Deriv1stVecMap& d_unit_normal, const bool reset = true) const;

      void deriv1st_jacobian(const Core::LinAlg::Matrix<probdim, 1>& unit_normal,
          const Deriv1stVecMap& d_non_unit_normal, Deriv1stMap& d_jac) const;

      /// @}

      /// @name Second derivatives
      /// @{

      void Deriv2nd_Jacobian(const Deriv1stVecMap& d_unit_normal,
          const Deriv1stVecMap& d_non_unit_normal,
          const Core::LinAlg::Matrix<probdim, 1>& unit_normal,
          const Deriv2ndVecMap& dd_non_unit_normal, Deriv2ndMap& dd_jac) const;

      void deriv2nd_non_unit_slave_element_normal(const Mortar::Element& sele,
          const Core::LinAlg::Matrix<probdim, SLAVENUMNODE, int>& nodal_dofs,
          const Core::LinAlg::Matrix<SLAVEDIM, SLAVENUMNODE>& deriv,
          Deriv2ndVecMap& dd_non_unit_normal) const;

      void deriv2nd_unit_slave_element_normal(const Core::LinAlg::Matrix<probdim, 1>& unit_normal,
          const double length_n_inv, const Deriv1stVecMap& d_non_unit_normal,
          const Deriv1stVecMap& d_unit_normal, const Deriv2ndVecMap& dd_non_unit_normal,
          Deriv2ndVecMap& dd_unit_normal) const;

      /// @}

      void AveragedNormalAtXi(Mortar::Element& sele,
          const Core::LinAlg::Matrix<SLAVENUMNODE, 1>& sval,
          Core::LinAlg::Matrix<probdim, 1>& snormal) const;

      /// complete data in nodal data container after successful integration
      void CompleteNodeData(Mortar::Element& sele) const;

     private:
      void inner_product_of_vector_and_deriv1st_vector(const Core::LinAlg::Matrix<probdim, 1>& vec,
          const Deriv1stVecMap& d_vec, Deriv1stMap& dvec_vec) const;

      void projection_into_tangential_plain(const Core::LinAlg::Matrix<probdim, 1>& unit_normal,
          Core::LinAlg::Matrix<probdim, probdim>& tproj_mat) const;

     protected:
      /// current Gauss point ID
      int gp_id_ = -1;
      mutable TimeMonitor<CONTACT::Aug::TimeID> timer_ = TimeMonitor<CONTACT::Aug::TimeID>();
    };

    /*--------------------------------------------------------------------------*/
    template <unsigned probdim, Core::FE::CellType slavetype, Core::FE::CellType mastertype>
    class BaseIntPolicy : public BaseSlaveIntPolicy<probdim, slavetype>
    {
      typedef BaseSlaveIntPolicy<probdim, slavetype> my;

     protected:
      static constexpr unsigned MASTERDIM = Core::FE::dim<mastertype>;
      static constexpr unsigned MASTERNUMNODE = Core::FE::num_nodes<mastertype>;

     public:
      /// constructor
      BaseIntPolicy(){/* empty */};


      void LMatrixInverse(const Core::LinAlg::Matrix<3, 2>& mtau,
          const Core::LinAlg::Matrix<probdim, 1>& snormal,
          Core::LinAlg::Matrix<probdim, probdim>& lmat_inv) const;

      void deriv1st_m_xi_gp(const Core::LinAlg::Matrix<probdim, probdim>& lmat_inv,
          Mortar::Element& sele, Mortar::Element& mele,
          const Core::LinAlg::Matrix<my::SLAVENUMNODE, 1>& sval,
          const Core::LinAlg::Matrix<MASTERNUMNODE, 1>& mval, const double alpha,
          Deriv1stVecMap& d_mxi, Deriv1stMap& d_alpha) const;

      void add_deriv2nd_ma_displ(const Core::LinAlg::Matrix<probdim, probdim>& lmat_inv,
          const Core::LinAlg::Matrix<probdim, MASTERNUMNODE, int>& mnodal_dofs,
          const Core::LinAlg::Matrix<MASTERDIM, MASTERNUMNODE>& mderiv,
          const Deriv1stVecMap& d_mxigp, Deriv2ndVecMap& dd_mxigp) const;

      void add_deriv1st_ma_metric(const Core::LinAlg::Matrix<probdim, probdim>& lmat_inv,
          const Core::LinAlg::Matrix<probdim, MASTERNUMNODE, int>& mnodal_dofs,
          const Core::LinAlg::Matrix<3, 2>& mtau,
          const Core::LinAlg::Matrix<MASTERDIM, MASTERNUMNODE>& mderiv,
          const Core::LinAlg::Matrix<3, MASTERNUMNODE>& mderiv2nd,
          const Core::LinAlg::Matrix<3, MASTERNUMNODE>& mcoord, const Deriv1stVecMap& d_mxigp,
          Deriv2ndVecMap& dd_mxigp) const;

      void add_deriv1st_alpha_deriv1st_normal(
          const Core::LinAlg::Matrix<probdim, probdim>& lmat_inv, const Deriv1stMap& d_alpha,
          const Mortar::Element& sele, const Core::LinAlg::Matrix<my::SLAVENUMNODE, 1>& sval,
          Deriv2ndVecMap& dd_mxigp) const;

      void add_alpha_deriv2nd_normal(const Core::LinAlg::Matrix<probdim, probdim>& lmat_inv,
          const double alpha, const Mortar::Element& sele,
          const Core::LinAlg::Matrix<my::SLAVENUMNODE, 1>& sval, Deriv2ndVecMap& dd_mxigp) const;

      void Deriv1st_GapN_Sl(const Core::Nodes::Node* const* snodes,
          const Core::LinAlg::Matrix<my::SLAVENUMNODE, 1>& sval, const double* gpn,
          Deriv1stMap& deriv_gapn_sl) const;

      void Deriv1st_GapN_Ma(const Core::Nodes::Node* const* mnodes,
          const Core::LinAlg::Matrix<MASTERNUMNODE, 1>& mval, const double* gpn,
          const Core::LinAlg::Matrix<3, 2>& mtau, const Deriv1stVecMap& d_mxigp,
          Deriv1stMap& deriv_gapn_ma) const;

      void add_deriv1st_gap_n_contributions(Core::Nodes::Node* const* snodes, const double scale,
          const Core::LinAlg::Matrix<my::SLAVENUMNODE, 1>& lmval, const Deriv1stMap& d_gapn_sl,
          const Deriv1stMap& d_gapn_ma) const;

      void add_deriv1st_jacobian_contributions(Core::Nodes::Node* const* snodes, const double wgt,
          const Core::LinAlg::Matrix<my::SLAVENUMNODE, 1>& lmval, const double gapn_sl,
          const double gapn_ma, const Deriv1stMap& d_jac) const;

      void add_var_gap_n_lin_jac(Mortar::Element& sele,
          const Core::LinAlg::Matrix<my::SLAVENUMNODE, 1>& lmval, const double wgt,
          const Deriv1stMap& d_gapn_sl, const Deriv1stMap& d_gapn_ma,
          const Deriv1stMap& d_jac) const;

      void add_var_jac_lin_gap_n(Mortar::Element& sele,
          const Core::LinAlg::Matrix<my::SLAVENUMNODE, 1>& lmval, const double wgt,
          const Deriv1stMap& d_gapn_sl, const Deriv1stMap& d_gapn_ma,
          const Deriv1stMap& d_jac) const;

      void add_gap_n_deriv2nd_jac(Mortar::Element& sele,
          const Core::LinAlg::Matrix<my::SLAVENUMNODE, 1>& lmval, const double wgt,
          const double gapn_sl, const double gapn_ma, const Deriv2ndMap& dd_jac) const;

      /// @name empty functions, if no debug policy is in use
      /// @{

      void Get_Debug(...) const { return; }

      void Get_Deriv1st_Debug(...) const { return; }

      void Get_Deriv2nd_Debug(...) const { return; }

      /// @}
    };

    /*--------------------------------------------------------------------------*/
    template <unsigned probdim, Core::FE::CellType slavetype, Core::FE::CellType mastertype>
    class IncompleteIntPolicy : public BaseIntPolicy<probdim, slavetype, mastertype>
    {
      typedef BaseIntPolicy<probdim, slavetype, mastertype> my;

     public:
      /// constructor
      IncompleteIntPolicy(){/* empty */};


      void get_deriv2nd_jacobian(const Mortar::Element& sele,
          const Core::LinAlg::Matrix<probdim, my::SLAVENUMNODE, int>& nodal_dofs,
          const Core::LinAlg::Matrix<my::SLAVEDIM, my::SLAVENUMNODE>& deriv,
          const Core::LinAlg::Matrix<probdim, 1>& unit_normal, const double length_n_inv,
          const Deriv1stVecMap& d_non_unit_normal, Deriv2ndMap& dd_jac) const;

      void Get_Deriv2nd_MXiGP(const Core::LinAlg::Matrix<probdim, probdim>& lmat_inv,
          Mortar::Element& sele, Mortar::Element& mele,
          const Core::LinAlg::Matrix<my::SLAVENUMNODE, 1>& sval,
          const Core::LinAlg::Matrix<my::MASTERNUMNODE, 1>& mval,
          const Core::LinAlg::Matrix<my::MASTERDIM, my::MASTERNUMNODE>& mderiv,
          const Core::LinAlg::Matrix<3, my::MASTERNUMNODE>& mderiv2,
          const Core::LinAlg::Matrix<3, 2>& mtau, const double* mxi, const double alpha,
          const Deriv1stVecMap& d_mxigp, const Deriv1stMap& d_alpha, Deriv2ndVecMap& dd_mxigp) const
      {
        // do nothing
        Core::Gen::reset(2, 0, dd_mxigp);
      }

      void Get_Deriv1st_GapN(Mortar::Element& sele, Mortar::Element& mele,
          const Core::LinAlg::Matrix<my::SLAVENUMNODE, 1>& sval,
          const Core::LinAlg::Matrix<my::MASTERNUMNODE, 1>& mval, const double* gpn,
          const Core::LinAlg::Matrix<3, 2>& mtau, const Deriv1stVecMap& d_mxigp,
          Deriv1stMap& deriv_gapn_sl, Deriv1stMap& deriv_gapn_ma) const;

      void Get_Deriv1st_WGap(Mortar::Element& sele,
          const Core::LinAlg::Matrix<my::SLAVENUMNODE, 1>& lmval, const double gapn_sl,
          const double gapn_ma, const double wgt, const double jac, const Deriv1stMap& d_jac,
          const Deriv1stMap& d_gapn_sl, const Deriv1stMap& d_gapn_ma) const;

      void get_deriv1st_w_gap_complete(const int linsize, Mortar::Element& sele,
          Mortar::Element& mele, const Core::LinAlg::Matrix<my::SLAVENUMNODE, 1>& sval,
          const Core::LinAlg::Matrix<my::MASTERNUMNODE, 1>& mval,
          const Core::LinAlg::Matrix<my::SLAVENUMNODE, 1>& lmval, const double* gpn,
          const Core::LinAlg::Matrix<3, 2>& mtau, const Deriv1stVecMap& d_mxigp,
          const double gapn_sl, const double gapn_ma, const double wgt, const double jac,
          const Deriv1stMap& d_jac) const;

      void Get_Deriv2nd_WGap(Mortar::Element& sele, Mortar::Element& mele,
          const Core::LinAlg::Matrix<my::SLAVENUMNODE, 1>& sval,
          const Core::LinAlg::Matrix<my::MASTERNUMNODE, 1>& mval,
          const Core::LinAlg::Matrix<my::SLAVENUMNODE, 1>& lmval,
          const Core::LinAlg::Matrix<my::MASTERDIM, my::MASTERNUMNODE>& mderiv,
          const Core::LinAlg::Matrix<3, my::MASTERNUMNODE>& mderiv2nd,
          const Core::LinAlg::Matrix<3, 2>& mtau, const double* gpn, const double wgt,
          const double gapn_sl, const double gapn_ma, const double jac, const Deriv1stMap& d_jac,
          const Deriv2ndMap& dd_jac, const Deriv1stVecMap& d_mxigp, const Deriv2ndVecMap& dd_mxigp,
          const Deriv1stVecMap& d_n_unit, const Deriv2ndVecMap& dd_n_unit,
          const Deriv1stMap& d_gapn_sl, const Deriv1stMap& d_gapn_ma) const;

      void get_deriv1st_w_gap_n_error(const Mortar::Element& sele,
          const std::vector<unsigned>& active_nlids,
          const Core::LinAlg::Matrix<my::SLAVENUMNODE, 1>& lmval, const double* gpn,
          const double gapn_sl, const double gapn_ma, const double wgt, const double jacslave,
          const Deriv1stMap& d_jac, const Core::LinAlg::Matrix<3, 2>& mtau,
          const Deriv1stVecMap& d_mxigp, Deriv1stMap& d_gapn_ma,
          std::unordered_map<int, Deriv1stMap>& error_ma,
          std::unordered_map<int, Deriv1stMap>& error_jac) const;

     protected:
      void add_deriv1st_gap_n_deriv1st_jac(Mortar::Element& sele,
          const Core::LinAlg::Matrix<my::SLAVENUMNODE, 1>& lmval, const double wgt,
          const Deriv1stMap& d_gapn_sl, const Deriv1stMap& d_gapn_ma,
          const Deriv1stMap& d_jac) const;

      void add_jac_deriv2nd_gap_n(Mortar::Element& sele, Mortar::Element& mele,
          const Core::LinAlg::Matrix<my::SLAVENUMNODE, 1>& sval,
          const Core::LinAlg::Matrix<my::MASTERNUMNODE, 1>& mval,
          const Core::LinAlg::Matrix<my::SLAVENUMNODE, 1>& lmval,
          const Core::LinAlg::Matrix<my::MASTERDIM, my::MASTERNUMNODE>& mderiv,
          const Core::LinAlg::Matrix<3, 2>& mtau, const double* gpn, const double wgt,
          const double jac, const Deriv1stVecMap& d_mxigp, const Deriv1stVecMap& d_n_unit,
          const Deriv2ndVecMap& dd_n_unit) const;
    };

    /*--------------------------------------------------------------------------*/
    template <unsigned probdim, Core::FE::CellType slavetype, Core::FE::CellType mastertype>
    class CompleteIntPolicy : public BaseIntPolicy<probdim, slavetype, mastertype>
    {
      typedef BaseIntPolicy<probdim, slavetype, mastertype> my;

     public:
      /// constructor
      CompleteIntPolicy(){/* empty */};


      void get_deriv2nd_jacobian(const Mortar::Element& sele,
          const Core::LinAlg::Matrix<probdim, my::SLAVENUMNODE, int>& nodal_dofs,
          const Core::LinAlg::Matrix<my::SLAVEDIM, my::SLAVENUMNODE>& deriv,
          const Core::LinAlg::Matrix<probdim, 1>& unit_normal, const double length_n_inv,
          const Deriv1stVecMap& d_non_unit_normal, Deriv2ndMap& dd_jac) const;

      void Get_Deriv2nd_MXiGP(const Core::LinAlg::Matrix<probdim, probdim>& lmat_inv,
          Mortar::Element& sele, Mortar::Element& mele,
          const Core::LinAlg::Matrix<my::SLAVENUMNODE, 1>& sval,
          const Core::LinAlg::Matrix<my::MASTERNUMNODE, 1>& mval,
          const Core::LinAlg::Matrix<my::MASTERDIM, my::MASTERNUMNODE>& mderiv,
          const Core::LinAlg::Matrix<3, my::MASTERNUMNODE>& mderiv2,
          const Core::LinAlg::Matrix<3, 2>& mtau, const double* mxi, const double alpha,
          const Deriv1stVecMap& d_mxigp, const Deriv1stMap& d_alpha,
          Deriv2ndVecMap& dd_mxigp) const;

      void Get_Deriv1st_GapN(Mortar::Element& sele, Mortar::Element& mele,
          const Core::LinAlg::Matrix<my::SLAVENUMNODE, 1>& sval,
          const Core::LinAlg::Matrix<my::MASTERNUMNODE, 1>& mval, const double* gpn,
          const Core::LinAlg::Matrix<3, 2>& mtau, const Deriv1stVecMap& d_mxigp,
          Deriv1stMap& deriv_gapn_sl, Deriv1stMap& deriv_gapn_ma) const;

      void Get_Deriv1st_WGap(Mortar::Element& sele,
          const Core::LinAlg::Matrix<my::SLAVENUMNODE, 1>& lmval, const double gapn_sl,
          const double gapn_ma, const double wgt, const double jac, const Deriv1stMap& d_jac,
          const Deriv1stMap& d_gapn_sl, const Deriv1stMap& d_gapn_ma) const;

      void get_deriv1st_w_gap_complete(const int linsize, Mortar::Element& sele,
          Mortar::Element& mele, const Core::LinAlg::Matrix<my::SLAVENUMNODE, 1>& sval,
          const Core::LinAlg::Matrix<my::MASTERNUMNODE, 1>& mval,
          const Core::LinAlg::Matrix<my::SLAVENUMNODE, 1>& lmval, const double* gpn,
          const Core::LinAlg::Matrix<3, 2>& mtau, const Deriv1stVecMap& d_mxigp,
          const double gapn_sl, const double gapn_ma, const double wgt, const double jac,
          const Deriv1stMap& d_jac) const;

      void Get_Deriv2nd_WGap(Mortar::Element& sele, Mortar::Element& mele,
          const Core::LinAlg::Matrix<my::SLAVENUMNODE, 1>& sval,
          const Core::LinAlg::Matrix<my::MASTERNUMNODE, 1>& mval,
          const Core::LinAlg::Matrix<my::SLAVENUMNODE, 1>& lmval,
          const Core::LinAlg::Matrix<my::MASTERDIM, my::MASTERNUMNODE>& mderiv,
          const Core::LinAlg::Matrix<3, my::MASTERNUMNODE>& mderiv2nd,
          const Core::LinAlg::Matrix<3, 2>& mtau, const double* gpn, const double wgt,
          const double gapn_sl, const double gapn_ma, const double jac, const Deriv1stMap& d_jac,
          const Deriv2ndMap& dd_jac, const Deriv1stVecMap& d_mxigp, const Deriv2ndVecMap& dd_mxigp,
          const Deriv1stVecMap& d_n_unit, const Deriv2ndVecMap& dd_n_unit,
          const Deriv1stMap& d_gapn_sl, const Deriv1stMap& d_gapn_ma) const;

      void get_deriv1st_w_gap_n_error(const Mortar::Element& sele,
          const std::vector<unsigned>& active_nlids,
          const Core::LinAlg::Matrix<my::SLAVENUMNODE, 1>& lmval, const double* gpn,
          const double gapn_sl, const double gapn_ma, const double wgt, const double jacslave,
          const Deriv1stMap& d_jac, const Core::LinAlg::Matrix<3, 2>& mtau,
          const Deriv1stVecMap& d_mxigp, Deriv1stMap& d_gapn_ma,
          std::unordered_map<int, Deriv1stMap>& error_ma,
          std::unordered_map<int, Deriv1stMap>& error_jac) const;

     protected:
      void add_deriv1st_gap_n_deriv1st_jac(Mortar::Element& sele,
          const Core::LinAlg::Matrix<my::SLAVENUMNODE, 1>& lmval, const double wgt,
          const Deriv1stMap& d_gapn_sl, const Deriv1stMap& d_gapn_ma,
          const Deriv1stMap& d_jac) const;

      void add_jac_deriv2nd_gap_n(Mortar::Element& sele, Mortar::Element& mele,
          const Core::LinAlg::Matrix<my::SLAVENUMNODE, 1>& sval,
          const Core::LinAlg::Matrix<my::MASTERNUMNODE, 1>& mval,
          const Core::LinAlg::Matrix<my::SLAVENUMNODE, 1>& lmval,
          const Core::LinAlg::Matrix<my::MASTERDIM, my::MASTERNUMNODE>& mderiv,
          const Core::LinAlg::Matrix<3, my::MASTERNUMNODE>& mderiv2nd,
          const Core::LinAlg::Matrix<3, 2>& mtau, const double* gpn, const double wgt,
          const double jac, const Deriv1stVecMap& d_mxigp, const Deriv2ndVecMap& dd_mxigp,
          const Deriv1stVecMap& d_n_unit, const Deriv2ndVecMap& dd_n_unit) const;
    };

    /*--------------------------------------------------------------------------*/
    template <unsigned probdim, Core::FE::CellType slavetype, Core::FE::CellType mastertype>
    class DebugIncompleteIntPolicy : public IncompleteIntPolicy<probdim, slavetype, mastertype>
    {
      typedef IncompleteIntPolicy<probdim, slavetype, mastertype> my;

     public:
      /// constructor
      DebugIncompleteIntPolicy(){/* empty */};


      void Get_Debug(Mortar::Element& sele, const Core::LinAlg::Matrix<my::SLAVENUMNODE, 1>& lmval,
          const double gapn_sl, const double gapn_ma, const double wgt, const double jac,
          const double* gpn, const double* mxigp) const;

      void Get_Deriv1st_Debug(Mortar::Element& sele,
          const Core::LinAlg::Matrix<my::SLAVENUMNODE, 1>& lmval,
          const Core::LinAlg::Matrix<my::SLAVENUMNODE, 1>& sval,
          const Core::LinAlg::Matrix<my::SLAVEDIM, my::SLAVENUMNODE>& sderiv,
          const Core::LinAlg::Matrix<3, 2>& stau, const Deriv1stMap& d_jac,
          const Deriv1stVecMap& dmxigp, const Deriv1stVecMap& d_gpn, const Deriv1stMap& d_gap_sl,
          const double gapn_sl, const double wgt, const double jac) const;

      void Get_Deriv2nd_Debug(Mortar::Element& sele,
          const Core::LinAlg::Matrix<my::SLAVENUMNODE, 1>& lmval,
          const Core::LinAlg::Matrix<my::SLAVENUMNODE, 1>& sval,
          const Core::LinAlg::Matrix<my::SLAVEDIM, my::SLAVENUMNODE>& sderiv,
          const Core::LinAlg::Matrix<3, 2>& stau, const Deriv1stMap& d_jac,
          const Deriv1stMap& d_gapn_sl, const Deriv2ndMap& dd_jac, const Deriv2ndVecMap& ddmxigp,
          const Deriv1stVecMap& d_n_unit, const Deriv2ndVecMap& dd_n_unit, const double gapn_sl,
          const double wgt, const double jac) const;

     private:
      void debug_deriv1st_w_gap(Mortar::Element& sele) const;

      void debug_deriv2nd_w_gap(Mortar::Element& sele) const;
    };

    /*--------------------------------------------------------------------------*/
    template <unsigned probdim, Core::FE::CellType slavetype, Core::FE::CellType mastertype>
    class DebugCompleteIntPolicy : public CompleteIntPolicy<probdim, slavetype, mastertype>
    {
      typedef CompleteIntPolicy<probdim, slavetype, mastertype> my;

     public:
      /// constructor
      DebugCompleteIntPolicy(){/* empty */};


      void Get_Debug(Mortar::Element& sele, const Core::LinAlg::Matrix<my::SLAVENUMNODE, 1>& lmval,
          const double gapn_sl, const double gapn_ma, const double wgt, const double jac,
          const double* gpn, const double* mxigp) const;

      void Get_Deriv1st_Debug(Mortar::Element& sele,
          const Core::LinAlg::Matrix<my::SLAVENUMNODE, 1>& lmval,
          const Core::LinAlg::Matrix<my::SLAVENUMNODE, 1>& sval,
          const Core::LinAlg::Matrix<my::SLAVEDIM, my::SLAVENUMNODE>& sderiv,
          const Core::LinAlg::Matrix<3, 2>& stau, const Deriv1stMap& d_jac,
          const Deriv1stVecMap& dmxigp, const Deriv1stVecMap& d_gpn, const Deriv1stMap& d_gap_sl,
          const double gapn_sl, const double wgt, const double jac) const;

      void Get_Deriv2nd_Debug(Mortar::Element& sele,
          const Core::LinAlg::Matrix<my::SLAVENUMNODE, 1>& lmval,
          const Core::LinAlg::Matrix<my::SLAVENUMNODE, 1>& sval,
          const Core::LinAlg::Matrix<my::SLAVEDIM, my::SLAVENUMNODE>& sderiv,
          const Core::LinAlg::Matrix<3, 2>& stau, const Deriv1stMap& d_jac,
          const Deriv1stMap& d_gapn_sl, const Deriv2ndMap& dd_jac, const Deriv2ndVecMap& ddmxigp,
          const Deriv1stVecMap& d_n_unit, const Deriv2ndVecMap& dd_n_unit, const double gapn_sl,
          const double wgt, const double jac) const;

     private:
      void debug_kappa(Mortar::Element& sele,
          const Core::LinAlg::Matrix<my::SLAVENUMNODE, 1>& lmval, const double jac,
          const double wgt) const;

      void debug_deriv1st_kappa(Mortar::Element& sele,
          const Core::LinAlg::Matrix<my::SLAVENUMNODE, 1>& lmval, const Deriv1stMap& djac,
          const double wgt) const;

      void debug_gpn(Mortar::Element& sele, const Core::LinAlg::Matrix<my::SLAVENUMNODE, 1>& lmval,
          const double* gpn) const;

      void debug_deriv1st_gpn(Mortar::Element& sele,
          const Core::LinAlg::Matrix<my::SLAVENUMNODE, 1>& lmval,
          const Deriv1stVecMap& d_gpn) const;

      void debug_w_gap(Mortar::Element& sele,
          const Core::LinAlg::Matrix<my::SLAVENUMNODE, 1>& lmval, const double gapn_sl,
          const double gapn_ma, const double wgt, const double jac) const;

      void debug_w_gap_sl(Mortar::Element& sele,
          const Core::LinAlg::Matrix<my::SLAVENUMNODE, 1>& lmval, const double gapn_sl,
          const double wgt, const double jac) const;

      void debug_deriv1st_gap_sl(Mortar::Element& sele,
          const Core::LinAlg::Matrix<my::SLAVENUMNODE, 1>& lmval,
          const Core::LinAlg::Matrix<my::SLAVENUMNODE, 1>& sval, const Deriv1stVecMap& d_n_unit,
          const Deriv1stMap& d_gapn_sl, const Deriv1stMap& d_jac, const double gapn_sl,
          const double wgt, const double jac, Deriv1stMap& d_gapn_sl_complete) const;

      void debug_deriv1st_w_gap_sl(Mortar::Element& sele,
          const Core::LinAlg::Matrix<my::SLAVENUMNODE, 1>& lmval,
          const Core::LinAlg::Matrix<my::SLAVENUMNODE, 1>& sval, const Deriv1stVecMap& d_gpn,
          const Deriv1stMap& d_gapn_sl, const Deriv1stMap& d_jac, const double gapn_sl,
          const double wgt, const double jac) const;

      void debug_deriv2nd_w_gap_sl(Mortar::Element& sele,
          const Core::LinAlg::Matrix<my::SLAVENUMNODE, 1>& lmval,
          const Core::LinAlg::Matrix<my::SLAVENUMNODE, 1>& sval, const Deriv1stVecMap& d_n_unit,
          const Deriv2ndVecMap& dd_n_unit, const Deriv1stMap& d_jac,
          const Deriv1stMap& d_gapn_sl_complete, const Deriv2ndMap& dd_jac, const double gapn_sl,
          const double wgt, const double jac) const;

      void debug_deriv1st_w_gap(Mortar::Element& sele) const;

      void debug_deriv2nd_w_gap(Mortar::Element& sele) const;

      void debug_m_xi(Mortar::Element& sele, const Core::LinAlg::Matrix<my::SLAVENUMNODE, 1>& lmval,
          const double* mxigp) const;

      void debug_deriv1st_m_xi(Mortar::Element& sele,
          const Core::LinAlg::Matrix<my::SLAVENUMNODE, 1>& lmval,
          const Deriv1stVecMap& dmxigp) const;

      void debug_deriv2nd_m_xi(Mortar::Element& sele,
          const Core::LinAlg::Matrix<my::SLAVENUMNODE, 1>& lmval,
          const Deriv2ndVecMap& ddmxigp) const;

      void debug_deriv1st_smooth_unit_normal(Mortar::Element& sele) const;

      void debug_deriv2nd_smooth_unit_normal(Mortar::Element& sele) const;

      void debug_deriv1st_jac(Mortar::Element& sele,
          const Core::LinAlg::Matrix<my::SLAVENUMNODE, 1>& lmval, const Deriv1stMap& djac) const;

      void debug_deriv2nd_jac(Mortar::Element& sele,
          const Core::LinAlg::Matrix<my::SLAVENUMNODE, 1>& lmval, const Deriv2ndMap& ddjac) const;

      void debug_deriv1st_jacobian(Mortar::Element& sele,
          const Core::LinAlg::Matrix<my::SLAVENUMNODE, 1>& lmval,
          const Core::LinAlg::Matrix<my::SLAVEDIM, my::SLAVENUMNODE>& sderiv,
          const Core::LinAlg::Matrix<3, 2>& stau) const;

      void debug_deriv2nd_jacobian(Mortar::Element& sele,
          const Core::LinAlg::Matrix<my::SLAVENUMNODE, 1>& lmval,
          const Core::LinAlg::Matrix<my::SLAVEDIM, my::SLAVENUMNODE>& sderiv,
          const Core::LinAlg::Matrix<3, 2>& stau) const;

      void debug_deriv1st_non_unit_normal(Mortar::Element& sele,
          const Core::LinAlg::Matrix<my::SLAVENUMNODE, 1>& lmval,
          const Core::LinAlg::Matrix<my::SLAVEDIM, my::SLAVENUMNODE>& sderiv,
          const Core::LinAlg::Matrix<3, 2>& stau) const;

      void debug_deriv2nd_non_unit_normal(Mortar::Element& sele,
          const Core::LinAlg::Matrix<my::SLAVENUMNODE, 1>& lmval,
          const Core::LinAlg::Matrix<my::SLAVEDIM, my::SLAVENUMNODE>& sderiv) const;

     private:
      const Deriv1stMap& get_nodal_deriv1st(NodeDataContainer& data) const;

      const Deriv2ndMap& get_nodal_deriv2nd(NodeDataContainer& data) const;
    };
  }  // namespace Aug
}  // namespace CONTACT


FOUR_C_NAMESPACE_CLOSE

#endif
