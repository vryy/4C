/*----------------------------------------------------------------------------*/
/*! \file
\brief (augmented) contact integration policies

\level 3

*/
/*----------------------------------------------------------------------------*/

#ifndef FOUR_C_CONTACT_AUG_INTEGRATOR_POLICY_HPP
#define FOUR_C_CONTACT_AUG_INTEGRATOR_POLICY_HPP

#include "baci_config.hpp"

#include "baci_contact_aug_contact_integrator_utils.hpp"
#include "baci_contact_aug_timemonitor.hpp"
#include "baci_lib_node.hpp"

BACI_NAMESPACE_OPEN

namespace CONTACT
{
  namespace AUG
  {
    // forward declarations
    class NodeDataContainer;

    /*--------------------------------------------------------------------------*/
    template <unsigned probdim, CORE::FE::CellType slavetype>
    class BaseSlaveIntPolicy
    {
     protected:
      static constexpr unsigned SLAVEDIM = CORE::FE::dim<slavetype>;
      static constexpr unsigned SLAVENUMNODE = CORE::FE::num_nodes<slavetype>;

     public:
      /// constructor
      BaseSlaveIntPolicy(){/* empty */};


      double UnitSlaveElementNormal(const MORTAR::Element& sele,
          const CORE::LINALG::Matrix<3, 2>& tau,
          CORE::LINALG::Matrix<probdim, 1>& unit_normal) const;

      /// @name First derivatives
      /// @{

      void Deriv1st_NonUnitSlaveElementNormal(const MORTAR::Element& sele,
          const CORE::LINALG::Matrix<probdim, SLAVENUMNODE, int>& nodal_dofs,
          const CORE::LINALG::Matrix<SLAVEDIM, SLAVENUMNODE>& deriv,
          const CORE::LINALG::Matrix<3, 2>& tau, Deriv1stVecMap& d_non_unit_normal) const;

      void Deriv1st_UnitSlaveElementNormal(const CORE::LINALG::Matrix<probdim, 1>& unit_normal,
          const double length_n_inv, const Deriv1stVecMap& d_non_unit_normal,
          Deriv1stVecMap& d_unit_normal, const bool reset = true) const;

      void Deriv1st_Jacobian(const CORE::LINALG::Matrix<probdim, 1>& unit_normal,
          const Deriv1stVecMap& d_non_unit_normal, Deriv1stMap& d_jac) const;

      /// @}

      /// @name Second derivatives
      /// @{

      void Deriv2nd_Jacobian(const Deriv1stVecMap& d_unit_normal,
          const Deriv1stVecMap& d_non_unit_normal,
          const CORE::LINALG::Matrix<probdim, 1>& unit_normal,
          const Deriv2ndVecMap& dd_non_unit_normal, Deriv2ndMap& dd_jac) const;

      void Deriv2nd_NonUnitSlaveElementNormal(const MORTAR::Element& sele,
          const CORE::LINALG::Matrix<probdim, SLAVENUMNODE, int>& nodal_dofs,
          const CORE::LINALG::Matrix<SLAVEDIM, SLAVENUMNODE>& deriv,
          Deriv2ndVecMap& dd_non_unit_normal) const;

      void Deriv2nd_UnitSlaveElementNormal(const CORE::LINALG::Matrix<probdim, 1>& unit_normal,
          const double length_n_inv, const Deriv1stVecMap& d_non_unit_normal,
          const Deriv1stVecMap& d_unit_normal, const Deriv2ndVecMap& dd_non_unit_normal,
          Deriv2ndVecMap& dd_unit_normal) const;

      /// @}

      void AveragedNormalAtXi(MORTAR::Element& sele,
          const CORE::LINALG::Matrix<SLAVENUMNODE, 1>& sval,
          CORE::LINALG::Matrix<probdim, 1>& snormal) const;

      /// complete data in nodal data container after successful integration
      void CompleteNodeData(MORTAR::Element& sele) const;

     private:
      void InnerProductOfVectorAndDeriv1stVector(const CORE::LINALG::Matrix<probdim, 1>& vec,
          const Deriv1stVecMap& d_vec, Deriv1stMap& dvec_vec) const;

      void ProjectionIntoTangentialPlain(const CORE::LINALG::Matrix<probdim, 1>& unit_normal,
          CORE::LINALG::Matrix<probdim, probdim>& tproj_mat) const;

     protected:
      /// current Gauss point ID
      int gp_id_ = -1;
      mutable TimeMonitor<CONTACT::AUG::TimeID> timer_ = TimeMonitor<CONTACT::AUG::TimeID>();
    };

    /*--------------------------------------------------------------------------*/
    template <unsigned probdim, CORE::FE::CellType slavetype, CORE::FE::CellType mastertype>
    class BaseIntPolicy : public BaseSlaveIntPolicy<probdim, slavetype>
    {
      typedef BaseSlaveIntPolicy<probdim, slavetype> my;

     protected:
      static constexpr unsigned MASTERDIM = CORE::FE::dim<mastertype>;
      static constexpr unsigned MASTERNUMNODE = CORE::FE::num_nodes<mastertype>;

     public:
      /// constructor
      BaseIntPolicy(){/* empty */};


      void LMatrixInverse(const CORE::LINALG::Matrix<3, 2>& mtau,
          const CORE::LINALG::Matrix<probdim, 1>& snormal,
          CORE::LINALG::Matrix<probdim, probdim>& lmat_inv) const;

      void Deriv1st_MXiGP(const CORE::LINALG::Matrix<probdim, probdim>& lmat_inv,
          MORTAR::Element& sele, MORTAR::Element& mele,
          const CORE::LINALG::Matrix<my::SLAVENUMNODE, 1>& sval,
          const CORE::LINALG::Matrix<MASTERNUMNODE, 1>& mval, const double alpha,
          Deriv1stVecMap& d_mxi, Deriv1stMap& d_alpha) const;

      void Add_Deriv2nd_MaDispl(const CORE::LINALG::Matrix<probdim, probdim>& lmat_inv,
          const CORE::LINALG::Matrix<probdim, MASTERNUMNODE, int>& mnodal_dofs,
          const CORE::LINALG::Matrix<MASTERDIM, MASTERNUMNODE>& mderiv,
          const Deriv1stVecMap& d_mxigp, Deriv2ndVecMap& dd_mxigp) const;

      void Add_Deriv1st_MaMetric(const CORE::LINALG::Matrix<probdim, probdim>& lmat_inv,
          const CORE::LINALG::Matrix<probdim, MASTERNUMNODE, int>& mnodal_dofs,
          const CORE::LINALG::Matrix<3, 2>& mtau,
          const CORE::LINALG::Matrix<MASTERDIM, MASTERNUMNODE>& mderiv,
          const CORE::LINALG::Matrix<3, MASTERNUMNODE>& mderiv2nd,
          const CORE::LINALG::Matrix<3, MASTERNUMNODE>& mcoord, const Deriv1stVecMap& d_mxigp,
          Deriv2ndVecMap& dd_mxigp) const;

      void Add_Deriv1st_Alpha_Deriv1st_Normal(
          const CORE::LINALG::Matrix<probdim, probdim>& lmat_inv, const Deriv1stMap& d_alpha,
          const MORTAR::Element& sele, const CORE::LINALG::Matrix<my::SLAVENUMNODE, 1>& sval,
          Deriv2ndVecMap& dd_mxigp) const;

      void Add_Alpha_Deriv2nd_Normal(const CORE::LINALG::Matrix<probdim, probdim>& lmat_inv,
          const double alpha, const MORTAR::Element& sele,
          const CORE::LINALG::Matrix<my::SLAVENUMNODE, 1>& sval, Deriv2ndVecMap& dd_mxigp) const;

      void Deriv1st_GapN_Sl(const DRT::Node* const* snodes,
          const CORE::LINALG::Matrix<my::SLAVENUMNODE, 1>& sval, const double* gpn,
          Deriv1stMap& deriv_gapn_sl) const;

      void Deriv1st_GapN_Ma(const DRT::Node* const* mnodes,
          const CORE::LINALG::Matrix<MASTERNUMNODE, 1>& mval, const double* gpn,
          const CORE::LINALG::Matrix<3, 2>& mtau, const Deriv1stVecMap& d_mxigp,
          Deriv1stMap& deriv_gapn_ma) const;

      void Add_Deriv1stGapNContributions(DRT::Node* const* snodes, const double scale,
          const CORE::LINALG::Matrix<my::SLAVENUMNODE, 1>& lmval, const Deriv1stMap& d_gapn_sl,
          const Deriv1stMap& d_gapn_ma) const;

      void Add_Deriv1stJacobianContributions(DRT::Node* const* snodes, const double wgt,
          const CORE::LINALG::Matrix<my::SLAVENUMNODE, 1>& lmval, const double gapn_sl,
          const double gapn_ma, const Deriv1stMap& d_jac) const;

      void Add_Var_GapN_Lin_Jac(MORTAR::Element& sele,
          const CORE::LINALG::Matrix<my::SLAVENUMNODE, 1>& lmval, const double wgt,
          const Deriv1stMap& d_gapn_sl, const Deriv1stMap& d_gapn_ma,
          const Deriv1stMap& d_jac) const;

      void Add_Var_Jac_Lin_GapN(MORTAR::Element& sele,
          const CORE::LINALG::Matrix<my::SLAVENUMNODE, 1>& lmval, const double wgt,
          const Deriv1stMap& d_gapn_sl, const Deriv1stMap& d_gapn_ma,
          const Deriv1stMap& d_jac) const;

      void Add_GapN_Deriv2nd_Jac(MORTAR::Element& sele,
          const CORE::LINALG::Matrix<my::SLAVENUMNODE, 1>& lmval, const double wgt,
          const double gapn_sl, const double gapn_ma, const Deriv2ndMap& dd_jac) const;

      /// @name empty functions, if no debug policy is in use
      /// @{

      void Get_Debug(...) const { return; }

      void Get_Deriv1st_Debug(...) const { return; }

      void Get_Deriv2nd_Debug(...) const { return; }

      /// @}
    };

    /*--------------------------------------------------------------------------*/
    template <unsigned probdim, CORE::FE::CellType slavetype, CORE::FE::CellType mastertype>
    class IncompleteIntPolicy : public BaseIntPolicy<probdim, slavetype, mastertype>
    {
      typedef BaseIntPolicy<probdim, slavetype, mastertype> my;

     public:
      /// constructor
      IncompleteIntPolicy(){/* empty */};


      void Get_Deriv2nd_Jacobian(const MORTAR::Element& sele,
          const CORE::LINALG::Matrix<probdim, my::SLAVENUMNODE, int>& nodal_dofs,
          const CORE::LINALG::Matrix<my::SLAVEDIM, my::SLAVENUMNODE>& deriv,
          const CORE::LINALG::Matrix<probdim, 1>& unit_normal, const double length_n_inv,
          const Deriv1stVecMap& d_non_unit_normal, Deriv2ndMap& dd_jac) const;

      void Get_Deriv2nd_MXiGP(const CORE::LINALG::Matrix<probdim, probdim>& lmat_inv,
          MORTAR::Element& sele, MORTAR::Element& mele,
          const CORE::LINALG::Matrix<my::SLAVENUMNODE, 1>& sval,
          const CORE::LINALG::Matrix<my::MASTERNUMNODE, 1>& mval,
          const CORE::LINALG::Matrix<my::MASTERDIM, my::MASTERNUMNODE>& mderiv,
          const CORE::LINALG::Matrix<3, my::MASTERNUMNODE>& mderiv2,
          const CORE::LINALG::Matrix<3, 2>& mtau, const double* mxi, const double alpha,
          const Deriv1stVecMap& d_mxigp, const Deriv1stMap& d_alpha, Deriv2ndVecMap& dd_mxigp) const
      {
        // do nothing
        CORE::GEN::reset(2, 0, dd_mxigp);
      }

      void Get_Deriv1st_GapN(MORTAR::Element& sele, MORTAR::Element& mele,
          const CORE::LINALG::Matrix<my::SLAVENUMNODE, 1>& sval,
          const CORE::LINALG::Matrix<my::MASTERNUMNODE, 1>& mval, const double* gpn,
          const CORE::LINALG::Matrix<3, 2>& mtau, const Deriv1stVecMap& d_mxigp,
          Deriv1stMap& deriv_gapn_sl, Deriv1stMap& deriv_gapn_ma) const;

      void Get_Deriv1st_WGap(MORTAR::Element& sele,
          const CORE::LINALG::Matrix<my::SLAVENUMNODE, 1>& lmval, const double gapn_sl,
          const double gapn_ma, const double wgt, const double jac, const Deriv1stMap& d_jac,
          const Deriv1stMap& d_gapn_sl, const Deriv1stMap& d_gapn_ma) const;

      void Get_Deriv1st_WGap_Complete(const int linsize, MORTAR::Element& sele,
          MORTAR::Element& mele, const CORE::LINALG::Matrix<my::SLAVENUMNODE, 1>& sval,
          const CORE::LINALG::Matrix<my::MASTERNUMNODE, 1>& mval,
          const CORE::LINALG::Matrix<my::SLAVENUMNODE, 1>& lmval, const double* gpn,
          const CORE::LINALG::Matrix<3, 2>& mtau, const Deriv1stVecMap& d_mxigp,
          const double gapn_sl, const double gapn_ma, const double wgt, const double jac,
          const Deriv1stMap& d_jac) const;

      void Get_Deriv2nd_WGap(MORTAR::Element& sele, MORTAR::Element& mele,
          const CORE::LINALG::Matrix<my::SLAVENUMNODE, 1>& sval,
          const CORE::LINALG::Matrix<my::MASTERNUMNODE, 1>& mval,
          const CORE::LINALG::Matrix<my::SLAVENUMNODE, 1>& lmval,
          const CORE::LINALG::Matrix<my::MASTERDIM, my::MASTERNUMNODE>& mderiv,
          const CORE::LINALG::Matrix<3, my::MASTERNUMNODE>& mderiv2nd,
          const CORE::LINALG::Matrix<3, 2>& mtau, const double* gpn, const double wgt,
          const double gapn_sl, const double gapn_ma, const double jac, const Deriv1stMap& d_jac,
          const Deriv2ndMap& dd_jac, const Deriv1stVecMap& d_mxigp, const Deriv2ndVecMap& dd_mxigp,
          const Deriv1stVecMap& d_n_unit, const Deriv2ndVecMap& dd_n_unit,
          const Deriv1stMap& d_gapn_sl, const Deriv1stMap& d_gapn_ma) const;

      void Get_Deriv1st_WGapN_Error(const MORTAR::Element& sele,
          const std::vector<unsigned>& active_nlids,
          const CORE::LINALG::Matrix<my::SLAVENUMNODE, 1>& lmval, const double* gpn,
          const double gapn_sl, const double gapn_ma, const double wgt, const double jacslave,
          const Deriv1stMap& d_jac, const CORE::LINALG::Matrix<3, 2>& mtau,
          const Deriv1stVecMap& d_mxigp, Deriv1stMap& d_gapn_ma,
          std::unordered_map<int, Deriv1stMap>& error_ma,
          std::unordered_map<int, Deriv1stMap>& error_jac) const;

     protected:
      void Add_Deriv1st_GapN_Deriv1st_Jac(MORTAR::Element& sele,
          const CORE::LINALG::Matrix<my::SLAVENUMNODE, 1>& lmval, const double wgt,
          const Deriv1stMap& d_gapn_sl, const Deriv1stMap& d_gapn_ma,
          const Deriv1stMap& d_jac) const;

      void Add_Jac_Deriv2nd_GapN(MORTAR::Element& sele, MORTAR::Element& mele,
          const CORE::LINALG::Matrix<my::SLAVENUMNODE, 1>& sval,
          const CORE::LINALG::Matrix<my::MASTERNUMNODE, 1>& mval,
          const CORE::LINALG::Matrix<my::SLAVENUMNODE, 1>& lmval,
          const CORE::LINALG::Matrix<my::MASTERDIM, my::MASTERNUMNODE>& mderiv,
          const CORE::LINALG::Matrix<3, 2>& mtau, const double* gpn, const double wgt,
          const double jac, const Deriv1stVecMap& d_mxigp, const Deriv1stVecMap& d_n_unit,
          const Deriv2ndVecMap& dd_n_unit) const;
    };

    /*--------------------------------------------------------------------------*/
    template <unsigned probdim, CORE::FE::CellType slavetype, CORE::FE::CellType mastertype>
    class CompleteIntPolicy : public BaseIntPolicy<probdim, slavetype, mastertype>
    {
      typedef BaseIntPolicy<probdim, slavetype, mastertype> my;

     public:
      /// constructor
      CompleteIntPolicy(){/* empty */};


      void Get_Deriv2nd_Jacobian(const MORTAR::Element& sele,
          const CORE::LINALG::Matrix<probdim, my::SLAVENUMNODE, int>& nodal_dofs,
          const CORE::LINALG::Matrix<my::SLAVEDIM, my::SLAVENUMNODE>& deriv,
          const CORE::LINALG::Matrix<probdim, 1>& unit_normal, const double length_n_inv,
          const Deriv1stVecMap& d_non_unit_normal, Deriv2ndMap& dd_jac) const;

      void Get_Deriv2nd_MXiGP(const CORE::LINALG::Matrix<probdim, probdim>& lmat_inv,
          MORTAR::Element& sele, MORTAR::Element& mele,
          const CORE::LINALG::Matrix<my::SLAVENUMNODE, 1>& sval,
          const CORE::LINALG::Matrix<my::MASTERNUMNODE, 1>& mval,
          const CORE::LINALG::Matrix<my::MASTERDIM, my::MASTERNUMNODE>& mderiv,
          const CORE::LINALG::Matrix<3, my::MASTERNUMNODE>& mderiv2,
          const CORE::LINALG::Matrix<3, 2>& mtau, const double* mxi, const double alpha,
          const Deriv1stVecMap& d_mxigp, const Deriv1stMap& d_alpha,
          Deriv2ndVecMap& dd_mxigp) const;

      void Get_Deriv1st_GapN(MORTAR::Element& sele, MORTAR::Element& mele,
          const CORE::LINALG::Matrix<my::SLAVENUMNODE, 1>& sval,
          const CORE::LINALG::Matrix<my::MASTERNUMNODE, 1>& mval, const double* gpn,
          const CORE::LINALG::Matrix<3, 2>& mtau, const Deriv1stVecMap& d_mxigp,
          Deriv1stMap& deriv_gapn_sl, Deriv1stMap& deriv_gapn_ma) const;

      void Get_Deriv1st_WGap(MORTAR::Element& sele,
          const CORE::LINALG::Matrix<my::SLAVENUMNODE, 1>& lmval, const double gapn_sl,
          const double gapn_ma, const double wgt, const double jac, const Deriv1stMap& d_jac,
          const Deriv1stMap& d_gapn_sl, const Deriv1stMap& d_gapn_ma) const;

      void Get_Deriv1st_WGap_Complete(const int linsize, MORTAR::Element& sele,
          MORTAR::Element& mele, const CORE::LINALG::Matrix<my::SLAVENUMNODE, 1>& sval,
          const CORE::LINALG::Matrix<my::MASTERNUMNODE, 1>& mval,
          const CORE::LINALG::Matrix<my::SLAVENUMNODE, 1>& lmval, const double* gpn,
          const CORE::LINALG::Matrix<3, 2>& mtau, const Deriv1stVecMap& d_mxigp,
          const double gapn_sl, const double gapn_ma, const double wgt, const double jac,
          const Deriv1stMap& d_jac) const;

      void Get_Deriv2nd_WGap(MORTAR::Element& sele, MORTAR::Element& mele,
          const CORE::LINALG::Matrix<my::SLAVENUMNODE, 1>& sval,
          const CORE::LINALG::Matrix<my::MASTERNUMNODE, 1>& mval,
          const CORE::LINALG::Matrix<my::SLAVENUMNODE, 1>& lmval,
          const CORE::LINALG::Matrix<my::MASTERDIM, my::MASTERNUMNODE>& mderiv,
          const CORE::LINALG::Matrix<3, my::MASTERNUMNODE>& mderiv2nd,
          const CORE::LINALG::Matrix<3, 2>& mtau, const double* gpn, const double wgt,
          const double gapn_sl, const double gapn_ma, const double jac, const Deriv1stMap& d_jac,
          const Deriv2ndMap& dd_jac, const Deriv1stVecMap& d_mxigp, const Deriv2ndVecMap& dd_mxigp,
          const Deriv1stVecMap& d_n_unit, const Deriv2ndVecMap& dd_n_unit,
          const Deriv1stMap& d_gapn_sl, const Deriv1stMap& d_gapn_ma) const;

      void Get_Deriv1st_WGapN_Error(const MORTAR::Element& sele,
          const std::vector<unsigned>& active_nlids,
          const CORE::LINALG::Matrix<my::SLAVENUMNODE, 1>& lmval, const double* gpn,
          const double gapn_sl, const double gapn_ma, const double wgt, const double jacslave,
          const Deriv1stMap& d_jac, const CORE::LINALG::Matrix<3, 2>& mtau,
          const Deriv1stVecMap& d_mxigp, Deriv1stMap& d_gapn_ma,
          std::unordered_map<int, Deriv1stMap>& error_ma,
          std::unordered_map<int, Deriv1stMap>& error_jac) const;

     protected:
      void Add_Deriv1st_GapN_Deriv1st_Jac(MORTAR::Element& sele,
          const CORE::LINALG::Matrix<my::SLAVENUMNODE, 1>& lmval, const double wgt,
          const Deriv1stMap& d_gapn_sl, const Deriv1stMap& d_gapn_ma,
          const Deriv1stMap& d_jac) const;

      void Add_Jac_Deriv2nd_GapN(MORTAR::Element& sele, MORTAR::Element& mele,
          const CORE::LINALG::Matrix<my::SLAVENUMNODE, 1>& sval,
          const CORE::LINALG::Matrix<my::MASTERNUMNODE, 1>& mval,
          const CORE::LINALG::Matrix<my::SLAVENUMNODE, 1>& lmval,
          const CORE::LINALG::Matrix<my::MASTERDIM, my::MASTERNUMNODE>& mderiv,
          const CORE::LINALG::Matrix<3, my::MASTERNUMNODE>& mderiv2nd,
          const CORE::LINALG::Matrix<3, 2>& mtau, const double* gpn, const double wgt,
          const double jac, const Deriv1stVecMap& d_mxigp, const Deriv2ndVecMap& dd_mxigp,
          const Deriv1stVecMap& d_n_unit, const Deriv2ndVecMap& dd_n_unit) const;
    };

    /*--------------------------------------------------------------------------*/
    template <unsigned probdim, CORE::FE::CellType slavetype, CORE::FE::CellType mastertype>
    class DebugIncompleteIntPolicy : public IncompleteIntPolicy<probdim, slavetype, mastertype>
    {
      typedef IncompleteIntPolicy<probdim, slavetype, mastertype> my;

     public:
      /// constructor
      DebugIncompleteIntPolicy(){/* empty */};


      void Get_Debug(MORTAR::Element& sele, const CORE::LINALG::Matrix<my::SLAVENUMNODE, 1>& lmval,
          const double gapn_sl, const double gapn_ma, const double wgt, const double jac,
          const double* gpn, const double* mxigp) const;

      void Get_Deriv1st_Debug(MORTAR::Element& sele,
          const CORE::LINALG::Matrix<my::SLAVENUMNODE, 1>& lmval,
          const CORE::LINALG::Matrix<my::SLAVENUMNODE, 1>& sval,
          const CORE::LINALG::Matrix<my::SLAVEDIM, my::SLAVENUMNODE>& sderiv,
          const CORE::LINALG::Matrix<3, 2>& stau, const Deriv1stMap& d_jac,
          const Deriv1stVecMap& dmxigp, const Deriv1stVecMap& d_gpn, const Deriv1stMap& d_gap_sl,
          const double gapn_sl, const double wgt, const double jac) const;

      void Get_Deriv2nd_Debug(MORTAR::Element& sele,
          const CORE::LINALG::Matrix<my::SLAVENUMNODE, 1>& lmval,
          const CORE::LINALG::Matrix<my::SLAVENUMNODE, 1>& sval,
          const CORE::LINALG::Matrix<my::SLAVEDIM, my::SLAVENUMNODE>& sderiv,
          const CORE::LINALG::Matrix<3, 2>& stau, const Deriv1stMap& d_jac,
          const Deriv1stMap& d_gapn_sl, const Deriv2ndMap& dd_jac, const Deriv2ndVecMap& ddmxigp,
          const Deriv1stVecMap& d_n_unit, const Deriv2ndVecMap& dd_n_unit, const double gapn_sl,
          const double wgt, const double jac) const;

     private:
      void Debug_Deriv1st_WGap(MORTAR::Element& sele) const;

      void Debug_Deriv2nd_WGap(MORTAR::Element& sele) const;
    };

    /*--------------------------------------------------------------------------*/
    template <unsigned probdim, CORE::FE::CellType slavetype, CORE::FE::CellType mastertype>
    class DebugCompleteIntPolicy : public CompleteIntPolicy<probdim, slavetype, mastertype>
    {
      typedef CompleteIntPolicy<probdim, slavetype, mastertype> my;

     public:
      /// constructor
      DebugCompleteIntPolicy(){/* empty */};


      void Get_Debug(MORTAR::Element& sele, const CORE::LINALG::Matrix<my::SLAVENUMNODE, 1>& lmval,
          const double gapn_sl, const double gapn_ma, const double wgt, const double jac,
          const double* gpn, const double* mxigp) const;

      void Get_Deriv1st_Debug(MORTAR::Element& sele,
          const CORE::LINALG::Matrix<my::SLAVENUMNODE, 1>& lmval,
          const CORE::LINALG::Matrix<my::SLAVENUMNODE, 1>& sval,
          const CORE::LINALG::Matrix<my::SLAVEDIM, my::SLAVENUMNODE>& sderiv,
          const CORE::LINALG::Matrix<3, 2>& stau, const Deriv1stMap& d_jac,
          const Deriv1stVecMap& dmxigp, const Deriv1stVecMap& d_gpn, const Deriv1stMap& d_gap_sl,
          const double gapn_sl, const double wgt, const double jac) const;

      void Get_Deriv2nd_Debug(MORTAR::Element& sele,
          const CORE::LINALG::Matrix<my::SLAVENUMNODE, 1>& lmval,
          const CORE::LINALG::Matrix<my::SLAVENUMNODE, 1>& sval,
          const CORE::LINALG::Matrix<my::SLAVEDIM, my::SLAVENUMNODE>& sderiv,
          const CORE::LINALG::Matrix<3, 2>& stau, const Deriv1stMap& d_jac,
          const Deriv1stMap& d_gapn_sl, const Deriv2ndMap& dd_jac, const Deriv2ndVecMap& ddmxigp,
          const Deriv1stVecMap& d_n_unit, const Deriv2ndVecMap& dd_n_unit, const double gapn_sl,
          const double wgt, const double jac) const;

     private:
      void Debug_Kappa(MORTAR::Element& sele,
          const CORE::LINALG::Matrix<my::SLAVENUMNODE, 1>& lmval, const double jac,
          const double wgt) const;

      void Debug_Deriv1st_Kappa(MORTAR::Element& sele,
          const CORE::LINALG::Matrix<my::SLAVENUMNODE, 1>& lmval, const Deriv1stMap& djac,
          const double wgt) const;

      void Debug_GPN(MORTAR::Element& sele, const CORE::LINALG::Matrix<my::SLAVENUMNODE, 1>& lmval,
          const double* gpn) const;

      void Debug_Deriv1st_GPN(MORTAR::Element& sele,
          const CORE::LINALG::Matrix<my::SLAVENUMNODE, 1>& lmval,
          const Deriv1stVecMap& d_gpn) const;

      void Debug_WGap(MORTAR::Element& sele, const CORE::LINALG::Matrix<my::SLAVENUMNODE, 1>& lmval,
          const double gapn_sl, const double gapn_ma, const double wgt, const double jac) const;

      void Debug_WGap_Sl(MORTAR::Element& sele,
          const CORE::LINALG::Matrix<my::SLAVENUMNODE, 1>& lmval, const double gapn_sl,
          const double wgt, const double jac) const;

      void Debug_Deriv1st_Gap_Sl(MORTAR::Element& sele,
          const CORE::LINALG::Matrix<my::SLAVENUMNODE, 1>& lmval,
          const CORE::LINALG::Matrix<my::SLAVENUMNODE, 1>& sval, const Deriv1stVecMap& d_n_unit,
          const Deriv1stMap& d_gapn_sl, const Deriv1stMap& d_jac, const double gapn_sl,
          const double wgt, const double jac, Deriv1stMap& d_gapn_sl_complete) const;

      void Debug_Deriv1st_WGap_Sl(MORTAR::Element& sele,
          const CORE::LINALG::Matrix<my::SLAVENUMNODE, 1>& lmval,
          const CORE::LINALG::Matrix<my::SLAVENUMNODE, 1>& sval, const Deriv1stVecMap& d_gpn,
          const Deriv1stMap& d_gapn_sl, const Deriv1stMap& d_jac, const double gapn_sl,
          const double wgt, const double jac) const;

      void Debug_Deriv2nd_WGap_Sl(MORTAR::Element& sele,
          const CORE::LINALG::Matrix<my::SLAVENUMNODE, 1>& lmval,
          const CORE::LINALG::Matrix<my::SLAVENUMNODE, 1>& sval, const Deriv1stVecMap& d_n_unit,
          const Deriv2ndVecMap& dd_n_unit, const Deriv1stMap& d_jac,
          const Deriv1stMap& d_gapn_sl_complete, const Deriv2ndMap& dd_jac, const double gapn_sl,
          const double wgt, const double jac) const;

      void Debug_Deriv1st_WGap(MORTAR::Element& sele) const;

      void Debug_Deriv2nd_WGap(MORTAR::Element& sele) const;

      void Debug_MXi(MORTAR::Element& sele, const CORE::LINALG::Matrix<my::SLAVENUMNODE, 1>& lmval,
          const double* mxigp) const;

      void Debug_Deriv1st_MXi(MORTAR::Element& sele,
          const CORE::LINALG::Matrix<my::SLAVENUMNODE, 1>& lmval,
          const Deriv1stVecMap& dmxigp) const;

      void Debug_Deriv2nd_MXi(MORTAR::Element& sele,
          const CORE::LINALG::Matrix<my::SLAVENUMNODE, 1>& lmval,
          const Deriv2ndVecMap& ddmxigp) const;

      void Debug_Deriv1st_Smooth_Unit_Normal(MORTAR::Element& sele) const;

      void Debug_Deriv2nd_Smooth_Unit_Normal(MORTAR::Element& sele) const;

      void Debug_Deriv1st_Jac(MORTAR::Element& sele,
          const CORE::LINALG::Matrix<my::SLAVENUMNODE, 1>& lmval, const Deriv1stMap& djac) const;

      void Debug_Deriv2nd_Jac(MORTAR::Element& sele,
          const CORE::LINALG::Matrix<my::SLAVENUMNODE, 1>& lmval, const Deriv2ndMap& ddjac) const;

      void Debug_Deriv1st_Jacobian(MORTAR::Element& sele,
          const CORE::LINALG::Matrix<my::SLAVENUMNODE, 1>& lmval,
          const CORE::LINALG::Matrix<my::SLAVEDIM, my::SLAVENUMNODE>& sderiv,
          const CORE::LINALG::Matrix<3, 2>& stau) const;

      void Debug_Deriv2nd_Jacobian(MORTAR::Element& sele,
          const CORE::LINALG::Matrix<my::SLAVENUMNODE, 1>& lmval,
          const CORE::LINALG::Matrix<my::SLAVEDIM, my::SLAVENUMNODE>& sderiv,
          const CORE::LINALG::Matrix<3, 2>& stau) const;

      void Debug_Deriv1st_Non_Unit_Normal(MORTAR::Element& sele,
          const CORE::LINALG::Matrix<my::SLAVENUMNODE, 1>& lmval,
          const CORE::LINALG::Matrix<my::SLAVEDIM, my::SLAVENUMNODE>& sderiv,
          const CORE::LINALG::Matrix<3, 2>& stau) const;

      void Debug_Deriv2nd_Non_Unit_Normal(MORTAR::Element& sele,
          const CORE::LINALG::Matrix<my::SLAVENUMNODE, 1>& lmval,
          const CORE::LINALG::Matrix<my::SLAVEDIM, my::SLAVENUMNODE>& sderiv) const;

     private:
      const Deriv1stMap& GetNodalDeriv1st(NodeDataContainer& data) const;

      const Deriv2ndMap& GetNodalDeriv2nd(NodeDataContainer& data) const;
    };
  }  // namespace AUG
}  // namespace CONTACT


BACI_NAMESPACE_CLOSE

#endif  // CONTACT_AUG_INTEGRATOR_POLICY_H
