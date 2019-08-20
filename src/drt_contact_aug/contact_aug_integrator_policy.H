/*----------------------------------------------------------------------------*/
/*! \file
\brief (augmented) contact integration policies

\level 3

\maintainer Matthias Mayr
*/
/*----------------------------------------------------------------------------*/

#ifndef SRC_DRT_CONTACT_AUG_CONTACT_AUG_INTEGRATOR_POLICY_H_
#define SRC_DRT_CONTACT_AUG_CONTACT_AUG_INTEGRATOR_POLICY_H_

//#define SWITCH_2ND_DERIV_NORMAL_OFF

#include "contact_integrator_utils.H"
#include "contact_aug_timemonitor.H"

namespace CONTACT
{
  namespace AUG
  {
    // forward declarations
    class NodeDataContainer;

    /*--------------------------------------------------------------------------*/
    template <unsigned probdim, DRT::Element::DiscretizationType slavetype>
    class BaseSlaveIntPolicy
    {
     protected:
      static const unsigned SLAVEDIM = DRT::UTILS::DisTypeToDim<slavetype>::dim;
      static const unsigned SLAVENUMNODE =
          DRT::UTILS::DisTypeToNumNodePerEle<slavetype>::numNodePerElement;

     public:
      /// constructor
      BaseSlaveIntPolicy(){/* empty */};

      /// destructor
      ~BaseSlaveIntPolicy(){/* empty */};

      double UnitSlaveElementNormal(const MORTAR::MortarElement& sele,
          const LINALG::Matrix<3, 2>& tau, LINALG::Matrix<probdim, 1>& unit_normal) const;

      /// @name First derivatives
      /// @{

      void Deriv1st_NonUnitSlaveElementNormal(const MORTAR::MortarElement& sele,
          const LINALG::TMatrix<int, probdim, SLAVENUMNODE>& nodal_dofs,
          const LINALG::Matrix<SLAVEDIM, SLAVENUMNODE>& deriv, const LINALG::Matrix<3, 2>& tau,
          Deriv1stVecMap& d_non_unit_normal) const;

      void Deriv1st_UnitSlaveElementNormal(const LINALG::Matrix<probdim, 1>& unit_normal,
          const double length_n_inv, const Deriv1stVecMap& d_non_unit_normal,
          Deriv1stVecMap& d_unit_normal, const bool reset = true) const;

      void Deriv1st_Jacobian(const LINALG::Matrix<probdim, 1>& unit_normal,
          const Deriv1stVecMap& d_non_unit_normal, Deriv1stMap& d_jac) const;

      /// @}

      /// @name Second derivatives
      /// @{

      void Deriv2nd_Jacobian(const Deriv1stVecMap& d_unit_normal,
          const Deriv1stVecMap& d_non_unit_normal, const LINALG::Matrix<probdim, 1>& unit_normal,
          const Deriv2ndVecMap& dd_non_unit_normal, Deriv2ndMap& dd_jac) const;

      void Deriv2nd_NonUnitSlaveElementNormal(const MORTAR::MortarElement& sele,
          const LINALG::TMatrix<int, probdim, SLAVENUMNODE>& nodal_dofs,
          const LINALG::Matrix<SLAVEDIM, SLAVENUMNODE>& deriv,
          Deriv2ndVecMap& dd_non_unit_normal) const;

      void Deriv2nd_UnitSlaveElementNormal(const LINALG::Matrix<probdim, 1>& unit_normal,
          const double length_n_inv, const Deriv1stVecMap& d_non_unit_normal,
          const Deriv1stVecMap& d_unit_normal, const Deriv2ndVecMap& dd_non_unit_normal,
          Deriv2ndVecMap& dd_unit_normal) const;

      /// @}

      void AveragedNormalAtXi(MORTAR::MortarElement& sele,
          const LINALG::Matrix<SLAVENUMNODE, 1>& sval, LINALG::Matrix<probdim, 1>& snormal) const;

      /// complete data in nodal data container after successful integration
      void CompleteNodeData(MORTAR::MortarElement& sele) const;

     private:
      void InnerProductOfVectorAndDeriv1stVector(const LINALG::Matrix<probdim, 1>& vec,
          const Deriv1stVecMap& d_vec, Deriv1stMap& dvec_vec) const;

      void ProjectionIntoTangentialPlain(const LINALG::Matrix<probdim, 1>& unit_normal,
          LINALG::Matrix<probdim, probdim>& tproj_mat) const;

     protected:
      /// current Gauss point ID
      int gp_id_ = -1;
      mutable TimeMonitor<CONTACT::AUG::TimeID> timer_ = TimeMonitor<CONTACT::AUG::TimeID>();
    };

    /*--------------------------------------------------------------------------*/
    template <unsigned probdim, DRT::Element::DiscretizationType slavetype,
        DRT::Element::DiscretizationType mastertype>
    class BaseIntPolicy : public BaseSlaveIntPolicy<probdim, slavetype>
    {
      typedef BaseSlaveIntPolicy<probdim, slavetype> my;

     protected:
      static const unsigned MASTERDIM = DRT::UTILS::DisTypeToDim<mastertype>::dim;
      static const unsigned MASTERNUMNODE =
          DRT::UTILS::DisTypeToNumNodePerEle<mastertype>::numNodePerElement;

     public:
      /// constructor
      BaseIntPolicy(){/* empty */};

      /// destructor
      ~BaseIntPolicy(){/* empty */};

      void LMatrixInverse(const LINALG::Matrix<3, 2>& mtau,
          const LINALG::Matrix<probdim, 1>& snormal,
          LINALG::Matrix<probdim, probdim>& lmat_inv) const;

      void Deriv1st_MXiGP(const LINALG::Matrix<probdim, probdim>& lmat_inv,
          MORTAR::MortarElement& sele, MORTAR::MortarElement& mele,
          const LINALG::Matrix<my::SLAVENUMNODE, 1>& sval,
          const LINALG::Matrix<MASTERNUMNODE, 1>& mval, const double alpha, Deriv1stVecMap& d_mxi,
          Deriv1stMap& d_alpha) const;

      void Add_Deriv2nd_MaDispl(const LINALG::Matrix<probdim, probdim>& lmat_inv,
          const LINALG::TMatrix<int, probdim, MASTERNUMNODE>& mnodal_dofs,
          const LINALG::Matrix<MASTERDIM, MASTERNUMNODE>& mderiv, const Deriv1stVecMap& d_mxigp,
          Deriv2ndVecMap& dd_mxigp) const;

      void Add_Deriv1st_MaMetric(const LINALG::Matrix<probdim, probdim>& lmat_inv,
          const LINALG::TMatrix<int, probdim, MASTERNUMNODE>& mnodal_dofs,
          const LINALG::Matrix<3, 2>& mtau, const LINALG::Matrix<MASTERDIM, MASTERNUMNODE>& mderiv,
          const LINALG::Matrix<3, MASTERNUMNODE>& mderiv2nd,
          const LINALG::Matrix<3, MASTERNUMNODE>& mcoord, const Deriv1stVecMap& d_mxigp,
          Deriv2ndVecMap& dd_mxigp) const;

      void Add_Deriv1st_Alpha_Deriv1st_Normal(const LINALG::Matrix<probdim, probdim>& lmat_inv,
          const Deriv1stMap& d_alpha, const MORTAR::MortarElement& sele,
          const LINALG::Matrix<my::SLAVENUMNODE, 1>& sval, Deriv2ndVecMap& dd_mxigp) const;

      void Add_Alpha_Deriv2nd_Normal(const LINALG::Matrix<probdim, probdim>& lmat_inv,
          const double alpha, const MORTAR::MortarElement& sele,
          const LINALG::Matrix<my::SLAVENUMNODE, 1>& sval, Deriv2ndVecMap& dd_mxigp) const;

      void Deriv1st_GapN_Sl(const DRT::Node* const* snodes,
          const LINALG::Matrix<my::SLAVENUMNODE, 1>& sval, const double* gpn,
          Deriv1stMap& deriv_gapn_sl) const;

      void Deriv1st_GapN_Ma(const DRT::Node* const* mnodes,
          const LINALG::Matrix<MASTERNUMNODE, 1>& mval, const double* gpn,
          const LINALG::Matrix<3, 2>& mtau, const Deriv1stVecMap& d_mxigp,
          Deriv1stMap& deriv_gapn_ma) const;

      void Add_Deriv1stGapNContributions(DRT::Node* const* snodes, const double scale,
          const LINALG::Matrix<my::SLAVENUMNODE, 1>& lmval, const Deriv1stMap& d_gapn_sl,
          const Deriv1stMap& d_gapn_ma) const;

      void Add_Deriv1stJacobianContributions(DRT::Node* const* snodes, const double wgt,
          const LINALG::Matrix<my::SLAVENUMNODE, 1>& lmval, const double gapn_sl,
          const double gapn_ma, const Deriv1stMap& d_jac) const;

      void Add_Var_GapN_Lin_Jac(MORTAR::MortarElement& sele,
          const LINALG::Matrix<my::SLAVENUMNODE, 1>& lmval, const double wgt,
          const Deriv1stMap& d_gapn_sl, const Deriv1stMap& d_gapn_ma,
          const Deriv1stMap& d_jac) const;

      void Add_Var_Jac_Lin_GapN(MORTAR::MortarElement& sele,
          const LINALG::Matrix<my::SLAVENUMNODE, 1>& lmval, const double wgt,
          const Deriv1stMap& d_gapn_sl, const Deriv1stMap& d_gapn_ma,
          const Deriv1stMap& d_jac) const;

      void Add_GapN_Deriv2nd_Jac(MORTAR::MortarElement& sele,
          const LINALG::Matrix<my::SLAVENUMNODE, 1>& lmval, const double wgt, const double gapn_sl,
          const double gapn_ma, const Deriv2ndMap& dd_jac) const;

      /// @name empty functions, if no debug policy is in use
      /// @{

      void Get_Debug(...) const { return; }

      void Get_Deriv1st_Debug(...) const { return; }

      void Get_Deriv2nd_Debug(...) const { return; }

      /// @}
    };

    /*--------------------------------------------------------------------------*/
    template <unsigned probdim, DRT::Element::DiscretizationType slavetype,
        DRT::Element::DiscretizationType mastertype>
    class IncompleteIntPolicy : public BaseIntPolicy<probdim, slavetype, mastertype>
    {
      typedef BaseIntPolicy<probdim, slavetype, mastertype> my;

     public:
      /// constructor
      IncompleteIntPolicy(){/* empty */};

      /// destructor
      ~IncompleteIntPolicy(){/* empty */};

      void Get_Deriv2nd_Jacobian(const MORTAR::MortarElement& sele,
          const LINALG::TMatrix<int, probdim, my::SLAVENUMNODE>& nodal_dofs,
          const LINALG::Matrix<my::SLAVEDIM, my::SLAVENUMNODE>& deriv,
          const LINALG::Matrix<probdim, 1>& unit_normal, const double length_n_inv,
          const Deriv1stVecMap& d_non_unit_normal, Deriv2ndMap& dd_jac) const;

      void Get_Deriv2nd_MXiGP(const LINALG::Matrix<probdim, probdim>& lmat_inv,
          MORTAR::MortarElement& sele, MORTAR::MortarElement& mele,
          const LINALG::Matrix<my::SLAVENUMNODE, 1>& sval,
          const LINALG::Matrix<my::MASTERNUMNODE, 1>& mval,
          const LINALG::Matrix<my::MASTERDIM, my::MASTERNUMNODE>& mderiv,
          const LINALG::Matrix<3, my::MASTERNUMNODE>& mderiv2, const LINALG::Matrix<3, 2>& mtau,
          const double* mxi, const double alpha, const Deriv1stVecMap& d_mxigp,
          const Deriv1stMap& d_alpha, Deriv2ndVecMap& dd_mxigp) const
      {
        // do nothing
        GEN_DATA::reset(2, 0, dd_mxigp);
      }

      void Get_Deriv1st_GapN(MORTAR::MortarElement& sele, MORTAR::MortarElement& mele,
          const LINALG::Matrix<my::SLAVENUMNODE, 1>& sval,
          const LINALG::Matrix<my::MASTERNUMNODE, 1>& mval, const double* gpn,
          const LINALG::Matrix<3, 2>& mtau, const Deriv1stVecMap& d_mxigp,
          Deriv1stMap& deriv_gapn_sl, Deriv1stMap& deriv_gapn_ma) const;

      void Get_Deriv1st_WGap(MORTAR::MortarElement& sele,
          const LINALG::Matrix<my::SLAVENUMNODE, 1>& lmval, const double gapn_sl,
          const double gapn_ma, const double wgt, const double jac, const Deriv1stMap& d_jac,
          const Deriv1stMap& d_gapn_sl, const Deriv1stMap& d_gapn_ma) const;

      void Get_Deriv1st_WGap_Complete(const int linsize, MORTAR::MortarElement& sele,
          MORTAR::MortarElement& mele, const LINALG::Matrix<my::SLAVENUMNODE, 1>& sval,
          const LINALG::Matrix<my::MASTERNUMNODE, 1>& mval,
          const LINALG::Matrix<my::SLAVENUMNODE, 1>& lmval, const double* gpn,
          const LINALG::Matrix<3, 2>& mtau, const Deriv1stVecMap& d_mxigp, const double gapn_sl,
          const double gapn_ma, const double wgt, const double jac, const Deriv1stMap& d_jac) const;

      void Get_Deriv2nd_WGap(MORTAR::MortarElement& sele, MORTAR::MortarElement& mele,
          const LINALG::Matrix<my::SLAVENUMNODE, 1>& sval,
          const LINALG::Matrix<my::MASTERNUMNODE, 1>& mval,
          const LINALG::Matrix<my::SLAVENUMNODE, 1>& lmval,
          const LINALG::Matrix<my::MASTERDIM, my::MASTERNUMNODE>& mderiv,
          const LINALG::Matrix<3, my::MASTERNUMNODE>& mderiv2nd, const LINALG::Matrix<3, 2>& mtau,
          const double* gpn, const double wgt, const double gapn_sl, const double gapn_ma,
          const double jac, const Deriv1stMap& d_jac, const Deriv2ndMap& dd_jac,
          const Deriv1stVecMap& d_mxigp, const Deriv2ndVecMap& dd_mxigp,
          const Deriv1stVecMap& d_n_unit, const Deriv2ndVecMap& dd_n_unit,
          const Deriv1stMap& d_gapn_sl, const Deriv1stMap& d_gapn_ma) const;

      void Get_Deriv1st_WGapN_Error(const MORTAR::MortarElement& sele,
          const std::vector<unsigned>& active_nlids,
          const LINALG::Matrix<my::SLAVENUMNODE, 1>& lmval, const double* gpn, const double gapn_sl,
          const double gapn_ma, const double wgt, const double jacslave, const Deriv1stMap& d_jac,
          const LINALG::Matrix<3, 2>& mtau, const Deriv1stVecMap& d_mxigp, Deriv1stMap& d_gapn_ma,
          std::unordered_map<int, Deriv1stMap>& error_ma,
          std::unordered_map<int, Deriv1stMap>& error_jac) const;

     protected:
      void Add_Deriv1st_GapN_Deriv1st_Jac(MORTAR::MortarElement& sele,
          const LINALG::Matrix<my::SLAVENUMNODE, 1>& lmval, const double wgt,
          const Deriv1stMap& d_gapn_sl, const Deriv1stMap& d_gapn_ma,
          const Deriv1stMap& d_jac) const;

      void Add_Jac_Deriv2nd_GapN(MORTAR::MortarElement& sele, MORTAR::MortarElement& mele,
          const LINALG::Matrix<my::SLAVENUMNODE, 1>& sval,
          const LINALG::Matrix<my::MASTERNUMNODE, 1>& mval,
          const LINALG::Matrix<my::SLAVENUMNODE, 1>& lmval,
          const LINALG::Matrix<my::MASTERDIM, my::MASTERNUMNODE>& mderiv,
          const LINALG::Matrix<3, 2>& mtau, const double* gpn, const double wgt, const double jac,
          const Deriv1stVecMap& d_mxigp, const Deriv1stVecMap& d_n_unit,
          const Deriv2ndVecMap& dd_n_unit) const;
    };

    /*--------------------------------------------------------------------------*/
    template <unsigned probdim, DRT::Element::DiscretizationType slavetype,
        DRT::Element::DiscretizationType mastertype>
    class CompleteIntPolicy : public BaseIntPolicy<probdim, slavetype, mastertype>
    {
      typedef BaseIntPolicy<probdim, slavetype, mastertype> my;

     public:
      /// constructor
      CompleteIntPolicy(){/* empty */};

      /// destructor
      ~CompleteIntPolicy(){/* empty */};

      void Get_Deriv2nd_Jacobian(const MORTAR::MortarElement& sele,
          const LINALG::TMatrix<int, probdim, my::SLAVENUMNODE>& nodal_dofs,
          const LINALG::Matrix<my::SLAVEDIM, my::SLAVENUMNODE>& deriv,
          const LINALG::Matrix<probdim, 1>& unit_normal, const double length_n_inv,
          const Deriv1stVecMap& d_non_unit_normal, Deriv2ndMap& dd_jac) const;

      void Get_Deriv2nd_MXiGP(const LINALG::Matrix<probdim, probdim>& lmat_inv,
          MORTAR::MortarElement& sele, MORTAR::MortarElement& mele,
          const LINALG::Matrix<my::SLAVENUMNODE, 1>& sval,
          const LINALG::Matrix<my::MASTERNUMNODE, 1>& mval,
          const LINALG::Matrix<my::MASTERDIM, my::MASTERNUMNODE>& mderiv,
          const LINALG::Matrix<3, my::MASTERNUMNODE>& mderiv2, const LINALG::Matrix<3, 2>& mtau,
          const double* mxi, const double alpha, const Deriv1stVecMap& d_mxigp,
          const Deriv1stMap& d_alpha, Deriv2ndVecMap& dd_mxigp) const;

      void Get_Deriv1st_GapN(MORTAR::MortarElement& sele, MORTAR::MortarElement& mele,
          const LINALG::Matrix<my::SLAVENUMNODE, 1>& sval,
          const LINALG::Matrix<my::MASTERNUMNODE, 1>& mval, const double* gpn,
          const LINALG::Matrix<3, 2>& mtau, const Deriv1stVecMap& d_mxigp,
          Deriv1stMap& deriv_gapn_sl, Deriv1stMap& deriv_gapn_ma) const;

      void Get_Deriv1st_WGap(MORTAR::MortarElement& sele,
          const LINALG::Matrix<my::SLAVENUMNODE, 1>& lmval, const double gapn_sl,
          const double gapn_ma, const double wgt, const double jac, const Deriv1stMap& d_jac,
          const Deriv1stMap& d_gapn_sl, const Deriv1stMap& d_gapn_ma) const;

      void Get_Deriv1st_WGap_Complete(const int linsize, MORTAR::MortarElement& sele,
          MORTAR::MortarElement& mele, const LINALG::Matrix<my::SLAVENUMNODE, 1>& sval,
          const LINALG::Matrix<my::MASTERNUMNODE, 1>& mval,
          const LINALG::Matrix<my::SLAVENUMNODE, 1>& lmval, const double* gpn,
          const LINALG::Matrix<3, 2>& mtau, const Deriv1stVecMap& d_mxigp, const double gapn_sl,
          const double gapn_ma, const double wgt, const double jac, const Deriv1stMap& d_jac) const;

      void Get_Deriv2nd_WGap(MORTAR::MortarElement& sele, MORTAR::MortarElement& mele,
          const LINALG::Matrix<my::SLAVENUMNODE, 1>& sval,
          const LINALG::Matrix<my::MASTERNUMNODE, 1>& mval,
          const LINALG::Matrix<my::SLAVENUMNODE, 1>& lmval,
          const LINALG::Matrix<my::MASTERDIM, my::MASTERNUMNODE>& mderiv,
          const LINALG::Matrix<3, my::MASTERNUMNODE>& mderiv2nd, const LINALG::Matrix<3, 2>& mtau,
          const double* gpn, const double wgt, const double gapn_sl, const double gapn_ma,
          const double jac, const Deriv1stMap& d_jac, const Deriv2ndMap& dd_jac,
          const Deriv1stVecMap& d_mxigp, const Deriv2ndVecMap& dd_mxigp,
          const Deriv1stVecMap& d_n_unit, const Deriv2ndVecMap& dd_n_unit,
          const Deriv1stMap& d_gapn_sl, const Deriv1stMap& d_gapn_ma) const;

      void Get_Deriv1st_WGapN_Error(const MORTAR::MortarElement& sele,
          const std::vector<unsigned>& active_nlids,
          const LINALG::Matrix<my::SLAVENUMNODE, 1>& lmval, const double* gpn, const double gapn_sl,
          const double gapn_ma, const double wgt, const double jacslave, const Deriv1stMap& d_jac,
          const LINALG::Matrix<3, 2>& mtau, const Deriv1stVecMap& d_mxigp, Deriv1stMap& d_gapn_ma,
          std::unordered_map<int, Deriv1stMap>& error_ma,
          std::unordered_map<int, Deriv1stMap>& error_jac) const;

     protected:
      void Add_Deriv1st_GapN_Deriv1st_Jac(MORTAR::MortarElement& sele,
          const LINALG::Matrix<my::SLAVENUMNODE, 1>& lmval, const double wgt,
          const Deriv1stMap& d_gapn_sl, const Deriv1stMap& d_gapn_ma,
          const Deriv1stMap& d_jac) const;

      void Add_Jac_Deriv2nd_GapN(MORTAR::MortarElement& sele, MORTAR::MortarElement& mele,
          const LINALG::Matrix<my::SLAVENUMNODE, 1>& sval,
          const LINALG::Matrix<my::MASTERNUMNODE, 1>& mval,
          const LINALG::Matrix<my::SLAVENUMNODE, 1>& lmval,
          const LINALG::Matrix<my::MASTERDIM, my::MASTERNUMNODE>& mderiv,
          const LINALG::Matrix<3, my::MASTERNUMNODE>& mderiv2nd, const LINALG::Matrix<3, 2>& mtau,
          const double* gpn, const double wgt, const double jac, const Deriv1stVecMap& d_mxigp,
          const Deriv2ndVecMap& dd_mxigp, const Deriv1stVecMap& d_n_unit,
          const Deriv2ndVecMap& dd_n_unit) const;
    };

    /*--------------------------------------------------------------------------*/
    template <unsigned probdim, DRT::Element::DiscretizationType slavetype,
        DRT::Element::DiscretizationType mastertype>
    class DebugIncompleteIntPolicy : public IncompleteIntPolicy<probdim, slavetype, mastertype>
    {
      typedef IncompleteIntPolicy<probdim, slavetype, mastertype> my;

     public:
      /// constructor
      DebugIncompleteIntPolicy(){/* empty */};

      /// destructor
      ~DebugIncompleteIntPolicy(){/* empty */};

      void Get_Debug(MORTAR::MortarElement& sele, const LINALG::Matrix<my::SLAVENUMNODE, 1>& lmval,
          const double gapn_sl, const double gapn_ma, const double wgt, const double jac,
          const double* gpn, const double* mxigp) const;

      void Get_Deriv1st_Debug(MORTAR::MortarElement& sele,
          const LINALG::Matrix<my::SLAVENUMNODE, 1>& lmval,
          const LINALG::Matrix<my::SLAVENUMNODE, 1>& sval,
          const LINALG::Matrix<my::SLAVEDIM, my::SLAVENUMNODE>& sderiv,
          const LINALG::Matrix<3, 2>& stau, const Deriv1stMap& d_jac, const Deriv1stVecMap& dmxigp,
          const Deriv1stVecMap& d_gpn, const Deriv1stMap& d_gap_sl, const double gapn_sl,
          const double wgt, const double jac) const;

      void Get_Deriv2nd_Debug(MORTAR::MortarElement& sele,
          const LINALG::Matrix<my::SLAVENUMNODE, 1>& lmval,
          const LINALG::Matrix<my::SLAVENUMNODE, 1>& sval,
          const LINALG::Matrix<my::SLAVEDIM, my::SLAVENUMNODE>& sderiv,
          const LINALG::Matrix<3, 2>& stau, const Deriv1stMap& d_jac, const Deriv1stMap& d_gapn_sl,
          const Deriv2ndMap& dd_jac, const Deriv2ndVecMap& ddmxigp, const Deriv1stVecMap& d_n_unit,
          const Deriv2ndVecMap& dd_n_unit, const double gapn_sl, const double wgt,
          const double jac) const;

     private:
      void Debug_Deriv1st_WGap(MORTAR::MortarElement& sele) const;

      void Debug_Deriv2nd_WGap(MORTAR::MortarElement& sele) const;
    };

    /*--------------------------------------------------------------------------*/
    template <unsigned probdim, DRT::Element::DiscretizationType slavetype,
        DRT::Element::DiscretizationType mastertype>
    class DebugCompleteIntPolicy : public CompleteIntPolicy<probdim, slavetype, mastertype>
    {
      typedef CompleteIntPolicy<probdim, slavetype, mastertype> my;

     public:
      /// constructor
      DebugCompleteIntPolicy(){/* empty */};

      /// destructor
      ~DebugCompleteIntPolicy(){/* empty */};

      void Get_Debug(MORTAR::MortarElement& sele, const LINALG::Matrix<my::SLAVENUMNODE, 1>& lmval,
          const double gapn_sl, const double gapn_ma, const double wgt, const double jac,
          const double* gpn, const double* mxigp) const;

      void Get_Deriv1st_Debug(MORTAR::MortarElement& sele,
          const LINALG::Matrix<my::SLAVENUMNODE, 1>& lmval,
          const LINALG::Matrix<my::SLAVENUMNODE, 1>& sval,
          const LINALG::Matrix<my::SLAVEDIM, my::SLAVENUMNODE>& sderiv,
          const LINALG::Matrix<3, 2>& stau, const Deriv1stMap& d_jac, const Deriv1stVecMap& dmxigp,
          const Deriv1stVecMap& d_gpn, const Deriv1stMap& d_gap_sl, const double gapn_sl,
          const double wgt, const double jac) const;

      void Get_Deriv2nd_Debug(MORTAR::MortarElement& sele,
          const LINALG::Matrix<my::SLAVENUMNODE, 1>& lmval,
          const LINALG::Matrix<my::SLAVENUMNODE, 1>& sval,
          const LINALG::Matrix<my::SLAVEDIM, my::SLAVENUMNODE>& sderiv,
          const LINALG::Matrix<3, 2>& stau, const Deriv1stMap& d_jac, const Deriv1stMap& d_gapn_sl,
          const Deriv2ndMap& dd_jac, const Deriv2ndVecMap& ddmxigp, const Deriv1stVecMap& d_n_unit,
          const Deriv2ndVecMap& dd_n_unit, const double gapn_sl, const double wgt,
          const double jac) const;

     private:
      void Debug_Kappa(MORTAR::MortarElement& sele,
          const LINALG::Matrix<my::SLAVENUMNODE, 1>& lmval, const double jac,
          const double wgt) const;

      void Debug_Deriv1st_Kappa(MORTAR::MortarElement& sele,
          const LINALG::Matrix<my::SLAVENUMNODE, 1>& lmval, const Deriv1stMap& djac,
          const double wgt) const;

      void Debug_GPN(MORTAR::MortarElement& sele, const LINALG::Matrix<my::SLAVENUMNODE, 1>& lmval,
          const double* gpn) const;

      void Debug_Deriv1st_GPN(MORTAR::MortarElement& sele,
          const LINALG::Matrix<my::SLAVENUMNODE, 1>& lmval, const Deriv1stVecMap& d_gpn) const;

      void Debug_WGap(MORTAR::MortarElement& sele, const LINALG::Matrix<my::SLAVENUMNODE, 1>& lmval,
          const double gapn_sl, const double gapn_ma, const double wgt, const double jac) const;

      void Debug_WGap_Sl(MORTAR::MortarElement& sele,
          const LINALG::Matrix<my::SLAVENUMNODE, 1>& lmval, const double gapn_sl, const double wgt,
          const double jac) const;

      void Debug_Deriv1st_Gap_Sl(MORTAR::MortarElement& sele,
          const LINALG::Matrix<my::SLAVENUMNODE, 1>& lmval,
          const LINALG::Matrix<my::SLAVENUMNODE, 1>& sval, const Deriv1stVecMap& d_n_unit,
          const Deriv1stMap& d_gapn_sl, const Deriv1stMap& d_jac, const double gapn_sl,
          const double wgt, const double jac, Deriv1stMap& d_gapn_sl_complete) const;

      void Debug_Deriv1st_WGap_Sl(MORTAR::MortarElement& sele,
          const LINALG::Matrix<my::SLAVENUMNODE, 1>& lmval,
          const LINALG::Matrix<my::SLAVENUMNODE, 1>& sval, const Deriv1stVecMap& d_gpn,
          const Deriv1stMap& d_gapn_sl, const Deriv1stMap& d_jac, const double gapn_sl,
          const double wgt, const double jac) const;

      void Debug_Deriv2nd_WGap_Sl(MORTAR::MortarElement& sele,
          const LINALG::Matrix<my::SLAVENUMNODE, 1>& lmval,
          const LINALG::Matrix<my::SLAVENUMNODE, 1>& sval, const Deriv1stVecMap& d_n_unit,
          const Deriv2ndVecMap& dd_n_unit, const Deriv1stMap& d_jac, const Deriv1stMap& d_gap_sl,
          const Deriv2ndMap& dd_jac, const double gapn_sl, const double wgt,
          const double jac) const;

      void Debug_Deriv1st_WGap(MORTAR::MortarElement& sele) const;

      void Debug_Deriv2nd_WGap(MORTAR::MortarElement& sele) const;

      void Debug_MXi(MORTAR::MortarElement& sele, const LINALG::Matrix<my::SLAVENUMNODE, 1>& lmval,
          const double* mxigp) const;

      void Debug_Deriv1st_MXi(MORTAR::MortarElement& sele,
          const LINALG::Matrix<my::SLAVENUMNODE, 1>& lmval, const Deriv1stVecMap& dmxigp) const;

      void Debug_Deriv2nd_MXi(MORTAR::MortarElement& sele,
          const LINALG::Matrix<my::SLAVENUMNODE, 1>& lmval, const Deriv2ndVecMap& ddmxigp) const;

      void Debug_Deriv1st_Smooth_Unit_Normal(MORTAR::MortarElement& sele) const;

      void Debug_Deriv2nd_Smooth_Unit_Normal(MORTAR::MortarElement& sele) const;

      void Debug_Deriv1st_Jac(MORTAR::MortarElement& sele,
          const LINALG::Matrix<my::SLAVENUMNODE, 1>& lmval, const Deriv1stMap& djac) const;

      void Debug_Deriv2nd_Jac(MORTAR::MortarElement& sele,
          const LINALG::Matrix<my::SLAVENUMNODE, 1>& lmval, const Deriv2ndMap& ddjac) const;

      void Debug_Deriv1st_Jacobian(MORTAR::MortarElement& sele,
          const LINALG::Matrix<my::SLAVENUMNODE, 1>& lmval,
          const LINALG::Matrix<my::SLAVEDIM, my::SLAVENUMNODE>& sderiv,
          const LINALG::Matrix<3, 2>& stau) const;

      void Debug_Deriv2nd_Jacobian(MORTAR::MortarElement& sele,
          const LINALG::Matrix<my::SLAVENUMNODE, 1>& lmval,
          const LINALG::Matrix<my::SLAVEDIM, my::SLAVENUMNODE>& sderiv,
          const LINALG::Matrix<3, 2>& stau) const;

      void Debug_Deriv1st_Non_Unit_Normal(MORTAR::MortarElement& sele,
          const LINALG::Matrix<my::SLAVENUMNODE, 1>& lmval,
          const LINALG::Matrix<my::SLAVEDIM, my::SLAVENUMNODE>& sderiv,
          const LINALG::Matrix<3, 2>& stau) const;

      void Debug_Deriv2nd_Non_Unit_Normal(MORTAR::MortarElement& sele,
          const LINALG::Matrix<my::SLAVENUMNODE, 1>& lmval,
          const LINALG::Matrix<my::SLAVEDIM, my::SLAVENUMNODE>& sderiv) const;

     private:
      const Deriv1stMap& GetNodalDeriv1st(NodeDataContainer& data) const;

      const Deriv2ndMap& GetNodalDeriv2nd(NodeDataContainer& data) const;
    };
  }  // namespace AUG
}  // namespace CONTACT


#endif /* SRC_DRT_CONTACT_AUG_CONTACT_AUG_INTEGRATOR_POLICY_H_ */
