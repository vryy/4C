/*----------------------------------------------------------------------------*/
/*! \file
\brief Debugging contact integration policies

\level 3

*/
/*----------------------------------------------------------------------------*/

#include "baci_contact_aug_integrator_policy.hpp"
#include "baci_contact_node.hpp"
#include "baci_mortar_element.hpp"

FOUR_C_NAMESPACE_OPEN

// #define DD_NON_UNIT_NORMAL
// #define DD_JACOBIAN
// #define DD_SMOOTH_UNIT_NORMAL
// #define DD_MXI
// #define DD_WGAP
// #define DD_WGAP_SL
// #define D_MXI
// #define D_WGAP_SL
// #define D_WGAP
// #define D_GPN               // 1-st order derivative of the GP normal
// #define D_KAPPA

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template <unsigned probdim, CORE::FE::CellType slavetype, CORE::FE::CellType mastertype>
void CONTACT::AUG::DebugIncompleteIntPolicy<probdim, slavetype, mastertype>::Get_Debug(
    MORTAR::Element& sele, const CORE::LINALG::Matrix<my::SLAVENUMNODE, 1>& lmval,
    const double gapn_sl, const double gapn_ma, const double wgt, const double jac,
    const double* gpn, const double* mxigp) const
{
  /* empty */
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template <unsigned probdim, CORE::FE::CellType slavetype, CORE::FE::CellType mastertype>
void CONTACT::AUG::DebugIncompleteIntPolicy<probdim, slavetype, mastertype>::Get_Deriv1st_Debug(
    MORTAR::Element& sele, const CORE::LINALG::Matrix<my::SLAVENUMNODE, 1>& lmval,
    const CORE::LINALG::Matrix<my::SLAVENUMNODE, 1>& sval,
    const CORE::LINALG::Matrix<my::SLAVEDIM, my::SLAVENUMNODE>& sderiv,
    const CORE::LINALG::Matrix<3, 2>& stau, const Deriv1stMap& d_jac, const Deriv1stVecMap& dmxigp,
    const Deriv1stVecMap& d_gpn, const Deriv1stMap& d_gap_sl, const double gapn_sl,
    const double wgt, const double jac) const
{
  this->CompleteNodeData(sele);

#ifdef DD_WGAP
  Debug_Deriv1st_WGap(sele);
#endif
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template <unsigned probdim, CORE::FE::CellType slavetype, CORE::FE::CellType mastertype>
void CONTACT::AUG::DebugIncompleteIntPolicy<probdim, slavetype, mastertype>::Get_Deriv2nd_Debug(
    MORTAR::Element& sele, const CORE::LINALG::Matrix<my::SLAVENUMNODE, 1>& lmval,
    const CORE::LINALG::Matrix<my::SLAVENUMNODE, 1>& sval,
    const CORE::LINALG::Matrix<my::SLAVEDIM, my::SLAVENUMNODE>& sderiv,
    const CORE::LINALG::Matrix<3, 2>& stau, const Deriv1stMap& d_jac, const Deriv1stMap& d_gapn_sl,
    const Deriv2ndMap& dd_jac, const Deriv2ndVecMap& ddmxigp, const Deriv1stVecMap& d_n_unit,
    const Deriv2ndVecMap& dd_n_unit, const double gapn_sl, const double wgt, const double jac) const
{
  this->CompleteNodeData(sele);

#ifdef DD_WGAP
  Debug_Deriv2nd_WGap(sele);
#endif
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template <unsigned probdim, CORE::FE::CellType slavetype, CORE::FE::CellType mastertype>
void CONTACT::AUG::DebugIncompleteIntPolicy<probdim, slavetype, mastertype>::Debug_Deriv1st_WGap(
    MORTAR::Element& sele) const
{
  DRT::Node* const* snodes = sele.Nodes();

  for (unsigned i = 0; i < my::SLAVENUMNODE; ++i)
  {
    Node& cnode = dynamic_cast<Node&>(*snodes[i]);

    Deriv1stMap& d_debug = cnode.AugData().GetDeriv1st_Debug();
    CORE::GEN::copy(cnode.AugData().GetDeriv1st_WGapMa(), d_debug);

    const Deriv1stMap& d_wgap_sl = cnode.AugData().GetDeriv1st_WGapSl();

    for (auto& d_wgap_sl_var : d_wgap_sl)
    {
      CORE::GEN::increaseCapacity(d_debug);
      d_debug[d_wgap_sl_var.first] -= d_wgap_sl_var.second;
    }
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template <unsigned probdim, CORE::FE::CellType slavetype, CORE::FE::CellType mastertype>
void CONTACT::AUG::DebugIncompleteIntPolicy<probdim, slavetype, mastertype>::Debug_Deriv2nd_WGap(
    MORTAR::Element& sele) const
{
  DRT::Node* const* snodes = sele.Nodes();

  for (unsigned i = 0; i < my::SLAVENUMNODE; ++i)
  {
    Node& cnode = dynamic_cast<Node&>(*snodes[i]);

    Deriv2ndMap& dd_debug = cnode.AugData().GetDeriv2nd_Debug();
    CORE::GEN::copy(cnode.AugData().GetDeriv2nd_WGapMa(), dd_debug);

    const Deriv2ndMap& dd_wgap_sl = cnode.AugData().GetDeriv2nd_WGapSl();
    for (auto& dd_wgap_sl_var : dd_wgap_sl)
    {
      CORE::GEN::increaseCapacity(dd_debug);
      Deriv1stMap& dd_debug_var = dd_debug[dd_wgap_sl_var.first];

      for (auto& dd_wgap_sl_var_lin : dd_wgap_sl_var.second)
      {
        CORE::GEN::increaseCapacity(dd_debug_var);
        dd_debug_var[dd_wgap_sl_var_lin.first] -= dd_wgap_sl_var_lin.second;
      }
    }
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template <unsigned probdim, CORE::FE::CellType slavetype, CORE::FE::CellType mastertype>
void CONTACT::AUG::DebugCompleteIntPolicy<probdim, slavetype, mastertype>::Get_Debug(
    MORTAR::Element& sele, const CORE::LINALG::Matrix<my::SLAVENUMNODE, 1>& lmval,
    const double gapn_sl, const double gapn_ma, const double wgt, const double jac,
    const double* gpn, const double* mxigp) const
{
  this->CompleteNodeData(sele);

#ifdef D_MXI
  Debug_MXi(sele, lmval, mxigp);
#endif

#ifdef D_WGAP
  Debug_WGap(sele, lmval, gapn_sl, gapn_ma, wgt, jac);
#endif

#ifdef D_GPN
  Debug_GPN(sele, lmval, gpn);
#endif

#ifdef D_WGAP_SL
  Debug_WGap_Sl(sele, lmval, gapn_sl, wgt, jac);
#endif

#ifdef D_KAPPA
  Debug_Kappa(sele, lmval, jac, wgt);
#endif
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template <unsigned probdim, CORE::FE::CellType slavetype, CORE::FE::CellType mastertype>
void CONTACT::AUG::DebugCompleteIntPolicy<probdim, slavetype, mastertype>::Get_Deriv1st_Debug(
    MORTAR::Element& sele, const CORE::LINALG::Matrix<my::SLAVENUMNODE, 1>& lmval,
    const CORE::LINALG::Matrix<my::SLAVENUMNODE, 1>& sval,
    const CORE::LINALG::Matrix<my::SLAVEDIM, my::SLAVENUMNODE>& sderiv,
    const CORE::LINALG::Matrix<3, 2>& stau, const Deriv1stMap& d_jac, const Deriv1stVecMap& dmxigp,
    const Deriv1stVecMap& d_gpn, const Deriv1stMap& d_gap_sl, const double gapn_sl,
    const double wgt, const double jac) const
{
  this->CompleteNodeData(sele);

  // non-unit normal
#ifdef DD_NON_UNIT_NORMAL
  Debug_Deriv1st_Non_Unit_Normal(sele, lmval, sderiv, stau);
#endif

  // jacobian
#ifdef DD_JACOBIAN
  //  Debug_Deriv1st_Jacobian( sele, lmval, sderiv, stau );
  Debug_Deriv1st_Jac(sele, lmval, djac);
#endif

#ifdef DD_SMOOTH_UNIT_NORMAL
  Debug_Deriv1st_Smooth_Unit_Normal(sele);
#endif

#if defined(DD_MXI) or defined(D_MXI)
  Debug_Deriv1st_MXi(sele, lmval, dmxigp);
#endif

  // weighetd normal gap (slave contributions)
#if defined(DD_WGAP) or defined(D_WGAP)
  Debug_Deriv1st_WGap(sele);
#endif

#ifdef D_GPN
  Debug_Deriv1st_GPN(sele, lmval, d_gpn);
#endif

#if defined(DD_WGAP_SL) or defined(D_WGAP_SL)

  Deriv1stMap d_gapn_sl_complete;
  Debug_Deriv1st_Gap_Sl(
      sele, lmval, sval, d_gpn, d_gap_sl, djac, gapn_sl, wgt, jac, d_gapn_sl_complete);

  Debug_Deriv1st_WGap_Sl(sele, lmval, sval, d_gpn, d_gapn_sl_complete, djac, gapn_sl, wgt, jac);

#endif

#ifdef D_KAPPA
  Debug_Deriv1st_Kappa(sele, lmval, djac, wgt);
#endif
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template <unsigned probdim, CORE::FE::CellType slavetype, CORE::FE::CellType mastertype>
void CONTACT::AUG::DebugCompleteIntPolicy<probdim, slavetype, mastertype>::Get_Deriv2nd_Debug(
    MORTAR::Element& sele, const CORE::LINALG::Matrix<my::SLAVENUMNODE, 1>& lmval,
    const CORE::LINALG::Matrix<my::SLAVENUMNODE, 1>& sval,
    const CORE::LINALG::Matrix<my::SLAVEDIM, my::SLAVENUMNODE>& sderiv,
    const CORE::LINALG::Matrix<3, 2>& stau, const Deriv1stMap& d_jac, const Deriv1stMap& d_gapn_sl,
    const Deriv2ndMap& dd_jac, const Deriv2ndVecMap& ddmxigp, const Deriv1stVecMap& d_n_unit,
    const Deriv2ndVecMap& dd_n_unit, const double gapn_sl, const double wgt, const double jac) const
{
  this->CompleteNodeData(sele);

  // non-unit normal
#ifdef DD_NON_UNIT_NORMAL
  Debug_Deriv2nd_Non_Unit_Normal(sele, lmval, sderiv);
#endif

  // jacobian
#ifdef DD_JACOBIAN
  //  Debug_Deriv2nd_Jacobian( sele, lmval, sderiv, stau );
  Debug_Deriv2nd_Jac(sele, lmval, ddjac);
#endif

#ifdef DD_SMOOTH_UNIT_NORMAL
  Debug_Deriv2nd_Smooth_Unit_Normal(sele);
#endif

#ifdef DD_MXI
  Debug_Deriv2nd_MXi(sele, lmval, ddmxigp);
#endif

  // weighetd normal gap (slave contributions)
#ifdef DD_WGAP
  Debug_Deriv2nd_WGap(sele);
#endif

#ifdef DD_WGAP_SL
  Deriv1stMap d_gapn_sl_complete;
  Debug_Deriv1st_Gap_Sl(
      sele, lmval, sval, d_n_unit, dgap_sl, djac, gapn_sl, wgt, jac, d_gapn_sl_complete);

  Debug_Deriv2nd_WGap_Sl(
      sele, lmval, sval, d_n_unit, dd_n_unit, djac, d_gapn_sl_complete, ddjac, gapn_sl, wgt, jac);
#endif
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template <unsigned probdim, CORE::FE::CellType slavetype, CORE::FE::CellType mastertype>
const CONTACT::AUG::Deriv1stMap&
CONTACT::AUG::DebugCompleteIntPolicy<probdim, slavetype, mastertype>::GetNodalDeriv1st(
    NodeDataContainer& data) const
{
  FOUR_C_THROW("No active scalar first derivative!");
  exit(EXIT_FAILURE);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template <unsigned probdim, CORE::FE::CellType slavetype, CORE::FE::CellType mastertype>
const CONTACT::AUG::Deriv2ndMap&
CONTACT::AUG::DebugCompleteIntPolicy<probdim, slavetype, mastertype>::GetNodalDeriv2nd(
    NodeDataContainer& data) const
{
  FOUR_C_THROW("No active scalar second derivative!");
  exit(EXIT_FAILURE);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template <unsigned probdim, CORE::FE::CellType slavetype, CORE::FE::CellType mastertype>
void CONTACT::AUG::DebugCompleteIntPolicy<probdim, slavetype, mastertype>::Debug_Kappa(
    MORTAR::Element& sele, const CORE::LINALG::Matrix<my::SLAVENUMNODE, 1>& lmval, const double jac,
    const double wgt) const
{
  DRT::Node* const* snodes = sele.Nodes();

  for (unsigned i = 0; i < my::SLAVENUMNODE; ++i)
  {
    Node& cnode = dynamic_cast<Node&>(*snodes[i]);
    std::pair<int, double>& val = cnode.AugData().Get_Debug();

    val.first = cnode.Id();
    val.second = cnode.AugData().GetKappa();
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template <unsigned probdim, CORE::FE::CellType slavetype, CORE::FE::CellType mastertype>
void CONTACT::AUG::DebugCompleteIntPolicy<probdim, slavetype, mastertype>::Debug_Deriv1st_Kappa(
    MORTAR::Element& sele, const CORE::LINALG::Matrix<my::SLAVENUMNODE, 1>& lmval,
    const Deriv1stMap& djac, const double wgt) const
{
  DRT::Node* const* snodes = sele.Nodes();

  for (unsigned i = 0; i < my::SLAVENUMNODE; ++i)
  {
    Node& cnode = dynamic_cast<Node&>(*snodes[i]);
    Deriv1stMap& d_debug = cnode.AugData().GetDeriv1st_Debug();

    //    CORE::GEN::copy( cnode.AugData().GetDeriv1st_Kappa(), d_debug );
    for (auto& djac_var : djac)
    {
      CORE::GEN::increaseCapacity(d_debug);
      d_debug[djac_var.first] += wgt * lmval(i, 0) * djac_var.second;
    }
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template <unsigned probdim, CORE::FE::CellType slavetype, CORE::FE::CellType mastertype>
void CONTACT::AUG::DebugCompleteIntPolicy<probdim, slavetype, mastertype>::Debug_GPN(
    MORTAR::Element& sele, const CORE::LINALG::Matrix<my::SLAVENUMNODE, 1>& lmval,
    const double* gpn) const
{
  DRT::Node* const* snodes = sele.Nodes();

  for (unsigned i = 0; i < my::SLAVENUMNODE; ++i)
  {
    Node& cnode = dynamic_cast<Node&>(*snodes[i]);
    std::vector<std::pair<int, double>>& vals = cnode.AugData().Get_DebugVec();

    for (unsigned d = 0; d < probdim; ++d)
    {
      vals[d].first = cnode.Id();
      vals[d].second += lmval(i, 0) * gpn[d];
    }
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template <unsigned probdim, CORE::FE::CellType slavetype, CORE::FE::CellType mastertype>
void CONTACT::AUG::DebugCompleteIntPolicy<probdim, slavetype, mastertype>::Debug_Deriv1st_GPN(
    MORTAR::Element& sele, const CORE::LINALG::Matrix<my::SLAVENUMNODE, 1>& lmval,
    const Deriv1stVecMap& d_gpn) const
{
  DRT::Node* const* snodes = sele.Nodes();

  for (unsigned i = 0; i < my::SLAVENUMNODE; ++i)
  {
    Node& cnode = dynamic_cast<Node&>(*snodes[i]);

    Deriv1stVecMap& d_debug = cnode.AugData().GetDeriv1st_DebugVec();
    if (d_debug.size() < probdim)
    {
      d_debug.resize(probdim);
      CORE::GEN::reset(1, d_debug);
    }

    for (unsigned d = 0; d < probdim; ++d)
    {
      for (auto& d_gpn_var : d_gpn[d])
      {
        CORE::GEN::increaseCapacity(d_debug[d]);
        d_debug[d][d_gpn_var.first] += lmval(i, 0) * d_gpn_var.second;
      }
    }
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template <unsigned probdim, CORE::FE::CellType slavetype, CORE::FE::CellType mastertype>
void CONTACT::AUG::DebugCompleteIntPolicy<probdim, slavetype, mastertype>::Debug_WGap(
    MORTAR::Element& sele, const CORE::LINALG::Matrix<my::SLAVENUMNODE, 1>& lmval,
    const double gapn_sl, const double gapn_ma, const double wgt, const double jac) const
{
  DRT::Node* const* snodes = sele.Nodes();

  for (unsigned i = 0; i < my::SLAVENUMNODE; ++i)
  {
    Node& cnode = dynamic_cast<Node&>(*snodes[i]);
    std::pair<int, double>& scalar = cnode.AugData().Get_Debug();

    scalar.first = cnode.Id();

    scalar.second += wgt * lmval(i, 0) * jac * (gapn_ma - gapn_sl);
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template <unsigned probdim, CORE::FE::CellType slavetype, CORE::FE::CellType mastertype>
void CONTACT::AUG::DebugCompleteIntPolicy<probdim, slavetype, mastertype>::Debug_WGap_Sl(
    MORTAR::Element& sele, const CORE::LINALG::Matrix<my::SLAVENUMNODE, 1>& lmval,
    const double gapn_sl, const double wgt, const double jac) const
{
  DRT::Node* const* snodes = sele.Nodes();

  for (unsigned i = 0; i < my::SLAVENUMNODE; ++i)
  {
    Node& cnode = dynamic_cast<Node&>(*snodes[i]);
    std::pair<int, double>& scalar = cnode.AugData().Get_Debug();

    scalar.first = cnode.Id();

    scalar.second += wgt * lmval(i, 0) * jac * (gapn_sl);
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template <unsigned probdim, CORE::FE::CellType slavetype, CORE::FE::CellType mastertype>
void CONTACT::AUG::DebugCompleteIntPolicy<probdim, slavetype, mastertype>::Debug_Deriv1st_Gap_Sl(
    MORTAR::Element& sele, const CORE::LINALG::Matrix<my::SLAVENUMNODE, 1>& lmval,
    const CORE::LINALG::Matrix<my::SLAVENUMNODE, 1>& sval, const Deriv1stVecMap& d_n_unit,
    const Deriv1stMap& d_gapn_sl, const Deriv1stMap& d_jac, const double gapn_sl, const double wgt,
    const double jac, Deriv1stMap& d_gapn_sl_complete) const
{
  DRT::Node** snodes = sele.Nodes();

  // --- standard slave gap part --------------------------------------------
  for (auto& d_gapn_sl_var : d_gapn_sl)
  {
    CORE::GEN::increaseCapacity(d_gapn_sl_complete);
    d_gapn_sl_complete[d_gapn_sl_var.first] += d_gapn_sl_var.second;
  }

  // --- Additional part for debugging only ----------------------------------
  CORE::LINALG::Matrix<probdim, 1> xs_int(true);
  for (unsigned k = 0; k < my::SLAVENUMNODE; ++k)
  {
    const Node& snode = static_cast<const Node&>(*snodes[k]);
    const double* xs = snode.xspatial();

    for (unsigned j = 0; j < probdim; ++j)
    {
      xs_int(j, 0) += sval(k, 0) * xs[j];
    }
  }

  for (unsigned j = 0; j < probdim; ++j)
  {
    for (auto& d_n_unit_j_var : d_n_unit[j])
    {
      CORE::GEN::increaseCapacity(d_gapn_sl_complete);
      d_gapn_sl_complete[d_n_unit_j_var.first] += d_n_unit_j_var.second * xs_int(j, 0);
    }
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template <unsigned probdim, CORE::FE::CellType slavetype, CORE::FE::CellType mastertype>
void CONTACT::AUG::DebugCompleteIntPolicy<probdim, slavetype, mastertype>::Debug_Deriv1st_WGap_Sl(
    MORTAR::Element& sele, const CORE::LINALG::Matrix<my::SLAVENUMNODE, 1>& lmval,
    const CORE::LINALG::Matrix<my::SLAVENUMNODE, 1>& sval, const Deriv1stVecMap& d_gpn,
    const Deriv1stMap& d_gapn_sl_complete, const Deriv1stMap& d_jac, const double gapn_sl,
    const double wgt, const double jac) const
{
  DRT::Node* const* snodes = sele.Nodes();

  for (unsigned i = 0; i < my::SLAVENUMNODE; ++i)
  {
    const double tmp = wgt * jac * lmval(i, 0);

    Node& cnode = static_cast<Node&>(*snodes[i]);

    Deriv1stMap& d_debug = cnode.AugData().GetDeriv1st_Debug();
    if (d_debug.capacity() == 0)
    {
      CORE::GEN::reset(1, d_debug);
    }

    // --- slave gap part -----------------------------------------------------
    for (auto& d_gapn_sl_var : d_gapn_sl_complete)
    {
      CORE::GEN::increaseCapacity(d_debug);
      d_debug[d_gapn_sl_var.first] += tmp * d_gapn_sl_var.second;
    }

    // --- jacobian part ------------------------------------------------------
    const double tmp_sl = lmval(i, 0) * wgt * gapn_sl;

    for (auto& d_jac_var : d_jac)
    {
      CORE::GEN::increaseCapacity(d_debug);
      d_debug[d_jac_var.first] += tmp_sl * d_jac_var.second;
    }
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template <unsigned probdim, CORE::FE::CellType slavetype, CORE::FE::CellType mastertype>
void CONTACT::AUG::DebugCompleteIntPolicy<probdim, slavetype, mastertype>::Debug_Deriv2nd_WGap_Sl(
    MORTAR::Element& sele, const CORE::LINALG::Matrix<my::SLAVENUMNODE, 1>& lmval,
    const CORE::LINALG::Matrix<my::SLAVENUMNODE, 1>& sval, const Deriv1stVecMap& d_n_unit,
    const Deriv2ndVecMap& dd_n_unit, const Deriv1stMap& d_jac,
    const Deriv1stMap& d_gapn_sl_complete, const Deriv2ndMap& dd_jac, const double gapn_sl,
    const double wgt, const double jac) const
{
  DRT::Node* const* snodes = sele.Nodes();

  for (unsigned i = 0; i < my::SLAVENUMNODE; ++i)
  {
    const double tmp = wgt * jac * lmval(i, 0);

    Node& cnode = static_cast<Node&>(*snodes[i]);

    Deriv2ndMap& dd_debug = cnode.AugData().GetDeriv2nd_Debug();
    if (dd_debug.capacity() == 0)
    {
      CORE::GEN::reset(1, dd_debug);
    }

    // --- linearized normal multiplied by the varied slave position ----------
    for (unsigned k = 0; k < my::SLAVENUMNODE; ++k)
    {
      const Node& snode = static_cast<const Node&>(*snodes[k]);
      const int* sdof = snode.Dofs().data();

      // variation of the slave position multiplied with the linearized
      // smooth normal
      for (unsigned d = 0; d < probdim; ++d)
      {
        CORE::GEN::increaseCapacity(dd_debug);
        Deriv1stMap& dd_debug_var = dd_debug[sdof[d]];

        for (auto& d_n_unit_lin : d_n_unit[d])
        {
          CORE::GEN::increaseCapacity(dd_debug_var);
          dd_debug_var[d_n_unit_lin.first] += tmp * sval(k, 0) * d_n_unit_lin.second;
        }
      }
    }

    // --- varied slave gap multiplied by linearized jacobian  ----------------
    for (auto& d_gapn_sl_var : d_gapn_sl_complete)
    {
      CORE::GEN::increaseCapacity(dd_debug);
      Deriv1stMap& dd_debug_var = dd_debug[d_gapn_sl_var.first];

      const double val = wgt * lmval(i, 0) * d_gapn_sl_var.second;

      for (auto& d_jac_lin : d_jac)
      {
        CORE::GEN::increaseCapacity(dd_debug_var);
        dd_debug_var[d_jac_lin.first] += d_jac_lin.second * val;
      }
    }

    // --- varied jacobian multiplied by linearized slave gap  ----------------
    for (auto& d_jac_var : d_jac)
    {
      const double val = wgt * lmval(i, 0) * d_jac_var.second;

      // linearized slave part of the normal gap
      CORE::GEN::increaseCapacity(dd_debug);
      Deriv1stMap& dd_debug_var = dd_debug[d_jac_var.first];

      for (auto& d_gapn_sl_lin : d_gapn_sl_complete)
      {
        CORE::GEN::increaseCapacity(dd_debug_var);
        dd_debug_var[d_gapn_sl_lin.first] += d_gapn_sl_lin.second * val;
      }
    }

    // --- 2-nd order derivative of the jacobian ------------------------------
    const double tmp_sl = wgt * lmval(i, 0) * gapn_sl;
    for (auto& dd_jac_var : dd_jac)
    {
      const int gid_var = dd_jac_var.first;

      CORE::GEN::increaseCapacity(dd_debug);
      Deriv1stMap& dd_debug_var = dd_debug[gid_var];

      for (auto& dd_jac_var_lin : dd_jac_var.second)
      {
        const int gid_lin = dd_jac_var_lin.first;

        CORE::GEN::increaseCapacity(dd_debug_var);
        dd_debug_var[gid_lin] += tmp_sl * dd_jac_var_lin.second;
      }
    }


    // --- Additional parts for debugging only --------------------------------
    // slave position multiplied with the 2-nd order derivative of the smooth
    // unit normal
    CORE::LINALG::Matrix<probdim, 1> xs_int(true);
    for (unsigned k = 0; k < my::SLAVENUMNODE; ++k)
    {
      const Node& snode = static_cast<const Node&>(*snodes[k]);
      const double* xs = snode.xspatial();

      for (unsigned j = 0; j < probdim; ++j)
      {
        xs_int(j, 0) += sval(k, 0) * xs[j];
      }
    }

    for (unsigned j = 0; j < probdim; ++j)
    {
      const double val = tmp * xs_int(j, 0);

      for (auto& dd_sun_var : dd_n_unit[j])
      {
        CORE::GEN::increaseCapacity(dd_debug);
        Deriv1stMap& dd_debug_var = dd_debug[dd_sun_var.first];

        for (auto& dd_sun_var_lin : dd_sun_var.second)
        {
          CORE::GEN::increaseCapacity(dd_debug_var);
          dd_debug_var[dd_sun_var_lin.first] += val * dd_sun_var_lin.second;
        }
      }
    }

    // varied smooth unit normal multiplied with the linearized slave position
    for (unsigned j = 0; j < probdim; ++j)
    {
      for (auto& d_sun_j_var : d_n_unit[j])
      {
        CORE::GEN::increaseCapacity(dd_debug);
        Deriv1stMap& dd_debug_var = dd_debug[d_sun_j_var.first];

        for (unsigned k = 0; k < my::SLAVENUMNODE; ++k)
        {
          const Node& snode = static_cast<const Node&>(*snodes[k]);
          const int* sdof = snode.Dofs().data();

          CORE::GEN::increaseCapacity(dd_debug_var);
          dd_debug_var[sdof[j]] += sval(k, 0) * d_sun_j_var.second * tmp;
        }
      }
    }
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template <unsigned probdim, CORE::FE::CellType slavetype, CORE::FE::CellType mastertype>
void CONTACT::AUG::DebugCompleteIntPolicy<probdim, slavetype, mastertype>::Debug_Deriv1st_WGap(
    MORTAR::Element& sele) const
{
  DRT::Node* const* snodes = sele.Nodes();

  for (unsigned i = 0; i < my::SLAVENUMNODE; ++i)
  {
    Node& cnode = dynamic_cast<Node&>(*snodes[i]);

    Deriv1stMap& d_debug = cnode.AugData().GetDeriv1st_Debug();
    CORE::GEN::copy(cnode.AugData().GetDeriv1st_WGapMa(), d_debug);

    const Deriv1stMap& d_wgap_sl = cnode.AugData().GetDeriv1st_WGapSl();

    for (auto& d_wgap_sl_var : d_wgap_sl)
    {
      CORE::GEN::increaseCapacity(d_debug);
      d_debug[d_wgap_sl_var.first] -= d_wgap_sl_var.second;
    }
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template <unsigned probdim, CORE::FE::CellType slavetype, CORE::FE::CellType mastertype>
void CONTACT::AUG::DebugCompleteIntPolicy<probdim, slavetype, mastertype>::Debug_Deriv2nd_WGap(
    MORTAR::Element& sele) const
{
  DRT::Node* const* snodes = sele.Nodes();

  for (unsigned i = 0; i < my::SLAVENUMNODE; ++i)
  {
    Node& cnode = dynamic_cast<Node&>(*snodes[i]);

    Deriv2ndMap& dd_debug = cnode.AugData().GetDeriv2nd_Debug();
    CORE::GEN::copy(cnode.AugData().GetDeriv2nd_WGapMa(), dd_debug);

    const Deriv2ndMap& dd_wgap_sl = cnode.AugData().GetDeriv2nd_WGapSl();
    for (auto& dd_wgap_sl_var : dd_wgap_sl)
    {
      CORE::GEN::increaseCapacity(dd_debug);
      Deriv1stMap& dd_debug_var = dd_debug[dd_wgap_sl_var.first];

      for (auto& dd_wgap_sl_var_lin : dd_wgap_sl_var.second)
      {
        CORE::GEN::increaseCapacity(dd_debug_var);
        dd_debug_var[dd_wgap_sl_var_lin.first] -= dd_wgap_sl_var_lin.second;
      }
    }
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template <unsigned probdim, CORE::FE::CellType slavetype, CORE::FE::CellType mastertype>
void CONTACT::AUG::DebugCompleteIntPolicy<probdim, slavetype, mastertype>::Debug_MXi(
    MORTAR::Element& sele, const CORE::LINALG::Matrix<my::SLAVENUMNODE, 1>& lmval,
    const double* mxigp) const
{
  DRT::Node* const* snodes = sele.Nodes();

  for (unsigned i = 0; i < my::SLAVENUMNODE; ++i)
  {
    Node& cnode = dynamic_cast<Node&>(*snodes[i]);
    std::vector<std::pair<int, double>>& vals = cnode.AugData().Get_DebugVec();

    if (vals.size() < my::MASTERDIM)
      FOUR_C_THROW("DebugVec has the wrong dimension! ( size=%d)", vals.size());

    for (unsigned d = 0; d < my::MASTERDIM; ++d)
    {
      vals[d].first = cnode.Id();
      vals[d].second += lmval(i, 0) * mxigp[d];
    }
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template <unsigned probdim, CORE::FE::CellType slavetype, CORE::FE::CellType mastertype>
void CONTACT::AUG::DebugCompleteIntPolicy<probdim, slavetype, mastertype>::Debug_Deriv1st_MXi(
    MORTAR::Element& sele, const CORE::LINALG::Matrix<my::SLAVENUMNODE, 1>& lmval,
    const Deriv1stVecMap& dmxigp) const
{
  DRT::Node* const* snodes = sele.Nodes();

  for (unsigned i = 0; i < my::SLAVENUMNODE; ++i)
  {
    Node& cnode = dynamic_cast<Node&>(*snodes[i]);

    Deriv1stVecMap& d_debug = cnode.AugData().GetDeriv1st_DebugVec();
    if (d_debug.size() < my::MASTERDIM)
    {
      d_debug.resize(my::MASTERDIM);
      CORE::GEN::reset(1, d_debug);
    }

    for (unsigned d = 0; d < my::MASTERDIM; ++d)
    {
      for (auto& dmxigp_var : dmxigp[d])
      {
        CORE::GEN::increaseCapacity(d_debug[d]);
        d_debug[d][dmxigp_var.first] += lmval(i, 0) * dmxigp_var.second;
      }
    }
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template <unsigned probdim, CORE::FE::CellType slavetype, CORE::FE::CellType mastertype>
void CONTACT::AUG::DebugCompleteIntPolicy<probdim, slavetype, mastertype>::Debug_Deriv2nd_MXi(
    MORTAR::Element& sele, const CORE::LINALG::Matrix<my::SLAVENUMNODE, 1>& lmval,
    const Deriv2ndVecMap& ddmxigp) const
{
  DRT::Node* const* snodes = sele.Nodes();

  for (unsigned i = 0; i < my::SLAVENUMNODE; ++i)
  {
    Node& cnode = dynamic_cast<Node&>(*snodes[i]);

    Deriv2ndVecMap& dd_debug = cnode.AugData().GetDeriv2nd_DebugVec();
    if (dd_debug.size() < my::MASTERDIM)
    {
      dd_debug.resize(my::MASTERDIM);
      CORE::GEN::reset(1, dd_debug);
    }

    for (unsigned d = 0; d < my::MASTERDIM; ++d)
    {
      for (auto& ddmxigp_var : ddmxigp[d])
      {
        CORE::GEN::increaseCapacity(dd_debug[d]);
        Deriv1stMap& dd_debug_var = dd_debug[d][ddmxigp_var.first];

        for (auto& ddmxigp_var_lin : ddmxigp_var.second)
        {
          CORE::GEN::increaseCapacity(dd_debug_var);
          dd_debug_var[ddmxigp_var_lin.first] += ddmxigp_var_lin.second * lmval(i, 0);
        }
      }
    }
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template <unsigned probdim, CORE::FE::CellType slavetype, CORE::FE::CellType mastertype>
void CONTACT::AUG::DebugCompleteIntPolicy<probdim, slavetype,
    mastertype>::Debug_Deriv1st_Smooth_Unit_Normal(MORTAR::Element& sele) const
{
  DRT::Node* const* snodes = sele.Nodes();

  for (unsigned i = 0; i < my::SLAVENUMNODE; ++i)
  {
    Node& cnode = dynamic_cast<Node&>(*snodes[i]);

    Deriv1stVecMap& d_debug = cnode.AugData().GetDeriv1st_DebugVec();
    CORE::GEN::copy(cnode.AugData().GetDeriv1st_N(), d_debug);
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template <unsigned probdim, CORE::FE::CellType slavetype, CORE::FE::CellType mastertype>
void CONTACT::AUG::DebugCompleteIntPolicy<probdim, slavetype,
    mastertype>::Debug_Deriv2nd_Smooth_Unit_Normal(MORTAR::Element& sele) const
{
  DRT::Node* const* snodes = sele.Nodes();

  for (unsigned i = 0; i < my::SLAVENUMNODE; ++i)
  {
    Node& cnode = dynamic_cast<Node&>(*snodes[i]);

    Deriv2ndVecMap& dd_debug = cnode.AugData().GetDeriv2nd_DebugVec();
    CORE::GEN::copy(cnode.AugData().GetDeriv2nd_N(), dd_debug);
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template <unsigned probdim, CORE::FE::CellType slavetype, CORE::FE::CellType mastertype>
void CONTACT::AUG::DebugCompleteIntPolicy<probdim, slavetype, mastertype>::Debug_Deriv1st_Jac(
    MORTAR::Element& sele, const CORE::LINALG::Matrix<my::SLAVENUMNODE, 1>& lmval,
    const Deriv1stMap& djac) const
{
  DRT::Node* const* snodes = sele.Nodes();

  for (unsigned i = 0; i < my::SLAVENUMNODE; ++i)
  {
    Node& cnode = dynamic_cast<Node&>(*snodes[i]);

    Deriv1stMap& d_debug = cnode.AugData().GetDeriv1st_Debug();
    if (d_debug.capacity() == 0)
    {
      CORE::GEN::reset(1, d_debug);
    }

    for (auto& djac_var : djac)
    {
      CORE::GEN::increaseCapacity(d_debug);
      d_debug[djac_var.first] += lmval(i, 0) * djac_var.second;
    }
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template <unsigned probdim, CORE::FE::CellType slavetype, CORE::FE::CellType mastertype>
void CONTACT::AUG::DebugCompleteIntPolicy<probdim, slavetype, mastertype>::Debug_Deriv2nd_Jac(
    MORTAR::Element& sele, const CORE::LINALG::Matrix<my::SLAVENUMNODE, 1>& lmval,
    const Deriv2ndMap& ddjac) const
{
  DRT::Node* const* snodes = sele.Nodes();

  for (unsigned i = 0; i < my::SLAVENUMNODE; ++i)
  {
    Node& cnode = dynamic_cast<Node&>(*snodes[i]);

    Deriv2ndMap& dd_debug = cnode.AugData().GetDeriv2nd_Debug();
    if (dd_debug.capacity() == 0)
    {
      CORE::GEN::reset(1, dd_debug);
    }

    for (auto& ddjac_var : ddjac)
    {
      CORE::GEN::increaseCapacity(dd_debug);
      Deriv1stMap& dd_debug_var = dd_debug[ddjac_var.first];

      for (auto& ddjac_var_lin : ddjac_var.second)
      {
        CORE::GEN::increaseCapacity(dd_debug_var);
        dd_debug_var[ddjac_var_lin.first] += ddjac_var_lin.second * lmval(i, 0);
      }
    }
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template <unsigned probdim, CORE::FE::CellType slavetype, CORE::FE::CellType mastertype>
void CONTACT::AUG::DebugCompleteIntPolicy<probdim, slavetype, mastertype>::Debug_Deriv1st_Jacobian(
    MORTAR::Element& sele, const CORE::LINALG::Matrix<my::SLAVENUMNODE, 1>& lmval,
    const CORE::LINALG::Matrix<my::SLAVEDIM, my::SLAVENUMNODE>& sderiv,
    const CORE::LINALG::Matrix<3, 2>& stau) const
{
  DRT::Node* const* snodes = sele.Nodes();

  CORE::LINALG::Matrix<probdim, my::SLAVENUMNODE, int> nodal_dofs;
  CONTACT::INTEGRATOR::GetElementNodalDofs(sele, nodal_dofs);

  // evaluate the non-unit slave element normal and the inverse of its length
  CORE::LINALG::Matrix<probdim, 1> unit_normal;
  this->UnitSlaveElementNormal(sele, stau, unit_normal);

  /*--------------------------------------------------------------------------*/
  // non-unit normal vector: 1-st order derivative
  // 1-st int: vector index corresponds to the normal component index
  // 2-nd int: paired vector key corresponds to varied dof GID
  Deriv1stVecMap d_non_unit_normal(probdim);
  this->Deriv1st_NonUnitSlaveElementNormal(sele, nodal_dofs, sderiv, stau, d_non_unit_normal);

  /*--------------------------------------------------------------------------*/
  // jacobian determinant: 1-st order derivative
  // 1-st int: 1-st paired vector key corresponds to varied dof GID
  Deriv1stMap deriv1st_jac;
  this->Deriv1st_Jacobian(unit_normal, d_non_unit_normal, deriv1st_jac);

  for (unsigned i = 0; i < my::SLAVENUMNODE; ++i)
  {
    Node& cnode = dynamic_cast<Node&>(*snodes[i]);

    Deriv1stMap& d_debug = cnode.AugData().GetDeriv1st_Debug();
    if (d_debug.capacity() == 0) CORE::GEN::reset(1, d_debug);

    for (auto& deriv1st_jac_var : deriv1st_jac)
    {
      CORE::GEN::increaseCapacity(d_debug);
      d_debug[deriv1st_jac_var.first] += lmval(i, 0) * deriv1st_jac_var.second;
    }
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template <unsigned probdim, CORE::FE::CellType slavetype, CORE::FE::CellType mastertype>
void CONTACT::AUG::DebugCompleteIntPolicy<probdim, slavetype, mastertype>::Debug_Deriv2nd_Jacobian(
    MORTAR::Element& sele, const CORE::LINALG::Matrix<my::SLAVENUMNODE, 1>& lmval,
    const CORE::LINALG::Matrix<my::SLAVEDIM, my::SLAVENUMNODE>& sderiv,
    const CORE::LINALG::Matrix<3, 2>& stau) const
{
  DRT::Node* const* snodes = sele.Nodes();

  CORE::LINALG::Matrix<probdim, my::SLAVENUMNODE, int> nodal_dofs;
  CONTACT::INTEGRATOR::GetElementNodalDofs(sele, nodal_dofs);

  // evaluate the non-unit slave element normal and the inverse of its length
  CORE::LINALG::Matrix<probdim, 1> unit_normal;
  const double length_n_inv = this->UnitSlaveElementNormal(sele, stau, unit_normal);

  /*--------------------------------------------------------------------------*/
  // non-unit normal vector: 1-st order derivative
  // 1-st int: vector index corresponds to the normal component index
  // 2-nd int: paired vector key corresponds to varied dof GID
  Deriv1stVecMap d_non_unit_normal(probdim);
  this->Deriv1st_NonUnitSlaveElementNormal(sele, nodal_dofs, sderiv, stau, d_non_unit_normal);

  Deriv2ndMap deriv2nd_jac;
  this->Get_Deriv2nd_Jacobian(
      sele, nodal_dofs, sderiv, unit_normal, length_n_inv, d_non_unit_normal, deriv2nd_jac);

  for (unsigned i = 0; i < my::SLAVENUMNODE; ++i)
  {
    Node& cnode = dynamic_cast<Node&>(*snodes[i]);

    Deriv2ndMap& dd_debug = cnode.AugData().GetDeriv2nd_Debug();
    if (dd_debug.capacity() == 0)
    {
      CORE::GEN::reset(1, dd_debug);
    }

    for (auto& deriv2nd_jac_var : deriv2nd_jac)
    {
      CORE::GEN::increaseCapacity(dd_debug);
      Deriv1stMap& dd_debug_var = dd_debug[deriv2nd_jac_var.first];
      for (auto& deriv2nd_jac_var_lin : deriv2nd_jac_var.second)
      {
        CORE::GEN::increaseCapacity(dd_debug_var);
        dd_debug_var[deriv2nd_jac_var_lin.first] += lmval(i, 0) * deriv2nd_jac_var_lin.second;
      }
    }
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template <unsigned probdim, CORE::FE::CellType slavetype, CORE::FE::CellType mastertype>
void CONTACT::AUG::DebugCompleteIntPolicy<probdim, slavetype,
    mastertype>::Debug_Deriv1st_Non_Unit_Normal(MORTAR::Element& sele,
    const CORE::LINALG::Matrix<my::SLAVENUMNODE, 1>& lmval,
    const CORE::LINALG::Matrix<my::SLAVEDIM, my::SLAVENUMNODE>& sderiv,
    const CORE::LINALG::Matrix<3, 2>& stau) const
{
  DRT::Node* const* snodes = sele.Nodes();

  CORE::LINALG::Matrix<probdim, my::SLAVENUMNODE, int> nodal_dofs;
  CONTACT::INTEGRATOR::GetElementNodalDofs(sele, nodal_dofs);

  /*--------------------------------------------------------------------------*/
  // non-unit normal vector: 1-st order derivative
  // 1-st int: vector index corresponds to the normal component index
  // 2-nd int: paired vector key corresponds to varied dof GID
  Deriv1stVecMap d_non_unit_normal(probdim);
  this->Deriv1st_NonUnitSlaveElementNormal(sele, nodal_dofs, sderiv, stau, d_non_unit_normal);

  for (unsigned i = 0; i < my::SLAVENUMNODE; ++i)
  {
    Node& cnode = dynamic_cast<Node&>(*snodes[i]);

    Deriv1stVecMap& d_debug = cnode.AugData().GetDeriv1st_DebugVec();
    if (d_debug.size() < 3)
    {
      d_debug.resize(probdim);
      CORE::GEN::reset(1, d_debug);
    }

    for (unsigned d = 0; d < probdim; ++d)
    {
      for (auto& d_nun_var : d_non_unit_normal[d])
      {
        CORE::GEN::increaseCapacity(d_debug[d]);
        d_debug[d][d_nun_var.first] += lmval(i, 0) * d_nun_var.second;
      }
    }
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template <unsigned probdim, CORE::FE::CellType slavetype, CORE::FE::CellType mastertype>
void CONTACT::AUG::DebugCompleteIntPolicy<probdim, slavetype,
    mastertype>::Debug_Deriv2nd_Non_Unit_Normal(MORTAR::Element& sele,
    const CORE::LINALG::Matrix<my::SLAVENUMNODE, 1>& lmval,
    const CORE::LINALG::Matrix<my::SLAVEDIM, my::SLAVENUMNODE>& sderiv) const
{
  DRT::Node* const* snodes = sele.Nodes();

  CORE::LINALG::Matrix<probdim, my::SLAVENUMNODE, int> nodal_dofs;
  CONTACT::INTEGRATOR::GetElementNodalDofs(sele, nodal_dofs);

  Deriv2ndVecMap dd_non_unit_normal(probdim);
  this->Deriv2nd_NonUnitSlaveElementNormal(sele, nodal_dofs, sderiv, dd_non_unit_normal);

  for (unsigned i = 0; i < my::SLAVENUMNODE; ++i)
  {
    Node& cnode = dynamic_cast<Node&>(*snodes[i]);

    Deriv2ndVecMap& dd_debug = cnode.AugData().GetDeriv2nd_DebugVec();
    if (dd_debug.size() < 3)
    {
      dd_debug.resize(probdim);
      CORE::GEN::reset(1, dd_debug);
    }

    for (unsigned d = 0; d < probdim; ++d)
    {
      for (auto& dd_nun_var : dd_non_unit_normal[d])
      {
        CORE::GEN::increaseCapacity(dd_debug[d]);
        Deriv1stMap& dd_debug_var = dd_debug[d][dd_nun_var.first];

        for (auto& dd_nun_var_lin : dd_nun_var.second)
        {
          CORE::GEN::increaseCapacity(dd_debug_var);
          dd_debug_var[dd_nun_var_lin.first] += dd_nun_var_lin.second * lmval(i, 0);
        }
      }
    }
  }
}

/*----------------------------------------------------------------------------*/
template class CONTACT::AUG::DebugIncompleteIntPolicy<2, CORE::FE::CellType::line2,
    CORE::FE::CellType::line2>;
template class CONTACT::AUG::DebugIncompleteIntPolicy<2, CORE::FE::CellType::line2,
    CORE::FE::CellType::nurbs2>;
template class CONTACT::AUG::DebugIncompleteIntPolicy<2, CORE::FE::CellType::line2,
    CORE::FE::CellType::nurbs3>;

template class CONTACT::AUG::DebugIncompleteIntPolicy<2, CORE::FE::CellType::nurbs2,
    CORE::FE::CellType::nurbs2>;
template class CONTACT::AUG::DebugIncompleteIntPolicy<2, CORE::FE::CellType::nurbs2,
    CORE::FE::CellType::line2>;
template class CONTACT::AUG::DebugIncompleteIntPolicy<2, CORE::FE::CellType::nurbs2,
    CORE::FE::CellType::nurbs3>;

template class CONTACT::AUG::DebugIncompleteIntPolicy<2, CORE::FE::CellType::nurbs3,
    CORE::FE::CellType::nurbs3>;
template class CONTACT::AUG::DebugIncompleteIntPolicy<2, CORE::FE::CellType::nurbs3,
    CORE::FE::CellType::line2>;
template class CONTACT::AUG::DebugIncompleteIntPolicy<2, CORE::FE::CellType::nurbs3,
    CORE::FE::CellType::nurbs2>;

template class CONTACT::AUG::DebugIncompleteIntPolicy<3, CORE::FE::CellType::quad4,
    CORE::FE::CellType::quad4>;
template class CONTACT::AUG::DebugIncompleteIntPolicy<3, CORE::FE::CellType::quad4,
    CORE::FE::CellType::tri3>;
template class CONTACT::AUG::DebugIncompleteIntPolicy<3, CORE::FE::CellType::quad4,
    CORE::FE::CellType::nurbs4>;
template class CONTACT::AUG::DebugIncompleteIntPolicy<3, CORE::FE::CellType::quad4,
    CORE::FE::CellType::nurbs9>;

template class CONTACT::AUG::DebugIncompleteIntPolicy<3, CORE::FE::CellType::tri3,
    CORE::FE::CellType::quad4>;
template class CONTACT::AUG::DebugIncompleteIntPolicy<3, CORE::FE::CellType::tri3,
    CORE::FE::CellType::tri3>;
template class CONTACT::AUG::DebugIncompleteIntPolicy<3, CORE::FE::CellType::tri3,
    CORE::FE::CellType::nurbs4>;
template class CONTACT::AUG::DebugIncompleteIntPolicy<3, CORE::FE::CellType::tri3,
    CORE::FE::CellType::nurbs9>;

template class CONTACT::AUG::DebugIncompleteIntPolicy<3, CORE::FE::CellType::nurbs4,
    CORE::FE::CellType::nurbs4>;
template class CONTACT::AUG::DebugIncompleteIntPolicy<3, CORE::FE::CellType::nurbs4,
    CORE::FE::CellType::quad4>;
template class CONTACT::AUG::DebugIncompleteIntPolicy<3, CORE::FE::CellType::nurbs4,
    CORE::FE::CellType::tri3>;
template class CONTACT::AUG::DebugIncompleteIntPolicy<3, CORE::FE::CellType::nurbs4,
    CORE::FE::CellType::nurbs9>;

template class CONTACT::AUG::DebugIncompleteIntPolicy<3, CORE::FE::CellType::nurbs9,
    CORE::FE::CellType::nurbs9>;
template class CONTACT::AUG::DebugIncompleteIntPolicy<3, CORE::FE::CellType::nurbs9,
    CORE::FE::CellType::quad4>;
template class CONTACT::AUG::DebugIncompleteIntPolicy<3, CORE::FE::CellType::nurbs9,
    CORE::FE::CellType::tri3>;
template class CONTACT::AUG::DebugIncompleteIntPolicy<3, CORE::FE::CellType::nurbs9,
    CORE::FE::CellType::nurbs4>;

/*----------------------------------------------------------------------------*/
template class CONTACT::AUG::DebugCompleteIntPolicy<2, CORE::FE::CellType::line2,
    CORE::FE::CellType::line2>;
template class CONTACT::AUG::DebugCompleteIntPolicy<2, CORE::FE::CellType::line2,
    CORE::FE::CellType::nurbs2>;
template class CONTACT::AUG::DebugCompleteIntPolicy<2, CORE::FE::CellType::line2,
    CORE::FE::CellType::nurbs3>;

template class CONTACT::AUG::DebugCompleteIntPolicy<2, CORE::FE::CellType::nurbs2,
    CORE::FE::CellType::nurbs2>;
template class CONTACT::AUG::DebugCompleteIntPolicy<2, CORE::FE::CellType::nurbs2,
    CORE::FE::CellType::line2>;
template class CONTACT::AUG::DebugCompleteIntPolicy<2, CORE::FE::CellType::nurbs2,
    CORE::FE::CellType::nurbs3>;

template class CONTACT::AUG::DebugCompleteIntPolicy<2, CORE::FE::CellType::nurbs3,
    CORE::FE::CellType::nurbs3>;
template class CONTACT::AUG::DebugCompleteIntPolicy<2, CORE::FE::CellType::nurbs3,
    CORE::FE::CellType::line2>;
template class CONTACT::AUG::DebugCompleteIntPolicy<2, CORE::FE::CellType::nurbs3,
    CORE::FE::CellType::nurbs2>;

template class CONTACT::AUG::DebugCompleteIntPolicy<3, CORE::FE::CellType::quad4,
    CORE::FE::CellType::quad4>;
template class CONTACT::AUG::DebugCompleteIntPolicy<3, CORE::FE::CellType::quad4,
    CORE::FE::CellType::tri3>;
template class CONTACT::AUG::DebugCompleteIntPolicy<3, CORE::FE::CellType::quad4,
    CORE::FE::CellType::nurbs4>;
template class CONTACT::AUG::DebugCompleteIntPolicy<3, CORE::FE::CellType::quad4,
    CORE::FE::CellType::nurbs9>;

template class CONTACT::AUG::DebugCompleteIntPolicy<3, CORE::FE::CellType::tri3,
    CORE::FE::CellType::quad4>;
template class CONTACT::AUG::DebugCompleteIntPolicy<3, CORE::FE::CellType::tri3,
    CORE::FE::CellType::tri3>;
template class CONTACT::AUG::DebugCompleteIntPolicy<3, CORE::FE::CellType::tri3,
    CORE::FE::CellType::nurbs4>;
template class CONTACT::AUG::DebugCompleteIntPolicy<3, CORE::FE::CellType::tri3,
    CORE::FE::CellType::nurbs9>;

template class CONTACT::AUG::DebugCompleteIntPolicy<3, CORE::FE::CellType::nurbs4,
    CORE::FE::CellType::nurbs4>;
template class CONTACT::AUG::DebugCompleteIntPolicy<3, CORE::FE::CellType::nurbs4,
    CORE::FE::CellType::quad4>;
template class CONTACT::AUG::DebugCompleteIntPolicy<3, CORE::FE::CellType::nurbs4,
    CORE::FE::CellType::tri3>;
template class CONTACT::AUG::DebugCompleteIntPolicy<3, CORE::FE::CellType::nurbs4,
    CORE::FE::CellType::nurbs9>;

template class CONTACT::AUG::DebugCompleteIntPolicy<3, CORE::FE::CellType::nurbs9,
    CORE::FE::CellType::nurbs9>;
template class CONTACT::AUG::DebugCompleteIntPolicy<3, CORE::FE::CellType::nurbs9,
    CORE::FE::CellType::quad4>;
template class CONTACT::AUG::DebugCompleteIntPolicy<3, CORE::FE::CellType::nurbs9,
    CORE::FE::CellType::tri3>;
template class CONTACT::AUG::DebugCompleteIntPolicy<3, CORE::FE::CellType::nurbs9,
    CORE::FE::CellType::nurbs4>;

FOUR_C_NAMESPACE_CLOSE
