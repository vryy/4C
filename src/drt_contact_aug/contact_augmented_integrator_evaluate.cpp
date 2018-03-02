/*---------------------------------------------------------------------*/
/*!
\file contact_augmented_integrator_evaluate.cpp

\brief A class to perform integrations of Mortar matrices on the overlap
       of two MortarElements in 1D and 2D (derived version for
       augmented contact). This file contains only the evaluate routines.

\level 2

\maintainer Michael Hiermeier

\date Mar 8, 2017

*/
/*---------------------------------------------------------------------*/

#include "contact_augmented_integrator.H"
#include "../drt_contact/contact_element.H"
#include "../drt_contact/contact_node.H"

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template < unsigned probdim,
           DRT::Element::DiscretizationType slavetype,
           DRT::Element::DiscretizationType mastertype,
           class IntPolicy >
void CONTACT::AUG::Integrator<probdim,slavetype,mastertype,IntPolicy>::
EvaluatorDeriv1stOnly::Deriv_Jacobian(
    MORTAR::MortarElement& ele,
    const double* xi,
    const LINALG::Matrix<my::SLAVEDIM,my::SLAVENUMNODE>& sderiv,
    const LINALG::Matrix<3,2>& stau )
{
  LINALG::TMatrix<int,probdim,my::SLAVENUMNODE> nodal_dofs;

  // evaluate the non-unit slave element normal and the inverse of its length
  LINALG::Matrix<probdim,1> unit_normal;
  double length_n_inv = 0.0;

  Deriv1stVecMap d_non_unit_normal( probdim );

  Deriv1st_Jacobian( ele, xi, sderiv, stau, nodal_dofs, unit_normal,
      length_n_inv, d_non_unit_normal );
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template < unsigned probdim,
           DRT::Element::DiscretizationType slavetype,
           DRT::Element::DiscretizationType mastertype,
           class IntPolicy >
void CONTACT::AUG::Integrator<probdim,slavetype,mastertype,IntPolicy>::
EvaluatorDeriv1stOnly::Deriv1st_Jacobian(
    MORTAR::MortarElement& ele,
    const double* xi,
    const LINALG::Matrix<my::SLAVEDIM,my::SLAVENUMNODE>& sderiv,
    const LINALG::Matrix<3,2>& stau,
    LINALG::TMatrix<int,probdim,my::SLAVENUMNODE>& nodal_dofs,
    LINALG::Matrix<probdim,1>& unit_normal,
    double& length_n_inv,
    Deriv1stVecMap& d_non_unit_normal )
{
  CONTACT::INTEGRATOR::GetElementNodalDofs( ele, nodal_dofs );

  // evaluate the non-unit slave element normal and the inverse of its length
  length_n_inv = this->parent_.IntPolicy::UnitSlaveElementNormal( ele, stau, unit_normal);

  /*--------------------------------------------------------------------------*/
  // non-unit normal vector: 1-st order derivative
  // 1-st int: vector index corresponds to the normal component index
  // 2-nd int: paired vector key corresponds to varied dof GID
  this->parent_.IntPolicy::Deriv1st_NonUnitSlaveElementNormal( ele, nodal_dofs, sderiv,
      stau, d_non_unit_normal );

  /*--------------------------------------------------------------------------*/
  // jacobian determinant: 1-st order derivative
  // 1-st int: 1-st paired vector key corresponds to varied dof GID
  Deriv1stMap& deriv1st_jac = this->parent_.derivjac_;
  this->parent_.IntPolicy::Deriv1st_Jacobian( unit_normal, d_non_unit_normal, deriv1st_jac );
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template < unsigned probdim,
           DRT::Element::DiscretizationType slavetype,
           DRT::Element::DiscretizationType mastertype,
           class IntPolicy >
void CONTACT::AUG::Integrator<probdim,slavetype,mastertype,IntPolicy>::
EvaluatorFull::Deriv_Jacobian(
    MORTAR::MortarElement& ele,
    const double* xi,
    const LINALG::Matrix<my::SLAVEDIM,my::SLAVENUMNODE>& sderiv,
    const LINALG::Matrix<3,2>& stau )
{
  LINALG::TMatrix<int,probdim,my::SLAVENUMNODE> nodal_dofs;

  // evaluate the non-unit slave element normal and the inverse of its length
  LINALG::Matrix<probdim,1> unit_normal;
  double length_n_inv = 0.0;

  Deriv1stVecMap d_non_unit_normal( probdim );

  base_type::Deriv1st_Jacobian( ele, xi, sderiv, stau, nodal_dofs, unit_normal,
      length_n_inv, d_non_unit_normal );

  Deriv2ndMap& deriv2nd_jac = this->parent_.deriv2ndjac_;
  this->parent_.IntPolicy::Get_Deriv2nd_Jacobian( ele, nodal_dofs, sderiv, unit_normal,
      length_n_inv, d_non_unit_normal, deriv2nd_jac );
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template < unsigned probdim,
           DRT::Element::DiscretizationType slavetype,
           DRT::Element::DiscretizationType mastertype,
           class IntPolicy >
void CONTACT::AUG::Integrator<probdim,slavetype,mastertype,IntPolicy>::
EvaluatorDeriv1stOnly::Deriv_MXiGP(
    MORTAR::MortarElement& sele,
    MORTAR::MortarElement& mele,
    const double* sxi,
    const double* mxi,
    const double  alpha,
    const LINALG::Matrix<my::SLAVENUMNODE,1>& sval,
    const LINALG::Matrix<my::MASTERNUMNODE,1>& mval,
    const LINALG::Matrix<my::MASTERDIM,my::MASTERNUMNODE>& mderiv,
    const LINALG::Matrix<3,2>& mtau )
{
  LINALG::Matrix<probdim,probdim> lmat_inv( false );

  Deriv1st_MXiGP( sele, mele, sxi, mxi, alpha, sval, mval, mderiv, mtau,
      lmat_inv );
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template < unsigned probdim,
           DRT::Element::DiscretizationType slavetype,
           DRT::Element::DiscretizationType mastertype,
           class IntPolicy >
void CONTACT::AUG::Integrator<probdim,slavetype,mastertype,IntPolicy>::
EvaluatorDeriv1stOnly::Deriv1st_MXiGP(
    MORTAR::MortarElement& sele,
    MORTAR::MortarElement& mele,
    const double* sxi,
    const double* mxi,
    const double  alpha,
    const LINALG::Matrix<my::SLAVENUMNODE,1>& sval,
    const LINALG::Matrix<my::MASTERNUMNODE,1>& mval,
    const LINALG::Matrix<my::MASTERDIM,my::MASTERNUMNODE>& mderiv,
    const LINALG::Matrix<3,2>& mtau,
    LINALG::Matrix<probdim,probdim>& lmat_inv )
{
  LINALG::Matrix<probdim,1> snormal( false );
  this->parent_.IntPolicy::AveragedNormalAtXi( sele, sval, snormal );

  this->parent_.IntPolicy::LMatrixInverse( mtau, snormal, lmat_inv );

  Deriv1stVecMap& d_mxigp = this->parent_.dmxigp_;
  Deriv1stMap& d_alpha = this->parent_.dalpha_;
  this->parent_.IntPolicy::Deriv1st_MXiGP( lmat_inv, sele, mele, sval,
      mval, alpha, d_mxigp, d_alpha );
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template < unsigned probdim,
           DRT::Element::DiscretizationType slavetype,
           DRT::Element::DiscretizationType mastertype,
           class IntPolicy >
void CONTACT::AUG::Integrator<probdim,slavetype,mastertype,IntPolicy>::
EvaluatorFull::Deriv_MXiGP(
    MORTAR::MortarElement& sele,
    MORTAR::MortarElement& mele,
    const double* sxi,
    const double* mxi,
    const double  alpha,
    const LINALG::Matrix<my::SLAVENUMNODE,1>& sval,
    const LINALG::Matrix<my::MASTERNUMNODE,1>& mval,
    const LINALG::Matrix<my::MASTERDIM,my::MASTERNUMNODE>& mderiv,
    const LINALG::Matrix<3,2>& mtau )
{
  LINALG::Matrix<probdim,probdim> lmat_inv( false );

  base_type::Deriv1st_MXiGP( sele, mele, sxi, mxi, alpha, sval, mval, mderiv,
      mtau, lmat_inv );

  Deriv1stVecMap& d_mxigp = this->parent_.dmxigp_;
  Deriv1stMap& d_alpha = this->parent_.dalpha_;

  const LINALG::Matrix<3,my::MASTERNUMNODE>& mderiv2 = this->parent_.mderiv2nd_;
  Deriv2ndVecMap& dd_mxigp = this->parent_.ddmxigp_;

  this->parent_.IntPolicy::Get_Deriv2nd_MXiGP( lmat_inv, sele, mele, sval, mval,
      mderiv, mderiv2, mtau, mxi, alpha, d_mxigp, d_alpha, dd_mxigp );
}


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template < unsigned probdim,
           DRT::Element::DiscretizationType slavetype,
           DRT::Element::DiscretizationType mastertype,
           class IntPolicy >
void CONTACT::AUG::Integrator<probdim,slavetype,mastertype,IntPolicy>::
EvaluatorFull::Get_Deriv2nd_AugA(
        MORTAR::MortarElement& sele,
        const LINALG::Matrix<my::SLAVENUMNODE,1>& lmval,
        const double wgt,
        const Deriv2ndMap& dd_jac ) const
{
  DRT::Node* const * snodes = sele.Nodes();

  for ( unsigned i=0; i<my::SLAVENUMNODE; ++i )
  {
    CoNode& cnode = static_cast<CoNode&>( *snodes[i] );
    Deriv2ndMap& dd_a = cnode.AugData().GetDeriv2nd_A();

    const double tmp = wgt * lmval(i,0);

    for ( auto& dd_jac_var : dd_jac )
    {
      Deriv1stMap& dd_a_var = dd_a[ dd_jac_var.first ];

      for ( auto& dd_jac_var_lin : dd_jac_var.second )
      {
        dd_a_var( dd_jac_var_lin.first ) += tmp * dd_jac_var_lin.second;
      }
    }

  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template < unsigned probdim,
           DRT::Element::DiscretizationType slavetype,
           DRT::Element::DiscretizationType mastertype,
           class IntPolicy >
void CONTACT::AUG::Integrator<probdim,slavetype,mastertype,IntPolicy>::GP_kappa(
    MORTAR::MortarElement& sele,
    const LINALG::Matrix<my::SLAVENUMNODE,1>& lmval,
    double wgt,
    double jac) const
{
  // Get slave nodes
  DRT::Node** snodes = sele.Nodes();
  dsassert( snodes, "ERROR: AugmentedIntegrator::GP_2D_kappa: Null pointer!" );

//  // number of nodes (slave)
//  int nrow = sele.NumNode();

  // add to node
  double val = 0.0;
  for (unsigned j=0;j<my::SLAVENUMNODE;++j)
  {
    CONTACT::CoNode* cnode = static_cast<CONTACT::CoNode*>(snodes[j]);

    val = lmval(j,0) * jac * wgt;

    // add current Gauss point's contribution kappaseg
    cnode->AddKappaValue( val );
  }

  return;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template < unsigned probdim,
           DRT::Element::DiscretizationType slavetype,
           DRT::Element::DiscretizationType mastertype,
           class IntPolicy >
void CONTACT::AUG::Integrator<probdim,slavetype,mastertype,IntPolicy>::Get_Deriv1st_Kappa(
    MORTAR::MortarElement& sele,
    const LINALG::Matrix<my::SLAVENUMNODE,1>& lmval,
    double wgt,
    const Deriv1stMap& d_jac )
{
  // get slave element nodes themselves
  DRT::Node* const * snodes = sele.Nodes();

  for ( unsigned i=0; i<my::SLAVENUMNODE; ++i )
  {
    CoNode& cnode = static_cast<CONTACT::CoNode&>( *snodes[i] );
    Deriv1stMap& d_kappa = cnode.AugData().GetDeriv1st_Kappa();

    for ( auto& d_jac_var : d_jac )
    {
      d_kappa( d_jac_var.first ) += d_jac_var.second * wgt * lmval(i,0);
    }
  }

  return;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template < unsigned probdim,
           DRT::Element::DiscretizationType slavetype,
           DRT::Element::DiscretizationType mastertype,
           class IntPolicy >
void CONTACT::AUG::Integrator<probdim,slavetype,mastertype,IntPolicy>::GP_Normal_DerivNormal(
    MORTAR::MortarElement& sele,
    const LINALG::Matrix<my::SLAVENUMNODE,1>& sval,
    double* gpn,
    Deriv1stVecMap& dn_non_unit,
    Deriv2ndVecMap& ddn_non_unit,
    Deriv1stVecMap& dn_unit,
    Deriv2ndVecMap& ddn_unit )
{
  const DRT::Node* const* snodes = sele.Nodes();

  for ( unsigned i=0; i<my::SLAVENUMNODE; ++i )
  {
    const CoNode& snode = static_cast<const CoNode&>( *snodes[i] );
    const double* nodal_smooth_unit_normal = snode.MoData().n();

    const Deriv1stVecMap& nodal_d_sun = snode.AugData().GetDeriv1st_N();
    const Deriv2ndVecMap& nodal_dd_sun = snode.AugData().GetDeriv2nd_N();

    for ( unsigned j=0; j<probdim; ++j )
    {
      gpn[j] += sval(i,0) * nodal_smooth_unit_normal[j];

      // build 1-st derivative of the non-unit gauss point normal by
      // interpolating the 1-st derivatives of the unit smooth nodal normals
      Deriv1stMap& d_nun_j = dn_non_unit[j];
      for ( auto& nodal_d_sun_var : nodal_d_sun[j] )
        d_nun_j.repetitive_access( nodal_d_sun_var.first, my::gp_id_ ) +=
            nodal_d_sun_var.second * sval(i,0);

      // build 2-nd order derivative of the non-unit gauss point normal
      Deriv2ndMap& dd_nun_j = ddn_non_unit[j];
      for ( auto& nodal_dd_sun_var : nodal_dd_sun[j] )
      {
        Deriv1stMap& dd_nun_j_var =
            dd_nun_j.repetitive_access( nodal_dd_sun_var.first, my::gp_id_ );

        for ( auto& nodal_dd_sun_var_lin : nodal_dd_sun_var.second )
        {
          dd_nun_j_var.repetitive_access( nodal_dd_sun_var_lin.first, my::gp_id_ ) +=
              nodal_dd_sun_var_lin.second * sval(i,0);
        }
      }
    }
  }

  GEN_DATA::complete( dn_non_unit );
  GEN_DATA::complete( ddn_non_unit );

  // unify the gauss point normal
  LINALG::Matrix<probdim,1> unit_gp_normal( gpn, true );
  const double length_n_inv = 1.0 / unit_gp_normal.Norm2();
  unit_gp_normal.Scale( length_n_inv );

  // calculate the 1-st derivative of the unified smooth gauss point normal
  IntPolicy::Deriv1st_UnitSlaveElementNormal( unit_gp_normal, length_n_inv,
      dn_non_unit, dn_unit, false );

  // calculate the 2-nd order derivative of the unified smooth gauss point normal
  IntPolicy::Deriv2nd_UnitSlaveElementNormal( unit_gp_normal, length_n_inv,
      dn_non_unit, dn_unit, ddn_non_unit, ddn_unit );
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template < unsigned probdim,
           DRT::Element::DiscretizationType slavetype,
           DRT::Element::DiscretizationType mastertype,
           class IntPolicy >
void CONTACT::AUG::Integrator<probdim,slavetype,mastertype,IntPolicy>::GP_AugA(
    MORTAR::MortarElement& sele,
    const LINALG::Matrix<my::SLAVENUMNODE,1>& lmval,
    double wgt,
    double jac) const
{
  // Get slave nodes
  DRT::Node** snodes = sele.Nodes();
  dsassert( snodes, "ERROR: AugmentedIntegrator::GP_2D_kappa: Null pointer!" );

  for ( unsigned it=0; it<my::SLAVENUMNODE; ++it )
  {
    CoNode* cnode = static_cast<CoNode*>(snodes[it]);
    double& augA = cnode->AugData().GetAugA();

    augA += lmval(it,0)*jac*wgt;
  }
  return;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template < unsigned probdim,
           DRT::Element::DiscretizationType slavetype,
           DRT::Element::DiscretizationType mastertype,
           class IntPolicy >
void CONTACT::AUG::Integrator<probdim,slavetype,mastertype,IntPolicy>::Get_Deriv1st_AugA(
    MORTAR::MortarElement& sele,
    const LINALG::Matrix<my::SLAVENUMNODE,1>& lmval,
    double wgt,
    double jac,
    const Deriv1stMap& derivjac) const
{
  // Get slave nodes
  DRT::Node** snodes = sele.Nodes();
  dsassert( snodes, "ERROR: AugmentedIntegrator::GP_2D_kappa: Null pointer!" );

  // map iterator
  typedef Deriv1stMap::const_iterator CI;

  for ( unsigned i=0; i<my::SLAVENUMNODE; ++i )
  {
    CONTACT::CoNode* cnode = static_cast<CONTACT::CoNode*>(snodes[i]);
    Deriv1stMap& augALinMap = cnode->AugData().GetDeriv1st_A();

    double val = wgt*lmval(i);

    // slave GP Jacobian
    for (CI p=derivjac.begin();p!=derivjac.end();++p)
      augALinMap[p->first] += val*p->second;
  }

  return;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template < unsigned probdim,
           DRT::Element::DiscretizationType slavetype,
           DRT::Element::DiscretizationType mastertype,
           class IntPolicy >
void CONTACT::AUG::Integrator<probdim,slavetype,mastertype,IntPolicy>::Get_Deriv2nd_Kappa(
    MORTAR::MortarElement& sele,
    const LINALG::Matrix<my::SLAVENUMNODE,1>& lmval,
    const double wgt,
    const Deriv2ndMap& dd_jac ) const
{
  DRT::Node* const * snodes = sele.Nodes();

  for ( unsigned i=0; i<my::SLAVENUMNODE; ++i )
  {
    CoNode& cnode = static_cast<CoNode&>( *snodes[i] );
    Deriv2ndMap& dd_kappa = cnode.AugData().GetDeriv2nd_Kappa();

    const double tmp = wgt * lmval(i,0);

    for ( auto& dd_jac_var : dd_jac )
    {
      Deriv1stMap& dd_kappa_var = dd_kappa[ dd_jac_var.first ];

      for ( auto& dd_jac_var_lin : dd_jac_var.second )
      {
        dd_kappa_var( dd_jac_var_lin.first ) += tmp * dd_jac_var_lin.second;
      }
    }

  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template < unsigned probdim,
           DRT::Element::DiscretizationType slavetype,
           DRT::Element::DiscretizationType mastertype,
           class IntPolicy >
void CONTACT::AUG::Integrator<probdim,slavetype,mastertype,IntPolicy>::GapN(
    MORTAR::MortarElement& sele,
    MORTAR::MortarElement& mele,
    const LINALG::Matrix<my::SLAVENUMNODE,1>& sval,
    const LINALG::Matrix<my::MASTERNUMNODE,1>& mval,
    const double* gpn,
    double & gapn_sl,
    double & gapn_ma ) const
{
  const DRT::Node* const * snodes = sele.Nodes();
  const DRT::Node* const * mnodes = mele.Nodes();

  // slave contribution of the gauss point normal gap
  gapn_sl = 0.0;
  for ( unsigned i=0; i<my::SLAVENUMNODE; ++i )
  {
    const CoNode& snode = static_cast<const CoNode&>( *snodes[i] );
    const double* xs = snode.xspatial();

    for ( unsigned k=0; k<probdim; ++k )
      gapn_sl += sval(i,0) * gpn[k] * xs[k];
  }

  // master contribution of the gauss point normal gap
  gapn_ma = 0.0;
  for ( unsigned i=0; i<my::MASTERNUMNODE; ++i )
  {
    const CoNode& mnode = static_cast<const CoNode&>( *mnodes[i] );
    const double* xm = mnode.xspatial();

    for ( unsigned k=0; k<probdim; ++k )
      gapn_ma += mval(i,0) * gpn[k] * xm[k];
  }
}
/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template < unsigned probdim,
           DRT::Element::DiscretizationType slavetype,
           DRT::Element::DiscretizationType mastertype,
           class IntPolicy >
void CONTACT::AUG::Integrator<probdim,slavetype,mastertype,IntPolicy>::GP_WGap(
    MORTAR::MortarElement& sele,
    const LINALG::Matrix<my::SLAVENUMNODE,1>& lmval,
    const double gapn_sl,
    const double gapn_ma,
    const double wg,
    const double jac ) const
{
  DRT::Node* const * snodes = sele.Nodes();
  const double scale = wg * jac;

  for ( unsigned i=0; i<my::SLAVENUMNODE; ++i )
  {
    CoNode& cnode = static_cast<CoNode&>( *snodes[i] );

    const double val = lmval(i,0) * scale * ( gapn_ma - gapn_sl );

    cnode.AddWGapValue( val );
  }
}


#include "contact_augmented_integrator_list.H"
