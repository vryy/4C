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
           unsigned slavedim,
           unsigned slavenumnode,
           unsigned masterdim,
           unsigned masternumnode >
void CONTACT::AUG::Integrator<probdim,slavetype,mastertype,slavedim,slavenumnode,
    masterdim,masternumnode>::GP_kappa(
    MORTAR::MortarElement& sele,
    const LINALG::Matrix<slavenumnode,1>& lmval,
    double wgt,
    double jac) const
{
  // Get slave nodes
  DRT::Node** snodes = sele.Nodes();
  dsassert( snodes, "ERROR: AugmentedIntegrator::GP_2D_kappa: Null pointer!" );

//  // number of nodes (slave)
//  int nrow = sele.NumNode();

  // add to node
  for (unsigned j=0;j<slavenumnode;++j)
  {
    CONTACT::CoNode* cnode = static_cast<CONTACT::CoNode*>(snodes[j]);

    double val = 0.0;
    val = lmval(j)*jac*wgt;

    // add current Gauss point's contribution kappaseg
    cnode->AddKappaValue(val);
  }

  return;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template < unsigned probdim,
           DRT::Element::DiscretizationType slavetype,
           DRT::Element::DiscretizationType mastertype,
           unsigned slavedim,
           unsigned slavenumnode,
           unsigned masterdim,
           unsigned masternumnode >
void CONTACT::AUG::Integrator<probdim,slavetype,mastertype,slavedim,slavenumnode,
    masterdim,masternumnode>::GP_2D_kappa_Lin(
    unsigned iter,
    MORTAR::MortarElement& sele,
    const LINALG::Matrix<slavenumnode,1>& lmval,
    const LINALG::Matrix<slavenumnode,slavedim>& lmderiv,
    double dsxideta,
    double dxdsxi,
    double dxdsxidsxi,
    double wgt,
    const GEN::pairedvector<int,double>& dsxigp,
    const GEN::pairedvector<int,double>& derivjac,
    const std::vector<GEN::pairedvector<int,double> >& ximaps)
{
  // Get slave nodes
  DRT::Node** snodes = sele.Nodes();
  dsassert( snodes, "ERROR: AugmentedIntegrator::GP_2D_kappa: Null pointer!" );

  // map iterator
  typedef GEN::pairedvector<int,double>::const_iterator CI;

  CONTACT::CoNode* cnode = static_cast<CONTACT::CoNode*>(snodes[iter]);
  std::map<int,double>& kappaLinMap = cnode->AugData().GetKappaLin();

  double fac = 0.0;

  // (0) Lin(LmSlave) - slave GP coordinates
  fac  = lmderiv(iter,0)*dxdsxi;
  // (1) Lin(dxdsxi) - slave GP coordinates
  fac += lmval(iter)*dxdsxidsxi;
  fac *= wgt*dsxideta;
  for (CI p=dsxigp.begin();p!=dsxigp.end();++p)
    kappaLinMap[p->first] += fac*(p->second);

  // (2) Lin(dsxideta) - segment end coordinates
  fac = wgt*lmval(iter)*dxdsxi;
  for (CI p=ximaps[0].begin();p!=ximaps[0].end();++p)
    kappaLinMap[p->first] -= 0.5*fac*(p->second);
  for (CI p=ximaps[1].begin();p!=ximaps[1].end();++p)
    kappaLinMap[p->first] += 0.5*fac*(p->second);

  // (3) Lin(dxdsxi) - slave GP Jacobian
  fac = wgt*lmval(iter)*dsxideta;
  for (CI p=derivjac.begin();p!=derivjac.end();++p)
    kappaLinMap[p->first] += fac*(p->second);


  return;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template < unsigned probdim,
           DRT::Element::DiscretizationType slavetype,
           DRT::Element::DiscretizationType mastertype,
           unsigned slavedim,
           unsigned slavenumnode,
           unsigned masterdim,
           unsigned masternumnode >
void CONTACT::AUG::Integrator<probdim,slavetype,mastertype,slavedim,slavenumnode,
    masterdim,masternumnode>::GP_2D_kappa_Ele_Lin(
    unsigned iter,
    MORTAR::MortarElement& sele,
    const LINALG::Matrix<slavenumnode,1>& lmval,
    const LINALG::Matrix<slavenumnode,slavedim>& lmderiv,
    double dxdsxi,
    double wgt,
    const GEN::pairedvector<int,double>& derivjac)
{
  // Get slave nodes
  DRT::Node** snodes = sele.Nodes();
  dsassert( snodes, "ERROR: AugmentedIntegrator::GP_2D_kappa: Null pointer!" );

  // map iterator
  typedef GEN::pairedvector<int,double>::const_iterator CI;

  CONTACT::CoNode* cnode = static_cast<CONTACT::CoNode*>(snodes[iter]);
  std::map<int,double>& kappaLinMap = cnode->AugData().GetKappaLin();

  // (0) Lin(LmSlave) - slave GP coordinates --> 0
  // (1) Lin(dxdsxi) - slave GP coordinates --> 0
  // (2) Lin(dsxideta) - segment end coordinates --> 0
  // (3) Lin(dxdsxi) - slave GP Jacobian
  const double fac = wgt*lmval(iter);
  for (CI p=derivjac.begin();p!=derivjac.end();++p)
    kappaLinMap[p->first] += fac*(p->second);

  return;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template < unsigned probdim,
           DRT::Element::DiscretizationType slavetype,
           DRT::Element::DiscretizationType mastertype,
           unsigned slavedim,
           unsigned slavenumnode,
           unsigned masterdim,
           unsigned masternumnode >
void CONTACT::AUG::Integrator<probdim,slavetype,mastertype,slavedim,slavenumnode,
    masterdim,masternumnode>::GP_3D_kappa_Lin(
    unsigned iter,
    MORTAR::MortarElement& sele,
    const LINALG::Matrix<slavenumnode,1>& lmval,
    const LINALG::Matrix<slavenumnode,slavedim>& lmderiv,
    double wgt,
    double jac,
    const std::vector<GEN::pairedvector<int,double> >& dsxigp,
    const GEN::pairedvector<int,double>& jacintcellmap)
{
  // get slave element nodes themselves
  DRT::Node** snodes = sele.Nodes();

  // map iterator
  typedef GEN::pairedvector<int,double>::const_iterator CI;

  CONTACT::CoNode* cnode = static_cast<CONTACT::CoNode*>(snodes[iter]);
  std::map<int,double>& kappaLinMap = cnode->AugData().GetKappaLin();
  double fac = 0.0;

  // (1) Lin(Phi) - dual shape functions
  // this vanishes here since there are no deformation-dependent dual functions

  // (2) Lin(LmSlave) - slave GP coordinates
  for ( unsigned i=0; i< dsxigp.size(); ++i )
  {
    fac = wgt*lmderiv(iter,i)*jac;
    for ( CI p=dsxigp[i].begin(); p!=dsxigp[i].end(); ++p )
      kappaLinMap[p->first] += fac*(p->second);
  }

  // (3) Lin(dsxideta) - intcell GP Jacobian
  fac = wgt*lmval(iter);
  for (CI p=jacintcellmap.begin();p!=jacintcellmap.end();++p)
    kappaLinMap[p->first] += fac*(p->second);

  return;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template < unsigned probdim,
           DRT::Element::DiscretizationType slavetype,
           DRT::Element::DiscretizationType mastertype,
           unsigned slavedim,
           unsigned slavenumnode,
           unsigned masterdim,
           unsigned masternumnode >
void CONTACT::AUG::Integrator<probdim,slavetype,mastertype,slavedim,slavenumnode,
    masterdim,masternumnode>::GP_Normal_DerivNormal(
    MORTAR::MortarElement& sele,
    const MORTAR::MortarElement& mele,
    const LINALG::Matrix<slavenumnode,1>& sval,
    const LINALG::Matrix<slavenumnode,slavedim>& sderiv,
    const std::vector<GEN::pairedvector<int,double> >& dsxigp,
    double* gpn,
    std::vector<GEN::pairedvector<int,double> >& dnmap_unit,
    int linsize)
{
  // map iterator
  typedef GEN::pairedvector<int,double>::const_iterator CI;

  // get slave element nodes
  DRT::Node** snodes = sele.Nodes();
  dsassert( snodes, "ERROR: GP_2D_Normal_DerivNormal: Null Pointer!" );

  // number of slave nodes
//  const unsigned nrow = slavenumnode;
//  const unsigned ncol = masternumnode;
//  const unsigned ndof = probdim;

  for (unsigned i=0;i<slavenumnode;++i)
  {
    MORTAR::MortarNode* mymrtnode = static_cast<MORTAR::MortarNode*> (snodes[i]);
    gpn[0]+=sval(i)*mymrtnode->MoData().n()[0];
    gpn[1]+=sval(i)*mymrtnode->MoData().n()[1];
    gpn[2]+=sval(i)*mymrtnode->MoData().n()[2];
  }

  // normalize interpolated GP normal back to length 1.0 !!!
  const double lengthn = sqrt(gpn[0]*gpn[0]+gpn[1]*gpn[1]+gpn[2]*gpn[2]);
  if (lengthn<1.0e-12)
    dserror("ERROR: IntegrateAndDerivSegment: Divide by zero!");

  for (int i=0;i<3;++i)
    gpn[i]/=lengthn;

  // ******************************
  // Linearization of the gp-normal
  // ******************************
  switch ( probdim )
  {
    // *** 2-D case ***********************************************
    case 2:
    {
      // build directional derivative of slave GP normal (non-unit)
      INTEGRATOR::ResetPairedVector( masternumnode*probdim+linsize, dmap_nxsl_gp_ );
      INTEGRATOR::ResetPairedVector( masternumnode*probdim+linsize, dmap_nysl_gp_ );

      for (unsigned i=0;i<slavenumnode;++i)
      {
        CoNode* snode = static_cast<CoNode*> (snodes[i]);

        GEN::pairedvector<int,double>& dmap_nxsl_i =
            static_cast<CoNode*>(snodes[i])->CoData().GetDerivN()[0];
        GEN::pairedvector<int,double>& dmap_nysl_i =
            static_cast<CoNode*>(snodes[i])->CoData().GetDerivN()[1];

        for (CI p=dmap_nxsl_i.begin();p!=dmap_nxsl_i.end();++p)
          dmap_nxsl_gp_[p->first] += sval(i)*(p->second);
        for (CI p=dmap_nysl_i.begin();p!=dmap_nysl_i.end();++p)
          dmap_nysl_gp_[p->first] += sval(i)*(p->second);

        for (CI p=dsxigp[0].begin();p!=dsxigp[0].end();++p)
        {
          const double valx =  sderiv(i,0)*snode->MoData().n()[0];
          dmap_nxsl_gp_[p->first] += valx*(p->second);
          const double valy =  sderiv(i,0)*snode->MoData().n()[1];
          dmap_nysl_gp_[p->first] += valy*(p->second);
        }
      }

      // build directional derivative of slave GP normal (unit)
      const double ll = lengthn*lengthn;
      const double sxsx = gpn[0]*gpn[0]*ll; // gpn is the unit normal --> multiplication with ll
      const double sxsy = gpn[0]*gpn[1]*ll; // to get the non-unit normal
      const double sysy = gpn[1]*gpn[1]*ll;

      for (CI p=dmap_nxsl_gp_.begin();p!=dmap_nxsl_gp_.end();++p)
      {
        dnmap_unit[0][p->first] += 1/lengthn*(p->second);
        dnmap_unit[0][p->first] -= 1/(lengthn*lengthn*lengthn)*sxsx*(p->second);
        dnmap_unit[1][p->first] -= 1/(lengthn*lengthn*lengthn)*sxsy*(p->second);
      }

      for (CI p=dmap_nysl_gp_.begin();p!=dmap_nysl_gp_.end();++p)
      {
        dnmap_unit[1][p->first] += 1/lengthn*(p->second);
        dnmap_unit[1][p->first] -= 1/(lengthn*lengthn*lengthn)*sysy*(p->second);
        dnmap_unit[0][p->first] -= 1/(lengthn*lengthn*lengthn)*sxsy*(p->second);
      }
      break;
    }
    // *** 3-D case ***********************************************
    case 3:
    {
      // build directional derivative of slave GP normal (non-unit)
      INTEGRATOR::ResetPairedVector( linsize, dmap_nxsl_gp_ );
      INTEGRATOR::ResetPairedVector( linsize, dmap_nysl_gp_ );
      INTEGRATOR::ResetPairedVector( linsize, dmap_nzsl_gp_ );

      for (unsigned i=0;i<slavenumnode;++i)
      {
        CoNode* cnode = static_cast<CoNode*> (snodes[i]);

        GEN::pairedvector<int,double>& dmap_nxsl_i = cnode->CoData().GetDerivN()[0];
        GEN::pairedvector<int,double>& dmap_nysl_i = cnode->CoData().GetDerivN()[1];
        GEN::pairedvector<int,double>& dmap_nzsl_i = cnode->CoData().GetDerivN()[2];

        for (CI p=dmap_nxsl_i.begin();p!=dmap_nxsl_i.end();++p)
          dmap_nxsl_gp_[p->first] += sval(i)*(p->second);
        for (CI p=dmap_nysl_i.begin();p!=dmap_nysl_i.end();++p)
          dmap_nysl_gp_[p->first] += sval(i)*(p->second);
        for (CI p=dmap_nzsl_i.begin();p!=dmap_nzsl_i.end();++p)
          dmap_nzsl_gp_[p->first] += sval(i)*(p->second);

        for (CI p=dsxigp[0].begin();p!=dsxigp[0].end();++p)
        {
          const double valx =  sderiv(i,0)*cnode->MoData().n()[0];
          dmap_nxsl_gp_[p->first] += valx*(p->second);
          const double valy =  sderiv(i,0)*cnode->MoData().n()[1];
          dmap_nysl_gp_[p->first] += valy*(p->second);
          const double valz =  sderiv(i,0)*cnode->MoData().n()[2];
          dmap_nzsl_gp_[p->first] += valz*(p->second);
        }

        for (CI p=dsxigp[1].begin();p!=dsxigp[1].end();++p)
        {
          const double valx =  sderiv(i,1)*cnode->MoData().n()[0];
          dmap_nxsl_gp_[p->first] += valx*(p->second);
          const double valy =  sderiv(i,1)*cnode->MoData().n()[1];
          dmap_nysl_gp_[p->first] += valy*(p->second);
          const double valz =  sderiv(i,1)*cnode->MoData().n()[2];
          dmap_nzsl_gp_[p->first] += valz*(p->second);
        }
      }

      double ll = lengthn*lengthn;
      double sxsx = gpn[0]*gpn[0]*ll;
      double sxsy = gpn[0]*gpn[1]*ll;
      double sxsz = gpn[0]*gpn[2]*ll;
      double sysy = gpn[1]*gpn[1]*ll;
      double sysz = gpn[1]*gpn[2]*ll;
      double szsz = gpn[2]*gpn[2]*ll;

      for (CI p=dmap_nxsl_gp_.begin();p!=dmap_nxsl_gp_.end();++p)
      {
        dnmap_unit[0][p->first] += 1/lengthn*(p->second);
        dnmap_unit[0][p->first] -= 1/(lengthn*lengthn*lengthn)*sxsx*(p->second);
        dnmap_unit[1][p->first] -= 1/(lengthn*lengthn*lengthn)*sxsy*(p->second);
        dnmap_unit[2][p->first] -= 1/(lengthn*lengthn*lengthn)*sxsz*(p->second);
      }

      for (CI p=dmap_nysl_gp_.begin();p!=dmap_nysl_gp_.end();++p)
      {
        dnmap_unit[1][p->first] += 1/lengthn*(p->second);
        dnmap_unit[1][p->first] -= 1/(lengthn*lengthn*lengthn)*sysy*(p->second);
        dnmap_unit[0][p->first] -= 1/(lengthn*lengthn*lengthn)*sxsy*(p->second);
        dnmap_unit[2][p->first] -= 1/(lengthn*lengthn*lengthn)*sysz*(p->second);
      }

      for (CI p=dmap_nzsl_gp_.begin();p!=dmap_nzsl_gp_.end();++p)
      {
        dnmap_unit[2][p->first] += 1/lengthn*(p->second);
        dnmap_unit[2][p->first] -= 1/(lengthn*lengthn*lengthn)*szsz*(p->second);
        dnmap_unit[0][p->first] -= 1/(lengthn*lengthn*lengthn)*sxsz*(p->second);
        dnmap_unit[1][p->first] -= 1/(lengthn*lengthn*lengthn)*sysz*(p->second);
      }
      break;
    }
    default:
      dserror( "Unsupported problem dimension of %d!", probdim );
      exit( EXIT_FAILURE );
  }

  return;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template < unsigned probdim,
           DRT::Element::DiscretizationType slavetype,
           DRT::Element::DiscretizationType mastertype,
           unsigned slavedim,
           unsigned slavenumnode,
           unsigned masterdim,
           unsigned masternumnode >
void CONTACT::AUG::Integrator<probdim,slavetype,mastertype,slavedim,slavenumnode,
    masterdim,masternumnode>::GP_VarWGap(
    MORTAR::MortarElement& sele,
    MORTAR::MortarElement& mele,
    const LINALG::Matrix<slavenumnode,1>& sval,
    const LINALG::Matrix<masternumnode,1>& mval,
    const LINALG::Matrix<slavenumnode,1>& lmval,
    const double* gpn,
    double wgt,
    double jac) const
{
  // get slave element nodes
  DRT::Node** snodes = sele.Nodes();
  DRT::Node** mnodes = mele.Nodes();
  dsassert( snodes, "ERROR: IntegrateAndDerivSegment: Null pointer!" );
  dsassert( mnodes, "ERROR: IntegrateAndDerivSegment: Null pointer!" );

  // number of nodes (slave, master)
//  const unsigned nrow = slavenumnode;
//  const unsigned ncol = masternumnode;

  // loop over all possible active slave nodes
  for (unsigned i=0;i<slavenumnode;++i)
  {
    CoNode* cnode = static_cast<CoNode*>(snodes[i]);

    // *** slave node contributions ***
    // loop over all slave nodes
    for (unsigned k=0;k<slavenumnode;++k)
    {
      CoNode* snode = static_cast<CoNode*>(snodes[k]);
      for ( unsigned kdof=0;kdof<probdim; ++kdof )
      {
        const int sGid   = snode->Id();
        const int sDofID = snode->Dofs()[kdof];
        const double val = lmval(i) * sval(k) * gpn[kdof] * wgt * jac;
        cnode->AddVarWGapSl(sDofID,sGid,val);
      }
    }

    // *** master node contributions ***
    // loop over all master nodes
    for (unsigned j=0;j<masternumnode;++j)
    {
      CoNode* mnode = static_cast<CoNode*>(mnodes[j]);
      for ( unsigned jdof=0; jdof<probdim; ++jdof )
      {
        const int mGid   = mnode->Id();
        const int mDofID = mnode->Dofs()[jdof];
        const double val = lmval(i) * mval(j) * gpn[jdof] * wgt * jac;
        cnode->AddVarWGapMa(mDofID,mGid,val);
      }
    }
  }

  return;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template < unsigned probdim,
           DRT::Element::DiscretizationType slavetype,
           DRT::Element::DiscretizationType mastertype,
           unsigned slavedim,
           unsigned slavenumnode,
           unsigned masterdim,
           unsigned masternumnode >
void CONTACT::AUG::Integrator<probdim,slavetype,mastertype,slavedim,slavenumnode,
    masterdim,masternumnode>::GP_2D_VarWGap_Lin(
    unsigned iter,
    MORTAR::MortarElement& sele,
    MORTAR::MortarElement& mele,
    const LINALG::Matrix<slavenumnode,1>& sval,
    const LINALG::Matrix<masternumnode,1>& mval,
    const LINALG::Matrix<slavenumnode,1>& lmval,
    const double* gpn,
    const LINALG::Matrix<slavenumnode,slavedim>& sderiv,
    const LINALG::Matrix<masternumnode,masterdim>& mderiv,
    const LINALG::Matrix<slavenumnode,slavedim>& lmderiv,
    double dsxideta,
    double dxdsxi,
    double dxdsxidsxi,
    double wgt,
    const GEN::pairedvector<int,double>& dsxigp,
    const GEN::pairedvector<int,double>& dmxigp,
    const GEN::pairedvector<int,double>& derivjac,
    const std::vector<GEN::pairedvector<int,double> >& dnmap_unit,
    const std::vector<GEN::pairedvector<int,double> >& ximaps) const
{
//  const unsigned nrow = slavenumnode;
//  const unsigned ncol = masternumnode;

  // get slave element nodes themselves
  DRT::Node** snodes = sele.Nodes();
  // get master element nodes themselves
  DRT::Node** mnodes = mele.Nodes();

  // map iterator
  typedef GEN::pairedvector<int,double>::const_iterator CI;

  CoNode* cnode = static_cast<CoNode*>(snodes[iter]);

  // *** integrate lin varWGapSl ****************************************
  for (unsigned k=0;k<slavenumnode;++k)
  {
    CoNode* snode = static_cast<CoNode*>(snodes[k]);

    double val[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

    // (0) Lin(LmSlave) - slave GP coordinates
    val[0] = wgt*lmderiv(iter,0)*sval(k)*dsxideta*dxdsxi;
    // (1) Lin(n-Slave) - normal direction
    val[1] = wgt*lmval(iter)*sval(k)*dsxideta*dxdsxi;
    // (2) Lin(NSlave) - slave GP coordinates
    val[2] = wgt*lmval(iter)*sderiv(k,0)*dsxideta*dxdsxi;
    // (3) Lin(dsxideta) - segment end coordinates
    val[3] = wgt*lmval(iter)*sval(k)*dxdsxi;
    // (4) Lin(dxdsxi) - slave GP Jacobian
    val[4] = wgt*lmval(iter)*sval(k)*dsxideta;
    // (5) Lin(dxdsxi) - slave GP coordinates
    val[5] = wgt*lmval(iter)*sval(k)*dsxideta*dxdsxidsxi;

    for (int kdof=0;kdof<snode->NumDof();++kdof)
    {
      int sDofId = snode->Dofs()[kdof];
      std::map<int,double>& varWGapLinSlMap = cnode->AugData().GetVarWGapLinSl()[sDofId];

      // (0,2,5) - slave GP coordinates
      for (CI p=dsxigp.begin();p!=dsxigp.end();++p)
        varWGapLinSlMap[p->first] += (val[0]+val[2]+val[5])*p->second*gpn[kdof];
      // (1) - normal direction
      for (CI p=dnmap_unit[kdof].begin();p!=dnmap_unit[kdof].end();++p)
        varWGapLinSlMap[p->first] += val[1]*p->second;
      // (3) - segment end coordinates
      for (CI p=ximaps[0].begin();p!=ximaps[0].end();++p)
        varWGapLinSlMap[p->first] -= 0.5*val[3]*p->second*gpn[kdof];
      for (CI p=ximaps[1].begin();p!=ximaps[1].end();++p)
        varWGapLinSlMap[p->first] += 0.5*val[3]*p->second*gpn[kdof];
      // (4) - Slave GP Jacobian
      for (CI p=derivjac.begin();p!=derivjac.end();++p)
        varWGapLinSlMap[p->first] += val[4]*p->second*gpn[kdof];
    }
  }
  // *** integrate lin varWGapMa ****************************************
  for (unsigned l=0;l<masternumnode;++l)
  {
    CoNode* mnode = static_cast<CoNode*>(mnodes[l]);

    double val[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

    // (0) Lin(LmSlave) - slave GP coordinates
    val[0] = wgt*lmderiv(iter,0)*mval(l)*dsxideta*dxdsxi;
    // (1) Lin(n-Slave) - normal direction
    val[1] = wgt*lmval(iter)*mval(l)*dsxideta*dxdsxi;
    // (2) Lin(NMaster) - master GP coordinates
    val[2] = wgt*lmval(iter)*mderiv(l,0)*dsxideta*dxdsxi;
    // (3) Lin(dsxideta) - segment end coordinates
    val[3] = wgt*lmval(iter)*mval(l)*dxdsxi;
    // (4) Lin(dxdsxi) - slave GP Jacobian
    val[4] = wgt*lmval(iter)*mval(l)*dsxideta;
    // (5) Lin(dxdsxi) - slave GP coordinates
    val[5] = wgt*lmval(iter)*mval(l)*dsxideta*dxdsxidsxi;

    for (int ldof=0;ldof<mnode->NumDof();++ldof)
    {
      int mDofId = mnode->Dofs()[ldof];
      std::map<int,double>& varWGapLinMaMap = cnode->AugData().GetVarWGapLinMa()[mDofId];

      // (0,5) - slave GP coordinates
      for (CI p=dsxigp.begin();p!=dsxigp.end();++p)
        varWGapLinMaMap[p->first] += (val[0]+val[5])*p->second*gpn[ldof];
      // (1) - normal direction
      for (CI p=dnmap_unit[ldof].begin();p!=dnmap_unit[ldof].end();++p)
        varWGapLinMaMap[p->first] += val[1]*p->second;
      // (2) - master GP coordinates
      for (CI p=dmxigp.begin();p!=dmxigp.end();++p)
        varWGapLinMaMap[p->first] += val[2]*p->second*gpn[ldof];
      // (3) - segement end coordinates
      for (CI p=ximaps[0].begin();p!=ximaps[0].end();++p)
        varWGapLinMaMap[p->first] -= 0.5*val[3]*p->second*gpn[ldof];
      for (CI p=ximaps[1].begin();p!=ximaps[1].end();++p)
        varWGapLinMaMap[p->first] += 0.5*val[3]*p->second*gpn[ldof];
      // (4) - slave GP Jacobian
      for (CI p=derivjac.begin();p!=derivjac.end();++p)
        varWGapLinMaMap[p->first] += val[4]*p->second*gpn[ldof];
    }
  }

  return;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template < unsigned probdim,
           DRT::Element::DiscretizationType slavetype,
           DRT::Element::DiscretizationType mastertype,
           unsigned slavedim,
           unsigned slavenumnode,
           unsigned masterdim,
           unsigned masternumnode >
void CONTACT::AUG::Integrator<probdim,slavetype,mastertype,slavedim,slavenumnode,
    masterdim,masternumnode>::GP_2D_VarWGap_Ele_Lin(
    unsigned iter,
    MORTAR::MortarElement& sele,
    MORTAR::MortarElement& mele,
    const LINALG::Matrix<slavenumnode,1>& sval,
    const LINALG::Matrix<masternumnode,1>& mval,
    const LINALG::Matrix<slavenumnode,1>& lmval,
    const double* gpn,
    const LINALG::Matrix<slavenumnode,slavedim>& sderiv,
    const LINALG::Matrix<masternumnode,masterdim>& mderiv,
    const LINALG::Matrix<slavenumnode,slavedim>& lmderiv,
    double dxdsxi,
    double wgt,
    const GEN::pairedvector<int,double>& dmxigp,
    const GEN::pairedvector<int,double>&derivjac,
    const std::vector<GEN::pairedvector<int,double> >& dnmap_unit) const
{
//  int nrow = sele.NumNode();
//  int ncol = mele.NumNode();

  // get slave element nodes themselves
  DRT::Node** snodes = sele.Nodes();
  // get master element nodes themselves
  DRT::Node** mnodes = mele.Nodes();

  // map iterator
  typedef GEN::pairedvector<int,double>::const_iterator CI;

  CoNode* cnode = static_cast<CoNode*>(snodes[iter]);

  // *** integrate lin varWGapSl ****************************************
  for (unsigned k=0;k<slavenumnode;++k)
  {
    CoNode* snode = static_cast<CoNode*>(snodes[k]);

    double val[2] = {0.0, 0.0};

    // (0) Lin(LmSlave) - slave GP coordinates --> 0

    // (1) Lin(n-Slave) - normal direction
    val[0] = wgt*lmval(iter)*sval(k)*dxdsxi;
    // (2) Lin(NSlave) - slave GP coordinates --> 0

    // (3) Lin(dsxideta) - segment end coordinates --> 0

    // (4) Lin(dxdsxi) - slave GP Jacobian
    val[1] = wgt*lmval(iter)*sval(k);
    // (5) Lin(dxdsxi) - slave GP coordinates --> 0

    for (int kdof=0;kdof<snode->NumDof();++kdof)
    {
      int sDofId = snode->Dofs()[kdof];
      std::map<int,double>& varWGapLinSlMap = cnode->AugData().GetVarWGapLinSl()[sDofId];

      // (1) - normal direction
      for (CI p=dnmap_unit[kdof].begin();p!=dnmap_unit[kdof].end();++p)
        varWGapLinSlMap[p->first] += val[0]*p->second;
      // (4) - Slave GP Jacobian
      for (CI p=derivjac.begin();p!=derivjac.end();++p)
        varWGapLinSlMap[p->first] += val[1]*p->second*gpn[kdof];
    }
  }
  // *** integrate lin varWGapMa ****************************************
  for (unsigned l=0;l<masternumnode;++l)
  {
    CoNode* mnode = static_cast<CoNode*>(mnodes[l]);

    double val[3] = {0.0, 0.0, 0.0};

    // (0) Lin(LmSlave) - slave GP coordinates --> 0

    // (1) Lin(n-Slave) - normal direction
    val[0] = wgt*lmval(iter)*mval(l)*dxdsxi;
    // (2) Lin(NMaster) - master GP coordinates
    val[1] = wgt*lmval(iter)*mderiv(l,0)*dxdsxi;
    // (3) Lin(dsxideta) - segment end coordinates --> 0

    // (4) Lin(dxdsxi) - slave GP Jacobian
    val[2] = wgt*lmval(iter)*mval(l);
    // (5) Lin(dxdsxi) - slave GP coordinates --> 0

    for (int ldof=0;ldof<mnode->NumDof();++ldof)
    {
      int mDofId = mnode->Dofs()[ldof];
      std::map<int,double>& varWGapLinMaMap = cnode->AugData().GetVarWGapLinMa()[mDofId];

      // (1) - normal direction
      for (CI p=dnmap_unit[ldof].begin();p!=dnmap_unit[ldof].end();++p)
        varWGapLinMaMap[p->first] += val[0]*p->second;
      // (2) - master GP coordinates
      for (CI p=dmxigp.begin();p!=dmxigp.end();++p)
        varWGapLinMaMap[p->first] += val[1]*p->second*gpn[ldof];
      // (4) - slave GP Jacobian
      for (CI p=derivjac.begin();p!=derivjac.end();++p)
        varWGapLinMaMap[p->first] += val[2]*p->second*gpn[ldof];
    }
  }

  return;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template < unsigned probdim,
           DRT::Element::DiscretizationType slavetype,
           DRT::Element::DiscretizationType mastertype,
           unsigned slavedim,
           unsigned slavenumnode,
           unsigned masterdim,
           unsigned masternumnode >
void CONTACT::AUG::Integrator<probdim,slavetype,mastertype,slavedim,slavenumnode,
    masterdim,masternumnode>::GP_3D_VarWGap_Lin(
    unsigned iter,
    MORTAR::MortarElement& sele,
    MORTAR::MortarElement& mele,
    const LINALG::Matrix<slavenumnode,1>& sval,
    const LINALG::Matrix<masternumnode,1>& mval,
    const LINALG::Matrix<slavenumnode,1>& lmval,
    const double* gpn,
    const LINALG::Matrix<slavenumnode,slavedim>& sderiv,
    const LINALG::Matrix<masternumnode,masterdim>& mderiv,
    const LINALG::Matrix<slavenumnode,slavedim>& lmderiv,
    double wgt,
    double jac,
    const std::vector<GEN::pairedvector<int,double> >& dsxigp,
    const std::vector<GEN::pairedvector<int,double> >& dmxigp,
    const GEN::pairedvector<int,double>& jacintcellmap,
    const std::vector<GEN::pairedvector<int,double> >& dnmap_unit) const
{
//  const unsigned nrow = slavenumnode;
//  const unsigned ncol = masternumnode;

  // get slave element nodes themselves
  DRT::Node** snodes = sele.Nodes();
  // get master element nodes themselves
  DRT::Node** mnodes = mele.Nodes();

  // map iterator
  typedef GEN::pairedvector<int,double>::const_iterator CI;

  CoNode* cnode = static_cast<CoNode*>(snodes[iter]);

  // *** integrate lin varWGapSl ****************************************
  for (unsigned k=0;k<slavenumnode;++k)
  {
    CoNode* snode = static_cast<CoNode*>(snodes[k]);
    double val[4] = {0.0,0.0,0.0,0.0};

    // (1) Lin(Phi) - dual shape functions
    // this vanishes here since there are no deformation-dependent dual functions
    // (2) Lin(LMShape) & Lin(NSlave) - 1st slave GP coordinate
    val[0]  = wgt*lmderiv(iter,0)*sval(k)*jac;
    val[0] += wgt*lmval(iter)*sderiv(k,0)*jac;
    // (3) Lin(LMShape) & Lin(NSlave) - 2nd slave GP coordinate
    val[1]  = wgt*lmderiv(iter,1)*sval(k)*jac;
    val[1] += wgt*lmval(iter)*sderiv(k,1)*jac;
    // (4) Lin(dsxideta) - intcell GP Jacobian
    val[2]  = wgt*lmval(iter)*sval(k);
    // (5) Lin(n-Slave) - normal direction
    val[3]  = wgt*lmval(iter)*sval(k)*jac;

    for (int kdof=0;kdof<snode->NumDof();++kdof)
    {
      int sDofId = snode->Dofs()[kdof];
      std::map<int,double>& varWGapLinSlMap = cnode->AugData().GetVarWGapLinSl()[sDofId];

      for (int i=0; i<(int) dsxigp.size();++i)
        for (CI p=dsxigp[i].begin(); p!=dsxigp[i].end(); ++p)
          varWGapLinSlMap[p->first] += val[i]*(p->second)*gpn[kdof];

      for (CI p=jacintcellmap.begin(); p!=jacintcellmap.end(); ++p)
        varWGapLinSlMap[p->first] += val[2]*(p->second)*gpn[kdof];

      for (CI p=dnmap_unit[kdof].begin();p!=dnmap_unit[kdof].end();++p)
        varWGapLinSlMap[p->first] += val[3]*(p->second);
    }
  }
  // *** integrate lin varWGapMa ****************************************
  for (unsigned l=0;l<masternumnode;++l)
  {
    CoNode* mnode = static_cast<CoNode*>(mnodes[l]);

    double val[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

    // (1) Lin(Phi) - dual shape functions
    // this vanishes here since there are no deformation-dependent dual functions
    // (2) Lin(LMShape) - 1st slave GP coordinate
    val[0] = wgt*lmderiv(iter,0)*mval(l)*jac;
    // (3) Lin(LMShape) - 2nd slave GP coordinate
    val[1] = wgt*lmderiv(iter,1)*mval(l)*jac;
    // (4) Lin(NMaster) - 1st master GP coordinate
    val[2] = wgt*lmval(iter)*mderiv(l,0)*jac;
    // (5) Lin(NMaster) - 1st master GP coordinate
    val[3] += wgt*lmval(iter)*mderiv(l,1)*jac;
    // (6) Lin(dsxideta) - intcell GP Jacobian
    val[4]  = wgt*lmval(iter)*mval(l);
    // (7) Lin(n-Slave) - normal direction
    val[5]  = wgt*lmval(iter)*mval(l)*jac;

    for (int ldof=0;ldof<mnode->NumDof();++ldof)
    {
      int mDofId = mnode->Dofs()[ldof];
      std::map<int,double>& varWGapLinMaMap = cnode->AugData().GetVarWGapLinMa()[mDofId];

      for (int i=0;i<(int) dsxigp.size();++i)
        for (CI p=dsxigp[i].begin(); p!=dsxigp[i].end(); ++p)
          varWGapLinMaMap[p->first] += val[i]*(p->second)*gpn[ldof];

      for (int i=0;i<(int) dmxigp.size();++i)
        for (CI p=dmxigp[i].begin(); p!=dmxigp[i].end(); ++p)
          varWGapLinMaMap[p->first] += val[i+2]*(p->second)*gpn[ldof];

      for (CI p=jacintcellmap.begin(); p!=jacintcellmap.end(); ++p)
        varWGapLinMaMap[p->first] += val[4]*(p->second)*gpn[ldof];

      for (CI p=dnmap_unit[ldof].begin();p!=dnmap_unit[ldof].end();++p)
        varWGapLinMaMap[p->first] += val[5]*(p->second);
    }
  }

  return;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template < unsigned probdim,
           DRT::Element::DiscretizationType slavetype,
           DRT::Element::DiscretizationType mastertype,
           unsigned slavedim,
           unsigned slavenumnode,
           unsigned masterdim,
           unsigned masternumnode >
void CONTACT::AUG::Integrator<probdim,slavetype,mastertype,slavedim,slavenumnode,
    masterdim,masternumnode>::GP_AugA(
    int it,
    MORTAR::MortarElement& sele,
    const LINALG::Matrix<slavenumnode,1>& lmval,
    double wgt,
    double jac) const
{
  // Get slave nodes
  DRT::Node** snodes = sele.Nodes();
  dsassert( snodes, "ERROR: AugmentedIntegrator::GP_2D_kappa: Null pointer!" );

  CoNode* cnode = static_cast<CoNode*>(snodes[it]);
  double& augA = cnode->AugData().GetAugA();

  augA += lmval(it)*jac*wgt;

  return;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template < unsigned probdim,
           DRT::Element::DiscretizationType slavetype,
           DRT::Element::DiscretizationType mastertype,
           unsigned slavedim,
           unsigned slavenumnode,
           unsigned masterdim,
           unsigned masternumnode >
void CONTACT::AUG::Integrator<probdim,slavetype,mastertype,slavedim,slavenumnode,
    masterdim,masternumnode>::GP_AugA_Lin(
    unsigned iter,
    MORTAR::MortarElement& sele,
    const LINALG::Matrix<slavenumnode,1>& lmval,
    double wgt,
    double jac,
    const GEN::pairedvector<int,double>& derivjac) const
{
  // Get slave nodes
  DRT::Node** snodes = sele.Nodes();
  dsassert( snodes, "ERROR: AugmentedIntegrator::GP_2D_kappa: Null pointer!" );

  // map iterator
  typedef GEN::pairedvector<int,double>::const_iterator CI;

  CONTACT::CoNode* cnode = static_cast<CONTACT::CoNode*>(snodes[iter]);
  GEN::pairedvector<int,double>& augALinMap = cnode->AugData().GetAugALin();

  // Lin(dxdsxi) - slave GP Jacobian
  double val = wgt*lmval(iter);

  // (2) - slave GP Jacobian
  for (CI p=derivjac.begin();p!=derivjac.end();++p)
    augALinMap[p->first] += val*p->second;

  return;
}

#include "contact_augmented_integrator_list.H"
