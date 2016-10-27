/*---------------------------------------------------------------------*/
/*!
\file contact_integrator_nitsche.cpp

\brief A class to perform integrations of nitsche related terms

\level 3

\maintainer Alexander Seitz

*/
/*---------------------------------------------------------------------*/
#include <Teuchos_StandardParameterEntryValidators.hpp>
#include "contact_integrator_nitsche.H"
#include "contact_node.H"
#include "contact_element.H"
#include "contact_defines.H"
#include "contact_paramsinterface.H"
#include "../drt_mortar/mortar_defines.H"
#include "../drt_inpar/inpar_contact.H"

#include "../drt_fem_general/drt_utils_boundary_integration.H"

#include "../drt_so3/so_base.H"

#include "../drt_mat/elasthyper.H"
#include <Epetra_FEVector.h>
#include <Epetra_CrsMatrix.h>
#include "../linalg/linalg_utils.H"



/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CONTACT::CoIntegratorNitsche::IntegrateGP_3D(
    MORTAR::MortarElement& sele,
    MORTAR::MortarElement& mele,
    LINALG::SerialDenseVector& sval,
    LINALG::SerialDenseVector& lmval,
    LINALG::SerialDenseVector& mval,
    LINALG::SerialDenseMatrix& sderiv,
    LINALG::SerialDenseMatrix& mderiv,
    LINALG::SerialDenseMatrix& lmderiv,
    GEN::pairedvector<int,Epetra_SerialDenseMatrix>& dualmap,
    double& wgt,
    double& jac,
    GEN::pairedvector<int, double>& derivjac,
    double* normal,
    std::vector<GEN::pairedvector<int, double> >& dnmap_unit,
    double& gap,
    GEN::pairedvector<int, double>& deriv_gap,
    double* sxi,
    double* mxi,
    std::vector<GEN::pairedvector<int,double> >& derivsxi,
    std::vector<GEN::pairedvector<int,double> >& derivmxi
    )
{
    GPTS_forces<3>(sele,mele,sval,sderiv,derivsxi,mval,mderiv,derivmxi,
        jac,derivjac,wgt,gap,deriv_gap,normal,dnmap_unit,sxi,mxi);

  return;
}



/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CONTACT::CoIntegratorNitsche::IntegrateGP_2D(
    MORTAR::MortarElement& sele,
    MORTAR::MortarElement& mele,
    LINALG::SerialDenseVector& sval,
    LINALG::SerialDenseVector& lmval,
    LINALG::SerialDenseVector& mval,
    LINALG::SerialDenseMatrix& sderiv,
    LINALG::SerialDenseMatrix& mderiv,
    LINALG::SerialDenseMatrix& lmderiv,
    GEN::pairedvector<int,Epetra_SerialDenseMatrix>& dualmap,
    double& wgt,
    double& jac,
    GEN::pairedvector<int, double>& derivjac,
    double* normal,
    std::vector<GEN::pairedvector<int, double> >& dnmap_unit,
    double& gap,
    GEN::pairedvector<int, double>& deriv_gap,
    double* sxi,
    double* mxi,
    std::vector<GEN::pairedvector<int,double> >& derivsxi,
    std::vector<GEN::pairedvector<int,double> >& derivmxi
    )
{
    GPTS_forces<2>(sele,mele,sval,sderiv,derivsxi,mval,mderiv,derivmxi,
        jac,derivjac,wgt,gap,deriv_gap,normal,dnmap_unit,sxi,mxi);

  return;
}



/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <int dim>
void CONTACT::CoIntegratorNitsche::GPTS_forces(
    MORTAR::MortarElement& sele, MORTAR::MortarElement& mele,
    const LINALG::SerialDenseVector& sval, const LINALG::SerialDenseMatrix& sderiv,const std::vector<GEN::pairedvector<int,double> >& dsxi,
    const LINALG::SerialDenseVector& mval, const LINALG::SerialDenseMatrix& mderiv,const std::vector<GEN::pairedvector<int,double> >& dmxi,
    const double jac,const GEN::pairedvector<int,double>& jacintcellmap, const double wgt,
    const double gap, const GEN::pairedvector<int,double>& dgapgp,
    double* gpn, std::vector<GEN::pairedvector<int,double> >& dnmap_unit,double* sxi,double* mxi)
{
  if (sele.Owner()!=Comm_.MyPID())
    return;

  if (fc_==Teuchos::null || kc_==Teuchos::null)
    dserror("matrix or vector not provided for Nitsche contact integration");

  if (dim!=Dim())
    dserror("dimension inconsistency");

  double pen = ppn_;

  double gap_plus_cauchy_nn=gap;
  std::map<int,double> gap_plus_cauchy_nn_deriv;
  for (GEN::pairedvector<int,double>::const_iterator p=dgapgp.begin();p!=dgapgp.end();++p)
    gap_plus_cauchy_nn_deriv[p->first]+=p->second;
  const LINALG::Matrix<dim,1> normal(gpn,true);

  if (stype_==INPAR::CONTACT::solution_nitsche)
  {
    pen*=std::max(dynamic_cast<CONTACT::CoElement&>(sele).TraceH(),dynamic_cast<CONTACT::CoElement&>(mele).TraceH());
    MAT::ElastHyper* mmat = dynamic_cast<MAT::ElastHyper*>(&*mele.ParentElement()->Material());
    MAT::ElastHyper* smat = dynamic_cast<MAT::ElastHyper*>(&*sele.ParentElement()->Material());
    if (smat==NULL || mmat==NULL)
      dserror("Nitsche contact only for elast hyper material");
    double ws=1.*dynamic_cast<CONTACT::CoElement&>(mele).TraceH()
        *mmat->GetYoung();
    double wm=1.*dynamic_cast<CONTACT::CoElement&>(sele).TraceH()
        *smat->GetYoung();
    ws/=(ws+wm);
    wm=1.-ws;

//    ws=1.;wm=0.; // fixme
    SoEleCauchy<dim>(sele,sxi,dsxi,wgt,gpn,dnmap_unit,false,ws/pen,gap_plus_cauchy_nn,gap_plus_cauchy_nn_deriv);
    SoEleCauchy<dim>(mele,mxi,dmxi,wgt,gpn,dnmap_unit,false,wm/pen,gap_plus_cauchy_nn,gap_plus_cauchy_nn_deriv);
  }
  else if (stype_==INPAR::CONTACT::solution_penalty)
    {} // do nothing
else
    dserror("unknown algorithm");

  if (gap_plus_cauchy_nn>=0.)
  {

    return;
  }
double val=0.;
  for (int d=0;d<Dim();++d)
  {
    double gp_force = pen * jac*wgt*gap_plus_cauchy_nn*gpn[d];

    for (int s=0;s<sele.NumNode();++s)
    {
      val = gp_force*sval(s);
      fc_->SumIntoGlobalValues(1,&(dynamic_cast<CONTACT::CoNode*>(sele.Nodes()[s])->Dofs()[d]),&val);
    }
    for (int m=0;m<sele.NumNode();++m)
    {
      val = -gp_force*mval(m);
      fc_->SumIntoGlobalValues(1,&(dynamic_cast<CONTACT::CoNode*>(mele.Nodes()[m])->Dofs()[d]),&val);
    }

    std::map<int,double> deriv_gp_force;

    for (GEN::pairedvector<int,double>::const_iterator p=jacintcellmap.begin();p!=jacintcellmap.end();++p)
      deriv_gp_force[p->first]+=pen * p->second * wgt*gap_plus_cauchy_nn*gpn[d];
    for (std::map<int,double>::const_iterator p=gap_plus_cauchy_nn_deriv.begin();p!=gap_plus_cauchy_nn_deriv.end();++p)
      deriv_gp_force[p->first]+=pen*jac*wgt*gpn[d]*p->second;
    for (GEN::pairedvector<int,double>::const_iterator p=dnmap_unit[d].begin();p!=dnmap_unit[d].end();++p)
      deriv_gp_force[p->first]+=pen * jac*wgt*gap_plus_cauchy_nn*p->second;

    for (std::map<int,double>::const_iterator p=deriv_gp_force.begin();p!=deriv_gp_force.end();++p)
    {
      for (int s=0;s<sele.NumNode();++s)
        kc_->FEAssemble(-p->second*sval(s),dynamic_cast<CONTACT::CoNode*>(sele.Nodes()[s])->Dofs()[d],p->first);
      for (int m=0;m<mele.NumNode();++m)
        kc_->FEAssemble(p->second*mval(m),dynamic_cast<CONTACT::CoNode*>(mele.Nodes()[m])->Dofs()[d],p->first);
    }

    for (int e=0;e<Dim()-1;++e)
      for (GEN::pairedvector<int,double>::const_iterator p=dsxi[e].begin();p!=dsxi[e].end();++p)
        for (int s=0;s<sele.NumNode();++s)
          kc_->FEAssemble(-gp_force*sderiv(s,e)*p->second,dynamic_cast<CONTACT::CoNode*>(sele.Nodes()[s])->Dofs()[d],p->first);

    for (int e=0;e<Dim()-1;++e)
      for (GEN::pairedvector<int,double>::const_iterator p=dmxi[e].begin();p!=dmxi[e].end();++p)
        for (int m=0;m<mele.NumNode();++m)
          kc_->FEAssemble(gp_force*mderiv(m,e)*p->second,dynamic_cast<CONTACT::CoNode*>(mele.Nodes()[m])->Dofs()[d],p->first);

  }

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType parentdistype, int dim>
void inline CONTACT::CoIntegratorNitsche::SoEleGP(
    MORTAR::MortarElement& sele,
    const double wgt,
    const double* gpcoord,
    LINALG::Matrix<dim,1>& pxsi,
    LINALG::Matrix<dim,dim>& derivtrafo
)
{
  DRT::UTILS::CollectedGaussPoints intpoints = DRT::UTILS::CollectedGaussPoints(1); //reserve just for 1 entry ...
  intpoints.Append(gpcoord[0], gpcoord[1],0.0, wgt);

  // get coordinates of gauss point w.r.t. local parent coordinate system
  LINALG::SerialDenseMatrix pqxg(1,dim);
  derivtrafo.Clear();

  DRT::UTILS::BoundaryGPToParentGP<dim>( pqxg,
      derivtrafo,
      intpoints ,
      sele.ParentElement()->Shape(),
      sele.Shape(),
      sele.FaceParentNumber());

  // coordinates of the current integration point in parent coordinate system
    for (int idim=0;idim<dim ;idim++)
  {
    pxsi(idim) = pqxg(0,idim);
  }
}


template <int dim>
void CONTACT::CoIntegratorNitsche::SoEleCauchy(
    MORTAR::MortarElement& moEle,
    double* boundary_gpcoord,
    std::vector<GEN::pairedvector<int,double> > boundary_gpcoord_lin,
    const double gp_wgt,
    const double* contactN,
    std::vector<GEN::pairedvector<int,double> >& contactN_deriv,
    const bool useEleN,
    const double w,
    double& cauchy_nn,
    std::map<int,double>& deriv_sigma_nn)
{
  LINALG::Matrix<dim,1> pxsi(true);
  LINALG::Matrix<dim,dim> derivtravo_slave;
  DRT::Element::DiscretizationType distype = moEle.ParentElement()->Shape();
  switch (distype)
  {
  case DRT::Element::hex8:
    SoEleGP<DRT::Element::hex8,dim>(moEle,gp_wgt,boundary_gpcoord,pxsi,derivtravo_slave);
    break;
  case DRT::Element::tet4:
    SoEleGP<DRT::Element::tet4,dim>(moEle,gp_wgt,boundary_gpcoord,pxsi,derivtravo_slave);
    break;
  case DRT::Element::quad4:
    SoEleGP<DRT::Element::quad4,dim>(moEle,gp_wgt,boundary_gpcoord,pxsi,derivtravo_slave);
    break;
  case DRT::Element::quad9:
    SoEleGP<DRT::Element::quad9,dim>(moEle,gp_wgt,boundary_gpcoord,pxsi,derivtravo_slave);
    break;
  case DRT::Element::tri3:
    SoEleGP<DRT::Element::tri3,dim>(moEle,gp_wgt,boundary_gpcoord,pxsi,derivtravo_slave);
    break;
  default:
    dserror("Nitsche contact not implemented for used (bulk) elements");
  }

  LINALG::Matrix<dim,1> normal(contactN,false);

  LINALG::Matrix<dim,dim> cauchy;
  Epetra_SerialDenseMatrix dsdd;
  LINALG::Matrix<DimToNumStr<dim>::NumStr,dim> dsdpxi;
  dynamic_cast<DRT::ELEMENTS::So_base*>(moEle.ParentElement())->GetCauchyAtXi(
      pxsi,moEle.MoData().ParentDisp(),cauchy,dsdd,dsdpxi);

  LINALG::Matrix<dim,1> eleN;
  std::vector<GEN::pairedvector<int,double> > deriv_eleN(dim,moEle.NumNode()*dim);
  if (useEleN)
  {
    dynamic_cast<CoElement&>(moEle).ComputeUnitNormalAtXi(boundary_gpcoord,eleN.A());
    dynamic_cast<CoElement&>(moEle).DerivUnitNormalAtXi(boundary_gpcoord,deriv_eleN);
  }
  else
  {
    eleN=normal;
    deriv_eleN=contactN_deriv;
  }

  LINALG::Matrix<dim,1> trac_eleN;
  trac_eleN.Multiply(cauchy,eleN);
  LINALG::Matrix<dim,1> trac_n;
  trac_n.Multiply(cauchy,normal);
  cauchy_nn += w*trac_eleN.Dot(normal);

  LINALG::Matrix<dim,dim> nn;
  nn.MultiplyNT(eleN,normal);
  Epetra_SerialDenseMatrix nn_vec(DimToNumStr<dim>::NumStr,1);
  for (int i=0;i<dim;i++)nn_vec(i,0)=nn(i,i);
  if (dim==3)
  {
    nn_vec(3,0)=nn(0,1)+nn(1,0);
    nn_vec(4,0)=nn(2,1)+nn(1,2);
    nn_vec(5,0)=nn(0,2)+nn(2,0);
  }
  else
    nn_vec(2,0)=nn(0,1)+nn(1,0);

  LINALG::Matrix<dim,1> dsnndpxi;
  dsnndpxi.MultiplyTN(dsdpxi,LINALG::Matrix<DimToNumStr<dim>::NumStr,1>(nn_vec.A(),true));

  Epetra_SerialDenseMatrix dsndd(moEle.ParentElement()->NumNode()*dim,1);
  if(dsdd.Multiply(true,nn_vec,dsndd)!=0) dserror("multiply failed");

  for (int i=0;i<moEle.ParentElement()->NumNode()*dim;++i)
    deriv_sigma_nn[moEle.MoData().ParentDof().at(i)] += w*dsndd(i,0);

  for (int i=0;i<Dim()-1;++i)
    for (GEN::pairedvector<int,double>::const_iterator p=boundary_gpcoord_lin[i].begin();p!=boundary_gpcoord_lin[i].end();++p)
    {
      double& ref=deriv_sigma_nn[p->first];
      for (int k=0;k<3;++k)
        ref+=dsnndpxi(k)*derivtravo_slave(k,i)*p->second*w;
    }

  for (int d=0;d<Dim();++d)
  {
    for (GEN::pairedvector<int,double>::const_iterator p=contactN_deriv[d].begin();p!=contactN_deriv[d].end();++p)
      deriv_sigma_nn[p->first]+=trac_eleN(d)*p->second*w;
    for (GEN::pairedvector<int,double>::const_iterator p=deriv_eleN[d].begin();p!=deriv_eleN[d].end();++p)
      deriv_sigma_nn[p->first]+=trac_n(d)*p->second*w;
  }
  return;
}
