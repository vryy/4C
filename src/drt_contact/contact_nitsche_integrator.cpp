/*---------------------------------------------------------------------*/
/*!
\file contact_nitsche_integrator.cpp

\brief A class to perform integrations of nitsche related terms

\level 3

\maintainer Alexander Seitz

*/
/*---------------------------------------------------------------------*/
#include "contact_nitsche_integrator.H"

#include <Teuchos_StandardParameterEntryValidators.hpp>
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
#include "contact_nitsche_utils.H"
#include "../drt_fem_general/drt_utils_local_connectivity_matrices.H"



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

  if (dim!=Dim())
    dserror("dimension inconsistency");

  double pen = ppn_;

  const LINALG::Matrix<dim,1> normal(gpn,true);

  if (stype_==INPAR::CONTACT::solution_nitsche)
  {
    double cauchy_nn_weighted_average=0.;
    GEN::pairedvector<int,double> cauchy_nn_weighted_average_deriv(sele.NumNode()*3*12+sele.MoData().ParentDisp().size()+mele.MoData().ParentDisp().size());

    LINALG::SerialDenseVector normal_adjoint_test_slave(sele.MoData().ParentDof().size());
    GEN::pairedvector<int,LINALG::SerialDenseVector> deriv_normal_adjoint_test_slave(
        sele.MoData().ParentDof().size()
        +dnmap_unit[0].size()
        +dsxi[0].size(),
        -1,
        LINALG::SerialDenseVector(sele.MoData().ParentDof().size(),true));

    LINALG::SerialDenseVector normal_adjoint_test_master(mele.MoData().ParentDof().size());
    GEN::pairedvector<int,LINALG::SerialDenseVector> deriv_normal_adjoint_test_master(
        mele.MoData().ParentDof().size()
        +dnmap_unit[0].size()
        +dmxi[0].size(),
        -1,
        LINALG::SerialDenseVector(mele.MoData().ParentDof().size(),true));

    MAT::ElastHyper* mmat = dynamic_cast<MAT::ElastHyper*>(&*mele.ParentElement()->Material());
    MAT::ElastHyper* smat = dynamic_cast<MAT::ElastHyper*>(&*sele.ParentElement()->Material());
    if (smat==NULL || mmat==NULL)
      dserror("Nitsche contact only for elast hyper material");

    double ws=0.;
    double wm=0.;
    switch(nit_wgt_)
    {
    case INPAR::CONTACT::NitWgt_slave :
      ws=1.;wm=0.;
      pen/=dynamic_cast<CONTACT::CoElement&>(sele).TraceHE();
      break;
    case INPAR::CONTACT::NitWgt_master:
      ws=0.;wm=1.;
      pen/=dynamic_cast<CONTACT::CoElement&>(mele).TraceHE();
      break;
    case INPAR::CONTACT::NitWgt_harmonic:
      ws=1./dynamic_cast<CONTACT::CoElement&>(mele).TraceHE();
      wm=1./dynamic_cast<CONTACT::CoElement&>(sele).TraceHE();
      ws/=(ws+wm);
      wm=1.-ws;
      pen/=ws*dynamic_cast<CONTACT::CoElement&>(sele).TraceHE()+wm*dynamic_cast<CONTACT::CoElement&>(mele).TraceHE();
      break;
    default: dserror("unknown Nitsche weighting"); break;
    }

    SoEleCauchy<dim>(sele,sxi,dsxi,wgt,gpn,dnmap_unit,false,ws,cauchy_nn_weighted_average,cauchy_nn_weighted_average_deriv,normal_adjoint_test_slave ,deriv_normal_adjoint_test_slave );
    SoEleCauchy<dim>(mele,mxi,dmxi,wgt,gpn,dnmap_unit,false,wm,cauchy_nn_weighted_average,cauchy_nn_weighted_average_deriv,normal_adjoint_test_master,deriv_normal_adjoint_test_master);

    if (gap+cauchy_nn_weighted_average/pen>=0.)
    {
       if (abs(theta_)>1.e-12)
         IntegrateAdjointTest<dim>(-theta_/pen,jac,jacintcellmap,wgt,cauchy_nn_weighted_average,cauchy_nn_weighted_average_deriv,sele,normal_adjoint_test_slave ,deriv_normal_adjoint_test_slave );
      if (abs(theta_)>1.e-12)
        IntegrateAdjointTest<dim>(-theta_/pen,jac,jacintcellmap,wgt,cauchy_nn_weighted_average,cauchy_nn_weighted_average_deriv,mele,normal_adjoint_test_master,deriv_normal_adjoint_test_master);
    }

    else
    {
      // test in normal contact direction
      IntegrateTest<dim>(+pen,sele,sval,sderiv,dsxi,jac,jacintcellmap,wgt,gap,dgapgp,normal,dnmap_unit);
      IntegrateTest<dim>(-pen,mele,mval,mderiv,dmxi,jac,jacintcellmap,wgt,gap,dgapgp,normal,dnmap_unit);

      IntegrateTest<dim>(+1.,sele,sval,sderiv,dsxi,jac,jacintcellmap,wgt,cauchy_nn_weighted_average,cauchy_nn_weighted_average_deriv,normal,dnmap_unit);
      IntegrateTest<dim>(-1.,mele,mval,mderiv,dmxi,jac,jacintcellmap,wgt,cauchy_nn_weighted_average,cauchy_nn_weighted_average_deriv,normal,dnmap_unit);

      if (abs(theta_)>1.e-12)
        IntegrateAdjointTest<dim>(theta_,jac,jacintcellmap,wgt,gap,dgapgp,sele,normal_adjoint_test_slave ,deriv_normal_adjoint_test_slave );
      if (abs(theta_)>1.e-12)
        IntegrateAdjointTest<dim>(theta_,jac,jacintcellmap,wgt,gap,dgapgp,mele,normal_adjoint_test_master,deriv_normal_adjoint_test_master);
    }
  }
  else if (stype_==INPAR::CONTACT::solution_penalty)
  {
    if (gap<0.)
    {
      IntegrateTest<dim>(+pen,sele,sval,sderiv,dsxi,jac,jacintcellmap,wgt,gap,dgapgp,normal,dnmap_unit);
      IntegrateTest<dim>(-pen,mele,mval,mderiv,dmxi,jac,jacintcellmap,wgt,gap,dgapgp,normal,dnmap_unit);
    }
  }
  else
    dserror("unknown algorithm");




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
    pxsi(idim) = pqxg(0,idim);
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
    GEN::pairedvector<int,double>& deriv_sigma_nn,
    LINALG::SerialDenseVector& adjoint_test,
    GEN::pairedvector<int,LINALG::SerialDenseVector>& deriv_adjoint_test)
{
  if (useEleN)
    dserror("not supported for now");

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

  double sigma_nn;
  Epetra_SerialDenseMatrix dsdd;
  Epetra_SerialDenseMatrix dsnndd , d2snndd2 , d2snnDdDn, d2snnDdDpxi;
  LINALG::Matrix<dim,1> dsnndn,dsnndpxi;
  dynamic_cast<DRT::ELEMENTS::So_base*>(moEle.ParentElement())->GetCauchyAtXi(
      pxsi,moEle.MoData().ParentDisp(),normal,normal,sigma_nn,&dsnndd,&d2snndd2,&d2snnDdDn,NULL,&d2snnDdDpxi,&dsnndn,NULL,&dsnndpxi
      );
  d2snnDdDn.Scale(2.);
  dsnndn.Scale(2.);

  cauchy_nn += w*sigma_nn;

  for (int i=0;i<moEle.ParentElement()->NumNode()*dim;++i)
    deriv_sigma_nn[moEle.MoData().ParentDof().at(i)] += w*dsnndd(i,0);

  for (int i=0;i<dim-1;++i)
    for (GEN::pairedvector<int,double>::const_iterator p=boundary_gpcoord_lin[i].begin();p!=boundary_gpcoord_lin[i].end();++p)
    {
      double& ref=deriv_sigma_nn[p->first];
      for (int k=0;k<dim;++k)
        ref+=dsnndpxi(k)*derivtravo_slave(k,i)*p->second*w;
    }

  for (int d=0;d<dim;++d)
    for (GEN::pairedvector<int,double>::const_iterator p=contactN_deriv[d].begin();p!=contactN_deriv[d].end();++p)
      deriv_sigma_nn[p->first]+=dsnndn(d)*p->second*w;

  if (abs(theta_)>1.e-12)
    BuildNormalAdjointTest<dim>(moEle,w,dsnndd,d2snndd2,d2snnDdDn,d2snnDdDpxi,boundary_gpcoord_lin,derivtravo_slave,contactN_deriv,adjoint_test,deriv_adjoint_test);

  return;
}

template <int dim>
void CONTACT::CoIntegratorNitsche::IntegrateTest(
    const double fac,
    MORTAR::MortarElement& ele,
    const LINALG::SerialDenseVector& shape,
    const LINALG::SerialDenseMatrix& deriv,
    const std::vector<GEN::pairedvector<int,double> >& dxi,
    const double jac,const GEN::pairedvector<int,double>& jacintcellmap, const double wgt,
    const double test_val, const GEN::pairedvector<int,double>& test_deriv,
    const LINALG::Matrix<dim,1>& test_dir, const std::vector<GEN::pairedvector<int,double> >& test_dir_deriv
    )
{
  for (int d=0;d<Dim();++d)
  {
    double val = fac*jac*wgt*test_val*test_dir(d);

    for (int s=0;s<ele.NumNode();++s)
      *(ele.GetNitscheContainer().rhs(
          DRT::UTILS::getParentNodeNumberFromFaceNodeNumber(ele.ParentElement()->Shape(),ele.FaceParentNumber(),s)*dim+d))+=val*shape(s);


    std::map<int,double> val_deriv;

    for (GEN::pairedvector<int,double>::const_iterator p=jacintcellmap.begin();p!=jacintcellmap.end();++p)
      val_deriv[p->first]+=fac * p->second * wgt*test_val*test_dir(d);
    for (GEN::pairedvector<int,double>::const_iterator p=test_deriv.begin();p!=test_deriv.end();++p)
      val_deriv[p->first]+=fac*jac*wgt*test_dir(d)*p->second;
    for (GEN::pairedvector<int,double>::const_iterator p=test_dir_deriv[d].begin();p!=test_dir_deriv[d].end();++p)
      val_deriv[p->first]+=fac * jac*wgt*test_val*p->second;

    for (std::map<int,double>::const_iterator p=val_deriv.begin();p!=val_deriv.end();++p)
    {
      double* row = ele.GetNitscheContainer().k(p->first);
      for (int s=0;s<ele.NumNode();++s)
        row[DRT::UTILS::getParentNodeNumberFromFaceNodeNumber(ele.ParentElement()->Shape(),ele.FaceParentNumber(),s)*dim+d]
             +=-p->second*shape(s);
    }

    for (int e=0;e<Dim()-1;++e)
      for (GEN::pairedvector<int,double>::const_iterator p=dxi[e].begin();p!=dxi[e].end();++p)
      {
        double* row = ele.GetNitscheContainer().k(p->first);
        for (int s=0;s<ele.NumNode();++s)
          row[DRT::UTILS::getParentNodeNumberFromFaceNodeNumber(ele.ParentElement()->Shape(),ele.FaceParentNumber(),s)*dim+d]
               +=-val*deriv(s,e)*p->second;
      }
  }
}

template <int dim>
void CONTACT::CoIntegratorNitsche::BuildNormalAdjointTest(
    MORTAR::MortarElement& moEle,
    const double fac,
    const Epetra_SerialDenseMatrix& dsnndd,
    const Epetra_SerialDenseMatrix& d2snndd2,
    const Epetra_SerialDenseMatrix& d2snnDdDn,
    const Epetra_SerialDenseMatrix& d2snnDdDpxi,
    const std::vector<GEN::pairedvector<int,double> > boundary_gpcoord_lin,
    LINALG::Matrix<dim,dim> derivtravo_slave,
    const std::vector<GEN::pairedvector<int,double> >& contactN_deriv,
    LINALG::SerialDenseVector& adjoint_test,
    GEN::pairedvector<int,LINALG::SerialDenseVector>& deriv_adjoint_test)
{
  for (int i=0;i<moEle.ParentElement()->NumNode()*dim;++i)
  {
    adjoint_test(i) = fac*dsnndd(i,0);
    LINALG::SerialDenseVector& at=deriv_adjoint_test[moEle.MoData().ParentDof().at(i)];
    for (int j=0;j<moEle.ParentElement()->NumNode()*dim;++j)
      at(j)+=fac*d2snndd2(i,j);
  }

  for (int d=0;d<dim;++d)
    for (GEN::pairedvector<int,double>::const_iterator p=contactN_deriv[d].begin();p!=contactN_deriv[d].end();++p)
    {
      LINALG::SerialDenseVector& at=deriv_adjoint_test[p->first];
      for (int i=0;i<moEle.ParentElement()->NumNode()*dim;++i)
        at(i)+=fac*d2snnDdDn(i,d)*p->second;
    }

  Epetra_SerialDenseMatrix tmp(moEle.ParentElement()->NumNode()*dim,dim,false);
  Epetra_SerialDenseMatrix deriv_trafo(::View,derivtravo_slave.A(),
      derivtravo_slave.Rows(),derivtravo_slave.Rows(),derivtravo_slave.Columns());
  if (tmp.Multiply('N','N',1.,d2snnDdDpxi,deriv_trafo,0.)) dserror("multiply failed");
  for (int d=0;d<dim-1;++d)
    for (GEN::pairedvector<int,double>::const_iterator p=boundary_gpcoord_lin[d].begin();p!=boundary_gpcoord_lin[d].end();++p)
    {
      LINALG::SerialDenseVector& at=deriv_adjoint_test[p->first];
        for (int i=0;i<moEle.ParentElement()->NumNode()*dim;++i)
          at(i)+=fac*tmp(i,d)*p->second;
    }
  return;
}


template <int dim>
void CONTACT::CoIntegratorNitsche::IntegrateAdjointTest(
        const double fac,
        const double jac, const GEN::pairedvector<int,double>& jacintcellmap, const double wgt,
        const double test, const GEN::pairedvector<int,double>& deriv_test,
        MORTAR::MortarElement& moEle,
        LINALG::SerialDenseVector& adjoint_test,
        GEN::pairedvector<int,LINALG::SerialDenseVector>& deriv_adjoint_test
        )
{
  LINALG::SerialDenseVector(View,moEle.GetNitscheContainer().rhs(),moEle.MoData().ParentDof().size()).Update(-fac*jac*wgt*test,adjoint_test,1.);

  for(GEN::pairedvector<int,LINALG::SerialDenseVector>::const_iterator p=deriv_adjoint_test.begin();p!=deriv_adjoint_test.end();++p)
    LINALG::SerialDenseVector(View,moEle.GetNitscheContainer().k(p->first),moEle.MoData().ParentDof().size()).Update(fac*jac*wgt*test,p->second,1.);

  for(GEN::pairedvector<int,double>::const_iterator p=jacintcellmap.begin();p!=jacintcellmap.end();++p)
    LINALG::SerialDenseVector(View,moEle.GetNitscheContainer().k(p->first),moEle.MoData().ParentDof().size()).Update(fac*p->second*wgt*test,adjoint_test,1.);

  for(GEN::pairedvector<int,double>::const_iterator p=deriv_test.begin();p!=deriv_test.end();++p)
    LINALG::SerialDenseVector(View,moEle.GetNitscheContainer().k(p->first),moEle.MoData().ParentDof().size()).Update(fac*jac*wgt*p->second,adjoint_test,1.);

  return;
}