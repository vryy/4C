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

  if (frtype_!=INPAR::CONTACT::friction_none && dim!=3)
    dserror("only 3D friction");
  if (frtype_!=INPAR::CONTACT::friction_none
      && frtype_!=INPAR::CONTACT::friction_coulomb
      && frtype_!=INPAR::CONTACT::friction_tresca)
    dserror("only coulomb or tresca friction");
  if (frtype_==INPAR::CONTACT::friction_coulomb && frcoeff_<0.)
    dserror("negative coulomb friction coefficient");
  if (frtype_==INPAR::CONTACT::friction_tresca && frbound_<0.)
    dserror("negative tresca friction bound");

  double pen = ppn_;
  double pet = ppt_;
  typedef GEN::pairedvector<int,double>::const_iterator _CI;

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

    double ws=0.;
    double wm=0.;
    switch(nit_wgt_)
    {
    case INPAR::CONTACT::NitWgt_slave :
      ws=1.;wm=0.;
      pen/=dynamic_cast<CONTACT::CoElement&>(sele).TraceHE();
      pet/=dynamic_cast<CONTACT::CoElement&>(sele).TraceHE();
      break;
    case INPAR::CONTACT::NitWgt_master:
      ws=0.;wm=1.;
      pen/=dynamic_cast<CONTACT::CoElement&>(mele).TraceHE();
      pet/=dynamic_cast<CONTACT::CoElement&>(mele).TraceHE();
      break;
    case INPAR::CONTACT::NitWgt_harmonic:
      ws=1./dynamic_cast<CONTACT::CoElement&>(mele).TraceHE();
      wm=1./dynamic_cast<CONTACT::CoElement&>(sele).TraceHE();
      ws/=(ws+wm);
      wm=1.-ws;
      pen/=ws*dynamic_cast<CONTACT::CoElement&>(sele).TraceHE()+wm*dynamic_cast<CONTACT::CoElement&>(mele).TraceHE();
      pet/=ws*dynamic_cast<CONTACT::CoElement&>(sele).TraceHE()+wm*dynamic_cast<CONTACT::CoElement&>(mele).TraceHE();
      break;
    default: dserror("unknown Nitsche weighting"); break;
    }
    const double characteristic_timescale=1.;
    pet*=characteristic_timescale/dt_;

    // variables for friction (declaration only)
    LINALG::Matrix<dim,1> t1, t2;
    std::vector<GEN::pairedvector<int,double> > dt1,dt2;
    LINALG::Matrix<dim,1> relVel;
    std::vector<GEN::pairedvector<int,double> > relVel_deriv(dim,sele.NumNode()*dim+mele.NumNode()*dim+dsxi[0].size()+dmxi[0].size());
    double vt1,vt2;
    GEN::pairedvector<int,double> dvt1(0);
    GEN::pairedvector<int,double> dvt2(0);
    double cauchy_nt1_weighted_average=0.;
    GEN::pairedvector<int,double> cauchy_nt1_weighted_average_deriv(sele.NumNode()*3*12+sele.MoData().ParentDisp().size()+mele.MoData().ParentDisp().size());
    LINALG::SerialDenseVector t1_adjoint_test_slave(sele.MoData().ParentDof().size());
    GEN::pairedvector<int,LINALG::SerialDenseVector> deriv_t1_adjoint_test_slave(
        sele.MoData().ParentDof().size()
        +dnmap_unit[0].size()
        +dsxi[0].size(),
        -1,
        LINALG::SerialDenseVector(sele.MoData().ParentDof().size(),true));
    LINALG::SerialDenseVector t1_adjoint_test_master(mele.MoData().ParentDof().size());
    GEN::pairedvector<int,LINALG::SerialDenseVector> deriv_t1_adjoint_test_master(
        mele.MoData().ParentDof().size()
        +dnmap_unit[0].size()
        +dmxi[0].size(),
        -1,
        LINALG::SerialDenseVector(mele.MoData().ParentDof().size(),true));
    double cauchy_nt2_weighted_average=0.;
    GEN::pairedvector<int,double> cauchy_nt2_weighted_average_deriv(sele.NumNode()*3*12+sele.MoData().ParentDisp().size()+mele.MoData().ParentDisp().size());
    LINALG::SerialDenseVector t2_adjoint_test_slave(sele.MoData().ParentDof().size());
    GEN::pairedvector<int,LINALG::SerialDenseVector> deriv_t2_adjoint_test_slave(
        sele.MoData().ParentDof().size()
        +dnmap_unit[0].size()
        +dsxi[0].size(),
        -1,
        LINALG::SerialDenseVector(sele.MoData().ParentDof().size(),true));
    LINALG::SerialDenseVector t2_adjoint_test_master(mele.MoData().ParentDof().size());
    GEN::pairedvector<int,LINALG::SerialDenseVector> deriv_t2_adjoint_test_master(
        mele.MoData().ParentDof().size()
        +dnmap_unit[0].size()
        +dmxi[0].size(),
        -1,
        LINALG::SerialDenseVector(mele.MoData().ParentDof().size(),true));
    double sigma_nt1_pen_vt1=0;
    double sigma_nt2_pen_vt2=0;
    GEN::pairedvector<int,double> d_sigma_nt1_pen_vt1(
        dgapgp.capacity()
        +cauchy_nn_weighted_average_deriv.capacity()
        +cauchy_nt1_weighted_average_deriv.capacity()
        +dvt1.capacity(),0,0);
    GEN::pairedvector<int,double> d_sigma_nt2_pen_vt2(
        dgapgp.capacity()
        +cauchy_nn_weighted_average_deriv.capacity()
        +cauchy_nt2_weighted_average_deriv.capacity()
        +dvt2.capacity(),0,0);
    // variables for friction (end)

    SoEleCauchy<dim>(sele,sxi,dsxi,wgt,normal,dnmap_unit,normal,dnmap_unit,ws,
        cauchy_nn_weighted_average,cauchy_nn_weighted_average_deriv,
        normal_adjoint_test_slave ,deriv_normal_adjoint_test_slave );
    SoEleCauchy<dim>(mele,mxi,dmxi,wgt,normal,dnmap_unit,normal,dnmap_unit,wm,
        cauchy_nn_weighted_average,cauchy_nn_weighted_average_deriv,
        normal_adjoint_test_master,deriv_normal_adjoint_test_master);

    const double snn_av_pen_gap = cauchy_nn_weighted_average+pen*gap;
    GEN::pairedvector<int,double> d_snn_av_pen_gap(
        cauchy_nn_weighted_average_deriv.size()+dgapgp.size());
    for (_CI p=cauchy_nn_weighted_average_deriv.begin();p!=cauchy_nn_weighted_average_deriv.end();++p)
      d_snn_av_pen_gap[p->first]+=p->second;
    for (_CI p=dgapgp.begin();p!=dgapgp.end();++p)
      d_snn_av_pen_gap[p->first]+=pen*p->second;

    // evaluation of tangential stuff
    if (frtype_)
    {
      CONTACT::UTILS::BuildTangentVectors<dim>(normal.A(),dnmap_unit,t1.A(),dt1,t2.A(),dt2);
      RelVelInvariant<dim>(sele,sxi,dsxi,sval,sderiv,mele,mxi,dmxi,mval,mderiv,gap,dgapgp,relVel,relVel_deriv);
      CONTACT::UTILS::VectorScalarProduct<dim>(t1,dt1,relVel,relVel_deriv,vt1,dvt1);
      CONTACT::UTILS::VectorScalarProduct<dim>(t2,dt2,relVel,relVel_deriv,vt2,dvt2);

      SoEleCauchy<dim>(sele,sxi,dsxi,wgt,normal,dnmap_unit,t1,dt1,ws,
          cauchy_nt1_weighted_average,cauchy_nt1_weighted_average_deriv,
          t1_adjoint_test_slave ,deriv_t1_adjoint_test_slave);
      SoEleCauchy<dim>(mele,mxi,dmxi,wgt,normal,dnmap_unit,t1,dt1,wm,
          cauchy_nt1_weighted_average,cauchy_nt1_weighted_average_deriv,
          t1_adjoint_test_master,deriv_t1_adjoint_test_master);

      SoEleCauchy<dim>(sele,sxi,dsxi,wgt,normal,dnmap_unit,t2,dt2,ws,
          cauchy_nt2_weighted_average,cauchy_nt2_weighted_average_deriv,
          t2_adjoint_test_slave ,deriv_t2_adjoint_test_slave);
      SoEleCauchy<dim>(mele,mxi,dmxi,wgt,normal,dnmap_unit,t2,dt2,wm,
          cauchy_nt2_weighted_average,cauchy_nt2_weighted_average_deriv,
          t2_adjoint_test_master,deriv_t2_adjoint_test_master);
    }// evaluation of tangential stuff

    if (frtype_)
    {
      IntegrateAdjointTest<dim>(-theta_/pet,jac,jacintcellmap,wgt,
          cauchy_nt1_weighted_average,cauchy_nt1_weighted_average_deriv,
          sele,t1_adjoint_test_slave,deriv_t1_adjoint_test_slave);
      IntegrateAdjointTest<dim>(-theta_/pet,jac,jacintcellmap,wgt,
          cauchy_nt2_weighted_average,cauchy_nt2_weighted_average_deriv,
          sele,t2_adjoint_test_slave,deriv_t2_adjoint_test_slave);

      IntegrateAdjointTest<dim>(-theta_/pet,jac,jacintcellmap,wgt,
          cauchy_nt1_weighted_average,cauchy_nt1_weighted_average_deriv,
          mele,t1_adjoint_test_master,deriv_t1_adjoint_test_master);
      IntegrateAdjointTest<dim>(-theta_/pet,jac,jacintcellmap,wgt,
          cauchy_nt2_weighted_average,cauchy_nt2_weighted_average_deriv,
          mele,t2_adjoint_test_master,deriv_t2_adjoint_test_master);
    }

    if (snn_av_pen_gap>=0.)
    {
      IntegrateAdjointTest<dim>(-theta_/pen,jac,jacintcellmap,wgt,
          cauchy_nn_weighted_average,cauchy_nn_weighted_average_deriv,
          sele,normal_adjoint_test_slave ,deriv_normal_adjoint_test_slave );
      IntegrateAdjointTest<dim>(-theta_/pen,jac,jacintcellmap,wgt,
          cauchy_nn_weighted_average,cauchy_nn_weighted_average_deriv,
          mele,normal_adjoint_test_master,deriv_normal_adjoint_test_master);
    }
    else
    {
      // test in normal contact direction
      IntegrateTest<dim>(-1.,sele,sval,sderiv,dsxi,jac,jacintcellmap,wgt,snn_av_pen_gap,d_snn_av_pen_gap,normal,dnmap_unit);
      IntegrateTest<dim>(+1.,mele,mval,mderiv,dmxi,jac,jacintcellmap,wgt,snn_av_pen_gap,d_snn_av_pen_gap,normal,dnmap_unit);

      IntegrateAdjointTest<dim>(theta_,jac,jacintcellmap,wgt,gap,dgapgp,sele,normal_adjoint_test_slave ,deriv_normal_adjoint_test_slave );
      IntegrateAdjointTest<dim>(theta_,jac,jacintcellmap,wgt,gap,dgapgp,mele,normal_adjoint_test_master,deriv_normal_adjoint_test_master);

      if (frtype_)
      {
        double fr =0.;
        switch (frtype_)
        {
        case INPAR::CONTACT::friction_coulomb: fr=frcoeff_*(-1.)*(snn_av_pen_gap); break;
        case INPAR::CONTACT::friction_tresca:  fr=frbound_; break;
        default: fr=0.; dserror("why are you here???"); break;
        }

        double tan_tr =            sqrt(
            (cauchy_nt1_weighted_average-pet*vt1)*(cauchy_nt1_weighted_average-pet*vt1)
            +
            (cauchy_nt2_weighted_average-pet*vt2)*(cauchy_nt2_weighted_average-pet*vt2));

        // stick
        if (tan_tr<fr)
        {
          sigma_nt1_pen_vt1=cauchy_nt1_weighted_average-pet*vt1;
          for (_CI p=dvt1.begin(); p!=dvt1.end();++p)                                                          d_sigma_nt1_pen_vt1[p->first]-=pet*p->second;
          for (_CI p=cauchy_nt1_weighted_average_deriv.begin();p!=cauchy_nt1_weighted_average_deriv.end();++p) d_sigma_nt1_pen_vt1[p->first]+=p->second;

          sigma_nt2_pen_vt2=cauchy_nt2_weighted_average-pet*vt2;
          for (_CI p=dvt2.begin(); p!=dvt2.end();++p)                                                          d_sigma_nt2_pen_vt2[p->first]-=pet*p->second;
          for (_CI p=cauchy_nt2_weighted_average_deriv.begin();p!=cauchy_nt2_weighted_average_deriv.end();++p) d_sigma_nt2_pen_vt2[p->first]+=p->second;
        }
        // slip
        else
        {
          GEN::pairedvector<int,double> tmp_d(
              dgapgp.size()
              +cauchy_nn_weighted_average_deriv.size()
              +cauchy_nt1_weighted_average_deriv.size()
              +dvt1.size(),0,0);
          if (frtype_==INPAR::CONTACT::friction_coulomb)
            for (_CI p=d_snn_av_pen_gap.begin();p!=d_snn_av_pen_gap.end();++p)
              tmp_d[p->first]+=-frcoeff_/tan_tr*p->second;
          for(_CI  p=cauchy_nt1_weighted_average_deriv.begin();p!=cauchy_nt1_weighted_average_deriv.end();++p)
            tmp_d[p->first]+=-fr/(tan_tr*tan_tr*tan_tr)*(cauchy_nt1_weighted_average-pet*vt1)*p->second;
          for(_CI p=dvt1.begin();p!=dvt1.end();++p)
            tmp_d[p->first]+=-fr/(tan_tr*tan_tr*tan_tr)*(cauchy_nt1_weighted_average-pet*vt1)*(-pet)*p->second;

          for(_CI  p=cauchy_nt2_weighted_average_deriv.begin();p!=cauchy_nt2_weighted_average_deriv.end();++p)
            tmp_d[p->first]+=-fr/(tan_tr*tan_tr*tan_tr)*(cauchy_nt2_weighted_average-pet*vt2)*p->second;
          for(_CI p=dvt2.begin();p!=dvt2.end();++p)
            tmp_d[p->first]+=-fr/(tan_tr*tan_tr*tan_tr)*(cauchy_nt2_weighted_average-pet*vt2)*(-pet)*p->second;

          sigma_nt1_pen_vt1=fr/tan_tr*(cauchy_nt1_weighted_average-pet*vt1);
          for(_CI p=tmp_d.begin();p!=tmp_d.end();++p)
            d_sigma_nt1_pen_vt1[p->first]+=p->second*(cauchy_nt1_weighted_average-pet*vt1);
          for(_CI p=cauchy_nt1_weighted_average_deriv.begin();p!=cauchy_nt1_weighted_average_deriv.end();++p)
            d_sigma_nt1_pen_vt1[p->first]+=fr/tan_tr*p->second;
          for(_CI p=dvt1.begin();p!=dvt1.end();++p)
            d_sigma_nt1_pen_vt1[p->first]+=-fr/tan_tr*pet*p->second;

          sigma_nt2_pen_vt2=fr/tan_tr*(cauchy_nt2_weighted_average-pet*vt2);
          for(_CI p=tmp_d.begin();p!=tmp_d.end();++p)
            d_sigma_nt2_pen_vt2[p->first]+=p->second*(cauchy_nt2_weighted_average-pet*vt2);
          for(_CI p=cauchy_nt2_weighted_average_deriv.begin();p!=cauchy_nt2_weighted_average_deriv.end();++p)
            d_sigma_nt2_pen_vt2[p->first]+=fr/tan_tr*p->second;
          for(_CI p=dvt2.begin();p!=dvt2.end();++p)
            d_sigma_nt2_pen_vt2[p->first]+=-fr/tan_tr*pet*p->second;
        }

        IntegrateTest<dim>(-1.,sele,sval,sderiv,dsxi,jac,jacintcellmap,wgt,sigma_nt1_pen_vt1,d_sigma_nt1_pen_vt1,t1,dt1);
        IntegrateTest<dim>(+1.,mele,mval,mderiv,dmxi,jac,jacintcellmap,wgt,sigma_nt1_pen_vt1,d_sigma_nt1_pen_vt1,t1,dt1);
        IntegrateTest<dim>(-1.,sele,sval,sderiv,dsxi,jac,jacintcellmap,wgt,sigma_nt2_pen_vt2,d_sigma_nt2_pen_vt2,t2,dt2);
        IntegrateTest<dim>(+1.,mele,mval,mderiv,dmxi,jac,jacintcellmap,wgt,sigma_nt2_pen_vt2,d_sigma_nt2_pen_vt2,t2,dt2);

        IntegrateAdjointTest<dim>(theta_/pet,jac,jacintcellmap,wgt,sigma_nt1_pen_vt1,d_sigma_nt1_pen_vt1,sele,t1_adjoint_test_slave ,deriv_t1_adjoint_test_slave);
        IntegrateAdjointTest<dim>(theta_/pet,jac,jacintcellmap,wgt,sigma_nt1_pen_vt1,d_sigma_nt1_pen_vt1,mele,t1_adjoint_test_master,deriv_t1_adjoint_test_master);
        IntegrateAdjointTest<dim>(theta_/pet,jac,jacintcellmap,wgt,sigma_nt2_pen_vt2,d_sigma_nt2_pen_vt2,sele,t2_adjoint_test_slave ,deriv_t2_adjoint_test_slave);
        IntegrateAdjointTest<dim>(theta_/pet,jac,jacintcellmap,wgt,sigma_nt2_pen_vt2,d_sigma_nt2_pen_vt2,mele,t2_adjoint_test_master,deriv_t2_adjoint_test_master);
      }
    }
  }
  else if (stype_==INPAR::CONTACT::solution_penalty)
  {
    if (gap<0.)
    {
      IntegrateTest<dim>(-pen,sele,sval,sderiv,dsxi,jac,jacintcellmap,wgt,gap,dgapgp,normal,dnmap_unit);
      IntegrateTest<dim>(+pen,mele,mval,mderiv,dmxi,jac,jacintcellmap,wgt,gap,dgapgp,normal,dnmap_unit);
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
    const LINALG::Matrix<dim,1>& normal,
    std::vector<GEN::pairedvector<int,double> >& normal_deriv,
    const LINALG::Matrix<dim,1>& direction,
    std::vector<GEN::pairedvector<int,double> >& direction_deriv,
    const double w,
    double& cauchy_nt,
    GEN::pairedvector<int,double>& deriv_sigma_nt,
    LINALG::SerialDenseVector& adjoint_test,
    GEN::pairedvector<int,LINALG::SerialDenseVector>& deriv_adjoint_test)
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

  double sigma_nt;
  Epetra_SerialDenseMatrix dsdd;
  Epetra_SerialDenseMatrix dsntdd , d2sntdd2 , d2sntDdDn, d2sntDdDt, d2sntDdDpxi;
  LINALG::Matrix<dim,1> dsntdn,dsntdt,dsntdpxi;
  dynamic_cast<DRT::ELEMENTS::So_base*>(moEle.ParentElement())->GetCauchyAtXi(
      pxsi,moEle.MoData().ParentDisp(),normal,direction,sigma_nt,&dsntdd,&d2sntdd2,&d2sntDdDn,&d2sntDdDt,&d2sntDdDpxi,&dsntdn,&dsntdt,&dsntdpxi
      );

  cauchy_nt += w*sigma_nt;

  for (int i=0;i<moEle.ParentElement()->NumNode()*dim;++i)
    deriv_sigma_nt[moEle.MoData().ParentDof().at(i)] += w*dsntdd(i,0);

  for (int i=0;i<dim-1;++i)
    for (GEN::pairedvector<int,double>::const_iterator p=boundary_gpcoord_lin[i].begin();p!=boundary_gpcoord_lin[i].end();++p)
    {
      double& ref=deriv_sigma_nt[p->first];
      for (int k=0;k<dim;++k)
        ref+=dsntdpxi(k)*derivtravo_slave(k,i)*p->second*w;
    }

  for (int d=0;d<dim;++d)
    for (GEN::pairedvector<int,double>::const_iterator p=normal_deriv[d].begin();p!=normal_deriv[d].end();++p)
      deriv_sigma_nt[p->first]+=dsntdn(d)*p->second*w;

  for (int d=0;d<dim;++d)
    for (GEN::pairedvector<int,double>::const_iterator p=direction_deriv[d].begin();p!=direction_deriv[d].end();++p)
      deriv_sigma_nt[p->first]+=dsntdt(d)*p->second*w;

  if (abs(theta_)>1.e-12)
    BuildAdjointTest<dim>(moEle,w,dsntdd,d2sntdd2,d2sntDdDn,d2sntDdDt,d2sntDdDpxi,boundary_gpcoord_lin,derivtravo_slave,normal_deriv,direction_deriv,adjoint_test,deriv_adjoint_test);

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
             +=p->second*shape(s);
    }

    for (int e=0;e<Dim()-1;++e)
      for (GEN::pairedvector<int,double>::const_iterator p=dxi[e].begin();p!=dxi[e].end();++p)
      {
        double* row = ele.GetNitscheContainer().k(p->first);
        for (int s=0;s<ele.NumNode();++s)
          row[DRT::UTILS::getParentNodeNumberFromFaceNodeNumber(ele.ParentElement()->Shape(),ele.FaceParentNumber(),s)*dim+d]
               +=val*deriv(s,e)*p->second;
      }
  }
}

template <int dim>
void CONTACT::CoIntegratorNitsche::BuildAdjointTest(
    MORTAR::MortarElement& moEle,
    const double fac,
    const Epetra_SerialDenseMatrix& dsntdd,
    const Epetra_SerialDenseMatrix& d2sntdd2,
    const Epetra_SerialDenseMatrix& d2sntDdDn,
    const Epetra_SerialDenseMatrix& d2sntDdDt,
    const Epetra_SerialDenseMatrix& d2sntDdDpxi,
    const std::vector<GEN::pairedvector<int,double> > boundary_gpcoord_lin,
    LINALG::Matrix<dim,dim> derivtravo_slave,
    const std::vector<GEN::pairedvector<int,double> >& normal_deriv,
    const std::vector<GEN::pairedvector<int,double> >& direction_deriv,
    LINALG::SerialDenseVector& adjoint_test,
    GEN::pairedvector<int,LINALG::SerialDenseVector>& deriv_adjoint_test)
{
  for (int i=0;i<moEle.ParentElement()->NumNode()*dim;++i)
  {
    adjoint_test(i) = fac*dsntdd(i,0);
    LINALG::SerialDenseVector& at=deriv_adjoint_test[moEle.MoData().ParentDof().at(i)];
    for (int j=0;j<moEle.ParentElement()->NumNode()*dim;++j)
      at(j)+=fac*d2sntdd2(i,j);
  }

  for (int d=0;d<dim;++d)
    for (GEN::pairedvector<int,double>::const_iterator p=normal_deriv[d].begin();p!=normal_deriv[d].end();++p)
    {
      LINALG::SerialDenseVector& at=deriv_adjoint_test[p->first];
      for (int i=0;i<moEle.ParentElement()->NumNode()*dim;++i)
        at(i)+=fac*d2sntDdDn(i,d)*p->second;
    }

  for (int d=0;d<dim;++d)
    for (GEN::pairedvector<int,double>::const_iterator p=direction_deriv[d].begin();p!=direction_deriv[d].end();++p)
    {
      LINALG::SerialDenseVector& at=deriv_adjoint_test[p->first];
      for (int i=0;i<moEle.ParentElement()->NumNode()*dim;++i)
        at(i)+=fac*d2sntDdDt(i,d)*p->second;
    }

  Epetra_SerialDenseMatrix tmp(moEle.ParentElement()->NumNode()*dim,dim,false);
  Epetra_SerialDenseMatrix deriv_trafo(::View,derivtravo_slave.A(),
      derivtravo_slave.Rows(),derivtravo_slave.Rows(),derivtravo_slave.Columns());
  if (tmp.Multiply('N','N',1.,d2sntDdDpxi,deriv_trafo,0.)) dserror("multiply failed");
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
  if (abs(theta_)<1.e-12) return;

  LINALG::SerialDenseVector(View,moEle.GetNitscheContainer().rhs(),moEle.MoData().ParentDof().size()).Update(fac*jac*wgt*test,adjoint_test,1.);

  for(GEN::pairedvector<int,LINALG::SerialDenseVector>::const_iterator p=deriv_adjoint_test.begin();p!=deriv_adjoint_test.end();++p)
    LINALG::SerialDenseVector(View,moEle.GetNitscheContainer().k(p->first),moEle.MoData().ParentDof().size()).Update(fac*jac*wgt*test,p->second,1.);

  for(GEN::pairedvector<int,double>::const_iterator p=jacintcellmap.begin();p!=jacintcellmap.end();++p)
    LINALG::SerialDenseVector(View,moEle.GetNitscheContainer().k(p->first),moEle.MoData().ParentDof().size()).Update(fac*p->second*wgt*test,adjoint_test,1.);

  for(GEN::pairedvector<int,double>::const_iterator p=deriv_test.begin();p!=deriv_test.end();++p)
    LINALG::SerialDenseVector(View,moEle.GetNitscheContainer().k(p->first),moEle.MoData().ParentDof().size()).Update(fac*jac*wgt*p->second,adjoint_test,1.);

  return;
}


template <int dim>
void CONTACT::CoIntegratorNitsche::RelVel(
    MORTAR::MortarElement& ele,
    const LINALG::SerialDenseVector& shape,
    const LINALG::SerialDenseMatrix& deriv,
    const std::vector<GEN::pairedvector<int,double> >& dxi,
    const double fac,
    LINALG::Matrix<dim,1>& relVel,
    std::vector<GEN::pairedvector<int,double> >& relVel_deriv
    )
{
  for (int n=0;n<ele.NumNode();++n)
    for (int d=0;d<dim;++d)
      {
      relVel(d)+=fac*shape(n)*(ele.GetNodalCoords(d,n)-ele.GetNodalCoordsOld(d,n));
      relVel_deriv[d][dynamic_cast<MORTAR::MortarNode*>(ele.Nodes()[n])->Dofs()[d]]+=fac*shape(n);

      for (int sd=0;sd<dim-1;++sd)
        for(GEN::pairedvector<int,double>::const_iterator p=dxi[sd].begin();p!=dxi[sd].end();++p)
          relVel_deriv[d][p->first]+=fac*(ele.GetNodalCoords(d,n)-ele.GetNodalCoordsOld(d,n))*deriv(n,sd)*p->second;
      }
}


template <int dim>
void CONTACT::CoIntegratorNitsche::RelVelInvariant(
    MORTAR::MortarElement& sele,
    double* sxi,
    const std::vector<GEN::pairedvector<int,double> >& derivsxi,
    const LINALG::SerialDenseVector& sval,
    const LINALG::SerialDenseMatrix& sderiv,
    MORTAR::MortarElement& mele,
    double* mxi,
    const std::vector<GEN::pairedvector<int,double> >& derivmxi,
    const LINALG::SerialDenseVector& mval,
    const LINALG::SerialDenseMatrix& mderiv,
    const double& gap,
    const GEN::pairedvector<int, double>& deriv_gap,
    LINALG::Matrix<dim,1>& relVel,
    std::vector<GEN::pairedvector<int,double> >& relVel_deriv
)
{
  LINALG::Matrix<3,1> n_old;
  LINALG::Matrix<3,2> d_n_old_dxi;
  dynamic_cast<CONTACT::CoElement&>(sele).OldUnitNormalAtXi(sxi,n_old,d_n_old_dxi);
  for (int i=0;i<sele.NumNode();++i)
    for (int d=0;d<dim;++d)
    {
      relVel(d)-=sele.GetNodalCoordsOld(d,i)*sval(i);

      for (int e=0;e<dim-1;++e)
        for (GEN::pairedvector<int,double>::const_iterator p=derivsxi[e].begin();p!=derivsxi[e].end();++p)
          relVel_deriv[d][p->first]-=sele.GetNodalCoordsOld(d,i)*sderiv(i,e)*p->second;
    }

  for (int i=0;i<mele.NumNode();++i)
    for (int d=0;d<dim;++d)
    {
      relVel(d)+=mele.GetNodalCoordsOld(d,i)*mval(i);

      for (int e=0;e<dim-1;++e)
        for (GEN::pairedvector<int,double>::const_iterator p=derivmxi[e].begin();p!=derivmxi[e].end();++p)
          relVel_deriv[d][p->first]+=mele.GetNodalCoordsOld(d,i)*mderiv(i,e)*p->second;
    }

  for (int d=0;d<dim;++d)
  {
    relVel(d)+=n_old(d)*gap;

    for (int e=0;e<dim-1;++e)
      for (GEN::pairedvector<int,double>::const_iterator p=derivsxi[e].begin();p!=derivsxi[e].end();++p)
        relVel_deriv[d][p->first]+=gap*d_n_old_dxi(d,e)*p->second;
    for (GEN::pairedvector<int,double>::const_iterator p=deriv_gap.begin();p!=deriv_gap.end();++p)
      relVel_deriv[d][p->first]+=n_old(d)*p->second;
  }
}

template <int dim>
void CONTACT::UTILS::VectorScalarProduct(
    const LINALG::Matrix<dim,1>& v1,
    const std::vector<GEN::pairedvector<int,double> >& v1d,
    const LINALG::Matrix<dim,1>& v2,
    const std::vector<GEN::pairedvector<int,double> >& v2d,
    double& val,
    GEN::pairedvector<int,double>& val_deriv
)
{
  val=v1.Dot(v2);
  val_deriv.resize(v1d[0].size()+v2d[0].size());
  for (int d=0;d<dim;++d)
  {
    for (GEN::pairedvector<int,double>::const_iterator p=v1d[d].begin();p!=v1d[d].end();++p)
      val_deriv[p->first]+=v2(d)*p->second;
    for (GEN::pairedvector<int,double>::const_iterator p=v2d[d].begin();p!=v2d[d].end();++p)
      val_deriv[p->first]+=v1(d)*p->second;
  }
}

void CONTACT::UTILS::BuildTangentVectors_3D(
    const double* np,
    const std::vector<GEN::pairedvector<int,double> >& dn,
    double* t1p,
    std::vector<GEN::pairedvector<int,double> >& dt1,
    double* t2p,
    std::vector<GEN::pairedvector<int,double> >& dt2
    )
{
  const LINALG::Matrix<3,1> n(np,false);
  LINALG::Matrix<3,1> t1(t1p,true);
  LINALG::Matrix<3,1> t2(t2p,true);

  typedef GEN::pairedvector<int,double>::const_iterator _CI;
  bool z=true;
  LINALG::Matrix<3,1> tmp;
  tmp(2)=1.;
  if (abs(tmp.Dot(n))>1.-1.e-4)
    {tmp(0)=1.;tmp(2)=0.; z=false;}

  t1.CrossProduct(tmp,n);
  dt1.resize(3,
      std::max(dn[0].size(),
          std::max(dn[1].size(),dn[2].size())));
  dt2.resize(3,
      std::max(dn[0].size(),
          std::max(dn[1].size(),dn[2].size())));

  const double lt1=t1.Norm2();
  t1.Scale(1./lt1);
  LINALG::Matrix<3,3> p;
  for (int i=0;i<3;++i)
    p(i,i)=1.;
  p.MultiplyNT(-1.,t1,t1,1.);
  p.Scale(1./lt1);
  if (z)
  {
    for (_CI i=dn[1].begin();i!=dn[1].end();++i)
        for (int d=0;d<3;++d)
      dt1[d][i->first]-=p(d,0)*i->second;

    for (_CI i=dn[0].begin();i!=dn[0].end();++i)
        for (int d=0;d<3;++d)
      dt1[d][i->first]+=p(d,1)*i->second;
  }
  else
  {
    for (_CI i=dn[2].begin();i!=dn[2].end();++i)
        for (int d=0;d<3;++d)
      dt1[d][i->first]-=p(d,1)*i->second;

    for (_CI i=dn[1].begin();i!=dn[1].end();++i)
        for (int d=0;d<3;++d)
      dt1[d][i->first]+=p(d,2)*i->second;
  }

  t2.CrossProduct(n,t1);
  if (abs(t2.Norm2()-1.)>1.e-10)
    dserror("this should already form an orthonormal basis");

  for (_CI p=dn[0].begin();p!=dn[0].end();++p)
  {
    dt2[1][p->first] -= t1(2)*(p->second);
    dt2[2][p->first] += t1(1)*(p->second);
  }
  for (_CI p=dn[1].begin();p!=dn[1].end();++p)
  {
    dt2[0][p->first] += t1(2)*(p->second);
    dt2[2][p->first] -= t1(0)*(p->second);
  }
  for (_CI p=dn[2].begin();p!=dn[2].end();++p)
  {
    dt2[0][p->first] -= t1(1)*(p->second);
    dt2[1][p->first] += t1(0)*(p->second);
  }
  for (_CI p=dt1[0].begin();p!=dt1[0].end();++p)
  {
    dt2[1][p->first] += n(2)*(p->second);
    dt2[2][p->first] -= n(1)*(p->second);
  }
  for (_CI p=dt1[1].begin();p!=dt1[1].end();++p)
  {
    dt2[0][p->first] -= n(2)*(p->second);
    dt2[2][p->first] += n(0)*(p->second);
  }
  for (_CI p=dt1[2].begin();p!=dt1[2].end();++p)
  {
    dt2[0][p->first] += n(1)*(p->second);
    dt2[1][p->first] -= n(0)*(p->second);
  }

  return;
}
template <int dim>
void CONTACT::UTILS::BuildTangentVectors(
    const double* np,
    const std::vector<GEN::pairedvector<int,double> >& dn,
    double* t1p,
    std::vector<GEN::pairedvector<int,double> >& dt1,
    double* t2p,
    std::vector<GEN::pairedvector<int,double> >& dt2
    )
{
  if(dim==3) BuildTangentVectors_3D(np,dn,t1p,dt1,t2p,dt2);
  else dserror("not implemented");
}

template void CONTACT::UTILS::BuildTangentVectors<2>(
    const double* ,
    const std::vector<GEN::pairedvector<int,double> >& ,
    double* ,
    std::vector<GEN::pairedvector<int,double> >& ,
    double* ,
    std::vector<GEN::pairedvector<int,double> >&
    );

template void CONTACT::UTILS::BuildTangentVectors<3>(
    const double* ,
    const std::vector<GEN::pairedvector<int,double> >& ,
    double* ,
    std::vector<GEN::pairedvector<int,double> >& ,
    double* ,
    std::vector<GEN::pairedvector<int,double> >&
    );



template
void CONTACT::CoIntegratorNitsche::IntegrateTest<2>(
    const double,
    MORTAR::MortarElement& ,
    const LINALG::SerialDenseVector& ,
    const LINALG::SerialDenseMatrix& ,
    const std::vector<GEN::pairedvector<int,double> >& i,
    const double ,const GEN::pairedvector<int,double>& , const double ,
    const double , const GEN::pairedvector<int,double>& ,
    const LINALG::Matrix<2,1>& test_dir, const std::vector<GEN::pairedvector<int,double> >& test_dir_deriv
    );
template
void CONTACT::CoIntegratorNitsche::IntegrateTest<3>(
    const double,
    MORTAR::MortarElement& ,
    const LINALG::SerialDenseVector& ,
    const LINALG::SerialDenseMatrix& ,
    const std::vector<GEN::pairedvector<int,double> >& i,
    const double ,const GEN::pairedvector<int,double>& , const double ,
    const double , const GEN::pairedvector<int,double>& ,
    const LINALG::Matrix<3,1>& test_dir, const std::vector<GEN::pairedvector<int,double> >& test_dir_deriv
    );

template
void CONTACT::CoIntegratorNitsche::IntegrateAdjointTest<2>(
    const double ,
    const double , const GEN::pairedvector<int,double>& , const double ,
    const double , const GEN::pairedvector<int,double>& ,
    MORTAR::MortarElement& ,
    LINALG::SerialDenseVector& ,
    GEN::pairedvector<int,LINALG::SerialDenseVector>&
    );

template
void CONTACT::CoIntegratorNitsche::IntegrateAdjointTest<3>(
    const double ,
    const double , const GEN::pairedvector<int,double>& , const double ,
    const double , const GEN::pairedvector<int,double>& ,
    MORTAR::MortarElement& ,
    LINALG::SerialDenseVector& ,
    GEN::pairedvector<int,LINALG::SerialDenseVector>&
    );

template
void CONTACT::CoIntegratorNitsche::BuildAdjointTest<2>(
    MORTAR::MortarElement& ,
    const double ,
    const Epetra_SerialDenseMatrix& ,
    const Epetra_SerialDenseMatrix& ,
    const Epetra_SerialDenseMatrix& ,
    const Epetra_SerialDenseMatrix& ,
    const Epetra_SerialDenseMatrix& ,
    const std::vector<GEN::pairedvector<int,double> > ,
    LINALG::Matrix<2,2> ,
    const std::vector<GEN::pairedvector<int,double> >& ,
    const std::vector<GEN::pairedvector<int,double> >& ,
    LINALG::SerialDenseVector& ,
    GEN::pairedvector<int,LINALG::SerialDenseVector>& );

template
void CONTACT::CoIntegratorNitsche::BuildAdjointTest<3>(
    MORTAR::MortarElement& ,
    const double ,
    const Epetra_SerialDenseMatrix& ,
    const Epetra_SerialDenseMatrix& ,
    const Epetra_SerialDenseMatrix& ,
    const Epetra_SerialDenseMatrix& ,
    const Epetra_SerialDenseMatrix& ,
    const std::vector<GEN::pairedvector<int,double> > ,
    LINALG::Matrix<3,3> ,
    const std::vector<GEN::pairedvector<int,double> >& ,
    const std::vector<GEN::pairedvector<int,double> >& ,
    LINALG::SerialDenseVector& ,
    GEN::pairedvector<int,LINALG::SerialDenseVector>& );


template
void CONTACT::CoIntegratorNitsche::RelVel<2>(
    MORTAR::MortarElement& ,
    const LINALG::SerialDenseVector& ,
    const LINALG::SerialDenseMatrix& ,
    const std::vector<GEN::pairedvector<int,double> >& ,
    const double ,
    LINALG::Matrix<2,1>& ,
    std::vector<GEN::pairedvector<int,double> >&
    );

template
void CONTACT::CoIntegratorNitsche::RelVel<3>(
    MORTAR::MortarElement& ,
    const LINALG::SerialDenseVector& ,
    const LINALG::SerialDenseMatrix& ,
    const std::vector<GEN::pairedvector<int,double> >& ,
    const double ,
    LINALG::Matrix<3,1>& ,
    std::vector<GEN::pairedvector<int,double> >&
    );

template
void CONTACT::UTILS::VectorScalarProduct<2>(
    const LINALG::Matrix<2,1>& ,
    const std::vector<GEN::pairedvector<int,double> >& ,
    const LINALG::Matrix<2,1>& ,
    const std::vector<GEN::pairedvector<int,double> >& ,
    double& ,
    GEN::pairedvector<int,double>&
);
template
void CONTACT::UTILS::VectorScalarProduct<3>(
    const LINALG::Matrix<3,1>& ,
    const std::vector<GEN::pairedvector<int,double> >& ,
    const LINALG::Matrix<3,1>& ,
    const std::vector<GEN::pairedvector<int,double> >& ,
    double& ,
    GEN::pairedvector<int,double>&
);
