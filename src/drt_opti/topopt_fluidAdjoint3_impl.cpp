/*!------------------------------------------------------------------------------------------------*
\file topopt_fluidAdjoint3_impl.cpp

\brief Functionality of the element level of the fluid adjoint equations

<pre>
Maintainer: Martin Winklmaier
            winklmaier@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15241
</pre>
 *------------------------------------------------------------------------------------------------*/


#include "topopt_fluidAdjoint3_impl.H"
#include "topopt_fluidAdjoint3_impl_parameter.H"

#include "../drt_fluid_ele/fluid_ele.H"
#include "../drt_fluid_ele/fluid_ele_utils.H"
#include "../drt_fem_general/drt_utils_gder2.H"
#include "../drt_geometry/position_array.H"
#include "../drt_inpar/inpar_topopt.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_condition_utils.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_lib/standardtypes_cpp.H"
#include "../drt_mat/newtonianfluid.H"



#include "../linalg/linalg_utils.H"



//----------------------------------------------------------------------*
//----------------------------------------------------------------------*/
DRT::ELEMENTS::FluidAdjoint3ImplInterface* DRT::ELEMENTS::FluidAdjoint3ImplInterface::Impl(DRT::Element::DiscretizationType distype)
{
  switch(distype)
  {
  case DRT::Element::hex8:
  {
    return FluidAdjoint3Impl<DRT::Element::hex8>::Instance();
  }
  case DRT::Element::hex20:
  {
    return FluidAdjoint3Impl<DRT::Element::hex20>::Instance();
  }
  case DRT::Element::hex27:
  {
    return FluidAdjoint3Impl<DRT::Element::hex27>::Instance();
  }
//  case DRT::Element::tet4:
//  {
//    return FluidAdjoint3Impl<DRT::Element::tet4>::Instance();
//  }
//  case DRT::Element::tet10:
//  {
//    return FluidAdjoint3Impl<DRT::Element::tet10>::Instance();
//  }
  case DRT::Element::quad4:
  {
    return FluidAdjoint3Impl<DRT::Element::quad4>::Instance();
  }
  case DRT::Element::quad8:
  {
    return FluidAdjoint3Impl<DRT::Element::quad8>::Instance();
  }
  case DRT::Element::quad9:
  {
    return FluidAdjoint3Impl<DRT::Element::quad9>::Instance();
  }
//  case DRT::Element::tri3:
//  {
//    return FluidAdjoint3Impl<DRT::Element::tri3>::Instance();
//  }
//  case DRT::Element::tri6:
//  {
//    return FluidAdjoint3Impl<DRT::Element::tri6>::Instance();
//  }
  // no 1D elements
  default:
  {
    dserror("Element shape %s not activated. Just do it.",DRT::DistypeToString(distype).c_str());
    break;
  }
  }
  return NULL;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::FluidAdjoint3Impl<distype> * DRT::ELEMENTS::FluidAdjoint3Impl<distype>::Instance( bool create )
{
  static FluidAdjoint3Impl<distype> * instance;
  if ( create )
  {
    if ( instance==NULL )
      instance = new FluidAdjoint3Impl<distype>();
  }
  else
  {
    if ( instance!=NULL )
      delete instance;
    instance = NULL;
  }
  return instance;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::FluidAdjoint3Impl<distype>::Done()
{
  // delete this pointer! Afterwards we have to go! But since this is a
  // cleanup call, we can do it this way.
  Instance( false );
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::FluidAdjoint3Impl<distype>::FluidAdjoint3Impl()
  : xyze_(true),
    funct_(true),
    deriv_(true),
    deriv2_(true),
    xjm_(true),
    xji_(true),
    derxy_(true),
    derxy2_(true),
    velint_(true),
    vderxy_(true),
    gradp_(true),
    fluidvelint_(true),
    fluidvelxy_(true),
    bodyforce_(true),
    contforce_(0.0),
    vdiv_(0.0),
    conv1_(true),
    conv2_(true),
    velint_old_(true),
    vderxy_old_(true),
    gradp_old_(true),
    fluidvelint_old_(true),
    fluidvelxy_old_(true),
    bodyforce_old_(true),
    contforce_old_(0.0),
    vdiv_old_(true),
    conv1_old_(true),
    conv2_old_(true),
    tau_(true),
    intpoints_( distype ),
    xsi_(true),
    det_(0.0),
    fac_(0.0),
    reacoeff_(0.0),
    is_higher_order_ele_(false)
{
  // pointer to class FluidImplParameter (access to the general parameter)
  fluidAdjoint3Parameter_ = DRT::ELEMENTS::FluidAdjoint3ImplParameter::Instance();
}



/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
int DRT::ELEMENTS::FluidAdjoint3Impl<distype>::Evaluate(DRT::ELEMENTS::Fluid*    ele,
                                                 DRT::Discretization & discretization,
                                                 const std::vector<int> & lm,
                                                 Teuchos::ParameterList&    params,
                                                 Teuchos::RCP<MAT::Material> & mat,
                                                 Epetra_SerialDenseMatrix&  elesysmat,
                                                 Epetra_SerialDenseVector&  elerhs
)
{
  return Evaluate( ele, discretization, lm, params, mat,
                   elesysmat, elerhs, intpoints_ );
}



/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
int DRT::ELEMENTS::FluidAdjoint3Impl<distype>::Evaluate(DRT::ELEMENTS::Fluid*    ele,
                                                 DRT::Discretization & discretization,
                                                 const std::vector<int> & lm,
                                                 Teuchos::ParameterList&    params,
                                                 Teuchos::RCP<MAT::Material> & mat,
                                                 Epetra_SerialDenseMatrix&  elesysmat,
                                                 Epetra_SerialDenseVector&  elerhs,
                                                 const DRT::UTILS::GaussIntegration & intpoints
)
{
  // construct views
  LINALG::Matrix<(nsd_+1)*nen_,(nsd_+1)*nen_> elemat(elesysmat,true);
  LINALG::Matrix<(nsd_+1)*nen_,            1> elevec(elerhs,true);

  // ---------------------------------------------------------------------
  // get all general state vectors: fluid/adjoint velocity/pressure
  // velocity/pressure values are at time n/n+1
  // ---------------------------------------------------------------------
  // fill the local element vector/matrix with the global values
  // ost:         velocity/pressure at time n+1
  LINALG::Matrix<nsd_,nen_> eveln(true);
  LINALG::Matrix<nen_,1>    epren(true);
  ExtractValuesFromGlobalVector(discretization,lm, &eveln, &epren,"veln");

  LINALG::Matrix<nsd_,nen_> evelnp(true);
  LINALG::Matrix<nen_,1>    eprenp(true);
  ExtractValuesFromGlobalVector(discretization,lm, &evelnp, &eprenp,"velnp");

  LINALG::Matrix<nsd_,nen_> efluidveln(true);
  ExtractValuesFromGlobalVector(discretization,lm, &efluidveln, NULL,"fluidveln");

  LINALG::Matrix<nsd_,nen_> efluidvelnp(true);
  ExtractValuesFromGlobalVector(discretization,lm, &efluidvelnp, NULL,"fluidvelnp");

  // evaluate nodal porosities
  LINALG::Matrix<nen_,1> eporo(true);
  {
    // read nodal values from global vector
    RCP<const Epetra_Vector> topopt_porosity = params.get<RCP<const Epetra_Vector> >("topopt_porosity");
    for (int nn=0;nn<nen_;++nn)
    {
      int lid = (ele->Nodes()[nn])->LID();
      eporo(nn,0) = (*topopt_porosity)[lid];
    }
  }

  // get node coordinates and number of elements per node
  GEO::fillInitialPositionArray<distype,nsd_,LINALG::Matrix<nsd_,nen_> >(ele,xyze_);

  // flag for higher order elements
  is_higher_order_ele_ = IsHigherOrder<distype>::ishigherorder;

  // overrule higher_order_ele if input-parameter is set
  // this might be interesting for fast (but slightly
  // less accurate) computations
  if (fluidAdjoint3Parameter_->IsInconsistent() == true) is_higher_order_ele_ = false;
  // TODO deactivate this maybe?!

  // ---------------------------------------------------------------------
  // call routine for calculating element matrix and right hand side
  // ---------------------------------------------------------------------
  Sysmat(ele->Id(),
         eveln,
         evelnp,
         epren,
         eprenp,
         efluidveln,
         efluidvelnp,
         elemat,
         elevec,
         eporo,
         mat,
         intpoints);

  return 0;
}


/*----------------------------------------------------------------------*
 |  calculate element matrix and right hand side        winklmaier 03/12|
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::FluidAdjoint3Impl<distype>::Sysmat(
  int                                           eid,
  const LINALG::Matrix<nsd_,nen_>&              eveln,
  const LINALG::Matrix<nsd_,nen_>&              evelnp,
  const LINALG::Matrix<nen_,1>&                 epren,
  const LINALG::Matrix<nen_,1>&                 eprenp,
  const LINALG::Matrix<nsd_,nen_> &             efluidveln,
  const LINALG::Matrix<nsd_,nen_> &             efluidvelnp,
  LINALG::Matrix<(nsd_+1)*nen_,(nsd_+1)*nen_>&  estif,
  LINALG::Matrix<(nsd_+1)*nen_,1>&              eforce,
  const LINALG::Matrix<nen_,1> &                eporo,
  Teuchos::RCP<const MAT::Material>             material,
  const DRT::UTILS::GaussIntegration &          intpoints
  )
{
  //------------------------------------------------------------------------
  //  preliminary definitions and evaluations
  //------------------------------------------------------------------------
  // definition of matrices
  LINALG::Matrix<nen_*nsd_,nen_*nsd_>  estif_w_v(true);
  LINALG::Matrix<nen_*nsd_,nen_>       estif_w_q(true);
  LINALG::Matrix<nen_, nen_*nsd_>      estif_r_v(true);
  LINALG::Matrix<nen_,nen_>            estif_r_q(true);

  // definition of vectors
  LINALG::Matrix<nen_,1>     preforce(true);
  LINALG::Matrix<nsd_,nen_>  velforce(true);

  // definition of velocity-based momentum residual vectors
  LINALG::Matrix<nsd_*nsd_,nen_>  lin_resM_Du(true);
  LINALG::Matrix<nsd_,1>          resM_Du(true);

  // evaluate shape functions and derivatives at element center
  EvalShapeFuncAndDerivsAtEleCenter(eid);

  // calculate subgrid viscosity and/or stabilization parameter at element center
  if (not fluidAdjoint3Parameter_->EvalTauAtGP())
  {
    // get velocity at element center
    velint_.Multiply(evelnp,funct_);

    // calculate stabilization parameters at element center
    CalcStabParameter(fac_);
  }

int gp = 0;
  //------------------------------------------------------------------------
  //  start loop over integration points
  //------------------------------------------------------------------------
  for ( DRT::UTILS::GaussIntegration::const_iterator iquad=intpoints.begin(); iquad!=intpoints.end(); ++iquad )
  {gp++;
    // evaluate shape functions and derivatives at integration point
    EvalShapeFuncAndDerivsAtIntPoint(iquad,eid);

    //----------------------------------------------------------------------
    //  evaluation of various values at integration point:
    //  1) velocity (including derivatives)
    //  2) fluid velocity (including derivatives)
    //  3) pressure (including derivatives)
    //  4) body-force vector
    //  5) "history" vector for momentum equation
    //----------------------------------------------------------------------
    // get velocity at integration point
    // 1) t^n=last iteration 2) t^n+1 = last time step
    velint_.Multiply(evelnp,funct_);
    velint_old_.Multiply(eveln,funct_);

    // get velocity derivatives at integration point
    // 1) t^n=last iteration 2) t^n+1 = last time step
    vderxy_.MultiplyNT(evelnp,derxy_);
    vderxy_old_.MultiplyNT(eveln,derxy_);

    // get fluid velocity at integration point
    // 1) t^n=last iteration 2) t^n+1 = last time step
    fluidvelint_.Multiply(efluidvelnp,funct_);
    fluidvelint_old_.Multiply(efluidveln,funct_);

    // get fluid velocity derivatives at integration point
    // 1) t^n=last iteration 2) t^n+1 = last time step
    fluidvelxy_.MultiplyNT(efluidvelnp,derxy_);
    fluidvelxy_old_.MultiplyNT(efluidveln,derxy_);

    // get pressure at integration point
    // 1) t^n=last iteration 2) t^n+1 = last time step
    double press = funct_.Dot(eprenp);
    double press_old = funct_.Dot(epren);

    // get pressure gradient at integration point
    // 1) t^n=last iteration 2) t^n+1 = last time step
    gradp_.Multiply(derxy_,eprenp);
    gradp_old_.Multiply(derxy_,epren);

    // get reaction coefficient due to porosity for topology optimization
    // !do this only at gauss point!
    // TODO does it make problems to evaluate at element center? (i think it should, winklmaier)
    reacoeff_ = funct_.Dot(eporo);

    // calculate stabilization parameter at integration point
    if (fluidAdjoint3Parameter_->EvalTauAtGP())
      CalcStabParameter(fac_);

    BodyForce(efluidveln,efluidvelnp);
    ContForce();

    // get first convective value at integration point
    // 1) t^n=last iteration 2) t^n+1 = last time step
    conv1_.Multiply(vderxy_,fluidvelint_);
    conv1_old_.Multiply(vderxy_old_,fluidvelint_old_);

    // get second convective value at integration point
    // 1) t^n=last iteration 2) t^n+1 = last time step
    conv2_.MultiplyTN(fluidvelxy_,velint_);
    conv2_old_.MultiplyTN(fluidvelxy_old_,velint_old_);

    // get divergence at integration point
    // 1) t^n=last iteration 2) t^n+1 = last time step
    vdiv_ = vdiv_old_ = 0.0;
    for (int idim = 0; idim <nsd_; ++idim)
    {
      vdiv_ += vderxy_(idim,idim);
      vdiv_old_ += vderxy_old_(idim,idim);
    }

    //----------------------------------------------------------------------
    // set time-integration factors for left- and right-hand side
    //----------------------------------------------------------------------
    const double timefacfac       = fluidAdjoint3Parameter_->Timefac()      * fac_;
    const double timefacfacrhs    = fluidAdjoint3Parameter_->TimefacRhs()   * fac_;

    const double timefacfacpre    = fluidAdjoint3Parameter_->TimefacPre()   * fac_;
    const double timefacfacprerhs = fluidAdjoint3Parameter_->TimefacPreRhs()* fac_;

    const double timefacfacdiv    = fluidAdjoint3Parameter_->TimefacDiv()   * fac_;
    const double timefacfacdivrhs = fluidAdjoint3Parameter_->TimefacDivRhs()* fac_;

    /* ------------------------------------------------------------------------ *
    * standard terms                                                            *
    * ------------------------------------------------------------------------- */

    // 1) mass matrix + reactive term
    MassReactionGalPart(estif_w_v,velforce,timefacfac,timefacfacrhs);

    // 2) convection terms
    ConvectionGalPart(estif_w_v,velforce,timefacfac,timefacfacrhs);

    // 3) viscous terms
    ViscousGalPart(estif_w_v,velforce,timefacfac,timefacfacrhs);

    Epetra_SerialDenseMatrix est(nen_*nsd_,nen_);
    for (int i=0;i<nen_*nsd_;i++){for (int j=0;j<nen_;j++) est(i,j) = -estif(i,j);}
    // 4) pressure term
    PressureGalPart(estif_w_q,velforce,timefacfacpre,timefacfacprerhs,press,press_old);
//    for (int i=0;i<nen_*nsd_;i++){for (int j=0;j<nen_;j++) est(i,j) = -estif(i,j);}
//    if (eid==0)
//    {
//      std::ostringstream filename;
//      filename << "tests/adjointsysmat_part_ele0_gp" << gp << ".mtl";
//      LINALG::PrintSerialDenseMatrixInMatlabFormat(filename.str(),est);
//    }

    // 5) continuity term
    ContinuityGalPart(estif_r_v,preforce,timefacfacdiv,timefacfacdivrhs);

    // 6) standard Galerkin bodyforce term on right-hand side
    BodyForceGalPart(velforce,timefacfac,timefacfacrhs);

    // 7) standard right-hand side term of continuity equation forces
    ContForceGalPart(preforce,timefacfacdiv,timefacfacdivrhs);
    /* ------------------------------------------------------------------------ *
    * standard terms done                                                       *
    * ------------------------------------------------------------------------- */



    /* ------------------------------------------------------------------------ *
    * stabilization part                                                        *
    * ------------------------------------------------------------------------- */

    if ((fluidAdjoint3Parameter_->PSPG() == INPAR::FLUID::pstab_use_pspg) or
        (fluidAdjoint3Parameter_->SUPG() == INPAR::FLUID::convective_stab_supg))
    {
      if (fluidAdjoint3Parameter_->AdjointType() == INPAR::TOPOPT::discrete_adjoint)
      {
        // prework for supg/psgp - stabilization: evaluate strong residual

        /* order of the derivatives in GalMomResnU is:
         * from 1 to nsd:       col-dim = x, row-dim = 1-nsd
         * from nsd+1 to 2*nsd: col-dim = y, row-dim = 1-nsd
         * and so on. so the outer loop is the column dimension
         * and the inner loop the row dimension */
        LINALG::Matrix<nsd_*nsd_,nen_> GalMomTestStat(true);

        DiscreteGalMom(GalMomTestStat,
            timefacfac,
            timefacfacrhs,
            timefacfacpre,
            timefacfacprerhs,
            eveln,
            evelnp,
            efluidveln,
            efluidvelnp);

        // 8) PSPG term
        if (fluidAdjoint3Parameter_->PSPG() == INPAR::FLUID::pstab_use_pspg)
        {
          LINALG::Matrix<nen_,nen_> estif(true);estif-=estif_r_q;

          DiscretePSPG(estif_w_q,
              estif_r_q,
              velforce,
              preforce,
              GalMomTestStat,
              timefacfac,
              timefacfacrhs,
              timefacfacpre,
              timefacfacprerhs);

          estif+=estif_r_q;
//          if (eid==0) cout << "entry is " << estif << endl;
        }

        // 9) SUPG term
        if (fluidAdjoint3Parameter_->SUPG() == INPAR::FLUID::convective_stab_supg)
        {
          DiscreteSUPG(estif_w_v,
              estif_w_q,
              velforce,
              preforce,
              GalMomTestStat,
              timefacfac,
              timefacfacrhs,
              timefacfacpre,
              timefacfacprerhs);
        }
      }
      else if (fluidAdjoint3Parameter_->AdjointType() == INPAR::TOPOPT::cont_adjoint)
      {
        // prework for supg/psgp - stabilization: evaluate strong residual

        /* order of the derivatives in GalMomResnU is:
         * from 1 to nsd:       col-dim = x, row-dim = 1-nsd
         * from nsd+1 to 2*nsd: col-dim = y, row-dim = 1-nsd
         * and so on. so the outer loop is the column dimension
         * and the inner loop the row dimension */
        LINALG::Matrix<nsd_*nsd_,nen_> GalMomResnU(true);

        // strong residual of momentum equation of last iteration, scaled with fac*dt/rho
        LINALG::Matrix<nsd_,1> StrongResMomScaled(true);

        MomRes(GalMomResnU,
            StrongResMomScaled,
            timefacfac,
            timefacfacrhs,
            timefacfacpre,
            timefacfacprerhs,
            eveln,
            evelnp,
            efluidveln,
            efluidvelnp);

        // 8) PSPG term
        if (fluidAdjoint3Parameter_->PSPG() == INPAR::FLUID::pstab_use_pspg)
        {
          PSPG(estif_r_v,
              estif_r_q,
              preforce,
              GalMomResnU,
              StrongResMomScaled,
              timefacfac,
              timefacfacrhs,
              timefacfacpre,
              timefacfacprerhs);
        }

        // 9) SUPG term
        if (fluidAdjoint3Parameter_->SUPG() == INPAR::FLUID::convective_stab_supg)
        {
          SUPG(estif_w_v,
              estif_w_q,
              velforce,
              GalMomResnU,
              StrongResMomScaled,
              timefacfac,
              timefacfacrhs,
              timefacfacpre,
              timefacfacprerhs);
        }
      }
    }

    // 10) continuity stabilization
    if (fluidAdjoint3Parameter_->CStab() == INPAR::FLUID::continuity_stab_yes)
    {
      ContStab(estif_w_v,
          velforce,
          timefacfacdiv,
          timefacfacdivrhs);
    }
  }
  //------------------------------------------------------------------------
  //  end loop over integration points
  //------------------------------------------------------------------------


  Epetra_SerialDenseMatrix test(nen_*(nsd_+1),nen_*(nsd_+1));
  if (eid==0)
    for (int i=0;i<nen_*(nsd_+1);i++){for (int j=0;j<nen_*(nsd_+1);j++) test(i,j) = -estif(i,j);}



  //------------------------------------------------------------------------
  //  add contributions to element matrix and right-hand-side vector
  //------------------------------------------------------------------------
  // add pressure part to right-hand-side vector
  for (int vi=0; vi<nen_; ++vi)
  {
    eforce(numdofpernode_*vi+nsd_)+=preforce(vi);
  }

  // add velocity part to right-hand-side vector
  for (int vi=0; vi<nen_; ++vi)
  {
    for (int idim=0; idim<nsd_; ++idim)
    {
      eforce(numdofpernode_*vi+idim)+=velforce(idim,vi);
    }
  }

  // add pressure-pressure part to matrix
  for (int ui=0; ui<nen_; ++ui)
  {
    const int fuippp = numdofpernode_*ui+nsd_;

    for (int vi=0; vi<nen_; ++vi)
    {
      const int numdof_vi_p_nsd = numdofpernode_*vi+nsd_;

      estif(numdof_vi_p_nsd,fuippp)+=estif_r_q(vi,ui);
    }
  }

  // add velocity-velocity part to matrix
  for (int ui=0; ui<nen_; ++ui)
  {
    const int numdof_ui = numdofpernode_*ui;
    const int nsd_ui = nsd_*ui;

    for (int jdim=0; jdim < nsd_;++jdim)
    {
      const int numdof_ui_jdim = numdof_ui+jdim;
      const int nsd_ui_jdim = nsd_ui+jdim;

      for (int vi=0; vi<nen_; ++vi)
      {
        const int numdof_vi = numdofpernode_*vi;
        const int nsd_vi = nsd_*vi;

        for (int idim=0; idim <nsd_; ++idim)
        {
          estif(numdof_vi+idim, numdof_ui_jdim) += estif_w_v(nsd_vi+idim, nsd_ui_jdim);
        }
      }
    }
  }

  // add velocity-pressure part to matrix
  for (int ui=0; ui<nen_; ++ui)
  {
    const int numdof_ui_nsd = numdofpernode_*ui + nsd_;

    for (int vi=0; vi<nen_; ++vi)
    {
      const int nsd_vi = nsd_*vi;
      const int numdof_vi = numdofpernode_*vi;

      for (int idim=0; idim <nsd_; ++idim)
      {
        estif(numdof_vi+idim, numdof_ui_nsd) += estif_w_q(nsd_vi+idim, ui);
      }
    }
  }

  // add pressure-velocity part to matrix
  for (int ui=0; ui<nen_; ++ui)
  {
    const int numdof_ui = numdofpernode_*ui;
    const int nsd_ui = nsd_*ui;

    for (int jdim=0; jdim < nsd_;++jdim)
    {
      const int numdof_ui_jdim = numdof_ui+jdim;
      const int nsd_ui_jdim = nsd_ui+jdim;

      for (int vi=0; vi<nen_; ++vi)
        estif(numdofpernode_*vi+nsd_, numdof_ui_jdim) += estif_r_v(vi, nsd_ui_jdim);
    }
  }



//  if (eid==0)
//  {
//    for (int i=0;i<nen_*(nsd_+1);i++){for (int j=0;j<nen_*(nsd_+1);j++) test(i,j) += estif(i,j);}
//    std::ostringstream filename;
//    filename << "tests/adjointsysmat_ele0.mtl";
//    LINALG::PrintSerialDenseMatrixInMatlabFormat(filename.str(),test);
//  }



  return;
}



/*---------------------------------------------------------------------------------*
 | evaluate shape functions and derivatives at element center     winklmaier 03/12 |
 *---------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::FluidAdjoint3Impl<distype>::EvalShapeFuncAndDerivsAtEleCenter(
  const int  eleid
)
{
  // use one-point Gauss rule
  DRT::UTILS::IntPointsAndWeights<nsd_> intpoints_stab(DRT::ELEMENTS::DisTypeToStabGaussRule<distype>::rule);

  // coordinates of the current integration point
  const double* gpcoord = (intpoints_stab.IP().qxg)[0];
  for (int idim=0;idim<nsd_;idim++)
  {
    xsi_(idim) = gpcoord[idim];
  }
  const double wquad = intpoints_stab.IP().qwgt[0];

  // shape functions and their first derivatives
  DRT::UTILS::shape_function<distype>(xsi_,funct_);
  DRT::UTILS::shape_function_deriv1<distype>(xsi_,deriv_);
  if (is_higher_order_ele_)
  {
    // get the second derivatives of standard element at current GP
    DRT::UTILS::shape_function_deriv2<distype>(xsi_,deriv2_);
  }


  // compute Jacobian matrix and determinant
  // actually compute its transpose....
  /*
    +-            -+ T      +-            -+
    | dx   dx   dx |        | dx   dy   dz |
    | --   --   -- |        | --   --   -- |
    | dr   ds   dt |        | dr   dr   dr |
    |              |        |              |
    | dy   dy   dy |        | dx   dy   dz |
    | --   --   -- |   =    | --   --   -- |
    | dr   ds   dt |        | ds   ds   ds |
    |              |        |              |
    | dz   dz   dz |        | dx   dy   dz |
    | --   --   -- |        | --   --   -- |
    | dr   ds   dt |        | dt   dt   dt |
    +-            -+        +-            -+
   */

  // get Jacobian matrix and determinant
  xjm_.MultiplyNT(deriv_,xyze_);
  det_ = xji_.Invert(xjm_);

  // check for degenerated elements
  if (det_ < 1E-16)
    dserror("GLOBAL ELEMENT NO.%i\nZERO OR NEGATIVE JACOBIAN DETERMINANT: %f", eleid, det_);

  // compute integration factor
  fac_ = wquad*det_;

  // compute global first derivates
  derxy_.Multiply(xji_,deriv_);

  //--------------------------------------------------------------
  //             compute global second derivatives
  //--------------------------------------------------------------
  if (is_higher_order_ele_)
  {
    DRT::UTILS::gder2<distype>(xjm_,derxy_,deriv2_,xyze_,derxy2_);
  }
  else derxy2_.Clear();

  return;
}



/*-------------------------------------------------------------------------------*
 | evaluate shape functions and derivatives at integr. point    winklmaier 02/12 |
 *-------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::FluidAdjoint3Impl<distype>::EvalShapeFuncAndDerivsAtIntPoint(
    DRT::UTILS::GaussIntegration::iterator & iquad,       // actual integration point
    const int                              eleid        // element ID
)
{
  // coordinates of the current integration point
  const double* gpcoord = iquad.Point();
  for (int idim=0;idim<nsd_;idim++)
  {
     xsi_(idim) = gpcoord[idim];
  }

  // shape functions and their first derivatives
  DRT::UTILS::shape_function<distype>(xsi_,funct_);
  DRT::UTILS::shape_function_deriv1<distype>(xsi_,deriv_);
  derxy2_.Clear();
  if (is_higher_order_ele_)
  {
    // get the second derivatives of standard element at current GP
    DRT::UTILS::shape_function_deriv2<distype>(xsi_,deriv2_);
  }

  // get Jacobian matrix and determinant
  // actually compute its transpose....
  /*
    +-            -+ T      +-            -+
    | dx   dx   dx |        | dx   dy   dz |
    | --   --   -- |        | --   --   -- |
    | dr   ds   dt |        | dr   dr   dr |
    |              |        |              |
    | dy   dy   dy |        | dx   dy   dz |
    | --   --   -- |   =    | --   --   -- |
    | dr   ds   dt |        | ds   ds   ds |
    |              |        |              |
    | dz   dz   dz |        | dx   dy   dz |
    | --   --   -- |        | --   --   -- |
    | dr   ds   dt |        | dt   dt   dt |
    +-            -+        +-            -+
  */
  xjm_.MultiplyNT(deriv_,xyze_);
  det_ = xji_.Invert(xjm_);

  if (det_ < 1E-16)
    dserror("GLOBAL ELEMENT NO.%i\nZERO OR NEGATIVE JACOBIAN DETERMINANT: %f", eleid, det_);

  // compute integration factor
  fac_ = iquad.Weight()*det_;

  // compute global first derivates
  derxy_.Multiply(xji_,deriv_);

  //--------------------------------------------------------------
  //             compute global second derivatives
  //--------------------------------------------------------------
  if (is_higher_order_ele_)
  {
    DRT::UTILS::gder2<distype>(xjm_,derxy_,deriv2_,xyze_,derxy2_);
  }
  else derxy2_.Clear();

  return;
}



/*----------------------------------------------------------------------*
 |  calculation of stabilization parameter             winklmaier 03/12 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::FluidAdjoint3Impl<distype>::CalcStabParameter(const double vol)
{
  //---------------------------------------------------------------------
  // preliminary definition of values which will already be computed for
  // tau_M and later be used for tau_C again by some of the subsequent
  // stabilization parameter definitions
  //---------------------------------------------------------------------
  double traceG = 0.0;
  double Gnormu = 0.0;
  double Gvisc  = 0.0;

  double strle    = 0.0;
  double hk       = 0.0;
  double fluidvel_norm = 0.0;
  double re12     = 0.0;
  double c3       = 0.0;

  // material parameters
  const double dens = fluidAdjoint3Parameter_->Density();
  const double visc = fluidAdjoint3Parameter_->Viscosity();

  //---------------------------------------------------------------------
  // first step: computation of tau_M with the following options
  // (both with or without inclusion of dt-part):
  // A) definition according to Taylor et al. (1998)
  //    -> see also Gravemeier and Wall (2010) for version for
  //       variable-density flow at low Mach number
  // B) combined definition according to Franca and Valentin (2000) as
  //    well as Barrenechea and Valentin (2002)
  //    -> differentiating tau_Mu and tau_Mp for this definition
  // C) definition according to Shakib (1989) / Shakib and Hughes (1991)
  //    -> differentiating tau_Mu and tau_Mp for this definition
  // D) definition according to Codina (1998)
  //    -> differentiating tau_Mu and tau_Mp for this definition
  // E) definition according to Franca et al. (2005) as well as Badia
  //    and Codina (2010)
  //    -> only for Darcy or Darcy-Stokes/Brinkman flow, hence only
  //       tau_Mp for this definition
  //---------------------------------------------------------------------
  // get element-type constant for tau
  const double mk = DRT::ELEMENTS::MK<distype>();

  // computation depending on which parameter definition is used
  switch (fluidAdjoint3Parameter_->TauType())
  {
  case INPAR::FLUID::tau_taylor_hughes_zarins:
  case INPAR::FLUID::tau_taylor_hughes_zarins_wo_dt:
  case INPAR::FLUID::tau_taylor_hughes_zarins_whiting_jansen:
  case INPAR::FLUID::tau_taylor_hughes_zarins_whiting_jansen_wo_dt:
  case INPAR::FLUID::tau_taylor_hughes_zarins_scaled:
  case INPAR::FLUID::tau_taylor_hughes_zarins_scaled_wo_dt:
  {
    /*

    literature:
    1) C.A. Taylor, T.J.R. Hughes, C.K. Zarins, Finite element modeling
       of blood flow in arteries, Comput. Methods Appl. Mech. Engrg. 158
       (1998) 155-196.
    2) V. Gravemeier, W.A. Wall, An algebraic variational multiscale-
       multigrid method for large-eddy simulation of turbulent variable-
       density flow at low Mach number, J. Comput. Phys. 229 (2010)
       6047-6070.
       -> version for variable-density low-Mach-number flow as implemented
          here, which corresponds to version for incompressible flow as
          given in the previous publications when density is constant

                                                                           1
                     +-                                               -+ - -
                     |        2                                        |   2
                     | c_1*rho                                  2      |
          tau  = C * | -------   +  c_2*rho*u*G*rho*u  +  c_3*mu *G:G  |
             M       |     2                                           |
                     |   dt                                            |
                     +-                                               -+

          with the constants and covariant metric tensor defined as follows:

          C   = 1.0 (not explicitly defined here),
          c_1 = 4.0 (for version with dt), 0.0 (for version without dt),
          c_2 = 1.0 (not explicitly defined here),
          c_3 = 12.0/m_k (36.0 for linear and 144.0 for quadratic elements)

                  +-           -+   +-           -+   +-           -+
                  |             |   |             |   |             |
                  |  dr    dr   |   |  ds    ds   |   |  dt    dt   |
            G   = |  --- * ---  | + |  --- * ---  | + |  --- * ---  |
             ij   |  dx    dx   |   |  dx    dx   |   |  dx    dx   |
                  |    i     j  |   |    i     j  |   |    i     j  |
                  +-           -+   +-           -+   +-           -+

                  +----
                   \
          G : G =   +   G   * G
                   /     ij    ij
                  +----
                   i,j
                             +----
                             \
          rho*u*G*rho*u  =   +   rho*u * G  *rho*u
                             /        i   ij      j
                            +----
                              i,j
    */

    // total reaction coefficient sigma_tot: sum of "artificial" reaction
    // due to time factor and reaction coefficient (reaction coefficient
    // ensured to remain zero in GetMaterialParams for non-reactive material)
    double sigma_tot = reacoeff_;
    if (fluidAdjoint3Parameter_->TauType() == INPAR::FLUID::tau_taylor_hughes_zarins or
        fluidAdjoint3Parameter_->TauType() == INPAR::FLUID::tau_taylor_hughes_zarins_whiting_jansen or
        fluidAdjoint3Parameter_->TauType() == INPAR::FLUID::tau_taylor_hughes_zarins_scaled)
      sigma_tot += 1.0/fluidAdjoint3Parameter_->Dt();

    // definition of constants as described above
    const double c1 = 4.0;
    c3 = 12.0/mk;

    // computation of various values derived from covariant metric tensor
    // (trace of covariant metric tensor required for computation of tau_C below)
    double G;
    double normG = 0.0;
    const double dens_sqr = dens*dens;
    for (int nn=0;nn<nsd_;++nn)
    {
      const double dens_sqr_velint_nn = dens_sqr*fluidvelint_(nn);
      for (int mm=0; mm<nsd_; ++mm)
      {
        traceG += xji_(nn,mm)*xji_(nn,mm);
      }
      for (int rr=0;rr<nsd_;++rr)
      {
        G = xji_(nn,0)*xji_(rr,0);
        for (int mm=1; mm<nsd_; ++mm)
        {
          G += xji_(nn,mm)*xji_(rr,mm);
        }
        normG  += G*G;
        Gnormu += dens_sqr_velint_nn*G*fluidvelint_(rr);
      }
    }

    // compute viscous part
    Gvisc = c3*visc*visc*normG;

    // computation of stabilization parameters tau_Mu and tau_Mp
    // -> identical for the present definitions
    tau_(0) = 1.0/(sqrt(c1*dens_sqr*DSQR(sigma_tot) + Gnormu + Gvisc));
    tau_(1) = tau_(0);
    break;
  }

  case INPAR::FLUID::tau_franca_barrenechea_valentin_frey_wall:
  {
    /*

    literature:
    1) L.P. Franca, F. Valentin, On an improved unusual stabilized
       finite element method for the advective-reactive-diffusive
       equation, Comput. Methods Appl. Mech. Engrg. 190 (2000) 1785-1800.
    2) G.R. Barrenechea, F. Valentin, An unusual stabilized finite
       element method for a generalized Stokes problem, Numer. Math.
       92 (2002) 652-677.


                  xi1,xi2 ^
                          |      /
                          |     /
                          |    /
                        1 +---+
                          |
                          |
                          |
                          +--------------> re1,re2
                              1

    */
    // get velocity norm
    fluidvel_norm = fluidvelint_.Norm2();

    // total reaction coefficient sigma_tot: sum of "artificial" reaction
    // due to time factor and reaction coefficient (reaction coefficient
    // ensured to remain zero in GetMaterialParams for non-reactive material)
    const double sigma_tot = 1.0/fluidAdjoint3Parameter_->Timefac() + reacoeff_;

    // calculate characteristic element length
    CalcCharEleLength(vol,fluidvel_norm,strle,hk);

    // various parameter computations for case with dt:
    // relating viscous to reactive part (re01: tau_Mu, re11: tau_Mp)
    const double re01 = 4.0 * visc / (mk * dens * sigma_tot * DSQR(strle));
    const double re11 = 4.0 * visc / (mk * dens * sigma_tot * DSQR(hk));

    // relating convective to viscous part (re02: tau_Mu, re12: tau_Mp)
    const double re02 = mk * dens * fluidvel_norm * strle / (2.0 * visc);
                 re12 = mk * dens * fluidvel_norm * hk / (2.0 * visc);

    // respective "switching" parameters
    const double xi01 = std::max(re01,1.0);
    const double xi11 = std::max(re11,1.0);
    const double xi02 = std::max(re02,1.0);
    const double xi12 = std::max(re12,1.0);

    tau_(0) = DSQR(strle)/(DSQR(strle)*dens*sigma_tot*xi01+(4.0*visc/mk)*xi02);
    tau_(1) = DSQR(hk)/(DSQR(hk)*dens*sigma_tot*xi11+(4.0*visc/mk)*xi12);
    break;
  }

  case INPAR::FLUID::tau_franca_barrenechea_valentin_frey_wall_wo_dt:
  {
    /*

     stabilization parameter as above without inclusion of dt-part

    */
    // get velocity norm
    fluidvel_norm = fluidvelint_.Norm2();

    // calculate characteristic element length
    CalcCharEleLength(vol,fluidvel_norm,strle,hk);

    // various parameter computations for case without dt:
    // relating viscous to reactive part (re01: tau_Mu, re11: tau_Mp)
    double re01 = 4.0 * visc / (mk * dens * reacoeff_ * DSQR(strle));
    double re11 = 4.0 * visc / (mk * dens * reacoeff_ * DSQR(hk));

    // relating convective to viscous part (re02: tau_Mu, re12: tau_Mp)
    const double re02 = mk * dens * fluidvel_norm * strle / (2.0 * visc);
                 re12 = mk * dens * fluidvel_norm * hk / (2.0 * visc);

    // respective "switching" parameters
    const double xi01 = std::max(re01,1.0);
    const double xi11 = std::max(re11,1.0);
    const double xi02 = std::max(re02,1.0);
    const double xi12 = std::max(re12,1.0);

    tau_(0) = DSQR(strle)/(DSQR(strle)*dens*reacoeff_*xi01+(4.0*visc/mk)*xi02);
    tau_(1) = DSQR(hk)/(DSQR(hk)*dens*reacoeff_*xi11+(4.0*visc/mk)*xi12);
    break;
  }

  case INPAR::FLUID::tau_shakib_hughes_codina:
  case INPAR::FLUID::tau_shakib_hughes_codina_wo_dt:
  {
    /*

    literature:
    1) F. Shakib, Finite element analysis of the compressible Euler and
       Navier-Stokes equations, PhD thesis, Division of Applied Mechanics,
       Stanford University, Stanford, CA, USA, 1989.
    2) F. Shakib, T.J.R. Hughes, A new finite element formulation for
       computational fluid dynamics: IX. Fourier analysis of space-time
       Galerkin/least-squares algorithms, Comput. Methods Appl. Mech.
       Engrg. 87 (1991) 35-58.
    3) R. Codina, Stabilized finite element approximation of transient
       incompressible flows using orthogonal subscales, Comput. Methods
       Appl. Mech. Engrg. 191 (2002) 4295-4321.

       constants defined as in Shakib (1989) / Shakib and Hughes (1991),
       merely slightly different with respect to c_3:

       c_1 = 4.0 (for version with dt), 0.0 (for version without dt),
       c_2 = 4.0,
       c_3 = 4.0/(m_k*m_k) (36.0 for linear, 576.0 for quadratic ele.)

       Codina (2002) proposed present version without dt and explicit
       definition of constants
       (condition for constants as defined here: c_2 <= sqrt(c_3)).

    */
    // get velocity norm
    fluidvel_norm = fluidvelint_.Norm2();

    // calculate characteristic element length
    CalcCharEleLength(vol,fluidvel_norm,strle,hk);

    // total reaction coefficient sigma_tot: sum of "artificial" reaction
    // due to time factor and reaction coefficient (reaction coefficient
    // ensured to remain zero in GetMaterialParams for non-reactive material)
    double sigma_tot = reacoeff_;
    if (fluidAdjoint3Parameter_->TauType() == INPAR::FLUID::tau_shakib_hughes_codina)
      sigma_tot += 1.0/fluidAdjoint3Parameter_->Dt();

    // definition of constants as described above
    const double c1 = 4.0;
    const double c2 = 4.0;
    c3 = 4.0/(mk*mk);
    // alternative value as proposed in Shakib (1989): c3 = 16.0/(mk*mk);

    tau_(0) = 1.0/(sqrt(c1*DSQR(dens)*DSQR(sigma_tot)
                      + c2*DSQR(dens)*DSQR(fluidvel_norm)/DSQR(strle)
                      + c3*DSQR(visc)/(DSQR(strle)*DSQR(strle))));
    tau_(1) = 1.0/(sqrt(c1*DSQR(dens)*DSQR(sigma_tot)
                      + c2*DSQR(dens)*DSQR(fluidvel_norm)/DSQR(hk)
                      + c3*DSQR(visc)/(DSQR(hk)*DSQR(hk))));
    break;
  }
  case INPAR::FLUID::tau_codina:
  case INPAR::FLUID::tau_codina_wo_dt:
  {
    /*

      literature:
         R. Codina, Comparison of some finite element methods for solving
         the diffusion-convection-reaction equation, Comput. Methods
         Appl. Mech. Engrg. 156 (1998) 185-210.

         constants:
         c_1 = 1.0 (for version with dt), 0.0 (for version without dt),
         c_2 = 2.0,
         c_3 = 4.0/m_k (12.0 for linear, 48.0 for quadratic elements)

         Codina (1998) proposed present version without dt.

    */
    // get velocity norm
    fluidvel_norm = fluidvelint_.Norm2();

    // calculate characteristic element length
    CalcCharEleLength(vol,fluidvel_norm,strle,hk);

    // total reaction coefficient sigma_tot: sum of "artificial" reaction
    // due to time factor and reaction coefficient (reaction coefficient
    // ensured to remain zero in GetMaterialParams for non-reactive material)
    double sigma_tot = reacoeff_;
    if (fluidAdjoint3Parameter_->TauType() == INPAR::FLUID::tau_codina)
      sigma_tot += 1.0/fluidAdjoint3Parameter_->Dt();

    // definition of constants as described above
    const double c1 = 1.0;
    const double c2 = 2.0;
    c3 = 4.0/mk;

    tau_(0) = 1.0/(sqrt(c1*dens*sigma_tot
                      + c2*dens*fluidvel_norm/strle
                      + c3*visc/DSQR(strle)));
    tau_(1) = 1.0/(sqrt(c1*dens*sigma_tot
                      + c2*dens*fluidvel_norm/hk
                      + c3*visc/DSQR(hk)));
    break;
  }
  default:
  {
    dserror("unknown definition for tau_M\n %i  ", fluidAdjoint3Parameter_->TauType());
    break;
  }
  }  // end switch (fluidAdjoint3Parameter_->whichtau_)


  //---------------------------------------------------------------------
  // second step: computation of tau_C with the following options:
  // A) definition according to Taylor et al. (1998)
  // B) definition according to Whiting (1999)/Whiting and Jansen (2001)
  // C) scaled version of definition according to Taylor et al. (1998)
  // D) definition according to Wall (1999)
  // E) definition according to Codina (2002)
  // F) definition according to Badia and Codina (2010)
  //    (only for Darcy or Darcy-Stokes/Brinkman flow)
  //---------------------------------------------------------------------
  // computation depending on which parameter definition is used
  switch (fluidAdjoint3Parameter_->TauType())
  {
  case INPAR::FLUID::tau_taylor_hughes_zarins:
  case INPAR::FLUID::tau_taylor_hughes_zarins_wo_dt:
  {
    /*

    literature:
       C.A. Taylor, T.J.R. Hughes, C.K. Zarins, Finite element modeling
       of blood flow in arteries, Comput. Methods Appl. Mech. Engrg. 158
       (1998) 155-196.

                                              1/2
                           (c_2*rho*u*G*rho*u)
                    tau  = -------------------
                       C       trace (G)


       -> see respective definitions for computation of tau_M above

    */

    tau_(2) = sqrt(Gnormu)/traceG;
  }
  break;

  case INPAR::FLUID::tau_taylor_hughes_zarins_whiting_jansen:
  case INPAR::FLUID::tau_taylor_hughes_zarins_whiting_jansen_wo_dt:
  {
    /*

    literature:
    1) C.H. Whiting, Stabilized finite element methods for fluid dynamics
       using a hierarchical basis, PhD thesis, Rensselaer Polytechnic
       Institute, Troy, NY, USA, 1999.
    2) C.H. Whiting, K.E. Jansen, A stabilized finite element method for
       the incompressible Navier-Stokes equations using a hierarchical
       basis, Int. J. Numer. Meth. Fluids 35 (2001) 93-116.

                                  1.0
                    tau  = ------------------
                       C    tau  * trace (G)
                               M

       -> see respective definitions for computation of tau_M above

    */

    tau_(2) = 1.0/(tau_(0)*traceG);
  }
  break;

  case INPAR::FLUID::tau_taylor_hughes_zarins_scaled:
  case INPAR::FLUID::tau_taylor_hughes_zarins_scaled_wo_dt:
  {
    /*

      Caution: This is an experimental version of a stabilization
               parameter definition which scales the definition
               for tau_C by Taylor et al. (1998) in a similar
               way as proposed below by Franca and Frey (1992)
               and Wall (1999) by appropriately defining an
               element Reynolds number based on the covariant
               metric tensor.

                  /                        1/2    \
                  |  /                    \       |                       1/2
                  | |  c_2*rho*u*G*rho*u  |       |    (c_2*rho*u*G*rho*u)
      tau  =  MIN | | ------------------- | | 1.0 | *  -------------------
         C        | |          2          |       |         trace (G)
                  | \    c_3*mu *G:G      /       |
                  \                               /
                    |                     |
                    -----------------------
                    element Reynolds number
                      based on covariant
                        metric tensor

       -> see respective definitions for computation of tau_M above

    */

    // element Reynolds number based on covariant metric tensor
    const double reG = std::sqrt(Gnormu/Gvisc);

    // "switching" parameter
    const double xi_tau_c = std::min(reG,1.0);

    tau_(2) = xi_tau_c*sqrt(Gnormu)/traceG;
  }
  break;

  case INPAR::FLUID::tau_franca_barrenechea_valentin_frey_wall:
  case INPAR::FLUID::tau_franca_barrenechea_valentin_frey_wall_wo_dt:
  {
    /*

    literature:
    1) L.P. Franca, S.L. Frey, Stabilized finite element methods:
       II. The incompressible Navier-Stokes equations, Comput. Methods
       Appl. Mech. Engrg. 99 (1992) 209-293.
    2) W.A. Wall, Fluid-Struktur-Interaktion mit stabilisierten Finiten
       Elementen, Dissertation, Universitaet Stuttgart, 1999.

                 xi_tau_c ^
                          |
                        1 |   +-----------
                          |  /
                          | /
                          |/
                          +--------------> re12
                              1

       -> see respective definitions for computation of tau_M above

    */

    // "switching" parameter
    const double xi_tau_c = std::min(re12,1.0);

    tau_(2) = 0.5 * dens * fluidvel_norm * hk * xi_tau_c;
  }
  break;

  case INPAR::FLUID::tau_shakib_hughes_codina:
  case INPAR::FLUID::tau_shakib_hughes_codina_wo_dt:
  case INPAR::FLUID::tau_codina:
  case INPAR::FLUID::tau_codina_wo_dt:
  {
    /*

    literature:
       R. Codina, Stabilized finite element approximations of transient
       incompressible flows using orthogonal subscales, Comput. Methods
       Appl. Mech. Engrg. 191 (2002) 4295-4321.

       -> see respective definitions for computation of tau_M above

    */

    tau_(2) = DSQR(hk)/(sqrt(c3)*tau_(1));
  }
  break;
  default:
  {
    dserror("unknown definition for tau_C\n %i  ", fluidAdjoint3Parameter_->TauType());
    break;
  }
  }  // end switch (fluidAdjoint3Parameter_->whichtau_)

  return;
}


/*----------------------------------------------------------------------*
 |  calculation of characteristic element length       winklmaier 03/12 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::FluidAdjoint3Impl<distype>::CalcCharEleLength(
    const double  vol,
    const double  fluidvel_norm,
    double&       strle,
    double&       hk
) const
{
  // cast dimension to a double varibale -> pow()
  const double dim = double (nsd_);

  //! direction of flow (normed velocity vector)
  LINALG::Matrix<nsd_,1> fluidvelino;

  //---------------------------------------------------------------------
  // various definitions for characteristic element length for tau_Mu
  //---------------------------------------------------------------------
  // a) streamlength due to Tezduyar et al. (1992) -> default
  // normed velocity vector
  if (fluidvel_norm>=1e-6) fluidvelino.Update(1.0/fluidvel_norm,fluidvelint_);
  else
  {
    fluidvelino.Clear();
    fluidvelino(0,0) = 1.0;
  }

  LINALG::Matrix<nen_,1> tmp;
  tmp.MultiplyTN(derxy_,fluidvelino);
  const double val = tmp.Norm1();
  strle = 2.0/val;

  // b) volume-equivalent diameter (warning: 3-D formula!)
  //strle = std::pow((6.*vol/M_PI),(1.0/3.0))/sqrt(3.0);

  // c) cubic/square root of element volume/area
  //strle = std::pow(vol,1/dim);

  //---------------------------------------------------------------------
  // various definitions for characteristic element length for tau_Mp
  //---------------------------------------------------------------------
  // a) volume-equivalent diameter -> default for 3-D computations
  if (nsd_==3) hk = std::pow((6.*vol/M_PI),(1.0/3.0))/sqrt(3.0);

  // b) square root of element area -> default for 2-D computations,
  // may also alternatively be used for 3-D computations
  else if (nsd_==2) hk = std::pow(vol,1/dim);
  // check for potential 1-D computations
  else dserror("element length calculation not implemented for 1-D computation!");

  return;
}



/*---------------------------------------------------------------------------------*
 | compute bodyforce of momentum equation                         winklmaier 03/12 |
 *---------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::FluidAdjoint3Impl<distype>::BodyForce(
           const LINALG::Matrix<nsd_,nen_>&           efluidveln,
           const LINALG::Matrix<nsd_,nen_>&           efluidvelnp
)
{
  bodyforce_.Clear();
  bodyforce_old_.Clear();

  if (fluidAdjoint3Parameter_->TestCase() == INPAR::TOPOPT::adjointtest_no)
  {
    if (fluidAdjoint3Parameter_->ObjDissipationTerm())
    {
      const double dissipation = fluidAdjoint3Parameter_->ObjDissipationFac();

      /* ------------------------------------------------------------------------ *
       * 1) evaluate bodyforce at new time step                                   *
       * ------------------------------------------------------------------------ */

      // dissipation term due to reaction
      bodyforce_.Update(-2*dissipation*reacoeff_,fluidvelint_);

      // dissipation term due to viscosity
      if (is_higher_order_ele_) // TODO check this
      {
        LINALG::Matrix<nsd_,numderiv2_> fluidvelxy2(true);
        fluidvelxy2.MultiplyNT(efluidvelnp,derxy2_);

        LINALG::Matrix<nsd_,1> laplaceU(true);
        for (int idim=0;idim<nsd_;++idim)
        {
          for (int jdim=0;jdim<nsd_;++jdim)
            laplaceU(idim) += fluidvelxy2(idim,jdim);
        }

        bodyforce_.Update(dissipation*fluidAdjoint3Parameter_->Viscosity(),laplaceU,1.0);
      }


      /* ------------------------------------------------------------------------ *
       * 2) evaluate bodyforce at old time step in instationary case              *
       * ------------------------------------------------------------------------ */
      if (not fluidAdjoint3Parameter_->IsStationary())
      {
        bodyforce_old_.Update(-2*dissipation*reacoeff_,fluidvelint_old_);

        // dissipation term due to viscosity
        if (is_higher_order_ele_) // TODO check this
        {
          LINALG::Matrix<nsd_,numderiv2_> fluidvelxy2_old(true);
          fluidvelxy2_old.MultiplyNT(efluidveln,derxy2_);

          LINALG::Matrix<nsd_,1> laplaceU_old(true);
          for (int idim=0;idim<nsd_;++idim)
          {
            for (int jdim=0;jdim<nsd_;++jdim)
              laplaceU_old(idim) += fluidvelxy2_old(idim,jdim);
          }

          bodyforce_old_.Update(dissipation*fluidAdjoint3Parameter_->Viscosity(),laplaceU_old,1.0);
        }
      }
    }
  }
  else // special cases
  {
    // get global coordinates of gauss point
    double x = 0.0;
    double y = 0.0;
    LINALG::Matrix<nsd_,1> coords(true);
    coords.Multiply(xyze_,funct_);
    x = coords(0);
    y = coords(1); // z-component currently not required in tests


    switch (fluidAdjoint3Parameter_->TestCase())
    {
    case INPAR::TOPOPT::adjointtest_stat_const_vel_lin_pres:
    {
      break;
    }
    case INPAR::TOPOPT::adjointtest_stat_lin_vel_quad_pres:
    {
      bodyforce_(0) = bodyforce_old_(0) = 4*y;
      bodyforce_(1) = bodyforce_old_(1) = 24994*x;
      break;
    }
    case INPAR::TOPOPT::adjointtest_stat_quad_vel_lin_pres:
    {
      bodyforce_(0) = bodyforce_old_(0) = 24999*x*x - 4*x*y;
      bodyforce_(1) = bodyforce_old_(1) = -74989*x*x + 50002*y*y
                                          + 12*x*y;
      break;
    }
    case INPAR::TOPOPT::adjointtest_stat_all_terms_all_constants:
    {
      bodyforce_(0) = bodyforce_old_(0) = 24998*x*x + 24994*x*y - 4*y*y;
      bodyforce_(1) = bodyforce_old_(1) = -74978*x*x + 50004*y*y + 28*x*y;
      break;
    }
    case INPAR::TOPOPT::adjointtest_instat_varying_theta:
    {
      double t = fluidAdjoint3Parameter_->Time();
      bodyforce_(0) = 5*x;
      bodyforce_(1) = -5*y;

      t += fluidAdjoint3Parameter_->Dt(); // old time = t + dt
      bodyforce_old_(0) = 5*x;
      bodyforce_old_(1) = -5*y;
      break;
    }
    case INPAR::TOPOPT::adjointtest_instat_all_terms_all_constants:
    {
      double t = fluidAdjoint3Parameter_->Time();
      bodyforce_(0) = -10*x*y + 5*x*x - 4*y*y*t + 9*x*y*t;
      bodyforce_(1) = -8*y*y*t + x*x + 24*x*y + 18*y*y*t*t + 4*x*y*t;

      t += fluidAdjoint3Parameter_->Dt(); // old time = t + dt
      bodyforce_old_(0) = -10*x*y + 5*x*x - 4*y*y*t + 9*x*y*t;
      bodyforce_old_(1) = -8*y*y*t + x*x + 24*x*y + 18*y*y*t*t + 4*x*y*t;
      break;
    }
    case INPAR::TOPOPT::adjointtest_instat_primal_and_dual:
    {
      double t = fluidAdjoint3Parameter_->Time();
      bodyforce_(0) = -3*x*y - 3*y*t*t - 6*y*y*t - 6*y + 8*x*y*t;
      bodyforce_(1) = 9*x + x*t + 6*x*y*t;

      t += fluidAdjoint3Parameter_->Dt(); // old time = t + dt
      bodyforce_old_(0) = -3*x*y - 3*y*t*t - 6*y*y*t - 6*y + 8*x*y*t;
      bodyforce_old_(1) = 9*x + x*t + 6*x*y*t;
      break;
    }
    case INPAR::TOPOPT::adjointtest_primal:
      break;
    default:
    {
      dserror("no dirichlet condition implemented for special test case");
      break;
    }
    }
  }
}



/*---------------------------------------------------------------------------------*
 | compute continuity force                                       winklmaier 03/12 |
 *---------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::FluidAdjoint3Impl<distype>::ContForce(
)
{
  contforce_ = 0.0;
  contforce_old_ = 0.0;

  if (fluidAdjoint3Parameter_->TestCase() == INPAR::TOPOPT::adjointtest_no)
  {
    ; // currently no domain pressure entries in objective -> no entry here
  }
  else
  {
    // get global coordinates of gauss point
    double x = 0.0;
    double y = 0.0;
    LINALG::Matrix<nsd_,1> coords(true);
    coords.Multiply(xyze_,funct_);
    x = coords(0);
    y = coords(1); // z-component currently not required in tests

    switch (fluidAdjoint3Parameter_->TestCase())
    {
    case INPAR::TOPOPT::adjointtest_stat_const_vel_lin_pres:
    {
      break;
    }
    case INPAR::TOPOPT::adjointtest_stat_lin_vel_quad_pres:
    {
      contforce_ = contforce_old_ = 12.0;
      break;
    }
    case INPAR::TOPOPT::adjointtest_stat_quad_vel_lin_pres:
    {
      contforce_ = contforce_old_ = 2*x + 4*y;
      break;
    }
    case INPAR::TOPOPT::adjointtest_stat_all_terms_all_constants:
    {
      contforce_ = contforce_old_ = 2*x + 5*y;
      break;
    }
    case INPAR::TOPOPT::adjointtest_instat_varying_theta:
    {
      contforce_ = contforce_old_ = 0;
      break;
    }
    case INPAR::TOPOPT::adjointtest_instat_all_terms_all_constants:
    {
      double t = fluidAdjoint3Parameter_->Time();
      contforce_ = 2*x + y*t + 4*y*t*t;

      t += fluidAdjoint3Parameter_->Dt(); // old time = t + dt
      contforce_old_ = 2*x + y*t + 4*y*t*t;
      break;
    }
    case INPAR::TOPOPT::adjointtest_instat_primal_and_dual:
    {
      double t = fluidAdjoint3Parameter_->Time();
      contforce_ = y*t + 1 + t;

      t += fluidAdjoint3Parameter_->Dt(); // old time = t + dt
      contforce_old_ = y*t + 1 + t;
      break;
    }
    case INPAR::TOPOPT::adjointtest_primal:
      break;
    default:
    {
      dserror("no dirichlet condition implemented for special test case");
      break;
    }
    }
  }
}



/*---------------------------------------------------------------------------------*
 | compute mass and reaction terms of galerkin part               winklmaier 03/12 |
 *---------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::FluidAdjoint3Impl<distype>::MassReactionGalPart(
    LINALG::Matrix<nen_*nsd_,nen_*nsd_> &   estif_w_v,
    LINALG::Matrix<nsd_,nen_> &             velforce,
    const double &                          timefacfac,
    const double &                          timefacfacrhs
) const
{
  /* inertia (contribution to mass matrix) if not is_stationary */
  /*
            /              \
           |                |
           |    rho*Dv , w  |
           |                |
            \              /
  */
  /*  reaction */
  /*
            /                \
           |                  |
           |    sigma*Dv , w  |
           |                  |
            \                /
  */
  double massreacfac = 0.0; // factor summing up coefficients of reactive term and mass-matrix
  if (fluidAdjoint3Parameter_->IsStationary())
    massreacfac = reacoeff_*timefacfac;
  else
    massreacfac = fluidAdjoint3Parameter_->Density()*fac_+reacoeff_*timefacfac; // fac -> mass matrix // reac*timefacfac/dens -> reactive

  for (int ui=0; ui<nen_; ++ui)
  {
    const int fui   = nsd_*ui;

    const double uifunct = massreacfac*funct_(ui);

    for (int vi=0; vi<nen_; ++vi)
    {
      const int fvi   = nsd_*vi;

      const double value = funct_(vi)*uifunct;

      for (int idim = 0; idim <nsd_; ++idim)
      {
        estif_w_v(fvi+idim,fui+idim) += value;
      } // end for (idim)
    } //vi
  } // ui

  // rhs at new time step
  LINALG::Matrix<nsd_,1> scaled_vel(true);
  scaled_vel.Update(massreacfac,velint_);
  for (int vi=0;vi<nen_;++vi)
  {
    for (int jdim=0;jdim<nsd_;++jdim)
    {
      velforce(jdim,vi)-=funct_(vi)*scaled_vel(jdim);
    }
  }

  // rhs at old time step
  if (not fluidAdjoint3Parameter_->IsStationary())
  {
    double massreacfacrhs = -fluidAdjoint3Parameter_->Density()*fac_+reacoeff_*timefacfacrhs; // fac -> mass matrix // reac*timefacfac/dens -> reactive
    scaled_vel.Update(massreacfacrhs,velint_old_);

    for (int vi=0;vi<nen_;++vi)
    {
      for (int jdim=0;jdim<nsd_;++jdim)
      {
        velforce(jdim,vi)-=funct_(vi)*scaled_vel(jdim);
      }
    }
  }

  return;
}



/*---------------------------------------------------------------------------------*
 | compute (two!) convection terms of galerkin part               winklmaier 03/12 |
 *---------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::FluidAdjoint3Impl<distype>::ConvectionGalPart(
    LINALG::Matrix<nen_*nsd_,nen_*nsd_> &   estif_w_v,
    LINALG::Matrix<nsd_,nen_> &             velforce,
    const double &                          timefacfac,
    const double &                          timefacfacrhs
) const
{
  /* convection, part 1 */
  /*
            /                          \
           |       /  n         \       |
     rho * | Dv , |  u   o nabla | w    |
           |       \            /   (i) |
            \                          /
  */

  double value = 0.0; // helper

  for (int vi=0; vi<nen_; ++vi)
  {
    const int fvi   = nsd_*vi;

    value = 0.0;

    for (int dim=0;dim<nsd_;++dim)
      value += derxy_(dim,vi)*fluidvelint_(dim);

    value *= fluidAdjoint3Parameter_->Density()*timefacfac;

    for (int ui=0; ui<nen_; ++ui)
    {
      const int fui   = nsd_*ui;

      for (int idim = 0; idim <nsd_; ++idim)
      {
        estif_w_v(fvi+idim,fui+idim) += value*funct_(ui);
      } // end for (idim)
    } //ui
  } // vi

  // rhs at new time step
  for (int vi=0; vi<nen_; ++vi)
  {
    value = 0.0;

    for (int dim=0;dim<nsd_;++dim)
      value += derxy_(dim,vi)*fluidvelint_(dim);

    value *= fluidAdjoint3Parameter_->Density()*timefacfac;

    for (int jdim=0;jdim<nsd_;++jdim)
      velforce(jdim,vi) -= value*velint_(jdim);
  } // vi

  // rhs at old time step
  if (not fluidAdjoint3Parameter_->IsStationary())
  {
    for (int vi=0; vi<nen_; ++vi)
    {
      value = 0.0;

      for (int dim=0;dim<nsd_;++dim)
        value += derxy_(dim,vi)*fluidvelint_old_(dim);

      value *= fluidAdjoint3Parameter_->Density()*timefacfacrhs;

      for (int jdim=0;jdim<nsd_;++jdim)
        velforce(jdim,vi) -= value*velint_old_(jdim);
    } // vi
  }

  /*  convection, part 2 */
  /*
            /                           \
           |  /               \   n      |
     rho * | |  Dv o nabla     | u  , w  |
           |  \           (i) /          |
            \                           /
  */
  for (int ui=0;ui<nen_;++ui)
  {
    value = fluidAdjoint3Parameter_->Density()*timefacfac*funct_(ui);

    for (int idim=0;idim<nsd_;++idim)
    {
      const int fui = nsd_*ui+idim;

      for (int vi=0;vi<nen_;++vi)
      {
        const int fvi = nsd_*vi;

        for (int jdim=0;jdim<nsd_;++jdim)
        {
          estif_w_v(fvi+jdim,fui) += funct_(vi)*fluidvelxy_(idim,jdim)*value;
        }
      }
    }
  }

  // rhs at new time step
  for (int jdim=0;jdim<nsd_;++jdim)
  {
    value = 0.0; // product of fluid and adjoint velocity

    for (int dim=0;dim<nsd_;++dim)
      value+=fluidvelxy_(dim,jdim)*velint_(dim);

    value*=fluidAdjoint3Parameter_->Density()*timefacfac;

    for (int vi=0;vi<nen_;++vi)
      velforce(jdim,vi) -= funct_(vi)*value;
  }

  // rhs at new and old time step
  if (not fluidAdjoint3Parameter_->IsStationary())
  {
    for (int jdim=0;jdim<nsd_;++jdim)
    {
      value = 0.0; // product of fluid and adjoint velocity

      for (int dim=0;dim<nsd_;++dim)
        value+=fluidvelxy_old_(dim,jdim)*velint_old_(dim);

      value*=fluidAdjoint3Parameter_->Density()*timefacfacrhs;

      for (int vi=0;vi<nen_;++vi)
        velforce(jdim,vi) -= funct_(vi)*value;
    }
  }
}



/*---------------------------------------------------------------------------------*
 | compute viscous terms of galerkin part                         winklmaier 03/12 |
 *---------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::FluidAdjoint3Impl<distype>::ViscousGalPart(
    LINALG::Matrix<nen_*nsd_,nen_*nsd_> &   estif_w_v,
    LINALG::Matrix<nsd_,nen_> &             velforce,
    const double &                          timefacfac,
    const double &                          timefacfacrhs
) const
{
  /* viscosity term: overall */
  /*
                   /                        \
                  |       /  \         / \   |
            2 mu  |  eps | Dv | , eps | w |  |
                  |       \  /         \ /   |
                   \                        /
  */

  /* split in two parts: */

  /* viscosity term, part 1 */
  /*
                   /                                    \
                  |         /      \           /     \   |
              mu  |  nabla | Dv     | , nabla | w     |  |
                  |         \  (i) /           \ (i) /   |
                   \                                    /
  */

  /* viscosity term, part 2 */
  /*
                   /                                  \
                  |           /  \           /     \   |
              mu  |  nabla   | Dv | , nabla | w     |  |
                  |       (i) \  /           \ (i) /   |
                   \                                  /
  */

  double value = 0.0; // helper
  const double viscdenstimefac = timefacfac*fluidAdjoint3Parameter_->Viscosity();

  for (int ui=0; ui<nen_; ++ui)
  {
    const int fui   = nsd_*ui;

    for (int vi=0; vi<nen_; ++vi)
    {
      const int fvi   = nsd_*vi;

      value = 0.0;

      for (int dim=0;dim<nsd_;++dim)
        value += derxy_(dim,vi)*derxy_(dim,ui);

      value *= viscdenstimefac;

      for (int idim = 0; idim <nsd_; ++idim)
      {
        estif_w_v(fvi+idim,fui+idim) += value;
      } // end for (idim)
    } //vi
  } // ui

  for (int ui=0; ui<nen_; ++ui)
  {
    const int fui   = nsd_*ui; // shp fcn of u known, derivative not

    for (int jdim=0;jdim<nsd_;++jdim) // derivative of u known by jdim
    {
      value = viscdenstimefac*derxy_(jdim,ui);

      for (int vi=0; vi<nen_; ++vi)
      {
        const int fvi   = nsd_*vi; // shp fcn of v known, derivative not

        for (int idim=0;idim<nsd_;++idim)
        {
          estif_w_v(fvi+jdim,fui+idim) += derxy_(idim,vi)*value;
        }
      } //vi
    }
  } // ui

  // viscosity at new and old time step (if instationary)
  LINALG::Matrix<nsd_,nsd_> viscstress(true);

  double viscdenstimefacrhs = timefacfacrhs*fluidAdjoint3Parameter_->Viscosity(); // for instationary problems
  for (int jdim = 0; jdim < nsd_; ++jdim)
  {
    for (int idim = 0; idim < nsd_; ++idim)
    {
      viscstress(idim,jdim)+=viscdenstimefac*(vderxy_(jdim,idim)+vderxy_(idim,jdim));

      if (not fluidAdjoint3Parameter_->IsStationary())
        viscstress(idim,jdim)+=viscdenstimefacrhs*(vderxy_old_(jdim,idim)+vderxy_old_(idim,jdim));
    }
  }


  // computation of right-hand-side viscosity term
  for (int vi=0; vi<nen_; ++vi)
  {
    for (int idim = 0; idim < nsd_; ++idim)
    {
      for (int jdim = 0; jdim < nsd_; ++jdim)
        velforce(idim,vi) -= viscstress(idim,jdim)*derxy_(jdim,vi);
    }
  }

  return;
}



/*---------------------------------------------------------------------------------*
 | compute pressure term of galerkin part                         winklmaier 03/12 |
 *---------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::FluidAdjoint3Impl<distype>::PressureGalPart(
    LINALG::Matrix<nen_*nsd_,nen_> &          estif_w_q,
    LINALG::Matrix<nsd_,nen_> &               velforce,
    const double &                            timefacfacpre,
    const double &                            timefacfacprerhs,
    const double &                            press,
    const double &                            press_old
) const
{
  /* pressure term */
  /*
       /                \
      |                  |
    + | Dq  , nabla o w  |
      |                  |
       \                /
  */

  double value = 0.0; // helper

  for (int ui=0;ui<nen_;++ui)
  {
    value = timefacfacpre*funct_(ui);

    for (int vi=0;vi<nen_;++vi)
    {
      const int fvi = vi*nsd_;

      for (int jdim=0;jdim<nsd_;++jdim)
        estif_w_q(fvi+jdim,ui) += derxy_(jdim,vi)*value;
    } // vi
  } // ui

// rhs at new time step
  value = timefacfacpre*press;
  for (int vi=0; vi<nen_; ++vi)
  {
    for (int jdim=0;jdim<nsd_;++jdim) // derivative of u known by jdim
      velforce(jdim,vi) -= derxy_(jdim,vi)*value;
  } // vi

  // rhs at old time step
  if (not fluidAdjoint3Parameter_->IsStationary())
  {
    value = timefacfacprerhs*press_old;
    for (int vi=0; vi<nen_; ++vi)
    {
      for (int jdim=0;jdim<nsd_;++jdim) // derivative of u known by jdim
        velforce(jdim,vi) -= derxy_(jdim,vi)*value;
    } // vi
  }

  return;
}



/*---------------------------------------------------------------------------------*
 | compute continuity term of galerkin part                       winklmaier 03/12 |
 *---------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::FluidAdjoint3Impl<distype>::ContinuityGalPart(
    LINALG::Matrix<nen_, nen_*nsd_> &         estif_r_v,
    LINALG::Matrix<nen_,1> &                  preforce,
    const double &                            timefacfacdiv,
    const double &                            timefacfacdivrhs
) const
{
  double value = 0.0; // helper
  /* continuity term */
  /*
       /                \
      |                  |
    - | nabla o Dv  , r  |
      |                  |
       \                /
  */
  for (int vi=0;vi<nen_;++vi)
  {
    value = timefacfacdiv*funct_(vi);

    for (int ui=0;ui<nen_;++ui)
    {
      const int fui = ui*nsd_;

      for (int idim=0;idim<nsd_;++idim)
        estif_r_v(vi,fui+idim) -= value*derxy_(idim,ui);
    } // vi
  } // ui

  // rhs at new and old time step
  value = timefacfacdiv*vdiv_;
  if (not fluidAdjoint3Parameter_->IsStationary())
    value += timefacfacdivrhs*vdiv_old_;

  for (int vi=0; vi<nen_; ++vi)
    preforce(vi) += funct_(vi)*value;


  return;
}


/*---------------------------------------------------------------------------------*
 | compute body force term of galerkin part                        winklmaier 03/12 |
 *---------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::FluidAdjoint3Impl<distype>::BodyForceGalPart(
    LINALG::Matrix<nsd_,nen_> &               velforce,
    const double &                            timefacfac,
    const double &                            timefacfacrhs
) const
{
  double value = 0.0;

  if (fluidAdjoint3Parameter_->TestCase() == INPAR::TOPOPT::adjointtest_no)
  {
    if (fluidAdjoint3Parameter_->ObjDissipationTerm())
    {
      const double dissipation = fluidAdjoint3Parameter_->ObjDissipationFac();

      /*
       *  d   /             \                    /         \
       * --- |   reac*u*u   | (w)   =   2*reac |   u , w   |
       *  du  \             /                    \         /
       */
      for (int idim = 0; idim <nsd_; ++idim)
      {
        value = 2*timefacfac*dissipation*reacoeff_*fluidvelint_(idim);

        if (not fluidAdjoint3Parameter_->IsStationary())
          value += timefacfacrhs*2*dissipation*reacoeff_*fluidvelint_old_(idim);

        for (int vi=0; vi<nen_; ++vi)
        {
          velforce(idim,vi)-=value*funct_(vi);
        }
      }  // end for(idim)


      /*
       *  d   /                    \                  /                    \
       * --- |  2*mu*eps(u)*eps(u)  | (w)   =   4*mu |   eps(u) , nabla w   |
       *  du  \                    /                  \                    /
       */
      for (int idim = 0; idim <nsd_; ++idim)
      {
        for (int jdim = 0; jdim<nsd_; ++jdim)
        {
          value = 2*timefacfac*fluidAdjoint3Parameter_->Viscosity()*(fluidvelxy_(idim,jdim)+fluidvelxy_(jdim,idim));

          if (not fluidAdjoint3Parameter_->IsStationary())
            value = 2*timefacfacrhs*fluidAdjoint3Parameter_->Viscosity()*(fluidvelxy_old_(idim,jdim)+fluidvelxy_old_(jdim,idim));

          for (int vi=0;vi<nen_; ++vi)
          {
            velforce(jdim,vi)+= derxy_(idim,vi)*value;
          }
        }
      }  // end for(idim)
    }
  }
  else // special cases -> no partial integration
  {
    for (int idim = 0; idim <nsd_; ++idim)
    {
      value = timefacfac*bodyforce_(idim);

      if (not fluidAdjoint3Parameter_->IsStationary())
        value += timefacfacrhs*bodyforce_old_(idim);

      for (int vi=0; vi<nen_; ++vi)
      {
        velforce(idim,vi)+=value*funct_(vi);
      }
    }  // end for(idim)
  }

  return;
}



/*---------------------------------------------------------------------------------*
 | compute continuity force term of galerkin part                 winklmaier 03/12 |
 *---------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::FluidAdjoint3Impl<distype>::ContForceGalPart(
    LINALG::Matrix<nen_,1> &                  preforce,
    const double &                            timefacfacdiv,
    const double &                            timefacfacdivrhs
) const
{
  double value = timefacfacdiv*contforce_;

  if (not fluidAdjoint3Parameter_->IsStationary())
    value += timefacfacdivrhs*contforce_old_;

  for (int vi=0; vi<nen_; ++vi)
  {
    preforce(vi,0)-=value*funct_(vi);
  }


  return;
}



/*---------------------------------------------------------------------------------*
 | compute momentum residuum                                      winklmaier 03/12 |
 *---------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::FluidAdjoint3Impl<distype>::DiscreteGalMom(
    LINALG::Matrix<nsd_*nsd_,nen_> &    GalMomTestStat,
    const double &                      timefacfac,
    const double &                      timefacfacrhs,
    const double &                      timefacfacpre,
    const double &                      timefacfacprerhs,
    const LINALG::Matrix<nsd_,nen_>&    eveln,
    const LINALG::Matrix<nsd_,nen_>&    evelnp,
    const LINALG::Matrix<nsd_,nen_>&    efluidveln,
    const LINALG::Matrix<nsd_,nen_>&    efluidvelnp
) const
{
  /*
//      Left hand side terms of Galerkin part for PSPG/SUPG with Dv
//
//    instationary + reactive         /
//                                   |
//      (rho + alpha) w +  dt*Theta  |
//                                   |
//                                    \
//
//                convective term 1             convective term 2 TODO hier anders!!!
//     VZ!!!     /  n             \           /               n    \
//      + rho * |  u o nabla w     | + rho * |   w o nabla   u      |
//               \            (i) /           \               (i)  /
//
//                     viscous term            \
//     VZ!!!    /                         \     |
//      + 2 mu |  nabla o epsilon   ( w )  |    |
//              \                (i)      /     |
//                                             /
  */
  int idim_nsd_p_idim[nsd_];

  for (int idim = 0; idim <nsd_; ++idim)
  {
    idim_nsd_p_idim[idim]=idim*nsd_+idim;
  }

  const double timefacfac_densaf=timefacfac*fluidAdjoint3Parameter_->Density();

  for (int ui=0; ui<nen_; ++ui)
  {
    double value = 0.0;

    for (int dim=0;dim<nsd_;++dim)
      value += fluidvelint_(dim)*derxy_(dim,ui);

    value *= fluidAdjoint3Parameter_->Density()*timefacfac;

    for (int idim = 0; idim <nsd_; ++idim)
    {
      GalMomTestStat(idim_nsd_p_idim[idim],ui)+=value;
    }
  }


// dr_j   d    /    du_j \          du_j         dN_B
// ----= ---- | u_i*----  | = N_B * ---- + u_i * ---- * d_jk
// du_k  du_k  \    dx_i /          dx_k         dx_i

  for (int ui=0; ui<nen_; ++ui)
  {
    const double temp=timefacfac_densaf*funct_(ui);

    for (int idim = 0; idim <nsd_; ++idim)
    {
      const int idim_nsd=idim*nsd_;

      for(int jdim=0;jdim<nsd_;++jdim)
      {
        GalMomTestStat(idim_nsd+jdim,ui)+=temp*fluidvelxy_(idim,jdim);
      }
    }
  }


  const double fac_reac=timefacfac*reacoeff_;

  for (int ui=0; ui<nen_; ++ui)
  {
    const double v=fac_reac*funct_(ui);

    for (int idim = 0; idim <nsd_; ++idim)
    {
      GalMomTestStat(idim_nsd_p_idim[idim],ui)+=v;
    }
  }

  // viscous
  LINALG::Matrix<nsd_,1> viscs(true);
  LINALG::Matrix<nsd_,1> viscs_old(true);
  if (is_higher_order_ele_)
  {
    // prework: evaluate div(eps(v))
    LINALG::Matrix<nsd_*nsd_,nen_> visc_shp(true);
    CalcDivEps(eveln,evelnp,viscs,viscs_old,visc_shp);

    const double v = -2.0*fluidAdjoint3Parameter_->Viscosity()*timefacfac;
    for (int jdim = 0; jdim <nsd_; ++jdim)
    {
      const int nsd_idim=nsd_*jdim;

      for(int idim=0;idim<nsd_;++idim)
      {
        const int nsd_idim_p_jdim=nsd_idim+idim;

        for (int ui=0; ui<nen_; ++ui)
        {
          GalMomTestStat(nsd_idim_p_jdim,ui)+=v*visc_shp(nsd_idim_p_jdim, ui);
        }
      }
    }
  }

  return;
}



template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::FluidAdjoint3Impl<distype>::MomRes(
    LINALG::Matrix<nsd_*nsd_,nen_> &    GalMomResnU,
    LINALG::Matrix<nsd_,1> &            StrongResMomScaled,
    const double &                      timefacfac,
    const double &                      timefacfacrhs,
    const double &                      timefacfacpre,
    const double &                      timefacfacprerhs,
    const LINALG::Matrix<nsd_,nen_>&    eveln,
    const LINALG::Matrix<nsd_,nen_>&    evelnp,
    const LINALG::Matrix<nsd_,nen_>&    efluidveln,
    const LINALG::Matrix<nsd_,nen_>&    efluidvelnp
) const
{
  /*
      Left hand side terms of Galerkin part for PSPG/SUPG with Dv

    instationary + reactive         /
                                   |
      (rho + alpha) Du +  dt*Theta |
                                   |
                                    \

                convective term 1             convective term 2
               /  n              \           /               n \
      - rho * |  u o nabla Dv     | + rho * |  Dv o nabla   u   |
               \             (i) /           \           (i)   /

                     viscous term            \
              /                         \     |
      - 2 mu |  nabla o epsilon   ( Dv ) |    |
              \                (i)      /     |
                                             /
  */

  GalMomResnU.Clear();

  // mass matrix + reaction
  double massreacfac = 0.0; // factor summing up coefficients of reactive term and mass-matrix
  if (fluidAdjoint3Parameter_->IsStationary())
    massreacfac = reacoeff_*timefacfac;
  else
    massreacfac = fluidAdjoint3Parameter_->Density()*fac_+reacoeff_*timefacfac; // fac -> mass matrix // reac*timefacfac/dens -> reactive

  for (int ui=0; ui<nen_; ++ui)
  {
    const double uifunct = massreacfac*funct_(ui);

    for (int idim=0; idim<nsd_; ++idim)
      GalMomResnU(idim*nsd_+idim,ui) += uifunct;
  } // ui

//   convection
  for (int ui=0; ui<nen_; ++ui)
  {
    double value = 0.0;

    for (int dim=0;dim<nsd_;++dim)
      value += fluidvelint_(dim)*derxy_(dim,ui);

    value *= fluidAdjoint3Parameter_->Density()*timefacfac;

    for (int idim = 0; idim <nsd_; ++idim)
      GalMomResnU(idim*nsd_+idim,ui) -= value;
  } //ui

  for (int ui=0;ui<nen_;++ui)
  {
    const double uifunct = fluidAdjoint3Parameter_->Density()*timefacfac*funct_(ui);

    for (int idim=0;idim<nsd_;++idim)
    {
      for (int jdim=0;jdim<nsd_;++jdim)
      {
        GalMomResnU(jdim+idim*nsd_,ui) += fluidvelxy_(jdim,idim)*uifunct;
      }
    }
  }

  // viscous
  LINALG::Matrix<nsd_,1> viscs(true);
  LINALG::Matrix<nsd_,1> viscs_old(true);
  if (is_higher_order_ele_)
  {
    // prework: evaluate div(eps(v))
    LINALG::Matrix<nsd_*nsd_,nen_> visc_shp(true);
    CalcDivEps(eveln,evelnp,viscs,viscs_old,visc_shp);

    // add viscous part
    GalMomResnU.Update(-2.0*timefacfac*fluidAdjoint3Parameter_->Viscosity(),visc_shp,1.0);
  }


  // residuum of momentum equation in strong form
  if (not fluidAdjoint3Parameter_->IsStationary())
  {
    for (int idim=0;idim<nsd_;++idim)
    {
      StrongResMomScaled(idim) = fluidAdjoint3Parameter_->Density()*fac_*(velint_(idim)-velint_old_(idim)) // mass term last iteration
          +timefacfac* // velocity part of last iteration (at t^n) coming
            (fluidAdjoint3Parameter_->Density()*(-conv1_(idim)+conv2_(idim))-2*fluidAdjoint3Parameter_->Viscosity()*viscs(idim)
            +reacoeff_*velint_(idim)-bodyforce_(idim))
          -timefacfacpre*gradp_(idim) // pressure part of last iteration (at t^n)
          +timefacfacrhs* // last time step (= t^n+1) coming
            (fluidAdjoint3Parameter_->Density()*(-conv1_old_(idim)+conv2_old_(idim))-2*fluidAdjoint3Parameter_->Viscosity()*viscs_old(idim)
            +reacoeff_*velint_old_(idim)-bodyforce_old_(idim))
          -timefacfacprerhs*gradp_old_(idim);
    }
  }
  else
  {
    for (int idim=0;idim<nsd_;++idim)
    {
      StrongResMomScaled(idim) = timefacfac*(fluidAdjoint3Parameter_->Density()*(-conv1_(idim)+conv2_(idim))-2*fluidAdjoint3Parameter_->Viscosity()*viscs(idim)
                      +reacoeff_*velint_(idim)-gradp_(idim)-bodyforce_(idim));
    }
  }

  return;
}



/*---------------------------------------------------------------------------------*
 | compute momentum residuum                                      winklmaier 03/12 |
 *---------------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------------*
 | compute divergence of epsilon of v                             winklmaier 03/12 |
 *---------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::FluidAdjoint3Impl<distype>::CalcDivEps(
    const LINALG::Matrix<nsd_,nen_>&      eveln,
    const LINALG::Matrix<nsd_,nen_>&      evelnp,
    LINALG::Matrix<nsd_,1>&               viscs,
    LINALG::Matrix<nsd_,1>&               viscs_old,
    LINALG::Matrix<nsd_*nsd_,nen_>&       visc_shp
) const
{
  /*--- viscous term: div(epsilon(u)) --------------------------------*/
  /*   /                                                \
       |  2 N_x,xx + N_x,yy + N_y,xy + N_x,zz + N_z,xz  |
     1 |                                                |
     - |  N_y,xx + N_x,yx + 2 N_y,yy + N_z,yz + N_y,zz  |
     2 |                                                |
       |  N_z,xx + N_x,zx + N_y,zy + N_z,yy + 2 N_z,zz  |
       \                                                /

       with N_x .. x-line of N
       N_y .. y-line of N                                             */

  if (nsd_==3)
  {
    for (int inode=0; inode<nen_; ++inode)
    {
      double sum = (derxy2_(0,inode)+derxy2_(1,inode)+derxy2_(2,inode));
      visc_shp(0,inode) = 0.5 * (sum + derxy2_(0,inode));
      visc_shp(1,inode) = 0.5 *  derxy2_(3,inode);
      visc_shp(2,inode) = 0.5 *  derxy2_(4,inode);
      visc_shp(3,inode) = 0.5 *  derxy2_(3,inode);
      visc_shp(4,inode) = 0.5 * (sum + derxy2_(1,inode));
      visc_shp(5,inode) = 0.5 *  derxy2_(5,inode);
      visc_shp(6,inode) = 0.5 *  derxy2_(4,inode);
      visc_shp(7,inode) = 0.5 *  derxy2_(5,inode);
      visc_shp(8,inode) = 0.5 * (sum + derxy2_(2,inode));
    }
  }
  else if (nsd_==2)
  {
    for (int inode=0; inode<nen_; ++inode)
    {
      double sum = (derxy2_(0,inode)+derxy2_(1,inode));
      visc_shp(0,inode) = 0.5 * (sum + derxy2_(0,inode));
      visc_shp(1,inode) = 0.5 * derxy2_(2,inode);
      visc_shp(2,inode) = 0.5 * derxy2_(2,inode);
      visc_shp(3,inode) = 0.5 * (sum + derxy2_(1,inode));
    }
  }
  else dserror("Epsilon(N) is not implemented for the 1D case");

  for (int inode=0; inode<nen_; ++inode)
  {
    for (int idim=0; idim<nsd_; ++idim)
    {
      const int nsd_idim = idim*nsd_;

      for (int jdim=0; jdim<nsd_; ++jdim)
        viscs(idim) += visc_shp(nsd_idim+jdim,inode)*evelnp(jdim,inode);
    }
  }

  if (not fluidAdjoint3Parameter_->IsStationary())
  {
    for (int inode=0; inode<nen_; ++inode)
    {
      for (int idim=0; idim<nsd_; ++idim)
      {
        const int nsd_idim = idim*nsd_;

        for (int jdim=0; jdim<nsd_; ++jdim)
          viscs_old(idim) += visc_shp(nsd_idim+jdim,inode)*eveln(jdim,inode);
      }
    }
  }

  return;
}



/*---------------------------------------------------------------------------------*
 | compute PSPG stabilization terms                               winklmaier 03/12 |
 *---------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::FluidAdjoint3Impl<distype>::DiscretePSPG(
    LINALG::Matrix<nen_*nsd_, nen_> &         estif_w_q,
    LINALG::Matrix<nen_,nen_> &               estif_r_q,
    LINALG::Matrix<nsd_,nen_> &               velforce,
    LINALG::Matrix<nen_,1> &                  preforce,
    const LINALG::Matrix<nsd_*nsd_,nen_> &    GalMomTestStat,
    const double &                            timefacfac,
    const double &                            timefacfacrhs,
    const double &                            timefacfacpre,
    const double &                            timefacfacprerhs
) const
{
  const double tau=tau_(1);

  /* pressure stabilisation: inertia if not stationary*/
  /*
              /                  \
             |                    |
             |  rho*Du , nabla q  |
             |                    |
              \                  /
   */
  /* pressure stabilisation: convection, convective part */
  /*
              /                                   \
             |  /       n+1       \                |
             | |   rho*u   o nabla | Du , nabla q  |
             |  \      (i)        /                |
              \                                   /
   */
  /* pressure stabilisation: convection, reactive part if Newton */
  /*
              /                                   \
             |  /                \   n+1           |
             | |   rho*Du o nabla | u     , grad q |
             |  \                /   (i)           |
              \                                   /
   */
  /* pressure stabilisation: reaction if included */
  /*
              /                     \
             |                      |
             |  sigma*Du , nabla q  |
             |                      |
              \                    /
   */
  /* pressure stabilisation: viscosity (-L_visc_u) */
  /*
              /                              \
             |               /  \             |
         mu  |  nabla o eps | Du | , nabla q  |
             |               \  /             |
              \                              /
   */

  // stationary part
  for(int jdim=0;jdim<nsd_;++jdim)
  {
    for (int ui=0; ui<nen_; ++ui)
    {
      const int fui_p_jdim   = nsd_*ui + jdim;

      for(int idim=0;idim<nsd_;++idim)
      {
        const int nsd_idim=nsd_*idim;

        velforce(jdim,ui) -= tau*GalMomTestStat(nsd_idim+jdim,ui)*gradp_(idim);
        if (not fluidAdjoint3Parameter_->IsStationary())
          velforce(jdim,ui) -= tau*timefacfacrhs/timefacfac*GalMomTestStat(nsd_idim+jdim,ui)*gradp_old_(idim);

        for (int vi=0; vi<nen_; ++vi)
        {
          const double temp_vi_idim=derxy_(idim,vi)*tau;

          estif_w_q(fui_p_jdim,vi) += GalMomTestStat(nsd_idim+jdim,ui)*temp_vi_idim;

        } // jdim
      } // vi
    } // ui
  } //idim


  // instationary part
    if (not fluidAdjoint3Parameter_->IsStationary())
    {
      double fac = fluidAdjoint3Parameter_->Density()*fac_; // fac -> mass matrix

      for (int ui=0; ui<nen_; ++ui)
      {
        const double uifunct = fac*funct_(ui);

        for (int idim=0; idim<nsd_; ++idim)
        {
          velforce(idim,ui) -= uifunct*(gradp_(idim)+gradp_old_(idim));

          for (int vi=0; vi<nen_;vi++)
          {
            estif_w_q(ui*nsd_+idim,ui) += uifunct*derxy_(idim,vi);
          }
        }
      } // ui
    }



  /* pressure stabilisation: pressure( L_pres_p) */
  /*
               /                    \
              |                      |
              |  nabla Dp , nabla q  |
              |                      |
               \                    /
   */
  for (int ui=0; ui<nen_; ++ui)
  {
    for (int idim = 0; idim <nsd_; ++idim)
    {
      const double v=timefacfacpre*derxy_(idim,ui)*tau;

      preforce(ui) -= v*gradp_(idim);
      if (not fluidAdjoint3Parameter_->IsStationary())
        preforce(ui) -= v*timefacfacprerhs/timefacfacpre*gradp_old_(idim);

      for (int vi=0; vi<nen_; ++vi)
      {
        estif_r_q(ui,vi)+=v*derxy_(idim,vi);
      } // vi
    } // end for(idim)
  }  // ui



  return;
}



/*---------------------------------------------------------------------------------*
 | compute PSPG stabilization terms                               winklmaier 03/12 |
 *---------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::FluidAdjoint3Impl<distype>::DiscreteSUPG(
    LINALG::Matrix<nen_*nsd_, nen_*nsd_> &    estif_w_v,
    LINALG::Matrix<nen_*nsd_,nen_> &          estif_r_v,
    LINALG::Matrix<nsd_,nen_> &               velforce,
    LINALG::Matrix<nen_,1> &                  preforce,
    const LINALG::Matrix<nsd_*nsd_,nen_> &    GalMomTestStat,
    const double &                            timefacfac,
    const double &                            timefacfacrhs,
    const double &                            timefacfacpre,
    const double &                            timefacfacprerhs
) const
{
  const double tau=tau_(1);
  dserror("not working");
  /* pressure stabilisation: inertia if not stationary*/
  /*
              /                  \
             |                    |
             |  rho*Du , nabla q  |
             |                    |
              \                  /
   */
  /* pressure stabilisation: convection, convective part */
  /*
              /                                   \
             |  /       n+1       \                |
             | |   rho*u   o nabla | Du , nabla q  |
             |  \      (i)        /                |
              \                                   /
   */
  /* pressure stabilisation: convection, reactive part if Newton */
  /*
              /                                   \
             |  /                \   n+1           |
             | |   rho*Du o nabla | u     , grad q |
             |  \                /   (i)           |
              \                                   /
   */
  /* pressure stabilisation: reaction if included */
  /*
              /                     \
             |                      |
             |  sigma*Du , nabla q  |
             |                      |
              \                    /
   */
  /* pressure stabilisation: viscosity (-L_visc_u) */
  /*
              /                              \
             |               /  \             |
         mu  |  nabla o eps | Du | , nabla q  |
             |               \  /             |
              \                              /
   */

  // stationary part
  for(int jdim=0;jdim<nsd_;++jdim)
  {
    for (int ui=0; ui<nen_; ++ui)
    {
      const int fui_p_jdim   = nsd_*ui + jdim;

      for(int idim=0;idim<nsd_;++idim)
      {
        const int nsd_idim=nsd_*idim;

        velforce(jdim,ui) -= tau*GalMomTestStat(nsd_idim+jdim,ui)*gradp_(idim);
        if (not fluidAdjoint3Parameter_->IsStationary())
          velforce(jdim,ui) -= tau*timefacfacrhs/timefacfac*GalMomTestStat(nsd_idim+jdim,ui)*gradp_old_(idim);

        for (int vi=0; vi<nen_; ++vi)
        {
          const double temp_vi_idim=derxy_(idim,vi)*tau;

          estif_w_v(fui_p_jdim,vi) += GalMomTestStat(nsd_idim+jdim,ui)*temp_vi_idim;

        } // jdim
      } // vi
    } // ui
  } //idim


  // instationary part
    if (not fluidAdjoint3Parameter_->IsStationary())
    {
      double fac = fluidAdjoint3Parameter_->Density()*fac_; // fac -> mass matrix

      for (int ui=0; ui<nen_; ++ui)
      {
        const double uifunct = fac*funct_(ui);

        for (int idim=0; idim<nsd_; ++idim)
        {
          velforce(idim,ui) -= uifunct*(gradp_(idim)+gradp_old_(idim));

          for (int vi=0; vi<nen_;vi++)
          {
            estif_r_v(ui*nsd_+idim,ui) += uifunct*derxy_(idim,vi);
          }
        }
      } // ui
    }



  /* pressure stabilisation: pressure( L_pres_p) */
  /*
               /                    \
              |                      |
              |  nabla Dp , nabla q  |
              |                      |
               \                    /
   */
//  for (int ui=0; ui<nen_; ++ui)
//  {
//    for (int idim = 0; idim <nsd_; ++idim)
//    {
//      const double v=timefacfacpre*derxy_(idim,ui)*tau;
//
//      preforce(ui) -= v*gradp_(idim);
//      if (not fluidAdjoint3Parameter_->IsStationary())
//        preforce(ui) -= v*timefacfacprerhs/timefacfacpre*gradp_old_(idim);
//
//      for (int vi=0; vi<nen_; ++vi)
//      {
//        estif_r_q(ui,vi)+=v*derxy_(idim,vi);
//      } // vi
//    } // end for(idim)
//  }  // ui



  return;
}



/*---------------------------------------------------------------------------------*
 | compute PSPG stabilization terms                               winklmaier 03/12 |
 *---------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::FluidAdjoint3Impl<distype>::PSPG(
    LINALG::Matrix<nen_, nen_*nsd_> &         estif_r_v,
    LINALG::Matrix<nen_,nen_> &               estif_r_q,
    LINALG::Matrix<nen_,1> &                  preforce,
    const LINALG::Matrix<nsd_*nsd_,nen_> &    GalMomResnU,
    const LINALG::Matrix<nsd_,1> &            StrongResMomScaled,
    const double &                            timefacfac,
    const double &                            timefacfacrhs,
    const double &                            timefacfacpre,
    const double &                            timefacfacprerhs
) const
{
  const double tau=tau_(1);

  /*
      pressure stabilization

        instationary + reactive
     /                            \
    |                              |
  - |  (rho + alpha) Dv , nabla r  |
    |                              |
     \                            /

                   convective term 1                      convective term 2
             /                          \           /                            \
            |   n                        |         |                 n            |
      + rho |  u o nabla Dv   , nabla r  | - rho * |  Dv o nabla    u  , nabla r  |
            |              (i)           |         |            (i)               |
             \                          /           \                            /

              /      viscous term                  \
             |                                      |
      + 2 mu |  nabla o epsilon   ( Dv ) , nabla r  |
             |                 (i)                  |
              \                                    /
  */

  for(int jdim=0;jdim<nsd_;++jdim)
  {
    for (int ui=0; ui<nen_; ++ui)
    {
      const int fui_p_jdim   = nsd_*ui + jdim;

      for(int idim=0;idim<nsd_;++idim)
      {
        const int nsd_idim=nsd_*idim;

        for (int vi=0; vi<nen_; ++vi)
        {
          estif_r_v(vi,fui_p_jdim) -= tau*derxy_(idim,vi)*GalMomResnU(nsd_idim+jdim,ui);
        } // jdim
      } // vi
    } // ui
  } //idim

  for (int ui=0; ui<nen_; ++ui)
  {
    for (int idim = 0; idim <nsd_; ++idim)
    {
      const double v=tau*timefacfacpre*derxy_(idim,ui);

      for (int vi=0; vi<nen_; ++vi)
      {
        /* pressure stabilisation: pressure( L_pres_p) */
        /*
               /                    \
              |                      |
            + |  nabla Dq , nabla r  |
              |                      |
               \                    /
         */
        estif_r_q(vi,ui)+=v*derxy_(idim,vi);
      } // vi
    } // end for(idim)
  }  // ui

  // rhs for new and old time step
  for (int idim = 0; idim <nsd_; ++idim)
  {
    const double resmom_scaled = -tau*StrongResMomScaled(idim);

    for (int vi=0; vi<nen_; ++vi)
    {
      // pressure stabilisation
      preforce(vi) -= derxy_(idim, vi)*resmom_scaled;
    }
  } // end for(idim)
  return;
}



/*---------------------------------------------------------------------------------*
 | compute SUPG stabilization terms                               winklmaier 03/12 |
 *---------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::FluidAdjoint3Impl<distype>::SUPG(
    LINALG::Matrix<nen_*nsd_,nen_*nsd_> &     estif_w_v,
    LINALG::Matrix<nen_*nsd_,nen_> &          estif_w_q,
    LINALG::Matrix<nsd_,nen_> &               velforce,
    const LINALG::Matrix<nsd_*nsd_,nen_> &    GalMomResnU,
    const LINALG::Matrix<nsd_,1> &            StrongResMomScaled,
    const double &                            timefacfac,
    const double &                            timefacfacrhs,
    const double &                            timefacfacpre,
    const double &                            timefacfacprerhs
) const
{
  /*
     test function
                     /  n             \
  supg_test =   rho |  u o nabla w     |
                     \            (i) /
   */

  double supgfac=-fluidAdjoint3Parameter_->Density()*tau_(0);

  LINALG::Matrix<nen_,1> supg_test(true);
  for (int vi=0; vi<nen_; ++vi)
  {
    for (int dim=0;dim<nsd_;++dim)
    {
      supg_test(vi)+=supgfac*derxy_(dim,vi)*fluidvelint_(dim);
    }
  }

  /*
      supg stabilization

        instationary + reactive
     /                              \
    |                                |
    |  (rho + alpha) Dv , supg_test  |
    |                                |
     \                              /

                   convective term 1                      convective term 2
             /                            \           /                              \
            |   n                          |         |                 n              |
      - rho |  u o nabla Dv   , supg_test  | + rho * |  Dv o nabla    u  , supg_test  |
            |              (i)             |         |            (i)                 |
             \                            /           \                              /

              /      viscous term                    \
             |                                        |
      - 2 mu |  nabla o epsilon   ( Dv ) , supg_test  |
             |                 (i)                    |
              \                                      /
  */

  for (int vi=0; vi<nen_; ++vi)
  {
    for(int idim=0;idim<nsd_;++idim)
    {
      const int nsd_idim=nsd_*idim;

      const int fvi_p_idim = nsd_*vi+idim;

      for(int jdim=0;jdim<nsd_;++jdim)
      {
        const int nsd_idim_p_jdim=nsd_idim+jdim;
        for (int ui=0; ui<nen_; ++ui)
        {
          const int fui_p_jdim   = nsd_*ui + jdim;

          estif_w_v(fvi_p_idim,fui_p_jdim) += supg_test(vi)*GalMomResnU(nsd_idim_p_jdim,ui);
        } // jdim
      } // vi
    } // ui
  } //idim

  /* supg stabilisation: pressure part  ( L_pres_p) */
  /*
              /                      \
             |                        |
           - |  nabla Dq , supg_test  |
             |                        |
              \                      /
   */
  for (int vi=0; vi<nen_; ++vi)
  {
    const double v = timefacfacpre*supg_test(vi);

    for (int idim = 0; idim <nsd_; ++idim)
    {
      const int fvi   = nsd_*vi + idim;

      for (int ui=0; ui<nen_; ++ui)
      {
        estif_w_q(fvi,ui) -= v*derxy_(idim, ui);
      }
    }
  }  // end for(idim)


  // rhs for new and old time step
  for (int idim = 0; idim <nsd_; ++idim)
  {
    for (int vi=0; vi<nen_; ++vi)
      velforce(idim,vi) -= supg_test(vi)*StrongResMomScaled(idim);
  }  // end for(idim)
  return;
}



/*---------------------------------------------------------------------------------*
 | compute residual of continuity equation                        winklmaier 03/12 |
 *---------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::FluidAdjoint3Impl<distype>::ContRes(
    double &                                  StrongResContScaled,
    const double &                            timefacfacdiv,
    const double &                            timefacfacdivrhs
) const
{
  StrongResContScaled = timefacfacdiv*(vdiv_-contforce_);

  if (not fluidAdjoint3Parameter_->IsStationary())
    StrongResContScaled += timefacfacdivrhs*(vdiv_old_-contforce_old_);
}



/*---------------------------------------------------------------------------------*
 | compute div-grad (=continuity) stabilization term              winklmaier 03/12 |
 *---------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::FluidAdjoint3Impl<distype>::ContStab(
    LINALG::Matrix<nen_*nsd_,nen_*nsd_> &     estif_w_v,
    LINALG::Matrix<nsd_,nen_> &               velforce,
    const double &                            timefacfacdiv,
    const double &                            timefacfacdivrhs
) const
{
  double cstabfac = timefacfacdiv*tau_(2);
  double value = 0.0;

  /* continuity stabilisation on left hand side */
  /*
              /                        \
             |                          |
        tauC | nabla o Dv  , nabla o w  |
             |                          |
              \                        /
  */

  for (int ui=0; ui<nen_; ++ui)
  {
    const int fui = nsd_*ui;

    for (int idim = 0; idim <nsd_; ++idim)
    {
      const int fui_p_idim = fui+idim;

      value = cstabfac*derxy_(idim,ui);

      for (int vi=0; vi<nen_; ++vi)
      {
        const int fvi = nsd_*vi;

        for(int jdim=0;jdim<nsd_;++jdim)
        {
          estif_w_v(fvi+jdim,fui_p_idim) += value*derxy_(jdim, vi);
        }
      }
    } // end for(idim)
  }

  double StrongResContScaled = 0.0;

  ContRes(StrongResContScaled,
      timefacfacdiv,
      timefacfacdivrhs);

  // computation of rhs viscosity term at new time step
  for (int vi=0; vi<nen_; ++vi)
  {
    for (int idim = 0; idim < nsd_; ++idim)
    {
      /* viscosity term on right-hand side */
      velforce(idim,vi) -= tau_(2)*derxy_(idim,vi)*StrongResContScaled;
    }
  }

  return;
}



/*---------------------------------------------------------------------------------*
 | extract element data from global vector                        winklmaier 03/12 |
 *---------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::FluidAdjoint3Impl<distype>::ExtractValuesFromGlobalVector(
    const DRT::Discretization&   discretization, ///< discretization
    const std::vector<int>&      lm,             ///<
    LINALG::Matrix<nsd_,nen_> *  matrixtofill,   ///< vector field
    LINALG::Matrix<nen_,1> *     vectortofill,   ///< scalar field
    const std::string            state          ///< state of the global vector
) const
{
  // get state of the global vector
  Teuchos::RCP<const Epetra_Vector> matrix_state = discretization.GetState(state);
  if(matrix_state == Teuchos::null)
    dserror("Cannot get state vector %s", state.c_str());

  // extract local values of the global vectors
  std::vector<double> mymatrix(lm.size());
  DRT::UTILS::ExtractMyValues(*matrix_state,mymatrix,lm);

  for (int inode=0; inode<nen_; ++inode)  // number of nodes
  {
    // fill a vector field via a pointer
    if (matrixtofill != NULL)
    {
      for(int idim=0; idim<nsd_; ++idim) // number of dimensions
      {
        (*matrixtofill)(idim,inode) = mymatrix[idim+(inode*numdofpernode_)];
      }  // end for(idim)
    }
    // fill a scalar field via a pointer
    if (vectortofill != NULL)
      (*vectortofill)(inode,0) = mymatrix[nsd_+(inode*numdofpernode_)];
  }
}
