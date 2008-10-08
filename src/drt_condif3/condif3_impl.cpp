/*----------------------------------------------------------------------*/
/*!
\file condif3_impl.cpp

\brief Internal implementation of Condif3 element

<pre>
Maintainer: Georg Bauer
            bauer@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15252
</pre>
*/
/*----------------------------------------------------------------------*/

#ifdef D_FLUID3
#ifdef CCADISCRET

#include "condif3_impl.H"
#include "condif3_utils.H"
#include "../drt_mat/convecdiffus.H"
#include "../drt_mat/matlist.H"
#include "../drt_lib/drt_timecurve.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_fem_general/drt_utils_fem_shapefunctions.H"
#include "../drt_lib/drt_globalproblem.H"
#include <Epetra_SerialDenseSolver.h>


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Condif3Impl* DRT::ELEMENTS::Condif3Impl::Impl(DRT::ELEMENTS::Condif3* c3)
{
  // we assume here, that numdofpernode is equal for every node within
  // the discretization and does not change during the computations
  const int numdofpernode = c3->NumDofPerNode(*(c3->Nodes()[0]));
  int numscal = numdofpernode;
  if (DRT::Problem::Instance()->ProblemType() == "elch")
    numscal -= 1;

  switch (c3->NumNode())
  {
  case 8:
  {
    static Condif3Impl* f8;
    if (f8==NULL)
      f8 = new Condif3Impl(8,numdofpernode,numscal);
    return f8;
  }
  case 20:
  {
    static Condif3Impl* f20;
    if (f20==NULL)
      f20 = new Condif3Impl(20,numdofpernode,numscal);
    return f20;
  }
  case 27:
  {
    static Condif3Impl* f27;
    if (f27==NULL)
      f27 = new Condif3Impl(27,numdofpernode,numscal);
    return f27;
  }
  case 4:
  {
    static Condif3Impl* f4;
    if (f4==NULL)
      f4 = new Condif3Impl(4,numdofpernode,numscal);
    return f4;
  }
  case 10:
  {
    static Condif3Impl* f10;
    if (f10==NULL)
      f10 = new Condif3Impl(10,numdofpernode,numscal);
    return f10;
  }
  case 6:
  {
    static Condif3Impl* f6;
    if (f6==NULL)
      f6 = new Condif3Impl(6,numdofpernode,numscal);
    return f6;
  }
  case 15:
  {
    static Condif3Impl* f15;
    if (f15==NULL)
      f15 = new Condif3Impl(15,numdofpernode,numscal);
    return f15;
  }
  case 5:
  {
    static Condif3Impl* f5;
    if (f5==NULL)
      f5 = new Condif3Impl(5,numdofpernode,numscal);
    return f5;
  }

  default:
    dserror("node number %d not supported", c3->NumNode());
  }
  return NULL;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Condif3Impl::Condif3Impl(int iel, int numdofpernode, int numscal)
  : iel_(iel),
    numdofpernode_(numdofpernode),
    numscal_(numscal),
    xyze_(3,iel_),
    bodyforce_(iel_*numdofpernode_),
    diffus_(numscal_),
    valence_(numscal_),
    shcacp_(0),
    funct_(iel_),
    densfunct_(iel_),
    deriv_(3,iel_),
    deriv2_(6,iel_),
    xjm_(3,3),
    xij_(3,3),
    derxy_(3,iel_),
    derxy2_(6,iel_),
    rhs_(numdofpernode_),
    hist_(numdofpernode_),
    velint_(3),
    migvelint_(3),
    tau_(numscal_),
    kart_(numscal_),
    xder2_(6,3),
    fac_(0),
    conv_(iel_),
    diff_(iel_),
    mig_(iel_),
    gradpot_(3),
    conint_(numscal_)
{
  return;
}


/*----------------------------------------------------------------------*
 |  calculate system matrix and rhs (public)                 g.bau 08/08|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Condif3Impl::Sysmat(
    const DRT::ELEMENTS::Condif3*   ele, ///< the element those matrix is calculated
    const vector<double>&           ephinp, ///< current scalar field
    const vector<double>&           ehist, ///< rhs from beginning of time step
    const vector<double>&           edens, ///< density*shc
    Epetra_SerialDenseMatrix&       sys_mat,///< element matrix to calculate
    Epetra_SerialDenseMatrix&       sys_mat_sd, ///< subgrid-diff. matrix
    Epetra_SerialDenseVector&       residual, ///< element rhs to calculate
    Epetra_SerialDenseVector&       sugrvisc, ///< subgrid-diff. vector
    const struct _MATERIAL*         material, ///< material pointer
    const double                    time, ///< current simulation time
    const double                    timefac, ///< time discretization factor
    const Epetra_SerialDenseVector& evel, ///< nodal velocities at n+1
    bool                            temperature, ///< temperature flag
    string                          fssgd, ///< subgrid-diff. flag
    const bool                      is_stationary, ///< flag indicating stationary formulation
    const double                    frt ///< factor F/RT needed for ELCH calculations
)
{
  const DRT::Element::DiscretizationType distype = ele->Shape();

  // get node coordinates
  for (int i=0;i<iel_;i++)
  {
    xyze_(0,i)=ele->Nodes()[i]->X()[0];
    xyze_(1,i)=ele->Nodes()[i]->X()[1];
    xyze_(2,i)=ele->Nodes()[i]->X()[2];
  }

  // dead load in element nodes
  BodyForce(ele,time);

  // get diffusivity / diffusivities
  if (material->mattyp == m_matlist)
  {
    for (int k = 0;k<numscal_;++k)
    {
      const int matid = material->m.matlist->matids[k];
      const _MATERIAL& singlemat =  DRT::Problem::Instance()->Material(matid-1);

      if (singlemat.mattyp == m_ion)
      {
        valence_[k]= singlemat.m.ion->valence;
        diffus_[k]= singlemat.m.ion->diffusivity;
      }
      else if (singlemat.mattyp == m_condif)
        diffus_[k]= singlemat.m.condif->diffusivity;
      else
        dserror("material type is not allowed");
    }
    // set specific heat capacity at constant pressure to 1.0
    shcacp_ = 1.0;
  }
  else if (material->mattyp == m_condif)
  {
    dsassert(numdofpernode_==1,"more than 1 dof per node for condif material");

    // in case of a temperature equation, we get thermal conductivity instead of
    // diffusivity and have to divide by the specific heat capacity at constant
    // pressure; otherwise, it is the "usual" diffusivity
    if (temperature)
    {
      shcacp_ = material->m.condif->shc;
      diffus_[0] = material->m.condif->diffusivity/shcacp_;
    }
    else
    {
      // set specific heat capacity at constant pressure to 1.0, get diffusivity
      shcacp_ = 1.0;
      diffus_[0] = material->m.condif->diffusivity;
    }
  }
  else
    dserror("Material type is not supported");

  /*----------------------------------------------------------------------*/
  // el. potential at element nodes (ELCH)
  /*----------------------------------------------------------------------*/
  Epetra_SerialDenseVector epot(iel_);
  if (numdofpernode_-numscal_==1) // only true for ELCH problems
    for (int inode=0;inode<iel_;++inode)
    {
      epot[inode] += ephinp[inode*numdofpernode_+numscal_];
    }

  /*----------------------------------------------------------------------*/
  // calculation of stabilization parameter(s) tau
  /*----------------------------------------------------------------------*/
  CalTau(ele,sugrvisc,evel,epot,distype,timefac,fssgd,is_stationary,frt);

  /*----------------------------------------------------------------------*/
  // integration loop for one condif3 element
  /*----------------------------------------------------------------------*/

  // flag for higher order elements
  const bool higher_order_ele = SCATRA::is3DHigherOrderElement(distype);

  // gaussian points
  const DRT::UTILS::IntegrationPoints3D intpoints(SCATRA::get3DOptimalGaussrule(distype));

  // integration loop
  for (int iquad=0; iquad<intpoints.nquad; ++iquad)
  {
    EvalShapeFuncAndDerivsAtIntPoint(intpoints,iquad,distype,higher_order_ele,ele);

    // density*specific heat capacity-weighted shape functions
    for (int j=0; j<iel_; j++)
    {
      densfunct_[j] = funct_[j]*edens[j];
    }

    // get (density*specific heat capacity-weighted) velocity at element center
    for (int i=0;i<3;i++)
    {
      velint_[i]=0.0;
      for (int j=0;j<iel_;j++)
      {
        velint_[i] += funct_[j]*evel[i+(3*j)];
      }
    }

    /*------------ get values of variables at integration point */
    for (int k = 0;k<numdofpernode_;++k)     // loop of each transported sclar
    {
      // get history data at integration point (weighted by density)
      hist_[k] = 0;
      for (int j=0;j<iel_;j++)
      {
        hist_[k] += densfunct_[j]*ehist[j*numdofpernode_+k];
      }

      // get bodyforce in gausspoint (divided by shcacp for temperature eq.)
      rhs_[k] = 0;
      for (int inode=0;inode<iel_;inode++)
      {
        rhs_[k]+= (1.0/shcacp_)*bodyforce_[inode*numdofpernode_+k]*funct_[inode];
      }
    }

    /*-------------- perform integration for entire matrix and rhs ---*/
    if (numdofpernode_-numscal_== 0) // 'standard' scalar transport
    {
      for (int k=0;k<numscal_;++k) // deal with a system of transported scalars
      {
        if(is_stationary==false)
          CalMat(sys_mat,sys_mat_sd,residual,higher_order_ele,fssgd,timefac,k);
        else
          CalMatStationary(sys_mat,sys_mat_sd,residual,higher_order_ele,fssgd,k);
      } // loop over each scalar
    }
    else  // ELCH problems
      CalMatElch(sys_mat,residual,ephinp,higher_order_ele,frt,is_stationary,timefac);

  } // integration loop

  return;
}


/*----------------------------------------------------------------------*
 |  get the body force  (private)                              gjb 06/08|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Condif3Impl::BodyForce(const DRT::ELEMENTS::Condif3* ele, const double time)
{
  vector<DRT::Condition*> myneumcond;

  // check whether all nodes have a unique VolumeNeumann condition
  DRT::UTILS::FindElementConditions(ele, "VolumeNeumann", myneumcond);

  if (myneumcond.size()>1)
    dserror("more than one VolumeNeumann cond on one node");

  if (myneumcond.size()==1)
  {
    // find out whether we will use a time curve
    const vector<int>* curve  = myneumcond[0]->Get<vector<int> >("curve");
    int curvenum = -1;

    if (curve) curvenum = (*curve)[0];

    // initialisation
    double curvefac    = 0.0;

    if (curvenum >= 0) // yes, we have a timecurve
    {
      // time factor for the intermediate step
      if(time >= 0.0)
      {
        curvefac = DRT::UTILS::TimeCurveManager::Instance().Curve(curvenum).f(time);
      }
      else
      {
        // A negative time value indicates an error.
        dserror("Negative time value in body force calculation: time = %f",time);
      }
    }
    else // we do not have a timecurve --- timefactors are constant equal 1
    {
      curvefac = 1.0;
    }

    // get values and switches from the condition
    const vector<int>*    onoff = myneumcond[0]->Get<vector<int> >   ("onoff");
    const vector<double>* val   = myneumcond[0]->Get<vector<double> >("val"  );

    // set this condition to the bodyforce array
    for (int jnode=0; jnode<iel_; jnode++)
    {
      for(int idof=0;idof<numdofpernode_;idof++)
      {
        bodyforce_(jnode*numdofpernode_+idof) = (*onoff)[idof]*(*val)[idof]*curvefac;
      }
    }
  }
  else
  {
    for (int jnode=0; jnode<iel_; jnode++)
    {
      for(int idof=0;idof<numdofpernode_;idof++)
      {
        // we have no dead load
        bodyforce_(jnode*numdofpernode_+idof) = 0.0;
      }
    }
  }

  return;

} //Condif3Impl::BodyForce


/*----------------------------------------------------------------------*
 |  calculate stabilization parameter  (private)              gjb 06/08 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Condif3Impl::CalTau(
    const DRT::ELEMENTS::Condif3*&          ele,
    Epetra_SerialDenseVector&               sugrvisc,
    const Epetra_SerialDenseVector&         evel,
    const Epetra_SerialDenseVector&         epot,
    const DRT::Element::DiscretizationType& distype,
    const double&                           timefac,
    string                                  fssgd,
    const bool&                             is_stationary,
    const double&                           frt
  )
{
  /*------------------------------------------------------- initialize ---*/
  // use one point gauss rule to calculate tau at element center
  DRT::UTILS::GaussRule3D integrationrule_stabili = DRT::UTILS::intrule3D_undefined;
  switch (distype)
  {
  case DRT::Element::hex8:
  case DRT::Element::hex20:
  case DRT::Element::hex27:
    integrationrule_stabili = DRT::UTILS::intrule_hex_1point;
    break;
  case DRT::Element::tet4:
  case DRT::Element::tet10:
    integrationrule_stabili = DRT::UTILS::intrule_tet_1point;
    break;
  case DRT::Element::wedge6:
  case DRT::Element::wedge15:
    integrationrule_stabili = DRT::UTILS::intrule_wedge_1point;
    break;
  case DRT::Element::pyramid5:
    integrationrule_stabili = DRT::UTILS::intrule_pyramid_1point;
    break;
  default:
    dserror("invalid discretization type for fluid3");
  }

  // gaussian points
  const DRT::UTILS::IntegrationPoints3D  intpoints_tau = getIntegrationPoints3D(integrationrule_stabili);

  // prepare the standard FE stuff for this single integration point
  // we do not need second derivatives for the calculation of tau
  // EvalShapeFuncAndDerivsAtIntPoint(intpoints_tau,0,distype,false,ele);

  // shape functions and derivs at element center
  const double e1    = intpoints_tau.qxg[0][0];
  const double e2    = intpoints_tau.qxg[0][1];
  const double e3    = intpoints_tau.qxg[0][2];
  const double wquad = intpoints_tau.qwgt[0];

  // shape functions and their derivatives
  DRT::UTILS::shape_function_3D(funct_,e1,e2,e3,distype);
  DRT::UTILS::shape_function_3D_deriv1(deriv_,e1,e2,e3,distype);

  /*------------------------------ get Jacobian matrix and determinant ---*/
  /*----------------------------- determine Jacobi Matrix at point r,s ---*/
  double dum;
  for (int i=0; i<3; i++)
  {
    for (int j=0; j<3; j++)
    {
      dum=0.0;
      for (int l=0; l<iel_; l++)
      {
        dum += deriv_(i,l)*xyze_(j,l);
      }
      xjm_(i,j)=dum;
    } /* end of loop j */
  } /* end of loop i */
  // ------------------------------------- calculate Jacobi determinant
  const double det = xjm_(0,0)*xjm_(1,1)*xjm_(2,2)+
                     xjm_(0,1)*xjm_(1,2)*xjm_(2,0)+
                     xjm_(0,2)*xjm_(1,0)*xjm_(2,1)-
                     xjm_(0,2)*xjm_(1,1)*xjm_(2,0)-
                     xjm_(0,0)*xjm_(1,2)*xjm_(2,1)-
                     xjm_(0,1)*xjm_(1,0)*xjm_(2,2);
  if (det < 0.0)
    dserror("GLOBAL ELEMENT NO.%i\nNEGATIVE JACOBIAN DETERMINANT: %f", ele->Id(), det);
  if (abs(det) < 1E-16)
    dserror("GLOBAL ELEMENT NO.%i\nZERO JACOBIAN DETERMINANT: %f", ele->Id(), det);

  const double vol = wquad*det;

  // There exist different definitions for 'the' characteristic element length hk:
  // 1) get element length for tau_Mp/tau_C: volume-equival. diameter
  // const double hk = pow((6.*vol/PI),(1.0/3.0));

  // 2) streamlength (based on velocity vector at element centre)
  if (numdofpernode_-numscal_== 1) // ELCH
  {
    /*------------------------------------------------------------------*/
    /* compute global derivates (only necessary for streamlength & ELCH)*/
    /*------------------------------------------------------------------*/
    /*------------------------------------------------- initialization */
    for(int k=0;k<iel_;k++)
    {
      derxy_(0,k)=0.0;
      derxy_(1,k)=0.0;
      derxy_(2,k)=0.0;
    }

    // ---------------------------------------inverse of transposed jacobian
    double idet = 1./det;
    xij_(0,0) = (  xjm_(1,1)*xjm_(2,2) - xjm_(2,1)*xjm_(1,2))*idet;
    xij_(1,0) = (- xjm_(1,0)*xjm_(2,2) + xjm_(2,0)*xjm_(1,2))*idet;
    xij_(2,0) = (  xjm_(1,0)*xjm_(2,1) - xjm_(2,0)*xjm_(1,1))*idet;
    xij_(0,1) = (- xjm_(0,1)*xjm_(2,2) + xjm_(2,1)*xjm_(0,2))*idet;
    xij_(1,1) = (  xjm_(0,0)*xjm_(2,2) - xjm_(2,0)*xjm_(0,2))*idet;
    xij_(2,1) = (- xjm_(0,0)*xjm_(2,1) + xjm_(2,0)*xjm_(0,1))*idet;
    xij_(0,2) = (  xjm_(0,1)*xjm_(1,2) - xjm_(1,1)*xjm_(0,2))*idet;
    xij_(1,2) = (- xjm_(0,0)*xjm_(1,2) + xjm_(1,0)*xjm_(0,2))*idet;
    xij_(2,2) = (  xjm_(0,0)*xjm_(1,1) - xjm_(1,0)*xjm_(0,1))*idet;

    //---------------------------------------- calculate global derivatives
    for (int k=0;k<iel_;k++)
    {
      derxy_(0,k) +=   xij_(0,0) * deriv_(0,k) + xij_(0,1) * deriv_(1,k) + xij_(0,2) * deriv_(2,k);
      derxy_(1,k) +=   xij_(1,0) * deriv_(0,k) + xij_(1,1) * deriv_(1,k) + xij_(1,2) * deriv_(2,k);
      derxy_(2,k) +=   xij_(2,0) * deriv_(0,k) + xij_(2,1) * deriv_(1,k) + xij_(2,2) * deriv_(2,k);
    }

    // get "migration velocity" divided by D_k*z_k at element center
    for (int i=0;i<3;i++)
    {
      migvelint_[i] = 0.0;
      for (int j=0;j<iel_;j++)
      {
        migvelint_[i] += (-frt)*derxy_(i,j)*epot[j];
      }
    }
  } // if ELCH

/*  double strle = 0.0;
  if (vel_norm>1e-6)
  {
    double val = 0;
    for (int i=0;i<iel_;++i)
    {
      double sum = 0;
      for (int j=0;j<3;++j)
      {
        sum += velint_[j]*derxy_(j,i);
      }
      val+= abs(sum);
    }
    strle = 2.0*vel_norm/val; //this formula is not working in 3D in case of HEX8 elements!!
  }
  else
  //case: 'zero' velocity vector => tau will be very small in diffusion-dominated regions
  // => usage of arbitrary vector velint = (1 0 0)^T in above formula possible
  {
     double val = 0;
     for (int i=0;i<iel_;++i)
     {
       val+=abs(derxy_(0,i));
     }
     strle = 2.0/val;
   }
   //const double hk = strle;*/

  // 3) use cubic root of the element volume as characteristic length
  const double hk = pow(vol,(1.0/3.0));

  // get element type constant for tau
  double mk=0.0;
  switch (distype)
  {
  case DRT::Element::tet4:
  case DRT::Element::pyramid5:
  case DRT::Element::hex8:
  case DRT::Element::wedge6:
    mk = 0.333333333333333333333;
    break;
  case DRT::Element::hex20:
  case DRT::Element::hex27:
  case DRT::Element::tet10:
  case DRT::Element::wedge15:
    mk = 0.083333333333333333333;
    break;
  default:
    dserror("type unknown!\n");
  }

  // get (density*specific heat capacity-weighted) velocity at element center
  for (int i=0;i<3;i++)
  {
    velint_[i]=0.0;
    for (int j=0;j<iel_;j++)
    {
      velint_[i] += funct_[j]*evel[i+(3*j)];
    }
  }

  // some necessary parameter definitions
  double vel_norm, epe1, epe2, xi1, xi2;

  for(int k = 0;k<numscal_;++k) // loop over all transported scalars
  {
    if (numdofpernode_ - numscal_ == 1) // ELCH
    {
      const double Dkzk = diffus_[k]*valence_[k];
      // get Euclidean norm of effective velocity at element center:
      // (weighted) convective velocity + individual migration velocity
      vel_norm = sqrt(DSQR(velint_[0]+Dkzk*migvelint_[0]) + DSQR(velint_[1]+Dkzk*migvelint_[1]) + DSQR(velint_[2]+Dkzk*migvelint_[2]));
    }
    else
    {
      // get Euclidean norm of (weighted) velocity at element center
      vel_norm = sqrt(DSQR(velint_[0]) + DSQR(velint_[1]) + DSQR(velint_[2]));
    }

  // stabilization parameter definition according to Franca and Valentin (2000)
  if (is_stationary == false)
  {
      /* parameter relating diffusive : reactive forces */
      epe1 = 2.0 * timefac * diffus_[k] / (mk * DSQR(hk)); 
      if (diffus_[k] == 0.0) dserror("diffusivity is zero: Preventing division by zero at evaluation of stabilization parameter");
      /* parameter relating convective : diffusive forces */
      epe2 = mk * vel_norm * hk / diffus_[k];
      xi1 = DMAX(epe1,1.0);
      xi2 = DMAX(epe2,1.0);

      tau_[k] = DSQR(hk)/((DSQR(hk)*xi1)/timefac + (2.0*diffus_[k]/mk)*xi2);
  }
  else
  {
      if (diffus_[k] == 0.0) dserror("diffusivity is zero: Preventing division by zero at evaluation of stabilization parameter");
      /* parameter relating convective : diffusive forces */
      epe2 = mk * vel_norm * hk / diffus_[k];
      xi2 = DMAX(epe2,1.0);

      tau_[k] = (DSQR(hk)*mk)/(2.0*diffus_[k]*xi2);
  }

  // compute artificial diffusivity kappa_art_[k]
  if (fssgd == "artificial_all")
  {
      if (diffus_[k] == 0.0) dserror("diffusivity is zero: Preventing division by zero at evaluation of stabilization parameter");
      /* parameter relating convective : diffusive forces */
      epe2 = mk * vel_norm * hk / diffus_[k];
      xi2 = DMAX(epe2,1.0);

      kart_[k] = (DSQR(hk)*mk*DSQR(vel_norm))/(2.0*diffus_[k]*xi2);

      for (int vi=0; vi<iel_; ++vi)
      {
        sugrvisc(vi) = kart_[k]/ele->Nodes()[vi]->NumElement();
      }
  }
  else if (fssgd == "artificial_small" || fssgd == "Smagorinsky_all" ||
           fssgd == "Smagorinsky_small" || fssgd == "scale_similarity" ||
           fssgd == "mixed_Smagorinsky_all" || fssgd == "mixed_Smagorinsky_small")
    dserror("only all-scale artficial diffusivity for convection-diffusion problems possible so far!\n");

  } // loop over scalars

#if 0
   // some debug output
   cout<<"hk (volume equiv diam)            = "<<pow((6.*vol/PI),(1.0/3.0))<<endl;
   cout<<"strle                             = "<<strle<<endl;
   cout<<"hk (cubic root of element volume) = "<<pow(vol,(1.0/3.0))<<endl;
   cout<<"tau = "<<tau<<endl;
#endif

  return;
} //Condif3Impl::Caltau


/*----------------------------------------------------------------------*
 | evaluate shape functions and derivatives at int. point     gjb 08/08 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Condif3Impl::EvalShapeFuncAndDerivsAtIntPoint(
    const DRT::UTILS::IntegrationPoints3D& intpoints, ///< integration points
    const int& iquad,                                 ///< id of current Gauss point
    const DRT::Element::DiscretizationType& distype,  ///< distinguish between DiscretizationType
    const bool& higher_order_ele,                     ///< are second derivatives needed?
    const DRT::ELEMENTS::Condif3*   ele               ///< the element
)
{
  // coordinates of the current integration point
  const double e1 = intpoints.qxg[iquad][0];
  const double e2 = intpoints.qxg[iquad][1];
  const double e3 = intpoints.qxg[iquad][2];

  // shape functions and their first derivatives
  DRT::UTILS::shape_function_3D(funct_,e1,e2,e3,distype);
  DRT::UTILS::shape_function_3D_deriv1(deriv_,e1,e2,e3,distype);

  /*----------------------------------------- compute Jacobian matrix */
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
  double dum;
  /*-------------------------------- determine jacobian at point r,s ---*/
  for (int i=0; i<3; i++)
  {
    for (int j=0; j<3; j++)
    {
      dum=0.0;
      for (int l=0; l<iel_; l++)
      {
        dum += deriv_(i,l)*xyze_(j,l);
      }
      xjm_(i,j)=dum;
    } // end of loop j
  } // end of loop i
  // ---------------------------------------- calculate determinant
  const double det =  xjm_(0,0)*xjm_(1,1)*xjm_(2,2)+
                      xjm_(0,1)*xjm_(1,2)*xjm_(2,0)+
                      xjm_(0,2)*xjm_(1,0)*xjm_(2,1)-
                      xjm_(0,2)*xjm_(1,1)*xjm_(2,0)-
                      xjm_(0,0)*xjm_(1,2)*xjm_(2,1)-
                      xjm_(0,1)*xjm_(1,0)*xjm_(2,2);

  if (det < 0.0)
    dserror("GLOBAL ELEMENT NO.%i\nNEGATIVE JACOBIAN DETERMINANT: %f", ele->Id(), det);
  if (abs(det) < 1E-16)
    dserror("GLOBAL ELEMENT NO.%i\nZERO JACOBIAN DETERMINANT: %f", ele->Id(), det);

  fac_ = intpoints.qwgt[iquad]*det; // Gauss weight * det(J)

  /*------------------------------------------------------------------*/
  /*                                         compute global derivates */
  /*------------------------------------------------------------------*/

  // ---------------------------------------inverse of transposed jacobian
  double idet = 1./det;
  xij_(0,0) = (  xjm_(1,1)*xjm_(2,2) - xjm_(2,1)*xjm_(1,2))*idet;
  xij_(1,0) = (- xjm_(1,0)*xjm_(2,2) + xjm_(2,0)*xjm_(1,2))*idet;
  xij_(2,0) = (  xjm_(1,0)*xjm_(2,1) - xjm_(2,0)*xjm_(1,1))*idet;
  xij_(0,1) = (- xjm_(0,1)*xjm_(2,2) + xjm_(2,1)*xjm_(0,2))*idet;
  xij_(1,1) = (  xjm_(0,0)*xjm_(2,2) - xjm_(2,0)*xjm_(0,2))*idet;
  xij_(2,1) = (- xjm_(0,0)*xjm_(2,1) + xjm_(2,0)*xjm_(0,1))*idet;
  xij_(0,2) = (  xjm_(0,1)*xjm_(1,2) - xjm_(1,1)*xjm_(0,2))*idet;
  xij_(1,2) = (- xjm_(0,0)*xjm_(1,2) + xjm_(1,0)*xjm_(0,2))*idet;
  xij_(2,2) = (  xjm_(0,0)*xjm_(1,1) - xjm_(1,0)*xjm_(0,1))*idet;

  /*------------------------------------------------- initialization */
  for(int k=0;k<iel_;k++)
  {
    derxy_(0,k)=0.0;
    derxy_(1,k)=0.0;
    derxy_(2,k)=0.0;
  }

  // ---------------------------------------- calculate global derivatives
  for (int k=0;k<iel_;k++)
  {
    derxy_(0,k) +=   xij_(0,0) * deriv_(0,k) + xij_(0,1) * deriv_(1,k) + xij_(0,2) * deriv_(2,k);
    derxy_(1,k) +=   xij_(1,0) * deriv_(0,k) + xij_(1,1) * deriv_(1,k) + xij_(1,2) * deriv_(2,k);
    derxy_(2,k) +=   xij_(2,0) * deriv_(0,k) + xij_(2,1) * deriv_(1,k) + xij_(2,2) * deriv_(2,k);
  }

  // ------------------------------------ compute second global derivatives
  if (higher_order_ele) CalSecondDeriv(e1,e2,e3,distype);

  // say goodbye
  return;
}


/*----------------------------------------------------------------------*
 |  calculate second global derivatives w.r.t. x,y,z at point r,s,t
 |                                            (private)      gammi 07/07
 |
 | From the six equations
 |
 |              +-                     -+
 |  d^2N     d  | dx dN   dy dN   dz dN |
 |  ----   = -- | --*-- + --*-- + --*-- |
 |  dr^2     dr | dr dx   dr dy   dr dz |
 |              +-                     -+
 |
 |              +-                     -+
 |  d^2N     d  | dx dN   dy dN   dz dN |
 |  ------ = -- | --*-- + --*-- + --*-- |
 |  ds^2     ds | ds dx   ds dy   ds dz |
 |              +-                     -+
 |
 |              +-                     -+
 |  d^2N     d  | dx dN   dy dN   dz dN |
 |  ----   = -- | --*-- + --*-- + --*-- |
 |  dt^2     dt | dt dx   dt dy   dt dz |
 |              +-                     -+
 |
 |              +-                     -+
 |  d^2N     d  | dx dN   dy dN   dz dN |
 | -----   = -- | --*-- + --*-- + --*-- |
 | ds dr     ds | dr dx   dr dy   dr dz |
 |              +-                     -+
 |
 |              +-                     -+
 |  d^2N     d  | dx dN   dy dN   dz dN |
 | -----   = -- | --*-- + --*-- + --*-- |
 | dt dr     dt | dr dx   dr dy   dr dz |
 |              +-                     -+
 |
 |              +-                     -+
 |  d^2N     d  | dx dN   dy dN   dz dN |
 | -----   = -- | --*-- + --*-- + --*-- |
 | ds dt     ds | dt dx   dt dy   dt dz |
 |              +-                     -+
 |
 | the matrix (jacobian-bar matrix) system
 |
 | +-                                                                                         -+   +-    -+
 | |   /dx\^2        /dy\^2         /dz\^2           dy dx           dz dx           dy dz     |   | d^2N |
 | |  | -- |        | ---|         | ---|          2*--*--         2*--*--         2*--*--     |   | ---- |
 | |   \dr/          \dr/           \dr/             dr dr           dr dr           dr dr     |   | dx^2 |
 | |                                                                                           |   |      |
 | |   /dx\^2        /dy\^2         /dz\^2           dy dx           dz dx           dy dz     |   | d^2N |
 | |  | -- |        | ---|         | ---|          2*--*--         2*--*--         2*--*--     |   | ---- |
 | |   \ds/          \ds/           \ds/             ds ds           ds ds           ds ds     |   | dy^2 |
 | |                                                                                           |   |      |
 | |   /dx\^2        /dy\^2         /dz\^2           dy dx           dz dx           dy dz     |   | d^2N |
 | |  | -- |        | ---|         | ---|          2*--*--         2*--*--         2*--*--     |   | ---- |
 | |   \dt/          \dt/           \dt/             dt dt           dt dt           dt dt     |   | dz^2 |
 | |                                                                                           | * |      |
 | |   dx dx         dy dy          dz dz        dx dy   dx dy   dx dz   dx dz  dy dz   dy dz  |   | d^2N |
 | |   --*--         --*--          --*--        --*-- + --*--   --*-- + --*--  --*-- + --*--  |   | ---- |
 | |   dr ds         dr ds          dr ds        dr ds   ds dr   dr ds   ds dr  dr ds   ds dr  |   | dxdy |
 | |                                                                                           |   |      |
 | |   dx dx         dy dy          dz dz        dx dy   dx dy   dx dz   dx dz  dy dz   dy dz  |   | d^2N |
 | |   --*--         --*--          --*--        --*-- + --*--   --*-- + --*--  --*-- + --*--  |   | ---- |
 | |   dr dt         dr dt          dr dt        dr dt   dt dr   dr dt   dt dr  dr dt   dt dr  |   | dxdz |
 | |                                                                                           |   |      |
 | |   dx dx         dy dy          dz dz        dx dy   dx dy   dx dz   dx dz  dy dz   dy dz  |   | d^2N |
 | |   --*--         --*--          --*--        --*-- + --*--   --*-- + --*--  --*-- + --*--  |   | ---- |
 | |   dt ds         dt ds          dt ds        dt ds   ds dt   dt ds   ds dt  dt ds   ds dt  |   | dydz |
 | +-                                                                                         -+   +-    -+
 |
 |                  +-    -+     +-                           -+
 |                  | d^2N |     | d^2x dN   d^2y dN   d^2y dN |
 |                  | ---- |     | ----*-- + ----*-- + ----*-- |
 |                  | dr^2 |     | dr^2 dx   dr^2 dy   dr^2 dz |
 |                  |      |     |                             |
 |                  | d^2N |     | d^2x dN   d^2y dN   d^2y dN |
 |                  | ---- |     | ----*-- + ----*-- + ----*-- |
 |                  | ds^2 |     | ds^2 dx   ds^2 dy   ds^2 dz |
 |                  |      |     |                             |
 |                  | d^2N |     | d^2x dN   d^2y dN   d^2y dN |
 |                  | ---- |     | ----*-- + ----*-- + ----*-- |
 |                  | dt^2 |     | dt^2 dx   dt^2 dy   dt^2 dz |
 |              =   |      |  -  |                             |
 |                  | d^2N |     | d^2x dN   d^2y dN   d^2y dN |
 |                  | ---- |     | ----*-- + ----*-- + ----*-- |
 |                  | drds |     | drds dx   drds dy   drds dz |
 |                  |      |     |                             |
 |                  | d^2N |     | d^2x dN   d^2y dN   d^2y dN |
 |                  | ---- |     | ----*-- + ----*-- + ----*-- |
 |                  | drdt |     | drdt dx   drdt dy   drdt dz |
 |                  |      |     |                             |
 |                  | d^2N |     | d^2x dN   d^2y dN   d^2z dN |
 |                  | ---- |     | ----*-- + ----*-- + ----*-- |
 |                  | dtds |     | dtds dx   dtds dy   dtds dz |
 |                  +-    -+     +-                           -+
 |
 |
 | is derived. This is solved for the unknown global derivatives.
 |
 |
 |             jacobian_bar * derxy2 = deriv2 - xder2 * derxy
 |                                              |           |
 |                                              +-----------+
 |                                              'chainrulerhs'
 |                                     |                    |
 |                                     +--------------------+
 |                                          'chainrulerhs'
 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Condif3Impl::CalSecondDeriv(
    const double&                            e1,      ///< first coordinate of GP
    const double&                            e2,      ///< second coordinate of GP
    const double&                            e3,      ///< third coordinate of GP
    const DRT::Element::DiscretizationType&  distype  ///< distinguish between DiscretizationType
    )
{
  /*--- get the second derivatives of standard element at current GP */
  DRT::UTILS::shape_function_3D_deriv2(deriv2_,e1,e2,e3,distype);


  /*----------- now we have to compute the second global derivatives */
  // initialize and zero out everything
  static Epetra_SerialDenseMatrix bm(6,6);
  
  /*------------------------------------------------- initialization */
  for(int k=0;k<iel_;k++)
  {
    derxy2_(0,k)=0.0;
    derxy2_(1,k)=0.0;
    derxy2_(2,k)=0.0;
    derxy2_(3,k)=0.0;
    derxy2_(4,k)=0.0;
    derxy2_(5,k)=0.0;
  } /* end of loop over k */

  // calculate elements of jacobian_bar matrix
  bm(0,0) = xjm_(0,0)*xjm_(0,0);
  bm(1,0) = xjm_(1,0)*xjm_(1,0);
  bm(2,0) = xjm_(2,0)*xjm_(2,0);
  bm(3,0) = xjm_(0,0)*xjm_(1,0);
  bm(4,0) = xjm_(0,0)*xjm_(2,0);
  bm(5,0) = xjm_(2,0)*xjm_(1,0);

  bm(0,1) = xjm_(0,1)*xjm_(0,1);
  bm(1,1) = xjm_(1,1)*xjm_(1,1);
  bm(2,1) = xjm_(2,1)*xjm_(2,1);
  bm(3,1) = xjm_(0,1)*xjm_(1,1);
  bm(4,1) = xjm_(0,1)*xjm_(2,1);
  bm(5,1) = xjm_(2,1)*xjm_(1,1);

  bm(0,2) = xjm_(0,2)*xjm_(0,2);
  bm(1,2) = xjm_(1,2)*xjm_(1,2);
  bm(2,2) = xjm_(2,2)*xjm_(2,2);
  bm(3,2) = xjm_(0,2)*xjm_(1,2);
  bm(4,2) = xjm_(0,2)*xjm_(2,2);
  bm(5,2) = xjm_(2,2)*xjm_(1,2);

  bm(0,3) = 2.*xjm_(0,0)*xjm_(0,1);
  bm(1,3) = 2.*xjm_(1,0)*xjm_(1,1);
  bm(2,3) = 2.*xjm_(2,0)*xjm_(2,1);
  bm(3,3) = xjm_(0,0)*xjm_(1,1)+xjm_(1,0)*xjm_(0,1);
  bm(4,3) = xjm_(0,0)*xjm_(2,1)+xjm_(2,0)*xjm_(0,1);
  bm(5,3) = xjm_(1,0)*xjm_(2,1)+xjm_(2,0)*xjm_(1,1);

  bm(0,4) = 2.*xjm_(0,0)*xjm_(0,2);
  bm(1,4) = 2.*xjm_(1,0)*xjm_(1,2);
  bm(2,4) = 2.*xjm_(2,0)*xjm_(2,2);
  bm(3,4) = xjm_(0,0)*xjm_(1,2)+xjm_(1,0)*xjm_(0,2);
  bm(4,4) = xjm_(0,0)*xjm_(2,2)+xjm_(2,0)*xjm_(0,2);
  bm(5,4) = xjm_(1,0)*xjm_(2,2)+xjm_(2,0)*xjm_(1,2);

  bm(0,5) = 2.*xjm_(0,1)*xjm_(0,2);
  bm(1,5) = 2.*xjm_(1,1)*xjm_(1,2);
  bm(2,5) = 2.*xjm_(2,1)*xjm_(2,2);
  bm(3,5) = xjm_(0,1)*xjm_(1,2)+xjm_(1,1)*xjm_(0,2);
  bm(4,5) = xjm_(0,1)*xjm_(2,2)+xjm_(2,1)*xjm_(0,2);
  bm(5,5) = xjm_(1,1)*xjm_(2,2)+xjm_(2,1)*xjm_(1,2);

  /*------------------ determine 2nd derivatives of coord.-functions */

  /*
  |
  |         0 1 2              0...iel-1
  |        +-+-+-+             +-+-+-+-+        0 1 2
  |        | | | | 0           | | | | | 0     +-+-+-+
  |        +-+-+-+             +-+-+-+-+       | | | | 0
  |        | | | | 1           | | | | | 1   * +-+-+-+ .
  |        +-+-+-+             +-+-+-+-+       | | | | .
  |        | | | | 2           | | | | | 2     +-+-+-+
  |        +-+-+-+       =     +-+-+-+-+       | | | | .
  |        | | | | 3           | | | | | 3     +-+-+-+ .
  |        +-+-+-+             +-+-+-+-+       | | | | .
  |        | | | | 4           | | | | | 4   * +-+-+-+ .
  |        +-+-+-+             +-+-+-+-+       | | | | .
  |        | | | | 5           | | | | | 5     +-+-+-+
  |        +-+-+-+             +-+-+-+-+       | | | | iel-1
  |                                            +-+-+-+
  |
  |        xder2               deriv2          xyze^T
  |
  |
  |                                     +-                  -+
  |                                     | d^2x   d^2y   d^2z |
  |                                     | ----   ----   ---- |
  |                                     | dr^2   dr^2   dr^2 |
  |                                     |                    |
  |                                     | d^2x   d^2y   d^2z |
  |                                     | ----   ----   ---- |
  |                                     | ds^2   ds^2   ds^2 |
  |                                     |                    |
  |                                     | d^2x   d^2y   d^2z |
  |                                     | ----   ----   ---- |
  |                                     | dt^2   dt^2   dt^2 |
  |               yields    xder2  =    |                    |
  |                                     | d^2x   d^2y   d^2z |
  |                                     | ----   ----   ---- |
  |                                     | drds   drds   drds |
  |                                     |                    |
  |                                     | d^2x   d^2y   d^2z |
  |                                     | ----   ----   ---- |
  |                                     | drdt   drdt   drdt |
  |                                     |                    |
  |                                     | d^2x   d^2y   d^2z |
  |                                     | ----   ----   ---- |
  |                                     | dsdt   dsdt   dsdt |
  |                                     +-                  -+
  |
  |
  */

  //xder2_ = blitz::sum(deriv2_(i,k)*xyze_(j,k),k);
  for (int i = 0; i < 6; ++i)
  {
      for (int j = 0; j < 3; ++j)
      {
          for (int k = 0; k < iel_; ++k)
          {
              xder2_(i,j) += deriv2_(i,k)*xyze_(j,k);
          }
      }
  }

  /*
  |        0...iel-1             0 1 2
  |        +-+-+-+-+            +-+-+-+
  |        | | | | | 0          | | | | 0
  |        +-+-+-+-+            +-+-+-+            0...iel-1
  |        | | | | | 1          | | | | 1         +-+-+-+-+
  |        +-+-+-+-+            +-+-+-+           | | | | | 0
  |        | | | | | 2          | | | | 2         +-+-+-+-+
  |        +-+-+-+-+       =    +-+-+-+       *   | | | | | 1 * (-1)
  |        | | | | | 3          | | | | 3         +-+-+-+-+
  |        +-+-+-+-+            +-+-+-+           | | | | | 2
  |        | | | | | 4          | | | | 4         +-+-+-+-+
  |        +-+-+-+-+            +-+-+-+
  |        | | | | | 5          | | | | 5          derxy
  |        +-+-+-+-+            +-+-+-+
  |
  |       chainrulerhs          xder2
  */

  //derxy2_ = -blitz::sum(xder2_(i,k)*derxy_(k,j),k);
  //derxy2_ = deriv2 - blitz::sum(xder2(i,k)*derxy(k,j),k);
  for (int i = 0; i < 6; ++i)
  {
      for (int j = 0; j < iel_; ++j)
      {
          derxy2_(i,j) += deriv2_(i,j);
          for (int k = 0; k < 3; ++k)
          {
              derxy2_(i,j) -= xder2_(i,k)*derxy_(k,j);
          }
      }
  }

  /*
  |        0...iel-1            0...iel-1         0...iel-1
  |        +-+-+-+-+            +-+-+-+-+         +-+-+-+-+
  |        | | | | | 0          | | | | | 0       | | | | | 0
  |        +-+-+-+-+            +-+-+-+-+         +-+-+-+-+
  |        | | | | | 1          | | | | | 1       | | | | | 1
  |        +-+-+-+-+            +-+-+-+-+         +-+-+-+-+
  |        | | | | | 2          | | | | | 2       | | | | | 2
  |        +-+-+-+-+       =    +-+-+-+-+    +    +-+-+-+-+
  |        | | | | | 3          | | | | | 3       | | | | | 3
  |        +-+-+-+-+            +-+-+-+-+         +-+-+-+-+
  |        | | | | | 4          | | | | | 4       | | | | | 4
  |        +-+-+-+-+            +-+-+-+-+         +-+-+-+-+
  |        | | | | | 5          | | | | | 5       | | | | | 5
  |        +-+-+-+-+            +-+-+-+-+         +-+-+-+-+
  |
  |       chainrulerhs         chainrulerhs        deriv2
  */

  //derxy2_ += deriv2_;

  /* make LR decomposition and solve system for all right hand sides
   * (i.e. the components of chainrulerhs)
  |
  |          0  1  2  3  4  5         i        i
  |        +--+--+--+--+--+--+       +-+      +-+
  |        |  |  |  |  |  |  | 0     | | 0    | | 0
  |        +--+--+--+--+--+--+       +-+      +-+
  |        |  |  |  |  |  |  | 1     | | 1    | | 1
  |        +--+--+--+--+--+--+       +-+      +-+
  |        |  |  |  |  |  |  | 2     | | 2    | | 2
  |        +--+--+--+--+--+--+    *  +-+   =  +-+      for i=0...iel-1
  |        |  |  |  |  |  |  | 3     | | 3    | | 3
  |        +--+--+--+--+--+--+       +-+      +-+
  |        |  |  |  |  |  |  | 4     | | 4    | | 4
  |        +--+--+--+--+--+--+       +-+      +-+
  |        |  |  |  |  |  |  | 5     | | 5    | | 5
  |        +--+--+--+--+--+--+       +-+      +-+
  |                                   |        |
  |                                   |        |
  |                                   derxy2[i]|
  |                                    |
  |                                    chainrulerhs[i]
  |
  |   yields
  |
  |                      0...iel-1
  |                      +-+-+-+-+
  |                      | | | | | 0 = drdr
  |                      +-+-+-+-+
  |                      | | | | | 1 = dsds
  |                      +-+-+-+-+
  |                      | | | | | 2 = dtdt
  |            derxy2 =  +-+-+-+-+
  |                      | | | | | 3 = drds
  |                      +-+-+-+-+
  |                      | | | | | 4 = drdt
  |                      +-+-+-+-+
  |                      | | | | | 5 = dsdt
  |                  +-+-+-+-+
  */

  Epetra_SerialDenseSolver solver;
  solver.SetMatrix(bm);

  // No need for a separate rhs. We assemble the rhs to the solution
  // vector. The solver will destroy the rhs and return the solution.
  solver.SetVectors(derxy2_,derxy2_);
  solver.Solve();

  return;
} //Condif3Impl::CalSecondDeriv


/*----------------------------------------------------------------------*
 |  evaluate instationary convection-diffusion matrix (private)gjb 06/08|
 *----------------------------------------------------------------------*/

/*
In this routine the Gauss point contributions to the elemental coefficient
matrix of a stabilized condif3 element are calculated for the instationary
case. The procedure is based on the Rothe method of first discretizing in
time. Hence the resulting terms include coefficients containing time
integration variables such as theta or delta t which are represented by
'timefac'.

The stabilization is based on the residuum:

R = rho * c_p * phi + timefac * rho * c_p * u * grad(phi)
                    - timefac * diffus * laplace(phi) - rhsint

The corresponding weighting operators are
L = timefac * rho * c_p * u * grad(w) +/- timefac * diffus * laplace(w)

'+': USFEM (default)
'-': GLS


The calculation proceeds as follows.
1) obtain single operators of R and L
2) build Galerkin terms from them
3) build stabilizing terms from them
4) build Galerkin and stabilizing terms of RHS

NOTE: Galerkin and stabilization matrices are calculated within one
      routine.


for further comments see comment lines within code.

</pre>
\param **estif      DOUBLE        (o)   ele stiffness matrix
\param  *eforce     DOUBLE        (o)   ele force vector
\return void
------------------------------------------------------------------------*/

void DRT::ELEMENTS::Condif3Impl::CalMat(
    Epetra_SerialDenseMatrix& estif,
    Epetra_SerialDenseMatrix& esd,
    Epetra_SerialDenseVector& eforce,
    const bool                higher_order_ele,
    string                    fssgd,
    const double&             timefac,
    const int&                dofindex
    )
{
static double             rhsint;           /* rhs at int. point     */

// stabilization parameter
const double taufac = tau_[dofindex]*fac_;

// integration factors and coefficients of single terms
const double timefacfac  = timefac * fac_;
const double timetaufac  = timefac * taufac;

/*-------------------------------- evaluate rhs at integration point ---*/
rhsint = hist_[dofindex] + rhs_[dofindex]*timefac;

for (int i=0; i<iel_; i++)
{
   /* convective part */
   /* rho * c_p * u_x * N,x  +  rho * c_p * u_y * N,y +  rho * c_p * u_z * N,z
      with  N .. form function matrix */
   conv_[i] = velint_[0] * derxy_(0,i) + velint_[1] * derxy_(1,i) + velint_[2] * derxy_(2,i);
}

if (higher_order_ele)
{
  for (int i=0; i<iel_; i++)
  {
    /* diffusive part */
    /* diffus * ( N,xx  +  N,yy +  N,zz ) */
    diff_[i] = diffus_[dofindex] * (derxy2_(0,i) + derxy2_(1,i) + derxy2_(2,i));
  }
}

/*--------------------------------- now build single stiffness terms ---*/
const int numdof =numdofpernode_;
// -------------------------------------------System matrix
for (int vi=0; vi<iel_; ++vi)
{
  for (int ui=0; ui<iel_; ++ui)
  {
    /* Standard Galerkin terms: */
    /* transient term */
    estif(vi*numdof+dofindex, ui*numdof+dofindex) += fac_*funct_[vi]*densfunct_[ui] ;

    /* convective term */
    estif(vi*numdof+dofindex, ui*numdof+dofindex) += timefacfac*funct_[vi]*conv_[ui] ;

    /* diffusive term */
    estif(vi*numdof+dofindex, ui*numdof+dofindex) += timefacfac*diffus_[dofindex]*(derxy_(0, ui)*derxy_(0, vi) + derxy_(1, ui)*derxy_(1, vi) + derxy_(2, ui)*derxy_(2, vi)) ;

    /* Stabilization terms: */
    /* 1) transient stabilization (USFEM assumed here, sign change necessary for GLS) */
    /* transient term */
    //estif(vi, ui) += -taufac*densfunct_[vi]*densfunct_[ui] ;

    /* convective term */
    //estif(vi, ui) += -timetaufac*densfunct_[vi]*conv_[ui] ;

    /* diffusive term */
    //if (higher_order_ele) estif(vi, ui) += timetaufac*densfunct_[vi]*diff[ui] ;

    /* 2) convective stabilization */
    /* transient term */
    estif(vi*numdof+dofindex, ui*numdof+dofindex) += taufac*conv_[vi]*densfunct_[ui];

    /* convective term */
    estif(vi*numdof+dofindex, ui*numdof+dofindex) += timetaufac*conv_[vi]*conv_[ui] ;

    if (higher_order_ele)
    {
      /* diffusive term */
      estif(vi*numdof+dofindex, ui*numdof+dofindex) += -timetaufac*conv_[vi]*diff_[ui] ;

      /* 2) diffusive stabilization (USFEM assumed here, sign change necessary for GLS) */
      /* transient term */
      estif(vi*numdof+dofindex, ui*numdof+dofindex) += taufac*diff_[vi]*densfunct_[ui] ;

      /* convective term */
      estif(vi*numdof+dofindex, ui*numdof+dofindex) += timetaufac*diff_[vi]*conv_[ui] ;

      /* diffusive term */
      estif(vi*numdof+dofindex, ui*numdof+dofindex) += -timetaufac*diff_[vi]*diff_[ui] ;
    }
  }
}

// ----------------------------------------------RHS
for (int vi=0; vi<iel_; ++vi)
{
  /* RHS source term */
  eforce[vi*numdof+dofindex] += fac_*funct_[vi]*rhsint ;

  /* transient stabilization of RHS source term */
  //eforce[vi] += -taufac*densfunct_[vi]*rhsint ;

  /* convective stabilization of RHS source term */
  eforce[vi*numdof+dofindex] += taufac*conv_[vi]*rhsint ;

  /* diffusive stabilization of RHS source term */
  if (higher_order_ele) eforce[vi*numdof+dofindex] += taufac*diff_[vi]*rhsint ;
}

// ---------------------------------------artifical-diffusivity matrix
if (fssgd != "No")
{
  // parameter for artificial diffusivity (scaled to one currently)
  const double kartfac = timefacfac;
  //const double kartfac = kart_[dofindex]*timefacfac;

  for (int vi=0; vi<iel_; ++vi)
  {
    for (int ui=0; ui<iel_; ++ui)
    {
      /* artificial diffusivity term */
      esd(vi*numdof+dofindex, ui*numdof+dofindex) += kartfac*(derxy_(0, vi)*derxy_(0, ui) + derxy_(1, vi)*derxy_(1, ui) + derxy_(2, vi)*derxy_(2, ui)) ;

      /*subtract SUPG term */
      //esd(vi*numdof+dofindex, ui*numdof+dofindex) -= taufac*conv[vi]*conv[ui] ;
    }
  }
}

return;
} //Condif3Impl::Condif3CalMat


/*----------------------------------------------------------------------*
 |  evaluate stationary convection-diffusion matrix (private)  gjb 06/08|
 *----------------------------------------------------------------------*/

/*
In this routine the Gauss point contributions to the elemental coefficient
matrix of a stabilized condif3 element are calculated for the stationary
case.

The stabilization is based on the residuum:

R = rho * c_p * u * grad(phi) - diffus *  laplace(phi) - rhsint

The corresponding weighting operators are
L = rho * c_p * u * grad(w) +/- diffus *  laplace(w)

'+': USFEM (default)
'-': GLS


The calculation proceeds as follows.
1) obtain single operators of R and L
2) build Galerkin terms from them
3) build stabilizing terms from them
4) build Galerkin and stabilizing terms of RHS

NOTE: Galerkin and stabilization matrices are calculated within one
      routine.


for further comments see comment lines within code.

</pre>
\param **estif      DOUBLE        (o)   ele stiffness matrix
\param  *eforce     DOUBLE        (o)   ele force vector
\return void
------------------------------------------------------------------------*/

void DRT::ELEMENTS::Condif3Impl::CalMatStationary(
    Epetra_SerialDenseMatrix& estif,
    Epetra_SerialDenseMatrix& esd,
    Epetra_SerialDenseVector& eforce,
    const bool                higher_order_ele,
    string                    fssgd,
    const int&                dofindex
    )
{
static double rhsint;           /* rhs at int. point     */
const double  fac_diffus = fac_*diffus_[dofindex];

// stabilization parameter
const double taufac = tau_[dofindex]*fac_;

/*------------------------------------- set rhs at integration point ---*/
rhsint = rhs_[dofindex];

for (int i=0; i<iel_; i++)
{
   /* convective part */
   /* rho * c_p * u_x * N,x  +  rho * c_p * u_y * N,y +  rho * c_p * u_z * N,z
      with  N .. form function matrix */
   conv_[i] = velint_[0] * derxy_(0,i) + velint_[1] * derxy_(1,i) + velint_[2] * derxy_(2,i);
}

if (higher_order_ele)
{
  for (int i=0; i<iel_; i++)
  {
    /* diffusive part */
    /* diffus * ( N,xx  +  N,yy +  N,zz ) */
    diff_[i] = diffus_[dofindex] * (derxy2_(0,i) + derxy2_(1,i) + derxy2_(2,i));
  }
}

/*--------------------------------- now build single stiffness terms ---*/
const int numdof = numdofpernode_;
// -------------------------------------------System matrix
for (int vi=0; vi<iel_; ++vi)
{
  for (int ui=0; ui<iel_; ++ui)
  {
    /* Standard Galerkin terms: */
    /* convective term */
    estif(vi*numdof+dofindex, ui*numdof+dofindex) += fac_*funct_[vi]*conv_[ui] ;

    /* diffusive term */
    estif(vi*numdof+dofindex, ui*numdof+dofindex) += fac_diffus*(derxy_(0, ui)*derxy_(0, vi) + derxy_(1, ui)*derxy_(1, vi)+ derxy_(2, ui)*derxy_(2, vi));

    /* Stabilization term: */
    /* 1) convective stabilization */

    /* convective term */
    estif(vi*numdof+dofindex, ui*numdof+dofindex) += taufac*conv_[vi]*conv_[ui] ;

    if (higher_order_ele)
    {
      /* diffusive term */
	  estif(vi*numdof+dofindex, ui*numdof+dofindex) += -taufac*conv_[vi]*diff_[ui] ;

      /* 2) diffusive stabilization (USFEM assumed here, sign change necessary for GLS) */

      /* convective term */
      estif(vi*numdof+dofindex, ui*numdof+dofindex) += taufac*diff_[vi]*conv_[ui] ;

      /* diffusive term */
      estif(vi*numdof+dofindex, ui*numdof+dofindex) -= taufac*diff_[vi]*diff_[ui] ;
    }
  }
}

// ----------------------------------------------RHS
for (int vi=0; vi<iel_; ++vi)
{
  /* RHS source term */
  eforce[vi*numdof+dofindex] += fac_*funct_[vi]*rhsint ;

  /* convective stabilization of RHS source term */
  eforce[vi*numdof+dofindex] += taufac*conv_[vi]*rhsint ;

  /* diffusive stabilization of RHS source term */
  if (higher_order_ele) eforce[vi*numdof+dofindex] += taufac*diff_[vi]*rhsint ;
}

// ---------------------------------------artifical-diffusivity matrix
if (fssgd != "No")
{
  // parameter for artificial diffusivity (scaled to one currently)
  const double kartfac = fac_;
  //const double kartfac = kart_[dofindex]*fac_;

  for (int vi=0; vi<iel_; ++vi)
  {
    for (int ui=0; ui<iel_; ++ui)
    {
      /* artificial diffusivity term */
      esd(vi*numdof+dofindex, ui*numdof+dofindex) += kartfac*(derxy_(0, vi)*derxy_(0, ui) + derxy_(1, vi)*derxy_(1, ui) + derxy_(2, vi)*derxy_(2, ui)) ;

      /*subtract SUPG term */
      //esd(vi*numdof+dofindex, ui*numdof+dofindex) -= taufac*conv[vi]*conv[ui] ;
    }
  }
}

return;
} //Condif3Impl::Condif3CalMatStationary


/*----------------------------------------------------------------------*
 | calculate mass matrix and rhs for initializing OST          gjb 08/08|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Condif3Impl::InitializeOST(
    const DRT::ELEMENTS::Condif3*   ele,
    const vector<double>&           ephi0,
    const vector<double>&           edens,
    Epetra_SerialDenseMatrix&       massmat,
    Epetra_SerialDenseVector&       rhs,
    Epetra_SerialDenseVector&       sugrvisc,
    const struct _MATERIAL*         material,
    const double                    time,
    const double                    timefac,
    const Epetra_SerialDenseVector& evel,
    bool                            temperature,
    string                          fssgd
    )
{
  const DRT::Element::DiscretizationType distype = ele->Shape();

  // get node coordinates
  for (int i=0;i<iel_;i++)
  {
    xyze_(0,i)=ele->Nodes()[i]->X()[0];
    xyze_(1,i)=ele->Nodes()[i]->X()[1];
    xyze_(2,i)=ele->Nodes()[i]->X()[2];
  }

  // dead load in element nodes
  BodyForce(ele,time);

  // get diffusivity / diffusivities
  if (material->mattyp == m_matlist)
  {
    for (int k = 0;k<numscal_;++k)
    {
      const int matid = material->m.matlist->matids[k];
      const _MATERIAL& singlemat =  DRT::Problem::Instance()->Material(matid-1);

      if (singlemat.mattyp == m_ion)
      {
        valence_[k]= singlemat.m.ion->valence;
        diffus_[k]= singlemat.m.ion->diffusivity;
      }
      else if (singlemat.mattyp == m_condif)
        diffus_[k]= singlemat.m.condif->diffusivity;
      else
        dserror("material type is not allowed");
#if 0
      cout<<"MatId: "<<material->m.matlist->matids[k]<<"diffusivity["<<k<<"] = "<<diffus[k]<<endl;
#endif
    }
    // set specific heat capacity at constant pressure to 1.0
    shcacp_ = 1.0;
  }
  else if (material->mattyp == m_condif)
  {
    dsassert(numdofpernode_==1,"more than 1 dof per node for condif material");

    // in case of a temperature equation, we get thermal conductivity instead of
    // diffusivity and have to divide by the specific heat capacity at constant
    // pressure; otherwise, it is the "usual" diffusivity
    if (temperature)
    {
      shcacp_ = material->m.condif->shc;
      diffus_[0] = material->m.condif->diffusivity/shcacp_;
    }
    else 
    {
      shcacp_ = 1.0;
      diffus_[0] = material->m.condif->diffusivity;
    }
  }
  else
    dserror("Material type is not supported");

  /*----------------------------------------------------------------------*/
  // el. potential at element nodes
  /*----------------------------------------------------------------------*/
  Epetra_SerialDenseVector epot(iel_); //dummy
  double frt = 1.0; //dummy
  
  /*----------------------------------------------------------------------*/
  // calculation of instationary(!) stabilization parameter(s)
  /*----------------------------------------------------------------------*/
  CalTau(ele,sugrvisc,evel,epot,distype,timefac,fssgd,false,frt);

  /*----------------------------------------------------------------------*/
  // integration loop for one condif3 element
  /*----------------------------------------------------------------------*/

  // flag for higher order elements
  const bool higher_order_ele = SCATRA::is3DHigherOrderElement(distype);

  // gaussian points
  const DRT::UTILS::IntegrationPoints3D intpoints(SCATRA::get3DOptimalGaussrule(distype));

  // integration loop
  for (int iquad=0; iquad<intpoints.nquad; ++iquad)
  {
    EvalShapeFuncAndDerivsAtIntPoint(intpoints,iquad,distype,higher_order_ele,ele);

    // density*specific heat capacity-weighted shape functions
    for (int j=0; j<iel_; j++)
    {
      densfunct_[j] = funct_[j]*edens[j];
    }

    // get (density*specific heat capacity-weighted) velocity at element center
    for (int i=0;i<3;i++)
    {
      velint_[i]=0.0;
      for (int j=0;j<iel_;j++)
      {
        velint_[i] += funct_[j]*evel[i+(3*j)];
      }
    }

    vector<double> phi0(numdofpernode_);

    /*------------ get values of variables at integration point */
    for (int k = 0;k<numscal_;++k)     // loop of each transported sclar
    {
      // get bodyforce in gausspoint (divided by shcacp for temperature eq.)
      rhs_[k] = 0;
      for (int inode=0;inode<iel_;inode++)
      {
        rhs_[k]+= (1.0/shcacp_)*bodyforce_[inode*numscal_+k]*funct_[inode];
        // note: bodyforce calculation isn't filled with functionality yet.
        // -> this line has no effect since bodyforce is always zero.
      }
    }

    /*-------------- perform integration for entire matrix and rhs ---*/
    for (int dofindex=0;dofindex<numscal_;++dofindex) // deal with a system of transported scalars
    {
      static double             rhsint;            /* rhs at int. point     */

      // stabilization parameter
      const double taufac = tau_[dofindex]*fac_;

      /*-------------------------------- evaluate rhs at integration point ---*/
      rhsint = rhs_[dofindex];

      for (int i=0; i<iel_; i++) /* loop over nodes of element */
      {
         /* convective part */
         /* rho * c_p * u_x * N,x  +  rho * c_p * u_y * N,y
            with  N .. form function matrix */
         conv_[i] = velint_[0] * derxy_(0,i) + velint_[1] * derxy_(1,i) + velint_[2] * derxy_(2,i);
      } // end of loop over nodes of element

      if (higher_order_ele)
      {
        for (int i=0; i<iel_; i++) /* loop over nodes of element */
        {
           /* diffusive part */
           /* diffus * ( N,xx  +  N,yy +  N,zz ) */
           diff_[i] = diffus_[dofindex] * (derxy2_(0,i) + derxy2_(1,i) + derxy2_(2,i));
        } // end of loop over nodes of element
      }

      /*--------------------------------- now build single stiffness terms ---*/
      const int numdof = numdofpernode_;
      // -------------------------------------------System matrix
      for (int vi=0; vi<iel_; ++vi)
      {
          for (int ui=0; ui<iel_; ++ui)
          {
          /* Standard Galerkin terms: */
          /* transient term */
          massmat(vi*numdof+dofindex, ui*numdof+dofindex) += fac_*funct_[vi]*densfunct_[ui] ;

          /* convective term */
          rhs[vi*numdof+dofindex] += -(fac_*funct_[vi]*conv_[ui]*ephi0[ui*numdof+dofindex]) ;

          /* diffusive term */
          rhs[vi*numdof+dofindex] += -(fac_*diffus_[dofindex]*(derxy_(0, ui)*derxy_(0, vi) + derxy_(1, ui)*derxy_(1, vi) + derxy_(2, ui)*derxy_(2, vi))*ephi0[ui*numdof+dofindex]);

          /* Stabilization terms: */
          /* 1) transient stabilization (USFEM assumed here, sign change necessary for GLS) */
          /* transient term */
          //massmat(vi, ui) += -taufac*densfunct_[vi]*densfunct_[ui] ;

          /* convective term */
          //massmat(vi, ui) += -timetaufac*densfunct_[vi]*conv[ui] ;

          /* diffusive term */
          //massmat(vi, ui) += timetaufac*densfunct_[vi]*diff[ui] ;

          /* 2) convective stabilization */
          /* transient term */
          massmat(vi*numdof+dofindex, ui*numdof+dofindex) += taufac*conv_[vi]*densfunct_[ui]*ephi0[ui*numdof+dofindex];

          /* convective term */
          rhs[vi*numdof+dofindex] += -(taufac*conv_[vi]*conv_[ui]*ephi0[ui*numdof+dofindex] );

          if (higher_order_ele)
          {
            /* diffusive term */
            rhs[vi*numdof+dofindex] += -(-taufac*conv_[vi]*diff_[ui]*ephi0[ui*numdof+dofindex] );

            /* 2) diffusive stabilization (USFEM assumed here, sign change necessary for GLS) */
            /* transient term */
            massmat(vi*numdof+dofindex, ui*numdof+dofindex) += taufac*diff_[vi]*densfunct_[ui];

            /* convective term */
            rhs[vi*numdof+dofindex] += -(taufac*diff_[vi]*conv_[ui]*ephi0[ui*numdof+dofindex]); 

            /* diffusive term */
            rhs[vi*numdof+dofindex] += -(-taufac*diff_[vi]*diff_[ui]*ephi0[ui*numdof+dofindex]);
          }
        }
      }

      // ----------------------------------------------RHS
      for (int vi=0; vi<iel_; ++vi)
      {
        /* RHS source term */
        rhs[vi*numdof+dofindex] += fac_*funct_[vi]*rhsint ;

        /* transient stabilization of RHS source term */
        //eforce[vi] += -taufac*densfunct[vi]*rhsint ;

        /* convective stabilization of RHS source term */
        rhs[vi*numdof+dofindex] += taufac*conv_[vi]*rhsint ;

        /* diffusive stabilization of RHS source term */
        rhs[vi*numdof+dofindex] += taufac*diff_[vi]*rhsint ;
      }

      // ------------------------------artifical-diffusivity matrix
      if (fssgd != "No")
      {
        // parameter for artificial diffusivity (scaled to one currently)
        const double kartfac = fac_;
        //const double kartfac = kart_[dofindex]*fac_;

        for (int vi=0; vi<iel_; ++vi)
        {
          for (int ui=0; ui<iel_; ++ui)
          {
            /* artificial diffusivity term */
            rhs[vi*numdof+dofindex] += -(kartfac*(derxy_(0, vi)*derxy_(0, ui) + derxy_(1, vi)*derxy_(1, ui))) ;

            /*subtract SUPG term */
            //rhs[vi*numdof+dofindex] -= -(taufac*conv[vi]*conv[ui]) ;
          }
        }
      }

    } // loop over each scalar

    if (numdofpernode_-numscal_== 1) // ELCH
    {
      // testing: set lower-right block to identity matrix:
      for (int vi=0; vi<iel_; ++vi)
      {
          //fac_funct_vi_densfunct_ui = fac_*funct_[vi]*densfunct_[ui];
          massmat(vi*numdofpernode_+numscal_, vi*numdofpernode_+numscal_) += 1.0;
      }
    }
    
  } // integration loop

  return;
} // Condif3Impl::InitializeOST


/*----------------------------------------------------------------------*
 |  add electrochemistry specific terms (private)              gjb 10/08|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Condif3Impl::CalMatElch(
    Epetra_SerialDenseMatrix& emat,
    Epetra_SerialDenseVector& erhs,
    const vector<double>&     ephinp,
    const bool&               higher_order_ele,
    const double&             frt,
    const bool&               is_stationary,
    const double&             timefac
)
{
  // get values of all transported scalars at integration point
  for (int k=0; k<numscal_; ++k)
  {
    conint_[k] = 0.0;
    for (int ui=0; ui<iel_; ++ui)
    {
      conint_[k] += funct_[ui]*ephinp[ui*numdofpernode_+k];
    }
    // when concentration becomes zero, the coupling terms in the system matrix get lost!
    if (abs(conint_[k])<1e-18) dserror("concentration is nearly singular: %g",conint_[k]);
  }

  // get gradient of el. potential at integration point
  static double pot_ui(0.0);
  gradpot_[0] = 0.0;
  gradpot_[1] = 0.0;
  gradpot_[2] = 0.0;
  for (int ui=0; ui<iel_; ++ui)
  {
    pot_ui = ephinp[ui*numdofpernode_+numscal_];
    gradpot_[0] += derxy_(0,ui)*pot_ui;
    gradpot_[1] += derxy_(1,ui)*pot_ui;
    gradpot_[2] += derxy_(2,ui)*pot_ui;
  }

  for (int i=0; i<iel_; i++)
  {
     /* convective part */
     /* rho * c_p * u_x * N,x  +  rho * c_p * u_y * N,y +  rho * c_p * u_z * N,z
        with  N .. form function matrix */
     conv_[i] = velint_[0] * derxy_(0,i) + velint_[1] * derxy_(1,i) + velint_[2] * derxy_(2,i);

     /* migration part */
     mig_[i] = (-frt)*(gradpot_[0] * derxy_(0,i) + gradpot_[1] * derxy_(1,i) + gradpot_[2] * derxy_(2,i));
  }

#if 0
  // DEBUG output
  cout<<endl<<"values at GP:"<<endl;
  cout<<"factor F/RT = "<<frt<<endl;
  for (int k=0;k<numscal_;++k)
  {cout<<"conint_["<<k<<"] = "<<conint_[k]<<endl;}
  for (int k=0;k<3;++k)
  {cout<<"gradpot_["<<k<<"] = "<<gradpot_[k]<<endl;}
#endif
  // some 'working doubles'
  static double phinp_k_ui(0.0);
  static double diffus_valence_k(0.0);
   double rhsint(0.0);  // rhs at int. point
  // integration factors and coefficients of single terms
  static double timefacfac(0.0);
  static double timetaufac(0.0);


  for (int k = 0; k < numscal_;++k) // loop over all transported sclars
  {
    // stabilization parameters
    const double taufac = tau_[k]*fac_;

    // factor D_k * z_k
    diffus_valence_k = diffus_[k]*valence_[k];

    if (is_stationary)
    {
      timefacfac  = fac_;
      timetaufac  = taufac;
      rhsint = rhs_[k];     // set rhs at integration point
    }
    else
    {
      timefacfac  = timefac * fac_;
      timetaufac  = timefac * taufac;
      rhsint = hist_[k] + rhs_[k]*timefac;     // set rhs at integration point
    }

    if (higher_order_ele)
    {
      for (int i=0; i<iel_; i++)
      {
        /* diffusive part */
        /* diffus * ( N,xx  +  N,yy +  N,zz ) */
        diff_[i] = diffus_[k] * (derxy2_(0,i) + derxy2_(1,i) + derxy2_(2,i));

        /* reactive part of migration*/
        /* diffus * ( N,xx  +  N,yy +  N,zz ) */
        //migr_[i] = diffus_[k] * funct_[i] * (derxy2_(0,i) + derxy2_(1,i) + derxy2_(2,i));
      }
    }

    // ----------------------------------------matrix entries
    for (int vi=0; vi<iel_; ++vi)
    {
      for (int ui=0; ui<iel_; ++ui)
      {
        /* Standard Galerkin terms: */
        /* convective term */
        emat(vi*numdofpernode_+k, ui*numdofpernode_+k) += timefacfac*funct_[vi]*conv_[ui] ;

        /* diffusive term */
        emat(vi*numdofpernode_+k, ui*numdofpernode_+k) += timefacfac*diffus_[k]*(derxy_(0, vi)*derxy_(0, ui) + derxy_(1, vi)*derxy_(1, ui) + derxy_(2, vi)*derxy_(2, ui));

        /* migration term (directional derivatives) */
        emat(vi*numdofpernode_+k, ui*numdofpernode_+k) -= timefacfac*diffus_valence_k*mig_[vi]*funct_[ui];
        emat(vi*numdofpernode_+k,ui*numdofpernode_+numscal_) += frt*timefacfac*diffus_valence_k*conint_[k]*(derxy_(0, vi)*derxy_(0, ui) + derxy_(1, vi)*derxy_(1, ui) + derxy_(2, vi)*derxy_(2, ui));

        /* electroneutrality condition */
        emat(vi*numdofpernode_+numscal_, ui*numdofpernode_+k) += valence_[k]*fac_*funct_[vi]*densfunct_[ui];

        /* Stabilization term: */
        /* 0) transient stabilization */
              // not implemented

        /* 1) convective stabilization */

        /* convective term */
        emat(vi*numdofpernode_+k, ui*numdofpernode_+k) += timetaufac*(conv_[vi]+diffus_valence_k*mig_[vi])*(conv_[ui]+diffus_valence_k*mig_[ui]);

        if (higher_order_ele)
        {
          /* diffusive term */
          emat(vi*numdofpernode_+k, ui*numdofpernode_+k) += -timetaufac*conv_[vi]*diff_[ui] ;

          /* 2) diffusive stabilization (USFEM assumed here, sign change necessary for GLS) */

          /* convective term */
          emat(vi*numdofpernode_+k, ui*numdofpernode_+k) += timetaufac*diff_[vi]*conv_[ui] ;

          /* diffusive term */
          emat(vi*numdofpernode_+k, ui*numdofpernode_+k) -= timetaufac*diff_[vi]*diff_[ui] ;
        }

      } // for ui
    } // for vi

    // ----------------------------------------------RHS
    for (int vi=0; vi<iel_; ++vi)
    {
      /* RHS source term */
      erhs[vi*numdofpernode_+k] += fac_*funct_[vi]*rhsint ;

      /* transient stabilization of RHS source term */
      // not implemented

      /* convective stabilization of RHS source term */
      erhs[vi*numdofpernode_+k] += taufac*conv_[vi]*rhsint ;

      /* diffusive stabilization of RHS source term */
      if (higher_order_ele) erhs[vi*numdofpernode_+k] += taufac*diff_[vi]*rhsint ;

      //-----------------------------------residual formulation
      // nonlinear migration term
      erhs[vi*numdofpernode_+k] += conint_[k]*timefacfac*diffus_valence_k*mig_[vi];

      for (int ui=0; ui<iel_; ++ui)
      {
        phinp_k_ui = ephinp[ui*numdofpernode_+k];
        // convective term
        erhs[vi*numdofpernode_+k] -= timefacfac*funct_[vi]*conv_[ui]*phinp_k_ui;
        // diffusive term
        erhs[vi*numdofpernode_+k] -= timefacfac*diffus_[k]*phinp_k_ui*(derxy_(0, ui)*derxy_(0, vi) + derxy_(1, ui)*derxy_(1, vi)+ derxy_(2, ui)*derxy_(2, vi));

        /* Stabilization term: */
        /* 1) convective stabilization */

        /* convective term */
        erhs[vi*numdofpernode_+k] -= timetaufac*(conv_[vi]+diffus_valence_k*mig_[vi])*(conv_[ui]+diffus_valence_k*mig_[ui])*phinp_k_ui ;

        if (higher_order_ele)
        {
          dserror("higher order terms not yet finished");
          /* diffusive term */
          erhs[vi*numdofpernode_+k] -= -timetaufac*conv_[vi]*diff_[ui]*phinp_k_ui ;

          /* 2) diffusive stabilization (USFEM assumed here, sign change necessary for GLS) */

          /* convective term */
          erhs[vi*numdofpernode_+k] -= timetaufac*diff_[vi]*conv_[ui]*phinp_k_ui ;

          /* diffusive term */
          erhs[vi*numdofpernode_+k] += timetaufac*diff_[vi]*diff_[ui]*phinp_k_ui ;
        }

        // electroneutrality condition
        // for incremental formulation, there is the residuum on the rhs! : 0-ENC*phi_i
        erhs[vi*numdofpernode_+numscal_] -= valence_[k]*fac_*funct_[vi]*densfunct_[ui]*phinp_k_ui;

      } // for ui
    } // for vi

    // -----------------------------------INSTATIONARY TERMS
    if (!is_stationary)
    {
      for (int vi=0; vi<iel_; ++vi)
      {
        for (int ui=0; ui<iel_; ++ui)
        {
          /* Standard Galerkin terms: */
          /* transient term */
          emat(vi*numdofpernode_+k, ui*numdofpernode_+k) += fac_*funct_[vi]*densfunct_[ui] ;

          /* 1) convective stabilization */
          /* transient term */
          emat(vi*numdofpernode_+k, ui*numdofpernode_+k) += taufac*(conv_[vi]+diffus_valence_k*mig_[vi])*densfunct_[ui];

          if (higher_order_ele)
          {
            /* 2) diffusive stabilization (USFEM assumed here, sign change necessary for GLS) */
            /* transient term */
            emat(vi*numdofpernode_+k, ui*numdofpernode_+k) += taufac*diff_[vi]*densfunct_[ui] ;
          }
        } // for ui

        // residuum on RHS:

        /* Standard Galerkin terms: */
        /* transient term */
        erhs[vi*numdofpernode_+k] -= fac_*funct_[vi]*conint_[k];

        /* 1) convective stabilization */
        /* transient term */
        erhs[vi*numdofpernode_+k] -= taufac*(conv_[vi]+diffus_valence_k*mig_[vi])*conint_[k];

        if (higher_order_ele)
        {
          /* 2) diffusive stabilization (USFEM assumed here, sign change necessary for GLS) */
          /* transient term */
          erhs[vi*numdofpernode_+k] -= taufac*diff_[vi]*conint_[k];
        }
      } // for vi
    } // instationary case

  } // loop over scalars

  return;
} // Condif3Impl::AddElchTerms

#endif
#endif
