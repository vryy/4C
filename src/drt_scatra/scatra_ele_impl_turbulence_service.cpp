/*----------------------------------------------------------------------*/
/*!
\file scatra_ele_impl_turbulence_service.cpp

\brief Internal implementation of ScaTra element

<pre>
Maintainer: Ursula Rasthofer
            rasthofer@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15236
</pre>
*/
/*----------------------------------------------------------------------*/

#include "scatra_ele_impl.H"

#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_lib/standardtypes_cpp.H"

#include "../drt_mat/newtonianfluid.H"
#include "../drt_mat/sutherland.H"

// include define flags for turbulence models under development
#include "../drt_fluid/fluid_turbulence_defines.H"


/*-----------------------------------------------------------------------------*
 | calculate filtered quantities for dynamic Smagorinsky model  rasthofer 08/12|
 *-----------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraImpl<distype>::scatra_apply_box_filter(
  const double               thermpress,
  double&                    dens_hat,
  double&                    temp_hat,
  double&                    dens_temp_hat,
  RCP<vector<double> >       vel_hat,
  RCP<vector<double> >       densvel_hat,
  RCP<vector<double> >       densveltemp_hat,
  RCP<vector<double> >       densstraintemp_hat,
  double&                    volume,
  const DRT::Element*        ele)
{
  // do preparations first
  // ---------------------------------------------

  // use one-point Gauss rule to do calculations at the element center
  DRT::UTILS::IntPointsAndWeights<nsd_> intpoints(SCATRA::DisTypeToStabGaussRule<distype>::rule);

  volume = EvalShapeFuncAndDerivsAtIntPoint(intpoints,0,0);

  // get phi at integration point
  phi_[0] = funct_.Dot(ephinp_[0]);

  // get material
  RefCountPtr<MAT::Material> material = ele->Material();
  if (material->MaterialType() == INPAR::MAT::m_scatra)
  {
    //access fluid discretization
    RCP<DRT::Discretization> fluiddis = null;
    fluiddis = DRT::Problem::Instance()->GetDis("fluid");
    //get corresponding fluid element (it has the same global ID as the scatra element)
    DRT::Element* fluidele = fluiddis->gElement(ele->Id());
    if (fluidele == NULL)
      dserror("Fluid element %i not on local processor", ele->Id());

    // get fluid material
    RCP<MAT::Material> fluidmat = fluidele->Material();
    if(fluidmat->MaterialType() != INPAR::MAT::m_fluid)
      dserror("Invalid fluid material for passive scalar transport in turbulent flow!");

    const Teuchos::RCP<const MAT::NewtonianFluid>& actfluidmat = Teuchos::rcp_dynamic_cast<const MAT::NewtonianFluid>(fluidmat);

    densnp_[0] = actfluidmat->Density();
     if (densnp_[0] != 1.0)
       dserror("Check your diffusivity! Dynamic diffusivity required!");
  }
  else if (material->MaterialType() == INPAR::MAT::m_sutherland)
  {
    const Teuchos::RCP<const MAT::Sutherland>& actmat
      = Teuchos::rcp_dynamic_cast<const MAT::Sutherland>(material);

    densnp_[0] = actmat->ComputeDensity(phi_[0],thermpress);
  }
  else dserror("material not supported");


  // get velocities (n+alpha_F/1,i) at integration point
  convelint_.Multiply(evelnp_,funct_);

  // compute rate of strain
  double rateofstrain = -1.0e30;
  rateofstrain = GetStrainRate(evelnp_);

  // gradient of scalar value
  gradphi_.Multiply(derxy_,ephinp_[0]);

  // perform integrations, i.e., convolution
  // ---------------------------------------------

  for (int rr=0;rr<nsd_;++rr)
  {
    double tmp=convelint_(rr)*volume;

    // add contribution to integral over velocities
    (*vel_hat)[rr] += tmp;

    // add contribution to integral over dens times velocity
    (*densvel_hat)[rr] += densnp_[0]*tmp;

    // add contribution to integral over dens times temperature times velocity
    (*densveltemp_hat)[rr] += densnp_[0]*phi_[0]*tmp;
  }

  for (int rr=0;rr<nsd_;++rr)
  {
    double tmp=gradphi_(rr)*volume;
    // add contribution to integral over dens times rate of strain times phi gradient
    (*densstraintemp_hat)[rr] += densnp_[0]*rateofstrain*tmp;
  }

  // add additional scalar quantities
  // i.e., filtered density, filtered density times scalar (i.e., temperature) and scalar
  dens_hat = densnp_[0]*volume;
  dens_temp_hat = densnp_[0]*phi_[0]*volume;
  temp_hat = phi_[0]*volume;

  return;
} //ScaTraImpl::scatra_apply_box_filter


/*----------------------------------------------------------------------------------*
 | calculate turbulent Prandtl number for dynamic Smagorinsky model  rasthofer 08/12|
 *----------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraImpl<distype>::scatra_calc_smag_const_LkMk_and_MkMk(
        RCP<Epetra_MultiVector>&  col_filtered_vel,
        RCP<Epetra_MultiVector>&  col_filtered_dens_vel,
        RCP<Epetra_MultiVector>&  col_filtered_dens_vel_temp,
        RCP<Epetra_MultiVector>&  col_filtered_dens_rateofstrain_temp,
        RCP<Epetra_Vector>&       col_filtered_temp,
        RCP<Epetra_Vector>&       col_filtered_dens,
        RCP<Epetra_Vector>&       col_filtered_dens_temp,
        double&                   LkMk,
        double&                   MkMk,
        double&                   xcenter,
        double&                   ycenter,
        double&                   zcenter,
        const DRT::Element*       ele
)
{
  LINALG::Matrix<nsd_,nen_> evel_hat;
  LINALG::Matrix<nsd_,nen_> edensvel_hat;
  LINALG::Matrix<nsd_,nen_> edensveltemp_hat;
  LINALG::Matrix<nsd_,nen_> edensstraintemp_hat;
  LINALG::Matrix<1,nen_> etemp_hat;
  LINALG::Matrix<1,nen_> edens_hat;
  LINALG::Matrix<1,nen_> edenstemp_hat;
  // extract required (node-based) filtered quantities
  DRT::UTILS::ExtractMyNodeBasedValues(ele,evel_hat,col_filtered_vel,nsd_);
  DRT::UTILS::ExtractMyNodeBasedValues(ele,edensvel_hat,col_filtered_dens_vel,nsd_);
  DRT::UTILS::ExtractMyNodeBasedValues(ele,edensveltemp_hat,col_filtered_dens_vel_temp,nsd_);
  DRT::UTILS::ExtractMyNodeBasedValues(ele,edensstraintemp_hat,col_filtered_dens_rateofstrain_temp,nsd_);
  DRT::UTILS::ExtractMyNodeBasedValues(ele,etemp_hat,col_filtered_temp,1);
  DRT::UTILS::ExtractMyNodeBasedValues(ele,edens_hat,col_filtered_dens,1);
  DRT::UTILS::ExtractMyNodeBasedValues(ele,edenstemp_hat,col_filtered_dens_temp,1);

  // get center coordinates of element
  xcenter = 0.0;
  ycenter = 0.0;
  zcenter = 0.0;
  for(int inode=0;inode<nen_;inode++)
  {
    // xyze_ has been initialized at the beginning of Impl()
    xcenter+=xyze_(0,inode);
    ycenter+=xyze_(1,inode);
    zcenter+=xyze_(2,inode);
  }
  xcenter/=nen_;
  ycenter/=nen_;
  zcenter/=nen_;

  // use one-point Gauss rule to do calculations at the element center
  DRT::UTILS::IntPointsAndWeights<nsd_> intpoints(SCATRA::DisTypeToStabGaussRule<distype>::rule);

  EvalShapeFuncAndDerivsAtIntPoint(intpoints,0,0);

  // get filtered dens * velocities (n+alpha_F/1,i) at integration point
  //
  //                     +-----
  //        ^ n+af/1     \                  ^ n+af/1
  //    rho * vel (x) =   +      N (x) * rho*vel
  //                     /        j             j
  //                     +-----
  //                    node j
  //
  LINALG::Matrix<nsd_,1> densvelint_hat;
  densvelint_hat.Multiply(edensvel_hat,funct_);

  // get filtered dens * temp * velocities (n+alpha_F/1,i) at integration point
  //
  //                           +-----
  //        ^ n+af/1            \                    ^ n+af/1
  //    rho * vel * temp (x) =   +      N (x) * rho*temp*vel
  //                            /        j                  j
  //                           +-----
  //                           node j
  //
  LINALG::Matrix<nsd_,1> densveltempint_hat;
  densveltempint_hat.Multiply(edensveltemp_hat,funct_);

  // get filtered dens * rate-of-strain * grad T (n+alpha_F/1,i) at integration point
  //
  //                                        +-----
  //                ^ n+af/1                 \                    ^ n+af/1
  //    rho * rate-of-strain * grad T (x) =   +      N (x) * rho * rate-of-strain * grad T (x)
  //                                         /        j                                       j
  //                                        +-----
  //                                         node j
  //
  LINALG::Matrix<nsd_,1> densstraintempint_hat;
  densstraintempint_hat.Multiply(edensstraintemp_hat,funct_);

  // get filtered density at integration point
  //
  //         +-----
  //    ^     \              ^
  //   rho =   +    N (x) * rho
  //          /      k         k
  //         +-----
  //         node k
  //
  //
  double densint_hat = 0.0;
  // get filtered density times temperature at integration point
  //
  //            +-----
  //     ^       \                 ^
  //  rho * T =   +    N (x) * rho * T
  //             /      k             k
  //            +-----
  //              node k
  //
  //
  double denstempint_hat = 0.0;
//  double tempint_hat = 0.0;
  for (int mm=0;mm<nen_;++mm)
  {
    densint_hat += funct_(mm,0)*edens_hat(0,mm);
    denstempint_hat += funct_(mm,0)*edenstemp_hat(0,mm);
  }

  // compute rate of strain
  double rateofstrain = -1.0e30;
  rateofstrain = GetStrainRate(evel_hat);

  // gradient of scalar value (i.e., of filtered temperature)
  LINALG::Matrix<nsd_,1> gradtemp_hat;
  gradtemp_hat.MultiplyNT(derxy_,etemp_hat);

  // this is sqrt(3)
  const double filterwidthratio = 1.73;

  // calculate L_k and M_k
  LINALG::Matrix<nsd_,1> L_k;
  LINALG::Matrix<nsd_,1> M_k;

  for(int rr=0;rr<nsd_;rr++)
  {
    L_k(rr,0) = densveltempint_hat(rr,0) - densvelint_hat(rr,0) * denstempint_hat / densint_hat;
    M_k(rr,0) = densstraintempint_hat(rr,0) - filterwidthratio * filterwidthratio * densint_hat * rateofstrain * gradtemp_hat(rr,0);
  }

  // perform contraction via dot product
  LkMk =0.0;
  MkMk =0.0;
  for(int rr=0;rr<nsd_;rr++)
  {
    LkMk += L_k(rr,0) * M_k(rr,0);
    MkMk += M_k(rr,0) * M_k(rr,0);
  }

  return;
}//ScaTraImpl::scatra_calc_smag_const_LkMk_and_MkMk


/*----------------------------------------------------------------------------------*
 | calculate mean turbulent Prandtl number                           rasthofer 08/12|
 *----------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraImpl<distype>::GetMeanPrtOfHomogenousDirection(
  ParameterList&             turbmodelparams,
  double&                    inv_Prt,
  int&                       nlayer
)
{
  if(nsd_ != 3)
    dserror("turbulence and 3D flow !");

  if (turbmodelparams.get<string>("HOMDIR","not_specified") !=  "not_specified")
  {
    RCP<vector<double> > averaged_LkMk = turbmodelparams.get<RCP<vector<double> > >("averaged_LkMk_");
    RCP<vector<double> > averaged_MkMk = turbmodelparams.get<RCP<vector<double> > >("averaged_MkMk_");

    // get homogeneous direction
    string homdir = turbmodelparams.get<string>("HOMDIR","not_specified");

    // here, the layer is determined via the element center in order to get the correct
    // averaged value from the vector of averaged (M/L)kMk
    double xcenter = 0.0;
    double ycenter = 0.0;
    double zcenter = 0.0;
    for(int inode=0;inode<nen_;inode++)
    {
      xcenter += xyze_( 0, inode );
      ycenter += xyze_( 1, inode );
      zcenter += xyze_( 2, inode );
    }
    xcenter/=nen_;
    ycenter/=nen_;
    zcenter/=nen_;

    // determine the layer
    if (homdir == "xy" or homdir == "xz" or homdir == "yz")
    {
      RCP<vector<double> > planecoords = turbmodelparams.get<RCP<vector<double> > >("planecoords_");
      // get center
      double center = 0.0;
      if (homdir == "xy")
        center = zcenter;
      else if (homdir == "xz")
        center = ycenter;
      else if (homdir == "yz")
        center = xcenter;

      bool found = false;
      for (nlayer=0;nlayer < static_cast<int>((*planecoords).size()-1);)
      {
        if(center<(*planecoords)[nlayer+1])
        {
          found = true;
          break;
        }
        nlayer++;
      }
      if (found ==false)
      {
        dserror("could not determine element layer");
      }
    }
    else if (homdir == "x" or homdir == "y" or homdir == "z")
    {
      RCP<vector<double> > dir1coords = turbmodelparams.get<RCP<vector<double> > >("dir1coords_");
      RCP<vector<double> > dir2coords = turbmodelparams.get<RCP<vector<double> > >("dir2coords_");
      // get center
      double dim1_center = 0.0;
      double dim2_center = 0.0;
      if (homdir == "x")
      {
        dim1_center = ycenter;
        dim2_center = zcenter;
      }
      else if (homdir == "y")
      {
        dim1_center = xcenter;
        dim2_center = zcenter;
      }
      else if (homdir == "z")
      {
        dim1_center = xcenter;
        dim2_center = ycenter;
      }

      int  n1layer = 0;
      int  n2layer = 0;
      bool dir1found = false;
      bool dir2found = false;
      for (n1layer=0;n1layer<(int)(*dir1coords).size()-1;)
      {
        if(dim1_center<(*dir1coords)[n1layer+1])
        {
          dir1found = true;
          break;
        }
        n1layer++;
      }
      if (dir1found ==false)
      {
        dserror("could not determine element layer");
      }
      for (n2layer=0;n2layer<(int)(*dir2coords).size()-1;)
      {
        if(dim2_center<(*dir2coords)[n2layer+1])
        {
          dir2found = true;
          break;
        }
        n2layer++;
      }
      if (dir2found ==false)
      {
        dserror("could not determine element layer");
      }

      const int numdir1layer = (int)(*dir1coords).size()-1;
      nlayer = numdir1layer * n2layer + n1layer;
    }
    else
      dserror("More than two homogeneous directions not supported!");

    // (Cs*Delta)^2/Prt is set by the averaged quantities
    if ((*averaged_MkMk)[nlayer] < 1E-16)
    {
    //  std::cout << "warning: abs(averaged_MkMk) < 1E-16 -> set inverse of turbulent Prandtl number to zero!"  << std::endl;
      inv_Prt = 0.0;
    }
    else
      inv_Prt = (*averaged_LkMk)[nlayer]/(*averaged_MkMk)[nlayer] ;
    // clipping to get algorithm stable
    if (inv_Prt<0.0)
    {
      inv_Prt=0.0;
    }

    }

  return;
}//ScaTraImpl::GetMeanOfHomogenousDirection


/*----------------------------------------------------------------------*
  |  calculate all-scale art. subgrid diffusivity (private)     vg 10/09 |
  *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraImpl<distype>::CalcSubgrDiff(
  const double                          dt,
  const double                          timefac,
  const enum INPAR::SCATRA::AssgdType   whichassgd,
  const bool                            assgd,
  const double                          Cs,
  const double                          tpn,
  const double                          vol,
  const int                             k
  )
{
  // get number of dimensions
  const double dim = (double) nsd_;

  // get characteristic element length as cubic root of element volume
  // (2D: square root of element area, 1D: element length)
  const double h = std::pow(vol,(1.0/dim));

  // artficial all-scale subgrid diffusivity
  if (assgd )
  {
    // classical linear artificial all-scale subgrid diffusivity
    if (whichassgd == INPAR::SCATRA::assgd_artificial)
    {
      // get element-type constant
      const double mk = SCATRA::MK<distype>();

      // velocity norm
      const double vel_norm = convelint_.Norm2();

      // parameter relating convective and diffusive forces + respective switch
      const double epe = mk * densnp_[k] * vel_norm * h / diffus_[k];
      const double xi = std::max(epe,1.0);

      // compute subgrid diffusivity
      sgdiff_[k] = (DSQR(h)*mk*DSQR(vel_norm)*DSQR(densnp_[k]))/(2.0*diffus_[k]*xi);
    }
    else
    {
      // gradient of current scalar value
      gradphi_.Multiply(derxy_,ephinp_[k]);

      // gradient norm
      const double grad_norm = gradphi_.Norm2();

      if (grad_norm > EPS10)
      {
        // compute residual of scalar transport equation
        // (subgrid-scale part of scalar, which is also computed, not required)
        CalcResidualAndSubgrScalar(dt,timefac,k);

        // for the present definitions, sigma and a specific term (either
        // residual or convective term) are different
        double sigma = 0.0;
        double specific_term = 0.0;
        switch (whichassgd)
        {
        case INPAR::SCATRA::assgd_hughes:
        {
          // get norm of velocity vector b_h^par
          const double vel_norm_bhpar = abs(conv_phi_[k]/grad_norm);

          // compute stabilization parameter based on b_h^par
          // (so far, only exact formula for stationary 1-D implemented)
          // element Peclet number relating convective and diffusive forces
          double epe = 0.5 * vel_norm_bhpar * h / diffus_[k];
          const double pp = exp(epe);
          const double pm = exp(-epe);
          double xi = 0.0;
          double tau_bhpar = 0.0;
          if (epe >= 700.0) tau_bhpar = 0.5*h/vel_norm_bhpar;
          else if (epe < 700.0 and epe > EPS15)
          {
            xi = (((pp+pm)/(pp-pm))-(1.0/epe)); // xi = coth(epe) - 1/epe
            // compute optimal stabilization parameter
            tau_bhpar = 0.5*h*xi/vel_norm_bhpar;
          }

          // compute sigma
          sigma = std::max(0.0,tau_bhpar-tau_[k]);

          // set specific term to convective term
          specific_term = conv_phi_[k];
        }
        break;
        case INPAR::SCATRA::assgd_tezduyar:
        {
          // velocity norm
          const double vel_norm = convelint_.Norm2();

          // get norm of velocity vector b_h^par
          const double vel_norm_bhpar = abs(conv_phi_[k]/grad_norm);

          // compute stabilization parameter based on b_h^par
          // (so far, only exact formula for stationary 1-D implemented)

          // compute sigma (version 1 according to John and Knobloch (2007))
          //sigma = (h/vel_norm)*(1.0-(vel_norm_bhpar/vel_norm));

          // compute sigma (version 2 according to John and Knobloch (2007))
          // setting scaling phi_0=1.0 as in John and Knobloch (2007)
          const double phi0 = 1.0;
          sigma = (h*h*grad_norm/(vel_norm*phi0))*(1.0-(vel_norm_bhpar/vel_norm));

          // set specific term to convective term
          specific_term = conv_phi_[k];
        }
        break;
        case INPAR::SCATRA::assgd_docarmo:
        case INPAR::SCATRA::assgd_almeida:
        {
          // velocity norm
          const double vel_norm = convelint_.Norm2();

          // get norm of velocity vector z_h
          const double vel_norm_zh = abs(scatrares_[k]/grad_norm);

          // parameter zeta differentiating approaches by doCarmo and Galeao (1991)
          // and Almeida and Silva (1997)
          double zeta = 0.0;
          if (whichassgd == INPAR::SCATRA::assgd_docarmo)
            zeta = 1.0;
          else zeta = std::max(1.0,(conv_phi_[k]/scatrares_[k]));

          // compute sigma
          sigma = tau_[k]*std::max(0.0,(vel_norm/vel_norm_zh)-zeta);

          // set specific term to residual
          specific_term = scatrares_[k];
        }
        break;
        default: dserror("unknown type of all-scale subgrid diffusivity\n");
        } //switch (whichassgd)

        // computation of subgrid diffusivity
        sgdiff_[k] = sigma*scatrares_[k]*specific_term/(grad_norm*grad_norm);
      }
      else sgdiff_[k] = 0.0;
    }
  }
  // all-scale subgrid diffusivity due to Smagorinsky model divided by
  // turbulent Prandtl number
  else if (turbmodel_ == INPAR::FLUID::smagorinsky)
  {
    //
    // SMAGORINSKY MODEL
    // -----------------
    //                                   +-                                 -+ 1
    //                               2   |          / h \           / h \    | -
    //    visc          = dens * lmix  * | 2 * eps | u   |   * eps | u   |   | 2
    //        turbulent           |      |          \   / ij        \   / ij |
    //                            |      +-                                 -+
    //                            |
    //                            |      |                                   |
    //                            |      +-----------------------------------+
    //                            |           'resolved' rate of strain
    //                    mixing length
    // -> either provided by dynamic modeling procedure and stored in Cs_delta_sq
    // -> or computed based on fixed Smagorinsky constant Cs:
    //             Cs = 0.17   (Lilly --- Determined from filter
    //                          analysis of Kolmogorov spectrum of
    //                          isotropic turbulence)
    //             0.1 < Cs < 0.24 (depending on the flow)
    //
    //                       dens * visc
    //                                  turbulent
    //    kappa           =  ---------------------
    //         turbulent         Pr
    //                             turbulent
    // -> Prt prescribed in input file or estimated dynamically
    //

    // compute (all-scale) rate of strain
    double rateofstrain = -1.0e30;
    rateofstrain = GetStrainRate(econvelnp_);

    // subgrid diffusivity = subgrid viscosity / turbulent Prandtl number
    sgdiff_[k] = densnp_[k] * Cs * Cs * h * h * rateofstrain / tpn;

    // add subgrid viscosity to physical viscosity for computation
    // of subgrid-scale velocity when turbulence model is applied
    if (sgvel_) visc_ += sgdiff_[k]*tpn;
  }
  else if (turbmodel_ == INPAR::FLUID::dynamic_smagorinsky)
  {
    // compute (all-scale) rate of strain
    double rateofstrain = -1.0e30;
    rateofstrain = GetStrainRate(econvelnp_);

    // subgrid diffusivity = subgrid viscosity / turbulent Prandtl number
    // remark: for dynamic estimation, tpn corresponds to (Cs*h)^2 / Pr_t
    sgdiff_[k] = densnp_[k] * rateofstrain * tpn;
  }

  // compute sum of physical and all-scale subgrid diffusivity
  // -> set internal variable for use when calculating matrix and rhs
  diffus_[k] += sgdiff_[k];
  if (update_mat_)
    dserror("Material update will overwrite effective diffusivity due to eddy-diffusivity model");

  return;
} //ScaTraImpl::CalcSubgrDiff


/*----------------------------------------------------------------------*
  |  calculate fine-scale art. subgrid diffusivity (private)    vg 10/09 |
  *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraImpl<distype>::CalcFineScaleSubgrDiff(
  DRT::Element*                         ele,
  Epetra_SerialDenseVector&             subgrdiff,
  const enum INPAR::SCATRA::FSSUGRDIFF  whichfssgd,
  const double                          Cs,
  const double                          tpn,
  const double                          vol,
  const int                             k
  )
{
  // get number of dimensions
  const double dim = (double) nsd_;

  // get characteristic element length as cubic root of element volume
  // (2D: square root of element area, 1D: element length)
  const double h = std::pow(vol,(1.0/dim));

  //----------------------------------------------------------------------
  // computation of fine-scale subgrid diffusivity for non-incremental
  // solver -> only artificial subgrid diffusivity
  // (values are stored in subgrid-diffusivity-scaling vector)
  //----------------------------------------------------------------------
  if (not is_incremental_)
  {
    // get element-type constant
    const double mk = SCATRA::MK<distype>();

    // velocity norm
    const double vel_norm = convelint_.Norm2();

    // parameter relating convective and diffusive forces + respective switch
    const double epe = mk * densnp_[k] * vel_norm * h / diffus_[k];
    const double xi = std::max(epe,1.0);

    // compute artificial subgrid diffusivity
    sgdiff_[k] = (DSQR(h)*mk*DSQR(vel_norm)*DSQR(densnp_[k]))/(2.0*diffus_[k]*xi);

    // compute entries of (fine-scale) subgrid-diffusivity-scaling vector
    for (int vi=0; vi<nen_; ++vi)
    {
      subgrdiff(vi) = sgdiff_[k]/ele->Nodes()[vi]->NumElement();
    }
  }
  //----------------------------------------------------------------------
  // computation of fine-scale subgrid diffusivity for incremental solver
  // -> only all-scale Smagorinsky model
  //----------------------------------------------------------------------
  else
  {
    if (whichfssgd == INPAR::SCATRA::fssugrdiff_smagorinsky_all)
    {
      //
      // ALL-SCALE SMAGORINSKY MODEL
      // ---------------------------
      //                                      +-                                 -+ 1
      //                                  2   |          / h \           / h \    | -
      //    visc          = dens * (C_S*h)  * | 2 * eps | u   |   * eps | u   |   | 2
      //        turbulent                     |          \   / ij        \   / ij |
      //                                      +-                                 -+
      //                                      |                                   |
      //                                      +-----------------------------------+
      //                                            'resolved' rate of strain
      //

      // compute (all-scale) rate of strain
      double rateofstrain = -1.0e30;
      rateofstrain = GetStrainRate(econvelnp_);

      // subgrid diffusivity = subgrid viscosity / turbulent Prandtl number
      sgdiff_[k] = densnp_[k] * Cs * Cs * h * h * rateofstrain / tpn;
    }
    else if (whichfssgd == INPAR::SCATRA::fssugrdiff_smagorinsky_small)
    {
      //
      // FINE-SCALE SMAGORINSKY MODEL
      // ----------------------------
      //                                      +-                                 -+ 1
      //                                  2   |          /    \          /   \    | -
      //    visc          = dens * (C_S*h)  * | 2 * eps | fsu |   * eps | fsu |   | 2
      //        turbulent                     |          \   / ij        \   / ij |
      //                                      +-                                 -+
      //                                      |                                   |
      //                                      +-----------------------------------+
      //                                           'fine-scale' rate of strain
      //

      // fine-scale rate of strain
      double fsrateofstrain = -1.0e30;
      fsrateofstrain = GetStrainRate(efsvel_);

      // subgrid diffusivity = subgrid viscosity / turbulent Prandtl number
      sgdiff_[k] = densnp_[k] * Cs * Cs * h * h * fsrateofstrain / tpn;
      }
  }

  return;
} //ScaTraImpl::CalcFineScaleSubgrDiff


/*----------------------------------------------------------------------*
 | calculation of coefficients B and D for multifractal subgrid-scales  |
 |                                                      rasthofer 12/11 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraImpl<distype>::CalcBAndDForMultifracSubgridScales(
    LINALG::Matrix<nsd_,1>&                     B_mfs, ///< coefficient for fine-scale velocity (will be filled)
    double &                                    D_mfs, ///< coefficient for fine-scale scalar (will be filled)
    const double                                Csgs_sgvel, ///< parameter of multifractal subgrid-scales (velocity)
    const double                                alpha, ///< grid-filter to test-filter ratio
    const bool                                  calc_N, ///< flag to activate calculation of N
    const double                                N_vel, ///< value for N if not calculated
    const enum INPAR::FLUID::RefVelocity        refvel, ///< reference velocity
    const enum INPAR::FLUID::RefLength          reflength, ///< reference length
    const double                                c_nu, ///< scaling for Re
    const bool                                  nwl, ///< flag to activate near-wall limit
    const double                                Csgs_sgphi, ///< parameter of multifractal subgrid-scales (phi)
    const double                                c_diff, ///< scaling for Re*Pr
    const double                                vol, ///< volume of element
    const int                                   k
)
{
  //----------------------------------------------------------------
  // calculation of B for fine-scale velocity
  //----------------------------------------------------------------

  // STEP1: determine N and Csgs

  // allocate vector for parameter N
  // N may depend on the direction -> currently unused
  vector<double> Nvel (3);
  // variable for final (corrected) Csgs_vel
  double Csgs_vel_nw = Csgs_sgvel;

  // potential calculation of Re to determine N
  double Re_ele = -1.0;
  // characteristic element length
  double hk = 1.0e+10;
  double strainnorm = 0.0;
  // ratio of viscous scale to element length
  double scale_ratio = 0.0;

  // get velocity at element center
  // convelint_.Multiply(econvelnp_,funct_);
  // get norm of velocity
  const double vel_norm = convelint_.Norm2();
  // also for fine-scale velocity
  // fsvelint_.Multiply(efsvel_,funct_);
  const double fsvel_norm = fsvelint_.Norm2();

  // do we have a fixed parameter N
  if (not calc_N)
  {
    // yes, store value
    for (int rr=1;rr<3;rr++)
      Nvel[rr] = N_vel;
  }
  else //no, so we calculate N from Re
  {
    // calculate characteristic element length
    double hk = CalcRefLength(reflength,vol);

    // warning: k=0, this first scalar is taken!
    // multifractal subgrid-scale model is for passive and active
    // scalar transport
    // therefore, we need the density of the fluid here
    switch (refvel)
    {
      case INPAR::FLUID::resolved:
      {
        Re_ele = vel_norm * hk * densnp_[0] / visc_;
        break;
      }
      case INPAR::FLUID::fine_scale:
      {
        Re_ele = fsvel_norm * hk * densnp_[0] / visc_;
        break;
      }
      case INPAR::FLUID::strainrate:
      {
        strainnorm = GetStrainRate(econvelnp_);
        strainnorm /= sqrt(2.0);
        Re_ele = strainnorm * hk * hk * densnp_[0] / visc_;
        break;
      }
      default:
        dserror("Unknown velocity!");
    }
    if (Re_ele < 0.0)
      dserror("Something went wrong!");
    // clip Re to prevent negative N
    if (Re_ele < 1.0)
       Re_ele = 1.0;

    //
    //   Delta
    //  ---------  ~ Re^(3/4)
    //  lambda_nu
    //
    scale_ratio = c_nu * pow(Re_ele,3.0/4.0);
    // scale_ratio < 1.0 leads to N < 0
    // therefore, we clip once more
    if (scale_ratio < 1.0)
      scale_ratio = 1.0;

    //         |   Delta     |
    //  N =log | ----------- |
    //        2|  lambda_nu  |
   double N_re = log(scale_ratio)/log(2.0);
   if (N_re < 0.0)
      dserror("Something went wrong when calculating N!");

    // store calculated N
    for (int i=0; i<nsd_; i++)
      Nvel[i] = N_re;
  }

  // calculate near-wall correction
  if (nwl)
  {
    // not yet calculated, estimate norm of strain rate
    if (not calc_N or refvel != INPAR::FLUID::strainrate)
    {
      strainnorm = GetStrainRate(econvelnp_);
      strainnorm /= sqrt(2.0); //cf. Burton & Dahm 2008
    }
    // and reference length
    if (not calc_N)
      hk = CalcRefLength(reflength,vol);

    // get Re from strain rate
    double Re_ele_str = strainnorm * hk * hk * densnp_[0] / visc_;
    if (Re_ele_str < 0.0)
      dserror("Something went wrong!");
    // ensure positive values
    if (Re_ele_str < 1.0)
       Re_ele_str = 1.0;

    // calculate corrected Csgs
    //           -3/16
    //  *(1 - (Re)   )
    //
    Csgs_vel_nw *= (1-pow(Re_ele_str,-3.0/16.0));
  }

  // STEP 2: calculate B

  //                                  1
  //          |       1              |2
  //  kappa = | -------------------- |
  //          |  1 - alpha ^ (-4/3)  |
  //
  double kappa = 1.0/(1.0-pow(alpha,-4.0/3.0));

  //                                                       1
  //                                  |                   |2
  //  B = Csgs * kappa * 2 ^ (-2*N/3) * | 2 ^ (4*N/3) - 1 |
  //                                  |                   |
  //
  for (int dim=0; dim<nsd_; dim++)
  {
    B_mfs(dim,0) = Csgs_vel_nw *sqrt(kappa) * pow(2.0,-2.0*Nvel[dim]/3.0) * sqrt((pow(2.0,4.0*Nvel[dim]/3.0)-1));
//    if (eid_ == 100)
//     std::cout << "B  " << std::setprecision (10) << B_mfs(dim,0) << std::endl;
  }

  //----------------------------------------------------------------
  // calculation of D for fine-scale scalar
  //----------------------------------------------------------------

  // STEP 1: determine N
  //         currently constant C_sgs for D assumed

  // calculate prandtl number or schmidt number (passive scalar)
  const double Pr = visc_/diffus_[k];

  // allocate vector for parameter N
  double Nphi = 0.0;
  // ratio of dissipation scale to element length
  double scale_ratio_phi = 0.0;

  if (calc_N)
  {
    //
    //   Delta
    //  ---------  ~ Re^(3/4)*Pr^(p)
    //  lambda_diff
    //
    // Pr <= 1: p=3/4
    // Pr >> 1: p=1/2
    double p = 3.0/4.0;
    if (Pr>1.0) p =1.0/2.0;

    scale_ratio_phi = c_diff * pow(Re_ele,3.0/4.0) * pow(Pr,p);
    // scale_ratio < 1.0 leads to N < 0
    // therefore, we clip again
    if (scale_ratio_phi < 1.0)
      scale_ratio_phi = 1.0;

    //         |   Delta     |
    //  N =log | ----------- |
    //        2|  lambda_nu  |
   Nphi = log(scale_ratio_phi)/log(2.0);
   if (Nphi < 0.0)
      dserror("Something went wrong when calculating N!");
  }
  else
   dserror("Multifractal subgrid-scales for loma with calculation of N, only!");

  // STEP 2: calculate D

  // here, we have to distinguish three different cases:
  // Pr ~ 1 : fluid and scalar field have the nearly the same cutoff (usual case)
  //          k^(-5/3) scaling -> gamma = 4/3
  // Pr >> 1: (i)  cutoff in the inertial-convective range (Nvel>0, tricky!)
  //               k^(-5/3) scaling in the inertial-convective range
  //               k^(-1) scaling in the viscous-convective range
  //          (ii) cutoff in the viscous-convective range (fluid field fully resolved, easier)
  //               k^(-1) scaling -> gamma = 2
  // rare:
  // Pr << 1: scatra field could be fully resolved, not necessary
  //          k^(-5/3) scaling -> gamma = 4/3
  // Remark: case 2.(i) not implemented, yet

#ifndef TESTING
  double gamma = 0.0;
  if (Pr < 2.0) // Pr <= 1, i.e., case 1 and 3
    gamma = 4.0/3.0;
  else if (Pr > 2.0 and Nvel[0]<1.0) // Pr >> 1, i.e., case 2 (ii)
    gamma = 2.0;
  else if (Pr > 2.0 and (Nvel[0]>=1.0 and Nvel[0]<Nphi))
  {
    gamma = 2.0;
//    std::cout << "Pr:" << Pr << std::endl;
//    std::cout << "Nvel:" << Nvel[0] << "  Nphi  " << Nphi << std::endl;
//    dserror("Inertial-convective and viscous-convective range?");
  }
  else
    dserror("Could not determine gamma!");

  //
  //   Phi    |       1                |
  //  kappa = | ---------------------- |
  //          |  1 - alpha ^ (-gamma)  |
  //
  double kappa_phi = 1.0/(1.0-pow(alpha,-gamma));

  //                                                             1
  //       Phi    Phi                       |                   |2
  //  D = Csgs * kappa * 2 ^ (-gamma*N/2) * | 2 ^ (gamma*N) - 1 |
  //                                        |                   |
  //
  D_mfs = Csgs_sgphi *sqrt(kappa_phi) * pow(2.0,-gamma*Nphi/2.0) * sqrt((pow(2.0,gamma*Nphi)-1));
//    if (eid_ == 100)
//     std::cout << "D  " << std::setprecision(10) << D_mfs << std::endl;
#endif

  // second implementation for tests on cluster
# ifdef TESTING
  double fac = 1.0;
# if 1
  // calculate near-wall correction
  if (nwl)
  {
    // not yet calculated, estimate norm of strain rate
    if (calc_N or refvel != INPAR::FLUID::strainrate)
    {
      strainnorm = GetStrainRate(econvelnp_);
      strainnorm /= sqrt(2.0);
    }

    // get Re from strain rate
    double Re_ele_str = strainnorm * hk * hk * densnp_[0] / visc_;
    if (Re_ele_str < 0.0)
      dserror("Something went wrong!");
    // ensure positive values
    if (Re_ele_str < 1.0)
       Re_ele_str = 1.0;

    // calculate corrected Csgs
    //           -3/16
    //  *(1 - (Re)   )
    //
    fac = (1-pow(Re_ele_str,-3.0/16.0)); //*pow(Pr,-1.0/8.0));
  }
#endif

// Pr <= 1
# if 1
  double gamma = 0.0;
  gamma = 4.0/3.0;
  double kappa_phi = 1.0/(1.0-pow(alpha,-gamma));
  D_mfs = Csgs_sgphi *sqrt(kappa_phi) * pow(2.0,-gamma*Nphi/2.0) * sqrt((pow(2.0,gamma*Nphi)-1)) * fac;
#endif

// Pr >> 1: cutoff viscous-convective
# if 0
  double gamma = 0.0;
  gamma = 2.0;
  double kappa_phi = 1.0/(1.0-pow(alpha,-gamma));
  D_mfs = Csgs_sgphi *sqrt(kappa_phi) * pow(2.0,-gamma*Nphi/2.0) * sqrt((pow(2.0,gamma*Nphi)-1)) * fac;
#endif

// Pr >> 1: cutoff inertial-convective
#if 0
  double gamma1 = 0.0;
  gamma1 = 4.0/3.0;
  double gamma2 = 0.0;
  gamma2 = 2.0;
  double kappa_phi = 1.0/(1.0-pow(alpha,-gamma1));
  D_mfs = Csgs_sgphi * sqrt(kappa_phi) * pow(2.0,-gamma2*Nphi/2.0) * sqrt((pow(2.0,gamma1*Nvel[dim])-1)+4.0/3.0*(PI/hk)*(pow(2.0,gamma2*Nphi)-pow(2.0,gamma2*Nvel[dim]))) * fac;
#endif

#endif

  return;
}


/*----------------------------------------------------------------------*
 | calculate reference length for multifractal subgrid-scales           |
 |                                                      rasthofer 09/12 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
double DRT::ELEMENTS::ScaTraImpl<distype>::CalcRefLength(
        const enum INPAR::FLUID::RefLength reflength,
        const double vol)
{
  // calculate characteristic element length
  double hk = 1.0e+10;
  // cf. stabilization parameters
  switch (reflength)
  {
    case INPAR::FLUID::streamlength:
    {
      // a) streamlength due to Tezduyar et al. (1992)
      // get norm of velocity
      const double vel_norm = convelint_.Norm2();
      // normed velocity vector
      LINALG::Matrix<nsd_,1> velino(true);
      if (vel_norm>=1e-6) velino.Update(1.0/vel_norm,convelint_);
      else
      {
        velino.Clear();
        velino(0,0) = 1.0;
      }
      LINALG::Matrix<nen_,1> tmp;
      tmp.MultiplyTN(derxy_,velino);
      const double val = tmp.Norm1();
      hk = 2.0/val;

      break;
    }
    case INPAR::FLUID::sphere_diameter:
    {
      // b) volume-equivalent diameter
      hk = std::pow((6.*vol/M_PI),(1.0/3.0))/sqrt(3.0);

      break;
    }
    case INPAR::FLUID::cube_edge:
    {
      // c) qubic element length
      hk = std::pow(vol,(1.0/nsd_));

      break;
    }
    case INPAR::FLUID::metric_tensor:
    {
      /*          +-           -+   +-           -+   +-           -+
                |             |   |             |   |             |
                |  dr    dr   |   |  ds    ds   |   |  dt    dt   |
          G   = |  --- * ---  | + |  --- * ---  | + |  --- * ---  |
           ij   |  dx    dx   |   |  dx    dx   |   |  dx    dx   |
                |    i     j  |   |    i     j  |   |    i     j  |
                +-           -+   +-           -+   +-           -+
      */
      LINALG::Matrix<3,3> G;

      for (int nn=0;nn<3;++nn)
      {
        for (int rr=0;rr<3;++rr)
        {
          G(nn,rr) = xij_(nn,0)*xij_(rr,0);
          for (int mm=1;mm<3;++mm)
          {
            G(nn,rr) += xij_(nn,mm)*xij_(rr,mm);
          }
        }
      }

      /*          +----
                 \
        G : G =   +   G   * G
        -   -    /     ij    ij
        -   -   +----
                 i,j
     */
     double normG = 0;
     for (int nn=0;nn<3;++nn)
     {
       for (int rr=0;rr<3;++rr)
       {
         normG+=G(nn,rr)*G(nn,rr);
       }
     }
     hk = std::pow(normG,-0.25);

     break;
  }
  case INPAR::FLUID::gradient_based:
  {
    LINALG::Matrix<nsd_,nsd_> convderxy;
    convderxy.MultiplyNT(econvelnp_,derxy_);
    LINALG::Matrix<3,1> normed_velgrad;

    for (int rr=0;rr<3;++rr)
    {
      normed_velgrad(rr)=sqrt(convderxy(0,rr)*convderxy(0,rr)
                              +
                              convderxy(1,rr)*convderxy(1,rr)
                              +
                              convderxy(2,rr)*convderxy(2,rr));
    }
    double norm=normed_velgrad.Norm2();

    // normed gradient
    if (norm>1e-6)
    {
      for (int rr=0;rr<3;++rr)
      {
        normed_velgrad(rr)/=norm;
      }
    }
    else
    {
      normed_velgrad(0) = 1.;
      for (int rr=1;rr<3;++rr)
      {
        normed_velgrad(rr)=0.0;
      }
    }

    // get length in this direction
    double val = 0.0;
    for (int rr=0;rr<nen_;++rr) /* loop element nodes */
    {
      val += fabs( normed_velgrad(0)*derxy_(0,rr)
                  +normed_velgrad(1)*derxy_(1,rr)
                  +normed_velgrad(2)*derxy_(2,rr));
    } /* end of loop over element nodes */

    hk = 2.0/val;

    break;
  }
  default:
    dserror("Unknown length");
  } // switch reflength
  if (hk == 1.0e+10)
   dserror("Something went wrong!");

  return hk;
}


/*----------------------------------------------------------------------*
 | output of model parameters                           rasthofer 09/12 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraImpl<distype>::StoreModelParametersForOutput(
  const bool                            isowned,
  ParameterList&                        turbulencelist,
  const int                             nlayer,
  const double                          tpn)
{
  if (isowned)
  {
    if (turbmodel_ == INPAR::FLUID::dynamic_smagorinsky)
    {
      if (turbulencelist.get<string>("TURBULENCE_APPROACH", "none") == "CLASSICAL_LES")
      {
        if (turbulencelist.get<string>("CANONICAL_FLOW","no")
            =="channel_flow_of_height_2" or
            turbulencelist.get<string>("CANONICAL_FLOW","no")
            =="scatra_channel_flow_of_height_2" or
            turbulencelist.get<string>("CANONICAL_FLOW","no")
            =="loma_channel_flow_of_height_2")
        {
           // calculate Prt form (Cs*h)^2/Prt
          if (tpn>1.0E-16)
          {
            // get dynamically estimated Smagorinsky constant from fluid element, i.e., (Cs*h)^2
            double Cs_delta_sq = (*(turbulencelist.get<RCP<vector<double> > >("global_Cs_delta_sq_sum")))[nlayer];
            // since Cs_delta_sq contains the sum over all elements of this layer,
            // we have to divide by the number of elements of this layer
            int numele_layer = turbulencelist.get<int>("numele_layer");
            (*(turbulencelist.get<RCP<vector<double> > >("local_Prt_sum")))[nlayer]+=(Cs_delta_sq/numele_layer)/tpn;
          }
          else
            (*(turbulencelist.get<RCP<vector<double> > >("local_Prt_sum")))         [nlayer]+=0.0;

          // set (Cs*h)^2/Prt and diffeff for output
          (*(turbulencelist.get<RCP<vector<double> > >("local_Cs_delta_sq_Prt_sum")))[nlayer]+=tpn;
          if (numscal_>1) dserror("One scalar assumed for dynamic Smagorinsky model!");
          (*(turbulencelist.get<RCP<vector<double> > >("local_diffeff_sum")))    [nlayer]+=diffus_[0];
        }
      }
    }
  }

  return;
}


template class DRT::ELEMENTS::ScaTraImpl<DRT::Element::hex8>;
//template class DRT::ELEMENTS::ScaTraImpl<DRT::Element::hex20>;
template class DRT::ELEMENTS::ScaTraImpl<DRT::Element::hex27>;
template class DRT::ELEMENTS::ScaTraImpl<DRT::Element::tet4>;
//template class DRT::ELEMENTS::ScaTraImpl<DRT::Element::tet10>;
template class DRT::ELEMENTS::ScaTraImpl<DRT::Element::wedge6>;
template class DRT::ELEMENTS::ScaTraImpl<DRT::Element::pyramid5>;
template class DRT::ELEMENTS::ScaTraImpl<DRT::Element::quad4>;
//template class DRT::ELEMENTS::ScaTraImpl<DRT::Element::quad8>;
template class DRT::ELEMENTS::ScaTraImpl<DRT::Element::quad9>;
template class DRT::ELEMENTS::ScaTraImpl<DRT::Element::tri3>;
//template class DRT::ELEMENTS::ScaTraImpl<DRT::Element::tri6>;
template class DRT::ELEMENTS::ScaTraImpl<DRT::Element::line2>;
template class DRT::ELEMENTS::ScaTraImpl<DRT::Element::line3>;
template class DRT::ELEMENTS::ScaTraImpl<DRT::Element::nurbs9>;
template class DRT::ELEMENTS::ScaTraImpl<DRT::Element::nurbs27>;
