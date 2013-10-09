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


/*-----------------------------------------------------------------------------*
 | calculate filtered quantities for dynamic Smagorinsky model  rasthofer 08/12|
 *-----------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraImpl<distype>::scatra_apply_box_filter(
  const double               thermpress,
  double&                    dens_hat,
  double&                    temp_hat,
  double&                    dens_temp_hat,
  double&                    phi2_hat,
  double&                    phiexpression_hat,
  Teuchos::RCP<std::vector<double> >       vel_hat,
  Teuchos::RCP<std::vector<double> >       densvel_hat,
  Teuchos::RCP<std::vector<double> >       densveltemp_hat,
  Teuchos::RCP<std::vector<double> >       densstraintemp_hat,
  Teuchos::RCP<std::vector<double> >       phi_hat,
  Teuchos::RCP<std::vector<std::vector<double> > > alphaijsc_hat,
  double&                    volume,
  const DRT::Element*        ele)
{

  // do preparations first
  // ---------------------------------------------
  LINALG::Matrix<nsd_,nsd_> vderxy (true);
  double alpha2 = 0.0;
  // use one-point Gauss rule to do calculations at the element center
  DRT::UTILS::IntPointsAndWeights<nsd_> intpoints(SCATRA::DisTypeToStabGaussRule<distype>::rule);

  volume = EvalShapeFuncAndDerivsAtIntPoint(intpoints,0,0);

  // get phi at integration point
  phi_[0] = funct_.Dot(ephinp_[0]);

  // get material
  RCP<MAT::Material> material = ele->Material();
  if (material->MaterialType() == INPAR::MAT::m_scatra)
  {
    //access fluid discretization
    RCP<DRT::Discretization> fluiddis = Teuchos::null;
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
    (*phi_hat)[rr] = tmp;
    phi2_hat += tmp*gradphi_(rr);
  }

  //calculate vreman part
  if (turbmodel_ == INPAR::FLUID::dynamic_vreman)
  {
    double beta00;
    double beta11;
    double beta22;
    double beta01;
    double beta02;
    double beta12;
    double bbeta;
    double hk2=pow(volume,(2.0/3.0));


    for (int nn=0;nn<nsd_;++nn)
    {
      for (int rr=0;rr<nsd_;++rr)
      {
        vderxy(nn,rr)=derxy_(rr,0)*evelnp_(nn,0);
        for (int mm=1;mm<nen_;++mm)
        {
          vderxy(nn,rr)+=derxy_(rr,mm)*evelnp_(nn,mm);
        }
        (*alphaijsc_hat)[rr][nn]=vderxy(nn,rr); //change indices to make compatible to paper
        alpha2 += vderxy(nn,rr)*vderxy(nn,rr);
      }
    }

    beta00=hk2 * (*alphaijsc_hat)[0][0] * (*alphaijsc_hat)[0][0] + hk2 * (*alphaijsc_hat)[1][0] * (*alphaijsc_hat)[1][0] + hk2 * (*alphaijsc_hat)[2][0] * (*alphaijsc_hat)[2][0];
    beta11=hk2 * (*alphaijsc_hat)[0][1] * (*alphaijsc_hat)[0][1] + hk2 * (*alphaijsc_hat)[1][1] * (*alphaijsc_hat)[1][1] + hk2 * (*alphaijsc_hat)[2][1] * (*alphaijsc_hat)[2][1];
    beta22=hk2 * (*alphaijsc_hat)[0][2] * (*alphaijsc_hat)[0][2] + hk2 * (*alphaijsc_hat)[1][2] * (*alphaijsc_hat)[1][2] + hk2 * (*alphaijsc_hat)[2][2] * (*alphaijsc_hat)[2][2];
    beta01=hk2 * (*alphaijsc_hat)[0][0] * (*alphaijsc_hat)[0][1] + hk2 * (*alphaijsc_hat)[1][0] * (*alphaijsc_hat)[1][1] + hk2 * (*alphaijsc_hat)[2][0] * (*alphaijsc_hat)[2][1];
    beta02=hk2 * (*alphaijsc_hat)[0][0] * (*alphaijsc_hat)[0][2] + hk2 * (*alphaijsc_hat)[1][0] * (*alphaijsc_hat)[1][2] + hk2 * (*alphaijsc_hat)[2][0] * (*alphaijsc_hat)[2][2];
    beta12=hk2 * (*alphaijsc_hat)[0][1] * (*alphaijsc_hat)[0][2] + hk2 * (*alphaijsc_hat)[1][1] * (*alphaijsc_hat)[1][2] + hk2 * (*alphaijsc_hat)[2][1] * (*alphaijsc_hat)[2][2];

    bbeta = beta00 * beta11 - beta01 * beta01
          + beta00 * beta22 - beta02 * beta02
          + beta11 * beta22 - beta12 * beta12;
    if (alpha2 < 1.0e-12)
      (phiexpression_hat)=0.0;
    else
      (phiexpression_hat) = (phi2_hat) * sqrt(bbeta/alpha2);

    for (int nn=0;nn<nsd_;++nn)
    {
      for (int rr=0;rr<nsd_;++rr)
      {
        (*alphaijsc_hat)[rr][nn]*=volume;
      }
    }

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


template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraImpl<distype>::        scatra_calc_vreman_dt(
  RCP<Epetra_MultiVector>& col_filtered_phi,
  RCP<Epetra_Vector>& col_filtered_phi2              ,
  RCP<Epetra_Vector>&   col_filtered_phiexpression         ,
  RCP<Epetra_MultiVector>& col_filtered_alphaijsc,
  double& dt_numerator,
  double& dt_denominator,
  const DRT::Element*       ele)
  {
  double phi_hat2=0.0;
  LINALG::Matrix<9,nen_> ealphaijsc_hat(true)                 ;
  LINALG::Matrix<nsd_,nen_> ephi_hat(true)                            ;
  LINALG::Matrix<1,nen_> ephi2_hat(true)                           ;
  LINALG::Matrix<1,nen_> ephiexpression_hat(true)                 ;
  LINALG::Matrix<nsd_,nsd_> alphaijsc_hat(true);
  LINALG::Matrix<nsd_,1> phi_hat(true);
  LINALG::Matrix<1,1> phi2_hat(true);
  LINALG::Matrix<1,1> phiexpression_hat(true);

  DRT::UTILS::ExtractMyNodeBasedValues(ele,ephi_hat,col_filtered_phi,nsd_);
  DRT::UTILS::ExtractMyNodeBasedValues(ele,ephi2_hat,col_filtered_phi2,1);
  DRT::UTILS::ExtractMyNodeBasedValues(ele,ephiexpression_hat,col_filtered_phiexpression,1);
  DRT::UTILS::ExtractMyNodeBasedValues(ele,ealphaijsc_hat,col_filtered_alphaijsc,nsd_*nsd_);
  // use one-point Gauss rule to do calculations at the element center
  DRT::UTILS::IntPointsAndWeights<nsd_> intpoints(SCATRA::DisTypeToStabGaussRule<distype>::rule);
  double volume = EvalShapeFuncAndDerivsAtIntPoint(intpoints,0,0);

  phi_hat.Multiply(ephi_hat,funct_);
  phi2_hat.Multiply(ephi2_hat,funct_);
  phiexpression_hat.Multiply(ephiexpression_hat,funct_);
  for (int nn=0;nn<3;++nn)
  {
    for (int rr=0;rr<3;++rr)
    {
      int index = 3*nn+rr;
      alphaijsc_hat(nn,rr)=funct_(0)*ealphaijsc_hat(index,0);

      for (int mm=1;mm<nen_;++mm)
      {
        alphaijsc_hat(nn,rr)+=funct_(mm)*ealphaijsc_hat(index,mm);
      }
    }
  }

  for (int rr=0;rr<nsd_;++rr)
  {
    phi_hat2+=phi_hat(rr,0)*phi_hat(rr,0);
  }

  dt_denominator=phi2_hat(0,0)-phi_hat2;

  //calculate vreman part
  {

    double beta00;
    double beta11;
    double beta22;
    double beta01;
    double beta02;
    double beta12;
    double bbeta;
    double hk2=3*pow(volume,(2.0/3.0));
    double alpha2=0.0;
    double phiexpressionf_hat=0.0;

    for (int nn=0;nn<nsd_;++nn)
    {
      for (int rr=0;rr<nsd_;++rr)
      {
        alpha2 += alphaijsc_hat(nn,rr)*alphaijsc_hat(nn,rr);
      }
    }

    beta00=hk2 * alphaijsc_hat(0,0) * alphaijsc_hat(0,0) + hk2 * alphaijsc_hat(1,0) * alphaijsc_hat(1,0) + hk2 * alphaijsc_hat(2,0) * alphaijsc_hat(2,0);
    beta11=hk2 * alphaijsc_hat(0,1) * alphaijsc_hat(0,1) + hk2 * alphaijsc_hat(1,1) * alphaijsc_hat(1,1) + hk2 * alphaijsc_hat(2,1) * alphaijsc_hat(2,1);
    beta22=hk2 * alphaijsc_hat(0,2) * alphaijsc_hat(0,2) + hk2 * alphaijsc_hat(1,2) * alphaijsc_hat(1,2) + hk2 * alphaijsc_hat(2,2) * alphaijsc_hat(2,2);
    beta01=hk2 * alphaijsc_hat(0,0) * alphaijsc_hat(0,1) + hk2 * alphaijsc_hat(1,0) * alphaijsc_hat(1,1) + hk2 * alphaijsc_hat(2,0) * alphaijsc_hat(2,1);
    beta02=hk2 * alphaijsc_hat(0,0) * alphaijsc_hat(0,2) + hk2 * alphaijsc_hat(1,0) * alphaijsc_hat(1,2) + hk2 * alphaijsc_hat(2,0) * alphaijsc_hat(2,2);
    beta12=hk2 * alphaijsc_hat(0,1) * alphaijsc_hat(0,2) + hk2 * alphaijsc_hat(1,1) * alphaijsc_hat(1,2) + hk2 * alphaijsc_hat(2,1) * alphaijsc_hat(2,2);

    bbeta = beta00 * beta11 - beta01 * beta01
          + beta00 * beta22 - beta02 * beta02
          + beta11 * beta22 - beta12 * beta12;
    if (alpha2 < 1.0e-12)
      phiexpressionf_hat=0.0;
    else
      phiexpressionf_hat = (phi_hat2) * sqrt(bbeta/alpha2);

    dt_numerator=phiexpressionf_hat-phiexpression_hat(0,0);
//    std::cout << "scatra_ele_impl_turb_service.cpp  dt_numerator  " << dt_numerator << std::endl;
//    std::cout << "scatra_ele_impl_turb_service.cpp  dt_denominator  " << dt_denominator << std::endl;
//    std::cout << "scatra_ele_impl_turb_service.cpp  phi_hat2  " << phi_hat2 << std::endl;
//    std::cout << "scatra_ele_impl_turb_service.cpp  phi2_hat  " << phi2_hat << std::endl;
//    std::cout << "scatra_ele_impl_turb_service.cpp  alpha2  " << alpha2 << std::endl;
  }


    return;
  }


/*----------------------------------------------------------------------------------*
 | calculate mean turbulent Prandtl number                           rasthofer 08/12|
 *----------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraImpl<distype>::GetMeanPrtOfHomogenousDirection(
  Teuchos::ParameterList&    turbmodelparams,
  double&                    inv_Prt,
  int&                       nlayer
)
{
  if(nsd_ != 3)
    dserror("turbulence and 3D flow !");

  if (turbmodelparams.get<std::string>("HOMDIR","not_specified") !=  "not_specified")
  {
    RCP<std::vector<double> > averaged_LkMk = turbmodelparams.get<RCP<std::vector<double> > >("averaged_LkMk_");
    RCP<std::vector<double> > averaged_MkMk = turbmodelparams.get<RCP<std::vector<double> > >("averaged_MkMk_");

    // get homogeneous direction
    std::string homdir = turbmodelparams.get<std::string>("HOMDIR","not_specified");

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
    if (homdir == "xyz")
    {
      nlayer = 0;
    }
    else if (homdir == "xy" or homdir == "xz" or homdir == "yz")
    {
      RCP<std::vector<double> > planecoords = turbmodelparams.get<RCP<std::vector<double> > >("planecoords_");
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
      RCP<std::vector<double> > dir1coords = turbmodelparams.get<RCP<std::vector<double> > >("dir1coords_");
      RCP<std::vector<double> > dir2coords = turbmodelparams.get<RCP<std::vector<double> > >("dir2coords_");
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
      dserror("Homogeneous directions not supported!");

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
  if (assgd)
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
	  if (eid_ == 0)
          {
            std::cout << "WARNING: Nonlinear isotropic artificial diffusion according to Hughes et al. (1986)\n";
            std::cout << "         is implemented based on the exact tau for 1D stationary problems!" << std::endl; 
          }
          // remark on this warning:
          // 1. Here, tau is calculated based on the exact formula for 1-D stationary problems. This is inconsistent
          //    if the problem is not 1-D and/or other definitions than the ecaxt tau are chosen in the input file.
          //    Consistently, one has to use the same definition here.
          // 2. Instead of using sigma = tau_bhbar, Hughes et al. suggested to use sigma = tau_bhbar - tau to not double
          //    the SUPG stabilization. This is another inconsitent aspect on this implementation. To have the right tau
          //    here (i.e, not the one the last gauss point or even last step), one has to calaculate tau first. Then,
          //    sigma and, hence, the addition diffusion is computed based on this tau. Next, tau is recomputed with the
          //    diffusivity repaced by the original (physical) diffusivity plus the estimated artificial diffusivity. This
          //    is pobably not a good choice, because, first, tau is considered in the estimation of the artificial diffusion
          //    and then this artificial diffusion is incorporated into tau. This would reduce the effect. Perhaps, one 
          //    should either consider tau in sigma or the artificial diffusion in tau. When changing this aspect, be aware
          //    that tau has to be computed after the subgrid-scale velocity has been calcuated since, for this calculation
          //    tau is overwritten by its value in the fluid field. Note that similar considerations may also hold for
          //    the methods by do Carmo and Almeida.

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
        case INPAR::SCATRA::assgd_tezduyar_wo_phizero:
        {
          // velocity norm
          const double vel_norm = convelint_.Norm2();
	  
          // calculate stream length
          // according to John and Knobloch stream length in direction of b_h^par should be used
          //const double h_stream = CalcCharEleLength(vol,vel_norm); 

          // get norm of velocity vector b_h^par
          const double vel_norm_bhpar = abs(conv_phi_[k]/grad_norm);

          // compute stabilization parameter based on b_h^par
          // (so far, only exact formula for stationary 1-D implemented)

          // compute sigma (version 1 according to John and Knobloch (2007))
          if (whichassgd == INPAR::SCATRA::assgd_tezduyar_wo_phizero)
          {
            if (vel_norm > EPS10)
              sigma = (h/vel_norm)*(1.0-(vel_norm_bhpar/vel_norm));
          }
          else
          {
            // compute sigma (version 2 according to John and Knobloch (2007))
            // setting scaling phi_0=1.0 as in John and Knobloch (2007)
            const double phi0 = 1.0;
            if (vel_norm > EPS10)
              sigma = (h*h*grad_norm/(vel_norm*phi0))*(1.0-(vel_norm_bhpar/vel_norm));
          }

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
          else
          {
            if (abs(scatrares_[k]) > EPS10)
              zeta = std::max(1.0,(conv_phi_[k]/scatrares_[k]));
          }

          // compute sigma
          if (vel_norm_zh > EPS10)
            sigma = tau_[k]*std::max(0.0,(vel_norm/vel_norm_zh)-zeta);

          // set specific term to residual
          specific_term = scatrares_[k];
        }
        break;
        default: dserror("unknown type of all-scale subgrid diffusivity\n"); break;
        } //switch (whichassgd)

        // computation of subgrid diffusivity
        sgdiff_[k] = sigma*scatrares_[k]*specific_term/(grad_norm*grad_norm);
        if (sgdiff_[k] < 0.0)
        {
          std::cout << "WARNING: isotropic artificial diffusion sgdiff < 0.0\n";
          std::cout << "         -> set sgdiff to abs(sgdiff)!" << std::endl;
          sgdiff_[k] = abs(sigma*scatrares_[k]*specific_term/(grad_norm*grad_norm));
        }
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
  else if (turbmodel_ == INPAR::FLUID::dynamic_vreman)
  {
    double beta00;
    double beta11;
    double beta22;
    double beta01;
    double beta02;
    double beta12;
    double bbeta;
    double alphavreman;
    double hkxpow2;
    double hkypow2;
    double hkzpow2;
    double sgviscwocv=0.0;
    double Dt=Cs;

    LINALG::Matrix<nsd_,nsd_> velderxy(true);

    velderxy.MultiplyNT(econvelnp_,derxy_);

    //- cube root of element volume

    hkxpow2=pow(vol,(2.0/3.0));
    hkypow2=hkxpow2;
    hkzpow2=hkxpow2;

    beta00=hkxpow2 * velderxy(0,0) * velderxy(0,0) + hkypow2 * velderxy(0,1) * velderxy(0,1) + hkzpow2 * velderxy(0,2) * velderxy(0,2);
    beta11=hkxpow2 * velderxy(1,0) * velderxy(1,0) + hkypow2 * velderxy(1,1) * velderxy(1,1) + hkzpow2 * velderxy(1,2) * velderxy(1,2);
    beta22=hkxpow2 * velderxy(2,0) * velderxy(2,0) + hkypow2 * velderxy(2,1) * velderxy(2,1) + hkzpow2 * velderxy(2,2) * velderxy(2,2);
    beta01=hkxpow2 * velderxy(0,0) * velderxy(1,0) + hkypow2 * velderxy(0,1) * velderxy(1,1) + hkzpow2 * velderxy(0,2) * velderxy(1,2);
    beta02=hkxpow2 * velderxy(0,0) * velderxy(2,0) + hkypow2 * velderxy(0,1) * velderxy(2,1) + hkzpow2 * velderxy(0,2) * velderxy(2,2);
    beta12=hkxpow2 * velderxy(1,0) * velderxy(2,0) + hkypow2 * velderxy(1,1) * velderxy(2,1) + hkzpow2 * velderxy(1,2) * velderxy(2,2);

    bbeta = beta00 * beta11 - beta01 * beta01
          + beta00 * beta22 - beta02 * beta02
          + beta11 * beta22 - beta12 * beta12;

    alphavreman = velderxy(0,0) * velderxy(0,0)
                + velderxy(0,1) * velderxy(0,1)
                + velderxy(0,2) * velderxy(0,2)
                + velderxy(1,0) * velderxy(1,0)
                + velderxy(1,1) * velderxy(1,1)
                + velderxy(1,2) * velderxy(1,2)
                + velderxy(2,0) * velderxy(2,0)
                + velderxy(2,1) * velderxy(2,1)
                + velderxy(2,2) * velderxy(2,2);

    if(alphavreman<1.0E-12)
      sgviscwocv=0.0;
    else
      sgviscwocv=   sqrt(bbeta / alphavreman);


    //remark: Cs corresponds to Dt, calculated in the vreman class
    //        The vreman constant Cv is not required here, since it cancelles out with the
    //        vreman constant omitted during the calculation of D_t
    if (Dt<=1.0E-12) sgdiff_[k]=0.0;
    else
      sgdiff_[k] = densnp_[k] * sgviscwocv / (Dt);

  }


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
    const bool                                  nwl_scatra, ///< flag to activate near-wall limit for scalar field
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
  std::vector<double> Nvel (3);
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
    hk = CalcRefLength(reflength,vol);

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
        strainnorm /= sqrt(2.0); //cf. Burton & Dahm 2008
        Re_ele = strainnorm * hk * hk * densnp_[0] / visc_;
        break;
      }
      default:
        dserror("Unknown velocity!");
        break;
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
  double Cai_phi = 0.0;
  if (nwl)
  {
    // not yet calculated, estimate norm of strain rate
    if ((not calc_N) or (refvel != INPAR::FLUID::strainrate))
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

    // store Cai for application to scalar field
    Cai_phi = (1-pow(Re_ele_str,-3.0/16.0));
  }

  // STEP 2: calculate B

  //       1                           1
  //       2   |       1              |2
  //  kappa  = | -------------------- |
  //           |  1 - alpha ^ (-4/3)  |
  //
  double kappa = 1.0/(1.0-pow(alpha,-4.0/3.0));

  //                  1                                     1
  //                  2                  |                 |2
  //  B = Csgs * kappa  * 2 ^ (-2*N/3) * | 2 ^ (4*N/3) - 1 |
  //                                     |                 |
  //
  for (int dim=0; dim<nsd_; dim++)
  {
    B_mfs(dim,0) = Csgs_vel_nw * sqrt(kappa) * pow(2.0,-2.0*Nvel[dim]/3.0) * sqrt((pow(2.0,4.0*Nvel[dim]/3.0)-1));
//    if (eid_ == 100)
//     std::cout << "B  " << std::setprecision (10) << B_mfs(dim,0) << std::endl;
  }

  //----------------------------------------------------------------
  // calculation of D for fine-scale scalar
  //----------------------------------------------------------------

  // STEP 1: determine N

  // calculate Prandtl number or Schmidt number (passive scalar)
  const double Pr = visc_/diffus_[k];

  // since there are differences in the physical behavior between low and high
  // Prandtl/Schmidt number regime, we define a limit 
  // to distinguish between the low and high Prandtl/Schmidt number regime
  // note: there is no clear definition of the ranges
  const double Pr_limit = 2.0;

  // allocate vector for parameter N
  double Nphi = 0.0;
  // ratio of diffusive scale to element length
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
    if (Pr>Pr_limit) p =1.0/2.0;

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

  // caution: compared to the mfs-loma paper, gamma denotes gamma+1 here
  double gamma = 0.0;
  // special option for case 2 (i)
  bool two_ranges = false;
  if (Pr < Pr_limit){ // Pr <= 1, i.e., case 1 and 3
    gamma = 4.0/3.0;
  }
  else // Pr >> 1
  {
    if (Nvel[0]<1.0){ // Pr >> 1 and fluid fully resolved, i.e., case 2 (ii)
      gamma = 2.0;
    }
    else // Pr >> 1 and fluid not fully resolved, i.e., case 2 (i)
    {
      if (Nvel[0] > Nphi)
      {
        std::cout << "Nvel   " << Nvel[0] << std::endl;
        std::cout << "Nphi   " << Nphi << std::endl;
        dserror("Nvel < Nphi expected!");
      }
      // here different options are possible
      // 1) we assume k^(-5/3) for the complete range
      gamma = 4.0/3.0;
#if 0
      // 2) we assume k^(-1) for the complete range
      gamma = 2.0;
      // 3) we take both ranges into account
      two_ranges = true;
#endif
    }
  }

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
  if (not two_ranges) // usual case
    D_mfs = Csgs_sgphi *sqrt(kappa_phi) * pow(2.0,-gamma*Nphi/2.0) * sqrt((pow(2.0,gamma*Nphi)-1));
  else
  {
    double gamma1 = 4.0/3.0;
    double gamma2 = 2.0;
    kappa_phi = 1.0/(1.0-pow(alpha,-gamma1));
    D_mfs = Csgs_sgphi * sqrt(kappa_phi) * pow(2.0,-gamma2*Nphi/2.0) * sqrt((pow(2.0,gamma1*Nvel[0])-1)+2.0/3.0*pow((M_PI/hk),2.0/3.0)*(pow(2.0,gamma2*Nphi)-pow(2.0,gamma2*Nvel[0])));
  }

  // apply near-wall limit if required
  if (nwl_scatra and nwl)
  {
    D_mfs *= Cai_phi;
  }

//  if (eid_ == 100){
////    std::cout << "sqrt(kappa_phi)  " << std::setprecision(10) << sqrt(kappa_phi) << std::endl;
////    std::cout << "pow(2.0,-gamma*Nphi/2.0)  " << std::setprecision(10) << pow(2.0,-gamma*Nphi/2.0) << std::endl;
////    std::cout << "sqrt((pow(2.0,gamma*Nphi)-1))  " << std::setprecision(10) << sqrt((pow(2.0,gamma*Nphi)-1)) << std::endl;
//    std::cout << "D  " << std::setprecision(10) << D_mfs << std::endl;
//    std::cout << "B  " << std::setprecision(10) << B_mfs(0,0) << "  " << B_mfs(1,0) << "  " << B_mfs(2,0) << "  " << std::endl;
//    if (nwl_scatra and nwl)
//     std::cout << "CsgsD  " << std::setprecision(10) << Csgs_sgphi*Cai_phi << std::endl;
//    else
//     std::cout << "CsgsD  " << std::setprecision(10) << Csgs_sgphi << std::endl;
//    std::cout << "CsgsB  " << std::setprecision(10) << Csgs_vel_nw << "  " << Csgs_sgvel << " " << Cai_phi << std::endl;
//  }

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
    break;
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
  Teuchos::ParameterList&               turbulencelist,
  const int                             nlayer,
  const double                          tpn)
{
  if (isowned)
  {
    if (turbmodel_ == INPAR::FLUID::dynamic_smagorinsky)
    {
      if (turbulencelist.get<std::string>("TURBULENCE_APPROACH", "none") == "CLASSICAL_LES")
      {
        if (turbulencelist.get<std::string>("CANONICAL_FLOW","no")
            =="channel_flow_of_height_2" or
            turbulencelist.get<std::string>("CANONICAL_FLOW","no")
            =="scatra_channel_flow_of_height_2" or
            turbulencelist.get<std::string>("CANONICAL_FLOW","no")
            =="loma_channel_flow_of_height_2")
        {
           // calculate Prt form (Cs*h)^2/Prt
          if (tpn>1.0E-16)
          {
            // get dynamically estimated Smagorinsky constant from fluid element, i.e., (Cs*h)^2
            double Cs_delta_sq = (*(turbulencelist.get<RCP<std::vector<double> > >("global_Cs_delta_sq_sum")))[nlayer];
            // since Cs_delta_sq contains the sum over all elements of this layer,
            // we have to divide by the number of elements of this layer
            int numele_layer = turbulencelist.get<int>("numele_layer");
            (*(turbulencelist.get<RCP<std::vector<double> > >("local_Prt_sum")))[nlayer]+=(Cs_delta_sq/numele_layer)/tpn;
          }
          else
            (*(turbulencelist.get<RCP<std::vector<double> > >("local_Prt_sum")))         [nlayer]+=0.0;

          // set (Cs*h)^2/Prt and diffeff for output
          (*(turbulencelist.get<RCP<std::vector<double> > >("local_Cs_delta_sq_Prt_sum")))[nlayer]+=tpn;
          if (numscal_>1) dserror("One scalar assumed for dynamic Smagorinsky model!");
          (*(turbulencelist.get<RCP<std::vector<double> > >("local_diffeff_sum")))    [nlayer]+=diffus_[0];
        }
      }
    }
  }

  return;
}


/*----------------------------------------------------------------------*
 | additional output for turbulent channel flow         rasthofer 11/12 |
 | dissipation introduced by stabilization and turbulence models        |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraImpl<distype>::CalcDissipation(
     Teuchos::ParameterList&               params,
     DRT::Element*                         ele,
     const enum INPAR::SCATRA::ScaTraType  scatratype,
     DRT::Discretization&                  discretization,
     std::vector<int>&                     lm)
{
  // do some checks first
  if (numscal_!=1 or numdofpernode_!=1)
    dserror("CalcDissipation only for one scalar field!");

  //----------------------------------------------------------------------
  // preliminary set-up of parameters
  // ---------------------------------------------------------------------

  // set thermodynamic pressure and its time derivative as well as
  // flag for turbulence model if required
  turbmodel_ = INPAR::FLUID::no_model;
  Teuchos::ParameterList& turbulencelist = params.sublist("TURBULENCE MODEL");
  Teuchos::ParameterList& sgvisclist = params.sublist("SUBGRID VISCOSITY");
  Teuchos::ParameterList& mfslist = params.sublist("MULTIFRACTAL SUBGRID SCALES");

  if (scatratype == INPAR::SCATRA::scatratype_loma)
  {
    thermpressnp_ = params.get<double>("thermodynamic pressure");
    thermpressdt_ = params.get<double>("time derivative of thermodynamic pressure");
    if (is_genalpha_)
      thermpressam_ = params.get<double>("thermodynamic pressure at n+alpha_M");

    // update material with subgrid-scale scalar
    update_mat_ = params.get<bool>("update material");
  }

  if (scatratype == INPAR::SCATRA::scatratype_loma or
      scatratype == INPAR::SCATRA::scatratype_turbpassivesca)
  {
    // set flag for turbulence model
    if (turbulencelist.get<std::string>("PHYSICAL_MODEL") == "Smagorinsky")
      turbmodel_ = INPAR::FLUID::smagorinsky;
    if (turbulencelist.get<std::string>("PHYSICAL_MODEL") == "Dynamic_Smagorinsky")
      turbmodel_ = INPAR::FLUID::dynamic_smagorinsky;
    if (turbulencelist.get<std::string>("PHYSICAL_MODEL") == "Multifractal_Subgrid_Scales")
      turbmodel_ = INPAR::FLUID::multifractal_subgrid_scales;
    if (turbulencelist.get<std::string>("PHYSICAL_MODEL") == "Dynamic_Vreman")
      turbmodel_ = INPAR::FLUID::dynamic_vreman;
    // as the scalar field is constant in the turbulent inflow section
    // we do not need any turbulence model
    if (params.get<bool>("turbulent inflow",false))
    {
      if (SCATRA::InflowElement(ele))
        turbmodel_ = INPAR::FLUID::no_model;
    }
  }

  // set time integration
  Teuchos::ParameterList& timeintlist = params.sublist("TIME INTEGRATION");
  is_incremental_ = true; // be careful
  is_stationary_  = timeintlist.get<bool>("using stationary formulation"); //turbulence is instationary!
  is_genalpha_    = timeintlist.get<bool>("using generalized-alpha time integration");
  // get current time and time-step length
  const double time = timeintlist.get<double>("total time");
  const double dt   = params.get<double>("time-step length");
  // get time factor and alpha_F if required
  // one-step-Theta:    timefac = theta*dt
  // BDF2:              timefac = 2/3 * dt
  // generalized-alpha: timefac = alphaF * (gamma*/alpha_M) * dt
  double timefac = 1.0;
  double alphaF  = 1.0;
  if (not is_stationary_)
  {
    timefac = timeintlist.get<double>("time factor");
    if (is_genalpha_)
    {
      alphaF = timeintlist.get<double>("alpha_F");
      timefac *= alphaF;
    }
    if (timefac < 0.0) dserror("time factor is negative.");
  }
  else dserror("Turbulence is instationary!");

  // set parameters for stabilization
  Teuchos::ParameterList& stablist = params.sublist("STABILIZATION");
  // get definition for stabilization parameter tau
  whichtau_ = DRT::INPUT::IntegralValue<INPAR::SCATRA::TauType>(stablist,"DEFINITION_TAU");
  // do one check
  if (whichtau_ == INPAR::SCATRA::tau_exact_1d)
    dserror("exact stabilization parameter only available for stationary case");
  // get type of stabilization
  const INPAR::SCATRA::StabType stabtype = DRT::INPUT::IntegralValue<INPAR::SCATRA::StabType>(stablist,"STABTYPE");

  // set flags for subgrid-scale velocity and all-scale subgrid-diffusivity term
  // (default: "false" for both flags)
  const bool sgvel(DRT::INPUT::IntegralValue<int>(stablist,"SUGRVEL"));
  sgvel_ = sgvel;
  // select type of all-scale subgrid diffusivity if included
  // this is just to have a dummy value here
  const INPAR::SCATRA::AssgdType whichassgd
    = DRT::INPUT::IntegralValue<INPAR::SCATRA::AssgdType>(stablist,"DEFINITION_ASSGD");
  const bool assgd(DRT::INPUT::IntegralValue<int>(stablist,"ASSUGRDIFF"));
  if (assgd) dserror("All-scale subgrid-diffusivity approach not supported!");

  // set flags for potential evaluation of tau and material law at int. point
  const INPAR::SCATRA::EvalTau tauloc = DRT::INPUT::IntegralValue<INPAR::SCATRA::EvalTau>(stablist,"EVALUATION_TAU");
  tau_gp_ = (tauloc == INPAR::SCATRA::evaltau_integration_point); // set true/false
  const INPAR::SCATRA::EvalMat matloc = DRT::INPUT::IntegralValue<INPAR::SCATRA::EvalMat>(stablist,"EVALUATION_MAT");
  mat_gp_ = (matloc == INPAR::SCATRA::evalmat_integration_point); // set true/false

  // set flag for fine-scale subgrid diffusivity and perform some checks
  bool fssgd = false; //default
  const INPAR::SCATRA::FSSUGRDIFF whichfssgd = DRT::INPUT::get<INPAR::SCATRA::FSSUGRDIFF>(params, "fs subgrid diffusivity");
  if (whichfssgd == INPAR::SCATRA::fssugrdiff_artificial)
  {
    fssgd = true;
    // checks are removed since they have already been done in the usual element call for mat and rhs
    // check for solver type
    //if (is_incremental_) dserror("Artificial fine-scale subgrid-diffusivity approach only in combination with non-incremental solver so far!");
  }
  else if (whichfssgd == INPAR::SCATRA::fssugrdiff_smagorinsky_all or whichfssgd == INPAR::SCATRA::fssugrdiff_smagorinsky_small)
  {
    fssgd = true;
    // checks are removed since they have already been done in the usual element call for mat and rhs
    // check for solver type
    //if (not is_incremental_) dserror("Fine-scale subgrid-diffusivity approach using all/small-scale Smagorinsky model only in combination with incremental solver so far!");
  }
  // check for combination of all-scale and fine-scale subgrid diffusivity
  if (assgd and fssgd) dserror("No combination of all-scale and fine-scale subgrid-diffusivity approach currently possible!");

  is_reactive_ = false;

  // parameters for subgrid-diffusivity models
  double Cs(0.0);
  double tpn(1.0);
  // parameters for multifractal subgrid-scale modeling
  double Csgs_sgvel = 0.0;
  double alpha = 0.0;
  bool calc_N = true;
  double N_vel = 1.0;
  INPAR::FLUID::RefVelocity refvel = INPAR::FLUID::strainrate;
  INPAR::FLUID::RefLength reflength = INPAR::FLUID::cube_edge;
  double c_nu = 1.0;
  bool nwl = false;
  bool nwl_scatra = false;
  bool beta = 0.0;
  bool BD_gp = false;
  double Csgs_sgphi = 0.0;
  double c_diff = 1.0;
  // parameter for averaging (dynamic Smagorinsky)
  if (turbmodel_!=INPAR::FLUID::no_model or fssgd)
  {
    // get Smagorinsky constant and turbulent Prandtl number
    Cs  = sgvisclist.get<double>("C_SMAGORINSKY");
    tpn = sgvisclist.get<double>("C_TURBPRANDTL");
    if (tpn <= 1.0E-16) dserror("Turbulent Prandtl number should be larger than zero!");

    if (turbmodel_ == INPAR::FLUID::dynamic_smagorinsky)
    {
      // remark: for dynamic estimation, this returns (Cs*h)^2 / Pr_t
      Teuchos::RCP<Epetra_Vector> ele_prt = turbulencelist.get<Teuchos::RCP<Epetra_Vector> >("col_ele_Prt");
      const int id = ele->LID();
      tpn = (*ele_prt)[id];

      // when no averaging was done, we just keep the calculated (clipped) value
      int dummy=0;
      if (DRT::INPUT::IntegralValue<int>(sgvisclist,"C_SMAGORINSKY_AVERAGED"))
        GetMeanPrtOfHomogenousDirection(params.sublist("TURBULENCE MODEL"),tpn,dummy);
    }

    // get model parameters
    if (turbmodel_ == INPAR::FLUID::multifractal_subgrid_scales)
    {
      // necessary parameters for subgrid-scale velocity estimation
      Csgs_sgvel = mfslist.get<double>("CSGS");
      if (mfslist.get<std::string>("SCALE_SEPARATION") == "algebraic_multigrid_operator")
       alpha = 3.0;
      else dserror("Scale-separation method not supported!");
      calc_N = DRT::INPUT::IntegralValue<int>(mfslist,"CALC_N");
      N_vel = mfslist.get<double>("N");
      if (mfslist.get<std::string>("REF_VELOCITY") == "strainrate")
       refvel = INPAR::FLUID::strainrate;
      else if (mfslist.get<std::string>("REF_VELOCITY") == "resolved")
       refvel = INPAR::FLUID::resolved;
      else if (mfslist.get<std::string>("REF_VELOCITY") == "fine_scale")
       refvel = INPAR::FLUID::fine_scale;
      else
       dserror("Unknown velocity!");
      if (mfslist.get<std::string>("REF_LENGTH") == "cube_edge")
       reflength = INPAR::FLUID::cube_edge;
      else if (mfslist.get<std::string>("REF_LENGTH") == "sphere_diameter")
       reflength = INPAR::FLUID::sphere_diameter;
      else if (mfslist.get<std::string>("REF_LENGTH") == "streamlength")
       reflength = INPAR::FLUID::streamlength;
      else if (mfslist.get<std::string>("REF_LENGTH") == "gradient_based")
       reflength = INPAR::FLUID::gradient_based;
      else if (mfslist.get<std::string>("REF_LENGTH") == "metric_tensor")
       reflength = INPAR::FLUID::metric_tensor;
      else
       dserror("Unknown length!");
      c_nu = mfslist.get<double>("C_NU");
      nwl = DRT::INPUT::IntegralValue<int>(mfslist,"NEAR_WALL_LIMIT");
      // necessary parameters for subgrid-scale scalar estimation
      Csgs_sgphi = mfslist.get<double>("CSGS_PHI");
      c_diff = mfslist.get<double>("C_DIFF");
      if (DRT::INPUT::IntegralValue<int>(mfslist,"ADAPT_CSGS_PHI") and nwl)
      {
        double meanCai = mfslist.get<double>("meanCai");
        Csgs_sgphi *= meanCai;
      }
      nwl_scatra = DRT::INPUT::IntegralValue<int>(mfslist,"NEAR_WALL_LIMIT_CSGS_PHI");
      // general parameters
      beta = mfslist.get<double>("BETA");
      if (beta!=0.0) dserror("Lhs terms for mfs not included! Fixed-point iteration only!");
      if (mfslist.get<std::string>("EVALUATION_B") == "element_center")
      BD_gp = false;
      else if (mfslist.get<std::string>("EVALUATION_B") == "integration_point")
      BD_gp = true;
      else
        dserror("Unknown evaluation point!");
      if (mfslist.get<std::string>("CONVFORM") == "convective")
      mfs_conservative_ = false;
      else if (mfslist.get<std::string>("CONVFORM") == "conservative")
      mfs_conservative_ = true;
      else
        dserror("Unknown form of convective term!");
    }
  }


  //----------------------------------------------------------------------
  // get all nodal values
  // ---------------------------------------------------------------------

  // get velocity at nodes
  const RCP<Epetra_MultiVector> convelocity = params.get< RCP<Epetra_MultiVector> >("convective velocity field");
  // set econvelnp_ equal to evelnp_ since ale is not supported
  DRT::UTILS::ExtractMyNodeBasedValues(ele,econvelnp_,convelocity,nsd_);
  DRT::UTILS::ExtractMyNodeBasedValues(ele,evelnp_,convelocity,nsd_);

  // get data required for subgrid-scale velocity: acceleration and pressure
  if (sgvel_)
  {
    // check for matching flags
    if (not mat_gp_ or not tau_gp_)
     dserror("Evaluation of material and stabilization parameters need to be done at the integration points if subgrid-scale velocity is included!");

    const RCP<Epetra_MultiVector> accpre = params.get< RCP<Epetra_MultiVector> >("acceleration/pressure field");
    LINALG::Matrix<nsd_+1,nen_> eaccprenp;
    DRT::UTILS::ExtractMyNodeBasedValues(ele,eaccprenp,accpre,nsd_+1);

    // split acceleration and pressure values
    for (int i=0;i<nen_;++i)
    {
      for (int j=0;j<nsd_;++j)
      {
        eaccnp_(j,i) = eaccprenp(j,i);
      }
      eprenp_(i) = eaccprenp(nsd_,i);
    }
  }

  // extract local values from the global vectors
  Teuchos::RCP<const Epetra_Vector> hist = discretization.GetState("hist");
  Teuchos::RCP<const Epetra_Vector> phinp = discretization.GetState("phinp");
  if (hist==Teuchos::null || phinp==Teuchos::null)
    dserror("Cannot get state vector 'hist' and/or 'phinp'");
  std::vector<double> myhist(lm.size());
  std::vector<double> myphinp(lm.size());
  DRT::UTILS::ExtractMyValues(*hist,myhist,lm);
  DRT::UTILS::ExtractMyValues(*phinp,myphinp,lm);

  // fill all element arrays
  for (int i=0;i<nen_;++i)
  {
    for (int k = 0; k< numscal_; ++k)
    {
      // split for each transported scalar, insert into element arrays
      ephinp_[k](i,0) = myphinp[k+(i*numdofpernode_)];
      ephin_[k](i,0) = 0.0; //reset to zero; used in GetMaterialParams if not incremental
                            //-> used to calculate densn which is not required here
    }
    for (int k = 0; k< numscal_; ++k)
    {
      // the history vectors contains information of time step t_n
      ehist_[k](i,0) = myhist[k+(i*numdofpernode_)];
    }
  } // for i

  if ((scatratype == INPAR::SCATRA::scatratype_loma) and is_genalpha_)
  {
    // extract additional local values from global vector
      Teuchos::RCP<const Epetra_Vector> phiam = discretization.GetState("phiam");
    if (phiam==Teuchos::null) dserror("Cannot get state vector 'phiam'");
    std::vector<double> myphiam(lm.size());
    DRT::UTILS::ExtractMyValues(*phiam,myphiam,lm);

    // fill element array
    for (int i=0;i<nen_;++i)
    {
      for (int k = 0; k< numscal_; ++k)
      {
        // split for each transported scalar, insert into element arrays
        ephiam_[k](i,0) = myphiam[k+(i*numdofpernode_)];
      }
    } // for i
  }

  // get fine-scale values
  if (whichfssgd == INPAR::SCATRA::fssugrdiff_smagorinsky_small
      or turbmodel_ == INPAR::FLUID::multifractal_subgrid_scales)
  {
    // get fine scale scalar field
    Teuchos::RCP<const Epetra_Vector> gfsphinp = discretization.GetState("fsphinp");
    if (gfsphinp==Teuchos::null) dserror("Cannot get state vector 'fsphinp'");

    std::vector<double> myfsphinp(lm.size());
    DRT::UTILS::ExtractMyValues(*gfsphinp,myfsphinp,lm);

    for (int i=0;i<nen_;++i)
    {
      for (int k = 0; k< numscal_; ++k)
      {
        // split for each transported scalar, insert into element arrays
        fsphinp_[k](i,0) = myfsphinp[k+(i*numdofpernode_)];
      }
    }

    // get fine-scale velocity at nodes
    const Teuchos::RCP<Epetra_MultiVector> fsvelocity = params.get< Teuchos::RCP<Epetra_MultiVector> >("fine-scale velocity field");
    DRT::UTILS::ExtractMyNodeBasedValues(ele,efsvel_,fsvelocity,nsd_);
  }


  //----------------------------------------------------------------------
  // prepare mean values
  // ---------------------------------------------------------------------

  // the coordinates of the element layers in the channel
  // planecoords are named nodeplanes in turbulence_statistics_channel!
  Teuchos::RCP<std::vector<double> > planecoords  = params.get<Teuchos::RCP<std::vector<double> > >("planecoords_",Teuchos::null);
  if(planecoords==Teuchos::null)
    dserror("planecoords is null, but need channel_flow_of_height_2\n");

  //this will be the y-coordinate of a point in the element interior
  double center = 0.0;
  // get node coordinates of element
  for(int inode=0;inode<nen_;inode++)
    center+=xyze_(1,inode);

  center/=nen_;

  // working arrays for the quantities we want to compute
  double vol             = 0.0;

  double averaged_tauS   = 0.0;

  double mean_resS       = 0.0;
  double mean_resS_sq    = 0.0;

  double eps_visc        = 0.0;
  double eps_conv        = 0.0;
  double eps_smag        = 0.0;
  double eps_avm3        = 0.0;
  double eps_mfs         = 0.0;
  double eps_mfscross    = 0.0;
  double eps_mfsrey      = 0.0;
  double eps_supg        = 0.0;
  double eps_cross       = 0.0;
  double eps_rey         = 0.0;

  //----------------------------------------------------------------------
  // calculation of element volume both for tau at ele. cent. and int. pt.
  //----------------------------------------------------------------------
  // use one-point Gauss rule to do calculations at the element center
  DRT::UTILS::IntPointsAndWeights<nsd_> intpoints_tau(SCATRA::DisTypeToStabGaussRule<distype>::rule);

  // volume of the element (2D: element surface area; 1D: element length)
  // (Integration of f(x) = 1 gives exactly the volume/surface/length of element)
  vol = EvalShapeFuncAndDerivsAtIntPoint(intpoints_tau,0,ele->Id());

  //----------------------------------------------------------------------
  // get material parameters (evaluation at element center)
  //----------------------------------------------------------------------
  // remark: last parameter is dt required for MAT::m_myocard, which is not supported for
  //         turbulence modeling -> dt=-1.0 (since negative time should hopefully provoke
  //         a dserror() or a segmentation fault
  if (not mat_gp_ or not tau_gp_) GetMaterialParams(ele,scatratype,-1.0);

  //----------------------------------------------------------------------
  // calculation of subgrid diffusivity and stabilization parameter(s)
  // at element center
  //----------------------------------------------------------------------
  if (not tau_gp_)
  {
    // get velocity at element center
    velint_.Multiply(evelnp_,funct_);
    convelint_.Multiply(econvelnp_,funct_);

    // calculation of all-scale subgrid diffusivity (artificial or due to
    // constant-coefficient Smagorinsky model) at element center
    if (assgd or turbmodel_ == INPAR::FLUID::smagorinsky or turbmodel_ == INPAR::FLUID::dynamic_smagorinsky or turbmodel_ == INPAR::FLUID::dynamic_vreman)
      CalcSubgrDiff(dt,timefac,whichassgd,assgd,Cs,tpn,vol,0);

    // calculation of fine-scale artificial subgrid diffusivity at element center
    if (fssgd)
    {
      // we have to set a vector, which is however only required for one
      // special case (computation of fine-scale subgrid diffusivity for non-incremental
      // solver -> only artificial subgrid diffusivity) not considered here
      Epetra_SerialDenseVector  elevec1_epetra_subgrdiff_dummy;
      CalcFineScaleSubgrDiff(ele,elevec1_epetra_subgrdiff_dummy,whichfssgd,Cs,tpn,vol,0);
    }

    // calculation of stabilization parameter at element center
    CalTau(ele,diffus_[0],dt,timefac,vol,0,0.0,false);
  }

  //----------------------------------------------------------------------
  // prepare multifractal subgrid-scale modeling
  // calculation of model coefficients B (velocity) and D (scalar)
  // at element center
  //----------------------------------------------------------------------
  // coefficient B of fine-scale velocity
  LINALG::Matrix<nsd_,1> B_mfs(true);
  // coefficient D of fine-scale scalar
  double D_mfs = 0.0;
  if (turbmodel_ == INPAR::FLUID::multifractal_subgrid_scales)
  {
    if (not BD_gp)
    {
      // make sure to get material parameters at element center
      // hence, determine them if not yet available
      if (mat_gp_) GetMaterialParams(ele,scatratype,dt);
      // provide necessary velocities and gradients at element center
      convelint_.Multiply(econvelnp_,funct_);
      fsvelint_.Multiply(efsvel_,funct_);
      // calculate model coefficients
      CalcBAndDForMultifracSubgridScales(B_mfs,D_mfs,Csgs_sgvel,alpha,calc_N,N_vel,refvel,reflength,c_nu,nwl,nwl_scatra,Csgs_sgphi,c_diff,vol,0);
      // and clear them
      convelint_.Clear();
      fsvelint_.Clear();
    }
  }

  // get body force
  BodyForce(ele,time);

  //----------------------------------------------------------------------
  //                       INTEGRATION LOOP
  //----------------------------------------------------------------------
  // integrations points and weights
  DRT::UTILS::IntPointsAndWeights<nsd_> intpoints(SCATRA::DisTypeToOptGaussRule<distype>::rule);

  for (int iquad=0; iquad<intpoints.IP().nquad; ++iquad)
  {
    //---------------------------------------------------------------
    // evaluate shape functions and derivatives at integration point
    //---------------------------------------------------------------
    const double fac = EvalShapeFuncAndDerivsAtIntPoint(intpoints,iquad,ele->Id());

    //----------------------------------------------------------------------
    // get material parameters (evaluation at integration point)
    //----------------------------------------------------------------------
    if (mat_gp_) GetMaterialParams(ele,scatratype,-1.0);

    // scalar at integration point
    phi_[0]=funct_.Dot(ephinp_[0]);
    // gradient of current scalar value
    gradphi_.Multiply(derxy_,ephinp_[0]);

    // get velocity at integration point
    velint_.Multiply(evelnp_,funct_);
    convelint_.Multiply(econvelnp_,funct_);

    // convective term using current scalar value
    conv_phi_[0] = convelint_.Dot(gradphi_);

    // diffusive term using current scalar value for higher-order elements
    if (use2ndderiv_)
    {
      // diffusive part:  diffus * ( N,xx  +  N,yy +  N,zz )
      //GetLaplacianStrongForm(diff_, derxy2_);
      diff_.Clear();
      // compute N,xx  +  N,yy +  N,zz for each shape function at integration point
      for (int i=0; i<nen_; ++i)
      {
        for (int j = 0; j<nsd_; ++j)
        {
          diff_(i) += derxy2_(j,i);
        }
      }
      diff_.Scale(diffus_[0]);
      diff_phi_[0] = diff_.Dot(ephinp_[0]);
    }

    // reactive term using current scalar value
    if (is_reactive_) //is_reactive_ set in GetMaterial!
    {
      std::cout << "Warning: Reaction!" << std::endl;
      rea_phi_[0] = densnp_[0]*reacterm_[0]; //reacterm_ set in GetMaterial!
    }
    else rea_phi_[0] = 0.0;

    // get fine-scale velocity and its derivatives at integration point
    if (turbmodel_ == INPAR::FLUID::multifractal_subgrid_scales)
      fsvelint_.Multiply(efsvel_,funct_);
    else
      fsvelint_.Clear();

    // compute gradient of fine-scale part of scalar value
    if (fssgd)
      fsgradphi_.Multiply(derxy_,fsphinp_[0]);
    else
      fsgradphi_.Clear();

    // get history data (or acceleration)
    hist_[0] = funct_.Dot(ehist_[0]);

    // compute rhs containing bodyforce (divided by specific heat capacity) and,
    // for temperature equation, the time derivative of thermodynamic pressure,
    // if not constant, and for temperature equation of a reactive
    // equation system, the reaction-rate term
    rhs_[0] = bodyforce_[0].Dot(funct_)/shc_;
    rhs_[0] += thermpressdt_/shc_;
    if (reatemprhs_[0]!=0.0) std::cout << "Warning: Reaction!" << std::endl;
    rhs_[0] += densnp_[0]*reatemprhs_[0]; //reatemprhs_ set in GetMatarialParams!

    //--------------------------------------------------------------------
    // calculation of (fine-scale) subgrid diffusivity, subgrid-scale
    // velocity and stabilization parameter(s) at integration point
    //--------------------------------------------------------------------
    // ensure that subgrid-scale velocity and subgrid-scale convective part
    // are zero if not computed below
    sgvelint_.Clear();
    if (tau_gp_)
    {
      // calculation of all-scale subgrid diffusivity (artificial or due to
      // constant-coefficient Smagorinsky model) at integration point
      if (assgd or turbmodel_ == INPAR::FLUID::smagorinsky or turbmodel_ == INPAR::FLUID::dynamic_smagorinsky)
        CalcSubgrDiff(dt,timefac,whichassgd,assgd,Cs,tpn,vol,0);

      // calculation of fine-scale artificial subgrid diffusivity
      // at integration point
      if (fssgd)
      {
        // we have to set a vector, which is however only required for one
        // special case (computation of fine-scale subgrid diffusivity for non-incremental
        // solver -> only artificial subgrid diffusivity) not considered here
        Epetra_SerialDenseVector  elevec1_epetra_subgrdiff_dummy;
        CalcFineScaleSubgrDiff(ele,elevec1_epetra_subgrdiff_dummy,whichfssgd,Cs,tpn,vol,0);
      }

      // calculation of subgrid-scale velocity at integration point if required
      // this has to be done here (before CalTau() for scatra is called), otherwise
      // the stabilization parameter for scatra would be overwritten
      if (sgvel_)
      {
        // calculation of stabilization parameter related to fluid momentum
        // equation at integration point
        CalTau(ele,visc_,dt,timefac,vol,0,0.0,false);

        if (scatratype != INPAR::SCATRA::scatratype_levelset)
          CalcSubgrVelocity(ele,time,dt,timefac,0,scatratype);
        else dserror("CalcSubgrVelocityLevelSet not available anymore");
      }

      // calculation of stabilization parameter at integration point
      CalTau(ele,diffus_[0],dt,timefac,vol,0,0.0,false);
    }

    // prepare multifractal subgrid-scale modeling
    // calculation of model coefficients B (velocity) and D (scalar)
    // at element gauss point
    if (turbmodel_ == INPAR::FLUID::multifractal_subgrid_scales)
    {
      if (BD_gp)
      {
        // make sure to get material parameters at element center
        // hence, determine them if not yet available
        if (not mat_gp_)
        {
          // GetMaterialParams(ele,scatratype,dt); would overwrite materials
          // at the element center, hence BD_gp should always be combined with
          // mat_gp_
          dserror("evaluation of B and D at gauss-point should always be combined with evaluation of material at gauss-point!");
        }
        // calculate model coefficients
        CalcBAndDForMultifracSubgridScales(B_mfs,D_mfs,Csgs_sgvel,alpha,calc_N,N_vel,refvel,reflength,c_nu,nwl,nwl_scatra,Csgs_sgphi,c_diff,vol,0);
      }

      // calculate fine-scale velocity for multifractal subgrid-scale modeling
      for (int idim=0; idim<nsd_; idim++)
        mfsgvelint_(idim,0) = fsvelint_(idim,0) * B_mfs(idim,0);

      // calculate fine-scale scalar for multifractal subgrid-scale modeling
      mfssgphi_[0] = D_mfs * funct_.Dot(fsphinp_[0]);
    }
    else
    {
      mfsgvelint_.Clear();
      mfssgphi_[0] = 0.0;
    }

    // get residual of convection-diffusion equation and residual-based subgrid-scale scalar
    CalcResidualAndSubgrScalar(dt,timefac,0);

    // update material parameters based on inclusion of subgrid-scale
    // part of scalar (active only for mixture fraction,
    // Sutherland law and progress variable, for the time being)
    if (update_mat_)
    {
      if (turbmodel_ == INPAR::FLUID::multifractal_subgrid_scales)
        UpdateMaterialParams(ele,mfssgphi_[0],0);
      else
        UpdateMaterialParams(ele,sgphi_[0],0);
    }
    // this yields updated material parameters (dens, visc, diffus)
    // scatrares_ and sgphi_ are not updated, since they are used
    // for the contributions of the stabilizations, for which we
    // do not use updated material parameters

//    if (ele->Id()==100)
//    {
//      std::cout << "densnp_[0] " << densnp_[0] << std::endl;
//      std::cout << "sgphi_[0] " << sgphi_[0] << std::endl;
//      std::cout << "convelint_ " << convelint_ << std::endl;
//      std::cout << "gradphi_ " << gradphi_ << std::endl;
//      std::cout << "phi_[0] " << phi_[0] << std::endl;
//      std::cout << "sgphi_[0] " << sgphi_[0] << std::endl;
//      std::cout << "sgvelint_ " << sgvelint_ << std::endl;
//      std::cout << "mfssgphi_[0] " << mfssgphi_[0] << std::endl;
//      std::cout << "mfsgvelint_ " << mfsgvelint_ << std::endl;
//      std::cout << "sgdiff_[0] " << sgdiff_[0] << std::endl;
//      std::cout << "fsgradphi_ " << fsgradphi_ << std::endl;
//      std::cout << "diffus_[0] " << diffus_[0] << std::endl;
//      std::cout << "tau_[0] " << tau_[0] << std::endl;
//      std::cout << "scatrares_[0] " << scatrares_[0] << std::endl;
//    }

    //---------------------------------------------------------------
    // element average dissipation and production rates
    //---------------------------------------------------------------

    //---------------------------------------------------------------
    // residual-based subgrid-scale modeling terms
    //---------------------------------------------------------------

    // dissipation by supg-stabilization
    if (stabtype == INPAR::SCATRA::stabtype_SUPG)
    {
      eps_supg -= densnp_[0] * sgphi_[0] * convelint_.Dot(gradphi_) * fac; //sgphi_ is negative

    }
    else if (stabtype == INPAR::SCATRA::stabtype_no_stabilization)
    {
       // nothing to do
    }
    else
      dserror("Stabtype not yet supported!");

    // dissipation by cross-stress-stabilization
    // dissipation by reynolds-stress-stabilization
    if (sgvel_)
    {
      eps_cross -= densnp_[0] * phi_[0] * sgvelint_.Dot(gradphi_) * fac;
      eps_rey -= densnp_[0] * sgphi_[0] * sgvelint_.Dot(gradphi_) * fac;
    }

    //---------------------------------------------------------------
    // multifractal subgrid-scale modeling terms
    //---------------------------------------------------------------

    // dissipation multifractal subgrid-scales
    if (turbmodel_ == INPAR::FLUID::multifractal_subgrid_scales)
    {
      eps_mfs -= densnp_[0] * (mfssgphi_[0] * convelint_.Dot(gradphi_)
                              + phi_[0] * mfsgvelint_.Dot(gradphi_)
                              + mfssgphi_[0] * mfsgvelint_.Dot(gradphi_)) * fac;
      eps_mfscross -= densnp_[0] * (mfssgphi_[0] * convelint_.Dot(gradphi_)
                                   + phi_[0] * mfsgvelint_.Dot(gradphi_)) * fac;
      eps_mfsrey -= densnp_[0] * mfssgphi_[0] * mfsgvelint_.Dot(gradphi_) * fac;
    }

    //---------------------------------------------------------------
    // small-scale subgrid-viscosity subgrid-scale modeling terms
    //---------------------------------------------------------------

    // dissipation AVM3
    if (fssgd)
    {
      eps_avm3 += sgdiff_[0] * fsgradphi_.Dot(fsgradphi_) * fac;
    }

    //---------------------------------------------------------------
    // Smagorinsky model
    //---------------------------------------------------------------

    // dissipation (Smagorinsky)
    if (assgd or turbmodel_ == INPAR::FLUID::smagorinsky or turbmodel_ == INPAR::FLUID::dynamic_smagorinsky)
    {
      eps_smag += sgdiff_[0] * gradphi_.Dot(gradphi_) * fac;
    }

    //---------------------------------------------------------------
    // standard Galerkin terms
    //---------------------------------------------------------------

    // convective (Galerkin)
    eps_conv -= densnp_[0] * phi_[0] * convelint_.Dot(gradphi_) * fac;

    // dissipation (Galerkin)
    eps_visc += diffus_[0] * gradphi_.Dot(gradphi_) * fac;

    //---------------------------------------------------------------
    // element averages of tau_S and residual
    //---------------------------------------------------------------
    averaged_tauS += tau_[0] * fac;
    mean_resS    += scatrares_[0] * fac;
    mean_resS_sq += scatrares_[0] * scatrares_[0] * fac;

  } // end loop integration points

  mean_resS       /= vol;
  mean_resS_sq    /= vol;

  averaged_tauS   /= vol;

  eps_visc /= vol;
  eps_conv /= vol;
  eps_smag /= vol;
  eps_avm3 /= vol;
  eps_mfs /= vol;
  eps_mfscross /= vol;
  eps_mfsrey /= vol;
  eps_supg /= vol;
  eps_cross /= vol;
  eps_rey /= vol;


  Teuchos::RCP<std::vector<double> > incrvol           = params.get<Teuchos::RCP<std::vector<double> > >("incrvol"          );

  Teuchos::RCP<std::vector<double> > incr_eps_visc      = params.get<Teuchos::RCP<std::vector<double> > >("incr_scatra_eps_visc"    );
  Teuchos::RCP<std::vector<double> > incr_eps_conv      = params.get<Teuchos::RCP<std::vector<double> > >("incr_scatra_eps_conv"    );
  Teuchos::RCP<std::vector<double> > incr_eps_smag      = params.get<Teuchos::RCP<std::vector<double> > >("incr_scatra_eps_eddyvisc");
  Teuchos::RCP<std::vector<double> > incr_eps_avm3      = params.get<Teuchos::RCP<std::vector<double> > >("incr_scatra_eps_avm3"    );
  Teuchos::RCP<std::vector<double> > incr_eps_mfs       = params.get<Teuchos::RCP<std::vector<double> > >("incr_scatra_eps_mfs"     );
  Teuchos::RCP<std::vector<double> > incr_eps_mfscross  = params.get<Teuchos::RCP<std::vector<double> > >("incr_scatra_eps_mfscross");
  Teuchos::RCP<std::vector<double> > incr_eps_mfsrey    = params.get<Teuchos::RCP<std::vector<double> > >("incr_scatra_eps_mfsrey"  );
  Teuchos::RCP<std::vector<double> > incr_eps_supg      = params.get<Teuchos::RCP<std::vector<double> > >("incr_scatra_eps_supg"    );
  Teuchos::RCP<std::vector<double> > incr_eps_cross     = params.get<Teuchos::RCP<std::vector<double> > >("incr_scatra_eps_cross"   );
  Teuchos::RCP<std::vector<double> > incr_eps_rey       = params.get<Teuchos::RCP<std::vector<double> > >("incr_scatra_eps_rey"     );

  Teuchos::RCP<std::vector<double> > incrresS          = params.get<Teuchos::RCP<std::vector<double> > >("incrresS"         );
  Teuchos::RCP<std::vector<double> > incrresS_sq       = params.get<Teuchos::RCP<std::vector<double> > >("incrresS_sq"      );

  Teuchos::RCP<std::vector<double> > incrtauS          = params.get<Teuchos::RCP<std::vector<double> > >("incrtauS"         );

  bool found = false;

  int nlayer = 0;
  for (nlayer=0;nlayer<(int)(*planecoords).size()-1;)
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

  // collect layer volume
  (*incrvol      )[nlayer] += vol;

  // averages of stabilization parameters
  (*incrtauS     )[nlayer] += averaged_tauS;

  // averages residual
  (*incrresS         )[nlayer] += mean_resS      ;
  (*incrresS_sq      )[nlayer] += mean_resS_sq   ;

  // averages dissipation
  (*incr_eps_visc    )[nlayer] += eps_visc       ;
  (*incr_eps_conv    )[nlayer] += eps_conv       ;
  (*incr_eps_smag    )[nlayer] += eps_smag       ;
  (*incr_eps_avm3    )[nlayer] += eps_avm3       ;
  (*incr_eps_mfs     )[nlayer] += eps_mfs        ;
  (*incr_eps_mfscross)[nlayer] += eps_mfscross   ;
  (*incr_eps_mfsrey  )[nlayer] += eps_mfsrey     ;
  (*incr_eps_supg    )[nlayer] += eps_supg       ;
  (*incr_eps_cross   )[nlayer] += eps_cross      ;
  (*incr_eps_rey     )[nlayer] += eps_rey        ;

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
