/*----------------------------------------------------------------------*/
/*!
\file scatra_ele_calc_service_turbulence.cpp

\brief Internal implementation of ScaTra element

\level 2

<pre>
\maintainer Ursula Rasthofer
            rasthofer@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15236
</pre>
*/
/*----------------------------------------------------------------------*/

#include "scatra_ele_calc.H"
#include "scatra_ele.H"
#include "scatra_ele_parameter_std.H"
#include "scatra_ele_parameter_timint.H"
#include "scatra_ele_parameter_turbulence.H"
#include "../drt_inpar/inpar_fluid.H"

#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_lib/standardtypes_cpp.H"
#include "../drt_mat/newtonianfluid.H"
#include "../drt_mat/matlist.H"

#include "../drt_lib/drt_dserror.H"

#include "../drt_fluid/fluid_rotsym_periodicbc.H"


/*-----------------------------------------------------------------------------*
 | calculate filtered quantities for dynamic Smagorinsky model  rasthofer 08/12|
 *-----------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype,int probdim>
void DRT::ELEMENTS::ScaTraEleCalc<distype,probdim>::scatra_apply_box_filter(
    double&                                            dens_hat,
    double&                                            temp_hat,
    double&                                            dens_temp_hat,
    double&                                            phi2_hat,
    double&                                            phiexpression_hat,
    Teuchos::RCP<std::vector<double> >                 vel_hat,
    Teuchos::RCP<std::vector<double> >                 densvel_hat,
    Teuchos::RCP<std::vector<double> >                 densveltemp_hat,
    Teuchos::RCP<std::vector<double> >                 densstraintemp_hat,
    Teuchos::RCP<std::vector<double> >                 phi_hat,
    Teuchos::RCP<std::vector<std::vector<double> > >   alphaijsc_hat,
    double&                                            volume,
    const DRT::Element*                                ele,
    Teuchos::ParameterList&                            params
    )
{
  // do preparations first
  // ---------------------------------------------
  LINALG::Matrix<nsd_,nsd_> vderxy (true);
  double alpha2 = 0.0;
  // use one-point Gauss rule to do calculations at the element center
  volume = EvalShapeFuncAndDerivsAtEleCenter();

  // get material
  Teuchos::RCP<const MAT::Material> material = ele->Material();

  // get phi at integration point
  const double phinp = funct_.Dot(ephinp_[0]);

  // get temperature at integration point
  double tempnp = phinp;
  if (material->MaterialType() == INPAR::MAT::m_matlist)
  {
    const MAT::MatList* actmat = static_cast<const MAT::MatList*>(material.get());
    tempnp = funct_.Dot(ephinp_[actmat->NumMat()-1]);
  }

  // density at time n+1
  const double densnp = GetDensity(ele,material,params,tempnp);

  // get velocities (n+alpha_F/1,i) at integration point
  LINALG::Matrix<nsd_,1> convelint(true);
  convelint.Multiply(evelnp_,funct_);

  // compute rate of strain
  double rateofstrain = -1.0e30;
  rateofstrain = GetStrainRate(evelnp_);

  // gradient of scalar value
  LINALG::Matrix<nsd_,1> gradphi(true);
  gradphi.Multiply(derxy_,ephinp_[0]);

  // perform integrations, i.e., convolution
  // ---------------------------------------------
  for (int rr=0;rr<nsd_;++rr)
  {
    double tmp=convelint(rr)*volume;

    // add contribution to integral over velocities
    (*vel_hat)[rr] += tmp;

    // add contribution to integral over dens times velocity
    (*densvel_hat)[rr] += densnp*tmp;

    // add contribution to integral over dens times temperature times velocity
    (*densveltemp_hat)[rr] += densnp*phinp*tmp;
  }

  for (int rr=0;rr<nsd_;++rr)
  {
    double tmp=gradphi(rr)*volume;
    // add contribution to integral over dens times rate of strain times phi gradient
    (*densstraintemp_hat)[rr] += densnp*rateofstrain*tmp;
    (*phi_hat)[rr] = tmp;
    phi2_hat += tmp*gradphi(rr);
  }

  //calculate vreman part
  if (turbparams_->TurbModel() == INPAR::FLUID::dynamic_vreman)
  {
    if(nsd_==3)
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
    else
      dserror("Vreman model only for nsd_==3");
  }
  // add additional scalar quantities
  // i.e., filtered density, filtered density times scalar (i.e., temperature) and scalar
  dens_hat = densnp*volume;
  dens_temp_hat = densnp*phinp*volume;
  temp_hat = phinp*volume;

  return;
} //ScaTraEleCalc::scatra_apply_box_filter


/*-----------------------------------------------------------------------------*
 | get density at integration point                                 fang 02/15 |
 *-----------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype,int probdim>
double DRT::ELEMENTS::ScaTraEleCalc<distype,probdim>::GetDensity(
    const DRT::Element*                 ele,
    Teuchos::RCP<const MAT::Material>   material,
    Teuchos::ParameterList&             params,
    const double                        tempnp
    )
{
  // initialization
  double density(0.);

  if(material->MaterialType() == INPAR::MAT::m_scatra)
  {
    // access fluid discretization
    Teuchos::RCP<DRT::Discretization> fluiddis = Teuchos::null;
    fluiddis = DRT::Problem::Instance()->GetDis("fluid");
    // get corresponding fluid element (it has the same global ID as the scatra element)
    DRT::Element* fluidele = fluiddis->gElement(ele->Id());
    if(fluidele == NULL)
      dserror("Fluid element %i not on local processor", ele->Id());
    // get fluid material
    Teuchos::RCP<MAT::Material> fluidmat = fluidele->Material();
    if(fluidmat->MaterialType() != INPAR::MAT::m_fluid)
      dserror("Invalid fluid material for passive scalar transport in turbulent flow!");

    density = Teuchos::rcp_dynamic_cast<const MAT::NewtonianFluid>(fluidmat)->Density();
    if(density != 1.0)
      dserror("Check your diffusivity! Dynamic diffusivity required!");
  }

  else
    dserror("Invalid material type!");

  return density;
} // DRT::ELEMENTS::ScaTraEleCalc<distype,probdim>::GetDensity


/*----------------------------------------------------------------------------------*
 | calculate turbulent Prandtl number for dynamic Smagorinsky model  rasthofer 08/12|
 *----------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype,int probdim>
void DRT::ELEMENTS::ScaTraEleCalc<distype,probdim>::scatra_calc_smag_const_LkMk_and_MkMk(
        Teuchos::RCP<Epetra_MultiVector>&  col_filtered_vel,
        Teuchos::RCP<Epetra_MultiVector>&  col_filtered_dens_vel,
        Teuchos::RCP<Epetra_MultiVector>&  col_filtered_dens_vel_temp,
        Teuchos::RCP<Epetra_MultiVector>&  col_filtered_dens_rateofstrain_temp,
        Teuchos::RCP<Epetra_Vector>&       col_filtered_temp,
        Teuchos::RCP<Epetra_Vector>&       col_filtered_dens,
        Teuchos::RCP<Epetra_Vector>&       col_filtered_dens_temp,
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
  DRT::UTILS::IntPointsAndWeights<nsd_ele_> intpoints(SCATRA::DisTypeToStabGaussRule<distype>::rule);

  EvalShapeFuncAndDerivsAtIntPoint(intpoints,0);

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
}//ScaTraEleCalc::scatra_calc_smag_const_LkMk_and_MkMk


/*----------------------------------------------------------------------------------*
 | calculate vreman constant                                             krank 08/13|
 *----------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype,int probdim>
void DRT::ELEMENTS::ScaTraEleCalc<distype,probdim>::scatra_calc_vreman_dt(
  Teuchos::RCP<Epetra_MultiVector>& col_filtered_phi,
  Teuchos::RCP<Epetra_Vector>& col_filtered_phi2              ,
  Teuchos::RCP<Epetra_Vector>&   col_filtered_phiexpression         ,
  Teuchos::RCP<Epetra_MultiVector>& col_filtered_alphaijsc,
  double& dt_numerator,
  double& dt_denominator,
  const DRT::Element*       ele
  )
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
  DRT::UTILS::IntPointsAndWeights<nsd_ele_> intpoints(SCATRA::DisTypeToStabGaussRule<distype>::rule);
  double volume = EvalShapeFuncAndDerivsAtIntPoint(intpoints,0);

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
template <DRT::Element::DiscretizationType distype,int probdim>
void DRT::ELEMENTS::ScaTraEleCalc<distype,probdim>::GetMeanPrtOfHomogenousDirection(
  Teuchos::ParameterList&    turbmodelparams,
  int&                       nlayer
)
{
  // NOTE: we calculate the inverse of the turbulent Prandtl number here (i.e., (Cs*h)^2 / Pr_t)

  if(nsd_ != 3)
    dserror("turbulence and 3D flow !");

  if (turbmodelparams.get<std::string>("HOMDIR","not_specified") !=  "not_specified")
  {
    Teuchos::RCP<std::vector<double> > averaged_LkMk = turbmodelparams.get<Teuchos::RCP<std::vector<double> > >("averaged_LkMk_");
    Teuchos::RCP<std::vector<double> > averaged_MkMk = turbmodelparams.get<Teuchos::RCP<std::vector<double> > >("averaged_MkMk_");

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
      Teuchos::RCP<std::vector<double> > planecoords = turbmodelparams.get<Teuchos::RCP<std::vector<double> > >("planecoords_");
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
      Teuchos::RCP<std::vector<double> > dir1coords = turbmodelparams.get<Teuchos::RCP<std::vector<double> > >("dir1coords_");
      Teuchos::RCP<std::vector<double> > dir2coords = turbmodelparams.get<Teuchos::RCP<std::vector<double> > >("dir2coords_");
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
      tpn_ = 0.0;
    }
    else
      tpn_ = (*averaged_LkMk)[nlayer]/(*averaged_MkMk)[nlayer] ;
    // clipping to get algorithm stable
    if (tpn_<0.0)
    {
      tpn_=0.0;
    }

    }

  return;
}//ScaTraEleCalc::GetMeanOfHomogenousDirection


/*----------------------------------------------------------------------*
  |  calculate all-scale art. subgrid diffusivity (private)     vg 10/09 |
  *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype,int probdim>
void DRT::ELEMENTS::ScaTraEleCalc<distype,probdim>::CalcSubgrDiff(
  double&                               visc,
  const double                          vol,
  const int                             k,
  const double                          densnp
  )
{
  // get number of dimensions
  const double dim = (double) nsd_;

  // get characteristic element length as cubic root of element volume
  // (2D: square root of element area, 1D: element length)
  const double h = std::pow(vol,(1.0/dim));

  // subgrid-scale diffusivity to be computed and added to diffus
  double sgdiff(0.0);

  // all-scale subgrid diffusivity due to Smagorinsky model divided by
  // turbulent Prandtl number
  if (turbparams_->TurbModel() == INPAR::FLUID::smagorinsky)
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
    sgdiff = densnp * turbparams_->Cs() * turbparams_->Cs() * h * h * rateofstrain / tpn_;

    // add subgrid viscosity to physical viscosity for computation
    // of subgrid-scale velocity when turbulence model is applied
    if (scatrapara_->RBSubGrVel()) visc += sgdiff*tpn_;
  }
  else if (turbparams_->TurbModel() == INPAR::FLUID::dynamic_smagorinsky)
  {
    // compute (all-scale) rate of strain
    double rateofstrain = -1.0e30;
    rateofstrain = GetStrainRate(econvelnp_);

    // subgrid diffusivity = subgrid viscosity / turbulent Prandtl number
    // remark: for dynamic estimation, tpn corresponds to (Cs*h)^2 / Pr_t
    sgdiff = densnp * rateofstrain * tpn_;
  }
  else if (turbparams_->TurbModel() == INPAR::FLUID::dynamic_vreman)
  {
    if(nsd_==3)
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
        sgviscwocv=sqrt(bbeta / alphavreman);


      //remark: Cs corresponds to Dt, calculated in the vreman class
      //        The vreman constant Cv is not required here, since it cancelles out with the
      //        vreman constant omitted during the calculation of D_t
      if (turbparams_->Cs()<=1.0E-12) sgdiff=0.0;
      else
        sgdiff = densnp * sgviscwocv / turbparams_->Cs();
    }
    else
      dserror("Vreman model only for nsd_==3");

  }

  // update diffusivity
  diffmanager_->SetIsotropicSubGridDiff(sgdiff,k);

  return;
} //ScaTraEleCalc::CalcSubgrDiff


/*----------------------------------------------------------------------*
  |  calculate fine-scale art. subgrid diffusivity (private)    vg 10/09 |
  *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype,int probdim>
void DRT::ELEMENTS::ScaTraEleCalc<distype,probdim>::CalcFineScaleSubgrDiff(
  double&                               sgdiff,
  Epetra_SerialDenseVector&             subgrdiff,
  DRT::Element*                         ele,
  const double                          vol,
  const int                             k,
  const double                          densnp,
  const double                          diffus,
  const LINALG::Matrix<nsd_,1>          convelint
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
  if (not scatraparatimint_->IsIncremental())
  {
    // get element-type constant
    const double mk = SCATRA::MK<distype>();

    // velocity norm
    const double vel_norm = convelint.Norm2();

    // parameter relating convective and diffusive forces + respective switch
    const double epe = mk * densnp * vel_norm * h / diffus;
    const double xi = std::max(epe,1.0);

    // compute artificial subgrid diffusivity
    sgdiff = (DSQR(h)*mk*DSQR(vel_norm)*DSQR(densnp))/(2.0*diffus*xi);

    // compute entries of (fine-scale) subgrid-diffusivity-scaling vector
    for (int vi=0; vi<nen_; ++vi)
    {
      subgrdiff(vi) = sgdiff/ele->Nodes()[vi]->NumElement();
    }
  }
  //----------------------------------------------------------------------
  // computation of fine-scale subgrid diffusivity for incremental solver
  // -> only all-scale Smagorinsky model
  //----------------------------------------------------------------------
  else
  {
    if (turbparams_->WhichFssgd() == INPAR::SCATRA::fssugrdiff_smagorinsky_all)
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
      sgdiff = densnp * turbparams_->Cs() * turbparams_->Cs() * h * h * rateofstrain / tpn_;
    }
    else if (turbparams_->WhichFssgd() == INPAR::SCATRA::fssugrdiff_smagorinsky_small)
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
      sgdiff = densnp * turbparams_->Cs() * turbparams_->Cs() * h * h * fsrateofstrain / tpn_;
    }
  }

  return;
} //ScaTraEleCalcCalc::FineScaleSubgrDiff


/*----------------------------------------------------------------------*
 | calculation of coefficients B and D for multifractal subgrid-scales  |
 |                                                      rasthofer 12/11 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype,int probdim>
void DRT::ELEMENTS::ScaTraEleCalc<distype,probdim>::CalcBAndDForMultifracSubgridScales(
    LINALG::Matrix<nsd_,1>&                     B_mfs, ///< coefficient for fine-scale velocity (will be filled)
    double &                                    D_mfs, ///< coefficient for fine-scale scalar (will be filled)
    const double                                vol, ///< volume of element
    const int                                   k,
    const double                                densnp,
    const double                                diffus,
    const double                                visc,
    const LINALG::Matrix<nsd_,1>                convelint,
    const LINALG::Matrix<nsd_,1>                fsvelint
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
  double Csgs_vel_nw = turbparams_->Csgs_SgVel();

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
  const double vel_norm = convelint.Norm2();
  // also for fine-scale velocity
  // fsvelint_.Multiply(efsvel_,funct_);
  const double fsvel_norm = fsvelint.Norm2();

  // do we have a fixed parameter N
  if (not turbparams_->Calc_N())
  {
    // yes, store value
    for (int rr=1;rr<3;rr++)
      Nvel[rr] = turbparams_->N_Vel();
  }
  else //no, so we calculate N from Re
  {
    // calculate characteristic element length
    hk = CalcRefLength(vol,convelint);

    // warning: k=0, this first scalar is taken!
    // multifractal subgrid-scale model is for passive and active
    // scalar transport
    // therefore, we need the density of the fluid here
    switch (turbparams_->RefVel())
    {
      case INPAR::FLUID::resolved:
      {
        Re_ele = vel_norm * hk * densnp / visc;
        break;
      }
      case INPAR::FLUID::fine_scale:
      {
        Re_ele = fsvel_norm * hk * densnp / visc;
        break;
      }
      case INPAR::FLUID::strainrate:
      {
        strainnorm = GetStrainRate(econvelnp_);
        strainnorm /= sqrt(2.0); //cf. Burton & Dahm 2008
        Re_ele = strainnorm * hk * hk * densnp / visc;
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
    scale_ratio = turbparams_->C_Nu() * pow(Re_ele,0.75);
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
  if (turbparams_->Nwl())
  {
    // not yet calculated, estimate norm of strain rate
    if ((not turbparams_->Calc_N()) or (turbparams_->RefVel() != INPAR::FLUID::strainrate))
    {
      strainnorm = GetStrainRate(econvelnp_);
      strainnorm /= sqrt(2.0); //cf. Burton & Dahm 2008
    }
    // and reference length
    if (not turbparams_->Calc_N())
      hk = CalcRefLength(vol,convelint);

    // get Re from strain rate
    double Re_ele_str = strainnorm * hk * hk * densnp / visc;
    if (Re_ele_str < 0.0)
      dserror("Something went wrong!");
    // ensure positive values
    if (Re_ele_str < 1.0)
       Re_ele_str = 1.0;

    // calculate corrected Csgs
    //           -3/16
    //  *(1 - (Re)   )
    //
    Csgs_vel_nw *= (1.0-pow(Re_ele_str,-3.0/16.0));

    // store Cai for application to scalar field
    Cai_phi = (1.0-pow(Re_ele_str,-3.0/16.0));
  }

  // STEP 2: calculate B

  //       1                           1
  //       2   |       1              |2
  //  kappa  = | -------------------- |
  //           |  1 - alpha ^ (-4/3)  |
  //
  double kappa = 1.0/(1.0-pow(turbparams_->Alpha(),-4.0/3.0));

  //                  1                                     1
  //                  2                  |                 |2
  //  B = Csgs * kappa  * 2 ^ (-2*N/3) * | 2 ^ (4*N/3) - 1 |
  //                                     |                 |
  //
  for (int dim=0; dim<nsd_; dim++)
  {
    B_mfs(dim,0) = Csgs_vel_nw * sqrt(kappa) * pow(2.0,-2.0*Nvel[dim]/3.0) * sqrt((pow(2.0,4.0*Nvel[dim]/3.0)-1.0));
//    if (eid_ == 100)
//     std::cout << "B  " << std::setprecision (10) << B_mfs(dim,0) << std::endl;
  }

  //----------------------------------------------------------------
  // calculation of D for fine-scale scalar
  //----------------------------------------------------------------

  // STEP 1: determine N

  // calculate Prandtl number or Schmidt number (passive scalar)
  const double Pr = visc/diffus;

  // since there are differences in the physical behavior between low and high
  // Prandtl/Schmidt number regime, we define a limit
  // to distinguish between the low and high Prandtl/Schmidt number regime
  // note: there is no clear definition of the ranges
  const double Pr_limit = 2.0;

  // allocate vector for parameter N
  double Nphi = 0.0;
  // ratio of diffusive scale to element length
  double scale_ratio_phi = 0.0;

  if (turbparams_->Calc_N())
  {
    //
    //   Delta
    //  ---------  ~ Re^(3/4)*Pr^(p)
    //  lambda_diff
    //
    // Pr <= 1: p=3/4
    // Pr >> 1: p=1/2
    double p = 0.75;
    if (Pr>Pr_limit) p =0.5;

    scale_ratio_phi = turbparams_->C_Diff() * pow(Re_ele,0.75) * pow(Pr,p);
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
  double kappa_phi = 1.0/(1.0-pow(turbparams_->Alpha(),-gamma));

  //                                                             1
  //       Phi    Phi                       |                   |2
  //  D = Csgs * kappa * 2 ^ (-gamma*N/2) * | 2 ^ (gamma*N) - 1 |
  //                                        |                   |
  //
  if (not two_ranges) // usual case
    D_mfs = turbparams_->Csgs_SgPhi() *sqrt(kappa_phi) * pow(2.0,-gamma*Nphi/2.0) * sqrt((pow(2.0,gamma*Nphi)-1.0));
  else
  {
    double gamma1 = 4.0/3.0;
    double gamma2 = 2.0;
    kappa_phi = 1.0/(1.0-pow(turbparams_->Alpha(),-gamma1));
    D_mfs = turbparams_->Csgs_SgPhi() * sqrt(kappa_phi) * pow(2.0,-gamma2*Nphi/2.0) * sqrt((pow(2.0,gamma1*Nvel[0])-1.0)+2.0/3.0*pow((M_PI/hk),2.0/3.0)*(pow(2.0,gamma2*Nphi)-pow(2.0,gamma2*Nvel[0])));
  }

  // apply near-wall limit if required
  if (turbparams_->Nwl_ScaTra() and turbparams_->Nwl())
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
template <DRT::Element::DiscretizationType distype,int probdim>
double DRT::ELEMENTS::ScaTraEleCalc<distype,probdim>::CalcRefLength(
        const double                       vol,
        const LINALG::Matrix<nsd_,1>       convelint
  )
{
  // calculate characteristic element length
  double hk = 1.0e+10;
  // cf. stabilization parameters
  switch (turbparams_->RefLength())
  {
    case INPAR::FLUID::streamlength:
    {
      // a) streamlength due to Tezduyar et al. (1992)
      // get norm of velocity
      const double vel_norm = convelint.Norm2();
      // normed velocity vector
      LINALG::Matrix<nsd_,1> velino(true);
      if (vel_norm>=1e-6) velino.Update(1.0/vel_norm,convelint);
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
      // c) cubic element length
      hk = std::pow(vol,(1.0/(double (nsd_))));

      break;
    }
    case INPAR::FLUID::metric_tensor:
    {
      if (nsd_ != 3) dserror("Turbulence is 3d!");
      /*          +-           -+   +-           -+   +-           -+
                |             |   |             |   |             |
                |  dr    dr   |   |  ds    ds   |   |  dt    dt   |
          G   = |  --- * ---  | + |  --- * ---  | + |  --- * ---  |
           ij   |  dx    dx   |   |  dx    dx   |   |  dx    dx   |
                |    i     j  |   |    i     j  |   |    i     j  |
                +-           -+   +-           -+   +-           -+
      */
      LINALG::Matrix<nsd_,nsd_> G(true);

      for (int nn=0;nn<nsd_;++nn)
      {
        for (int rr=0;rr<nsd_;++rr)
        {
          G(nn,rr) = xij_(nn,0)*xij_(rr,0);
          for (int mm=1;mm<nsd_;++mm)
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
     double normG = 0.0;
     for (int nn=0;nn<nsd_;++nn)
     {
       for (int rr=0;rr<nsd_;++rr)
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
    if (nsd_ != 3) dserror("Turbulence is 3d!");
    LINALG::Matrix<nsd_,1> normed_velgrad;

    for (int rr=0;rr<nsd_;++rr)
    {
      double val = 0.0;
      for (int idim = 0; idim < nsd_; idim ++)
         val += convderxy(idim,rr)*convderxy(idim,rr);

      normed_velgrad(rr) = std::sqrt(val);

      //normed_velgrad(rr)=sqrt(convderxy(0,rr)*convderxy(0,rr)
      //                        +
      //                        convderxy(1,rr)*convderxy(1,rr)
      //                        +
      //                        convderxy(2,rr)*convderxy(2,rr));
    }
    double norm=normed_velgrad.Norm2();

    // normed gradient
    if (norm>1e-6)
    {
      for (int rr=0;rr<nsd_;++rr)
      {
        normed_velgrad(rr)/=norm;
      }
    }
    else
    {
      normed_velgrad(0) = 1.;
      for (int rr=1;rr<nsd_;++rr)
      {
        normed_velgrad(rr)=0.0;
      }
    }

    // get length in this direction
    double val = 0.0;
    for (int rr=0;rr<nen_;++rr) /* loop element nodes */
    {
      double loc = 0.0;
      for (int idim = 0; idim < nsd_; idim ++)
        loc += normed_velgrad(idim)*derxy_(idim,rr);

      val += fabs(loc);

      //val += fabs( normed_velgrad(0)*derxy_(0,rr)
      //            +normed_velgrad(1)*derxy_(1,rr)
      //            +normed_velgrad(2)*derxy_(2,rr));
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
template <DRT::Element::DiscretizationType distype,int probdim>
void DRT::ELEMENTS::ScaTraEleCalc<distype,probdim>::StoreModelParametersForOutput(
  const DRT::Element*                   ele,
  const bool                            isowned,
  Teuchos::ParameterList&               turbulencelist,
  const int                             nlayer)
{
  if (isowned)
  {
    if (turbparams_->TurbModel() == INPAR::FLUID::dynamic_smagorinsky)
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
           // calculate Prt from (Cs*h)^2/Prt
          if (tpn_>1.0E-16)
          {
            // get dynamically estimated Smagorinsky constant from fluid element, i.e., (Cs*h)^2
            double Cs_delta_sq = (*(turbulencelist.get<Teuchos::RCP<std::vector<double> > >("global_Cs_delta_sq_sum")))[nlayer];
            // since Cs_delta_sq contains the sum over all elements of this layer,
            // we have to divide by the number of elements of this layer
            int numele_layer = turbulencelist.get<int>("numele_layer");
            (*(turbulencelist.get<Teuchos::RCP<std::vector<double> > >("local_Prt_sum")))[nlayer]+=(Cs_delta_sq/numele_layer)/tpn_;

          }
          else
            (*(turbulencelist.get<Teuchos::RCP<std::vector<double> > >("local_Prt_sum")))         [nlayer]+=0.0;

          // set (Cs*h)^2/Prt and diffeff for output
          (*(turbulencelist.get<Teuchos::RCP<std::vector<double> > >("local_Cs_delta_sq_Prt_sum")))[nlayer]+=tpn_;
          if (numscal_>1) dserror("One scalar assumed for dynamic Smagorinsky model!");

          // calculation of effective diffusion coefficient
          const double vol=EvalShapeFuncAndDerivsAtEleCenter();

          // get material  at element center
          // density at t_(n)
          std::vector<double> densn(numscal_,1.0);
          // density at t_(n+1) or t_(n+alpha_F)
          std::vector<double> densnp(numscal_,1.0);
          // density at t_(n+alpha_M)
          std::vector<double> densam(numscal_,1.0);

          // fluid viscosity
          double visc(0.0);

          SetInternalVariablesForMatAndRHS();

          GetMaterialParams(ele,densn,densnp,densam,visc);

          CalcSubgrDiff(visc,vol,0,densnp[0]);

          (*(turbulencelist.get<Teuchos::RCP<std::vector<double> > >("local_diffeff_sum")))    [nlayer]+=diffmanager_->GetIsotropicDiff(0);
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
template <DRT::Element::DiscretizationType distype,int probdim>
void DRT::ELEMENTS::ScaTraEleCalc<distype,probdim>::CalcDissipation(
    Teuchos::ParameterList&               params,          //!< parameter list
    DRT::Element*                         ele,             //!< pointer to element
    DRT::Discretization&                  discretization,  //!< scatra discretization
    DRT::Element::LocationArray&          la               //!< location array
    )
{
  // do some checks first
  if (numscal_!=1 or numdofpernode_!=1)
    dserror("CalcDissipation only for one scalar field!");

  //----------------------------------------------------------------------
  // preliminary set-up of parameters
  // ---------------------------------------------------------------------

    // as the scalar field is constant in the turbulent inflow section
    // we do not need any turbulence model
    if (turbparams_->TurbInflow()) dserror("CalcDissipation in combination with inflow generation not supported!");
//    if (params.get<bool>("turbulent inflow",false))
//    {
//      if (SCATRA::InflowElement(ele))
//        turbmodel_ = INPAR::FLUID::no_model;
//    }

  // set time integration
  if (scatraparatimint_->IsStationary()) dserror("Turbulence is instationary!");

  // set turbulent Prandt number to value given in parameterlist
  tpn_ = turbparams_->TPN();

  // if we have a dynamic model,we overwrite this value by a local element-based one here
  if (turbparams_->TurbModel() == INPAR::FLUID::dynamic_smagorinsky)
  {
    Teuchos::ParameterList& turbulencelist = params.sublist("TURBULENCE MODEL");
    // remark: for dynamic estimation, this returns (Cs*h)^2 / Pr_t
    Teuchos::RCP<Epetra_Vector> ele_prt = turbulencelist.get<Teuchos::RCP<Epetra_Vector> >("col_ele_Prt");
    const int id = ele->LID();
    tpn_ = (*ele_prt)[id];

    int dummy=0;
    // when no averaging was done, we just keep the calculated (clipped) value
    if (turbparams_->CsAv())
      GetMeanPrtOfHomogenousDirection(params.sublist("TURBULENCE MODEL"),dummy);
  }

  //----------------------------------------------------------------------
  // get all nodal values
  // ---------------------------------------------------------------------

  // get number of dofset associated with velocity related dofs
  const int ndsvel = params.get<int>("ndsvel");

  // get velocity values at nodes
  const Teuchos::RCP<const Epetra_Vector> convel = discretization.GetState(ndsvel,"convective velocity field");

  // safety check
  if(convel == Teuchos::null)
    dserror("Cannot get state vector convective velocity");

  // determine number of velocity related dofs per node
  const int numveldofpernode = la[ndsvel].lm_.size()/nen_;

  // construct location vector for velocity related dofs
  std::vector<int> lmvel(nsd_*nen_,-1);
  for (int inode=0; inode<nen_; ++inode)
    for (int idim=0; idim<nsd_; ++idim)
      lmvel[inode*nsd_+idim] = la[ndsvel].lm_[inode*numveldofpernode+idim];

  // extract local values of convective velocity field from global state vector
  DRT::UTILS::ExtractMyValues<LINALG::Matrix<nsd_,nen_> >(*convel,econvelnp_,lmvel);

  // rotate the vector field in the case of rotationally symmetric boundary conditions
  rotsymmpbc_->RotateMyValuesIfNecessary(econvelnp_);

  // set econvelnp_ equal to evelnp_ since ale is not supported
  evelnp_ = econvelnp_;

  // get data required for subgrid-scale velocity: acceleration and pressure
  if (scatrapara_->RBSubGrVel())
  {
    // get acceleration values at nodes
    const Teuchos::RCP<const Epetra_Vector> acc = discretization.GetState(ndsvel,"acceleration field");
    if(acc == Teuchos::null)
      dserror("Cannot get state vector acceleration field");

    // extract local values of acceleration field from global state vector
    DRT::UTILS::ExtractMyValues<LINALG::Matrix<nsd_,nen_> >(*acc,eaccnp_,lmvel);

    // rotate the vector field in the case of rotationally symmetric boundary conditions
    rotsymmpbc_->RotateMyValuesIfNecessary(eaccnp_);

    // construct location vector for pressure dofs
    std::vector<int> lmpre(nen_,-1);
    for (int inode=0; inode<nen_; ++inode)
      lmpre[inode] = la[ndsvel].lm_[inode*numveldofpernode+nsd_];

    // extract local values of pressure field from global state vector
    DRT::UTILS::ExtractMyValues<LINALG::Matrix<nen_,1> >(*convel,eprenp_,lmpre);
  }

  // extract local values from the global vectors
  Teuchos::RCP<const Epetra_Vector> hist = discretization.GetState("hist");
  Teuchos::RCP<const Epetra_Vector> phinp = discretization.GetState("phinp");
  if (hist==Teuchos::null || phinp==Teuchos::null)
    dserror("Cannot get state vector 'hist' and/or 'phinp'");
  DRT::UTILS::ExtractMyValues<LINALG::Matrix<nen_,1> >(*hist,ehist_,la[0].lm_);
  DRT::UTILS::ExtractMyValues<LINALG::Matrix<nen_,1> >(*phinp,ephinp_,la[0].lm_);

  // reset to zero; used in GetMaterialParams if not incremental -> used to calculate densn which is not required here
  for (int k=0; k<numdofpernode_; ++k)
    ephin_[k].Clear();

  // get fine-scale values
  if (turbparams_->WhichFssgd() == INPAR::SCATRA::fssugrdiff_smagorinsky_small
      or turbparams_->TurbModel() == INPAR::FLUID::multifractal_subgrid_scales)
  {
    // get fine scale scalar field
    Teuchos::RCP<const Epetra_Vector> gfsphinp = discretization.GetState("fsphinp");
    if (gfsphinp==Teuchos::null) dserror("Cannot get state vector 'fsphinp'");

    DRT::UTILS::ExtractMyValues<LINALG::Matrix<nen_,1> >(*gfsphinp,fsphinp_,la[0].lm_);

    // get fine-scale velocity at nodes
    const Teuchos::RCP<const Epetra_Vector> fsvelocity = discretization.GetState(ndsvel,"fine-scale velocity field");
    if(fsvelocity == Teuchos::null)
      dserror("Cannot get fine-scale velocity field from scatra discretization!");

    // extract local values of fine-scale velocity field from global state vector
    DRT::UTILS::ExtractMyValues<LINALG::Matrix<nsd_,nen_> >(*fsvelocity,efsvel_,lmvel);

    // rotate the vector field in the case of rotationally symmetric boundary conditions
    rotsymmpbc_->RotateMyValuesIfNecessary(efsvel_);
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

  // volume of the element (2D: element surface area; 1D: element length)
  // (Integration of f(x) = 1 gives exactly the volume/surface/length of element)
  vol = EvalShapeFuncAndDerivsAtEleCenter();

  //----------------------------------------------------------------------
  // get material parameters (evaluation at element center)
  //----------------------------------------------------------------------

  // density at t_(n)
  std::vector<double> densn(numscal_,1.0);
  // density at t_(n+1) or t_(n+alpha_F)
  std::vector<double> densnp(numscal_,1.0);
  // density at t_(n+alpha_M)
  std::vector<double> densam(numscal_,1.0);

  // fluid viscosity
  double visc(0.0);

  // material parameter at the element center are also necessary
  // even if the stabilization parameter is evaluated at the element center
  if (not scatrapara_->MatGP() or scatrapara_->TauGP())
  {
    SetInternalVariablesForMatAndRHS();

    GetMaterialParams(ele,densn,densnp,densam,visc);
  }

  //----------------------------------------------------------------------
  // calculation of subgrid diffusivity and stabilization parameter(s)
  // at element center
  //----------------------------------------------------------------------

  // the stabilization parameters (one per transported scalar)
  std::vector<double> tau(numscal_,0.0);
  // subgrid-scale diffusion coefficient
  double sgdiff(0.0);

  if (not scatrapara_->TauGP())
  {
    // calculation of all-scale subgrid diffusivity (by, e.g.,
    // Smagorinsky model) at element center
    if (turbparams_->TurbModel() == INPAR::FLUID::smagorinsky
        or turbparams_->TurbModel() == INPAR::FLUID::dynamic_smagorinsky
         or turbparams_->TurbModel() == INPAR::FLUID::dynamic_vreman)
    {
      CalcSubgrDiff(visc,vol,0,densnp[0]);
    }

    // calculation of fine-scale artificial subgrid diffusivity at element center
    // we have to set a vector, which is however only required for one
    // special case (computation of fine-scale subgrid diffusivity for non-incremental
    // solver -> only artificial subgrid diffusivity) not considered here
    Epetra_SerialDenseVector  elevec1_epetra_subgrdiff_dummy;
    if (turbparams_->FSSGD()) CalcFineScaleSubgrDiff(sgdiff,elevec1_epetra_subgrdiff_dummy,ele,vol,0,densnp[0],diffmanager_->GetIsotropicDiff(0),scatravarmanager_->ConVel());

    // calculation of stabilization parameter at element center
    CalcTau(tau[0],diffmanager_->GetIsotropicDiff(0),reamanager_->GetReaCoeff(0),densnp[0],scatravarmanager_->ConVel(),vol);
  }


  // prepare multifractal subgrid-scale modeling
  // calculation of model coefficients B (velocity) and D (scalar)
  // at element center
  // coefficient B of fine-scale velocity
  LINALG::Matrix<nsd_,1> B_mfs(true);
  // coefficient D of fine-scale scalar
  double D_mfs = 0.0;
  if (turbparams_->TurbModel() == INPAR::FLUID::multifractal_subgrid_scales)
  {
    if (not turbparams_->BD_Gp())
    {
      // make sure to get material parameters at element center
      // hence, determine them if not yet available
      if (scatrapara_->MatGP())
      {
        SetInternalVariablesForMatAndRHS();

        GetMaterialParams(ele,densn,densnp,densam,visc);
      }
      // provide necessary velocities and gradients at element center
      // get velocity at element center
      LINALG::Matrix<nsd_,1> fsvelint(true);
      fsvelint.Multiply(efsvel_,funct_);

      // calculate model coefficients
      for (int k = 0;k<numscal_;++k) // loop of each transported scalar
        CalcBAndDForMultifracSubgridScales(B_mfs,D_mfs,vol,k,densnp[0],diffmanager_->GetIsotropicDiff(k),visc,scatravarmanager_->ConVel(),fsvelint);
    }
  }

  // get body force
  BodyForce(ele);

  //----------------------------------------------------------------------
  //                       INTEGRATION LOOP
  //----------------------------------------------------------------------
  // integration points and weights
  DRT::UTILS::IntPointsAndWeights<nsd_ele_> intpoints(SCATRA::DisTypeToOptGaussRule<distype>::rule);

  for (int iquad=0; iquad<intpoints.IP().nquad; ++iquad)
  {
    //---------------------------------------------------------------
    // evaluate shape functions and derivatives at integration point
    //---------------------------------------------------------------
    const double fac = EvalShapeFuncAndDerivsAtIntPoint(intpoints,iquad);

    //----------------------------------------------------------------------
    // get material parameters (evaluation at integration point)
    //----------------------------------------------------------------------

    SetInternalVariablesForMatAndRHS();

    if (scatrapara_->MatGP())
      GetMaterialParams(ele,densn,densnp,densam,visc);

    // get velocity at integration point
    const LINALG::Matrix<nsd_,1>& convelint = scatravarmanager_->ConVel();

    // scalar at integration point at time step n+1
    const double& phinp = scatravarmanager_->Phinp(0);

    // gradient of current scalar value at integration point
    LINALG::Matrix<nsd_,1> gradphi(true);
    gradphi.Multiply(derxy_,ephinp_[0]);

    // reactive part of the form: (reaction coefficient)*phi
    const double rea_phi = densnp[0]*phinp*reamanager_->GetReaCoeff(0);

    // get fine-scale velocity and its derivatives at integration point
    LINALG::Matrix<nsd_,1> fsvelint(true);
    if (turbparams_->TurbModel() == INPAR::FLUID::multifractal_subgrid_scales)
      fsvelint.Multiply(efsvel_,funct_);

    // compute gradient of fine-scale part of scalar value
    LINALG::Matrix<nsd_,1> fsgradphi(true);
    if (turbparams_->FSSGD())
      fsgradphi.Multiply(derxy_,fsphinp_[0]);

    double rhsint(0.0);
    GetRhsInt(rhsint,densnp[0],0);

    //--------------------------------------------------------------------
    // calculation of (fine-scale) subgrid diffusivity, subgrid-scale
    // velocity and stabilization parameter(s) at integration point
    //--------------------------------------------------------------------

    // subgrid-scale convective term
    LINALG::Matrix<nen_,1> sgconv(true);
    // subgrid-scale velocity vector in gausspoint
    LINALG::Matrix<nsd_,1> sgvelint(true);

    if (scatrapara_->TauGP())
    {
      // artificial diffusion / shock capturing: adaption of diffusion coefficient
      if (scatrapara_->ASSGD())
       dserror("Not supported");

      // calculation of all-scale subgrid diffusivity (by, e.g.,
      // Smagorinsky model) at element center
      if (turbparams_->TurbModel() == INPAR::FLUID::smagorinsky
          or turbparams_->TurbModel() == INPAR::FLUID::dynamic_smagorinsky
          or turbparams_->TurbModel() == INPAR::FLUID::dynamic_vreman)
      {
        CalcSubgrDiff(visc,vol,0,densnp[0]);
      }

      // calculation of fine-scale artificial subgrid diffusivity at element center
      // we have to set a vector, which is however only required for one
      // special case (computation of fine-scale subgrid diffusivity for non-incremental
      // solver -> only artificial subgrid diffusivity) not considered here
      Epetra_SerialDenseVector  elevec1_epetra_subgrdiff_dummy;
      if (turbparams_->FSSGD()) CalcFineScaleSubgrDiff(sgdiff,elevec1_epetra_subgrdiff_dummy,ele,vol,0,densnp[0],diffmanager_->GetIsotropicDiff(0),scatravarmanager_->ConVel());

      // calculation of subgrid-scale velocity at integration point if required
      if (scatrapara_->RBSubGrVel())
      {
        // calculation of stabilization parameter related to fluid momentum
        // equation at integration point
        CalcTau(tau[0],visc,0.0,densnp[0],convelint,vol);
        // calculation of residual-based subgrid-scale velocity
        CalcSubgrVelocity(ele,sgvelint,densam[0],densnp[0],visc,convelint,tau[0]);

        // calculation of subgrid-scale convective part
        sgconv.MultiplyTN(derxy_,sgvelint);
      }

      // calculation of stabilization parameter at integration point
      CalcTau(tau[0],diffmanager_->GetIsotropicDiff(0),reamanager_->GetReaCoeff(0),densnp[0],convelint,vol);
    }

    // prepare multifractal subgrid-scale modeling
    // calculation of model coefficients B (velocity) and D (scalar)
    // at Gauss point as well as calculation
    // of multifractal subgrid-scale quantities
    LINALG::Matrix<nsd_,1> mfsgvelint(true);
    double mfsvdiv(0.0);
    double mfssgphi(0.0);
    LINALG::Matrix<nsd_,1> mfsggradphi(true);
    if (turbparams_->TurbModel() == INPAR::FLUID::multifractal_subgrid_scales)
    {
      if (turbparams_->BD_Gp())
        // calculate model coefficients
        CalcBAndDForMultifracSubgridScales(B_mfs,D_mfs,vol,0,densnp[0],diffmanager_->GetIsotropicDiff(0),visc,convelint,fsvelint);

      // calculate fine-scale velocity, its derivative and divergence for multifractal subgrid-scale modeling
      for (int idim=0; idim<nsd_; idim++)
        mfsgvelint(idim,0) = fsvelint(idim,0) * B_mfs(idim,0);
      // required for conservative formulation in the context of passive scalar transport
      if (turbparams_->MfsConservative() or scatrapara_->IsConservative())
      {
        // get divergence of subgrid-scale velocity
        LINALG::Matrix<nsd_,nsd_> mfsvderxy;
        mfsvderxy.MultiplyNT(efsvel_,derxy_);
        for (int idim = 0; idim<nsd_; idim++)
          mfsvdiv += mfsvderxy(idim,idim) * B_mfs(idim,0);
      }

      // calculate fine-scale scalar and its derivative for multifractal subgrid-scale modeling
      mfssgphi = D_mfs * funct_.Dot(fsphinp_[0]);
      mfsggradphi.Multiply(derxy_,fsphinp_[0]);
      mfsggradphi.Scale(D_mfs);
    }

    // residual of convection-diffusion-reaction eq
    double scatrares(0.0);

    // compute residual of scalar transport equation and
    // subgrid-scale part of scalar
    CalcStrongResidual(0,scatrares,densam[0],densnp[0],rea_phi,rhsint,tau[0]);

    //--------------------------------------------------------------------
    // calculation of subgrid-scale part of scalar
    //--------------------------------------------------------------------
    double sgphi = -tau[0]*scatrares;

    // not supported anymore
    // update material parameters based on inclusion of subgrid-scale
    // part of scalar (active only for mixture fraction,
    // Sutherland law and progress variable, for the time being)
//    if (update_mat_)
//    {
//      if (turbmodel_ == INPAR::FLUID::multifractal_subgrid_scales)
//        UpdateMaterialParams(ele,mfssgphi_[0],0);
//      else
//        UpdateMaterialParams(ele,sgphi_[0],0);
//    }
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
    if (scatrapara_->StabType() == INPAR::SCATRA::stabtype_SUPG)
    {
      eps_supg -= densnp[0] * sgphi * convelint.Dot(gradphi) * fac; //sgphi_ is negative

    }
    else if (scatrapara_->StabType() == INPAR::SCATRA::stabtype_no_stabilization)
    {
       // nothing to do
    }
    else
      dserror("Stabtype not yet supported!");

    // dissipation by cross-stress-stabilization
    // dissipation by reynolds-stress-stabilization
    if (scatrapara_->RBSubGrVel())
    {
      eps_cross -= densnp[0] * phinp * sgvelint.Dot(gradphi) * fac;
      eps_rey -= densnp[0] * sgphi * sgvelint.Dot(gradphi) * fac;
    }

    //---------------------------------------------------------------
    // multifractal subgrid-scale modeling terms
    //---------------------------------------------------------------

    // dissipation multifractal subgrid-scales
    if (turbparams_->TurbModel() == INPAR::FLUID::multifractal_subgrid_scales)
    {
      eps_mfs -= densnp[0] * (mfssgphi * convelint.Dot(gradphi)
                              + phinp * mfsgvelint.Dot(gradphi)
                              + mfssgphi * mfsgvelint.Dot(gradphi)) * fac;
      eps_mfscross -= densnp[0] * (mfssgphi * convelint.Dot(gradphi)
                                   + phinp * mfsgvelint.Dot(gradphi)) * fac;
      eps_mfsrey -= densnp[0] * mfssgphi * mfsgvelint.Dot(gradphi) * fac;
    }

    //---------------------------------------------------------------
    // small-scale subgrid-viscosity subgrid-scale modeling terms
    //---------------------------------------------------------------

    // dissipation AVM3
    if (turbparams_->FSSGD())
    {
      eps_avm3 += sgdiff * fsgradphi.Dot(fsgradphi) * fac;
    }

    //---------------------------------------------------------------
    // Smagorinsky model
    //---------------------------------------------------------------

    // dissipation (Smagorinsky)
    if (scatrapara_->ASSGD()
        or turbparams_->TurbModel() == INPAR::FLUID::smagorinsky
        or turbparams_->TurbModel() == INPAR::FLUID::dynamic_smagorinsky)
    {
      eps_smag += diffmanager_->GetSubGrDiff(0) * gradphi.Dot(gradphi) * fac;
    }

    //---------------------------------------------------------------
    // standard Galerkin terms
    //---------------------------------------------------------------

    // convective (Galerkin)
    eps_conv -= densnp[0] * phinp * convelint.Dot(gradphi) * fac;

    // dissipation (Galerkin) (diffus is diffus+sgdiff here)
    eps_visc += (diffmanager_->GetIsotropicDiff(0) - diffmanager_->GetSubGrDiff(0)) * gradphi.Dot(gradphi) * fac;

    //---------------------------------------------------------------
    // element averages of tau_S and residual
    //---------------------------------------------------------------
    averaged_tauS += tau[0] * fac;
    mean_resS    += scatrares * fac;
    mean_resS_sq += scatrares * scatrares * fac;

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


// template classes

#include "scatra_ele_calc_fwd.hpp"
