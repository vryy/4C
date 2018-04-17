/*----------------------------------------------------------------------*/
/*!
\file fluid_ele_calc_turbulence_service.cpp

\brief Internal implementation of Fluid element

\level 2

<pre>
\maintainer Martin Kronbichler
            kronbichler@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15235
</pre>
*/
/*----------------------------------------------------------------------*/

#include "fluid_ele_calc.H"
#include "fluid_ele_parameter.H"

#include "../drt_inpar/inpar_turbulence.H"

#include "../drt_lib/drt_elementtype.H"


/*----------------------------------------------------------------------*
 |  compute turbulence parameters                       rasthofer 10/11 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype, DRT::ELEMENTS::Fluid::EnrichmentType enrtype>
void DRT::ELEMENTS::FluidEleCalc<distype,enrtype>::GetTurbulenceParams(
                               Teuchos::ParameterList&    turbmodelparams,
                               double&                    Cs_delta_sq,
                               double&                    Ci_delta_sq,
                               int&                       nlayer,
                               double                     CsDeltaSq,
                               double                     CiDeltaSq)
{
  if(fldpara_->TurbModAction() != INPAR::FLUID::no_model and nsd_ == 2)
    dserror("turbulence and 2D flow does not make any sense");

  // classical smagorinsky does only have constant parameter
  if (fldpara_->TurbModAction() == INPAR::FLUID::smagorinsky_with_van_Driest_damping)
  {
    // this will be the y-coordinate of a point in the element interior
    // we will determine the element layer in which he is contained to
    // be able to do the output of visceff etc.
    double center = 0.0;

    for(int inode=0;inode<nen_;inode++)
    {
      center += xyze_( 1, inode );
    }
    center/=nen_;

    // node coordinates of plane to the element layer
    Teuchos::RCP<std::vector<double> > planecoords
      =
      turbmodelparams.get<Teuchos::RCP<std::vector<double> > >("planecoords_");

    bool found = false;
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
  }
  // --------------------------------------------------
  // Smagorinsky model with dynamic Computation of Cs
  //else if (physical_turbulence_model == "Dynamic_Smagorinsky")
  else if (fldpara_->TurbModAction() == INPAR::FLUID::dynamic_smagorinsky)
  {
    //turb_mod_action_ = Fluid::dynamic_smagorinsky;

    // for homogeneous flow, use averaged quantities
    if (fldpara_->CsAveraged()==true){
    if (turbmodelparams.get<std::string>("HOMDIR","not_specified")
            !=
            "not_specified")
    {
      Teuchos::RCP<std::vector<double> > averaged_LijMij = turbmodelparams.get<Teuchos::RCP<std::vector<double> > >("averaged_LijMij_");
      Teuchos::RCP<std::vector<double> > averaged_MijMij = turbmodelparams.get<Teuchos::RCP<std::vector<double> > >("averaged_MijMij_");
      Teuchos::RCP<std::vector<double> > averaged_CI_numerator = Teuchos::null;
      Teuchos::RCP<std::vector<double> > averaged_CI_denominator = Teuchos::null;
      if (fldpara_->PhysicalType()==INPAR::FLUID::loma)
      {
        averaged_CI_numerator = turbmodelparams.get<Teuchos::RCP<std::vector<double> > >("averaged_CI_numerator_");
        averaged_CI_denominator = turbmodelparams.get<Teuchos::RCP<std::vector<double> > >("averaged_CI_denominator_");
      }

      // get homogeneous direction
      std::string homdir = turbmodelparams.get<std::string>("HOMDIR","not_specified");

      // here, the layer is determined in order to get the correct
      // averaged value from the vector of averaged (M/L)ijMij
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

      // Cs_delta_sq is set by the averaged quantities
      if ((*averaged_MijMij)[nlayer] > 1E-16)
        Cs_delta_sq = 0.5 * (*averaged_LijMij)[nlayer]/(*averaged_MijMij)[nlayer] ;
      else
        Cs_delta_sq = 0.0;

      // clipping to get algorithm stable
      if (Cs_delta_sq<0)
        Cs_delta_sq=0.0;

      // Ci_delta_sq is set by the averaged quantities
      if (fldpara_->PhysicalType()==INPAR::FLUID::loma and fldpara_->IncludeCi()==true)
      {
        if ((*averaged_CI_denominator)[nlayer] > 1E-16)
          Ci_delta_sq = 0.5 * (*averaged_CI_numerator)[nlayer]/(*averaged_CI_denominator)[nlayer] ;
        else
          Ci_delta_sq = 0.0;

        // clipping to get algorithm stable
        if (Ci_delta_sq<0.0)
          Ci_delta_sq=0.0;
      }

    }
    }
    else
    {
      // when no averaging was done, we just keep the calculated (clipped) value
      Cs_delta_sq = CsDeltaSq;
      if (fldpara_->PhysicalType()==INPAR::FLUID::loma and fldpara_->IncludeCi()==true)
        Ci_delta_sq = CiDeltaSq;
    }
  }

  return;
} // FluidEleCalc::GetTurbulenceParams


/*----------------------------------------------------------------------*
 |  calculation of (all-scale) subgrid viscosity               vg 09/09 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype, DRT::ELEMENTS::Fluid::EnrichmentType enrtype>
void DRT::ELEMENTS::FluidEleCalc<distype,enrtype>::CalcSubgrVisc(
  const LINALG::Matrix<nsd_,nen_>&        evelaf,
  const double                            vol,
  double&                                 Cs_delta_sq,
  double&                                 Ci_delta_sq
  )
{
  // cast dimension to a double varibale -> pow()
  const double dim = double (nsd_);
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

  // compute (all-scale) rate of strain
  double rateofstrain = -1.0e30;
  rateofstrain = GetStrainRate(evelaf);

  if (fldpara_->TurbModAction() == INPAR::FLUID::dynamic_smagorinsky)
  {
    // subgrid viscosity
    sgvisc_ = densaf_ * Cs_delta_sq * rateofstrain;

    // calculate isotropic part of subgrid-stress tensor (loma only)
    if (fldpara_->PhysicalType() == INPAR::FLUID::loma)
    {
      if (fldpara_->IncludeCi() and is_inflow_ele_ == false)
      {
        if (fldpara_->Ci()<0.0)
          q_sq_ = densaf_ * Ci_delta_sq * rateofstrain * rateofstrain; //remark: missing factor 2 is added in function ContStab()
        else
        {
          // get characteristic element length for Smagorinsky model for 2D and 3D
          // 3D: delta = V^1/3
          // 2D: delta = A^1/2
          const double delta = pow(vol,(1.0/dim));
          q_sq_ = densaf_ * fldpara_->Ci() * delta * delta * rateofstrain * rateofstrain;
        }
      }
      else
        q_sq_ = 0.0;
    }
  }
  //Vreman turbulence model according to
  //"An eddy-viscosity subgrid-scale model for turbulent shear flow: Algebraic theory and applications", 2004
  else if (fldpara_->TurbModAction() == INPAR::FLUID::vreman or fldpara_->TurbModAction() == INPAR::FLUID::dynamic_vreman)
  {
    if(nsd_==3)
    {
      const double cs = fldpara_->Cs();
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

      LINALG::Matrix<nsd_,nsd_> velderxy;

      velderxy.MultiplyNT(evelaf,derxy_);

      //calculate grid filter width: 3 options:
      //- direction dependent using the element length in x, y and z
      //- cube root of element volume
      //- minimum element length (of x, y, z)
      if (fldpara_->Vrfi()==INPAR::FLUID::dir_dep)
      {
        double xmin = 0.0;
        double ymin = 0.0;
        double zmin = 0.0;
        double xmax = 0.0;
        double ymax = 0.0;
        double zmax = 0.0;
        for (int inen=0; inen<nen_; inen++)
        {
          if (inen == 0)
          {
            xmin = xyze_(0,inen);
            xmax = xyze_(0,inen);
            ymin = xyze_(1,inen);
            ymax = xyze_(1,inen);
            zmin = xyze_(2,inen);
            zmax = xyze_(2,inen);
          }
          else
          {
            if(xyze_(0,inen)<xmin)
              xmin = xyze_(0,inen);
            if(xyze_(0,inen)>xmax)
              xmax = xyze_(0,inen);
            if(xyze_(1,inen)<ymin)
              ymin = xyze_(1,inen);
            if(xyze_(1,inen)>ymax)
              ymax = xyze_(1,inen);
            if(xyze_(2,inen)<zmin)
              zmin = xyze_(2,inen);
            if(xyze_(2,inen)>zmax)
              zmax = xyze_(2,inen);
          }
        }
        hkxpow2 =  (xmax - xmin) * (xmax - xmin); //filter width = 2 grid spacing?
        hkypow2 =  (ymax - ymin) * (ymax - ymin);
        hkzpow2 =  (zmax - zmin) * (zmax - zmin);
      }
      else if (fldpara_->Vrfi()==INPAR::FLUID::cuberootvol)
      {
        hkxpow2=pow(vol,(2.0/3.0));
        hkypow2=hkxpow2;
        hkzpow2=hkxpow2;
      }
      else //minimum element length
      {
        double hk=0.0;
        double xmin = 0.0;
        double ymin = 0.0;
        double zmin = 0.0;
        double xmax = 0.0;
        double ymax = 0.0;
        double zmax = 0.0;
        for (int inen=0; inen<nen_; inen++)
        {
          if (inen == 0)
          {
            xmin = xyze_(0,inen);
            xmax = xyze_(0,inen);
            ymin = xyze_(1,inen);
            ymax = xyze_(1,inen);
            zmin = xyze_(2,inen);
            zmax = xyze_(2,inen);
          }
          else
          {
            if(xyze_(0,inen)<xmin)
              xmin = xyze_(0,inen);
            if(xyze_(0,inen)>xmax)
              xmax = xyze_(0,inen);
            if(xyze_(1,inen)<ymin)
              ymin = xyze_(1,inen);
            if(xyze_(1,inen)>ymax)
              ymax = xyze_(1,inen);
            if(xyze_(2,inen)<zmin)
              zmin = xyze_(2,inen);
            if(xyze_(2,inen)>zmax)
              zmax = xyze_(2,inen);
          }
        }
        if ((xmax-xmin) < (ymax-ymin))
        {
          if ((xmax-xmin) < (zmax-zmin))
             hk = xmax-xmin;
        }
        else
        {
          if ((ymax-ymin) < (zmax-zmin))
             hk = ymax-ymin;
          else
             hk = zmax-zmin;
        }
        hkxpow2=hk*hk;
        hkypow2=hkxpow2;
        hkzpow2=hkxpow2;
      }


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
        sgvisc_=0.0;
      else
      {
        if (fldpara_->TurbModAction() == INPAR::FLUID::vreman)
          sgvisc_= densaf_ * cs * sqrt(bbeta / alphavreman); // c_vreman=2.5*(c_smagorinsky*c_smagorinsky)
        else
        {
          double Cv=Cs_delta_sq; //the variable of Cs_delta_sq has only been used to get the Vreman constant here.

          sgvisc_= densaf_ * Cv * sqrt(bbeta / alphavreman);
        }
      }


    }
    else
      dserror("Vreman model only for nsd_==3");
  }
  else
  {
    double van_Driest_damping = 1.0;
    if (fldpara_->TurbModAction() == INPAR::FLUID::smagorinsky_with_van_Driest_damping)
    {
      // since the Smagorinsky constant is only valid if hk is in the inertial
      // subrange of turbulent flows, the mixing length is damped in the
      // viscous near wall region using the van Driest damping function
      /*
                                     /         /   y+ \ \
                   lmix = Cs * hk * | 1 - exp | - ---- | |
                                     \         \   A+ / /
      */
      // A+ is a constant parameter, y+ the distance from the wall in wall
     // units
      const double A_plus = 26.0;
      double y_plus;

      // the integration point coordinate is defined by the isometric approach
      /*
                  +-----
                   \
              x =   +      N (x) * x
                   /        j       j
                  +-----
                  node j
      */

      LINALG::Matrix<nsd_,1> centernodecoord;
      centernodecoord.Multiply(xyze_,funct_);

      if (centernodecoord(1,0)>0) y_plus=(1.0-centernodecoord(1,0))/fldpara_->ltau();
      else                        y_plus=(1.0+centernodecoord(1,0))/fldpara_->ltau();

      // lmix *= (1.0-exp(-y_plus/A_plus));
      // multiply with van Driest damping function
      van_Driest_damping = (1.0-exp(-y_plus/A_plus));
      fldpara_->SetvanDriestdamping(van_Driest_damping);
    }

    // get characteristic element length for Smagorinsky model for 2D and 3D
    // 3D: hk = V^1/3
    // 2D: hk = A^1/2
    const double hk = std::pow(vol,(1.0/dim));

    // mixing length set proportional to grid width: lmix = Cs * hk
    double lmix = fldpara_->Cs() * van_Driest_damping * hk;

    Cs_delta_sq = lmix * lmix;

    // subgrid viscosity
    sgvisc_ = densaf_ * Cs_delta_sq * rateofstrain;

    // calculate isotropic part of subgrid-stress tensor (loma only)
    if (fldpara_->PhysicalType() == INPAR::FLUID::loma)
    {
      if (fldpara_->IncludeCi() and is_inflow_ele_ == false)
      {
        if (fldpara_->Ci()<0.0)
          dserror("Ci expected!");
        else
        {
          // get characteristic element length for Smagorinsky model for 2D and 3D
          // 3D: delta = V^1/3
          // 2D: delta = A^1/2
          const double delta = pow(vol,(1.0/dim));
          q_sq_ = densaf_ * fldpara_->Ci() * delta * delta * rateofstrain * rateofstrain;
        }
      }
      else
        q_sq_ = 0.0;
    }
  }
  return;
}


/*----------------------------------------------------------------------*
 |  calculation of fine-scale subgrid viscosity                vg 09/09 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype, DRT::ELEMENTS::Fluid::EnrichmentType enrtype>
void DRT::ELEMENTS::FluidEleCalc<distype,enrtype>::CalcFineScaleSubgrVisc(
  const LINALG::Matrix<nsd_,nen_>&        evelaf,
  const LINALG::Matrix<nsd_,nen_>&        fsevelaf,
  const double                            vol
  )
{
  // cast dimension to a double varibale -> pow()
  const double dim = double (nsd_);

  //     // get characteristic element length for Smagorinsky model for 2D and 3D
  // 3D: hk = V^1/3
  // 2D: hk = A^1/2
  const double hk = std::pow(vol,(1.0/dim));

  if (fldpara_->Fssgv() == INPAR::FLUID::smagorinsky_all)
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
    rateofstrain = GetStrainRate(evelaf);

    fssgvisc_ = densaf_ * fldpara_->Cs() * fldpara_->Cs() * hk * hk * rateofstrain;
  }
  else if (fldpara_->Fssgv() == INPAR::FLUID::smagorinsky_small)
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
    //                                            'resolved' rate of strain
    //

    // fine-scale rate of strain
    double fsrateofstrain = -1.0e30;
    fsrateofstrain = GetStrainRate(fsevelaf);

    fssgvisc_ = densaf_ * fldpara_->Cs() * fldpara_->Cs() * hk * hk * fsrateofstrain;
  }

  return;
}


/*----------------------------------------------------------------------*
 |  compute multifractal subgrid scales parameters    rasthofer 04/2011 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype, DRT::ELEMENTS::Fluid::EnrichmentType enrtype>
void DRT::ELEMENTS::FluidEleCalc<distype,enrtype>::PrepareMultifractalSubgrScales(
  LINALG::Matrix<nsd_,1>&           B_mfs,
  double &                          D_mfs,
  const LINALG::Matrix<nsd_,nen_>&  evelaf,
  const LINALG::Matrix<nsd_,nen_>&  fsevelaf,
  const double                      vol
)
{
    // set input parameters
    double Csgs = fldpara_->Csgs();
    double alpha = fldpara_->Alpha();

    // allocate vector for parameter N
    // N may depend on the direction
    std::vector<double> Nvel (3);

    // potential calculation of Re to determine N
    double Re_ele = -1.0;
    // characteristic element length
    double hk = 1.0e+10;
    double strainnorm = 0.0;
    // ratio of viscous scale to element length
    double scale_ratio = 0.0;

    // get norm
    const double vel_norm = velint_.Norm2();
    const double fsvel_norm = fsvelint_.Norm2();

    // do we have a fixed parameter N
    if (not fldpara_->CalcN())
    {
      for (int rr=1;rr<3;rr++)
        Nvel[rr] = fldpara_->N();
    }
    else //no, so we calculate N from Re
    {
      // calculate characteristic element length
      // cf. stabilization parameters
      switch (fldpara_->RefLength()){
      case INPAR::FLUID::streamlength:
      {
        // a) streamlength due to Tezduyar et al. (1992)
        // normed velocity vector
        LINALG::Matrix<nsd_,1> velino(true);
        if (vel_norm>=1e-6) velino.Update(1.0/vel_norm,velint_);
        else
        {
          velino.Clear();
          velino(0,0) = 1.0;
        }
        LINALG::Matrix<nen_,1> tmp;
        //enriched dofs are not interpolatory with respect to geometry
        if(enrtype == DRT::ELEMENTS::Fluid::xwall)
        {
          LINALG::Matrix<nsd_,nen_> derxy_copy(derxy_);
          for(int inode=1;inode<nen_;inode+=2)
          {
            for (int idim=0; idim<nsd_;idim++)
              derxy_copy(idim,inode)=0.0;
          }
          tmp.MultiplyTN(derxy_copy,velino);
        }
        else
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
            G(nn,rr) = xji_(nn,0)*xji_(rr,0);
            for (int mm=1;mm<nsd_;++mm)
            {
              G(nn,rr) += xji_(nn,mm)*xji_(rr,mm);
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
        if (nsd_ != 3) dserror("Turbulence is 3d!");
        LINALG::Matrix<nsd_,1> normed_velgrad;

        for (int rr=0;rr<nsd_;++rr)
        {
          double val = 0.0;
          for (int idim = 0; idim < nsd_; idim ++)
            val += vderxy_(idim,rr)*vderxy_(idim,rr);

          normed_velgrad(rr) = std::sqrt(val);

          //normed_velgrad(rr)=sqrt(vderxy_(0,rr)*vderxy_(0,rr)
          //                      +
          //                      vderxy_(1,rr)*vderxy_(1,rr)
          //                      +
          //                      vderxy_(2,rr)*vderxy_(2,rr));
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

          val += std::abs(loc);

          //val += abs( normed_velgrad(0)*derxy_(0,rr)
          //            +normed_velgrad(1)*derxy_(1,rr)
          //            +normed_velgrad(2)*derxy_(2,rr));
        } /* end of loop over element nodes */

        hk = 2.0/val;

        break;
      }
      default:
        dserror("Unknown length");
      }

// alternative length for comparison, currently not used
#ifdef HMIN // minimal element length
      double xmin = 0.0;
      double ymin = 0.0;
      double zmin = 0.0;
      double xmax = 0.0;
      double ymax = 0.0;
      double zmax = 0.0;
      for (int inen=0; inen<nen_; inen++)
      {
        if (inen == 0)
        {
          xmin = xyze_(0,inen);
          xmax = xyze_(0,inen);
          ymin = xyze_(1,inen);
          ymax = xyze_(1,inen);
          zmin = xyze_(2,inen);
          zmax = xyze_(2,inen);
        }
        else
        {
          if(xyze_(0,inen)<xmin)
            xmin = xyze_(0,inen);
          if(xyze_(0,inen)>xmax)
            xmax = xyze_(0,inen);
          if(xyze_(1,inen)<ymin)
            ymin = xyze_(1,inen);
          if(xyze_(1,inen)>ymax)
            ymax = xyze_(1,inen);
          if(xyze_(2,inen)<zmin)
            zmin = xyze_(2,inen);
          if(xyze_(2,inen)>zmax)
            zmax = xyze_(2,inen);
        }
      }
      if ((xmax-xmin) < (ymax-ymin))
      {
        if ((xmax-xmin) < (zmax-zmin))
           hk = xmax-xmin;
      }
      else
      {
        if ((ymax-ymin) < (zmax-zmin))
           hk = ymax-ymin;
        else
           hk = zmax-zmin;
      }
#endif
#ifdef HMAX // maximal element length
      double xmin = 0.0;
      double ymin = 0.0;
      double zmin = 0.0;
      double xmax = 0.0;
      double ymax = 0.0;
      double zmax = 0.0;
      for (int inen=0; inen<nen_; inen++)
      {
        if (inen == 0)
        {
          xmin = xyze_(0,inen);
          xmax = xyze_(0,inen);
          ymin = xyze_(1,inen);
          ymax = xyze_(1,inen);
          zmin = xyze_(2,inen);
          zmax = xyze_(2,inen);
        }
        else
        {
          if(xyze_(0,inen)<xmin)
            xmin = xyze_(0,inen);
          if(xyze_(0,inen)>xmax)
            xmax = xyze_(0,inen);
          if(xyze_(1,inen)<ymin)
            ymin = xyze_(1,inen);
          if(xyze_(1,inen)>ymax)
            ymax = xyze_(1,inen);
          if(xyze_(2,inen)<zmin)
            zmin = xyze_(2,inen);
          if(xyze_(2,inen)>zmax)
            zmax = xyze_(2,inen);
        }
      }
      if ((xmax-xmin) > (ymax-ymin))
      {
        if ((xmax-xmin) > (zmax-zmin))
           hk = xmax-xmin;
      }
      else
      {
        if ((ymax-ymin) > (zmax-zmin))
           hk = ymax-ymin;
        else
           hk = zmax-zmin;
      }
#endif

      if (hk == 1.0e+10)
        dserror("Something went wrong!");

      switch (fldpara_->RefVel()){
      case INPAR::FLUID::resolved:
      {
        Re_ele = vel_norm * hk * densaf_ / visc_;
        break;
      }
      case INPAR::FLUID::fine_scale:
      {
        Re_ele = fsvel_norm * hk * densaf_ / visc_;
        break;
      }
      case INPAR::FLUID::strainrate:
      {
        strainnorm = GetStrainRate(evelaf);
        strainnorm /= sqrt(2.0); //cf. Burton & Dahm 2005
        Re_ele = strainnorm * hk * hk * densaf_ / visc_;
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
      scale_ratio = fldpara_->CNu() * pow(Re_ele,0.75);
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
    if (fldpara_->NearWallLimit())
    {
      // if not yet calculated, estimate norm of strain rate
      if ((not fldpara_->CalcN()) or (fldpara_->RefVel() != INPAR::FLUID::strainrate))
      {
        //strainnorm = GetNormStrain(evelaf,derxy_,vderxy_);
        strainnorm = GetStrainRate(evelaf);
        strainnorm /= sqrt(2.0); //cf. Burton & Dahm 2005
      }
      // and reference length
      if (not fldpara_->CalcN()) dserror("hk not yet calculated"); // solution see scatra

      // get Re from strain rate
      double Re_ele_str = strainnorm * hk * hk * densaf_ / visc_;
      if (Re_ele_str < 0.0)
        dserror("Something went wrong!");
      // ensure positive values
      if (Re_ele_str < 1.0)
         Re_ele_str = 1.0;

      // calculate corrected Csgs
      //           -3/16
      //  *(1 - (Re)   )
      //
      Csgs *= (1.0-pow(Re_ele_str,-3.0/16.0));

      // store Cai for application to scalar field
      Cai_phi = (1.0-pow(Re_ele_str,-3.0/16.0));
    }

    // call function to compute coefficient B
    CalcMultiFracSubgridVelCoef(Csgs,alpha,Nvel,B_mfs);

    // prepare calculation of subgrid-scalar coefficient for loma
    // required if further subgrid-scale terms of cross- and Reynolds-stress
    // type arising in the continuity equation should be included
    if (fldpara_->PhysicalType() == INPAR::FLUID::loma)
    {
      // set input parameters
      double Csgs_phi = fldpara_->CsgsPhi();

      // calculate prandtl number
      double Pr = visc_/diffus_;

      // since there are differences in the physical behavior between low and high
      // Prandtl/Schmidt number regime, we define a limit
      // to distinguish between the low and high Prandtl/Schmidt number regime
      // note: there is no clear definition of the ranges
      const double Pr_limit = 2.0;

      // allocate double for parameter N
      double Nphi = 0.0;
      // ratio of dissipation scale to element length
      double scale_ratio_phi = 0.0;

      if (fldpara_->CalcN())
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

        scale_ratio_phi = fldpara_->CDiff() * pow(Re_ele,0.75) * pow(Pr,p);
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

      // call function to compute coefficient D
      if (not fldpara_->NearWallLimitScatra())
        CalcMultiFracSubgridScaCoef(Csgs_phi,alpha,Pr,Pr_limit,Nvel,Nphi,D_mfs);
      else
      {
        if (not fldpara_->NearWallLimit())
          dserror("Near-wall limit expected!");

        CalcMultiFracSubgridScaCoef(Csgs_phi*Cai_phi,alpha,Pr,Pr_limit,Nvel,Nphi,D_mfs);
      }
    }
}



/*-------------------------------------------------------------------------------*
 |calculation parameter for multifractal subgrid scale modeling  rasthofer 03/11 |
 *-------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype, DRT::ELEMENTS::Fluid::EnrichmentType enrtype>
void DRT::ELEMENTS::FluidEleCalc<distype,enrtype>::CalcMultiFracSubgridVelCoef(
  const double            Csgs,
  const double            alpha,
  const std::vector<double> Nvel,
  LINALG::Matrix<nsd_,1>& B_mfs
  )
{
  //
  //          |       1              |
  //  kappa = | -------------------- |
  //          |  1 - alpha ^ (-4/3)  |
  //
  double kappa = 1.0/(1.0-pow(alpha,-4.0/3.0));

  //                  1                                    1
  //                  2                 |                 |2
  //  B = Csgs * kappa * 2 ^ (-2*N/3) * | 2 ^ (4*N/3) - 1 |
  //                                    |                 |
  //
  for (int dim=0; dim<nsd_; dim++)
  {
    B_mfs(dim,0) = Csgs *sqrt(kappa) * pow(2.0,-2.0*Nvel[dim]/3.0) * sqrt((pow(2.0,4.0*Nvel[dim]/3.0)-1.0));
  }

//  if (eid_ == 100){
//    std::cout << "B  " << std::setprecision(10) << B_mfs(0,0) << "  " << B_mfs(1,0) << "  " << B_mfs(2,0) << "  " << std::endl;
//    std::cout << "CsgsB  " << std::setprecision(10) << Csgs << std::endl;
//  }

  return;
}


/*-------------------------------------------------------------------------------*
 |calculation parameter for multifractal subgrid scale modeling  rasthofer 02/12 |
 |subgrid-scale scalar for loma                                                  |
 *-------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype, DRT::ELEMENTS::Fluid::EnrichmentType enrtype>
void DRT::ELEMENTS::FluidEleCalc<distype,enrtype>::CalcMultiFracSubgridScaCoef(
  const double            Csgs,
  const double            alpha,
  const double            Pr,
  const double            Pr_limit,
  const std::vector<double>    Nvel,
  double                  Nphi,
  double &                D_mfs
  )
{
  // here, we have to distinguish tree different cases:
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
  if (Pr < Pr_limit) // Pr <= 1, i.e., case 1 and 3
    gamma = 4.0/3.0;
  else // Pr >> 1
  {
    dserror("Loma with Pr>>1?");
    if (Nvel[0]<1.0) // Pr >> 1 and fluid fully resolved, i.e., case 2 (ii)
      gamma = 2.0;
    else // Pr >> 1 and fluid not fully resolved, i.e., case 2 (i)
    {
      if (Nvel[0]>Nphi)
       dserror("Nvel < Nphi expected!");
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
    D_mfs = Csgs *sqrt(kappa_phi) * pow(2.0,-gamma*Nphi/2.0) * sqrt((pow(2.0,gamma*Nphi)-1.0));
  else
  {
    dserror("Special option for passive scalars only!");
//    double gamma1 = 4.0/3.0;
//    double gamma2 = 2.0;
//    kappa_phi = 1.0/(1.0-pow(alpha,-gamma1));
//    D_mfs = Csgs * sqrt(kappa_phi) * pow(2.0,-gamma2*Nphi/2.0) * sqrt((pow(2.0,gamma1*Nvel[0])-1)+4.0/3.0*(M_PI/hk)*(pow(2.0,gamma2*Nphi)-pow(2.0,gamma2*Nvel[0])));
  }

//  if (eid_ == 100){
//    std::cout << "D  " << std::setprecision(10) << D_mfs << std::endl;
//    std::cout << "CsgsD  " << std::setprecision(10) << Csgs << std::endl;
//  }

  return;
}


template <DRT::Element::DiscretizationType distype, DRT::ELEMENTS::Fluid::EnrichmentType enrtype>
void DRT::ELEMENTS::FluidEleCalc<distype,enrtype>::FineScaleSubGridViscosityTerm(
    LINALG::Matrix<nsd_,nen_> &             velforce,
    const double &                          fssgviscfac)
{
  if (nsd_ == 2)
  {
    for (int vi=0; vi<nen_; ++vi)
    {
      /* fine-scale subgrid-viscosity term on right hand side */
      /*
                          /                          \
                         |       /    \         / \   |
         - mu_art(fsu) * |  eps | Dfsu | , eps | v |  |
                         |       \    /         \ /   |
                          \                          /
      */
      velforce(0, vi) -= fssgviscfac*(2.0*derxy_(0, vi)*fsvderxy_(0, 0)
                                     +    derxy_(1, vi)*fsvderxy_(0, 1)
                                     +    derxy_(1, vi)*fsvderxy_(1, 0)) ;
      velforce(1, vi) -= fssgviscfac*(    derxy_(0, vi)*fsvderxy_(0, 1)
                                     +    derxy_(0, vi)*fsvderxy_(1, 0)
                                     +2.0*derxy_(1, vi)*fsvderxy_(1, 1)) ;
    }
  }
  else if(nsd_ == 3)
  {
    for (int vi=0; vi<nen_; ++vi)
    {
      /* fine-scale subgrid-viscosity term on right hand side */
      /*
                            /                          \
                           |       /    \         / \   |
           - mu_art(fsu) * |  eps | Dfsu | , eps | v |  |
                           |       \    /         \ /   |
                            \                          /
      */
      velforce(0, vi) -= fssgviscfac*(2.0*derxy_(0, vi)*fsvderxy_(0, 0)
                                     +    derxy_(1, vi)*fsvderxy_(0, 1)
                                     +    derxy_(1, vi)*fsvderxy_(1, 0)
                                     +    derxy_(2, vi)*fsvderxy_(0, 2)
                                     +    derxy_(2, vi)*fsvderxy_(2, 0)) ;
      velforce(1, vi) -= fssgviscfac*(    derxy_(0, vi)*fsvderxy_(0, 1)
                                     +    derxy_(0, vi)*fsvderxy_(1, 0)
                                     +2.0*derxy_(1, vi)*fsvderxy_(1, 1)
                                     +    derxy_(2, vi)*fsvderxy_(1, 2)
                                     +    derxy_(2, vi)*fsvderxy_(2, 1)) ;
      velforce(2, vi) -= fssgviscfac*(    derxy_(0, vi)*fsvderxy_(0, 2)
                                     +    derxy_(0, vi)*fsvderxy_(2, 0)
                                     +    derxy_(1, vi)*fsvderxy_(1, 2)
                                     +    derxy_(1, vi)*fsvderxy_(2, 1)
                                     +2.0*derxy_(2, vi)*fsvderxy_(2, 2)) ;
    }
  }
  else dserror("fine-scale subgrid viscosity not implemented for 1-D problems!");

  return;
}


//----------------------------------------------------------------------
// Cross-stress terms: multifractal subgrid-scales       rasthofer 06/11
//----------------------------------------------------------------------
template <DRT::Element::DiscretizationType distype, DRT::ELEMENTS::Fluid::EnrichmentType enrtype>
void DRT::ELEMENTS::FluidEleCalc<distype,enrtype>::MultfracSubGridScalesCross(
    LINALG::Matrix<nen_*nsd_,nen_*nsd_> &   estif_u,
    LINALG::Matrix<nsd_,nen_> &             velforce,
    const double &                          timefacfac,
    const double &                          rhsfac)
{
  //--------------------------------------------------------------------
  // rhs contribution
  //--------------------------------------------------------------------
  if (nsd_ == 3)
  {
    for (int vi=0; vi<nen_; ++vi)
    {
      /* cross-stress term on right hand side */
      /*
               /                                      \
              |                                        |
              | ( du o nabla u + u o nabla du ) ,  v   |
              |                                        |
               \                                      /
      */
      velforce(0,vi) -= rhsfac * densaf_ * funct_(vi,0)
                      * (convvelint_(0,0) * mffsvderxy_(0,0)
                        +convvelint_(1,0) * mffsvderxy_(0,1)
                        +convvelint_(2,0) * mffsvderxy_(0,2)
                        +mffsvelint_(0,0) * vderxy_(0,0)
                        +mffsvelint_(1,0) * vderxy_(0,1)
                        +mffsvelint_(2,0) * vderxy_(0,2));
      velforce(1,vi) -= rhsfac * densaf_ * funct_(vi,0)
                      * (convvelint_(0,0) * mffsvderxy_(1,0)
                        +convvelint_(1,0) * mffsvderxy_(1,1)
                        +convvelint_(2,0) * mffsvderxy_(1,2)
                        +mffsvelint_(0,0) * vderxy_(1,0)
                        +mffsvelint_(1,0) * vderxy_(1,1)
                        +mffsvelint_(2,0) * vderxy_(1,2));
      velforce(2,vi) -= rhsfac * densaf_ * funct_(vi,0)
                      * (convvelint_(0,0) * mffsvderxy_(2,0)
                        +convvelint_(1,0) * mffsvderxy_(2,1)
                        +convvelint_(2,0) * mffsvderxy_(2,2)
                        +mffsvelint_(0,0) * vderxy_(2,0)
                        +mffsvelint_(1,0) * vderxy_(2,1)
                        +mffsvelint_(2,0) * vderxy_(2,2));

      /* cross-stress term on right hand side */
      /* additional terms conservative form */
      /*
               /                                         \
              |                                           |
              | ( du (nabla o u) + u (nabla o du ) ,  v   |
              |                                           |
               \                                         /
      */
      if (fldpara_->MfsIsConservative() or fldpara_->IsConservative())
      {
        velforce(0,vi) -= rhsfac * densaf_ * funct_(vi,0)
                        * (mffsvelint_(0,0) * vdiv_
                          +convvelint_(0,0) * mffsvdiv_);
        velforce(1,vi) -= rhsfac * densaf_ * funct_(vi,0)
                        * (mffsvelint_(1,0) * vdiv_
                          +convvelint_(1,0) * mffsvdiv_);
        velforce(2,vi) -= rhsfac * densaf_ * funct_(vi,0)
                        * (mffsvelint_(2,0) * vdiv_
                          +convvelint_(2,0) * mffsvdiv_);

        if (fldpara_->PhysicalType() == INPAR::FLUID::loma)
        {
          dserror("Conservative formulation not supported for loma!");
        }
        if(gridvelint_.Norm2()>1.0e-9)
          dserror("Conservative formulation not supported with ale!");
      }
    }
  }
  else
    dserror("Multifractal subgrid-scale modeling model for 3D-problems only!");

  //--------------------------------------------------------------------
  // lhs contribution
  //--------------------------------------------------------------------
  // linearized as far as possible due to the filter

  // turn left-hand-side contribution on
  double beta = fldpara_->Beta();

  if (beta > 1.0e-9)
  {
    if(gridvelint_.Norm2()>1.0e-9)
      dserror("left hand side terms of MFS not supported with ale");
    LINALG::Matrix<nen_,1> mfconv_c(true);
    mfconv_c.MultiplyTN(derxy_,mffsvelint_);
    // convective part
    for (int ui=0; ui<nen_; ui++)
    {
      for (int idim=0; idim<nsd_; idim++)
      {
        int fui = ui * nsd_ + idim;
        for (int vi=0; vi<nen_; vi++)
        {
          for (int jdim=0; jdim<nsd_; jdim++)
          {
            int fvi = vi * nsd_ + jdim;
            /*
                    /                             \
                   |  /                 \          |
                   | |   rho*Du  o nabla | du , v  |
                   |  \                 /          |
                    \                             /
            */
            estif_u(fvi,fui) += beta * timefacfac * densaf_ * funct_(vi)
                              * funct_(ui) * mffsvderxy_(jdim,idim);
            /*
                    /                             \
                   |  /                 \          |
                   | |   rho*du  o nabla | Du , v  |
                   |  \                 /          |
                    \                             /
            */
            if (jdim == idim)
            {
              estif_u(fvi,fui) += beta * timefacfac * densaf_ * funct_(vi)
                                * mfconv_c(ui);
            }

            // additional terms conservative part
            if (fldpara_->MfsIsConservative() or fldpara_->IsConservative())
            {
              /*
                   /                                     \
                   |      /               \       \      |
                   |  du | rho*nabla o Du  | , v   |     |
                   |      \               /       /      |
                   \                                     /
              */
              estif_u(fvi,fui) += beta * timefacfac * densaf_ * funct_(vi)
                                * mffsvelint_(jdim) * derxy_(idim, ui);
              /*
                    /                                     \
                    |      /               \       \      |
                    |  Du | rho*nabla o du  | , v   |     |
                    |      \               /       /      |
                    \                                     /
              */
              if (jdim == idim)
              {
                estif_u(fvi,fui) += beta * timefacfac * densaf_
                                  * funct_(vi) * funct_(ui) * mffsvdiv_;
              }

              if (fldpara_->PhysicalType() == INPAR::FLUID::loma)
              {
                dserror("Conservative formulation not supported for loma!");
              }
            }
          }
        }
      }
    }
  }

  return;
}


//----------------------------------------------------------------------
// Reynolds-stress terms: multifractal subgrid-scales    rasthofer 06/11
//----------------------------------------------------------------------
template <DRT::Element::DiscretizationType distype, DRT::ELEMENTS::Fluid::EnrichmentType enrtype>
void DRT::ELEMENTS::FluidEleCalc<distype,enrtype>::MultfracSubGridScalesReynolds(
    LINALG::Matrix<nen_*nsd_,nen_*nsd_> &   estif_u,
    LINALG::Matrix<nsd_,nen_> &             velforce,
    const double &                          timefacfac,
    const double &                          rhsfac)
{
  //--------------------------------------------------------------------
  // rhs contribution
  //--------------------------------------------------------------------
  if (nsd_ == 3)
  {
    for (int vi=0; vi<nen_; ++vi)
    {
      /* reynolds-stress term on right hand side */
      /*
               /                       \
              |                         |
              | ( du o nabla du) ,  v   |
              |                         |
               \                       /
      */
      velforce(0,vi) -= rhsfac * densaf_ * funct_(vi,0)
                      * (mffsvelint_(0,0) * mffsvderxy_(0,0)
                        +mffsvelint_(1,0) * mffsvderxy_(0,1)
                        +mffsvelint_(2,0) * mffsvderxy_(0,2));
      velforce(1,vi) -= rhsfac * densaf_ * funct_(vi,0)
                      * (mffsvelint_(0,0) * mffsvderxy_(1,0)
                        +mffsvelint_(1,0) * mffsvderxy_(1,1)
                        +mffsvelint_(2,0) * mffsvderxy_(1,2));
      velforce(2,vi) -= rhsfac * densaf_ * funct_(vi,0)
                      * (mffsvelint_(0,0) * mffsvderxy_(2,0)
                        +mffsvelint_(1,0) * mffsvderxy_(2,1)
                        +mffsvelint_(2,0) * mffsvderxy_(2,2));

      /* reynolds-stress term on right hand side */
      /* additional terms conservative form */
      /*
               /                       \
              |                         |
              |   du (nabla o du),  v   |
              |                         |
               \                       /
      */
      if (fldpara_->MfsIsConservative() or fldpara_->IsConservative())
      {
        velforce(0,vi) -= rhsfac * densaf_ * funct_(vi,0)
                        * (mffsvelint_(0,0) * mffsvdiv_);
        velforce(1,vi) -= rhsfac * densaf_ * funct_(vi,0)
                        * (mffsvelint_(1,0) * mffsvdiv_);
        velforce(2,vi) -= rhsfac * densaf_ * funct_(vi,0)
                        * (mffsvelint_(2,0) * mffsvdiv_);

        if (fldpara_->PhysicalType() == INPAR::FLUID::loma)
        {
          dserror("Conservative formulation not supported for loma!");
        }
      }
    }
  }
  else
    dserror("Multifractal subgrid-scale modeling for 3D-problems only!");

  //--------------------------------------------------------------------
  // lhs contribution
  //--------------------------------------------------------------------
  // no contribution, due to necessary linearization of filter

  return;
}


//----------------------------------------------------------------------
// concistency terms for residual-based stabilization           bk 10/14
//----------------------------------------------------------------------
template <DRT::Element::DiscretizationType distype, DRT::ELEMENTS::Fluid::EnrichmentType enrtype>
void DRT::ELEMENTS::FluidEleCalc<distype,enrtype>::MultfracSubGridScalesConsistentResidual()
{

  if (fldpara_->TurbModAction() == INPAR::FLUID::multifractal_subgrid_scales && fldpara_->ConsistentMFSResidual())
  {
    if(not fldparatimint_->IsGenalpha())
      dserror("please check the implementation of the consistent residual for mfs for your time-integrator!");

    LINALG::Matrix<nsd_,1>              velforce(true);

      /* cross-stress term of residual */
      /*
               /                                      \
              |                                        |
              | ( du o nabla u + u o nabla du )        |
              |                                        |
               \                                      /
      */
      velforce(0,0) += densaf_
                      * (convvelint_(0,0) * mffsvderxy_(0,0)
                        +convvelint_(1,0) * mffsvderxy_(0,1)
                        +convvelint_(2,0) * mffsvderxy_(0,2)
                        +mffsvelint_(0,0) * vderxy_(0,0)
                        +mffsvelint_(1,0) * vderxy_(0,1)
                        +mffsvelint_(2,0) * vderxy_(0,2));
      velforce(1,0) += densaf_
                      * (convvelint_(0,0) * mffsvderxy_(1,0)
                        +convvelint_(1,0) * mffsvderxy_(1,1)
                        +convvelint_(2,0) * mffsvderxy_(1,2)
                        +mffsvelint_(0,0) * vderxy_(1,0)
                        +mffsvelint_(1,0) * vderxy_(1,1)
                        +mffsvelint_(2,0) * vderxy_(1,2));
      velforce(2,0) += densaf_
                      * (convvelint_(0,0) * mffsvderxy_(2,0)
                        +convvelint_(1,0) * mffsvderxy_(2,1)
                        +convvelint_(2,0) * mffsvderxy_(2,2)
                        +mffsvelint_(0,0) * vderxy_(2,0)
                        +mffsvelint_(1,0) * vderxy_(2,1)
                        +mffsvelint_(2,0) * vderxy_(2,2));


      /* reynolds-stress term of residual */
      /*
               /                       \
              |                         |
              | ( du o nabla du)        |
              |                         |
               \                       /
      */
      velforce(0,0) += densaf_
                      * (mffsvelint_(0,0) * mffsvderxy_(0,0)
                        +mffsvelint_(1,0) * mffsvderxy_(0,1)
                        +mffsvelint_(2,0) * mffsvderxy_(0,2));
      velforce(1,0) += densaf_
                      * (mffsvelint_(0,0) * mffsvderxy_(1,0)
                        +mffsvelint_(1,0) * mffsvderxy_(1,1)
                        +mffsvelint_(2,0) * mffsvderxy_(1,2));
      velforce(2,0) += densaf_
                      * (mffsvelint_(0,0) * mffsvderxy_(2,0)
                        +mffsvelint_(1,0) * mffsvderxy_(2,1)
                        +mffsvelint_(2,0) * mffsvderxy_(2,2));

    for (int rr=0;rr<nsd_;++rr)
    {
      momres_old_(rr)+=velforce(rr);
    }
  }
  return;
}

//----------------------------------------------------------------------
// outpu for statistics of dynamic Smagorinsky           rasthofer 09/12
//----------------------------------------------------------------------
template <DRT::Element::DiscretizationType distype, DRT::ELEMENTS::Fluid::EnrichmentType enrtype>
void DRT::ELEMENTS::FluidEleCalc<distype,enrtype>::StoreModelParametersForOutput(
   const double Cs_delta_sq,
   const double Ci_delta_sq,
   const int    nlayer,
   const bool   isowned,
   Teuchos::ParameterList&  turbmodelparams)
{
  // do the fastest test first
  if (isowned)
  {
    if (fldpara_->TurbModAction() == INPAR::FLUID::dynamic_smagorinsky or
        fldpara_->TurbModAction() == INPAR::FLUID::smagorinsky_with_van_Driest_damping)
    {
      if (turbmodelparams.get<std::string>("TURBULENCE_APPROACH", "none") == "CLASSICAL_LES")
      {
        if (turbmodelparams.get<std::string>("CANONICAL_FLOW","no") == "channel_flow_of_height_2" or
            turbmodelparams.get<std::string>("CANONICAL_FLOW","no") == "loma_channel_flow_of_height_2" or
            turbmodelparams.get<std::string>("CANONICAL_FLOW","no") == "scatra_channel_flow_of_height_2")
        {
          // recompute delta = pow((vol),(1.0/3.0))
          // evaluate shape functions and derivatives at element center
          EvalShapeFuncAndDerivsAtEleCenter();
          // set element volume
          const double vol = fac_;

          // to compare it with the standard Smagorinsky Cs
          if (fldpara_->TurbModAction() == INPAR::FLUID::dynamic_smagorinsky)
            (*(turbmodelparams.get<Teuchos::RCP<std::vector<double> > >("local_Cs_sum")))         [nlayer]+=sqrt(Cs_delta_sq)/pow((vol),(1.0/3.0));
          else if (fldpara_->TurbModAction() == INPAR::FLUID::smagorinsky_with_van_Driest_damping)
            (*(turbmodelparams.get<Teuchos::RCP<std::vector<double> > >("local_Cs_sum")))         [nlayer]+=fldpara_->Cs() * fldpara_->VanDriestdamping();
          else dserror("Dynamic Smagorinsky or Smagorisnsky with van Driest damping expected!");
          (*(turbmodelparams.get<Teuchos::RCP<std::vector<double> > >("local_Cs_delta_sq_sum")))[nlayer]+=Cs_delta_sq;
          (*(turbmodelparams.get<Teuchos::RCP<std::vector<double> > >("local_visceff_sum")))    [nlayer]+=visceff_;

          if (turbmodelparams.get<std::string>("CANONICAL_FLOW","no") == "loma_channel_flow_of_height_2")
          {
            (*(turbmodelparams.get<Teuchos::RCP<std::vector<double> > >("local_Ci_sum")))         [nlayer]+=sqrt(Ci_delta_sq)/pow((vol),(1.0/3.0));
            (*(turbmodelparams.get<Teuchos::RCP<std::vector<double> > >("local_Ci_delta_sq_sum")))[nlayer]+=Ci_delta_sq;
          }
        }
      }
    }
  }

  return;
}


template class DRT::ELEMENTS::FluidEleCalc<DRT::Element::hex8,DRT::ELEMENTS::Fluid::none>;
template class DRT::ELEMENTS::FluidEleCalc<DRT::Element::hex8,DRT::ELEMENTS::Fluid::xwall>;
template class DRT::ELEMENTS::FluidEleCalc<DRT::Element::tet4,DRT::ELEMENTS::Fluid::xwall>;
template class DRT::ELEMENTS::FluidEleCalc<DRT::Element::hex20,DRT::ELEMENTS::Fluid::none>;
template class DRT::ELEMENTS::FluidEleCalc<DRT::Element::hex27,DRT::ELEMENTS::Fluid::none>;
template class DRT::ELEMENTS::FluidEleCalc<DRT::Element::tet4,DRT::ELEMENTS::Fluid::none>;
template class DRT::ELEMENTS::FluidEleCalc<DRT::Element::tet10,DRT::ELEMENTS::Fluid::none>;
template class DRT::ELEMENTS::FluidEleCalc<DRT::Element::wedge6,DRT::ELEMENTS::Fluid::none>;
template class DRT::ELEMENTS::FluidEleCalc<DRT::Element::wedge15,DRT::ELEMENTS::Fluid::none>;
template class DRT::ELEMENTS::FluidEleCalc<DRT::Element::pyramid5,DRT::ELEMENTS::Fluid::none>;
template class DRT::ELEMENTS::FluidEleCalc<DRT::Element::quad4,DRT::ELEMENTS::Fluid::none>;
template class DRT::ELEMENTS::FluidEleCalc<DRT::Element::quad8,DRT::ELEMENTS::Fluid::none>;
template class DRT::ELEMENTS::FluidEleCalc<DRT::Element::quad9,DRT::ELEMENTS::Fluid::none>;
template class DRT::ELEMENTS::FluidEleCalc<DRT::Element::tri3,DRT::ELEMENTS::Fluid::none>;
template class DRT::ELEMENTS::FluidEleCalc<DRT::Element::tri6,DRT::ELEMENTS::Fluid::none>;
template class DRT::ELEMENTS::FluidEleCalc<DRT::Element::nurbs9,DRT::ELEMENTS::Fluid::none>;
template class DRT::ELEMENTS::FluidEleCalc<DRT::Element::nurbs27,DRT::ELEMENTS::Fluid::none>;
