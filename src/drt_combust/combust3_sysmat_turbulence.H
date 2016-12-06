//*----------------------------------------------------------------------*/
/*!
\file combust3_sysmat_turbulence.H

\brief  and rhs term of turbulence models as well as helper functions

\level 2

<pre>
\maintainer Benedikt Schott
            schott@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15241
</pre>
*/
/*----------------------------------------------------------------------*/

#ifndef COMBUST3_SYSMAT_TURBULENCE_H_
#define COMBUST3_SYSMAT_TURBULENCE_H_

#include "../drt_inpar/inpar_turbulence.H"
#include "../drt_inpar/inpar_fluid.H"

#include "combust3_utils.H"


namespace COMBUST
{

template <DRT::Element::DiscretizationType DISTYPE,
          XFEM::AssemblyType ASSTYPE,
          size_t NUMDOF,
          size_t shpVecSize>
void BuildTurbulenceDomainIntegrals(
    const INPAR::FLUID::FineSubgridVisc    fssgv,
    const INPAR::FLUID::TurbModelAction    turbmodel,
    COMBUST::LocalAssembler<DISTYPE,ASSTYPE,NUMDOF>& assembler,
    const XFEM::ApproxFunc<2,shpVecSize>& shpvel,      // shape functions and their derivatives (velocity)
    const LINALG::Matrix<3,3>& fsvderxy,               // velocity derivative at GP at n+1
    const LINALG::Matrix<3,1>& velint,
    const LINALG::Matrix<3,1>& mffsvelint,
    const LINALG::Matrix<3,3>& vderxy,
    const LINALG::Matrix<3,3>& mffsvderxy,
    const double vdiv,
    const double mffsvdiv,
    const bool   genalpha,                             // factor for integration in space
    const double fac,                                  // factor for integration in space
    const double timefac,                              // factor for integration in time
    const double rhsfac,                               // overall integration factor for the rhs
    const double dens,                                 // density at n+1 (OST), alpha_f (GenAlpha)
    const double fssgvisc,                             // finescale subgrid viscosity
    const bool   mfsconservative
)
{
  const double rhsdensfac    = dens * rhsfac;
//  if (genalpha)
//    return rhsfac;
//  else
//    return rhsfac*timefac;
  const double altrhsdensfac = (genalpha)?(rhsdensfac):(rhsdensfac*timefac);

  //---------------------------------------------------------------------
  //      fine-scale subgrid-viscosity term on right hand side
  //---------------------------------------------------------------------
  if (fssgv != INPAR::FLUID::no_fssgv)
  {
    double fssgviscfac=fssgvisc*altrhsdensfac;
    /*
                            /                          \
                           |       /    \         / \   |
           - mu_art(fsu) * |  eps | Dfsu | , eps | v |  |
                           |       \    /         \ /   |
                            \                          /
      */

    assembler.template Vector<XFEM::PHYSICS::Velx>(shpvel.dx,-fssgviscfac*(fsvderxy(0, 0)+fsvderxy(0, 0)));
    assembler.template Vector<XFEM::PHYSICS::Velx>(shpvel.dy,-fssgviscfac*(fsvderxy(0, 1)+fsvderxy(1, 0)));
    assembler.template Vector<XFEM::PHYSICS::Velx>(shpvel.dz,-fssgviscfac*(fsvderxy(0, 2)+fsvderxy(2, 0)));

    assembler.template Vector<XFEM::PHYSICS::Vely>(shpvel.dx,-fssgviscfac*(fsvderxy(0, 1)+fsvderxy(1, 0)));
    assembler.template Vector<XFEM::PHYSICS::Vely>(shpvel.dy,-fssgviscfac*(fsvderxy(1, 1)+fsvderxy(1, 1)));
    assembler.template Vector<XFEM::PHYSICS::Vely>(shpvel.dz,-fssgviscfac*(fsvderxy(1, 2)+fsvderxy(2, 1)));

    assembler.template Vector<XFEM::PHYSICS::Velz>(shpvel.dx,-fssgviscfac*(fsvderxy(0, 2)+fsvderxy(2, 0)));
    assembler.template Vector<XFEM::PHYSICS::Velz>(shpvel.dy,-fssgviscfac*(fsvderxy(1, 2)+fsvderxy(2, 1)));
    assembler.template Vector<XFEM::PHYSICS::Velz>(shpvel.dz,-fssgviscfac*(fsvderxy(2, 2)+fsvderxy(2, 2)));
  }

  //---------------------------------------------------------------------
  //    multifractal subgrid-scale modeling terms on right hand side
  //---------------------------------------------------------------------
  if (turbmodel == INPAR::FLUID::multifractal_subgrid_scales)
  {
    /* cross-stress term on right hand side */
    /*
               /                                      \
              |                                        |
       - dens | ( du o nabla u + u o nabla du ) ,  v   |
              |                                        |
               \                                      /
    */

    assembler.template Vector<XFEM::PHYSICS::Velx>(shpvel.d0,-altrhsdensfac * (velint(0,0) * mffsvderxy(0,0)
                                                                           +velint(1,0) * mffsvderxy(0,1)
                                                                           +velint(2,0) * mffsvderxy(0,2)
                                                                           +mffsvelint(0,0) * vderxy(0,0)
                                                                           +mffsvelint(1,0) * vderxy(0,1)
                                                                           +mffsvelint(2,0) * vderxy(0,2)));
    assembler.template Vector<XFEM::PHYSICS::Vely>(shpvel.d0,-altrhsdensfac * (velint(0,0) * mffsvderxy(1,0)
                                                                           +velint(1,0) * mffsvderxy(1,1)
                                                                           +velint(2,0) * mffsvderxy(1,2)
                                                                           +mffsvelint(0,0) * vderxy(1,0)
                                                                           +mffsvelint(1,0) * vderxy(1,1)
                                                                           +mffsvelint(2,0) * vderxy(1,2)));
    assembler.template Vector<XFEM::PHYSICS::Velz>(shpvel.d0,-altrhsdensfac * (velint(0,0) * mffsvderxy(2,0)
                                                                           +velint(1,0) * mffsvderxy(2,1)
                                                                           +velint(2,0) * mffsvderxy(2,2)
                                                                           +mffsvelint(0,0) * vderxy(2,0)
                                                                           +mffsvelint(1,0) * vderxy(2,1)
                                                                           +mffsvelint(2,0) * vderxy(2,2)));

    /* cross-stress term on right hand side */
    /* additional terms conservative form */
    /*
                 /                                         \
                |                                           |
         - dens | ( du (nabla o u) + u (nabla o du ) ,  v   |
                |                                           |
                 \                                         /
    */

    if (mfsconservative) // no conservative form for combustion
    {
      assembler.template Vector<XFEM::PHYSICS::Velx>(shpvel.d0,-altrhsdensfac * (mffsvelint(0,0) * vdiv + velint(0,0) * mffsvdiv));
      assembler.template Vector<XFEM::PHYSICS::Vely>(shpvel.d0,-altrhsdensfac * (mffsvelint(1,0) * vdiv + velint(1,0) * mffsvdiv));
      assembler.template Vector<XFEM::PHYSICS::Velz>(shpvel.d0,-altrhsdensfac * (mffsvelint(2,0) * vdiv + velint(2,0) * mffsvdiv));
    }

    /* reynolds-stress term on right hand side */
    /*
             /                       \
            |                         |
     - dens | ( du o nabla du) ,  v   |
            |                         |
             \                       /
    */

    assembler.template Vector<XFEM::PHYSICS::Velx>(shpvel.d0,-altrhsdensfac * (mffsvelint(0,0) * mffsvderxy(0,0)
                                                                           +mffsvelint(1,0) * mffsvderxy(0,1)
                                                                           +mffsvelint(2,0) * mffsvderxy(0,2)));
    assembler.template Vector<XFEM::PHYSICS::Vely>(shpvel.d0,-altrhsdensfac * (mffsvelint(0,0) * mffsvderxy(1,0)
                                                                           +mffsvelint(1,0) * mffsvderxy(1,1)
                                                                           +mffsvelint(2,0) * mffsvderxy(1,2)));
    assembler.template Vector<XFEM::PHYSICS::Velz>(shpvel.d0,-altrhsdensfac * (mffsvelint(0,0) * mffsvderxy(2,0)
                                                                           +mffsvelint(1,0) * mffsvderxy(2,1)
                                                                           +mffsvelint(2,0) * mffsvderxy(2,2)));

    /* Reynolds-stress term on right hand side */
    /* additional terms conservative form */
    /*
             /                       \
            |                         |
     - dens |   du (nabla o du),  v   |
            |                         |
             \                       /
    */

    if (mfsconservative) // no conservative form for combustion
    {
      assembler.template Vector<XFEM::PHYSICS::Velx>(shpvel.d0,-altrhsdensfac * (mffsvelint(0,0) * mffsvdiv));
      assembler.template Vector<XFEM::PHYSICS::Vely>(shpvel.d0,-altrhsdensfac * (mffsvelint(1,0) * mffsvdiv));
      assembler.template Vector<XFEM::PHYSICS::Velz>(shpvel.d0,-altrhsdensfac * (mffsvelint(2,0) * mffsvdiv));
    }
  }

  return;
}


/*----------------------------------------------------------------------*
 |  calculation of fine-scale subgrid viscosity                vg 09/09 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType DISTYPE>
  double CalcFineScaleSubgrVisc(
  const LINALG::Matrix<3,3>&      vderxy,   // velocity derivative at GP at n+1
  const LINALG::Matrix<3,3>&      fsvderxy, // fine-scale velocity derivative at GP at n+1
  const double                    vol,
  const double                    Cs,
  const INPAR::FLUID::FineSubgridVisc fssgrv
  )
{
  // kinematic fine-scale viscosity
  double fssgvisc=0.0;

  // get characteristic element length for Smagorinsky model for 2D and 3D
  // NOTICE: here really the element volume is taken (not the volume of the integration cell!)
  // 3D: hk = V^1/3
  double hk= pow(vol,(1.0/3.0));

  if (fssgrv==INPAR::FLUID::smagorinsky_all)
  {
    // ALL-SCALE SMAGORINSKY MODEL
    // ---------------------------
    //                                   +-                                 -+ 1
    //                               2   |          / h \           / h \    | -
    //    kin-visc          = (C_S*h)  * | 2 * eps | u   |   * eps | u   |   | 2
    //            turbulent              |          \   / ij        \   / ij |
    //                                   +-                                 -+
    //                                   |                                   |
    //                                   +-----------------------------------+
    //                                           'resolved' rate of strain
    //

    // compute (all-scale) rate of strain
    double rateofstrain = COMBUST::GetStrainRate<DISTYPE>(vderxy);
    fssgvisc = Cs * Cs * hk * hk * rateofstrain;

  }
  else if (fssgrv==INPAR::FLUID::smagorinsky_small)
  {
    // FINE-SCALE SMAGORINSKY MODEL
    // ----------------------------
    //                                   +-                                 -+ 1
    //                               2   |          /    \          /   \    | -
    //    kin-visc          = (C_S*h)  * | 2 * eps | fsu |   * eps | fsu |   | 2
    //            turbulent              |          \   / ij        \   / ij |
    //                                   +-                                 -+
    //                                   |                                   |
    //                                   +-----------------------------------+
    //                                            'resolved' rate of strain
    //

    // fine-scale rate of strain
    double fsrateofstrain=COMBUST::GetStrainRate<DISTYPE>(fsvderxy);
    fssgvisc =  Cs * Cs * hk * hk * fsrateofstrain;
  }
  else
    dserror("Unknown approach for avm3!");

  return fssgvisc;
}


/*----------------------------------------------------------------------*
 |  compute multifractal subgrid scales parameters    rasthofer 04/2011 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType DISTYPE>
void PrepareMultifractalSubgrScales(
  double&                           B_mfs,
  const LINALG::Matrix<3,1>&        velint,
  const LINALG::Matrix<3,1>&        fsvelint,
  const LINALG::Matrix<3,3>&        vderxy,
  const double                      vol,
  const double                      Csgs,
  const double                      alpha,
  const bool                        calcN,
  const double                      N,
  const INPAR::FLUID::RefLength     reflength,
  const INPAR::FLUID::RefVelocity   refvel,
  const double                      cnu,
  const bool                        nwl,
  const double                      dens,
  const double                      visc
)
{
  // cascade steps
  double Nvel = 0.0;
  // corrected parameter
  double Csgs_ai = Csgs;

  // potential calculation of Re to determine N
  double Re_ele = -1.0;
  // characteristic element length
  double hk = 1.0e+10;
  double strainnorm = 0.0;
  // ratio of viscous scale to element length
  double scale_ratio = 0.0;

  // get norm
  const double vel_norm = velint.Norm2();
  const double fsvel_norm = fsvelint.Norm2();

  // do we have a fixed parameter N
  if (not calcN)
    Nvel = N;
  else //no, so we calculate N from Re
  {
    // calculate characteristic element length
    // cf. stabilization parameters
    switch (reflength)
    {
      case INPAR::FLUID::streamlength:
      {
        dserror("Reflength streamlength not yet supported for xfem!");
#if 0
        // a) streamlength due to Tezduyar et al. (1992)
        // normed velocity vector
        LINALG::Matrix<3,1> velino(true);
        if (vel_norm>=1e-6) velino.Update(1.0/vel_norm,velint);
        else
        {
          velino.Clear();
          velino(0,0) = 1.0;
        }
        LINALG::Matrix<nen,1> tmp;
        tmp.MultiplyTN(derxy,velino);
        const double val = tmp.Norm1();
        hk = 2.0/val;

        break;
#endif
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
        hk = std::pow(vol,(1.0/3.0));
        break;
      }
      case INPAR::FLUID::metric_tensor:
      {
        dserror("Reflength metric_tensor not yet supported for xfem!");
#if 0
        /*          +-           -+   +-           -+   +-           -+
                    |             |   |             |   |             |
                    |  dr    dr   |   |  ds    ds   |   |  dt    dt   |
              G   = |  --- * ---  | + |  --- * ---  | + |  --- * ---  |
               ij   |  dx    dx   |   |  dx    dx   |   |  dx    dx   |
                    |    i     j  |   |    i     j  |   |    i     j  |
                    +-           -+   +-           -+   +-           -+
        */
        LINALG::Matrix<3,3> G(true);

        for (int nn=0;nn<3;++nn)
        {
          for (int rr=0;rr<3;++rr)
          {
            G(nn,rr) = xji(nn,0)*xji(rr,0);
            for (int mm=1;mm<3;++mm)
            {
              G(nn,rr) += xji(nn,mm)*xji(rr,mm);
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
        for (int nn=0;nn<3;++nn)
        {
          for (int rr=0;rr<3;++rr)
          {
            normG+=G(nn,rr)*G(nn,rr);
          }
        }
        hk = std::pow(normG,-0.25);
#endif
        break;
      }
      case INPAR::FLUID::gradient_based:
      {
        dserror("Reflength gradient_based not yet supported for xfem!");
#if 0
        LINALG::Matrix<3,1> normed_velgrad;

        for (int rr=0;rr<3;++rr)
        {
          double val = 0.0;
          for (int idim = 0; idim < 3; idim ++)
            val += vderxy(idim,rr)*vderxy(idim,rr);

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
        for (int rr=0;rr<3;++rr) /* loop element nodes */
        {
          double loc = 0.0;
          for (int idim = 0; idim < 3; idim ++)
            loc += normed_velgrad(idim)*derxy(idim,rr);

          val += std::abs(loc);

          //val += abs( normed_velgrad(0)*derxy_(0,rr)
          //            +normed_velgrad(1)*derxy_(1,rr)
          //            +normed_velgrad(2)*derxy_(2,rr));
        } /* end of loop over element nodes */

        hk = 2.0/val;
#endif
        break;
      }
      default:
        dserror("Unknown length");
      }

      if (hk == 1.0e+10)
        dserror("Something went wrong!");

      switch (refvel){
      case INPAR::FLUID::resolved:
      {
        Re_ele = vel_norm * hk * dens / visc;
        break;
      }
      case INPAR::FLUID::fine_scale:
      {
        Re_ele = fsvel_norm * hk * dens / visc;
        break;
      }
      case INPAR::FLUID::strainrate:
      {
        strainnorm = COMBUST::GetStrainRate<DISTYPE>(vderxy);
        strainnorm /= sqrt(2.0); //cf. Burton & Dahm 2005
        Re_ele = strainnorm * hk * hk * dens / visc;
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
    scale_ratio = cnu * pow(Re_ele,0.75);
    // scale_ratio < 1.0 leads to N < 0
    // therefore, we clip once more
    if (scale_ratio < 1.0)
      scale_ratio = 1.0;

    //         |   Delta     |
    //  N =log | ----------- |
    //        2|  lambda_nu  |
    Nvel = log(scale_ratio)/log(2.0);
    if (Nvel < 0.0)
      dserror("Something went wrong when calculating N!");
  }

  // calculate near-wall correction
  if (nwl)
  {
    // if not yet calculated, estimate norm of strain rate
    if ((not calcN) or (refvel != INPAR::FLUID::strainrate))
    {
      strainnorm = COMBUST::GetStrainRate<DISTYPE>(vderxy);
      strainnorm /= sqrt(2.0); //cf. Burton & Dahm 2005
    }
    // and reference length
    if (not calcN) dserror("hk not yet calculated"); // solution see scatra

    // get Re from strain rate
    double Re_ele_str = strainnorm * hk * hk * dens / visc;
    if (Re_ele_str < 0.0)
      dserror("Something went wrong!");
    // ensure positive values
    if (Re_ele_str < 1.0)
      Re_ele_str = 1.0;

    // calculate corrected Csgs_ai
    //           -3/16
    //  *(1 - (Re)   )
    //
    Csgs_ai *= (1.0-pow(Re_ele_str,-3.0/16.0));
  }

  // compute coefficient B

  //
  //          |       1              |
  //  kappa = | -------------------- |
  //          |  1 - alpha ^ (-4/3)  |
  //
  double kappa = 1.0/(1.0-pow(alpha,-4.0/3.0));

  //                  1                                    1
  //                  2                 |                 |2
  //  B = Csgs_ai * kappa * 2 ^ (-2*N/3) * | 2 ^ (4*N/3) - 1 |
  //                                    |                 |
  //

  B_mfs = Csgs_ai *sqrt(kappa) * pow(2.0,-2.0*Nvel/3.0) * sqrt((pow(2.0,4.0*Nvel/3.0)-1.0));

    //if (eid_ == 100)
//    {
//      std::cout << "B  " << std::setprecision(10) << B_mfs << std::endl;
//      std::cout << "CsgsB  " << std::setprecision(10) << Csgs_ai << std::endl;
//    }

  return;
}

}
#endif /* COMBUST3_SYSMAT_TURBULENCE_H_ */
