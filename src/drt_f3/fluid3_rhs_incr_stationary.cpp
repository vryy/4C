// This file is currently not used by any of the fluid implementations
// Fluid3Impl, Fluid3GenalphaResVMM or Fluid3Stationary. They use
// rearranged versions of this code (slightly adopted to the different
// approaches) --- they are expected to be faster and more readable.
// 
// Axel is still referencing these loops in his xfem routines, so they
// should not be deleted until he declares them obsolete...

for (int vi=0; vi<iel; ++vi)
{
  /* Konvektionsterm */
  eforce_(vi*4)     += -(fac*(velint_(0)*conv_r_(0, 0, vi) + velint_(1)*conv_r_(0, 1, vi) + velint_(2)*conv_r_(0, 2, vi))) ;
  eforce_(vi*4 + 1) += -(fac*(velint_(0)*conv_r_(1, 0, vi) + velint_(1)*conv_r_(1, 1, vi) + velint_(2)*conv_r_(1, 2, vi))) ;
  eforce_(vi*4 + 2) += -(fac*(velint_(0)*conv_r_(2, 0, vi) + velint_(1)*conv_r_(2, 1, vi) + velint_(2)*conv_r_(2, 2, vi))) ;

  /* Viskositätsterm */
  eforce_(vi*4)     += -(nu_*fac*(2.0*derxyz_(0, vi)*vderxyz_(0, 0) + derxyz_(1, vi)*vderxyz_(0, 1) + derxyz_(1, vi)*vderxyz_(1, 0) + derxyz_(2, vi)*vderxyz_(0, 2) + derxyz_(2, vi)*vderxyz_(2, 0))) ;
  eforce_(vi*4 + 1) += -(nu_*fac*(derxyz_(0, vi)*vderxyz_(0, 1) + derxyz_(0, vi)*vderxyz_(1, 0) + 2.0*derxyz_(1, vi)*vderxyz_(1, 1) + derxyz_(2, vi)*vderxyz_(1, 2) + derxyz_(2, vi)*vderxyz_(2, 1))) ;
  eforce_(vi*4 + 2) += -(nu_*fac*(derxyz_(0, vi)*vderxyz_(0, 2) + derxyz_(0, vi)*vderxyz_(2, 0) + derxyz_(1, vi)*vderxyz_(1, 2) + derxyz_(1, vi)*vderxyz_(2, 1) + 2.0*derxyz_(2, vi)*vderxyz_(2, 2))) ;

  /* Druckterm */
  eforce_(vi*4)     += press*fac*derxyz_(0, vi) ;
  eforce_(vi*4 + 1) += press*fac*derxyz_(1, vi) ;
  eforce_(vi*4 + 2) += press*fac*derxyz_(2, vi) ;
  
  /* Divergenzfreiheit */
  eforce_(vi*4 + 3) += -(fac*(conv_r_(0, 0, vi) + conv_r_(1, 1, vi) + conv_r_(2, 2, vi))) ;

  /* konvektive Stabilisierung: */
  /* Stabilisierung der Konvektion ( L_conv_u) */
  eforce_(vi*4)     += -(tau_M*conv_c_(vi)*conv_old_(0)) ;
  eforce_(vi*4 + 1) += -(tau_M*conv_c_(vi)*conv_old_(1)) ;
  eforce_(vi*4 + 2) += -(tau_M*conv_c_(vi)*conv_old_(2)) ;
  
    /* Stabilisierung der Konvektion (-L_visc_u) */
  eforce_(vi*4)     += 2.0*nu_*tau_M*conv_c_(vi)*visc_old_(0) ;
  eforce_(vi*4 + 1) += 2.0*nu_*tau_M*conv_c_(vi)*visc_old_(1) ;
  eforce_(vi*4 + 2) += 2.0*nu_*tau_M*conv_c_(vi)*visc_old_(2) ;

  /* Stabilisierung der Konvektion ( L_pres_p) */
  eforce_(vi*4)     += -(tau_M*conv_c_(vi)*gradp_(0)) ;
  eforce_(vi*4 + 1) += -(tau_M*conv_c_(vi)*gradp_(1)) ;
  eforce_(vi*4 + 2) += -(tau_M*conv_c_(vi)*gradp_(2)) ;
  
  /* Druck-Stabilisierung: */
  /* Stabilisierung des Drucks ( L_conv_u) */
  eforce_(vi*4 + 3) += -(tau_Mp*(conv_old_(0)*derxyz_(0, vi) + conv_old_(1)*derxyz_(1, vi) + conv_old_(2)*derxyz_(2, vi))) ;

  /* Stabilisierung des Drucks (-L_visc_u) */
  eforce_(vi*4 + 3) += 2.0*nu_*tau_Mp*(visc_old_(0)*derxyz_(0, vi) + visc_old_(1)*derxyz_(1, vi) + visc_old_(2)*derxyz_(2, vi)) ;

  /* Stabilisierung des Drucks ( L_pres_p) */
  eforce_(vi*4 + 3) += -(tau_Mp*(gradp_(0)*derxyz_(0, vi) + gradp_(1)*derxyz_(1, vi) + gradp_(2)*derxyz_(2, vi))) ;

  /* Kontinuitätsstabilisierung */
  eforce_(vi*4)     += -(tau_C*derxyz_(0, vi)*(vderxyz_(0, 0) + vderxyz_(1, 1) + vderxyz_(2, 2))) ;
  eforce_(vi*4 + 1) += -(tau_C*derxyz_(1, vi)*(vderxyz_(0, 0) + vderxyz_(1, 1) + vderxyz_(2, 2))) ;
  eforce_(vi*4 + 2) += -(tau_C*derxyz_(2, vi)*(vderxyz_(0, 0) + vderxyz_(1, 1) + vderxyz_(2, 2))) ;
}
