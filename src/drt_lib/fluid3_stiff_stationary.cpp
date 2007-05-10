for (int vi=0; vi<iel; ++vi)    
{
  for (int ui=0; ui<iel; ++ui)
  {
    /* Konvektionsterm */
    estif_(vi*4, ui*4)         += fac*funct_(vi)*(conv_c_(ui) + conv_r_(0, 0, ui)) ;
    estif_(vi*4, ui*4 + 1)     += fac*funct_(vi)*conv_r_(0, 1, ui) ;
    estif_(vi*4, ui*4 + 2)     += fac*funct_(vi)*conv_r_(0, 2, ui) ;
    estif_(vi*4 + 1, ui*4)     += fac*funct_(vi)*conv_r_(1, 0, ui) ;
    estif_(vi*4 + 1, ui*4 + 1) += fac*funct_(vi)*(conv_c_(ui) + conv_r_(1, 1, ui)) ;
    estif_(vi*4 + 1, ui*4 + 2) += fac*funct_(vi)*conv_r_(1, 2, ui) ;
    estif_(vi*4 + 2, ui*4)     += fac*funct_(vi)*conv_r_(2, 0, ui) ;
    estif_(vi*4 + 2, ui*4 + 1) += fac*funct_(vi)*conv_r_(2, 1, ui) ;
    estif_(vi*4 + 2, ui*4 + 2) += fac*funct_(vi)*(conv_c_(ui) + conv_r_(2, 2, ui)) ;

    /* Viskositätsterm */
    estif_(vi*4, ui*4)         += nu_*fac*(2.0*derxyz_(0, ui)*derxyz_(0, vi) + derxyz_(1, ui)*derxyz_(1, vi) + derxyz_(2, ui)*derxyz_(2, vi)) ;
    estif_(vi*4, ui*4 + 1)     += nu_*fac*derxyz_(0, ui)*derxyz_(1, vi) ;
    estif_(vi*4, ui*4 + 2)     += nu_*fac*derxyz_(0, ui)*derxyz_(2, vi) ;
    estif_(vi*4 + 1, ui*4)     += nu_*fac*derxyz_(0, vi)*derxyz_(1, ui) ;
    estif_(vi*4 + 1, ui*4 + 1) += nu_*fac*(derxyz_(0, ui)*derxyz_(0, vi) + 2.0*derxyz_(1, ui)*derxyz_(1, vi) + derxyz_(2, ui)*derxyz_(2, vi)) ;
    estif_(vi*4 + 1, ui*4 + 2) += nu_*fac*derxyz_(1, ui)*derxyz_(2, vi) ;
    estif_(vi*4 + 2, ui*4)     += nu_*fac*derxyz_(0, vi)*derxyz_(2, ui) ;
    estif_(vi*4 + 2, ui*4 + 1) += nu_*fac*derxyz_(1, vi)*derxyz_(2, ui) ;
    estif_(vi*4 + 2, ui*4 + 2) += nu_*fac*(derxyz_(0, ui)*derxyz_(0, vi) + derxyz_(1, ui)*derxyz_(1, vi) + 2.0*derxyz_(2, ui)*derxyz_(2, vi)) ;
    
    /* Druckterm */
    estif_(vi*4, ui*4 + 3)     += -(fac*funct_(ui)*derxyz_(0, vi)) ;
    estif_(vi*4 + 1, ui*4 + 3) += -(fac*funct_(ui)*derxyz_(1, vi)) ;
    estif_(vi*4 + 2, ui*4 + 3) += -(fac*funct_(ui)*derxyz_(2, vi)) ;
    
    /* Divergenzfreiheit */
    estif_(vi*4 + 3, ui*4)     += fac*funct_(vi)*derxyz_(0, ui) ;
    estif_(vi*4 + 3, ui*4 + 1) += fac*funct_(vi)*derxyz_(1, ui) ;
    estif_(vi*4 + 3, ui*4 + 2) += fac*funct_(vi)*derxyz_(2, ui) ;
    
    /* konvektive Stabilisierung: */
    /* Stabilisierung der Konvektion ( L_conv_u) */
    estif_(vi*4, ui*4)         += tau_M*(conv_c_(ui)*conv_c_(vi) + conv_c_(vi)*conv_r_(0, 0, ui) + velint_(0)*derxyz_(0, vi)*conv_r_(0, 0, ui) + velint_(1)*derxyz_(0, vi)*conv_r_(0, 1, ui) + velint_(2)*derxyz_(0, vi)*conv_r_(0, 2, ui)) ;
    estif_(vi*4, ui*4 + 1)     += tau_M*(conv_c_(vi)*conv_r_(0, 1, ui) + velint_(0)*derxyz_(1, vi)*conv_r_(0, 0, ui) + velint_(1)*derxyz_(1, vi)*conv_r_(0, 1, ui) + velint_(2)*derxyz_(1, vi)*conv_r_(0, 2, ui)) ;
    estif_(vi*4, ui*4 + 2)     += tau_M*(conv_c_(vi)*conv_r_(0, 2, ui) + velint_(0)*derxyz_(2, vi)*conv_r_(0, 0, ui) + velint_(1)*derxyz_(2, vi)*conv_r_(0, 1, ui) + velint_(2)*derxyz_(2, vi)*conv_r_(0, 2, ui)) ;
    estif_(vi*4 + 1, ui*4)     += tau_M*(conv_c_(vi)*conv_r_(1, 0, ui) + velint_(0)*derxyz_(0, vi)*conv_r_(1, 0, ui) + velint_(1)*derxyz_(0, vi)*conv_r_(1, 1, ui) + velint_(2)*derxyz_(0, vi)*conv_r_(1, 2, ui)) ;
    estif_(vi*4 + 1, ui*4 + 1) += tau_M*(conv_c_(ui)*conv_c_(vi) + conv_c_(vi)*conv_r_(1, 1, ui) + velint_(0)*derxyz_(1, vi)*conv_r_(1, 0, ui) + velint_(1)*derxyz_(1, vi)*conv_r_(1, 1, ui) + velint_(2)*derxyz_(1, vi)*conv_r_(1, 2, ui)) ;
    estif_(vi*4 + 1, ui*4 + 2) += tau_M*(conv_c_(vi)*conv_r_(1, 2, ui) + velint_(0)*derxyz_(2, vi)*conv_r_(1, 0, ui) + velint_(1)*derxyz_(2, vi)*conv_r_(1, 1, ui) + velint_(2)*derxyz_(2, vi)*conv_r_(1, 2, ui)) ;
    estif_(vi*4 + 2, ui*4)     += tau_M*(conv_c_(vi)*conv_r_(2, 0, ui) + velint_(0)*derxyz_(0, vi)*conv_r_(2, 0, ui) + velint_(1)*derxyz_(0, vi)*conv_r_(2, 1, ui) + velint_(2)*derxyz_(0, vi)*conv_r_(2, 2, ui)) ;
    estif_(vi*4 + 2, ui*4 + 1) += tau_M*(conv_c_(vi)*conv_r_(2, 1, ui) + velint_(0)*derxyz_(1, vi)*conv_r_(2, 0, ui) + velint_(1)*derxyz_(1, vi)*conv_r_(2, 1, ui) + velint_(2)*derxyz_(1, vi)*conv_r_(2, 2, ui)) ;
    estif_(vi*4 + 2, ui*4 + 2) += tau_M*(conv_c_(ui)*conv_c_(vi) + conv_c_(vi)*conv_r_(2, 2, ui) + velint_(0)*derxyz_(2, vi)*conv_r_(2, 0, ui) + velint_(1)*derxyz_(2, vi)*conv_r_(2, 1, ui) + velint_(2)*derxyz_(2, vi)*conv_r_(2, 2, ui)) ;

    /* Stabilisierung der Konvektion (-L_visc_u) */
    estif_(vi*4, ui*4)         += -2.0*nu_*tau_M*(-(conv_c_(vi)*viscs2_(0, 0, ui)) + funct_(ui)*visc_old_(0)*derxyz_(0, vi)) ;
    estif_(vi*4, ui*4 + 1)     += -2.0*nu_*tau_M*(-(conv_c_(vi)*viscs2_(0, 1, ui)) + funct_(ui)*visc_old_(0)*derxyz_(1, vi)) ;
    estif_(vi*4, ui*4 + 2)     += -2.0*nu_*tau_M*(-(conv_c_(vi)*viscs2_(0, 2, ui)) + funct_(ui)*visc_old_(0)*derxyz_(2, vi)) ;
    estif_(vi*4 + 1, ui*4)     += -2.0*nu_*tau_M*(-(conv_c_(vi)*viscs2_(0, 1, ui)) + funct_(ui)*visc_old_(1)*derxyz_(0, vi)) ;
    estif_(vi*4 + 1, ui*4 + 1) += -2.0*nu_*tau_M*(-(conv_c_(vi)*viscs2_(1, 1, ui)) + funct_(ui)*visc_old_(1)*derxyz_(1, vi)) ;
    estif_(vi*4 + 1, ui*4 + 2) += -2.0*nu_*tau_M*(-(conv_c_(vi)*viscs2_(1, 2, ui)) + funct_(ui)*visc_old_(1)*derxyz_(2, vi)) ;
    estif_(vi*4 + 2, ui*4)     += -2.0*nu_*tau_M*(-(conv_c_(vi)*viscs2_(0, 2, ui)) + funct_(ui)*visc_old_(2)*derxyz_(0, vi)) ;
    estif_(vi*4 + 2, ui*4 + 1) += -2.0*nu_*tau_M*(-(conv_c_(vi)*viscs2_(1, 2, ui)) + funct_(ui)*visc_old_(2)*derxyz_(1, vi)) ;
    estif_(vi*4 + 2, ui*4 + 2) += -2.0*nu_*tau_M*(-(conv_c_(vi)*viscs2_(2, 2, ui)) + funct_(ui)*visc_old_(2)*derxyz_(2, vi)) ;

    /* Stabilisierung der Konvektion ( L_pres_p) */
    estif_(vi*4, ui*4)         += tau_M*funct_(ui)*gradp_(0)*derxyz_(0, vi) ;
    estif_(vi*4, ui*4 + 1)     += tau_M*funct_(ui)*gradp_(0)*derxyz_(1, vi) ;
    estif_(vi*4, ui*4 + 2)     += tau_M*funct_(ui)*gradp_(0)*derxyz_(2, vi) ;
    estif_(vi*4, ui*4 + 3)     += tau_M*conv_c_(vi)*derxyz_(0, ui) ;
    estif_(vi*4 + 1, ui*4)     += tau_M*funct_(ui)*gradp_(1)*derxyz_(0, vi) ;
    estif_(vi*4 + 1, ui*4 + 1) += tau_M*funct_(ui)*gradp_(1)*derxyz_(1, vi) ;
    estif_(vi*4 + 1, ui*4 + 2) += tau_M*funct_(ui)*gradp_(1)*derxyz_(2, vi) ;
    estif_(vi*4 + 1, ui*4 + 3) += tau_M*conv_c_(vi)*derxyz_(1, ui) ;
    estif_(vi*4 + 2, ui*4)     += tau_M*funct_(ui)*gradp_(2)*derxyz_(0, vi) ;
    estif_(vi*4 + 2, ui*4 + 1) += tau_M*funct_(ui)*gradp_(2)*derxyz_(1, vi) ;
    estif_(vi*4 + 2, ui*4 + 2) += tau_M*funct_(ui)*gradp_(2)*derxyz_(2, vi) ;
    estif_(vi*4 + 2, ui*4 + 3) += tau_M*conv_c_(vi)*derxyz_(2, ui) ;

    /* Druck-Stabilisierung: */
    /* Stabilisierung des Drucks ( L_conv_u) */
    estif_(vi*4 + 3, ui*4)     += tau_Mp*(conv_c_(ui)*derxyz_(0, vi) + derxyz_(0, vi)*conv_r_(0, 0, ui) + derxyz_(1, vi)*conv_r_(1, 0, ui) + derxyz_(2, vi)*conv_r_(2, 0, ui)) ;
    estif_(vi*4 + 3, ui*4 + 1) += tau_Mp*(conv_c_(ui)*derxyz_(1, vi) + derxyz_(0, vi)*conv_r_(0, 1, ui) + derxyz_(1, vi)*conv_r_(1, 1, ui) + derxyz_(2, vi)*conv_r_(2, 1, ui)) ;
    estif_(vi*4 + 3, ui*4 + 2) += tau_Mp*(conv_c_(ui)*derxyz_(2, vi) + derxyz_(0, vi)*conv_r_(0, 2, ui) + derxyz_(1, vi)*conv_r_(1, 2, ui) + derxyz_(2, vi)*conv_r_(2, 2, ui)) ;

    /* Stabilisierung des Drucks (-L_visc_u) */
    estif_(vi*4 + 3, ui*4)     += 2.0*nu_*tau_Mp*(derxyz_(0, vi)*viscs2_(0, 0, ui) + derxyz_(1, vi)*viscs2_(0, 1, ui) + derxyz_(2, vi)*viscs2_(0, 2, ui)) ;
    estif_(vi*4 + 3, ui*4 + 1) += 2.0*nu_*tau_Mp*(derxyz_(0, vi)*viscs2_(0, 1, ui) + derxyz_(1, vi)*viscs2_(1, 1, ui) + derxyz_(2, vi)*viscs2_(1, 2, ui)) ;
    estif_(vi*4 + 3, ui*4 + 2) += 2.0*nu_*tau_Mp*(derxyz_(0, vi)*viscs2_(0, 2, ui) + derxyz_(1, vi)*viscs2_(1, 2, ui) + derxyz_(2, vi)*viscs2_(2, 2, ui)) ;

    /* Stabilisierung des Drucks ( L_pres_p) */
    estif_(vi*4 + 3, ui*4 + 3) += tau_Mp*(derxyz_(0, ui)*derxyz_(0, vi) + derxyz_(1, ui)*derxyz_(1, vi) + derxyz_(2, ui)*derxyz_(2, vi)) ;

    /* Kontinuitätsstabilisierung */
    estif_(vi*4, ui*4)         += tau_C*derxyz_(0, ui)*derxyz_(0, vi) ;
    estif_(vi*4, ui*4 + 1)     += tau_C*derxyz_(0, vi)*derxyz_(1, ui) ;
    estif_(vi*4, ui*4 + 2)     += tau_C*derxyz_(0, vi)*derxyz_(2, ui) ;
    estif_(vi*4 + 1, ui*4)     += tau_C*derxyz_(0, ui)*derxyz_(1, vi) ;
    estif_(vi*4 + 1, ui*4 + 1) += tau_C*derxyz_(1, ui)*derxyz_(1, vi) ;
    estif_(vi*4 + 1, ui*4 + 2) += tau_C*derxyz_(1, vi)*derxyz_(2, ui) ;
    estif_(vi*4 + 2, ui*4)     += tau_C*derxyz_(0, ui)*derxyz_(2, vi) ;
    estif_(vi*4 + 2, ui*4 + 1) += tau_C*derxyz_(1, ui)*derxyz_(2, vi) ;
    estif_(vi*4 + 2, ui*4 + 2) += tau_C*derxyz_(2, ui)*derxyz_(2, vi) ;
  }
}
