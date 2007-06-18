for (vi=0; vi<iel; ++vi)
{
  for (ui=0; ui<iel; ++ui)
  {
    /* Konvektionsterm */
    estif_(vi*4, ui*4)         += timefacfac*funct_(vi)*(conv_c_(ui) + conv_g_(ui) + conv_r_(0, 0, ui)) ;
    estif_(vi*4, ui*4 + 1)     += timefacfac*funct_(vi)*conv_r_(0, 1, ui) ;
    estif_(vi*4, ui*4 + 2)     += timefacfac*funct_(vi)*conv_r_(0, 2, ui) ;
    estif_(vi*4 + 1, ui*4)     += timefacfac*funct_(vi)*conv_r_(1, 0, ui) ;
    estif_(vi*4 + 1, ui*4 + 1) += timefacfac*funct_(vi)*(conv_c_(ui) + conv_g_(ui) + conv_r_(1, 1, ui)) ;
    estif_(vi*4 + 1, ui*4 + 2) += timefacfac*funct_(vi)*conv_r_(1, 2, ui) ;
    estif_(vi*4 + 2, ui*4)     += timefacfac*funct_(vi)*conv_r_(2, 0, ui) ;
    estif_(vi*4 + 2, ui*4 + 1) += timefacfac*funct_(vi)*conv_r_(2, 1, ui) ;
    estif_(vi*4 + 2, ui*4 + 2) += timefacfac*funct_(vi)*(conv_c_(ui) + conv_g_(ui) + conv_r_(2, 2, ui)) ;

    /* Stabilisierung der Konvektion ( L_conv_u) */
    estif_(vi*4, ui*4)         += ttimetauM*(conv_c_(ui)*conv_c_(vi) + conv_c_(ui)*conv_g_(vi) + conv_c_(vi)*conv_g_(ui) + conv_g_(ui)*conv_g_(vi) + conv_c_(vi)*conv_r_(0, 0, ui) + conv_g_(vi)*conv_r_(0, 0, ui) + velint_(0)*derxyz_(0, vi)*conv_r_(0, 0, ui) + velint_(1)*derxyz_(0, vi)*conv_r_(0, 1, ui) + velint_(2)*derxyz_(0, vi)*conv_r_(0, 2, ui) - gridvint_(0)*derxyz_(0, vi)*conv_r_(0, 0, ui) - gridvint_(1)*derxyz_(0, vi)*conv_r_(0, 1, ui) - gridvint_(2)*derxyz_(0, vi)*conv_r_(0, 2, ui)) ;
    estif_(vi*4, ui*4 + 1)     += ttimetauM*(conv_c_(vi)*conv_r_(0, 1, ui) + conv_g_(vi)*conv_r_(0, 1, ui) + velint_(0)*derxyz_(1, vi)*conv_r_(0, 0, ui) + velint_(1)*derxyz_(1, vi)*conv_r_(0, 1, ui) + velint_(2)*derxyz_(1, vi)*conv_r_(0, 2, ui) - gridvint_(0)*derxyz_(1, vi)*conv_r_(0, 0, ui) - gridvint_(1)*derxyz_(1, vi)*conv_r_(0, 1, ui) - gridvint_(2)*derxyz_(1, vi)*conv_r_(0, 2, ui)) ;
    estif_(vi*4, ui*4 + 2)     += ttimetauM*(conv_c_(vi)*conv_r_(0, 2, ui) + conv_g_(vi)*conv_r_(0, 2, ui) + velint_(0)*derxyz_(2, vi)*conv_r_(0, 0, ui) + velint_(1)*derxyz_(2, vi)*conv_r_(0, 1, ui) + velint_(2)*derxyz_(2, vi)*conv_r_(0, 2, ui) - gridvint_(0)*derxyz_(2, vi)*conv_r_(0, 0, ui) - gridvint_(1)*derxyz_(2, vi)*conv_r_(0, 1, ui) - gridvint_(2)*derxyz_(2, vi)*conv_r_(0, 2, ui)) ;
    estif_(vi*4 + 1, ui*4)     += ttimetauM*(conv_c_(vi)*conv_r_(1, 0, ui) + conv_g_(vi)*conv_r_(1, 0, ui) + velint_(0)*derxyz_(0, vi)*conv_r_(1, 0, ui) + velint_(1)*derxyz_(0, vi)*conv_r_(1, 1, ui) + velint_(2)*derxyz_(0, vi)*conv_r_(1, 2, ui) - gridvint_(0)*derxyz_(0, vi)*conv_r_(1, 0, ui) - gridvint_(1)*derxyz_(0, vi)*conv_r_(1, 1, ui) - gridvint_(2)*derxyz_(0, vi)*conv_r_(1, 2, ui)) ;
    estif_(vi*4 + 1, ui*4 + 1) += ttimetauM*(conv_c_(ui)*conv_c_(vi) + conv_c_(ui)*conv_g_(vi) + conv_c_(vi)*conv_g_(ui) + conv_g_(ui)*conv_g_(vi) + conv_c_(vi)*conv_r_(1, 1, ui) + conv_g_(vi)*conv_r_(1, 1, ui) + velint_(0)*derxyz_(1, vi)*conv_r_(1, 0, ui) + velint_(1)*derxyz_(1, vi)*conv_r_(1, 1, ui) + velint_(2)*derxyz_(1, vi)*conv_r_(1, 2, ui) - gridvint_(0)*derxyz_(1, vi)*conv_r_(1, 0, ui) - gridvint_(1)*derxyz_(1, vi)*conv_r_(1, 1, ui) - gridvint_(2)*derxyz_(1, vi)*conv_r_(1, 2, ui)) ;
    estif_(vi*4 + 1, ui*4 + 2) += ttimetauM*(conv_c_(vi)*conv_r_(1, 2, ui) + conv_g_(vi)*conv_r_(1, 2, ui) + velint_(0)*derxyz_(2, vi)*conv_r_(1, 0, ui) + velint_(1)*derxyz_(2, vi)*conv_r_(1, 1, ui) + velint_(2)*derxyz_(2, vi)*conv_r_(1, 2, ui) - gridvint_(0)*derxyz_(2, vi)*conv_r_(1, 0, ui) - gridvint_(1)*derxyz_(2, vi)*conv_r_(1, 1, ui) - gridvint_(2)*derxyz_(2, vi)*conv_r_(1, 2, ui)) ;
    estif_(vi*4 + 2, ui*4)     += ttimetauM*(conv_c_(vi)*conv_r_(2, 0, ui) + conv_g_(vi)*conv_r_(2, 0, ui) + velint_(0)*derxyz_(0, vi)*conv_r_(2, 0, ui) + velint_(1)*derxyz_(0, vi)*conv_r_(2, 1, ui) + velint_(2)*derxyz_(0, vi)*conv_r_(2, 2, ui) - gridvint_(0)*derxyz_(0, vi)*conv_r_(2, 0, ui) - gridvint_(1)*derxyz_(0, vi)*conv_r_(2, 1, ui) - gridvint_(2)*derxyz_(0, vi)*conv_r_(2, 2, ui)) ;
    estif_(vi*4 + 2, ui*4 + 1) += ttimetauM*(conv_c_(vi)*conv_r_(2, 1, ui) + conv_g_(vi)*conv_r_(2, 1, ui) + velint_(0)*derxyz_(1, vi)*conv_r_(2, 0, ui) + velint_(1)*derxyz_(1, vi)*conv_r_(2, 1, ui) + velint_(2)*derxyz_(1, vi)*conv_r_(2, 2, ui) - gridvint_(0)*derxyz_(1, vi)*conv_r_(2, 0, ui) - gridvint_(1)*derxyz_(1, vi)*conv_r_(2, 1, ui) - gridvint_(2)*derxyz_(1, vi)*conv_r_(2, 2, ui)) ;
    estif_(vi*4 + 2, ui*4 + 2) += ttimetauM*(conv_c_(ui)*conv_c_(vi) + conv_c_(ui)*conv_g_(vi) + conv_c_(vi)*conv_g_(ui) + conv_g_(ui)*conv_g_(vi) + conv_c_(vi)*conv_r_(2, 2, ui) + conv_g_(vi)*conv_r_(2, 2, ui) + velint_(0)*derxyz_(2, vi)*conv_r_(2, 0, ui) + velint_(1)*derxyz_(2, vi)*conv_r_(2, 1, ui) + velint_(2)*derxyz_(2, vi)*conv_r_(2, 2, ui) - gridvint_(0)*derxyz_(2, vi)*conv_r_(2, 0, ui) - gridvint_(1)*derxyz_(2, vi)*conv_r_(2, 1, ui) - gridvint_(2)*derxyz_(2, vi)*conv_r_(2, 2, ui)) ;

    /* Stabilisierung der Konvektion (-L_visc_u) */
    estif_(vi*4, ui*4)         += -2.0*nu_*ttimetauM*(-(conv_c_(vi)*viscs2_(0, 0, ui)) - conv_g_(vi)*viscs2_(0, 0, ui) + funct_(ui)*visc_old_(0)*derxyz_(0, vi)) ;
    estif_(vi*4, ui*4 + 1)     += -2.0*nu_*ttimetauM*(-(conv_c_(vi)*viscs2_(0, 1, ui)) - conv_g_(vi)*viscs2_(0, 1, ui) + funct_(ui)*visc_old_(0)*derxyz_(1, vi)) ;
    estif_(vi*4, ui*4 + 2)     += -2.0*nu_*ttimetauM*(-(conv_c_(vi)*viscs2_(0, 2, ui)) - conv_g_(vi)*viscs2_(0, 2, ui) + funct_(ui)*visc_old_(0)*derxyz_(2, vi)) ;
    estif_(vi*4 + 1, ui*4)     += -2.0*nu_*ttimetauM*(-(conv_c_(vi)*viscs2_(0, 1, ui)) - conv_g_(vi)*viscs2_(0, 1, ui) + funct_(ui)*visc_old_(1)*derxyz_(0, vi)) ;
    estif_(vi*4 + 1, ui*4 + 1) += -2.0*nu_*ttimetauM*(-(conv_c_(vi)*viscs2_(1, 1, ui)) - conv_g_(vi)*viscs2_(1, 1, ui) + funct_(ui)*visc_old_(1)*derxyz_(1, vi)) ;
    estif_(vi*4 + 1, ui*4 + 2) += -2.0*nu_*ttimetauM*(-(conv_c_(vi)*viscs2_(1, 2, ui)) - conv_g_(vi)*viscs2_(1, 2, ui) + funct_(ui)*visc_old_(1)*derxyz_(2, vi)) ;
    estif_(vi*4 + 2, ui*4)     += -2.0*nu_*ttimetauM*(-(conv_c_(vi)*viscs2_(0, 2, ui)) - conv_g_(vi)*viscs2_(0, 2, ui) + funct_(ui)*visc_old_(2)*derxyz_(0, vi)) ;
    estif_(vi*4 + 2, ui*4 + 1) += -2.0*nu_*ttimetauM*(-(conv_c_(vi)*viscs2_(1, 2, ui)) - conv_g_(vi)*viscs2_(1, 2, ui) + funct_(ui)*visc_old_(2)*derxyz_(1, vi)) ;
    estif_(vi*4 + 2, ui*4 + 2) += -2.0*nu_*ttimetauM*(-(conv_c_(vi)*viscs2_(2, 2, ui)) - conv_g_(vi)*viscs2_(2, 2, ui) + funct_(ui)*visc_old_(2)*derxyz_(2, vi)) ;

    /* Stabilisierung der Konvektion ( L_pres_p) */
    estif_(vi*4, ui*4)         += ttimetauM*funct_(ui)*gradp_(0)*derxyz_(0, vi) ;
    estif_(vi*4, ui*4 + 1)     += ttimetauM*funct_(ui)*gradp_(0)*derxyz_(1, vi) ;
    estif_(vi*4, ui*4 + 2)     += ttimetauM*funct_(ui)*gradp_(0)*derxyz_(2, vi) ;
    estif_(vi*4, ui*4 + 3)     += ttimetauM*(conv_c_(vi) + conv_g_(vi))*derxyz_(0, ui) ;
    estif_(vi*4 + 1, ui*4)     += ttimetauM*funct_(ui)*gradp_(1)*derxyz_(0, vi) ;
    estif_(vi*4 + 1, ui*4 + 1) += ttimetauM*funct_(ui)*gradp_(1)*derxyz_(1, vi) ;
    estif_(vi*4 + 1, ui*4 + 2) += ttimetauM*funct_(ui)*gradp_(1)*derxyz_(2, vi) ;
    estif_(vi*4 + 1, ui*4 + 3) += ttimetauM*(conv_c_(vi) + conv_g_(vi))*derxyz_(1, ui) ;
    estif_(vi*4 + 2, ui*4)     += ttimetauM*funct_(ui)*gradp_(2)*derxyz_(0, vi) ;
    estif_(vi*4 + 2, ui*4 + 1) += ttimetauM*funct_(ui)*gradp_(2)*derxyz_(1, vi) ;
    estif_(vi*4 + 2, ui*4 + 2) += ttimetauM*funct_(ui)*gradp_(2)*derxyz_(2, vi) ;
    estif_(vi*4 + 2, ui*4 + 3) += ttimetauM*(conv_c_(vi) + conv_g_(vi))*derxyz_(2, ui) ;

    /* Viskositätsterm */
    estif_(vi*4, ui*4)         += nu_*timefacfac*(2.0*derxyz_(0, ui)*derxyz_(0, vi) + derxyz_(1, ui)*derxyz_(1, vi) + derxyz_(2, ui)*derxyz_(2, vi)) ;
    estif_(vi*4, ui*4 + 1)     += nu_*timefacfac*derxyz_(0, ui)*derxyz_(1, vi) ;
    estif_(vi*4, ui*4 + 2)     += nu_*timefacfac*derxyz_(0, ui)*derxyz_(2, vi) ;
    estif_(vi*4 + 1, ui*4)     += nu_*timefacfac*derxyz_(0, vi)*derxyz_(1, ui) ;
    estif_(vi*4 + 1, ui*4 + 1) += nu_*timefacfac*(derxyz_(0, ui)*derxyz_(0, vi) + 2.0*derxyz_(1, ui)*derxyz_(1, vi) + derxyz_(2, ui)*derxyz_(2, vi)) ;
    estif_(vi*4 + 1, ui*4 + 2) += nu_*timefacfac*derxyz_(1, ui)*derxyz_(2, vi) ;
    estif_(vi*4 + 2, ui*4)     += nu_*timefacfac*derxyz_(0, vi)*derxyz_(2, ui) ;
    estif_(vi*4 + 2, ui*4 + 1) += nu_*timefacfac*derxyz_(1, vi)*derxyz_(2, ui) ;
    estif_(vi*4 + 2, ui*4 + 2) += nu_*timefacfac*(derxyz_(0, ui)*derxyz_(0, vi) + derxyz_(1, ui)*derxyz_(1, vi) + 2.0*derxyz_(2, ui)*derxyz_(2, vi)) ;

    /* Stabilisierung der Viskosität ( L_conv_u) */
    estif_(vi*4, ui*4)         += 2.0*nu_*ttimetauMp*(conv_c_(ui)*viscs2_(0, 0, vi) + conv_g_(ui)*viscs2_(0, 0, vi) + viscs2_(0, 0, vi)*conv_r_(0, 0, ui) + viscs2_(0, 1, vi)*conv_r_(1, 0, ui) + viscs2_(0, 2, vi)*conv_r_(2, 0, ui)) ;
    estif_(vi*4, ui*4 + 1)     += 2.0*nu_*ttimetauMp*(conv_c_(ui)*viscs2_(0, 1, vi) + conv_g_(ui)*viscs2_(0, 1, vi) + viscs2_(0, 0, vi)*conv_r_(0, 1, ui) + viscs2_(0, 1, vi)*conv_r_(1, 1, ui) + viscs2_(0, 2, vi)*conv_r_(2, 1, ui)) ;
    estif_(vi*4, ui*4 + 2)     += 2.0*nu_*ttimetauMp*(conv_c_(ui)*viscs2_(0, 2, vi) + conv_g_(ui)*viscs2_(0, 2, vi) + viscs2_(0, 0, vi)*conv_r_(0, 2, ui) + viscs2_(0, 1, vi)*conv_r_(1, 2, ui) + viscs2_(0, 2, vi)*conv_r_(2, 2, ui)) ;
    estif_(vi*4 + 1, ui*4)     += 2.0*nu_*ttimetauMp*(conv_c_(ui)*viscs2_(0, 1, vi) + conv_g_(ui)*viscs2_(0, 1, vi) + viscs2_(0, 1, vi)*conv_r_(0, 0, ui) + viscs2_(1, 1, vi)*conv_r_(1, 0, ui) + viscs2_(1, 2, vi)*conv_r_(2, 0, ui)) ;
    estif_(vi*4 + 1, ui*4 + 1) += 2.0*nu_*ttimetauMp*(conv_c_(ui)*viscs2_(1, 1, vi) + conv_g_(ui)*viscs2_(1, 1, vi) + viscs2_(0, 1, vi)*conv_r_(0, 1, ui) + viscs2_(1, 1, vi)*conv_r_(1, 1, ui) + viscs2_(1, 2, vi)*conv_r_(2, 1, ui)) ;
    estif_(vi*4 + 1, ui*4 + 2) += 2.0*nu_*ttimetauMp*(conv_c_(ui)*viscs2_(1, 2, vi) + conv_g_(ui)*viscs2_(1, 2, vi) + viscs2_(0, 1, vi)*conv_r_(0, 2, ui) + viscs2_(1, 1, vi)*conv_r_(1, 2, ui) + viscs2_(1, 2, vi)*conv_r_(2, 2, ui)) ;
    estif_(vi*4 + 2, ui*4)     += 2.0*nu_*ttimetauMp*(conv_c_(ui)*viscs2_(0, 2, vi) + conv_g_(ui)*viscs2_(0, 2, vi) + viscs2_(0, 2, vi)*conv_r_(0, 0, ui) + viscs2_(1, 2, vi)*conv_r_(1, 0, ui) + viscs2_(2, 2, vi)*conv_r_(2, 0, ui)) ;
    estif_(vi*4 + 2, ui*4 + 1) += 2.0*nu_*ttimetauMp*(conv_c_(ui)*viscs2_(1, 2, vi) + conv_g_(ui)*viscs2_(1, 2, vi) + viscs2_(0, 2, vi)*conv_r_(0, 1, ui) + viscs2_(1, 2, vi)*conv_r_(1, 1, ui) + viscs2_(2, 2, vi)*conv_r_(2, 1, ui)) ;
    estif_(vi*4 + 2, ui*4 + 2) += 2.0*nu_*ttimetauMp*(conv_c_(ui)*viscs2_(2, 2, vi) + conv_g_(ui)*viscs2_(2, 2, vi) + viscs2_(0, 2, vi)*conv_r_(0, 2, ui) + viscs2_(1, 2, vi)*conv_r_(1, 2, ui) + viscs2_(2, 2, vi)*conv_r_(2, 2, ui)) ;

    /* Stabilisierung der Viskosität (-L_visc_u) */
    estif_(vi*4, ui*4)         += 4.0*(nu_*nu_)*ttimetauMp*(viscs2_(0, 0, ui)*viscs2_(0, 0, vi) + viscs2_(0, 1, ui)*viscs2_(0, 1, vi) + viscs2_(0, 2, ui)*viscs2_(0, 2, vi)) ;
    estif_(vi*4, ui*4 + 1)     += 4.0*(nu_*nu_)*ttimetauMp*(viscs2_(0, 0, vi)*viscs2_(0, 1, ui) + viscs2_(0, 1, vi)*viscs2_(1, 1, ui) + viscs2_(0, 2, vi)*viscs2_(1, 2, ui)) ;
    estif_(vi*4, ui*4 + 2)     += 4.0*(nu_*nu_)*ttimetauMp*(viscs2_(0, 0, vi)*viscs2_(0, 2, ui) + viscs2_(0, 1, vi)*viscs2_(1, 2, ui) + viscs2_(0, 2, vi)*viscs2_(2, 2, ui)) ;
    estif_(vi*4 + 1, ui*4)     += 4.0*(nu_*nu_)*ttimetauMp*(viscs2_(0, 0, ui)*viscs2_(0, 1, vi) + viscs2_(0, 1, ui)*viscs2_(1, 1, vi) + viscs2_(0, 2, ui)*viscs2_(1, 2, vi)) ;
    estif_(vi*4 + 1, ui*4 + 1) += 4.0*(nu_*nu_)*ttimetauMp*(viscs2_(0, 1, ui)*viscs2_(0, 1, vi) + viscs2_(1, 1, ui)*viscs2_(1, 1, vi) + viscs2_(1, 2, ui)*viscs2_(1, 2, vi)) ;
    estif_(vi*4 + 1, ui*4 + 2) += 4.0*(nu_*nu_)*ttimetauMp*(viscs2_(0, 1, vi)*viscs2_(0, 2, ui) + viscs2_(1, 1, vi)*viscs2_(1, 2, ui) + viscs2_(1, 2, vi)*viscs2_(2, 2, ui)) ;
    estif_(vi*4 + 2, ui*4)     += 4.0*(nu_*nu_)*ttimetauMp*(viscs2_(0, 0, ui)*viscs2_(0, 2, vi) + viscs2_(0, 1, ui)*viscs2_(1, 2, vi) + viscs2_(0, 2, ui)*viscs2_(2, 2, vi)) ;
    estif_(vi*4 + 2, ui*4 + 1) += 4.0*(nu_*nu_)*ttimetauMp*(viscs2_(0, 1, ui)*viscs2_(0, 2, vi) + viscs2_(1, 1, ui)*viscs2_(1, 2, vi) + viscs2_(1, 2, ui)*viscs2_(2, 2, vi)) ;
    estif_(vi*4 + 2, ui*4 + 2) += 4.0*(nu_*nu_)*ttimetauMp*(viscs2_(0, 2, ui)*viscs2_(0, 2, vi) + viscs2_(1, 2, ui)*viscs2_(1, 2, vi) + viscs2_(2, 2, ui)*viscs2_(2, 2, vi)) ;

    /* Stabilisierung der Viskosität ( L_pres_p) */
    estif_(vi*4, ui*4 + 3)     += 2.0*nu_*ttimetauMp*(derxyz_(0, ui)*viscs2_(0, 0, vi) + derxyz_(1, ui)*viscs2_(0, 1, vi) + derxyz_(2, ui)*viscs2_(0, 2, vi)) ;
    estif_(vi*4 + 1, ui*4 + 3) += 2.0*nu_*ttimetauMp*(derxyz_(0, ui)*viscs2_(0, 1, vi) + derxyz_(1, ui)*viscs2_(1, 1, vi) + derxyz_(2, ui)*viscs2_(1, 2, vi)) ;
    estif_(vi*4 + 2, ui*4 + 3) += 2.0*nu_*ttimetauMp*(derxyz_(0, ui)*viscs2_(0, 2, vi) + derxyz_(1, ui)*viscs2_(1, 2, vi) + derxyz_(2, ui)*viscs2_(2, 2, vi)) ;

    /* Druckterm */
    estif_(vi*4, ui*4 + 3)     += -(timefacfac*funct_(ui)*derxyz_(0, vi)) ;
    estif_(vi*4 + 1, ui*4 + 3) += -(timefacfac*funct_(ui)*derxyz_(1, vi)) ;
    estif_(vi*4 + 2, ui*4 + 3) += -(timefacfac*funct_(ui)*derxyz_(2, vi)) ;

    /* Stabilisierung des Drucks ( L_conv_u) */
    estif_(vi*4 + 3, ui*4)     += ttimetauMp*(conv_c_(ui)*derxyz_(0, vi) + conv_g_(ui)*derxyz_(0, vi) + derxyz_(0, vi)*conv_r_(0, 0, ui) + derxyz_(1, vi)*conv_r_(1, 0, ui) + derxyz_(2, vi)*conv_r_(2, 0, ui)) ;
    estif_(vi*4 + 3, ui*4 + 1) += ttimetauMp*(conv_c_(ui)*derxyz_(1, vi) + conv_g_(ui)*derxyz_(1, vi) + derxyz_(0, vi)*conv_r_(0, 1, ui) + derxyz_(1, vi)*conv_r_(1, 1, ui) + derxyz_(2, vi)*conv_r_(2, 1, ui)) ;
    estif_(vi*4 + 3, ui*4 + 2) += ttimetauMp*(conv_c_(ui)*derxyz_(2, vi) + conv_g_(ui)*derxyz_(2, vi) + derxyz_(0, vi)*conv_r_(0, 2, ui) + derxyz_(1, vi)*conv_r_(1, 2, ui) + derxyz_(2, vi)*conv_r_(2, 2, ui)) ;

    /* Stabilisierung des Drucks (-L_visc_u) */
    estif_(vi*4 + 3, ui*4)     += 2.0*nu_*ttimetauMp*(derxyz_(0, vi)*viscs2_(0, 0, ui) + derxyz_(1, vi)*viscs2_(0, 1, ui) + derxyz_(2, vi)*viscs2_(0, 2, ui)) ;
    estif_(vi*4 + 3, ui*4 + 1) += 2.0*nu_*ttimetauMp*(derxyz_(0, vi)*viscs2_(0, 1, ui) + derxyz_(1, vi)*viscs2_(1, 1, ui) + derxyz_(2, vi)*viscs2_(1, 2, ui)) ;
    estif_(vi*4 + 3, ui*4 + 2) += 2.0*nu_*ttimetauMp*(derxyz_(0, vi)*viscs2_(0, 2, ui) + derxyz_(1, vi)*viscs2_(1, 2, ui) + derxyz_(2, vi)*viscs2_(2, 2, ui)) ;

    /* Stabilisierung des Drucks ( L_pres_p) */
    estif_(vi*4 + 3, ui*4 + 3) += ttimetauMp*(derxyz_(0, ui)*derxyz_(0, vi) + derxyz_(1, ui)*derxyz_(1, vi) + derxyz_(2, ui)*derxyz_(2, vi)) ;

    /* Divergenzfreiheit */
    estif_(vi*4 + 3, ui*4)     += timefacfac*funct_(vi)*derxyz_(0, ui) ;
    estif_(vi*4 + 3, ui*4 + 1) += timefacfac*funct_(vi)*derxyz_(1, ui) ;
    estif_(vi*4 + 3, ui*4 + 2) += timefacfac*funct_(vi)*derxyz_(2, ui) ;

    /* Kontinuitätsstabilisierung */
    estif_(vi*4, ui*4)         += (thsl*thsl)*tau_C*derxyz_(0, ui)*derxyz_(0, vi) ;
    estif_(vi*4, ui*4 + 1)     += (thsl*thsl)*tau_C*derxyz_(0, vi)*derxyz_(1, ui) ;
    estif_(vi*4, ui*4 + 2)     += (thsl*thsl)*tau_C*derxyz_(0, vi)*derxyz_(2, ui) ;
    estif_(vi*4 + 1, ui*4)     += (thsl*thsl)*tau_C*derxyz_(0, ui)*derxyz_(1, vi) ;
    estif_(vi*4 + 1, ui*4 + 1) += (thsl*thsl)*tau_C*derxyz_(1, ui)*derxyz_(1, vi) ;
    estif_(vi*4 + 1, ui*4 + 2) += (thsl*thsl)*tau_C*derxyz_(1, vi)*derxyz_(2, ui) ;
    estif_(vi*4 + 2, ui*4)     += (thsl*thsl)*tau_C*derxyz_(0, ui)*derxyz_(2, vi) ;
    estif_(vi*4 + 2, ui*4 + 1) += (thsl*thsl)*tau_C*derxyz_(1, ui)*derxyz_(2, vi) ;
    estif_(vi*4 + 2, ui*4 + 2) += (thsl*thsl)*tau_C*derxyz_(2, ui)*derxyz_(2, vi) ;

    /* Massenterm */
    estif_(vi*4, ui*4)         += fac*funct_(ui)*funct_(vi) ;
    estif_(vi*4 + 1, ui*4 + 1) += fac*funct_(ui)*funct_(vi) ;
    estif_(vi*4 + 2, ui*4 + 2) += fac*funct_(ui)*funct_(vi) ;

    /* Konvektionsstabilisierung */
    estif_(vi*4, ui*4)         += timetauM*funct_(ui)*(conv_g_(vi) + 2.0*velint_(0)*derxyz_(0, vi) + velint_(1)*derxyz_(1, vi) + velint_(2)*derxyz_(2, vi)) ;
    estif_(vi*4, ui*4 + 1)     += timetauM*funct_(ui)*velint_(0)*derxyz_(1, vi) ;
    estif_(vi*4, ui*4 + 2)     += timetauM*funct_(ui)*velint_(0)*derxyz_(2, vi) ;
    estif_(vi*4 + 1, ui*4)     += timetauM*funct_(ui)*velint_(1)*derxyz_(0, vi) ;
    estif_(vi*4 + 1, ui*4 + 1) += timetauM*funct_(ui)*(conv_g_(vi) + velint_(0)*derxyz_(0, vi) + 2.0*velint_(1)*derxyz_(1, vi) + velint_(2)*derxyz_(2, vi)) ;
    estif_(vi*4 + 1, ui*4 + 2) += timetauM*funct_(ui)*velint_(1)*derxyz_(2, vi) ;
    estif_(vi*4 + 2, ui*4)     += timetauM*funct_(ui)*velint_(2)*derxyz_(0, vi) ;
    estif_(vi*4 + 2, ui*4 + 1) += timetauM*funct_(ui)*velint_(2)*derxyz_(1, vi) ;
    estif_(vi*4 + 2, ui*4 + 2) += timetauM*funct_(ui)*(conv_g_(vi) + velint_(0)*derxyz_(0, vi) + velint_(1)*derxyz_(1, vi) + 2.0*velint_(2)*derxyz_(2, vi)) ;

    /* Viskositätsstabilisierung */
    estif_(vi*4, ui*4)         += 2.0*nu_*timetauMp*funct_(ui)*viscs2_(0, 0, vi) ;
    estif_(vi*4, ui*4 + 1)     += 2.0*nu_*timetauMp*funct_(ui)*viscs2_(0, 1, vi) ;
    estif_(vi*4, ui*4 + 2)     += 2.0*nu_*timetauMp*funct_(ui)*viscs2_(0, 2, vi) ;
    estif_(vi*4 + 1, ui*4)     += 2.0*nu_*timetauMp*funct_(ui)*viscs2_(0, 1, vi) ;
    estif_(vi*4 + 1, ui*4 + 1) += 2.0*nu_*timetauMp*funct_(ui)*viscs2_(1, 1, vi) ;
    estif_(vi*4 + 1, ui*4 + 2) += 2.0*nu_*timetauMp*funct_(ui)*viscs2_(1, 2, vi) ;
    estif_(vi*4 + 2, ui*4)     += 2.0*nu_*timetauMp*funct_(ui)*viscs2_(0, 2, vi) ;
    estif_(vi*4 + 2, ui*4 + 1) += 2.0*nu_*timetauMp*funct_(ui)*viscs2_(1, 2, vi) ;
    estif_(vi*4 + 2, ui*4 + 2) += 2.0*nu_*timetauMp*funct_(ui)*viscs2_(2, 2, vi) ;

    /* Stabilisierung der Druckgleichung */
    estif_(vi*4 + 3, ui*4)     += timetauMp*funct_(ui)*derxyz_(0, vi) ;
    estif_(vi*4 + 3, ui*4 + 1) += timetauMp*funct_(ui)*derxyz_(1, vi) ;
    estif_(vi*4 + 3, ui*4 + 2) += timetauMp*funct_(ui)*derxyz_(2, vi) ;

    /* Quellterm der rechten Seite */

    /* Konvektionsstabilisierung */
    estif_(vi*4, ui*4)         += -(timetauM*funct_(ui)*rhsint_(0)*derxyz_(0, vi)) ;
    estif_(vi*4, ui*4 + 1)     += -(timetauM*funct_(ui)*rhsint_(0)*derxyz_(1, vi)) ;
    estif_(vi*4, ui*4 + 2)     += -(timetauM*funct_(ui)*rhsint_(0)*derxyz_(2, vi)) ;
    estif_(vi*4 + 1, ui*4)     += -(timetauM*funct_(ui)*rhsint_(1)*derxyz_(0, vi)) ;
    estif_(vi*4 + 1, ui*4 + 1) += -(timetauM*funct_(ui)*rhsint_(1)*derxyz_(1, vi)) ;
    estif_(vi*4 + 1, ui*4 + 2) += -(timetauM*funct_(ui)*rhsint_(1)*derxyz_(2, vi)) ;
    estif_(vi*4 + 2, ui*4)     += -(timetauM*funct_(ui)*rhsint_(2)*derxyz_(0, vi)) ;
    estif_(vi*4 + 2, ui*4 + 1) += -(timetauM*funct_(ui)*rhsint_(2)*derxyz_(1, vi)) ;
    estif_(vi*4 + 2, ui*4 + 2) += -(timetauM*funct_(ui)*rhsint_(2)*derxyz_(2, vi)) ;

    /* Viskositätsstabilisierung */

    /* Stabilisierung der Druckgleichung */

  }
}
