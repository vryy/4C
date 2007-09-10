#ifndef CCADISCRET
for (vi=0; vi<iel; ++vi)
{
  /* Konvektionsterm */
  eforce_(vi*4)     += timefacfac*(-(velint_(0)*conv_r_(0, 0, vi)) - velint_(1)*conv_r_(0, 1, vi) - velint_(2)*conv_r_(0, 2, vi) + gridvint_(0)*conv_r_(0, 0, vi) + gridvint_(1)*conv_r_(0, 1, vi) + gridvint_(2)*conv_r_(0, 2, vi)) ;
  eforce_(vi*4 + 1) += timefacfac*(-(velint_(0)*conv_r_(1, 0, vi)) - velint_(1)*conv_r_(1, 1, vi) - velint_(2)*conv_r_(1, 2, vi) + gridvint_(0)*conv_r_(1, 0, vi) + gridvint_(1)*conv_r_(1, 1, vi) + gridvint_(2)*conv_r_(1, 2, vi)) ;
  eforce_(vi*4 + 2) += -(timefacfac*(velint_(0)*conv_r_(2, 0, vi) + velint_(1)*conv_r_(2, 1, vi) + velint_(2)*conv_r_(2, 2, vi) - gridvint_(0)*conv_r_(2, 0, vi) - gridvint_(1)*conv_r_(2, 1, vi) - gridvint_(2)*conv_r_(2, 2, vi))) ;

  /* Stabilisierung der Konvektion ( L_conv_u) */
  eforce_(vi*4)     += ttimetauM*(-(conv_c_(vi)*conv_old_(0)) + conv_g_(vi)*gridvint_(0)*vderxyz_(0, 0) - (gridvint_(1)*gridvint_(1))*derxyz_(1, vi)*vderxyz_(0, 1) - (gridvint_(2)*gridvint_(2))*derxyz_(2, vi)*vderxyz_(0, 2) + 2.0*velint_(0)*gridvint_(0)*derxyz_(0, vi)*vderxyz_(0, 0) + velint_(0)*gridvint_(1)*derxyz_(0, vi)*vderxyz_(0, 1) + velint_(0)*gridvint_(1)*derxyz_(1, vi)*vderxyz_(0, 0) + velint_(1)*gridvint_(0)*derxyz_(0, vi)*vderxyz_(0, 1) + velint_(1)*gridvint_(0)*derxyz_(1, vi)*vderxyz_(0, 0) + velint_(0)*gridvint_(2)*derxyz_(0, vi)*vderxyz_(0, 2) + velint_(0)*gridvint_(2)*derxyz_(2, vi)*vderxyz_(0, 0) + velint_(2)*gridvint_(0)*derxyz_(0, vi)*vderxyz_(0, 2) + velint_(2)*gridvint_(0)*derxyz_(2, vi)*vderxyz_(0, 0) + 2.0*velint_(1)*gridvint_(1)*derxyz_(1, vi)*vderxyz_(0, 1) + velint_(1)*gridvint_(2)*derxyz_(1, vi)*vderxyz_(0, 2) + velint_(1)*gridvint_(2)*derxyz_(2, vi)*vderxyz_(0, 1) + velint_(2)*gridvint_(1)*derxyz_(1, vi)*vderxyz_(0, 2) + velint_(2)*gridvint_(1)*derxyz_(2, vi)*vderxyz_(0, 1) + 2.0*velint_(2)*gridvint_(2)*derxyz_(2, vi)*vderxyz_(0, 2) - gridvint_(0)*gridvint_(1)*derxyz_(0, vi)*vderxyz_(0, 1) - gridvint_(0)*gridvint_(2)*derxyz_(0, vi)*vderxyz_(0, 2) - gridvint_(1)*gridvint_(2)*derxyz_(1, vi)*vderxyz_(0, 2) - gridvint_(1)*gridvint_(2)*derxyz_(2, vi)*vderxyz_(0, 1)) ;
  eforce_(vi*4 + 1) += ttimetauM*(-(conv_c_(vi)*conv_old_(1)) + conv_g_(vi)*gridvint_(0)*vderxyz_(1, 0) - (gridvint_(1)*gridvint_(1))*derxyz_(1, vi)*vderxyz_(1, 1) - (gridvint_(2)*gridvint_(2))*derxyz_(2, vi)*vderxyz_(1, 2) + 2.0*velint_(0)*gridvint_(0)*derxyz_(0, vi)*vderxyz_(1, 0) + velint_(0)*gridvint_(1)*derxyz_(0, vi)*vderxyz_(1, 1) + velint_(0)*gridvint_(1)*derxyz_(1, vi)*vderxyz_(1, 0) + velint_(1)*gridvint_(0)*derxyz_(0, vi)*vderxyz_(1, 1) + velint_(1)*gridvint_(0)*derxyz_(1, vi)*vderxyz_(1, 0) + velint_(0)*gridvint_(2)*derxyz_(0, vi)*vderxyz_(1, 2) + velint_(0)*gridvint_(2)*derxyz_(2, vi)*vderxyz_(1, 0) + velint_(2)*gridvint_(0)*derxyz_(0, vi)*vderxyz_(1, 2) + velint_(2)*gridvint_(0)*derxyz_(2, vi)*vderxyz_(1, 0) + 2.0*velint_(1)*gridvint_(1)*derxyz_(1, vi)*vderxyz_(1, 1) + velint_(1)*gridvint_(2)*derxyz_(1, vi)*vderxyz_(1, 2) + velint_(1)*gridvint_(2)*derxyz_(2, vi)*vderxyz_(1, 1) + velint_(2)*gridvint_(1)*derxyz_(1, vi)*vderxyz_(1, 2) + velint_(2)*gridvint_(1)*derxyz_(2, vi)*vderxyz_(1, 1) + 2.0*velint_(2)*gridvint_(2)*derxyz_(2, vi)*vderxyz_(1, 2) - gridvint_(0)*gridvint_(1)*derxyz_(0, vi)*vderxyz_(1, 1) - gridvint_(0)*gridvint_(2)*derxyz_(0, vi)*vderxyz_(1, 2) - gridvint_(1)*gridvint_(2)*derxyz_(1, vi)*vderxyz_(1, 2) - gridvint_(1)*gridvint_(2)*derxyz_(2, vi)*vderxyz_(1, 1)) ;
  eforce_(vi*4 + 2) += ttimetauM*(-(conv_c_(vi)*conv_old_(2)) + conv_g_(vi)*gridvint_(0)*vderxyz_(2, 0) - (gridvint_(1)*gridvint_(1))*derxyz_(1, vi)*vderxyz_(2, 1) - (gridvint_(2)*gridvint_(2))*derxyz_(2, vi)*vderxyz_(2, 2) + 2.0*velint_(0)*gridvint_(0)*derxyz_(0, vi)*vderxyz_(2, 0) + velint_(0)*gridvint_(1)*derxyz_(0, vi)*vderxyz_(2, 1) + velint_(0)*gridvint_(1)*derxyz_(1, vi)*vderxyz_(2, 0) + velint_(1)*gridvint_(0)*derxyz_(0, vi)*vderxyz_(2, 1) + velint_(1)*gridvint_(0)*derxyz_(1, vi)*vderxyz_(2, 0) + velint_(0)*gridvint_(2)*derxyz_(0, vi)*vderxyz_(2, 2) + velint_(0)*gridvint_(2)*derxyz_(2, vi)*vderxyz_(2, 0) + velint_(2)*gridvint_(0)*derxyz_(0, vi)*vderxyz_(2, 2) + velint_(2)*gridvint_(0)*derxyz_(2, vi)*vderxyz_(2, 0) + 2.0*velint_(1)*gridvint_(1)*derxyz_(1, vi)*vderxyz_(2, 1) + velint_(1)*gridvint_(2)*derxyz_(1, vi)*vderxyz_(2, 2) + velint_(1)*gridvint_(2)*derxyz_(2, vi)*vderxyz_(2, 1) + velint_(2)*gridvint_(1)*derxyz_(1, vi)*vderxyz_(2, 2) + velint_(2)*gridvint_(1)*derxyz_(2, vi)*vderxyz_(2, 1) + 2.0*velint_(2)*gridvint_(2)*derxyz_(2, vi)*vderxyz_(2, 2) - gridvint_(0)*gridvint_(1)*derxyz_(0, vi)*vderxyz_(2, 1) - gridvint_(0)*gridvint_(2)*derxyz_(0, vi)*vderxyz_(2, 2) - gridvint_(1)*gridvint_(2)*derxyz_(1, vi)*vderxyz_(2, 2) - gridvint_(1)*gridvint_(2)*derxyz_(2, vi)*vderxyz_(2, 1)) ;

  /* Stabilisierung der Konvektion (-L_visc_u) */
  eforce_(vi*4)     += 2.0*nu_*ttimetauM*(conv_c_(vi) + conv_g_(vi))*visc_old_(0) ;
  eforce_(vi*4 + 1) += 2.0*nu_*ttimetauM*(conv_c_(vi) + conv_g_(vi))*visc_old_(1) ;
  eforce_(vi*4 + 2) += 2.0*nu_*ttimetauM*(conv_c_(vi) + conv_g_(vi))*visc_old_(2) ;

  /* Stabilisierung der Konvektion ( L_pres_p) */
  eforce_(vi*4)     += -(ttimetauM*(conv_c_(vi) + conv_g_(vi))*gradp_(0)) ;
  eforce_(vi*4 + 1) += -(ttimetauM*(conv_c_(vi) + conv_g_(vi))*gradp_(1)) ;
  eforce_(vi*4 + 2) += -(ttimetauM*(conv_c_(vi) + conv_g_(vi))*gradp_(2)) ;

  /* Viskositätsterm */
  eforce_(vi*4)     += -(nu_*timefacfac*(2.0*derxyz_(0, vi)*vderxyz_(0, 0) + derxyz_(1, vi)*vderxyz_(0, 1) + derxyz_(1, vi)*vderxyz_(1, 0) + derxyz_(2, vi)*vderxyz_(0, 2) + derxyz_(2, vi)*vderxyz_(2, 0))) ;
  eforce_(vi*4 + 1) += -(nu_*timefacfac*(derxyz_(0, vi)*vderxyz_(0, 1) + derxyz_(0, vi)*vderxyz_(1, 0) + 2.0*derxyz_(1, vi)*vderxyz_(1, 1) + derxyz_(2, vi)*vderxyz_(1, 2) + derxyz_(2, vi)*vderxyz_(2, 1))) ;
  eforce_(vi*4 + 2) += -(nu_*timefacfac*(derxyz_(0, vi)*vderxyz_(0, 2) + derxyz_(0, vi)*vderxyz_(2, 0) + derxyz_(1, vi)*vderxyz_(1, 2) + derxyz_(1, vi)*vderxyz_(2, 1) + 2.0*derxyz_(2, vi)*vderxyz_(2, 2))) ;

  /* Stabilisierung der Viskosität ( L_conv_u) */
  eforce_(vi*4)     += -2.0*nu_*ttimetauMp*(conv_old_(0)*viscs2_(0, 0, vi) + conv_old_(1)*viscs2_(0, 1, vi) + conv_old_(2)*viscs2_(0, 2, vi) - gridvint_(0)*vderxyz_(0, 0)*viscs2_(0, 0, vi) - gridvint_(0)*vderxyz_(1, 0)*viscs2_(0, 1, vi) - gridvint_(1)*vderxyz_(0, 1)*viscs2_(0, 0, vi) - gridvint_(0)*vderxyz_(2, 0)*viscs2_(0, 2, vi) - gridvint_(2)*vderxyz_(0, 2)*viscs2_(0, 0, vi) - gridvint_(1)*vderxyz_(1, 1)*viscs2_(0, 1, vi) - gridvint_(1)*vderxyz_(2, 1)*viscs2_(0, 2, vi) - gridvint_(2)*vderxyz_(1, 2)*viscs2_(0, 1, vi) - gridvint_(2)*vderxyz_(2, 2)*viscs2_(0, 2, vi)) ;
  eforce_(vi*4 + 1) += 2.0*nu_*ttimetauMp*(-(conv_old_(0)*viscs2_(0, 1, vi)) - conv_old_(1)*viscs2_(1, 1, vi) - conv_old_(2)*viscs2_(1, 2, vi) + gridvint_(0)*vderxyz_(0, 0)*viscs2_(0, 1, vi) + gridvint_(0)*vderxyz_(1, 0)*viscs2_(1, 1, vi) + gridvint_(1)*vderxyz_(0, 1)*viscs2_(0, 1, vi) + gridvint_(0)*vderxyz_(2, 0)*viscs2_(1, 2, vi) + gridvint_(2)*vderxyz_(0, 2)*viscs2_(0, 1, vi) + gridvint_(1)*vderxyz_(1, 1)*viscs2_(1, 1, vi) + gridvint_(1)*vderxyz_(2, 1)*viscs2_(1, 2, vi) + gridvint_(2)*vderxyz_(1, 2)*viscs2_(1, 1, vi) + gridvint_(2)*vderxyz_(2, 2)*viscs2_(1, 2, vi)) ;
  eforce_(vi*4 + 2) += 2.0*nu_*ttimetauMp*(-(conv_old_(0)*viscs2_(0, 2, vi)) - conv_old_(1)*viscs2_(1, 2, vi) - conv_old_(2)*viscs2_(2, 2, vi) + gridvint_(0)*vderxyz_(0, 0)*viscs2_(0, 2, vi) + gridvint_(0)*vderxyz_(1, 0)*viscs2_(1, 2, vi) + gridvint_(1)*vderxyz_(0, 1)*viscs2_(0, 2, vi) + gridvint_(0)*vderxyz_(2, 0)*viscs2_(2, 2, vi) + gridvint_(2)*vderxyz_(0, 2)*viscs2_(0, 2, vi) + gridvint_(1)*vderxyz_(1, 1)*viscs2_(1, 2, vi) + gridvint_(1)*vderxyz_(2, 1)*viscs2_(2, 2, vi) + gridvint_(2)*vderxyz_(1, 2)*viscs2_(1, 2, vi) + gridvint_(2)*vderxyz_(2, 2)*viscs2_(2, 2, vi)) ;

  /* Stabilisierung der Viskosität (-L_visc_u) */
  eforce_(vi*4)     += 4.0*(nu_*nu_)*ttimetauMp*(visc_old_(0)*viscs2_(0, 0, vi) + visc_old_(1)*viscs2_(0, 1, vi) + visc_old_(2)*viscs2_(0, 2, vi)) ;
  eforce_(vi*4 + 1) += 4.0*(nu_*nu_)*ttimetauMp*(visc_old_(0)*viscs2_(0, 1, vi) + visc_old_(1)*viscs2_(1, 1, vi) + visc_old_(2)*viscs2_(1, 2, vi)) ;
  eforce_(vi*4 + 2) += 4.0*(nu_*nu_)*ttimetauMp*(visc_old_(0)*viscs2_(0, 2, vi) + visc_old_(1)*viscs2_(1, 2, vi) + visc_old_(2)*viscs2_(2, 2, vi)) ;

  /* Stabilisierung der Viskosität ( L_pres_p) */
  eforce_(vi*4)     += -2.0*nu_*ttimetauMp*(gradp_(0)*viscs2_(0, 0, vi) + gradp_(1)*viscs2_(0, 1, vi) + gradp_(2)*viscs2_(0, 2, vi)) ;
  eforce_(vi*4 + 1) += -2.0*nu_*ttimetauMp*(gradp_(0)*viscs2_(0, 1, vi) + gradp_(1)*viscs2_(1, 1, vi) + gradp_(2)*viscs2_(1, 2, vi)) ;
  eforce_(vi*4 + 2) += -2.0*nu_*ttimetauMp*(gradp_(0)*viscs2_(0, 2, vi) + gradp_(1)*viscs2_(1, 2, vi) + gradp_(2)*viscs2_(2, 2, vi)) ;

  /* Druckterm */
  eforce_(vi*4)     += press*timefacfac*derxyz_(0, vi) ;
  eforce_(vi*4 + 1) += press*timefacfac*derxyz_(1, vi) ;
  eforce_(vi*4 + 2) += press*timefacfac*derxyz_(2, vi) ;

  /* Stabilisierung des Drucks ( L_conv_u) */
  eforce_(vi*4 + 3) += ttimetauMp*(-(conv_old_(0)*derxyz_(0, vi)) - conv_old_(1)*derxyz_(1, vi) - conv_old_(2)*derxyz_(2, vi) + gridvint_(0)*derxyz_(0, vi)*vderxyz_(0, 0) + gridvint_(0)*derxyz_(1, vi)*vderxyz_(1, 0) + gridvint_(1)*derxyz_(0, vi)*vderxyz_(0, 1) + gridvint_(0)*derxyz_(2, vi)*vderxyz_(2, 0) + gridvint_(2)*derxyz_(0, vi)*vderxyz_(0, 2) + gridvint_(1)*derxyz_(1, vi)*vderxyz_(1, 1) + gridvint_(1)*derxyz_(2, vi)*vderxyz_(2, 1) + gridvint_(2)*derxyz_(1, vi)*vderxyz_(1, 2) + gridvint_(2)*derxyz_(2, vi)*vderxyz_(2, 2)) ;

  /* Stabilisierung des Drucks (-L_visc_u) */
  eforce_(vi*4 + 3) += 2.0*nu_*ttimetauMp*(visc_old_(0)*derxyz_(0, vi) + visc_old_(1)*derxyz_(1, vi) + visc_old_(2)*derxyz_(2, vi)) ;

  /* Stabilisierung des Drucks ( L_pres_p) */
  eforce_(vi*4 + 3) += -(ttimetauMp*(gradp_(0)*derxyz_(0, vi) + gradp_(1)*derxyz_(1, vi) + gradp_(2)*derxyz_(2, vi))) ;

  /* Divergenzfreiheit */
  eforce_(vi*4 + 3) += -(timefacfac*(conv_r_(0, 0, vi) + conv_r_(1, 1, vi) + conv_r_(2, 2, vi))) ;

  /* Kontinuitätsstabilisierung */
  eforce_(vi*4)     += -((thsl*thsl)*tau_C*derxyz_(0, vi)*(vderxyz_(0, 0) + vderxyz_(1, 1) + vderxyz_(2, 2))) ;
  eforce_(vi*4 + 1) += -((thsl*thsl)*tau_C*derxyz_(1, vi)*(vderxyz_(0, 0) + vderxyz_(1, 1) + vderxyz_(2, 2))) ;
  eforce_(vi*4 + 2) += -((thsl*thsl)*tau_C*derxyz_(2, vi)*(vderxyz_(0, 0) + vderxyz_(1, 1) + vderxyz_(2, 2))) ;

  /* Massenterm */
  eforce_(vi*4)     += -(fac*funct_(vi)*velint_(0)) ;
  eforce_(vi*4 + 1) += -(fac*funct_(vi)*velint_(1)) ;
  eforce_(vi*4 + 2) += -(fac*funct_(vi)*velint_(2)) ;

  /* Konvektionsstabilisierung */
  eforce_(vi*4)     += -(timetauM*(conv_c_(vi) + conv_g_(vi))*velint_(0)) ;
  eforce_(vi*4 + 1) += -(timetauM*(conv_c_(vi) + conv_g_(vi))*velint_(1)) ;
  eforce_(vi*4 + 2) += -(timetauM*(conv_c_(vi) + conv_g_(vi))*velint_(2)) ;

  /* Viskositätsstabilisierung */
  eforce_(vi*4)     += -2.0*nu_*timetauMp*(velint_(0)*viscs2_(0, 0, vi) + velint_(1)*viscs2_(0, 1, vi) + velint_(2)*viscs2_(0, 2, vi)) ;
  eforce_(vi*4 + 1) += -2.0*nu_*timetauMp*(velint_(0)*viscs2_(0, 1, vi) + velint_(1)*viscs2_(1, 1, vi) + velint_(2)*viscs2_(1, 2, vi)) ;
  eforce_(vi*4 + 2) += -2.0*nu_*timetauMp*(velint_(0)*viscs2_(0, 2, vi) + velint_(1)*viscs2_(1, 2, vi) + velint_(2)*viscs2_(2, 2, vi)) ;

  /* Stabilisierung der Druckgleichung */
  eforce_(vi*4 + 3) += -(timetauMp*conv_c_(vi)) ;

  /* Quellterm der rechten Seite */
  eforce_(vi*4)     += fac*funct_(vi)*rhsint_(0) ;
  eforce_(vi*4 + 1) += fac*funct_(vi)*rhsint_(1) ;
  eforce_(vi*4 + 2) += fac*funct_(vi)*rhsint_(2) ;

  /* Konvektionsstabilisierung */
  eforce_(vi*4)     += timetauM*(conv_c_(vi) + conv_g_(vi))*rhsint_(0) ;
  eforce_(vi*4 + 1) += timetauM*(conv_c_(vi) + conv_g_(vi))*rhsint_(1) ;
  eforce_(vi*4 + 2) += timetauM*(conv_c_(vi) + conv_g_(vi))*rhsint_(2) ;

  /* Viskositätsstabilisierung */
  eforce_(vi*4)     += 2.0*nu_*timetauMp*(rhsint_(0)*viscs2_(0, 0, vi) + rhsint_(1)*viscs2_(0, 1, vi) + rhsint_(2)*viscs2_(0, 2, vi)) ;
  eforce_(vi*4 + 1) += 2.0*nu_*timetauMp*(rhsint_(0)*viscs2_(0, 1, vi) + rhsint_(1)*viscs2_(1, 1, vi) + rhsint_(2)*viscs2_(1, 2, vi)) ;
  eforce_(vi*4 + 2) += 2.0*nu_*timetauMp*(rhsint_(0)*viscs2_(0, 2, vi) + rhsint_(1)*viscs2_(1, 2, vi) + rhsint_(2)*viscs2_(2, 2, vi)) ;

  /* Stabilisierung der Druckgleichung */
  eforce_(vi*4 + 3) += timetauMp*(rhsint_(0)*derxyz_(0, vi) + rhsint_(1)*derxyz_(1, vi) + rhsint_(2)*derxyz_(2, vi)) ;

}
#endif
