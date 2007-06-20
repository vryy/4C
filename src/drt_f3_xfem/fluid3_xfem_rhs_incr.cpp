for (int vi=0; vi<iel; ++vi)
{
  const int dUDF = vi*4;
  const int dVDF = dUDF + 1;
  const int dWDF = dUDF + 2;
  const int dPDF = dUDF + 3;
  
  /* Konvektionsterm */
  eforce_(dUDF) += -(timefacfac*(velint_(0)*conv_r_(0, 0, vi) + velint_(1)*conv_r_(0, 1, vi) + velint_(2)*conv_r_(0, 2, vi))) ;
  eforce_(dVDF) += -(timefacfac*(velint_(0)*conv_r_(1, 0, vi) + velint_(1)*conv_r_(1, 1, vi) + velint_(2)*conv_r_(1, 2, vi))) ;
  eforce_(dWDF) += -(timefacfac*(velint_(0)*conv_r_(2, 0, vi) + velint_(1)*conv_r_(2, 1, vi) + velint_(2)*conv_r_(2, 2, vi))) ;

  /* Stabilisierung der Konvektion ( L_conv_u) */
  eforce_(dUDF) += -(ttimetauM*conv_c_(vi)*conv_old_(0)) ;
  eforce_(dVDF) += -(ttimetauM*conv_c_(vi)*conv_old_(1)) ;
  eforce_(dWDF) += -(ttimetauM*conv_c_(vi)*conv_old_(2)) ;

  /* Stabilisierung der Konvektion (-L_visc_u) */
  eforce_(dUDF) += 2.0*nu_*ttimetauM*conv_c_(vi)*visc_old_(0) ;
  eforce_(dVDF) += 2.0*nu_*ttimetauM*conv_c_(vi)*visc_old_(1) ;
  eforce_(dWDF) += 2.0*nu_*ttimetauM*conv_c_(vi)*visc_old_(2) ;

  /* Stabilisierung der Konvektion ( L_pres_p) */
  eforce_(dUDF) += -(ttimetauM*conv_c_(vi)*gradp_(0)) ;
  eforce_(dVDF) += -(ttimetauM*conv_c_(vi)*gradp_(1)) ;
  eforce_(dWDF) += -(ttimetauM*conv_c_(vi)*gradp_(2)) ;

  /* Viskosit�tsterm */
  eforce_(dUDF) += -(nu_*timefacfac*(2.0*derxyz_(0, vi)*vderxyz_(0, 0) + derxyz_(1, vi)*vderxyz_(0, 1) + derxyz_(1, vi)*vderxyz_(1, 0) + derxyz_(2, vi)*vderxyz_(0, 2) + derxyz_(2, vi)*vderxyz_(2, 0))) ;
  eforce_(dVDF) += -(nu_*timefacfac*(derxyz_(0, vi)*vderxyz_(0, 1) + derxyz_(0, vi)*vderxyz_(1, 0) + 2.0*derxyz_(1, vi)*vderxyz_(1, 1) + derxyz_(2, vi)*vderxyz_(1, 2) + derxyz_(2, vi)*vderxyz_(2, 1))) ;
  eforce_(dWDF) += -(nu_*timefacfac*(derxyz_(0, vi)*vderxyz_(0, 2) + derxyz_(0, vi)*vderxyz_(2, 0) + derxyz_(1, vi)*vderxyz_(1, 2) + derxyz_(1, vi)*vderxyz_(2, 1) + 2.0*derxyz_(2, vi)*vderxyz_(2, 2))) ;

  /* Stabilisierung der Viskosit�t ( L_conv_u) */
  eforce_(dUDF) += -2.0*nu_*ttimetauMp*(conv_old_(0)*viscs2_(0, 0, vi) + conv_old_(1)*viscs2_(0, 1, vi) + conv_old_(2)*viscs2_(0, 2, vi)) ;
  eforce_(dVDF) += -2.0*nu_*ttimetauMp*(conv_old_(0)*viscs2_(0, 1, vi) + conv_old_(1)*viscs2_(1, 1, vi) + conv_old_(2)*viscs2_(1, 2, vi)) ;
  eforce_(dWDF) += -2.0*nu_*ttimetauMp*(conv_old_(0)*viscs2_(0, 2, vi) + conv_old_(1)*viscs2_(1, 2, vi) + conv_old_(2)*viscs2_(2, 2, vi)) ;

  /* Stabilisierung der Viskosit�t (-L_visc_u) */
  eforce_(dUDF) += 4.0*(nu_*nu_)*ttimetauMp*(visc_old_(0)*viscs2_(0, 0, vi) + visc_old_(1)*viscs2_(0, 1, vi) + visc_old_(2)*viscs2_(0, 2, vi)) ;
  eforce_(dVDF) += 4.0*(nu_*nu_)*ttimetauMp*(visc_old_(0)*viscs2_(0, 1, vi) + visc_old_(1)*viscs2_(1, 1, vi) + visc_old_(2)*viscs2_(1, 2, vi)) ;
  eforce_(dWDF) += 4.0*(nu_*nu_)*ttimetauMp*(visc_old_(0)*viscs2_(0, 2, vi) + visc_old_(1)*viscs2_(1, 2, vi) + visc_old_(2)*viscs2_(2, 2, vi)) ;

  /* Stabilisierung der Viskosit�t ( L_pres_p) */
  eforce_(dUDF) += -2.0*nu_*ttimetauMp*(gradp_(0)*viscs2_(0, 0, vi) + gradp_(1)*viscs2_(0, 1, vi) + gradp_(2)*viscs2_(0, 2, vi)) ;
  eforce_(dVDF) += -2.0*nu_*ttimetauMp*(gradp_(0)*viscs2_(0, 1, vi) + gradp_(1)*viscs2_(1, 1, vi) + gradp_(2)*viscs2_(1, 2, vi)) ;
  eforce_(dWDF) += -2.0*nu_*ttimetauMp*(gradp_(0)*viscs2_(0, 2, vi) + gradp_(1)*viscs2_(1, 2, vi) + gradp_(2)*viscs2_(2, 2, vi)) ;

  /* Druckterm */
  eforce_(dUDF) += press*timefacfac*derxyz_(0, vi) ;
  eforce_(dVDF) += press*timefacfac*derxyz_(1, vi) ;
  eforce_(dWDF) += press*timefacfac*derxyz_(2, vi) ;

  /* Stabilisierung des Drucks ( L_conv_u) */
  eforce_(dPDF) += -(ttimetauMp*(conv_old_(0)*derxyz_(0, vi) + conv_old_(1)*derxyz_(1, vi) + conv_old_(2)*derxyz_(2, vi))) ;

  /* Stabilisierung des Drucks (-L_visc_u) */
  eforce_(dPDF) += 2.0*nu_*ttimetauMp*(visc_old_(0)*derxyz_(0, vi) + visc_old_(1)*derxyz_(1, vi) + visc_old_(2)*derxyz_(2, vi)) ;

  /* Stabilisierung des Drucks ( L_pres_p) */
  eforce_(dPDF) += -(ttimetauMp*(gradp_(0)*derxyz_(0, vi) + gradp_(1)*derxyz_(1, vi) + gradp_(2)*derxyz_(2, vi))) ;

  /* Divergenzfreiheit */
  eforce_(dPDF) += -(timefacfac*(conv_r_(0, 0, vi) + conv_r_(1, 1, vi) + conv_r_(2, 2, vi))) ;

  /* Kontinuit�tsstabilisierung */
  eforce_(dUDF) += -((thsl*thsl)*tau_C*derxyz_(0, vi)*(vderxyz_(0, 0) + vderxyz_(1, 1) + vderxyz_(2, 2))) ;
  eforce_(dVDF) += -((thsl*thsl)*tau_C*derxyz_(1, vi)*(vderxyz_(0, 0) + vderxyz_(1, 1) + vderxyz_(2, 2))) ;
  eforce_(dWDF) += -((thsl*thsl)*tau_C*derxyz_(2, vi)*(vderxyz_(0, 0) + vderxyz_(1, 1) + vderxyz_(2, 2))) ;

  /* Massenterm */
  eforce_(dUDF) += -(fac*funct_(vi)*velint_(0)) ;
  eforce_(dVDF) += -(fac*funct_(vi)*velint_(1)) ;
  eforce_(dWDF) += -(fac*funct_(vi)*velint_(2)) ;

  /* Konvektionsstabilisierung */
  eforce_(dUDF) += -(timetauM*conv_c_(vi)*velint_(0)) ;
  eforce_(dVDF) += -(timetauM*conv_c_(vi)*velint_(1)) ;
  eforce_(dWDF) += -(timetauM*conv_c_(vi)*velint_(2)) ;

  /* Viskosit�tsstabilisierung */
  eforce_(dUDF) += -2.0*nu_*timetauMp*(velint_(0)*viscs2_(0, 0, vi) + velint_(1)*viscs2_(0, 1, vi) + velint_(2)*viscs2_(0, 2, vi)) ;
  eforce_(dVDF) += -2.0*nu_*timetauMp*(velint_(0)*viscs2_(0, 1, vi) + velint_(1)*viscs2_(1, 1, vi) + velint_(2)*viscs2_(1, 2, vi)) ;
  eforce_(dWDF) += -2.0*nu_*timetauMp*(velint_(0)*viscs2_(0, 2, vi) + velint_(1)*viscs2_(1, 2, vi) + velint_(2)*viscs2_(2, 2, vi)) ;

  /* Stabilisierung der Druckgleichung */
  eforce_(dPDF) += -(timetauMp*conv_c_(vi)) ;

  /* Quellterm der rechten Seite */
  eforce_(dUDF) += fac*funct_(vi)*rhsint_(0) ;
  eforce_(dVDF) += fac*funct_(vi)*rhsint_(1) ;
  eforce_(dWDF) += fac*funct_(vi)*rhsint_(2) ;

  /* Konvektionsstabilisierung */
  eforce_(dUDF) += timetauM*conv_c_(vi)*rhsint_(0) ;
  eforce_(dVDF) += timetauM*conv_c_(vi)*rhsint_(1) ;
  eforce_(dWDF) += timetauM*conv_c_(vi)*rhsint_(2) ;

  /* Viskosit�tsstabilisierung */
  eforce_(dUDF) += 2.0*nu_*timetauMp*(rhsint_(0)*viscs2_(0, 0, vi) + rhsint_(1)*viscs2_(0, 1, vi) + rhsint_(2)*viscs2_(0, 2, vi)) ;
  eforce_(dVDF) += 2.0*nu_*timetauMp*(rhsint_(0)*viscs2_(0, 1, vi) + rhsint_(1)*viscs2_(1, 1, vi) + rhsint_(2)*viscs2_(1, 2, vi)) ;
  eforce_(dWDF) += 2.0*nu_*timetauMp*(rhsint_(0)*viscs2_(0, 2, vi) + rhsint_(1)*viscs2_(1, 2, vi) + rhsint_(2)*viscs2_(2, 2, vi)) ;

  /* Stabilisierung der Druckgleichung */
  eforce_(dPDF) += timetauMp*(rhsint_(0)*derxyz_(0, vi) + rhsint_(1)*derxyz_(1, vi) + rhsint_(2)*derxyz_(2, vi)) ;

}
