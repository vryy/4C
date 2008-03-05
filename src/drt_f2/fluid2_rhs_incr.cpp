for (int vi=0; vi<iel; ++vi)
{
  /* standard Galerkin terms: */
  /* convective term */
  eforce_(vi*3)     += -(timefacfac*funct_(vi)*conv_old_(0)) ;
  eforce_(vi*3 + 1) += -(timefacfac*funct_(vi)*conv_old_(1)) ;

  /* viscous term */
  eforce_(vi*3)     += -(nu_*timefacfac*(2.0*derxy_(0, vi)*vderxy_(0, 0) + derxy_(1, vi)*vderxy_(0, 1) + derxy_(1, vi)*vderxy_(1, 0))) ;
  eforce_(vi*3 + 1) += -(nu_*timefacfac*(derxy_(0, vi)*vderxy_(0, 1) + derxy_(0, vi)*vderxy_(1, 0) + 2.0*derxy_(1, vi)*vderxy_(1, 1))) ;

  /* pressure term */
  eforce_(vi*3)     += press*timefacfac*derxy_(0, vi) ;
  eforce_(vi*3 + 1) += press*timefacfac*derxy_(1, vi) ;

  /* transient term */
  eforce_(vi*3)     += -(fac*funct_(vi)*velint_(0)) ;
  eforce_(vi*3 + 1) += -(fac*funct_(vi)*velint_(1)) ;

  /* continuity term */
  eforce_(vi*3 + 2) += -(timefacfac*(conv_r_(0, 0, vi) + conv_r_(1, 1, vi))) ;

  /* stabilization terms: */
  /* convective stabilization: */
  /* convective term */
  eforce_(vi*3)     += -(ttimetauM*conv_c_(vi)*conv_old_(0)) ;
  eforce_(vi*3 + 1) += -(ttimetauM*conv_c_(vi)*conv_old_(1)) ;

  /* viscous term */
  eforce_(vi*3)     += 2.0*nu_*ttimetauM*conv_c_(vi)*visc_old_(0) ;
  eforce_(vi*3 + 1) += 2.0*nu_*ttimetauM*conv_c_(vi)*visc_old_(1) ;

  /* pressure term */
  eforce_(vi*3)     += -(ttimetauM*conv_c_(vi)*gradp_(0)) ;
  eforce_(vi*3 + 1) += -(ttimetauM*conv_c_(vi)*gradp_(1)) ;

  /* transient term */
  eforce_(vi*3)     += -(timetauM*conv_c_(vi)*velint_(0)) ;
  eforce_(vi*3 + 1) += -(timetauM*conv_c_(vi)*velint_(1)) ;

  /* viscous stabilization: */
  /* convective term */
  /*eforce_(vi*3)     += -2.0*nu_*ttimetauMp*(conv_old_(0)*viscs2_(0, 0, vi) + conv_old_(1)*viscs2_(0, 1, vi)) ;
  eforce_(vi*3 + 1) += -2.0*nu_*ttimetauMp*(conv_old_(0)*viscs2_(0, 1, vi) + conv_old_(1)*viscs2_(1, 1, vi)) ;*/

  /* viscous term */
  /*eforce_(vi*3)     += 4.0*(nu_*nu_)*ttimetauMp*(visc_old_(0)*viscs2_(0, 0, vi) + visc_old_(1)*viscs2_(0, 1, vi)) ;
  eforce_(vi*3 + 1) += 4.0*(nu_*nu_)*ttimetauMp*(visc_old_(0)*viscs2_(0, 1, vi) + visc_old_(1)*viscs2_(1, 1, vi)) ;*/

  /* pressure term */
  /*eforce_(vi*3)     += -2.0*nu_*ttimetauMp*(gradp_(0)*viscs2_(0, 0, vi) + gradp_(1)*viscs2_(0, 1, vi)) ;
  eforce_(vi*3 + 1) += -2.0*nu_*ttimetauMp*(gradp_(0)*viscs2_(0, 1, vi) + gradp_(1)*viscs2_(1, 1, vi)) ;*/

  /* transient term */
  /*eforce_(vi*3)     += -2.0*nu_*timetauMp*(velint_(0)*viscs2_(0, 0, vi) + velint_(1)*viscs2_(0, 1, vi)) ;
  eforce_(vi*3 + 1) += -2.0*nu_*timetauMp*(velint_(0)*viscs2_(0, 1, vi) + velint_(1)*viscs2_(1, 1, vi)) ;*/

  /* pressure stabilization: */
  /* convective term */
  eforce_(vi*3 + 2) += -(ttimetauMp*(derxy_(0, vi)*conv_old_(0) + derxy_(1, vi)*conv_old_(1))) ;

  /* viscous term */
  eforce_(vi*3 + 2) += 2.0*nu_*ttimetauMp*(derxy_(0, vi)*visc_old_(0) + derxy_(1, vi)*visc_old_(1)) ;

  /* pressure term */
  eforce_(vi*3 + 2) += -(ttimetauMp*(gradp_(0)*derxy_(0, vi) + gradp_(1)*derxy_(1, vi))) ;

  /* transient term(?) */
  eforce_(vi*3 + 2) += -(timetauMp*conv_c_(vi)) ;

  /* continuity stabilization: */
  eforce_(vi*3)     += -((thsl*thsl)*tau_C*derxy_(0, vi)*(vderxy_(0, 0) + vderxy_(1, 1))) ;
  eforce_(vi*3 + 1) += -((thsl*thsl)*tau_C*derxy_(1, vi)*(vderxy_(0, 0) + vderxy_(1, 1))) ;

  /* right-hand-side term */
  eforce_(vi*3)     += fac*funct_(vi)*rhsint_(0) ;
  eforce_(vi*3 + 1) += fac*funct_(vi)*rhsint_(1) ;

  /* stabilization terms: */
  /* convective stabilization: */
  eforce_(vi*3)     += timetauM*conv_c_(vi)*rhsint_(0) ;
  eforce_(vi*3 + 1) += timetauM*conv_c_(vi)*rhsint_(1) ;

  /* viscous stabilization: */
  eforce_(vi*3)     += tau_Mp*time2nue*(rhsint_(0)*viscs2_(0, 0, vi) + rhsint_(1)*viscs2_(0, 1, vi)) ;
  eforce_(vi*3 + 1) += tau_Mp*time2nue*(rhsint_(0)*viscs2_(0, 1, vi) + rhsint_(1)*viscs2_(1, 1, vi)) ;

  /* pressure stabilization: */
  eforce_(vi*3 + 2) += timetauMp*(rhsint_(0)*derxy_(0, vi) + rhsint_(1)*derxy_(1, vi)) ;

  if (fssgv != "No" && fssgv != "scale_similarity")
  {
    /* viscous term */
    eforce_(vi*3)     -= vartfac*(2.0*derxy_(0, vi)*fsvderxy_(0, 0) + derxy_(1, vi)*fsvderxy_(0, 1) + derxy_(1, vi)*fsvderxy_(1, 0)) ;
    eforce_(vi*3 + 1) -= vartfac*(derxy_(0, vi)*fsvderxy_(0, 1) + derxy_(0, vi)*fsvderxy_(1, 0) + 2.0*derxy_(1, vi)*fsvderxy_(1, 1)) ;
  }
}
