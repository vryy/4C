for (vi=0; vi<iel; ++vi)
{
    /* Konvektionsterm (u*grad(u),v) */
    eforce_(vi*2)     += -(timefacfac*funct_(vi)*conv_old_(0)) ;
    eforce_(vi*2 + 1) += -(timefacfac*funct_(vi)*conv_old_(1)) ;

    /* Stabilisierung der Konvektion (u*grad(u),u*grad(v)) */
    eforce_(vi*2)     += -(ttimetauM*conv_c_(vi)*conv_old_(0)) ;
    eforce_(vi*2 + 1) += -(ttimetauM*conv_c_(vi)*conv_old_(1)) ;

    /* Stabilisierung der Konvektion (-2*nu*div(epsilon((u))),u*grad(v)) */
    eforce_(vi*2)     += 2.0*visc_*ttimetauM*conv_c_(vi)*visc_old_(0) ;
    eforce_(vi*2 + 1) += 2.0*visc_*ttimetauM*conv_c_(vi)*visc_old_(1) ;

    /* Stabilisierung der Konvektion (grad(p),u*grad(v)) */
    eforce_(vi*2)     += -(ttimetauM*conv_c_(vi)*gradp_(0)) ;
    eforce_(vi*2 + 1) += -(ttimetauM*conv_c_(vi)*gradp_(1)) ;

    /* Viskositätsterm (2*nu*epsilon(u),epsilon(v)) */
    eforce_(vi*2)     += -(visc_*timefacfac*(2.0*derxy_(0, vi)*vderxy_(0, 0) + derxy_(1, vi)*vderxy_(0, 1) + derxy_(1, vi)*vderxy_(1, 0))) ;
    eforce_(vi*2 + 1) += -(visc_*timefacfac*(derxy_(0, vi)*vderxy_(0, 1) + derxy_(0, vi)*vderxy_(1, 0) + 2.0*derxy_(1, vi)*vderxy_(1, 1))) ;

    /* Stabilisierung der Viskosität (u*grad(u),-2*nu*div(epsilon((v)))) */
    eforce_(vi*2)     += -2.0*visc_*ttimetauMp*(conv_old_(0)*viscs2_(0, 0, vi) + conv_old_(1)*viscs2_(0, 1, vi)) ;
    eforce_(vi*2 + 1) += -2.0*visc_*ttimetauMp*(conv_old_(0)*viscs2_(0, 1, vi) + conv_old_(1)*viscs2_(1, 1, vi)) ;

    /* Stabilisierung der Viskosität (-2*nu*div(epsilon((u))),-2*nu*div(epsilon((v)))) */
    eforce_(vi*2)     += 4.0*(visc_*visc_)*ttimetauMp*(visc_old_(0)*viscs2_(0, 0, vi) + visc_old_(1)*viscs2_(0, 1, vi)) ;
    eforce_(vi*2 + 1) += 4.0*(visc_*visc_)*ttimetauMp*(visc_old_(0)*viscs2_(0, 1, vi) + visc_old_(1)*viscs2_(1, 1, vi)) ;

    /* Stabilisierung der Viskosität (grad(p),-2*nu*div(epsilon((v)))) */
    eforce_(vi*2)     += -2.0*visc_*ttimetauMp*(gradp_(0)*viscs2_(0, 0, vi) + gradp_(1)*viscs2_(0, 1, vi)) ;
    eforce_(vi*2 + 1) += -2.0*visc_*ttimetauMp*(gradp_(0)*viscs2_(0, 1, vi) + gradp_(1)*viscs2_(1, 1, vi)) ;

    /* Druckterm */
    eforce_(vi*2)     += press*timefacfac*derxy_(0, vi) ;
    eforce_(vi*2 + 1) += press*timefacfac*derxy_(1, vi) ;

    /* Kontinuitätsstabilisierung */
    eforce_(vi*2)     += -((thsl*thsl)*tau_C*derxy_(0, vi)*(vderxy_(0, 0) + vderxy_(1, 1))) ;
    eforce_(vi*2 + 1) += -((thsl*thsl)*tau_C*derxy_(1, vi)*(vderxy_(0, 0) + vderxy_(1, 1))) ;

    /* Massenterm (u,v) */
    eforce_(vi*2)     += -(fac*funct_(vi)*velint_(0)) ;
    eforce_(vi*2 + 1) += -(fac*funct_(vi)*velint_(1)) ;

    /* Konvektionsstabilisierung (u,u*grad(v)) */
    eforce_(vi*2)     += -(timetauM*conv_c_(vi)*velint_(0)) ;
    eforce_(vi*2 + 1) += -(timetauM*conv_c_(vi)*velint_(1)) ;

    /* Viskositätsstabilisierung (u,-2*nu*div(epsilon((v)))) */
    eforce_(vi*2)     += -2.0*visc_*timetauMp*(velint_(0)*viscs2_(0, 0, vi) + velint_(1)*viscs2_(0, 1, vi)) ;
    eforce_(vi*2 + 1) += -2.0*visc_*timetauMp*(velint_(0)*viscs2_(0, 1, vi) + velint_(1)*viscs2_(1, 1, vi)) ;

    /* Quellterm der rechten Seite (b,v) */
    eforce_(vi*2)     += fac*funct_(vi)*rhsint_(0) ;
    eforce_(vi*2 + 1) += fac*funct_(vi)*rhsint_(1) ;

    /* Konvektionsstabilisierung (b,u*grad(v)) */
    eforce_(vi*2)     += timetauM*conv_c_(vi)*rhsint_(0) ;
    eforce_(vi*2 + 1) += timetauM*conv_c_(vi)*rhsint_(1) ;

    /* Viskositätsstabilisierung (b,-2*nu*div(epsilon((v))) */
    eforce_(vi*2)     += tau_Mp*time2nue*(rhsint_(0)*viscs2_(0, 0, vi) + rhsint_(1)*viscs2_(0, 1, vi)) ;
    eforce_(vi*2 + 1) += tau_Mp*time2nue*(rhsint_(0)*viscs2_(0, 1, vi) + rhsint_(1)*viscs2_(1, 1, vi)) ;

}
