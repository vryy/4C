for (vi=0; vi<iel; ++vi)
{
    /* Konvektionsterm (u*grad(u),v) */
    eforce_(vi*2)     += timefacfac*funct_(vi)*conv_old_(0) ;
    eforce_(vi*2 + 1) += timefacfac*funct_(vi)*conv_old_(1) ;

    /* Stabilisierung der Konvektion (u*grad(u),u*grad(v)) */
    eforce_(vi*2)     += 2.0*ttimetauM*conv_c_(vi)*conv_old_(0) ;
    eforce_(vi*2 + 1) += 2.0*ttimetauM*conv_c_(vi)*conv_old_(1) ;

    /* Stabilisierung der Konvektion (-2*nu*div(epsilon((u))),u*grad(v)) */
    eforce_(vi*2)     += -2.0*visc_*ttimetauM*conv_c_(vi)*visc_old_(0) ;
    eforce_(vi*2 + 1) += -2.0*visc_*ttimetauM*conv_c_(vi)*visc_old_(1) ;

    /* Stabilisierung der Konvektion (grad(p),u*grad(v)) */

    /* Viskositätsterm (2*nu*epsilon(u),epsilon(v)) */

    /* Stabilisierung der Viskosität (u*grad(u),-2*nu*div(epsilon((v)))) */
    eforce_(vi*2)     += 2.0*visc_*ttimetauMp*(conv_old_(0)*viscs2_(0, 0, vi) + conv_old_(1)*viscs2_(0, 1, vi)) ;
    eforce_(vi*2 + 1) += 2.0*visc_*ttimetauMp*(conv_old_(0)*viscs2_(0, 1, vi) + conv_old_(1)*viscs2_(1, 1, vi)) ;

    /* Stabilisierung der Viskosität (-2*nu*div(epsilon((u))),-2*nu*div(epsilon((v)))) */

    /* Stabilisierung der Viskosität (grad(p),-2*nu*div(epsilon((v)))) */
    eforce_(vi*2)     += -2.0*visc_*ttimetauMp*(gradp_(0)*viscs2_(0, 0, vi) + gradp_(1)*viscs2_(0, 1, vi)) ;
    eforce_(vi*2 + 1) += -2.0*visc_*ttimetauMp*(gradp_(0)*viscs2_(0, 1, vi) + gradp_(1)*viscs2_(1, 1, vi)) ;

    /* Kontinuitätsstabilisierung */

    /* Massenterm (u,v) */

    /* Konvektionsstabilisierung (u,u*grad(v)) */
    eforce_(vi*2)     += timetauM*conv_c_(vi)*velint_(0) ;
    eforce_(vi*2 + 1) += timetauM*conv_c_(vi)*velint_(1) ;

    /* Viskositätsstabilisierung (u,-2*nu*div(epsilon((v)))) */

    /* Quellterm der rechten Seite (b,v) */
    eforce_(vi*2)     += fac*funct_(vi)*rhsint_(0) ;
    eforce_(vi*2 + 1) += fac*funct_(vi)*rhsint_(1) ;

    /* Konvektionsstabilisierung (b,u*grad(v)) */

    /* Viskositätsstabilisierung (b,-2*nu*div(epsilon((v))) */
    eforce_(vi*2)     += tau_Mp*time2nue*(rhsint_(0)*viscs2_(0, 0, vi) + rhsint_(1)*viscs2_(0, 1, vi)) ;
    eforce_(vi*2 + 1) += tau_Mp*time2nue*(rhsint_(0)*viscs2_(0, 1, vi) + rhsint_(1)*viscs2_(1, 1, vi)) ;

}
