/*
<pre>
Maintainer: Ulrich Küttler
            kuettler@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/Members/kuettler
            089 - 289-15238
</pre>
 */

for (vi=0; vi<4; ++vi)
{
#ifdef FLUID2_IS_TERM1
    /* Konvektionsterm */
    eforce_(vi*3)     += timefacfac*(-(funct_(vi)*conv_old_(0)) + gridvint_(0)*conv_r_(0, 0, vi) + gridvint_(1)*conv_r_(0, 1, vi)) ;
    eforce_(vi*3 + 1) += -(timefacfac*(funct_(vi)*conv_old_(1) - gridvint_(0)*conv_r_(1, 0, vi) - gridvint_(1)*conv_r_(1, 1, vi))) ;
#endif

#ifdef FLUID2_IS_TERM2
    /* Stabilisierung der Konvektion ( L_conv_u) */
    eforce_(vi*3)     += ttimetauM*(-(conv_c_(vi)*conv_old_(0)) + conv_g_(vi)*gridvint_(0)*vderxy_(0, 0) - derxy_(1, vi)*(gridvint_(1)*gridvint_(1))*vderxy_(0, 1) + 2.0*velint_(0)*derxy_(0, vi)*gridvint_(0)*vderxy_(0, 0) + velint_(0)*derxy_(0, vi)*gridvint_(1)*vderxy_(0, 1) + velint_(0)*derxy_(1, vi)*gridvint_(1)*vderxy_(0, 0) + velint_(1)*derxy_(0, vi)*gridvint_(0)*vderxy_(0, 1) + velint_(1)*derxy_(1, vi)*gridvint_(0)*vderxy_(0, 0) + 2.0*velint_(1)*derxy_(1, vi)*gridvint_(1)*vderxy_(0, 1) - derxy_(0, vi)*gridvint_(0)*gridvint_(1)*vderxy_(0, 1)) ;
    eforce_(vi*3 + 1) += ttimetauM*(-(conv_c_(vi)*conv_old_(1)) + conv_g_(vi)*gridvint_(0)*vderxy_(1, 0) - derxy_(1, vi)*(gridvint_(1)*gridvint_(1))*vderxy_(1, 1) + 2.0*velint_(0)*derxy_(0, vi)*gridvint_(0)*vderxy_(1, 0) + velint_(0)*derxy_(0, vi)*gridvint_(1)*vderxy_(1, 1) + velint_(0)*derxy_(1, vi)*gridvint_(1)*vderxy_(1, 0) + velint_(1)*derxy_(0, vi)*gridvint_(0)*vderxy_(1, 1) + velint_(1)*derxy_(1, vi)*gridvint_(0)*vderxy_(1, 0) + 2.0*velint_(1)*derxy_(1, vi)*gridvint_(1)*vderxy_(1, 1) - derxy_(0, vi)*gridvint_(0)*gridvint_(1)*vderxy_(1, 1)) ;
#endif

#ifdef FLUID2_IS_TERM3
    /* Stabilisierung der Konvektion (-L_visc_u) */
    eforce_(vi*3)     += 2.0*nu_*ttimetauM*(conv_c_(vi) + conv_g_(vi))*visc_old_(0) ;
    eforce_(vi*3 + 1) += 2.0*nu_*ttimetauM*(conv_c_(vi) + conv_g_(vi))*visc_old_(1) ;
#endif

#ifdef FLUID2_IS_TERM4
    /* Viskositätsterm */
    eforce_(vi*3)     += -(nu_*timefacfac*(2.0*derxy_(0, vi)*vderxy_(0, 0) + derxy_(1, vi)*vderxy_(0, 1) + derxy_(1, vi)*vderxy_(1, 0))) ;
    eforce_(vi*3 + 1) += -(nu_*timefacfac*(derxy_(0, vi)*vderxy_(0, 1) + derxy_(0, vi)*vderxy_(1, 0) + 2.0*derxy_(1, vi)*vderxy_(1, 1))) ;
#endif

#ifdef FLUID2_IS_TERM5
    /* Stabilisierung der Viskosität ( L_conv_u) */
    eforce_(vi*3)     += 2.0*nu_*ttimetauMp*(-(conv_old_(0)*viscs2_(0, 0, vi)) - conv_old_(1)*viscs2_(0, 1, vi) + gridvint_(0)*vderxy_(0, 0)*viscs2_(0, 0, vi) + gridvint_(0)*vderxy_(1, 0)*viscs2_(0, 1, vi) + gridvint_(1)*vderxy_(0, 1)*viscs2_(0, 0, vi) + gridvint_(1)*vderxy_(1, 1)*viscs2_(0, 1, vi)) ;
    eforce_(vi*3 + 1) += 2.0*nu_*ttimetauMp*(-(conv_old_(0)*viscs2_(0, 1, vi)) - conv_old_(1)*viscs2_(1, 1, vi) + gridvint_(0)*vderxy_(0, 0)*viscs2_(0, 1, vi) + gridvint_(0)*vderxy_(1, 0)*viscs2_(1, 1, vi) + gridvint_(1)*vderxy_(0, 1)*viscs2_(0, 1, vi) + gridvint_(1)*vderxy_(1, 1)*viscs2_(1, 1, vi)) ;
#endif

#ifdef FLUID2_IS_TERM6
    /* Stabilisierung der Viskosität (-L_visc_u) */
    eforce_(vi*3)     += 4.0*(nu_*nu_)*ttimetauMp*(visc_old_(0)*viscs2_(0, 0, vi) + visc_old_(1)*viscs2_(0, 1, vi)) ;
    eforce_(vi*3 + 1) += 4.0*(nu_*nu_)*ttimetauMp*(visc_old_(0)*viscs2_(0, 1, vi) + visc_old_(1)*viscs2_(1, 1, vi)) ;
#endif

#ifdef FLUID2_IS_TERM7
    /* Druckterm */
    eforce_(vi*3)     += press*timefacfac*derxy_(0, vi) ;
    eforce_(vi*3 + 1) += press*timefacfac*derxy_(1, vi) ;
#endif

#ifdef FLUID2_IS_TERM8
    /* Divergenzfreiheit */
    eforce_(vi*3 + 2) += -(timefacfac*funct_p_(vi)*(vderxy_(0, 0) + vderxy_(1, 1))) ;
#endif

#ifdef FLUID2_IS_TERM9
    /* Kontinuitätsstabilisierung */
    eforce_(vi*3)     += -((thsl*thsl)*tau_C*derxy_(0, vi)*(vderxy_(0, 0) + vderxy_(1, 1))) ;
    eforce_(vi*3 + 1) += -((thsl*thsl)*tau_C*derxy_(1, vi)*(vderxy_(0, 0) + vderxy_(1, 1))) ;
#endif

#ifdef FLUID2_IS_TERM10
    /* Massenterm */
    eforce_(vi*3)     += -(fac*funct_(vi)*velint_(0)) ;
    eforce_(vi*3 + 1) += -(fac*funct_(vi)*velint_(1)) ;
#endif

#ifdef FLUID2_IS_TERM11
    /* Konvektionsstabilisierung */
    eforce_(vi*3)     += -(timetauM*(conv_c_(vi) + conv_g_(vi))*velint_(0)) ;
    eforce_(vi*3 + 1) += -(timetauM*(conv_c_(vi) + conv_g_(vi))*velint_(1)) ;
#endif

#ifdef FLUID2_IS_TERM12
    /* Viskositätsstabilisierung */
    eforce_(vi*3)     += -2.0*nu_*timetauMp*(velint_(0)*viscs2_(0, 0, vi) + velint_(1)*viscs2_(0, 1, vi)) ;
    eforce_(vi*3 + 1) += -2.0*nu_*timetauMp*(velint_(0)*viscs2_(0, 1, vi) + velint_(1)*viscs2_(1, 1, vi)) ;
#endif

#ifdef FLUID2_IS_TERM13
    /* Quellterm der rechten Seite */
    eforce_(vi*3)     += fac*funct_(vi)*rhsint_(0) ;
    eforce_(vi*3 + 1) += fac*funct_(vi)*rhsint_(1) ;
#endif

#ifdef FLUID2_IS_TERM14
    /* Konvektionsstabilisierung */
    eforce_(vi*3)     += timetauM*(conv_c_(vi) + conv_g_(vi))*rhsint_(0) ;
    eforce_(vi*3 + 1) += timetauM*(conv_c_(vi) + conv_g_(vi))*rhsint_(1) ;
#endif

#ifdef FLUID2_IS_TERM15
    /* Viskositätsstabilisierung */
    eforce_(vi*3)     += tau_Mp*time2nue*(rhsint_(0)*viscs2_(0, 0, vi) + rhsint_(1)*viscs2_(0, 1, vi)) ;
    eforce_(vi*3 + 1) += tau_Mp*time2nue*(rhsint_(0)*viscs2_(0, 1, vi) + rhsint_(1)*viscs2_(1, 1, vi)) ;
#endif

}

for (vi=4; vi<iel; ++vi)
{
#ifdef FLUID2_IS_TERM1
    /* Konvektionsterm */
    eforce_(vi*2 + 4) += timefacfac*(-(funct_(vi)*conv_old_(0)) + gridvint_(0)*conv_r_(0, 0, vi) + gridvint_(1)*conv_r_(0, 1, vi)) ;
    eforce_(vi*2 + 5) += -(timefacfac*(funct_(vi)*conv_old_(1) - gridvint_(0)*conv_r_(1, 0, vi) - gridvint_(1)*conv_r_(1, 1, vi))) ;
#endif

#ifdef FLUID2_IS_TERM2
    /* Stabilisierung der Konvektion ( L_conv_u) */
    eforce_(vi*2 + 4) += ttimetauM*(-(conv_c_(vi)*conv_old_(0)) + conv_g_(vi)*gridvint_(0)*vderxy_(0, 0) - derxy_(1, vi)*(gridvint_(1)*gridvint_(1))*vderxy_(0, 1) + 2.0*velint_(0)*derxy_(0, vi)*gridvint_(0)*vderxy_(0, 0) + velint_(0)*derxy_(0, vi)*gridvint_(1)*vderxy_(0, 1) + velint_(0)*derxy_(1, vi)*gridvint_(1)*vderxy_(0, 0) + velint_(1)*derxy_(0, vi)*gridvint_(0)*vderxy_(0, 1) + velint_(1)*derxy_(1, vi)*gridvint_(0)*vderxy_(0, 0) + 2.0*velint_(1)*derxy_(1, vi)*gridvint_(1)*vderxy_(0, 1) - derxy_(0, vi)*gridvint_(0)*gridvint_(1)*vderxy_(0, 1)) ;
    eforce_(vi*2 + 5) += ttimetauM*(-(conv_c_(vi)*conv_old_(1)) + conv_g_(vi)*gridvint_(0)*vderxy_(1, 0) - derxy_(1, vi)*(gridvint_(1)*gridvint_(1))*vderxy_(1, 1) + 2.0*velint_(0)*derxy_(0, vi)*gridvint_(0)*vderxy_(1, 0) + velint_(0)*derxy_(0, vi)*gridvint_(1)*vderxy_(1, 1) + velint_(0)*derxy_(1, vi)*gridvint_(1)*vderxy_(1, 0) + velint_(1)*derxy_(0, vi)*gridvint_(0)*vderxy_(1, 1) + velint_(1)*derxy_(1, vi)*gridvint_(0)*vderxy_(1, 0) + 2.0*velint_(1)*derxy_(1, vi)*gridvint_(1)*vderxy_(1, 1) - derxy_(0, vi)*gridvint_(0)*gridvint_(1)*vderxy_(1, 1)) ;
#endif

#ifdef FLUID2_IS_TERM3
    /* Stabilisierung der Konvektion (-L_visc_u) */
    eforce_(vi*2 + 4) += 2.0*nu_*ttimetauM*(conv_c_(vi) + conv_g_(vi))*visc_old_(0) ;
    eforce_(vi*2 + 5) += 2.0*nu_*ttimetauM*(conv_c_(vi) + conv_g_(vi))*visc_old_(1) ;
#endif

#ifdef FLUID2_IS_TERM4
    /* Viskositätsterm */
    eforce_(vi*2 + 4) += -(nu_*timefacfac*(2.0*derxy_(0, vi)*vderxy_(0, 0) + derxy_(1, vi)*vderxy_(0, 1) + derxy_(1, vi)*vderxy_(1, 0))) ;
    eforce_(vi*2 + 5) += -(nu_*timefacfac*(derxy_(0, vi)*vderxy_(0, 1) + derxy_(0, vi)*vderxy_(1, 0) + 2.0*derxy_(1, vi)*vderxy_(1, 1))) ;
#endif

#ifdef FLUID2_IS_TERM5
    /* Stabilisierung der Viskosität ( L_conv_u) */
    eforce_(vi*2 + 4) += 2.0*nu_*ttimetauMp*(-(conv_old_(0)*viscs2_(0, 0, vi)) - conv_old_(1)*viscs2_(0, 1, vi) + gridvint_(0)*vderxy_(0, 0)*viscs2_(0, 0, vi) + gridvint_(0)*vderxy_(1, 0)*viscs2_(0, 1, vi) + gridvint_(1)*vderxy_(0, 1)*viscs2_(0, 0, vi) + gridvint_(1)*vderxy_(1, 1)*viscs2_(0, 1, vi)) ;
    eforce_(vi*2 + 5) += 2.0*nu_*ttimetauMp*(-(conv_old_(0)*viscs2_(0, 1, vi)) - conv_old_(1)*viscs2_(1, 1, vi) + gridvint_(0)*vderxy_(0, 0)*viscs2_(0, 1, vi) + gridvint_(0)*vderxy_(1, 0)*viscs2_(1, 1, vi) + gridvint_(1)*vderxy_(0, 1)*viscs2_(0, 1, vi) + gridvint_(1)*vderxy_(1, 1)*viscs2_(1, 1, vi)) ;
#endif

#ifdef FLUID2_IS_TERM6
    /* Stabilisierung der Viskosität (-L_visc_u) */
    eforce_(vi*2 + 4) += 4.0*(nu_*nu_)*ttimetauMp*(visc_old_(0)*viscs2_(0, 0, vi) + visc_old_(1)*viscs2_(0, 1, vi)) ;
    eforce_(vi*2 + 5) += 4.0*(nu_*nu_)*ttimetauMp*(visc_old_(0)*viscs2_(0, 1, vi) + visc_old_(1)*viscs2_(1, 1, vi)) ;
#endif

#ifdef FLUID2_IS_TERM7
    /* Druckterm */
    eforce_(vi*2 + 4) += press*timefacfac*derxy_(0, vi) ;
    eforce_(vi*2 + 5) += press*timefacfac*derxy_(1, vi) ;
#endif

#ifdef FLUID2_IS_TERM8
    /* Divergenzfreiheit */
#endif

#ifdef FLUID2_IS_TERM9
    /* Kontinuitätsstabilisierung */
    eforce_(vi*2 + 4) += -((thsl*thsl)*tau_C*derxy_(0, vi)*(vderxy_(0, 0) + vderxy_(1, 1))) ;
    eforce_(vi*2 + 5) += -((thsl*thsl)*tau_C*derxy_(1, vi)*(vderxy_(0, 0) + vderxy_(1, 1))) ;
#endif

#ifdef FLUID2_IS_TERM10
    /* Massenterm */
    eforce_(vi*2 + 4) += -(fac*funct_(vi)*velint_(0)) ;
    eforce_(vi*2 + 5) += -(fac*funct_(vi)*velint_(1)) ;
#endif

#ifdef FLUID2_IS_TERM11
    /* Konvektionsstabilisierung */
    eforce_(vi*2 + 4) += -(timetauM*(conv_c_(vi) + conv_g_(vi))*velint_(0)) ;
    eforce_(vi*2 + 5) += -(timetauM*(conv_c_(vi) + conv_g_(vi))*velint_(1)) ;
#endif

#ifdef FLUID2_IS_TERM12
    /* Viskositätsstabilisierung */
    eforce_(vi*2 + 4) += -2.0*nu_*timetauMp*(velint_(0)*viscs2_(0, 0, vi) + velint_(1)*viscs2_(0, 1, vi)) ;
    eforce_(vi*2 + 5) += -2.0*nu_*timetauMp*(velint_(0)*viscs2_(0, 1, vi) + velint_(1)*viscs2_(1, 1, vi)) ;
#endif

#ifdef FLUID2_IS_TERM13
    /* Quellterm der rechten Seite */
    eforce_(vi*2 + 4) += fac*funct_(vi)*rhsint_(0) ;
    eforce_(vi*2 + 5) += fac*funct_(vi)*rhsint_(1) ;
#endif

#ifdef FLUID2_IS_TERM14
    /* Konvektionsstabilisierung */
    eforce_(vi*2 + 4) += timetauM*(conv_c_(vi) + conv_g_(vi))*rhsint_(0) ;
    eforce_(vi*2 + 5) += timetauM*(conv_c_(vi) + conv_g_(vi))*rhsint_(1) ;
#endif

#ifdef FLUID2_IS_TERM15
    /* Viskositätsstabilisierung */
    eforce_(vi*2 + 4) += tau_Mp*time2nue*(rhsint_(0)*viscs2_(0, 0, vi) + rhsint_(1)*viscs2_(0, 1, vi)) ;
    eforce_(vi*2 + 5) += tau_Mp*time2nue*(rhsint_(0)*viscs2_(0, 1, vi) + rhsint_(1)*viscs2_(1, 1, vi)) ;
#endif

}
