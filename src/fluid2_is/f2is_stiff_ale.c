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
  for (ui=0; ui<4; ++ui)
  {
#ifdef FLUID2_IS_TERM1
    /* Konvektionsterm */
    estif_(vi*3, ui*3)         += timefacfac*funct_(vi)*(conv_c_(ui) + conv_g_(ui) + conv_r_(0, 0, ui)) ;
    estif_(vi*3, ui*3 + 1)     += timefacfac*funct_(vi)*conv_r_(0, 1, ui) ;
    estif_(vi*3 + 1, ui*3)     += timefacfac*funct_(vi)*conv_r_(1, 0, ui) ;
    estif_(vi*3 + 1, ui*3 + 1) += timefacfac*funct_(vi)*(conv_c_(ui) + conv_g_(ui) + conv_r_(1, 1, ui)) ;
#endif

#ifdef FLUID2_IS_TERM2
    /* Stabilisierung der Konvektion ( L_conv_u) */
    estif_(vi*3, ui*3)         += ttimetauM*(conv_c_(ui)*conv_c_(vi) - conv_g_(ui)*derxy_(0, vi)*gridvint_(0) + 2.0*velint_(0)*derxy_(0, vi)*conv_r_(0, 0, ui) + velint_(1)*derxy_(0, vi)*conv_r_(0, 1, ui) + velint_(1)*derxy_(1, vi)*conv_r_(0, 0, ui) - 2.0*derxy_(0, vi)*gridvint_(0)*conv_r_(0, 0, ui) - derxy_(0, vi)*gridvint_(1)*conv_r_(0, 1, ui) - derxy_(1, vi)*gridvint_(1)*conv_r_(0, 0, ui) + derxy_(1, ui)*derxy_(1, vi)*(gridvint_(1)*gridvint_(1)) - 2.0*velint_(0)*derxy_(0, ui)*derxy_(0, vi)*gridvint_(0) - velint_(0)*derxy_(0, ui)*derxy_(1, vi)*gridvint_(1) - velint_(0)*derxy_(0, vi)*derxy_(1, ui)*gridvint_(1) - velint_(1)*derxy_(0, ui)*derxy_(1, vi)*gridvint_(0) - velint_(1)*derxy_(0, vi)*derxy_(1, ui)*gridvint_(0) - 2.0*velint_(1)*derxy_(1, ui)*derxy_(1, vi)*gridvint_(1) + derxy_(0, ui)*derxy_(1, vi)*gridvint_(0)*gridvint_(1)) ;
    estif_(vi*3, ui*3 + 1)     += ttimetauM*(conv_c_(vi)*conv_r_(0, 1, ui) + velint_(0)*derxy_(1, vi)*conv_r_(0, 0, ui) + velint_(1)*derxy_(1, vi)*conv_r_(0, 1, ui) - derxy_(0, vi)*gridvint_(0)*conv_r_(0, 1, ui) - derxy_(1, vi)*gridvint_(0)*conv_r_(0, 0, ui) - 2.0*derxy_(1, vi)*gridvint_(1)*conv_r_(0, 1, ui)) ;
    estif_(vi*3 + 1, ui*3)     += ttimetauM*(conv_c_(vi)*conv_r_(1, 0, ui) + velint_(0)*derxy_(0, vi)*conv_r_(1, 0, ui) + velint_(1)*derxy_(0, vi)*conv_r_(1, 1, ui) - 2.0*derxy_(0, vi)*gridvint_(0)*conv_r_(1, 0, ui) - derxy_(0, vi)*gridvint_(1)*conv_r_(1, 1, ui) - derxy_(1, vi)*gridvint_(1)*conv_r_(1, 0, ui)) ;
    estif_(vi*3 + 1, ui*3 + 1) += ttimetauM*(conv_c_(ui)*conv_c_(vi) - conv_g_(ui)*derxy_(0, vi)*gridvint_(0) + velint_(0)*derxy_(0, vi)*conv_r_(1, 1, ui) + velint_(0)*derxy_(1, vi)*conv_r_(1, 0, ui) + 2.0*velint_(1)*derxy_(1, vi)*conv_r_(1, 1, ui) - derxy_(0, vi)*gridvint_(0)*conv_r_(1, 1, ui) - derxy_(1, vi)*gridvint_(0)*conv_r_(1, 0, ui) - 2.0*derxy_(1, vi)*gridvint_(1)*conv_r_(1, 1, ui) + derxy_(1, ui)*derxy_(1, vi)*(gridvint_(1)*gridvint_(1)) - 2.0*velint_(0)*derxy_(0, ui)*derxy_(0, vi)*gridvint_(0) - velint_(0)*derxy_(0, ui)*derxy_(1, vi)*gridvint_(1) - velint_(0)*derxy_(0, vi)*derxy_(1, ui)*gridvint_(1) - velint_(1)*derxy_(0, ui)*derxy_(1, vi)*gridvint_(0) - velint_(1)*derxy_(0, vi)*derxy_(1, ui)*gridvint_(0) - 2.0*velint_(1)*derxy_(1, ui)*derxy_(1, vi)*gridvint_(1) + derxy_(0, ui)*derxy_(1, vi)*gridvint_(0)*gridvint_(1)) ;
#endif

#ifdef FLUID2_IS_TERM3
    /* Stabilisierung der Konvektion (-L_visc_u) */
    estif_(vi*3, ui*3)         += -2.0*nu_*ttimetauM*(-(conv_c_(vi)*viscs2_(0, 0, ui)) - conv_g_(vi)*viscs2_(0, 0, ui) + funct_(ui)*derxy_(0, vi)*visc_old_(0)) ;
    estif_(vi*3, ui*3 + 1)     += -2.0*nu_*ttimetauM*(-(conv_c_(vi)*viscs2_(0, 1, ui)) - conv_g_(vi)*viscs2_(0, 1, ui) + funct_(ui)*derxy_(1, vi)*visc_old_(0)) ;
    estif_(vi*3 + 1, ui*3)     += -2.0*nu_*ttimetauM*(-(conv_c_(vi)*viscs2_(0, 1, ui)) - conv_g_(vi)*viscs2_(0, 1, ui) + funct_(ui)*derxy_(0, vi)*visc_old_(1)) ;
    estif_(vi*3 + 1, ui*3 + 1) += -2.0*nu_*ttimetauM*(-(conv_c_(vi)*viscs2_(1, 1, ui)) - conv_g_(vi)*viscs2_(1, 1, ui) + funct_(ui)*derxy_(1, vi)*visc_old_(1)) ;
#endif

#ifdef FLUID2_IS_TERM4
    /* Viskositätsterm */
    estif_(vi*3, ui*3)         += fac*time2nue*derxy_(0, ui)*derxy_(0, vi) + nu_*timefacfac*derxy_(1, ui)*derxy_(1, vi) ;
    estif_(vi*3, ui*3 + 1)     += nu_*timefacfac*derxy_(0, ui)*derxy_(1, vi) ;
    estif_(vi*3 + 1, ui*3)     += nu_*timefacfac*derxy_(0, vi)*derxy_(1, ui) ;
    estif_(vi*3 + 1, ui*3 + 1) += fac*time2nue*derxy_(1, ui)*derxy_(1, vi) + nu_*timefacfac*derxy_(0, ui)*derxy_(0, vi) ;
#endif

#ifdef FLUID2_IS_TERM5
    /* Stabilisierung der Viskosität ( L_conv_u) */
    estif_(vi*3, ui*3)         += 2.0*nu_*ttimetauMp*(conv_c_(ui)*viscs2_(0, 0, vi) + conv_g_(ui)*viscs2_(0, 0, vi) + viscs2_(0, 0, vi)*conv_r_(0, 0, ui) + viscs2_(0, 1, vi)*conv_r_(1, 0, ui)) ;
    estif_(vi*3, ui*3 + 1)     += 2.0*nu_*ttimetauMp*(conv_c_(ui)*viscs2_(0, 1, vi) + conv_g_(ui)*viscs2_(0, 1, vi) + viscs2_(0, 0, vi)*conv_r_(0, 1, ui) + viscs2_(0, 1, vi)*conv_r_(1, 1, ui)) ;
    estif_(vi*3 + 1, ui*3)     += 2.0*nu_*ttimetauMp*(conv_c_(ui)*viscs2_(0, 1, vi) + conv_g_(ui)*viscs2_(0, 1, vi) + viscs2_(0, 1, vi)*conv_r_(0, 0, ui) + viscs2_(1, 1, vi)*conv_r_(1, 0, ui)) ;
    estif_(vi*3 + 1, ui*3 + 1) += 2.0*nu_*ttimetauMp*(conv_c_(ui)*viscs2_(1, 1, vi) + conv_g_(ui)*viscs2_(1, 1, vi) + viscs2_(0, 1, vi)*conv_r_(0, 1, ui) + viscs2_(1, 1, vi)*conv_r_(1, 1, ui)) ;
#endif

#ifdef FLUID2_IS_TERM6
    /* Stabilisierung der Viskosität (-L_visc_u) */
    estif_(vi*3, ui*3)         += 4.0*(nu_*nu_)*ttimetauMp*(viscs2_(0, 0, ui)*viscs2_(0, 0, vi) + viscs2_(0, 1, ui)*viscs2_(0, 1, vi)) ;
    estif_(vi*3, ui*3 + 1)     += 4.0*(nu_*nu_)*ttimetauMp*(viscs2_(0, 0, vi)*viscs2_(0, 1, ui) + viscs2_(0, 1, vi)*viscs2_(1, 1, ui)) ;
    estif_(vi*3 + 1, ui*3)     += 4.0*(nu_*nu_)*ttimetauMp*(viscs2_(0, 0, ui)*viscs2_(0, 1, vi) + viscs2_(0, 1, ui)*viscs2_(1, 1, vi)) ;
    estif_(vi*3 + 1, ui*3 + 1) += 4.0*(nu_*nu_)*ttimetauMp*(viscs2_(0, 1, ui)*viscs2_(0, 1, vi) + viscs2_(1, 1, ui)*viscs2_(1, 1, vi)) ;
#endif

#ifdef FLUID2_IS_TERM7
    /* Druckterm */
    estif_(vi*3, ui*3 + 2)     += -(timefacfac*funct_p_(ui)*derxy_(0, vi)) ;
    estif_(vi*3 + 1, ui*3 + 2) += -(timefacfac*funct_p_(ui)*derxy_(1, vi)) ;
#endif

#ifdef FLUID2_IS_TERM8
    /* Divergenzfreiheit */
    estif_(vi*3 + 2, ui*3)     += timefacfac*funct_p_(vi)*derxy_(0, ui) ;
    estif_(vi*3 + 2, ui*3 + 1) += timefacfac*funct_p_(vi)*derxy_(1, ui) ;
#endif

#ifdef FLUID2_IS_TERM9
    /* Kontinuitätsstabilisierung */
    estif_(vi*3, ui*3)         += (thsl*thsl)*tau_C*derxy_(0, ui)*derxy_(0, vi) ;
    estif_(vi*3, ui*3 + 1)     += (thsl*thsl)*tau_C*derxy_(0, vi)*derxy_(1, ui) ;
    estif_(vi*3 + 1, ui*3)     += (thsl*thsl)*tau_C*derxy_(0, ui)*derxy_(1, vi) ;
    estif_(vi*3 + 1, ui*3 + 1) += (thsl*thsl)*tau_C*derxy_(1, ui)*derxy_(1, vi) ;
#endif

#ifdef FLUID2_IS_TERM10
    /* Massenterm */
    estif_(vi*3, ui*3)         += fac*funct_(ui)*funct_(vi) ;
    estif_(vi*3 + 1, ui*3 + 1) += fac*funct_(ui)*funct_(vi) ;
#endif

#ifdef FLUID2_IS_TERM11
    /* Konvektionsstabilisierung */
    estif_(vi*3, ui*3)         += timetauM*funct_(ui)*(conv_c_(vi) + conv_g_(vi) + velint_(0)*derxy_(0, vi)) ;
    estif_(vi*3, ui*3 + 1)     += timetauM*funct_(ui)*velint_(0)*derxy_(1, vi) ;
    estif_(vi*3 + 1, ui*3)     += timetauM*funct_(ui)*velint_(1)*derxy_(0, vi) ;
    estif_(vi*3 + 1, ui*3 + 1) += timetauM*funct_(ui)*(conv_c_(vi) + conv_g_(vi) + velint_(1)*derxy_(1, vi)) ;
#endif

#ifdef FLUID2_IS_TERM12
    /* Viskositätsstabilisierung */
    estif_(vi*3, ui*3)         += tau_Mp*time2nue*funct_(ui)*viscs2_(0, 0, vi) ;
    estif_(vi*3, ui*3 + 1)     += tau_Mp*time2nue*funct_(ui)*viscs2_(0, 1, vi) ;
    estif_(vi*3 + 1, ui*3)     += tau_Mp*time2nue*funct_(ui)*viscs2_(0, 1, vi) ;
    estif_(vi*3 + 1, ui*3 + 1) += tau_Mp*time2nue*funct_(ui)*viscs2_(1, 1, vi) ;
#endif

#ifdef FLUID2_IS_TERM13
    /* Quellterm der rechten Seite */
#endif

#ifdef FLUID2_IS_TERM14
    /* Konvektionsstabilisierung */
    estif_(vi*3, ui*3)         += -(timetauM*funct_(ui)*rhsint_(0)*derxy_(0, vi)) ;
    estif_(vi*3, ui*3 + 1)     += -(timetauM*funct_(ui)*rhsint_(0)*derxy_(1, vi)) ;
    estif_(vi*3 + 1, ui*3)     += -(timetauM*funct_(ui)*rhsint_(1)*derxy_(0, vi)) ;
    estif_(vi*3 + 1, ui*3 + 1) += -(timetauM*funct_(ui)*rhsint_(1)*derxy_(1, vi)) ;
#endif

#ifdef FLUID2_IS_TERM15
    /* Viskositätsstabilisierung */
#endif

  }
}

for (vi=0; vi<4; ++vi)
{
  for (ui=4; ui<9; ++ui)
  {
#ifdef FLUID2_IS_TERM1
    /* Konvektionsterm */
    estif_(vi*3, ui*2 + 4)     += timefacfac*funct_(vi)*(conv_c_(ui) + conv_g_(ui) + conv_r_(0, 0, ui)) ;
    estif_(vi*3, ui*2 + 5)     += timefacfac*funct_(vi)*conv_r_(0, 1, ui) ;
    estif_(vi*3 + 1, ui*2 + 4) += timefacfac*funct_(vi)*conv_r_(1, 0, ui) ;
    estif_(vi*3 + 1, ui*2 + 5) += timefacfac*funct_(vi)*(conv_c_(ui) + conv_g_(ui) + conv_r_(1, 1, ui)) ;
#endif

#ifdef FLUID2_IS_TERM2
    /* Stabilisierung der Konvektion ( L_conv_u) */
    estif_(vi*3, ui*2 + 4)     += ttimetauM*(conv_c_(ui)*conv_c_(vi) - conv_g_(ui)*derxy_(0, vi)*gridvint_(0) + 2.0*velint_(0)*derxy_(0, vi)*conv_r_(0, 0, ui) + velint_(1)*derxy_(0, vi)*conv_r_(0, 1, ui) + velint_(1)*derxy_(1, vi)*conv_r_(0, 0, ui) - 2.0*derxy_(0, vi)*gridvint_(0)*conv_r_(0, 0, ui) - derxy_(0, vi)*gridvint_(1)*conv_r_(0, 1, ui) - derxy_(1, vi)*gridvint_(1)*conv_r_(0, 0, ui) + derxy_(1, ui)*derxy_(1, vi)*(gridvint_(1)*gridvint_(1)) - 2.0*velint_(0)*derxy_(0, ui)*derxy_(0, vi)*gridvint_(0) - velint_(0)*derxy_(0, ui)*derxy_(1, vi)*gridvint_(1) - velint_(0)*derxy_(0, vi)*derxy_(1, ui)*gridvint_(1) - velint_(1)*derxy_(0, ui)*derxy_(1, vi)*gridvint_(0) - velint_(1)*derxy_(0, vi)*derxy_(1, ui)*gridvint_(0) - 2.0*velint_(1)*derxy_(1, ui)*derxy_(1, vi)*gridvint_(1) + derxy_(0, ui)*derxy_(1, vi)*gridvint_(0)*gridvint_(1)) ;
    estif_(vi*3, ui*2 + 5)     += ttimetauM*(conv_c_(vi)*conv_r_(0, 1, ui) + velint_(0)*derxy_(1, vi)*conv_r_(0, 0, ui) + velint_(1)*derxy_(1, vi)*conv_r_(0, 1, ui) - derxy_(0, vi)*gridvint_(0)*conv_r_(0, 1, ui) - derxy_(1, vi)*gridvint_(0)*conv_r_(0, 0, ui) - 2.0*derxy_(1, vi)*gridvint_(1)*conv_r_(0, 1, ui)) ;
    estif_(vi*3 + 1, ui*2 + 4) += ttimetauM*(conv_c_(vi)*conv_r_(1, 0, ui) + velint_(0)*derxy_(0, vi)*conv_r_(1, 0, ui) + velint_(1)*derxy_(0, vi)*conv_r_(1, 1, ui) - 2.0*derxy_(0, vi)*gridvint_(0)*conv_r_(1, 0, ui) - derxy_(0, vi)*gridvint_(1)*conv_r_(1, 1, ui) - derxy_(1, vi)*gridvint_(1)*conv_r_(1, 0, ui)) ;
    estif_(vi*3 + 1, ui*2 + 5) += ttimetauM*(conv_c_(ui)*conv_c_(vi) - conv_g_(ui)*derxy_(0, vi)*gridvint_(0) + velint_(0)*derxy_(0, vi)*conv_r_(1, 1, ui) + velint_(0)*derxy_(1, vi)*conv_r_(1, 0, ui) + 2.0*velint_(1)*derxy_(1, vi)*conv_r_(1, 1, ui) - derxy_(0, vi)*gridvint_(0)*conv_r_(1, 1, ui) - derxy_(1, vi)*gridvint_(0)*conv_r_(1, 0, ui) - 2.0*derxy_(1, vi)*gridvint_(1)*conv_r_(1, 1, ui) + derxy_(1, ui)*derxy_(1, vi)*(gridvint_(1)*gridvint_(1)) - 2.0*velint_(0)*derxy_(0, ui)*derxy_(0, vi)*gridvint_(0) - velint_(0)*derxy_(0, ui)*derxy_(1, vi)*gridvint_(1) - velint_(0)*derxy_(0, vi)*derxy_(1, ui)*gridvint_(1) - velint_(1)*derxy_(0, ui)*derxy_(1, vi)*gridvint_(0) - velint_(1)*derxy_(0, vi)*derxy_(1, ui)*gridvint_(0) - 2.0*velint_(1)*derxy_(1, ui)*derxy_(1, vi)*gridvint_(1) + derxy_(0, ui)*derxy_(1, vi)*gridvint_(0)*gridvint_(1)) ;
#endif

#ifdef FLUID2_IS_TERM3
    /* Stabilisierung der Konvektion (-L_visc_u) */
    estif_(vi*3, ui*2 + 4)     += -2.0*nu_*ttimetauM*(-(conv_c_(vi)*viscs2_(0, 0, ui)) - conv_g_(vi)*viscs2_(0, 0, ui) + funct_(ui)*derxy_(0, vi)*visc_old_(0)) ;
    estif_(vi*3, ui*2 + 5)     += -2.0*nu_*ttimetauM*(-(conv_c_(vi)*viscs2_(0, 1, ui)) - conv_g_(vi)*viscs2_(0, 1, ui) + funct_(ui)*derxy_(1, vi)*visc_old_(0)) ;
    estif_(vi*3 + 1, ui*2 + 4) += -2.0*nu_*ttimetauM*(-(conv_c_(vi)*viscs2_(0, 1, ui)) - conv_g_(vi)*viscs2_(0, 1, ui) + funct_(ui)*derxy_(0, vi)*visc_old_(1)) ;
    estif_(vi*3 + 1, ui*2 + 5) += -2.0*nu_*ttimetauM*(-(conv_c_(vi)*viscs2_(1, 1, ui)) - conv_g_(vi)*viscs2_(1, 1, ui) + funct_(ui)*derxy_(1, vi)*visc_old_(1)) ;
#endif

#ifdef FLUID2_IS_TERM4
    /* Viskositätsterm */
    estif_(vi*3, ui*2 + 4)     += fac*time2nue*derxy_(0, ui)*derxy_(0, vi) + nu_*timefacfac*derxy_(1, ui)*derxy_(1, vi) ;
    estif_(vi*3, ui*2 + 5)     += nu_*timefacfac*derxy_(0, ui)*derxy_(1, vi) ;
    estif_(vi*3 + 1, ui*2 + 4) += nu_*timefacfac*derxy_(0, vi)*derxy_(1, ui) ;
    estif_(vi*3 + 1, ui*2 + 5) += fac*time2nue*derxy_(1, ui)*derxy_(1, vi) + nu_*timefacfac*derxy_(0, ui)*derxy_(0, vi) ;
#endif

#ifdef FLUID2_IS_TERM5
    /* Stabilisierung der Viskosität ( L_conv_u) */
    estif_(vi*3, ui*2 + 4)     += 2.0*nu_*ttimetauMp*(conv_c_(ui)*viscs2_(0, 0, vi) + conv_g_(ui)*viscs2_(0, 0, vi) + viscs2_(0, 0, vi)*conv_r_(0, 0, ui) + viscs2_(0, 1, vi)*conv_r_(1, 0, ui)) ;
    estif_(vi*3, ui*2 + 5)     += 2.0*nu_*ttimetauMp*(conv_c_(ui)*viscs2_(0, 1, vi) + conv_g_(ui)*viscs2_(0, 1, vi) + viscs2_(0, 0, vi)*conv_r_(0, 1, ui) + viscs2_(0, 1, vi)*conv_r_(1, 1, ui)) ;
    estif_(vi*3 + 1, ui*2 + 4) += 2.0*nu_*ttimetauMp*(conv_c_(ui)*viscs2_(0, 1, vi) + conv_g_(ui)*viscs2_(0, 1, vi) + viscs2_(0, 1, vi)*conv_r_(0, 0, ui) + viscs2_(1, 1, vi)*conv_r_(1, 0, ui)) ;
    estif_(vi*3 + 1, ui*2 + 5) += 2.0*nu_*ttimetauMp*(conv_c_(ui)*viscs2_(1, 1, vi) + conv_g_(ui)*viscs2_(1, 1, vi) + viscs2_(0, 1, vi)*conv_r_(0, 1, ui) + viscs2_(1, 1, vi)*conv_r_(1, 1, ui)) ;
#endif

#ifdef FLUID2_IS_TERM6
    /* Stabilisierung der Viskosität (-L_visc_u) */
    estif_(vi*3, ui*2 + 4)     += 4.0*(nu_*nu_)*ttimetauMp*(viscs2_(0, 0, ui)*viscs2_(0, 0, vi) + viscs2_(0, 1, ui)*viscs2_(0, 1, vi)) ;
    estif_(vi*3, ui*2 + 5)     += 4.0*(nu_*nu_)*ttimetauMp*(viscs2_(0, 0, vi)*viscs2_(0, 1, ui) + viscs2_(0, 1, vi)*viscs2_(1, 1, ui)) ;
    estif_(vi*3 + 1, ui*2 + 4) += 4.0*(nu_*nu_)*ttimetauMp*(viscs2_(0, 0, ui)*viscs2_(0, 1, vi) + viscs2_(0, 1, ui)*viscs2_(1, 1, vi)) ;
    estif_(vi*3 + 1, ui*2 + 5) += 4.0*(nu_*nu_)*ttimetauMp*(viscs2_(0, 1, ui)*viscs2_(0, 1, vi) + viscs2_(1, 1, ui)*viscs2_(1, 1, vi)) ;
#endif

#ifdef FLUID2_IS_TERM7
    /* Druckterm */
#endif

#ifdef FLUID2_IS_TERM8
    /* Divergenzfreiheit */
    estif_(vi*3 + 2, ui*2 + 4) += timefacfac*funct_p_(vi)*derxy_(0, ui) ;
    estif_(vi*3 + 2, ui*2 + 5) += timefacfac*funct_p_(vi)*derxy_(1, ui) ;
#endif

#ifdef FLUID2_IS_TERM9
    /* Kontinuitätsstabilisierung */
    estif_(vi*3, ui*2 + 4)     += (thsl*thsl)*tau_C*derxy_(0, ui)*derxy_(0, vi) ;
    estif_(vi*3, ui*2 + 5)     += (thsl*thsl)*tau_C*derxy_(0, vi)*derxy_(1, ui) ;
    estif_(vi*3 + 1, ui*2 + 4) += (thsl*thsl)*tau_C*derxy_(0, ui)*derxy_(1, vi) ;
    estif_(vi*3 + 1, ui*2 + 5) += (thsl*thsl)*tau_C*derxy_(1, ui)*derxy_(1, vi) ;
#endif

#ifdef FLUID2_IS_TERM10
    /* Massenterm */
    estif_(vi*3, ui*2 + 4)     += fac*funct_(ui)*funct_(vi) ;
    estif_(vi*3 + 1, ui*2 + 5) += fac*funct_(ui)*funct_(vi) ;
#endif

#ifdef FLUID2_IS_TERM11
    /* Konvektionsstabilisierung */
    estif_(vi*3, ui*2 + 4)     += timetauM*funct_(ui)*(conv_c_(vi) + conv_g_(vi) + velint_(0)*derxy_(0, vi)) ;
    estif_(vi*3, ui*2 + 5)     += timetauM*funct_(ui)*velint_(0)*derxy_(1, vi) ;
    estif_(vi*3 + 1, ui*2 + 4) += timetauM*funct_(ui)*velint_(1)*derxy_(0, vi) ;
    estif_(vi*3 + 1, ui*2 + 5) += timetauM*funct_(ui)*(conv_c_(vi) + conv_g_(vi) + velint_(1)*derxy_(1, vi)) ;
#endif

#ifdef FLUID2_IS_TERM12
    /* Viskositätsstabilisierung */
    estif_(vi*3, ui*2 + 4)     += tau_Mp*time2nue*funct_(ui)*viscs2_(0, 0, vi) ;
    estif_(vi*3, ui*2 + 5)     += tau_Mp*time2nue*funct_(ui)*viscs2_(0, 1, vi) ;
    estif_(vi*3 + 1, ui*2 + 4) += tau_Mp*time2nue*funct_(ui)*viscs2_(0, 1, vi) ;
    estif_(vi*3 + 1, ui*2 + 5) += tau_Mp*time2nue*funct_(ui)*viscs2_(1, 1, vi) ;
#endif

#ifdef FLUID2_IS_TERM13
    /* Quellterm der rechten Seite */
#endif

#ifdef FLUID2_IS_TERM14
    /* Konvektionsstabilisierung */
    estif_(vi*3, ui*2 + 4)     += -(timetauM*funct_(ui)*rhsint_(0)*derxy_(0, vi)) ;
    estif_(vi*3, ui*2 + 5)     += -(timetauM*funct_(ui)*rhsint_(0)*derxy_(1, vi)) ;
    estif_(vi*3 + 1, ui*2 + 4) += -(timetauM*funct_(ui)*rhsint_(1)*derxy_(0, vi)) ;
    estif_(vi*3 + 1, ui*2 + 5) += -(timetauM*funct_(ui)*rhsint_(1)*derxy_(1, vi)) ;
#endif

#ifdef FLUID2_IS_TERM15
    /* Viskositätsstabilisierung */
#endif

  }
}

for (vi=4; vi<9; ++vi)
{
  for (ui=0; ui<4; ++ui)
  {
#ifdef FLUID2_IS_TERM1
    /* Konvektionsterm */
    estif_(vi*2 + 4, ui*3)     += timefacfac*funct_(vi)*(conv_c_(ui) + conv_g_(ui) + conv_r_(0, 0, ui)) ;
    estif_(vi*2 + 4, ui*3 + 1) += timefacfac*funct_(vi)*conv_r_(0, 1, ui) ;
    estif_(vi*2 + 5, ui*3)     += timefacfac*funct_(vi)*conv_r_(1, 0, ui) ;
    estif_(vi*2 + 5, ui*3 + 1) += timefacfac*funct_(vi)*(conv_c_(ui) + conv_g_(ui) + conv_r_(1, 1, ui)) ;
#endif

#ifdef FLUID2_IS_TERM2
    /* Stabilisierung der Konvektion ( L_conv_u) */
    estif_(vi*2 + 4, ui*3)     += ttimetauM*(conv_c_(ui)*conv_c_(vi) - conv_g_(ui)*derxy_(0, vi)*gridvint_(0) + 2.0*velint_(0)*derxy_(0, vi)*conv_r_(0, 0, ui) + velint_(1)*derxy_(0, vi)*conv_r_(0, 1, ui) + velint_(1)*derxy_(1, vi)*conv_r_(0, 0, ui) - 2.0*derxy_(0, vi)*gridvint_(0)*conv_r_(0, 0, ui) - derxy_(0, vi)*gridvint_(1)*conv_r_(0, 1, ui) - derxy_(1, vi)*gridvint_(1)*conv_r_(0, 0, ui) + derxy_(1, ui)*derxy_(1, vi)*(gridvint_(1)*gridvint_(1)) - 2.0*velint_(0)*derxy_(0, ui)*derxy_(0, vi)*gridvint_(0) - velint_(0)*derxy_(0, ui)*derxy_(1, vi)*gridvint_(1) - velint_(0)*derxy_(0, vi)*derxy_(1, ui)*gridvint_(1) - velint_(1)*derxy_(0, ui)*derxy_(1, vi)*gridvint_(0) - velint_(1)*derxy_(0, vi)*derxy_(1, ui)*gridvint_(0) - 2.0*velint_(1)*derxy_(1, ui)*derxy_(1, vi)*gridvint_(1) + derxy_(0, ui)*derxy_(1, vi)*gridvint_(0)*gridvint_(1)) ;
    estif_(vi*2 + 4, ui*3 + 1) += ttimetauM*(conv_c_(vi)*conv_r_(0, 1, ui) + velint_(0)*derxy_(1, vi)*conv_r_(0, 0, ui) + velint_(1)*derxy_(1, vi)*conv_r_(0, 1, ui) - derxy_(0, vi)*gridvint_(0)*conv_r_(0, 1, ui) - derxy_(1, vi)*gridvint_(0)*conv_r_(0, 0, ui) - 2.0*derxy_(1, vi)*gridvint_(1)*conv_r_(0, 1, ui)) ;
    estif_(vi*2 + 5, ui*3)     += ttimetauM*(conv_c_(vi)*conv_r_(1, 0, ui) + velint_(0)*derxy_(0, vi)*conv_r_(1, 0, ui) + velint_(1)*derxy_(0, vi)*conv_r_(1, 1, ui) - 2.0*derxy_(0, vi)*gridvint_(0)*conv_r_(1, 0, ui) - derxy_(0, vi)*gridvint_(1)*conv_r_(1, 1, ui) - derxy_(1, vi)*gridvint_(1)*conv_r_(1, 0, ui)) ;
    estif_(vi*2 + 5, ui*3 + 1) += ttimetauM*(conv_c_(ui)*conv_c_(vi) - conv_g_(ui)*derxy_(0, vi)*gridvint_(0) + velint_(0)*derxy_(0, vi)*conv_r_(1, 1, ui) + velint_(0)*derxy_(1, vi)*conv_r_(1, 0, ui) + 2.0*velint_(1)*derxy_(1, vi)*conv_r_(1, 1, ui) - derxy_(0, vi)*gridvint_(0)*conv_r_(1, 1, ui) - derxy_(1, vi)*gridvint_(0)*conv_r_(1, 0, ui) - 2.0*derxy_(1, vi)*gridvint_(1)*conv_r_(1, 1, ui) + derxy_(1, ui)*derxy_(1, vi)*(gridvint_(1)*gridvint_(1)) - 2.0*velint_(0)*derxy_(0, ui)*derxy_(0, vi)*gridvint_(0) - velint_(0)*derxy_(0, ui)*derxy_(1, vi)*gridvint_(1) - velint_(0)*derxy_(0, vi)*derxy_(1, ui)*gridvint_(1) - velint_(1)*derxy_(0, ui)*derxy_(1, vi)*gridvint_(0) - velint_(1)*derxy_(0, vi)*derxy_(1, ui)*gridvint_(0) - 2.0*velint_(1)*derxy_(1, ui)*derxy_(1, vi)*gridvint_(1) + derxy_(0, ui)*derxy_(1, vi)*gridvint_(0)*gridvint_(1)) ;
#endif

#ifdef FLUID2_IS_TERM3
    /* Stabilisierung der Konvektion (-L_visc_u) */
    estif_(vi*2 + 4, ui*3)     += -2.0*nu_*ttimetauM*(-(conv_c_(vi)*viscs2_(0, 0, ui)) - conv_g_(vi)*viscs2_(0, 0, ui) + funct_(ui)*derxy_(0, vi)*visc_old_(0)) ;
    estif_(vi*2 + 4, ui*3 + 1) += -2.0*nu_*ttimetauM*(-(conv_c_(vi)*viscs2_(0, 1, ui)) - conv_g_(vi)*viscs2_(0, 1, ui) + funct_(ui)*derxy_(1, vi)*visc_old_(0)) ;
    estif_(vi*2 + 5, ui*3)     += -2.0*nu_*ttimetauM*(-(conv_c_(vi)*viscs2_(0, 1, ui)) - conv_g_(vi)*viscs2_(0, 1, ui) + funct_(ui)*derxy_(0, vi)*visc_old_(1)) ;
    estif_(vi*2 + 5, ui*3 + 1) += -2.0*nu_*ttimetauM*(-(conv_c_(vi)*viscs2_(1, 1, ui)) - conv_g_(vi)*viscs2_(1, 1, ui) + funct_(ui)*derxy_(1, vi)*visc_old_(1)) ;
#endif

#ifdef FLUID2_IS_TERM4
    /* Viskositätsterm */
    estif_(vi*2 + 4, ui*3)     += fac*time2nue*derxy_(0, ui)*derxy_(0, vi) + nu_*timefacfac*derxy_(1, ui)*derxy_(1, vi) ;
    estif_(vi*2 + 4, ui*3 + 1) += nu_*timefacfac*derxy_(0, ui)*derxy_(1, vi) ;
    estif_(vi*2 + 5, ui*3)     += nu_*timefacfac*derxy_(0, vi)*derxy_(1, ui) ;
    estif_(vi*2 + 5, ui*3 + 1) += fac*time2nue*derxy_(1, ui)*derxy_(1, vi) + nu_*timefacfac*derxy_(0, ui)*derxy_(0, vi) ;
#endif

#ifdef FLUID2_IS_TERM5
    /* Stabilisierung der Viskosität ( L_conv_u) */
    estif_(vi*2 + 4, ui*3)     += 2.0*nu_*ttimetauMp*(conv_c_(ui)*viscs2_(0, 0, vi) + conv_g_(ui)*viscs2_(0, 0, vi) + viscs2_(0, 0, vi)*conv_r_(0, 0, ui) + viscs2_(0, 1, vi)*conv_r_(1, 0, ui)) ;
    estif_(vi*2 + 4, ui*3 + 1) += 2.0*nu_*ttimetauMp*(conv_c_(ui)*viscs2_(0, 1, vi) + conv_g_(ui)*viscs2_(0, 1, vi) + viscs2_(0, 0, vi)*conv_r_(0, 1, ui) + viscs2_(0, 1, vi)*conv_r_(1, 1, ui)) ;
    estif_(vi*2 + 5, ui*3)     += 2.0*nu_*ttimetauMp*(conv_c_(ui)*viscs2_(0, 1, vi) + conv_g_(ui)*viscs2_(0, 1, vi) + viscs2_(0, 1, vi)*conv_r_(0, 0, ui) + viscs2_(1, 1, vi)*conv_r_(1, 0, ui)) ;
    estif_(vi*2 + 5, ui*3 + 1) += 2.0*nu_*ttimetauMp*(conv_c_(ui)*viscs2_(1, 1, vi) + conv_g_(ui)*viscs2_(1, 1, vi) + viscs2_(0, 1, vi)*conv_r_(0, 1, ui) + viscs2_(1, 1, vi)*conv_r_(1, 1, ui)) ;
#endif

#ifdef FLUID2_IS_TERM6
    /* Stabilisierung der Viskosität (-L_visc_u) */
    estif_(vi*2 + 4, ui*3)     += 4.0*(nu_*nu_)*ttimetauMp*(viscs2_(0, 0, ui)*viscs2_(0, 0, vi) + viscs2_(0, 1, ui)*viscs2_(0, 1, vi)) ;
    estif_(vi*2 + 4, ui*3 + 1) += 4.0*(nu_*nu_)*ttimetauMp*(viscs2_(0, 0, vi)*viscs2_(0, 1, ui) + viscs2_(0, 1, vi)*viscs2_(1, 1, ui)) ;
    estif_(vi*2 + 5, ui*3)     += 4.0*(nu_*nu_)*ttimetauMp*(viscs2_(0, 0, ui)*viscs2_(0, 1, vi) + viscs2_(0, 1, ui)*viscs2_(1, 1, vi)) ;
    estif_(vi*2 + 5, ui*3 + 1) += 4.0*(nu_*nu_)*ttimetauMp*(viscs2_(0, 1, ui)*viscs2_(0, 1, vi) + viscs2_(1, 1, ui)*viscs2_(1, 1, vi)) ;
#endif

#ifdef FLUID2_IS_TERM7
    /* Druckterm */
    estif_(vi*2 + 4, ui*3 + 2) += -(timefacfac*funct_p_(ui)*derxy_(0, vi)) ;
    estif_(vi*2 + 5, ui*3 + 2) += -(timefacfac*funct_p_(ui)*derxy_(1, vi)) ;
#endif

#ifdef FLUID2_IS_TERM8
    /* Divergenzfreiheit */
#endif

#ifdef FLUID2_IS_TERM9
    /* Kontinuitätsstabilisierung */
    estif_(vi*2 + 4, ui*3)     += (thsl*thsl)*tau_C*derxy_(0, ui)*derxy_(0, vi) ;
    estif_(vi*2 + 4, ui*3 + 1) += (thsl*thsl)*tau_C*derxy_(0, vi)*derxy_(1, ui) ;
    estif_(vi*2 + 5, ui*3)     += (thsl*thsl)*tau_C*derxy_(0, ui)*derxy_(1, vi) ;
    estif_(vi*2 + 5, ui*3 + 1) += (thsl*thsl)*tau_C*derxy_(1, ui)*derxy_(1, vi) ;
#endif

#ifdef FLUID2_IS_TERM10
    /* Massenterm */
    estif_(vi*2 + 4, ui*3)     += fac*funct_(ui)*funct_(vi) ;
    estif_(vi*2 + 5, ui*3 + 1) += fac*funct_(ui)*funct_(vi) ;
#endif

#ifdef FLUID2_IS_TERM11
    /* Konvektionsstabilisierung */
    estif_(vi*2 + 4, ui*3)     += timetauM*funct_(ui)*(conv_c_(vi) + conv_g_(vi) + velint_(0)*derxy_(0, vi)) ;
    estif_(vi*2 + 4, ui*3 + 1) += timetauM*funct_(ui)*velint_(0)*derxy_(1, vi) ;
    estif_(vi*2 + 5, ui*3)     += timetauM*funct_(ui)*velint_(1)*derxy_(0, vi) ;
    estif_(vi*2 + 5, ui*3 + 1) += timetauM*funct_(ui)*(conv_c_(vi) + conv_g_(vi) + velint_(1)*derxy_(1, vi)) ;
#endif

#ifdef FLUID2_IS_TERM12
    /* Viskositätsstabilisierung */
    estif_(vi*2 + 4, ui*3)     += tau_Mp*time2nue*funct_(ui)*viscs2_(0, 0, vi) ;
    estif_(vi*2 + 4, ui*3 + 1) += tau_Mp*time2nue*funct_(ui)*viscs2_(0, 1, vi) ;
    estif_(vi*2 + 5, ui*3)     += tau_Mp*time2nue*funct_(ui)*viscs2_(0, 1, vi) ;
    estif_(vi*2 + 5, ui*3 + 1) += tau_Mp*time2nue*funct_(ui)*viscs2_(1, 1, vi) ;
#endif

#ifdef FLUID2_IS_TERM13
    /* Quellterm der rechten Seite */
#endif

#ifdef FLUID2_IS_TERM14
    /* Konvektionsstabilisierung */
    estif_(vi*2 + 4, ui*3)     += -(timetauM*funct_(ui)*rhsint_(0)*derxy_(0, vi)) ;
    estif_(vi*2 + 4, ui*3 + 1) += -(timetauM*funct_(ui)*rhsint_(0)*derxy_(1, vi)) ;
    estif_(vi*2 + 5, ui*3)     += -(timetauM*funct_(ui)*rhsint_(1)*derxy_(0, vi)) ;
    estif_(vi*2 + 5, ui*3 + 1) += -(timetauM*funct_(ui)*rhsint_(1)*derxy_(1, vi)) ;
#endif

#ifdef FLUID2_IS_TERM15
    /* Viskositätsstabilisierung */
#endif

  }
}

for (vi=4; vi<9; ++vi)
{
  for (ui=4; ui<9; ++ui)
  {
#ifdef FLUID2_IS_TERM1
    /* Konvektionsterm */
    estif_(vi*2 + 4, ui*2 + 4) += timefacfac*funct_(vi)*(conv_c_(ui) + conv_g_(ui) + conv_r_(0, 0, ui)) ;
    estif_(vi*2 + 4, ui*2 + 5) += timefacfac*funct_(vi)*conv_r_(0, 1, ui) ;
    estif_(vi*2 + 5, ui*2 + 4) += timefacfac*funct_(vi)*conv_r_(1, 0, ui) ;
    estif_(vi*2 + 5, ui*2 + 5) += timefacfac*funct_(vi)*(conv_c_(ui) + conv_g_(ui) + conv_r_(1, 1, ui)) ;
#endif

#ifdef FLUID2_IS_TERM2
    /* Stabilisierung der Konvektion ( L_conv_u) */
    estif_(vi*2 + 4, ui*2 + 4) += ttimetauM*(conv_c_(ui)*conv_c_(vi) - conv_g_(ui)*derxy_(0, vi)*gridvint_(0) + 2.0*velint_(0)*derxy_(0, vi)*conv_r_(0, 0, ui) + velint_(1)*derxy_(0, vi)*conv_r_(0, 1, ui) + velint_(1)*derxy_(1, vi)*conv_r_(0, 0, ui) - 2.0*derxy_(0, vi)*gridvint_(0)*conv_r_(0, 0, ui) - derxy_(0, vi)*gridvint_(1)*conv_r_(0, 1, ui) - derxy_(1, vi)*gridvint_(1)*conv_r_(0, 0, ui) + derxy_(1, ui)*derxy_(1, vi)*(gridvint_(1)*gridvint_(1)) - 2.0*velint_(0)*derxy_(0, ui)*derxy_(0, vi)*gridvint_(0) - velint_(0)*derxy_(0, ui)*derxy_(1, vi)*gridvint_(1) - velint_(0)*derxy_(0, vi)*derxy_(1, ui)*gridvint_(1) - velint_(1)*derxy_(0, ui)*derxy_(1, vi)*gridvint_(0) - velint_(1)*derxy_(0, vi)*derxy_(1, ui)*gridvint_(0) - 2.0*velint_(1)*derxy_(1, ui)*derxy_(1, vi)*gridvint_(1) + derxy_(0, ui)*derxy_(1, vi)*gridvint_(0)*gridvint_(1)) ;
    estif_(vi*2 + 4, ui*2 + 5) += ttimetauM*(conv_c_(vi)*conv_r_(0, 1, ui) + velint_(0)*derxy_(1, vi)*conv_r_(0, 0, ui) + velint_(1)*derxy_(1, vi)*conv_r_(0, 1, ui) - derxy_(0, vi)*gridvint_(0)*conv_r_(0, 1, ui) - derxy_(1, vi)*gridvint_(0)*conv_r_(0, 0, ui) - 2.0*derxy_(1, vi)*gridvint_(1)*conv_r_(0, 1, ui)) ;
    estif_(vi*2 + 5, ui*2 + 4) += ttimetauM*(conv_c_(vi)*conv_r_(1, 0, ui) + velint_(0)*derxy_(0, vi)*conv_r_(1, 0, ui) + velint_(1)*derxy_(0, vi)*conv_r_(1, 1, ui) - 2.0*derxy_(0, vi)*gridvint_(0)*conv_r_(1, 0, ui) - derxy_(0, vi)*gridvint_(1)*conv_r_(1, 1, ui) - derxy_(1, vi)*gridvint_(1)*conv_r_(1, 0, ui)) ;
    estif_(vi*2 + 5, ui*2 + 5) += ttimetauM*(conv_c_(ui)*conv_c_(vi) - conv_g_(ui)*derxy_(0, vi)*gridvint_(0) + velint_(0)*derxy_(0, vi)*conv_r_(1, 1, ui) + velint_(0)*derxy_(1, vi)*conv_r_(1, 0, ui) + 2.0*velint_(1)*derxy_(1, vi)*conv_r_(1, 1, ui) - derxy_(0, vi)*gridvint_(0)*conv_r_(1, 1, ui) - derxy_(1, vi)*gridvint_(0)*conv_r_(1, 0, ui) - 2.0*derxy_(1, vi)*gridvint_(1)*conv_r_(1, 1, ui) + derxy_(1, ui)*derxy_(1, vi)*(gridvint_(1)*gridvint_(1)) - 2.0*velint_(0)*derxy_(0, ui)*derxy_(0, vi)*gridvint_(0) - velint_(0)*derxy_(0, ui)*derxy_(1, vi)*gridvint_(1) - velint_(0)*derxy_(0, vi)*derxy_(1, ui)*gridvint_(1) - velint_(1)*derxy_(0, ui)*derxy_(1, vi)*gridvint_(0) - velint_(1)*derxy_(0, vi)*derxy_(1, ui)*gridvint_(0) - 2.0*velint_(1)*derxy_(1, ui)*derxy_(1, vi)*gridvint_(1) + derxy_(0, ui)*derxy_(1, vi)*gridvint_(0)*gridvint_(1)) ;
#endif

#ifdef FLUID2_IS_TERM3
    /* Stabilisierung der Konvektion (-L_visc_u) */
    estif_(vi*2 + 4, ui*2 + 4) += -2.0*nu_*ttimetauM*(-(conv_c_(vi)*viscs2_(0, 0, ui)) - conv_g_(vi)*viscs2_(0, 0, ui) + funct_(ui)*derxy_(0, vi)*visc_old_(0)) ;
    estif_(vi*2 + 4, ui*2 + 5) += -2.0*nu_*ttimetauM*(-(conv_c_(vi)*viscs2_(0, 1, ui)) - conv_g_(vi)*viscs2_(0, 1, ui) + funct_(ui)*derxy_(1, vi)*visc_old_(0)) ;
    estif_(vi*2 + 5, ui*2 + 4) += -2.0*nu_*ttimetauM*(-(conv_c_(vi)*viscs2_(0, 1, ui)) - conv_g_(vi)*viscs2_(0, 1, ui) + funct_(ui)*derxy_(0, vi)*visc_old_(1)) ;
    estif_(vi*2 + 5, ui*2 + 5) += -2.0*nu_*ttimetauM*(-(conv_c_(vi)*viscs2_(1, 1, ui)) - conv_g_(vi)*viscs2_(1, 1, ui) + funct_(ui)*derxy_(1, vi)*visc_old_(1)) ;
#endif

#ifdef FLUID2_IS_TERM4
    /* Viskositätsterm */
    estif_(vi*2 + 4, ui*2 + 4) += fac*time2nue*derxy_(0, ui)*derxy_(0, vi) + nu_*timefacfac*derxy_(1, ui)*derxy_(1, vi) ;
    estif_(vi*2 + 4, ui*2 + 5) += nu_*timefacfac*derxy_(0, ui)*derxy_(1, vi) ;
    estif_(vi*2 + 5, ui*2 + 4) += nu_*timefacfac*derxy_(0, vi)*derxy_(1, ui) ;
    estif_(vi*2 + 5, ui*2 + 5) += fac*time2nue*derxy_(1, ui)*derxy_(1, vi) + nu_*timefacfac*derxy_(0, ui)*derxy_(0, vi) ;
#endif

#ifdef FLUID2_IS_TERM5
    /* Stabilisierung der Viskosität ( L_conv_u) */
    estif_(vi*2 + 4, ui*2 + 4) += 2.0*nu_*ttimetauMp*(conv_c_(ui)*viscs2_(0, 0, vi) + conv_g_(ui)*viscs2_(0, 0, vi) + viscs2_(0, 0, vi)*conv_r_(0, 0, ui) + viscs2_(0, 1, vi)*conv_r_(1, 0, ui)) ;
    estif_(vi*2 + 4, ui*2 + 5) += 2.0*nu_*ttimetauMp*(conv_c_(ui)*viscs2_(0, 1, vi) + conv_g_(ui)*viscs2_(0, 1, vi) + viscs2_(0, 0, vi)*conv_r_(0, 1, ui) + viscs2_(0, 1, vi)*conv_r_(1, 1, ui)) ;
    estif_(vi*2 + 5, ui*2 + 4) += 2.0*nu_*ttimetauMp*(conv_c_(ui)*viscs2_(0, 1, vi) + conv_g_(ui)*viscs2_(0, 1, vi) + viscs2_(0, 1, vi)*conv_r_(0, 0, ui) + viscs2_(1, 1, vi)*conv_r_(1, 0, ui)) ;
    estif_(vi*2 + 5, ui*2 + 5) += 2.0*nu_*ttimetauMp*(conv_c_(ui)*viscs2_(1, 1, vi) + conv_g_(ui)*viscs2_(1, 1, vi) + viscs2_(0, 1, vi)*conv_r_(0, 1, ui) + viscs2_(1, 1, vi)*conv_r_(1, 1, ui)) ;
#endif

#ifdef FLUID2_IS_TERM6
    /* Stabilisierung der Viskosität (-L_visc_u) */
    estif_(vi*2 + 4, ui*2 + 4) += 4.0*(nu_*nu_)*ttimetauMp*(viscs2_(0, 0, ui)*viscs2_(0, 0, vi) + viscs2_(0, 1, ui)*viscs2_(0, 1, vi)) ;
    estif_(vi*2 + 4, ui*2 + 5) += 4.0*(nu_*nu_)*ttimetauMp*(viscs2_(0, 0, vi)*viscs2_(0, 1, ui) + viscs2_(0, 1, vi)*viscs2_(1, 1, ui)) ;
    estif_(vi*2 + 5, ui*2 + 4) += 4.0*(nu_*nu_)*ttimetauMp*(viscs2_(0, 0, ui)*viscs2_(0, 1, vi) + viscs2_(0, 1, ui)*viscs2_(1, 1, vi)) ;
    estif_(vi*2 + 5, ui*2 + 5) += 4.0*(nu_*nu_)*ttimetauMp*(viscs2_(0, 1, ui)*viscs2_(0, 1, vi) + viscs2_(1, 1, ui)*viscs2_(1, 1, vi)) ;
#endif

#ifdef FLUID2_IS_TERM7
    /* Druckterm */
#endif

#ifdef FLUID2_IS_TERM8
    /* Divergenzfreiheit */
#endif

#ifdef FLUID2_IS_TERM9
    /* Kontinuitätsstabilisierung */
    estif_(vi*2 + 4, ui*2 + 4) += (thsl*thsl)*tau_C*derxy_(0, ui)*derxy_(0, vi) ;
    estif_(vi*2 + 4, ui*2 + 5) += (thsl*thsl)*tau_C*derxy_(0, vi)*derxy_(1, ui) ;
    estif_(vi*2 + 5, ui*2 + 4) += (thsl*thsl)*tau_C*derxy_(0, ui)*derxy_(1, vi) ;
    estif_(vi*2 + 5, ui*2 + 5) += (thsl*thsl)*tau_C*derxy_(1, ui)*derxy_(1, vi) ;
#endif

#ifdef FLUID2_IS_TERM10
    /* Massenterm */
    estif_(vi*2 + 4, ui*2 + 4) += fac*funct_(ui)*funct_(vi) ;
    estif_(vi*2 + 5, ui*2 + 5) += fac*funct_(ui)*funct_(vi) ;
#endif

#ifdef FLUID2_IS_TERM11
    /* Konvektionsstabilisierung */
    estif_(vi*2 + 4, ui*2 + 4) += timetauM*funct_(ui)*(conv_c_(vi) + conv_g_(vi) + velint_(0)*derxy_(0, vi)) ;
    estif_(vi*2 + 4, ui*2 + 5) += timetauM*funct_(ui)*velint_(0)*derxy_(1, vi) ;
    estif_(vi*2 + 5, ui*2 + 4) += timetauM*funct_(ui)*velint_(1)*derxy_(0, vi) ;
    estif_(vi*2 + 5, ui*2 + 5) += timetauM*funct_(ui)*(conv_c_(vi) + conv_g_(vi) + velint_(1)*derxy_(1, vi)) ;
#endif

#ifdef FLUID2_IS_TERM12
    /* Viskositätsstabilisierung */
    estif_(vi*2 + 4, ui*2 + 4) += tau_Mp*time2nue*funct_(ui)*viscs2_(0, 0, vi) ;
    estif_(vi*2 + 4, ui*2 + 5) += tau_Mp*time2nue*funct_(ui)*viscs2_(0, 1, vi) ;
    estif_(vi*2 + 5, ui*2 + 4) += tau_Mp*time2nue*funct_(ui)*viscs2_(0, 1, vi) ;
    estif_(vi*2 + 5, ui*2 + 5) += tau_Mp*time2nue*funct_(ui)*viscs2_(1, 1, vi) ;
#endif

#ifdef FLUID2_IS_TERM13
    /* Quellterm der rechten Seite */
#endif

#ifdef FLUID2_IS_TERM14
    /* Konvektionsstabilisierung */
    estif_(vi*2 + 4, ui*2 + 4) += -(timetauM*funct_(ui)*rhsint_(0)*derxy_(0, vi)) ;
    estif_(vi*2 + 4, ui*2 + 5) += -(timetauM*funct_(ui)*rhsint_(0)*derxy_(1, vi)) ;
    estif_(vi*2 + 5, ui*2 + 4) += -(timetauM*funct_(ui)*rhsint_(1)*derxy_(0, vi)) ;
    estif_(vi*2 + 5, ui*2 + 5) += -(timetauM*funct_(ui)*rhsint_(1)*derxy_(1, vi)) ;
#endif

#ifdef FLUID2_IS_TERM15
    /* Viskositätsstabilisierung */
#endif

  }
}
