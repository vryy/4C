/*
<pre>
Maintainer: Ulrich Küttler
            kuettler@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/Members/kuettler
            089 - 289-15238
</pre>
*/

for (vi=0; vi<8; ++vi)
{
  for (ui=0; ui<8; ++ui)
  {
#ifdef FLUID3_IS_TERM1
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
#endif

#ifdef FLUID3_IS_TERM2
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
#endif

#ifdef FLUID3_IS_TERM3
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
#endif

#ifdef FLUID3_IS_TERM4
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
#endif

#ifdef FLUID3_IS_TERM5
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
#endif

#ifdef FLUID3_IS_TERM6
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
#endif

#ifdef FLUID3_IS_TERM7
    /* Druckterm */
    estif_(vi*4, ui*4 + 3)     += -(timefacfac*funct_p_(ui)*derxyz_(0, vi)) ;
    estif_(vi*4 + 1, ui*4 + 3) += -(timefacfac*funct_p_(ui)*derxyz_(1, vi)) ;
    estif_(vi*4 + 2, ui*4 + 3) += -(timefacfac*funct_p_(ui)*derxyz_(2, vi)) ;
#endif

#ifdef FLUID3_IS_TERM8
    /* Divergenzfreiheit */
    estif_(vi*4 + 3, ui*4)     += timefacfac*funct_p_(vi)*derxyz_(0, ui) ;
    estif_(vi*4 + 3, ui*4 + 1) += timefacfac*funct_p_(vi)*derxyz_(1, ui) ;
    estif_(vi*4 + 3, ui*4 + 2) += timefacfac*funct_p_(vi)*derxyz_(2, ui) ;
#endif

#ifdef FLUID3_IS_TERM9
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
#endif

#ifdef FLUID3_IS_TERM10
    /* Massenterm */
    estif_(vi*4, ui*4)         += fac*funct_(ui)*funct_(vi) ;
    estif_(vi*4 + 1, ui*4 + 1) += fac*funct_(ui)*funct_(vi) ;
    estif_(vi*4 + 2, ui*4 + 2) += fac*funct_(ui)*funct_(vi) ;
#endif

#ifdef FLUID3_IS_TERM11
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
#endif

#ifdef FLUID3_IS_TERM12
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
#endif

#ifdef FLUID3_IS_TERM13
    /* Quellterm der rechten Seite */
#endif

#ifdef FLUID3_IS_TERM14
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
#endif

#ifdef FLUID3_IS_TERM15
    /* Viskositätsstabilisierung */
#endif

  }
}

for (vi=0; vi<8; ++vi)
{
  for (ui=8; ui<iel; ++ui)
  {
#ifdef FLUID3_IS_TERM1
    /* Konvektionsterm */
    estif_(vi*4, ui*3 + 8)     += timefacfac*funct_(vi)*(conv_c_(ui) + conv_g_(ui) + conv_r_(0, 0, ui)) ;
    estif_(vi*4, ui*3 + 9)     += timefacfac*funct_(vi)*conv_r_(0, 1, ui) ;
    estif_(vi*4, ui*3 + 10)     += timefacfac*funct_(vi)*conv_r_(0, 2, ui) ;
    estif_(vi*4 + 1, ui*3 + 8) += timefacfac*funct_(vi)*conv_r_(1, 0, ui) ;
    estif_(vi*4 + 1, ui*3 + 9) += timefacfac*funct_(vi)*(conv_c_(ui) + conv_g_(ui) + conv_r_(1, 1, ui)) ;
    estif_(vi*4 + 1, ui*3 + 10) += timefacfac*funct_(vi)*conv_r_(1, 2, ui) ;
    estif_(vi*4 + 2, ui*3 + 8) += timefacfac*funct_(vi)*conv_r_(2, 0, ui) ;
    estif_(vi*4 + 2, ui*3 + 9) += timefacfac*funct_(vi)*conv_r_(2, 1, ui) ;
    estif_(vi*4 + 2, ui*3 + 10) += timefacfac*funct_(vi)*(conv_c_(ui) + conv_g_(ui) + conv_r_(2, 2, ui)) ;
#endif

#ifdef FLUID3_IS_TERM2
    /* Stabilisierung der Konvektion ( L_conv_u) */
    estif_(vi*4, ui*3 + 8)     += ttimetauM*(conv_c_(ui)*conv_c_(vi) + conv_c_(ui)*conv_g_(vi) + conv_c_(vi)*conv_g_(ui) + conv_g_(ui)*conv_g_(vi) + conv_c_(vi)*conv_r_(0, 0, ui) + conv_g_(vi)*conv_r_(0, 0, ui) + velint_(0)*derxyz_(0, vi)*conv_r_(0, 0, ui) + velint_(1)*derxyz_(0, vi)*conv_r_(0, 1, ui) + velint_(2)*derxyz_(0, vi)*conv_r_(0, 2, ui) - gridvint_(0)*derxyz_(0, vi)*conv_r_(0, 0, ui) - gridvint_(1)*derxyz_(0, vi)*conv_r_(0, 1, ui) - gridvint_(2)*derxyz_(0, vi)*conv_r_(0, 2, ui)) ;
    estif_(vi*4, ui*3 + 9)     += ttimetauM*(conv_c_(vi)*conv_r_(0, 1, ui) + conv_g_(vi)*conv_r_(0, 1, ui) + velint_(0)*derxyz_(1, vi)*conv_r_(0, 0, ui) + velint_(1)*derxyz_(1, vi)*conv_r_(0, 1, ui) + velint_(2)*derxyz_(1, vi)*conv_r_(0, 2, ui) - gridvint_(0)*derxyz_(1, vi)*conv_r_(0, 0, ui) - gridvint_(1)*derxyz_(1, vi)*conv_r_(0, 1, ui) - gridvint_(2)*derxyz_(1, vi)*conv_r_(0, 2, ui)) ;
    estif_(vi*4, ui*3 + 10)     += ttimetauM*(conv_c_(vi)*conv_r_(0, 2, ui) + conv_g_(vi)*conv_r_(0, 2, ui) + velint_(0)*derxyz_(2, vi)*conv_r_(0, 0, ui) + velint_(1)*derxyz_(2, vi)*conv_r_(0, 1, ui) + velint_(2)*derxyz_(2, vi)*conv_r_(0, 2, ui) - gridvint_(0)*derxyz_(2, vi)*conv_r_(0, 0, ui) - gridvint_(1)*derxyz_(2, vi)*conv_r_(0, 1, ui) - gridvint_(2)*derxyz_(2, vi)*conv_r_(0, 2, ui)) ;
    estif_(vi*4 + 1, ui*3 + 8) += ttimetauM*(conv_c_(vi)*conv_r_(1, 0, ui) + conv_g_(vi)*conv_r_(1, 0, ui) + velint_(0)*derxyz_(0, vi)*conv_r_(1, 0, ui) + velint_(1)*derxyz_(0, vi)*conv_r_(1, 1, ui) + velint_(2)*derxyz_(0, vi)*conv_r_(1, 2, ui) - gridvint_(0)*derxyz_(0, vi)*conv_r_(1, 0, ui) - gridvint_(1)*derxyz_(0, vi)*conv_r_(1, 1, ui) - gridvint_(2)*derxyz_(0, vi)*conv_r_(1, 2, ui)) ;
    estif_(vi*4 + 1, ui*3 + 9) += ttimetauM*(conv_c_(ui)*conv_c_(vi) + conv_c_(ui)*conv_g_(vi) + conv_c_(vi)*conv_g_(ui) + conv_g_(ui)*conv_g_(vi) + conv_c_(vi)*conv_r_(1, 1, ui) + conv_g_(vi)*conv_r_(1, 1, ui) + velint_(0)*derxyz_(1, vi)*conv_r_(1, 0, ui) + velint_(1)*derxyz_(1, vi)*conv_r_(1, 1, ui) + velint_(2)*derxyz_(1, vi)*conv_r_(1, 2, ui) - gridvint_(0)*derxyz_(1, vi)*conv_r_(1, 0, ui) - gridvint_(1)*derxyz_(1, vi)*conv_r_(1, 1, ui) - gridvint_(2)*derxyz_(1, vi)*conv_r_(1, 2, ui)) ;
    estif_(vi*4 + 1, ui*3 + 10) += ttimetauM*(conv_c_(vi)*conv_r_(1, 2, ui) + conv_g_(vi)*conv_r_(1, 2, ui) + velint_(0)*derxyz_(2, vi)*conv_r_(1, 0, ui) + velint_(1)*derxyz_(2, vi)*conv_r_(1, 1, ui) + velint_(2)*derxyz_(2, vi)*conv_r_(1, 2, ui) - gridvint_(0)*derxyz_(2, vi)*conv_r_(1, 0, ui) - gridvint_(1)*derxyz_(2, vi)*conv_r_(1, 1, ui) - gridvint_(2)*derxyz_(2, vi)*conv_r_(1, 2, ui)) ;
    estif_(vi*4 + 2, ui*3 + 8) += ttimetauM*(conv_c_(vi)*conv_r_(2, 0, ui) + conv_g_(vi)*conv_r_(2, 0, ui) + velint_(0)*derxyz_(0, vi)*conv_r_(2, 0, ui) + velint_(1)*derxyz_(0, vi)*conv_r_(2, 1, ui) + velint_(2)*derxyz_(0, vi)*conv_r_(2, 2, ui) - gridvint_(0)*derxyz_(0, vi)*conv_r_(2, 0, ui) - gridvint_(1)*derxyz_(0, vi)*conv_r_(2, 1, ui) - gridvint_(2)*derxyz_(0, vi)*conv_r_(2, 2, ui)) ;
    estif_(vi*4 + 2, ui*3 + 9) += ttimetauM*(conv_c_(vi)*conv_r_(2, 1, ui) + conv_g_(vi)*conv_r_(2, 1, ui) + velint_(0)*derxyz_(1, vi)*conv_r_(2, 0, ui) + velint_(1)*derxyz_(1, vi)*conv_r_(2, 1, ui) + velint_(2)*derxyz_(1, vi)*conv_r_(2, 2, ui) - gridvint_(0)*derxyz_(1, vi)*conv_r_(2, 0, ui) - gridvint_(1)*derxyz_(1, vi)*conv_r_(2, 1, ui) - gridvint_(2)*derxyz_(1, vi)*conv_r_(2, 2, ui)) ;
    estif_(vi*4 + 2, ui*3 + 10) += ttimetauM*(conv_c_(ui)*conv_c_(vi) + conv_c_(ui)*conv_g_(vi) + conv_c_(vi)*conv_g_(ui) + conv_g_(ui)*conv_g_(vi) + conv_c_(vi)*conv_r_(2, 2, ui) + conv_g_(vi)*conv_r_(2, 2, ui) + velint_(0)*derxyz_(2, vi)*conv_r_(2, 0, ui) + velint_(1)*derxyz_(2, vi)*conv_r_(2, 1, ui) + velint_(2)*derxyz_(2, vi)*conv_r_(2, 2, ui) - gridvint_(0)*derxyz_(2, vi)*conv_r_(2, 0, ui) - gridvint_(1)*derxyz_(2, vi)*conv_r_(2, 1, ui) - gridvint_(2)*derxyz_(2, vi)*conv_r_(2, 2, ui)) ;
#endif

#ifdef FLUID3_IS_TERM3
    /* Stabilisierung der Konvektion (-L_visc_u) */
    estif_(vi*4, ui*3 + 8)     += -2.0*nu_*ttimetauM*(-(conv_c_(vi)*viscs2_(0, 0, ui)) - conv_g_(vi)*viscs2_(0, 0, ui) + funct_(ui)*visc_old_(0)*derxyz_(0, vi)) ;
    estif_(vi*4, ui*3 + 9)     += -2.0*nu_*ttimetauM*(-(conv_c_(vi)*viscs2_(0, 1, ui)) - conv_g_(vi)*viscs2_(0, 1, ui) + funct_(ui)*visc_old_(0)*derxyz_(1, vi)) ;
    estif_(vi*4, ui*3 + 10)     += -2.0*nu_*ttimetauM*(-(conv_c_(vi)*viscs2_(0, 2, ui)) - conv_g_(vi)*viscs2_(0, 2, ui) + funct_(ui)*visc_old_(0)*derxyz_(2, vi)) ;
    estif_(vi*4 + 1, ui*3 + 8) += -2.0*nu_*ttimetauM*(-(conv_c_(vi)*viscs2_(0, 1, ui)) - conv_g_(vi)*viscs2_(0, 1, ui) + funct_(ui)*visc_old_(1)*derxyz_(0, vi)) ;
    estif_(vi*4 + 1, ui*3 + 9) += -2.0*nu_*ttimetauM*(-(conv_c_(vi)*viscs2_(1, 1, ui)) - conv_g_(vi)*viscs2_(1, 1, ui) + funct_(ui)*visc_old_(1)*derxyz_(1, vi)) ;
    estif_(vi*4 + 1, ui*3 + 10) += -2.0*nu_*ttimetauM*(-(conv_c_(vi)*viscs2_(1, 2, ui)) - conv_g_(vi)*viscs2_(1, 2, ui) + funct_(ui)*visc_old_(1)*derxyz_(2, vi)) ;
    estif_(vi*4 + 2, ui*3 + 8) += -2.0*nu_*ttimetauM*(-(conv_c_(vi)*viscs2_(0, 2, ui)) - conv_g_(vi)*viscs2_(0, 2, ui) + funct_(ui)*visc_old_(2)*derxyz_(0, vi)) ;
    estif_(vi*4 + 2, ui*3 + 9) += -2.0*nu_*ttimetauM*(-(conv_c_(vi)*viscs2_(1, 2, ui)) - conv_g_(vi)*viscs2_(1, 2, ui) + funct_(ui)*visc_old_(2)*derxyz_(1, vi)) ;
    estif_(vi*4 + 2, ui*3 + 10) += -2.0*nu_*ttimetauM*(-(conv_c_(vi)*viscs2_(2, 2, ui)) - conv_g_(vi)*viscs2_(2, 2, ui) + funct_(ui)*visc_old_(2)*derxyz_(2, vi)) ;
#endif

#ifdef FLUID3_IS_TERM4
    /* Viskositätsterm */
    estif_(vi*4, ui*3 + 8)     += nu_*timefacfac*(2.0*derxyz_(0, ui)*derxyz_(0, vi) + derxyz_(1, ui)*derxyz_(1, vi) + derxyz_(2, ui)*derxyz_(2, vi)) ;
    estif_(vi*4, ui*3 + 9)     += nu_*timefacfac*derxyz_(0, ui)*derxyz_(1, vi) ;
    estif_(vi*4, ui*3 + 10)     += nu_*timefacfac*derxyz_(0, ui)*derxyz_(2, vi) ;
    estif_(vi*4 + 1, ui*3 + 8) += nu_*timefacfac*derxyz_(0, vi)*derxyz_(1, ui) ;
    estif_(vi*4 + 1, ui*3 + 9) += nu_*timefacfac*(derxyz_(0, ui)*derxyz_(0, vi) + 2.0*derxyz_(1, ui)*derxyz_(1, vi) + derxyz_(2, ui)*derxyz_(2, vi)) ;
    estif_(vi*4 + 1, ui*3 + 10) += nu_*timefacfac*derxyz_(1, ui)*derxyz_(2, vi) ;
    estif_(vi*4 + 2, ui*3 + 8) += nu_*timefacfac*derxyz_(0, vi)*derxyz_(2, ui) ;
    estif_(vi*4 + 2, ui*3 + 9) += nu_*timefacfac*derxyz_(1, vi)*derxyz_(2, ui) ;
    estif_(vi*4 + 2, ui*3 + 10) += nu_*timefacfac*(derxyz_(0, ui)*derxyz_(0, vi) + derxyz_(1, ui)*derxyz_(1, vi) + 2.0*derxyz_(2, ui)*derxyz_(2, vi)) ;
#endif

#ifdef FLUID3_IS_TERM5
    /* Stabilisierung der Viskosität ( L_conv_u) */
    estif_(vi*4, ui*3 + 8)     += 2.0*nu_*ttimetauMp*(conv_c_(ui)*viscs2_(0, 0, vi) + conv_g_(ui)*viscs2_(0, 0, vi) + viscs2_(0, 0, vi)*conv_r_(0, 0, ui) + viscs2_(0, 1, vi)*conv_r_(1, 0, ui) + viscs2_(0, 2, vi)*conv_r_(2, 0, ui)) ;
    estif_(vi*4, ui*3 + 9)     += 2.0*nu_*ttimetauMp*(conv_c_(ui)*viscs2_(0, 1, vi) + conv_g_(ui)*viscs2_(0, 1, vi) + viscs2_(0, 0, vi)*conv_r_(0, 1, ui) + viscs2_(0, 1, vi)*conv_r_(1, 1, ui) + viscs2_(0, 2, vi)*conv_r_(2, 1, ui)) ;
    estif_(vi*4, ui*3 + 10)     += 2.0*nu_*ttimetauMp*(conv_c_(ui)*viscs2_(0, 2, vi) + conv_g_(ui)*viscs2_(0, 2, vi) + viscs2_(0, 0, vi)*conv_r_(0, 2, ui) + viscs2_(0, 1, vi)*conv_r_(1, 2, ui) + viscs2_(0, 2, vi)*conv_r_(2, 2, ui)) ;
    estif_(vi*4 + 1, ui*3 + 8) += 2.0*nu_*ttimetauMp*(conv_c_(ui)*viscs2_(0, 1, vi) + conv_g_(ui)*viscs2_(0, 1, vi) + viscs2_(0, 1, vi)*conv_r_(0, 0, ui) + viscs2_(1, 1, vi)*conv_r_(1, 0, ui) + viscs2_(1, 2, vi)*conv_r_(2, 0, ui)) ;
    estif_(vi*4 + 1, ui*3 + 9) += 2.0*nu_*ttimetauMp*(conv_c_(ui)*viscs2_(1, 1, vi) + conv_g_(ui)*viscs2_(1, 1, vi) + viscs2_(0, 1, vi)*conv_r_(0, 1, ui) + viscs2_(1, 1, vi)*conv_r_(1, 1, ui) + viscs2_(1, 2, vi)*conv_r_(2, 1, ui)) ;
    estif_(vi*4 + 1, ui*3 + 10) += 2.0*nu_*ttimetauMp*(conv_c_(ui)*viscs2_(1, 2, vi) + conv_g_(ui)*viscs2_(1, 2, vi) + viscs2_(0, 1, vi)*conv_r_(0, 2, ui) + viscs2_(1, 1, vi)*conv_r_(1, 2, ui) + viscs2_(1, 2, vi)*conv_r_(2, 2, ui)) ;
    estif_(vi*4 + 2, ui*3 + 8) += 2.0*nu_*ttimetauMp*(conv_c_(ui)*viscs2_(0, 2, vi) + conv_g_(ui)*viscs2_(0, 2, vi) + viscs2_(0, 2, vi)*conv_r_(0, 0, ui) + viscs2_(1, 2, vi)*conv_r_(1, 0, ui) + viscs2_(2, 2, vi)*conv_r_(2, 0, ui)) ;
    estif_(vi*4 + 2, ui*3 + 9) += 2.0*nu_*ttimetauMp*(conv_c_(ui)*viscs2_(1, 2, vi) + conv_g_(ui)*viscs2_(1, 2, vi) + viscs2_(0, 2, vi)*conv_r_(0, 1, ui) + viscs2_(1, 2, vi)*conv_r_(1, 1, ui) + viscs2_(2, 2, vi)*conv_r_(2, 1, ui)) ;
    estif_(vi*4 + 2, ui*3 + 10) += 2.0*nu_*ttimetauMp*(conv_c_(ui)*viscs2_(2, 2, vi) + conv_g_(ui)*viscs2_(2, 2, vi) + viscs2_(0, 2, vi)*conv_r_(0, 2, ui) + viscs2_(1, 2, vi)*conv_r_(1, 2, ui) + viscs2_(2, 2, vi)*conv_r_(2, 2, ui)) ;
#endif

#ifdef FLUID3_IS_TERM6
    /* Stabilisierung der Viskosität (-L_visc_u) */
    estif_(vi*4, ui*3 + 8)     += 4.0*(nu_*nu_)*ttimetauMp*(viscs2_(0, 0, ui)*viscs2_(0, 0, vi) + viscs2_(0, 1, ui)*viscs2_(0, 1, vi) + viscs2_(0, 2, ui)*viscs2_(0, 2, vi)) ;
    estif_(vi*4, ui*3 + 9)     += 4.0*(nu_*nu_)*ttimetauMp*(viscs2_(0, 0, vi)*viscs2_(0, 1, ui) + viscs2_(0, 1, vi)*viscs2_(1, 1, ui) + viscs2_(0, 2, vi)*viscs2_(1, 2, ui)) ;
    estif_(vi*4, ui*3 + 10)     += 4.0*(nu_*nu_)*ttimetauMp*(viscs2_(0, 0, vi)*viscs2_(0, 2, ui) + viscs2_(0, 1, vi)*viscs2_(1, 2, ui) + viscs2_(0, 2, vi)*viscs2_(2, 2, ui)) ;
    estif_(vi*4 + 1, ui*3 + 8) += 4.0*(nu_*nu_)*ttimetauMp*(viscs2_(0, 0, ui)*viscs2_(0, 1, vi) + viscs2_(0, 1, ui)*viscs2_(1, 1, vi) + viscs2_(0, 2, ui)*viscs2_(1, 2, vi)) ;
    estif_(vi*4 + 1, ui*3 + 9) += 4.0*(nu_*nu_)*ttimetauMp*(viscs2_(0, 1, ui)*viscs2_(0, 1, vi) + viscs2_(1, 1, ui)*viscs2_(1, 1, vi) + viscs2_(1, 2, ui)*viscs2_(1, 2, vi)) ;
    estif_(vi*4 + 1, ui*3 + 10) += 4.0*(nu_*nu_)*ttimetauMp*(viscs2_(0, 1, vi)*viscs2_(0, 2, ui) + viscs2_(1, 1, vi)*viscs2_(1, 2, ui) + viscs2_(1, 2, vi)*viscs2_(2, 2, ui)) ;
    estif_(vi*4 + 2, ui*3 + 8) += 4.0*(nu_*nu_)*ttimetauMp*(viscs2_(0, 0, ui)*viscs2_(0, 2, vi) + viscs2_(0, 1, ui)*viscs2_(1, 2, vi) + viscs2_(0, 2, ui)*viscs2_(2, 2, vi)) ;
    estif_(vi*4 + 2, ui*3 + 9) += 4.0*(nu_*nu_)*ttimetauMp*(viscs2_(0, 1, ui)*viscs2_(0, 2, vi) + viscs2_(1, 1, ui)*viscs2_(1, 2, vi) + viscs2_(1, 2, ui)*viscs2_(2, 2, vi)) ;
    estif_(vi*4 + 2, ui*3 + 10) += 4.0*(nu_*nu_)*ttimetauMp*(viscs2_(0, 2, ui)*viscs2_(0, 2, vi) + viscs2_(1, 2, ui)*viscs2_(1, 2, vi) + viscs2_(2, 2, ui)*viscs2_(2, 2, vi)) ;
#endif

#ifdef FLUID3_IS_TERM7
    /* Druckterm */
#endif

#ifdef FLUID3_IS_TERM8
    /* Divergenzfreiheit */
    estif_(vi*4 + 3, ui*3 + 8) += timefacfac*funct_p_(vi)*derxyz_(0, ui) ;
    estif_(vi*4 + 3, ui*3 + 9) += timefacfac*funct_p_(vi)*derxyz_(1, ui) ;
    estif_(vi*4 + 3, ui*3 + 10) += timefacfac*funct_p_(vi)*derxyz_(2, ui) ;
#endif

#ifdef FLUID3_IS_TERM9
    /* Kontinuitätsstabilisierung */
    estif_(vi*4, ui*3 + 8)     += (thsl*thsl)*tau_C*derxyz_(0, ui)*derxyz_(0, vi) ;
    estif_(vi*4, ui*3 + 9)     += (thsl*thsl)*tau_C*derxyz_(0, vi)*derxyz_(1, ui) ;
    estif_(vi*4, ui*3 + 10)     += (thsl*thsl)*tau_C*derxyz_(0, vi)*derxyz_(2, ui) ;
    estif_(vi*4 + 1, ui*3 + 8) += (thsl*thsl)*tau_C*derxyz_(0, ui)*derxyz_(1, vi) ;
    estif_(vi*4 + 1, ui*3 + 9) += (thsl*thsl)*tau_C*derxyz_(1, ui)*derxyz_(1, vi) ;
    estif_(vi*4 + 1, ui*3 + 10) += (thsl*thsl)*tau_C*derxyz_(1, vi)*derxyz_(2, ui) ;
    estif_(vi*4 + 2, ui*3 + 8) += (thsl*thsl)*tau_C*derxyz_(0, ui)*derxyz_(2, vi) ;
    estif_(vi*4 + 2, ui*3 + 9) += (thsl*thsl)*tau_C*derxyz_(1, ui)*derxyz_(2, vi) ;
    estif_(vi*4 + 2, ui*3 + 10) += (thsl*thsl)*tau_C*derxyz_(2, ui)*derxyz_(2, vi) ;
#endif

#ifdef FLUID3_IS_TERM10
    /* Massenterm */
    estif_(vi*4, ui*3 + 8)     += fac*funct_(ui)*funct_(vi) ;
    estif_(vi*4 + 1, ui*3 + 9) += fac*funct_(ui)*funct_(vi) ;
    estif_(vi*4 + 2, ui*3 + 10) += fac*funct_(ui)*funct_(vi) ;
#endif

#ifdef FLUID3_IS_TERM11
    /* Konvektionsstabilisierung */
    estif_(vi*4, ui*3 + 8)     += timetauM*funct_(ui)*(conv_g_(vi) + 2.0*velint_(0)*derxyz_(0, vi) + velint_(1)*derxyz_(1, vi) + velint_(2)*derxyz_(2, vi)) ;
    estif_(vi*4, ui*3 + 9)     += timetauM*funct_(ui)*velint_(0)*derxyz_(1, vi) ;
    estif_(vi*4, ui*3 + 10)     += timetauM*funct_(ui)*velint_(0)*derxyz_(2, vi) ;
    estif_(vi*4 + 1, ui*3 + 8) += timetauM*funct_(ui)*velint_(1)*derxyz_(0, vi) ;
    estif_(vi*4 + 1, ui*3 + 9) += timetauM*funct_(ui)*(conv_g_(vi) + velint_(0)*derxyz_(0, vi) + 2.0*velint_(1)*derxyz_(1, vi) + velint_(2)*derxyz_(2, vi)) ;
    estif_(vi*4 + 1, ui*3 + 10) += timetauM*funct_(ui)*velint_(1)*derxyz_(2, vi) ;
    estif_(vi*4 + 2, ui*3 + 8) += timetauM*funct_(ui)*velint_(2)*derxyz_(0, vi) ;
    estif_(vi*4 + 2, ui*3 + 9) += timetauM*funct_(ui)*velint_(2)*derxyz_(1, vi) ;
    estif_(vi*4 + 2, ui*3 + 10) += timetauM*funct_(ui)*(conv_g_(vi) + velint_(0)*derxyz_(0, vi) + velint_(1)*derxyz_(1, vi) + 2.0*velint_(2)*derxyz_(2, vi)) ;
#endif

#ifdef FLUID3_IS_TERM12
    /* Viskositätsstabilisierung */
    estif_(vi*4, ui*3 + 8)     += 2.0*nu_*timetauMp*funct_(ui)*viscs2_(0, 0, vi) ;
    estif_(vi*4, ui*3 + 9)     += 2.0*nu_*timetauMp*funct_(ui)*viscs2_(0, 1, vi) ;
    estif_(vi*4, ui*3 + 10)     += 2.0*nu_*timetauMp*funct_(ui)*viscs2_(0, 2, vi) ;
    estif_(vi*4 + 1, ui*3 + 8) += 2.0*nu_*timetauMp*funct_(ui)*viscs2_(0, 1, vi) ;
    estif_(vi*4 + 1, ui*3 + 9) += 2.0*nu_*timetauMp*funct_(ui)*viscs2_(1, 1, vi) ;
    estif_(vi*4 + 1, ui*3 + 10) += 2.0*nu_*timetauMp*funct_(ui)*viscs2_(1, 2, vi) ;
    estif_(vi*4 + 2, ui*3 + 8) += 2.0*nu_*timetauMp*funct_(ui)*viscs2_(0, 2, vi) ;
    estif_(vi*4 + 2, ui*3 + 9) += 2.0*nu_*timetauMp*funct_(ui)*viscs2_(1, 2, vi) ;
    estif_(vi*4 + 2, ui*3 + 10) += 2.0*nu_*timetauMp*funct_(ui)*viscs2_(2, 2, vi) ;
#endif

#ifdef FLUID3_IS_TERM13
    /* Quellterm der rechten Seite */
#endif

#ifdef FLUID3_IS_TERM14
    /* Konvektionsstabilisierung */
    estif_(vi*4, ui*3 + 8)     += -(timetauM*funct_(ui)*rhsint_(0)*derxyz_(0, vi)) ;
    estif_(vi*4, ui*3 + 9)     += -(timetauM*funct_(ui)*rhsint_(0)*derxyz_(1, vi)) ;
    estif_(vi*4, ui*3 + 10)     += -(timetauM*funct_(ui)*rhsint_(0)*derxyz_(2, vi)) ;
    estif_(vi*4 + 1, ui*3 + 8) += -(timetauM*funct_(ui)*rhsint_(1)*derxyz_(0, vi)) ;
    estif_(vi*4 + 1, ui*3 + 9) += -(timetauM*funct_(ui)*rhsint_(1)*derxyz_(1, vi)) ;
    estif_(vi*4 + 1, ui*3 + 10) += -(timetauM*funct_(ui)*rhsint_(1)*derxyz_(2, vi)) ;
    estif_(vi*4 + 2, ui*3 + 8) += -(timetauM*funct_(ui)*rhsint_(2)*derxyz_(0, vi)) ;
    estif_(vi*4 + 2, ui*3 + 9) += -(timetauM*funct_(ui)*rhsint_(2)*derxyz_(1, vi)) ;
    estif_(vi*4 + 2, ui*3 + 10) += -(timetauM*funct_(ui)*rhsint_(2)*derxyz_(2, vi)) ;
#endif

#ifdef FLUID3_IS_TERM15
    /* Viskositätsstabilisierung */
#endif

  }
}

for (vi=8; vi<iel; ++vi)
{
  for (ui=0; ui<8; ++ui)
  {
#ifdef FLUID3_IS_TERM1
    /* Konvektionsterm */
    estif_(vi*3 + 8, ui*4)     += timefacfac*funct_(vi)*(conv_c_(ui) + conv_g_(ui) + conv_r_(0, 0, ui)) ;
    estif_(vi*3 + 8, ui*4 + 1) += timefacfac*funct_(vi)*conv_r_(0, 1, ui) ;
    estif_(vi*3 + 8, ui*4 + 2) += timefacfac*funct_(vi)*conv_r_(0, 2, ui) ;
    estif_(vi*3 + 9, ui*4)     += timefacfac*funct_(vi)*conv_r_(1, 0, ui) ;
    estif_(vi*3 + 9, ui*4 + 1) += timefacfac*funct_(vi)*(conv_c_(ui) + conv_g_(ui) + conv_r_(1, 1, ui)) ;
    estif_(vi*3 + 9, ui*4 + 2) += timefacfac*funct_(vi)*conv_r_(1, 2, ui) ;
    estif_(vi*3 + 10, ui*4)     += timefacfac*funct_(vi)*conv_r_(2, 0, ui) ;
    estif_(vi*3 + 10, ui*4 + 1) += timefacfac*funct_(vi)*conv_r_(2, 1, ui) ;
    estif_(vi*3 + 10, ui*4 + 2) += timefacfac*funct_(vi)*(conv_c_(ui) + conv_g_(ui) + conv_r_(2, 2, ui)) ;
#endif

#ifdef FLUID3_IS_TERM2
    /* Stabilisierung der Konvektion ( L_conv_u) */
    estif_(vi*3 + 8, ui*4)     += ttimetauM*(conv_c_(ui)*conv_c_(vi) + conv_c_(ui)*conv_g_(vi) + conv_c_(vi)*conv_g_(ui) + conv_g_(ui)*conv_g_(vi) + conv_c_(vi)*conv_r_(0, 0, ui) + conv_g_(vi)*conv_r_(0, 0, ui) + velint_(0)*derxyz_(0, vi)*conv_r_(0, 0, ui) + velint_(1)*derxyz_(0, vi)*conv_r_(0, 1, ui) + velint_(2)*derxyz_(0, vi)*conv_r_(0, 2, ui) - gridvint_(0)*derxyz_(0, vi)*conv_r_(0, 0, ui) - gridvint_(1)*derxyz_(0, vi)*conv_r_(0, 1, ui) - gridvint_(2)*derxyz_(0, vi)*conv_r_(0, 2, ui)) ;
    estif_(vi*3 + 8, ui*4 + 1) += ttimetauM*(conv_c_(vi)*conv_r_(0, 1, ui) + conv_g_(vi)*conv_r_(0, 1, ui) + velint_(0)*derxyz_(1, vi)*conv_r_(0, 0, ui) + velint_(1)*derxyz_(1, vi)*conv_r_(0, 1, ui) + velint_(2)*derxyz_(1, vi)*conv_r_(0, 2, ui) - gridvint_(0)*derxyz_(1, vi)*conv_r_(0, 0, ui) - gridvint_(1)*derxyz_(1, vi)*conv_r_(0, 1, ui) - gridvint_(2)*derxyz_(1, vi)*conv_r_(0, 2, ui)) ;
    estif_(vi*3 + 8, ui*4 + 2) += ttimetauM*(conv_c_(vi)*conv_r_(0, 2, ui) + conv_g_(vi)*conv_r_(0, 2, ui) + velint_(0)*derxyz_(2, vi)*conv_r_(0, 0, ui) + velint_(1)*derxyz_(2, vi)*conv_r_(0, 1, ui) + velint_(2)*derxyz_(2, vi)*conv_r_(0, 2, ui) - gridvint_(0)*derxyz_(2, vi)*conv_r_(0, 0, ui) - gridvint_(1)*derxyz_(2, vi)*conv_r_(0, 1, ui) - gridvint_(2)*derxyz_(2, vi)*conv_r_(0, 2, ui)) ;
    estif_(vi*3 + 9, ui*4)     += ttimetauM*(conv_c_(vi)*conv_r_(1, 0, ui) + conv_g_(vi)*conv_r_(1, 0, ui) + velint_(0)*derxyz_(0, vi)*conv_r_(1, 0, ui) + velint_(1)*derxyz_(0, vi)*conv_r_(1, 1, ui) + velint_(2)*derxyz_(0, vi)*conv_r_(1, 2, ui) - gridvint_(0)*derxyz_(0, vi)*conv_r_(1, 0, ui) - gridvint_(1)*derxyz_(0, vi)*conv_r_(1, 1, ui) - gridvint_(2)*derxyz_(0, vi)*conv_r_(1, 2, ui)) ;
    estif_(vi*3 + 9, ui*4 + 1) += ttimetauM*(conv_c_(ui)*conv_c_(vi) + conv_c_(ui)*conv_g_(vi) + conv_c_(vi)*conv_g_(ui) + conv_g_(ui)*conv_g_(vi) + conv_c_(vi)*conv_r_(1, 1, ui) + conv_g_(vi)*conv_r_(1, 1, ui) + velint_(0)*derxyz_(1, vi)*conv_r_(1, 0, ui) + velint_(1)*derxyz_(1, vi)*conv_r_(1, 1, ui) + velint_(2)*derxyz_(1, vi)*conv_r_(1, 2, ui) - gridvint_(0)*derxyz_(1, vi)*conv_r_(1, 0, ui) - gridvint_(1)*derxyz_(1, vi)*conv_r_(1, 1, ui) - gridvint_(2)*derxyz_(1, vi)*conv_r_(1, 2, ui)) ;
    estif_(vi*3 + 9, ui*4 + 2) += ttimetauM*(conv_c_(vi)*conv_r_(1, 2, ui) + conv_g_(vi)*conv_r_(1, 2, ui) + velint_(0)*derxyz_(2, vi)*conv_r_(1, 0, ui) + velint_(1)*derxyz_(2, vi)*conv_r_(1, 1, ui) + velint_(2)*derxyz_(2, vi)*conv_r_(1, 2, ui) - gridvint_(0)*derxyz_(2, vi)*conv_r_(1, 0, ui) - gridvint_(1)*derxyz_(2, vi)*conv_r_(1, 1, ui) - gridvint_(2)*derxyz_(2, vi)*conv_r_(1, 2, ui)) ;
    estif_(vi*3 + 10, ui*4)     += ttimetauM*(conv_c_(vi)*conv_r_(2, 0, ui) + conv_g_(vi)*conv_r_(2, 0, ui) + velint_(0)*derxyz_(0, vi)*conv_r_(2, 0, ui) + velint_(1)*derxyz_(0, vi)*conv_r_(2, 1, ui) + velint_(2)*derxyz_(0, vi)*conv_r_(2, 2, ui) - gridvint_(0)*derxyz_(0, vi)*conv_r_(2, 0, ui) - gridvint_(1)*derxyz_(0, vi)*conv_r_(2, 1, ui) - gridvint_(2)*derxyz_(0, vi)*conv_r_(2, 2, ui)) ;
    estif_(vi*3 + 10, ui*4 + 1) += ttimetauM*(conv_c_(vi)*conv_r_(2, 1, ui) + conv_g_(vi)*conv_r_(2, 1, ui) + velint_(0)*derxyz_(1, vi)*conv_r_(2, 0, ui) + velint_(1)*derxyz_(1, vi)*conv_r_(2, 1, ui) + velint_(2)*derxyz_(1, vi)*conv_r_(2, 2, ui) - gridvint_(0)*derxyz_(1, vi)*conv_r_(2, 0, ui) - gridvint_(1)*derxyz_(1, vi)*conv_r_(2, 1, ui) - gridvint_(2)*derxyz_(1, vi)*conv_r_(2, 2, ui)) ;
    estif_(vi*3 + 10, ui*4 + 2) += ttimetauM*(conv_c_(ui)*conv_c_(vi) + conv_c_(ui)*conv_g_(vi) + conv_c_(vi)*conv_g_(ui) + conv_g_(ui)*conv_g_(vi) + conv_c_(vi)*conv_r_(2, 2, ui) + conv_g_(vi)*conv_r_(2, 2, ui) + velint_(0)*derxyz_(2, vi)*conv_r_(2, 0, ui) + velint_(1)*derxyz_(2, vi)*conv_r_(2, 1, ui) + velint_(2)*derxyz_(2, vi)*conv_r_(2, 2, ui) - gridvint_(0)*derxyz_(2, vi)*conv_r_(2, 0, ui) - gridvint_(1)*derxyz_(2, vi)*conv_r_(2, 1, ui) - gridvint_(2)*derxyz_(2, vi)*conv_r_(2, 2, ui)) ;
#endif

#ifdef FLUID3_IS_TERM3
    /* Stabilisierung der Konvektion (-L_visc_u) */
    estif_(vi*3 + 8, ui*4)     += -2.0*nu_*ttimetauM*(-(conv_c_(vi)*viscs2_(0, 0, ui)) - conv_g_(vi)*viscs2_(0, 0, ui) + funct_(ui)*visc_old_(0)*derxyz_(0, vi)) ;
    estif_(vi*3 + 8, ui*4 + 1) += -2.0*nu_*ttimetauM*(-(conv_c_(vi)*viscs2_(0, 1, ui)) - conv_g_(vi)*viscs2_(0, 1, ui) + funct_(ui)*visc_old_(0)*derxyz_(1, vi)) ;
    estif_(vi*3 + 8, ui*4 + 2) += -2.0*nu_*ttimetauM*(-(conv_c_(vi)*viscs2_(0, 2, ui)) - conv_g_(vi)*viscs2_(0, 2, ui) + funct_(ui)*visc_old_(0)*derxyz_(2, vi)) ;
    estif_(vi*3 + 9, ui*4)     += -2.0*nu_*ttimetauM*(-(conv_c_(vi)*viscs2_(0, 1, ui)) - conv_g_(vi)*viscs2_(0, 1, ui) + funct_(ui)*visc_old_(1)*derxyz_(0, vi)) ;
    estif_(vi*3 + 9, ui*4 + 1) += -2.0*nu_*ttimetauM*(-(conv_c_(vi)*viscs2_(1, 1, ui)) - conv_g_(vi)*viscs2_(1, 1, ui) + funct_(ui)*visc_old_(1)*derxyz_(1, vi)) ;
    estif_(vi*3 + 9, ui*4 + 2) += -2.0*nu_*ttimetauM*(-(conv_c_(vi)*viscs2_(1, 2, ui)) - conv_g_(vi)*viscs2_(1, 2, ui) + funct_(ui)*visc_old_(1)*derxyz_(2, vi)) ;
    estif_(vi*3 + 10, ui*4)     += -2.0*nu_*ttimetauM*(-(conv_c_(vi)*viscs2_(0, 2, ui)) - conv_g_(vi)*viscs2_(0, 2, ui) + funct_(ui)*visc_old_(2)*derxyz_(0, vi)) ;
    estif_(vi*3 + 10, ui*4 + 1) += -2.0*nu_*ttimetauM*(-(conv_c_(vi)*viscs2_(1, 2, ui)) - conv_g_(vi)*viscs2_(1, 2, ui) + funct_(ui)*visc_old_(2)*derxyz_(1, vi)) ;
    estif_(vi*3 + 10, ui*4 + 2) += -2.0*nu_*ttimetauM*(-(conv_c_(vi)*viscs2_(2, 2, ui)) - conv_g_(vi)*viscs2_(2, 2, ui) + funct_(ui)*visc_old_(2)*derxyz_(2, vi)) ;
#endif

#ifdef FLUID3_IS_TERM4
    /* Viskositätsterm */
    estif_(vi*3 + 8, ui*4)     += nu_*timefacfac*(2.0*derxyz_(0, ui)*derxyz_(0, vi) + derxyz_(1, ui)*derxyz_(1, vi) + derxyz_(2, ui)*derxyz_(2, vi)) ;
    estif_(vi*3 + 8, ui*4 + 1) += nu_*timefacfac*derxyz_(0, ui)*derxyz_(1, vi) ;
    estif_(vi*3 + 8, ui*4 + 2) += nu_*timefacfac*derxyz_(0, ui)*derxyz_(2, vi) ;
    estif_(vi*3 + 9, ui*4)     += nu_*timefacfac*derxyz_(0, vi)*derxyz_(1, ui) ;
    estif_(vi*3 + 9, ui*4 + 1) += nu_*timefacfac*(derxyz_(0, ui)*derxyz_(0, vi) + 2.0*derxyz_(1, ui)*derxyz_(1, vi) + derxyz_(2, ui)*derxyz_(2, vi)) ;
    estif_(vi*3 + 9, ui*4 + 2) += nu_*timefacfac*derxyz_(1, ui)*derxyz_(2, vi) ;
    estif_(vi*3 + 10, ui*4)     += nu_*timefacfac*derxyz_(0, vi)*derxyz_(2, ui) ;
    estif_(vi*3 + 10, ui*4 + 1) += nu_*timefacfac*derxyz_(1, vi)*derxyz_(2, ui) ;
    estif_(vi*3 + 10, ui*4 + 2) += nu_*timefacfac*(derxyz_(0, ui)*derxyz_(0, vi) + derxyz_(1, ui)*derxyz_(1, vi) + 2.0*derxyz_(2, ui)*derxyz_(2, vi)) ;
#endif

#ifdef FLUID3_IS_TERM5
    /* Stabilisierung der Viskosität ( L_conv_u) */
    estif_(vi*3 + 8, ui*4)     += 2.0*nu_*ttimetauMp*(conv_c_(ui)*viscs2_(0, 0, vi) + conv_g_(ui)*viscs2_(0, 0, vi) + viscs2_(0, 0, vi)*conv_r_(0, 0, ui) + viscs2_(0, 1, vi)*conv_r_(1, 0, ui) + viscs2_(0, 2, vi)*conv_r_(2, 0, ui)) ;
    estif_(vi*3 + 8, ui*4 + 1) += 2.0*nu_*ttimetauMp*(conv_c_(ui)*viscs2_(0, 1, vi) + conv_g_(ui)*viscs2_(0, 1, vi) + viscs2_(0, 0, vi)*conv_r_(0, 1, ui) + viscs2_(0, 1, vi)*conv_r_(1, 1, ui) + viscs2_(0, 2, vi)*conv_r_(2, 1, ui)) ;
    estif_(vi*3 + 8, ui*4 + 2) += 2.0*nu_*ttimetauMp*(conv_c_(ui)*viscs2_(0, 2, vi) + conv_g_(ui)*viscs2_(0, 2, vi) + viscs2_(0, 0, vi)*conv_r_(0, 2, ui) + viscs2_(0, 1, vi)*conv_r_(1, 2, ui) + viscs2_(0, 2, vi)*conv_r_(2, 2, ui)) ;
    estif_(vi*3 + 9, ui*4)     += 2.0*nu_*ttimetauMp*(conv_c_(ui)*viscs2_(0, 1, vi) + conv_g_(ui)*viscs2_(0, 1, vi) + viscs2_(0, 1, vi)*conv_r_(0, 0, ui) + viscs2_(1, 1, vi)*conv_r_(1, 0, ui) + viscs2_(1, 2, vi)*conv_r_(2, 0, ui)) ;
    estif_(vi*3 + 9, ui*4 + 1) += 2.0*nu_*ttimetauMp*(conv_c_(ui)*viscs2_(1, 1, vi) + conv_g_(ui)*viscs2_(1, 1, vi) + viscs2_(0, 1, vi)*conv_r_(0, 1, ui) + viscs2_(1, 1, vi)*conv_r_(1, 1, ui) + viscs2_(1, 2, vi)*conv_r_(2, 1, ui)) ;
    estif_(vi*3 + 9, ui*4 + 2) += 2.0*nu_*ttimetauMp*(conv_c_(ui)*viscs2_(1, 2, vi) + conv_g_(ui)*viscs2_(1, 2, vi) + viscs2_(0, 1, vi)*conv_r_(0, 2, ui) + viscs2_(1, 1, vi)*conv_r_(1, 2, ui) + viscs2_(1, 2, vi)*conv_r_(2, 2, ui)) ;
    estif_(vi*3 + 10, ui*4)     += 2.0*nu_*ttimetauMp*(conv_c_(ui)*viscs2_(0, 2, vi) + conv_g_(ui)*viscs2_(0, 2, vi) + viscs2_(0, 2, vi)*conv_r_(0, 0, ui) + viscs2_(1, 2, vi)*conv_r_(1, 0, ui) + viscs2_(2, 2, vi)*conv_r_(2, 0, ui)) ;
    estif_(vi*3 + 10, ui*4 + 1) += 2.0*nu_*ttimetauMp*(conv_c_(ui)*viscs2_(1, 2, vi) + conv_g_(ui)*viscs2_(1, 2, vi) + viscs2_(0, 2, vi)*conv_r_(0, 1, ui) + viscs2_(1, 2, vi)*conv_r_(1, 1, ui) + viscs2_(2, 2, vi)*conv_r_(2, 1, ui)) ;
    estif_(vi*3 + 10, ui*4 + 2) += 2.0*nu_*ttimetauMp*(conv_c_(ui)*viscs2_(2, 2, vi) + conv_g_(ui)*viscs2_(2, 2, vi) + viscs2_(0, 2, vi)*conv_r_(0, 2, ui) + viscs2_(1, 2, vi)*conv_r_(1, 2, ui) + viscs2_(2, 2, vi)*conv_r_(2, 2, ui)) ;
#endif

#ifdef FLUID3_IS_TERM6
    /* Stabilisierung der Viskosität (-L_visc_u) */
    estif_(vi*3 + 8, ui*4)     += 4.0*(nu_*nu_)*ttimetauMp*(viscs2_(0, 0, ui)*viscs2_(0, 0, vi) + viscs2_(0, 1, ui)*viscs2_(0, 1, vi) + viscs2_(0, 2, ui)*viscs2_(0, 2, vi)) ;
    estif_(vi*3 + 8, ui*4 + 1) += 4.0*(nu_*nu_)*ttimetauMp*(viscs2_(0, 0, vi)*viscs2_(0, 1, ui) + viscs2_(0, 1, vi)*viscs2_(1, 1, ui) + viscs2_(0, 2, vi)*viscs2_(1, 2, ui)) ;
    estif_(vi*3 + 8, ui*4 + 2) += 4.0*(nu_*nu_)*ttimetauMp*(viscs2_(0, 0, vi)*viscs2_(0, 2, ui) + viscs2_(0, 1, vi)*viscs2_(1, 2, ui) + viscs2_(0, 2, vi)*viscs2_(2, 2, ui)) ;
    estif_(vi*3 + 9, ui*4)     += 4.0*(nu_*nu_)*ttimetauMp*(viscs2_(0, 0, ui)*viscs2_(0, 1, vi) + viscs2_(0, 1, ui)*viscs2_(1, 1, vi) + viscs2_(0, 2, ui)*viscs2_(1, 2, vi)) ;
    estif_(vi*3 + 9, ui*4 + 1) += 4.0*(nu_*nu_)*ttimetauMp*(viscs2_(0, 1, ui)*viscs2_(0, 1, vi) + viscs2_(1, 1, ui)*viscs2_(1, 1, vi) + viscs2_(1, 2, ui)*viscs2_(1, 2, vi)) ;
    estif_(vi*3 + 9, ui*4 + 2) += 4.0*(nu_*nu_)*ttimetauMp*(viscs2_(0, 1, vi)*viscs2_(0, 2, ui) + viscs2_(1, 1, vi)*viscs2_(1, 2, ui) + viscs2_(1, 2, vi)*viscs2_(2, 2, ui)) ;
    estif_(vi*3 + 10, ui*4)     += 4.0*(nu_*nu_)*ttimetauMp*(viscs2_(0, 0, ui)*viscs2_(0, 2, vi) + viscs2_(0, 1, ui)*viscs2_(1, 2, vi) + viscs2_(0, 2, ui)*viscs2_(2, 2, vi)) ;
    estif_(vi*3 + 10, ui*4 + 1) += 4.0*(nu_*nu_)*ttimetauMp*(viscs2_(0, 1, ui)*viscs2_(0, 2, vi) + viscs2_(1, 1, ui)*viscs2_(1, 2, vi) + viscs2_(1, 2, ui)*viscs2_(2, 2, vi)) ;
    estif_(vi*3 + 10, ui*4 + 2) += 4.0*(nu_*nu_)*ttimetauMp*(viscs2_(0, 2, ui)*viscs2_(0, 2, vi) + viscs2_(1, 2, ui)*viscs2_(1, 2, vi) + viscs2_(2, 2, ui)*viscs2_(2, 2, vi)) ;
#endif

#ifdef FLUID3_IS_TERM7
    /* Druckterm */
    estif_(vi*3 + 8, ui*4 + 3) += -(timefacfac*funct_p_(ui)*derxyz_(0, vi)) ;
    estif_(vi*3 + 9, ui*4 + 3) += -(timefacfac*funct_p_(ui)*derxyz_(1, vi)) ;
    estif_(vi*3 + 10, ui*4 + 3) += -(timefacfac*funct_p_(ui)*derxyz_(2, vi)) ;
#endif

#ifdef FLUID3_IS_TERM8
    /* Divergenzfreiheit */
#endif

#ifdef FLUID3_IS_TERM9
    /* Kontinuitätsstabilisierung */
    estif_(vi*3 + 8, ui*4)     += (thsl*thsl)*tau_C*derxyz_(0, ui)*derxyz_(0, vi) ;
    estif_(vi*3 + 8, ui*4 + 1) += (thsl*thsl)*tau_C*derxyz_(0, vi)*derxyz_(1, ui) ;
    estif_(vi*3 + 8, ui*4 + 2) += (thsl*thsl)*tau_C*derxyz_(0, vi)*derxyz_(2, ui) ;
    estif_(vi*3 + 9, ui*4)     += (thsl*thsl)*tau_C*derxyz_(0, ui)*derxyz_(1, vi) ;
    estif_(vi*3 + 9, ui*4 + 1) += (thsl*thsl)*tau_C*derxyz_(1, ui)*derxyz_(1, vi) ;
    estif_(vi*3 + 9, ui*4 + 2) += (thsl*thsl)*tau_C*derxyz_(1, vi)*derxyz_(2, ui) ;
    estif_(vi*3 + 10, ui*4)     += (thsl*thsl)*tau_C*derxyz_(0, ui)*derxyz_(2, vi) ;
    estif_(vi*3 + 10, ui*4 + 1) += (thsl*thsl)*tau_C*derxyz_(1, ui)*derxyz_(2, vi) ;
    estif_(vi*3 + 10, ui*4 + 2) += (thsl*thsl)*tau_C*derxyz_(2, ui)*derxyz_(2, vi) ;
#endif

#ifdef FLUID3_IS_TERM10
    /* Massenterm */
    estif_(vi*3 + 8, ui*4)     += fac*funct_(ui)*funct_(vi) ;
    estif_(vi*3 + 9, ui*4 + 1) += fac*funct_(ui)*funct_(vi) ;
    estif_(vi*3 + 10, ui*4 + 2) += fac*funct_(ui)*funct_(vi) ;
#endif

#ifdef FLUID3_IS_TERM11
    /* Konvektionsstabilisierung */
    estif_(vi*3 + 8, ui*4)     += timetauM*funct_(ui)*(conv_g_(vi) + 2.0*velint_(0)*derxyz_(0, vi) + velint_(1)*derxyz_(1, vi) + velint_(2)*derxyz_(2, vi)) ;
    estif_(vi*3 + 8, ui*4 + 1) += timetauM*funct_(ui)*velint_(0)*derxyz_(1, vi) ;
    estif_(vi*3 + 8, ui*4 + 2) += timetauM*funct_(ui)*velint_(0)*derxyz_(2, vi) ;
    estif_(vi*3 + 9, ui*4)     += timetauM*funct_(ui)*velint_(1)*derxyz_(0, vi) ;
    estif_(vi*3 + 9, ui*4 + 1) += timetauM*funct_(ui)*(conv_g_(vi) + velint_(0)*derxyz_(0, vi) + 2.0*velint_(1)*derxyz_(1, vi) + velint_(2)*derxyz_(2, vi)) ;
    estif_(vi*3 + 9, ui*4 + 2) += timetauM*funct_(ui)*velint_(1)*derxyz_(2, vi) ;
    estif_(vi*3 + 10, ui*4)     += timetauM*funct_(ui)*velint_(2)*derxyz_(0, vi) ;
    estif_(vi*3 + 10, ui*4 + 1) += timetauM*funct_(ui)*velint_(2)*derxyz_(1, vi) ;
    estif_(vi*3 + 10, ui*4 + 2) += timetauM*funct_(ui)*(conv_g_(vi) + velint_(0)*derxyz_(0, vi) + velint_(1)*derxyz_(1, vi) + 2.0*velint_(2)*derxyz_(2, vi)) ;
#endif

#ifdef FLUID3_IS_TERM12
    /* Viskositätsstabilisierung */
    estif_(vi*3 + 8, ui*4)     += 2.0*nu_*timetauMp*funct_(ui)*viscs2_(0, 0, vi) ;
    estif_(vi*3 + 8, ui*4 + 1) += 2.0*nu_*timetauMp*funct_(ui)*viscs2_(0, 1, vi) ;
    estif_(vi*3 + 8, ui*4 + 2) += 2.0*nu_*timetauMp*funct_(ui)*viscs2_(0, 2, vi) ;
    estif_(vi*3 + 9, ui*4)     += 2.0*nu_*timetauMp*funct_(ui)*viscs2_(0, 1, vi) ;
    estif_(vi*3 + 9, ui*4 + 1) += 2.0*nu_*timetauMp*funct_(ui)*viscs2_(1, 1, vi) ;
    estif_(vi*3 + 9, ui*4 + 2) += 2.0*nu_*timetauMp*funct_(ui)*viscs2_(1, 2, vi) ;
    estif_(vi*3 + 10, ui*4)     += 2.0*nu_*timetauMp*funct_(ui)*viscs2_(0, 2, vi) ;
    estif_(vi*3 + 10, ui*4 + 1) += 2.0*nu_*timetauMp*funct_(ui)*viscs2_(1, 2, vi) ;
    estif_(vi*3 + 10, ui*4 + 2) += 2.0*nu_*timetauMp*funct_(ui)*viscs2_(2, 2, vi) ;
#endif

#ifdef FLUID3_IS_TERM13
    /* Quellterm der rechten Seite */
#endif

#ifdef FLUID3_IS_TERM14
    /* Konvektionsstabilisierung */
    estif_(vi*3 + 8, ui*4)     += -(timetauM*funct_(ui)*rhsint_(0)*derxyz_(0, vi)) ;
    estif_(vi*3 + 8, ui*4 + 1) += -(timetauM*funct_(ui)*rhsint_(0)*derxyz_(1, vi)) ;
    estif_(vi*3 + 8, ui*4 + 2) += -(timetauM*funct_(ui)*rhsint_(0)*derxyz_(2, vi)) ;
    estif_(vi*3 + 9, ui*4)     += -(timetauM*funct_(ui)*rhsint_(1)*derxyz_(0, vi)) ;
    estif_(vi*3 + 9, ui*4 + 1) += -(timetauM*funct_(ui)*rhsint_(1)*derxyz_(1, vi)) ;
    estif_(vi*3 + 9, ui*4 + 2) += -(timetauM*funct_(ui)*rhsint_(1)*derxyz_(2, vi)) ;
    estif_(vi*3 + 10, ui*4)     += -(timetauM*funct_(ui)*rhsint_(2)*derxyz_(0, vi)) ;
    estif_(vi*3 + 10, ui*4 + 1) += -(timetauM*funct_(ui)*rhsint_(2)*derxyz_(1, vi)) ;
    estif_(vi*3 + 10, ui*4 + 2) += -(timetauM*funct_(ui)*rhsint_(2)*derxyz_(2, vi)) ;
#endif

#ifdef FLUID3_IS_TERM15
    /* Viskositätsstabilisierung */
#endif

  }
}

for (vi=8; vi<iel; ++vi)
{
  for (ui=8; ui<iel; ++ui)
  {
#ifdef FLUID3_IS_TERM1
    /* Konvektionsterm */
    estif_(vi*3 + 8, ui*3 + 8) += timefacfac*funct_(vi)*(conv_c_(ui) + conv_g_(ui) + conv_r_(0, 0, ui)) ;
    estif_(vi*3 + 8, ui*3 + 9) += timefacfac*funct_(vi)*conv_r_(0, 1, ui) ;
    estif_(vi*3 + 8, ui*3 + 10) += timefacfac*funct_(vi)*conv_r_(0, 2, ui) ;
    estif_(vi*3 + 9, ui*3 + 8) += timefacfac*funct_(vi)*conv_r_(1, 0, ui) ;
    estif_(vi*3 + 9, ui*3 + 9) += timefacfac*funct_(vi)*(conv_c_(ui) + conv_g_(ui) + conv_r_(1, 1, ui)) ;
    estif_(vi*3 + 9, ui*3 + 10) += timefacfac*funct_(vi)*conv_r_(1, 2, ui) ;
    estif_(vi*3 + 10, ui*3 + 8) += timefacfac*funct_(vi)*conv_r_(2, 0, ui) ;
    estif_(vi*3 + 10, ui*3 + 9) += timefacfac*funct_(vi)*conv_r_(2, 1, ui) ;
    estif_(vi*3 + 10, ui*3 + 10) += timefacfac*funct_(vi)*(conv_c_(ui) + conv_g_(ui) + conv_r_(2, 2, ui)) ;
#endif

#ifdef FLUID3_IS_TERM2
    /* Stabilisierung der Konvektion ( L_conv_u) */
    estif_(vi*3 + 8, ui*3 + 8) += ttimetauM*(conv_c_(ui)*conv_c_(vi) + conv_c_(ui)*conv_g_(vi) + conv_c_(vi)*conv_g_(ui) + conv_g_(ui)*conv_g_(vi) + conv_c_(vi)*conv_r_(0, 0, ui) + conv_g_(vi)*conv_r_(0, 0, ui) + velint_(0)*derxyz_(0, vi)*conv_r_(0, 0, ui) + velint_(1)*derxyz_(0, vi)*conv_r_(0, 1, ui) + velint_(2)*derxyz_(0, vi)*conv_r_(0, 2, ui) - gridvint_(0)*derxyz_(0, vi)*conv_r_(0, 0, ui) - gridvint_(1)*derxyz_(0, vi)*conv_r_(0, 1, ui) - gridvint_(2)*derxyz_(0, vi)*conv_r_(0, 2, ui)) ;
    estif_(vi*3 + 8, ui*3 + 9) += ttimetauM*(conv_c_(vi)*conv_r_(0, 1, ui) + conv_g_(vi)*conv_r_(0, 1, ui) + velint_(0)*derxyz_(1, vi)*conv_r_(0, 0, ui) + velint_(1)*derxyz_(1, vi)*conv_r_(0, 1, ui) + velint_(2)*derxyz_(1, vi)*conv_r_(0, 2, ui) - gridvint_(0)*derxyz_(1, vi)*conv_r_(0, 0, ui) - gridvint_(1)*derxyz_(1, vi)*conv_r_(0, 1, ui) - gridvint_(2)*derxyz_(1, vi)*conv_r_(0, 2, ui)) ;
    estif_(vi*3 + 8, ui*3 + 10) += ttimetauM*(conv_c_(vi)*conv_r_(0, 2, ui) + conv_g_(vi)*conv_r_(0, 2, ui) + velint_(0)*derxyz_(2, vi)*conv_r_(0, 0, ui) + velint_(1)*derxyz_(2, vi)*conv_r_(0, 1, ui) + velint_(2)*derxyz_(2, vi)*conv_r_(0, 2, ui) - gridvint_(0)*derxyz_(2, vi)*conv_r_(0, 0, ui) - gridvint_(1)*derxyz_(2, vi)*conv_r_(0, 1, ui) - gridvint_(2)*derxyz_(2, vi)*conv_r_(0, 2, ui)) ;
    estif_(vi*3 + 9, ui*3 + 8) += ttimetauM*(conv_c_(vi)*conv_r_(1, 0, ui) + conv_g_(vi)*conv_r_(1, 0, ui) + velint_(0)*derxyz_(0, vi)*conv_r_(1, 0, ui) + velint_(1)*derxyz_(0, vi)*conv_r_(1, 1, ui) + velint_(2)*derxyz_(0, vi)*conv_r_(1, 2, ui) - gridvint_(0)*derxyz_(0, vi)*conv_r_(1, 0, ui) - gridvint_(1)*derxyz_(0, vi)*conv_r_(1, 1, ui) - gridvint_(2)*derxyz_(0, vi)*conv_r_(1, 2, ui)) ;
    estif_(vi*3 + 9, ui*3 + 9) += ttimetauM*(conv_c_(ui)*conv_c_(vi) + conv_c_(ui)*conv_g_(vi) + conv_c_(vi)*conv_g_(ui) + conv_g_(ui)*conv_g_(vi) + conv_c_(vi)*conv_r_(1, 1, ui) + conv_g_(vi)*conv_r_(1, 1, ui) + velint_(0)*derxyz_(1, vi)*conv_r_(1, 0, ui) + velint_(1)*derxyz_(1, vi)*conv_r_(1, 1, ui) + velint_(2)*derxyz_(1, vi)*conv_r_(1, 2, ui) - gridvint_(0)*derxyz_(1, vi)*conv_r_(1, 0, ui) - gridvint_(1)*derxyz_(1, vi)*conv_r_(1, 1, ui) - gridvint_(2)*derxyz_(1, vi)*conv_r_(1, 2, ui)) ;
    estif_(vi*3 + 9, ui*3 + 10) += ttimetauM*(conv_c_(vi)*conv_r_(1, 2, ui) + conv_g_(vi)*conv_r_(1, 2, ui) + velint_(0)*derxyz_(2, vi)*conv_r_(1, 0, ui) + velint_(1)*derxyz_(2, vi)*conv_r_(1, 1, ui) + velint_(2)*derxyz_(2, vi)*conv_r_(1, 2, ui) - gridvint_(0)*derxyz_(2, vi)*conv_r_(1, 0, ui) - gridvint_(1)*derxyz_(2, vi)*conv_r_(1, 1, ui) - gridvint_(2)*derxyz_(2, vi)*conv_r_(1, 2, ui)) ;
    estif_(vi*3 + 10, ui*3 + 8) += ttimetauM*(conv_c_(vi)*conv_r_(2, 0, ui) + conv_g_(vi)*conv_r_(2, 0, ui) + velint_(0)*derxyz_(0, vi)*conv_r_(2, 0, ui) + velint_(1)*derxyz_(0, vi)*conv_r_(2, 1, ui) + velint_(2)*derxyz_(0, vi)*conv_r_(2, 2, ui) - gridvint_(0)*derxyz_(0, vi)*conv_r_(2, 0, ui) - gridvint_(1)*derxyz_(0, vi)*conv_r_(2, 1, ui) - gridvint_(2)*derxyz_(0, vi)*conv_r_(2, 2, ui)) ;
    estif_(vi*3 + 10, ui*3 + 9) += ttimetauM*(conv_c_(vi)*conv_r_(2, 1, ui) + conv_g_(vi)*conv_r_(2, 1, ui) + velint_(0)*derxyz_(1, vi)*conv_r_(2, 0, ui) + velint_(1)*derxyz_(1, vi)*conv_r_(2, 1, ui) + velint_(2)*derxyz_(1, vi)*conv_r_(2, 2, ui) - gridvint_(0)*derxyz_(1, vi)*conv_r_(2, 0, ui) - gridvint_(1)*derxyz_(1, vi)*conv_r_(2, 1, ui) - gridvint_(2)*derxyz_(1, vi)*conv_r_(2, 2, ui)) ;
    estif_(vi*3 + 10, ui*3 + 10) += ttimetauM*(conv_c_(ui)*conv_c_(vi) + conv_c_(ui)*conv_g_(vi) + conv_c_(vi)*conv_g_(ui) + conv_g_(ui)*conv_g_(vi) + conv_c_(vi)*conv_r_(2, 2, ui) + conv_g_(vi)*conv_r_(2, 2, ui) + velint_(0)*derxyz_(2, vi)*conv_r_(2, 0, ui) + velint_(1)*derxyz_(2, vi)*conv_r_(2, 1, ui) + velint_(2)*derxyz_(2, vi)*conv_r_(2, 2, ui) - gridvint_(0)*derxyz_(2, vi)*conv_r_(2, 0, ui) - gridvint_(1)*derxyz_(2, vi)*conv_r_(2, 1, ui) - gridvint_(2)*derxyz_(2, vi)*conv_r_(2, 2, ui)) ;
#endif

#ifdef FLUID3_IS_TERM3
    /* Stabilisierung der Konvektion (-L_visc_u) */
    estif_(vi*3 + 8, ui*3 + 8) += -2.0*nu_*ttimetauM*(-(conv_c_(vi)*viscs2_(0, 0, ui)) - conv_g_(vi)*viscs2_(0, 0, ui) + funct_(ui)*visc_old_(0)*derxyz_(0, vi)) ;
    estif_(vi*3 + 8, ui*3 + 9) += -2.0*nu_*ttimetauM*(-(conv_c_(vi)*viscs2_(0, 1, ui)) - conv_g_(vi)*viscs2_(0, 1, ui) + funct_(ui)*visc_old_(0)*derxyz_(1, vi)) ;
    estif_(vi*3 + 8, ui*3 + 10) += -2.0*nu_*ttimetauM*(-(conv_c_(vi)*viscs2_(0, 2, ui)) - conv_g_(vi)*viscs2_(0, 2, ui) + funct_(ui)*visc_old_(0)*derxyz_(2, vi)) ;
    estif_(vi*3 + 9, ui*3 + 8) += -2.0*nu_*ttimetauM*(-(conv_c_(vi)*viscs2_(0, 1, ui)) - conv_g_(vi)*viscs2_(0, 1, ui) + funct_(ui)*visc_old_(1)*derxyz_(0, vi)) ;
    estif_(vi*3 + 9, ui*3 + 9) += -2.0*nu_*ttimetauM*(-(conv_c_(vi)*viscs2_(1, 1, ui)) - conv_g_(vi)*viscs2_(1, 1, ui) + funct_(ui)*visc_old_(1)*derxyz_(1, vi)) ;
    estif_(vi*3 + 9, ui*3 + 10) += -2.0*nu_*ttimetauM*(-(conv_c_(vi)*viscs2_(1, 2, ui)) - conv_g_(vi)*viscs2_(1, 2, ui) + funct_(ui)*visc_old_(1)*derxyz_(2, vi)) ;
    estif_(vi*3 + 10, ui*3 + 8) += -2.0*nu_*ttimetauM*(-(conv_c_(vi)*viscs2_(0, 2, ui)) - conv_g_(vi)*viscs2_(0, 2, ui) + funct_(ui)*visc_old_(2)*derxyz_(0, vi)) ;
    estif_(vi*3 + 10, ui*3 + 9) += -2.0*nu_*ttimetauM*(-(conv_c_(vi)*viscs2_(1, 2, ui)) - conv_g_(vi)*viscs2_(1, 2, ui) + funct_(ui)*visc_old_(2)*derxyz_(1, vi)) ;
    estif_(vi*3 + 10, ui*3 + 10) += -2.0*nu_*ttimetauM*(-(conv_c_(vi)*viscs2_(2, 2, ui)) - conv_g_(vi)*viscs2_(2, 2, ui) + funct_(ui)*visc_old_(2)*derxyz_(2, vi)) ;
#endif

#ifdef FLUID3_IS_TERM4
    /* Viskositätsterm */
    estif_(vi*3 + 8, ui*3 + 8) += nu_*timefacfac*(2.0*derxyz_(0, ui)*derxyz_(0, vi) + derxyz_(1, ui)*derxyz_(1, vi) + derxyz_(2, ui)*derxyz_(2, vi)) ;
    estif_(vi*3 + 8, ui*3 + 9) += nu_*timefacfac*derxyz_(0, ui)*derxyz_(1, vi) ;
    estif_(vi*3 + 8, ui*3 + 10) += nu_*timefacfac*derxyz_(0, ui)*derxyz_(2, vi) ;
    estif_(vi*3 + 9, ui*3 + 8) += nu_*timefacfac*derxyz_(0, vi)*derxyz_(1, ui) ;
    estif_(vi*3 + 9, ui*3 + 9) += nu_*timefacfac*(derxyz_(0, ui)*derxyz_(0, vi) + 2.0*derxyz_(1, ui)*derxyz_(1, vi) + derxyz_(2, ui)*derxyz_(2, vi)) ;
    estif_(vi*3 + 9, ui*3 + 10) += nu_*timefacfac*derxyz_(1, ui)*derxyz_(2, vi) ;
    estif_(vi*3 + 10, ui*3 + 8) += nu_*timefacfac*derxyz_(0, vi)*derxyz_(2, ui) ;
    estif_(vi*3 + 10, ui*3 + 9) += nu_*timefacfac*derxyz_(1, vi)*derxyz_(2, ui) ;
    estif_(vi*3 + 10, ui*3 + 10) += nu_*timefacfac*(derxyz_(0, ui)*derxyz_(0, vi) + derxyz_(1, ui)*derxyz_(1, vi) + 2.0*derxyz_(2, ui)*derxyz_(2, vi)) ;
#endif

#ifdef FLUID3_IS_TERM5
    /* Stabilisierung der Viskosität ( L_conv_u) */
    estif_(vi*3 + 8, ui*3 + 8) += 2.0*nu_*ttimetauMp*(conv_c_(ui)*viscs2_(0, 0, vi) + conv_g_(ui)*viscs2_(0, 0, vi) + viscs2_(0, 0, vi)*conv_r_(0, 0, ui) + viscs2_(0, 1, vi)*conv_r_(1, 0, ui) + viscs2_(0, 2, vi)*conv_r_(2, 0, ui)) ;
    estif_(vi*3 + 8, ui*3 + 9) += 2.0*nu_*ttimetauMp*(conv_c_(ui)*viscs2_(0, 1, vi) + conv_g_(ui)*viscs2_(0, 1, vi) + viscs2_(0, 0, vi)*conv_r_(0, 1, ui) + viscs2_(0, 1, vi)*conv_r_(1, 1, ui) + viscs2_(0, 2, vi)*conv_r_(2, 1, ui)) ;
    estif_(vi*3 + 8, ui*3 + 10) += 2.0*nu_*ttimetauMp*(conv_c_(ui)*viscs2_(0, 2, vi) + conv_g_(ui)*viscs2_(0, 2, vi) + viscs2_(0, 0, vi)*conv_r_(0, 2, ui) + viscs2_(0, 1, vi)*conv_r_(1, 2, ui) + viscs2_(0, 2, vi)*conv_r_(2, 2, ui)) ;
    estif_(vi*3 + 9, ui*3 + 8) += 2.0*nu_*ttimetauMp*(conv_c_(ui)*viscs2_(0, 1, vi) + conv_g_(ui)*viscs2_(0, 1, vi) + viscs2_(0, 1, vi)*conv_r_(0, 0, ui) + viscs2_(1, 1, vi)*conv_r_(1, 0, ui) + viscs2_(1, 2, vi)*conv_r_(2, 0, ui)) ;
    estif_(vi*3 + 9, ui*3 + 9) += 2.0*nu_*ttimetauMp*(conv_c_(ui)*viscs2_(1, 1, vi) + conv_g_(ui)*viscs2_(1, 1, vi) + viscs2_(0, 1, vi)*conv_r_(0, 1, ui) + viscs2_(1, 1, vi)*conv_r_(1, 1, ui) + viscs2_(1, 2, vi)*conv_r_(2, 1, ui)) ;
    estif_(vi*3 + 9, ui*3 + 10) += 2.0*nu_*ttimetauMp*(conv_c_(ui)*viscs2_(1, 2, vi) + conv_g_(ui)*viscs2_(1, 2, vi) + viscs2_(0, 1, vi)*conv_r_(0, 2, ui) + viscs2_(1, 1, vi)*conv_r_(1, 2, ui) + viscs2_(1, 2, vi)*conv_r_(2, 2, ui)) ;
    estif_(vi*3 + 10, ui*3 + 8) += 2.0*nu_*ttimetauMp*(conv_c_(ui)*viscs2_(0, 2, vi) + conv_g_(ui)*viscs2_(0, 2, vi) + viscs2_(0, 2, vi)*conv_r_(0, 0, ui) + viscs2_(1, 2, vi)*conv_r_(1, 0, ui) + viscs2_(2, 2, vi)*conv_r_(2, 0, ui)) ;
    estif_(vi*3 + 10, ui*3 + 9) += 2.0*nu_*ttimetauMp*(conv_c_(ui)*viscs2_(1, 2, vi) + conv_g_(ui)*viscs2_(1, 2, vi) + viscs2_(0, 2, vi)*conv_r_(0, 1, ui) + viscs2_(1, 2, vi)*conv_r_(1, 1, ui) + viscs2_(2, 2, vi)*conv_r_(2, 1, ui)) ;
    estif_(vi*3 + 10, ui*3 + 10) += 2.0*nu_*ttimetauMp*(conv_c_(ui)*viscs2_(2, 2, vi) + conv_g_(ui)*viscs2_(2, 2, vi) + viscs2_(0, 2, vi)*conv_r_(0, 2, ui) + viscs2_(1, 2, vi)*conv_r_(1, 2, ui) + viscs2_(2, 2, vi)*conv_r_(2, 2, ui)) ;
#endif

#ifdef FLUID3_IS_TERM6
    /* Stabilisierung der Viskosität (-L_visc_u) */
    estif_(vi*3 + 8, ui*3 + 8) += 4.0*(nu_*nu_)*ttimetauMp*(viscs2_(0, 0, ui)*viscs2_(0, 0, vi) + viscs2_(0, 1, ui)*viscs2_(0, 1, vi) + viscs2_(0, 2, ui)*viscs2_(0, 2, vi)) ;
    estif_(vi*3 + 8, ui*3 + 9) += 4.0*(nu_*nu_)*ttimetauMp*(viscs2_(0, 0, vi)*viscs2_(0, 1, ui) + viscs2_(0, 1, vi)*viscs2_(1, 1, ui) + viscs2_(0, 2, vi)*viscs2_(1, 2, ui)) ;
    estif_(vi*3 + 8, ui*3 + 10) += 4.0*(nu_*nu_)*ttimetauMp*(viscs2_(0, 0, vi)*viscs2_(0, 2, ui) + viscs2_(0, 1, vi)*viscs2_(1, 2, ui) + viscs2_(0, 2, vi)*viscs2_(2, 2, ui)) ;
    estif_(vi*3 + 9, ui*3 + 8) += 4.0*(nu_*nu_)*ttimetauMp*(viscs2_(0, 0, ui)*viscs2_(0, 1, vi) + viscs2_(0, 1, ui)*viscs2_(1, 1, vi) + viscs2_(0, 2, ui)*viscs2_(1, 2, vi)) ;
    estif_(vi*3 + 9, ui*3 + 9) += 4.0*(nu_*nu_)*ttimetauMp*(viscs2_(0, 1, ui)*viscs2_(0, 1, vi) + viscs2_(1, 1, ui)*viscs2_(1, 1, vi) + viscs2_(1, 2, ui)*viscs2_(1, 2, vi)) ;
    estif_(vi*3 + 9, ui*3 + 10) += 4.0*(nu_*nu_)*ttimetauMp*(viscs2_(0, 1, vi)*viscs2_(0, 2, ui) + viscs2_(1, 1, vi)*viscs2_(1, 2, ui) + viscs2_(1, 2, vi)*viscs2_(2, 2, ui)) ;
    estif_(vi*3 + 10, ui*3 + 8) += 4.0*(nu_*nu_)*ttimetauMp*(viscs2_(0, 0, ui)*viscs2_(0, 2, vi) + viscs2_(0, 1, ui)*viscs2_(1, 2, vi) + viscs2_(0, 2, ui)*viscs2_(2, 2, vi)) ;
    estif_(vi*3 + 10, ui*3 + 9) += 4.0*(nu_*nu_)*ttimetauMp*(viscs2_(0, 1, ui)*viscs2_(0, 2, vi) + viscs2_(1, 1, ui)*viscs2_(1, 2, vi) + viscs2_(1, 2, ui)*viscs2_(2, 2, vi)) ;
    estif_(vi*3 + 10, ui*3 + 10) += 4.0*(nu_*nu_)*ttimetauMp*(viscs2_(0, 2, ui)*viscs2_(0, 2, vi) + viscs2_(1, 2, ui)*viscs2_(1, 2, vi) + viscs2_(2, 2, ui)*viscs2_(2, 2, vi)) ;
#endif

#ifdef FLUID3_IS_TERM7
    /* Druckterm */
#endif

#ifdef FLUID3_IS_TERM8
    /* Divergenzfreiheit */
#endif

#ifdef FLUID3_IS_TERM9
    /* Kontinuitätsstabilisierung */
    estif_(vi*3 + 8, ui*3 + 8) += (thsl*thsl)*tau_C*derxyz_(0, ui)*derxyz_(0, vi) ;
    estif_(vi*3 + 8, ui*3 + 9) += (thsl*thsl)*tau_C*derxyz_(0, vi)*derxyz_(1, ui) ;
    estif_(vi*3 + 8, ui*3 + 10) += (thsl*thsl)*tau_C*derxyz_(0, vi)*derxyz_(2, ui) ;
    estif_(vi*3 + 9, ui*3 + 8) += (thsl*thsl)*tau_C*derxyz_(0, ui)*derxyz_(1, vi) ;
    estif_(vi*3 + 9, ui*3 + 9) += (thsl*thsl)*tau_C*derxyz_(1, ui)*derxyz_(1, vi) ;
    estif_(vi*3 + 9, ui*3 + 10) += (thsl*thsl)*tau_C*derxyz_(1, vi)*derxyz_(2, ui) ;
    estif_(vi*3 + 10, ui*3 + 8) += (thsl*thsl)*tau_C*derxyz_(0, ui)*derxyz_(2, vi) ;
    estif_(vi*3 + 10, ui*3 + 9) += (thsl*thsl)*tau_C*derxyz_(1, ui)*derxyz_(2, vi) ;
    estif_(vi*3 + 10, ui*3 + 10) += (thsl*thsl)*tau_C*derxyz_(2, ui)*derxyz_(2, vi) ;
#endif

#ifdef FLUID3_IS_TERM10
    /* Massenterm */
    estif_(vi*3 + 8, ui*3 + 8) += fac*funct_(ui)*funct_(vi) ;
    estif_(vi*3 + 9, ui*3 + 9) += fac*funct_(ui)*funct_(vi) ;
    estif_(vi*3 + 10, ui*3 + 10) += fac*funct_(ui)*funct_(vi) ;
#endif

#ifdef FLUID3_IS_TERM11
    /* Konvektionsstabilisierung */
    estif_(vi*3 + 8, ui*3 + 8) += timetauM*funct_(ui)*(conv_g_(vi) + 2.0*velint_(0)*derxyz_(0, vi) + velint_(1)*derxyz_(1, vi) + velint_(2)*derxyz_(2, vi)) ;
    estif_(vi*3 + 8, ui*3 + 9) += timetauM*funct_(ui)*velint_(0)*derxyz_(1, vi) ;
    estif_(vi*3 + 8, ui*3 + 10) += timetauM*funct_(ui)*velint_(0)*derxyz_(2, vi) ;
    estif_(vi*3 + 9, ui*3 + 8) += timetauM*funct_(ui)*velint_(1)*derxyz_(0, vi) ;
    estif_(vi*3 + 9, ui*3 + 9) += timetauM*funct_(ui)*(conv_g_(vi) + velint_(0)*derxyz_(0, vi) + 2.0*velint_(1)*derxyz_(1, vi) + velint_(2)*derxyz_(2, vi)) ;
    estif_(vi*3 + 9, ui*3 + 10) += timetauM*funct_(ui)*velint_(1)*derxyz_(2, vi) ;
    estif_(vi*3 + 10, ui*3 + 8) += timetauM*funct_(ui)*velint_(2)*derxyz_(0, vi) ;
    estif_(vi*3 + 10, ui*3 + 9) += timetauM*funct_(ui)*velint_(2)*derxyz_(1, vi) ;
    estif_(vi*3 + 10, ui*3 + 10) += timetauM*funct_(ui)*(conv_g_(vi) + velint_(0)*derxyz_(0, vi) + velint_(1)*derxyz_(1, vi) + 2.0*velint_(2)*derxyz_(2, vi)) ;
#endif

#ifdef FLUID3_IS_TERM12
    /* Viskositätsstabilisierung */
    estif_(vi*3 + 8, ui*3 + 8) += 2.0*nu_*timetauMp*funct_(ui)*viscs2_(0, 0, vi) ;
    estif_(vi*3 + 8, ui*3 + 9) += 2.0*nu_*timetauMp*funct_(ui)*viscs2_(0, 1, vi) ;
    estif_(vi*3 + 8, ui*3 + 10) += 2.0*nu_*timetauMp*funct_(ui)*viscs2_(0, 2, vi) ;
    estif_(vi*3 + 9, ui*3 + 8) += 2.0*nu_*timetauMp*funct_(ui)*viscs2_(0, 1, vi) ;
    estif_(vi*3 + 9, ui*3 + 9) += 2.0*nu_*timetauMp*funct_(ui)*viscs2_(1, 1, vi) ;
    estif_(vi*3 + 9, ui*3 + 10) += 2.0*nu_*timetauMp*funct_(ui)*viscs2_(1, 2, vi) ;
    estif_(vi*3 + 10, ui*3 + 8) += 2.0*nu_*timetauMp*funct_(ui)*viscs2_(0, 2, vi) ;
    estif_(vi*3 + 10, ui*3 + 9) += 2.0*nu_*timetauMp*funct_(ui)*viscs2_(1, 2, vi) ;
    estif_(vi*3 + 10, ui*3 + 10) += 2.0*nu_*timetauMp*funct_(ui)*viscs2_(2, 2, vi) ;
#endif

#ifdef FLUID3_IS_TERM13
    /* Quellterm der rechten Seite */
#endif

#ifdef FLUID3_IS_TERM14
    /* Konvektionsstabilisierung */
    estif_(vi*3 + 8, ui*3 + 8) += -(timetauM*funct_(ui)*rhsint_(0)*derxyz_(0, vi)) ;
    estif_(vi*3 + 8, ui*3 + 9) += -(timetauM*funct_(ui)*rhsint_(0)*derxyz_(1, vi)) ;
    estif_(vi*3 + 8, ui*3 + 10) += -(timetauM*funct_(ui)*rhsint_(0)*derxyz_(2, vi)) ;
    estif_(vi*3 + 9, ui*3 + 8) += -(timetauM*funct_(ui)*rhsint_(1)*derxyz_(0, vi)) ;
    estif_(vi*3 + 9, ui*3 + 9) += -(timetauM*funct_(ui)*rhsint_(1)*derxyz_(1, vi)) ;
    estif_(vi*3 + 9, ui*3 + 10) += -(timetauM*funct_(ui)*rhsint_(1)*derxyz_(2, vi)) ;
    estif_(vi*3 + 10, ui*3 + 8) += -(timetauM*funct_(ui)*rhsint_(2)*derxyz_(0, vi)) ;
    estif_(vi*3 + 10, ui*3 + 9) += -(timetauM*funct_(ui)*rhsint_(2)*derxyz_(1, vi)) ;
    estif_(vi*3 + 10, ui*3 + 10) += -(timetauM*funct_(ui)*rhsint_(2)*derxyz_(2, vi)) ;
#endif

#ifdef FLUID3_IS_TERM15
    /* Viskositätsstabilisierung */
#endif

  }
}
