
for (vi=0; vi<iel; ++vi)
{
  for (ui=0; ui<iel; ++ui)
  {
    /* Konvektionsterm (u*grad(u),v) */
    estif_(vi*2, ui*2)         += timefacfac*funct_(vi)*(conv_c_(ui) + conv_r_(0, 0, ui)) ;
    estif_(vi*2, ui*2 + 1)     += timefacfac*funct_(vi)*conv_r_(0, 1, ui) ;
    estif_(vi*2 + 1, ui*2)     += timefacfac*funct_(vi)*conv_r_(1, 0, ui) ;
    estif_(vi*2 + 1, ui*2 + 1) += timefacfac*funct_(vi)*(conv_c_(ui) + conv_r_(1, 1, ui)) ;

    /* Stabilisierung der Konvektion (u*grad(u),u*grad(v)) */
    estif_(vi*2, ui*2)         += ttimetauM*(conv_c_(ui)*conv_c_(vi) + conv_c_(vi)*conv_r_(0, 0, ui) + derxy_(0, vi)*vconv_r_(0, ui)) ;
    estif_(vi*2, ui*2 + 1)     += ttimetauM*(conv_c_(vi)*conv_r_(0, 1, ui) + derxy_(1, vi)*vconv_r_(0, ui)) ;
    estif_(vi*2 + 1, ui*2)     += ttimetauM*(conv_c_(vi)*conv_r_(1, 0, ui) + derxy_(0, vi)*vconv_r_(1, ui)) ;
    estif_(vi*2 + 1, ui*2 + 1) += ttimetauM*(conv_c_(ui)*conv_c_(vi) + conv_c_(vi)*conv_r_(1, 1, ui) + derxy_(1, vi)*vconv_r_(1, ui)) ;

    /* Stabilisierung der Konvektion (-2*nu*div(epsilon((u))),u*grad(v)) */
    estif_(vi*2, ui*2)         += -2.0*visc_*ttimetauM*(-(conv_c_(vi)*viscs2_(0, 0, ui)) + funct_(ui)*derxy_(0, vi)*visc_old_(0)) ;
    estif_(vi*2, ui*2 + 1)     += -2.0*visc_*ttimetauM*(-(conv_c_(vi)*viscs2_(0, 1, ui)) + funct_(ui)*derxy_(1, vi)*visc_old_(0)) ;
    estif_(vi*2 + 1, ui*2)     += -2.0*visc_*ttimetauM*(-(conv_c_(vi)*viscs2_(0, 1, ui)) + funct_(ui)*derxy_(0, vi)*visc_old_(1)) ;
    estif_(vi*2 + 1, ui*2 + 1) += -2.0*visc_*ttimetauM*(-(conv_c_(vi)*viscs2_(1, 1, ui)) + funct_(ui)*derxy_(1, vi)*visc_old_(1)) ;

    /* Stabilisierung der Konvektion (grad(p),u*grad(v)) */
    estif_(vi*2, ui*2)         += ttimetauM*funct_(ui)*gradp_(0)*derxy_(0, vi) ;
    estif_(vi*2, ui*2 + 1)     += ttimetauM*funct_(ui)*gradp_(0)*derxy_(1, vi) ;
    estif_(vi*2 + 1, ui*2)     += ttimetauM*funct_(ui)*gradp_(1)*derxy_(0, vi) ;
    estif_(vi*2 + 1, ui*2 + 1) += ttimetauM*funct_(ui)*gradp_(1)*derxy_(1, vi) ;

    /* Viskositätsterm (2*nu*epsilon(u),epsilon(v)) */
    estif_(vi*2, ui*2)         += fac*time2nue*derxy_(0, ui)*derxy_(0, vi) + visc_*timefacfac*derxy_(1, ui)*derxy_(1, vi) ;
    estif_(vi*2, ui*2 + 1)     += visc_*timefacfac*derxy_(0, ui)*derxy_(1, vi) ;
    estif_(vi*2 + 1, ui*2)     += visc_*timefacfac*derxy_(0, vi)*derxy_(1, ui) ;
    estif_(vi*2 + 1, ui*2 + 1) += fac*time2nue*derxy_(1, ui)*derxy_(1, vi) + visc_*timefacfac*derxy_(0, ui)*derxy_(0, vi) ;

    /* Stabilisierung der Viskosität (u*grad(u),-2*nu*div(epsilon((v)))) */
    estif_(vi*2, ui*2)         += 2.0*visc_*ttimetauM*(conv_c_(ui)*viscs2_(0, 0, vi) + viscs2_(0, 0, vi)*conv_r_(0, 0, ui) + viscs2_(0, 1, vi)*conv_r_(1, 0, ui)) ;
    estif_(vi*2, ui*2 + 1)     += 2.0*visc_*ttimetauM*(conv_c_(ui)*viscs2_(0, 1, vi) + viscs2_(0, 0, vi)*conv_r_(0, 1, ui) + viscs2_(0, 1, vi)*conv_r_(1, 1, ui)) ;
    estif_(vi*2 + 1, ui*2)     += 2.0*visc_*ttimetauM*(conv_c_(ui)*viscs2_(0, 1, vi) + viscs2_(0, 1, vi)*conv_r_(0, 0, ui) + viscs2_(1, 1, vi)*conv_r_(1, 0, ui)) ;
    estif_(vi*2 + 1, ui*2 + 1) += 2.0*visc_*ttimetauM*(conv_c_(ui)*viscs2_(1, 1, vi) + viscs2_(0, 1, vi)*conv_r_(0, 1, ui) + viscs2_(1, 1, vi)*conv_r_(1, 1, ui)) ;

    /* Stabilisierung der Viskosität (-2*nu*div(epsilon((u))),-2*nu*div(epsilon((v)))) */
    estif_(vi*2, ui*2)         += 4.0*(visc_*visc_)*ttimetauM*(viscs2_(0, 0, ui)*viscs2_(0, 0, vi) + viscs2_(0, 1, ui)*viscs2_(0, 1, vi)) ;
    estif_(vi*2, ui*2 + 1)     += 4.0*(visc_*visc_)*ttimetauM*(viscs2_(0, 0, vi)*viscs2_(0, 1, ui) + viscs2_(0, 1, vi)*viscs2_(1, 1, ui)) ;
    estif_(vi*2 + 1, ui*2)     += 4.0*(visc_*visc_)*ttimetauM*(viscs2_(0, 0, ui)*viscs2_(0, 1, vi) + viscs2_(0, 1, ui)*viscs2_(1, 1, vi)) ;
    estif_(vi*2 + 1, ui*2 + 1) += 4.0*(visc_*visc_)*ttimetauM*(viscs2_(0, 1, ui)*viscs2_(0, 1, vi) + viscs2_(1, 1, ui)*viscs2_(1, 1, vi)) ;

    /* Stabilisierung der Viskosität (grad(p),-2*nu*div(epsilon((v)))) */

    /* Massenterm (u,v) */
    estif_(vi*2, ui*2)         += fac*funct_(ui)*funct_(vi) ;
    estif_(vi*2 + 1, ui*2 + 1) += fac*funct_(ui)*funct_(vi) ;

    /* Konvektionsstabilisierung (u,u*grad(v)) */
    estif_(vi*2, ui*2)         += timetauM*funct_(ui)*(conv_c_(vi) + velint_(0)*derxy_(0, vi)) ;
    estif_(vi*2, ui*2 + 1)     += timetauM*funct_(ui)*velint_(0)*derxy_(1, vi) ;
    estif_(vi*2 + 1, ui*2)     += timetauM*funct_(ui)*velint_(1)*derxy_(0, vi) ;
    estif_(vi*2 + 1, ui*2 + 1) += timetauM*funct_(ui)*(conv_c_(vi) + velint_(1)*derxy_(1, vi)) ;

    /* Viskositätsstabilisierung (u,-2*nu*div(epsilon((v)))) */
    estif_(vi*2, ui*2)         += tau_M*time2nue*funct_(ui)*viscs2_(0, 0, vi) ;
    estif_(vi*2, ui*2 + 1)     += tau_M*time2nue*funct_(ui)*viscs2_(0, 1, vi) ;
    estif_(vi*2 + 1, ui*2)     += tau_M*time2nue*funct_(ui)*viscs2_(0, 1, vi) ;
    estif_(vi*2 + 1, ui*2 + 1) += tau_M*time2nue*funct_(ui)*viscs2_(1, 1, vi) ;

    /* Quellterm der rechten Seite (b,v) */

    /* Konvektionsstabilisierung (b,u*grad(v)) */
    estif_(vi*2, ui*2)         += -(timetauM*funct_(ui)*rhsint_(0)*derxy_(0, vi)) ;
    estif_(vi*2, ui*2 + 1)     += -(timetauM*funct_(ui)*rhsint_(0)*derxy_(1, vi)) ;
    estif_(vi*2 + 1, ui*2)     += -(timetauM*funct_(ui)*rhsint_(1)*derxy_(0, vi)) ;
    estif_(vi*2 + 1, ui*2 + 1) += -(timetauM*funct_(ui)*rhsint_(1)*derxy_(1, vi)) ;

    /* Viskositätsstabilisierung (b,-2*nu*div(epsilon((v))) */

  }
}
