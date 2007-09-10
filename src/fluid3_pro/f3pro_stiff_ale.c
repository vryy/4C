#ifndef CCADISCRET
for (vi=0; vi<iel; ++vi)
{
  for (ui=0; ui<iel; ++ui)
  {
    /* Konvektionsterm (u*grad(u),v) */
    estif_(vi*3, ui*3)         += timefacfac*funct_(vi)*(conv_c_(ui) + conv_g_(ui) + conv_r_(0, 0, ui)) ;
    estif_(vi*3, ui*3 + 1)     += timefacfac*funct_(vi)*conv_r_(0, 1, ui) ;
    estif_(vi*3, ui*3 + 2)     += timefacfac*funct_(vi)*conv_r_(0, 2, ui) ;
    estif_(vi*3 + 1, ui*3)     += timefacfac*funct_(vi)*conv_r_(1, 0, ui) ;
    estif_(vi*3 + 1, ui*3 + 1) += timefacfac*funct_(vi)*(conv_c_(ui) + conv_g_(ui) + conv_r_(1, 1, ui)) ;
    estif_(vi*3 + 1, ui*3 + 2) += timefacfac*funct_(vi)*conv_r_(1, 2, ui) ;
    estif_(vi*3 + 2, ui*3)     += timefacfac*funct_(vi)*conv_r_(2, 0, ui) ;
    estif_(vi*3 + 2, ui*3 + 1) += timefacfac*funct_(vi)*conv_r_(2, 1, ui) ;
    estif_(vi*3 + 2, ui*3 + 2) += timefacfac*funct_(vi)*(conv_c_(ui) + conv_g_(ui) + conv_r_(2, 2, ui)) ;

    /* Stabilisierung der Konvektion (u*grad(u),u*grad(v)) */
    estif_(vi*3, ui*3)         += ttimetauM*(conv_c_(ui)*conv_c_(vi) - conv_g_(ui)*gridvint_(0)*derxyz_(0, vi) + 2.0*velint_(0)*derxyz_(0, vi)*conv_r_(0, 0, ui) + velint_(1)*derxyz_(0, vi)*conv_r_(0, 1, ui) + velint_(1)*derxyz_(1, vi)*conv_r_(0, 0, ui) + velint_(2)*derxyz_(0, vi)*conv_r_(0, 2, ui) + velint_(2)*derxyz_(2, vi)*conv_r_(0, 0, ui) - 2.0*gridvint_(0)*derxyz_(0, vi)*conv_r_(0, 0, ui) - gridvint_(1)*derxyz_(0, vi)*conv_r_(0, 1, ui) - gridvint_(1)*derxyz_(1, vi)*conv_r_(0, 0, ui) - gridvint_(2)*derxyz_(0, vi)*conv_r_(0, 2, ui) - gridvint_(2)*derxyz_(2, vi)*conv_r_(0, 0, ui) + (gridvint_(1)*gridvint_(1))*derxyz_(1, ui)*derxyz_(1, vi) + (gridvint_(2)*gridvint_(2))*derxyz_(2, ui)*derxyz_(2, vi) - 2.0*velint_(0)*gridvint_(0)*derxyz_(0, ui)*derxyz_(0, vi) - velint_(0)*gridvint_(1)*derxyz_(0, ui)*derxyz_(1, vi) - velint_(0)*gridvint_(1)*derxyz_(0, vi)*derxyz_(1, ui) - velint_(1)*gridvint_(0)*derxyz_(0, ui)*derxyz_(1, vi) - velint_(1)*gridvint_(0)*derxyz_(0, vi)*derxyz_(1, ui) - velint_(0)*gridvint_(2)*derxyz_(0, ui)*derxyz_(2, vi) - velint_(0)*gridvint_(2)*derxyz_(0, vi)*derxyz_(2, ui) - velint_(2)*gridvint_(0)*derxyz_(0, ui)*derxyz_(2, vi) - velint_(2)*gridvint_(0)*derxyz_(0, vi)*derxyz_(2, ui) - 2.0*velint_(1)*gridvint_(1)*derxyz_(1, ui)*derxyz_(1, vi) - velint_(1)*gridvint_(2)*derxyz_(1, ui)*derxyz_(2, vi) - velint_(1)*gridvint_(2)*derxyz_(1, vi)*derxyz_(2, ui) - velint_(2)*gridvint_(1)*derxyz_(1, ui)*derxyz_(2, vi) - velint_(2)*gridvint_(1)*derxyz_(1, vi)*derxyz_(2, ui) - 2.0*velint_(2)*gridvint_(2)*derxyz_(2, ui)*derxyz_(2, vi) + gridvint_(0)*gridvint_(1)*derxyz_(0, ui)*derxyz_(1, vi) + gridvint_(0)*gridvint_(2)*derxyz_(0, ui)*derxyz_(2, vi) + gridvint_(1)*gridvint_(2)*derxyz_(1, ui)*derxyz_(2, vi) + gridvint_(1)*gridvint_(2)*derxyz_(1, vi)*derxyz_(2, ui)) ;
    estif_(vi*3, ui*3 + 1)     += ttimetauM*(conv_c_(vi)*conv_r_(0, 1, ui) + derxyz_(1, vi)*vconv_r_(0, ui) - gridvint_(0)*derxyz_(0, vi)*conv_r_(0, 1, ui) - gridvint_(0)*derxyz_(1, vi)*conv_r_(0, 0, ui) - 2.0*gridvint_(1)*derxyz_(1, vi)*conv_r_(0, 1, ui) - gridvint_(2)*derxyz_(1, vi)*conv_r_(0, 2, ui) - gridvint_(2)*derxyz_(2, vi)*conv_r_(0, 1, ui)) ;
    estif_(vi*3, ui*3 + 2)     += ttimetauM*(conv_c_(vi)*conv_r_(0, 2, ui) + derxyz_(2, vi)*vconv_r_(0, ui) - gridvint_(0)*derxyz_(0, vi)*conv_r_(0, 2, ui) - gridvint_(0)*derxyz_(2, vi)*conv_r_(0, 0, ui) - gridvint_(1)*derxyz_(1, vi)*conv_r_(0, 2, ui) - gridvint_(1)*derxyz_(2, vi)*conv_r_(0, 1, ui) - 2.0*gridvint_(2)*derxyz_(2, vi)*conv_r_(0, 2, ui)) ;
    estif_(vi*3 + 1, ui*3)     += ttimetauM*(conv_c_(vi)*conv_r_(1, 0, ui) + derxyz_(0, vi)*vconv_r_(1, ui) - 2.0*gridvint_(0)*derxyz_(0, vi)*conv_r_(1, 0, ui) - gridvint_(1)*derxyz_(0, vi)*conv_r_(1, 1, ui) - gridvint_(1)*derxyz_(1, vi)*conv_r_(1, 0, ui) - gridvint_(2)*derxyz_(0, vi)*conv_r_(1, 2, ui) - gridvint_(2)*derxyz_(2, vi)*conv_r_(1, 0, ui)) ;
    estif_(vi*3 + 1, ui*3 + 1) += ttimetauM*(conv_c_(ui)*conv_c_(vi) - conv_g_(ui)*gridvint_(0)*derxyz_(0, vi) + velint_(0)*derxyz_(0, vi)*conv_r_(1, 1, ui) + velint_(0)*derxyz_(1, vi)*conv_r_(1, 0, ui) + 2.0*velint_(1)*derxyz_(1, vi)*conv_r_(1, 1, ui) + velint_(2)*derxyz_(1, vi)*conv_r_(1, 2, ui) + velint_(2)*derxyz_(2, vi)*conv_r_(1, 1, ui) - gridvint_(0)*derxyz_(0, vi)*conv_r_(1, 1, ui) - gridvint_(0)*derxyz_(1, vi)*conv_r_(1, 0, ui) - 2.0*gridvint_(1)*derxyz_(1, vi)*conv_r_(1, 1, ui) - gridvint_(2)*derxyz_(1, vi)*conv_r_(1, 2, ui) - gridvint_(2)*derxyz_(2, vi)*conv_r_(1, 1, ui) + (gridvint_(1)*gridvint_(1))*derxyz_(1, ui)*derxyz_(1, vi) + (gridvint_(2)*gridvint_(2))*derxyz_(2, ui)*derxyz_(2, vi) - 2.0*velint_(0)*gridvint_(0)*derxyz_(0, ui)*derxyz_(0, vi) - velint_(0)*gridvint_(1)*derxyz_(0, ui)*derxyz_(1, vi) - velint_(0)*gridvint_(1)*derxyz_(0, vi)*derxyz_(1, ui) - velint_(1)*gridvint_(0)*derxyz_(0, ui)*derxyz_(1, vi) - velint_(1)*gridvint_(0)*derxyz_(0, vi)*derxyz_(1, ui) - velint_(0)*gridvint_(2)*derxyz_(0, ui)*derxyz_(2, vi) - velint_(0)*gridvint_(2)*derxyz_(0, vi)*derxyz_(2, ui) - velint_(2)*gridvint_(0)*derxyz_(0, ui)*derxyz_(2, vi) - velint_(2)*gridvint_(0)*derxyz_(0, vi)*derxyz_(2, ui) - 2.0*velint_(1)*gridvint_(1)*derxyz_(1, ui)*derxyz_(1, vi) - velint_(1)*gridvint_(2)*derxyz_(1, ui)*derxyz_(2, vi) - velint_(1)*gridvint_(2)*derxyz_(1, vi)*derxyz_(2, ui) - velint_(2)*gridvint_(1)*derxyz_(1, ui)*derxyz_(2, vi) - velint_(2)*gridvint_(1)*derxyz_(1, vi)*derxyz_(2, ui) - 2.0*velint_(2)*gridvint_(2)*derxyz_(2, ui)*derxyz_(2, vi) + gridvint_(0)*gridvint_(1)*derxyz_(0, ui)*derxyz_(1, vi) + gridvint_(0)*gridvint_(2)*derxyz_(0, ui)*derxyz_(2, vi) + gridvint_(1)*gridvint_(2)*derxyz_(1, ui)*derxyz_(2, vi) + gridvint_(1)*gridvint_(2)*derxyz_(1, vi)*derxyz_(2, ui)) ;
    estif_(vi*3 + 1, ui*3 + 2) += ttimetauM*(conv_c_(vi)*conv_r_(1, 2, ui) + derxyz_(2, vi)*vconv_r_(1, ui) - gridvint_(0)*derxyz_(0, vi)*conv_r_(1, 2, ui) - gridvint_(0)*derxyz_(2, vi)*conv_r_(1, 0, ui) - gridvint_(1)*derxyz_(1, vi)*conv_r_(1, 2, ui) - gridvint_(1)*derxyz_(2, vi)*conv_r_(1, 1, ui) - 2.0*gridvint_(2)*derxyz_(2, vi)*conv_r_(1, 2, ui)) ;
    estif_(vi*3 + 2, ui*3)     += ttimetauM*(conv_c_(vi)*conv_r_(2, 0, ui) + derxyz_(0, vi)*vconv_r_(2, ui) - 2.0*gridvint_(0)*derxyz_(0, vi)*conv_r_(2, 0, ui) - gridvint_(1)*derxyz_(0, vi)*conv_r_(2, 1, ui) - gridvint_(1)*derxyz_(1, vi)*conv_r_(2, 0, ui) - gridvint_(2)*derxyz_(0, vi)*conv_r_(2, 2, ui) - gridvint_(2)*derxyz_(2, vi)*conv_r_(2, 0, ui)) ;
    estif_(vi*3 + 2, ui*3 + 1) += ttimetauM*(conv_c_(vi)*conv_r_(2, 1, ui) + derxyz_(1, vi)*vconv_r_(2, ui) - gridvint_(0)*derxyz_(0, vi)*conv_r_(2, 1, ui) - gridvint_(0)*derxyz_(1, vi)*conv_r_(2, 0, ui) - 2.0*gridvint_(1)*derxyz_(1, vi)*conv_r_(2, 1, ui) - gridvint_(2)*derxyz_(1, vi)*conv_r_(2, 2, ui) - gridvint_(2)*derxyz_(2, vi)*conv_r_(2, 1, ui)) ;
    estif_(vi*3 + 2, ui*3 + 2) += ttimetauM*(conv_c_(ui)*conv_c_(vi) - conv_g_(ui)*gridvint_(0)*derxyz_(0, vi) + velint_(0)*derxyz_(0, vi)*conv_r_(2, 2, ui) + velint_(0)*derxyz_(2, vi)*conv_r_(2, 0, ui) + velint_(1)*derxyz_(1, vi)*conv_r_(2, 2, ui) + velint_(1)*derxyz_(2, vi)*conv_r_(2, 1, ui) + 2.0*velint_(2)*derxyz_(2, vi)*conv_r_(2, 2, ui) - gridvint_(0)*derxyz_(0, vi)*conv_r_(2, 2, ui) - gridvint_(0)*derxyz_(2, vi)*conv_r_(2, 0, ui) - gridvint_(1)*derxyz_(1, vi)*conv_r_(2, 2, ui) - gridvint_(1)*derxyz_(2, vi)*conv_r_(2, 1, ui) - 2.0*gridvint_(2)*derxyz_(2, vi)*conv_r_(2, 2, ui) + (gridvint_(1)*gridvint_(1))*derxyz_(1, ui)*derxyz_(1, vi) + (gridvint_(2)*gridvint_(2))*derxyz_(2, ui)*derxyz_(2, vi) - 2.0*velint_(0)*gridvint_(0)*derxyz_(0, ui)*derxyz_(0, vi) - velint_(0)*gridvint_(1)*derxyz_(0, ui)*derxyz_(1, vi) - velint_(0)*gridvint_(1)*derxyz_(0, vi)*derxyz_(1, ui) - velint_(1)*gridvint_(0)*derxyz_(0, ui)*derxyz_(1, vi) - velint_(1)*gridvint_(0)*derxyz_(0, vi)*derxyz_(1, ui) - velint_(0)*gridvint_(2)*derxyz_(0, ui)*derxyz_(2, vi) - velint_(0)*gridvint_(2)*derxyz_(0, vi)*derxyz_(2, ui) - velint_(2)*gridvint_(0)*derxyz_(0, ui)*derxyz_(2, vi) - velint_(2)*gridvint_(0)*derxyz_(0, vi)*derxyz_(2, ui) - 2.0*velint_(1)*gridvint_(1)*derxyz_(1, ui)*derxyz_(1, vi) - velint_(1)*gridvint_(2)*derxyz_(1, ui)*derxyz_(2, vi) - velint_(1)*gridvint_(2)*derxyz_(1, vi)*derxyz_(2, ui) - velint_(2)*gridvint_(1)*derxyz_(1, ui)*derxyz_(2, vi) - velint_(2)*gridvint_(1)*derxyz_(1, vi)*derxyz_(2, ui) - 2.0*velint_(2)*gridvint_(2)*derxyz_(2, ui)*derxyz_(2, vi) + gridvint_(0)*gridvint_(1)*derxyz_(0, ui)*derxyz_(1, vi) + gridvint_(0)*gridvint_(2)*derxyz_(0, ui)*derxyz_(2, vi) + gridvint_(1)*gridvint_(2)*derxyz_(1, ui)*derxyz_(2, vi) + gridvint_(1)*gridvint_(2)*derxyz_(1, vi)*derxyz_(2, ui)) ;

    /* Stabilisierung der Konvektion (-2*nu*div(epsilon((u))),u*grad(v)) */
    estif_(vi*3, ui*3)         += -2.0*visc_*ttimetauM*(-(conv_c_(vi)*viscs2_(0, 0, ui)) - conv_g_(vi)*viscs2_(0, 0, ui) + funct_(ui)*visc_old_(0)*derxyz_(0, vi)) ;
    estif_(vi*3, ui*3 + 1)     += -2.0*visc_*ttimetauM*(-(conv_c_(vi)*viscs2_(0, 1, ui)) - conv_g_(vi)*viscs2_(0, 1, ui) + funct_(ui)*visc_old_(0)*derxyz_(1, vi)) ;
    estif_(vi*3, ui*3 + 2)     += -2.0*visc_*ttimetauM*(-(conv_c_(vi)*viscs2_(0, 2, ui)) - conv_g_(vi)*viscs2_(0, 2, ui) + funct_(ui)*visc_old_(0)*derxyz_(2, vi)) ;
    estif_(vi*3 + 1, ui*3)     += -2.0*visc_*ttimetauM*(-(conv_c_(vi)*viscs2_(0, 1, ui)) - conv_g_(vi)*viscs2_(0, 1, ui) + funct_(ui)*visc_old_(1)*derxyz_(0, vi)) ;
    estif_(vi*3 + 1, ui*3 + 1) += -2.0*visc_*ttimetauM*(-(conv_c_(vi)*viscs2_(1, 1, ui)) - conv_g_(vi)*viscs2_(1, 1, ui) + funct_(ui)*visc_old_(1)*derxyz_(1, vi)) ;
    estif_(vi*3 + 1, ui*3 + 2) += -2.0*visc_*ttimetauM*(-(conv_c_(vi)*viscs2_(1, 2, ui)) - conv_g_(vi)*viscs2_(1, 2, ui) + funct_(ui)*visc_old_(1)*derxyz_(2, vi)) ;
    estif_(vi*3 + 2, ui*3)     += -2.0*visc_*ttimetauM*(-(conv_c_(vi)*viscs2_(0, 2, ui)) - conv_g_(vi)*viscs2_(0, 2, ui) + funct_(ui)*visc_old_(2)*derxyz_(0, vi)) ;
    estif_(vi*3 + 2, ui*3 + 1) += -2.0*visc_*ttimetauM*(-(conv_c_(vi)*viscs2_(1, 2, ui)) - conv_g_(vi)*viscs2_(1, 2, ui) + funct_(ui)*visc_old_(2)*derxyz_(1, vi)) ;
    estif_(vi*3 + 2, ui*3 + 2) += -2.0*visc_*ttimetauM*(-(conv_c_(vi)*viscs2_(2, 2, ui)) - conv_g_(vi)*viscs2_(2, 2, ui) + funct_(ui)*visc_old_(2)*derxyz_(2, vi)) ;

    /* Stabilisierung der Konvektion (grad(p),u*grad(v)) */
    estif_(vi*3, ui*3)         += ttimetauM*funct_(ui)*gradp_(0)*derxyz_(0, vi) ;
    estif_(vi*3, ui*3 + 1)     += ttimetauM*funct_(ui)*gradp_(0)*derxyz_(1, vi) ;
    estif_(vi*3, ui*3 + 2)     += ttimetauM*funct_(ui)*gradp_(0)*derxyz_(2, vi) ;
    estif_(vi*3 + 1, ui*3)     += ttimetauM*funct_(ui)*gradp_(1)*derxyz_(0, vi) ;
    estif_(vi*3 + 1, ui*3 + 1) += ttimetauM*funct_(ui)*gradp_(1)*derxyz_(1, vi) ;
    estif_(vi*3 + 1, ui*3 + 2) += ttimetauM*funct_(ui)*gradp_(1)*derxyz_(2, vi) ;
    estif_(vi*3 + 2, ui*3)     += ttimetauM*funct_(ui)*gradp_(2)*derxyz_(0, vi) ;
    estif_(vi*3 + 2, ui*3 + 1) += ttimetauM*funct_(ui)*gradp_(2)*derxyz_(1, vi) ;
    estif_(vi*3 + 2, ui*3 + 2) += ttimetauM*funct_(ui)*gradp_(2)*derxyz_(2, vi) ;

    /* Viskositätsterm (2*nu*epsilon(u),epsilon(v)) */
    estif_(vi*3, ui*3)         += fac*time2nue*derxyz_(0, ui)*derxyz_(0, vi) + visc_*timefacfac*derxyz_(1, ui)*derxyz_(1, vi) + visc_*timefacfac*derxyz_(2, ui)*derxyz_(2, vi) ;
    estif_(vi*3, ui*3 + 1)     += visc_*timefacfac*derxyz_(0, ui)*derxyz_(1, vi) ;
    estif_(vi*3, ui*3 + 2)     += visc_*timefacfac*derxyz_(0, ui)*derxyz_(2, vi) ;
    estif_(vi*3 + 1, ui*3)     += visc_*timefacfac*derxyz_(0, vi)*derxyz_(1, ui) ;
    estif_(vi*3 + 1, ui*3 + 1) += fac*time2nue*derxyz_(1, ui)*derxyz_(1, vi) + visc_*timefacfac*derxyz_(0, ui)*derxyz_(0, vi) + visc_*timefacfac*derxyz_(2, ui)*derxyz_(2, vi) ;
    estif_(vi*3 + 1, ui*3 + 2) += visc_*timefacfac*derxyz_(1, ui)*derxyz_(2, vi) ;
    estif_(vi*3 + 2, ui*3)     += visc_*timefacfac*derxyz_(0, vi)*derxyz_(2, ui) ;
    estif_(vi*3 + 2, ui*3 + 1) += visc_*timefacfac*derxyz_(1, vi)*derxyz_(2, ui) ;
    estif_(vi*3 + 2, ui*3 + 2) += fac*time2nue*derxyz_(2, ui)*derxyz_(2, vi) + visc_*timefacfac*derxyz_(0, ui)*derxyz_(0, vi) + visc_*timefacfac*derxyz_(1, ui)*derxyz_(1, vi) ;

    /* Stabilisierung der Viskosität (u*grad(u),-2*nu*div(epsilon((v)))) */
    estif_(vi*3, ui*3)         += 2.0*visc_*ttimetauMp*(conv_c_(ui)*viscs2_(0, 0, vi) + conv_g_(ui)*viscs2_(0, 0, vi) + viscs2_(0, 0, vi)*conv_r_(0, 0, ui) + viscs2_(0, 1, vi)*conv_r_(1, 0, ui) + viscs2_(0, 2, vi)*conv_r_(2, 0, ui)) ;
    estif_(vi*3, ui*3 + 1)     += 2.0*visc_*ttimetauMp*(conv_c_(ui)*viscs2_(0, 1, vi) + conv_g_(ui)*viscs2_(0, 1, vi) + viscs2_(0, 0, vi)*conv_r_(0, 1, ui) + viscs2_(0, 1, vi)*conv_r_(1, 1, ui) + viscs2_(0, 2, vi)*conv_r_(2, 1, ui)) ;
    estif_(vi*3, ui*3 + 2)     += 2.0*visc_*ttimetauMp*(conv_c_(ui)*viscs2_(0, 2, vi) + conv_g_(ui)*viscs2_(0, 2, vi) + viscs2_(0, 0, vi)*conv_r_(0, 2, ui) + viscs2_(0, 1, vi)*conv_r_(1, 2, ui) + viscs2_(0, 2, vi)*conv_r_(2, 2, ui)) ;
    estif_(vi*3 + 1, ui*3)     += 2.0*visc_*ttimetauMp*(conv_c_(ui)*viscs2_(0, 1, vi) + conv_g_(ui)*viscs2_(0, 1, vi) + viscs2_(0, 1, vi)*conv_r_(0, 0, ui) + viscs2_(1, 1, vi)*conv_r_(1, 0, ui) + viscs2_(1, 2, vi)*conv_r_(2, 0, ui)) ;
    estif_(vi*3 + 1, ui*3 + 1) += 2.0*visc_*ttimetauMp*(conv_c_(ui)*viscs2_(1, 1, vi) + conv_g_(ui)*viscs2_(1, 1, vi) + viscs2_(0, 1, vi)*conv_r_(0, 1, ui) + viscs2_(1, 1, vi)*conv_r_(1, 1, ui) + viscs2_(1, 2, vi)*conv_r_(2, 1, ui)) ;
    estif_(vi*3 + 1, ui*3 + 2) += 2.0*visc_*ttimetauMp*(conv_c_(ui)*viscs2_(1, 2, vi) + conv_g_(ui)*viscs2_(1, 2, vi) + viscs2_(0, 1, vi)*conv_r_(0, 2, ui) + viscs2_(1, 1, vi)*conv_r_(1, 2, ui) + viscs2_(1, 2, vi)*conv_r_(2, 2, ui)) ;
    estif_(vi*3 + 2, ui*3)     += 2.0*visc_*ttimetauMp*(conv_c_(ui)*viscs2_(0, 2, vi) + conv_g_(ui)*viscs2_(0, 2, vi) + viscs2_(0, 2, vi)*conv_r_(0, 0, ui) + viscs2_(1, 2, vi)*conv_r_(1, 0, ui) + viscs2_(2, 2, vi)*conv_r_(2, 0, ui)) ;
    estif_(vi*3 + 2, ui*3 + 1) += 2.0*visc_*ttimetauMp*(conv_c_(ui)*viscs2_(1, 2, vi) + conv_g_(ui)*viscs2_(1, 2, vi) + viscs2_(0, 2, vi)*conv_r_(0, 1, ui) + viscs2_(1, 2, vi)*conv_r_(1, 1, ui) + viscs2_(2, 2, vi)*conv_r_(2, 1, ui)) ;
    estif_(vi*3 + 2, ui*3 + 2) += 2.0*visc_*ttimetauMp*(conv_c_(ui)*viscs2_(2, 2, vi) + conv_g_(ui)*viscs2_(2, 2, vi) + viscs2_(0, 2, vi)*conv_r_(0, 2, ui) + viscs2_(1, 2, vi)*conv_r_(1, 2, ui) + viscs2_(2, 2, vi)*conv_r_(2, 2, ui)) ;

    /* Stabilisierung der Viskosität (-2*nu*div(epsilon((u))),-2*nu*div(epsilon((v)))) */
    estif_(vi*3, ui*3)         += 4.0*(visc_*visc_)*ttimetauMp*(viscs2_(0, 0, ui)*viscs2_(0, 0, vi) + viscs2_(0, 1, ui)*viscs2_(0, 1, vi) + viscs2_(0, 2, ui)*viscs2_(0, 2, vi)) ;
    estif_(vi*3, ui*3 + 1)     += 4.0*(visc_*visc_)*ttimetauMp*(viscs2_(0, 0, vi)*viscs2_(0, 1, ui) + viscs2_(0, 1, vi)*viscs2_(1, 1, ui) + viscs2_(0, 2, vi)*viscs2_(1, 2, ui)) ;
    estif_(vi*3, ui*3 + 2)     += 4.0*(visc_*visc_)*ttimetauMp*(viscs2_(0, 0, vi)*viscs2_(0, 2, ui) + viscs2_(0, 1, vi)*viscs2_(1, 2, ui) + viscs2_(0, 2, vi)*viscs2_(2, 2, ui)) ;
    estif_(vi*3 + 1, ui*3)     += 4.0*(visc_*visc_)*ttimetauMp*(viscs2_(0, 0, ui)*viscs2_(0, 1, vi) + viscs2_(0, 1, ui)*viscs2_(1, 1, vi) + viscs2_(0, 2, ui)*viscs2_(1, 2, vi)) ;
    estif_(vi*3 + 1, ui*3 + 1) += 4.0*(visc_*visc_)*ttimetauMp*(viscs2_(0, 1, ui)*viscs2_(0, 1, vi) + viscs2_(1, 1, ui)*viscs2_(1, 1, vi) + viscs2_(1, 2, ui)*viscs2_(1, 2, vi)) ;
    estif_(vi*3 + 1, ui*3 + 2) += 4.0*(visc_*visc_)*ttimetauMp*(viscs2_(0, 1, vi)*viscs2_(0, 2, ui) + viscs2_(1, 1, vi)*viscs2_(1, 2, ui) + viscs2_(1, 2, vi)*viscs2_(2, 2, ui)) ;
    estif_(vi*3 + 2, ui*3)     += 4.0*(visc_*visc_)*ttimetauMp*(viscs2_(0, 0, ui)*viscs2_(0, 2, vi) + viscs2_(0, 1, ui)*viscs2_(1, 2, vi) + viscs2_(0, 2, ui)*viscs2_(2, 2, vi)) ;
    estif_(vi*3 + 2, ui*3 + 1) += 4.0*(visc_*visc_)*ttimetauMp*(viscs2_(0, 1, ui)*viscs2_(0, 2, vi) + viscs2_(1, 1, ui)*viscs2_(1, 2, vi) + viscs2_(1, 2, ui)*viscs2_(2, 2, vi)) ;
    estif_(vi*3 + 2, ui*3 + 2) += 4.0*(visc_*visc_)*ttimetauMp*(viscs2_(0, 2, ui)*viscs2_(0, 2, vi) + viscs2_(1, 2, ui)*viscs2_(1, 2, vi) + viscs2_(2, 2, ui)*viscs2_(2, 2, vi)) ;

    /* Stabilisierung der Viskosität (grad(p),-2*nu*div(epsilon((v)))) */

    /* Kontinuitätsstabilisierung */
    estif_(vi*3, ui*3)         += (thsl*thsl)*tau_C*derxyz_(0, ui)*derxyz_(0, vi) ;
    estif_(vi*3, ui*3 + 1)     += (thsl*thsl)*tau_C*derxyz_(0, vi)*derxyz_(1, ui) ;
    estif_(vi*3, ui*3 + 2)     += (thsl*thsl)*tau_C*derxyz_(0, vi)*derxyz_(2, ui) ;
    estif_(vi*3 + 1, ui*3)     += (thsl*thsl)*tau_C*derxyz_(0, ui)*derxyz_(1, vi) ;
    estif_(vi*3 + 1, ui*3 + 1) += (thsl*thsl)*tau_C*derxyz_(1, ui)*derxyz_(1, vi) ;
    estif_(vi*3 + 1, ui*3 + 2) += (thsl*thsl)*tau_C*derxyz_(1, vi)*derxyz_(2, ui) ;
    estif_(vi*3 + 2, ui*3)     += (thsl*thsl)*tau_C*derxyz_(0, ui)*derxyz_(2, vi) ;
    estif_(vi*3 + 2, ui*3 + 1) += (thsl*thsl)*tau_C*derxyz_(1, ui)*derxyz_(2, vi) ;
    estif_(vi*3 + 2, ui*3 + 2) += (thsl*thsl)*tau_C*derxyz_(2, ui)*derxyz_(2, vi) ;

    /* Massenterm (u,v) */
    estif_(vi*3, ui*3)         += fac*funct_(ui)*funct_(vi) ;
    estif_(vi*3 + 1, ui*3 + 1) += fac*funct_(ui)*funct_(vi) ;
    estif_(vi*3 + 2, ui*3 + 2) += fac*funct_(ui)*funct_(vi) ;

    /* Konvektionsstabilisierung (u,u*grad(v)) */
    estif_(vi*3, ui*3)         += timetauM*funct_(ui)*(conv_c_(vi) + conv_g_(vi) + velint_(0)*derxyz_(0, vi)) ;
    estif_(vi*3, ui*3 + 1)     += timetauM*funct_(ui)*velint_(0)*derxyz_(1, vi) ;
    estif_(vi*3, ui*3 + 2)     += timetauM*funct_(ui)*velint_(0)*derxyz_(2, vi) ;
    estif_(vi*3 + 1, ui*3)     += timetauM*funct_(ui)*velint_(1)*derxyz_(0, vi) ;
    estif_(vi*3 + 1, ui*3 + 1) += timetauM*funct_(ui)*(conv_c_(vi) + conv_g_(vi) + velint_(1)*derxyz_(1, vi)) ;
    estif_(vi*3 + 1, ui*3 + 2) += timetauM*funct_(ui)*velint_(1)*derxyz_(2, vi) ;
    estif_(vi*3 + 2, ui*3)     += timetauM*funct_(ui)*velint_(2)*derxyz_(0, vi) ;
    estif_(vi*3 + 2, ui*3 + 1) += timetauM*funct_(ui)*velint_(2)*derxyz_(1, vi) ;
    estif_(vi*3 + 2, ui*3 + 2) += timetauM*funct_(ui)*(conv_c_(vi) + conv_g_(vi) + velint_(2)*derxyz_(2, vi)) ;

    /* Viskositätsstabilisierung (u,-2*nu*div(epsilon((v)))) */
    estif_(vi*3, ui*3)         += tau_Mp*time2nue*funct_(ui)*viscs2_(0, 0, vi) ;
    estif_(vi*3, ui*3 + 1)     += tau_Mp*time2nue*funct_(ui)*viscs2_(0, 1, vi) ;
    estif_(vi*3, ui*3 + 2)     += tau_Mp*time2nue*funct_(ui)*viscs2_(0, 2, vi) ;
    estif_(vi*3 + 1, ui*3)     += tau_Mp*time2nue*funct_(ui)*viscs2_(0, 1, vi) ;
    estif_(vi*3 + 1, ui*3 + 1) += tau_Mp*time2nue*funct_(ui)*viscs2_(1, 1, vi) ;
    estif_(vi*3 + 1, ui*3 + 2) += tau_Mp*time2nue*funct_(ui)*viscs2_(1, 2, vi) ;
    estif_(vi*3 + 2, ui*3)     += tau_Mp*time2nue*funct_(ui)*viscs2_(0, 2, vi) ;
    estif_(vi*3 + 2, ui*3 + 1) += tau_Mp*time2nue*funct_(ui)*viscs2_(1, 2, vi) ;
    estif_(vi*3 + 2, ui*3 + 2) += tau_Mp*time2nue*funct_(ui)*viscs2_(2, 2, vi) ;

    /* Quellterm der rechten Seite (b,v) */

    /* Konvektionsstabilisierung (b,u*grad(v)) */
    estif_(vi*3, ui*3)         += -(timetauM*funct_(ui)*rhsint_(0)*derxyz_(0, vi)) ;
    estif_(vi*3, ui*3 + 1)     += -(timetauM*funct_(ui)*rhsint_(0)*derxyz_(1, vi)) ;
    estif_(vi*3, ui*3 + 2)     += -(timetauM*funct_(ui)*rhsint_(0)*derxyz_(2, vi)) ;
    estif_(vi*3 + 1, ui*3)     += -(timetauM*funct_(ui)*rhsint_(1)*derxyz_(0, vi)) ;
    estif_(vi*3 + 1, ui*3 + 1) += -(timetauM*funct_(ui)*rhsint_(1)*derxyz_(1, vi)) ;
    estif_(vi*3 + 1, ui*3 + 2) += -(timetauM*funct_(ui)*rhsint_(1)*derxyz_(2, vi)) ;
    estif_(vi*3 + 2, ui*3)     += -(timetauM*funct_(ui)*rhsint_(2)*derxyz_(0, vi)) ;
    estif_(vi*3 + 2, ui*3 + 1) += -(timetauM*funct_(ui)*rhsint_(2)*derxyz_(1, vi)) ;
    estif_(vi*3 + 2, ui*3 + 2) += -(timetauM*funct_(ui)*rhsint_(2)*derxyz_(2, vi)) ;

    /* Viskositätsstabilisierung (b,-2*nu*div(epsilon((v))) */

  }
}
#endif
