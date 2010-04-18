#ifndef CCADISCRET
for (vi=0; vi<iel; ++vi)
{
  for (ui=0; ui<iel; ++ui)
  {
#if 0
    /* Konvektionsterm (u*grad(u),v) */
    estif_(vi*3, ui*3)         += timefacfac*funct_(vi)*(conv_c_(ui) + conv_r_(0, 0, ui)) ;
    estif_(vi*3, ui*3 + 1)     += timefacfac*funct_(vi)*conv_r_(0, 1, ui) ;
    estif_(vi*3, ui*3 + 2)     += timefacfac*funct_(vi)*conv_r_(0, 2, ui) ;
    estif_(vi*3 + 1, ui*3)     += timefacfac*funct_(vi)*conv_r_(1, 0, ui) ;
    estif_(vi*3 + 1, ui*3 + 1) += timefacfac*funct_(vi)*(conv_c_(ui) + conv_r_(1, 1, ui)) ;
    estif_(vi*3 + 1, ui*3 + 2) += timefacfac*funct_(vi)*conv_r_(1, 2, ui) ;
    estif_(vi*3 + 2, ui*3)     += timefacfac*funct_(vi)*conv_r_(2, 0, ui) ;
    estif_(vi*3 + 2, ui*3 + 1) += timefacfac*funct_(vi)*conv_r_(2, 1, ui) ;
    estif_(vi*3 + 2, ui*3 + 2) += timefacfac*funct_(vi)*(conv_c_(ui) + conv_r_(2, 2, ui)) ;
#endif

    /* Viskositaetsterm (2*nu*epsilon(u),epsilon(v)) */
    estif_(vi*3, ui*3)         += fac*time2nue*derxyz_(0, ui)*derxyz_(0, vi) + visc_*timefacfac*derxyz_(1, ui)*derxyz_(1, vi) + visc_*timefacfac*derxyz_(2, ui)*derxyz_(2, vi) ;
    estif_(vi*3, ui*3 + 1)     += visc_*timefacfac*derxyz_(0, ui)*derxyz_(1, vi) ;
    estif_(vi*3, ui*3 + 2)     += visc_*timefacfac*derxyz_(0, ui)*derxyz_(2, vi) ;
    estif_(vi*3 + 1, ui*3)     += visc_*timefacfac*derxyz_(0, vi)*derxyz_(1, ui) ;
    estif_(vi*3 + 1, ui*3 + 1) += fac*time2nue*derxyz_(1, ui)*derxyz_(1, vi) + visc_*timefacfac*derxyz_(0, ui)*derxyz_(0, vi) + visc_*timefacfac*derxyz_(2, ui)*derxyz_(2, vi) ;
    estif_(vi*3 + 1, ui*3 + 2) += visc_*timefacfac*derxyz_(1, ui)*derxyz_(2, vi) ;
    estif_(vi*3 + 2, ui*3)     += visc_*timefacfac*derxyz_(0, vi)*derxyz_(2, ui) ;
    estif_(vi*3 + 2, ui*3 + 1) += visc_*timefacfac*derxyz_(1, vi)*derxyz_(2, ui) ;
    estif_(vi*3 + 2, ui*3 + 2) += fac*time2nue*derxyz_(2, ui)*derxyz_(2, vi) + visc_*timefacfac*derxyz_(0, ui)*derxyz_(0, vi) + visc_*timefacfac*derxyz_(1, ui)*derxyz_(1, vi) ;

    /* Massenterm (u,v) */
    estif_(vi*3, ui*3)         += fac*funct_(ui)*funct_(vi) ;
    estif_(vi*3 + 1, ui*3 + 1) += fac*funct_(ui)*funct_(vi) ;
    estif_(vi*3 + 2, ui*3 + 2) += fac*funct_(ui)*funct_(vi) ;
  }
}
#endif
