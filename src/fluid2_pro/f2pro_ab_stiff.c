for (vi=0; vi<iel; ++vi)
{
  for (ui=0; ui<iel; ++ui)
  {
#if 0
    /* Konvektionsterm (u*grad(u),v) */
    estif_(vi*2, ui*2)         += timefacfac*funct_(vi)*(conv_c_(ui) + conv_r_(0, 0, ui)) ;
    estif_(vi*2, ui*2 + 1)     += timefacfac*funct_(vi)*conv_r_(0, 1, ui) ;
    estif_(vi*2 + 1, ui*2)     += timefacfac*funct_(vi)*conv_r_(1, 0, ui) ;
    estif_(vi*2 + 1, ui*2 + 1) += timefacfac*funct_(vi)*(conv_c_(ui) + conv_r_(1, 1, ui)) ;
#endif

    /* Viskositaetsterm (2*nu*epsilon(u),epsilon(v)) */
    estif_(vi*2, ui*2)         += fac*time2nue*derxy_(0, ui)*derxy_(0, vi) + visc_*timefacfac*derxy_(1, ui)*derxy_(1, vi) ;
    estif_(vi*2, ui*2 + 1)     += visc_*timefacfac*derxy_(0, ui)*derxy_(1, vi) ;
    estif_(vi*2 + 1, ui*2)     += visc_*timefacfac*derxy_(0, vi)*derxy_(1, ui) ;
    estif_(vi*2 + 1, ui*2 + 1) += fac*time2nue*derxy_(1, ui)*derxy_(1, vi) + visc_*timefacfac*derxy_(0, ui)*derxy_(0, vi) ;

    /* Massenterm (u,v) */
    estif_(vi*2, ui*2)         += fac*funct_(ui)*funct_(vi) ;
    estif_(vi*2 + 1, ui*2 + 1) += fac*funct_(ui)*funct_(vi) ;
  }
}
