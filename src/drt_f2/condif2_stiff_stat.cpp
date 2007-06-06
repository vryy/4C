for (int vi=0; vi<iel; ++vi)
{
    for (int ui=0; ui<iel; ++ui)
    {
   
    /* Standard Galerkin terms: */
    /* convective term */
    estif_(vi, ui) += fac*funct_(vi)*conv_(ui) ;

    /* diffusive term */
    estif_(vi, ui) += fac*diffus*(derxy_(0, ui)*derxy_(0, vi) + derxy_(1, ui)*derxy_(1, vi)) ;

    /* Stabilization terms: */
    /* 1) convective stabilization */
    /* convective term */
    estif_(vi, ui) += taufac*conv_(vi)*conv_(ui) ;

    /* diffusive term */
    estif_(vi, ui) += -taufac*conv_(vi)*diff_(ui) ;

    /* 2) diffusive stabilization (USFEM assumed here, sign change necessary for GLS) */
    /* convective term */
    estif_(vi, ui) += taufac*diff_(vi)*conv_(ui) ;

    /* diffusive term */
    estif_(vi, ui) += -taufac*diff_(vi)*diff_(ui) ;

  }
}
