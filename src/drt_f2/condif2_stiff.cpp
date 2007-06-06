for (int vi=0; vi<iel; ++vi)
{
    for (int ui=0; ui<iel; ++ui)
    {
   
    /* Standard Galerkin terms: */
    /* transient term */
    estif_(vi, ui) += timefacl*funct_(vi)*funct_(ui) ;

    /* convective term */
    estif_(vi, ui) += fac*funct_(vi)*conv_(ui) ;

    /* diffusive term */
    estif_(vi, ui) += fac*diffus*(derxy_(0, ui)*derxy_(0, vi) + derxy_(1, ui)*derxy_(1, vi)) ;

    /* Stabilization terms: */
    /* 1) transient stabilization (USFEM assumed here, sign change necessary for GLS) */
    /* transient term */
    estif_(vi, ui) += -ttimetaufacl*funct_(vi)*funct_(ui) ;

    /* convective term */
    estif_(vi, ui) += -timetaufacl*funct_(vi)*conv_(ui) ;

    /* diffusive term */
    estif_(vi, ui) += timetaufacl*funct_(vi)*diff_(ui) ;

    /* 2) convective stabilization */
    /* transient term */
    estif_(vi, ui) += timetaufacl*conv_(vi)*funct_(ui) ;

    /* convective term */
    estif_(vi, ui) += taufac*conv_(vi)*conv_(ui) ;

    /* diffusive term */
    estif_(vi, ui) += -taufac*conv_(vi)*diff_(ui) ;

    /* 3) diffusive stabilization (USFEM assumed here, sign change necessary for GLS) */
    /* transient term */
    estif_(vi, ui) += timetaufacl*diff_(vi)*funct_(ui) ;

    /* convective term */
    estif_(vi, ui) += taufac*diff_(vi)*conv_(ui) ;

    /* diffusive term */
    estif_(vi, ui) += -taufac*diff_(vi)*diff_(ui) ;

  }
}
