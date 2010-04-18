#ifndef CCADISCRET
for (vi=0; vi<iel; ++vi)
{
#if 0
    /* no neumann or body forces right now! */

    /* Quellterm der rechten Seite (b,v) */
    eforce_(vi*2)     += fac*funct_(vi)*rhsint_(0) ;
    eforce_(vi*2 + 1) += fac*funct_(vi)*rhsint_(1) ;
#endif

    /*
     * Here we have to add terms from the time step beginning because
     * we do not have a mass rhs. I might be doing something stupid. */

    /* Viskositaetsterm (2*nu*epsilon(u),epsilon(v)) */
    eforce_(vi*2)     += (visc_*fdyn->thsr*fac*(2.0*derxy_(0, vi)*vderxy_n_(0, 0) +
						    derxy_(1, vi)*vderxy_n_(0, 1) +
						    derxy_(1, vi)*vderxy_n_(1, 0))) ;
    eforce_(vi*2 + 1) += (visc_*fdyn->thsr*fac*(    derxy_(0, vi)*vderxy_n_(0, 1) +
						    derxy_(0, vi)*vderxy_n_(1, 0) +
						2.0*derxy_(1, vi)*vderxy_n_(1, 1))) ;

    /* Massenterm (u,v) */
    eforce_(vi*2)     += (fac*funct_(vi)*velint_n_(0)) ;
    eforce_(vi*2 + 1) += (fac*funct_(vi)*velint_n_(1)) ;

#if 0
    /* no neumann or body forces right now! */

    /* Quellterm der rechten Seite (b,v) */
    eforce_(vi*2)     += fac*funct_(vi)*rhsint_n_(0) ;
    eforce_(vi*2 + 1) += fac*funct_(vi)*rhsint_n_(1) ;
#endif

    eforce_(vi*2    ) += -press_n*fdyn->thsr*fac*derxy_(0, vi) ;
    eforce_(vi*2 + 1) += -press_n*fdyn->thsr*fac*derxy_(1, vi) ;

    /* Konvektionsterm -3/2 dt (u*grad(u),v) at (n) */
    eforce_(vi*2)     += -3./2.*fdyn->dta*fac*funct_(vi)*conv_old_n_(0) ;
    eforce_(vi*2 + 1) += -3./2.*fdyn->dta*fac*funct_(vi)*conv_old_n_(1) ;

    /* Konvektionsterm  1/2 dt (u*grad(u),v) at (n-1) */
    eforce_(vi*2)     +=  1./2.*fdyn->dta*fac*funct_(vi)*conv_old_nm_(0) ;
    eforce_(vi*2 + 1) +=  1./2.*fdyn->dta*fac*funct_(vi)*conv_old_nm_(1) ;
}
#endif
