for (vi=0; vi<iel; ++vi)
{
#if 0
    /* no neumann or body forces right now! */

    /* Quellterm der rechten Seite (b,v) */
    eforce_(vi*3)     += fac*funct_(vi)*rhsint_(0) ;
    eforce_(vi*3 + 1) += fac*funct_(vi)*rhsint_(1) ;
    eforce_(vi*3 + 2) += fac*funct_(vi)*rhsint_(2) ;
#endif

    /*
     * Here we have to add terms from the time step beginning because
     * we do not have a mass rhs. I might be doing something stupid. */

    /* Viskositätsterm (2*nu*epsilon(u),epsilon(v)) */
    eforce_(vi*3)     += (visc_*fdyn->thsr*fac*(2.0*derxyz_(0, vi)*vderxyz_n_(0, 0) +
                                                    derxyz_(1, vi)*vderxyz_n_(0, 1) +
                                                    derxyz_(1, vi)*vderxyz_n_(1, 0) +
                                                    derxyz_(2, vi)*vderxyz_n_(0, 2) +
                                                    derxyz_(2, vi)*vderxyz_n_(2, 0))) ;
    eforce_(vi*3 + 1) += (visc_*fdyn->thsr*fac*(    derxyz_(0, vi)*vderxyz_n_(0, 1) +
                                                    derxyz_(0, vi)*vderxyz_n_(1, 0) +
                                                2.0*derxyz_(1, vi)*vderxyz_n_(1, 1) +
                                                    derxyz_(2, vi)*vderxyz_n_(1, 2) +
                                                    derxyz_(2, vi)*vderxyz_n_(2, 1))) ;
    eforce_(vi*3 + 2) += (visc_*fdyn->thsr*fac*(    derxyz_(0, vi)*vderxyz_n_(0, 2) +
                                                    derxyz_(0, vi)*vderxyz_n_(2, 0) +
                                                    derxyz_(1, vi)*vderxyz_n_(1, 2) +
                                                    derxyz_(1, vi)*vderxyz_n_(2, 1) +
                                                2.0*derxyz_(2, vi)*vderxyz_n_(2, 2))) ;

    /* Massenterm (u,v) */
    eforce_(vi*3)     += (fac*funct_(vi)*velint_n_(0)) ;
    eforce_(vi*3 + 1) += (fac*funct_(vi)*velint_n_(1)) ;
    eforce_(vi*3 + 2) += (fac*funct_(vi)*velint_n_(2)) ;

#if 0
    /* no neumann or body forces right now! */

    /* Quellterm der rechten Seite (b,v) */
    eforce_(vi*3)     += fac*funct_(vi)*rhsint_n_(0) ;
    eforce_(vi*3 + 1) += fac*funct_(vi)*rhsint_n_(1) ;
    eforce_(vi*3 + 2) += fac*funct_(vi)*rhsint_n_(2) ;
#endif

    /* Druckterm */
    eforce_(vi*3)     += -press_n*fdyn->thsr*fac*derxyz_(0, vi) ;
    eforce_(vi*3 + 1) += -press_n*fdyn->thsr*fac*derxyz_(1, vi) ;
    eforce_(vi*3 + 2) += -press_n*fdyn->thsr*fac*derxyz_(2, vi) ;

    /* Konvektionsterm -3/2 dt (u*grad(u),v) at (n) */
    eforce_(vi*3)     += -3./2.*fdyn->dta*fac*funct_(vi)*conv_old_n_(0) ;
    eforce_(vi*3 + 1) += -3./2.*fdyn->dta*fac*funct_(vi)*conv_old_n_(1) ;
    eforce_(vi*3 + 2) += -3./2.*fdyn->dta*fac*funct_(vi)*conv_old_n_(2) ;

    /* Konvektionsterm  1/2 dt (u*grad(u),v) at (n-1) */
    eforce_(vi*3)     +=  1./2.*fdyn->dta*fac*funct_(vi)*conv_old_nm_(0) ;
    eforce_(vi*3 + 1) +=  1./2.*fdyn->dta*fac*funct_(vi)*conv_old_nm_(1) ;
    eforce_(vi*3 + 2) +=  1./2.*fdyn->dta*fac*funct_(vi)*conv_old_nm_(2) ;
}
