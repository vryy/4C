for (int vi=0; vi<iel; ++vi)
{
    for (int ui=0; ui<iel; ++ui)
    {

    /* artificial viscosity term on left hand side (stress-divergence form) */
    esv_(vi*3, ui*3)         += vartfac*(2.0*derxy_(0, vi)*derxy_(0, ui) + derxy_(1, vi)*derxy_(1, ui)) ;
    esv_(vi*3, ui*3 + 1)     += vartfac*derxy_(1, vi)*derxy_(0, ui) ;
    esv_(vi*3 + 1, ui*3)     += vartfac*derxy_(0, vi)*derxy_(1, ui) ;
    esv_(vi*3 + 1, ui*3 + 1) += vartfac*(2.0*derxy_(1, vi)*derxy_(1, ui) + derxy_(0, vi)*derxy_(0, ui)) ;

    /* artificial viscosity term on left hand side (conventional form) */
    //esv_(vi*3, ui*3)         += vartfac*derxy_(0, vi)*derxy_(0, ui)  ;
    //esv_(vi*3, ui*3 + 1)     += vartfac*derxy_(1, vi)*derxy_(0, ui)  ;
    //esv_(vi*3 + 1, ui*3)     += vartfac*derxy_(0, vi)*derxy_(1, ui)  ;
    //esv_(vi*3 + 1, ui*3 + 1) += vartfac*derxy_(1, vi)*derxy_(1, ui)  ;

    /*subtracted SUPG term on left hand side */
    esv_(vi*3, ui*3)         -= taumfac*conv_c_(vi)*conv_c_(ui) ;
    esv_(vi*3 + 1, ui*3 + 1) -= taumfac*conv_c_(vi)*conv_c_(ui) ;

    }

}


