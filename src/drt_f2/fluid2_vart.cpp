for (int vi=0; vi<iel; ++vi)
{
    for (int ui=0; ui<iel; ++ui)
    {

    /* artificial viscosity term */
    edc_(vi*3, ui*3)         += vartfac*derxy_(0, vi)*derxy_(0, ui)  ;
    edc_(vi*3, ui*3 + 1)     += vartfac*derxy_(1, vi)*derxy_(0, ui)  ;
    edc_(vi*3 + 1, ui*3)     += vartfac*derxy_(0, vi)*derxy_(1, ui)  ;
    edc_(vi*3 + 1, ui*3 + 1) += vartfac*derxy_(1, vi)*derxy_(1, ui)  ;

    /*subtract SUPG term */
    edc_(vi*3, ui*3)         -= taumfac*conv_c_(vi)*conv_c_(ui) ;
    edc_(vi*3 + 1, ui*3 + 1) -= taumfac*conv_c_(vi)*conv_c_(ui) ;

    }
}
