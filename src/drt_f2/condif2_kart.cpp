for (int vi=0; vi<iel; ++vi)
{
    for (int ui=0; ui<iel; ++ui)
    {

    /* artificial diffusivity term */
    edc_(vi, ui) += kartfac*(derxy_(0, ui)*derxy_(0, vi) + derxy_(1, ui)*derxy_(1, vi)) ;

    }
}
