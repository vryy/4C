for (int vi=0; vi<iel; ++vi)
{
    for (int ui=0; ui<iel; ++ui)
    {
    /* Konvektionsterm */
    estif_(vi*3, ui*3)         += fac*funct_(vi)*(conv_c_(ui) + conv_r_(0, 0, ui)) ;
    estif_(vi*3, ui*3 + 1)     += fac*funct_(vi)*conv_r_(0, 1, ui) ;
    estif_(vi*3 + 1, ui*3)     += fac*funct_(vi)*conv_r_(1, 0, ui) ;
    estif_(vi*3 + 1, ui*3 + 1) += fac*funct_(vi)*(conv_c_(ui) + conv_r_(1, 1, ui)) ;

    /* Viskositaetsterm */
    estif_(vi*3, ui*3)         += fac*2.0*nu_*derxy_(0, ui)*derxy_(0, vi) +
	                          nu_*fac*derxy_(1, ui)*derxy_(1, vi) ;
    estif_(vi*3, ui*3 + 1)     += nu_*fac*derxy_(0, ui)*derxy_(1, vi) ;
    estif_(vi*3 + 1, ui*3)     += nu_*fac*derxy_(0, vi)*derxy_(1, ui) ;
    estif_(vi*3 + 1, ui*3 + 1) += fac*2.0*nu_*derxy_(1, ui)*derxy_(1, vi) +
	                          nu_*fac*derxy_(0, ui)*derxy_(0, vi) ;

    /* Druckterm */
    estif_(vi*3, ui*3 + 2)     += -(fac*funct_(ui)*derxy_(0, vi)) ;
    estif_(vi*3 + 1, ui*3 + 2) += -(fac*funct_(ui)*derxy_(1, vi)) ;

    /* Divergenzfreiheit */
    estif_(vi*3 + 2, ui*3)     += fac*funct_(vi)*derxy_(0, ui) ;
    estif_(vi*3 + 2, ui*3 + 1) += fac*funct_(vi)*derxy_(1, ui) ;

    /* konvektive Stabilisierung: */
    /* konvektive Stabilisierung der Konvektion */
    estif_(vi*3, ui*3)         += tau_M*(conv_c_(ui)*conv_c_(vi) +
					     conv_c_(vi)*conv_r_(0, 0, ui) +
					     velint_(0)*derxy_(0, vi)*conv_r_(0, 0, ui) +
					     velint_(1)*derxy_(0, vi)*conv_r_(0, 1, ui)) ;
    estif_(vi*3, ui*3 + 1)     += tau_M*(conv_c_(vi)*conv_r_(0, 1, ui) +
					     velint_(0)*derxy_(1, vi)*conv_r_(0, 0, ui) +
					     velint_(1)*derxy_(1, vi)*conv_r_(0, 1, ui)) ;
    estif_(vi*3 + 1, ui*3)     += tau_M*(conv_c_(vi)*conv_r_(1, 0, ui) +
					     velint_(0)*derxy_(0, vi)*conv_r_(1, 0, ui) +
					     velint_(1)*derxy_(0, vi)*conv_r_(1, 1, ui)) ;
    estif_(vi*3 + 1, ui*3 + 1) += tau_M*(conv_c_(ui)*conv_c_(vi) +
					     conv_c_(vi)*conv_r_(1, 1, ui) +
					     velint_(0)*derxy_(1, vi)*conv_r_(1, 0, ui) +
					     velint_(1)*derxy_(1, vi)*conv_r_(1, 1, ui)) ;

    /* konvektive Stabilisierung der Viskositaet */
    estif_(vi*3, ui*3)         += -2.0*nu_*tau_M*(-(conv_c_(vi)*viscs2_(0, 0, ui)) +
						      funct_(ui)*derxy_(0, vi)*visc_old_(0)) ;
    estif_(vi*3, ui*3 + 1)     += -2.0*nu_*tau_M*(-(conv_c_(vi)*viscs2_(0, 1, ui)) +
						      funct_(ui)*derxy_(1, vi)*visc_old_(0)) ;
    estif_(vi*3 + 1, ui*3)     += -2.0*nu_*tau_M*(-(conv_c_(vi)*viscs2_(0, 1, ui)) +
						      funct_(ui)*derxy_(0, vi)*visc_old_(1)) ;
    estif_(vi*3 + 1, ui*3 + 1) += -2.0*nu_*tau_M*(-(conv_c_(vi)*viscs2_(1, 1, ui)) +
						      funct_(ui)*derxy_(1, vi)*visc_old_(1)) ;

    /* konvektive Stabilisierung des Drucks */
    estif_(vi*3, ui*3)         += tau_M*funct_(ui)*gradp_(0)*derxy_(0, vi) ;
    estif_(vi*3, ui*3 + 1)     += tau_M*funct_(ui)*gradp_(0)*derxy_(1, vi) ;
    estif_(vi*3, ui*3 + 2)     += tau_M*conv_c_(vi)*derxy_(0, ui) ;
    estif_(vi*3 + 1, ui*3)     += tau_M*funct_(ui)*gradp_(1)*derxy_(0, vi) ;
    estif_(vi*3 + 1, ui*3 + 1) += tau_M*funct_(ui)*gradp_(1)*derxy_(1, vi) ;
    estif_(vi*3 + 1, ui*3 + 2) += tau_M*conv_c_(vi)*derxy_(1, ui) ;

    /* Druck-Stabilisierung: */
    /* Druck-Stabilisierung der Konvektion */
    estif_(vi*3 + 2, ui*3)     += tau_Mp*(conv_c_(ui)*derxy_(0, vi) +
					      derxy_(0, vi)*conv_r_(0, 0, ui) +
					      derxy_(1, vi)*conv_r_(1, 0, ui)) ;
    estif_(vi*3 + 2, ui*3 + 1) += tau_Mp*(conv_c_(ui)*derxy_(1, vi) +
					      derxy_(0, vi)*conv_r_(0, 1, ui) +
					      derxy_(1, vi)*conv_r_(1, 1, ui)) ;

    /* Druck-Stabilisierung der Viskositaet */
    estif_(vi*3 + 2, ui*3)     += 2.0*nu_*tau_Mp*(derxy_(0, vi)*viscs2_(0, 0, ui) +
						      derxy_(1, vi)*viscs2_(0, 1, ui)) ;
    estif_(vi*3 + 2, ui*3 + 1) += 2.0*nu_*tau_Mp*(derxy_(0, vi)*viscs2_(0, 1, ui) +
						      derxy_(1, vi)*viscs2_(1, 1, ui)) ;

    /* Druck-Stabilisierung des Drucks */
    estif_(vi*3 + 2, ui*3 + 2) += tau_Mp*(derxy_(0, ui)*derxy_(0, vi) +
					      derxy_(1, ui)*derxy_(1, vi)) ;

    /* Kontinuitaetsstabilisierung */
    estif_(vi*3, ui*3)         += tau_C*derxy_(0, ui)*derxy_(0, vi) ;
    estif_(vi*3, ui*3 + 1)     += tau_C*derxy_(0, vi)*derxy_(1, ui) ;
    estif_(vi*3 + 1, ui*3)     += tau_C*derxy_(0, ui)*derxy_(1, vi) ;
    estif_(vi*3 + 1, ui*3 + 1) += tau_C*derxy_(1, ui)*derxy_(1, vi) ;

  }
}
