for (int vi=0; vi<iel; ++vi)
{
  /* Konvektionsterm */
  eforce_(vi*3)     += -(fac*funct_(vi)*conv_old_(0)) ;
  eforce_(vi*3 + 1) += -(fac*funct_(vi)*conv_old_(1)) ;

  /* Viskositaetsterm */
  eforce_(vi*3)     += -(nu_*fac*(2.0*derxy_(0, vi)*vderxy_(0, 0) + derxy_(1, vi)*vderxy_(0, 1) + derxy_(1, vi)*vderxy_(1, 0))) ;
  eforce_(vi*3 + 1) += -(nu_*fac*(derxy_(0, vi)*vderxy_(0, 1) + derxy_(0, vi)*vderxy_(1, 0) + 2.0*derxy_(1, vi)*vderxy_(1, 1))) ;

  /* Druckterm */
  eforce_(vi*3)     += press*fac*derxy_(0, vi) ;
  eforce_(vi*3 + 1) += press*fac*derxy_(1, vi) ;

  /* Divergenzfreiheit */
  eforce_(vi*3 + 2) += -(fac*(conv_r_(0, 0, vi) + conv_r_(1, 1, vi))) ;

  /* konvektive Stabilisierung: */
  /* konvektive Stabilisierung der Konvektion */
  eforce_(vi*3)     += -(tau_M*conv_c_(vi)*conv_old_(0)) ;
  eforce_(vi*3 + 1) += -(tau_M*conv_c_(vi)*conv_old_(1)) ;

  /* konvektive Stabilisierung der Viskositaet */
  eforce_(vi*3)     += 2.0*nu_*tau_M*conv_c_(vi)*visc_old_(0) ;
  eforce_(vi*3 + 1) += 2.0*nu_*tau_M*conv_c_(vi)*visc_old_(1) ;

  /* konvektive Stabilisierung des Drucks */
  eforce_(vi*3)     += -(tau_M*conv_c_(vi)*gradp_(0)) ;
  eforce_(vi*3 + 1) += -(tau_M*conv_c_(vi)*gradp_(1)) ;

  /* Druck-Stabilisierung: */
  /* Druck-Stabilisierung der Konvektion */
  eforce_(vi*3 + 2) += -(tau_Mp*(derxy_(0, vi)*conv_old_(0) + derxy_(1, vi)*conv_old_(1))) ;

  /* Druck-Stabilisierung der Viskositaet */
  eforce_(vi*3 + 2) += 2.0*nu_*tau_Mp*(derxy_(0, vi)*visc_old_(0) + derxy_(1, vi)*visc_old_(1)) ;

  /* Druck-Stabilisierung des Drucks */
  eforce_(vi*3 + 2) += -(tau_Mp*(gradp_(0)*derxy_(0, vi) + gradp_(1)*derxy_(1, vi))) ;

  /* Kontinuitaetsstabilisierung */
  eforce_(vi*3)     += -(tau_C*derxy_(0, vi)*(vderxy_(0, 0) + vderxy_(1, 1))) ;
  eforce_(vi*3 + 1) += -(tau_C*derxy_(1, vi)*(vderxy_(0, 0) + vderxy_(1, 1))) ;

}
