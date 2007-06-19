for (int ui=0; ui<iel; ++ui)    // exchanged ui- and vi-loop => slight speedup!?   g.bau 03/07
{
  const int UDF = ui*4;
  const int VDF = UDF + 1;
  const int WDF = UDF + 2;
  const int PDF = UDF + 3;
  
  for (int vi=0; vi<iel; ++vi)
  {
    const int dUDF = vi*4;
    const int dVDF = dUDF + 1;
    const int dWDF = dUDF + 2;
    const int dPDF = dUDF + 3;
    
    /* Konvektionsterm */
    estif_(dUDF, UDF) += timefacfac*funct_(vi)*(conv_c_(ui) + conv_r_(0, 0, ui)) ;
    estif_(dUDF, VDF) += timefacfac*funct_(vi)*conv_r_(0, 1, ui) ;
    estif_(dUDF, WDF) += timefacfac*funct_(vi)*conv_r_(0, 2, ui) ;
    estif_(dVDF, UDF) += timefacfac*funct_(vi)*conv_r_(1, 0, ui) ;
    estif_(dVDF, VDF) += timefacfac*funct_(vi)*(conv_c_(ui) + conv_r_(1, 1, ui)) ;
    estif_(dVDF, WDF) += timefacfac*funct_(vi)*conv_r_(1, 2, ui) ;
    estif_(dWDF, UDF) += timefacfac*funct_(vi)*conv_r_(2, 0, ui) ;
    estif_(dWDF, VDF) += timefacfac*funct_(vi)*conv_r_(2, 1, ui) ;
    estif_(dWDF, WDF) += timefacfac*funct_(vi)*(conv_c_(ui) + conv_r_(2, 2, ui)) ;

    /* Stabilisierung der Konvektion ( L_conv_u) */
    estif_(dUDF, UDF) += ttimetauM*(conv_c_(ui)*conv_c_(vi) + conv_c_(vi)*conv_r_(0, 0, ui) + velint_(0)*derxyz_(0, vi)*conv_r_(0, 0, ui) + velint_(1)*derxyz_(0, vi)*conv_r_(0, 1, ui) + velint_(2)*derxyz_(0, vi)*conv_r_(0, 2, ui)) ;
    estif_(dUDF, VDF) += ttimetauM*(conv_c_(vi)*conv_r_(0, 1, ui) + velint_(0)*derxyz_(1, vi)*conv_r_(0, 0, ui) + velint_(1)*derxyz_(1, vi)*conv_r_(0, 1, ui) + velint_(2)*derxyz_(1, vi)*conv_r_(0, 2, ui)) ;
    estif_(dUDF, WDF) += ttimetauM*(conv_c_(vi)*conv_r_(0, 2, ui) + velint_(0)*derxyz_(2, vi)*conv_r_(0, 0, ui) + velint_(1)*derxyz_(2, vi)*conv_r_(0, 1, ui) + velint_(2)*derxyz_(2, vi)*conv_r_(0, 2, ui)) ;
    estif_(dVDF, UDF) += ttimetauM*(conv_c_(vi)*conv_r_(1, 0, ui) + velint_(0)*derxyz_(0, vi)*conv_r_(1, 0, ui) + velint_(1)*derxyz_(0, vi)*conv_r_(1, 1, ui) + velint_(2)*derxyz_(0, vi)*conv_r_(1, 2, ui)) ;
    estif_(dVDF, VDF) += ttimetauM*(conv_c_(ui)*conv_c_(vi) + conv_c_(vi)*conv_r_(1, 1, ui) + velint_(0)*derxyz_(1, vi)*conv_r_(1, 0, ui) + velint_(1)*derxyz_(1, vi)*conv_r_(1, 1, ui) + velint_(2)*derxyz_(1, vi)*conv_r_(1, 2, ui)) ;
    estif_(dVDF, WDF) += ttimetauM*(conv_c_(vi)*conv_r_(1, 2, ui) + velint_(0)*derxyz_(2, vi)*conv_r_(1, 0, ui) + velint_(1)*derxyz_(2, vi)*conv_r_(1, 1, ui) + velint_(2)*derxyz_(2, vi)*conv_r_(1, 2, ui)) ;
    estif_(dWDF, UDF) += ttimetauM*(conv_c_(vi)*conv_r_(2, 0, ui) + velint_(0)*derxyz_(0, vi)*conv_r_(2, 0, ui) + velint_(1)*derxyz_(0, vi)*conv_r_(2, 1, ui) + velint_(2)*derxyz_(0, vi)*conv_r_(2, 2, ui)) ;
    estif_(dWDF, VDF) += ttimetauM*(conv_c_(vi)*conv_r_(2, 1, ui) + velint_(0)*derxyz_(1, vi)*conv_r_(2, 0, ui) + velint_(1)*derxyz_(1, vi)*conv_r_(2, 1, ui) + velint_(2)*derxyz_(1, vi)*conv_r_(2, 2, ui)) ;
    estif_(dWDF, WDF) += ttimetauM*(conv_c_(ui)*conv_c_(vi) + conv_c_(vi)*conv_r_(2, 2, ui) + velint_(0)*derxyz_(2, vi)*conv_r_(2, 0, ui) + velint_(1)*derxyz_(2, vi)*conv_r_(2, 1, ui) + velint_(2)*derxyz_(2, vi)*conv_r_(2, 2, ui)) ;

    /* Stabilisierung der Konvektion (-L_visc_u) */
    estif_(dUDF, UDF) += -2.0*nu_*ttimetauM*(-(conv_c_(vi)*viscs2_(0, 0, ui)) + funct_(ui)*visc_old_(0)*derxyz_(0, vi)) ;
    estif_(dUDF, VDF) += -2.0*nu_*ttimetauM*(-(conv_c_(vi)*viscs2_(0, 1, ui)) + funct_(ui)*visc_old_(0)*derxyz_(1, vi)) ;
    estif_(dUDF, WDF) += -2.0*nu_*ttimetauM*(-(conv_c_(vi)*viscs2_(0, 2, ui)) + funct_(ui)*visc_old_(0)*derxyz_(2, vi)) ;
    estif_(dVDF, UDF) += -2.0*nu_*ttimetauM*(-(conv_c_(vi)*viscs2_(0, 1, ui)) + funct_(ui)*visc_old_(1)*derxyz_(0, vi)) ;
    estif_(dVDF, VDF) += -2.0*nu_*ttimetauM*(-(conv_c_(vi)*viscs2_(1, 1, ui)) + funct_(ui)*visc_old_(1)*derxyz_(1, vi)) ;
    estif_(dVDF, WDF) += -2.0*nu_*ttimetauM*(-(conv_c_(vi)*viscs2_(1, 2, ui)) + funct_(ui)*visc_old_(1)*derxyz_(2, vi)) ;
    estif_(dWDF, UDF) += -2.0*nu_*ttimetauM*(-(conv_c_(vi)*viscs2_(0, 2, ui)) + funct_(ui)*visc_old_(2)*derxyz_(0, vi)) ;
    estif_(dWDF, VDF) += -2.0*nu_*ttimetauM*(-(conv_c_(vi)*viscs2_(1, 2, ui)) + funct_(ui)*visc_old_(2)*derxyz_(1, vi)) ;
    estif_(dWDF, WDF) += -2.0*nu_*ttimetauM*(-(conv_c_(vi)*viscs2_(2, 2, ui)) + funct_(ui)*visc_old_(2)*derxyz_(2, vi)) ;

    /* Stabilisierung der Konvektion ( L_pres_p) */
    estif_(dUDF, UDF) += ttimetauM*funct_(ui)*gradp_(0)*derxyz_(0, vi) ;
    estif_(dUDF, VDF) += ttimetauM*funct_(ui)*gradp_(0)*derxyz_(1, vi) ;
    estif_(dUDF, WDF) += ttimetauM*funct_(ui)*gradp_(0)*derxyz_(2, vi) ;
    estif_(dUDF, PDF) += ttimetauM*conv_c_(vi)*derxyz_(0, ui) ;
    estif_(dVDF, UDF) += ttimetauM*funct_(ui)*gradp_(1)*derxyz_(0, vi) ;
    estif_(dVDF, VDF) += ttimetauM*funct_(ui)*gradp_(1)*derxyz_(1, vi) ;
    estif_(dVDF, WDF) += ttimetauM*funct_(ui)*gradp_(1)*derxyz_(2, vi) ;
    estif_(dVDF, PDF) += ttimetauM*conv_c_(vi)*derxyz_(1, ui) ;
    estif_(dWDF, UDF) += ttimetauM*funct_(ui)*gradp_(2)*derxyz_(0, vi) ;
    estif_(dWDF, VDF) += ttimetauM*funct_(ui)*gradp_(2)*derxyz_(1, vi) ;
    estif_(dWDF, WDF) += ttimetauM*funct_(ui)*gradp_(2)*derxyz_(2, vi) ;
    estif_(dWDF, PDF) += ttimetauM*conv_c_(vi)*derxyz_(2, ui) ;

    /* Viskosit�tsterm */
    estif_(dUDF, UDF)         += nu_*timefacfac*(2.0*derxyz_(0, ui)*derxyz_(0, vi) + derxyz_(1, ui)*derxyz_(1, vi) + derxyz_(2, ui)*derxyz_(2, vi)) ;
    estif_(dUDF, VDF)     += nu_*timefacfac*derxyz_(0, ui)*derxyz_(1, vi) ;
    estif_(dUDF, WDF)     += nu_*timefacfac*derxyz_(0, ui)*derxyz_(2, vi) ;
    estif_(dVDF, UDF)     += nu_*timefacfac*derxyz_(0, vi)*derxyz_(1, ui) ;
    estif_(dVDF, VDF) += nu_*timefacfac*(derxyz_(0, ui)*derxyz_(0, vi) + 2.0*derxyz_(1, ui)*derxyz_(1, vi) + derxyz_(2, ui)*derxyz_(2, vi)) ;
    estif_(dVDF, WDF) += nu_*timefacfac*derxyz_(1, ui)*derxyz_(2, vi) ;
    estif_(dWDF, UDF)     += nu_*timefacfac*derxyz_(0, vi)*derxyz_(2, ui) ;
    estif_(dWDF, VDF) += nu_*timefacfac*derxyz_(1, vi)*derxyz_(2, ui) ;
    estif_(dWDF, WDF) += nu_*timefacfac*(derxyz_(0, ui)*derxyz_(0, vi) + derxyz_(1, ui)*derxyz_(1, vi) + 2.0*derxyz_(2, ui)*derxyz_(2, vi)) ;

    /* Stabilisierung der Viskosit�t ( L_conv_u) */
    estif_(dUDF, UDF)         += 2.0*nu_*ttimetauMp*(conv_c_(ui)*viscs2_(0, 0, vi) + viscs2_(0, 0, vi)*conv_r_(0, 0, ui) + viscs2_(0, 1, vi)*conv_r_(1, 0, ui) + viscs2_(0, 2, vi)*conv_r_(2, 0, ui)) ;
    estif_(dUDF, VDF)     += 2.0*nu_*ttimetauMp*(conv_c_(ui)*viscs2_(0, 1, vi) + viscs2_(0, 0, vi)*conv_r_(0, 1, ui) + viscs2_(0, 1, vi)*conv_r_(1, 1, ui) + viscs2_(0, 2, vi)*conv_r_(2, 1, ui)) ;
    estif_(dUDF, WDF)     += 2.0*nu_*ttimetauMp*(conv_c_(ui)*viscs2_(0, 2, vi) + viscs2_(0, 0, vi)*conv_r_(0, 2, ui) + viscs2_(0, 1, vi)*conv_r_(1, 2, ui) + viscs2_(0, 2, vi)*conv_r_(2, 2, ui)) ;
    estif_(dVDF, UDF)     += 2.0*nu_*ttimetauMp*(conv_c_(ui)*viscs2_(0, 1, vi) + viscs2_(0, 1, vi)*conv_r_(0, 0, ui) + viscs2_(1, 1, vi)*conv_r_(1, 0, ui) + viscs2_(1, 2, vi)*conv_r_(2, 0, ui)) ;
    estif_(dVDF, VDF) += 2.0*nu_*ttimetauMp*(conv_c_(ui)*viscs2_(1, 1, vi) + viscs2_(0, 1, vi)*conv_r_(0, 1, ui) + viscs2_(1, 1, vi)*conv_r_(1, 1, ui) + viscs2_(1, 2, vi)*conv_r_(2, 1, ui)) ;
    estif_(dVDF, WDF) += 2.0*nu_*ttimetauMp*(conv_c_(ui)*viscs2_(1, 2, vi) + viscs2_(0, 1, vi)*conv_r_(0, 2, ui) + viscs2_(1, 1, vi)*conv_r_(1, 2, ui) + viscs2_(1, 2, vi)*conv_r_(2, 2, ui)) ;
    estif_(dWDF, UDF)     += 2.0*nu_*ttimetauMp*(conv_c_(ui)*viscs2_(0, 2, vi) + viscs2_(0, 2, vi)*conv_r_(0, 0, ui) + viscs2_(1, 2, vi)*conv_r_(1, 0, ui) + viscs2_(2, 2, vi)*conv_r_(2, 0, ui)) ;
    estif_(dWDF, VDF) += 2.0*nu_*ttimetauMp*(conv_c_(ui)*viscs2_(1, 2, vi) + viscs2_(0, 2, vi)*conv_r_(0, 1, ui) + viscs2_(1, 2, vi)*conv_r_(1, 1, ui) + viscs2_(2, 2, vi)*conv_r_(2, 1, ui)) ;
    estif_(dWDF, WDF) += 2.0*nu_*ttimetauMp*(conv_c_(ui)*viscs2_(2, 2, vi) + viscs2_(0, 2, vi)*conv_r_(0, 2, ui) + viscs2_(1, 2, vi)*conv_r_(1, 2, ui) + viscs2_(2, 2, vi)*conv_r_(2, 2, ui)) ;

    /* Stabilisierung der Viskosit�t (-L_visc_u) */
    estif_(dUDF, UDF)         += 4.0*(nu_*nu_)*ttimetauMp*(viscs2_(0, 0, ui)*viscs2_(0, 0, vi) + viscs2_(0, 1, ui)*viscs2_(0, 1, vi) + viscs2_(0, 2, ui)*viscs2_(0, 2, vi)) ;
    estif_(dUDF, VDF)     += 4.0*(nu_*nu_)*ttimetauMp*(viscs2_(0, 0, vi)*viscs2_(0, 1, ui) + viscs2_(0, 1, vi)*viscs2_(1, 1, ui) + viscs2_(0, 2, vi)*viscs2_(1, 2, ui)) ;
    estif_(dUDF, WDF)     += 4.0*(nu_*nu_)*ttimetauMp*(viscs2_(0, 0, vi)*viscs2_(0, 2, ui) + viscs2_(0, 1, vi)*viscs2_(1, 2, ui) + viscs2_(0, 2, vi)*viscs2_(2, 2, ui)) ;
    estif_(dVDF, UDF)     += 4.0*(nu_*nu_)*ttimetauMp*(viscs2_(0, 0, ui)*viscs2_(0, 1, vi) + viscs2_(0, 1, ui)*viscs2_(1, 1, vi) + viscs2_(0, 2, ui)*viscs2_(1, 2, vi)) ;
    estif_(dVDF, VDF) += 4.0*(nu_*nu_)*ttimetauMp*(viscs2_(0, 1, ui)*viscs2_(0, 1, vi) + viscs2_(1, 1, ui)*viscs2_(1, 1, vi) + viscs2_(1, 2, ui)*viscs2_(1, 2, vi)) ;
    estif_(dVDF, WDF) += 4.0*(nu_*nu_)*ttimetauMp*(viscs2_(0, 1, vi)*viscs2_(0, 2, ui) + viscs2_(1, 1, vi)*viscs2_(1, 2, ui) + viscs2_(1, 2, vi)*viscs2_(2, 2, ui)) ;
    estif_(dWDF, UDF)     += 4.0*(nu_*nu_)*ttimetauMp*(viscs2_(0, 0, ui)*viscs2_(0, 2, vi) + viscs2_(0, 1, ui)*viscs2_(1, 2, vi) + viscs2_(0, 2, ui)*viscs2_(2, 2, vi)) ;
    estif_(dWDF, VDF) += 4.0*(nu_*nu_)*ttimetauMp*(viscs2_(0, 1, ui)*viscs2_(0, 2, vi) + viscs2_(1, 1, ui)*viscs2_(1, 2, vi) + viscs2_(1, 2, ui)*viscs2_(2, 2, vi)) ;
    estif_(dWDF, WDF) += 4.0*(nu_*nu_)*ttimetauMp*(viscs2_(0, 2, ui)*viscs2_(0, 2, vi) + viscs2_(1, 2, ui)*viscs2_(1, 2, vi) + viscs2_(2, 2, ui)*viscs2_(2, 2, vi)) ;

    /* Stabilisierung der Viskosit�t ( L_pres_p) */
    estif_(dUDF, PDF)     += 2.0*nu_*ttimetauMp*(derxyz_(0, ui)*viscs2_(0, 0, vi) + derxyz_(1, ui)*viscs2_(0, 1, vi) + derxyz_(2, ui)*viscs2_(0, 2, vi)) ;
    estif_(dVDF, PDF) += 2.0*nu_*ttimetauMp*(derxyz_(0, ui)*viscs2_(0, 1, vi) + derxyz_(1, ui)*viscs2_(1, 1, vi) + derxyz_(2, ui)*viscs2_(1, 2, vi)) ;
    estif_(dWDF, PDF) += 2.0*nu_*ttimetauMp*(derxyz_(0, ui)*viscs2_(0, 2, vi) + derxyz_(1, ui)*viscs2_(1, 2, vi) + derxyz_(2, ui)*viscs2_(2, 2, vi)) ;

    /* Druckterm */
    estif_(dUDF, PDF)     += -(timefacfac*funct_(ui)*derxyz_(0, vi)) ;
    estif_(dVDF, PDF) += -(timefacfac*funct_(ui)*derxyz_(1, vi)) ;
    estif_(dWDF, PDF) += -(timefacfac*funct_(ui)*derxyz_(2, vi)) ;

    /* Stabilisierung des Drucks ( L_conv_u) */
    estif_(dPDF, UDF)     += ttimetauMp*(conv_c_(ui)*derxyz_(0, vi) + derxyz_(0, vi)*conv_r_(0, 0, ui) + derxyz_(1, vi)*conv_r_(1, 0, ui) + derxyz_(2, vi)*conv_r_(2, 0, ui)) ;
    estif_(dPDF, VDF) += ttimetauMp*(conv_c_(ui)*derxyz_(1, vi) + derxyz_(0, vi)*conv_r_(0, 1, ui) + derxyz_(1, vi)*conv_r_(1, 1, ui) + derxyz_(2, vi)*conv_r_(2, 1, ui)) ;
    estif_(dPDF, WDF) += ttimetauMp*(conv_c_(ui)*derxyz_(2, vi) + derxyz_(0, vi)*conv_r_(0, 2, ui) + derxyz_(1, vi)*conv_r_(1, 2, ui) + derxyz_(2, vi)*conv_r_(2, 2, ui)) ;

    /* Stabilisierung des Drucks (-L_visc_u) */
    estif_(dPDF, UDF)     += 2.0*nu_*ttimetauMp*(derxyz_(0, vi)*viscs2_(0, 0, ui) + derxyz_(1, vi)*viscs2_(0, 1, ui) + derxyz_(2, vi)*viscs2_(0, 2, ui)) ;
    estif_(dPDF, VDF) += 2.0*nu_*ttimetauMp*(derxyz_(0, vi)*viscs2_(0, 1, ui) + derxyz_(1, vi)*viscs2_(1, 1, ui) + derxyz_(2, vi)*viscs2_(1, 2, ui)) ;
    estif_(dPDF, WDF) += 2.0*nu_*ttimetauMp*(derxyz_(0, vi)*viscs2_(0, 2, ui) + derxyz_(1, vi)*viscs2_(1, 2, ui) + derxyz_(2, vi)*viscs2_(2, 2, ui)) ;

    /* Stabilisierung des Drucks ( L_pres_p) */
    estif_(dPDF, PDF) += ttimetauMp*(derxyz_(0, ui)*derxyz_(0, vi) + derxyz_(1, ui)*derxyz_(1, vi) + derxyz_(2, ui)*derxyz_(2, vi)) ;

    /* Divergenzfreiheit */
    estif_(dPDF, UDF) += timefacfac*funct_(vi)*derxyz_(0, ui) ;
    estif_(dPDF, VDF) += timefacfac*funct_(vi)*derxyz_(1, ui) ;
    estif_(dPDF, WDF) += timefacfac*funct_(vi)*derxyz_(2, ui) ;

    /* Kontinuit�tsstabilisierung */
    estif_(dUDF, UDF) += (thsl*thsl)*tau_C*derxyz_(0, ui)*derxyz_(0, vi) ;
    estif_(dUDF, VDF) += (thsl*thsl)*tau_C*derxyz_(0, vi)*derxyz_(1, ui) ;
    estif_(dUDF, WDF) += (thsl*thsl)*tau_C*derxyz_(0, vi)*derxyz_(2, ui) ;
    estif_(dVDF, UDF) += (thsl*thsl)*tau_C*derxyz_(0, ui)*derxyz_(1, vi) ;
    estif_(dVDF, VDF) += (thsl*thsl)*tau_C*derxyz_(1, ui)*derxyz_(1, vi) ;
    estif_(dVDF, WDF) += (thsl*thsl)*tau_C*derxyz_(1, vi)*derxyz_(2, ui) ;
    estif_(dWDF, UDF) += (thsl*thsl)*tau_C*derxyz_(0, ui)*derxyz_(2, vi) ;
    estif_(dWDF, VDF) += (thsl*thsl)*tau_C*derxyz_(1, ui)*derxyz_(2, vi) ;
    estif_(dWDF, WDF) += (thsl*thsl)*tau_C*derxyz_(2, ui)*derxyz_(2, vi) ;

    /* Massenterm */
    estif_(dUDF, UDF)         += fac*funct_(ui)*funct_(vi) ;
    estif_(dVDF, VDF) += fac*funct_(ui)*funct_(vi) ;
    estif_(dWDF, WDF) += fac*funct_(ui)*funct_(vi) ;

    /* Konvektionsstabilisierung */
    estif_(dUDF, UDF)         += timetauM*funct_(ui)*(2.0*velint_(0)*derxyz_(0, vi) + velint_(1)*derxyz_(1, vi) + velint_(2)*derxyz_(2, vi)) ;
    estif_(dUDF, VDF)     += timetauM*funct_(ui)*velint_(0)*derxyz_(1, vi) ;
    estif_(dUDF, WDF)     += timetauM*funct_(ui)*velint_(0)*derxyz_(2, vi) ;
    estif_(dVDF, UDF)     += timetauM*funct_(ui)*velint_(1)*derxyz_(0, vi) ;
    estif_(dVDF, VDF) += timetauM*funct_(ui)*(velint_(0)*derxyz_(0, vi) + 2.0*velint_(1)*derxyz_(1, vi) + velint_(2)*derxyz_(2, vi)) ;
    estif_(dVDF, WDF) += timetauM*funct_(ui)*velint_(1)*derxyz_(2, vi) ;
    estif_(dWDF, UDF)     += timetauM*funct_(ui)*velint_(2)*derxyz_(0, vi) ;
    estif_(dWDF, VDF) += timetauM*funct_(ui)*velint_(2)*derxyz_(1, vi) ;
    estif_(dWDF, WDF) += timetauM*funct_(ui)*(velint_(0)*derxyz_(0, vi) + velint_(1)*derxyz_(1, vi) + 2.0*velint_(2)*derxyz_(2, vi)) ;

    /* Viskosit�tsstabilisierung */
    estif_(dUDF, UDF)         += 2.0*nu_*timetauMp*funct_(ui)*viscs2_(0, 0, vi) ;
    estif_(dUDF, VDF)     += 2.0*nu_*timetauMp*funct_(ui)*viscs2_(0, 1, vi) ;
    estif_(dUDF, WDF)     += 2.0*nu_*timetauMp*funct_(ui)*viscs2_(0, 2, vi) ;
    estif_(dVDF, UDF)     += 2.0*nu_*timetauMp*funct_(ui)*viscs2_(0, 1, vi) ;
    estif_(dVDF, VDF) += 2.0*nu_*timetauMp*funct_(ui)*viscs2_(1, 1, vi) ;
    estif_(dVDF, WDF) += 2.0*nu_*timetauMp*funct_(ui)*viscs2_(1, 2, vi) ;
    estif_(dWDF, UDF)     += 2.0*nu_*timetauMp*funct_(ui)*viscs2_(0, 2, vi) ;
    estif_(dWDF, VDF) += 2.0*nu_*timetauMp*funct_(ui)*viscs2_(1, 2, vi) ;
    estif_(dWDF, WDF) += 2.0*nu_*timetauMp*funct_(ui)*viscs2_(2, 2, vi) ;

    /* Stabilisierung der Druckgleichung */
    estif_(dPDF, UDF)     += timetauMp*funct_(ui)*derxyz_(0, vi) ;
    estif_(dPDF, VDF) += timetauMp*funct_(ui)*derxyz_(1, vi) ;
    estif_(dPDF, WDF) += timetauMp*funct_(ui)*derxyz_(2, vi) ;

    /* Quellterm der rechten Seite */

    /* Konvektionsstabilisierung */
    estif_(dUDF, UDF)         += -(timetauM*funct_(ui)*rhsint_(0)*derxyz_(0, vi)) ;
    estif_(dUDF, VDF)     += -(timetauM*funct_(ui)*rhsint_(0)*derxyz_(1, vi)) ;
    estif_(dUDF, WDF)     += -(timetauM*funct_(ui)*rhsint_(0)*derxyz_(2, vi)) ;
    estif_(dVDF, UDF)     += -(timetauM*funct_(ui)*rhsint_(1)*derxyz_(0, vi)) ;
    estif_(dVDF, VDF) += -(timetauM*funct_(ui)*rhsint_(1)*derxyz_(1, vi)) ;
    estif_(dVDF, WDF) += -(timetauM*funct_(ui)*rhsint_(1)*derxyz_(2, vi)) ;
    estif_(dWDF, UDF)     += -(timetauM*funct_(ui)*rhsint_(2)*derxyz_(0, vi)) ;
    estif_(dWDF, VDF) += -(timetauM*funct_(ui)*rhsint_(2)*derxyz_(1, vi)) ;
    estif_(dWDF, WDF) += -(timetauM*funct_(ui)*rhsint_(2)*derxyz_(2, vi)) ;

    /* Viskosit�tsstabilisierung */

    /* Stabilisierung der Druckgleichung */

  }
}
