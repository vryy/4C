// This file is currently not used by any of the fluid implementations
// Fluid3Impl, Fluid3GenalphaResVMM or Fluid3Stationary. They use
// rearranged versions of this code (slightly adopted to the different
// approaches) --- they are expected to be faster and more readable.
// 
// Axel is still referencing these loops in his xfem routines, so they
// should not be deleted until he declares them obsolete...

for (int ui=0; ui<iel; ++ui)    // exchanged ui- and vi-loop => slight speedup!?   g.bau 03/07
{
  for (int vi=0; vi<iel; ++vi)
  {

    /* Konvektionsterm */
#if 0
    /* 

                 /                                              \
                |  / n+1       \        /          \   n+1       |
                | | u   o nabla | Du + | Du o nabla | u     , v  |
                |  \ (i)       /        \          /   (i)       |
                 \                                              /
    */
    estif_(vi*4, ui*4)         += timefacfac*funct_(vi)*(conv_c_(ui) + conv_r_(0, 0, ui)) ;
    estif_(vi*4, ui*4 + 1)     += timefacfac*funct_(vi)*conv_r_(0, 1, ui) ;
    estif_(vi*4, ui*4 + 2)     += timefacfac*funct_(vi)*conv_r_(0, 2, ui) ;
    estif_(vi*4 + 1, ui*4)     += timefacfac*funct_(vi)*conv_r_(1, 0, ui) ;
    estif_(vi*4 + 1, ui*4 + 1) += timefacfac*funct_(vi)*(conv_c_(ui) + conv_r_(1, 1, ui)) ;
    estif_(vi*4 + 1, ui*4 + 2) += timefacfac*funct_(vi)*conv_r_(1, 2, ui) ;
    estif_(vi*4 + 2, ui*4)     += timefacfac*funct_(vi)*conv_r_(2, 0, ui) ;
    estif_(vi*4 + 2, ui*4 + 1) += timefacfac*funct_(vi)*conv_r_(2, 1, ui) ;
    estif_(vi*4 + 2, ui*4 + 2) += timefacfac*funct_(vi)*(conv_c_(ui) + conv_r_(2, 2, ui)) ;
#else
    /* 
         
                          /                       \
                         |  / n+1       \          |
                         | | u   o nabla | Du , v  |
                         |  \ (i)       /          |
                          \                       /
    */   
    estif_(vi*4, ui*4)         += timefacfac*funct_(vi)*conv_c_(ui) ;
    estif_(vi*4 + 1, ui*4 + 1) += timefacfac*funct_(vi)*conv_c_(ui) ;
    estif_(vi*4 + 2, ui*4 + 2) += timefacfac*funct_(vi)*conv_c_(ui) ;
#endif    
    /* Stabilisierung der Konvektion ( L_conv_u) */
#if 0
    /*
                /                                           \
               |    / n+1        \   n+1    /          \     |
               |   | u    o nabla | u    , | Du o nabla | v  |
               |    \ (i)        /   (i)    \          /     |
                \                                           /
   
                /                                           \
               |    / n+1        \        / n+1        \     |
               |   | u    o nabla | Du , | u    o nabla | v  |
               |    \ (i)        /        \ (i)        /     |
                \                                           /
   
                /                                           \
               |    /          \   n+1    / n+1        \     |
               |   | Du o nabla | u    , | u    o nabla | v  |
               |    \          /   (i)    \ (i)        /     |
                \                                           /
    */ 
    estif_(vi*4, ui*4)         += ttimetauM*(conv_c_(ui)*conv_c_(vi) + conv_c_(vi)*conv_r_(0, 0, ui) + velint_(0)*derxyz_(0, vi)*conv_r_(0, 0, ui) + velint_(1)*derxyz_(0, vi)*conv_r_(0, 1, ui) + velint_(2)*derxyz_(0, vi)*conv_r_(0, 2, ui)) ;
    estif_(vi*4, ui*4 + 1)     += ttimetauM*(conv_c_(vi)*conv_r_(0, 1, ui) + velint_(0)*derxyz_(1, vi)*conv_r_(0, 0, ui) + velint_(1)*derxyz_(1, vi)*conv_r_(0, 1, ui) + velint_(2)*derxyz_(1, vi)*conv_r_(0, 2, ui)) ;
    estif_(vi*4, ui*4 + 2)     += ttimetauM*(conv_c_(vi)*conv_r_(0, 2, ui) + velint_(0)*derxyz_(2, vi)*conv_r_(0, 0, ui) + velint_(1)*derxyz_(2, vi)*conv_r_(0, 1, ui) + velint_(2)*derxyz_(2, vi)*conv_r_(0, 2, ui)) ;
    estif_(vi*4 + 1, ui*4)     += ttimetauM*(conv_c_(vi)*conv_r_(1, 0, ui) + velint_(0)*derxyz_(0, vi)*conv_r_(1, 0, ui) + velint_(1)*derxyz_(0, vi)*conv_r_(1, 1, ui) + velint_(2)*derxyz_(0, vi)*conv_r_(1, 2, ui)) ;
    estif_(vi*4 + 1, ui*4 + 1) += ttimetauM*(conv_c_(ui)*conv_c_(vi) + conv_c_(vi)*conv_r_(1, 1, ui) + velint_(0)*derxyz_(1, vi)*conv_r_(1, 0, ui) + velint_(1)*derxyz_(1, vi)*conv_r_(1, 1, ui) + velint_(2)*derxyz_(1, vi)*conv_r_(1, 2, ui)) ;
    estif_(vi*4 + 1, ui*4 + 2) += ttimetauM*(conv_c_(vi)*conv_r_(1, 2, ui) + velint_(0)*derxyz_(2, vi)*conv_r_(1, 0, ui) + velint_(1)*derxyz_(2, vi)*conv_r_(1, 1, ui) + velint_(2)*derxyz_(2, vi)*conv_r_(1, 2, ui)) ;
    estif_(vi*4 + 2, ui*4)     += ttimetauM*(conv_c_(vi)*conv_r_(2, 0, ui) + velint_(0)*derxyz_(0, vi)*conv_r_(2, 0, ui) + velint_(1)*derxyz_(0, vi)*conv_r_(2, 1, ui) + velint_(2)*derxyz_(0, vi)*conv_r_(2, 2, ui)) ;
    estif_(vi*4 + 2, ui*4 + 1) += ttimetauM*(conv_c_(vi)*conv_r_(2, 1, ui) + velint_(0)*derxyz_(1, vi)*conv_r_(2, 0, ui) + velint_(1)*derxyz_(1, vi)*conv_r_(2, 1, ui) + velint_(2)*derxyz_(1, vi)*conv_r_(2, 2, ui)) ;
    estif_(vi*4 + 2, ui*4 + 2) += ttimetauM*(conv_c_(ui)*conv_c_(vi) + conv_c_(vi)*conv_r_(2, 2, ui) + velint_(0)*derxyz_(2, vi)*conv_r_(2, 0, ui) + velint_(1)*derxyz_(2, vi)*conv_r_(2, 1, ui) + velint_(2)*derxyz_(2, vi)*conv_r_(2, 2, ui)) ;
#else
    /*
                    /                                         \
                   |    / n+1       \        / n+1       \     |
                   |   | u   o nabla | Du , | u   o nabla | v  |
                   |    \ (i)       /        \ (i)       /     |
                    \                                         /

    */ 
    estif_(vi*4, ui*4)         += ttimetauM*conv_c_(ui)*conv_c_(vi) ;
    estif_(vi*4 + 1, ui*4 + 1) += ttimetauM*conv_c_(ui)*conv_c_(vi) ;
    estif_(vi*4 + 2, ui*4 + 2) += ttimetauM*conv_c_(ui)*conv_c_(vi) ;
#endif
    /* Stabilisierung der Konvektion (-L_visc_u) */
#if 0
    /*   
                    /                                        \
                   |               / n+1 \    /          \     |
                   |  nabla o eps | u     |, | Du o nabla | v  |
                   |               \ (i) /    \          /     |
                    \                                         /
   
   
                    /                                        \
                   |               /  \    / n+1        \     |
                   |  nabla o eps | Du |, | u    o nabla | v  |
                   |               \  /    \ (i)        /     |
                    \                                        /
    */
    estif_(vi*4, ui*4)         += -2.0*nu_*ttimetauM*(-(conv_c_(vi)*viscs2_(0, 0, ui)) + funct_(ui)*visc_old_(0)*derxyz_(0, vi)) ;
    estif_(vi*4, ui*4 + 1)     += -2.0*nu_*ttimetauM*(-(conv_c_(vi)*viscs2_(0, 1, ui)) + funct_(ui)*visc_old_(0)*derxyz_(1, vi)) ;
    estif_(vi*4, ui*4 + 2)     += -2.0*nu_*ttimetauM*(-(conv_c_(vi)*viscs2_(0, 2, ui)) + funct_(ui)*visc_old_(0)*derxyz_(2, vi)) ;
    estif_(vi*4 + 1, ui*4)     += -2.0*nu_*ttimetauM*(-(conv_c_(vi)*viscs2_(0, 1, ui)) + funct_(ui)*visc_old_(1)*derxyz_(0, vi)) ;
    estif_(vi*4 + 1, ui*4 + 1) += -2.0*nu_*ttimetauM*(-(conv_c_(vi)*viscs2_(1, 1, ui)) + funct_(ui)*visc_old_(1)*derxyz_(1, vi)) ;
    estif_(vi*4 + 1, ui*4 + 2) += -2.0*nu_*ttimetauM*(-(conv_c_(vi)*viscs2_(1, 2, ui)) + funct_(ui)*visc_old_(1)*derxyz_(2, vi)) ;
    estif_(vi*4 + 2, ui*4)     += -2.0*nu_*ttimetauM*(-(conv_c_(vi)*viscs2_(0, 2, ui)) + funct_(ui)*visc_old_(2)*derxyz_(0, vi)) ;
    estif_(vi*4 + 2, ui*4 + 1) += -2.0*nu_*ttimetauM*(-(conv_c_(vi)*viscs2_(1, 2, ui)) + funct_(ui)*visc_old_(2)*derxyz_(1, vi)) ;
    estif_(vi*4 + 2, ui*4 + 2) += -2.0*nu_*ttimetauM*(-(conv_c_(vi)*viscs2_(2, 2, ui)) + funct_(ui)*visc_old_(2)*derxyz_(2, vi)) ;
#else
    /*
                  /                                        \
                 |               /  \    / n+1        \     |
                 |  nabla o eps | Du |, | u    o nabla | v  |
                 |               \  /    \ (i)        /     |
                  \                                        /
    */
    estif_(vi*4, ui*4)         += 2.0*nu_*ttimetauM*conv_c_(vi)*viscs2_(0, 0, ui) ;
    estif_(vi*4, ui*4 + 1)     += 2.0*nu_*ttimetauM*conv_c_(vi)*viscs2_(0, 1, ui) ;
    estif_(vi*4, ui*4 + 2)     += 2.0*nu_*ttimetauM*conv_c_(vi)*viscs2_(0, 2, ui) ;
    estif_(vi*4 + 1, ui*4)     += 2.0*nu_*ttimetauM*conv_c_(vi)*viscs2_(0, 1, ui) ;
    estif_(vi*4 + 1, ui*4 + 1) += 2.0*nu_*ttimetauM*conv_c_(vi)*viscs2_(1, 1, ui) ;
    estif_(vi*4 + 1, ui*4 + 2) += 2.0*nu_*ttimetauM*conv_c_(vi)*viscs2_(1, 2, ui) ;
    estif_(vi*4 + 2, ui*4)     += 2.0*nu_*ttimetauM*conv_c_(vi)*viscs2_(0, 2, ui) ;
    estif_(vi*4 + 2, ui*4 + 1) += 2.0*nu_*ttimetauM*conv_c_(vi)*viscs2_(1, 2, ui) ;
    estif_(vi*4 + 2, ui*4 + 2) += 2.0*nu_*ttimetauM*conv_c_(vi)*viscs2_(2, 2, ui) ;
#endif  
    /* Stabilisierung der Konvektion ( L_pres_p) */
#if 0
    /*
   
                    /                               \
                   |         n+1    /          \     |
                   |  nabla p    , | Du o nabla | v  |
                   |         (i)    \          /     |
                    \                               /
 
   
                    /                              \
                   |              / n+1       \     |
                   |  nabla Dp , | u   o nabla | v  |
                   |              \ (i)       /     |
                    \                              /
    */ 
    estif_(vi*4, ui*4)         += ttimetauM*funct_(ui)*gradp_(0)*derxyz_(0, vi) ;
    estif_(vi*4, ui*4 + 1)     += ttimetauM*funct_(ui)*gradp_(0)*derxyz_(1, vi) ;
    estif_(vi*4, ui*4 + 2)     += ttimetauM*funct_(ui)*gradp_(0)*derxyz_(2, vi) ;
    estif_(vi*4, ui*4 + 3)     += ttimetauM*conv_c_(vi)*derxyz_(0, ui) ;
    estif_(vi*4 + 1, ui*4)     += ttimetauM*funct_(ui)*gradp_(1)*derxyz_(0, vi) ;
    estif_(vi*4 + 1, ui*4 + 1) += ttimetauM*funct_(ui)*gradp_(1)*derxyz_(1, vi) ;
    estif_(vi*4 + 1, ui*4 + 2) += ttimetauM*funct_(ui)*gradp_(1)*derxyz_(2, vi) ;
    estif_(vi*4 + 1, ui*4 + 3) += ttimetauM*conv_c_(vi)*derxyz_(1, ui) ;
    estif_(vi*4 + 2, ui*4)     += ttimetauM*funct_(ui)*gradp_(2)*derxyz_(0, vi) ;
    estif_(vi*4 + 2, ui*4 + 1) += ttimetauM*funct_(ui)*gradp_(2)*derxyz_(1, vi) ;
    estif_(vi*4 + 2, ui*4 + 2) += ttimetauM*funct_(ui)*gradp_(2)*derxyz_(2, vi) ;
    estif_(vi*4 + 2, ui*4 + 3) += ttimetauM*conv_c_(vi)*derxyz_(2, ui) ;
#else
    /*
   
                     /                              \
                    |              / n+1       \     |
                    |  nabla Dp , | u   o nabla | v  |
                    |              \ (i)       /     |
                     \                              /
    */ 
    estif_(vi*4, ui*4 + 3)     += ttimetauM*conv_c_(vi)*derxyz_(0, ui) ;
    estif_(vi*4 + 1, ui*4 + 3) += ttimetauM*conv_c_(vi)*derxyz_(1, ui) ;
    estif_(vi*4 + 2, ui*4 + 3) += ttimetauM*conv_c_(vi)*derxyz_(2, ui) ;
#endif
    /* Viskositätsterm */
    /*
                         /                        \
                        |       /  \         / \   |
                        |  eps | Du | , eps | v |  |
                        |       \  /         \ /   |
                         \                        /
    */
    estif_(vi*4, ui*4)         += nu_*timefacfac*(2.0*derxyz_(0, ui)*derxyz_(0, vi) + derxyz_(1, ui)*derxyz_(1, vi) + derxyz_(2, ui)*derxyz_(2, vi)) ;
    estif_(vi*4, ui*4 + 1)     += nu_*timefacfac*derxyz_(0, ui)*derxyz_(1, vi) ;
    estif_(vi*4, ui*4 + 2)     += nu_*timefacfac*derxyz_(0, ui)*derxyz_(2, vi) ;
    estif_(vi*4 + 1, ui*4)     += nu_*timefacfac*derxyz_(0, vi)*derxyz_(1, ui) ;
    estif_(vi*4 + 1, ui*4 + 1) += nu_*timefacfac*(derxyz_(0, ui)*derxyz_(0, vi) + 2.0*derxyz_(1, ui)*derxyz_(1, vi) + derxyz_(2, ui)*derxyz_(2, vi)) ;
    estif_(vi*4 + 1, ui*4 + 2) += nu_*timefacfac*derxyz_(1, ui)*derxyz_(2, vi) ;
    estif_(vi*4 + 2, ui*4)     += nu_*timefacfac*derxyz_(0, vi)*derxyz_(2, ui) ;
    estif_(vi*4 + 2, ui*4 + 1) += nu_*timefacfac*derxyz_(1, vi)*derxyz_(2, ui) ;
    estif_(vi*4 + 2, ui*4 + 2) += nu_*timefacfac*(derxyz_(0, ui)*derxyz_(0, vi) + derxyz_(1, ui)*derxyz_(1, vi) + 2.0*derxyz_(2, ui)*derxyz_(2, vi)) ;

    /* Stabilisierung der Viskosität ( L_conv_u) */
#if 0
    /*    
            /                                                      \
           |  / n+1       \        /          \   n+1               |
           | | u   o nabla | Du + | Du o nabla | u    , div eps (v) |
           |  \ (i)       /        \          /   (i)               |
            \                                                      /
    */  
    estif_(vi*4, ui*4)         += 2.0*nu_*ttimetauMp*(conv_c_(ui)*viscs2_(0, 0, vi) + viscs2_(0, 0, vi)*conv_r_(0, 0, ui) + viscs2_(0, 1, vi)*conv_r_(1, 0, ui) + viscs2_(0, 2, vi)*conv_r_(2, 0, ui)) ;
    estif_(vi*4, ui*4 + 1)     += 2.0*nu_*ttimetauMp*(conv_c_(ui)*viscs2_(0, 1, vi) + viscs2_(0, 0, vi)*conv_r_(0, 1, ui) + viscs2_(0, 1, vi)*conv_r_(1, 1, ui) + viscs2_(0, 2, vi)*conv_r_(2, 1, ui)) ;
    estif_(vi*4, ui*4 + 2)     += 2.0*nu_*ttimetauMp*(conv_c_(ui)*viscs2_(0, 2, vi) + viscs2_(0, 0, vi)*conv_r_(0, 2, ui) + viscs2_(0, 1, vi)*conv_r_(1, 2, ui) + viscs2_(0, 2, vi)*conv_r_(2, 2, ui)) ;
    estif_(vi*4 + 1, ui*4)     += 2.0*nu_*ttimetauMp*(conv_c_(ui)*viscs2_(0, 1, vi) + viscs2_(0, 1, vi)*conv_r_(0, 0, ui) + viscs2_(1, 1, vi)*conv_r_(1, 0, ui) + viscs2_(1, 2, vi)*conv_r_(2, 0, ui)) ;
    estif_(vi*4 + 1, ui*4 + 1) += 2.0*nu_*ttimetauMp*(conv_c_(ui)*viscs2_(1, 1, vi) + viscs2_(0, 1, vi)*conv_r_(0, 1, ui) + viscs2_(1, 1, vi)*conv_r_(1, 1, ui) + viscs2_(1, 2, vi)*conv_r_(2, 1, ui)) ;
    estif_(vi*4 + 1, ui*4 + 2) += 2.0*nu_*ttimetauMp*(conv_c_(ui)*viscs2_(1, 2, vi) + viscs2_(0, 1, vi)*conv_r_(0, 2, ui) + viscs2_(1, 1, vi)*conv_r_(1, 2, ui) + viscs2_(1, 2, vi)*conv_r_(2, 2, ui)) ;
    estif_(vi*4 + 2, ui*4)     += 2.0*nu_*ttimetauMp*(conv_c_(ui)*viscs2_(0, 2, vi) + viscs2_(0, 2, vi)*conv_r_(0, 0, ui) + viscs2_(1, 2, vi)*conv_r_(1, 0, ui) + viscs2_(2, 2, vi)*conv_r_(2, 0, ui)) ;
    estif_(vi*4 + 2, ui*4 + 1) += 2.0*nu_*ttimetauMp*(conv_c_(ui)*viscs2_(1, 2, vi) + viscs2_(0, 2, vi)*conv_r_(0, 1, ui) + viscs2_(1, 2, vi)*conv_r_(1, 1, ui) + viscs2_(2, 2, vi)*conv_r_(2, 1, ui)) ;
    estif_(vi*4 + 2, ui*4 + 2) += 2.0*nu_*ttimetauMp*(conv_c_(ui)*viscs2_(2, 2, vi) + viscs2_(0, 2, vi)*conv_r_(0, 2, ui) + viscs2_(1, 2, vi)*conv_r_(1, 2, ui) + viscs2_(2, 2, vi)*conv_r_(2, 2, ui)) ;
#else
    /*    
                    /                                \
                   |  / n+1       \                   |
                   | | u   o nabla | Du , div eps (v) |
                   |  \ (i)       /                   |
                    \                                /
    */
    estif_(vi*4, ui*4)         += 2.0*nu_*ttimetauMp*conv_c_(ui)*viscs2_(0, 0, vi) ;
    estif_(vi*4, ui*4 + 1)     += 2.0*nu_*ttimetauMp*conv_c_(ui)*viscs2_(0, 1, vi) ;
    estif_(vi*4, ui*4 + 2)     += 2.0*nu_*ttimetauMp*conv_c_(ui)*viscs2_(0, 2, vi) ;
    estif_(vi*4 + 1, ui*4)     += 2.0*nu_*ttimetauMp*conv_c_(ui)*viscs2_(0, 1, vi) ;
    estif_(vi*4 + 1, ui*4 + 1) += 2.0*nu_*ttimetauMp*conv_c_(ui)*viscs2_(1, 1, vi) ;
    estif_(vi*4 + 1, ui*4 + 2) += 2.0*nu_*ttimetauMp*conv_c_(ui)*viscs2_(1, 2, vi) ;
    estif_(vi*4 + 2, ui*4)     += 2.0*nu_*ttimetauMp*conv_c_(ui)*viscs2_(0, 2, vi) ;
    estif_(vi*4 + 2, ui*4 + 1) += 2.0*nu_*ttimetauMp*conv_c_(ui)*viscs2_(1, 2, vi) ;
    estif_(vi*4 + 2, ui*4 + 2) += 2.0*nu_*ttimetauMp*conv_c_(ui)*viscs2_(2, 2, vi) ;
#endif
    /* Stabilisierung der Viskosität (-L_visc_u) */
    /*   
                    /                                 \
                   |               /  \                |
                   |  nabla o eps | Du | , div eps (v) |
                   |               \  /                |
                    \                                 /
    */
    estif_(vi*4, ui*4)         += 4.0*(nu_*nu_)*ttimetauMp*(viscs2_(0, 0, ui)*viscs2_(0, 0, vi) + viscs2_(0, 1, ui)*viscs2_(0, 1, vi) + viscs2_(0, 2, ui)*viscs2_(0, 2, vi)) ;
    estif_(vi*4, ui*4 + 1)     += 4.0*(nu_*nu_)*ttimetauMp*(viscs2_(0, 0, vi)*viscs2_(0, 1, ui) + viscs2_(0, 1, vi)*viscs2_(1, 1, ui) + viscs2_(0, 2, vi)*viscs2_(1, 2, ui)) ;
    estif_(vi*4, ui*4 + 2)     += 4.0*(nu_*nu_)*ttimetauMp*(viscs2_(0, 0, vi)*viscs2_(0, 2, ui) + viscs2_(0, 1, vi)*viscs2_(1, 2, ui) + viscs2_(0, 2, vi)*viscs2_(2, 2, ui)) ;
    estif_(vi*4 + 1, ui*4)     += 4.0*(nu_*nu_)*ttimetauMp*(viscs2_(0, 0, ui)*viscs2_(0, 1, vi) + viscs2_(0, 1, ui)*viscs2_(1, 1, vi) + viscs2_(0, 2, ui)*viscs2_(1, 2, vi)) ;
    estif_(vi*4 + 1, ui*4 + 1) += 4.0*(nu_*nu_)*ttimetauMp*(viscs2_(0, 1, ui)*viscs2_(0, 1, vi) + viscs2_(1, 1, ui)*viscs2_(1, 1, vi) + viscs2_(1, 2, ui)*viscs2_(1, 2, vi)) ;
    estif_(vi*4 + 1, ui*4 + 2) += 4.0*(nu_*nu_)*ttimetauMp*(viscs2_(0, 1, vi)*viscs2_(0, 2, ui) + viscs2_(1, 1, vi)*viscs2_(1, 2, ui) + viscs2_(1, 2, vi)*viscs2_(2, 2, ui)) ;
    estif_(vi*4 + 2, ui*4)     += 4.0*(nu_*nu_)*ttimetauMp*(viscs2_(0, 0, ui)*viscs2_(0, 2, vi) + viscs2_(0, 1, ui)*viscs2_(1, 2, vi) + viscs2_(0, 2, ui)*viscs2_(2, 2, vi)) ;
    estif_(vi*4 + 2, ui*4 + 1) += 4.0*(nu_*nu_)*ttimetauMp*(viscs2_(0, 1, ui)*viscs2_(0, 2, vi) + viscs2_(1, 1, ui)*viscs2_(1, 2, vi) + viscs2_(1, 2, ui)*viscs2_(2, 2, vi)) ;
    estif_(vi*4 + 2, ui*4 + 2) += 4.0*(nu_*nu_)*ttimetauMp*(viscs2_(0, 2, ui)*viscs2_(0, 2, vi) + viscs2_(1, 2, ui)*viscs2_(1, 2, vi) + viscs2_(2, 2, ui)*viscs2_(2, 2, vi)) ;

    /* Stabilisierung der Viskosität ( L_pres_p) */
    /*    
                       /                        \
                      |                          |
                      |  nabla Dp , div eps (v)  |
                      |                          |
                       \                        /
    */
    estif_(vi*4, ui*4 + 3)     += 2.0*nu_*ttimetauMp*(derxyz_(0, ui)*viscs2_(0, 0, vi) + derxyz_(1, ui)*viscs2_(0, 1, vi) + derxyz_(2, ui)*viscs2_(0, 2, vi)) ;
    estif_(vi*4 + 1, ui*4 + 3) += 2.0*nu_*ttimetauMp*(derxyz_(0, ui)*viscs2_(0, 1, vi) + derxyz_(1, ui)*viscs2_(1, 1, vi) + derxyz_(2, ui)*viscs2_(1, 2, vi)) ;
    estif_(vi*4 + 2, ui*4 + 3) += 2.0*nu_*ttimetauMp*(derxyz_(0, ui)*viscs2_(0, 2, vi) + derxyz_(1, ui)*viscs2_(1, 2, vi) + derxyz_(2, ui)*viscs2_(2, 2, vi)) ;

    /* Druckterm */
    /*  
        
                          /                \
                         |                  |
                         |  Dp , nabla o v  |
                         |                  |
                          \                /
    */

    estif_(vi*4, ui*4 + 3)     += -(timefacfac*funct_(ui)*derxyz_(0, vi)) ;
    estif_(vi*4 + 1, ui*4 + 3) += -(timefacfac*funct_(ui)*derxyz_(1, vi)) ;
    estif_(vi*4 + 2, ui*4 + 3) += -(timefacfac*funct_(ui)*derxyz_(2, vi)) ;

    /* Stabilisierung des Drucks ( L_conv_u) */
#if 0
    /*
                    /                               \
                   |         n+1    /          \     |
                   |  nabla p    , | Du o nabla | v  |
                   |         (i)    \          /     |
                    \                               /


                    /                              \
                   |              / n+1       \     |
                   |  nabla Dp , | u   o nabla | v  |
                   |              \ (i)       /     |
                    \                              /
    */
    estif_(vi*4 + 3, ui*4)     += ttimetauMp*(conv_c_(ui)*derxyz_(0, vi) + derxyz_(0, vi)*conv_r_(0, 0, ui) + derxyz_(1, vi)*conv_r_(1, 0, ui) + derxyz_(2, vi)*conv_r_(2, 0, ui)) ;
    estif_(vi*4 + 3, ui*4 + 1) += ttimetauMp*(conv_c_(ui)*derxyz_(1, vi) + derxyz_(0, vi)*conv_r_(0, 1, ui) + derxyz_(1, vi)*conv_r_(1, 1, ui) + derxyz_(2, vi)*conv_r_(2, 1, ui)) ;
    estif_(vi*4 + 3, ui*4 + 2) += ttimetauMp*(conv_c_(ui)*derxyz_(2, vi) + derxyz_(0, vi)*conv_r_(0, 2, ui) + derxyz_(1, vi)*conv_r_(1, 2, ui) + derxyz_(2, vi)*conv_r_(2, 2, ui)) ;
#else
    /*    
                    /                              \
                   |              / n+1       \     |
                   |  nabla Dp , | u   o nabla | v  |
                   |              \ (i)       /     |
                    \                              /
    */
    estif_(vi*4 + 3, ui*4)     += ttimetauMp*(conv_c_(ui)*derxyz_(0, vi)) ;
    estif_(vi*4 + 3, ui*4 + 1) += ttimetauMp*(conv_c_(ui)*derxyz_(1, vi)) ;
    estif_(vi*4 + 3, ui*4 + 2) += ttimetauMp*(conv_c_(ui)*derxyz_(2, vi)) ;
#endif
    
    /* Stabilisierung des Drucks (-L_visc_u) */
    /* 
                    /                              \
                   |               /  \             |
                   |  nabla o eps | Du | , nabla q  |
                   |               \  /             |
                    \                              /
    */
    estif_(vi*4 + 3, ui*4)     += 2.0*nu_*ttimetauMp*(derxyz_(0, vi)*viscs2_(0, 0, ui) + derxyz_(1, vi)*viscs2_(0, 1, ui) + derxyz_(2, vi)*viscs2_(0, 2, ui)) ;
    estif_(vi*4 + 3, ui*4 + 1) += 2.0*nu_*ttimetauMp*(derxyz_(0, vi)*viscs2_(0, 1, ui) + derxyz_(1, vi)*viscs2_(1, 1, ui) + derxyz_(2, vi)*viscs2_(1, 2, ui)) ;
    estif_(vi*4 + 3, ui*4 + 2) += 2.0*nu_*ttimetauMp*(derxyz_(0, vi)*viscs2_(0, 2, ui) + derxyz_(1, vi)*viscs2_(1, 2, ui) + derxyz_(2, vi)*viscs2_(2, 2, ui)) ;

    /* Stabilisierung des Drucks ( L_pres_p) */
    /*   
                          /                    \
                         |                      |
                         |  nabla Dp , nabla q  |
                         |                      |
                          \                    /
    */
    estif_(vi*4 + 3, ui*4 + 3) += ttimetauMp*(derxyz_(0, ui)*derxyz_(0, vi) + derxyz_(1, ui)*derxyz_(1, vi) + derxyz_(2, ui)*derxyz_(2, vi)) ;

    /* Divergenzfreiheit */
    /*
                            /                \
                           |                  |
                           | nabla o Du  , q  |
                           |                  |
                            \                /
    */
    estif_(vi*4 + 3, ui*4)     += timefacfac*funct_(vi)*derxyz_(0, ui) ;
    estif_(vi*4 + 3, ui*4 + 1) += timefacfac*funct_(vi)*derxyz_(1, ui) ;
    estif_(vi*4 + 3, ui*4 + 2) += timefacfac*funct_(vi)*derxyz_(2, ui) ;

    /* Kontinuitätsstabilisierung */
    /*   
                        /                        \
                       |                          |
                       | nabla o Du  , nabla o v  |
                       |                          |
                        \                        /
    */
    estif_(vi*4, ui*4)         += (thsl*thsl)*tau_C*derxyz_(0, ui)*derxyz_(0, vi) ;
    estif_(vi*4, ui*4 + 1)     += (thsl*thsl)*tau_C*derxyz_(0, vi)*derxyz_(1, ui) ;
    estif_(vi*4, ui*4 + 2)     += (thsl*thsl)*tau_C*derxyz_(0, vi)*derxyz_(2, ui) ;
    estif_(vi*4 + 1, ui*4)     += (thsl*thsl)*tau_C*derxyz_(0, ui)*derxyz_(1, vi) ;
    estif_(vi*4 + 1, ui*4 + 1) += (thsl*thsl)*tau_C*derxyz_(1, ui)*derxyz_(1, vi) ;
    estif_(vi*4 + 1, ui*4 + 2) += (thsl*thsl)*tau_C*derxyz_(1, vi)*derxyz_(2, ui) ;
    estif_(vi*4 + 2, ui*4)     += (thsl*thsl)*tau_C*derxyz_(0, ui)*derxyz_(2, vi) ;
    estif_(vi*4 + 2, ui*4 + 1) += (thsl*thsl)*tau_C*derxyz_(1, ui)*derxyz_(2, vi) ;
    estif_(vi*4 + 2, ui*4 + 2) += (thsl*thsl)*tau_C*derxyz_(2, ui)*derxyz_(2, vi) ;

    /* Massenterm */
    /*

                             /        \
                            |          |
                            |  Du , v  |
                            |          |
                             \        /
    */
    estif_(vi*4, ui*4)         += fac*funct_(ui)*funct_(vi) ;
    estif_(vi*4 + 1, ui*4 + 1) += fac*funct_(ui)*funct_(vi) ;
    estif_(vi*4 + 2, ui*4 + 2) += fac*funct_(ui)*funct_(vi) ;

    /* Konvektionsstabilisierung */
#if 0
    /*   
                    /                           \
                   |   n+1      /          \     |
                   |  u      , | Du o nabla | v  |
                   |   (i)      \          /     |
                    \                           /

   
                     /                        \
                    |        / n+1       \     |
                    |  Du , | u   o nabla | v  |
                    |        \ (i)       /     |
                     \                        /
    */
    estif_(vi*4, ui*4)         += timetauM*funct_(ui)*(2.0*velint_(0)*derxyz_(0, vi) + velint_(1)*derxyz_(1, vi) + velint_(2)*derxyz_(2, vi)) ;
    estif_(vi*4, ui*4 + 1)     += timetauM*funct_(ui)*velint_(0)*derxyz_(1, vi) ;
    estif_(vi*4, ui*4 + 2)     += timetauM*funct_(ui)*velint_(0)*derxyz_(2, vi) ;
    estif_(vi*4 + 1, ui*4)     += timetauM*funct_(ui)*velint_(1)*derxyz_(0, vi) ;
    estif_(vi*4 + 1, ui*4 + 1) += timetauM*funct_(ui)*(velint_(0)*derxyz_(0, vi) + 2.0*velint_(1)*derxyz_(1, vi) + velint_(2)*derxyz_(2, vi)) ;
    estif_(vi*4 + 1, ui*4 + 2) += timetauM*funct_(ui)*velint_(1)*derxyz_(2, vi) ;
    estif_(vi*4 + 2, ui*4)     += timetauM*funct_(ui)*velint_(2)*derxyz_(0, vi) ;
    estif_(vi*4 + 2, ui*4 + 1) += timetauM*funct_(ui)*velint_(2)*derxyz_(1, vi) ;
    estif_(vi*4 + 2, ui*4 + 2) += timetauM*funct_(ui)*(velint_(0)*derxyz_(0, vi) + velint_(1)*derxyz_(1, vi) + 2.0*velint_(2)*derxyz_(2, vi)) ;
#else
    /*   
                    /                        \
                   |        / n+1       \     |
                   |  Du , | u   o nabla | v  |
                   |        \ (i)       /     |
                    \                        /
    */
    estif_(vi*4, ui*4)         += timetauM*funct_(ui)*conv_c_(vi);
    estif_(vi*4 + 1, ui*4 + 1) += timetauM*funct_(ui)*conv_c_(vi);
    estif_(vi*4 + 2, ui*4 + 2) += timetauM*funct_(ui)*conv_c_(vi);
#endif

    
    /* Viskositätsstabilisierung */
    /*           
                         /                  \
                        |                    |
                        |  Du , div eps (v)  |
                        |                    |
                         \                  /
    */
    estif_(vi*4, ui*4)         += 2.0*nu_*timetauMp*funct_(ui)*viscs2_(0, 0, vi) ;
    estif_(vi*4, ui*4 + 1)     += 2.0*nu_*timetauMp*funct_(ui)*viscs2_(0, 1, vi) ;
    estif_(vi*4, ui*4 + 2)     += 2.0*nu_*timetauMp*funct_(ui)*viscs2_(0, 2, vi) ;
    estif_(vi*4 + 1, ui*4)     += 2.0*nu_*timetauMp*funct_(ui)*viscs2_(0, 1, vi) ;
    estif_(vi*4 + 1, ui*4 + 1) += 2.0*nu_*timetauMp*funct_(ui)*viscs2_(1, 1, vi) ;
    estif_(vi*4 + 1, ui*4 + 2) += 2.0*nu_*timetauMp*funct_(ui)*viscs2_(1, 2, vi) ;
    estif_(vi*4 + 2, ui*4)     += 2.0*nu_*timetauMp*funct_(ui)*viscs2_(0, 2, vi) ;
    estif_(vi*4 + 2, ui*4 + 1) += 2.0*nu_*timetauMp*funct_(ui)*viscs2_(1, 2, vi) ;
    estif_(vi*4 + 2, ui*4 + 2) += 2.0*nu_*timetauMp*funct_(ui)*viscs2_(2, 2, vi) ;

    /* Stabilisierung der Druckgleichung */
    /*                
                           /              \
                          |                |
                          |  Du , nabla q  |
                          |                |
                           \              /
    */

    estif_(vi*4 + 3, ui*4)     += timetauMp*funct_(ui)*derxyz_(0, vi) ;
    estif_(vi*4 + 3, ui*4 + 1) += timetauMp*funct_(ui)*derxyz_(1, vi) ;
    estif_(vi*4 + 3, ui*4 + 2) += timetauMp*funct_(ui)*derxyz_(2, vi) ;

    /* Konvektionsstabilisierung */
    /*   
                   /                             \
                  |              /          \     |
                  |  rhsint   , | Du o nabla | v  |
                  |              \          /     |
                   \                             /

    */
    estif_(vi*4    , ui*4)     += -(timetauM*funct_(ui)*derxyz_(0, vi)*rhsint_(0)) ;
    estif_(vi*4    , ui*4 + 1) += -(timetauM*funct_(ui)*derxyz_(1, vi)*rhsint_(0)) ;
    estif_(vi*4    , ui*4 + 2) += -(timetauM*funct_(ui)*derxyz_(2, vi)*rhsint_(0)) ;
    estif_(vi*4 + 1, ui*4)     += -(timetauM*funct_(ui)*derxyz_(0, vi)*rhsint_(1)) ;
    estif_(vi*4 + 1, ui*4 + 1) += -(timetauM*funct_(ui)*derxyz_(1, vi)*rhsint_(1)) ;
    estif_(vi*4 + 1, ui*4 + 2) += -(timetauM*funct_(ui)*derxyz_(2, vi)*rhsint_(1)) ;
    estif_(vi*4 + 2, ui*4)     += -(timetauM*funct_(ui)*derxyz_(0, vi)*rhsint_(2)) ;
    estif_(vi*4 + 2, ui*4 + 1) += -(timetauM*funct_(ui)*derxyz_(1, vi)*rhsint_(2)) ;
    estif_(vi*4 + 2, ui*4 + 2) += -(timetauM*funct_(ui)*derxyz_(2, vi)*rhsint_(2)) ;


  }
}
