for (int vi=0; vi<iel; ++vi)
{
  /* RHS source term */
  eforce_(vi) += fac*funct_(vi)*rhsint_ ;

  /* time right hand side */
  /* transient part from last time-step */ 
  eforce_(vi) += timefacl*funct_(vi)*phin_ ;

  /* convective and RHS part from last time-step */ 
  eforce_(vi) += timefacr*funct_(vi)*(convn_-rhsint_) ;

  /* diffusive part from last time-step */ 
  eforce_(vi) += timefacr*diffus*(derxy_(0, vi)*phidern_(0) + derxy_(1, vi)*phidern_(1)) ;

  /* Stabilization terms: */
  /* 1) transient stabilization */
  /* RHS source term */
  eforce_(vi) += -timetaufacl*funct_(vi)*rhsint_ ;

  /* transient part from last time-step */ 
  eforce_(vi) += -ttimetaufacl*funct_(vi)*phin_ ;

  /* convective and RHS part from last time-step */ 
  eforce_(vi) += -timetaufaclr*funct_(vi)*(convn_-rhsint_) ;

  /* diffusive part from last time-step */ 
  eforce_(vi) += timetaufaclr*funct_(vi)*diffn_ ;

  /* 2) convective stabilization */
  /* RHS source term */
  eforce_(vi) += taufac*conv_(vi)*rhsint_ ;

  /* transient part from last time-step */ 
  eforce_(vi) += timetaufacl*conv_(vi)*phin_ ;

  /* convective and RHS part from last time-step */ 
  eforce_(vi) += timetaufacr*conv_(vi)*(convn_-rhsint_) ;

  /* diffusive part from last time-step */ 
  eforce_(vi) += timetaufacr*conv_(vi)*diffn_ ;

  /* 3) diffusive stabilization (USFEM assumed here, sign change necessary for GLS) */
  /* RHS source term */
  eforce_(vi) += taufac*diff_(vi)*rhsint_ ;

  /* transient part from last time-step */ 
  eforce_(vi) += timetaufacl*diff_(vi)*phin_ ;

  /* convective and RHS part from last time-step */ 
  eforce_(vi) += timetaufacr*diff_(vi)*(convn_-rhsint_) ;

  /* diffusive part from last time-step */ 
  eforce_(vi) += timetaufacr*diff_(vi)*diffn_ ;

}
