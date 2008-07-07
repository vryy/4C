for (int vi=0; vi<iel; ++vi)
{
  /* RHS source term */
  eforce_(vi) += fac*funct_(vi)*rhsint_ ;

  /* transient stabilization of RHS source term */
  /*eforce_(vi) += -taufac*funct_(vi)*rhsint_ ;*/

  /* convective stabilization of RHS source term */
  eforce_(vi) += taufac*conv_(vi)*rhsint_ ;

  /* diffusive stabilization of RHS source term */
  eforce_(vi) += taufac*diff_(vi)*rhsint_ ;

}
