/*---------------------------------------------------------------------*/
/*! \file

\brief for an exact gradient, we need fully linearized adjoint equations with
 respect to the primal variables (velocity, pressure). Therefore, the stabilization
 parameter has to be derived with respect to velocity, which is not done
 currently. In order to test the adjoints and the gradient computation, velocity
 independent stabilization parameters are used by these defines. You have to
 recompute the stabilization parameters in the fluid field, in the adjoint field
 and in the optimizers computation of gradient and target functional. Moreover,
 the new, consistent One-Step-Theta scheme has to be used as well as a gauss-point
 based computation of the bodyforce.
 I am currently setting the parameter for stationary testing constant and for
 instationary testing spatially constant and temporally instationary so that
 the code at the three files looks similar like
 #ifdef TAU_CONST_STAT
      tau_(0) = 1.0e-04;
      tau_(1) = 1.0e-04;
      tau_(2) = 1.0;
      tau_old_(0) = 1.0e-04;
      tau_old_(1) = 1.0e-04;
      tau_old_(2) = 1.0;
 #endif
 #ifdef TAU_CONST_INSTAT
      tau_(0) = 7.0e-05+1.0e-03*fldAdPara_->Time();
      tau_(1) = 7.0e-05+1.0e-03*fldAdPara_->Time();
      tau_(2) = 1.0+10*fldAdPara_->Time();
      tau_old_(0) = 7.0e-05+1.0e-03*(fldAdPara_->Time()+fldAdPara_->Dt());
      tau_old_(1) = 7.0e-05+1.0e-03*(fldAdPara_->Time()+fldAdPara_->Dt());
      tau_old_(2) = 1.0+10*(fldAdPara_->Time()+fldAdPara_->Dt());
 #endif

\maintainer Martin Kronbichler

\level 3

*/
/*---------------------------------------------------------------------*/


#ifndef OPTI_DEFINES_H_
#define OPTI_DEFINES_H_

//#define TAU_CONST_STAT
//#define TAU_CONST_INSTAT

#endif /* OPTI_DEFINES_H_ */
