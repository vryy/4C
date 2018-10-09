/*-----------------------------------------------------------*/
/*!
\file str_nln_solver_ptc.H

\brief pseudo transient solution method

\maintainer Michael Hiermeier

\level 3

*/
/*-----------------------------------------------------------*/

#ifndef STR_NLN_SOLVER_PTC_H_
#define STR_NLN_SOLVER_PTC_H_

#include "str_nln_solver_nox.H"  // base class

namespace STR
{
  namespace NLN
  {
    namespace SOLVER
    {
      /*! \brief Interface to pseudo transient continuation (PTC) in structural dynamics
       *
       *  The nonlinear solve is done using NOX. PTC is realized by setting specific NOX parameters.
       *
       *  For updating the pseudo-time-step-size, only the "SER" mechanism as available for now.
       *
       *  <h3> References </h3>
       *
       *  - Kelley CT, Keyes DE:
       *    Convergence analysis of pseudo-transient continuation,
       *    SIAM Journal on Numerical Analysis, 35:508--523 (1998)
       *  - Coffey TS, Kelley CT, Keyes DE:
       *    Pseudo-transient continuation and differential-algebraic equations,
       *    SIAM Journal on Scientific Computing, 25:553--569 (2003)
       *  - Kelley CT, Liao L-Z, Qi L, Chu MT, Reese JP, Winton C:
       *    Projected pseudo-transient continuation,
       *    SIAM Journal on Numerical Analysis, 46(6):3071--3083 (2007)
       *  - Gee MW, Kelley CT, Lehoucq RB:
       *    Pseudo-transient continuation for nonlinear transient elasticity,
       *    International Journal for Numerical Methods in Engineering, 78(10):1209--1219 (2009)
       */
      class PseudoTransient : public Nox
      {
       public:
        //! constructor
        PseudoTransient();

        //! destructor
        virtual ~PseudoTransient(){};

        //! derived from the base class
        virtual void Setup();

       protected:
        //! Insert PTC configuration into NOX parameter list
        void SetPseudoTransientParams();

      };  // class PTC
    }     // namespace SOLVER
  }       // namespace NLN
}  // namespace STR


#endif /* STR_NLN_SOLVER_PTC_H_ */
