/*-----------------------------------------------------------*/
/*!
\file str_nln_solver_ptc.H

\brief pseudo transient solution method

\maintainer Michael Hiermeier

\date Aug 13, 2015

\level 3

*/
/*-----------------------------------------------------------*/

#ifndef STR_NLN_SOLVER_PTC_H_
#define STR_NLN_SOLVER_PTC_H_

#include "str_nln_solver_nox.H"   // base class

namespace STR {
namespace NLN {
namespace SOLVER {
  class PseudoTransient : public Nox
  {
  public:
    //! constructor
    PseudoTransient();

    //! destructor
    virtual ~PseudoTransient() {};

    // derived from the base class
    virtual void Setup();

  protected:
    void SetPseudoTransientParams();

  };  // class PTC
} // namespace SOLVER
} // namespace NLN
} // namespace STR


#endif /* STR_NLN_SOLVER_PTC_H_ */
