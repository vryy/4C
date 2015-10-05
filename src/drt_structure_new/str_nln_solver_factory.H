/*
 * str_nln_solver_factory.H
 *
 *  Created on: Aug 13, 2015
 *      Author: farah
 */

#ifndef STR_NLN_SOLVER_FACTORY_H_
#define STR_NLN_SOLVER_FACTORY_H_

#include "../drt_inpar/inpar_structure.H"

// forward declaration
namespace STR {
namespace NLN {
namespace SOLVER {
  class Generic;

  class Factory
  {
  public:
    //! constructor
    Factory();

    //! destructor
    virtual ~Factory() {};

    //! build the specific nonlinear solver
    Teuchos::RCP<STR::NLN::SOLVER::Generic>  BuildNlnSolver(
        const enum INPAR::STR::NonlinSolTech& nlnSolType);
  };

  //! non-member function (please call this one!)
  Teuchos::RCP<STR::NLN::SOLVER::Generic> BuildNlnSolver(
      const enum INPAR::STR::NonlinSolTech& nlnSolType);
} // namespace SOLVER
} // namespace NLN
} // namespace STR


#endif /* STR_NLN_SOLVER_FACTORY_H_ */
