/*
 * str_solver_factory.H
 *
 *  Created on: Sep 9, 2015
 *      Author: hiermeier
 */

#ifndef STR_SOLVER_FACTORY_H_
#define STR_SOLVER_FACTORY_H_

#include <Teuchos_RCP.hpp>
#include "../drt_inpar/inpar_structure.H"

namespace LINALG{
  class Solver;
}
namespace DRT{
class Discretization;
}
namespace Teuchos{
class ParameterList;
}
namespace STR {
namespace SOLVER {
  class Factory
  {
  public:
    //! constructor
    Factory();

    //! destructor
    virtual ~Factory(){};

    //! build the desired linear solvers
    Teuchos::RCP<std::map<enum INPAR::STR::ModelType, Teuchos::RCP<LINALG::Solver> > >
        BuildLinSolvers(const std::set<enum INPAR::STR::ModelType>& modeltypes,
                        const Teuchos::ParameterList& sdyn,
                        DRT::Discretization& actdis) const;

  private:
    //! create the structural linear solver (should be called by default)
    Teuchos::RCP<LINALG::Solver> BuildStructureLinSolver(
        const Teuchos::ParameterList& sdyn,
        DRT::Discretization& actdis) const;

    //! create the meshtying/contact linear solver
    Teuchos::RCP<LINALG::Solver> BuildMeshtyingContactLinSolver(
        const Teuchos::ParameterList& sdyn,
        DRT::Discretization& actdis) const;

    //! create the Windkessel linear solver
    Teuchos::RCP<LINALG::Solver> BuildLagPenConstraintLinSolver(
        const Teuchos::ParameterList& sdyn,
        DRT::Discretization& actdis) const;

    //! create the Lagrange/penalty enforced constraint linear solver
    Teuchos::RCP<LINALG::Solver> BuildWindkesselLinSolver(
        const Teuchos::ParameterList& sdyn,
        DRT::Discretization& actdis) const;

    //! create the spring dashpot linear solver
    Teuchos::RCP<LINALG::Solver> BuildSpringDashpotLinSolver(
        const Teuchos::ParameterList& sdyn,
        DRT::Discretization& actdis) const;
  }; // class Factory

  //! non-member function
  Teuchos::RCP<std::map<enum INPAR::STR::ModelType, Teuchos::RCP<LINALG::Solver> > >
      BuildLinSolvers(const std::set<enum INPAR::STR::ModelType>& modeltypes,
                      const Teuchos::ParameterList& sdyn,
                      DRT::Discretization& actdis);

} // namespace SOLVER
} // namespace STR

#endif /* STR_SOLVER_FACTORY_H_ */
