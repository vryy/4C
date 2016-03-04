/*-----------------------------------------------------------*/
/*!
\file str_solver_factory.H

\maintainer Michael Hiermeier

\date Sep 9, 2015

\level 3

*/
/*-----------------------------------------------------------*/


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
  /*! \brief Factory to build the desired linear solver std::map corresponding
   *  to the active model types.
   *
   *  \author Michael Hiermeier */
  class Factory
  {
  private:
    typedef std::map<enum INPAR::STR::ModelType, Teuchos::RCP<LINALG::Solver> > LinSolMap;

  public:
    //! constructor
    Factory();

    //! destructor
    virtual ~Factory(){};

    //! build the desired linear solvers
    Teuchos::RCP<LinSolMap> BuildLinSolvers(
        const std::set<enum INPAR::STR::ModelType>& modeltypes,
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

  }; // class Factory

  /*! Non-member function, which relates to the STR::SOLVER::Factory class
   *  Please call this method from outside! */
  Teuchos::RCP<std::map<enum INPAR::STR::ModelType, Teuchos::RCP<LINALG::Solver> > >
      BuildLinSolvers(const std::set<enum INPAR::STR::ModelType>& modeltypes,
                      const Teuchos::ParameterList& sdyn,
                      DRT::Discretization& actdis);
} // namespace SOLVER
} // namespace STR

#endif /* STR_SOLVER_FACTORY_H_ */
