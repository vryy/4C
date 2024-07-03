/*-----------------------------------------------------------*/
/*! \file

\brief Factory to build the desired linear solver std::map corresponding to the active model types


\level 3

*/
/*-----------------------------------------------------------*/


#ifndef FOUR_C_STRUCTURE_NEW_SOLVER_FACTORY_HPP
#define FOUR_C_STRUCTURE_NEW_SOLVER_FACTORY_HPP

#include "4C_config.hpp"

#include <Teuchos_RCP.hpp>

#include <set>

namespace Teuchos
{
  class ParameterList;
}

FOUR_C_NAMESPACE_OPEN

// forward declarations...
namespace Inpar
{
  namespace Solid
  {
    enum ModelType : int;
  }
  namespace CONTACT
  {
    enum SolvingStrategy : int;
    enum SystemType : int;
  }  // namespace CONTACT
}  // namespace Inpar
namespace Core::LinAlg
{
  class Solver;
}
namespace Core::FE
{
  class Discretization;
}  // namespace Core::FE
namespace Solid
{
  namespace SOLVER
  {
    /*! \brief Factory to build the desired linear solver std::map corresponding
     *  to the active model types.
     *
     *  \author Michael Hiermeier */
    class Factory
    {
     private:
      typedef std::map<enum Inpar::Solid::ModelType, Teuchos::RCP<Core::LinAlg::Solver>> LinSolMap;

     public:
      //! constructor
      Factory();

      //! destructor
      virtual ~Factory() = default;

      //! build the desired linear solvers
      Teuchos::RCP<LinSolMap> build_lin_solvers(
          const std::set<enum Inpar::Solid::ModelType>& modeltypes,
          const Teuchos::ParameterList& sdyn, Core::FE::Discretization& actdis) const;

      //! create the meshtying/contact linear solver
      static Teuchos::RCP<Core::LinAlg::Solver> build_meshtying_contact_lin_solver(
          Core::FE::Discretization& actdis, enum Inpar::CONTACT::SolvingStrategy sol_type,
          enum Inpar::CONTACT::SystemType sys_type, const int lin_solver_id);

     private:
      //! create the structural linear solver (should be called by default)
      Teuchos::RCP<Core::LinAlg::Solver> build_structure_lin_solver(
          const Teuchos::ParameterList& sdyn, Core::FE::Discretization& actdis) const;

      //! create the meshtying/contact linear solver
      Teuchos::RCP<Core::LinAlg::Solver> build_meshtying_contact_lin_solver(
          Core::FE::Discretization& actdis) const;

      //! create the Lagrange/penalty enforced constraint linear solver
      Teuchos::RCP<Core::LinAlg::Solver> build_lag_pen_constraint_lin_solver(
          const Teuchos::ParameterList& sdyn, Core::FE::Discretization& actdis) const;

      //! create the Windkessel linear solver
      Teuchos::RCP<Core::LinAlg::Solver> build_cardiovascular0_d_lin_solver(
          const Teuchos::ParameterList& sdyn, Core::FE::Discretization& actdis) const;

    };  // class Factory

    /*! Non-member function, which relates to the Solid::SOLVER::Factory class
     *  Please call this method from outside! */
    Teuchos::RCP<std::map<enum Inpar::Solid::ModelType, Teuchos::RCP<Core::LinAlg::Solver>>>
    build_lin_solvers(const std::set<enum Inpar::Solid::ModelType>& modeltypes,
        const Teuchos::ParameterList& sdyn, Core::FE::Discretization& actdis);
  }  // namespace SOLVER
}  // namespace Solid

FOUR_C_NAMESPACE_CLOSE

#endif
