/*----------------------------------------------------------------------*/
/*! \file
 \brief utils methods for for porous multiphase flow through elastic medium problems

   \level 3

 *----------------------------------------------------------------------*/

#ifndef FOUR_C_POROMULTIPHASE_UTILS_HPP
#define FOUR_C_POROMULTIPHASE_UTILS_HPP

#include "4C_config.hpp"

#include "4C_inpar_poromultiphase.hpp"

#include <Epetra_Vector.h>
#include <Teuchos_RCP.hpp>

#include <set>

FOUR_C_NAMESPACE_OPEN

namespace Core::FE
{
  class Discretization;
}  // namespace Core::FE

namespace Adapter
{
  class PoroMultiPhase;
}

namespace POROMULTIPHASE
{
  namespace UTILS
  {
    /// setup discretizations and dofsets
    std::map<int, std::set<int>> SetupDiscretizationsAndFieldCoupling(const Epetra_Comm& comm,
        const std::string& struct_disname, const std::string& fluid_disname, int& nds_disp,
        int& nds_vel, int& nds_solidpressure);

    //! exchange material pointers of both discretizations
    void assign_material_pointers(
        const std::string& struct_disname, const std::string& fluid_disname);

    /// create solution algorithm depending on input file
    Teuchos::RCP<Adapter::PoroMultiPhase> CreatePoroMultiPhaseAlgorithm(
        Inpar::POROMULTIPHASE::SolutionSchemeOverFields
            solscheme,                             //!< solution scheme to build (i)
        const Teuchos::ParameterList& timeparams,  //!< problem parameters (i)
        const Epetra_Comm& comm                    //!< communicator(i)
    );

    //! Determine norm of vector
    double calculate_vector_norm(
        const enum Inpar::POROMULTIPHASE::VectorNorm norm,  //!< norm to use
        const Teuchos::RCP<const Epetra_Vector> vect        //!< the vector of interest
    );

  }  // namespace UTILS
  // Print the logo
  void PrintLogo();
}  // namespace POROMULTIPHASE



FOUR_C_NAMESPACE_CLOSE

#endif
