/*----------------------------------------------------------------------*/
/*! \file
 \brief utils methods for for porous multiphase flow through elastic medium problems

   \level 3

 *----------------------------------------------------------------------*/

#ifndef FOUR_C_POROMULTIPHASE_UTILS_HPP
#define FOUR_C_POROMULTIPHASE_UTILS_HPP

#include "baci_config.hpp"

#include "baci_inpar_poromultiphase.hpp"

#include <Epetra_Vector.h>
#include <Teuchos_RCP.hpp>

#include <set>

BACI_NAMESPACE_OPEN

namespace DRT
{
  class Discretization;
}

namespace ADAPTER
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
    void AssignMaterialPointers(
        const std::string& struct_disname, const std::string& fluid_disname);

    /// create solution algorithm depending on input file
    Teuchos::RCP<ADAPTER::PoroMultiPhase> CreatePoroMultiPhaseAlgorithm(
        INPAR::POROMULTIPHASE::SolutionSchemeOverFields
            solscheme,                             //!< solution scheme to build (i)
        const Teuchos::ParameterList& timeparams,  //!< problem parameters (i)
        const Epetra_Comm& comm                    //!< communicator(i)
    );

    //! Determine norm of vector
    double CalculateVectorNorm(const enum INPAR::POROMULTIPHASE::VectorNorm norm,  //!< norm to use
        const Teuchos::RCP<const Epetra_Vector> vect  //!< the vector of interest
    );

  }  // namespace UTILS
  // Print the logo
  void PrintLogo();
}  // namespace POROMULTIPHASE



BACI_NAMESPACE_CLOSE

#endif  // POROMULTIPHASE_UTILS_H
