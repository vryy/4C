/*----------------------------------------------------------------------*/
/*! \file
 \brief utils methods for for porous multiphase flow through elastic medium problems

   \level 3

 *----------------------------------------------------------------------*/

#ifndef BACI_POROMULTIPHASE_UTILS_H
#define BACI_POROMULTIPHASE_UTILS_H

#include "baci_inpar_poromultiphase.H"

#include <Epetra_Vector.h>
#include <Teuchos_RCP.hpp>

#include <set>

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



#endif  // POROMULTIPHASE_UTILS_H
