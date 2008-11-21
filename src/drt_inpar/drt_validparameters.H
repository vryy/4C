/*----------------------------------------------------------------------*/
/*!
\file drt_validparameters.H

\brief Setup of the list of valid input parameters

<pre>
Maintainer: Ulrich Kuettler
            kuettler@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15238
</pre>
*/
/*----------------------------------------------------------------------*/

#ifdef CCADISCRET

#ifndef DRT_VALIDPARAMETERS_H
#define DRT_VALIDPARAMETERS_H

#ifdef __cplusplus /* stuff which is c++ specific */

#include <string>
#include <iostream>

#include <Teuchos_ParameterList.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_StandardParameterEntryValidators.hpp>


/// Namespace for Baci enums of input parameters
///
/// This namespace should one day contain all the enums now still
/// defined in the various old C-headers (enums.h, dynamic.h,...)
/// The enums are used to represent the different values of string
/// input variables within Baci. The idea behind this new namespace
/// is twofold: first, of course, the old C-headers should go away
/// one day and we have to establish an alternative for them.
/// Second, it makes sense to define the valid parameters on a
/// global level and not within the different types of problems.
/// Of course, then drt_validparameters.H has to be included
/// wherever we want to work with the enums!
/// (a.popp 04/08)
/// \date 11/08
namespace INPAR {}

/*!
\brief INPUTPARAMS: Namespace for baci enums of input parameters

// This namespace should one day contain all the enums now still
// defined in the various old c-headers (enums.h, dynamic.h,...)
// The enums are used to represent the different values of string
// input variables within baci. The idea behind this new namespace
// is twofold: first, of course, the old c-headers should go away
// one day and we have to establish an alternative for them.
// Second, it makes sense to define the valid parameters on a
// global level and not within the different types of problems.
// Of course, then drt_validparameters.H has to be included
// wherever we want to work with the enums!
// (a.popp 04/08)

*/
namespace INPUTPARAMS
{
//-------------------------------------------STATISTICAL MECHANICS

  /// Type of thermal bath affecting elements
  /// (this enum represents the input file parameter THERMALBATH)
  enum ThermalBathType
  {
    thermalbath_none,
    thermalbath_uniform,
    thermalbath_shearflow
  };
  /// (this enum represents the input file parameter STAT_OUTPUT)
  enum StatOutput
  {
    statout_none,
    statout_endtoendlog,
    statout_endtoendergodicity,
    statout_viscoelasticity
  };
  //-------------------------------------------STRUCTURAL CONTACT

  /// Type of structural contact problem
  /// (this enum represents the input file parameter CONTACT)
  enum ContactType
  {
    contact_none,
    contact_normal,
    contact_frictional,
    contact_meshtying
  };

  /// Type of contact friction law
  /// (this enum represents the input file parameter FRICTION)
  enum ContactFrictionType
  {
    friction_none,
    friction_stick,
    friction_tresca,
    friction_coulomb
  };

  //----------------------------------------------------------FSI

  /// Type of partitioned coupling for FSI problems
  enum FSIPartitionedCouplingMethod
  {
    fsi_DirichletNeumann,
    fsi_RobinNeumann,
    fsi_DirichletRobin,
    fsi_RobinRobin
  };

  /// Linear solver algorithm for monolithic block system in monolithic FSI
  enum FSILinearBlockSolver
  {
    fsi_PreconditionedKrylov,
    fsi_PartitionedAitken,
    fsi_PartitionedVectorExtrapolation,
    fsi_PartitionedJacobianFreeNewtonKrylov,
    fsi_BGSAitken,
    fsi_BGSVectorExtrapolation,
    fsi_BGSJacobianFreeNewtonKrylov
  };

  //--------------------------------------SCALAR TRANSPORT DYNAMIC

  /// time integration schemes for scalar transport problems
  enum ScaTraTimeIntegrationScheme
  {
    timeint_stationary,
    timeint_one_step_theta,
    timeint_bdf2,
    timeint_gen_alpha
  };
  
  //--------------------------------------STRUCTURAL DYNAMIC
  
  /// possible constraint solvers
  enum ConSolveAlgo
  {
    consolve_direct,
    consolve_iterative
  };
}

namespace DRT {
namespace INPUT {

  /// local wrapper for Teuchos::setIntParameter() that allows only integers
  void IntParameter(std::string const& paramName,
                    int const value,
                    std::string const& docString,
                    Teuchos::ParameterList* paramList);

  /// local wrapper for Teuchos::setDoubleParameter() that allows only doubles
  void DoubleParameter(std::string const& paramName,
                       double const& value,
                       std::string const& docString,
                       Teuchos::ParameterList* paramList);

  /// construct list with all parameters and documentation
  Teuchos::RCP<const Teuchos::ParameterList> ValidParameters();

  /// set parameters read by our solvers
  void SetValidSolverParameters(Teuchos::ParameterList& list);

  /// set parameters for time adaptivity
  void SetValidTimeAdaptivityParameters(Teuchos::ParameterList& list);

  /// print all parameters that have a default value
  void PrintDefaultParameters(std::ostream& stream, const Teuchos::ParameterList& list);

  /// print flag sections of dat file with given list
  void PrintDatHeader(std::ostream& stream, const Teuchos::ParameterList& list, std::string parentname="", bool color=false);
}
}

#endif

/* headers and prototypes of functions also acessible from c */
#ifdef __cplusplus
extern "C"
{
#endif

  /*! print list of valid parameters with documentation */
  void PrintValidParameters();

  /*! print flag sections of dat file with default flags */
  void PrintDefaultDatHeader();

#ifdef __cplusplus
}
#endif

#endif
#endif
