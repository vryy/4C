/*!-----------------------------------------------------------------------------------------------*
\file levelset_algorithm.H

  \brief level-set algorithm

<pre>
Maintainer: Ursula Rasthofer
            rasthofer@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15236
</pre>
 *------------------------------------------------------------------------------------------------*/

#ifndef LEVELSET_ALGORITHM_H
#define LEVELSET_ALGORITHM_H

#include <Epetra_MpiComm.h>
#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>


/*----------------------------------------------------------------------*
 | forward declarations                                     ghamm 11/12 |
 *----------------------------------------------------------------------*/
namespace SCATRA
{
  class ScaTraTimIntImpl;
}
namespace PARTICLE
{
  class ScatraParticleCoupling;
}


/*----------------------------------------------------------------------*
 | cavitation algorithm                                     ghamm 11/12 |
 *----------------------------------------------------------------------*/
/// SCATRA
namespace SCATRA
{

  class LevelSetAlgorithm
  {
  public:
    LevelSetAlgorithm(
      const Epetra_Comm& comm,              ///< communicator
      const Teuchos::ParameterList& prbdyn, ///< problem-specific parameters
      bool isale,  ///< do we need an ALE formulation of the fields?
      const std::string disname, ///< scatra discretization name
      const Teuchos::ParameterList& solverparams
      );

    /// virtual destructor to support polymorph destruction
    virtual ~LevelSetAlgorithm();

  protected:


  private:
    Teuchos::RCP<SCATRA::ScaTraTimIntImpl> scatra_;

    Teuchos::RCP<PARTICLE::ScatraParticleCoupling> particle_;
  };

}

#endif // #ifndef LEVELSET_ALGORITHM_H
