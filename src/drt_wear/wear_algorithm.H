/*----------------------------------------------------------------------*/
/*!
\file wear_algorithm.H

\brief Basis of all WEAR algorithms that perform a coupling between the
       structural field equation and ALE field equations

\level 2

\maintainer Alexander Popp
*/


/*----------------------------------------------------------------------*
 | definitions                                              farah 11/13 |
 *----------------------------------------------------------------------*/
#ifndef WEAR_ALGORITHM_H_
#define WEAR_ALGORITHM_H_


/*----------------------------------------------------------------------*
 | headers                                                  farah 11/13 |
 *----------------------------------------------------------------------*/
#include "../drt_adapter/adapter_algorithmbase.H"

#include <Epetra_Vector.h>


/*----------------------------------------------------------------------*
 | forward declarations                                     farah 11/13 |
 *----------------------------------------------------------------------*/
namespace DRT
{
  class Discretization;
}

namespace ADAPTER
{
  class Structure;
  class FSIStructureWrapper;
  class AleWearWrapper;
}

namespace ALE
{
  class Ale;
}

namespace MORTAR
{
  class ManagerBase;
}

namespace CONTACT
{
  class CoInterface;
}

/*----------------------------------------------------------------------*
 |                                                          farah 11/13 |
 *----------------------------------------------------------------------*/
//! WEAR: Structure-ALE Interaction for wear
namespace WEAR
{
  class Algorithm : public ADAPTER::AlgorithmBase
  {

  public:

    //! create using a Epetra_Comm
    explicit Algorithm(const Epetra_Comm& comm);

    //! virtual destructor to support polymorph destruction
    virtual ~Algorithm();

    //! outer level time loop (to be implemented by deriving classes)
    virtual void TimeLoop() = 0;

    //! read restart data
    virtual void ReadRestart(
      int step  //!< step number where the calculation is continued
      ) = 0;

    //! access to structural field
    Teuchos::RCP<ADAPTER::FSIStructureWrapper> StructureField(){return structure_;}

    //! access to ALE field
    ADAPTER::AleWearWrapper& AleField() { return *ale_; }

  private:

    //! check compatibility if input parameters
    void CheckInput();

    //! create mortar interfaces for material conf.
    void CreateMaterialInterface();

    //! @name Underlying fields
    Teuchos::RCP< ::ADAPTER::FSIStructureWrapper>      structure_;  //! underlying structure
    Teuchos::RCP<ADAPTER::AleWearWrapper>              ale_;        //! underlying ALE

  protected:
    int                                                dim_;        //! problem dimension
    //@}

    Teuchos::RCP<MORTAR::ManagerBase>                  cmtman_;        // contact manager
    std::vector<Teuchos::RCP<CONTACT::CoInterface> >   interfaces_;    // contact/wear interfaces
    std::vector<Teuchos::RCP<CONTACT::CoInterface> >   interfacesMat_; // contact interfaces in Mat.

  };  // Algorithm
}  // namespace WEAR


/*----------------------------------------------------------------------*/
#endif /* WEAR_ALGORITHM_H_ */
