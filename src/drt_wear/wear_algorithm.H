/*----------------------------------------------------------------------*/
/*!
\file wear_algorithm.H

\brief Basis of all WEAR algorithms that perform a coupling between the
       structural field equation and ALE field equations

<pre>
Maintainer: Philipp Farah
            farah@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15257
</pre>
*/


/*----------------------------------------------------------------------*
 | definitions                                              farah 11/13 |
 *----------------------------------------------------------------------*/
#ifndef WEAR_ALGORITHM_H_
#define WEAR_ALGORITHM_H_


/*----------------------------------------------------------------------*
 | headers                                                  farah 11/13 |
 *----------------------------------------------------------------------*/
#include "../drt_adapter/ad_str_fsiwrapper.H"
#include "../drt_adapter/adapter_algorithmbase.H"
#include "../drt_ale/ale.H"

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
  //class Structure;
  class coupling;
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

    //! create mortar interfaces for material conf.
    void CreateMaterialInterface();

    //! outer level time loop (to be implemented by deriving classes)
    virtual void TimeLoop() = 0;

    //! check compatibility if input parameters
    void CheckInput();

    //! read restart data
    virtual void ReadRestart(
      int step  //!< step number where the calculation is continued
      ) = 0;

    //! access to structural field
    Teuchos::RCP<ADAPTER::FSIStructureWrapper> StructureField(){return structure_;}

    //! access to ALE field
    ALE::Ale& AleField() { return *ale_; }

  protected:

    int                                                dim_;        //! problem dimension

    //! @name Underlying fields
    Teuchos::RCP< ::ADAPTER::FSIStructureWrapper>      structure_;  //! underlying structure
    Teuchos::RCP<ALE::Ale>                             ale_;        //! underlying ALE
    //@}

    Teuchos::RCP<MORTAR::ManagerBase>                  cmtman_;        // contact manager
    std::vector<Teuchos::RCP<CONTACT::CoInterface> >   interfaces_;    // contact/wear interfaces
    std::vector<Teuchos::RCP<CONTACT::CoInterface> >   interfacesMat_; // contact interfaces in Mat.

  };  // Algorithm
}  // namespace WEAR


/*----------------------------------------------------------------------*/
#endif /* WEAR_ALGORITHM_H_ */
