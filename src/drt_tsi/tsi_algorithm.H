/*----------------------------------------------------------------------*/
/*!
\file tsi_algorithm.H

\brief Basis of all TSI algorithms that perform a coupling between the
       structural field equation and temperature field equations
\level 2
<pre>
  \maintainer  Alexander Seitz
               seitz@lnm.mw.tum.de
               http://www.lnm.mw.tum.de
               089 - 289-15271
</pre>
*/


/*----------------------------------------------------------------------*
 | definitions                                               dano 12/09 |
 *----------------------------------------------------------------------*/
#ifndef TSI_ALGORITHM_H
#define TSI_ALGORITHM_H


/*----------------------------------------------------------------------*
 | headers                                                   dano 12/09 |
 *----------------------------------------------------------------------*/
#include "../drt_adapter/adapter_algorithmbase.H"
#include <Epetra_Vector.h>


/*----------------------------------------------------------------------*
 | forward declarations                                      dano 02/12 |
 *----------------------------------------------------------------------*/
namespace DRT
{
  class Discretization;
}

namespace CONTACT
{
  class CoTSILagrangeStrategy;
  class CoNitscheStrategyTsi;
}  // namespace CONTACT

namespace ADAPTER
{
  class Thermo;
  class Structure;
  class MortarVolCoupl;
  class Coupling;
}  // namespace ADAPTER

namespace MORTAR
{
  class MultiFieldCoupling;
}


/*----------------------------------------------------------------------*
 |                                                           dano 12/09 |
 *----------------------------------------------------------------------*/
//! TSI: Thermo-Structure Interaction
namespace TSI
{
  //! TSI algorithm base
  //!
  //!
  //!  Base class of TSI algorithms. Derives from StructureBaseAlgorithm and
  //!  ThermoBaseAlgorithm with temperature field.
  //!  There can (and will) be different subclasses that implement different
  //!  coupling schemes.
  //!
  //!  \warning The order of calling the two BaseAlgorithm-constructors (that
  //!  is the order in which we list the base classes) is important here! In the
  //!  constructors control file entries are written. And these entries define
  //!  the order in which the filters handle the Discretizations, which in turn
  //!  defines the dof number ordering of the Discretizations... Don't get
  //!  confused. Just always list structure, thermo. In that order.
  //!
  //!  \author u.kue
  //!  \date 02/08
  class Algorithm : public ADAPTER::AlgorithmBase
  {
   public:
    //! create using a Epetra_Comm
    explicit Algorithm(const Epetra_Comm& comm);

    //! virtual destructor to support polymorph destruction
    virtual ~Algorithm();

    //! outer level time loop (to be implemented by deriving classes)
    virtual void TimeLoop() = 0;

    /// initialise TSI system
    virtual void SetupSystem() = 0;

    /// non-linear solve, i.e. (multiple) corrector
    virtual void Solve() = 0;

    //! read restart data
    virtual void ReadRestart(int step  //!< step number where the calculation is continued
        ) = 0;

    //! access to structural field
    const Teuchos::RCP<ADAPTER::Structure>& StructureField() { return structure_; }

    //! access to thermal field
    const Teuchos::RCP<::ADAPTER::Thermo>& ThermoField() { return thermo_; }

   protected:
    //! @name Time loop building blocks

    //! start a new time step
    virtual void PrepareTimeStep() = 0;

    //! calculate stresses, strains, energies
    virtual void PrepareOutput() = 0;

    //! take current results for converged and save for next time step
    virtual void Update();

    //! write output
    virtual void Output(bool forced_writerestart = false);

    //! communicate displacement vector to thermal field to enable their
    //! visualisation on the deformed body
    void OutputDeformationInThr(
        Teuchos::RCP<const Epetra_Vector> dispnp, Teuchos::RCP<DRT::Discretization> structdis);

    //@}

    //! @name Transfer methods

    //! apply temperature state on structure discretization
    virtual void ApplyThermoCouplingState(Teuchos::RCP<const Epetra_Vector> temp,
        Teuchos::RCP<const Epetra_Vector> temp_res = Teuchos::null);

    //! apply structural displacements and velocities on thermo discretization
    virtual void ApplyStructCouplingState(
        Teuchos::RCP<const Epetra_Vector> disp, Teuchos::RCP<const Epetra_Vector> vel);

    //! Get a ptr to the contact strategy from the structural field,
    //! store it in tsi and hand it to the thermal field
    void GetContactStrategy();

    //! Access to the dof coupling for matching grid TSI
    ADAPTER::Coupling& StructureThermoCoupling() { return *coupST_; }
    //@}

    //! @name Access methods

    //! velocity calculation given the displacements (like in FSI)
    Teuchos::RCP<const Epetra_Vector> CalcVelocity(Teuchos::RCP<const Epetra_Vector> dispnp);

    //! displacements at time n+1 for thermal output
    Teuchos::RCP<Epetra_MultiVector> dispnp_;

    //! temperatures at time n+1 for structure output
    //! introduced for non-matching discretizations
    Teuchos::RCP<Epetra_MultiVector> tempnp_;

    //@}


    //! @name Underlying fields

    //! underlying structure of the FSI problem
    Teuchos::RCP<::ADAPTER::Structure> structure_;

    //! underlying fluid of the FSI problem
    Teuchos::RCP<::ADAPTER::Thermo> thermo_;

    //! contact strategies
    Teuchos::RCP<CONTACT::CoTSILagrangeStrategy> contact_strategy_lagrange_;
    Teuchos::RCP<CONTACT::CoNitscheStrategyTsi> contact_strategy_nitsche_;

    //@}

    //! @name Volume Mortar stuff

    //! flag for matchinggrid
    const bool matchinggrid_;
    //! volume coupling (using mortar) adapter
    Teuchos::RCP<::ADAPTER::MortarVolCoupl> volcoupl_;

    Teuchos::RCP<ADAPTER::Coupling> coupST_;  // S: master, T: slave
    //@}


    //! @name Surface Mortar stuff
    Teuchos::RCP<MORTAR::MultiFieldCoupling> mortar_coupling_;

    //@}
  };  // Algorithm
}  // namespace TSI


/*----------------------------------------------------------------------*/
#endif  // TSI_ALGORITHM_H
