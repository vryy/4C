/*----------------------------------------------------------------------*/
/*!
\file wear_partitioned.H

\brief  Basis of partitioned wear algorithm
        (Lagrangian step followed by Eulerian step )
<pre>
Maintainer: Philipp Farah
            farah@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15257
</pre>
*/

/*----------------------------------------------------------------------*
 | definitions                                               farah 11/13 |
 *----------------------------------------------------------------------*/

#ifndef WEAR_PARTITIONED_H
#define WEAR_PARTITIONED_H

/*----------------------------------------------------------------------*
 | headers                                                   farah 11/13 |
 *----------------------------------------------------------------------*/
//#include <Epetra_Comm.h>
//
#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_StandardParameterEntryValidators.hpp>

#include "wear_algorithm.H"

#include "../drt_adapter/ad_str_structure.H"
#include "../drt_adapter/ad_str_fsiwrapper.H"
#include "../drt_lib/drt_element.H"
#include "../drt_ale/ale.H"

/*----------------------------------------------------------------------*
 | forward declarations                                     farah 11/13 |
 *----------------------------------------------------------------------*/
namespace MORTAR
{
  class ManagerBase;
}
namespace ADAPTER
{
  class Coupling;
}
namespace DRT
{
  class LocationArray;
}
/*----------------------------------------------------------------------*
 |                                                          farah 11/13 |
 *----------------------------------------------------------------------*/
/// STRUE_ALE: Structure with ale
namespace WEAR
{
  /// WEAR stru_ale partitioned algorithm

  class Partitioned : public Algorithm
  {
  public:

    /// create using a Epetra_Comm
    explicit Partitioned(const Epetra_Comm& comm);

    /// virtual destructor to support polymorph destruction
    virtual ~Partitioned();

    // do ale step
    void AleStep(Teuchos::RCP<Epetra_Vector> idisale_global);

    // transform from ale to structure dofs
    virtual Teuchos::RCP<Epetra_Vector> AleToStructure(Teuchos::RCP<Epetra_Vector> vec) const;

    // transform from ale to structure dofs
    virtual Teuchos::RCP<Epetra_Vector> AleToStructure(Teuchos::RCP<const Epetra_Vector> vec) const;

    /// Application of mesh displacements (frictional contact)
    void ApplyMeshDisplacement(bool iterated=false);

    /// material coordinates evaluated from spatial ones
    void AdvectionMap(double* XMat,
                      double* XMesh,
                      DRT::Element** ElementPtr,
                      int numelements);

    // check convergence
    bool ConvergenceCheck(int iter);

    // Dof Coupling
    void DispCoupling (Teuchos::RCP<Epetra_Vector>& disinterface);

    /// Interface displacements (frictional contact)
    void InterfaceDisp (Teuchos::RCP<Epetra_Vector>& disinterface_s,
                        Teuchos::RCP<Epetra_Vector>& disinterface_m);

    // Merge wear from slave and master surface to one wear vector
    void MergeWear (Teuchos::RCP<Epetra_Vector>& disinterface_s,
                    Teuchos::RCP<Epetra_Vector>& disinterface_m,
                    Teuchos::RCP<Epetra_Vector>& disinterface_g);

    // write output for ale and structure
    void Output();

    // ale parameter list
    virtual Teuchos::ParameterList& ParamsAle() {return alepara_;};

    // prepare time step for ale and structure
    void PrepareTimeStep();

    /// read restart data
    virtual void ReadRestart(
      int step  ///< step number where the calculation is continued
      );

    // redistribute material interfaces according to current interfaces
    void RedistributeMatInterfaces();

    // transform from ale to structure dofs
    virtual Teuchos::RCP<Epetra_Vector> StructureToAle(Teuchos::RCP<Epetra_Vector> vec) const;

    // transform from ale to structure dofs
    virtual Teuchos::RCP<Epetra_Vector> StructureToAle(Teuchos::RCP<const Epetra_Vector> vec) const;

    /// general time loop
    virtual void TimeLoop();

    // time loop for staggered coupling
    void TimeLoopStagg(bool alestep);

    //time loop for iterative stagered coupling
    void TimeLoopIterStagg();

    // update ale and structure
    void Update();

    // update spatial displacements due to mat. displ
    void UpdateDispnp();

    /// pull-back operation for wear from current to material conf.
    void WearPullBackSlave(Teuchos::RCP<Epetra_Vector>& disinterface_s);

    /// pull-back operation for wear from current to material conf.
    void WearPullBackMaster(Teuchos::RCP<Epetra_Vector>& disinterface_m);

    /// wear in sp conf.
    void WearSpatialMaster(Teuchos::RCP<Epetra_Vector>& disinterface_m);

    /// wear in sp conf. with mortar map
    void WearSpatialMasterMap(Teuchos::RCP<Epetra_Vector>& disinterface_s,
                              Teuchos::RCP<Epetra_Vector>& disinterface_m);

    /// wear in sp conf.
    void WearSpatialSlave(Teuchos::RCP<Epetra_Vector>& disinterface_s);

  protected:
    //nothing

  private:

    Teuchos::RCP<Epetra_Vector>                      wearnp_i_;     // wear in timestep n+1 and nonlin iter i
    Teuchos::RCP<Epetra_Vector>                      wearnp_ip_;    // wear in timestep n+1 and nonlin iter i+1
    Teuchos::RCP<Epetra_Vector>                      wearincr_;     // wear incr between wearnp_i_ and wearnp_ip_

    Teuchos::RCP<Epetra_Vector>                      delta_ale_;

    Teuchos::RCP< ::ADAPTER::Coupling>               coupalestru_;   // ale struct cpupling
    Teuchos::RCP< ::ADAPTER::Coupling>               coupstrualei_;  // ale struct coupling on ale interface

    Teuchos::ParameterList                           alepara_;       // ale parameter list


 };  // Algorithm

}  // namespace WEAR

/*----------------------------------------------------------------------*/
#endif  // WEAR_PARTITIONED_H
