/*----------------------------------------------------------------------*/
/*!
\file wear_partitioned.H

\brief  Basis of partitioned wear algorithm
        (Lagrangian step followed by Eulerian step)

\level 2

\maintainer Alexander Popp
*/

/*----------------------------------------------------------------------*
 | definitions                                              farah 11/13 |
 *----------------------------------------------------------------------*/

#ifndef WEAR_PARTITIONED_H
#define WEAR_PARTITIONED_H

/*----------------------------------------------------------------------*
 | headers                                                  farah 11/13 |
 *----------------------------------------------------------------------*/
#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_StandardParameterEntryValidators.hpp>

#include "wear_algorithm.H"

#include "../drt_adapter/ad_str_fsiwrapper.H"

/*----------------------------------------------------------------------*
 | forward declarations                                     farah 11/13 |
 *----------------------------------------------------------------------*/
namespace MORTAR
{
  class ManagerBase;
}
namespace ADAPTER
{
  class CouplingBase;
  class Coupling;
}
namespace DRT
{
  class LocationArray;
  class Element;
}
namespace ALE
{
  class Ale;
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

    // write output for ale and structure
    void Output();

    /// read restart data
    virtual void ReadRestart(
      int step  ///< step number where the calculation is continued
      );

    /// general time loop
    virtual void TimeLoop();

    // update ale and structure
    void Update();

  protected:
    //nothing

  private:
    // do ale step
    void AleStep(Teuchos::RCP<Epetra_Vector> idisale_global);

    // transform from ale to structure dofs
    virtual Teuchos::RCP<Epetra_Vector> AleToStructure(Teuchos::RCP<Epetra_Vector> vec) const;

    // transform from ale to structure dofs
    virtual Teuchos::RCP<Epetra_Vector> AleToStructure(Teuchos::RCP<const Epetra_Vector> vec) const;

    /// Application of mesh displacements (frictional contact) to material conf
    void UpdateMatConf();

    /*!
    \brief parameter space mapping between configurations

    \param Xtarget           (out) : new material coordinate of considered node (spatialtomaterial=true), or
                                     new spatial  coordinate of considered node (spatialtomaterial=false)
    \param Xsource           (in)  : current mesh coordinate of configuration that is NOT supposed to be transported
    \param ElementPtr        (in)  : pointer to elements adjacent to considered node
    \param numelements       (in)  : number of elements adjacent to considered node
    \param spatialtomaterial (in)  : true if considered node already has correct spatial coordinate
                                     and corresponding material coordinate is to be determined
     */
    void AdvectionMap(
        double* Xtarget,            // out
        double* Xsource,            // in
        DRT::Element** ElementPtr,  // in
        int numelements,            // in
        bool spatialtomaterial);    // in

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

    // ale parameter list
    virtual Teuchos::ParameterList& ParamsAle() {return alepara_;};

    // prepare time step for ale and structure
    void PrepareTimeStep();

    // redistribute material interfaces according to current interfaces
    void RedistributeMatInterfaces();

    // transform from ale to structure dofs
    virtual Teuchos::RCP<Epetra_Vector> StructureToAle(Teuchos::RCP<Epetra_Vector> vec) const;

    // transform from ale to structure dofs
    virtual Teuchos::RCP<Epetra_Vector> StructureToAle(Teuchos::RCP<const Epetra_Vector> vec) const;

    // time loop for staggered coupling
    void TimeLoopStagg(bool alestep);

    //time loop for iterative stagered coupling
    void TimeLoopIterStagg();

    // update spatial displacements due to mat. displ
    void UpdateSpatConf();

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

    Teuchos::RCP<Epetra_Vector>                      wearnp_i_;     // wear in timestep n+1 and nonlin iter i
    Teuchos::RCP<Epetra_Vector>                      wearnp_ip_;    // wear in timestep n+1 and nonlin iter i+1
    Teuchos::RCP<Epetra_Vector>                      wearincr_;     // wear incr between wearnp_i_ and wearnp_ip_

    Teuchos::RCP<Epetra_Vector>                      delta_ale_;
    Teuchos::RCP<Epetra_Vector>                      ale_i_;

    Teuchos::RCP< ::ADAPTER::Coupling>               coupstrualei_;  // ale struct coupling on ale interface
    Teuchos::RCP< ::ADAPTER::CouplingBase>           coupalestru_;   // ale struct cpupling

    Teuchos::ParameterList                           alepara_;       // ale parameter list

 };  // Algorithm

}  // namespace WEAR

/*----------------------------------------------------------------------*/
#endif  // WEAR_PARTITIONED_H
