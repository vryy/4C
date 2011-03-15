/*!
\file dof_distribution_switcher.H

\brief provides the dofmanager classes

<pre>
Maintainer: Axel Gerstenberger
            gerstenberger@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15236
</pre>
*/
#ifdef CCADISCRET
#ifndef DOF_DISTRIBUTION_SWITCHER_H
#define DOF_DISTRIBUTION_SWITCHER_H

#include "dofkey.H"
#include "../drt_cut/cut_boundingbox.H"
#include "../drt_lib/drt_exporter.H"

namespace XFEM
{
  class InterfaceHandle;

  /*!
  \brief this class is used to switch from one distribution to a second

  */
  class DofDistributionSwitcher
  {
  public:

    //! constructor
    explicit DofDistributionSwitcher(
        const Teuchos::RCP<InterfaceHandle>       ih,
        const Teuchos::RCP<DofManager>            dofman,
        const Epetra_Map&                         olddofrowmap,
        const Epetra_Map&                         newdofrowmap,
        const std::map<DofKey<onNode>, DofGID>&   oldNodalDofDistrib,
        const std::map<DofKey<onNode>, DofGID>&   newNodalDofDistrib,
        const std::map<DofKey<onElem>, DofGID>&   oldElementalDofDistrib,
        const std::map<DofKey<onElem>, DofGID>&   newElementalDofDistrib
    ) :
      ih_(ih),
      dofman_(dofman),
      olddofrowmap_(olddofrowmap),
      newdofrowmap_(newdofrowmap),
      oldNodalDofDistrib_(oldNodalDofDistrib),
      newNodalDofDistrib_(newNodalDofDistrib),
      oldElementalDofDistrib_(oldElementalDofDistrib),
      newElementalDofDistrib_(newElementalDofDistrib)
    {
      extractDofKeysForInitialization(unknownFieldEnr_);
      cout << "Number of nodes with missing nodal unknowns: " << unknownFieldEnr_.size() << endl;
      return;
    }

    //! constructor for FluidXFluid. Here all the nodes which were located in
    //xfem-fluid are unknown nodes.
    explicit DofDistributionSwitcher(
        const Teuchos::RCP<InterfaceHandle>       ih,
        const Teuchos::RCP<DofManager>            dofman,
        const Epetra_Map&                         olddofrowmap,
        const Epetra_Map&                         newdofrowmap,
        const std::map<DofKey<onNode>, DofGID>&   oldNodalDofDistrib,
        const std::map<DofKey<onNode>, DofGID>&   newNodalDofDistrib,
        const std::map<DofKey<onElem>, DofGID>&   oldElementalDofDistrib,
        const std::map<DofKey<onElem>, DofGID>&   newElementalDofDistrib,
        const std::map<int,int>                   nodelabeln
    ) :
      ih_(ih),
      dofman_(dofman),
      olddofrowmap_(olddofrowmap),
      newdofrowmap_(newdofrowmap),
      oldNodalDofDistrib_(oldNodalDofDistrib),
      newNodalDofDistrib_(newNodalDofDistrib),
      oldElementalDofDistrib_(oldElementalDofDistrib),
      newElementalDofDistrib_(newElementalDofDistrib)
    {
      extractDofKeysForInitializationFluidXFluid(unknownFieldEnr_,nodelabeln);
      cout << "Number of nodes with missing nodal unknowns: " << unknownFieldEnr_.size() << endl;
      return;
    }


    //! destructor
    ~DofDistributionSwitcher()
    {
      return;
    }


    //! switch from old to new dof distribution
    void mapVectorToNewDofDistribution(
        Teuchos::RCP<Epetra_Vector>&     vector
    ) const;

    //! switch from old to new dof distribution for jump/kink enrichment problems
    void mapVectorToNewDofDistributionCombust(
        Teuchos::RCP<Epetra_Vector>&     vector,
        const bool                       quasi_static_enr
    ) const;

    //! switch from old to new dof distribution for fluidxfluid Coupling
    void mapVectorToNewDofDistributionFluidXFluid(
      RCP<Epetra_Vector>&    vector,
      std::map<int,int>      nodelabeln
    ) const;

    void extractDofKeysForInitialization(
        std::map<int, std::set<FieldEnr> >& unknownFieldEnr
    ) const;

    void extractDofKeysForInitializationFluidXFluid(
      std::map<int, std::set<FieldEnr> >& unknownFieldEnr,
      std::map<int,int>      nodelabeln
    ) const;

    void extrapolateOldTimeStepValues(
        const Teuchos::RCP<DRT::Discretization>    bdis,
        const std::map<int,LINALG::Matrix<3,1> >&  cutterposnp,
        const Teuchos::RCP<const Epetra_Vector>    ivector,
        Teuchos::RCP<Epetra_Vector>                state_vector
    ) const;

    // for fluid-fluid Coupling
    // patchdis   (in): patchdis
    // patchboxes (in): map of element id of the fluid-patch element to it's patchbox
    void projectEmbeddedtoBackgroundfluid(
      const RCP<DRT::Discretization>                patchdis,
      const std::map<int,GEO::CUT::BoundingBox>&    patchboxes,
      const RCP<const Epetra_Vector>                fluidstate_vector_n,
      const RCP<Epetra_Vector>                      xfluidstate_vector_n,
      const RCP<const Epetra_Vector>                fluiddispnm,
      int                                           step
    ) const;

    // Transforms a backgroundFluid node in current coordinates to a point of patch fluid
    // in element coordinates due to the given patch-element. After this the
    // value of the node is interpolated from the embedded fluid.
    // The nonlinar system is: F(xsi)=X_patch(xsi)-X_xfem=N(xsi)X_I_patch-X_xfem=0
    // This nonlinear system of equation is solved with help of the Newton-method.
    // The following linear system of equation is solved in every Newton-step:
    // J(xsi)*delta_xsi=-F(xsi)
    // pele                 (in)   :patch_element
    // patchdis             (in)   :embedded
    // x                    (in)   :backgr_node in current coordinates (x, y, z)
    // xsi                  (inout):node in element coordinates
    bool ComputeSpacialToElementCoordAndProject(
      const DRT::Element*                                 pele, //patch_element
      const RCP<DRT::Discretization>                      patchdis, // patch dis
      const LINALG::Matrix<3,1>&                          x, // background node's coordinates (x,y,z)
      const RCP<const Epetra_Vector>                      fluidstate_vector_n, //source
      const RCP<Epetra_Vector>                            xfluidstate_vector_n, //the vector with unknowns
      std::map<int, set<XFEM::FieldEnr> >::const_iterator unknownnodaldof,
      const RCP<const Epetra_Vector>                      fluiddispnm  // embedded fluid position
    ) const;

    //! switch from old to new dof distribution
    void generateTransferInformation(
        const Teuchos::RCP<Epetra_Vector>&     vector
    ) const;

    // needed for Round Robin
    void PackValuestoSend(
      const RCP<Epetra_Vector>                   state_vector,
      std::map<int, set<XFEM::FieldEnr> >        unknownFieldEnr,
      DRT::PackBuffer                          & sblock
    ) const;

    void UnpackValues(
      RCP<Epetra_Vector>                         state_vector,
      std::map<int, set<XFEM::FieldEnr> >         unknownFieldEnr,
      vector<char>                             &  rblock
    ) const;

#ifdef PARALLEL
    //! receive a block in the round robin communication pattern
    void ReceiveBlock(
      vector<char>                              & rblock,
      DRT::Exporter                             & exporter,
      MPI_Request                               & request
    ) const;
#endif

#ifdef PARALLEL
    //! send a block in the round robin communication pattern
    void SendBlock(
      vector<char>                             & sblock,
      DRT::Exporter                            & exporter,
      MPI_Request                              & request
    ) const;
#endif

    void GmshOutput(
      int                               step,
      const RCP<DRT::Discretization>    patchdis,
      std::map<int,int>                mapforoutput
    ) const;

  private:

    //! disabled copy constructor
    explicit DofDistributionSwitcher(const DofDistributionSwitcher& dofdist);

    //! disabled assignment operator
    DofDistributionSwitcher operator = (const DofDistributionSwitcher& old);

    //! disabled default constructor
    explicit DofDistributionSwitcher();

    const Teuchos::RCP<InterfaceHandle>       ih_;                      ///< interface information
    const Teuchos::RCP<DofManager>            dofman_;                  ///< dofmanager for this proc
    const Epetra_Map                          olddofrowmap_;            ///< dofrowmap of the old state vector
    const Epetra_Map                          newdofrowmap_;            ///< dofrowmap of the new state vector
    const std::map<DofKey<onNode>, DofGID>    oldNodalDofDistrib_;      ///< old positions of dof in Epetra Vector per nodal dofkey
    const std::map<DofKey<onNode>, DofGID>    newNodalDofDistrib_;      ///< new positions of dof in Epetra Vector per nodal dofkey
    const std::map<DofKey<onElem>, DofGID>    oldElementalDofDistrib_;  ///< old positions of dof in Epetra Vector per element dofkey
    const std::map<DofKey<onElem>, DofGID>    newElementalDofDistrib_;  ///< new positions of dof in Epetra Vector per element dofkey

    std::map<int, std::set<FieldEnr> >  unknownFieldEnr_;
  };

  //! try to find another enrichment for this physical field
  Enrichment genAlternativeEnrichment(
      const int                             gnodeid,
      const PHYSICS::Field            oldphysvar,
      const Teuchos::RCP<DofManager>& dofman
  );


} // namespace XFEM

#endif
#endif  // #ifdef CCADISCRET
