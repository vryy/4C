/*!----------------------------------------------------------------------
\file ssi_partitioned_2wc_protrusionformation.H

\brief specialization of ssi2wc, including "structale"-surface growth

\level 2

\maintainer Andreas Rauch
            rauch@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289 -15240

*----------------------------------------------------------------------*/
#ifndef SSI_PARTITIONED2_WC_PROTRUSIONFORMATION_H_
#define SSI_PARTITIONED2_WC_PROTRUSIONFORMATION_H_

#include "../drt_ssi/ssi_partitioned_2wc.H"

/*==========================================================================*/
// forward declarations
/*==========================================================================*/
namespace ADAPTER
{
  class FSIStructureWrapper;
  class AleFluidWrapper;
  class Coupling;
}

namespace DRT
{
  class Element;
  class ImmersedFieldExchangeManager;
}


namespace SSI
{

/*!
\brief Specialized class to deal with partitioned 2WC SSI during protrusion formation in cell migration.

Protrusion formation, means actin polymerization in the cytosol, pushing the cell membrane forward. This
is a Biochemo-Mechano coupled model, meaning, that there is a biochemical model and a mechanical model
acting at the same time. These two submodels have to be in equilibrium. This equilibrium is established by
a two way coupled partitioned Structure-Scalar Interaction (SSI) approach, implemented in this class.

<h3>Biochemical Part</h3>

The biochemical model consists of reaction networks in the cytosol. They can be defined in the input file.
Note, that there are a volume scatra discretization, modeling transport in the cytosol, AND a surface scatra
discretization modeling transport at the cell membrane. We may also model binding/unbinding at the membrane
by reactions in the input file. This is important because some species, like the Rho-GTPases cycle between
inactivated cytosol-bound and activated membrane-bound form.
Furthermore, we model source/sink terms at the transport boundary of the cytosolic volume. Note, that actin
polymerization is happening in the volume, however directly underneath the membrane. Actin polymerization
introduces a sink on the actin monomer concentration. At the same time new branches may form. This is modeled by
a source on the branch concentration. Branches start at the side of actin filaments. The branch concentration
determines to how many actin filaments the polymerized actin monomers are distributed. At a fixed amount of polymerized
actin monomers in a given time step, the more branches there are, the lower the growth value is, because the actin
monomers a distributed to more filaments. Each actin monomer polymerizing below the membrane, may displace the membrane
by a value, proportional to the actin monomer size.

<h3>Mechanical Part</h3>

<h3>Implementation</h3>

For the implementation of the actual growth equation and the corresponding source/sink terms, it is important to note,
that this class relies on the \ref ssicoupling_ object in SSI::SSI_Base::SSI_Base. In this coupling object we build
\ref DRT::DofSetMappedProxy::DofSetMappedProxy . See also the documentation for this class. The cell discretization
wears three dofsets: 1) the cell structural dofs, 2) the scatra volume dofs, and 3) the scatra surface dofs. By
scatra surface we mean actual surface elements that are given in the input file. The coupling between the surface
discretization and the volume discretization (actually they are both combined in one large discretization) is managed
in \ref SCATRA::HeterogeneousReactionStrategy::HeterogeneousReactionStrategy .


<h3>Usage</h3>

For an example on how to use this model see test case 'immersed_cell_biochemo_mechano_pureProtrusion_h8.dat'


<h3>Result Filter</h3>

If you filter your results using post_drt_ensight, the maximum dof number of all nodes in the Transport discretization is determined. This gives
the amount of phis which are written. If your membrane wears e.g. 2 dofs and the cell volume only 1, fields phi_1 and phi_2 are written.
For the cell volume phi_2 then has no meaning, i.e. it should just be zero everywhere in the domain.

You can split cell membrane discretization and volume discretization by applying the Connectivity Filter to the scatra.case. If you apply the Threshold
filter subsequently, you get either the volume or the surface depending on your threshold value (0 or 1).

\date 08/2016

\author rauch (rauch@lnm.mw.tum.de)
*/
class SSI_Part2WC_PROTRUSIONFORMATION : public SSI_Part2WC
{
public:

  //! constructor (ctor)
  SSI_Part2WC_PROTRUSIONFORMATION(const Epetra_Comm& comm,
      const Teuchos::ParameterList& globaltimeparams);

  //! virtual destructor
  virtual ~SSI_Part2WC_PROTRUSIONFORMATION(){};

  /*!
  \brief Setup this object

   Initializes members and performs problem specific setup.

  \note Must only be called after parallel (re-)distribution of discretizations is finished !
        Otherwise, vectors may have wrong maps.

  \warning none
  \return void
  \date 08/16
  \author rauch
  */
  virtual void Setup(const Epetra_Comm& comm,
      const Teuchos::ParameterList& globaltimeparams,
      const Teuchos::ParameterList& scatraparams,
      const Teuchos::ParameterList& structparams,
      const std::string struct_disname,
      const std::string scatra_disname);

  /*!
  \brief Setup discretizations and dofsets

   Fill discretizations and setup coupling object \ref ssicoupling_

  \note Must only be called after parallel (re-)distribution of discretizations is finished !
        Otherwise, vectors may have wrong maps.

  \warning Parallel redistribution inside !

  \return void
  \date 08/16
  \author rauch
  */
  virtual void SetupDiscretizationsAndFieldCoupling(
      const Epetra_Comm& comm,
      const std::string& struct_disname,
      const std::string& scatra_disname);

  //! access to ale field
  const Teuchos::RCP<ADAPTER::AleFluidWrapper> AleField(){return ale_;};

  //! update time step and print to screen
  virtual void UpdateAndOutput();


protected:

  //! perform iteration step of ale field
  void DoAleStep(Teuchos::RCP<Epetra_Vector> growthincrement);

  //! solve the structure partition and set the new states in the scatra field
  virtual void DoStructStep();

  //! solve the scatra partition and set the new state in the structure field
  virtual void DoScatraStep();

  //! access to scatra field
  const Teuchos::RCP< ::ADAPTER::FSIStructureWrapper>&  SpecStructureField(){return specialized_structure_;};

  //! return coupling object
  Teuchos::RCP< ::ADAPTER::Coupling> AleStruCoupling(){return coupalestru_; };

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

  //! evaluate sources with respect to growth
  void EvaluateSources();

  //! evaluate cell growth with respect to sources
  void EvaluateGrowth();

  //! update material configuration
  void UpdateMatConf();

  //! update spatial configuration
  void UpdateSpatConf();

  /*!
    \brief build map containing row dofs of conditioned nodes

    \param condition  (in)  : condition to be considered
    \param dis        (in)  : discretization wearing the condition
    \param conddofmap (out) : dofrowmap for conditioned nodes

   */
  void BuildConditionDofRowMap(const DRT::Condition* condition,
                               const Teuchos::RCP<const DRT::Discretization> dis,
                               Teuchos::RCP<Epetra_Map>& conddofmap);

  //! transform from structure to ale dofs
  virtual Teuchos::RCP<Epetra_Vector> StructureToAle(Teuchos::RCP<Epetra_Vector> vec);

  //! transform from structure to ale dofs
  virtual Teuchos::RCP<Epetra_Vector> StructureToAle(Teuchos::RCP<const Epetra_Vector> vec);

  //! transform from ale to structure dofs
  virtual Teuchos::RCP<Epetra_Vector> AleToStructure(Teuchos::RCP<Epetra_Vector> vec);

  //! transform from ale to structure dofs
  virtual Teuchos::RCP<Epetra_Vector> AleToStructure(Teuchos::RCP<const Epetra_Vector> vec);

  //! convergence check of outer loop
  virtual bool ConvergenceCheck(int itnum);

  //! update the current states in every iteration
  //! states are set to the last solutions obtained
  virtual void IterUpdateStates();


private:

  //! number of dofset associated with displacement dofs
  int nds_disp_;

  //! relaxation parameter
  double omega_;

  //! @name Pointers to Discretizations
  //@{
  Teuchos::RCP<DRT::Discretization> backgroundstructuredis_; //!< background ecm discretisation
  Teuchos::RCP<DRT::Discretization> backgroundfluiddis_;     //!< background interstitial flow discretisation
  Teuchos::RCP<DRT::Discretization> celldis_;                //!< immersed cell discretisation
  Teuchos::RCP<DRT::Discretization> cellscatradis_;          //!< immersed cell scatra discretisation
  //@}

  //! @name Data Management
  //@{
  DRT::ImmersedFieldExchangeManager* exchange_manager_ ; //!< manager for data exchange between fields
  //@}

  //! @name Information for parallelism
  //@{
  int myrank_; //!< processor id
  //@}

  //! at each scatra node the numdof_actin_'th dof is the actin dof
  int numdof_actin_;

  //! growth state
  Teuchos::RCP<Epetra_Vector> x_;

  //! growth at t_n
  Teuchos::RCP<Epetra_Vector> growth_n_;

  //! growth at t_{n+1}
  Teuchos::RCP<Epetra_Vector> growth_np_;

  //! growth step increment (increment containing growth from t_n -> t_{n+1})
  Teuchos::RCP<Epetra_Vector> growth_step_;

  //! growth iteration increment of the growth loop (increment containing growth from iteration i -> i+1)
  Teuchos::RCP<Epetra_Vector> growthinc_;

  //! source values for actin monomers, branches, and filament barbed ends (to be added to scatra rhs)
  Teuchos::RCP<Epetra_Vector> sources_;

  //! yet unknown vector
  Teuchos::RCP<Epetra_Vector> delta_ale_;

  //! @name coupling objects declaration
  //@{
  Teuchos::RCP< ::ADAPTER::Coupling> coupalestru_;    //!< ale struct coupling
  //@}

  //! underlying ALE
  Teuchos::RCP<ADAPTER::AleFluidWrapper> ale_;

  //! underlying structure of the SSI problem
  Teuchos::RCP< ::ADAPTER::FSIStructureWrapper> specialized_structure_;

  //! map associating structure to scatra (bdry) element ids
  Teuchos::RCP<std::map<int,int> > structscatraelemap_;

  //! map containing the dofs of conditioned growth nodes
  Teuchos::RCP<Epetra_Map> conditiondofrowmap_;

  //! matched master nodemap (with coupled node IDs)
  Teuchos::RCP<Epetra_Map> structnodemap_;

  //! matched slave nodemap (with coupled node IDs)
  Teuchos::RCP<Epetra_Map> scatranodemap_;

}; // class SSI_Part2WC_PROTRUSIONFORMATION


} // namespace SSI

#endif /* SSI_PARTITIONED2_WC__PROTRUSIONFORMATION_H_ */
