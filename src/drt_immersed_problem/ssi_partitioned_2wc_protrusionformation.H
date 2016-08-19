/*!----------------------------------------------------------------------
\file ssi_partitioned_2wc_protrusionformation.H

\brief specialization of ssi2wc, including "structale"-surface growth

\level 3

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

//! specialized class to deal with partitioned 2WC SSI during protrusion formation in cell migration.
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
    \brief match conditioned nodes of master and slave fields and fill maps

    \param masterdis             (in)  : master discretization conditioned with 'condition' whose nodes should be matched
    \param slavedis              (in)  : slave discretization conditioned with 'condition' that is matched to the masterdis
    \param masternodemap_matched (out) : map containing ids of matched master nodes
    \param slavenodemap_matched  (out) : map containing ids of matched slave nodes
    \param condname              (in)  : condition to be matched

    \remark The condition has to build a geometry.

    \note Currently it is checked whether 4 noded elements are conditioned. This is just for current debugging purposes.
          This check can be removed any time and the method works for general purposes.

   */
  void MatchNodes(
      DRT::Discretization& masterdis,
      DRT::Discretization& slavedis,
      Teuchos::RCP<Epetra_Map>& masternodemap_matched,
      Teuchos::RCP<Epetra_Map>& slavenodemap_matched,
      const std::string& condname);

  //! match nodes of master and slave field that have specified condition
  void MatchConditionedNodes(
      DRT::Discretization& masterdis,
      DRT::Discretization& slavedis,
      Teuchos::RCP<Epetra_Map>& masternodemap_matched,
      Teuchos::RCP<Epetra_Map>& slavenodemap_matched,
      const std::string& condname);

  /*!
    \brief build map associating master elements to matching slaveelements

    \param masterdis             (in)  : master discretization conditioned with 'condition' whose nodes should be matched
    \param slavedis              (in)  : slave discretization conditioned with 'condition' that is matched to the masterdis
    \param masternodemap_matched (in)  : map containing ids of matched master nodes
    \param slavenodemap_matched  (in)  : map containing ids of matched slave nodes
    \param maptofill             (out) : map that is filled master id (key) -> slave id (value)
    \param condname              (in)  : condition to be matched

    \remark The condition has to build a geometry.

    \note Currently it is checked whether 4-node elements are conditioned. This is just for current debugging purposes.
          This check can be removed any time and the method works for general purposes.

   */
  void BuildMasterGeometryToSlaveDisEleMap(
      const DRT::Discretization& masterdis,
      const DRT::Discretization& slavedis,
      const Teuchos::RCP<Epetra_Map>& masternodemap_matched,
      const Teuchos::RCP<Epetra_Map>& slavenodemap_matched,
      Teuchos::RCP<std::map<int,int> > maptofill,
      const std::string& condname);

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

  //! setup discretizations and dofsets
  virtual void SetupDiscretizations(const Epetra_Comm& comm, const std::string struct_disname, const std::string scatra_disname);

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
  Teuchos::RCP< ::ADAPTER::Coupling> coupstrualei_;   //!< ale struct coupling on ale interface
  Teuchos::RCP< ::ADAPTER::Coupling> coupscatrastrui_;//!< scatra struct coupling on whole interface
  //@}

  //! underlying ALE
  Teuchos::RCP<ADAPTER::AleFluidWrapper> ale_;

  //! underlying structure of the SSI problem
  Teuchos::RCP< ::ADAPTER::FSIStructureWrapper> specialized_structure_;

  //! map associating structure to scatra element ids
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
