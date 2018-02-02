/*-----------------------------------------------------------*/
/*!
\file beaminteraction_submodel_evaluator_crosslinking.H

\brief class for submodel crosslinking

\maintainer Jonas Eichinger, Maximilian Grill

\level 3

*/
/*-----------------------------------------------------------*/


#ifndef BEAMINTERACTION_SUBMODEL_EVALUATOR_CROSSLINKING_H_
#define BEAMINTERACTION_SUBMODEL_EVALUATOR_CROSSLINKING_H_

#include "../drt_binstrategy/binning_strategy_utils.H"
#include "../drt_inpar/inpar_beaminteraction.H"
#include "../linalg/linalg_fixedsizematrix.H"
#include <Epetra_MpiComm.h>
#include "beaminteraction_submodel_evaluator_generic.H"


// forward declaration
class DiscretizationRuntimeVtpWriter;

namespace MAT {
  class CrosslinkerMat;
}
namespace DRT {
  class Exporter;
  class Element;
  class PackBuffer;
  class Node;
}
namespace CROSSLINKING {
  class CrosslinkerNode;
}
namespace BEAMINTERACTION {
namespace DATA {
  struct CrosslinkerData;
  struct BeamData;
  struct BindEventData;
  struct UnBindEventData;
}
class CrosslinkingParams;
class BeamLink;

namespace SUBMODELEVALUATOR {

  class Crosslinking : public Generic
  {
  public:

    //! constructor
    Crosslinking();

    //! destructor
    virtual ~Crosslinking() {};

    //! setup class variables
    virtual void Setup();

    //! derive
    virtual bool PostPartitionProblem();

    //! derive
    virtual void PostSetup();

    //! Returns the type of the current model evaluator
    INPAR::BEAMINTERACTION::SubModelType Type() const
    { return INPAR::BEAMINTERACTION::submodel_crosslinking; }

    //! derived
    virtual void Reset();

    //! derived
    virtual bool EvaluateForce();

    //! derived
    virtual bool EvaluateStiff();

    //! derived
    virtual bool EvaluateForceStiff();

    //! derived
    virtual void UpdateStepState(const double& timefac_n);

    //! derived
    virtual bool PreUpdateStepElement( bool beam_redist );

    //! derived
    virtual void UpdateStepElement( bool repartition_was_done );

    //! derived
    virtual void PostUpdateStepElement();

    //! derived
    virtual void OutputStepState(IO::DiscretizationWriter& iowriter) const;

    //! derived
    virtual void RuntimeOutputStepState() const;

    //! derived
    virtual void ResetStepState();

    //! derived
    virtual void WriteRestart(
        IO::DiscretizationWriter & ia_writer,
        IO::DiscretizationWriter & bin_writer) const;

    //! derived
    virtual void ReadRestart(
        IO::DiscretizationReader & ia_reader,
        IO::DiscretizationReader & bin_reader);

    //! derived
    virtual void PostReadRestart();

    //! derived
    virtual void RunPostIterate( const NOX::Solver::Generic& solver )
    { /*empty*/ };

    //! derived
    virtual void InitSubmodelDependencies(Teuchos::RCP<STR::MODELEVALUATOR::BeamInteraction::Map> const submodelmap);

    //! derived
    virtual void AddBinsToBinColMap( std::set< int >& colbins);

    //! derived
    virtual void AddBinsWithRelevantContentForIaDiscretColMap( std::set< int >& colbins) const;

    //! derived
    virtual void GetHalfInteractionDistance( double & half_interaction_distance );

    //! @}

    //!@name routines that are not derived and handle crosslinking
    //! @{

    /// unbind all crosslinker residing in bingid and its neighborhood
    virtual void UnbindCrosslinkerInBinsAndNeighborhood(
        std::set<int> const& bingids,
        bool const doubleunbind);

    /// unbind all crosslinker residing in bingid and its neighborhood
    virtual void UnbindCrosslinkerInBinsAndNeighborhood(
        std::set<int> const& bingids);

    /// determine which proc is responsible for forced crosslinker unbinding in certain bins
    /// by checking bin ownership
    virtual void DetermineResponsilbeProcsForForcedCrosslinkerUnbinding(
        std::set<int> const& bingids,
        std::set< int >& binsonmyrank) const;

    virtual void CommunicateBinIds(
        std::map< int, std::vector< int > > const& binstosend,
        std::set< int >& binsonmyrank) const;

    /// set all double bonds in bind and neighborhood
    virtual void DoubleBindCrosslinkerInBinsAndNeighborhood(
        std::set<int> const& bingids);

    //! @}

private:

    //! writes output for discretization structure in VTP format
    void WriteOutputRuntimeVtpStructure() const;

    //! init output for discretization structure in VTP format
    void InitOutputRuntimeVtpStructure();

    //! @name Small data structs for faster access and better structure
    //! @{

    //! struct that stores all necessary data to handle the crosslinking between two elements on each proc
    struct NewDoubleBonds {
        int                                   id;          // gid of crosslinker
        std::vector< std::pair<int, int> >    eleids;      // elegid and local binding spot number of first element
        std::vector< LINALG::Matrix< 3, 1 > > bspotposs;
        std::vector< LINALG::Matrix< 3, 3 > > bspottriads;
      };

    //! @}


    //!@name routines that are not derived and handle crosslinking
    //! @{

    /// add crosslinker to bin discretization initially
    virtual void AddCrosslinkerToBinDiscretization();

    /// set binding spot positions on filament initially
    virtual void SetFilamentBindingSpotPositions();

    /// set filament types on elements
    virtual void SetFilamentTypes();

    /// extent ghosting so that each proc has all nodes and elemens of a filament containing
    /// at least one of its nodes
    virtual void ExtendGhostingForFilamentBspotSetup(
        std::set<int>& relevantfilaments);

    /// determine nodes with relevant element cloud
    virtual void DetermineOffMyRankNodesWithRelevantEleCloudForFilamentBspotSetup(
        std::set<int>& examinedfilaments,
        std::set<int>& setofrequirednodes) const;

    // compute total filament length using its elements reference length and sort elements
    virtual void ComputeFilamentLengthAndSortItsElements(
        std::vector< DRT::Element* >& sortedfilamenteles,
        std::vector<int> const* nodeids,
        double & filreflength) const;

    /// set binding spot positions on element of a filament
    virtual void SetBindingSpotsPositionsOnFilament(
        std::vector< DRT::Element* >& sortedfilamenteles,
        double const start,
        int const    numbspot,
        double const filamentbspotinterval,
        double const tol);


    /// set double bonded linker between all binding spots that match certain
    /// neigbhoring criteria
    virtual void SetAllPossibleInitialDoubleBondedCrosslinker(
        std::vector< Teuchos::RCP<BEAMINTERACTION::DATA::CrosslinkerData> > & newlinker,
        std::map< int, NewDoubleBonds > & mynewdbondcl);

    /// diffuse crosslinker depending on number of bonds they have
    virtual void DiffuseCrosslinker();

    /// diffuse unbound crosslinker according to brownian dynamics
    virtual void DiffuseUnboundCrosslinker(
        DRT::Node * crosslinker,
        BEAMINTERACTION::DATA::CrosslinkerData * cldata_i);

    /// get binding spot of crosslinker that is currently occupied
    virtual int GetSingleOccupiedClBspot(
        std::vector<std::pair<int, int> > const& clbspots
    ) const;

    void SetPositionOfDoubleBondedCrosslinkerPBCconsistent(
        LINALG::Matrix<3,1>& clpos,
        LINALG::Matrix<3,1> const& bspot1pos,
        LINALG::Matrix<3,1> const& bspot2pos
    ) const;

    /// new position after transition from single to not bonded
    virtual void SetPositionOfNewlyFreeCrosslinker(
        CROSSLINKING::CrosslinkerNode* crosslinker,
        BEAMINTERACTION::DATA::CrosslinkerData * cldata
    );

    /// new position after transition from double to single bonded
    virtual void  SetPositionOfNewlySingleBondedCrosslinker(
        DRT::Node* crosslinker,
        BEAMINTERACTION::DATA::CrosslinkerData * cldata,
        int const stayoccpotid);

    /// fill epetar vectors to write vtp output
    void FillStateDataVectorsForOutput(
        Teuchos::RCP<Epetra_Vector> displacement,
        Teuchos::RCP<Epetra_Vector> orientation,
        Teuchos::RCP<Epetra_Vector> numberofbonds,
        Teuchos::RCP<Epetra_Vector> owner,
        Teuchos::RCP<Epetra_Vector> force,
        Teuchos::RCP<Epetra_Vector> force_norm
    ) const;

    /// update maps
    virtual void StoreMapsPriorRedistribution();

    /// get crosslink data before interaction evaluation
    virtual void UpdateAndExportCrosslinkerData();

    /// get beam data before interaction evaluation
    virtual void UpdateAndExportBeamData( bool update_states = true );

    /// bind and unbind crosslinker
    virtual void BindAndUnbindCrosslinker();

    /// bind crossslinker
    virtual void BindCrosslinker();

    /// search and set crosslinker
    /*! -------------------------------------------------------------------------
     note: only the owner of a beam element is allowed to change the status of
     of a binding spot. Therefore we utilize the one layer ghosting around bins
     containing a crosslinker and the ghosting around bins that are touched
     by a row element (this can lead to two layer ghosting) of a proc. Thus we
     exclude the binding of two crosslinker on different procs on the same
     binding spot without loosing any potential interaction.
     To ensure that no crosslinker is bonded to often but still totally random over
     all procs, each binding event of a col crosslinker to a row element needs to
     be communicated to the crosslinker owner, he randomly decides who is allowed
     to bind, sets the according stuff for the cl and  informs back the
     requesting procs so they can set the stuff for the elements.
     As no proc on his own can decide whether a crosslink should be set, two
     binding events for one crosslinker in one time step are excluded (for this
     the proc must be sure that a crosslink is set as the binding range is
     different for a single bonded crosslinker compared to a free one)
    *  \author J. Eichinger
     -------------------------------------------------------------------------*/
    virtual void FindPotentialBindingEvents(
        std::map<int, Teuchos::RCP< BEAMINTERACTION::DATA::BindEventData > > & mybonds,
        std::map<int, std::vector< Teuchos::RCP< BEAMINTERACTION::DATA::BindEventData > > > & undecidedbonds
    );

    /// find potential binding events in one bin
    virtual void FindPotentialBindingEventsInBinAndNeighborhood(
        DRT::Element* bin,
        std::map<int, Teuchos::RCP< BEAMINTERACTION::DATA::BindEventData > > & mybonds,
        std::map<int, std::vector< Teuchos::RCP< BEAMINTERACTION::DATA::BindEventData > > >& undecidedbonds,
        std::vector< std::map< int, std::set<int> > >& intendedbeambonds,
        bool const checklinkingprop
    );

    /// search for binding events on each proc separately (i.e. pretending myrank is alone)
    /// communication to ensure correct binding over all procs is done afterwards
    virtual void PrepareBinding(
        DRT::Node * node_i,
        std::set<DRT::Element*> const & neighboring_beams,
        std::map< int, Teuchos::RCP< BEAMINTERACTION::DATA::BindEventData > > & mybonds,
        std::map< int, std::vector< Teuchos::RCP< BEAMINTERACTION::DATA::BindEventData > > > & undecidedbonds,
        std::vector< std::map< int, std::set<int> > > & intendedbeambonds,
        bool const checklinkingprop
    );

    /// check criteria if binding event is feasible
    virtual bool CheckBindEventCriteria(
        CROSSLINKING::CrosslinkerNode const * const crosslinker_i,
        DRT::Element const * const potbeampartner,
        BEAMINTERACTION::DATA::CrosslinkerData * cldata_i,
        BEAMINTERACTION::DATA::BeamData const * beamdata_i,
        int locnbspot,
        std::vector< std::map< int, std::set<int> > >& intendedbeambonds,
        bool checklinkingprop
    ) const;

    // check if identical bond alread exists
    virtual bool ReturnFalseIfIdenticalBondAlreadyExists(
        BEAMINTERACTION::DATA::CrosslinkerData * cldata_i,
        std::vector< std::map< int, std::set<int> > > & intendedbeambonds,
        BEAMINTERACTION::DATA::BeamData const * beamdata_i,
        int locnbspot,
        int potbeampartnerrowlid
    ) const;

    /// check if crosslinke and filament type are compatible
    virtual bool CheckLinkerAndFilamentTypeCompatibility(
        INPAR::BEAMINTERACTION::CrosslinkerType linkertype,
        INPAR::BEAMINTERACTION::FilamentType filamenttype) const;

    /// if crosslinker is singly bound, we fetch the orientation vector of the
    /// filament axis at the already occupied binding spot for the orientation
    /// criterion (enclosed angle) to be checked later on
    virtual void GetOccupiedClBSpotBeamTangent(
        BEAMINTERACTION::DATA::CrosslinkerData * cldata_i,
        LINALG::Matrix<3,1>& occ_bindingspot_beam_tangent,
        int const clgid
    ) const;

    /// decide by asking other procs who is allowed to set specific crosslinker,
    /// this is necessary to avoid setting crosslinker more than once per time step
    virtual void ManageBindingInParallel(
        std::map< int, Teuchos::RCP< BEAMINTERACTION::DATA::BindEventData> > & mysetbonds,  // clgid to cldata
        std::map< int, std::vector< Teuchos::RCP< BEAMINTERACTION::DATA::BindEventData > > > & undecidedbonds,  // owner of cldatas in vector to be requested
        std::map< int, Teuchos::RCP< BEAMINTERACTION::DATA::BindEventData > > & myelebonds
    ) const;

    /// communicate requests
    virtual void CommunicateUndecidedBonds(
        std::map< int, std::vector< Teuchos::RCP< BEAMINTERACTION::DATA::BindEventData > > > & undecidedbonds,
        std::map< int, std::vector< Teuchos::RCP< BEAMINTERACTION::DATA::BindEventData > > > & requestedcl
    ) const;

    /*
     now myrank needs to decide which proc is allowed to set the requested
     link, add it to his own list as row owner of cl sets stuff for cls, send
     back the answers to the row ele owner and receive the decisions made for
     its own requests:
     - if only one proc is requesting, the link can be set
     - if two procs are requesting or the current proc wants to set a link with
       a requested crosslinker, a random decision who is allowed to set the link
       has to be made
     -------------------------------------------------------------------------*/
    virtual void DecideBindingInParallel(
        std::map< int, std::vector< Teuchos::RCP< BEAMINTERACTION::DATA::BindEventData > > >& requestedcl,
        std::map< int, Teuchos::RCP< BEAMINTERACTION::DATA::BindEventData > >&              mybonds,
        std::map< int, std::vector< Teuchos::RCP< BEAMINTERACTION::DATA::BindEventData > > >& decidedbonds
    ) const;

    /// communicate decisions for binding events
    virtual void CommunicateDecidedBonds(
        std::map< int, std::vector< Teuchos::RCP< BEAMINTERACTION::DATA::BindEventData > > >& decidedbonds,
        std::map< int, Teuchos::RCP< BEAMINTERACTION::DATA::BindEventData >  >&              myelebonds
    ) const;

    /* now have two distinct maps of binding events on each proc, depending
       on ownership of crosslinker and elements myrank has different tasks:
        - mybonds: myrank takes care of crosslinker and (most) elements
        - myelebonds: myrank takes care of elements
       within those maps, different treatment is necessary for free and single
       bonded linker
                                                                            */
    virtual void UpdateMyCrosslinkerAndElementBindingStates(
        std::map< int, Teuchos::RCP< BEAMINTERACTION::DATA::BindEventData > > & mybonds,
        std::map< int, Teuchos::RCP< BEAMINTERACTION::DATA::BindEventData > >& myelbonds
    );

    /// bind row linker of myrank
    virtual void UpdateMyCrosslinkerBindingStates(
        const std::map< int, Teuchos::RCP< BEAMINTERACTION::DATA::BindEventData > >& mybonds,
        std::map< int, NewDoubleBonds>&       mynewdbondcl
    );

    /// bind row elements of myrank
    virtual void UpdateMyElementBindingStates(
        std::map< int, Teuchos::RCP< BEAMINTERACTION::DATA::BindEventData > > const& myelebonds
    );

    /// setup new double bonds
    virtual void CreateNewDoubleBondedCrosslinkerElementPairs( const std::map<int, NewDoubleBonds>& mynewdbondcl);

    /// unbind crosslinker if criteria are met
    virtual void UnBindCrosslinker();

    /// calclulate force dependent unbind probability for double bonded crosslinker
    /// according to Bell's equation (Howard, eq 5.10, p.89)
    void CalcBellsForceDependentUnbindProbability(
        CROSSLINKING::CrosslinkerNode* linker,
        Teuchos::RCP< BEAMINTERACTION::BeamLink > const& elepairptr,
        std::vector< double >& punlinkforcedependent
    ) const;

    /// communicate crosslinker unbinding event data
    virtual void CommunicateCrosslinkerUnbinding(
        std::map< int, std::vector< Teuchos::RCP< BEAMINTERACTION::DATA::UnBindEventData > > > & sendunbindevent,
        std::vector< Teuchos::RCP< BEAMINTERACTION::DATA::UnBindEventData > > &                 myrankunbindevent
    ) const;

    /// update binding status of beams after unbinding
    virtual void UpdateBeamBindingStatusAfterUnbinding(
        std::vector< Teuchos::RCP< BEAMINTERACTION::DATA::UnBindEventData > > const & unbindevent
    );

    /// -------------------------------------------------------------------------
    /// in case we have double bonded crosslinker on myrank we have to check if
    /// myrank is still owner of all its crosslinker (if not, set up double bond on
    /// other proc that is now responsible)
    /// -------------------------------------------------------------------------
    virtual void UpdateMyDoubleBondsAfterRedistribution();

    /// in case char vector containing double bonds is read by proc != proc
    /// that has written
    virtual void UpdateMyDoubleBondsRemoteIdList();

    /// dissolve certain bonds
    virtual void DissolveBond(
        DRT::Node* linker,
        int const freedbspotid,
        int const numbondsold,
        std::map< int, std::vector< Teuchos::RCP< BEAMINTERACTION::DATA::UnBindEventData > > > & sendunbindevents,
        std::vector< Teuchos::RCP< BEAMINTERACTION::DATA::UnBindEventData > > & myrankunbindevents);

    /// send double bonds to new owner if crosslinker ownership change
    /// in the course of redistribution
    virtual void CommunicateBeamLinkAfterRedistribution(
        std::map<int, std::vector<Teuchos::RCP<BEAMINTERACTION::BeamLink> > >& dbondcltosend
    );

    /// send data T to rank= mapkey
    template<typename T>
    void ISend(
        DRT::Exporter & exporter,
        std::vector<MPI_Request> & request,
        std::map< int, std::vector< Teuchos::RCP<T> > > const & send
    ) const;

    /// get number of request for each proc
    template<typename T>
    void PrepareReceivingProcs(
        std::map<int, std::vector< Teuchos::RCP<T> > > const & datasenttorank,
        std::vector<int> & summedtargets
    ) const;

    /// recieve "receivesize" number of T and store in vector recv
    template<typename T>
    void RecvAny(
        DRT::Exporter & exporter,
        int const receivesize,
        std::vector< Teuchos::RCP<T> > & recv
    ) const;

    /// unblocking send and blocking recvany
    template<typename T>
    void ISendRecvAny(
        std::map< int, std::vector< Teuchos::RCP<T> > > const & send,
        std::vector< Teuchos::RCP<T> > & recv
    ) const;

    // wait for all communication to finish
    virtual void Wait(
        DRT::Exporter & exporter,
        std::vector<MPI_Request> & request,
        int length
    ) const;

    /// debug feature to check bindevent structs
    virtual void PrintAndCheckBindEventData ( Teuchos::RCP< BEAMINTERACTION::DATA::BindEventData > bindeventdata ) const;

    //! @}

  private:

    //!@name crosslinking member variables
    //! @{

    //! data container holding all beam contact related parameters
    Teuchos::RCP<BEAMINTERACTION::CrosslinkingParams> crosslinking_params_ptr_;

    //! temporary storage for all relevant crosslinker data
    //! (vector key is col lid of crosslinker)
    std::vector< Teuchos::RCP< BEAMINTERACTION::DATA::CrosslinkerData > > crosslinker_data_;

    //! crosslinker exporter for crosslinker data container
    Teuchos::RCP<DRT::Exporter> cl_exporter_;

    //! beam exporter for beam data container
    Teuchos::RCP<DRT::Exporter> beam_exporter_;

    //! temporary storage for all relevant beam data during crosslinking
    //  (vector index is col lid of beamele)
    std::vector< Teuchos::RCP< BEAMINTERACTION::DATA::BeamData > > beam_data_;

    //! double bonded crosslinker that exert forces on network (map key is crosslinker gid)
    std::map<int, Teuchos::RCP<BEAMINTERACTION::BeamLink> > doublebondcl_;

    //! particle, i.e. crosslinker molecule discretization runtime vtp writer
    Teuchos::RCP<DiscretizationRuntimeVtpWriter> vtp_writer_ptr_;

    //! current linker displacement
    Teuchos::RCP<Epetra_Vector> linker_disnp_;

    //! summarized displacement of nodes since last redistribution
    Teuchos::RCP<Epetra_Vector> dis_at_last_redistr_;

    //! half interaction distance considering largest linker + tolerance
    double half_interaction_distance_;

    //! store node row map before current redistribution
    Teuchos::RCP<Epetra_Map> cl_noderowmap_prior_redistr_;

    //! store node row map before current redistribution
    Teuchos::RCP<Epetra_Map> cl_nodecolmap_prior_redistr_;

    //! store node row map before current redistribution
    Teuchos::RCP<Epetra_Map> beam_elerowmap_prior_redistr_;

    //! store node row map before current redistribution
    Teuchos::RCP<Epetra_Map> beam_elecolmap_prior_redistr_;

    //! @}

  };

}
}

#endif
