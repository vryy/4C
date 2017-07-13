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
class CrosslinkingParams;
class BeamToBeamLinkage;

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
    virtual void PreUpdateStepElement();

    //! derived
    virtual void UpdateStepElement();

    //! derived
    virtual void PostUpdateStepElement();

    //! derived
    virtual void OutputStepState(IO::DiscretizationWriter& iowriter) const;
    virtual void OutputStepState( bool const restart ) const;

    //! derived
    virtual void RuntimeOutputStepState() const;

    //! derived
    virtual void ResetStepState();

    //! derived
    virtual void WriteRestart(
        IO::DiscretizationWriter& iowriter,
        const bool& forced_writerestart) const;

    //! derived
    virtual void ReadRestart( IO::DiscretizationReader& ioreader );

    //! derived
    virtual void PostReadRestart();

    //! derived
    virtual void InitSubmodelDependencies(Teuchos::RCP<STR::MODELEVALUATOR::BeamInteraction::Map> const submodelmap);

    //! derived
    virtual void AddBinsToBinColMap( std::set< int >& colbins);

    //! derived
    virtual void AddBinsWithRelevantContentForIaDiscretColMap( std::set< int >& colbins) const;

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

    //! struct to store crosslinker data for faster access
    struct CrosslinkerData {
        LINALG::Matrix<3,1> clpos;                        // current position of crosslinker
        std::vector<std::pair<int, int> > clbspots;       // gid of element to local number of bspot, [0] and [1] first and second bspot
        int clnumbond;                                    // number of active bonds
        Teuchos::RCP<MAT::CrosslinkerMat> clmat;          // material attached to crosslinker node
        int clowner;                                      // owner of crosslinker
    };

    //! struct to store beam data for faster access
    struct BeamData {
        std::map<int, LINALG::Matrix<3,1> > bbspotpos;      // current position at bindingspots (xi) (key is local number of binding spot)
        std::map<int, LINALG::Matrix<3,3> > bbspottriad;    // current triad at bindingspots (xi) (key is local number of binding spot)
        std::map<int, int> bbspotstatus;                    // key is locn of bspot, holds gid of crosslinker if binded, if free = -1
        int bowner;                                         // owner of current beam
        INPAR::BEAMINTERACTION::FilamentType filamenttype;  // type of filament
//        int filnumber;
    };

    //! struct that stores all necessary data to handle the crosslining between two elements on each proc
    struct NewDoubleBonds {
        int id;                                           // gid of crosslinker
        std::vector< std::pair<int, int> > eleids;         // elegid and local binding spot number of first element
      };

    //! @}

    //! @name Small data structs that need to be communicated
    //! @{
    //! struct that stores binding event data that needs to be be communicated
    struct BindEventData {
        int clgid;            // gid of crosslinker
        int elegid;           // ele gid crosslinker wants to bind to
        int bspotlocn;        // loc number of bspot on ele cl wants to bind to
        int requestproc;      // myrank, processor that is requesting
        int permission;       // permission/veto, if crosslinker is allowed to bind
    };

    //! struct that stores unbinding event data that needs to be be communicated
    struct UnBindEventData {
      std::pair<int, int> eletoupdate;   // element gid (first) that needs to be updated at local binding (second)
    };

    //! @}

    //!@name routines that are not derived and handle crosslinking
    //! @{

    /// add crosslinker to bin discretization initially
    virtual void AddCrosslinkerToBinDiscretization();

    /// set binding spot positions on filament initially
    virtual void SetFilamentBindingSpotPositions();

    /// set filament types on elements
    virtual void SetFilamentFilamentTypes();

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

    /// diffuse crosslinker depending on number of bonds they have
    virtual void DiffuseCrosslinker();

    /// diffuse unbound crosslinker according to brownian dynamics
    virtual void DiffuseUnboundCrosslinker(CROSSLINKING::CrosslinkerNode* crosslinker) const;

    /// get binding spot of crosslinker that is currently occupied
    virtual int GetSingleOccupiedClBspot(
        std::vector<std::pair<int, int> > const& clbspots
    ) const;

    void SetPositionOfDoubleBondedCrosslinkerPBCconsistent(
        DRT::Node* crosslinker,
        LINALG::Matrix<3,1>& clpos,
        LINALG::Matrix<3,1> const& bspot1pos,
        LINALG::Matrix<3,1> const& bspot2pos
    ) const;

    /// set crosslinker position (change X)
    virtual void SetCrosslinkerPosition(
        DRT::Node* crosslinker,
        LINALG::Matrix<3,1> const& newclpos
    ) const;

    /// new position after transition from single to not bonded
    virtual void SetPositionOfNewlyFreeCrosslinker(
        CROSSLINKING::CrosslinkerNode* crosslinker,
        LINALG::Matrix<3,1>& clpos
    ) const;

    /// new position after transition from double to single bonded
    virtual void  SetPositionOfNewlySingleBondedCrosslinker(
        DRT::Node* crosslinker,
        CrosslinkerData& cldata,
        int const stayoccpotid);

    /// check if crosslink of adjacent elements on same filament
    virtual bool CheckCrosslinkOfAdjacentElements(
        DRT::Element const * const nbbeam,
        CrosslinkerData const& cldata_i
    ) const;

    /// check whether distance of two points is out of given range
    bool IsDistanceOutOfBindingRange(
        LINALG::Matrix<3,1> const& pos1,
        LINALG::Matrix<3,1> const& pos2,
        double const lowerbound,
        double const upperbound
    ) const;

    /// check whether enclosed angle of two beam centerline tangents is out of given range
    bool IsEnclosedAngleOfBSpotTangentsOutOfRange(
        LINALG::Matrix<3,1> const& direction1,
        LINALG::Matrix<3,1> const& direction2,
        double const lowerbound,
        double const upperbound
    ) const;

    /// fill epetar vectors to write vtp output
    void FillStateDataVectorsForOutput(
        Teuchos::RCP<Epetra_Vector> displacement,
        Teuchos::RCP<Epetra_Vector> orientation,
        Teuchos::RCP<Epetra_Vector> numberofbonds,
        Teuchos::RCP<Epetra_Vector> owner
    ) const;

    // -------------------------------------------------------------------------
    // loop over all column crosslinker and pre compute their data that is needed
    // multiple times in the following, therefore this gives faster access
    // note: we get references here, i.e. changes in the crosslinker status will
    // be done with the variables created here (only row owner change something,
    // this is ensured in the actual algorithm)
    // -------------------------------------------------------------------------
    /// get crosslink data before interaction evaluation
    virtual void PreComputeCrosslinkerData();

    /// get beam data before interaction evaluation
    virtual void PreComputeBeamData();

    /// get crosslink and beam data before interaction evaluation
    virtual void PreComputeCrosslinkerAndBeamData();

    /// get crosslink data before evaluation
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
        std::map<int, BindEventData >&              mybonds,
        std::map<int, std::vector<BindEventData> >& undecidedbonds
    );

    /// find potential binding events in one bin
    virtual void FindPotentialBindingEventsInBinAndNeighborhood(
        DRT::Element* bin,
        std::map<int, BindEventData >& mybonds,
        std::map<int, std::vector<BindEventData> >& undecidedbonds,
        std::vector< std::set <int> >& intendedbeambonds,
        bool const checklinkingprop
    );

    /// search for binding events on each proc separately (i.e. pretending myran is alone)
    /// communication to ensure correct binding over all procs is done afterwards
    virtual void PrepareBinding(
        DRT::Node*                                  crosslinker_i,
        std::set<DRT::Element*> const&              neighboring_beams,
        std::map<int, BindEventData >&              mybonds,
        std::map<int, std::vector<BindEventData> >& undecidedbonds,
        std::vector< std::set <int> >&              intendedbeambonds,
        bool const checklinkingprop
    );

    /// check criteria if binding event is feasible
    virtual bool CheckBindEventCriteria(
        CROSSLINKING::CrosslinkerNode const * const crosslinker_i,
        DRT::Element const * const potbeampartner,
        CrosslinkerData const& cldata_i,
        BeamData const& beamdata_i,
        int locnbspot,
        std::vector< std::set <int> >& intendedbeambonds,
        bool checklinkingprop
    ) const;

    /// check if crosslinke and filament type are compatible
    virtual bool CheckLinkerAndFilamentTypeCompatibility(
        CROSSLINKING::CrosslinkerNode const * const crosslinker_i,
        BeamData const& beamdata) const;

    /// if crosslinker is singly bound, we fetch the orientation vector of the
    /// filament axis at the already occupied binding spot for the orientation
    /// criterion (enclosed angle) to be checked later on
    virtual void GetOccupiedClBSpotBeamTangent(
        CrosslinkerData const& cldata_i,
        LINALG::Matrix<3,1>& occ_bindingspot_beam_tangent,
        int const clgid
    ) const;

    /// decide by asking other procs who is allowed to set specific crosslinker,
    /// this is necessary to avoid setting crosslinker more than once per time step
    virtual void ManageBindingInParallel(
        std::map<int, BindEventData>&               mysetbonds,        // clgid to cldata
        std::map<int, std::vector<BindEventData> >& undecidedbonds,    // owner of cldatas in vector to be requested
        std::map<int, BindEventData >&              myelebonds
    ) const;

    /// communicate requests
    virtual void CommunicateUndecidedBonds(
        DRT::Exporter & exporter,
        std::map<int, std::vector<BindEventData> >& undecidedbonds,
        int& numrecrequest,
        std::map<int, std::vector<BindEventData> >& requestedcl
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
        std::map<int, std::vector<BindEventData> >& requestedcl,
        std::map<int, BindEventData >&              mybonds,
        std::map<int, std::vector<BindEventData> >& decidedbonds
    ) const;

    /// communicate decisions for binding events
    virtual void CommunicateDecidedBonds(
        DRT::Exporter& exporter,
        std::map<int, std::vector<BindEventData> >& decidedbonds,
        std::map<int, BindEventData >&              myelebonds,
        int const numrecrequest,
        int const answersize
    ) const;

    /* now have two distinct maps of binding events on each proc, depending
       on ownership of crosslinker and elements myrank has different tasks:
        - mybonds: myrank takes care of crosslinker and (most) elements
        - myelebonds: myrank takes care of elements
       within those maps, different treatment is necessary for free and single
       bonded linker
                                                                            */
    virtual void UpdateMyCrosslinkerAndElementBindingStates(
        std::map<int, BindEventData >& mybonds,
        std::map<int, BindEventData >& myelbonds
    );

    /// bind row linker of myrank
    virtual void UpdateMyCrosslinkerBindingStates(
        const std::map<int, BindEventData >& mybonds,
        std::map<int, NewDoubleBonds>&       mynewdbondcl
    );

    /// bind row elements of myrank
    virtual void UpdateMyElementBindingStates(
        std::map<int, BindEventData > const& myelebonds
    );

    /// setup new double bonds
    virtual void CreateNewDoubleBondedCrosslinkerElementPairs( const std::map<int, NewDoubleBonds>& mynewdbondcl);

    /// unbind crosslinker if criteria are met
    virtual void UnBindCrosslinker();

    /// calclulate force dependent unbind probability for double bonded crosslinker
    /// according to Bell's equation (Howard, eq 5.10, p.89)
    void CalcBellsForceDependentUnbindProbability(
        CROSSLINKING::CrosslinkerNode* linker,
        Teuchos::RCP< BEAMINTERACTION::BeamToBeamLinkage > const& elepairptr,
        std::vector< double >& punlinkforcedependent
    ) const;

    /// communicate crosslinker unbinding event data
    virtual void CommunicateCrosslinkerUnbinding(
        std::map<int, std::vector<UnBindEventData> >& sendunbindevent,
        std::vector<UnBindEventData>&                 myrankunbindevent
    ) const;

    /// update binding status of beams after unbinding
    virtual void UpdateBeamBindingStatusAfterUnbinding(
        const std::vector<UnBindEventData>& unbindevent
    );

    /// -------------------------------------------------------------------------
    /// in case we have double bonded crosslinker on myrank we have to check if
    /// myrank is still owner of all its crosslinker (if not, set up double bond on
    /// other proc that is now responsible)
    /// -------------------------------------------------------------------------
    virtual void UpdateMyDoubleBondsAfterRedistribution();

    /// in case char vector containing double bonds is read by proc != proc
    /// that has written
    virtual void UpdateMyDoubleBondsAfterRestartRemoteIdList();

    /// dissolve certain bonds
    virtual void DissolveBond(
        DRT::Node* linker,
        int const freedbspotid,
        int const numbondsold,
        std::map< int, std::vector< UnBindEventData > >& sendunbindevents,
        std::vector< UnBindEventData >& myrankunbindevents);

    /// send double bonds to new owner if crosslinker ownership change
    /// in the course of redistribution
    virtual void CommunicateBeamToBeamLinkageAfterRedistribution(
        std::map<int, std::vector<Teuchos::RCP<BEAMINTERACTION::BeamToBeamLinkage> > >& dbondcltosend
    );

    /// send data T to rank= mapkey
    template<typename T>
    void ISend(
        DRT::Exporter& exporter,
        std::vector<MPI_Request>& request,
        std::map<int, std::vector<T> > const& send
    ) const;

    /// get number of request for each proc
    template<typename T>
    void PrepareReceivingProcs(
        std::map<int, std::vector<T> > const& datasenttorank,
        std::vector<int>& summedtargets
    ) const;

    /// recieve "receivesize" number of T and store in vector recv
    template<typename T>
    void RecvAny(
        DRT::Exporter&  exporter,
        int const receivesize,
        std::vector<T>& recv
    ) const;

    /// unblocking send and blocking recvany
    template<typename T>
    void ISendRecvAny(
        std::map<int, std::vector<T> > const& send,
        std::vector<T>&               recv
    ) const;

    // wait for all communication to finish
    virtual void Wait(
        DRT::Exporter& exporter,
        std::vector<MPI_Request>& request,
        int const length
    ) const;

    /// pack binding event data
    virtual void Pack(
        DRT::PackBuffer&     data,
        BindEventData const& bindeventdata
    ) const;

    /// pack unbinding event data
    virtual void Pack(
        DRT::PackBuffer&       data,
        UnBindEventData const& unbindeventdata
    ) const;

    /// unpack binding event data
    virtual void UnPack(
        std::vector<char>::size_type& position,
        std::vector<char> data,
        BindEventData&    bindeventdata
    ) const;

    /// unpack unbinding event data
    virtual void UnPack(
        std::vector<char>::size_type& position,
        std::vector<char> data,
        UnBindEventData&  unbindeventdata
    ) const;

    /// debug feature to check bindevent structs
    virtual void PrintAndCheckBindEventData ( BindEventData const& bindeventdata ) const;

    //! @}

  private:

    //!@name crosslinking member variables
    //! @{

    //! data container holding all beam contact related parameters
    Teuchos::RCP<BEAMINTERACTION::CrosslinkingParams> crosslinking_params_ptr_;

    //! type of eles in bins
    BINSTRATEGY::UTILS::BinContentType bin_beamcontent_;

    //! temporary storage for all relevant crosslinker data
    //! (vector key is col lid of crosslinker)
    std::vector<CrosslinkerData> crosslinker_data_;

    //! temporary storage for all relevant beam data during crosslinking
    //  (vector index is col lid of beamele)
    std::vector<BeamData> beam_data_;

    //! double bonded crosslinker that exert forces on network (map key is crosslinker gid)
    std::map<int, Teuchos::RCP<BEAMINTERACTION::BeamToBeamLinkage> > doublebondcl_;

    //! particle, i.e. crosslinker molecule discretization runtime vtp writer
    Teuchos::RCP<DiscretizationRuntimeVtpWriter> vtp_writer_ptr_;
    //! @}

  };

}
}

#endif
