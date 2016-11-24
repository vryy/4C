/*-----------------------------------------------------------*/
/*!
\file str_model_evaluator_crosslinking.H

\brief model evaluator for crosslinking in biopolymer networks

\maintainer Jonas Eichinger

\date May, 2016

\level 3

*/
/*-----------------------------------------------------------*/


#ifndef STR_MODEL_EVALUATOR_CROSSLINKING_H_
#define STR_MODEL_EVALUATOR_CROSSLINKING_H_

#include "str_model_evaluator_generic.H"        // base class
#include <Teuchos_RCP.hpp>
#include <Epetra_SerialDenseVector.h>
#include "str_elements_paramsinterface.H"       // interface to the element evaluation
#include "../drt_inpar/inpar_binstrategy.H"
#include "../linalg/linalg_fixedsizematrix.H"
#include <Epetra_Comm.h>
#include <Epetra_MpiComm.h>


//forward declaration ...
namespace ADAPTER
{
  class Coupling;
} // namespace ADAPTER
namespace FSI{
namespace UTILS{
  class MatrixRowTransform;
} // namespace FSI
} // namespace UTILS
namespace DRT
{
  class Exporter;
  class Element;
  class Node;
  class PackBuffer;
} // namespace DRT
namespace CROSSLINKING
{
 class CrosslinkerNode;
} // namespace CROSSLINKING
namespace MAT
{
  class CrosslinkerMat;
} //namespace MATERIAL
namespace LINALG {
  class SparseMatrix;
} // namespace LINALG
namespace PARTICLE
{
  class Algorithm;
} // namespace PARTICLE
namespace BEAMINTERACTION {
class BeamToBeamLinkage;
}

namespace STR{
namespace MODELEVALUATOR{

  // forward declaration
  class StatMechData;

  class Crosslinking : public Generic
  {
  public:
    //! constructor
    Crosslinking();

    //! destructor
    virtual ~Crosslinking(){};

    virtual void Setup();

    //! @name Derived public STR::MODELEVALUATOR::Generic methods
    //! @{
    //! derived

    //! derived
    INPAR::STR::ModelType Type() const
    { return INPAR::STR::model_crosslinking; }

    //! derived
    virtual bool EvaluateForce();

    //! derived
    virtual bool EvaluateStiff();

    //! derived
    virtual bool EvaluateForceStiff();

    //! derived
    virtual bool AssembleForce(Epetra_Vector& f,
        const double & timefac_np) const;

    //! derived
    virtual bool AssembleJacobian(LINALG::SparseOperator& jac,
        const double & timefac_np) const;

    //! derived
    virtual void WriteRestart(
        IO::DiscretizationWriter& iowriter,
        const bool& forced_writerestart) const;

    //! derived
    virtual void ReadRestart(
        IO::DiscretizationReader& ioreader);

    //! derived
    virtual void RecoverState(
            const Epetra_Vector& xold,
            const Epetra_Vector& dir,
            const Epetra_Vector& xnew);

    //! derived
    virtual void UpdateStepState(const double& timefac_n);

    //! derived
    virtual void UpdateStepElement();

    //! derived
    void DetermineStressStrain();

    //! derived
    void DetermineEnergy();

    //! derived
    virtual void OutputStepState(IO::DiscretizationWriter& iowriter) const;

    //! derived
    virtual Teuchos::RCP<const Epetra_Map> GetBlockDofRowMapPtr() const;

    //! derived
    virtual Teuchos::RCP<const Epetra_Vector> GetCurrentSolutionPtr() const;

    //! derived
    virtual Teuchos::RCP<const Epetra_Vector> GetLastTimeStepSolutionPtr() const;

    //! derived
    virtual void PostOutput();

    //! derived
    virtual void ResetStepState();
    //! @}

    //! @name Small data structs for faster access, better structure and communication
    //! @{

    //! struct to store crosslinker data for faster access
    struct CrosslinkerData {
        LINALG::Matrix<3,1> clpos;                      // current position of crosslinker
        std::vector<std::pair<int, int> > clbspots;     // gid of element to local number of bspot, [0] and [1] first and second bspot
        std::vector<std::pair<int, int> > bnodegids_;   // gid of first and second node of binding partner element, [0] and [1] first and second bspot
        int clnumbond;                                  // number of active bonds
        Teuchos::RCP<MAT::CrosslinkerMat> clmat;        // material attached to crosslinker node
        int clowner;                                    // owner of crosslinker
    };

    //! struct to store beam data for faster access
    struct BeamData {
        std::map<int, LINALG::Matrix<3,1> > bbspotpos;    // current position at bindingspots (xi) (key is local number of binding spot)
        std::map<int, LINALG::Matrix<3,3> > bbspottriad;  // current triad at bindingspots (xi) (key is local number of binding spot)
        std::map<int, int> bbspotstatus;                  // key is locn of bspot, holds gid of crosslinker if binded, if free = -1
//      int type;
//      int filnumber;
        int bowner;                                                      // owner of current beam
    };

    //! struct that stores binding event data that needs to be be communicated
    struct BindEventData {
        int clgid;            // gid of crosslinker
        int elegid;           // ele gid crosslinker wants to bind to
        int bspotlocn;        // loc number of bspot on ele cl wants to bind to
        int requestproc;      // myrank, processor that is requesting
        int permission;       // permission/veto, if crosslinker is allowed to bind
    };

    //! struct that stores all necessary data to handle the crosslining between two elements on each proc
    struct NewDoubleBonds {
        int id;                                    // gid of crosslinker
        std::vector<std::pair<int, int> > eleids;  // elegid and local binding spot number of first element
      };

    //! struct that stores unbinding event data that needs to be be communicated
    struct UnBindEventData {
      std::pair<int, int> eletoupdate;   // element gid (first) that needs to be updated at
    };

    //! struct that stores force and stiff contributions that need to be communicated between procs
    struct CommForceStiff {
      int elegid1;
      int elegid2;
      Epetra_SerialDenseVector ele1force;
      Epetra_SerialDenseVector ele2force;
      Epetra_SerialDenseMatrix ele11stiff;
      Epetra_SerialDenseMatrix ele12stiff;
      Epetra_SerialDenseMatrix ele21stiff;
      Epetra_SerialDenseMatrix ele22stiff;
    };

    //! @}


  protected:
    //! derived
    virtual void Reset(const Epetra_Vector& x);

  private:
    //!@name routines that manage crosslinker as particles as well as
    //       the whole corresponding binstrategy
    //! @{

    virtual void ResetStateOfElementPairs();

    /// set random initial position of crosslinker
    virtual void InitializeBinDiscret();

    /// update states based on bindis after its redistribution
    virtual void UpdateMaps();

    /// build ele to bin map
    virtual void BuildEleToBinMap();

    /// build ele to bin map
    virtual void UpdateMyDoubleBondsAfterRedistribution();

    /// diffuse crosslinker depending on number of bonds they have
    virtual void DiffuseCrosslinker();

    /// diffuse unbound crosslinker according to brownian dynamics
    virtual void DiffuseUnboundCrosslinker(CROSSLINKING::CrosslinkerNode* crosslinker) const;

    /// set crosslinker position (change X)
    virtual void SetCrosslinkerPosition(
        DRT::Node* crosslinker,
        const LINALG::Matrix<3,1>& newclpos
    ) const;

    virtual void SetPositionOfNewlyFreeCrosslinker(
        DRT::Node* crosslinker,
        LINALG::Matrix<3,1>& clpos
    ) const;

    /// update bin, crosslinker and beam distribution
    virtual void UpdateBinStrategy(
        bool transfer          = true,
        bool partition         = false,
        bool repartition       = false,
        bool createxaabb       = false,
        bool setcutoff         = false,
        bool cltoclinteraction = false
    );

    /// search and set crosslinker
    virtual void UpdateCrosslinking();

    /// search and set crosslinker
    /* -------------------------------------------------------------------------
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
    virtual void LookForBindingEvents(
        std::map<int, BindEventData >&              mybonds,
        std::map<int, std::vector<BindEventData> >& undecidedbonds
    );

    /// search for binding events on each proc separately (i.e. pretending myran is alone)
    /// communication to ensure correct binding over all procs is done afterwards
    virtual void PrepareBinding(
        DRT::Node*                                  crosslinker_i,
        const std::set<DRT::Element*>&              neighboring_beams,
        std::map<int, BindEventData >&              mybonds,
        std::map<int, std::vector<BindEventData> >& undecidedbonds
    );

    /// check if crosslink of adjacent elements on same filament
    virtual bool CheckCrosslinkOfAdjacentElements(
        DRT::Element* nbbeam,
        const std::vector<std::pair<int, int> >& bnodegids
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

    /// communicate decisions for binding events
    virtual void CommunicateDecidedBonds(
        DRT::Exporter& exporter,
        std::map<int, std::vector<BindEventData> >& decidedbonds,
        std::map<int, BindEventData >&              myelebonds,
        const int& numrecrequest,
        const int& answersize
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

    /* now have two distinct maps of binding events on each proc, depending
       on ownership of crosslinker and elements myrank has different tasks:
        - mybonds: myrank takes care of crosslinker and (most) elements
        - myelebonds: myrank takes care of elements
       within those maps, different treatment is necessary for free and single
       bonded linker
                                                                            */
    virtual void BindMyCrosslinkerAndElements(
        std::map<int, BindEventData >& mybonds,
        std::map<int, BindEventData >& myelbonds
    );

    /// bind row linker of myrank
    virtual void BindMyCrosslinker(
        const std::map<int, BindEventData >& mybonds,
        std::map<int, NewDoubleBonds>&       mynewdbondcl
    );

    /// bind row elements of myrank
    virtual void BindMyElements(
        const std::map<int, BindEventData >& myelebonds
    );

    /// setup new double bonds
    virtual void SetupNewDoubleBonds(
        const std::map<int, NewDoubleBonds>& mynewdbondcl,
        bool precomputed = false
    );

    /// search and set crosslinker
    virtual void UnBindCrosslinker();

    /// update binding status of beams after unbinding
    virtual void GetProcsInvolvedInMyDoubleBonds();

    /// update binding status of beams after unbinding
    virtual void UpdateBeamBindingStatusAfterUnbinding(
        const std::vector<UnBindEventData>& unbindevent
    );

    // update list of double bonds on myrank after unbinding
    virtual void UpdateDoubleBondsAfterUbinding(
        const std::vector<int>& doublebondstodelete
    );

    // -------------------------------------------------------------------------
    // loop over all column crosslinker and pre compute their data that is needed
    // multiple times in the following, therefore this gives faster access
    // note: we get references here, i.e. changes in the crosslinker status will
    // be done with the variables created here (only row owner change something,
    // this is ensured in the actual algorithm)
    // -------------------------------------------------------------------------
    /// get crosslink data before evaluation
    virtual void PreComputeCrosslinkerData(
        const int numcolcl
    );

    /// get crosslink data before evaluation
    virtual void PreComputeBeamData(
        const int numcolbeams
    );

    /// get current element displacements
    virtual void GetCurrentElementDis(
        const DRT::Element* ele,
        const Teuchos::RCP<Epetra_Vector> ia_discolnp,
        std::vector<double>& eledisp
    ) const;

    /// get position and triad at certain binding spot of element
    virtual void GetPosAndTriadOfBindingSpot(
        DRT::Element* ele,
        const Teuchos::RCP<Epetra_Vector> ia_discolnp,
        const int& locbspotnum,
        LINALG::Matrix<3,1>& bspotpos,
        LINALG::Matrix<3,3>& bspottriad
    ) const;

    /// get position and triad at certain binding spot of element with pre computed element displacements
    virtual void GetPosAndTriadOfBindingSpot(
        DRT::Element* ele,
        const Teuchos::RCP<Epetra_Vector> ia_discolnp,
        const int& locbspotnum,
        LINALG::Matrix<3,1>& bspotpos,
        LINALG::Matrix<3,3>& bspottriad,
        std::vector<double>& eledisp
    ) const;

    /// get binding spot of crosslinker that is currently occupied
    virtual void GetOccupiedClBspot(
        int& occbspotid,
        const std::vector<std::pair<int, int> >& clbspots
    ) const;

    /// send double bonds to new owner if crosslinker ownership change
    /// in the course of redistribution
    virtual void CommunicateBeamToBeamLinkageAfterRedistribution(
        std::map<int, std::vector<Teuchos::RCP<BEAMINTERACTION::BeamToBeamLinkage> > >& dbondcltosend
    );

    /// communicate force and stiff contribution data
    virtual void CommunicateForceStiff(
        std::map<int, std::vector<CommForceStiff> >& sendforcestiff,
        std::vector<CommForceStiff>&                 recvforcestiff
    ) const;

    /// communicate crosslinker unbinding event data
    virtual void CommunicateCrosslinkerUnbinding(
        std::map<int, std::vector<UnBindEventData> >& sendunbindevent,
        std::vector<UnBindEventData>&                 myrankunbindevent
    ) const;

    /// update map of vector (based on dof map) after redistribution
    virtual void UpdateDofMapOfVector(
        Teuchos::RCP<DRT::Discretization> discret,
        Teuchos::RCP<Epetra_Vector>&      dofmapvec,
        Teuchos::RCP<Epetra_Vector> old = Teuchos::null
    );

    /// assemble element force vector and stiffness matrix into system vector and matrix
    virtual void AssembleEleForceStiffIntoSystemVectorMatrix(
        const DRT::Discretization&       discret,
        const int                        elegid1,
        const int                        elegid2,
        const Epetra_SerialDenseVector&  elevec1,
        const Epetra_SerialDenseVector&  elevec2,
        const Epetra_SerialDenseMatrix&  elemat11,
        const Epetra_SerialDenseMatrix&  elemat12,
        const Epetra_SerialDenseMatrix&  elemat21,
        const Epetra_SerialDenseMatrix&  elemat22,
        Teuchos::RCP<Epetra_Vector>                   sysvec,
        Teuchos::RCP<LINALG::SparseMatrix>            sysmat
      ) const;

    /// assemble element force vector and stiffness matrix received from other procs into system vector and matrix
    virtual void AssembleRecvEleForceStiffIntoSystemVectorMatrix(
      const std::vector<CommForceStiff>& recvforcestiff,
      Teuchos::RCP<Epetra_Vector>        sysvec,
      Teuchos::RCP<LINALG::SparseMatrix> sysmat
    ) const;

    /// transform force vector and stiffness matrix from Discret to ia_discret_
    virtual void TransformForceAndStiff(
        bool force = true,
        bool stiff = true
      );

    /// get neighbouring eles in discret
    virtual void GetNeighboringEles();

    /// exclude certein neighbors from interaction evaluation
    virtual void VerifyNeighbors(
        DRT::Element* currele,
        std::set<DRT::Element*>& neighbors
    ) const;

    /// send data T to rank= mapkey
    template<typename T>
    void ISend(
        DRT::Exporter& exporter,
        std::vector<MPI_Request>& request,
        const std::map<int, std::vector<T> >& send
    ) const;

    /// get number of request for each proc
    template<typename T>
    void PrepareReceivingProcs(
        const std::map<int, std::vector<T> >& datasenttorank,
        std::vector<int>& summedtargets
    ) const;

    /// recieve "receivesize" number of T and store in vector recv
    template<typename T>
    void RecvAny(
        DRT::Exporter&  exporter,
        const int& receivesize,
        std::vector<T>& recv
    ) const;

    /// unblocking send and blocking recvany
    template<typename T>
    void ISendRecvAny(
        const std::map<int, std::vector<T> >& send,
        std::vector<T>&               recv
    ) const;

    // wait for all communication to finish
    virtual void Wait(
        DRT::Exporter& exporter,
        std::vector<MPI_Request>& request,
        const int& length
    ) const;

    /// pack binding event data
    virtual void Pack(
        DRT::PackBuffer&     data,
        const BindEventData& bindeventdata
    ) const;

    /// pack unbinding event data
    virtual void Pack(
        DRT::PackBuffer&       data,
        const UnBindEventData& unbindeventdata
    ) const;

    /// pack force and stiff contribution data
    virtual void Pack(
        DRT::PackBuffer&       data,
        const CommForceStiff&  forcestiffdata
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

    /// unpack force and stiff contribution data
    virtual void UnPack(
        std::vector<char>::size_type& position,
        std::vector<char> data,
        CommForceStiff&   forcestiffdata
    ) const;

    /// print welcome to biopolymer network simulation
    virtual void Logo() const;

    //! @}

  private:

    //! statmech evaluation data container
    Teuchos::RCP<STR::MODELEVALUATOR::StatMechData> eval_statmech_ptr_;

    //! myrank
    int myrank_;

    //! number of procs
    int numproc_;

    //! interaction discretization handling all interactions (e.g. crosslinker to beam,
    //! beam to beam, potential ...)
    Teuchos::RCP<DRT::Discretization> ia_discret_;

    //! coupling adapter to transfer vectors and matrices between Discret() and intactids_
    Teuchos::RCP<ADAPTER::Coupling> coupsia_;

    //! transform object for structure stiffness matrix
    Teuchos::RCP<FSI::UTILS::MatrixRowTransform> siatransform_;

    //! global force based on Discret() at \f$t_{n+1}\f$
    Teuchos::RCP<Epetra_Vector> force_crosslink_;

    //! structural stiffness matrix based on Discret()
    Teuchos::RCP<LINALG::SparseMatrix> stiff_crosslink_;

    //! global force based on intactdis_ at \f$t_{n+1}\f$
    Teuchos::RCP<Epetra_Vector> ia_force_crosslink_;

    //! structural stiffness matrix based on intactdis_
    Teuchos::RCP<LINALG::SparseMatrix> ia_stiff_crosslink_;

    //! global displacements \f${D}_{n+1}\f$ at \f$t_{n+1}\f$ based on maps of interaction discretization
    Teuchos::RCP<Epetra_Vector> ia_disnp_;

    //! crosslinker management
    Teuchos::RCP<PARTICLE::Algorithm> binning_;

    //! crosslinker and bin discretization
    Teuchos::RCP<DRT::Discretization> bindis_;

    //! elerowmap of bindis
    Teuchos::RCP<Epetra_Map> rowbins_;

    //! extended bin to ele map
    std::map<int, std::set<int> > extbintoelemap_;

    //! extended ele to bin map
    std::map<int, std::set<int> > exteletobinmap_;

    //! type of eles in bins
    INPAR::BINSTRATEGY::BinContent bin_beamcontent_;

    //! temporary storage for all relevant crosslinker data
    //! (vector key is col lid of crosslinker)
    std::vector<CrosslinkerData> crosslinker_data_;

    //! temporary storage for all relevant beam data during crosslinking
    //  (vector key is col lid of beamele)
    std::vector<BeamData> beam_data_;

    //! double bonded crosslinker that exert forces on network
    std::map<int, Teuchos::RCP<BEAMINTERACTION::BeamToBeamLinkage> > doublebondcl_;



  }; // class Crosslink

} // namespace ModelEvaluator
} // namespace STR



#endif /* STR_MODEL_EVALUATOR_CROSSLINKING_H_ */
