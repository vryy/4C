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
#include "str_elements_paramsinterface.H"       // interface to the element evaluation
#include "../drt_inpar/inpar_binstrategy.H"
#include "../linalg/linalg_fixedsizematrix.H"


//forward declaration ...
namespace DRT
{
  class Element;
  class Node;
  class PackBuffer;
}
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
      int clnumbond;                                  // number of active bonds
      Teuchos::RCP<MAT::CrosslinkerMat> clmat;        // material attached to crosslinker node
      int clowner;                                    // owner of crosslinker
    };

    //! struct to store beam data for faster access
    struct BeamData {
      std::map<int, std::vector<LINALG::Matrix<3,1> > > bbspotdofs;    // currrent values of dofs at bindingspots (xi) (key is local number of binding spot)
      std::map<int, int> bbspotstatus;                                 // key is locn of bspot, holds gid of crosslinker if binded, if free = -1
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
    struct DoubleBondedCl {
      int id;                          // gid of crosslinker
      std::pair<int, int> eleone;      // elegid and local binding spot number of first element
      std::pair<int, int> eletwo;      // elegid and local binding spot number of second element
    };

    //! struct that stores unbinding event data that needs to be be communicated
    struct UnBindEventData {
      int id;                            // this id represents the clgid that needs to be deleted from the double bond list
      std::pair<int, int> eletoupdate;   // element gid (first) that needs to be updated at
    };

    //! @}


  protected:
    //! derived
    virtual void Reset(const Epetra_Vector& x);

  private:
    //! apply Browian (stochastic and damping forces)
    bool ApplyForceBrownian();

    //! apply statmech specific neumann conditions
    bool ApplyForceExternal();

    //! apply Browian (stochastic and damping forces)
    bool ApplyForceStiffBrownian();

    //! apply statmech specific neumann conditions
    bool ApplyForceStiffExternal();

    //! evaluate statmech specific neumann conditions
    void EvaluateNeumannStatMech(
        Teuchos::RCP<Epetra_Vector> eval_vec,
        Teuchos::RCP<LINALG::SparseOperator> eval_mat);

    //! evaluate Browian (stochastic and damping forces)
    void EvaluateBrownian(
        Teuchos::RCP<LINALG::SparseOperator>* eval_mat,
        Teuchos::RCP<Epetra_Vector>* eval_vec);

    //! evaluate Browian (stochastic and damping forces)
    void EvaluateBrownian(
        Teuchos::ParameterList& p,
        Teuchos::RCP<LINALG::SparseOperator>* eval_mat,
        Teuchos::RCP<Epetra_Vector>* eval_vec);


    //!@name routines that manage crosslinker as particles as well as
    //       the whole corresponding binstrategy
    //! @{

    /// set random initial position of crosslinker
    virtual void InitializeBinDiscret();

    /// update states based on bindis after its redistribution
    virtual void UpdateMaps();

    /// build ele to bin map
    virtual void BuildEleToBinMap();

    /// update states based on bindis after its redistribution
    virtual void DiffuseCrosslinker();

    /// update bin, crosslinker and beam distribution
    virtual void UpdateBinStrategy(
        bool transfer    = true,
        bool partition   = false,
        bool repartition = false,
        bool createxaabb = false,
        bool setcutoff   = false
    );

    /// search and set crosslinker
    virtual void UpdateCrosslinking();

    /// search and set crosslinker
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

    /// decide by asking other procs who is allowed to set specific crosslinker,
    /// this is necessary to avoid setting crosslinker more than once per time step
    virtual void ManageBindingInParallel(
      std::map<int, BindEventData>&               mysetbonds,        // clgid to cldata
      std::map<int, std::vector<BindEventData> >& undecidedbonds,    // owner of cldatas in vector to be requested
      std::map<int, BindEventData >&              myelebonds
    );

    /// search for binding events on each proc separately (i.e. pretending myran is alone)
    /// communication to ensure correct binding over all procs is done afterwards
    virtual void BindCrosslinker(
        std::map<int, BindEventData >& mybonds,
        std::map<int, BindEventData >& myelbonds
    );

    /// search and set crosslinker
    virtual void UnBindCrosslinker();

    /// get crosslink data before evaluation
    virtual void PreComputeCrosslinkerData(
      const int numcolcl
    );

    /// get crosslink data before evaluation
    virtual void PreComputeBeamData(
      const int numcolbeams
    );

    virtual void CommunicateDoubleBondedCl(
      std::map<int, DoubleBondedCl>&               mydbondcl,
      std::map<int, std::vector<DoubleBondedCl> >& dbondcltosend
    );

    /// communicate new crosslinker elements
    virtual void CommunicateCrosslinkerUnbinding(
      std::map<int, std::vector<UnBindEventData> >& sendunbindevent,
      std::vector<UnBindEventData>&                 myrankunbindevent
    );

    /// update map of vector (based on dof map) after redistribution
    virtual void UpdateDofMapOfVector(
      Teuchos::RCP<DRT::Discretization> discret,
      Teuchos::RCP<Epetra_Vector>&      dofmapvec,
      Teuchos::RCP<Epetra_Vector> old = Teuchos::null
    );

    /// get neighbouring eles in discret
    virtual void GetNeighboringEles();

//    ///
//    template<typename T>
//    void ISendRecv(
//      std::map<int, std::vector<T> >& send,
//      std::map<int, T>&               recv
//    );

    /// pack binding event data
    virtual void Pack(
      DRT::PackBuffer&     data,
      const BindEventData& bindeventdata
    );

    /// pack unbinding event data
    virtual void Pack(
      DRT::PackBuffer&      data,
      const DoubleBondedCl& dbondedcldata
    );

    /// pack unbinding event data
    virtual void Pack(
      DRT::PackBuffer&       data,
      const UnBindEventData& unbindeventdata
    );

    /// unpack binding event data
    virtual void UnPack(
      std::vector<char>::size_type& position,
      std::vector<char> data,
      BindEventData&    bindeventdata
    );

    /// unpack unbinding event data
    virtual void UnPack(
      std::vector<char>::size_type& position,
      std::vector<char> data,
      DoubleBondedCl&   dbondedcldata
    );

    /// unpack unbinding event data
    virtual void UnPack(
      std::vector<char>::size_type& position,
      std::vector<char> data,
      UnBindEventData&  unbindeventdata
    );

    /// print welcome to biopolymer network simulation
    virtual void Logo();

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
    Teuchos::RCP<DRT::Discretization> intactdis_;

    //! global displacements \f${D}_{n+1}\f$ at \f$t_{n+1}\f$ based on maps of interaction discretization
    Teuchos::RCP<Epetra_Vector> ia_disnp_;

    //! crosslinker management
    Teuchos::RCP<PARTICLE::Algorithm> particlealgo_;

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
    std::map<int, DoubleBondedCl> doublebondcl_;

  }; // class Crosslink

} // namespace ModelEvaluator
} // namespace STR



#endif /* STR_MODEL_EVALUATOR_CROSSLINKING_H_ */
