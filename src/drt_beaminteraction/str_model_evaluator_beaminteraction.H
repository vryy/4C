/*-----------------------------------------------------------*/
/*!
\file str_model_evaluator_beaminteraction.H

\brief Evaluation of all beam interaction terms

\maintainer Jonas Eichinger, Maximilian Grill

\level 3
*/
/*-----------------------------------------------------------*/


#ifndef STR_MODEL_EVALUATOR_BEAMINTERACTION_H
#define STR_MODEL_EVALUATOR_BEAMINTERACTION_H

#include "../drt_structure_new/str_model_evaluator_generic.H"        // base class
#include "../drt_inpar/inpar_beaminteraction.H"
#include "../linalg/linalg_mapextractor.H"


//forward declaration ...
namespace ADAPTER
{
  class Coupling;
}
namespace DRT {
  class Discretization;
}
namespace FSI{
namespace UTILS{
  class MatrixRowTransform;
}
}
namespace LINALG {
  class SparseMatrix;
}
namespace BINSTRATEGY{
  class BinningStrategy;
}
namespace PARTICLE
{
  class ParticleHandler;
}
namespace BEAMINTERACTION {
namespace SUBMODELEVALUATOR{
  class Generic;
}
}
namespace STR{
namespace MODELEVALUATOR{

  // forward declaration
  class BeamInteractionDataState;

  class BeamInteraction : public Generic
  {
  public:

    typedef std::map<enum INPAR::BEAMINTERACTION::SubModelType, Teuchos::RCP<BEAMINTERACTION::SUBMODELEVALUATOR::Generic> > Map;
    typedef std::vector<Teuchos::RCP<BEAMINTERACTION::SUBMODELEVALUATOR::Generic> > Vector;

    //! constructor
    BeamInteraction();

    //! destructor
    virtual ~BeamInteraction(){};

    virtual void Setup();

    /// print welcome to biopolymer network simulation
    virtual void Logo() const;

    //! @name Derived public STR::MODELEVALUATOR::Generic methods
    //! @{
    //! derived

    //! derived
    INPAR::STR::ModelType Type() const
    { return INPAR::STR::model_beaminteraction; }

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

    //! [derived]
    virtual void Predict(const INPAR::STR::PredEnum& pred_type)
    {return;};

    //! derived
    virtual void RunPreComputeX(
            const Epetra_Vector& xold,
            Epetra_Vector& dir_mutable,
            const NOX::NLN::Group& curr_grp )
    { return; };

    //! derived
    virtual void RecoverState(
            const Epetra_Vector& xold,
            const Epetra_Vector& dir,
            const Epetra_Vector& xnew);

    //! derived
    virtual void RunPostIterate( const NOX::Solver::Generic& solver )
    { return; };

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
    virtual void OutputStepState() const;

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

    //!@name routines for submodel management
    //! @{

    //! build, init and setup submodel evaluator
    void PartitionProblem();
    //! @}

  protected:
    //! derived
    virtual void Reset(const Epetra_Vector& x);

  private:

    //!@name routines for submodel management
    //! @{

    //! set beaminteraction sub models
    void SetSubModelTypes();

    /// check if the given model type is active.
    bool HaveSubModelType(INPAR::BEAMINTERACTION::SubModelType const& submodeltype) const;

    //! build, init and setup submodel evaluator
    void InitAndSetupSubModeEvaluators();

    //! give submodels a certain order in which they are evaluated
    virtual Teuchos::RCP< STR::MODELEVALUATOR::BeamInteraction::Vector > Sort(
        STR::MODELEVALUATOR::BeamInteraction::Map submodel_map,
        std::vector<INPAR::BEAMINTERACTION::SubModelType>& sorted_submodel_types
    ) const;

    //! @}

    //!@name routines that manage to discretizations with distinct parallel distribution
    //! @{

    /// update coupling adapter and matrix transformation object with new maps
    void UpdateCouplingAdapterAndMatrixTransformation();

    /// transform force vector from ia_discret_ to Discret()
    virtual void TransformForce();

    /// transform stiffness matrix from ia_discret_ to Discret()
    virtual void TransformStiff();

    /// transform force vector and stiffness matrix from ia_discret_ to Discret()
    virtual void TransformForceStiff();

    /// update states based on bindis after its redistribution
    virtual void UpdateMaps();

    /// update map of vector (based on dof map) after redistribution
    virtual void UpdateDofMapOfVector(
        Teuchos::RCP<DRT::Discretization> discret,
        Teuchos::RCP<Epetra_Vector>&      dofmapvec,
        Teuchos::RCP<Epetra_Vector> old = Teuchos::null
    );

    /// get eletype map extractor pointer
    Teuchos::RCP<LINALG::MultiMapExtractor>& EleTypeMapExtractorPtr()
    {
      CheckInit();
      eletypeextractor_->CheckForValidMapExtractor();
      return eletypeextractor_;
    }

    //! @}

    //!@name routines that manage binning strategy
    //! @{

    /// change parallel distribution of bindis and ia_discret and assign (beam) eles to bins
    virtual void ExtendGhosting( std::map< int, std::set< int > >& extbintoelemap );

    /// build ele to bin map
    virtual void BuildRowEleToBinMap();

    /// create new bins
    virtual void CreateNewBins(bool newxaabb, bool newcutoff);

    /// print some information about binning
    virtual void PrintBinningInfoToScreen() const;

    //! @}

  private:
    //! pointer to the problem discretization (cast of base class member)
    Teuchos::RCP<DRT::Discretization> discret_ptr_;

    //!@name data for submodel management
    //! @{
    /// current active model types for the model evaluator
    Teuchos::RCP<std::set<enum INPAR::BEAMINTERACTION::SubModelType> > submodeltypes_;

    Teuchos::RCP<STR::MODELEVALUATOR::BeamInteraction::Map> me_map_ptr_;

    Teuchos::RCP<STR::MODELEVALUATOR::BeamInteraction::Vector> me_vec_ptr_;
    //! @}

    //!@name data for handling two distinct parallel distributed discretizations
    //! @{
    //! myrank
    int myrank_;

    //! coupling adapter to transfer vectors and matrices between Discret() and intactids_
    Teuchos::RCP<ADAPTER::Coupling> coupsia_;

    //! transform object for structure stiffness matrix
    Teuchos::RCP<FSI::UTILS::MatrixRowTransform> siatransform_;
    //! @}


    //!@name data for beaminteraction with binning strategy
    //! @{
    //! interaction discretization handling all interactions (e.g. crosslinker to beam,
    //! beam to beam, potential ...)
    Teuchos::RCP<DRT::Discretization> ia_discret_;

    /// map extractor for split of different element types
    Teuchos::RCP<LINALG::MultiMapExtractor> eletypeextractor_;

    //! pointer to the global state data container
    Teuchos::RCP<STR::MODELEVALUATOR::BeamInteractionDataState> ia_state_ptr_;

    //! force based on ia_discret at \f$t_{n+1}\f$
    Teuchos::RCP<Epetra_Vector> ia_force_beaminteraction_;

    //! global force based on Discret() at \f$t_{n+1}\f$
    Teuchos::RCP<Epetra_Vector> force_beaminteraction_;

    //! structural stiffness matrix based on Discret()
    Teuchos::RCP<LINALG::SparseMatrix> stiff_beaminteraction_;

    //! particle handler and binning strategy
    Teuchos::RCP<PARTICLE::ParticleHandler> particlehandler_;

    //! particle handler and binning strategy
    Teuchos::RCP<BINSTRATEGY::BinningStrategy> binstrategy_;

    //! crosslinker and bin discretization
    Teuchos::RCP<DRT::Discretization> bindis_;

    //! elerowmap of bindis
    Teuchos::RCP<Epetra_Map> rowbins_;

//    Epetra_Map* colbins_;  // Todo kept line for future improvement

    //! @}

  };
}
}


#endif
