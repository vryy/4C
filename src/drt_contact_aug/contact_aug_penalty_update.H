/*----------------------------------------------------------------------------*/
/*! \file
\brief different strategies for the update/correction of the regularization
parameter cn

\level 3

\maintainer Matthias Mayr
*/
/*----------------------------------------------------------------------------*/

#ifndef SRC_DRT_CONTACT_AUG_CONTACT_AUG_PENALTY_UPDATE_H_
#define SRC_DRT_CONTACT_AUG_CONTACT_AUG_PENALTY_UPDATE_H_

#include <Teuchos_RCP.hpp>

// forward declarations
class Epetra_Vector;
class Epetra_Map;

namespace INPAR
{
  namespace STR
  {
    enum ModelType : int;
  }  // namespace STR
  namespace CONTACT
  {
    enum class PenaltyUpdate : char;
  }  // namespace CONTACT
}  // namespace INPAR

namespace LINALG
{
  class SparseMatrix;
}  // namespace LINALG

namespace Teuchos
{
  class ParameterList;
}  // namespace Teuchos

namespace CONTACT
{
  class ParamsInterface;
  namespace AUG
  {
    class Strategy;
    class DataContainer;

    /*--------------------------------------------------------------------------*/
    /** \brief Base class for all strategies concerning the cn parameter correction
     *
     *  These routines have a huge impact on the performance of the SteepestAscent
     *  contact strategies.
     *
     *  \author hiermeier */
    class PenaltyUpdate
    {
     protected:
      // internal forward declaration
      struct State;

      /// current status of the cN value
      enum class Status
      {
        unevaluated,  ///< the update is unevaluated
        unchanged,    ///< cn stays unchanged
        increased,    ///< cn has been increased
        decreased     ///< cn has been decreased
      };

      /// convert to string
      inline std::string Status2String(enum Status status) const
      {
        switch (status)
        {
          case Status::unevaluated:
            return "Status::unevaluated";
          case Status::unchanged:
            return "Status::unchanged";
          case Status::increased:
            return "Status::increased";
          case Status::decreased:
            return "Status::decreased";
          default:
            return "Status::UNKNOWN";
        }
      }

     public:
      /// create a new penalty update object
      static PenaltyUpdate* Create(const Teuchos::ParameterList& sa_params);

      /** create a new penalty update object and use the data of another
       *  penalty update object to set it up. */
      static PenaltyUpdate* Create(
          const INPAR::CONTACT::PenaltyUpdate update_type, const PenaltyUpdate* pu_src = NULL);

     public:
      /// default constructor
      PenaltyUpdate() : state_(*this){/*empty*/};

      /// use default copy constructor
      PenaltyUpdate(const PenaltyUpdate& pu) = default;

      /// destructor
      virtual ~PenaltyUpdate(){};

      /// return the type enum
      virtual INPAR::CONTACT::PenaltyUpdate Type() const = 0;

      /// initialize the penalty update object
      void Init(Strategy* const strategy, DataContainer* const data);

      /// update the contact regularization parameter (increase)
      void Update(const CONTACT::ParamsInterface& cparams);

      /// update the state in this object
      virtual void SetState(const CONTACT::ParamsInterface& cparams, const Epetra_Vector& xold,
          const Epetra_Vector& dir);

      /* Step length parameter: Variations are possible, see for example:
       * "Constrained optimization and Lagrange multiplier methods",
       * Dimitri P. Bertsekas, 1996, pp.125-133
       * (interpolation strategy on the pages 132 and 133). -- hiermeier 03/17 */
      virtual double ScaleDirection(Epetra_Vector& dir) { return 1.0; };

      /// perform a possible decrease of the regularization parameter
      void Decrease(const CONTACT::ParamsInterface& cparams);

      /// print the update
      void PrintUpdate(std::ostream& os) const;

      /// print info about the stored state
      void PrintInfo(std::ostream& os) const;

     protected:
      /// access the structural stiffness matrix
      Teuchos::RCP<const LINALG::SparseMatrix> GetStructuralStiffnessMatrix(
          const CONTACT::ParamsInterface& cparams) const;

      /// access the right hand side vector of the entire problem
      Teuchos::RCP<const Epetra_Vector> GetProblemRhs(const CONTACT::ParamsInterface& cparams,
          const std::vector<INPAR::STR::ModelType>* without_these_models) const;

      /// do stuff before the update
      virtual void PreUpdate(){/* empty */};

      /// execute the update/increase
      virtual Status Execute(const CONTACT::ParamsInterface& cparams) = 0;

      /// do stuff after the update
      virtual void PostUpdate();

      /// execute the decrease
      virtual Status ExecuteDecrease(const CONTACT::ParamsInterface& cparams);

      /// do stuff after a possible decrease
      virtual void PostDecrease();

      /// reset class members
      void Reset();

      /// Throw if Init() has not been called
      void ThrowIfNotInitialized() const;

      /// access the surrounding strategy
      AUG::Strategy& Strategy() { return *strategy_ptr_; };
      const AUG::Strategy& Strategy() const { return *strategy_ptr_; };

      /// access the data container of the surrounding strategy
      DataContainer& Data() { return *data_ptr_; };
      const DataContainer& Data() const { return *data_ptr_; };

      /** Evaluate the directional derivative of the consistently computed
       *  gradient of the weighted gap w.r.t. to the current displ. solution
       *  increment
       *
       *  \f[
       *      \langle \nabla_{\underline{d}} \underline{\tilde{g}}_{\mathrm{N}},
       *              \Delta \underline{d}_{\mathcal{S}\!\mathcal{M}}
       *  \f]
       */
      Teuchos::RCP<const Epetra_Vector> Get_DGapN(const Epetra_Vector& dincr_slma) const;

      /** Evaluate the directional derivative of the potentially
       *  inconsistently computed gradient of the weighted gap w.r.t. to the
       *  current displ. solution increment
       *
       *  \f[
       *      \langle \tilde{\nabla}_{\underline{d}} \underline{\tilde{g}}_{\mathrm{N}},
       *              \Delta \underline{d}_{\mathcal{S}\!\mathcal{M}}
       *  \f]
       */
      Teuchos::RCP<const Epetra_Vector> Get_inconsistent_DGapN(
          const Epetra_Vector& dincr_slma) const;

      /// access the state container
      const State& GetState() const { return state_; }

      /// access the ratio variable (cn_new / cn_old)
      double& Ratio() { return ratio_; };
      double Ratio() const { return ratio_; };

     protected:
      /*----------------------------------------------------------------------*/
      /// internal state container
      struct State
      {
        /// constructor
        State(PenaltyUpdate& pu) : pu_(pu){/* empty */};

        /// copy constructor
        State(const State& state)
            : full_direction_(Teuchos::rcp(state.full_direction_.get(), false)),
              xold_(Teuchos::rcp(state.xold_.get(), false)),
              wgap_(Teuchos::rcp(state.wgap_.get(), false)),
              tributary_area_active_(Teuchos::rcp(state.tributary_area_active_.get(), false)),
              tributary_area_inactive_(Teuchos::rcp(state.tributary_area_inactive_.get(), false)),
              gn_gn_(state.gn_gn_),
              gn_dgn_(state.gn_dgn_),
              pu_(state.pu_){/* empty */};

        /// set state
        void Set(const Epetra_Vector& xold, const Epetra_Vector& dir,
            const CONTACT::AUG::DataContainer& data);

        /// get direction
        const Epetra_Vector& GetDirection() const;

        /// access previously accepted state
        const Epetra_Vector& GetPreviouslyAcceptedState() const;

        /// access weighted gap vector
        const Epetra_Vector& GetWGap() const;

        /// access tributary area vector (active nodes)
        const Epetra_Vector& GetActiveTributaryArea() const;

        /// access tributary area vector (inactive nodes)
        const Epetra_Vector& GetInactiveTributaryArea() const;

        /// reset the state variables
        void Reset();

        /// print the state
        void Print(std::ostream& os) const;

        /// current full direction vector
        Teuchos::RCP<const Epetra_Vector> full_direction_ = Teuchos::null;

        /// old state vector
        Teuchos::RCP<const Epetra_Vector> xold_ = Teuchos::null;

        /// weighted gap vector
        Teuchos::RCP<const Epetra_Vector> wgap_ = Teuchos::null;

        /// tributary area vector (active nodes)
        Teuchos::RCP<const Epetra_Vector> tributary_area_active_ = Teuchos::null;

        /// tributary area vector (inactive nodes)
        Teuchos::RCP<const Epetra_Vector> tributary_area_inactive_ = Teuchos::null;

        /// squared l2-norm of the active weighted gap
        double gn_gn_ = 0.0;

        /// scalar product of the weighted gap and its directional derivative
        double gn_dgn_ = 0.0;

        /// call-back
        PenaltyUpdate& pu_;
      };

     private:
      /// has Init() been called?
      bool isinit_ = false;

      /// call-back to the surrounding strategy
      AUG::Strategy* strategy_ptr_ = NULL;

      /// container of the surrounding strategy
      DataContainer* data_ptr_ = NULL;

      /// internal state object (see container description)
      State state_;

      /// status of this updating strategy
      Status status_ = Status::unevaluated;

      /// l2-norm of the current direction
      double dir_norm2_ = 0.0;

      /// ratio of new and old cn value
      double ratio_ = 0.0;
    };

    /*--------------------------------------------------------------------------*/
    /// empty correction strategy
    class PenaltyUpdate_Empty : public PenaltyUpdate
    {
     public:
      PenaltyUpdate_Empty(){/* empty */};

      explicit PenaltyUpdate_Empty(const PenaltyUpdate& pu) : PenaltyUpdate(pu){/* empty */};

      PenaltyUpdate_Empty(const PenaltyUpdate_Empty& pu) = delete;

      ~PenaltyUpdate_Empty() override{};

      INPAR::CONTACT::PenaltyUpdate Type() const override;

     protected:
      Status Execute(const CONTACT::ParamsInterface& cparams) override
      {
        return Status::unchanged;
      };
    };

    /*--------------------------------------------------------------------------*/
    /** \brief sufficient linear reduction correction strategy
     *
     *  \author hiermeier */
    class PenaltyUpdate_SufficientLinReduction : public PenaltyUpdate
    {
     public:
      PenaltyUpdate_SufficientLinReduction() = default;

      explicit PenaltyUpdate_SufficientLinReduction(const PenaltyUpdate& pu)
          : PenaltyUpdate(pu){/* empty */};

      PenaltyUpdate_SufficientLinReduction(const PenaltyUpdate_SufficientLinReduction& pu) = delete;

      ~PenaltyUpdate_SufficientLinReduction() override{};

      INPAR::CONTACT::PenaltyUpdate Type() const override;

     protected:
      Status Execute(const CONTACT::ParamsInterface& cparams) override;

      Status ExecuteDecrease(const CONTACT::ParamsInterface& cparams) override;

      void SetState(const CONTACT::ParamsInterface& cparams, const Epetra_Vector& xold,
          const Epetra_Vector& dir) override;

     private:
      /// return the \f$\beta^{c_{\mathrm{N}}}_{\Theta}\f$ parameter set by the user
      double BetaTheta() const;

      /// return the \f$\beta^{c_{\mathrm{N}}}_{\Theta\mathrm{crit}}\f$ parameter set by the user
      double BetaThetaDecrease() const;
    };

    /*--------------------------------------------------------------------------*/
    /** \brief sufficient angle correction strategy
     *
     *  \author hiermeier */
    class PenaltyUpdate_SufficientAngle : public PenaltyUpdate
    {
     public:
      PenaltyUpdate_SufficientAngle() = default;

      explicit PenaltyUpdate_SufficientAngle(const PenaltyUpdate& pu)
          : PenaltyUpdate(pu){/* empty */};

      PenaltyUpdate_SufficientAngle(const PenaltyUpdate_SufficientAngle& pu) = delete;

      ~PenaltyUpdate_SufficientAngle() override{};

      INPAR::CONTACT::PenaltyUpdate Type() const override;

     protected:
      Status Execute(const CONTACT::ParamsInterface& cparams) override;

      void SetState(const CONTACT::ParamsInterface& cparams, const Epetra_Vector& xold,
          const Epetra_Vector& dir) override;

     private:
      /// return the \f$\beta^{c_{\mathrm{N}}}_{\varphi}\f$ parameter set by the user
      double BetaAngle() const;
    };

  }  // namespace AUG
}  // namespace CONTACT


#endif /* SRC_DRT_CONTACT_AUG_CONTACT_AUG_PENALTY_UPDATE_H_ */
