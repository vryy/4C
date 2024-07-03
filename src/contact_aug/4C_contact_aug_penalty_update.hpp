/*----------------------------------------------------------------------------*/
/*! \file
\brief different strategies for the update/correction of the regularization
parameter cn

\level 3

*/
/*----------------------------------------------------------------------------*/

#ifndef FOUR_C_CONTACT_AUG_PENALTY_UPDATE_HPP
#define FOUR_C_CONTACT_AUG_PENALTY_UPDATE_HPP

#include "4C_config.hpp"

#include <Teuchos_RCP.hpp>

#include <vector>

// forward declarations
class Epetra_Vector;
class Epetra_Map;
namespace Teuchos
{
  class ParameterList;
}  // namespace Teuchos

FOUR_C_NAMESPACE_OPEN

namespace Inpar
{
  namespace Solid
  {
    enum ModelType : int;
  }  // namespace Solid
  namespace CONTACT
  {
    enum class PenaltyUpdate : char;
  }  // namespace CONTACT
}  // namespace Inpar

namespace Core::LinAlg
{
  class SparseMatrix;
}  // namespace Core::LinAlg

namespace CONTACT
{
  class ParamsInterface;
  namespace Aug
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
      inline std::string status2_string(enum Status status) const
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
          const Inpar::CONTACT::PenaltyUpdate update_type, const PenaltyUpdate* pu_src = nullptr);

     public:
      /// default constructor
      PenaltyUpdate() : state_(*this){/*empty*/};

      /// use default copy constructor
      PenaltyUpdate(const PenaltyUpdate& pu) = default;

      /// destructor
      virtual ~PenaltyUpdate() = default;

      /// return the type enum
      virtual Inpar::CONTACT::PenaltyUpdate Type() const = 0;

      /// initialize the penalty update object
      void init(Strategy* const strategy, DataContainer* const data);

      /// update the contact regularization parameter (increase)
      void Update(const CONTACT::ParamsInterface& cparams);

      /// update the state in this object
      virtual void set_state(const CONTACT::ParamsInterface& cparams, const Epetra_Vector& xold,
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
      Teuchos::RCP<const Core::LinAlg::SparseMatrix> get_structural_stiffness_matrix(
          const CONTACT::ParamsInterface& cparams) const;

      /// access the right hand side vector of the entire problem
      Teuchos::RCP<const Epetra_Vector> get_problem_rhs(const CONTACT::ParamsInterface& cparams,
          const std::vector<Inpar::Solid::ModelType>* without_these_models) const;

      /// do stuff before the update
      virtual void pre_update(){/* empty */};

      /// execute the update/increase
      virtual Status execute(const CONTACT::ParamsInterface& cparams) = 0;

      /// do stuff after the update
      virtual void post_update();

      /// execute the decrease
      virtual Status execute_decrease(const CONTACT::ParamsInterface& cparams);

      /// do stuff after a possible decrease
      virtual void post_decrease();

      /// reset class members
      void reset();

      /// Throw if init() has not been called
      void throw_if_not_initialized() const;

      /// access the surrounding strategy
      Aug::Strategy& strategy() { return *strategy_ptr_; };
      const Aug::Strategy& strategy() const { return *strategy_ptr_; };

      /// access the data container of the surrounding strategy
      DataContainer& data() { return *data_ptr_; };
      const DataContainer& data() const { return *data_ptr_; };

      /** Evaluate the directional derivative of the consistently computed
       *  gradient of the weighted gap w.r.t. to the current displ. solution
       *  increment
       *
       *  \f[
       *      \langle \nabla_{\underline{d}} \underline{\tilde{g}}_{\mathrm{N}},
       *              \Delta \underline{d}_{\mathcal{S}\!\mathcal{M}}
       *  \f]
       */
      Teuchos::RCP<const Epetra_Vector> get_d_gap_n(const Epetra_Vector& dincr_slma) const;

      /** Evaluate the directional derivative of the potentially
       *  inconsistently computed gradient of the weighted gap w.r.t. to the
       *  current displ. solution increment
       *
       *  \f[
       *      \langle \tilde{\nabla}_{\underline{d}} \underline{\tilde{g}}_{\mathrm{N}},
       *              \Delta \underline{d}_{\mathcal{S}\!\mathcal{M}}
       *  \f]
       */
      Teuchos::RCP<const Epetra_Vector> get_inconsistent_d_gap_n(
          const Epetra_Vector& dincr_slma) const;

      /// access the state container
      const State& get_state() const { return state_; }

      /// access the ratio variable (cn_new / cn_old)
      double& ratio() { return ratio_; };
      double ratio() const { return ratio_; };

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
        void set(const Epetra_Vector& xold, const Epetra_Vector& dir,
            const CONTACT::Aug::DataContainer& data);

        /// get direction
        const Epetra_Vector& GetDirection() const;

        /// access previously accepted state
        const Epetra_Vector& get_previously_accepted_state() const;

        /// access weighted gap vector
        const Epetra_Vector& GetWGap() const;

        /// access tributary area vector (active nodes)
        const Epetra_Vector& get_active_tributary_area() const;

        /// access tributary area vector (inactive nodes)
        const Epetra_Vector& get_inactive_tributary_area() const;

        /// reset the state variables
        void reset();

        /// print the state
        void print(std::ostream& os) const;

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
      /// has init() been called?
      bool isinit_ = false;

      /// call-back to the surrounding strategy
      Aug::Strategy* strategy_ptr_ = nullptr;

      /// container of the surrounding strategy
      DataContainer* data_ptr_ = nullptr;

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
    class PenaltyUpdateEmpty : public PenaltyUpdate
    {
     public:
      PenaltyUpdateEmpty(){/* empty */};

      explicit PenaltyUpdateEmpty(const PenaltyUpdate& pu) : PenaltyUpdate(pu){/* empty */};

      PenaltyUpdateEmpty(const PenaltyUpdateEmpty& pu) = delete;

      Inpar::CONTACT::PenaltyUpdate Type() const override;

     protected:
      Status execute(const CONTACT::ParamsInterface& cparams) override
      {
        return Status::unchanged;
      };
    };

    /*--------------------------------------------------------------------------*/
    /** \brief sufficient linear reduction correction strategy
     *
     *  \author hiermeier */
    class PenaltyUpdateSufficientLinReduction : public PenaltyUpdate
    {
     public:
      PenaltyUpdateSufficientLinReduction() = default;

      explicit PenaltyUpdateSufficientLinReduction(const PenaltyUpdate& pu)
          : PenaltyUpdate(pu){/* empty */};

      PenaltyUpdateSufficientLinReduction(const PenaltyUpdateSufficientLinReduction& pu) = delete;

      Inpar::CONTACT::PenaltyUpdate Type() const override;

     protected:
      Status execute(const CONTACT::ParamsInterface& cparams) override;

      Status execute_decrease(const CONTACT::ParamsInterface& cparams) override;

      void set_state(const CONTACT::ParamsInterface& cparams, const Epetra_Vector& xold,
          const Epetra_Vector& dir) override;

     private:
      /// return the \f$\beta^{c_{\mathrm{N}}}_{\Theta}\f$ parameter set by the user
      double beta_theta() const;

      /// return the \f$\beta^{c_{\mathrm{N}}}_{\Theta\mathrm{crit}}\f$ parameter set by the user
      double beta_theta_decrease() const;
    };

    /*--------------------------------------------------------------------------*/
    /** \brief sufficient angle correction strategy
     *
     *  \author hiermeier */
    class PenaltyUpdateSufficientAngle : public PenaltyUpdate
    {
     public:
      PenaltyUpdateSufficientAngle() = default;

      explicit PenaltyUpdateSufficientAngle(const PenaltyUpdate& pu)
          : PenaltyUpdate(pu){/* empty */};

      PenaltyUpdateSufficientAngle(const PenaltyUpdateSufficientAngle& pu) = delete;

      Inpar::CONTACT::PenaltyUpdate Type() const override;

     protected:
      Status execute(const CONTACT::ParamsInterface& cparams) override;

      void set_state(const CONTACT::ParamsInterface& cparams, const Epetra_Vector& xold,
          const Epetra_Vector& dir) override;

     private:
      /// return the \f$\beta^{c_{\mathrm{N}}}_{\varphi}\f$ parameter set by the user
      double beta_angle() const;
    };

  }  // namespace Aug
}  // namespace CONTACT


FOUR_C_NAMESPACE_CLOSE

#endif
