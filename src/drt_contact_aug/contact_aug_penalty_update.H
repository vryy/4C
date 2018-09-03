/*----------------------------------------------------------------------------*/
/*!
\file contact_aug_penalty_update.H

\brief different strategies for the update/correction of the regularization
parameter cn

\level 3

\maintainer Michael Hiermeier
\date Jul 28, 2017

*/
/*----------------------------------------------------------------------------*/

#ifndef SRC_DRT_CONTACT_AUG_CONTACT_AUG_PENALTY_UPDATE_H_
#define SRC_DRT_CONTACT_AUG_CONTACT_AUG_PENALTY_UPDATE_H_

#include <Teuchos_RCP.hpp>

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
    class PenaltyUpdate
    {
     protected:
      // internal forward declaration
      struct State;

     public:
      static PenaltyUpdate* Create(const Teuchos::ParameterList& sa_params);
      static PenaltyUpdate* Create(
          const INPAR::CONTACT::PenaltyUpdate update_type, const PenaltyUpdate* pu_src = NULL);

     public:
      PenaltyUpdate() : state_(*this){/*empty*/};

      PenaltyUpdate(const PenaltyUpdate& pu) = default;

      virtual ~PenaltyUpdate(){};

      virtual INPAR::CONTACT::PenaltyUpdate Type() const = 0;

      void Init(Strategy* const strategy, DataContainer* const data);

      void Update(const CONTACT::ParamsInterface& cparams);

      virtual void SetState(const CONTACT::ParamsInterface& cparams, const Epetra_Vector& xold,
          const Epetra_Vector& dir);

      /* Step length parameter: Variations are possible, see for example:
       * "Constrained optimization and Lagrange multiplier methods",
       * Dimitri P. Bertsekas, 1996, pp.125-133
       * (interpolation strategy on the pages 132 and 133). -- hiermeier 03/17 */
      virtual double ScaleDirection(Epetra_Vector& dir) { return 1.0; };

      void PrintUpdate(std::ostream& os) const;

      void PrintInfo(std::ostream& os) const;

     protected:
      Teuchos::RCP<const LINALG::SparseMatrix> GetStructuralStiffnessMatrix(
          const CONTACT::ParamsInterface& cparams) const;

      Teuchos::RCP<const Epetra_Vector> GetProblemRhs(const CONTACT::ParamsInterface& cparams,
          const std::vector<INPAR::STR::ModelType>* without_these_models) const;

      virtual void PreUpdate(){/* empty */};

      virtual void Execute(const CONTACT::ParamsInterface& cparams) = 0;

      virtual void PostUpdate();

      void ThrowIfNotInitialized() const;

      AUG::Strategy& Strategy() { return *strategy_ptr_; };
      const AUG::Strategy& Strategy() const { return *strategy_ptr_; };
      DataContainer& Data() { return *data_ptr_; };
      const DataContainer& Data() const { return *data_ptr_; };

      /** \brief Evaluate the directional derivative of the consistently computed
       *  gradient of the weighted gap w.r.t. to the current displ. solution
       *  increment
       *
       *  \f[
       *      \langle \nabla_{\underline{d}} \underline{\tilde{g}}_{\mathrm{N}},
       *              \Delta \underline{d}_{\mathcal{S}\!\mathcal{M}}
       *  \f]
       */
      Teuchos::RCP<const Epetra_Vector> Get_DGapN(const Epetra_Vector& dincr_slma) const;

      /** \brief Evaluate the directional derivative of the potentially
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

      const State& GetState() const { return state_; }
      double& Ratio() { return ratio_; };
      double Ratio() const { return ratio_; };

     protected:
      struct State
      {
        State(PenaltyUpdate& pu) : pu_(pu){/* empty */};

        State(const State& state)
            : full_direction_(Teuchos::rcp(state.full_direction_.get(), false)),
              xold_(Teuchos::rcp(state.xold_.get(), false)),
              wgap_(Teuchos::rcp(state.wgap_.get(), false)),
              tributary_area_active_(Teuchos::rcp(state.tributary_area_active_.get(), false)),
              tributary_area_inactive_(Teuchos::rcp(state.tributary_area_inactive_.get(), false)),
              gn_gn_(state.gn_gn_),
              gn_dgn_(state.gn_dgn_),
              pu_(state.pu_){/* empty */};

        void Set(const Epetra_Vector& xold, const Epetra_Vector& dir,
            const CONTACT::AUG::DataContainer& data);

        const Epetra_Vector& GetDirection() const;

        const Epetra_Vector& GetPreviouslyAcceptedState() const;

        const Epetra_Vector& GetWGap() const;

        const Epetra_Vector& GetActiveTributaryArea() const;

        const Epetra_Vector& GetInactiveTributaryArea() const;

        void Reset();

        void Print(std::ostream& os) const;

        Teuchos::RCP<const Epetra_Vector> full_direction_ = Teuchos::null;
        Teuchos::RCP<const Epetra_Vector> xold_ = Teuchos::null;
        Teuchos::RCP<const Epetra_Vector> wgap_ = Teuchos::null;
        Teuchos::RCP<const Epetra_Vector> tributary_area_active_ = Teuchos::null;
        Teuchos::RCP<const Epetra_Vector> tributary_area_inactive_ = Teuchos::null;

        double gn_gn_ = 0.0;
        double gn_dgn_ = 0.0;

        PenaltyUpdate& pu_;
      };

     private:
      bool isinit_ = false;
      AUG::Strategy* strategy_ptr_ = NULL;
      DataContainer* data_ptr_ = NULL;

      State state_;

      double dir_norm2_ = 0.0;
      double ratio_ = 0.0;
    };

    /*--------------------------------------------------------------------------*/
    class PenaltyUpdate_Empty : public PenaltyUpdate
    {
     public:
      PenaltyUpdate_Empty(){/* empty */};

      explicit PenaltyUpdate_Empty(const PenaltyUpdate& pu) : PenaltyUpdate(pu){/* empty */};

      PenaltyUpdate_Empty(const PenaltyUpdate_Empty& pu) = delete;

      virtual ~PenaltyUpdate_Empty(){};

      virtual INPAR::CONTACT::PenaltyUpdate Type() const;

     protected:
      virtual void Execute(const CONTACT::ParamsInterface& cparams){/* empty */};
    };

    /*--------------------------------------------------------------------------*/
    class PenaltyUpdate_LagrMultiplierGapRatio : public PenaltyUpdate
    {
     public:
      PenaltyUpdate_LagrMultiplierGapRatio() = default;

      explicit PenaltyUpdate_LagrMultiplierGapRatio(const PenaltyUpdate& pu)
          : PenaltyUpdate(pu){/* empty */};

      PenaltyUpdate_LagrMultiplierGapRatio(const PenaltyUpdate_LagrMultiplierGapRatio& pu) = delete;

      virtual ~PenaltyUpdate_LagrMultiplierGapRatio(){};

      virtual INPAR::CONTACT::PenaltyUpdate Type() const;

     protected:
      virtual void Execute(const CONTACT::ParamsInterface& cparams);

     private:
      void Print2Screen(const double constr_violation, const double old_constr_violation,
          const double cn_max, const double cn_new) const;
    };

    /*--------------------------------------------------------------------------*/
    class PenaltyUpdate_Complementarity : public PenaltyUpdate
    {
     public:
      PenaltyUpdate_Complementarity() = default;

      explicit PenaltyUpdate_Complementarity(const PenaltyUpdate& pu)
          : PenaltyUpdate(pu){/* empty */};

      virtual ~PenaltyUpdate_Complementarity(){};

      virtual INPAR::CONTACT::PenaltyUpdate Type() const;

     protected:
      virtual void Execute(const CONTACT::ParamsInterface& cparams);

      virtual double ScaleDirection(Epetra_Vector& dir);

      virtual void SetState(const CONTACT::ParamsInterface& cparams, const Epetra_Vector& xold,
          const Epetra_Vector& dir);
    };

    /*--------------------------------------------------------------------------*/
    class PenaltyUpdate_SufficientLinReduction : public PenaltyUpdate
    {
     public:
      PenaltyUpdate_SufficientLinReduction() = default;

      explicit PenaltyUpdate_SufficientLinReduction(const PenaltyUpdate& pu)
          : PenaltyUpdate(pu){/* empty */};

      PenaltyUpdate_SufficientLinReduction(const PenaltyUpdate_SufficientLinReduction& pu) = delete;

      virtual ~PenaltyUpdate_SufficientLinReduction(){};

      virtual INPAR::CONTACT::PenaltyUpdate Type() const;

     protected:
      virtual void Execute(const CONTACT::ParamsInterface& cparams);

      virtual void SetState(const CONTACT::ParamsInterface& cparams, const Epetra_Vector& xold,
          const Epetra_Vector& dir);

     private:
      double GammaTheta() const;
    };

    /*--------------------------------------------------------------------------*/
    class PenaltyUpdate_SufficientAngle : public PenaltyUpdate
    {
     public:
      PenaltyUpdate_SufficientAngle() = default;

      explicit PenaltyUpdate_SufficientAngle(const PenaltyUpdate& pu)
          : PenaltyUpdate(pu){/* empty */};

      PenaltyUpdate_SufficientAngle(const PenaltyUpdate_SufficientAngle& pu) = delete;

      virtual ~PenaltyUpdate_SufficientAngle(){};

      virtual INPAR::CONTACT::PenaltyUpdate Type() const;

     protected:
      virtual void Execute(const CONTACT::ParamsInterface& cparams);

      virtual void SetState(const CONTACT::ParamsInterface& cparams, const Epetra_Vector& xold,
          const Epetra_Vector& dir);

     private:
      double GammaAngle() const;
    };

  }  // namespace AUG
}  // namespace CONTACT


#endif /* SRC_DRT_CONTACT_AUG_CONTACT_AUG_PENALTY_UPDATE_H_ */
