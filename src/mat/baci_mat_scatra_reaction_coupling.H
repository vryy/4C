/*----------------------------------------------------------------------*/
/*! \file
 \brief helper class encapsulating the reaction terms and its derivatives

   \level 3

 *----------------------------------------------------------------------*/

#ifndef BACI_MAT_SCATRA_REACTION_COUPLING_H
#define BACI_MAT_SCATRA_REACTION_COUPLING_H

#include "baci_comm_parobjectfactory.H"
#include "baci_mat_material.H"
#include "baci_mat_par_parameter.H"
#include "baci_mat_scatra_reaction_mat.H"

namespace DRT
{
  namespace UTILS
  {
    class FunctionOfAnything;
  }  // namespace UTILS
}  // namespace DRT

namespace MAT
{
  namespace PAR
  {
    namespace REACTIONCOUPLING
    {
      //! interface class for generic reaction coupling
      class ReactionInterface
      {
       public:
        /// factory method
        static Teuchos::RCP<ReactionInterface> CreateReaction(
            MAT::PAR::reaction_coupling couplingtype,  //!< coupling type definig reaction
            bool isreacstart,                          //!< flag for reaction start feature
            const std::vector<double>& reacstart       //!< reaction start vector
        );

        /// standard constructor
        ReactionInterface(){};

        /// destructor
        virtual ~ReactionInterface() = default;
        /// initialization (to be called by derived classes)
        virtual void Initialize(int numscal,     //!< number of scalars
            const std::vector<double>& couprole  //!< coupling role vector
            ) = 0;

        /// check for initialization
        virtual bool IsInit() const = 0;

        /// helper for calculating advanced reaction terms
        virtual double CalcReaBodyForceTerm(int k,  //!< current scalar id
            int numscal,                            //!< number of scalars
            const std::vector<double>& phinp,       //!< scalar values at t_(n+1)
            const std::vector<std::pair<std::string, double>>&
                constants,  //!< vector containing values which are independent of the scalars (e.g.
                            //!< t,x,y,z)
            const std::vector<double>& couprole,  //!< coupling role vector
            double scale_reac,  //!< scaling factor for reaction term (= reaction coefficient *
                                //!< stoichometry)
            double
                scale_phi  //!< scaling factor for scalar values (used for reference concentrations)
            ) = 0;

        /// helper for calculating advanced reaction term derivatives
        virtual void CalcReaBodyForceDeriv(int k,  //!< current scalar id
            int numscal,                           //!< number of scalars
            std::vector<double>& derivs,           //!< vector with derivatives (to be filled)
            const std::vector<double>& phinp,      //!< scalar values at t_(n+1)
            const std::vector<std::pair<std::string, double>>&
                constants,  //!< vector containing values which are independent of the scalars (e.g.
                            //!< t,x,y,z)
            const std::vector<double>& couprole,  //!< coupling role vector
            double scale_reac,  //!< scaling factor for reaction term (= reaction coefficient *
                                //!< stoichometry)
            double
                scale_phi  //!< scaling factor for scalar values (used for reference concentrations)
            ) = 0;

        /// add additional variables for by-function reaction
        virtual void AddAdditionalVariables(const int k,                   //!< current scalar id
            const std::vector<std::pair<std::string, double>>& variables,  //!< variables
            const std::vector<double>& couprole                            //!< coupling role vector
        )
        {
          // do nothing in this case --> only the by-function reaction will override this method
          // if you want to evaluate any coupling terms in your own reaction coupling just override
          // it in your own function
          dserror("Only the by-function coupling is capable of adding additional variables");
          return;
        }

        /// helper for calculating advanced reaction term derivatives w.r.t. additional variables
        virtual void CalcReaBodyForceDerivAddVariables(const int k,  //!< current scalar id
            std::vector<double>& derivs,  //!< vector with derivatives (to be filled)
            const std::vector<std::pair<std::string, double>>& variables,  //!< variables
            const std::vector<std::pair<std::string, double>>&
                constants,                        //!< constants (including scalar values phinp)
            const std::vector<double>& couprole,  //!< coupling role vector
            double scale_reac,  //!< scaling factor for reaction term (= reaction coefficient *
                                //!< stoichometry)
            double
                scale_phi  //!< scaling factor for scalar values (used for reference concentrations)
        )
        {
          // do nothing in this case --> only the by-function reaction will override this method
          // if you want to evaluate any coupling terms in your own reaction coupling just override
          // it in your own function
          dserror("Only the by-function coupling is capable of calculating additional derivatives");
          return;
        }
      };

      /*----------------------------------------------------------------------*
       * **********************************************************************
       *----------------------------------------------------------------------*/

      //! reaction start feature for reaction coupling
      //! it wraps another reaction class
      class ReacStart : public ReactionInterface
      {
       public:
        /// standard constructor
        ReacStart(Teuchos::RCP<ReactionInterface> reaction, const std::vector<double>& reacstart)
            : reaction_(reaction), reacstart_(reacstart){};

        /// initialization
        void Initialize(int numscal,             //!< number of scalars
            const std::vector<double>& couprole  //!< coupling role vector
            ) override;

        /// check for initialization
        bool IsInit() const override { return reaction_->IsInit(); };

        /// helper for calculating advanced reaction terms
        double CalcReaBodyForceTerm(int k,     //!< current scalar id
            int numscal,                       //!< number of scalars
            const std::vector<double>& phinp,  //!< scalar values at t_(n+1)
            const std::vector<std::pair<std::string, double>>&
                constants,  //!< vector containing values which are independent of the scalars (e.g.
                            //!< t,x,y,z)
            const std::vector<double>& couprole,  //!< coupling role vector
            double scale_reac,  //!< scaling factor for reaction term (= reaction coefficient *
                                //!< stoichometry)
            double
                scale_phi  //!< scaling factor for scalar values (used for reference concentrations)
            ) override;

        /// helper for calculating advanced reaction term derivatives
        void CalcReaBodyForceDeriv(int k,      //!< current scalar id
            int numscal,                       //!< number of scalars
            std::vector<double>& derivs,       //!< vector with derivatives (to be filled)
            const std::vector<double>& phinp,  //!< scalar values at t_(n+1)
            const std::vector<std::pair<std::string, double>>&
                constants,  //!< vector containing values which are independent of the scalars (e.g.
                            //!< t,x,y,z)
            const std::vector<double>& couprole,  //!< coupling role vector
            double scale_reac,  //!< scaling factor for reaction term (= reaction coefficient *
                                //!< stoichometry)
            double
                scale_phi  //!< scaling factor for scalar values (used for reference concentrations)
            ) override;

       private:
        /// helper for reaction start feature
        virtual std::vector<double> ModifyPhi(const std::vector<double>& phinp);

        /// actual reaction
        Teuchos::RCP<ReactionInterface> reaction_;
        /// reacstart vector
        const std::vector<double>& reacstart_;
      };

      /*----------------------------------------------------------------------*
       * **********************************************************************
       *----------------------------------------------------------------------*/

      //! wrapper class for reaction coupling with potential phi scaling
      //! it wraps another reaction class
      class ReactionWithPhiScaling : public ReactionInterface
      {
       public:
        /// standard constructor
        ReactionWithPhiScaling(Teuchos::RCP<ReactionInterface> reaction) : reaction_(reaction){};

        /// initialization (to be called by derived classes)
        void Initialize(int numscal,             //!< number of scalars
            const std::vector<double>& couprole  //!< coupling role vector
            ) override
        {
          reaction_->Initialize(numscal, couprole);
        };

        /// check for initialization
        bool IsInit() const override { return reaction_->IsInit(); };

        /// helper for calculating advanced reaction terms
        double CalcReaBodyForceTerm(int k,     //!< current scalar id
            int numscal,                       //!< number of scalars
            const std::vector<double>& phinp,  //!< scalar values at t_(n+1)
            const std::vector<std::pair<std::string, double>>&
                constants,  //!< vector containing values which are independent of the scalars (e.g.
                            //!< t,x,y,z)
            const std::vector<double>& couprole,  //!< coupling role vector
            double scale_reac,  //!< scaling factor for reaction term (= reaction coefficient *
                                //!< stoichometry)
            double
                scale_phi  //!< scaling factor for scalar values (used for reference concentrations)
            ) override;

        /// helper for calculating advanced reaction term derivatives
        void CalcReaBodyForceDeriv(int k,      //!< current scalar id
            int numscal,                       //!< number of scalars
            std::vector<double>& derivs,       //!< vector with derivatives (to be filled)
            const std::vector<double>& phinp,  //!< scalar values at t_(n+1)
            const std::vector<std::pair<std::string, double>>&
                constants,  //!< vector containing values which are independent of the scalars (e.g.
                            //!< t,x,y,z)
            const std::vector<double>& couprole,  //!< coupling role vector
            double scale_reac,  //!< scaling factor for reaction term (= reaction coefficient *
                                //!< stoichometry)
            double
                scale_phi  //!< scaling factor for scalar values (used for reference concentrations)
            ) override;

        /// helper for calculating advanced reaction term derivatives after additional variables
        void CalcReaBodyForceDerivAddVariables(const int k,  //!< current scalar id
            std::vector<double>& derivs,  //!< vector with derivatives (to be filled)
            const std::vector<std::pair<std::string, double>>& variables,  //!< variables
            const std::vector<std::pair<std::string, double>>&
                constants,                        //!< constants (including scalar values phinp)
            const std::vector<double>& couprole,  //!< coupling role vector
            double scale_reac,  //!< scaling factor for reaction term (= reaction coefficient *
                                //!< stoichometry)
            double
                scale_phi  //!< scaling factor for scalar values (used for reference concentrations)
            ) override;

        /// add additional variables for by-function reaction
        void AddAdditionalVariables(const int k,                           //!< current scalar id
            const std::vector<std::pair<std::string, double>>& variables,  //!< variables
            const std::vector<double>& couprole                            //!< coupling role vector
            ) override
        {
          reaction_->AddAdditionalVariables(k, variables, couprole);
        }

       protected:
        /// actual reaction
        Teuchos::RCP<ReactionInterface> reaction_;

        /// helper for scaling
        virtual std::vector<double> ModifyPhi(const std::vector<double>& phinp, double scale_phi);
      };

      /*----------------------------------------------------------------------*
       * **********************************************************************
       *----------------------------------------------------------------------*/

      //! base class for reaction coupling kinetics
      class ReactionBase : public ReactionInterface
      {
       public:
        /// standard constructor
        ReactionBase() : isinit_(false){};

        /// initialization (to be called by derived classes)
        void Initialize(int numscal,             //!< number of scalars
            const std::vector<double>& couprole  //!< coupling role vector
            ) override
        {
          isinit_ = true;
        };

        /// check for initialization
        bool IsInit() const override { return isinit_; };

        /// helper for calculating advanced reaction terms
        double CalcReaBodyForceTerm(int k,     //!< current scalar id
            int numscal,                       //!< number of scalars
            const std::vector<double>& phinp,  //!< scalar values at t_(n+1)
            const std::vector<std::pair<std::string, double>>&
                constants,  //!< vector containing values which are independent of the scalars (e.g.
                            //!< t,x,y,z)
            const std::vector<double>& couprole,  //!< coupling role vector
            double scale_reac,  //!< scaling factor for reaction term (= reaction coefficient *
                                //!< stoichometry)
            double
                scale_phi  //!< scaling factor for scalar values (used for reference concentrations)
            ) override
        {
          // check
          dsassert(IsInit(), "Reaction class has not been initialized!");

          // call the real evaluation (scale_phi should have been applied in wrapper class)
          return CalcReaBodyForceTerm(k, numscal, phinp, constants, couprole, scale_reac);
        };

        /// helper for calculating advanced reaction term derivatives
        void CalcReaBodyForceDeriv(int k,      //!< current scalar id
            int numscal,                       //!< number of scalars
            std::vector<double>& derivs,       //!< vector with derivatives (to be filled)
            const std::vector<double>& phinp,  //!< scalar values at t_(n+1)
            const std::vector<std::pair<std::string, double>>&
                constants,  //!< vector containing values which are independent of the scalars (e.g.
                            //!< t,x,y,z)
            const std::vector<double>& couprole,  //!< coupling role vector
            double scale_reac,  //!< scaling factor for reaction term (= reaction coefficient *
                                //!< stoichometry)
            double
                scale_phi  //!< scaling factor for scalar values (used for reference concentrations)
            ) override
        {
          // check
          dsassert(IsInit(), "Reaction class has not been initialized!");

          // call the real evaluation (scale_phi should have been applied in wrapper class)
          CalcReaBodyForceDeriv(k, numscal, derivs, phinp, constants, couprole, scale_reac);
          return;
        };

       private:
        /// helper for calculating advanced reaction terms
        virtual double CalcReaBodyForceTerm(int k,  //!< current scalar id
            int numscal,                            //!< number of scalars
            const std::vector<double>& phinp,       //!< scalar values at t_(n+1)
            const std::vector<std::pair<std::string, double>>&
                constants,  //!< vector containing values which are independent of the scalars (e.g.
                            //!< t,x,y,z)
            const std::vector<double>& couprole,  //!< coupling role vector
            double scale_reac  //!< scaling factor for reaction term (= reaction coefficient *
                               //!< stoichometry)
            ) = 0;

        /// helper for calculating advanced reaction term derivatives
        virtual void CalcReaBodyForceDeriv(int k,  //!< current scalar id
            int numscal,                           //!< number of scalars
            std::vector<double>& derivs,           //!< vector with derivatives (to be filled)
            const std::vector<double>& phinp,      //!< scalar values at t_(n+1)
            const std::vector<std::pair<std::string, double>>&
                constants,  //!< vector containing values which are independent of the scalars (e.g.
                            //!< t,x,y,z)
            const std::vector<double>& couprole,  //!< coupling role vector
            double scale_reac  //!< scaling factor for reaction term (= reaction coefficient *
                               //!< stoichometry)
            ) = 0;

       private:
        // initialization flag
        bool isinit_;
      };

      /*----------------------------------------------------------------------*
       * **********************************************************************
       *----------------------------------------------------------------------*/

      //! simple multiplicative reaction coupling
      class SimpleMultiplicative : public ReactionBase
      {
       public:
        /// standard constructor
        SimpleMultiplicative(){};

        /// initialization
        void Initialize(int numscal,             //!< number of scalars
            const std::vector<double>& couprole  //!< coupling role vector
            ) override;

       private:
        /// helper for calculating advanced reaction terms
        double CalcReaBodyForceTerm(int k,     //!< current scalar id
            int numscal,                       //!< number of scalars
            const std::vector<double>& phinp,  //!< scalar values at t_(n+1)
            const std::vector<std::pair<std::string, double>>&
                constants,  //!< vector containing values which are independent of the scalars (e.g.
                            //!< t,x,y,z)
            const std::vector<double>& couprole,  //!< coupling role vector
            double scale_reac  //!< scaling factor for reaction term (= reaction coefficient *
                               //!< stoichometry)
            ) override;

        /// helper for calculating advanced reaction term derivatives
        void CalcReaBodyForceDeriv(int k,      //!< current scalar id
            int numscal,                       //!< number of scalars
            std::vector<double>& derivs,       //!< vector with derivatives (to be filled)
            const std::vector<double>& phinp,  //!< scalar values at t_(n+1)
            const std::vector<std::pair<std::string, double>>&
                constants,  //!< vector containing values which are independent of the scalars (e.g.
                            //!< t,x,y,z)
            const std::vector<double>& couprole,  //!< coupling role vector
            double scale_reac  //!< scaling factor for reaction term (= reaction coefficient *
                               //!< stoichometry)
            ) override;
      };

      /*----------------------------------------------------------------------*
       * **********************************************************************
       *----------------------------------------------------------------------*/

      //! power multiplicative reaction coupling
      class PowerMultiplicative : public ReactionBase
      {
       public:
        /// standard constructor
        PowerMultiplicative(){};

        /// initialization
        void Initialize(int numscal,             //!< number of scalars
            const std::vector<double>& couprole  //!< coupling role vector
            ) override;

        /// helper for calculating advanced reaction terms
        double CalcReaBodyForceTerm(int k,     //!< current scalar id
            int numscal,                       //!< number of scalars
            const std::vector<double>& phinp,  //!< scalar values at t_(n+1)
            const std::vector<std::pair<std::string, double>>&
                constants,  //!< vector containing values which are independent of the scalars (e.g.
                            //!< t,x,y,z)
            const std::vector<double>& couprole,  //!< coupling role vector
            double scale_reac  //!< scaling factor for reaction term (= reaction coefficient *
                               //!< stoichometry)
            ) override;

        /// helper for calculating advanced reaction term derivatives
        void CalcReaBodyForceDeriv(int k,      //!< current scalar id
            int numscal,                       //!< number of scalars
            std::vector<double>& derivs,       //!< vector with derivatives (to be filled)
            const std::vector<double>& phinp,  //!< scalar values at t_(n+1)
            const std::vector<std::pair<std::string, double>>&
                constants,  //!< vector containing values which are independent of the scalars (e.g.
                            //!< t,x,y,z)
            const std::vector<double>& couprole,  //!< coupling role vector
            double scale_reac  //!< scaling factor for reaction term (= reaction coefficient *
                               //!< stoichometry)
            ) override;
      };

      /*----------------------------------------------------------------------*
       * **********************************************************************
       *----------------------------------------------------------------------*/

      //! simple multiplicative reaction coupling
      class Constant : public ReactionBase
      {
       public:
        /// standard constructor
        Constant(){};

        /// initialization
        void Initialize(int numscal,             //!< number of scalars
            const std::vector<double>& couprole  //!< coupling role vector
            ) override;

       private:
        /// helper for calculating advanced reaction terms
        double CalcReaBodyForceTerm(int k,     //!< current scalar id
            int numscal,                       //!< number of scalars
            const std::vector<double>& phinp,  //!< scalar values at t_(n+1)
            const std::vector<std::pair<std::string, double>>&
                constants,  //!< vector containing values which are independent of the scalars (e.g.
                            //!< t,x,y,z)
            const std::vector<double>& couprole,  //!< coupling role vector
            double scale_reac  //!< scaling factor for reaction term (= reaction coefficient *
                               //!< stoichometry)
            ) override;

        /// helper for calculating advanced reaction term derivatives
        void CalcReaBodyForceDeriv(const int k,  //!< current scalar id
            int numscal,                         //!< number of scalars
            std::vector<double>& derivs,         //!< vector with derivatives (to be filled)
            const std::vector<double>& phinp,    //!< scalar values at t_(n+1)
            const std::vector<std::pair<std::string, double>>&
                constants,  //!< vector containing values which are independent of the scalars (e.g.
                            //!< t,x,y,z)
            const std::vector<double>& couprole,  //!< coupling role vector
            double scale_reac  //!< scaling factor for reaction term (= reaction coefficient *
                               //!< stoichometry)
            ) override;
      };

      /*----------------------------------------------------------------------*
       * **********************************************************************
       *----------------------------------------------------------------------*/

      //! simple multiplicative reaction coupling
      class MichaelisMenten : public ReactionBase
      {
       public:
        /// standard constructor
        MichaelisMenten(){};

        /// initialization
        void Initialize(int numscal,             //!< number of scalars
            const std::vector<double>& couprole  //!< coupling role vector
            ) override;

       private:
        /// helper for calculating advanced reaction terms
        double CalcReaBodyForceTerm(const int k,  //!< current scalar id
            int numscal,                          //!< number of scalars
            const std::vector<double>& phinp,     //!< scalar values at t_(n+1)
            const std::vector<std::pair<std::string, double>>&
                constants,  //!< vector containing values which are independent of the scalars (e.g.
                            //!< t,x,y,z)
            const std::vector<double>& couprole,  //!< coupling role vector
            double scale_reac  //!< scaling factor for reaction term (= reaction coefficient *
                               //!< stoichometry)
            ) override;

        /// helper for calculating advanced reaction term derivatives
        void CalcReaBodyForceDeriv(const int k,  //!< current scalar id
            int numscal,                         //!< number of scalars
            std::vector<double>& derivs,         //!< vector with derivatives (to be filled)
            const std::vector<double>& phinp,    //!< scalar values at t_(n+1)
            const std::vector<std::pair<std::string, double>>&
                constants,  //!< vector containing values which are independent of the scalars (e.g.
                            //!< t,x,y,z)
            const std::vector<double>& couprole,  //!< coupling role vector
            double scale_reac  //!< scaling factor for reaction term (= reaction coefficient *
                               //!< stoichometry)
            ) override;
      };

      /*----------------------------------------------------------------------*
       * **********************************************************************
       *----------------------------------------------------------------------*/

      //! simple multiplicative reaction coupling
      class ByFunction : public ReactionBase
      {
       public:
        /// standard constructor
        ByFunction(){};

        /// initialization
        void Initialize(int numscal,             //!< number of scalars
            const std::vector<double>& couprole  //!< coupling role vector
            ) override;

        /// helper for calculating advanced reaction term derivatives after additional variables
        void CalcReaBodyForceDerivAddVariables(const int k,  //!< current scalar id
            std::vector<double>& derivs,  //!< vector with derivatives (to be filled)
            const std::vector<std::pair<std::string, double>>& variables,  //!< variables
            const std::vector<std::pair<std::string, double>>&
                constants,                        //!< constants (including scalar values phinp)
            const std::vector<double>& couprole,  //!< coupling role vector
            double scale_reac,  //!< scaling factor for reaction term (= reaction coefficient *
                                //!< stoichometry)
            double
                scale_phi  //!< scaling factor for scalar values (used for reference concentrations)
            ) override
        {
          // check
          dsassert(IsInit(), "Reaction class has not been initialized!");

          // call the real evaluation (scale_phi should have been applied in wrapper class)
          CalcReaBodyForceDerivAddVariables(k, derivs, variables, constants, couprole, scale_reac);
        }

        /// add additional variables for by-function reaction
        void AddAdditionalVariables(const int k,                           //!< current scalar id
            const std::vector<std::pair<std::string, double>>& variables,  //!< variables
            const std::vector<double>& couprole                            //!< coupling role vector
            ) override;

       private:
        /// helper for calculating advanced reaction terms
        double CalcReaBodyForceTerm(const int k,  //!< current scalar id
            int numscal,                          //!< number of scalars
            const std::vector<double>& phinp,     //!< scalar values at t_(n+1)
            const std::vector<std::pair<std::string, double>>&
                constants,  //!< vector containing values which are independent of the scalars (e.g.
                            //!< t,x,y,z)
            const std::vector<double>& couprole,  //!< coupling role vector
            double scale_reac  //!< scaling factor for reaction term (= reaction coefficient *
                               //!< stoichometry)
            ) override;

        /// helper for calculating advanced reaction term derivatives
        void CalcReaBodyForceDeriv(const int k,  //!< current scalar id
            int numscal,                         //!< number of scalars
            std::vector<double>& derivs,         //!< vector with derivatives (to be filled)
            const std::vector<double>& phinp,    //!< scalar values at t_(n+1)
            const std::vector<std::pair<std::string, double>>&
                constants,  //!< vector containing values which are independent of the scalars (e.g.
                            //!< t,x,y,z)
            const std::vector<double>& couprole,  //!< coupling role vector
            double scale_reac  //!< scaling factor for reaction term (= reaction coefficient *
                               //!< stoichometry)
            ) override;

        /// helper for calculating advanced reaction term derivatives after additional variables
        void CalcReaBodyForceDerivAddVariables(const int k,  //!< current scalar id
            std::vector<double>& derivs,  //!< vector with derivatives (to be filled)
            const std::vector<std::pair<std::string, double>>& variables,  //!< variables
            const std::vector<std::pair<std::string, double>>&
                constants,                        //!< constants (including scalar values phinp)
            const std::vector<double>& couprole,  //!< coupling role vector
            double scale_reac  //!< scaling factor for reaction term (= reaction coefficient *
                               //!< stoichometry)
        );

        /// helper for evaluation by function
        void BuildPhiVectorForFunction(
            const std::vector<double>& phinp_org,  //!< scalar values at t_(n+1)
            int numscal                            //!< number of scalars
        );


        //! templated internal AddAdditionalVariables implementation
        template <int dim>
        void AddAdditionalVariablesInternal(const int k,                   //!< current scalar id
            const std::vector<std::pair<std::string, double>>& variables,  //!< variables
            const std::vector<double>& couprole                            //!< coupling role vector
        );

        //! templated internal CalcReaBodyForceDerivAddVariables implementation
        template <int dim>
        void CalcReaBodyForceDerivAddVariablesInternal(const int k,  //!< current scalar id
            std::vector<double>& derivs,  //!< vector with derivatives (to be filled)
            const std::vector<std::pair<std::string, double>>& variables,  //!< variables
            const std::vector<std::pair<std::string, double>>&
                constants,                        //!< constants (including scalar values phinp)
            const std::vector<double>& couprole,  //!< coupling role vector
            double scale_reac  //!< scaling factor for reaction term (= reaction coefficient *
                               //!< stoichometry)
        );

        //! templated internal CalcReaBodyForceDeriv implementation
        template <int dim>
        void CalcReaBodyForceDerivInternal(int k,  //!< current scalar id
            int numscal,                           //!< number of scalars
            std::vector<double>& derivs,           //!< vector with derivatives (to be filled)
            const std::vector<double>& phinp,      //!< scalar values at t_(n+1)
            const std::vector<std::pair<std::string, double>>&
                constants,  //!< vector containing values which are independent of the scalars (e.g.
                            //!< t,x,y,z)
            const std::vector<double>& couprole,  //!< coupling role vector
            double scale_reac  //!< scaling factor for reaction term (= reaction coefficient *
                               //!< stoichometry)
        );

        //! templated internal CalcReaBodyForceTerm implementation
        template <int dim>
        double CalcReaBodyForceTermInternal(int k,  //!< current scalar id
            int numscal,                            //!< number of scalars
            const std::vector<double>& phinp,       //!< scalar values at t_(n+1)
            const std::vector<std::pair<std::string, double>>&
                constants,  //!< vector containing values which are independent of the scalars (e.g.
                            //!< t,x,y,z)
            const std::vector<double>& couprole,  //!< coupling role vector
            double scale_reac  //!< scaling factor for reaction term (= reaction coefficient *
                               //!< stoichometry)
        );

        //! templated internal Initialize implementation
        template <int dim>
        void InitializeInternal(int numscal,     //!< number of scalars
            const std::vector<double>& couprole  //!< coupling role vector
        );

        /// variable vector for function evaluation
        std::vector<std::pair<std::string, double>> variables_;
      };

    }  // namespace REACTIONCOUPLING

  }  // namespace PAR
}  // namespace MAT


#endif  // MAT_SCATRA_REACTION_COUPLING_H
