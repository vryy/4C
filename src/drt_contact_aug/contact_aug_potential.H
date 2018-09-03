/*---------------------------------------------------------------------*/
/*!
\file contact_aug_potential.H

\brief Class for the evaluation of the contact potential and its
       linearization.

\level 2

\maintainer Michael Hiermeier

\date Mar 14, 2017

*/
/*---------------------------------------------------------------------*/

#ifndef SRC_DRT_CONTACT_AUG_CONTACT_AUG_POTENTIAL_H_
#define SRC_DRT_CONTACT_AUG_CONTACT_AUG_POTENTIAL_H_

#include <Teuchos_RCP.hpp>
#include "../drt_lib/drt_dserror.H"

class Epetra_Vector;

namespace CONTACT
{
  namespace AUG
  {
    class Strategy;
    class DataContainer;
    namespace POTENTIAL
    {
      enum class Type : int
      {
        lagrangian,            ///< Lagrangian potential
        augmented_lagrangian,  ///< Augmented lagrangian potential
        infeasibility_measure  ///< infeasibility measure
      };

      enum class SetType : int
      {
        active,    ///< active set contributions
        inactive,  ///< inactive set contributions
        all        ///< all contributions
      };

      enum class LinTerm : int
      {
        wrt_d,        ///< linearization w.r.t. the displacements
        wrt_z,        ///< linearization w.r.t. the Lagrange multipliers
        wrt_d_and_z,  ///< linearization w.r.t. the displacements and the Lagrange multipliers
        wrt_z_and_z   ///< 2-nd derivative w.r.t. the Lagrange multipliers
      };
    }  // namespace POTENTIAL

    /** \brief Evaluate the contact potential
     *
     *  \author hiermeier \date 03/17 */
    class Potential
    {
     public:
      /// constructor
      Potential(const CONTACT::AUG::Strategy& strategy, const CONTACT::AUG::DataContainer& data);

      /// destructor
      virtual ~Potential(){};

      void Setup();

      void SetActiveInactiveState();

      /// compute the contact potential terms
      void Compute();

      /// compute the linearizations terms of the contact potential terms
      void ComputeLin(const Epetra_Vector& dir);

      double Get(enum POTENTIAL::Type pot_type, enum POTENTIAL::SetType pot_set) const;

      double GetLin(enum POTENTIAL::Type pot_type, enum POTENTIAL::SetType pot_set,
          enum POTENTIAL::LinTerm lin_term) const;

      inline const Epetra_Vector& GetZnActive() const
      {
        if (not isvalid_.state_) dserror("zn_active_ is not valid!");

        return *zn_active_;
      }

      inline const Epetra_Vector& GetZnInactive() const
      {
        if (not isvalid_.state_) dserror("zn_inactive_ is not valid!");

        return *zn_inactive_;
      }

     private:
      void ResetIsValid();

      void ComputeLinActive(const Epetra_Vector& dincrSlMa, const Epetra_Vector& znincr_active);
      void ComputeLinInactive(const Epetra_Vector& znincr_inactive);

      double GetLagrangian(enum POTENTIAL::SetType pot_set) const;
      double GetLagrangianLin(
          enum POTENTIAL::LinTerm lin_term, enum POTENTIAL::SetType pot_set) const;
      double GetActiveLagrangianLin(const enum POTENTIAL::LinTerm lin_term) const;
      double GetInactiveLagrangianLin(const enum POTENTIAL::LinTerm lin_term) const;

      double GetAugmentedLagrangian(enum POTENTIAL::SetType pot_set) const;
      double GetAugmentedLagrangianLin(
          enum POTENTIAL::LinTerm lin_term, enum POTENTIAL::SetType pot_set) const;
      double GetActiveAugmentedLagrangianLin(const enum POTENTIAL::LinTerm lin_term) const;
      double GetInactiveAugmentedLagrangianLin(const enum POTENTIAL::LinTerm lin_term) const;

      double GetInfeasibilityMeasure(enum POTENTIAL::SetType pot_set) const;
      double GetInfeasibilityMeasureLin(
          enum POTENTIAL::LinTerm lin_term, enum POTENTIAL::SetType pot_set) const;
      double GetActiveInfeasibilityMeasureLin(enum POTENTIAL::LinTerm lin_term) const;
      double GetInactiveInfeasibilityMeasureLin(enum POTENTIAL::LinTerm lin_term) const;

      void SetDirection(const Epetra_Vector& direction, Epetra_Vector& dincrSlMa,
          Epetra_Vector& znincr_active, Epetra_Vector& znincr_inactive);

     private:
      /// container for all isvalid flags
      class IsValid
      {
       public:
        IsValid() : potential_(false), linearization_(false), state_(false), dir_nrm2_(-1.0)
        { /* empty */
        }

        /** \brief return true if the same direction with the same step-length is
         *  repeatedly passed as input argument
         *
         *  This is supposed to reduce the computational cost of the directional
         *  derivative computation.
         *
         *  \param[in] dir current search direction
         *
         *  \author hiermeier \date 08/17 */
        bool isSameDirection(const Epetra_Vector& dir);

        /// true if the potential values are valid
        bool potential_;

        /// true if the directional derivative values are valid
        bool linearization_;

        /// true if the current state is valid
        bool state_;

       private:
        /** L2-norm of the step used for the last directional derivative
         * calculation */
        double dir_nrm2_;
      };

      IsValid isvalid_;

      bool issetup_;

      const CONTACT::AUG::Strategy& strategy_;
      const CONTACT::AUG::DataContainer& data_;

      Teuchos::RCP<Epetra_Vector> zn_active_;
      Teuchos::RCP<Epetra_Vector> zn_inactive_;

      Teuchos::RCP<Epetra_Vector> zt_active_;
      Teuchos::RCP<Epetra_Vector> zt_inactive_;

      /// container for the potential data
      struct PotData
      {
        PotData()
            : zn_gn_(0.0), gn_gn_(0.0), zn_zn_(0.0), zt_zt_(0.0), inf_gn_gn_(0.0), inf_zn_zn_(0.0)
        { /* empty */
        }

        void print(std::ostream& os) const;

        /// active part
        double zn_gn_;

        /// active augmented part
        double gn_gn_;

        /// inactive normal part
        double zn_zn_;

        /// tangential part
        double zt_zt_;

        /// infeasibility active square part
        double inf_gn_gn_;

        /// infeasibility inactive square part
        double inf_zn_zn_;
      };

      /// container for the linearization data
      struct LinData
      {
        LinData()
            : gn_dzn_(0.0),
              zn_dgn_(0.0),
              gn_dgn_(0.0),
              zn_dzn_(0.0),
              dzn_dgn_(0.0),
              dzn_dzn_(0.0),
              inf_gn_dgn_(0.0),
              inf_zn_dzn_(0.0)
        { /* empty */
        }

        void reset()
        {
          reset_active();
          reset_inactive();
        }

        void reset_active()
        {
          gn_dzn_ = 0.0;
          zn_dgn_ = 0.0;
          gn_dgn_ = 0.0;
          dzn_dgn_ = 0.0;
          inf_gn_dgn_ = 0.0;
        }

        void reset_inactive()
        {
          zn_dzn_ = 0.0;
          dzn_dzn_ = 0.0;
          inf_zn_dzn_ = 0.0;
        }

        void print(std::ostream& os) const;

        /// linearization of the active part w.r.t. zn
        double gn_dzn_;

        /// linearization of the active part w.r.t. d
        double zn_dgn_;

        /// linearization of the augmented part w.r.t. d
        double gn_dgn_;

        /// linearization of the inactive part w.r.t. zn
        double zn_dzn_;

        /// linearization of the active part w.r.t. d and zn
        double dzn_dgn_;

        /// 2-nd derivative of the inactive part w.r.t. zn
        double dzn_dzn_;

        /// linearization of the infeasibility active square part
        double inf_gn_dgn_;

        /// linearization of the infeasibility inactive square part
        double inf_zn_dzn_;
      };

      PotData potdata_;

      LinData lindata_;

    };  // class potential
  }     // namespace AUG
}  // namespace CONTACT


#endif /* SRC_DRT_CONTACT_AUG_CONTACT_AUG_POTENTIAL_H_ */
