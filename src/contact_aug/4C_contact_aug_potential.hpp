/*---------------------------------------------------------------------*/
/*! \file
\brief Class for the evaluation of the contact potential and its
       linearization.

\level 2

*/
/*---------------------------------------------------------------------*/

#ifndef FOUR_C_CONTACT_AUG_POTENTIAL_HPP
#define FOUR_C_CONTACT_AUG_POTENTIAL_HPP

#include "4C_config.hpp"

#include "4C_utils_exceptions.hpp"

#include <Teuchos_RCP.hpp>

class Epetra_Vector;

FOUR_C_NAMESPACE_OPEN

namespace CONTACT
{
  namespace Aug
  {
    class Strategy;
    class DataContainer;
    namespace POTENTIAL
    {
      /// supported potential/function types
      enum class Type : int
      {
        lagrangian,            ///< Lagrangian potential
        augmented_lagrangian,  ///< Augmented lagrangian potential
        infeasibility_measure  ///< infeasibility measure
      };

      /// set affiliation
      enum class SetType : int
      {
        active,    ///< active set contributions
        inactive,  ///< inactive set contributions
        all        ///< all contributions
      };

      /// linearization term type
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
     *  This class computes the objective/merit function contributions for a variety
     *  of difference functions. The computed contributions, either for the function
     *  value itself or for a linear model, are stored in respective data containers
     *  and can be accessed by different accessors.
     *
     *  \author hiermeier \date 03/17 */
    class Potential
    {
     public:
      /// constructor
      Potential(const CONTACT::Aug::Strategy& strategy, const CONTACT::Aug::DataContainer& data);

      /// destructor
      virtual ~Potential() = default;

      /// setup the potential class
      void setup();

      /// set active and inactive Lagrange multiplier state vectors
      void set_active_inactive_state();

      /// compute the contact potential terms
      void Compute();

      /// compute the linearizations terms of the contact potential terms
      void ComputeLin(const Epetra_Vector& dir);

      /** \brief Return the potential value of \c pot_type belonging to \c pot_set
       *
       *  \param pot_type  The desired potential type
       *  \param set_type  Specifies the desired set (active, inactive, all/both, ...)
       *
       *  \author hiermeier */
      double get(enum POTENTIAL::Type pot_type, enum POTENTIAL::SetType pot_set) const;

      /** \brief Return a part of the linear potential model
       *
       *  \param pot_type  potential type
       *  \param set_type  set type (active, inactive, ...)
       *  \param lin_term  linearization w.r.t. which quantity?
       *
       *  \author hiermeier */
      double GetLin(enum POTENTIAL::Type pot_type, enum POTENTIAL::SetType pot_set,
          enum POTENTIAL::LinTerm lin_term) const;

      /// return active Lagrange multiplier vector (normal part)
      inline const Epetra_Vector& GetZnActive() const
      {
        if (not isvalid_.state_) FOUR_C_THROW("zn_active_ is not valid!");

        return *zn_active_;
      }

      /// return inactive Lagrange multiplier vector (normal part)
      inline const Epetra_Vector& GetZnInactive() const
      {
        if (not isvalid_.state_) FOUR_C_THROW("zn_inactive_ is not valid!");

        return *zn_inactive_;
      }

     private:
      /// reset all internal isvalid flags
      void reset_is_valid();

      /// compute all terms for the active linear models
      void compute_lin_active(const Epetra_Vector& dincrSlMa, const Epetra_Vector& znincr_active);

      /// compute all terms for the inactive linear models
      void compute_lin_inactive(const Epetra_Vector& znincr_inactive);

      /// @name Lagrangian potential (contact part)
      /// @{

      /// return the Lagrangian constraint part belonging to \c pot_set
      double get_lagrangian(enum POTENTIAL::SetType pot_set) const;

      /// return contributions for the linear Lagrangian model
      double get_lagrangian_lin(
          enum POTENTIAL::LinTerm lin_term, enum POTENTIAL::SetType pot_set) const;

      /// return active contributions for the linear Lagrangian model
      double get_active_lagrangian_lin(const enum POTENTIAL::LinTerm lin_term) const;

      /// return inactive contributions for the linear Lagrangian model
      double get_inactive_lagrangian_lin(const enum POTENTIAL::LinTerm lin_term) const;

      /// @}

      /// @name Augmented Lagrangian potential (contact part)
      /// @{

      /// return the augmented Lagrangian constraint part belonging to \c pot_set
      double get_augmented_lagrangian(enum POTENTIAL::SetType pot_set) const;

      /// return contributions for the linear augmented Lagrangian model
      double get_augmented_lagrangian_lin(
          enum POTENTIAL::LinTerm lin_term, enum POTENTIAL::SetType pot_set) const;

      /// return active contributions for the linear aug. Lagrangian model
      double get_active_augmented_lagrangian_lin(const enum POTENTIAL::LinTerm lin_term) const;

      /// return inactive contributions for the linear aug. Lagrangian model
      double get_inactive_augmented_lagrangian_lin(const enum POTENTIAL::LinTerm lin_term) const;

      /// @}

      /// @name Infeasibility function (contact)
      /// @{

      /// return the infeasibility function value belonging to \c pot_set
      double get_infeasibility_measure(enum POTENTIAL::SetType pot_set) const;

      /// return contributions for the linear infeasibility model
      double get_infeasibility_measure_lin(
          enum POTENTIAL::LinTerm lin_term, enum POTENTIAL::SetType pot_set) const;

      /// return active contributions for the linear infeasibility model
      double get_active_infeasibility_measure_lin(enum POTENTIAL::LinTerm lin_term) const;

      /// return inactive contributions for the linear infeasibility model
      double get_inactive_infeasibility_measure_lin(enum POTENTIAL::LinTerm lin_term) const;

      /// @}

      /** Split the given global direction vector into its distinct parts and
       *  return these parts */
      void set_direction(const Epetra_Vector& direction, Epetra_Vector& dincrSlMa,
          Epetra_Vector& znincr_active, Epetra_Vector& znincr_inactive);

      /** apply necessary time integration scaling in dynamic simulations to
       *  current contributions, i.e. \f$t_{n+1}\f$ */
      double time_int_scale_np(const double static_part) const;

      /// return the time integration coefficient corresponding to \f$t_{n}\f$
      double get_time_integration_factor() const;

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

      /// isvalid container
      IsValid isvalid_;

      /// has setup been called?
      bool issetup_;

      /// call-back to the wrapping contact strategy
      const CONTACT::Aug::Strategy& strategy_;

      /// call-back to the data container of the surrounding contact strategy
      const CONTACT::Aug::DataContainer& data_;

      /// active Lagrange multiplier vector in normal direction
      Teuchos::RCP<Epetra_Vector> zn_active_;

      /// inactive Lagrange multiplier vector in normal direction
      Teuchos::RCP<Epetra_Vector> zn_inactive_;

      /// active Lagrange multiplier vector in tangential direction
      Teuchos::RCP<Epetra_Vector> zt_active_;

      /// inactive Lagrange multiplier vector in tangential direction
      Teuchos::RCP<Epetra_Vector> zt_inactive_;

      /// container for the potential data
      struct PotData
      {
        /// constructor
        PotData() = default;

        /// print content to screen
        void print(std::ostream& os, const Potential& pot) const;

        /// active part
        double zn_gn_ = 0.0;

        /// active augmented part
        double gn_gn_ = 0.0;

        /// inactive normal part
        double zn_zn_ = 0.0;

        /// tangential part
        double zt_zt_ = 0.0;

        /// infeasibility active square part
        double inf_gn_gn_ = 0.0;

        /// infeasibility inactive square part
        double inf_zn_zn_ = 0.0;
      };

      /// container for the linearization data
      struct LinData
      {
        /// constructor
        LinData() = default;

        /// reset all container variables
        void reset()
        {
          reset_active();
          reset_inactive();
        }

        /// reset all active variables
        void reset_active()
        {
          gn_dzn_ = 0.0;
          zn_dgn_ = 0.0;
          gn_dgn_ = 0.0;
          dzn_dgn_ = 0.0;
          inf_gn_dgn_ = 0.0;
        }

        /// reset all inactive variables
        void reset_inactive()
        {
          zn_dzn_ = 0.0;
          dzn_dzn_ = 0.0;
          inf_zn_dzn_ = 0.0;
        }

        /// print content to screen
        void print(std::ostream& os, const Potential& pot) const;

        /// linearization of the active part w.r.t. zn
        double gn_dzn_ = 0.0;

        /// linearization of the active part w.r.t. d
        double zn_dgn_ = 0.0;

        /// linearization of the augmented part w.r.t. d
        double gn_dgn_ = 0.0;

        /// linearization of the inactive part w.r.t. zn
        double zn_dzn_ = 0.0;

        /// linearization of the active part w.r.t. d and zn
        double dzn_dgn_ = 0.0;

        /// 2-nd derivative of the inactive part w.r.t. zn
        double dzn_dzn_ = 0.0;

        /// linearization of the infeasibility active square part
        double inf_gn_dgn_ = 0.0;

        /// linearization of the infeasibility inactive square part
        double inf_zn_dzn_ = 0.0;
      };

      /// instance of the potential data container
      PotData potdata_;

      /// instance of the linear model data container
      LinData lindata_;
    };
  }  // namespace Aug
}  // namespace CONTACT


FOUR_C_NAMESPACE_CLOSE

#endif
