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

namespace CONTACT {
namespace AUG {
  class Strategy;
  class DataContainer;

  /** \brief Evaluate the contact potential
   *
   *  \author hiermeier \date 03/17 */
  class Potential
  {
  public:
    enum PotentialType
    {
      type_lagrangian,             ///< Lagrangian potential
      type_augmented_lagrangian,   ///< Augmented lagrangian potential
      type_infeasibility_measure   ///< infeasibility measure
    };

    enum PotentialTerm
    {
      term_active,               ///< active potential terms
      term_augmented,            ///< augmented potential terms
      term_inactive,             ///< inactive potential terms
      term_all                   ///< entire contact potential
    };

    enum LinearizationTerm
    {
      lin_active_wrt_d,         ///< linearization of the active term w.r.t. the displacements
      lin_active_wrt_z,         ///< linearization of the active term w.r.t. the Lagrange multipliers
      lin_augmented_wrt_d,      ///< linearization of the augmented term w.r.t. the displacements
      lin_inactive_wrt_z,       ///< linearization of the inactive term w.r.t. the Lagrange multipliers
      lin_active_wrt_d_and_z,   ///< linearization of the active term w.r.t. the displacements and the Lagrange multipliers
      lin_inactive_wrt_z_and_z  ///< 2-nd derivative of the inactive term w.r.t. the Lagrange multipliers
    };

  public:
    /// constructor
    Potential( const CONTACT::AUG::Strategy& strategy,
        const CONTACT::AUG::DataContainer& data );

    /// destructor
    virtual ~Potential() { };

    void Setup( bool active_inactive_only = false );

    void ResetIsValid();

    void SetActiveInactiveState();

    void SetDirection( const Epetra_Vector& direction );

    /// compute the contact potential terms
    void Compute();

    /// compute the linearizations terms of the contact potential terms
    void ComputeLin();

    double Get( enum PotentialType pot_type, enum PotentialTerm pot_term ) const;

    double GetLin(  enum PotentialType pot_type, enum LinearizationTerm lin_term ) const;

    inline const Epetra_Vector& GetZnIncrActive() const
    {
      if ( not isvalidDirection_ )
        dserror( "znincr_active_ is not valid!" );

      return *znincr_active_;
    }

    inline const Epetra_Vector& GetZnIncrInactive() const
    {
      if ( not isvalidDirection_ )
        dserror( "znincr_inactive_ is not valid!" );

      return *znincr_inactive_;
    }

    inline const Epetra_Vector& GetZnActive() const
    {
      if ( not isvalidState_ )
        dserror( "zn_active_ is not valid!" );

      return *zn_active_;
    }

    inline const Epetra_Vector& GetZnInactive() const
    {
      if ( not isvalidState_ )
        dserror( "zn_inactive_ is not valid!" );

      return *zn_inactive_;
    }

  private:

    double GetLagrangian( enum PotentialTerm pot_term ) const;
    double GetLagrangianLin( enum LinearizationTerm lin_term ) const;

    double GetAugmentedLagrangian( enum PotentialTerm pot_term ) const;
    double GetAugmentedLagrangianLin( enum LinearizationTerm lin_term ) const;

    double GetInfeasibilityMeasure( enum PotentialTerm pot_term ) const;
    double GetInfeasibilityMeasureLin( enum LinearizationTerm lin_term ) const;

  private:
    bool isvalidPotential_;
    bool isvalidLinearization_;
    bool isvalidState_;
    bool isvalidDirection_;
    bool issetup_;

    const CONTACT::AUG::Strategy& strategy_;
    const CONTACT::AUG::DataContainer& data_;

    Teuchos::RCP<Epetra_Vector> zn_active_;
    Teuchos::RCP<Epetra_Vector> zn_inactive_;

    Teuchos::RCP<Epetra_Vector> zt_active_;
    Teuchos::RCP<Epetra_Vector> zt_inactive_;

    Teuchos::RCP<Epetra_Vector> dincrSlMa_;

    Teuchos::RCP<Epetra_Vector> znincr_active_;
    Teuchos::RCP<Epetra_Vector> znincr_inactive_;

    /// container for the potential data
    struct PotData
    {
      PotData() : zn_gn_( 0.0 ),
              gn_gn_( 0.0 ),
              zn_zn_( 0.0 ),
              zt_zt_( 0.0 ),
              inf_gn_gn_( 0.0 ),
              inf_zn_zn_( 0.0 )
      {

      }

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
      LinData() : gn_dzn_( 0.0 ),
              zn_dgn_( 0.0 ),
              gn_dgn_( 0.0 ),
              zn_dzn_( 0.0 ),
              dzn_dgn_( 0.0 ),
              dzn_dzn_( 0.0 ),
              inf_gn_dgn_( 0.0 ),
              inf_zn_dzn_( 0.0 )
      {

      }

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
}
}


#endif /* SRC_DRT_CONTACT_AUG_CONTACT_AUG_POTENTIAL_H_ */
