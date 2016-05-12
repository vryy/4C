/*-----------------------------------------------------------*/
/*!
\file drt_elements_paramsinterface.H

\brief Declaration
\level 3
\maintainer Michael Hiermeier

\date Mar 24, 2016

*/
/*-----------------------------------------------------------*/

#ifndef SRC_DRT_LIB_DRT_ELEMENTS_PARAMSINTERFACE_H_
#define SRC_DRT_LIB_DRT_ELEMENTS_PARAMSINTERFACE_H_

#include <Teuchos_RCP.hpp>

namespace DRT {
namespace ELEMENTS {
  /*! action which the element has to perform
   *
   * ToDo Please add ALL your action types here! */
  enum ActionType
  {
    none,
    struct_calc_linstiff,
    struct_calc_nlnstiff,                   //!< evaluate the tangential stiffness matrix and the internal force vector
    struct_calc_internalforce,              //!< evaluate only the internal forces (no need for the stiffness terms)
    struct_calc_linstiffmass,
    struct_calc_nlnstiffmass,               //!< evaluate the dynamic state: internal forces vector, stiffness and the default/nln mass matrix
    struct_calc_nlnstifflmass,              //!< evaluate the dynamic state: internal forces vector, stiffness and the lumped mass matrix
    struct_calc_nlnstiff_gemm,              //!< internal force, stiffness and mass for GEMM
    struct_calc_recover,                    //!< recover elementwise condensed internal variables
    struct_calc_stress,
    struct_calc_eleload,
    struct_calc_fsiload,
    struct_calc_update_istep,
    struct_calc_reset_istep,                //!< reset elementwise internal variables, during iteration to last converged state
    struct_calc_store_istep,                //!< store internal information in history
    struct_calc_recover_istep,              //!< recover internal information from history
    struct_calc_reset_all,                  //!< reset elementwise internal variables to state in the beginning of the computation
    struct_calc_energy,                     //!< compute internal energy
    struct_calc_errornorms,                 //!< compute error norms (L2,H1,energy)
    struct_postprocess_stress,              //!< postprocess stresses
    struct_update_prestress,
    inversedesign_update,
    inversedesign_switch,
    multi_readrestart,                      //!< multi-scale: read restart on microscale
    multi_init_eas,                         //!< multi-scale: initialize EAS parameters on microscale
    multi_set_eas,                          //!< multi-scale: set EAS parameters on microscale
    multi_calc_dens,                        //!< multi-scale: calculate homogenized density
    calc_potential_stiff,                   //!< volume interaction potential
    shell_calc_stc_matrix,                  //!< calculate scaled director matrix for thin shell structures
    shell_calc_stc_matrix_inverse,          //!< calculate inverse of scaled director matrix for thin shell structures
    struct_calc_ptcstiff,                   //!< calculate artificial stiffness due to PTC solution strategy
    struct_calc_stifftemp,                  //!< TSI specific: mechanical-thermal stiffness
    struct_calc_global_gpstresses_map,      //!< basically calc_struct_stress but with assembly of global gpstresses map
    struct_interpolate_velocity_to_point    //!< interpolate the structural velocity to a given point
  };  // enum ActionType

  //! Map action type enum to std::string
  static inline const std::string ActionType2String(const enum ActionType& type)
  {
    switch (type)
    {
      case none:
        return "none";
        break;
      case struct_calc_linstiff:
        return "struct_calc_linstiff";
        break;
      case struct_calc_nlnstiff:
        return "struct_calc_nlnstiff";
        break;
      case struct_calc_internalforce:
        return "struct_calc_internalforce";
        break;
      case struct_calc_linstiffmass:
        return "struct_calc_linstiffmass";
        break;
      case struct_calc_nlnstiffmass:
        return "struct_calc_nlnstiffmass";
        break;
      case struct_calc_nlnstifflmass:
        return "struct_calc_nlnstifflmass";
        break;
      case struct_calc_nlnstiff_gemm:
        return "struct_calc_nlnstiff_gemm";
        break;
      case struct_calc_recover:
        return "struct_calc_recover";
        break;
      case struct_calc_stress:
        return "struct_calc_stress";
        break;
      case struct_calc_eleload:
        return "struct_calc_eleload";
        break;
      case struct_calc_fsiload:
        return "struct_calc_fsiload";
        break;
      case struct_calc_update_istep:
        return "struct_calc_update_istep";
        break;
      case struct_calc_reset_istep:
        return "struct_calc_reset_istep";
        break;
      case struct_calc_store_istep:
        return "struct_calc_store_istep";
        break;
      case struct_calc_recover_istep:
        return "struct_calc_recover_istep";
        break;
      case struct_calc_reset_all:
        return "struct_calc_reset_all";
        break;
      case struct_calc_energy:
        return "struct_calc_energy";
        break;
      case struct_calc_errornorms:
        return "struct_calc_errornorms";
        break;
      case struct_postprocess_stress:
        return "struct_postprocess_stress";
        break;
      case struct_update_prestress:
        return "struct_update_prestress";
        break;
      case inversedesign_update:
        return "inversedesign_update";
        break;
      case inversedesign_switch:
        return "inversedesign_switch";
        break;
      case multi_readrestart:
        return "multi_readrestart";
        break;
      case multi_init_eas:
        return "multi_init_eas";
        break;
      case multi_set_eas:
        return "multi_set_eas";
        break;
      case multi_calc_dens:
        return "multi_calc_dens";
        break;
      case calc_potential_stiff:
        return "calc_potential_stiff";
        break;
      case shell_calc_stc_matrix:
        return "shell_calc_stc_matrix";
        break;
      case shell_calc_stc_matrix_inverse:
        return "shell_calc_stc_matrix_inverse";
        break;
        case struct_calc_ptcstiff:
        return "struct_calc_ptcstiff";
        break;
      case struct_calc_stifftemp:
        return "struct_calc_stifftemp";
        break;
      case struct_calc_global_gpstresses_map:
        return "struct_calc_global_gpstresses_map";
        break;
      case struct_interpolate_velocity_to_point:
        return "struct_interpolate_velocity_to_point";
        break;
      default:
        return "unknown";
        break;
    }
    return "";
  }; // ActionType2String

  /*! \brief Declaration Parameter interface for the element <--> time integrator data exchange
   *
   *  Pure virtual interface class. This class is supposed to replace the current
   *  tasks of the Teuchos::ParameterList.
   *  Please consider to derive a special interface class, if you need special parameters inside
   *  of your element. Keep the Evaluate call untouched and cast the interface object to the
   *  desired specification when and where you need it.
   *
   *  ToDo Currently we set the interface in the elements via the Teuchos::ParameterList.
   *  Theoretically, the Teuchos::ParameterList can be replaced by the interface itself!
   *
   *  \date 03/2016
   *  \author hiermeier */
  class ParamsInterface
  {
  public:
    //! destructor
    virtual ~ParamsInterface() {};

    //! @name Access general control parameters
    //! @{
    //! get the desired action type
    virtual enum ActionType GetActionType() const = 0;

    //! get the current total time for the evaluate call
    virtual double GetTotalTime() const = 0;

    //! get the current time step
    virtual double GetDeltaTime() const = 0;
    //! @}
  };  // class ParamsInterface
} // namespace ELEMENTS
} // namespace DRT


#endif /* SRC_DRT_LIB_DRT_ELEMENTS_PARAMSINTERFACE_H_ */
