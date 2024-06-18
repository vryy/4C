/*----------------------------------------------------------------------*/
/*! \file

\brief Monolithic coupling of 3D structure Cardiovascular0D models

\level 2


A closed-loop cardiovascular model with 0D bi-resistive cardiac valve laws and lumped parameter
models for the systemic and pulmonary circulation; heart can either be fully 3D, partly 3D (only
ventricles, atria 0D elastance models) or fully 0D (elastance models for ventricles and atria)
specification: DESIGN SURF CARDIOVASCULAR 0D SYS-PUL CIRCULATION CONDITIONS

detailed fully in Hirschvogel, Bassilious, Jagschies, Wildhirt, Gee, "A monolithic 3D-0D coupled
closed-loop model of the heart and the vascular system: Experiment-based parameter estimation for
patient-specific cardiac mechanics", IJNMBE, 2016

The residual vector reads:

      [d(p_at_l/E_at_l)/dt - q_ven_pul + q_vin_l  OR*  d(V_at_l)/dt - q_ven_pul + q_vin_l  ]   [ 0 ]
      [(p_at_l - p_v_l)/R_atv_l - q_vin_l                                                  ]   [ 0 ]
      [d(V_v_l)/dt - q_vin_l + q_vout_l  OR*  d(p_v_l/E_v_l)/dt - q_vin_l + q_vout_l       ]   [ 0 ]
      [(p_v_l - p_ar_sys)/R_arv_l - q_vout_l                                               ]   [ 0 ]
      [C_ar_sys * (d(p_ar_sys)/dt - Z_ar_sys * d(q_vout_l)/dt) - q_vout_l + q_ar_sys       ]   [ 0 ]
      [L_ar_sys/R_ar_sys + (p_ven_sys - p_ar_sys + Z_ar_sys * q_vout_l)/R_ar_sys + q_ar_sys]   [ 0 ]
      [C_ven_sys * d(p_ven_sys)/dt - q_ar_sys + q_ven_sys                                  ]   [ 0 ]
      [L_ven_sys/R_ven_sys + (p_at_r - p_ven_sys)/R_ven_sys + q_ven_sys                    ]   [ 0 ]
Res =                                                                                        =
      [d(p_at_r/E_at_r)/dt - q_ven_sys + q_vin_r  OR*  d(V_at_r)/dt - q_ven_sys + q_vin_r  ]   [ 0 ]
      [(p_at_r - p_v_r)/R_atv_r - q_vin_r                                                  ]   [ 0 ]
      [d(V_v_r)/dt - q_vin_r + q_vout_r  OR*  d(p_v_r/E_v_r)/dt - q_vin_r + q_vout_r       ]   [ 0 ]
      [(p_v_r - p_ar_pul)/R_arv_r - q_vout_r                                               ]   [ 0 ]
      [C_ar_pul * (d(p_ar_pul)/dt - Z_ar_pul * d(q_vout_r)/dt) - q_vout_r + q_ar_pul       ]   [ 0 ]
      [L_ar_pul/R_ar_pul + (p_ven_pul - p_ar_pul + Z_ar_pul * q_vout_r)/R_ar_pul + q_ar_pul]   [ 0 ]
      [C_ven_pul * d(p_ven_pul)/dt - q_ar_pul + q_ven_pul                                  ]   [ 0 ]
      [L_ven_pul/R_ven_pul + (p_at_l - p_ven_pul)/R_ven_pul + q_ven_pul                    ]   [ 0 ]

* depending on atrial/ventricular model (0D elastance vs. 3D structural)

*----------------------------------------------------------------------*/

#ifndef FOUR_C_CARDIOVASCULAR0D_SYSPULCIRCULATION_HPP
#define FOUR_C_CARDIOVASCULAR0D_SYSPULCIRCULATION_HPP

#include "4C_config.hpp"

#include "4C_cardiovascular0d.hpp"
#include "4C_fem_general_utils_integration.hpp"
#include "4C_inpar_cardiovascular0d.hpp"

#include <Epetra_FECrsMatrix.h>
#include <Epetra_Operator.h>
#include <Epetra_RowMatrix.h>
#include <Epetra_Vector.h>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_RCP.hpp>

FOUR_C_NAMESPACE_OPEN

// forward declarations
namespace Core::FE
{
  class Discretization;
}  // namespace Core::FE

namespace Core::LinAlg
{
  class SparseMatrix;
  class SparseOperator;
}  // namespace Core::LinAlg

namespace UTILS
{
  class Cardiovascular0DSysPulCirculation : public Cardiovascular0D

  {
   public:
    /*!
    \brief Constructor of a Cardiovascular0D based on conditions with a given name. It also
    takes care of the Cardiovascular0D IDs.
    */

    Cardiovascular0DSysPulCirculation(
        Teuchos::RCP<Core::FE::Discretization>
            discr,                         ///< discretization where Cardiovascular0D lives on
        const std::string& conditionname,  ///< Name of condition to create Cardiovascular0D from
        std::vector<int>& curID            ///< current ID
    );



    /// initialization routine called by the manager ctor to get correct reference base values and
    /// activating the right conditions at the beginning
    void Initialize(
        Teuchos::ParameterList&
            params,  ///< parameter list to communicate between elements and discretization
        Teuchos::RCP<Epetra_Vector> sysvec1,  ///< distributed vector that may be filled by assembly
                                              ///< of element contributions
        Teuchos::RCP<Epetra_Vector>
            sysvec2  ///< distributed vector that may be filled by assembly of element contributions
        ) override;

    //! Evaluate routine to call from outside. In here the right action is determined and the
    //! #EvaluateCardiovascular0D routine is called
    void evaluate(
        Teuchos::ParameterList&
            params,  ///< parameter list to communicate between elements and discretization
        Teuchos::RCP<Core::LinAlg::SparseMatrix> sysmat1,  ///< Cardiovascular0D stiffness matrix
        Teuchos::RCP<Core::LinAlg::SparseOperator>
            sysmat2,  ///< Cardiovascular0D offdiagonal matrix dV/dd
        Teuchos::RCP<Core::LinAlg::SparseOperator>
            sysmat3,                          ///< Cardiovascular0D offdiagonal matrix dfext/dp
        Teuchos::RCP<Epetra_Vector> sysvec1,  ///< distributed vectors that may be filled by
                                              ///< assembly of element contributions
        Teuchos::RCP<Epetra_Vector> sysvec2, Teuchos::RCP<Epetra_Vector> sysvec3,
        const Teuchos::RCP<Epetra_Vector> sysvec4, Teuchos::RCP<Epetra_Vector> sysvec5) override;

   private:
    // parameters
    double r_arvalve_max_l_;       ///< maximum aortic valve resistance
    double r_arvalve_min_l_;       ///< minimum aortic valve resistance
    double r_atvalve_max_l_;       ///< maximum mitral valve resistance
    double r_atvalve_min_l_;       ///< minimum mitral valve resistance
    double r_arvalve_max_r_;       ///< maximum pulmonary valve resistance
    double r_arvalve_min_r_;       ///< minimum pulmonary valve resistance
    double r_atvalve_max_r_;       ///< maximum tricuspid valve resistance
    double r_atvalve_min_r_;       ///< minimum tricuspid valve resistance
    int atrium_act_curve_l_;       ///< left atrial activation curve (ONLY for ATRIUM_MODEL "0D"!)
    int atrium_act_curve_r_;       ///< right atrial activation curve (ONLY for ATRIUM_MODEL "0D"!)
    int ventricle_act_curve_l_;    ///< left ventricular activation curve (ONLY for VENTRICLE_MODEL
                                   ///< "0D"!)
    int ventricle_act_curve_r_;    ///< right ventricular activation curve (ONLY for VENTRICLE_MODEL
                                   ///< "0D"!)
    int atrium_prescr_e_curve_l_;  ///< left atrial elastance prescription curve (ONLY for
                                   ///< ATRIUM_MODEL "prescribed"!)
    int atrium_prescr_e_curve_r_;  ///< right atrial elastance prescription curve (ONLY for
                                   ///< ATRIUM_MODEL "prescribed"!)
    int ventricle_prescr_e_curve_l_;  ///< left ventricular elastance prescription curve (ONLY for
                                      ///< VENTRICLE_MODEL "prescribed"!)
    int ventricle_prescr_e_curve_r_;  ///< right ventricular elastance prescription curve (ONLY for
                                      ///< VENTRICLE_MODEL "prescribed"!)
    double e_at_max_l_;               ///< maximum left atrial elastance (for 0D atria)
    double e_at_min_l_;               ///< minimum left atrial elastance (for 0D atria)
    double e_at_max_r_;               ///< maximum right atrial elastance (for 0D atria)
    double e_at_min_r_;               ///< minimum right atrial elastance (for 0D atria)
    double e_v_max_l_;                ///< maximum left ventricular elastance (for 0D ventricles)
    double e_v_min_l_;                ///< minimum left ventricular elastance (for 0D ventricles)
    double e_v_max_r_;                ///< maximum right ventricular elastance (for 0D ventricles)
    double e_v_min_r_;                ///< minimum right ventricular elastance (for 0D ventricles)
    double c_ar_sys_;                 ///< systemic arterial compliance
    double r_ar_sys_;                 ///< systemic arterial resistance
    double l_ar_sys_;                 ///< systemic arterial inertance
    double z_ar_sys_;                 ///< systemic arterial impedance
    double c_ar_pul_;                 ///< pulmonary arterial compliance
    double r_ar_pul_;                 ///< pulmonary arterial resistance
    double l_ar_pul_;                 ///< pulmonary arterial inertance
    double z_ar_pul_;                 ///< pulmonary arterial impedance
    double c_ven_sys_;                ///< systemic venous compliance
    double r_ven_sys_;                ///< systemic venous resistance
    double l_ven_sys_;                ///< systemic venous inertance
    double c_ven_pul_;                ///< pulmonary venous compliance
    double r_ven_pul_;                ///< pulmonary venous resistance
    double l_ven_pul_;                ///< pulmonary venous inertance

    // unstressed 0D volumes - for post-processing only, since all 0D equations are formulated in
    // terms of fluxes!
    double v_v_l_u_;      ///< unstressed left ventricular volume (for 0D ventricles)
    double v_at_l_u_;     ///< unstressed left atrial volume (for 0D atria)
    double v_ar_sys_u_;   ///< unstressed systemic arterial volume
    double v_ven_sys_u_;  ///< unstressed systemic venous volume
    double v_v_r_u_;      ///< dead right ventricular volume (at zero pressure) (for 0D ventricles)
    double v_at_r_u_;     ///< dead right atrial volume (at zero pressure) (for 0D atria)
    double v_ar_pul_u_;   ///< unstressed pulmonary arterial volume
    double v_ven_pul_u_;  ///< unstressed pulmonary venous volume

    // don't want = operator, cctor and destructor

    Cardiovascular0DSysPulCirculation operator=(const Cardiovascular0DSysPulCirculation& old);
    Cardiovascular0DSysPulCirculation(const Cardiovascular0DSysPulCirculation& old);



  };  // class
}  // namespace UTILS

FOUR_C_NAMESPACE_CLOSE

#endif
