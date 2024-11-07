// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_CARDIOVASCULAR0D_RESPIRATORY_SYSPULPERIPHCIRCULATION_HPP
#define FOUR_C_CARDIOVASCULAR0D_RESPIRATORY_SYSPULPERIPHCIRCULATION_HPP

#include "4C_config.hpp"

#include "4C_cardiovascular0d.hpp"
#include "4C_fem_general_utils_integration.hpp"
#include "4C_inpar_cardiovascular0d.hpp"
#include "4C_linalg_vector.hpp"
#include "4C_utils_parameter_list.fwd.hpp"

#include <Epetra_FECrsMatrix.h>
#include <Epetra_Operator.h>
#include <Epetra_RowMatrix.h>

#include <memory>

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

namespace Utils
{
  class CardiovascularRespiratory0DSysPulPeriphCirculation : public Cardiovascular0D

  {
   public:
    /*!
    \brief Constructor of a Cardiovascular0D based on conditions with a given name. It also
    takes care of the Cardiovascular0D IDs.
    */

    CardiovascularRespiratory0DSysPulPeriphCirculation(
        std::shared_ptr<Core::FE::Discretization>
            discr,                         ///< discretization where Cardiovascular0D lives on
        const std::string& conditionname,  ///< Name of condition to create Cardiovascular0D from
        std::vector<int>& curID            ///< current ID
    );



    /// initialization routine called by the manager ctor to get correct reference base values and
    /// activating the right conditions at the beginning
    void initialize(
        Teuchos::ParameterList&
            params,  ///< parameter list to communicate between elements and discretization
        std::shared_ptr<Core::LinAlg::Vector<double>>
            sysvec1,  ///< distributed vector that may be filled by
                      ///< assembly of element contributions
        std::shared_ptr<Core::LinAlg::Vector<double>>
            sysvec2  ///< distributed vector that may be filled by assembly of element contributions
        ) override;

    //! Evaluate routine to call from outside. In here the right action is determined and the
    //! #EvaluateCardiovascular0D routine is called
    void evaluate(
        Teuchos::ParameterList&
            params,  ///< parameter list to communicate between elements and discretization
        std::shared_ptr<Core::LinAlg::SparseMatrix> sysmat1,  ///< Cardiovascular0D stiffness matrix
        std::shared_ptr<Core::LinAlg::SparseOperator>
            sysmat2,  ///< Cardiovascular0D offdiagonal matrix dV/dd
        std::shared_ptr<Core::LinAlg::SparseOperator>
            sysmat3,  ///< Cardiovascular0D offdiagonal matrix dfext/dp
        std::shared_ptr<Core::LinAlg::Vector<double>>
            sysvec1,  ///< distributed vectors that may be filled by
                      ///< assembly of element contributions
        std::shared_ptr<Core::LinAlg::Vector<double>> sysvec2,
        std::shared_ptr<Core::LinAlg::Vector<double>> sysvec3,
        const std::shared_ptr<Core::LinAlg::Vector<double>> sysvec4,
        std::shared_ptr<Core::LinAlg::Vector<double>> sysvec5) override;

    // cbO2 and its derivatives
    double cb_o2(double ppCO2, double ppO2);
    // w.r.t. O2
    double dcb_o2_dpp_o2(double ppCO2, double ppO2);
    double d2cb_o2_dpp_o22(double ppCO2, double ppO2);
    // w.r.t. CO2
    double dcb_o2_dpp_c_o2(double ppCO2, double ppO2);
    double d2cb_o2_dpp_c_o22(double ppCO2, double ppO2);
    double d2cb_o2_dpp_o2dpp_c_o2(double ppCO2, double ppO2);

    // cbCO2 and its derivatives
    double cb_c_o2(double ppCO2, double ppO2);
    // w.r.t. CO2
    double dcb_c_o2_dpp_c_o2(double ppCO2, double ppO2);
    double d2cb_c_o2_dpp_c_o22(double ppCO2, double ppO2);
    // w.r.t. O2
    double dcb_c_o2_dpp_o2(double ppCO2, double ppO2);
    double d2cb_c_o2_dpp_o22(double ppCO2, double ppO2);
    double d2cb_c_o2_dpp_c_o2dpp_o2(double ppCO2, double ppO2);

    double ct_o2(double ppO2);
    double dct_o2_dpp_o2(double ppO2);
    double d2ct_o2_dpp_o22(double ppO2);

    double s_o2(double ppCO2, double ppO2);

    double ct_c_o2(double ppCO2);
    double dct_c_o2_dpp_c_o2(double ppCO2);
    double d2ct_c_o2_dpp_c_o22(double ppCO2);


    //! Evaluate routine to call from outside. In here the right action is determined and the
    //! #EvaluateCardiovascular0D routine is called
    virtual void evaluate_respiratory(Teuchos::ParameterList& params, std::vector<double>& df_np,
        std::vector<double>& f_np, Core::LinAlg::SerialDenseMatrix& wkstiff,
        std::shared_ptr<Core::LinAlg::Vector<double>> dofvec,
        std::shared_ptr<Core::LinAlg::Vector<double>> volvec, bool evalstiff);

   private:
    // number of degrees of freedom for submodels
    int num_dof_cardio_;
    int num_dof_respir_;

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
    // peripheral arterial compliances and resistances
    double c_arspl_sys_;   ///< systemic arterial splanchnic compliance
    double r_arspl_sys_;   ///< systemic arterial splanchnic resistance
    double c_arespl_sys_;  ///< systemic arterial extra-splanchnic compliance
    double r_arespl_sys_;  ///< systemic arterial extra-splanchnic resistance
    double c_armsc_sys_;   ///< systemic arterial muscular compliance
    double r_armsc_sys_;   ///< systemic arterial muscular resistance
    double c_arcer_sys_;   ///< systemic arterial cerebral compliance
    double r_arcer_sys_;   ///< systemic arterial cerebral resistance
    double c_arcor_sys_;   ///< systemic arterial coronary compliance
    double r_arcor_sys_;   ///< systemic arterial coronary resistance
    // peripheral venous compliances and resistances
    double c_venspl_sys_;   ///< systemic venous splanchnic compliance
    double r_venspl_sys_;   ///< systemic venous splanchnic resistance
    double c_venespl_sys_;  ///< systemic venous extra-splanchnic compliance
    double r_venespl_sys_;  ///< systemic venous extra-splanchnic resistance
    double c_venmsc_sys_;   ///< systemic venous muscular compliance
    double r_venmsc_sys_;   ///< systemic venous muscular resistance
    double c_vencer_sys_;   ///< systemic venous cerebral compliance
    double r_vencer_sys_;   ///< systemic venous cerebral resistance
    double c_vencor_sys_;   ///< systemic venous coronary compliance
    double r_vencor_sys_;   ///< systemic venous coronary resistance

    double c_ar_pul_;  ///< pulmonary arterial compliance
    double r_ar_pul_;  ///< pulmonary arterial resistance
    double l_ar_pul_;  ///< pulmonary arterial inertance
    double z_ar_pul_;  ///< pulmonary arterial impedance

    // pulmonary capillary compliance and resistance
    double c_cap_pul_;  ///< pulmonary capillary compliance
    double r_cap_pul_;  ///< pulmonary capillary resistance

    double c_ven_sys_;  ///< systemic venous compliance
    double r_ven_sys_;  ///< systemic venous resistance
    double l_ven_sys_;  ///< systemic venous inertance
    double c_ven_pul_;  ///< pulmonary venous compliance
    double r_ven_pul_;  ///< pulmonary venous resistance
    double l_ven_pul_;  ///< pulmonary venous inertance

    // unstressed 0D volumes
    double v_v_l_u_;     ///< dead left ventricular volume (at zero pressure) (for 0D ventricles)
    double v_at_l_u_;    ///< dead left atrial volume (at zero pressure) (for 0D atria)
    double v_ar_sys_u_;  ///< unstressed systemic arterial volume

    double v_arspl_sys_u_;    ///< unstressed systemic arterial splanchnic volume
    double v_arespl_sys_u_;   ///< unstressed systemic arterial extra-splanchnic volume
    double v_armsc_sys_u_;    ///< unstressed systemic arterial muscular volume
    double v_arcer_sys_u_;    ///< unstressed systemic arterial cerebral volume
    double v_arcor_sys_u_;    ///< unstressed systemic arterial coronary volume
    double v_venspl_sys_u_;   ///< unstressed systemic venous splanchnic volume
    double v_venespl_sys_u_;  ///< unstressed systemic venous extra-splanchnic volume
    double v_venmsc_sys_u_;   ///< unstressed systemic venous muscular volume
    double v_vencer_sys_u_;   ///< unstressed systemic venous cerebral volume
    double v_vencor_sys_u_;   ///< unstressed systemic venous coronary volume

    double v_ven_sys_u_;  ///< unstressed systemic venous volume
    double v_v_r_u_;      ///< unstressed  right ventricular volume (for 0D ventricles)
    double v_at_r_u_;     ///< unstressed  right atrial volume (for 0D atria)
    double v_ar_pul_u_;   ///< unstressed pulmonary arterial volume
    double v_cap_pul_u_;  ///< unstressed pulmonary capillary volume
    double v_ven_pul_u_;  ///< unstressed pulmonary venous volume

    // total number of degrees of freedom
    int num_dof_;

    // 0D lung
    double l_alv_;   ///< alveolar inertance
    double r_alv_;   ///< alveolar resistance
    double e_alv_;   ///< alveolar elastance
    int u_t_curve_;  ///< time-varying, prescribed pleural pressure curve driven by diaphragm
    double u_m_;     ///< in-breath pressure

    double v_lung_tidal_;  ///< tidal volume (the total volume of inspired air, in a single breath)
    double v_lung_dead_;   ///< dead space volume
    double v_lung_u_;  ///< unstressed lung volume (volume of the lung when it is fully collapsed
                       ///< outside the body)

    double f_c_o2_ext_;  ///< atmospheric CO2 gas fraction
    double f_o2_ext_;    ///< atmospheric O2 gas fraction

    // should be 22.4 liters per mol !
    // however we specify it as an input parameter since its decimal power depends on the system of
    // units your whole model is specified in! i.e. if you have kg - mm - s - mmol, it's 22.4e3 mm^3
    // / mmol
    double v_m_gas_;  ///< molar volume of an ideal gas

    // should be 47.1 mmHg = 6.279485 kPa !
    // however we specify it as an input parameter since its decimal power depends on the system of
    // units your whole model is specified in! i.e. if you have kg - mm - s - mmol, it's 6.279485
    // kPa
    double p_vap_water_37_;  ///< vapor pressure of water at 37  degrees celsius

    double kappa_c_o2_;  ///< diffusion coefficient for CO2 across the hemato-alveolar membrane, in
                         ///< molar value / (time * pressure)
    double kappa_o2_;    ///< diffusion coefficient for O2 across the hemato-alveolar membrane, in
                         ///< molar value / (time * pressure)

    double alpha_c_o2_;  ///< CO2 solubility constant, in molar value / (volume * pressure)
    double alpha_o2_;    ///< O2 solubility constant, in molar value / (volume * pressure)

    double c_hb_;  ///< hemoglobin concentration of the blood, in molar value / volume

    double m_c_o2_arspl_;   ///< splanchnic metabolic rate of CO2 production
    double m_o2_arspl_;     ///< splanchnic metabolic rate of O2 consumption
    double m_c_o2_arespl_;  ///< extra-splanchnic metabolic rate of CO2 production
    double m_o2_arespl_;    ///< extra-splanchnic metabolic rate of O2 consumption
    double m_c_o2_armsc_;   ///< muscular metabolic rate of CO2 production
    double m_o2_armsc_;     ///< muscular metabolic rate of O2 consumption
    double m_c_o2_arcer_;   ///< cerebral metabolic rate of CO2 production
    double m_o2_arcer_;     ///< cerebral metabolic rate of O2 consumption
    double m_c_o2_arcor_;   ///< coronary metabolic rate of CO2 production
    double m_o2_arcor_;     ///< coronary metabolic rate of O2 consumption

    // tissue columes
    double v_tissspl_;   ///< splanchnic tissue volume
    double v_tissespl_;  ///< extra-splanchnic tissue volume
    double v_tissmsc_;   ///< muscular tissue volume
    double v_tisscer_;   ///< cerebral tissue volume
    double v_tisscor_;   ///< coronary tissue volume



    // don't want = operator, cctor and destructor

    CardiovascularRespiratory0DSysPulPeriphCirculation operator=(
        const CardiovascularRespiratory0DSysPulPeriphCirculation& old);
    CardiovascularRespiratory0DSysPulPeriphCirculation(
        const CardiovascularRespiratory0DSysPulPeriphCirculation& old);



  };  // class
}  // namespace Utils

FOUR_C_NAMESPACE_CLOSE

#endif
