/*!----------------------------------------------------------------------

\brief Monolithic coupling of 3D structure Cardiovascular0D models

\level 3

\maintainer Amadeus Gebauer

A closed-loop cardiovascular model with 0D bi-resistive cardiac valve laws and lumped parameter
models for the systemic and pulmonary circulation; heart can either be fully 3D, partly 3D (only
ventricles, atria 0D elastance models) or fully 0D (elastance models for ventricles and atria)
specification: DESIGN SURF CARDIOVASCULAR RESPIRATORY 0D SYS-PUL CIRCULATION PERIPH CONDITIONS

Extension of the model in Hirschvogel, Bassilious, Jagschies, Wildhirt, Gee, "A monolithic 3D-0D
coupled closed-loop model of the heart and the vascular system: Experiment-based parameter
estimation for patient-specific cardiac mechanics", IJNMBE, 2016, extension according to Trenhago et
al., "An integrated mathematical model of the cardiovascular and respiratory systems", IJNMBE, 2016,
and Ursino and Magosso, "Acute cardiovascular response to isocapnic hypoxia. I. A mathematical
model", Am J Physiol Heart Circ Physiol, 279: H149–H165, 2000

The residual vector reads:

      [d(p_at_l/E_at_l)/dt - q_ven_pul + q_vin_l  OR*  d(V_at_l)/dt - q_ven_pul + q_vin_l ]   [ 0 ]
      [(p_at_l - p_v_l)/R_atv_l - q_vin_l ]   [ 0 ] [d(V_v_l)/dt - q_vin_l + q_vout_l  OR*
d(p_v_l/E_v_l)/dt - q_vin_l + q_vout_l                       ]   [ 0 ]
      [(p_v_l - p_ar_sys)/R_arv_l - q_vout_l ]   [ 0 ] [C_ar_sys * (d(p_ar_sys)/dt - Z_ar_sys *
d(q_vout_l)/dt) - q_vout_l + q_ar_sys                       ]   [ 0 ] [L_ar_sys/R_ar_sys +
(p_arperi_sys - p_ar_sys + Z_ar_sys * q_vout_l)/R_ar_sys + q_ar_sys             ]   [ 0 ]
      [\sum_{j} C_ar,j_sys * d(p_ar,peri_sys)/dt + \sum_{j} q_ar,j_sys - q_ar_sys,
j=spl,espl,msc,cer,cor]   [ 0 ]
      [(p_ven,spl_sys - p_ar,peri_sys)/R_ar,spl_sys + q_ar,spl_sys ]   [ 0 ]
      [(p_ven,espl_sys - p_ar,peri_sys)/R_ar,espl_sys + q_ar,espl_sys ]   [ 0 ]
      [(p_ven,msc_sys - p_ar,peri_sys)/R_ar,msc_sys + q_ar,msc_sys ]   [ 0 ]
      [(p_ven,cer_sys - p_ar,peri_sys)/R_ar,cer_sys + q_ar,cer_sys ]   [ 0 ]
      [(p_ven,cor_sys - p_ar,peri_sys)/R_ar,cor_sys + q_ar,cor_sys ]   [ 0 ] [C_ven,spl_sys *
d(p_ven,spl_sys)/dt + q_ven,spl_sys - q_ar,spl_sys                                  ]   [ 0 ]
      [(p_ven_sys - p_ven,spl_sys)/R_ven,spl_sys + q_ven,spl_sys ]   [ 0 ] [C_ven,espl_sys *
d(p_ven,espl_sys)/dt + q_ven,espl_sys - q_ar,espl_sys                              ]   [ 0 ]
      [(p_ven_sys - p_ven,espl_sys)/R_ven,espl_sys + q_ven,espl_sys ]   [ 0 ] [C_ven,msc_sys *
d(p_ven,msc_sys)/dt + q_ven,msc_sys - q_ar,msc_sys                                  ]   [ 0 ]
      [(p_ven_sys - p_ven,msc_sys)/R_ven,msc_sys + q_ven,msc_sys ]   [ 0 ] [C_ven,cer_sys *
d(p_ven,cer_sys)/dt + q_ven,cer_sys - q_ar,cer_sys                                  ]   [ 0 ]
      [(p_ven_sys - p_ven,cer_sys)/R_ven,cer_sys + q_ven,cer_sys ]   [ 0 ] [C_ven,cor_sys *
d(p_ven,cor_sys)/dt + q_ven,cor_sys - q_ar,cor_sys                                  ]   [ 0 ]
      [(p_ven_sys - p_ven,cor_sys)/R_ven,cor_sys + q_ven,cor_sys ]   [ 0 ] [C_ven_sys *
d(p_ven_sys)/dt + q_ven_sys - \sum_{j} q_ven,j_sys,   j=spl,espl,msc,cer,cor            ]   [ 0 ]
      [L_ven_sys/R_ven_sys + (p_at_r - p_ven_sys)/R_ven_sys + q_ven_sys ]   [ 0 ] Res = =
      [d(p_at_r/E_at_r)/dt - q_ven_sys + q_vin_r  OR*  d(V_at_r)/dt - q_ven_sys + q_vin_r ]   [ 0 ]
      [(p_at_r - p_v_r)/R_atv_r - q_vin_r ]   [ 0 ] [d(V_v_r)/dt - q_vin_r + q_vout_r  OR*
d(p_v_r/E_v_r)/dt - q_vin_r + q_vout_r                       ]   [ 0 ]
      [(p_v_r - p_ar_pul)/R_arv_r - q_vout_r ]   [ 0 ] [C_ar_pul * (d(p_ar_pul)/dt - Z_ar_pul *
d(q_vout_r)/dt) - q_vout_r + q_ar_pul                       ]   [ 0 ] [L_ar_pul/R_ar_pul +
(p_cap_pul - p_ar_pul + Z_ar_pul * q_vout_r)/R_ar_pul + q_ar_pul                ]   [ 0 ] [C_cap_pul
* d(p_cap_pul)/dt - q_ar_pul + q_cap_pul                                                  ]   [ 0 ]
      [(p_ven_pul - p_cap_pul)/R_cap_pul + q_cap_pul ]   [ 0 ] [C_ven_pul * d(p_ven_pul)/dt -
q_cap_pul + q_ven_pul                                                 ]   [ 0 ] [L_ven_pul/R_ven_pul
+ (p_at_l - p_ven_pul)/R_ven_pul + q_ven_pul                                    ]   [ 0 ]

* depending on atrial/ventricular model (0D elastance vs. 3D structural)

if we set RESPIRATORY_MODEL (in CARDIOVASCULAR 0D-STRUCTURE COUPLING/CARDIOVASCULAR RESPIRATORY 0D
PARAMETERS) from 'None' to 'Standard', we additionally model oxygen (O2) uptake and carbon dioxide
(CO2) release from a 0D lung model as well as O2 and CO2 transport through the vascular system, plus
potential O2 consumption and CO2 production in the systemic capillaries (differentiated in
splanchnic (spl), extra-splanchnic (espl), muscular (msc), cerebral (cer) and coronary (cor)),
according to Trenhago et al., "An integrated mathematical model of the cardiovascular and
respiratory systems", IJNMBE, 2016, while the lung model is according to Ben-Tal, "Simplified models
for gas exchange in the human lungs", J Theor biol (2006)

*----------------------------------------------------------------------*/

#ifndef CARDIOVASCULARRESPIRATORY0DSYSPULPERIPHCIRCULATION_H
#define CARDIOVASCULARRESPIRATORY0DSYSPULPERIPHCIRCULATION_H

#include <Teuchos_RCP.hpp>
#include <Epetra_Vector.h>
#include <Epetra_Operator.h>
#include <Epetra_RowMatrix.h>
#include <Teuchos_ParameterList.hpp>
#include "Epetra_FECrsMatrix.h"
#include "../drt_fem_general/drt_utils_integration.H"
#include "../drt_inpar/inpar_cardiovascular0d.H"
// class which is inherited from
#include "cardiovascular0d.H"


// forward declarations
namespace DRT
{
  class Condition;
  class Discretization;
}  // namespace DRT

namespace LINALG
{
  class SparseMatrix;
  class SparseOperator;
}  // namespace LINALG

namespace UTILS
{
  class CardiovascularRespiratory0DSysPulPeriphCirculation : public Cardiovascular0D

  {
   public:
    /*!
    \brief Constructor of a Cardiovascular0D based on conditions with a given name. It also
    takes care of the Cardiovascular0D IDs.
    */

    CardiovascularRespiratory0DSysPulPeriphCirculation(
        Teuchos::RCP<DRT::Discretization>
            discr,                    ///< Discretization where Cardiovascular0D lives on
        const std::string& CondName,  ///< Name of condition to create Cardiovascular0D from
        std::vector<int>& curID       ///< current ID
    );

    /*!
        \brief Destructor

     */
    virtual ~CardiovascularRespiratory0DSysPulPeriphCirculation() { return; };

    /// initialization routine called by the manager ctor to get correct reference base values and
    /// activating the right conditions at the beginning
    virtual void Initialize(
        Teuchos::ParameterList&
            params,  ///< parameter list to communicate between elements and discretization
        Teuchos::RCP<Epetra_Vector> sysvec1,  ///< distributed vector that may be filled by assembly
                                              ///< of element contributions
        Teuchos::RCP<Epetra_Vector>
            sysvec2  ///< distributed vector that may be filled by assembly of element contributions
    );

    //! Evaluate routine to call from outside. In here the right action is determined and the
    //! #EvaluateCardiovascular0D routine is called
    virtual void Evaluate(
        Teuchos::ParameterList&
            params,  ///< parameter list to communicate between elements and discretization
        Teuchos::RCP<LINALG::SparseMatrix> sysmat1,  ///< Cardiovascular0D stiffness matrix
        Teuchos::RCP<LINALG::SparseOperator>
            sysmat2,  ///< Cardiovascular0D offdiagonal matrix dV/dd
        Teuchos::RCP<LINALG::SparseOperator>
            sysmat3,                          ///< Cardiovascular0D offdiagonal matrix dfext/dp
        Teuchos::RCP<Epetra_Vector> sysvec1,  ///< distributed vectors that may be filled by
                                              ///< assembly of element contributions
        Teuchos::RCP<Epetra_Vector> sysvec2, Teuchos::RCP<Epetra_Vector> sysvec3,
        const Teuchos::RCP<Epetra_Vector> sysvec4, Teuchos::RCP<Epetra_Vector> sysvec5);

    // cbO2 and its derivatives
    double cbO2(double ppCO2, double ppO2);
    // w.r.t. O2
    double dcbO2_dppO2(double ppCO2, double ppO2);
    double d2cbO2_dppO22(double ppCO2, double ppO2);
    // w.r.t. CO2
    double dcbO2_dppCO2(double ppCO2, double ppO2);
    double d2cbO2_dppCO22(double ppCO2, double ppO2);
    double d2cbO2_dppO2dppCO2(double ppCO2, double ppO2);

    // cbCO2 and its derivatives
    double cbCO2(double ppCO2, double ppO2);
    // w.r.t. CO2
    double dcbCO2_dppCO2(double ppCO2, double ppO2);
    double d2cbCO2_dppCO22(double ppCO2, double ppO2);
    // w.r.t. O2
    double dcbCO2_dppO2(double ppCO2, double ppO2);
    double d2cbCO2_dppO22(double ppCO2, double ppO2);
    double d2cbCO2_dppCO2dppO2(double ppCO2, double ppO2);

    double ctO2(double ppO2);
    double dctO2_dppO2(double ppO2);
    double d2ctO2_dppO22(double ppO2);

    double SO2(double ppCO2, double ppO2);

    double ctCO2(double ppCO2);
    double dctCO2_dppCO2(double ppCO2);
    double d2ctCO2_dppCO22(double ppCO2);


    //! Evaluate routine to call from outside. In here the right action is determined and the
    //! #EvaluateCardiovascular0D routine is called
    virtual void EvaluateRespiratory(Teuchos::ParameterList& params, std::vector<double>& df_np,
        std::vector<double>& f_np, Epetra_SerialDenseMatrix& wkstiff,
        Teuchos::RCP<Epetra_Vector> dofvec, Teuchos::RCP<Epetra_Vector> volvec, bool evalstiff);

   private:
    // number of degrees of freedom for submodels
    int num_dof_cardio_;
    int num_dof_respir_;

    // parameters
    double R_arvalve_max_l_;       ///< maximum aortic valve resistance
    double R_arvalve_min_l_;       ///< minimum aortic valve resistance
    double R_atvalve_max_l_;       ///< maximum mitral valve resistance
    double R_atvalve_min_l_;       ///< minimum mitral valve resistance
    double R_arvalve_max_r_;       ///< maximum pulmonary valve resistance
    double R_arvalve_min_r_;       ///< minimum pulmonary valve resistance
    double R_atvalve_max_r_;       ///< maximum tricuspid valve resistance
    double R_atvalve_min_r_;       ///< minimum tricuspid valve resistance
    int Atrium_act_curve_l_;       ///< left atrial activation curve (ONLY for ATRIUM_MODEL "0D"!)
    int Atrium_act_curve_r_;       ///< right atrial activation curve (ONLY for ATRIUM_MODEL "0D"!)
    int Ventricle_act_curve_l_;    ///< left ventricular activation curve (ONLY for VENTRICLE_MODEL
                                   ///< "0D"!)
    int Ventricle_act_curve_r_;    ///< right ventricular activation curve (ONLY for VENTRICLE_MODEL
                                   ///< "0D"!)
    int Atrium_prescr_E_curve_l_;  ///< left atrial elastance prescription curve (ONLY for
                                   ///< ATRIUM_MODEL "prescribed"!)
    int Atrium_prescr_E_curve_r_;  ///< right atrial elastance prescription curve (ONLY for
                                   ///< ATRIUM_MODEL "prescribed"!)
    int Ventricle_prescr_E_curve_l_;  ///< left ventricular elastance prescription curve (ONLY for
                                      ///< VENTRICLE_MODEL "prescribed"!)
    int Ventricle_prescr_E_curve_r_;  ///< right ventricular elastance prescription curve (ONLY for
                                      ///< VENTRICLE_MODEL "prescribed"!)
    double E_at_max_l_;               ///< maximum left atrial elastance (for 0D atria)
    double E_at_min_l_;               ///< minimum left atrial elastance (for 0D atria)
    double E_at_max_r_;               ///< maximum right atrial elastance (for 0D atria)
    double E_at_min_r_;               ///< minimum right atrial elastance (for 0D atria)
    double E_v_max_l_;                ///< maximum left ventricular elastance (for 0D ventricles)
    double E_v_min_l_;                ///< minimum left ventricular elastance (for 0D ventricles)
    double E_v_max_r_;                ///< maximum right ventricular elastance (for 0D ventricles)
    double E_v_min_r_;                ///< minimum right ventricular elastance (for 0D ventricles)
    double C_ar_sys_;                 ///< systemic arterial compliance
    double R_ar_sys_;                 ///< systemic arterial resistance
    double L_ar_sys_;                 ///< systemic arterial inertance
    double Z_ar_sys_;                 ///< systemic arterial impedance
    // peripheral arterial compliances and resistances
    double C_arspl_sys_;   ///< systemic arterial splanchnic compliance
    double R_arspl_sys_;   ///< systemic arterial splanchnic resistance
    double C_arespl_sys_;  ///< systemic arterial extra-splanchnic compliance
    double R_arespl_sys_;  ///< systemic arterial extra-splanchnic resistance
    double C_armsc_sys_;   ///< systemic arterial muscular compliance
    double R_armsc_sys_;   ///< systemic arterial muscular resistance
    double C_arcer_sys_;   ///< systemic arterial cerebral compliance
    double R_arcer_sys_;   ///< systemic arterial cerebral resistance
    double C_arcor_sys_;   ///< systemic arterial coronary compliance
    double R_arcor_sys_;   ///< systemic arterial coronary resistance
    // peripheral venous compliances and resistances
    double C_venspl_sys_;   ///< systemic venous splanchnic compliance
    double R_venspl_sys_;   ///< systemic venous splanchnic resistance
    double C_venespl_sys_;  ///< systemic venous extra-splanchnic compliance
    double R_venespl_sys_;  ///< systemic venous extra-splanchnic resistance
    double C_venmsc_sys_;   ///< systemic venous muscular compliance
    double R_venmsc_sys_;   ///< systemic venous muscular resistance
    double C_vencer_sys_;   ///< systemic venous cerebral compliance
    double R_vencer_sys_;   ///< systemic venous cerebral resistance
    double C_vencor_sys_;   ///< systemic venous coronary compliance
    double R_vencor_sys_;   ///< systemic venous coronary resistance

    double C_ar_pul_;  ///< pulmonary arterial compliance
    double R_ar_pul_;  ///< pulmonary arterial resistance
    double L_ar_pul_;  ///< pulmonary arterial inertance
    double Z_ar_pul_;  ///< pulmonary arterial impedance

    // pulmonary capillary compliance and resistance
    double C_cap_pul_;  ///< pulmonary capillary compliance
    double R_cap_pul_;  ///< pulmonary capillary resistance

    double C_ven_sys_;  ///< systemic venous compliance
    double R_ven_sys_;  ///< systemic venous resistance
    double L_ven_sys_;  ///< systemic venous inertance
    double C_ven_pul_;  ///< pulmonary venous compliance
    double R_ven_pul_;  ///< pulmonary venous resistance
    double L_ven_pul_;  ///< pulmonary venous inertance

    // unstressed 0D volumes
    double V_v_l_u_;     ///< dead left ventricular volume (at zero pressure) (for 0D ventricles)
    double V_at_l_u_;    ///< dead left atrial volume (at zero pressure) (for 0D atria)
    double V_ar_sys_u_;  ///< unstressed systemic arterial volume

    double V_arspl_sys_u_;    ///< unstressed systemic arterial splanchnic volume
    double V_arespl_sys_u_;   ///< unstressed systemic arterial extra-splanchnic volume
    double V_armsc_sys_u_;    ///< unstressed systemic arterial muscular volume
    double V_arcer_sys_u_;    ///< unstressed systemic arterial cerebral volume
    double V_arcor_sys_u_;    ///< unstressed systemic arterial coronary volume
    double V_venspl_sys_u_;   ///< unstressed systemic venous splanchnic volume
    double V_venespl_sys_u_;  ///< unstressed systemic venous extra-splanchnic volume
    double V_venmsc_sys_u_;   ///< unstressed systemic venous muscular volume
    double V_vencer_sys_u_;   ///< unstressed systemic venous cerebral volume
    double V_vencor_sys_u_;   ///< unstressed systemic venous coronary volume

    double V_ven_sys_u_;  ///< unstressed systemic venous volume
    double V_v_r_u_;      ///< unstressed  right ventricular volume (for 0D ventricles)
    double V_at_r_u_;     ///< unstressed  right atrial volume (for 0D atria)
    double V_ar_pul_u_;   ///< unstressed pulmonary arterial volume
    double V_cap_pul_u_;  ///< unstressed pulmonary capillary volume
    double V_ven_pul_u_;  ///< unstressed pulmonary venous volume

    // total number of degrees of freedom
    int num_dof_;

    // 0D lung
    double L_alv_;   ///< alveolar inertance
    double R_alv_;   ///< alveolar resistance
    double E_alv_;   ///< alveolar elastance
    int U_t_curve_;  ///< time-varying, prescribed pleural pressure curve driven by diaphragm
    double U_m_;     ///< in-breath pressure

    double V_lung_tidal_;  ///< tidal volume (the total volume of inspired air, in a single breath)
    double V_lung_dead_;   ///< dead space volume
    double V_lung_u_;  ///< unstressed lung volume (volume of the lung when it is fully collapsed
                       ///< outside the body)

    double fCO2_ext_;  ///< atmospheric CO2 gas fraction
    double fO2_ext_;   ///< atmospheric O2 gas fraction

    // should be 22.4 liters per mol !
    // however we specify it as an input parameter since its decimal power depends on the system of
    // units your whole model is specified in! i.e. if you have kg - mm - s - mmol, it's 22.4e3 mm^3
    // / mmol
    double V_m_gas_;  ///< molar volume of an ideal gas

    // should be 47.1 mmHg = 6.279485 kPa !
    // however we specify it as an input parameter since its decimal power depends on the system of
    // units your whole model is specified in! i.e. if you have kg - mm - s - mmol, it's 6.279485
    // kPa
    double p_vap_water_37_;  ///< vapor pressure of water at 37 °C

    double kappa_CO2_;  ///< diffusion coefficient for CO2 across the hemato-alveolar membrane, in
                        ///< molar value / (time * pressure)
    double kappa_O2_;   ///< diffusion coefficient for O2 across the hemato-alveolar membrane, in
                        ///< molar value / (time * pressure)

    double alpha_CO2_;  ///< CO2 solubility constant, in molar value / (volume * pressure)
    double alpha_O2_;   ///< O2 solubility constant, in molar value / (volume * pressure)

    double c_Hb_;  ///< hemoglobin concentration of the blood, in molar value / volume

    double M_CO2_arspl_;   ///< splanchnic metabolic rate of CO2 production
    double M_O2_arspl_;    ///< splanchnic metabolic rate of O2 consumption
    double M_CO2_arespl_;  ///< extra-splanchnic metabolic rate of CO2 production
    double M_O2_arespl_;   ///< extra-splanchnic metabolic rate of O2 consumption
    double M_CO2_armsc_;   ///< muscular metabolic rate of CO2 production
    double M_O2_armsc_;    ///< muscular metabolic rate of O2 consumption
    double M_CO2_arcer_;   ///< cerebral metabolic rate of CO2 production
    double M_O2_arcer_;    ///< cerebral metabolic rate of O2 consumption
    double M_CO2_arcor_;   ///< coronary metabolic rate of CO2 production
    double M_O2_arcor_;    ///< coronary metabolic rate of O2 consumption

    // tissue columes
    double V_tissspl_;   ///< splanchnic tissue volume
    double V_tissespl_;  ///< extra-splanchnic tissue volume
    double V_tissmsc_;   ///< muscular tissue volume
    double V_tisscer_;   ///< cerebral tissue volume
    double V_tisscor_;   ///< coronary tissue volume



    // don't want = operator, cctor and destructor

    CardiovascularRespiratory0DSysPulPeriphCirculation operator=(
        const CardiovascularRespiratory0DSysPulPeriphCirculation& old);
    CardiovascularRespiratory0DSysPulPeriphCirculation(
        const CardiovascularRespiratory0DSysPulPeriphCirculation& old);



  };  // class
}  // namespace UTILS

#endif /*CARDIOVASCULARRESPIRATORY0DSYSPULPERIPHCIRCULATION_H*/
