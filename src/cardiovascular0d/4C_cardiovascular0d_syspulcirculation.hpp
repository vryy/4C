#ifndef FOUR_C_CARDIOVASCULAR0D_SYSPULCIRCULATION_HPP
#define FOUR_C_CARDIOVASCULAR0D_SYSPULCIRCULATION_HPP

#include "4C_config.hpp"

#include "4C_cardiovascular0d.hpp"
#include "4C_fem_general_utils_integration.hpp"
#include "4C_inpar_cardiovascular0d.hpp"
#include "4C_linalg_vector.hpp"
#include "4C_utils_parameter_list.fwd.hpp"

#include <Epetra_FECrsMatrix.h>
#include <Epetra_Operator.h>
#include <Epetra_RowMatrix.h>
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

namespace Utils
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
    void initialize(
        Teuchos::ParameterList&
            params,  ///< parameter list to communicate between elements and discretization
        Teuchos::RCP<Core::LinAlg::Vector<double>>
            sysvec1,  ///< distributed vector that may be filled by
                      ///< assembly of element contributions
        Teuchos::RCP<Core::LinAlg::Vector<double>>
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
            sysmat3,  ///< Cardiovascular0D offdiagonal matrix dfext/dp
        Teuchos::RCP<Core::LinAlg::Vector<double>>
            sysvec1,  ///< distributed vectors that may be filled by
                      ///< assembly of element contributions
        Teuchos::RCP<Core::LinAlg::Vector<double>> sysvec2,
        Teuchos::RCP<Core::LinAlg::Vector<double>> sysvec3,
        const Teuchos::RCP<Core::LinAlg::Vector<double>> sysvec4,
        Teuchos::RCP<Core::LinAlg::Vector<double>> sysvec5) override;

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
}  // namespace Utils

FOUR_C_NAMESPACE_CLOSE

#endif
