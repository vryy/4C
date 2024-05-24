/*--------------------------------------------------------------------------*/
/*! \file

\brief supplementary element calculation class providing general utility for thermodynamic scalar
transport


\level 2
*/
/*--------------------------------------------------------------------------*/
#ifndef FOUR_C_SCATRA_ELE_STI_THERMO_HPP
#define FOUR_C_SCATRA_ELE_STI_THERMO_HPP

#include "4C_config.hpp"

#include "4C_scatra_ele_calc.hpp"

FOUR_C_NAMESPACE_OPEN

namespace DRT
{
  namespace ELEMENTS
  {
    // forward declaration
    class ScaTraEleDiffManagerSTIThermo;

    // class implementation
    template <CORE::FE::CellType distype>
    class ScaTraEleSTIThermo
    {
     public:
      virtual ~ScaTraEleSTIThermo() = default;

     protected:
      using my = ScaTraEleCalc<distype>;

      //! number of element nodes
      static constexpr int nen_ = my::nen_;

      //! number of space dimensions
      static constexpr int nsd_ = my::nsd_;

      //! protected constructor for singletons
      ScaTraEleSTIThermo(const int& numscal  //!< number of transported scalars
      );

      //! extract quantities for element evaluation
      virtual void extract_element_and_node_values(DRT::Element* ele,  //!< current element
          Teuchos::ParameterList& params,                              //!< parameter list
          DRT::Discretization& discretization,                         //!< discretization
          DRT::Element::LocationArray& la                              //!< location array
      );

      //! provide element matrix with linearizations of Soret effect term in discrete scatra
      //! residuals w.r.t. scatra dofs
      void CalcMatSoret(CORE::LINALG::SerialDenseMatrix& emat,  //!< element matrix
          const double& timefacfac,  //!< domain integration factor times time integration factor
          const double& conc,        //!< concentration
          const double& diffcoeff,   //!< diffusion coefficient
          const double&
              diffcoeffderiv,  //!< derivative of diffusion coefficient w.r.t. concentration
          const double& temp,  //!< temperature
          const CORE::LINALG::Matrix<nsd_, 1>& gradtemp,  //!< gradient of temperature
          const CORE::LINALG::Matrix<nen_, 1>& funct,     //!< shape functions
          const CORE::LINALG::Matrix<nsd_, nen_>& derxy  //!< spatial derivatives of shape functions
      );

      //! provide element matrix with linearizations of Soret effect term in discrete scatra
      //! residuals w.r.t. thermo dofs
      void CalcMatSoretOD(CORE::LINALG::SerialDenseMatrix& emat,  //!< element matrix
          const double& timefacfac,     //!< time integration factor times domain integration factor
          const double& concentration,  //!< concentration
          const double& diffcoeff,      //!< diffusion coefficient
          const double& temperature,    //!< temperature
          const CORE::LINALG::Matrix<nsd_, 1>& gradtemp,  //!< gradient of temperature
          const CORE::LINALG::Matrix<nen_, 1>& funct,     //!< shape functions
          const CORE::LINALG::Matrix<nsd_, nen_>& derxy  //!< spatial derivatives of shape functions
      );

      //! provide element right-hand side vector with contributions of Soret effect term to discrete
      //! scatra residuals
      void CalcRHSSoret(CORE::LINALG::SerialDenseVector& erhs,  //!< element right-hand side vector
          const double& rhsfac,     //!< domain integration factor times time integration factor for
                                    //!< right-hand side vector
          const double& conc,       //!< concentration
          const double& diffcoeff,  //!< diffusion coefficient
          const double& temp,       //!< temperature
          const CORE::LINALG::Matrix<nsd_, 1>& gradtemp,  //!< gradient of temperature
          const CORE::LINALG::Matrix<nsd_, nen_>& derxy  //!< spatial derivatives of shape functions
      );

      //! provide element matrix with linearization of diffusion coefficient with respect to
      //! discrete temperature
      //!
      //! \param emat           element matrix
      //! \param numdofpernode  number of dof per node
      //! \param timefacfac     time integration factor times domain integration factor
      //! \param invF  1/F      with F being the faraday constant
      //! \param gradconc       spatial derivative of concentration
      //! \param gradpot        spatial derivative of electric potential
      //! \param tempderivisodiffcoef  derivative of isotropic diffusion coefficient wrt temperature
      //! \param tempderivcond         derivative of electric conductivity wrt temperature
      //! \param funct     shape functions
      //! \param derxy     spatial derivatives of shape functions
      //! \param scalefac  scaling factor for pot. contributions
      void CalcMatDiffThermoOD(CORE::LINALG::SerialDenseMatrix& emat, const int& numdofpernode,
          const double& timefacfac, const double& invF,
          const CORE::LINALG::Matrix<nsd_, 1>& gradconc,
          const CORE::LINALG::Matrix<nsd_, 1>& gradpot, const double& tempderivisodiffcoef,
          const double& tempderivcond, const CORE::LINALG::Matrix<nen_, 1>& funct,
          const CORE::LINALG::Matrix<nsd_, nen_>& derxy, const double& scalefac);

      //! evaluate Soret material
      void mat_soret(const Teuchos::RCP<const CORE::MAT::Material> material  //!< Soret material
      );

      //! compute gradient of test function times gradient of shape function
      void get_laplacian_weak_form(double& result,       //!< result variable
          const int& vi,                                 //!< index of test function
          const int& ui,                                 //!< index of shape function
          const CORE::LINALG::Matrix<nsd_, nen_>& derxy  //!< spatial derivatives of shape functions
      )
      {
        // initialize result variable
        result = 0.;

        // compute gradient of test function times gradient of shape function
        for (int idim = 0; idim < nsd_; ++idim) result += derxy(idim, vi) * derxy(idim, ui);
      };

      //! compute gradient of test function times given gradient
      void get_laplacian_weak_form_rhs(double& result,    //!< result variable
          const int& vi,                                  //!< index of test function
          const CORE::LINALG::Matrix<nsd_, 1>& gradient,  //!< given gradient
          const CORE::LINALG::Matrix<nsd_, nen_>& derxy  //!< spatial derivatives of shape functions
      )
      {
        // initialize result variable
        result = 0.;

        // compute gradient of test function times given gradient
        for (int idim = 0; idim < nsd_; ++idim) result += derxy(idim, vi) * gradient(idim);
      };

      //! local nodal values of temperature
      CORE::LINALG::Matrix<nen_, 1> etempnp_;

      //! thermo diffusion manager
      const Teuchos::RCP<ScaTraEleDiffManagerSTIThermo> diffmanagerstithermo_;
    };  // class ScaTraEleSTIThermo


    //! implementation of ScaTraEleDiffManagerSTIThermo
    class ScaTraEleDiffManagerSTIThermo : public ScaTraEleDiffManager
    {
     public:
      //! constructor
      ScaTraEleDiffManagerSTIThermo(const int& numscal  //!< number of transported scalars
          )
          :  // constructor of base class
            ScaTraEleDiffManager(numscal),

            // initialize internal member variable
            soret_(0.){};

      //! set Soret coefficient
      void SetSoret(const double& soret  //!< Soret coefficient
      )
      {
        soret_ = soret;
      };

      //! get Soret coefficient
      const double& GetSoret() const { return soret_; };

     protected:
      double soret_;
    };  // class ScaTraEleDiffManagerSTIThermo


    // implementation of ScaTraEleInternalVariableManagerSTIThermo
    template <int NSD, int NEN>
    class ScaTraEleInternalVariableManagerSTIThermo
    {
     public:
      //! constructor
      ScaTraEleInternalVariableManagerSTIThermo() : temp_(0.), gradtemp_(true){};

      //! set internal variables for element evaluation
      void set_internal_variables_sti_thermo(
          const CORE::LINALG::Matrix<NEN, 1>& funct,    //!< shape functions
          const CORE::LINALG::Matrix<NSD, NEN>& derxy,  //!< spatial derivatives of shape functions
          const CORE::LINALG::Matrix<NEN, 1>&
              etempnp  //!< nodal temperature values at time t_(n+1) or t_(n+alpha_F)
      )
      {
        // set local values of thermo variables at time t_(n+1) or t_(n+alpha_F)
        temp_ = funct.Dot(etempnp);  // temperature
        gradtemp_.Multiply(derxy, etempnp);
      };

      //! return temperature
      const double& Temp() const { return temp_; };

      //! return gradient of temperature
      const CORE::LINALG::Matrix<NSD, 1>& GradTemp() const { return gradtemp_; };

     protected:
      //! temperature
      double temp_;

      //! gradient of temperature
      CORE::LINALG::Matrix<NSD, 1> gradtemp_;
    };
  }  // namespace ELEMENTS
}  // namespace DRT
FOUR_C_NAMESPACE_CLOSE

#endif
