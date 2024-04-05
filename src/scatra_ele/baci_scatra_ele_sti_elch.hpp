/*--------------------------------------------------------------------------*/
/*! \file

\brief supplementary element calculation class providing general utility for evaluation of heat
transport within electrochemical substances


\level 2
*/
/*--------------------------------------------------------------------------*/
#ifndef FOUR_C_SCATRA_ELE_STI_ELCH_HPP
#define FOUR_C_SCATRA_ELE_STI_ELCH_HPP

#include "baci_config.hpp"

#include "baci_scatra_ele_calc.hpp"

BACI_NAMESPACE_OPEN

namespace DRT
{
  namespace ELEMENTS
  {
    // class implementation
    template <CORE::FE::CellType distype>
    class ScaTraEleSTIElch
    {
     public:
      //! Virtual destructor.
      virtual ~ScaTraEleSTIElch() = default;

     protected:
      //! number of element nodes
      static constexpr int nen_ = CORE::FE::num_nodes<distype>;

      //! protected constructor for singletons
      ScaTraEleSTIElch(const int numdofpernode, const int numscal, const std::string& disname);

      //! element matrix and right-hand side vector contributions arising from thermal source terms
      //! in discrete thermo residuals
      void CalcMatAndRhsSource(CORE::LINALG::SerialDenseMatrix& emat,  //!< element matrix
          CORE::LINALG::SerialDenseVector& erhs,  //!< element right-hand side vector
          const double& timefacfac,  //!< domain integration factor times time integration factor
          const double& rhsfac  //!< domain integration factor times time integration factor for
                                //!< right-hand side vector
      );

      //! provide element matrix with linearizations of source terms in discrete thermo residuals
      //! w.r.t. scatra dofs
      void CalcMatSourceOD(CORE::LINALG::SerialDenseMatrix& emat,  //!< element matrix
          const double& timefacfac  //!< domain integration factor times time integration factor
      );

      //! element matrix and right-hand side vector contributions arising from Joule's heat
      virtual void CalcMatAndRhsJoule(CORE::LINALG::SerialDenseMatrix& emat,  //!< element matrix
          CORE::LINALG::SerialDenseVector& erhs,  //!< element right-hand side vector
          const double& timefacfac,  //!< domain integration factor times time integration factor
          const double& rhsfac  //!< domain integration factor times time integration factor for
                                //!< right-hand side vector
          ) = 0;

      //! element matrix and right-hand side vector contributions arising from heat of mixing
      virtual void CalcMatAndRhsMixing(CORE::LINALG::SerialDenseMatrix& emat,  //!< element matrix
          CORE::LINALG::SerialDenseVector& erhs,  //!< element right-hand side vector
          const double& timefacfac,  //!< domain integration factor times time integration factor
          const double& rhsfac  //!< domain integration factor times time integration factor for
                                //!< right-hand side vector
          ) = 0;

      //! element matrix and right-hand side vector contributions arising from Soret effect
      virtual void CalcMatAndRhsSoret(CORE::LINALG::SerialDenseMatrix& emat,  //!< element matrix
          CORE::LINALG::SerialDenseVector& erhs,  //!< element right-hand side vector
          const double& timefacfac,  //!< domain integration factor times time integration factor
          const double& rhsfac  //!< domain integration factor times time integration factor for
                                //!< right-hand side vector
          ) = 0;

      //! provide element matrix with linearizations of Joule's heat term in discrete thermo
      //! residuals w.r.t. scatra dofs
      virtual void CalcMatJouleOD(CORE::LINALG::SerialDenseMatrix& emat,  //!< element matrix
          const double& timefacfac  //!< domain integration factor times time integration factor
          ) = 0;

      //! provide element matrix with linearizations of heat of mixing term in discrete thermo
      //! residuals w.r.t. scatra dofs
      virtual void CalcMatMixingOD(CORE::LINALG::SerialDenseMatrix& emat,  //!< element matrix
          const double& timefacfac  //!< domain integration factor times time integration factor
          ) = 0;

      //! provide element matrix with linearizations of Soret effect term in discrete thermo
      //! residuals w.r.t. scatra dofs
      virtual void CalcMatSoretOD(CORE::LINALG::SerialDenseMatrix& emat,  //!< element matrix
          const double& timefacfac  //!< domain integration factor times time integration factor
          ) = 0;

      //! extract quantities for element evaluation
      virtual void ExtractElementAndNodeValues(DRT::Element* ele,  //!< current element
          Teuchos::ParameterList& params,                          //!< parameter list
          DRT::Discretization& discretization,                     //!< discretization
          DRT::Element::LocationArray& la                          //!< location array
      );

      //! local nodal values of concentration
      CORE::LINALG::Matrix<nen_, 1> econcnp_;

      //! local nodal values of electric potential
      CORE::LINALG::Matrix<nen_, 1> epotnp_;
    };  // class ScaTraEleSTIElch


    //! implementation of ScaTraEleInternalVariableManagerSTIElch
    template <int NSD, int NEN>
    class ScaTraEleInternalVariableManagerSTIElch
        : public ScaTraEleInternalVariableManager<NSD, NEN>
    {
     public:
      //! abbreviation
      typedef ScaTraEleInternalVariableManager<NSD, NEN> vm;

      //! constructor
      ScaTraEleInternalVariableManagerSTIElch(const int& numscal)
          :  // call base class constructor
            ScaTraEleInternalVariableManager<NSD, NEN>(numscal),

            // initialize internal member variables
            conc_(0.),
            gradconc_(true),
            gradpot_(true)
      {
        return;
      };

      //! set internal variables for element evaluation
      void SetInternalVariablesSTIElch(
          const CORE::LINALG::Matrix<NEN, 1>& funct,    //!< shape functions
          const CORE::LINALG::Matrix<NSD, NEN>& derxy,  //!< spatial derivatives of shape functions
          const std::vector<CORE::LINALG::Matrix<NEN, 1>>&
              etempnp,  //!< nodal temperature values at time t_(n+1) or t_(n+alpha_F)
          const std::vector<CORE::LINALG::Matrix<NEN, 1>>&
              etempn,  //!< nodal temperature values at time t_(n)
          const CORE::LINALG::Matrix<NEN, 1>&
              econcnp,  //!< nodal concentration values at time t_(n+1) or t_(n+alpha_F)
          const CORE::LINALG::Matrix<NEN, 1>&
              epotnp,  //!< nodal electric potential values at time t_(n+1) or t_(n+alpha_F)
          const CORE::LINALG::Matrix<NSD, NEN>&
              econvelnp,  //!< nodal convective velocity values at time t_(n+1) or t_(n+alpha_F)
          const std::vector<CORE::LINALG::Matrix<NEN, 1>>& ehist  //!< nodal history values
      )
      {
        // call base class routine to set thermo variables
        const CORE::LINALG::Matrix<NSD, NEN> eforcevelocity(true);
        vm::SetInternalVariables(funct, derxy, etempnp, etempn, econvelnp, ehist, eforcevelocity);

        // set local values of scatra variables at time t_(n+1) or t_(n+alpha_F)
        conc_ = funct.Dot(econcnp);          // concentration
        gradconc_.Multiply(derxy, econcnp);  // gradient of concentration
        gradpot_.Multiply(derxy, epotnp);    // gradient of electric potential

        return;
      }

      //! return concentration
      const double& Conc() const { return conc_; };

      //! return gradient of concentration
      const CORE::LINALG::Matrix<NSD, 1>& GradConc() const { return gradconc_; };

      //! return gradient of electric potential
      const CORE::LINALG::Matrix<NSD, 1>& GradPot() const { return gradpot_; };

      //! concentration
      double conc_;

      //! gradient of concentration
      CORE::LINALG::Matrix<NSD, 1> gradconc_;

      //! gradient of electric potential
      CORE::LINALG::Matrix<NSD, 1> gradpot_;
    };  // class ScaTraEleInternalVariableManagerSTIElch
  }     // namespace ELEMENTS
}  // namespace DRT
BACI_NAMESPACE_CLOSE

#endif
