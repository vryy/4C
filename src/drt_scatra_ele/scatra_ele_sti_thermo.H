/*--------------------------------------------------------------------------*/
/*!
\file scatra_ele_sti_thermo.H

\brief supplementary element calculation class providing general utility for thermodynamic scalar
transport

<pre>
Maintainer: Rui Fang
            fang@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/
            089-289-15251
</pre>
*/
/*--------------------------------------------------------------------------*/
#ifndef SCATRA_ELE_STI_THERMO_H
#define SCATRA_ELE_STI_THERMO_H

#include "../drt_scatra_ele/scatra_ele_calc.H"

namespace DRT
{
  namespace ELEMENTS
  {
    // forward declaration
    class ScaTraEleDiffManagerSTIThermo;

    // class implementation
    template <DRT::Element::DiscretizationType distype>
    class ScaTraEleSTIThermo
    {
     public:
      //! singleton destruction
      virtual void Done() = 0;

     protected:
      //! number of element nodes
      static const int nen_ = DRT::UTILS::DisTypeToNumNodePerEle<distype>::numNodePerElement;

      //! number of space dimensions
      static const int nsd_ = DRT::UTILS::DisTypeToDim<distype>::dim;

      //! protected constructor for singletons
      ScaTraEleSTIThermo(const int& numscal  //!< number of transported scalars
      );

      //! destructor
      virtual ~ScaTraEleSTIThermo() { return; };

      //! extract quantities for element evaluation
      virtual void ExtractElementAndNodeValues(DRT::Element* ele,  //!< current element
          Teuchos::ParameterList& params,                          //!< parameter list
          DRT::Discretization& discretization,                     //!< discretization
          DRT::Element::LocationArray& la                          //!< location array
      );

      //! provide element matrix with linearizations of Soret effect term in discrete scatra
      //! residuals w.r.t. scatra dofs
      void CalcMatSoret(Epetra_SerialDenseMatrix& emat,  //!< element matrix
          const double& timefacfac,  //!< domain integration factor times time integration factor
          const double& conc,        //!< concentration
          const double& diffcoeff,   //!< diffusion coefficient
          const double&
              diffcoeffderiv,  //!< derivative of diffusion coefficient w.r.t. concentration
          const double& temp,  //!< temperature
          const LINALG::Matrix<nsd_, 1>& gradtemp,  //!< gradient of temperature
          const LINALG::Matrix<nen_, 1>& funct,     //!< shape functions
          const LINALG::Matrix<nsd_, nen_>& derxy   //!< spatial derivatives of shape functions
      );

      //! provide element matrix with linearizations of Soret effect term in discrete scatra
      //! residuals w.r.t. thermo dofs
      void CalcMatSoretOD(Epetra_SerialDenseMatrix& emat,  //!< element matrix
          const double& timefacfac,     //!< time integration factor times domain integration factor
          const double& concentration,  //!< concentration
          const double& diffcoeff,      //!< diffusion coefficient
          const double& temperature,    //!< temperature
          const LINALG::Matrix<nsd_, 1>& gradtemp,  //!< gradient of temperature
          const LINALG::Matrix<nen_, 1>& funct,     //!< shape functions
          const LINALG::Matrix<nsd_, nen_>& derxy   //!< spatial derivatives of shape functions
      );

      //! provide element right-hand side vector with contributions of Soret effect term to discrete
      //! scatra residuals
      void CalcRHSSoret(Epetra_SerialDenseVector& erhs,  //!< element right-hand side vector
          const double& rhsfac,     //!< domain integration factor times time integration factor for
                                    //!< right-hand side vector
          const double& conc,       //!< concentration
          const double& diffcoeff,  //!< diffusion coefficient
          const double& temp,       //!< temperature
          const LINALG::Matrix<nsd_, 1>& gradtemp,  //!< gradient of temperature
          const LINALG::Matrix<nsd_, nen_>& derxy   //!< spatial derivatives of shape functions
      );

      //! evaluate Soret material
      void MatSoret(const Teuchos::RCP<const MAT::Material> material  //!< Soret material
      );

      //! compute gradient of test function times gradient of shape function
      void GetLaplacianWeakForm(double& result,    //!< result variable
          const int& vi,                           //!< index of test function
          const int& ui,                           //!< index of shape function
          const LINALG::Matrix<nsd_, nen_>& derxy  //!< spatial derivatives of shape functions
      )
      {
        // initialize result variable
        result = 0.;

        // compute gradient of test function times gradient of shape function
        for (int idim = 0; idim < nsd_; ++idim) result += derxy(idim, vi) * derxy(idim, ui);

        return;
      };

      //! compute gradient of test function times given gradient
      void GetLaplacianWeakFormRHS(double& result,  //!< result variable
          const int& vi,                            //!< index of test function
          const LINALG::Matrix<nsd_, 1>& gradient,  //!< given gradient
          const LINALG::Matrix<nsd_, nen_>& derxy   //!< spatial derivatives of shape functions
      )
      {
        // initialize result variable
        result = 0.;

        // compute gradient of test function times given gradient
        for (int idim = 0; idim < nsd_; ++idim) result += derxy(idim, vi) * gradient(idim);

        return;
      };

      //! local nodal values of temperature
      LINALG::Matrix<nen_, 1> etempnp_;

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
            soret_(0.)
      {
        return;
      };

      //! destructor
      virtual ~ScaTraEleDiffManagerSTIThermo() { return; };

      //! set Soret coefficient
      void SetSoret(const double& soret  //!< Soret coefficient
      )
      {
        soret_ = soret;
        return;
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
      ScaTraEleInternalVariableManagerSTIThermo() : temp_(0.), gradtemp_(true) { return; };

      //! destructor
      virtual ~ScaTraEleInternalVariableManagerSTIThermo() { return; };

      //! set internal variables for element evaluation
      void SetInternalVariablesSTIThermo(const LINALG::Matrix<NEN, 1>& funct,  //!< shape functions
          const LINALG::Matrix<NSD, NEN>& derxy,  //!< spatial derivatives of shape functions
          const LINALG::Matrix<NEN, 1>&
              etempnp  //!< nodal temperature values at time t_(n+1) or t_(n+alpha_F)
      )
      {
        // set local values of thermo variables at time t_(n+1) or t_(n+alpha_F)
        temp_ = funct.Dot(etempnp);          // temperature
        gradtemp_.Multiply(derxy, etempnp);  // gradient of temperature

        return;
      };

      //! return temperature
      const double& Temp() const { return temp_; };

      //! return gradient of temperature
      const LINALG::Matrix<NSD, 1>& GradTemp() const { return gradtemp_; };

     protected:
      //! temperature
      double temp_;

      //! gradient of temperature
      LINALG::Matrix<NSD, 1> gradtemp_;
    };
  }  // namespace ELEMENTS
}  // namespace DRT
#endif
