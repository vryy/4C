/*--------------------------------------------------------------------------*/
/*!
\file scatra_ele_calc_elch_electrode_sti_thermo.H

\brief evaluation of scatra elements for conservation of mass concentration and electronic charge within thermodynamic electrodes

\level 2

<pre>
\maintainer Rui Fang
            fang@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/
            089-289-15251
</pre>
*/
/*--------------------------------------------------------------------------*/
#ifndef SCATRA_ELE_CALC_ELCH_ELECTRODE_STI_THERMO_H
#define SCATRA_ELE_CALC_ELCH_ELECTRODE_STI_THERMO_H

#include "scatra_ele_sti_thermo.H"

#include "../drt_scatra_ele/scatra_ele_calc_elch_electrode.H"

namespace DRT
{
  namespace ELEMENTS
  {
    // forward declaration
    template<int NSD, int NEN> class ScaTraEleInternalVariableManagerElchElectrodeSTIThermo;

    // class implementation
    template<DRT::Element::DiscretizationType distype>
    class ScaTraEleCalcElchElectrodeSTIThermo : public ScaTraEleCalcElchElectrode<distype>, public ScaTraEleSTIThermo<distype>
    {
      public:

        //! singleton access method
        static ScaTraEleCalcElchElectrodeSTIThermo<distype>* Instance(
            const int numdofpernode,
            const int numscal,
            const std::string& disname,
            const ScaTraEleCalcElchElectrodeSTIThermo* delete_me = NULL
            );

        //! singleton destruction
        void Done();

      private:

        //! abbreviations
        typedef ScaTraEleCalc<distype> my;
        typedef ScaTraEleCalcElch<distype> myelch;
        typedef ScaTraEleCalcElchElectrode<distype> myelectrode;
        typedef ScaTraEleSTIThermo<distype> mythermo;

        //! private constructor for singletons
        ScaTraEleCalcElchElectrodeSTIThermo(
            const int numdofpernode,
            const int numscal,
            const std::string& disname
            );

        //! destructor
        virtual ~ScaTraEleCalcElchElectrodeSTIThermo(){return;};

        //! evaluate action for off-diagonal system matrix block
        int EvaluateActionOD(
            DRT::Element*                  ele,              //!< current element
            Teuchos::ParameterList&        params,           //!< parameter list
            DRT::Discretization&           discretization,   //!< discretization
            const SCATRA::Action&          action,           //!< action parameter
            DRT::Element::LocationArray&   la,               //!< location array
            Epetra_SerialDenseMatrix&      elemat1_epetra,   //!< element matrix 1
            Epetra_SerialDenseMatrix&      elemat2_epetra,   //!< element matrix 2
            Epetra_SerialDenseVector&      elevec1_epetra,   //!< element right-hand side vector 1
            Epetra_SerialDenseVector&      elevec2_epetra,   //!< element right-hand side vector 2
            Epetra_SerialDenseVector&      elevec3_epetra    //!< element right-hand side vector 3
            );

        //! extract quantities for element evaluation
        void ExtractElementAndNodeValues(
            DRT::Element*                  ele,              //!< current element
            Teuchos::ParameterList&        params,           //!< parameter list
            DRT::Discretization&           discretization,   //!< discretization
            DRT::Element::LocationArray&   la                //!< location array
            );

        //! get material parameters
        void GetMaterialParams(
            const DRT::Element*   ele,           //!< current element
            std::vector<double>&  densn,         //!< density at t_(n)
            std::vector<double>&  densnp,        //!< density at t_(n+1) or t_(n+alpha_F)
            std::vector<double>&  densam,        //!< density at t_(n+alpha_M)
            double&               visc,          //!< fluid viscosity
            const int             iquad          //!< ID of current integration point
            );

        //! calculate element matrix and element right-hand side vector
        void CalcMatAndRhs(
            Epetra_SerialDenseMatrix&     emat,         //!< element matrix to calculate
            Epetra_SerialDenseVector&     erhs,         //!< element rhs to calculate+
            const int                     k,            //!< index of current scalar
            const double                  fac,          //!< domain-integration factor
            const double                  timefacfac,   //!< domain-integration factor times time-integration factor
            const double                  rhsfac,       //!< time-integration factor for rhs times domain-integration factor
            const double                  taufac,       //!< tau times domain-integration factor
            const double                  timetaufac,   //!< domain-integration factor times tau times time-integration factor
            const double                  rhstaufac,    //!< time-integration factor for rhs times tau times domain-integration factor
            LINALG::Matrix<my::nen_,1>&   tauderpot,    //!< derivatives of stabilization parameter w.r.t. electric potential
            double&                       rhsint        //!< rhs at Gauss point
            );

        //! fill element matrix with linearizations of discrete scatra residuals w.r.t. thermo dofs
        void SysmatODScatraThermo(
            DRT::Element*               ele,   //!< current element
            Epetra_SerialDenseMatrix&   emat   //!< element matrix
            );

        //! set internal variables for element evaluation
        void SetInternalVariablesForMatAndRHS();

        //! get internal variable manager for thermodynamic electrodes
        Teuchos::RCP<ScaTraEleInternalVariableManagerElchElectrodeSTIThermo<my::nsd_,my::nen_> > VarManager(){return Teuchos::rcp_static_cast<ScaTraEleInternalVariableManagerElchElectrodeSTIThermo<my::nsd_,my::nen_> >(my::scatravarmanager_);};
    }; // class ScaTraEleCalcElchElectrodeSTIThermo


    //! implementation of ScaTraEleInternalVariableManagerElchElectrodeSTIThermo
    template<int NSD,int NEN>
    class ScaTraEleInternalVariableManagerElchElectrodeSTIThermo : public ScaTraEleInternalVariableManagerElchElectrode<NSD,NEN>, public ScaTraEleInternalVariableManagerSTIThermo<NSD,NEN>
    {
      public:

        //! abbreviations
        typedef ScaTraEleInternalVariableManagerElch<NSD,NEN> vmelch;
        typedef ScaTraEleInternalVariableManagerElchElectrode<NSD,NEN> vmelchelectrode;
        typedef ScaTraEleInternalVariableManagerSTIThermo<NSD,NEN> vmthermo;

        //! constructor
        ScaTraEleInternalVariableManagerElchElectrodeSTIThermo(
            int numscal,
            const DRT::ELEMENTS::ScaTraEleParameterElch* elchpara
            ) :
            // call base class constructors
            ScaTraEleInternalVariableManagerElchElectrode<NSD,NEN>(numscal,elchpara),
            ScaTraEleInternalVariableManagerSTIThermo<NSD,NEN>()
        {return;};

        //! destructor
        virtual ~ScaTraEleInternalVariableManagerElchElectrodeSTIThermo(){return;};

        //! set internal variables for element evaluation
        void SetInternalVariables(
            const LINALG::Matrix<NEN,1>&                 funct,       //!< shape functions
            const LINALG::Matrix<NSD,NEN>&               derxy,       //!< spatial derivatives of shape functions
            const std::vector<LINALG::Matrix<NEN,1> >&   ephinp,      //!< nodal concentration and electric potential values at time t_(n+1) or t_(n+alpha_F)
            const std::vector<LINALG::Matrix<NEN,1> >&   ephin,       //!< nodal concentration and electric potential values at time t_(n)
            const LINALG::Matrix<NEN,1>&                 etempnp,     //!< nodal temperature values at time t_(n+1) or t_(n+alpha_F)
            const LINALG::Matrix<NSD,NEN>&               econvelnp,   //!< nodal convective velocity values at time t_(n+1) or t_(n+alpha_F)
            const std::vector<LINALG::Matrix<NEN,1> >&   ehist        //!< nodal history values
        )
        {
          // set thermo variables
          vmthermo::SetInternalVariablesSTIThermo(funct,derxy,etempnp);

          // set scatra variables
          // this requires the temperature to be already set
          vmelchelectrode::SetInternalVariablesElchElectrode(funct,derxy,ephinp,ephin,econvelnp,ehist);

          return;
        }

        //! set factor F/RT
        void SetFRT()
        {
          const double faraday = DRT::ELEMENTS::ScaTraEleParameterElch::Instance("scatra")->Faraday();
          const double gasconstant = DRT::ELEMENTS::ScaTraEleParameterElch::Instance("scatra")->GasConstant();

          vmelch::frt_ = faraday/(gasconstant*vmthermo::Temp());

          return;
        }
    }; // class ScaTraEleInternalVariableManagerElchElectrodeSTIThermo
  } // namespace ELEMENTS
} // namespace DRT
#endif
