/*--------------------------------------------------------------------------*/
/*!
\file scatra_ele_calc_elch_diffcond_sti_thermo.H

\brief evaluation of scatra elements for thermodynamic diffusion-conduction ion-transport equations

\level 2

<pre>
\maintainer Rui Fang
            fang@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/
            089-289-15251
</pre>
*/
/*--------------------------------------------------------------------------*/
#ifndef SCATRA_ELE_CALC_ELCH_DIFFCOND_STI_THERMO_H
#define SCATRA_ELE_CALC_ELCH_DIFFCOND_STI_THERMO_H

#include "scatra_ele_sti_thermo.H"

#include "../drt_scatra_ele/scatra_ele_calc_elch_diffcond.H"

namespace DRT
{
  namespace ELEMENTS
  {
    // forward declaration
    template<int NSD,int NEN> class ScaTraEleInternalVariableManagerElchDiffCondSTIThermo;

    // class implementation
    template<DRT::Element::DiscretizationType distype>
    class ScaTraEleCalcElchDiffCondSTIThermo : public ScaTraEleCalcElchDiffCond<distype>, public ScaTraEleSTIThermo<distype>
    {
      public:

        //! singleton access method
        static ScaTraEleCalcElchDiffCondSTIThermo<distype>* Instance(
            const int numdofpernode,
            const int numscal,
            const std::string& disname,
            const ScaTraEleCalcElchDiffCondSTIThermo* delete_me = NULL
            );

        //! singleton destruction
        virtual void Done();

      private:

        //! abbreviations
        typedef ScaTraEleCalc<distype> my;
        typedef ScaTraEleCalcElch<distype> myelch;
        typedef ScaTraEleCalcElchDiffCond<distype> mydiffcond;
        typedef ScaTraEleSTIThermo<distype> mythermo;

        //! private constructor for singletons
        ScaTraEleCalcElchDiffCondSTIThermo(const int numdofpernode,const int numscal,const std::string& disname);

        //! destructor
        virtual ~ScaTraEleCalcElchDiffCondSTIThermo(){return;}

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
            const int             iquad = -1     //!< ID of current integration point
            );

        //! calculate element matrix and element right-hand side vector
        void CalcMatAndRhs(
            Epetra_SerialDenseMatrix&     emat,         //!< element matrix
            Epetra_SerialDenseVector&     erhs,         //!< element right-hand side vector
            const int                     k,            //!< index of current scalar
            const double                  fac,          //!< domain integration factor
            const double                  timefacfac,   //!< domain integration factor times time integration factor
            const double                  rhsfac,       //!< domain integration factor times time integration factor for right-hand side vector
            const double                  taufac,       //!< domain integration factor times stabilization parameter
            const double                  timetaufac,   //!< domain integration factor times stabilization parameter times time integration factor
            const double                  rhstaufac,    //!< domain integration factor times stabilization parameter times time integration factor for right-hand side vector
            LINALG::Matrix<my::nen_,1>&   tauderpot,    //!< derivatives of stabilization parameter w.r.t. electric potential
            double&                       rhsint        //!< body force value
            );

        //! fill element matrix with linearizations of discrete scatra residuals w.r.t. thermo dofs
        void SysmatODScatraThermo(
            DRT::Element*               ele,   //!< current element
            Epetra_SerialDenseMatrix&   emat   //!< element matrix
            );

        //! set internal variables for element evaluation
        void SetInternalVariablesForMatAndRHS();

        //! get internal variable manager for thermodynamic diffusion-conduction formulation
        Teuchos::RCP<ScaTraEleInternalVariableManagerElchDiffCondSTIThermo<my::nsd_,my::nen_> > VarManager(){return Teuchos::rcp_static_cast<ScaTraEleInternalVariableManagerElchDiffCondSTIThermo<my::nsd_,my::nen_> >(my::scatravarmanager_);};
    }; // class ScaTraEleCalcElchDiffCondSTIThermo


    //! implementation of ScaTraEleInternalVariableManagerElchDiffCondSTIThermo
    template<int NSD,int NEN>
    class ScaTraEleInternalVariableManagerElchDiffCondSTIThermo : public ScaTraEleInternalVariableManagerElchDiffCond<NSD,NEN>, public ScaTraEleInternalVariableManagerSTIThermo<NSD,NEN>
    {
      public:

        //! abbreviations
        typedef ScaTraEleInternalVariableManagerElch<NSD,NEN> vmelch;
        typedef ScaTraEleInternalVariableManagerElchDiffCond<NSD,NEN> vmdiffcond;
        typedef ScaTraEleInternalVariableManagerSTIThermo<NSD,NEN> vmthermo;

        //! constructor
        ScaTraEleInternalVariableManagerElchDiffCondSTIThermo(
            int numscal,
            const DRT::ELEMENTS::ScaTraEleParameterElch* elchparams,
            const DRT::ELEMENTS::ScaTraEleParameterElchDiffCond* diffcondparams
            ) :
            // call base class constructors
            ScaTraEleInternalVariableManagerElchDiffCond<NSD,NEN>(numscal,elchparams,diffcondparams),
            ScaTraEleInternalVariableManagerSTIThermo<NSD,NEN>()
        {return;};

        //! destructor
        virtual ~ScaTraEleInternalVariableManagerElchDiffCondSTIThermo(){return;};

        //! set internal variables for element evaluation
        void SetInternalVariables(
            const LINALG::Matrix<NEN,1>&                 funct,       //!< shape functions
            const LINALG::Matrix<NSD,NEN>&               derxy,       //!< spatial derivatives of shape functions
            const LINALG::Matrix<NEN,1>&                 etempnp,     //!< nodal temperature values at time t_(n+1) or t_(n+alpha_F)
            const std::vector<LINALG::Matrix<NEN,1> >&   ephinp,      //!< nodal concentration and electric potential values at time t_(n+1) or t_(n+alpha_F)
            const std::vector<LINALG::Matrix<NEN,1> >&   ephin,       //!< nodal concentration and electric potential values at time t_(n)
            const LINALG::Matrix<NSD,NEN>&               econvelnp,   //!< nodal convective velocity values at time t_(n+1) or t_(n+alpha_F)
            const std::vector<LINALG::Matrix<NEN,1> >&   ehist        //!< nodal history values
            )
        {
          // set thermo variables
          vmthermo::SetInternalVariablesSTIThermo(funct,derxy,etempnp);

          // set scatra variables
          // this requires the temperature to be already set
          vmdiffcond::SetInternalVariablesElchDiffCond(funct,derxy,ephinp,ephin,econvelnp,ehist);

          return;
        }

        //! set factor F/RT
        void SetFRT()
        {
          vmelch::frt_ = DRT::ELEMENTS::ScaTraEleParameterElch::Instance("scatra")->Faraday()/(DRT::ELEMENTS::ScaTraEleParameterElch::Instance("scatra")->GasConstant()*vmthermo::Temp());

          return;
        }
    }; // class ScaTraEleInternalVariableManagerElchDiffCondSTIThermo
  } // namespace ELEMENTS
} // namespace DRT
#endif
