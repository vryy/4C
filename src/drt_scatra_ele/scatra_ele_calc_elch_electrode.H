/*--------------------------------------------------------------------------*/
/*!
\file scatra_ele_calc_elch_electrode.H

\brief evaluation of scatra elements for conservation of mass concentration and electronic charge within isothermal electrodes

\level 2

<pre>
\maintainer Rui Fang
            fang@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/
            089-289-15251
</pre>
*/
/*--------------------------------------------------------------------------*/
#ifndef SCATRA_ELE_CALC_ELCH_ELECTRODE_H
#define SCATRA_ELE_CALC_ELCH_ELECTRODE_H

#include "scatra_ele_calc_elch.H"
#include "scatra_ele_utils_elch_electrode.H"

namespace DRT
{
  namespace ELEMENTS
  {
    // forward declarations
    class ScaTraEleDiffManagerElchElectrode;
    template<int NSD, int NEN> class ScaTraEleInternalVariableManagerElchElectrode;
    template<DRT::Element::DiscretizationType distype> class ScaTraEleUtilsElchElectrode;

    // class implementation
    template<DRT::Element::DiscretizationType distype>
    class ScaTraEleCalcElchElectrode : public ScaTraEleCalcElch<distype>
    {
      public:

        //! singleton access method
        static ScaTraEleCalcElchElectrode<distype>* Instance(
            const int numdofpernode,
            const int numscal,
            const std::string& disname,
            const ScaTraEleCalcElchElectrode* delete_me = NULL
            );

        //! singleton destruction
        virtual void Done();

        //! evaluate action
        virtual int EvaluateAction(
            DRT::Element*                 ele,
            Teuchos::ParameterList&       params,
            DRT::Discretization&          discretization,
            const SCATRA::Action&         action,
            DRT::Element::LocationArray&  la,
            Epetra_SerialDenseMatrix&     elemat1_epetra,
            Epetra_SerialDenseMatrix&     elemat2_epetra,
            Epetra_SerialDenseVector&     elevec1_epetra,
            Epetra_SerialDenseVector&     elevec2_epetra,
            Epetra_SerialDenseVector&     elevec3_epetra
            );

      protected:

        //! protected constructor for singletons
        ScaTraEleCalcElchElectrode(const int numdofpernode,const int numscal,const std::string& disname);

        //! destructor
        virtual ~ScaTraEleCalcElchElectrode(){return;}

        //! abbreviations
        typedef ScaTraEleCalc<distype> my;
        typedef ScaTraEleCalcElch<distype> myelch;

        /*========================================================================*/
        //! @name general framework
        /*========================================================================*/

        //! calculate contributions to element matrix and residual (inside loop over all scalars)
        virtual void CalcMatAndRhs(
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

        //! calculate contributions to element matrix and residual (outside loop over all scalars)
        virtual void CalcMatAndRhsOutsideScalarLoop(
            Epetra_SerialDenseMatrix&   emat,         //!< element matrix to calculate
            Epetra_SerialDenseVector&   erhs,         //!< element rhs to calculate
            const double                fac,          //!< domain-integration factor
            const double                timefacfac,   //!< domain-integration factor times time-integration factor
            const double                rhsfac        //!< time-integration factor for rhs times domain-integration factor
            );

        //! compute additional flux terms across Dirichlet boundaries
        virtual void CorrectionForFluxAcrossDC(
          DRT::Discretization&        discretization,   //!< discretization
          const std::vector<int>&     lm,               //!< location vector
          Epetra_SerialDenseMatrix&   emat,             //!< element matrix to calculate
          Epetra_SerialDenseVector&   erhs              //!< element rhs to calculate
        ){return;};

        /*========================================================================*/
        //! @name material and related and related functions
        /*========================================================================*/

        //! get material parameters
        virtual void GetMaterialParams(
          const DRT::Element*                  ele,           //!< the element we are dealing with
          double&                              densn,         //!< density at t_(n)
          double&                              densnp,        //!< density at t_(n+1) or t_(n+alpha_F)
          double&                              densam,        //!< density at t_(n+alpha_M)
          double&                              visc,          //!< fluid viscosity
          const int                            iquad=-1       //!< id of current gauss point (default = -1)
        );

        /*========================================================================*/
        //! @name methods for evaluation of individual terms
        /*========================================================================*/

        //! CalcMat: linearization of diffusion coefficient in diffusion term
        void CalcMatDiffCoeffLin(
          Epetra_SerialDenseMatrix&                          emat,         //!< element matrix to be filled
          const int                                          k,            //!< index of current scalar
          const double                                       timefacfac,   //!< domain-integration factor times time-integration factor
          const LINALG::Matrix<my::nsd_,1>&                  gradphi,      //!< gradient of concentration at GP
          const double                                       scalar        //!< scaling factor for element matrix contributions
        );

        //! CalcMat: potential equation div i with inserted current - ohmic overpotential
        void CalcMatPotEquDiviOhm(
          Epetra_SerialDenseMatrix&                          emat,         //!< element matrix to be filled
          const double                                       timefacfac,   //!< domain-integration factor times time-integration factor
          const double                                       invf,         //!< 1/F
          const LINALG::Matrix<my::nsd_,1>&                  gradpot,      //!< gradient of potenial at GP
          const double                                       scalar        //!< scaling factor for element matrix contributions
        );

        //! CalcRhs: potential equation div i with inserted current - ohmic overpotential
        void CalcRhsPotEquDiviOhm(
          Epetra_SerialDenseVector&                          erhs,      //!< element vector to be filled
          const double                                       rhsfac,    //!< time-integration factor for rhs times domain-integration factor
          const double                                       invf,      //!< 1./F
          const LINALG::Matrix<my::nsd_,1>&                  gradpot,   //!< gradient of potenial at GP
          const double                                       scalar     //!< scaling factor for element residual contributions
        );

        //! calculate weighted current density
        virtual void CalculateCurrent(
            LINALG::Matrix<my::nsd_,1>&     q,          //!< flux of species k
            const INPAR::SCATRA::FluxType   fluxtype,   //!< type fo flux
            const double                    fac         //!< integration factor
            );

        /*========================================================================*/
        //! @name additional service routines
        /*========================================================================*/

        //! validity check with respect to input parameters, degrees of freedom, number of scalars etc.
        virtual void CheckElchElementParameter(
            DRT::Element*              ele  //!< current element
        );

        //! get conductivity
        void GetConductivity(
          const enum INPAR::ELCH::EquPot   equpot,      //!< type of closing equation for electric potential
          double&                          sigma_all,   //!< conductivity of electrolyte solution
          std::vector<double>&             sigma,       //!< conductivity or a single ion + overall electrolyte solution
          bool                             effCond      //!< flag if effective conductivity should be calculated
        );

        //! calculate electrode state of charge
        void CalculateElectrodeSOC(
            const DRT::Element*               ele,              //!< the element we are dealing with
            Teuchos::ParameterList&           params,           //!< parameter list
            DRT::Discretization&              discretization,   //!< discretization
            const std::vector<int>&           lm,               //!< location vector
            Epetra_SerialDenseVector&         scalars           //!< result vector for scalar integrals to be computed
        );

        //! calculate weighted mass flux (no reactive flux so far)
        virtual void CalculateFlux(
            LINALG::Matrix<my::nsd_,1>&     q,          //!< flux of species k
            const INPAR::SCATRA::FluxType   fluxtype,   //!< type fo flux
            const int                       k           //!< index of current scalar
        );

        //! calculate error of numerical solution with respect to analytical solution
        void CalErrorComparedToAnalytSolution(
          const DRT::Element*                   ele,        //!< the element we are dealing with
          Teuchos::ParameterList&               params,     //!< parameter list
          Epetra_SerialDenseVector&             errors      //!< vector containing L2-error norm
        ){return;};

        //! set internal variables for electrodes
        virtual void SetInternalVariablesForMatAndRHS();

        //! get diffusion manager for electrodes
        Teuchos::RCP<ScaTraEleDiffManagerElchElectrode> DiffManager(){return Teuchos::rcp_static_cast<ScaTraEleDiffManagerElchElectrode>(my::diffmanager_);};

      private:

        //! get internal variable manager for electrodes
        Teuchos::RCP<ScaTraEleInternalVariableManagerElchElectrode<my::nsd_,my::nen_> > VarManager(){return Teuchos::rcp_static_cast<ScaTraEleInternalVariableManagerElchElectrode<my::nsd_,my::nen_> >(my::scatravarmanager_);};

        //! get utility class supporting element evaluation for electrodes
        DRT::ELEMENTS::ScaTraEleUtilsElchElectrode<distype>* Utils(){return static_cast<DRT::ELEMENTS::ScaTraEleUtilsElchElectrode<distype>*>(myelch::utils_);};
    };


    //! ScaTraEleDiffManagerElchElectrode implementation
    /*!
      This class is derived from the standard ScaTraEleDiffManager and keeps all electrode-specific transport parameters.
    */
    class ScaTraEleDiffManagerElchElectrode : public ScaTraEleDiffManagerElch
    {
      public:

        ScaTraEleDiffManagerElchElectrode(int numscal) :
          ScaTraEleDiffManagerElch(numscal),
          derivdiff_(numscal,std::vector<double>(numscal,0.0)),
          cond_(0.0),
          derivcond_(numscal,0.0)
        {return;}

        virtual ~ScaTraEleDiffManagerElchElectrode(){return;};

        /*========================================================================*/
        //! @name access methods
        /*========================================================================*/

        //! Set derivative of diffusion coefficients with respect to concentrations
        void SetDerivIsoDiffCoef(const double derivdiff,const int k,const int iscal)
        {
          (derivdiff_[k])[iscal] = derivdiff;
          return;
        };

        //! Access routine for derivative of diffusion coefficients with respect to concentrations
        double GetDerivIsoDiffCoef(const int k,const int iscal){return (derivdiff_[k])[iscal];};

        //! Set conductivity of the electrolyte solution
        void SetCond(const double cond)
        {
          cond_ = cond;
          return;
        };

        //! Access routine for conductivity of the electrolyte solution
        double GetCond(){return cond_;};

        //! Set derivative of the conductivity with respect to concentrations
        void SetDerivCond(const double derivcond,const int k)
        {
          derivcond_[k] = derivcond;
          return;
        };

        //! Access routine for derivative of the conductivity with respect to concentrations
        double GetDerivCond(const int k){return derivcond_[k];};

        /*========================================================================*/
        //! @name output
        /*========================================================================*/

        //! print transport parameters to screen
        virtual void OutputTransportParams(const int numscal)
        {
          for(int k=0;k<numscal;++k)
            std::cout << "diffusion coefficient " << k << ":   " << diff_[k] << std::endl;

          for(int k=0;k<numscal;++k)
          {
            for(int iscal=0;iscal<numscal;++iscal)
            {
              std::cout << "derivative of diffusion coefficient (" << k << "," << iscal << "):  " << (derivdiff_[k])[iscal] << std::endl;
            }
          }
          std::cout << std::endl;

          std::cout << "conductivity:   " << cond_ << std::endl;

          for(int k=0;k<numscal;++k)
            std::cout << "derivative of conductivity  " << k << ":   " << derivcond_[k] << std::endl;
          std::cout << std::endl;

          return;
        };

      protected:

        /*========================================================================*/
        //! @name transport parameter
        /*========================================================================*/

        //! derivative of diffusion coefficients with respect to concentrations
        std::vector< std::vector<double> > derivdiff_;

        //! conductivity of the electrolyte solution
        double cond_;

        //! derivative of the conductivity with respect to concentrations
        std::vector<double> derivcond_;
    };


    //! ScaTraEleInternalVariableManagerElchElectrode implementation
    /*!
      This class keeps all internal electrode-specific variables.
    */
    template<int NSD, int NEN>
    class ScaTraEleInternalVariableManagerElchElectrode : public ScaTraEleInternalVariableManagerElch<NSD,NEN>
    {
      public:

        typedef ScaTraEleInternalVariableManager<NSD,NEN> vm;
        typedef ScaTraEleInternalVariableManagerElch<NSD,NEN> vmelch;

        //! constructor
        ScaTraEleInternalVariableManagerElchElectrode(int numscal,const DRT::ELEMENTS::ScaTraEleParameterElch* elchpara):
          ScaTraEleInternalVariableManagerElch<NSD,NEN>(numscal,elchpara),
          invf_(1./vmelch::parameters_->Faraday())
        {return;}

        //! destructor
        virtual ~ScaTraEleInternalVariableManagerElchElectrode(){return;};

        //! compute and set internal electrode-specific variables
        void SetInternalVariablesElchElectrode(
            const LINALG::Matrix<NEN,1>&                 funct,       //!< array for shape functions
            const LINALG::Matrix<NSD,NEN>&               derxy,       //!< global derivatives of shape functions w.r.t x,y,z
            const std::vector<LINALG::Matrix<NEN,1> >&   ephinp,      //!< nodal state variables at t_(n+1) or t_(n+alpha_F)
            const std::vector<LINALG::Matrix<NEN,1> >&   ephin,       //!< nodal state variables at t_(n)
            const LINALG::Matrix<NSD,NEN>&               econvelnp,   //!< nodal convective velocity values at t_(n+1) or t_(n+alpha_F)
            const std::vector<LINALG::Matrix<NEN,1> >&   ehist        //!< history vector of transported scalars
            )
        {
          // set internal variables in base variable manager
          vmelch::SetInternalVariablesElch(funct,derxy,ephinp,ephin,econvelnp,ehist);

          return;
        };

        //! return constant parameter 1./F
        double InvF(){return invf_;};

      protected:

        //! constant parameter 1./F
        const double invf_;
    }; // class ScaTraEleInternalVariableManagerElchElectrode
  } // namespace ELEMENTS
} // namespace DRT
#endif
