/*--------------------------------------------------------------------------*/
/*! \file

\brief evaluation of scatra elements for conservation of mass concentration and electronic charge
within isothermal electrodes

\level 2

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
    template <int NSD, int NEN>
    class ScaTraEleInternalVariableManagerElchElectrode;
    template <DRT::Element::DiscretizationType distype>
    class ScaTraEleUtilsElchElectrode;

    // class implementation
    template <DRT::Element::DiscretizationType distype,
        int probdim = DRT::UTILS::DisTypeToDim<distype>::dim>
    class ScaTraEleCalcElchElectrode : public ScaTraEleCalcElch<distype, probdim>
    {
     public:
      //! singleton access method
      static ScaTraEleCalcElchElectrode<distype, probdim>* Instance(
          const int numdofpernode, const int numscal, const std::string& disname);



      //! evaluate action
      virtual int EvaluateAction(DRT::Element* ele, Teuchos::ParameterList& params,
          DRT::Discretization& discretization, const SCATRA::Action& action,
          DRT::Element::LocationArray& la, Epetra_SerialDenseMatrix& elemat1_epetra,
          Epetra_SerialDenseMatrix& elemat2_epetra, Epetra_SerialDenseVector& elevec1_epetra,
          Epetra_SerialDenseVector& elevec2_epetra, Epetra_SerialDenseVector& elevec3_epetra);

     protected:
      //! protected constructor for singletons
      ScaTraEleCalcElchElectrode(
          const int numdofpernode, const int numscal, const std::string& disname);

      //! abbreviations
      using my = ScaTraEleCalc<distype, probdim>;
      using myelch = ScaTraEleCalcElch<distype, probdim>;

      /*========================================================================*/
      //! @name general framework
      /*========================================================================*/

      //! calculate contributions to element matrix and residual (inside loop over all scalars)
      virtual void CalcMatAndRhs(Epetra_SerialDenseMatrix& emat,  //!< element matrix to calculate
          Epetra_SerialDenseVector& erhs,                         //!< element rhs to calculate+
          const int k,                                            //!< index of current scalar
          const double fac,                                       //!< domain-integration factor
          const double timefacfac,  //!< domain-integration factor times time-integration factor
          const double rhsfac,  //!< time-integration factor for rhs times domain-integration factor
          const double taufac,  //!< tau times domain-integration factor
          const double
              timetaufac,  //!< domain-integration factor times tau times time-integration factor
          const double rhstaufac,  //!< time-integration factor for rhs times tau times
                                   //!< domain-integration factor
          LINALG::Matrix<my::nen_, 1>&
              tauderpot,  //!< derivatives of stabilization parameter w.r.t. electric potential
          double& rhsint  //!< rhs at Gauss point
      );

      //! calculate contributions to element matrix and residual (outside loop over all scalars)
      virtual void CalcMatAndRhsOutsideScalarLoop(
          Epetra_SerialDenseMatrix& emat,  //!< element matrix to calculate
          Epetra_SerialDenseVector& erhs,  //!< element rhs to calculate
          const double fac,                //!< domain-integration factor
          const double timefacfac,  //!< domain-integration factor times time-integration factor
          const double rhsfac  //!< time-integration factor for rhs times domain-integration factor
      );

      //! compute additional flux terms across Dirichlet boundaries
      virtual void CorrectionForFluxAcrossDC(
          DRT::Discretization& discretization,  //!< discretization
          const std::vector<int>& lm,           //!< location vector
          Epetra_SerialDenseMatrix& emat,       //!< element matrix to calculate
          Epetra_SerialDenseVector& erhs        //!< element rhs to calculate
      ){};

      /*========================================================================*/
      //! @name material and related and related functions
      /*========================================================================*/

      //! get material parameters
      virtual void GetMaterialParams(const DRT::Element* ele,  //!< the element we are dealing with
          std::vector<double>& densn,                          //!< density at t_(n)
          std::vector<double>& densnp,  //!< density at t_(n+1) or t_(n+alpha_F)
          std::vector<double>& densam,  //!< density at t_(n+alpha_M)
          double& visc,                 //!< fluid viscosity
          const int iquad = -1          //!< id of current gauss point (default = -1)
      );

      /*========================================================================*/
      //! @name methods for evaluation of individual terms
      /*========================================================================*/

      //! CalcMat: linearizations of diffusion term and Ohmic overpotential w.r.t. structural
      //! displacements
      void CalcDiffODMesh(Epetra_SerialDenseMatrix& emat,  //!< element matrix
          const int k,                                     //!< index of current scalar
          const int ndofpernodemesh,  //!< number of structural degrees of freedom per node
          const double diffcoeff,     //!< diffusion coefficient
          const double fac,           //!< domain-integration factor
          const double rhsfac,  //!< time-integration factor for rhs times domain-integration factor
          const double J,       //!< Jacobian determinant det(dx/ds)
          const LINALG::Matrix<my::nsd_, 1>& gradphi,    //!< gradient of current scalar
          const LINALG::Matrix<my::nsd_, 1>& convelint,  //!< convective velocity
          const LINALG::Matrix<1, my::nsd_ * my::nen_>&
              dJ_dmesh  //!< derivatives of Jacobian determinant det(dx/ds) w.r.t. structural
                        //!< displacements
      );

      //! CalcMat: linearization of diffusion coefficient in diffusion term
      void CalcMatDiffCoeffLin(Epetra_SerialDenseMatrix& emat,  //!< element matrix to be filled
          const int k,                                          //!< index of current scalar
          const double timefacfac,  //!< domain-integration factor times time-integration factor
          const LINALG::Matrix<my::nsd_, 1>& gradphi,  //!< gradient of concentration at GP
          const double scalar  //!< scaling factor for element matrix contributions
      );

      //! CalcMat: potential equation div i with inserted current - ohmic overpotential
      void CalcMatPotEquDiviOhm(Epetra_SerialDenseMatrix& emat,  //!< element matrix to be filled
          const double timefacfac,  //!< domain-integration factor times time-integration factor
          const double invf,        //!< 1/F
          const LINALG::Matrix<my::nsd_, 1>& gradpot,  //!< gradient of potential at GP
          const double scalar  //!< scaling factor for element matrix contributions
      );

      /*!
       * @brief calculate the conservative part of the convective term of the right-hand side
       * vector, due to e.g. deformation of the body
       *
       * @param[in,out] erhs  element rhs vector to be filled
       * @param[in] k         index of current scalar
       * @param[in] rhsfac    time-integration factor for rhs times domain-integration factor
       * @param[in] vdiv      divergence of velocity
       */
      void CalcRhsConservativePartOfConvectiveTerm(
          Epetra_SerialDenseVector& erhs, const int k, const double rhsfac, const double vdiv);

      //! CalcRhs: potential equation div i with inserted current - ohmic overpotential
      void CalcRhsPotEquDiviOhm(Epetra_SerialDenseVector& erhs,  //!< element vector to be filled
          const double rhsfac,  //!< time-integration factor for rhs times domain-integration factor
          const double invf,    //!< 1./F
          const LINALG::Matrix<my::nsd_, 1>& gradpot,  //!< gradient of potential at GP
          const double scalar  //!< scaling factor for element residual contributions
      );

      //! calculate weighted current density
      virtual void CalculateCurrent(LINALG::Matrix<my::nsd_, 1>& q,  //!< flux of species k
          const INPAR::SCATRA::FluxType fluxtype,                    //!< type fo flux
          const double fac                                           //!< integration factor
      );

      //! get utility class supporting element evaluation for electrodes
      DRT::ELEMENTS::ScaTraEleUtilsElchElectrode<distype>* Utils()
      {
        return static_cast<DRT::ELEMENTS::ScaTraEleUtilsElchElectrode<distype>*>(myelch::utils_);
      };

      /*========================================================================*/
      //! @name additional service routines
      /*========================================================================*/

      //! validity check with respect to input parameters, degrees of freedom, number of scalars
      //! etc.
      virtual void CheckElchElementParameter(DRT::Element* ele  //!< current element
      );

      //! get conductivity
      void GetConductivity(const enum INPAR::ELCH::EquPot
                               equpot,  //!< type of closing equation for electric potential
          double& sigma_all,            //!< conductivity of electrolyte solution
          std::vector<double>&
              sigma,    //!< conductivity or a single ion + overall electrolyte solution
          bool effCond  //!< flag if effective conductivity should be calculated
      );

      //! calculate electrode state of charge and C rate
      virtual void CalculateElectrodeSOCAndCRate(
          const DRT::Element* const& ele,             //!< the element we are dealing with
          const DRT::Discretization& discretization,  //!< discretization
          DRT::Element::LocationArray& la,            //!< location array
          Epetra_SerialDenseVector& scalars  //!< result vector for scalar integrals to be computed
      );

      //! calculate mean concentration of micro discretization
      //! \param ele              the element we are dealing with
      //! \param discretization   discretization
      //! \param la               location array
      //! \param conc             result vector for scalar integrals to be computed
      virtual void CalculateMeanElectrodeConcentration(const DRT::Element* const& ele,
          const DRT::Discretization& discretization, DRT::Element::LocationArray& la,
          Epetra_SerialDenseVector& conc);

      //! calculate weighted mass flux (no reactive flux so far)
      virtual void CalculateFlux(LINALG::Matrix<my::nsd_, 1>& q,  //!< flux of species k
          const INPAR::SCATRA::FluxType fluxtype,                 //!< type fo flux
          const int k                                             //!< index of current scalar
      );

      //! calculate error of numerical solution with respect to analytical solution
      virtual void CalErrorComparedToAnalytSolution(const DRT::Element* ele,  //!< element
          Teuchos::ParameterList& params,                                     //!< parameter list
          Epetra_SerialDenseVector& errors  //!< vector containing L2 and H1 error norms
      );

      //! set internal variables for electrodes
      virtual void SetInternalVariablesForMatAndRHS();

      //! get diffusion manager for electrodes
      Teuchos::RCP<ScaTraEleDiffManagerElchElectrode> DiffManager()
      {
        return Teuchos::rcp_static_cast<ScaTraEleDiffManagerElchElectrode>(my::diffmanager_);
      };

     private:
      //! get internal variable manager for electrodes
      Teuchos::RCP<ScaTraEleInternalVariableManagerElchElectrode<my::nsd_, my::nen_>> VarManager()
      {
        return Teuchos::rcp_static_cast<
            ScaTraEleInternalVariableManagerElchElectrode<my::nsd_, my::nen_>>(
            my::scatravarmanager_);
      };
    };


    //! ScaTraEleDiffManagerElchElectrode implementation
    /*!
      This class is derived from the standard ScaTraEleDiffManager and keeps all electrode-specific
      transport parameters.
    */
    class ScaTraEleDiffManagerElchElectrode : public ScaTraEleDiffManagerElch
    {
     public:
      ScaTraEleDiffManagerElchElectrode(int numscal)
          : ScaTraEleDiffManagerElch(numscal),
            concderivdiff_(numscal, std::vector<double>(numscal, 0.0)),
            tempderivdiff_(numscal, std::vector<double>(numscal, 0.0)),
            cond_(0.0),
            concderivcond_(numscal, 0.0),
            tempderivcond_(numscal, 0.0)
      {
      }

      ~ScaTraEleDiffManagerElchElectrode() override = default;
      ;

      /*========================================================================*/
      //! @name access methods
      /*========================================================================*/

      //! Set derivative of diffusion coefficients with respect to concentrations
      void SetConcDerivIsoDiffCoef(const double concderivdiff, const int k, const int iscal)
      {
        (concderivdiff_[k])[iscal] = concderivdiff;
      };

      //! Set derivative of diffusion coefficient with respect to concentrations
      void SetTempDerivIsoDiffCoef(const double tempderivdiff, const int k, const int iscal)
      {
        (tempderivdiff_[k])[iscal] = tempderivdiff;
      };

      //! Access routine for derivative of diffusion coefficients with respect to concentrations
      double GetConcDerivIsoDiffCoef(const int k, const int iscal)
      {
        return (concderivdiff_[k])[iscal];
      };

      //! Access routine for derivative of diffusion coefficients with respect to concentrations
      double GetTempDerivIsoDiffCoef(const int k, const int iscal)
      {
        return (tempderivdiff_[k])[iscal];
      };

      //! Set conductivity of the electrolyte solution and electrode
      void SetCond(const double cond) { cond_ = cond; };

      //! Access routine for conductivity of the electrolyte solution and electrode
      double GetCond() { return cond_; };

      //! Set derivative of the conductivity with respect to concentrations
      void SetConcDerivCond(const double concderivcond, const int k)
      {
        concderivcond_[k] = concderivcond;
      };

      //! Access routine for derivative of the conductivity with respect to concentrations
      double GetConcDerivCond(const int k) { return concderivcond_[k]; };

      //! Set derivative of the conductivity with respect to temperature
      void SetTempDerivCond(const double tempderivcond, const int k)
      {
        tempderivcond_[k] = tempderivcond;
      };

      //! Access routine for derivative of the conductivity with respect to temperature
      double GetTempDerivCond(const int k) { return tempderivcond_[k]; };

      /*========================================================================*/
      //! @name output
      /*========================================================================*/

      //! print transport parameters to screen
      virtual void OutputTransportParams(const int numscal)
      {
        for (int k = 0; k < numscal; ++k)
          std::cout << "diffusion coefficient " << k << ":   " << diff_[k] << std::endl;

        for (int k = 0; k < numscal; ++k)
        {
          for (int iscal = 0; iscal < numscal; ++iscal)
          {
            std::cout << "derivative of diffusion coefficient (" << k << "," << iscal
                      << "):  " << (concderivdiff_[k])[iscal] << std::endl;
            std::cout << "derivative of diffusion coefficient w.r.t. temperature (" << k << ","
                      << iscal << "):  " << (tempderivdiff_[k])[iscal] << std::endl;
          }
        }
        std::cout << std::endl;

        std::cout << "conductivity:   " << cond_ << std::endl;

        for (int k = 0; k < numscal; ++k)
        {
          std::cout << "derivative of conductivity " << k
                    << " w.r.t. concentration:   " << concderivcond_[k] << std::endl;
          std::cout << "derivative of conductivity " << k
                    << " w.r.t. temperature:   " << tempderivcond_[k] << std::endl;
        }
        std::cout << std::endl;
      };

     protected:
      /*========================================================================*/
      //! @name transport parameter
      /*========================================================================*/

      //! derivative of diffusion coefficients with respect to concentrations
      std::vector<std::vector<double>> concderivdiff_;

      //! derivative of diffusion coefficients with respect to temperature
      std::vector<std::vector<double>> tempderivdiff_;

      //! conductivity of the electrolyte solution
      double cond_;

      //! derivative of the conductivity with respect to concentrations
      std::vector<double> concderivcond_;

      //! derivative of the conductivity with respect to temperature
      std::vector<double> tempderivcond_;
    };

    //! ScaTraEleInternalVariableManagerElchElectrode implementation
    /*!
      This class keeps all internal electrode-specific variables.
    */
    template <int NSD, int NEN>
    class ScaTraEleInternalVariableManagerElchElectrode
        : public ScaTraEleInternalVariableManagerElch<NSD, NEN>
    {
     public:
      using vm = ScaTraEleInternalVariableManager<NSD, NEN>;
      using vmelch = ScaTraEleInternalVariableManagerElch<NSD, NEN>;

      //! constructor
      ScaTraEleInternalVariableManagerElchElectrode(
          int numscal, const DRT::ELEMENTS::ScaTraEleParameterElch* elchpara)
          : ScaTraEleInternalVariableManagerElch<NSD, NEN>(numscal, elchpara),
            invf_(1. / vmelch::parameters_->Faraday())
      {
      }

      //! destructor
      ~ScaTraEleInternalVariableManagerElchElectrode() override = default;

      //! compute and set internal electrode-specific variables
      void SetInternalVariablesElchElectrode(
          const LINALG::Matrix<NEN, 1>& funct,  //!< array for shape functions
          const LINALG::Matrix<NSD, NEN>&
              derxy,  //!< global derivatives of shape functions w.r.t x,y,z
          const std::vector<LINALG::Matrix<NEN, 1>>&
              ephinp,  //!< nodal state variables at t_(n+1) or t_(n+alpha_F)
          const std::vector<LINALG::Matrix<NEN, 1>>& ephin,  //!< nodal state variables at t_(n)
          const LINALG::Matrix<NSD, NEN>&
              econvelnp,  //!< nodal convective velocity values at t_(n+1) or t_(n+alpha_F)
          const std::vector<LINALG::Matrix<NEN, 1>>&
              ehist  //!< history vector of transported scalars
      )
      {
        // set internal variables in base variable manager
        vmelch::SetInternalVariablesElch(funct, derxy, ephinp, ephin, econvelnp, ehist);
      };

      //! compute and set internal electrode-specific variables for evaluation of SOC and Crate
      void SetInternalVariablesElchElectrodeSOCAndCRate(
          const LINALG::Matrix<NEN, 1>& funct,  //!< array for shape functions
          const LINALG::Matrix<NSD, NEN>&
              derxy,  //!< global derivatives of shape functions w.r.t x,y,z
          const std::vector<LINALG::Matrix<NEN, 1>>&
              ephinp,  //!< nodal state variables at t_(n+1) or t_(n+alpha_F)
          const std::vector<LINALG::Matrix<NEN, 1>>& ephin,  //!< nodal state variables at t_(n)
          const LINALG::Matrix<NSD, NEN>&
              econvelnp,  //!< nodal convective velocity values at t_(n+1) or t_(n+alpha_F)
          const std::vector<LINALG::Matrix<NEN, 1>>&
              ehist  //!< history vector of transported scalars
      )
      {
        // FRT not needed for calculation of SOC and C Rate. Temperature is not extracted in case of
        // non constant temperatures,
        const bool do_setfrt = false;
        // set internal variables in base variable manager
        vmelch::SetInternalVariablesElch(funct, derxy, ephinp, ephin, econvelnp, ehist, do_setfrt);
      };

      //! return constant parameter 1./F
      double InvF() { return invf_; };

     protected:
      //! constant parameter 1./F
      const double invf_;
    };  // class ScaTraEleInternalVariableManagerElchElectrode
  }     // namespace ELEMENTS
}  // namespace DRT
#endif
