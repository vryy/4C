/*--------------------------------------------------------------------------*/
/*!
\file scatra_ele_calc_variational.H

\brief evaluation of ScaTra elements for chemical diffusion under a variational framework
\Equations : Evolution eq: Conservation of mass
\      Constitutive Eq: Mass Action Law

\level 2

<pre>
\maintainer Jorge De Anda Salazar
            DeAnda@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/
            089-289-15251
</pre>
*/
/*--------------------------------------------------------------------------*/

#ifndef SCATRA_ELE_CALC_VARIATIONAL_H
#define SCATRA_ELE_CALC_VARIATIONAL_H

#include "../drt_fem_general/drt_utils_local_connectivity_matrices.H"  //Check if really needed

#include "scatra_ele_calc.H"
#include "scatra_ele_parameter_elch.H"
#include "../drt_mat/scatra_mat_var_chemdiffusion.H"


#include <Epetra_IntSerialDenseVector.h>  //Check if really needed

#include "scatra_ele_action.H"      //Check if really needed
#include "scatra_ele_calc_utils.H"  //Check if really needed
#include "scatra_ele_interface.H"   //Check if really needed

#include "../drt_fluid_ele/fluid_ele.H"  //Check if really needed

namespace DRT
{
  namespace ELEMENTS
  {
    // forward declarations
    template <int NSD, int NEN>
    class ScaTraEleDiffManagerVar;
    template <int NSD, int NEN>
    class ScaTraEleInternalVariableManagerVar;
    template <DRT::Element::DiscretizationType distype>
    class ScaTraEleUtilsElch;

    // class implementation
    template <DRT::Element::DiscretizationType distype>
    class ScaTraEleCalVariational : public ScaTraEleCalc<distype>
    {
     protected:
      /// protected constructor, since we are a singleton
      /// this constructor is called from a derived class
      /// -> therefore, it has to be protected instead of private
      ScaTraEleCalVariational(
          const int numdofpernode, const int numscal, const std::string& disname);

      typedef ScaTraEleCalc<distype> my;
      //  typedef ScaTraEleCalc<distype> myVar;

     public:
      virtual ~ScaTraEleCalVariational() {}

      /// called upon destruction
      virtual void Done();

      //! singleton access method
      static ScaTraEleCalVariational<distype>* Instance(const int numdofpernode, const int numscal,
          const std::string& disname, const ScaTraEleCalVariational* delete_me = NULL);

      //! evaluate the element
      int Evaluate(DRT::Element* ele, Teuchos::ParameterList& params,
          DRT::Discretization& discretization, DRT::Element::LocationArray& la,
          Epetra_SerialDenseMatrix& elemat1_epetra, Epetra_SerialDenseMatrix& elemat2_epetra,
          Epetra_SerialDenseVector& elevec1_epetra, Epetra_SerialDenseVector& elevec2_epetra,
          Epetra_SerialDenseVector& elevec3_epetra);

      //! get the element initial conditions
      void ExtractInitialNodeValues(DRT::Element* ele, Teuchos::ParameterList& params,
          DRT::Discretization& discretization, DRT::Element::LocationArray& la);

     protected:
      /*========================================================================*/
      //! @name general framework
      /*========================================================================*/

      //! Prepare everything what is needed in CallMatAndRhs() to calculate the sysmat and the rhs
      void Sysmat(DRT::Element* ele,           //!< the element we are dealing with
          Epetra_SerialDenseMatrix& emat,      //!< element matrix to calculate
          Epetra_SerialDenseVector& erhs,      //!< element rhs to calculate
          Epetra_SerialDenseVector& subgrdiff  //!< subgrid-diff.-scaling vector
      );

      //! calculate contributions to matrix and rhs (inside of loop over all scalars)
      void CalcMatAndRhs(Epetra_SerialDenseMatrix& emat,  //!< element matrix to calculate
          Epetra_SerialDenseVector& erhs,                 //!< element rhs to calculate+
          const int k,                                    //!< index of current scalar
          const double weight,      //!< quadrature weight at integration point
          const double fac,         //!< domain-integration factor
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

      //! calculation of generic element RHS
      void CalcRhsGeneric(Epetra_SerialDenseVector& erhs,  //!< element vector to be filled
          const int k,                                     //!< index of current scalar
          const double& discrFact,  //!< discrete factors(time &temporal discretization, weights,
                                    //!< Jacobian, etc...)
          const double& kernel,     //!< integral kernel (excluding shape functions)
          const LINALG::Matrix<my::nen_, 1>& sfunct  //!< first shape function
          ) const;

      //! calculation of generic element RHS
      void CalcRhsGeneric(Epetra_SerialDenseVector& erhs,  //!< element vector to be filled
          const int k,                                     //!< index of current scalar
          const double& discrFact,  //!< discrete factors(time &temporal discretization, weights,
                                    //!< Jacobian, etc...)
          const LINALG::Matrix<my::nsd_, 1>&
              kernel,  //!< integral kernel (excluding shape functions)
          const LINALG::Matrix<my::nsd_, my::nen_>& sfunct  //!< first shape function
          ) const;

      //! calculation of generic element matrix
      void CalcMatGeneric(Epetra_SerialDenseMatrix& emat,  //!< element matrix to be filled
          const int& row,                                  //!< starting row position index
          const int& col,                                  //!< starting column position index
          const double& discrFact,  //!< discrete factors(time &temporal discretization, weights,
                                    //!< Jacobian, etc...)
          const double& kernel,     //!< integral kernel (excluding shape functions)
          const LINALG::Matrix<my::nen_, 1>& sfunct,  //!< first shape function
          const LINALG::Matrix<my::nen_, 1>& tfunct   //!< second shape function
          ) const;

      //! calculation of generic element matrix
      void CalcMatGeneric(Epetra_SerialDenseMatrix& emat,  //!< element matrix to be filled
          const int& row,                                  //!< starting row position index
          const int& col,                                  //!< starting column position index
          const double& discrFact,  //!< discrete factors(time &temporal discretization, weights,
                                    //!< Jacobian, etc...)
          const LINALG::Matrix<my::nsd_, my::nsd_>&
              kernel,  //!< integral kernel (excluding shape functions)
          const LINALG::Matrix<my::nsd_, my::nen_>& sfunct,  //!< first shape function
          const LINALG::Matrix<my::nsd_, my::nen_>& tfunct   //!< second shape function
          ) const;

      //! get material parameters
      virtual void GetMaterialParams(const DRT::Element* ele,  //!< the element we are dealing with
          std::vector<double>& densn,                          //!< density at t_(n)
          std::vector<double>& densnp,  //!< density at t_(n+1) or t_(n+alpha_F)
          std::vector<double>& densam,  //!< density at t_(n+alpha_M)
          double& visc,                 //!< fluid viscosity
          const int iquad = -1          //!< id of current gauss point (default = -1)
      );

      //! evaluate material
      virtual void Materials(
          const Teuchos::RCP<const MAT::Material> material,  //!< pointer to current material
          const int k,                                       //!< id of current scalar
          double& densn,                                     //!< density at t_(n)
          double& densnp,       //!< density at t_(n+1) or t_(n+alpha_F)
          double& densam,       //!< density at t_(n+alpha_M)
          double& visc,         //!< fluid viscosity
          const int iquad = -1  //!< id of current gauss point (default = -1)
      );


      //! material ScaTra
      virtual void MatScaTra_Var_ChemDiff(
          const Teuchos::RCP<const MAT::Material> material,  //!< pointer to current material
          const int k,                                       //!< id of current scalar
          double& densn,                                     //!< density at t_(n)
          double& densnp,       //!< density at t_(n+1) or t_(n+alpha_F)
          double& densam,       //!< density at t_(n+alpha_M)
          double& visc,         //!< fluid viscosity
          const int iquad = -1  //!< id of current gauss point (default = -1)
      );

      //! get material fields interpolated
      virtual void GetMaterialInterpolatedFields(
          const DRT::Element* ele,  //!< the element we are dealing with
          const std::vector<LINALG::Matrix<my::nen_, 1>>&
              ephinp  //!< nodal state variables at t_(n+1) or t_(n+alpha_F)
      );

      //! set internal variables
      void SetInternalVariablesForMatAndRHS();

      //! finite difference check for debugging purposes
      void FDCheck(DRT::Element* ele,          //!< the element we are dealing with
          Epetra_SerialDenseMatrix& emat,      //!< element matrix to calculate
          Epetra_SerialDenseVector& erhs,      //!< element rhs to calculate
          Epetra_SerialDenseVector& subgrdiff  //!< subgrid-diff.-scaling vector
      );

      /*========================================================================*/
      //! @name additional service routines
      /*========================================================================*/
      //! calculate error of numerical solution with respect to analytical solution
      void CalErrorComparedToAnalytSolution(
          const DRT::Element* ele,          //!< the element we are dealing with
          Teuchos::ParameterList& params,   //!< parameter list
          Epetra_SerialDenseVector& errors  //!< vector containing L2-error norm
      );

      /*========================================================================*/
      //! @name scalar degrees of freedom and related
      /*========================================================================*/

      //! state variables at t_(0)
      std::vector<LINALG::Matrix<my::nen_, 1>> ephi0_;

      //! get var diffusion manager
      Teuchos::RCP<ScaTraEleDiffManagerVar<my::nsd_, my::nen_>> DiffManager()
      {
        return Teuchos::rcp_static_cast<ScaTraEleDiffManagerVar<my::nsd_, my::nen_>>(
            my::diffmanager_);
      };

      // Static since we only need to read it once
      bool IsSemImplicitFunctional_;  // TODO make static

     private:
      //! get internal variable manager for electrodes
      Teuchos::RCP<ScaTraEleInternalVariableManagerVar<my::nsd_, my::nen_>> VarManager()
      {
        return Teuchos::rcp_static_cast<ScaTraEleInternalVariableManagerVar<my::nsd_, my::nen_>>(
            my::scatravarmanager_);
      };
    };


    /// ScaTraEleDiffManagerVar implementation
    /*!
      This class keeps all transport parameter needed for the evaluation of an element under a
      variational formulation for chemical diffusion. The ScaTraEleDiffManagerVar is derived from
      the standard ScaTraEleDiffManager.
    */
    template <int NSD, int NEN>
    class ScaTraEleDiffManagerVar : public ScaTraEleDiffManager
    {
     public:
      ScaTraEleDiffManagerVar(int numscal)
          : ScaTraEleDiffManager(numscal)  //,
            ,
            mobility_(numscal, 0.0),
            refMu_(numscal, 0.0),
            refC_(numscal, 0.0),
            rt_(0.0),
            IntEnergy_(numscal, 0.0),
            IntEnergy1deriv_(numscal, 0.0),
            IntEnergy2deriv_(numscal, 0.0),
            DissipationPot_(numscal, 0.0),
            DissipationPot1deriv_(numscal, LINALG::Matrix<NSD, 1>(true)),
            DissipationPot2deriv_(numscal, LINALG::Matrix<NSD, NSD>(true))
      {
        return;
      }

      ~ScaTraEleDiffManagerVar(){};

      //! Set mobility of the single ionic species k
      void SetMobility(const double diffusivity, const double rt, const int k)
      {
        mobility_[k] = diffusivity / rt;
        return;
      };

      //! Set reference chemical potential of the single ionic species k
      void SetRefMu(const double refMu, const int k)
      {
        refMu_[k] = refMu;
        return;
      };

      //! Set reference concentration of the single ionic species k
      void SetRefC(const double refC, const int k)
      {
        refC_[k] = refC;
        return;
      };

      //! Set gas constant times temperature
      void SetRT(const double temp, const double gascon)
      {
        rt_ = temp * gascon;
        return;
      };

      //! Set Internal energy
      void SetInternalEnergy(const Teuchos::RCP<const MAT::ScatraMatVarChemDiffusion>
                                 material,  //!< pointer to current material
          const double concentration,       //!< Concentration at t_(n+1)
          const double refMu,               //!< Reference chemical potential
          const double refC,                //!< Reference concentration
          const double rt,                  //!< RT factor
          const int k                       //!< species number
      )
      {
        IntEnergy_[k] = material->ComputeInternalEnergy(concentration, refMu, refC, rt, 0);
        return;
      }

      //! Set Internal energy 1st derivative
      void SetInternalEnergy1deriv(const Teuchos::RCP<const MAT::ScatraMatVarChemDiffusion>
                                       material,  //!< pointer to current material
          const double concentration,             //!< Concentration at t_(n+1)
          const double refMu,                     //!< Reference chemical potential
          const double refC,                      //!< Reference concentration
          const double rt,                        //!< RT factor
          const int k                             //!< species number
      )
      {
        IntEnergy1deriv_[k] = material->ComputeInternalEnergy(concentration, refMu, refC, rt, 1);
        return;
      }

      //! Replace Internal energy 1st derivative
      /*  This is needed because sometimes we obtain local negative values for the concentration due
         to the interpolation, which gives error computing the internal energy due to the log... In
         this case the opted solution is to interpolate not the concentration but the internal
         energy itself and we use this function to do that change)*/
      void ChangeInternalEnergy1deriv(
          const double NewIntEnergy,  //!< New Internal energy (1derv) values
          const int k                 //!< species number
      )
      {
        IntEnergy1deriv_[k] = NewIntEnergy;
        return;
      }

      //! Set Internal energy 2nd derivative
      void SetInternalEnergy2deriv(const Teuchos::RCP<const MAT::ScatraMatVarChemDiffusion>
                                       material,  //!< pointer to current material
          const double concentration,             //!< Concentration at t_(n+1)
          const double refMu,                     //!< Reference chemical potential
          const double refC,                      //!< Reference concentration
          const double rt,                        //!< RT factor
          const int k                             //!< species number
      )
      {
        IntEnergy2deriv_[k] = material->ComputeInternalEnergy(concentration, refMu, refC, rt, 2);
        return;
      }

      //! Replace Internal energy 2nd derivative
      /*  This is needed because sometimes we obtain local negative values for the concentration due
         to the interpolation, which gives error computing the internal energy due to the log... In
         this case the opted solution is to interpolate not the concentration but the internal
         energy itself and we use this function to do that change)*/
      void ChangeInternalEnergy2deriv(
          const double NewIntEnergy,  //!< New Internal energy (1derv) values
          const int k                 //!< species number
      )
      {
        IntEnergy2deriv_[k] = NewIntEnergy;
        return;
      }

      //! Set Dissipation Potential
      void SetDissipationPot(const Teuchos::RCP<const MAT::ScatraMatVarChemDiffusion>
                                 material,             //!< pointer to current material
          double concentration,                        //!< Concentration at t_(n+1)
          const double refC,                           //!< Reference concentration
          const LINALG::Matrix<NSD, 1> chemicalfield,  //!< Chemical field at t_(n+1)
          const double mobility,                       //!< mobility
          const int k                                  //!< species number
      )
      {
        material->ComputeDissipationPot<NSD>(
            concentration, refC, chemicalfield, mobility, DissipationPot_[k]);
        return;
      }

      //! Set Dissipation Potential 1st derivative
      void SetDissipationPot1deriv(const Teuchos::RCP<const MAT::ScatraMatVarChemDiffusion>
                                       material,       //!< pointer to current material
          double concentration,                        //!< Concentration at t_(n+1)
          const double refC,                           //!< Reference concentration
          const LINALG::Matrix<NSD, 1> chemicalfield,  //!< Chemical field at t_(n+1)
          const double mobility,                       //!< mobility
          const int k                                  //!< species number
      )
      {
        material->ComputeDissipationPot_D1<NSD>(
            concentration, refC, chemicalfield, mobility, DissipationPot1deriv_[k]);
        return;
      }

      //! Set Dissipation Potential 2nd derivative
      void SetDissipationPot2deriv(const Teuchos::RCP<const MAT::ScatraMatVarChemDiffusion>
                                       material,       //!< pointer to current material
          double concentration,                        //!< Concentration at t_(n+1)
          const double refC,                           //!< Reference concentration
          const LINALG::Matrix<NSD, 1> chemicalfield,  //!< Chemical field at t_(n+1)
          const double mobility,                       //!< mobility
          const int k                                  //!< species number
      )
      {
        material->ComputeDissipationPot_D2<NSD>(
            concentration, refC, chemicalfield, mobility, DissipationPot2deriv_[k]);
        return;
      }


      //! Access routine for all (or of single ionic species k) mobility
      const std::vector<double> GetMobility() { return mobility_; };
      double GetMobility(const int k) { return mobility_[k]; };

      //! Access routine for all (or of single ionic species k) reference chemical potential
      const std::vector<double> GetRefMu() { return refMu_; };
      double GetRefMu(const int k) { return refMu_[k]; };

      //! Access routine for all (or of single ionic species k) reference concentration
      const std::vector<double> GetRefC() { return refC_; };
      double GetRefC(const int k) { return refC_[k]; };

      //! Access routine for RT factor
      double GetRT() { return rt_; };

      //! Access routine for Internal energy
      const std::vector<double> GetInternalEnergy() { return IntEnergy_; };
      double GetInternalEnergy(const int k) { return IntEnergy_[k]; };

      //! Access routine for Internal energy 1st derivatives
      const std::vector<LINALG::Matrix<NSD, 1>> GetInternalEnergy1deriv()
      {
        return IntEnergy1deriv_;
      };
      double GetInternalEnergy1deriv(const int k) { return IntEnergy1deriv_[k]; };

      //! Access routine for Internal energy 2nd derivatives
      const std::vector<double> GetInternalEnergy2deriv() { return IntEnergy2deriv_; };
      double GetInternalEnergy2deriv(const int k) { return IntEnergy2deriv_[k]; };

      //! Access routine for Dissipation Potential
      const std::vector<double> GetDissipationPot() { return DissipationPot_; };
      double GetDissipationPot(const int k) { return DissipationPot_[k]; };

      //! Access routine for Dissipation Potential 1st derivatives
      const std::vector<LINALG::Matrix<NSD, 1>> GetDissipationPot1deriv()
      {
        return DissipationPot1deriv_;
      };
      LINALG::Matrix<NSD, 1> GetDissipationPot1deriv(const int k)
      {
        return DissipationPot1deriv_[k];
      };

      //! Access routine for Dissipation Potential 2nd derivatives
      const std::vector<LINALG::Matrix<NSD, NSD>> GetDissipationPot2deriv()
      {
        return DissipationPot2deriv_;
      };
      LINALG::Matrix<NSD, NSD> GetDissipationPot2deriv(const int k)
      {
        return DissipationPot2deriv_[k];
      };

     protected:
      //! Reference mobility of the ionic species
      std::vector<double> mobility_;
      //! Reference chemical potential of the ionic species
      std::vector<double> refMu_;
      //! Reference concentration of the ionic species
      std::vector<double> refC_;
      //! Gas constant times temperature
      double rt_;
      //! Internal energy of the ionic species
      std::vector<double> IntEnergy_;
      //! Internal energy first derivative of the ionic species
      std::vector<double> IntEnergy1deriv_;
      //! Internal energy second derivative of the ionic species
      std::vector<double> IntEnergy2deriv_;
      //! Dissipation Potential of the ionic species
      std::vector<double> DissipationPot_;
      //! Dissipation Potential first derivative of the ionic species
      std::vector<LINALG::Matrix<NSD, 1>> DissipationPot1deriv_;
      //! Dissipation Potential second derivative of the  ionic species
      std::vector<LINALG::Matrix<NSD, NSD>> DissipationPot2deriv_;
    };


    /// ScaTraEleInternalVariableManagerVar implementation
    /*!
      This class manages all internal variables needed for the evaluation of an element.
      The internal variables stored in this class are used by the Nernst-Planck formulation
      as well as the diffusion-conduction formulation.
      All formulation-specific internal variables are stored and managed by a class derived from
      this class (class ScaTraEleInternalVariableManagerVarNP, class
      ScaTraEleInternalVariableManagerVarElectrode, and class
      ScaTraEleInternalVariableManagerVarDiffCond).
    */
    template <int NSD, int NEN>
    class ScaTraEleInternalVariableManagerVar : public ScaTraEleInternalVariableManager<NSD, NEN>
    {
      typedef ScaTraEleInternalVariableManager<NSD, NEN> my;
      typedef ScaTraEleInternalVariableManagerVar<NSD, NEN> myvar;

     public:
      ScaTraEleInternalVariableManagerVar(int numscal)
          : ScaTraEleInternalVariableManager<NSD, NEN>(numscal),
            conjugate_phinp_(numscal, 0.0),
            conjugate_phin_(numscal, 0.0),
            ConjugateField_(numscal),  // Chemical field
            phi0_(numscal, 0.0),
            conjugate_phi0_(numscal, 0.0)
      {
        return;
      }

      virtual ~ScaTraEleInternalVariableManagerVar(){};

      // compute and set internal variables used by chemical diffusion under a variational
      // fomrulation
      void SetInternalVariablesVariationalChemDiff(
          const LINALG::Matrix<NEN, 1>& funct,  //!< array for shape functions
          const LINALG::Matrix<NSD, NEN>&
              derxy,  //!< global derivatives of shape functions w.r.t x,y,z
          const std::vector<LINALG::Matrix<NEN, 1>>&
              ephinp,  //!< nodal state variables at t_(n+1) or t_(n+alpha_F)
          const std::vector<LINALG::Matrix<NEN, 1>>& ephin,  //!< nodal state variables at t_(n)
          const LINALG::Matrix<NSD, NEN>&
              econvelnp,  //!< nodal convective velocity values at t_(n+1) or t_(n+alpha_F)
          const std::vector<LINALG::Matrix<NEN, 1>>&
              ehist,  //!< history vector of transported scalars
          const std::vector<LINALG::Matrix<NEN, 1>>& ephi0  //!< initial value vector of d.o.f.
      )
      {
        // call base class (scatra)
        my::SetInternalVariables(funct, derxy, ephinp, ephin, econvelnp, ehist);

        const int ConjugateOffsetIndex = my::numscal_;

        for (int k = 0; k < ConjugateOffsetIndex; ++k)
        {
          // calculate initial values for scalar at t_(0)
          phi0_[k] = funct.Dot(ephi0[k]);
        }

        // This loop only concerns to the conjugate field
        for (int k = ConjugateOffsetIndex; k < 2 * my::numscal_; ++k)
        {
          // calculate conjugate scalar at t_(n+1) or t_(n+alpha_F)
          conjugate_phinp_[k - ConjugateOffsetIndex] = funct.Dot(ephinp[k]);
          // calculate conjugate scalar at t_(n)
          conjugate_phin_[k - ConjugateOffsetIndex] = funct.Dot(ephin[k]);
          // spatial negative gradient (field) of current conjugate scalar value
          ConjugateField_[k - ConjugateOffsetIndex].Multiply(-1.0, derxy, ephinp[k]);
          // calculate initial values for scalar at t_(0)
          conjugate_phi0_[k - ConjugateOffsetIndex] = funct.Dot(ephi0[k]);
        }

        return;
      };

      /*========================================================================*/
      //! @name return methods for internal variables
      /*========================================================================*/

      //! return conjugate scalar values at t_(n+1) or t_(n+alpha_F)
      const std::vector<double>& ConjugatePhinp() const { return conjugate_phinp_; };
      //! return conjugate scalar value at t_(n+1) or t_(n+alpha_F)
      const double& ConjugatePhinp(const int k) const { return conjugate_phinp_[k]; };
      //! return conjugate scalar values at t_(n)
      const std::vector<double>& ConjugatePhin() const { return conjugate_phin_; };
      //! return conjugate scalar value at t_(n)
      const double& ConjugatePhin(const int k) const { return conjugate_phin_[k]; };
      //! return spatial negative gradient of all conjugate scalar values
      const std::vector<LINALG::Matrix<NSD, 1>>& ConjugateField() const { return ConjugateField_; };
      //! return spatial negative gradient of current conjugate scalar value
      const LINALG::Matrix<NSD, 1>& ConjugateField(const int k) const
      {
        return ConjugateField_[k];
      };
      //! return initial scalar values at t_(0)
      const std::vector<double>& Phi0() const { return phi0_; };
      //! return initial scalar value at t_(0)
      const double& Phi0(const int k) const { return phi0_[k]; };
      //! return conjugate initial scalar values at t_(0)
      const std::vector<double>& ConjugatePhi0() const { return conjugate_phi0_; };
      //! return conjugate initial scalar value at t_(0)
      const double& ConjugatePhi0(const int k) const { return conjugate_phi0_[k]; };

     protected:
      /*========================================================================*/
      //! @name internal variables evaluated at element center or Gauss point
      /*========================================================================*/

      //! scalar conjugate at t_(n+1) or t_(n+alpha_F)
      std::vector<double> conjugate_phinp_;
      //! scalar conjugate at t_(n)
      std::vector<double> conjugate_phin_;
      //! spatial negative conjugate gradient of current scalar value (conjugate field)
      std::vector<LINALG::Matrix<NSD, 1>> ConjugateField_;
      //! scalar at t_(0) (initial values)
      std::vector<double> phi0_;
      //! scalar conjugate at t_(0) (initial values)
      std::vector<double> conjugate_phi0_;
    };
  }  // namespace ELEMENTS
}  // namespace DRT
#endif
