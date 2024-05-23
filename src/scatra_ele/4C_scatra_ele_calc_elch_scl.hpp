/*--------------------------------------------------------------------------*/
/*! \file

\brief evaluation of ScaTra elements with scls

\level 2

*/

#ifndef FOUR_C_SCATRA_ELE_CALC_ELCH_SCL_HPP
#define FOUR_C_SCATRA_ELE_CALC_ELCH_SCL_HPP
#include "4C_config.hpp"

#include "4C_scatra_ele_calc_elch_diffcond.hpp"
#include "4C_scatra_ele_parameter_elch_diffcond.hpp"

FOUR_C_NAMESPACE_OPEN

namespace DRT
{
  namespace ELEMENTS
  {
    // forward declarations
    class ScaTraEleDiffManagerElchScl;
    template <int NSD, int NEN>
    class ScaTraEleInternalVariableManagerElchScl;
    template <CORE::FE::CellType distype>
    class ScaTraEleUtilsElchScl;

    // class implementation
    template <CORE::FE::CellType distype, int probdim = CORE::FE::dim<distype>>
    class ScaTraEleCalcElchScl : public ScaTraEleCalcElchDiffCond<distype, probdim>
    {
     protected:
      /// protected constructor, since we are a Singleton, but a derived class exists
      ScaTraEleCalcElchScl(const int numdofpernode, const int numscal, const std::string& disname);

      using my = ScaTraEleCalc<distype, probdim>;
      using myelch = ScaTraEleCalcElch<distype, probdim>;
      using mydiffcond = ScaTraEleCalcElchDiffCond<distype, probdim>;

     public:
      /// Singleton access method
      static ScaTraEleCalcElchScl<distype, probdim>* Instance(
          const int numdofpernode, const int numscal, const std::string& disname);

     protected:
      /*========================================================================*/
      //! @name general framework
      /*========================================================================*/

      void CalcMatAndRhs(CORE::LINALG::SerialDenseMatrix& emat,
          CORE::LINALG::SerialDenseVector& erhs, const int k, const double fac,
          const double timefacfac, const double rhsfac, const double taufac,
          const double timetaufac, const double rhstaufac,
          CORE::LINALG::Matrix<my::nen_, 1>& tauderpot, double& rhsint) override;

      //! CalcMat: Coulomb's equation (no source)
      //!
      //! \param emat  element matrix to be filled
      //! \param timefacfac  domain-integration factor times time integration factor
      //! \param invf  1/Faraday
      //! \param scalefac  scaling factor for consistent micro-macro-coupling
      //! \param gradpot  spatial gradient of electric potential
      //! \param epsilon  dielectric permittivity of bulk electrolyte

      void CalcMatPotCoulomb(CORE::LINALG::SerialDenseMatrix& emat, double timefacfac, double invf,
          double scalefac, const CORE::LINALG::Matrix<my::nsd_, 1>& gradpot, double epsilon);

      //! CalcRhs: Coulomb's equation (no source)
      //!
      //! \param erhs  element rhs to be filled
      //! \param fac  domain-integration factor
      //! \param invf  1/Faraday
      //! \param cond_invperm  conductivity / permitivity to scale pot. equations for consistent
      //! micro-macro-coupling
      //! \param gradpot  gradient of electric potential
      //! \param epsilon dielectric permittitivity of bulk electrolyte
      void CalcRhsPotCoulomb(CORE::LINALG::SerialDenseVector& erhs, double fac, double invf,
          double cond_invperm, const CORE::LINALG::Matrix<my::nsd_, 1>& gradpot, double epsilon);

      //! CalcRhs: Source contribution to electric potential (free charge)
      //!
      //! \param emat  element matrix to calculate
      //! \param k  index of current scalar
      //! \param fac  domain-integration factor
      //! \param invf  1/Faraday
      //! \param cond_invperm  conductivity / permitivity to scale pot. equations for consistent
      //! micro-macro-coupling
      //! \param z_k_F  transference number times Faraday constant
      void CalcMatPotSrc(CORE::LINALG::SerialDenseMatrix& emat, int k, double fac, double invf,
          double cond_invperm, double z_k_F);

      //! CalcRhs: Source contribution to electric potential (free charge)
      //!
      //! \param erhs  element rhs to calculate
      //! \param k  index of current scalar
      //! \param fac  domain-integration factor
      //! \param invf  1/Faraday
      //! \param cond_invperm  conductivity / permitivity to scale pot. equations for consistent
      //! micro-macro-coupling
      //! \param q_f  free charge density
      void CalcRhsPotSrc(CORE::LINALG::SerialDenseVector& erhs, int k, double fac, double invf,
          double cond_invperm, double q_f);

      //! CalcRhs: Calculate diffusion contribution to current
      //!
      //! \param erhs  element rhs to calculate
      //! \param rhsfac  time-integration factor for rhs times domain-integration factor
      //! \param invfval  transference number time Faraday const. inverted
      //! \param gradphi  ector of scalar gradients at t_(n+1)
      void CalcRhsDiffCur(CORE::LINALG::SerialDenseVector& erhs, double rhsfac,
          const std::vector<double>& invfval,
          const std::vector<CORE::LINALG::Matrix<my::nsd_, 1>>& gradphi);

      //! CalcMat: Calculate diffusion contribution to current
      //!
      //! \param emat   element matrix to calculate
      //! \param timefacfac  domain-integration factor times time-integration factor
      //! \param invfval  transference number time Faraday const. inverted
      //! \param gradphi  vector of scalar gradients at t_(n+1)
      void CalcMatDiffCur(CORE::LINALG::SerialDenseMatrix& emat, double timefacfac,
          const std::vector<double>& invfval,
          const std::vector<CORE::LINALG::Matrix<my::nsd_, 1>>& gradphi);

      void calc_mat_and_rhs_outside_scalar_loop(CORE::LINALG::SerialDenseMatrix& emat,
          CORE::LINALG::SerialDenseVector& erhs, const double fac, const double timefacfac,
          const double rhsfac) override;

      /*========================================================================*/
      //! @name material and related and related functions
      /*========================================================================*/

      void GetMaterialParams(const DRT::Element* ele, std::vector<double>& densn,
          std::vector<double>& densnp, std::vector<double>& densam, double& visc,
          const int iquad = -1) override;

      /*========================================================================*/
      //! @name methods for evaluation of individual terms
      /*========================================================================*/

      //! Calculate free charge based on concentration
      double CalcFreeCharge(const double concentration);

      //! Calculate derivative of free charge w.r.t. concentration
      double calc_free_charge_der_conc();

      /*========================================================================*/
      //! @name additional service routines
      /*========================================================================*/

      //! get diffusion manager for diffusion-conduction formulation
      Teuchos::RCP<ScaTraEleDiffManagerElchScl> DiffManager()
      {
        return Teuchos::rcp_static_cast<ScaTraEleDiffManagerElchScl>(my::diffmanager_);
      }

      //! get internal variable manager for diffusion-conduction formulation
      Teuchos::RCP<ScaTraEleInternalVariableManagerElchScl<my::nsd_, my::nen_>> VarManager()
      {
        return Teuchos::rcp_static_cast<
            ScaTraEleInternalVariableManagerElchScl<my::nsd_, my::nen_>>(my::scatravarmanager_);
      }

      //! get utility class supporting element evaluation for diffusion-conduction formulation
      ScaTraEleUtilsElchScl<distype>* Utils()
      {
        return static_cast<ScaTraEleUtilsElchScl<distype>*>(myelch::utils_);
      }

      //! flag for used element formulation (material based)
      INPAR::ELCH::DiffCondMat diffcondmat_;

      //! parameter class for diffusion-conduction formulation
      const ScaTraEleParameterElchDiffCond* diffcondparams_;
    };

    /// ScaTraEleDiffManagerElchSclimplementation
    /*!
      This class keeps all Diffusion-Conduction-specific transport parameter needed for the
      evaluation of an element.
    */
    class ScaTraEleDiffManagerElchScl : public ScaTraEleDiffManagerElchDiffCond
    {
     public:
      using dmdiffcond = ScaTraEleDiffManagerElchDiffCond;

      ScaTraEleDiffManagerElchScl(int numscal)
          : ScaTraEleDiffManagerElchDiffCond(numscal),
            c_max_(-1.0),
            c_bulk_(-1.0),
            c_lim_(-1.0),
            chi_(-1.0),
            epsilon_(-1.0)
      {
      }

      /*========================================================================*/
      //! @name access methods for transport parameter
      /*========================================================================*/

      //! Set dieelectric susceptibility of medium
      void SetSusceptibility(const double chi) { chi_ = chi; }

      //! Set dieelectric permittivity of the medium
      void SetPermittivity(const double epsilon) { epsilon_ = epsilon; }

      //! Get dieelectric permittivity of medium
      double GetPermittivity() { return epsilon_; }

      //! Set Bulk concentration of cations
      void SetBulkConc(const double cbulk) { c_bulk_ = cbulk; }

      //! Get Bulk concentration of cations
      double GetBulkConc() { return c_bulk_; }

      void SetTransNum(const double transnum, const int k) override
      {
        if (transnum != 1.0)
          FOUR_C_THROW("Only transference number of 1.0 is allowed in the SCL-model.");
        dmdiffcond::SetTransNum(transnum, k);
      }

      void SetDerivTransNum(const double derivtransnum, const int k, const int iscal) override
      {
        if (derivtransnum != 0.0)
        {
          FOUR_C_THROW(
              "Only constant transference number (1.0) is allowed. Derivative w.r.t. concentration "
              "has to be zero.");
        }
        dmdiffcond::SetDerivTransNum(derivtransnum, k, iscal);
      }

      void SetPhasePoro(const double eps, const int phase) override
      {
        // Only porosity of 1.0 can be assigned to the SCL-material
        if (eps != 1.0) FOUR_C_THROW("Only constant porosity (1.0) is allowed for SCL-material.");
        dmdiffcond::SetPhasePoro(eps, phase);
      }

      void SetPhaseTort(const double tort, const int phase) override
      {
        // only tortuosity of 1.0 can be assigned to the SCL-material
        if (tort != 1.0)
          FOUR_C_THROW("Only constant tortuosity (1.0) is allowed for SCL-material.");
        dmdiffcond::SetPhaseTort(tort, phase);
      }

      /*========================================================================*/
      //! @name output
      /*========================================================================*/

      void output_transport_params(const int numscal) override
      {
        // call base class routine
        dmdiffcond::output_transport_params(numscal);

        // additional outputs
        std::cout << "maximum concentration:    " << c_max_ << std::endl;

        for (int k = 0; k < numscal; ++k)
          std::cout << "limiting concentration:    " << c_lim_ << std::endl;

        std::cout << "dieelectric susceptibility   " << chi_ << std::endl;
        std::cout << "dieelectric permittivity   " << epsilon_ << std::endl;
      }

     private:
      /*========================================================================*/
      //! @name transport parameter
      /*========================================================================*/

      //! maximum concentration of species
      double c_max_;

      //! bulk concentration (= anion concentration)
      double c_bulk_;

      //! limit concentration for extrapolation strategy of diffusion coefficient
      double c_lim_;

      //! dieelectric susceptibility of electrolyte material
      double chi_;

      //! dieelectric permittivity of electrolyte material
      double epsilon_;

      /*========================================================================*/
      //! @name geometrical parameters of the porous medium
      /*========================================================================*/
    };

    /// ScaTraEleInternalVariableManagerElchDiffCond implementation
    /*!
      This class keeps all internal variables needed for the diffusion-conduction formulation.
    */
    template <int NSD, int NEN>
    class ScaTraEleInternalVariableManagerElchScl
        : public ScaTraEleInternalVariableManagerElchDiffCond<NSD, NEN>
    {
     public:
      using vm = ScaTraEleInternalVariableManager<NSD, NEN>;

      ScaTraEleInternalVariableManagerElchScl(int numscal,
          const DRT::ELEMENTS::ScaTraEleParameterElch* elchparams,
          const DRT::ELEMENTS::ScaTraEleParameterElchDiffCond* diffcondparams)
          : ScaTraEleInternalVariableManagerElchDiffCond<NSD, NEN>(
                numscal, elchparams, diffcondparams),
            curint_(true)
      {
      }

     private:
      /*========================================================================*/
      //! @name internal variables evaluated at the Gauss point
      /*========================================================================*/

      //! current density at Gauss point
      CORE::LINALG::Matrix<NSD, 1> curint_;
    };
  }  // namespace ELEMENTS
}  // namespace DRT
FOUR_C_NAMESPACE_CLOSE

#endif
