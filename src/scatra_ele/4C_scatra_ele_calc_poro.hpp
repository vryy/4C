/*----------------------------------------------------------------------*/
/*! \file

 \brief evaluation class containing routines for calculation of scalar transport
        within porous medium

\level 2

 *----------------------------------------------------------------------*/


#ifndef FOUR_C_SCATRA_ELE_CALC_PORO_HPP
#define FOUR_C_SCATRA_ELE_CALC_PORO_HPP

#include "4C_config.hpp"

#include "4C_scatra_ele_calc.hpp"

FOUR_C_NAMESPACE_OPEN

// forward declarations
namespace Mat
{
  class ScatraMat;
}

namespace Discret
{
  namespace ELEMENTS
  {
    class ScaTraEleDiffManagerPoro;


    template <Core::FE::CellType distype>
    class ScaTraEleCalcPoro : public virtual ScaTraEleCalc<distype>
    {
     protected:
      /// (private) protected constructor, since we are a Singleton.
      ScaTraEleCalcPoro(const int numdofpernode, const int numscal, const std::string& disname);

     private:
      typedef ScaTraEleCalc<distype> my;
      using my::nen_;
      using my::nsd_;
      using my::nsd_ele_;

     public:
      /// Singleton access method
      static ScaTraEleCalcPoro<distype>* Instance(
          const int numdofpernode, const int numscal, const std::string& disname);

      /// Evaluate the element
      /*!
        Generic virtual interface function. Called via base pointer.
       */
      //   virtual int Evaluate(Core::Elements::Element*                 ele,
      //                        Teuchos::ParameterList&       params,
      //                        Discret::Discretization &         discretization,
      //                        const std::vector<int> &      lm,
      //                        Core::LinAlg::SerialDenseMatrix&     elemat1_epetra,
      //                        Core::LinAlg::SerialDenseMatrix&     elemat2_epetra,
      //                        Core::LinAlg::SerialDenseVector&     elevec1_epetra,
      //                        Core::LinAlg::SerialDenseVector&     elevec2_epetra,
      //                        Core::LinAlg::SerialDenseVector&     elevec3_epetra);

     protected:
      //! evaluate action
      int evaluate_action(Core::Elements::Element* ele, Teuchos::ParameterList& params,
          Discret::Discretization& discretization, const ScaTra::Action& action,
          Core::Elements::Element::LocationArray& la,
          Core::LinAlg::SerialDenseMatrix& elemat1_epetra,
          Core::LinAlg::SerialDenseMatrix& elemat2_epetra,
          Core::LinAlg::SerialDenseVector& elevec1_epetra,
          Core::LinAlg::SerialDenseVector& elevec2_epetra,
          Core::LinAlg::SerialDenseVector& elevec3_epetra) override;

      //   int EvaluateODMesh(  Core::Elements::Element*                 ele,
      //                        Teuchos::ParameterList&       params,
      //                        Discret::Discretization &         discretization,
      //                        const std::vector<int> &      lm,
      //                        Core::LinAlg::SerialDenseMatrix&     elemat1_epetra,
      //                        Core::LinAlg::SerialDenseMatrix&     elemat2_epetra,
      //                        Core::LinAlg::SerialDenseVector&     elevec1_epetra,
      //                        Core::LinAlg::SerialDenseVector&     elevec2_epetra,
      //                        Core::LinAlg::SerialDenseVector&     elevec3_epetra);
      //
      //   int EvaluateODFluid(  Core::Elements::Element*                 ele,
      //                        Teuchos::ParameterList&       params,
      //                        Discret::Discretization &         discretization,
      //                        const std::vector<int> &      lm,
      //                        Core::LinAlg::SerialDenseMatrix&     elemat1_epetra,
      //                        Core::LinAlg::SerialDenseMatrix&     elemat2_epetra,
      //                        Core::LinAlg::SerialDenseVector&     elevec1_epetra,
      //                        Core::LinAlg::SerialDenseVector&     elevec2_epetra,
      //                        Core::LinAlg::SerialDenseVector&     elevec3_epetra);
      //
      //   //! calculate matrix and rhs. Here the whole thing is hidden.
      //   virtual void sysmat_od_mesh(
      //     Core::Elements::Element*                         ele,       //!< the element we are
      //     dealing with Core::LinAlg::SerialDenseMatrix&             emat,      //!< element
      //     matrix to calculate const int                     numdofpernode
      //   );
      //
      //   //! calculate matrix and rhs. Here the whole thing is hidden.
      //   virtual void sysmat_od_fluid(
      //     Core::Elements::Element*                         ele,       //!< the element we are
      //     dealing with Core::LinAlg::SerialDenseMatrix&             emat,      //!< element
      //     matrix to calculate const int                     numdofpernode
      //   );

      //! read element coordinates
      void read_element_coordinates(const Core::Elements::Element* ele) override;

      //! extract element based or nodal values
      //  return extracted values of phinp
      void extract_element_and_node_values(Core::Elements::Element* ele,
          Teuchos::ParameterList& params, Discret::Discretization& discretization,
          Core::Elements::Element::LocationArray& la) override;

      //! extract element based or nodal values
      //  return extracted values of phinp
      virtual void extract_element_and_node_values_poro(Core::Elements::Element* ele,
          Teuchos::ParameterList& params, Discret::Discretization& discretization,
          Core::Elements::Element::LocationArray& la);

      //! get the material parameters
      void get_material_params(
          const Core::Elements::Element* ele,  //!< the element we are dealing with
          std::vector<double>& densn,          //!< density at t_(n)
          std::vector<double>& densnp,         //!< density at t_(n+1) or t_(n+alpha_F)
          std::vector<double>& densam,         //!< density at t_(n+alpha_M)
          double& visc,                        //!< fluid viscosity
          const int iquad = -1                 //!< id of current gauss point (default = -1)
          ) override;

      //! compute porosity based on solid, fluid and (potentially) scatra solution
      virtual void compute_porosity(
          const Core::Elements::Element* ele  //!< the element we are dealing with
      );

      //! compute pore pressure
      virtual double compute_pore_pressure();

      //! material ScaTra
      void mat_scatra(
          const Teuchos::RCP<const Core::Mat::Material> material,  //!< pointer to current material
          const int k,                                             //!< id of current scalar
          double& densn,                                           //!< density at t_(n)
          double& densnp,       //!< density at t_(n+1) or t_(n+alpha_F)
          double& densam,       //!< density at t_(n+alpha_M)
          double& visc,         //!< fluid viscosity
          const int iquad = -1  //!< id of current gauss point (default = -1)
          ) override;


      //! set diffusivity for poro scatra problem (i.e. scale by porosity)
      virtual void set_diffusivity(
          const Teuchos::RCP<const Mat::ScatraMat>& material, const int k, const double scale);

      //! set densisties for poro scatra problem (i.e. scale by porosity)
      virtual void set_densities(double porosity,
          double& densn,   //!< density at t_(n)
          double& densnp,  //!< density at t_(n+1) or t_(n+alpha_F)
          double& densam   //!< density at t_(n+alpha_M));
      );

      //! calculate scalar(s) and domain integral
      void calculate_scalars(const Core::Elements::Element* ele,
          Core::LinAlg::SerialDenseVector& scalars, bool inverting, bool calc_grad_phi) override;


      //! get poro diffusion manager
      Teuchos::RCP<ScaTraEleDiffManagerPoro> diff_manager()
      {
        return Teuchos::rcp_static_cast<ScaTraEleDiffManagerPoro>(my::diffmanager_);
      };

      /*========================================================================*/
      //! @name Galerkin approximation and related
      /*========================================================================*/

      //! initial node coordinates
      Core::LinAlg::Matrix<nsd_, nen_> xyze0_;

      //! nodal porosity values at t_(n+1)
      Core::LinAlg::Matrix<nen_, 1> eporosity_;

      //! flag indacting a node based porosity
      bool isnodalporosity_;
    };

    /// ScaTraEleDiffManagerPoro implementation
    /*!
      This class keeps all poro-specific transport parameter needed for the evaluation of an
      element. The ScaTraEleDiffManagerPoro is derived from the standard ScaTraEleDiffManager.
    */

    ////TODO: HACK
    // const int NO_CONVECTION_NR = 6;

    class ScaTraEleDiffManagerPoro : public ScaTraEleDiffManager
    {
     public:
      ScaTraEleDiffManagerPoro(int numscal) : ScaTraEleDiffManager(numscal), porosity_(0.0)
      {
        return;
      }

      void SetPorosity(double porosity) { porosity_ = porosity; }

      double GetPorosity(const int k) const
      {
        //      if(k<NO_CONVECTION_NR)
        return porosity_;
        //      else
        //        return 1.0;
      }

     protected:
      //! porosity at gauss point
      double porosity_;
    };

    template <int NSD, int NEN>
    class ScaTraEleInternalVariableManagerPoro : public ScaTraEleInternalVariableManager<NSD, NEN>
    {
      typedef ScaTraEleInternalVariableManager<NSD, NEN> my;

     public:
      ScaTraEleInternalVariableManagerPoro(int numscal)
          : ScaTraEleInternalVariableManager<NSD, NEN>(numscal),
            zeroconvel_(true),
            zeroconv_(true),
            zero_(0.0)
      {
        return;
      }

      virtual ~ScaTraEleInternalVariableManagerPoro() = default;

      /*========================================================================*/
      //! @name return methods for internal variables
      /*========================================================================*/

      //! return convective velocity
      virtual const Core::LinAlg::Matrix<NSD, 1>& ConVel(const int k) const {
          //    if(k<NO_CONVECTION_NR)
          {return my::convelint_;
    }
    //    else
    //      return zeroconvel_;
  };  // namespace ELEMENTS
  //! return convective part in convective form
  virtual const Core::LinAlg::Matrix<NEN, 1>& Conv(const int k) const
  {  //    if(k<NO_CONVECTION_NR)
    {
      return my::conv_;
    }  // namespace Discret
       //    else
       //      return zeroconv_;
  }

  //! return convective term of current scalar value
  virtual const double& ConvPhi(const int k) const
  {  //    if(k<NO_CONVECTION_NR)
    {
      return my::conv_phi_[k];
    }
    //    else
    //      return zero_;
  }

 private:
  Core::LinAlg::Matrix<NSD, 1> zeroconvel_;
  Core::LinAlg::Matrix<NEN, 1> zeroconv_;
  double zero_;
};  // namespace Discret
}
}


FOUR_C_NAMESPACE_CLOSE

#endif
