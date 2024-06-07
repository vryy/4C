/*----------------------------------------------------------------------*/
/*! \file

\brief main file containing routines for calculation of HDG cardiac monodomain element

\level 3


*----------------------------------------------------------------------*/

#ifndef FOUR_C_SCATRA_ELE_CALC_CARDIAC_MONODOMAIN_HDG_HPP
#define FOUR_C_SCATRA_ELE_CALC_CARDIAC_MONODOMAIN_HDG_HPP


#include "4C_config.hpp"

#include "4C_fem_general_fiber_node_holder.hpp"
#include "4C_fem_general_utils_polynomial.hpp"
#include "4C_scatra_ele_calc_hdg.hpp"
#include "4C_scatra_ele_hdg.hpp"

FOUR_C_NAMESPACE_OPEN

namespace Discret
{
  namespace ELEMENTS
  {
    /// Scatra HDG element implementation
    template <Core::FE::CellType distype, int probdim = Core::FE::dim<distype>>
    class ScaTraEleCalcHDGCardiacMonodomain : public ScaTraEleCalcHDG<distype, probdim>
    {
     protected:
     private:
      /// (private) protected constructor, since we are a Singleton.
      /// this constructor is called from a derived class
      /// -> therefore, it has to be protected instead of private
      ScaTraEleCalcHDGCardiacMonodomain(
          const int numdofpernode, const int numscal, const std::string& disname);

      //    typedef ScaTraEleCalc<distype,probdim> my;
      //    typedef ScaTraEleCalcAniso<distype,probdim> aniso;
      //    typedef ScaTraEleCalcAdvReac<distype,probdim> advreac;

      std::vector<Core::LinAlg::SerialDenseVector> values_mat_gp_all_;
      std::vector<double> gp_mat_alpha_;

     public:
      /// Singleton access method
      static ScaTraEleCalcHDGCardiacMonodomain<distype, probdim>* Instance(const int numdofpernode,
          const int numscal, const std::string& disname, bool create = true);

      //! evaluate the element
      int evaluate_action(Core::Elements::Element* ele, Teuchos::ParameterList& params,
          Discret::Discretization& discretization, const ScaTra::Action& action,
          Core::Elements::Element::LocationArray& la,
          Core::LinAlg::SerialDenseMatrix& elemat1_epetra,
          Core::LinAlg::SerialDenseMatrix& elemat2_epetra,
          Core::LinAlg::SerialDenseVector& elevec1_epetra,
          Core::LinAlg::SerialDenseVector& elevec2_epetra,
          Core::LinAlg::SerialDenseVector& elevec3_epetra)
      {
        return 0;
      };

     protected:
      /*========================================================================*/
      //! @name material and related and related functions
      /*========================================================================*/

      //! evaluate material
      void prepare_materials(Core::Elements::Element* ele,  //!< the element we are dealing with
          const Teuchos::RCP<const Core::Mat::Material> material,  //!< pointer to current material
          const int k,                                             //!< id of current scalar
          Teuchos::RCP<std::vector<Core::LinAlg::SerialDenseMatrix>> difftensor) override;

      //! evaluate material
      virtual void prepare_materials_all(
          Core::Elements::Element* ele,  //!< the element we are dealing with
          const Teuchos::RCP<const Core::Mat::Material> material,  //!< pointer to current material
          const int k,                                             //!< id of current scalar
          Teuchos::RCP<std::vector<Core::LinAlg::SerialDenseMatrix>> difftensor);

      //! evaluate material
      virtual void prepare_materials_tet(
          Core::Elements::Element* ele,  //!< the element we are dealing with
          const Teuchos::RCP<const Core::Mat::Material> material,  //!< pointer to current material
          const int k,                                             //!< id of current scalar
          Teuchos::RCP<std::vector<Core::LinAlg::SerialDenseMatrix>> difftensor);


      //! evaluate material
      void materials(
          const Teuchos::RCP<const Core::Mat::Material> material,  //!< pointer to current material
          const int k,                                             //!< id of current scalar
          Core::LinAlg::SerialDenseMatrix& difftensor, Core::LinAlg::SerialDenseVector& ivecn,
          Core::LinAlg::SerialDenseVector& ivecnp,
          Core::LinAlg::SerialDenseMatrix& ivecnpderiv) override;

      //! material ScaTra
      void mat_myocard(
          const Teuchos::RCP<const Core::Mat::Material> material,  //!< pointer to current material
          const int k,                                             //!< id of current scalar
          Core::LinAlg::SerialDenseMatrix& difftensor, Core::LinAlg::SerialDenseVector& ivecn,
          Core::LinAlg::SerialDenseVector& ivecnp, Core::LinAlg::SerialDenseMatrix& ivecnpderiv);

      //! update time dependent material
      void time_update_material(
          const Core::Elements::Element* ele  //!< the element we are dealing with
          ) override;

      //! get material internal state for output
      void get_material_internal_state(const Core::Elements::Element* ele,
          Teuchos::ParameterList& params, Discret::Discretization& discretization) override;

      //! set material internal state after restart
      void set_material_internal_state(const Core::Elements::Element* ele,
          Teuchos::ParameterList& params, Discret::Discretization& discretization) override;

      //! project material field
      int project_material_field(const Core::Elements::Element* ele) override;

      //! project material field
      int project_material_field_all(const Core::Elements::Element* ele);

      //! project material field for Tet elements, because quadrature not working for higher order
      //! polynomials with Intrepid
      int project_material_field_tet(const Core::Elements::Element* ele);

      /*!
       * @brief Setup cardiac fibers. If the fiber direction fiber1 is directly given, it is
       * returned. Otherwise, the fiber will be calculated with the RAD-AXI-CIR coordinate system
       * and the helix and transversal angles
       *
       * @tparam dim space dimension of the problem
       * @param fibers Fiber data projected
       * @param f Cardiac fiber direction setup from the fiber or coordinate system data
       */
      template <std::size_t dim>
      void setup_cardiac_fibers(const Core::Nodes::NodalFiberHolder& fibers,
          std::vector<Core::LinAlg::Matrix<dim, 1>>& f);

      /// polynomial space for element interior for various Gauss Points for the evaluation of the
      /// material
      Teuchos::RCP<Core::FE::PolynomialSpace<probdim>> polySpace_;
    };

  }  // namespace ELEMENTS
}  // namespace Discret

FOUR_C_NAMESPACE_CLOSE

#endif
