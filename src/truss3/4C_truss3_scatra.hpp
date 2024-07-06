/*----------------------------------------------------------------------------*/
/*! \file
\brief three dimensional total Lagrange truss element used for scalar transport coupling

\level 3

*/
/*---------------------------------------------------------------------------*/
#ifndef FOUR_C_TRUSS3_SCATRA_HPP
#define FOUR_C_TRUSS3_SCATRA_HPP

#include "4C_config.hpp"

#include "4C_inpar_scatra.hpp"
#include "4C_truss3.hpp"

FOUR_C_NAMESPACE_OPEN

namespace Discret
{
  namespace ELEMENTS
  {
    class Truss3ScatraType : public Truss3Type
    {
     public:
      Core::Communication::ParObject* create(const std::vector<char>& data) override;

      Teuchos::RCP<Core::Elements::Element> create(const std::string eletype,
          const std::string eledistype, const int id, const int owner) override;

      Teuchos::RCP<Core::Elements::Element> create(const int id, const int owner) override;

      static Truss3ScatraType& instance();

      std::string name() const override { return "Truss3ScatraType"; }

      void setup_element_definition(
          std::map<std::string, std::map<std::string, Input::LineDefinition>>& definitions)
          override;

     private:
      static Truss3ScatraType instance_;
    };

    /*!
     \brief three dimensional total Lagrange truss element for scalar transport coupling

     */
    class Truss3Scatra : public Truss3
    {
     public:
      /*!
       \brief Standard Constructor

       \param id    (in): A globally unique element id
       \param owner (in): owner processor of the element
       */
      Truss3Scatra(int id, int owner);

      Truss3Scatra(const Truss3Scatra& old);

      Core::Elements::Element* clone() const override;

      Core::Elements::ElementType& element_type() const override
      {
        return Truss3ScatraType::instance();
      }

      /// return ScaTra::ImplType
      const Inpar::ScaTra::ImplType& impl_type() const { return impltype_; };

      bool read_element(const std::string& eletype, const std::string& distype,
          Input::LineDefinition* linedef) override;

      int unique_par_object_id() const override
      {
        return Truss3ScatraType::instance().unique_par_object_id();
      }

      void calc_internal_force_stiff_tot_lag(
          const std::map<std::string, std::vector<double>>& ele_state,
          Core::LinAlg::SerialDenseVector& forcevec,
          Core::LinAlg::SerialDenseMatrix& stiffmat) override;

      void calc_gp_stresses(Teuchos::ParameterList& params,
          const std::map<std::string, std::vector<double>>& ele_state) override;

      void pack(Core::Communication::PackBuffer& data) const override;
      void unpack(const std::vector<char>& data) override;

     protected:
      void extract_elemental_variables(LocationArray& la,
          const Core::FE::Discretization& discretization, const Teuchos::ParameterList& params,
          std::map<std::string, std::vector<double>>& ele_state) override;

      void energy(const std::map<std::string, std::vector<double>>& ele_state,
          Teuchos::ParameterList& params, Core::LinAlg::SerialDenseVector& intenergy) override;

     private:
      //! scalar transport implementation type (physics)
      Inpar::ScaTra::ImplType impltype_;

      //! evaluate elemental specific values
      //!
      //! @param[in] ele_state              elemental states (depending on the instantiated element)
      //! @param[out] curr_nodal_coords     nodal displacement
      //! @param[out] dcurr_nodal_coords_du derivative of truss displacement w.r.t. global
      //! displacement
      //! @param[out] dN_dx               derivative of shape functions
      //! @param[out] nodal_concentration nodal concentrations
      void prep_calc_internal_force_stiff_tot_lag_sca_tra(
          const std::map<std::string, std::vector<double>>& ele_state,
          Core::LinAlg::Matrix<6, 1>& curr_nodal_coords,
          Core::LinAlg::Matrix<6, 6>& dcurr_nodal_coords_du, Core::LinAlg::Matrix<6, 1>& dN_dx,
          Core::LinAlg::Matrix<2, 1>& nodal_concentration);

      //! calculation of concentration at Gauss Points, given concentration at nodes
      double project_scalar_to_gauss_point(double xi, const Core::LinAlg::Matrix<2, 1>& c) const;

      // don't want = operator
      Truss3Scatra& operator=(const Truss3Scatra& old);
    };
  }  // namespace ELEMENTS
}  // namespace Discret

FOUR_C_NAMESPACE_CLOSE

#endif
