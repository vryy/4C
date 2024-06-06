/*----------------------------------------------------------------------*/
/*! \file

\brief evaluation of scatra boundary terms at integration points

\level 1

 */
/*----------------------------------------------------------------------*/

#ifndef FOUR_C_SCATRA_ELE_BOUNDARY_CALC_HPP
#define FOUR_C_SCATRA_ELE_BOUNDARY_CALC_HPP

#include "4C_config.hpp"

#include "4C_discretization_fem_general_utils_local_connectivity_matrices.hpp"
#include "4C_fluid_ele.hpp"
#include "4C_inpar_elch.hpp"
#include "4C_scatra_ele_action.hpp"
#include "4C_scatra_ele_boundary_interface.hpp"
#include "4C_scatra_ele_calc_utils.hpp"

FOUR_C_NAMESPACE_OPEN

namespace FLD
{
  template <Core::FE::CellType distype, int numdofpernode,
      Discret::ELEMENTS::Fluid::EnrichmentType enrtype>
  class RotationallySymmetricPeriodicBC;
}

namespace Discret
{
  namespace ELEMENTS
  {
    class ScaTraEleParameterStd;
    class ScaTraEleParameterTimInt;
    class ScaTraEleParameterBoundary;

    /*! \brief Internal Scalar transport element implementation

      This internal class keeps all the working arrays needed to
      calculate the transport element. Additionally the method Sysmat()
      provides a clean and fast element implementation.

      <h3>Purpose</h3>

      The idea is to separate the element maintenance (class Transport)
      from the mathematical contents (this class). Of course there are
      different implementations of the Transport element, this is just one
      such implementation.

      The Transport element will allocate exactly one object of this class
      for all transport elements with the same number of nodes in the mesh.
      This allows us to use exactly matching working arrays (and keep them
      around.)

      The code is meant to be as clean as possible. This is the only way
      to keep it fast. The number of working arrays has to be reduced to
      a minimum so that the element fits into the cache. (There might be
      room for improvements.)

      <h3>History</h3>

      The implementation here is the standard convection-diffusion element
      capable of dealing with systems of transported scalars.

      Right now we do not read any stabilization parameters from the
      input file but have a fixed version.

      \author gjb
      \date 08/08
    */
    template <Core::FE::CellType distype, int probdim = Core::FE::dim<distype> + 1>
    class ScaTraEleBoundaryCalc : public ScaTraBoundaryInterface
    {
     public:
      //! number of element nodes (nomenclature: T. Hughes, The finite element method)
      static constexpr int nen_ = Core::FE::num_nodes<distype>;

      //! number of boundary(!) space dimensions
      static constexpr int nsd_ = probdim;
      static constexpr int nsd_ele_ = Core::FE::dim<distype>;

      // schott
      /*----------------------------------------------------------------------*/
      /*----------------------------------------------------------------------*/
      //! compute largest element diameter for reinitialization pseudo time step size
      template <Core::FE::CellType pdistype, class M1>
      double GetEleDiameter(const M1& xyze)
      {
        double elediam = 0.0;

        // number of nodes of this element
        const size_t numnode = Core::FE::num_nodes<pdistype>;

        // check all possible connections between nodes of an element
        // node 1 to 2
        //    :
        // node 1 to 8 = numnode
        // node 2 to 3
        //    :
        // node 2 to 8
        //    :
        //    :
        //    :
        // node 7 to 8
        for (size_t i_start = 0; i_start < numnode - 2; ++i_start)
        {
          for (size_t i_end = i_start + 1; i_end < numnode - 1; ++i_end)
          {
            Core::LinAlg::Matrix<3, 1> direction;
            direction.Clear();
            direction(0) = xyze(0, i_start) - xyze(0, i_end);
            direction(1) = xyze(1, i_start) - xyze(1, i_end);
            direction(2) = xyze(2, i_start) - xyze(2, i_end);

            // update elediam
            if (direction.Norm2() > elediam) elediam = direction.Norm2();
          }
        }

        return elediam;
      }

      //! setup element evaluation
      int SetupCalc(Core::Elements::FaceElement* ele, Teuchos::ParameterList& params,
          Discret::Discretization& discretization);

      /*!
       * @brief This method extracts the displacement values if ALE is activated
       *
       * @param[in] ele             currently evaluated face element
       * @param[in] discretization  discretization
       * @param[in] la              location array
       */
      void extract_displacement_values(Core::Elements::FaceElement* ele,
          const Discret::Discretization& discretization,
          Core::Elements::Element::LocationArray& la);

      /*!
       * @brief This method extracts the displacement values if ALE is activated
       *
       * @tparam parentdistype      discretization type of the parent element of the evaluated face
       *                            element
       * @param[in] ele             currently evaluated face element
       * @param[in] discretization  discretization
       * @param[in] la              location array
       */
      template <Core::FE::CellType parentdistype>
      void extract_displacement_values(Core::Elements::FaceElement* ele,
          const Discret::Discretization& discretization,
          Core::Elements::Element::LocationArray& la);

      //! Evaluate the element (using location array)
      int Evaluate(Core::Elements::FaceElement* ele, Teuchos::ParameterList& params,
          Discret::Discretization& discretization, Core::Elements::Element::LocationArray& la,
          Core::LinAlg::SerialDenseMatrix& elemat1_epetra,
          Core::LinAlg::SerialDenseMatrix& elemat2_epetra,
          Core::LinAlg::SerialDenseVector& elevec1_epetra,
          Core::LinAlg::SerialDenseVector& elevec2_epetra,
          Core::LinAlg::SerialDenseVector& elevec3_epetra) override;

      //! evaluate action
      virtual int evaluate_action(Core::Elements::FaceElement* ele, Teuchos::ParameterList& params,
          Discret::Discretization& discretization, ScaTra::BoundaryAction action,
          Core::Elements::Element::LocationArray& la,
          Core::LinAlg::SerialDenseMatrix& elemat1_epetra,
          Core::LinAlg::SerialDenseMatrix& elemat2_epetra,
          Core::LinAlg::SerialDenseVector& elevec1_epetra,
          Core::LinAlg::SerialDenseVector& elevec2_epetra,
          Core::LinAlg::SerialDenseVector& elevec3_epetra);

      //! evaluate Neumann boundary condition
      int evaluate_neumann(Core::Elements::FaceElement* ele, Teuchos::ParameterList& params,
          Discret::Discretization& discretization, Core::Conditions::Condition& condition,
          Core::Elements::Element::LocationArray& la, Core::LinAlg::SerialDenseVector& elevec1,
          const double scalar) override;

      /*!
       * \brief evaluate scatra-scatra interface coupling condition at integration point
       *
       * \remark This is a static method as it is also called from
       * `scatra_timint_meshtying_strategy_s2i.cpp` for the mortar implementation.
       *
       * @param[in] eslavephinp     state variables at slave-side nodes
       * @param[in] emasterphinp    state variables at master-side nodes
       * @param[in] pseudo_contact_fac  factor, modeling pseudo contact by considering the
       *                                mechanical stress state at the interface (1.0 if under
       *                                compressive stresses, 0.0 if under tensile stresses)
       * @param[in] funct_slave     slave-side shape function values
       * @param[in] funct_master    master-side shape function values
       * @param[in] test_slave      slave-side test function values
       * @param[in] test_master     master-side test function values
       * @param[in] numscal         number of transported scalars
       * @param[in] scatra_parameter_boundary  interface parameter class
       * @param[in] timefacfac      time-integration factor times domain-integration factor
       * @param[in] timefacrhsfac   time-integration factor for right-hand side times
       *                            domain-integration factor
       * @param[out] k_ss           linearizations of slave-side residuals w.r.t. slave-side dofs
       * @param[out] k_sm           linearizations of slave-side residuals w.r.t. master-side dofs
       * @param[out] k_ms           linearizations of master-side residuals w.r.t. slave-side dofs
       * @param[out] k_mm           linearizations of master-side residuals w.r.t. master-side dofs
       * @param[out] r_s            slave-side residual vector
       * @param[out] r_m            master-side residual vector
       *
       * \tparam distype_master This method is templated on the master-side discretization type.
       */
      template <Core::FE::CellType distype_master>
      static void evaluate_s2_i_coupling_at_integration_point(
          const std::vector<Core::LinAlg::Matrix<nen_, 1>>& eslavephinp,
          const std::vector<Core::LinAlg::Matrix<Core::FE::num_nodes<distype_master>, 1>>&
              emasterphinp,
          double pseudo_contact_fac, const Core::LinAlg::Matrix<nen_, 1>& funct_slave,
          const Core::LinAlg::Matrix<Core::FE::num_nodes<distype_master>, 1>& funct_master,
          const Core::LinAlg::Matrix<nen_, 1>& test_slave,
          const Core::LinAlg::Matrix<Core::FE::num_nodes<distype_master>, 1>& test_master,
          int numscal,
          const Discret::ELEMENTS::ScaTraEleParameterBoundary* const scatra_parameter_boundary,
          double timefacfac, double timefacrhsfac, Core::LinAlg::SerialDenseMatrix& k_ss,
          Core::LinAlg::SerialDenseMatrix& k_sm, Core::LinAlg::SerialDenseMatrix& k_ms,
          Core::LinAlg::SerialDenseMatrix& k_mm, Core::LinAlg::SerialDenseVector& r_s,
          Core::LinAlg::SerialDenseVector& r_m);

      /*!
       * @brief This method returns the determinant of the deformation gradient of the parent
       * element at the current gauss point coordinates
       *
       * @param[in] faceele     pointer to current face element
       * @param[in] faceele_xsi coordinates of current gauss point in parameter space of the face
       *                        element
       * @return determinant of the deformation gradient of the parent element
       */
      double calculate_det_f_of_parent_element(
          const Core::Elements::FaceElement* faceele, const double* faceele_xsi);

      /*!
       * @brief This method returns the determinant of the deformation gradient of the parent
       * element at the current gauss point coordinates
       *
       * @param[in] faceele     pointer to current face element
       * @param[in] faceele_xsi coordinates of current gauss point in parameter space of the face
       *                        element
       * @return determinant of the deformation gradient of the parent element
       */
      template <Core::FE::CellType parentdistype>
      double calculate_det_f_of_parent_element(
          const Core::Elements::FaceElement* faceele, const double* faceele_xi);

      /*!
       * @brief Method that calculates the prefactor indicating if interface at current gauss point
       *        is under tensile or compressive stresses. If it is under tensile stresses no flux is
       *        allowed and the prefactor is set to zero, otherwise it is set to one.
       *
       * @param[in] is_pseudo_contact           flag indicating if pseudo contact interface model
       *                                        should be evaluated
       * @param[in] eslavestress_vector         vector containing the stress contributions at the
       *                                        current gauss point
       * @param[in] gp_normal                   normal at the current gauss point
       * @param[in] funct_slave                 slave-side shape function values
       * @return prefactor indicating if interface at current gauss point is under tensile (0.0) or
       *         compressive stresses (1.0)
       */
      double calculate_pseudo_contact_factor(bool is_pseudo_contact,
          const std::vector<Core::LinAlg::Matrix<nen_, 1>>& eslavestress_vector,
          const Core::LinAlg::Matrix<nsd_, 1>& gp_normal,
          const Core::LinAlg::Matrix<nen_, 1>& funct_slave);

     private:
      // add nodal displacements to point coordinates
      void update_node_coordinates() { xyze_ += edispnp_; };

     protected:
      //! Constructor
      ScaTraEleBoundaryCalc(const int numdofpernode, const int numscal, const std::string& disname);

      //! calculate normal vectors
      void calc_normal_vectors(Teuchos::ParameterList& params, Core::Elements::FaceElement* ele);

      /*!
       * @brief evaluate shape functions, derivatives and domain integration factor at integration
       * point
       *
       * @param[in] intpoints   integration points
       * @param[in] iquad       id of current Gauss point
       * @param[out] normalvec  normal vector at Gauss point(optional)
       * @return domain integration factor
       */
      double eval_shape_func_and_int_fac(const Core::FE::IntPointsAndWeights<nsd_ele_>& intpoints,
          const int iquad, Core::LinAlg::Matrix<nsd_, 1>* normalvec = nullptr);

      /*!
       * @brief evaluate shape functions and derivatives at integration point
       *
       * @param[in] intpoints  integration points
       * @param[in] iquad      id of current Gauss point
       */
      void evaluate_shape_func_and_derivative_at_int_point(
          const Core::FE::IntPointsAndWeights<nsd_ele_>& intpoints, const int iquad);

      /*!
       * @brief Calculate the derivative of the square root of the determinant of the metric tensor
       * w.r.t. spatial coordinates
       *
       * @param[in] intpoints  integration points
       * @param[in] iquad      id of current Gauss point
       * @param[out] dsqrtdetg_dd  derivative of the square root of the determinant of the metric
       *                           tensor w.r.t. spatial coordinates
       */
      void evaluate_spatial_derivative_of_area_integration_factor(
          const Core::FE::IntPointsAndWeights<nsd_ele_>& intpoints, int iquad,
          Core::LinAlg::Matrix<nsd_, nen_>& dsqrtdetg_dd);

      //! evaluate scatra-scatra interface coupling condition
      virtual void evaluate_s2_i_coupling(
          const Core::Elements::FaceElement* ele,          ///< current boundary element
          Teuchos::ParameterList& params,                  ///< parameter list
          Discret::Discretization& discretization,         ///< discretization
          Core::Elements::Element::LocationArray& la,      ///< location array
          Core::LinAlg::SerialDenseMatrix& eslavematrix,   ///< element matrix for slave side
          Core::LinAlg::SerialDenseMatrix& emastermatrix,  ///< element matrix for master side
          Core::LinAlg::SerialDenseVector& eslaveresidual  ///< element residual for slave side
      );

      //! evaluate size of element projected to node, i.e. int N dv
      virtual void evaluate_nodal_size(
          const Core::Elements::FaceElement* ele,      ///< current boundary element
          Teuchos::ParameterList& params,              ///< parameter list
          Discret::Discretization& discretization,     ///< discretization
          Core::Elements::Element::LocationArray& la,  ///< location array
          Core::LinAlg::SerialDenseVector& nodalsize   ///< elemental size projected to nodes
      );

      /*!
       * @brief evaluate flux due to capacitance of scatra-scatra interface coupling condition
       *
       * @param[in] discretization   discretization
       * @param[in] la               location array
       * @param[out] eslavematrix    element matrix for slave side, i.e. slave-side residuals w.r.t.
       *                             slave-side dofs
       * @param[out] emastermatrix   element matrix for master side, i.e. master-side residuals
       *                             w.r.t. slave-side dofs
       * @param[out] eslaveresidual  element residual due to capacitive flux at slave side
       * @param[out] emasterresidual element residual due to capacitive flux at master side
       */
      virtual void evaluate_s2_i_coupling_capacitance(const Discret::Discretization& discretization,
          Core::Elements::Element::LocationArray& la, Core::LinAlg::SerialDenseMatrix& eslavematrix,
          Core::LinAlg::SerialDenseMatrix& emastermatrix,
          Core::LinAlg::SerialDenseVector& eslaveresidual,
          Core::LinAlg::SerialDenseVector& emasterresidual)
      {
        FOUR_C_THROW("not yet implemented!");
      };

      //! evaluate off-diagonal system matrix contributions associated with scatra-scatra interface
      //! coupling condition
      virtual void evaluate_s2_i_coupling_od(
          const Core::Elements::FaceElement* ele,        ///< current boundary element
          Teuchos::ParameterList& params,                ///< parameter list
          Discret::Discretization& discretization,       ///< discretization
          Core::Elements::Element::LocationArray& la,    ///< location array
          Core::LinAlg::SerialDenseMatrix& eslavematrix  ///< element matrix for slave side
      );

      virtual void evaluate_s2_i_coupling_capacitance_od(Teuchos::ParameterList& params,
          Discret::Discretization& discretization, Core::Elements::Element::LocationArray& la,
          Core::LinAlg::SerialDenseMatrix& eslavematrix,
          Core::LinAlg::SerialDenseMatrix& emastermatrix)
      {
        FOUR_C_THROW("not yet implemented");
      };

      //! extract nodal state variables associated with boundary element
      virtual void extract_node_values(
          const Discret::Discretization& discretization,  //!< discretization
          Core::Elements::Element::LocationArray& la      //!< location array
      );

      /*!
       * @brief extract nodal state variables associated with boundary element
       *
       * @param estate          nodal state variables
       * @param discretization  discretization
       * @param la              location array
       * @param statename       name of relevant state
       * @param nds             number of relevant dofset
       */
      void extract_node_values(Core::LinAlg::Matrix<nen_, 1>& estate,
          const Discret::Discretization& discretization, Core::Elements::Element::LocationArray& la,
          const std::string& statename = "phinp", const int& nds = 0) const;

      /*!
       * @brief extract nodal state variables associated with boundary element
       *
       * @param estate          nodal state variables
       * @param discretization  discretization
       * @param la              location array
       * @param statename       name of relevant state
       * @param nds             number of relevant dofset
       */
      void extract_node_values(std::vector<Core::LinAlg::Matrix<nen_, 1>>& estate,
          const Discret::Discretization& discretization, Core::Elements::Element::LocationArray& la,
          const std::string& statename = "phinp", const int& nds = 0) const;

      //! calculate boundary integral, i.e., surface area of boundary element
      void calc_boundary_integral(
          const Core::Elements::FaceElement* ele,  //!< the element we are dealing with
          Core::LinAlg::SerialDenseVector&
              scalar  //!< result vector for scalar integral to be computed
      );

      //! calculate boundary mass matrix
      void calc_mat_mass(const Core::Elements::FaceElement* const element,  //!< boundary element
          Core::LinAlg::SerialDenseMatrix& massmatrix                       //!< element mass matrix
      );

      //! evaluate Kedem-Katchalsky interface condition
      void evaluate_kedem_katchalsky(
          const Core::Elements::FaceElement* ele,      ///< current boundary element
          Teuchos::ParameterList& params,              ///< parameter list
          Discret::Discretization& discretization,     ///< discretization
          Core::Elements::Element::LocationArray& la,  ///< location array
          Core::LinAlg::SerialDenseMatrix& elemat1,    ///< element matrix for slave side
          Core::LinAlg::SerialDenseVector& elevec1     ///< element residual for slave side
      );

      //! integral of shape functions over boundary surface
      void integrate_shape_functions(
          const Core::Elements::FaceElement* ele,    ///< the actual boundary element
          Teuchos::ParameterList& params,            ///< the parameter list
          Core::LinAlg::SerialDenseVector& elevec1,  ///< result vector (to be assembled)
          const bool addarea                         ///< flag for area calculation
      );

      //! evaluate scatra-scatra interface coupling condition
      void evaluate_surface_permeability(
          const Core::Elements::FaceElement* ele,      ///< current boundary element
          Teuchos::ParameterList& params,              ///< parameter list
          Discret::Discretization& discretization,     ///< discretization
          Core::Elements::Element::LocationArray& la,  ///< location array
          Core::LinAlg::SerialDenseMatrix& elemat1,    ///< element matrix for slave side
          Core::LinAlg::SerialDenseVector& elevec1     ///< element residual for slave side
      );

      //! computes the factor for the wall shear stress dependent interface flux
      //! @note: Calculated as in Calvez, V., "Mathematical and numerical modeling of early
      //!        atherosclerotic lesions",ESAIM: Proceedings. Vol. 30. EDP Sciences, 2010
      double ws_sinfluence(const Core::LinAlg::Matrix<nsd_, nen_>&
                               ewss,  ///< vector of euclidean norm of shear stress
                                      ///< tensor of each node of the element
          const bool wss_onoff,       ///< flag if WSS has influence on interface permeability
          const std::vector<double>* coeffs);  ///< coefficients of the log law to determine the
                                               ///< influence of WSS on concentration flux

      //! Factor needed for the calculation of reference concentrations
      virtual double fac_for_ref_conc(const int iquad,  ///< current boundary integration point
          const Core::Elements::FaceElement* bele,      ///< current boundary element
          Teuchos::ParameterList& params,               ///< parameter list
          Discret::Discretization& discretization       ///< discretization
      )
      {
        return 1.0;
      };

      //! compute integral of convective mass/heat flux over boundary surface
      virtual std::vector<double> calc_convective_flux(const Core::Elements::FaceElement* ele,
          const std::vector<Core::LinAlg::Matrix<nen_, 1>>& ephinp,
          const Core::LinAlg::Matrix<nsd_, nen_>& evelnp, Core::LinAlg::SerialDenseVector& erhs);

      //! Compute the normal vector for a 2D boundary element in a 3D space
      //!
      //! \param xyze   position of boundary nodes
      //! \return       normal vector of boundary element
      Core::LinAlg::Matrix<3, 1> get_const_normal(const Core::LinAlg::Matrix<3, nen_>& xyze);

      //! Compute the normal vector for a 1D boundary element in a 2D space
      //!
      //! \param xyze   position of boundary nodes
      //! \return       normal vector of boundary element
      Core::LinAlg::Matrix<2, 1> get_const_normal(const Core::LinAlg::Matrix<2, nen_>& xyze);

      //! Compute the normal vector for a 1D boundary element in a 3D space
      //!
      //! \param xyze   position of boundary nodes
      //! \param nodes_parent_ele position of 3 nodes in parent element to compute normal vector of
      //! parent element
      //! \return       normal vector of boundary element
      Core::LinAlg::Matrix<3, 1> get_const_normal(const Core::LinAlg::Matrix<3, nen_>& xyze,
          const Core::LinAlg::Matrix<3, 3>& nodes_parent_ele);

      //! calculate potential Neumann inflow terms
      virtual void neumann_inflow(const Core::Elements::FaceElement* ele,
          Teuchos::ParameterList& params, Discret::Discretization& discretization,
          Core::Elements::Element::LocationArray& la, Core::LinAlg::SerialDenseMatrix& emat,
          Core::LinAlg::SerialDenseVector& erhs);

      //! get density at integration point
      virtual double get_density(Teuchos::RCP<const Core::Mat::Material> material,
          const std::vector<Core::LinAlg::Matrix<nen_, 1>>& ephinp, const int k);

      //! calculate boundary condition due to convective heat transfer
      void convective_heat_transfer(const Core::Elements::FaceElement* ele,
          Teuchos::RCP<const Core::Mat::Material> material,
          const std::vector<Core::LinAlg::Matrix<nen_, 1>>& ephinp,
          Core::LinAlg::SerialDenseMatrix& emat, Core::LinAlg::SerialDenseVector& erhs,
          const double heatranscoeff, const double surtemp);

      /*!
      \brief Evaluate weak Dirichlet boundary conditions

      \param params (in)        : ParameterList for communication between control routine
      \param discretization (in): A reference to the underlying discretization
      \param material (in)      : material of this element
      \param elemat1 (out)      : matrix to be filled by element.
      \param elevec1 (out)      : vector to be filled by element.

      */
      template <Core::FE::CellType bdistype, Core::FE::CellType pdistype>
      void weak_dirichlet(Core::Elements::FaceElement* ele, Teuchos::ParameterList& params,
          Discret::Discretization& discretization, Teuchos::RCP<const Core::Mat::Material> material,
          Core::LinAlg::SerialDenseMatrix& elemat_epetra,
          Core::LinAlg::SerialDenseVector& elevec_epetra);


      //! calculate boundary conditions for impl. Characteristic Galerkin time integration, just for
      //! the reinitialization equation
      template <Core::FE::CellType bdistype,
          Core::FE::CellType pdistype>
      void reinit_characteristic_galerkin_boundary(
          Core::Elements::FaceElement* ele,                  //!< transport element
          Teuchos::ParameterList& params,                    //!< parameter list
          Discret::Discretization& discretization,           //!< discretization
          Teuchos::RCP<const Core::Mat::Material> material,  //!< material
          Core::LinAlg::SerialDenseMatrix& elemat_epetra,    //!< ele sysmat
          Core::LinAlg::SerialDenseVector& elevec_epetra     //!< ele rhs
      );

      //! evaluate Robin boundary condition
      virtual void calc_robin_boundary(Core::Elements::FaceElement* ele,
          Teuchos::ParameterList& params, Discret::Discretization& discretization,
          Core::Elements::Element::LocationArray& la,  ///< location array
          Core::LinAlg::SerialDenseMatrix& elemat1_epetra,
          Core::LinAlg::SerialDenseVector& elevec1_epetra, const double scalar);

      //! evaluate integral of all positive fluxes on s2i condition
      virtual void calc_s2_i_coupling_flux(const Core::Elements::FaceElement* ele,
          const Teuchos::ParameterList& params, Discret::Discretization& discretization,
          Core::Elements::Element::LocationArray& la, Core::LinAlg::SerialDenseVector& scalars)
      {
        FOUR_C_THROW("Evaluation of coupling flux not implemented in base class.");
      }

      //! pointer to parameter list for time integration
      Discret::ELEMENTS::ScaTraEleParameterTimInt* scatraparamstimint_;
      //! pointer to parameter list
      Discret::ELEMENTS::ScaTraEleParameterStd* scatraparams_;
      //! pointer to scatra boundary parameter list
      Discret::ELEMENTS::ScaTraEleParameterBoundary* scatraparamsboundary_;

      //! number of dof per node
      const int numdofpernode_;
      //! number of transported scalars (numscal_ <= numdofpernode_)
      const int numscal_;

      //! node coordinates
      Core::LinAlg::Matrix<nsd_, nen_> xyze_;
      //! weights for nurbs elements
      Core::LinAlg::Matrix<nen_, 1> weights_;
      //! knot vector for nurbs elements
      std::vector<Core::LinAlg::SerialDenseVector> myknots_;
      //! knot vector of corresponding parent element
      std::vector<Core::LinAlg::SerialDenseVector> mypknots_;
      //! for nurbs elements the normal vector is multiplied with normalfac_!
      double normalfac_;
      //! nodal state variables associated with time t_{n+1} or t_{n+alpha_f}
      std::vector<Core::LinAlg::Matrix<nen_, 1>> ephinp_;
      //! nodal displacement values for ALE
      Core::LinAlg::Matrix<nsd_, nen_> edispnp_;
      //! nodal displacement values of parent element for ALE
      std::vector<double> eparentdispnp_;
      //! diffusivity / diffusivities (in case of systems) / thermal conductivity
      std::vector<double> diffus_;
      //! specific heat capacity at constant pressure (in case of temperature eq.)
      double shcacp_;
      //! coordinates of current integration point in reference coordinates
      Core::LinAlg::Matrix<nsd_ele_, 1> xsi_;
      //! array for shape functions
      Core::LinAlg::Matrix<nen_, 1> funct_;
      //! array for shape function derivatives w.r.t r,s,t
      Core::LinAlg::Matrix<nsd_ele_, nen_> deriv_;
      //! global derivatives of shape functions w.r.t x,y,z
      Core::LinAlg::Matrix<nsd_, nen_> derxy_;
      //! unit normal vector at integration point
      Core::LinAlg::Matrix<nsd_, 1> normal_;
      //! velocity vector in gausspoint
      Core::LinAlg::Matrix<nsd_, 1> velint_;
      //! metric tensor at integration point
      Core::LinAlg::Matrix<nsd_ele_, nsd_ele_> metrictensor_;
      //! for the handling of rotationally symmetric periodic boundary conditions
      Teuchos::RCP<
          FLD::RotationallySymmetricPeriodicBC<distype, nsd_ + 1, Discret::ELEMENTS::Fluid::none>>
          rotsymmpbc_;
    };
  }  // namespace ELEMENTS
}  // namespace Discret

FOUR_C_NAMESPACE_CLOSE

#endif
