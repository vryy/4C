// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_FLUID_ELE_BOUNDARY_CALC_HPP
#define FOUR_C_FLUID_ELE_BOUNDARY_CALC_HPP

#include "4C_config.hpp"

#include "4C_fem_condition.hpp"
#include "4C_fem_general_element.hpp"
#include "4C_fem_general_utils_integration.hpp"
#include "4C_fem_general_utils_local_connectivity_matrices.hpp"
#include "4C_fluid_ele_boundary_interface.hpp"

FOUR_C_NAMESPACE_OPEN

namespace Core::FE
{
  class Discretization;
}  // namespace Core::FE

namespace Discret
{
  namespace Elements
  {
    class FluidBoundary;
    class FluidEleParameter;
    class FluidEleParameterTimInt;

    /// Internal FluidBoundary element implementation
    /*!
      This internal class keeps all the working arrays needed to
      calculate the FluidBoundary element. Additionally the method Sysmat()
      provides a clean and fast element implementation.

      <h3>Purpose</h3>

      The idea is to separate the element maintenance (class FluidBoundary)
      from the mathematical contents (this class). Of course there are
      different implementations of the FluidBoundary element, this is just one
      such implementation.

      The FluidBoundary element will allocate exactly one object of this class
      for all FluidBoundary elements with the same number of nodes in the mesh.
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
      \author ehrl
      \date 03/10
    */
    template <Core::FE::CellType distype>
    class FluidBoundaryImpl : public FluidBoundaryInterface
    {
      friend class FluidEleParameter;

     public:
      /// Constructor
      FluidBoundaryImpl();

      void evaluate_action(Discret::Elements::FluidBoundary* ele1, Teuchos::ParameterList& params,
          Core::FE::Discretization& discretization, std::vector<int>& lm,
          Core::LinAlg::SerialDenseMatrix& elemat1, Core::LinAlg::SerialDenseMatrix& elemat2,
          Core::LinAlg::SerialDenseVector& elevec1, Core::LinAlg::SerialDenseVector& elevec2,
          Core::LinAlg::SerialDenseVector& elevec3) override;

      //! number of element nodes
      static constexpr int bdrynen_ = Core::FE::num_nodes<distype>;

      //! number of spatial dimensions of boundary element
      static constexpr int bdrynsd_ = Core::FE::dim<distype>;

      //! number of spatial dimensions of parent element
      static constexpr int nsd_ = bdrynsd_ + 1;

      //! number of degrees of freedom per node
      static constexpr int numdofpernode_ = nsd_ + 1;

      //! Evaluate a Neumann boundary condition
      int evaluate_neumann(Discret::Elements::FluidBoundary* ele, Teuchos::ParameterList& params,
          Core::FE::Discretization& discretization, Core::Conditions::Condition& condition,
          std::vector<int>& lm, Core::LinAlg::SerialDenseVector& elevec1_epetra,
          Core::LinAlg::SerialDenseMatrix* elemat1) override;

     protected:
      /*!
      \brief integrate shapefunction over surface element

      \param params (in/out)    : ParameterList for communication between control routine
                                  and elements
      \param discretization (in): A reference to the underlying discretization
      \param lm (in)            : location vector of this element
      \param elevec1 (out)      : vector to be filled by element. If nullptr on input,
      */
      virtual void integrate_shape_function(Discret::Elements::FluidBoundary* ele,
          Teuchos::ParameterList& params, Core::FE::Discretization& discretization,
          std::vector<int>& lm, Core::LinAlg::SerialDenseVector& elevec1,
          const std::vector<double>& edispnp);

      /*!
        \brief Calculate mean curvature H. Interpolate the results to achieve better
        results in the surface tension algorithm (c0 field -> c1 field).

        \param elevec1  (out)     : Nodal values of mean curvature
        \param edispnp  (in)      : Displacement-vector
        \param enormals (in)      : Node normals
      */
      void element_mean_curvature(Discret::Elements::FluidBoundary* ele,
          Teuchos::ParameterList& params, Core::FE::Discretization& discretization,
          std::vector<int>& lm, Core::LinAlg::SerialDenseVector& elevec1,
          const std::vector<double>& edispnp, std::vector<double>& enormals);

      /*!
      brief integrate elemental areas over a surface

      \param params (in/out)    : ParameterList for communication between control routine
                                  and elements
      */
      void area_calculation(Discret::Elements::FluidBoundary* ele, Teuchos::ParameterList& params,
          Core::FE::Discretization& discretization, std::vector<int>& lm);

      /*!
       brief integrate elemental pressure over a surface

       \param params (in/out)    : ParameterList for communication between control routine
                                   and elements
       \param discretization (in): A reference to the underlying discretization
       \param lm (in)            : location vector of this element
       */
      void pressure_boundary_integral(Discret::Elements::FluidBoundary* ele,
          Teuchos::ParameterList& params, Core::FE::Discretization& discretization,
          std::vector<int>& lm);

      /*!
       brief center of mass of a surface

       \param params (in/out)    : ParameterList for communication between control routine
                                   and elements
       \param discretization (in): A reference to the underlying discretization
       \param lm (in)            : location vector of this element
       */
      void center_of_mass_calculation(Discret::Elements::FluidBoundary* ele,
          Teuchos::ParameterList& params, Core::FE::Discretization& discretization,
          std::vector<int>& lm);

      /*!
      \brief calculate the velocity component of the traction

      \param params (in/out)    : ParameterList for communication between control routine
                                  and elements
      \param discretization (in): A reference to the underlying discretization
      \param lm (in)            : location vector of this element
      \param elevec1 (out)      : vector to be filled by element. If nullptr on input,
      */
      void calc_traction_velocity_component(Discret::Elements::FluidBoundary* ele,
          Teuchos::ParameterList& params, Core::FE::Discretization& discretization,
          std::vector<int>& lm, Core::LinAlg::SerialDenseVector& elevec1);

      /*!
      \brief calculate the integral of the neumann inflow component

      \param params (in/out)    : ParameterList for communication between control routine
                                  and elements
      \param discretization (in): A reference to the underlying discretization
      \param lm (in)            : location vector of this element
      \param elevec1 (out)      : vector to be filled by element. If nullptr on input,
      */
      void compute_neumann_uv_integral(Discret::Elements::FluidBoundary* ele,
          Teuchos::ParameterList& params, Core::FE::Discretization& discretization,
          std::vector<int>& lm, Core::LinAlg::SerialDenseVector& elevec1);

      /*!
       brief integrate elemental flow rates over a surface

       \param params (in/out)    : ParameterList for communication between control routine
                                   and elements
       \param discretization (in): A reference to the underlying discretization
       \param lm (in)            : location vector of this element
       */
      void compute_flow_rate(Discret::Elements::FluidBoundary* ele, Teuchos::ParameterList& params,
          Core::FE::Discretization& discretization, std::vector<int>& lm,
          Core::LinAlg::SerialDenseVector& elevec1);

      /*!
      \brief determine elemental flow rate and its derivatives wrt dofs

      \param params (in/out)    : parameterList for communication between control routine
                                  and elements
      \param discretization (in): reference to the underlying discretization
      \param lm (in)            : location vector of this element
      \param elemat1 (out)      : 2nd derivative of flowrates wrt
                                  velocities and displacements
      \param elemat2 (out)      : 2nd derivative of flowrates wrt displacements
      \param elevec1 (out)      : derivative of flowrates wrt velocities
      \param elevec2 (out)      : derivative of flowrate wrt displacements
      \param elevec3 (out)      : flowrate
      */
      void flow_rate_deriv(Discret::Elements::FluidBoundary* ele, Teuchos::ParameterList& params,
          Core::FE::Discretization& discretization, std::vector<int>& lm,
          Core::LinAlg::SerialDenseMatrix& elemat1, Core::LinAlg::SerialDenseMatrix& elemat2,
          Core::LinAlg::SerialDenseVector& elevec1, Core::LinAlg::SerialDenseVector& elevec2,
          Core::LinAlg::SerialDenseVector& elevec3);

      /*!
      \brief apply impedance boundary condition (outlet pressures)

      \param params (in/out)    : ParameterList for communication between control routine
                                  and elements
      \param discretization (in): A reference to the underlying discretization
      \param lm (in)            : location vector of this element
      \param elevec1 (out)      : vector to be filled by element. If nullptr on input,
      */
      void impedance_integration(Discret::Elements::FluidBoundary* ele,
          Teuchos::ParameterList& params, Core::FE::Discretization& discretization,
          std::vector<int>& lm, Core::LinAlg::SerialDenseVector& elevec1);

      /*!
      \brief linearisation of flux w.r.t velocities on boundary elements

      \param params (in/out)    : ParameterList for communication between control routine
                                  and elements
      \param discretization (in): A reference to the underlying discretization
      \param lm (in)            : location vector of this element
      \param elevec1 (out)      : vector to be filled by element. If nullptr on input,
      */
      void d_qdu(Discret::Elements::FluidBoundary* ele, Teuchos::ParameterList& params,
          Core::FE::Discretization& discretization, std::vector<int>& lm,
          Core::LinAlg::SerialDenseVector& elevec1);

      /*!
      \brief compute potential Neumann inflow

      \param params (in)        : ParameterList for communication between control routine
      \param discretization (in): A reference to the underlying discretization
      \param lm (in)            : location vector of this element
      \param elemat1 (out)      : matrix to be filled by element. If nullptr on input,
                                  the controling method does not epxect the element to fill
                                  this matrix.                          and elements
      \param elevec1 (out)      : vector to be filled by element. If nullptr on input,

      */
      void neumann_inflow(Discret::Elements::FluidBoundary* ele, Teuchos::ParameterList& params,
          Core::FE::Discretization& discretization, std::vector<int>& lm,
          Core::LinAlg::SerialDenseMatrix& elemat1, Core::LinAlg::SerialDenseVector& elevec1);

     protected:
      void create_lines_tri(const int& nline, const int& nnode);

      void create_lines_quad(const int& nline, const int& nnode);

      //! get density
      void get_density(
          std::shared_ptr<const Core::Mat::Material> material,  ///< reference pointer to material
          const Core::LinAlg::Matrix<bdrynen_, 1>& escaaf,      ///< scalar at time n+alpha_f / n+1
          const double thermpressaf,  ///< thermodynamic pressure at time n+alpha_f / n+1
          const Core::LinAlg::Matrix<bdrynen_, 1>& epreaf);  ///< pressure at time n+alpha_f / n+1


      //! pointer to parameter list for time integration
      Discret::Elements::FluidEleParameterTimInt* fldparatimint_;
      //! pointer to parameter list
      Discret::Elements::FluidEleParameter* fldpara_;

      //! node coordinates for boundary element
      Core::LinAlg::Matrix<nsd_, bdrynen_> xyze_;
      //! node coordinates for boundary element at timestep n
      Core::LinAlg::Matrix<nsd_, bdrynen_> xyze_n_;
      //! coordinates of current integration point in reference coordinates
      Core::LinAlg::Matrix<bdrynsd_, 1> xsi_;
      //! array for shape functions for boundary element
      Core::LinAlg::Matrix<bdrynen_, 1> funct_;
      //! array for shape function derivatives for boundary element
      Core::LinAlg::Matrix<bdrynsd_, bdrynen_> deriv_;
      //! normal vector pointing out of the domain
      Core::LinAlg::Matrix<nsd_, 1> unitnormal_;
      //! normal vector pointing out of the domain at timestep n
      Core::LinAlg::Matrix<nsd_, 1> unitnormal_n_;
      //! velocity vector at integration point
      Core::LinAlg::Matrix<nsd_, 1> velint_;
      //! infinitesimal area element drs
      double drs_;
      //! integration factor
      double fac_;
      //! physical viscosity
      double visc_;
      //! density at t_(n+alpha_F) or t_(n+1)
      double densaf_;
      //! density factor for Neumann boundary conditions (set according to problem type)
      double densfac_;

    };  // class FluidBoundaryImpl

  }  // namespace Elements
}  // namespace Discret

FOUR_C_NAMESPACE_CLOSE

#endif
