/*----------------------------------------------------------------------------*/
/*! \file
\brief 2D wall element for structure part of porous medium.

\level 2


*/
/*---------------------------------------------------------------------------*/

#ifndef FOUR_C_W1_PORO_HPP
#define FOUR_C_W1_PORO_HPP

#include "4C_config.hpp"

#include "4C_discretization_fem_general_utils_gausspoints.hpp"
#include "4C_inpar_structure.hpp"
#include "4C_w1.hpp"
#include "4C_w1_poro_eletypes.hpp"

FOUR_C_NAMESPACE_OPEN

namespace MAT
{
  class StructPoro;
  class FluidPoro;
  class FluidPoroMultiPhase;
}  // namespace MAT

namespace DRT
{
  class Discretization;

  namespace ELEMENTS
  {
    /*!
    \brief A C++ version of a 2 dimensional solid element with modifications for porous media

    */
    template <CORE::FE::CellType distype>
    class Wall1Poro : public Wall1
    {
     public:
      //!@}
      //! @name Constructors and destructors and related methods

      /*!
      \brief Standard Constructor

      \param id : A unique global id
      \param owner : elements owner
      */
      Wall1Poro(int id, int owner);

      /*!
      \brief Copy Constructor

      Makes a deep copy of a Element

      */
      Wall1Poro(const Wall1Poro& old);

      //!@}

      //! number of element nodes (
      static constexpr int numnod_ = CORE::FE::num_nodes<distype>;

      //! number of strains per node
      static constexpr int numstr_ = 3;

      //! number or second derivatives
      static constexpr int numderiv2_ = CORE::FE::DisTypeToNumDeriv2<distype>::numderiv2;

      //! number of degrees of freedom of element
      static constexpr int numdof_ = numnod_ * noddof_;

      //! total gauss points per element
      int numgpt_;

      //! @name Access methods

      /*!
      \brief Deep copy this instance of Solid3 and return pointer to the copy

      The Clone() method is used from the virtual base class Element in cases
      where the type of the derived class is unknown and a copy-ctor is needed

      */
      DRT::Element* Clone() const override;

      /*!
      \brief Return unique ParObject id

      every class implementing ParObject needs a unique id defined at the
      top of this file.
      */
      int UniqueParObjectId() const override
      {
        switch (distype)
        {
          case CORE::FE::CellType::tri3:
            return DRT::ELEMENTS::WallTri3PoroType::Instance().UniqueParObjectId();
          case CORE::FE::CellType::quad4:
            return DRT::ELEMENTS::WallQuad4PoroType::Instance().UniqueParObjectId();
          case CORE::FE::CellType::quad9:
            return DRT::ELEMENTS::WallQuad9PoroType::Instance().UniqueParObjectId();
          case CORE::FE::CellType::nurbs4:
            return DRT::ELEMENTS::WallNurbs4PoroType::Instance().UniqueParObjectId();
          case CORE::FE::CellType::nurbs9:
            return DRT::ELEMENTS::WallNurbs9PoroType::Instance().UniqueParObjectId();
          default:
            FOUR_C_THROW("unknown element type");
            break;
        }
        return -1;
      }

      /*!
      \brief Pack this class so it can be communicated

      \ref Pack and \ref Unpack are used to communicate this element

      */
      void Pack(CORE::COMM::PackBuffer& data) const override;

      /*!
      \brief Unpack data from a char vector into this class

      \ref Pack and \ref Unpack are used to communicate this element

      */
      void Unpack(const std::vector<char>& data) override;

      //! Get vector of Teuchos::RCPs to the lines of this element
      std::vector<Teuchos::RCP<DRT::Element>> Lines() override;

      //! Get vector of Teuchos::RCPs to the surfaces of this element
      std::vector<Teuchos::RCP<DRT::Element>> Surfaces() override;

      //! @name Access methods

      /*!
      \brief Print this element
      */
      void Print(std::ostream& os) const override;

      DRT::ElementType& ElementType() const override
      {
        switch (distype)
        {
          case CORE::FE::CellType::tri3:
            return DRT::ELEMENTS::WallTri3PoroType::Instance();
          case CORE::FE::CellType::quad4:
            return DRT::ELEMENTS::WallQuad4PoroType::Instance();
          case CORE::FE::CellType::quad9:
            return DRT::ELEMENTS::WallQuad9PoroType::Instance();
          case CORE::FE::CellType::nurbs4:
            return DRT::ELEMENTS::WallNurbs4PoroType::Instance();
          case CORE::FE::CellType::nurbs9:
            return DRT::ELEMENTS::WallNurbs9PoroType::Instance();
          default:
            FOUR_C_THROW("unknown element type");
            break;
        }
        return DRT::ELEMENTS::WallQuad4PoroType::Instance();
      }

      //!@}

      //! @name Evaluation

      void pre_evaluate(
          Teuchos::ParameterList&
              params,  //!< ParameterList for communication between control routine and elements
          DRT::Discretization& discretization,  //!< pointer to discretization for de-assembly
          DRT::Element::LocationArray& la       //!< location array for de-assembly
      );

      /*!
      \brief Evaluate an element

      Evaluate element stiffness, mass, internal forces, etc.

      If nullptr on input, the controlling method does not expect the element
      to fill these matrices or vectors.

      \return 0 if successful, negative otherwise
      */
      int Evaluate(
          Teuchos::ParameterList&
              params,  //!< ParameterList for communication between control routine and elements
          DRT::Discretization& discretization,  //!< pointer to discretization for de-assembly
          DRT::Element::LocationArray& la,      //!< location array for de-assembly
          CORE::LINALG::SerialDenseMatrix&
              elemat1,  //!< (stiffness-)matrix to be filled by element.
          CORE::LINALG::SerialDenseMatrix& elemat2,  //!< (mass-)matrix to be filled by element.
          CORE::LINALG::SerialDenseVector&
              elevec1,  //!< (internal force-)vector to be filled by element
          CORE::LINALG::SerialDenseVector& elevec2,  //!< vector to be filled by element
          CORE::LINALG::SerialDenseVector& elevec3   //!< vector to be filled by element
          ) override;


      //! initialize the inverse of the jacobian and its determinant in the material configuration
      virtual void InitElement();

      //!@}

      //! @name Input and Creation

      /*!
      \brief Read input for this element
      */
      bool ReadElement(const std::string& eletype, const std::string& eledistype,
          INPUT::LineDefinition* linedef) override;

      //!@}

      /*!
      \brief Read input for this element
      */
      /*!
      \brief Query names of element data to be visualized using BINIO

      The element fills the provided map with key names of
      visualization data the element wants to visualize AT THE CENTER
      of the element geometry. The values is supposed to be dimension of the
      data to be visualized. It can either be 1 (scalar), 3 (vector), 6 (sym. tensor)
      or 9 (nonsym. tensor)

      Example:
      \code
        // Name of data is 'Owner', dimension is 1 (scalar value)
        names.insert(std::pair<string,int>("Owner",1));
        // Name of data is 'StressesXYZ', dimension is 6 (sym. tensor value)
        names.insert(std::pair<string,int>("StressesXYZ",6));
      \endcode

      \param names (out): On return, the derived class has filled names with
                          key names of data it wants to visualize and with int dimensions
                          of that data.
      */
      void VisNames(std::map<std::string, int>& names) override;

      /*!
      \brief Query data to be visualized using BINIO of a given name

      The method is supposed to call this base method to visualize the owner of
      the element.
      If the derived method recognizes a supported data name, it shall fill it
      with corresponding data.
      If it does NOT recognizes the name, it shall do nothing.

      \warning The method must not change size of data

      \param name (in):   Name of data that is currently processed for visualization
      \param data (out):  data to be filled by element if element recognizes the name
      */
      bool VisData(const std::string& name, std::vector<double>& data) override;


      //! return anisotropic permeability directions (used for cloning)
      const std::vector<std::vector<double>>& get_anisotropic_permeability_directions() const
      {
        return anisotropic_permeability_directions_;
      }

      //! return scaling coefficients for anisotropic permeability (used for cloning)
      const std::vector<std::vector<double>>& get_anisotropic_permeability_nodal_coeffs() const
      {
        return anisotropic_permeability_nodal_coeffs_;
      }

      //! don't want = operator
      Wall1Poro& operator=(const Wall1Poro& old) = delete;

     protected:
      /*!
      \brief Evaluate an element

      Evaluate Wall1_poro element stiffness, mass, internal forces, etc.
      Templated evaluate routine of element matrixes

      If nullptr on input, the controlling method does not expect the element
      to fill these matrices or vectors.

      \return 0 if successful, negative otherwise
      */
      virtual int my_evaluate(
          Teuchos::ParameterList&
              params,  //!< ParameterList for communication between control routine and elements
          DRT::Discretization& discretization,  //!< pointer to discretization for de-assembly
          DRT::Element::LocationArray& la,      //!< location array for de-assembly
          CORE::LINALG::SerialDenseMatrix&
              elemat1,  //!< (stiffness-)matrix to be filled by element.
          CORE::LINALG::SerialDenseMatrix& elemat2,  //!< (mass-)matrix to be filled by element.
          CORE::LINALG::SerialDenseVector&
              elevec1,  //!< (internal force-)vector to be filled by element
          CORE::LINALG::SerialDenseVector& elevec2,  //!< vector to be filled by element
          CORE::LINALG::SerialDenseVector& elevec3   //!< vector to be filled by element
      );

      //! vector of inverses of the jacobian in material frame
      std::vector<CORE::LINALG::Matrix<numdim_, numdim_>> invJ_;
      //! determinant of Jacobian in material frame
      std::vector<double> detJ_;
      //! vector of coordinates of current integration point in reference coordinates
      std::vector<CORE::LINALG::Matrix<numdim_, 1>> xsi_;

      //! Calculate nonlinear stiffness and internal force for poroelasticity problems
      virtual void nonlinear_stiffness_poroelast(std::vector<int>& lm,  //!< location matrix
          CORE::LINALG::Matrix<numdim_, numnod_>& disp,                 //!< current displacements
          CORE::LINALG::Matrix<numdim_, numnod_>& vel,                  //!< current velocities
          CORE::LINALG::Matrix<numdim_, numnod_>& evelnp,       //!< fluid velocity of element
          CORE::LINALG::Matrix<numnod_, 1>& epreaf,             //!< fluid pressure of element
          CORE::LINALG::Matrix<numdof_, numdof_>* stiffmatrix,  //!< element stiffness matrix
          CORE::LINALG::Matrix<numdof_, numdof_>* reamatrix,    //!< element reactive matrix
          CORE::LINALG::Matrix<numdof_, 1>* force,              //!< element internal force vector
          Teuchos::ParameterList& params  //!< algorithmic parameters e.g. time
      );

      //! Calculate nonlinear stiffness and internal force for poroelasticity problems (pressure
      //! based formulation)
      virtual void nonlinear_stiffness_poroelast_pressure_based(
          std::vector<int>& lm,                          //!< location matrix
          CORE::LINALG::Matrix<numdim_, numnod_>& disp,  //!< current displacements
          const std::vector<double>& ephi,  //!< current primary variable for poro-multiphase flow
          CORE::LINALG::Matrix<numdof_, numdof_>* stiffmatrix,  //!< element stiffness matrix
          CORE::LINALG::Matrix<numdof_, 1>* force,              //!< element internal force vector
          Teuchos::ParameterList& params  //!< algorithmic parameters e.g. time
      );

      //! Calculate coupling terms in nonlinear stiffness and internal force for poroelasticity
      //! problems
      virtual void coupling_poroelast(std::vector<int>& lm,  //!< location matrix
          CORE::LINALG::Matrix<numdim_, numnod_>& disp,      //!< current displacements
          CORE::LINALG::Matrix<numdim_, numnod_>& vel,       //!< current velocities
          CORE::LINALG::Matrix<numdim_, numnod_>& evelnp,    //!< fluid velocity of element
          CORE::LINALG::Matrix<numnod_, 1>& epreaf,          //!< fluid pressure of element
          CORE::LINALG::Matrix<numdof_, (numdim_ + 1) * numnod_>*
              stiffmatrix,  //!< element stiffness matrix
          CORE::LINALG::Matrix<numdof_, (numdim_ + 1) * numnod_>*
              reamatrix,                            //!< element reactive matrix
          CORE::LINALG::Matrix<numdof_, 1>* force,  //!< element internal force vector
          Teuchos::ParameterList& params            //!< algorithmic parameters e.g. time
      );

      //! Calculate coupling terms in nonlinear stiffness and internal force for poroelasticity
      //! problems (pressure based formulation)
      virtual void coupling_poroelast_pressure_based(std::vector<int>& lm,  //!< location matrix
          CORE::LINALG::Matrix<numdim_, numnod_>& disp,  //!< current displacements
          const std::vector<double>& ephi,  //! current primary variable for poro-multiphase flow
          CORE::LINALG::SerialDenseMatrix& couplmat,  //!< element stiffness matrix
          Teuchos::ParameterList& params              //!< algorithmic parameters e.g. time
      );

      //! Calculate coupling stress for poroelasticity problems
      virtual void coupling_stress_poroelast(
          CORE::LINALG::Matrix<numdim_, numnod_>& disp,    //!< current displacements
          CORE::LINALG::Matrix<numdim_, numnod_>& evelnp,  //!< fluid velocity of element
          CORE::LINALG::Matrix<numnod_, 1>& epreaf,        //!< fluid pressure of element
          CORE::LINALG::SerialDenseMatrix* elestress,      //!< stresses at GP
          CORE::LINALG::SerialDenseMatrix* elestrain,      //!< strains at GP
          Teuchos::ParameterList& params,                  //!< algorithmic parameters e.g. time
          const INPAR::STR::StressType iostress            //!< stress output option
      );

      //! Gauss Point Loop evaluating stiffness and force
      void gauss_point_loop(Teuchos::ParameterList& params,
          const CORE::LINALG::Matrix<numdim_, numnod_>& xrefe,
          const CORE::LINALG::Matrix<numdim_, numnod_>& xcurr,
          const CORE::LINALG::Matrix<numdim_, numnod_>& nodaldisp,
          const CORE::LINALG::Matrix<numdim_, numnod_>& nodalvel,
          const CORE::LINALG::Matrix<numdim_, numnod_>& evelnp,
          const CORE::LINALG::Matrix<numnod_, 1>& epreaf,
          const CORE::LINALG::Matrix<numnod_, 1>* porosity_dof,
          CORE::LINALG::Matrix<numdof_, numdof_>& erea_v,
          CORE::LINALG::Matrix<numdof_, numdof_>* stiffmatrix,
          CORE::LINALG::Matrix<numdof_, numdof_>* reamatrix,
          CORE::LINALG::Matrix<numdof_, 1>* force);

      //! Gauss Point Loop evaluating stiffness (off diagonal)
      void gauss_point_loop_od(Teuchos::ParameterList& params,
          const CORE::LINALG::Matrix<numdim_, numnod_>& xrefe,
          const CORE::LINALG::Matrix<numdim_, numnod_>& xcurr,
          const CORE::LINALG::Matrix<numdim_, numnod_>& nodaldisp,
          const CORE::LINALG::Matrix<numdim_, numnod_>& nodalvel,
          const CORE::LINALG::Matrix<numdim_, numnod_>& evelnp,
          const CORE::LINALG::Matrix<numnod_, 1>& epreaf,
          const CORE::LINALG::Matrix<numnod_, 1>* porosity_dof,
          CORE::LINALG::Matrix<numdof_, (numdim_ + 1) * numnod_>& ecoupl);

      //! Gauss Point Loop evaluating stiffness and force
      void gauss_point_loop_pressure_based(Teuchos::ParameterList& params,
          const CORE::LINALG::Matrix<numdim_, numnod_>& xrefe,
          const CORE::LINALG::Matrix<numdim_, numnod_>& xcurr,
          const CORE::LINALG::Matrix<numdim_, numnod_>& nodaldisp, const std::vector<double>& ephi,
          CORE::LINALG::Matrix<numdof_, numdof_>* stiffmatrix,
          CORE::LINALG::Matrix<numdof_, 1>* force);

      //! Gauss Point Loop evaluating stiffness (off diagonal)
      void gauss_point_loop_od_pressure_based(Teuchos::ParameterList& params,
          const CORE::LINALG::Matrix<numdim_, numnod_>& xrefe,
          const CORE::LINALG::Matrix<numdim_, numnod_>& xcurr,
          const CORE::LINALG::Matrix<numdim_, numnod_>& nodaldisp, const std::vector<double>& ephi,
          CORE::LINALG::SerialDenseMatrix& couplmat);

      //! compute porosity at gausspoint and linearization of porosity w.r.t. structural
      //! displacements
      virtual void compute_porosity_and_linearization(Teuchos::ParameterList& params,
          const double& press, const double& J, const int& gp,
          const CORE::LINALG::Matrix<numnod_, 1>& shapfct,
          const CORE::LINALG::Matrix<numnod_, 1>* myporosity,
          const CORE::LINALG::Matrix<1, numdof_>& dJ_dus, double& porosity,
          CORE::LINALG::Matrix<1, numdof_>& dphi_dus);

      //! compute porosity at gausspoint and linearization of porosity w.r.t. fluid pressure
      virtual void compute_porosity_and_linearization_od(Teuchos::ParameterList& params,
          const double& press, const double& J, const int& gp,
          const CORE::LINALG::Matrix<numnod_, 1>& shapfct,
          const CORE::LINALG::Matrix<numnod_, 1>* myporosity, double& porosity, double& dphi_dp);

      //! Compute Jacobian Determinant
      void compute_jacobian_determinant_volume_change_and_linearizations(double& J,
          double& volchange, CORE::LINALG::Matrix<1, numdof_>& dJ_dus,
          CORE::LINALG::Matrix<1, numdof_>& dvolchange_dus,
          const CORE::LINALG::Matrix<numdim_, numdim_>& defgrd,
          const CORE::LINALG::Matrix<numdim_, numdim_>& defgrd_inv,
          const CORE::LINALG::Matrix<numdim_, numnod_>& N_XYZ,
          const CORE::LINALG::Matrix<numdim_, numnod_>& nodaldisp);

      //! Compute Jacobian Determinant
      void compute_jacobian_determinant_volume_change(double& J, double& volchange,
          const CORE::LINALG::Matrix<numdim_, numdim_>& defgrd,
          const CORE::LINALG::Matrix<numdim_, numnod_>& N_XYZ,
          const CORE::LINALG::Matrix<numdim_, numnod_>& nodaldisp);

      //! fill stiffness matrix and rhs vector for darcy flow
      void fill_matrix_and_vectors(const int& gp, const CORE::LINALG::Matrix<numnod_, 1>& shapefct,
          const CORE::LINALG::Matrix<numdim_, numnod_>& N_XYZ, const double& J, const double& press,
          const double& porosity, const CORE::LINALG::Matrix<numdim_, 1>& velint,
          const CORE::LINALG::Matrix<numdim_, 1>& fvelint,
          const CORE::LINALG::Matrix<numdim_, numdim_>& fvelder,
          const CORE::LINALG::Matrix<numdim_, numdim_>& defgrd_inv,
          const CORE::LINALG::Matrix<numstr_, numdof_>& bop,
          const CORE::LINALG::Matrix<numdim_, numdim_>& C_inv,
          const CORE::LINALG::Matrix<numdim_, 1>& Finvgradp,
          const CORE::LINALG::Matrix<1, numdof_>& dphi_dus,
          const CORE::LINALG::Matrix<1, numdof_>& dJ_dus,
          const CORE::LINALG::Matrix<numstr_, numdof_>& dCinv_dus,
          const CORE::LINALG::Matrix<numdim_, numdof_>& dFinvdus_gradp,
          const CORE::LINALG::Matrix<numdim_ * numdim_, numdof_>& dFinvTdus,
          CORE::LINALG::Matrix<numdof_, numdof_>& erea_v,
          CORE::LINALG::Matrix<numdof_, numdof_>* stiffmatrix,
          CORE::LINALG::Matrix<numdof_, 1>* force, CORE::LINALG::Matrix<numstr_, 1>& fstress);

      //! fill stiffness matrix and rhs vector for pressure-based formulation
      void fill_matrix_and_vectors_pressure_based(const int& gp,
          const CORE::LINALG::Matrix<numnod_, 1>& shapefct,
          const CORE::LINALG::Matrix<numdim_, numnod_>& N_XYZ, const double& J, const double& press,
          const CORE::LINALG::Matrix<numstr_, numdof_>& bop,
          const CORE::LINALG::Matrix<numdim_, numdim_>& C_inv,
          const CORE::LINALG::Matrix<1, numdof_>& dJ_dus,
          const CORE::LINALG::Matrix<numstr_, numdof_>& dCinv_dus,
          const CORE::LINALG::Matrix<1, numdof_>& dps_dus,
          CORE::LINALG::Matrix<numdof_, numdof_>* stiffmatrix,
          CORE::LINALG::Matrix<numdof_, 1>* force);

      //! fill stiffness matrix and rhs vector for brinkman flow
      void fill_matrix_and_vectors_brinkman(const int& gp, const double& J, const double& porosity,
          const CORE::LINALG::Matrix<numdim_, numdim_>& fvelder,
          const CORE::LINALG::Matrix<numdim_, numdim_>& defgrd_inv,
          const CORE::LINALG::Matrix<numstr_, numdof_>& bop,
          const CORE::LINALG::Matrix<numdim_, numdim_>& C_inv,
          const CORE::LINALG::Matrix<1, numdof_>& dphi_dus,
          const CORE::LINALG::Matrix<1, numdof_>& dJ_dus,
          const CORE::LINALG::Matrix<numstr_, numdof_>& dCinv_dus,
          const CORE::LINALG::Matrix<numdim_ * numdim_, numdof_>& dFinvTdus,
          CORE::LINALG::Matrix<numdof_, numdof_>* stiffmatrix,
          CORE::LINALG::Matrix<numdof_, 1>* force, CORE::LINALG::Matrix<numstr_, 1>& fstress);

      //! fill stiffness matrix and rhs vector for darcy flow (off diagonal terms)
      void fill_matrix_and_vectors_od(const int& gp,
          const CORE::LINALG::Matrix<numnod_, 1>& shapefct,
          const CORE::LINALG::Matrix<numdim_, numnod_>& N_XYZ, const double& J,
          const double& porosity, const double& dphi_dp,
          const CORE::LINALG::Matrix<numdim_, 1>& velint,
          const CORE::LINALG::Matrix<numdim_, 1>& fvelint,
          const CORE::LINALG::Matrix<numdim_, numdim_>& defgrd_inv,
          const CORE::LINALG::Matrix<numdim_, 1>& Gradp,
          const CORE::LINALG::Matrix<numstr_, numdof_>& bop,
          const CORE::LINALG::Matrix<numdim_, numdim_>& C_inv,
          CORE::LINALG::Matrix<numdof_, (numdim_ + 1) * numnod_>& ecoupl);

      //! fill stiffness matrix (off diagonal terms) -- pressure-based formulation
      void fill_matrix_and_vectors_od_pressure_based(const int& gp,
          const CORE::LINALG::Matrix<numnod_, 1>& shapefct,
          const CORE::LINALG::Matrix<numdim_, numnod_>& N_XYZ, const double& J,
          const CORE::LINALG::Matrix<numstr_, numdof_>& bop,
          const CORE::LINALG::Matrix<numdim_, numdim_>& C_inv,
          const std::vector<double>& solpressderiv, CORE::LINALG::SerialDenseMatrix& couplmat);

      //! fill stiffness matrix and rhs vector for darcy brinkman flow (off diagonal terms)
      void fill_matrix_and_vectors_brinkman_od(const int& gp,
          const CORE::LINALG::Matrix<numnod_, 1>& shapefct,
          const CORE::LINALG::Matrix<numdim_, numnod_>& N_XYZ, const double& J,
          const double& porosity, const double& dphi_dp,
          const CORE::LINALG::Matrix<numdim_, numdim_>& fvelder,
          const CORE::LINALG::Matrix<numdim_, numdim_>& defgrd_inv,
          const CORE::LINALG::Matrix<numstr_, numdof_>& bop,
          const CORE::LINALG::Matrix<numdim_, numdim_>& C_inv,
          CORE::LINALG::Matrix<numdof_, (numdim_ + 1) * numnod_>& ecoupl);

      //! Compute  nonlinear b-operator
      void compute_b_operator(CORE::LINALG::Matrix<numstr_, numdof_>& bop,
          const CORE::LINALG::Matrix<numdim_, numdim_>& defgrd,
          const CORE::LINALG::Matrix<numdim_, numnod_>& N_XYZ);

      //! evaluate shape functions and their derivatives at gauss point
      void compute_shape_functions_and_derivatives(const int& gp,
          CORE::LINALG::Matrix<numnod_, 1>& shapefct, CORE::LINALG::Matrix<numdim_, numnod_>& deriv,
          CORE::LINALG::Matrix<numdim_, numnod_>& N_XYZ);

      //! Compute Jacobian Determinant
      double compute_jacobian_determinant(const int& gp,
          const CORE::LINALG::Matrix<numdim_, numnod_>& xcurr,
          const CORE::LINALG::Matrix<numdim_, numnod_>& deriv);

      //! Compute Linearization Of Jacobian
      void compute_linearization_of_jacobian(CORE::LINALG::Matrix<1, numdof_>& dJ_dus,
          const double& J, const CORE::LINALG::Matrix<numdim_, numnod_>& N_XYZ,
          const CORE::LINALG::Matrix<numdim_, numdim_>& defgrd_inv);

      //! helper functions to get element vectors from global vector
      void compute_auxiliary_values(const CORE::LINALG::Matrix<numdim_, numnod_>& N_XYZ,
          const CORE::LINALG::Matrix<numdim_, numdim_>& defgrd_inv,
          const CORE::LINALG::Matrix<numdim_, numdim_>& C_inv,
          const CORE::LINALG::Matrix<numdim_, 1>& Gradp,
          CORE::LINALG::Matrix<numdim_ * numdim_, numdof_>& dFinvTdus,
          CORE::LINALG::Matrix<numdim_, 1>& Finvgradp,
          CORE::LINALG::Matrix<numdim_, numdof_>& dFinvdus_gradp,
          CORE::LINALG::Matrix<numstr_, numdof_>& dCinv_dus);

      //! push forward of material stresses to the current, spatial configuration (for output only)
      void p_k2to_cauchy(CORE::LINALG::Matrix<Wall1::numstr_, 1>& stress,
          CORE::LINALG::Matrix<numdim_, numdim_>& defgrd,
          CORE::LINALG::Matrix<numdim_, numdim_>& cauchystress);

      //! Compute deformation gradient
      void compute_def_gradient(CORE::LINALG::Matrix<numdim_, numdim_>&
                                    defgrd,  //!<<    (i) deformation gradient at gausspoint
          const CORE::LINALG::Matrix<numdim_, numnod_>&
              N_XYZ,  //!<<    (i) derivatives of shape functions w.r.t. reference coordinates
          const CORE::LINALG::Matrix<numdim_, numnod_>&
              xcurr  //!<<    (i) current position of gausspoint
      );

      //! helper functions to get element vectors from global vector
      void extract_values_from_global_vector(
          const DRT::Discretization& discretization,             //!< discretization
          const int& dofset,                                     //!< number of dofset
          const std::vector<int>& lm,                            //!< location vector
          CORE::LINALG::Matrix<numdim_, numnod_>* matrixtofill,  //!< vector field
          CORE::LINALG::Matrix<numnod_, 1>* vectortofill,        //!< scalar field
          const std::string state                                //!< state of the global vector
      );

      //! Compute Solid-pressure derivative w.r.t. primary variable at GP
      void compute_sol_pressure_deriv(const std::vector<double>& phiAtGP,  //!<< primary variable
          const int numfluidphases,             //!<< number of fluid phases
          std::vector<double>& solidpressderiv  //!<< solid pressure derivative at GP
      );

      //! Compute solid pressure at GP
      double compute_sol_pressure_at_gp(
          const int totalnumdofpernode,       //!<< total number of multiphase dofs
          const int numfluidphases,           //!<< number of fluid phases
          const std::vector<double>& phiAtGP  //!<< primary variable
      );

      //! recalculate solid pressure at GP in case of volfracs
      double recalculate_sol_pressure_at_gp(double press, const double porosity,
          const int totalnumdofpernode, const int numfluidphases, const int numvolfrac,
          const std::vector<double>& phiAtGP);

      //! recalculate solid pressure derivative in case of volfracs
      void recalculate_sol_pressure_deriv(const std::vector<double>& phiAtGP,
          const int totalnumdofpernode, const int numfluidphases, const int numvolfrac,
          const double press, const double porosity, std::vector<double>& solidpressderiv);

      //! Compute primary variable for multiphase flow at GP
      void compute_primary_variable_at_gp(
          const std::vector<double>& ephi,                   //!<< primary variable at node
          const int totalnumdofpernode,                      //!<< total number of multiphase dofs
          const CORE::LINALG::Matrix<numnod_, 1>& shapefct,  //!<< shapefct
          std::vector<double>& phiAtGP                       //!<< primary variable at GP
      );

      //! Compute linearizaton of solid press w.r.t. displacements
      //! only needed if additional volume fractions are present and porosity depends on
      //! Jacobian of deformation gradient
      void compute_linearization_of_sol_press_wrt_disp(const double fluidpress,
          const double porosity, const int totalnumdofpernode, const int numfluidphases,
          const int numvolfrac, const std::vector<double>& phiAtGP,
          const CORE::LINALG::Matrix<1, numdof_>& dphi_dus,
          CORE::LINALG::Matrix<1, numdof_>& dps_dus);

      //! get materials (solid and fluid)
      void get_materials();

      //! get materials (solid and fluidmulti)
      void get_materials_pressure_based();

      //! anisotropic permeability directions in the element definition
      void read_anisotropic_permeability_directions_from_element_line_definition(
          INPUT::LineDefinition* linedef);

      //! read nodal anisotropic permeability scaling coefficients in the element definition
      void read_anisotropic_permeability_nodal_coeffs_from_element_line_definition(
          INPUT::LineDefinition* linedef);

      //! interpolate the anisotropic permeability coefficients at GP from nodal values
      std::vector<double> compute_anisotropic_permeability_coeffs_at_gp(
          const CORE::LINALG::Matrix<numnod_, 1>& shapefct) const;

      //! Gauss integration rule
      CORE::FE::GaussIntegration intpoints_;

      //! flag indicating initialization of element
      bool init_;

      //! flag for scatra coupling
      bool scatra_coupling_;

      //! corresponding fluid material
      Teuchos::RCP<MAT::FluidPoro> fluid_mat_;

      //! corresponding multiphase fluid material
      Teuchos::RCP<MAT::FluidPoroMultiPhase> fluidmulti_mat_;

      //! own poro structure material
      Teuchos::RCP<MAT::StructPoro> struct_mat_;

      //! weights for nurbs elements
      CORE::LINALG::Matrix<numnod_, 1> weights_;
      //! knot vector for nurbs elements
      std::vector<CORE::LINALG::SerialDenseVector> myknots_;

      //! directions for anisotropic permeability
      std::vector<std::vector<double>> anisotropic_permeability_directions_;

      //! scaling coefficients for nodal anisotropic permeability
      std::vector<std::vector<double>> anisotropic_permeability_nodal_coeffs_;
    };

  }  // namespace ELEMENTS
}  // namespace DRT

FOUR_C_NAMESPACE_CLOSE

#endif
