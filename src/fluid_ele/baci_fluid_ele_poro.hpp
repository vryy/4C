/*-----------------------------------------------------------*/
/*! \file

\brief A C++ wrapper for the poro fluid element

This file contains the element-specific service routines such as
Pack, Unpack, NumDofPerNode etc.


\level 2

*/
/*-----------------------------------------------------------*/


#ifndef FOUR_C_FLUID_ELE_PORO_HPP
#define FOUR_C_FLUID_ELE_PORO_HPP

#include "baci_config.hpp"

#include "baci_fluid_ele.hpp"
#include "baci_inpar_structure.hpp"

BACI_NAMESPACE_OPEN

namespace DRT
{
  class Discretization;

  namespace ELEMENTS
  {
    class FluidPoroEleType : public FluidType
    {
     public:
      std::string Name() const override { return "FluidPoroEleType"; }

      static FluidPoroEleType& Instance();

      CORE::COMM::ParObject* Create(const std::vector<char>& data) override;

      Teuchos::RCP<DRT::Element> Create(const std::string eletype, const std::string eledistype,
          const int id, const int owner) override;

      Teuchos::RCP<DRT::Element> Create(const int id, const int owner) override;

      void SetupElementDefinition(
          std::map<std::string, std::map<std::string, INPUT::LineDefinition>>& definitions)
          override;

      //! pre-evaluation
      void PreEvaluate(DRT::Discretization& dis, Teuchos::ParameterList& p,
          Teuchos::RCP<CORE::LINALG::SparseOperator> systemmatrix1,
          Teuchos::RCP<CORE::LINALG::SparseOperator> systemmatrix2,
          Teuchos::RCP<Epetra_Vector> systemvector1, Teuchos::RCP<Epetra_Vector> systemvector2,
          Teuchos::RCP<Epetra_Vector> systemvector3) override;

     private:
      static FluidPoroEleType instance_;
    };

    /*!
    \brief A C++ wrapper for the fluid element
    */
    class FluidPoro : public virtual Fluid
    {
     public:
      //!@}
      //! @name constructors and destructors and related methods

      /*!
      \brief standard constructor
      */
      FluidPoro(int id,  //!< A unique global id
          int owner      //!< number of owning processor
      );

      /*!
      \brief Copy Constructor

      Makes a deep copy of a Element

      */
      FluidPoro(const FluidPoro& old);

      /*!
      \brief Deep copy this instance of fluid and return pointer to the copy

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
        return FluidPoroEleType::Instance().UniqueParObjectId();
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

      /*!
      \brief Get vector of Teuchos::RCPs to the lines of this element
      */
      std::vector<Teuchos::RCP<DRT::Element>> Lines() override;

      /*!
      \brief Get vector of Teuchos::RCPs to the surfaces of this element
      */
      std::vector<Teuchos::RCP<DRT::Element>> Surfaces() override;

      //!@}

      //! @name Acess methods

      /*!
      \brief Print this element
      */
      void Print(std::ostream& os) const override;

      DRT::ElementType& ElementType() const override { return FluidPoroEleType::Instance(); }


      void SetKinematicType(INPAR::STR::KinemType kintype) { kintype_ = kintype; }

      INPAR::STR::KinemType KinematicType() const { return kintype_; }

      //! Set anisotropic permeability directions obtained from structure element during cloning
      void SetAnisotropicPermeabilityDirections(
          const std::vector<std::vector<double>>& ref_anisotropic_permeability_directions)
      {
        anisotropic_permeability_directions_ = ref_anisotropic_permeability_directions;
      }

      //! Set anisotropic permeability coefficients obtained from structure element during cloning
      void SetAnisotropicPermeabilityNodalCoeffs(
          const std::vector<std::vector<double>>& ref_anisotropic_permeability_nodal_coeffs)
      {
        anisotropic_permeability_nodal_coeffs_ = ref_anisotropic_permeability_nodal_coeffs;
      }

      //! Provide the anisotropic permeability directions of the element
      const std::vector<std::vector<double>>& GetAnisotropicPermeabilityDirections() const
      {
        return anisotropic_permeability_directions_;
      }

      //! Provide the nodal anisotropic permeability coefficients of the element
      const std::vector<std::vector<double>>& GetAnisotropicPermeabilityNodalCoeffs() const
      {
        return anisotropic_permeability_nodal_coeffs_;
      }

      //!@}

      //! @name Evaluation

      /*!
      \brief Evaluate an element, that is, call the element routines to evaluate fluid
      element matrices and vectors or evaluate errors, statistics or updates etc. directly.

      \param params (in/out): ParameterList for communication between control routine
                              and elements
      \param elemat1 (out)  : matrix to be filled by element. If nullptr on input,
                              the controling method does not epxect the element to fill
                              this matrix.
      \param elemat2 (out)  : matrix to be filled by element. If nullptr on input,
                              the controling method does not epxect the element to fill
                              this matrix.
      \param elevec1 (out)  : vector to be filled by element. If nullptr on input,
                              the controlling method does not epxect the element
                              to fill this vector
      \param elevec2 (out)  : vector to be filled by element. If nullptr on input,
                              the controlling method does not epxect the element
                              to fill this vector
      \param elevec3 (out)  : vector to be filled by element. If nullptr on input,
                              the controlling method does not epxect the element
                              to fill this vector
      \return 0 if successful, negative otherwise
      */
      int Evaluate(Teuchos::ParameterList& params, DRT::Discretization& discretization,
          std::vector<int>& lm, CORE::LINALG::SerialDenseMatrix& elemat1,
          CORE::LINALG::SerialDenseMatrix& elemat2, CORE::LINALG::SerialDenseVector& elevec1,
          CORE::LINALG::SerialDenseVector& elevec2,
          CORE::LINALG::SerialDenseVector& elevec3) override;

      //!@}

     private:
      //! don't want = operator
      FluidPoro& operator=(const FluidPoro& old);

      //! kinematic type
      INPAR::STR::KinemType kintype_;

      //! directions for anisotropic permeability
      std::vector<std::vector<double>> anisotropic_permeability_directions_;

      //! nodal coefficients for anisotropic permeability
      std::vector<std::vector<double>> anisotropic_permeability_nodal_coeffs_;
    };

    /*!
    \brief An element representing a boundary element of a fluid element

    */

    class FluidPoroBoundaryType : public FluidBoundaryType
    {
     public:
      std::string Name() const override { return "FluidPoroBoundaryType"; }

      Teuchos::RCP<DRT::Element> Create(const int id, const int owner) override;

      static FluidPoroBoundaryType& Instance();

     private:
      static FluidPoroBoundaryType instance_;
    };

    class FluidPoroBoundary : public FluidBoundary
    {
     public:
      //! @name Constructors and destructors and related methods

      //! number of space dimensions
      /*!
      \brief Standard Constructor

      \param id : A unique global id
      \param owner: Processor owning this surface
      \param nnode: Number of nodes attached to this element
      \param nodeids: global ids of nodes attached to this element
      \param nodes: the discretizations map of nodes to build ptrs to nodes from
      \param parent: The parent fluid element of this surface
      \param lsurface: the local surface number of this surface w.r.t. the parent element
      */
      FluidPoroBoundary(int id, int owner, int nnode, const int* nodeids, DRT::Node** nodes,
          DRT::ELEMENTS::Fluid* parent, const int lsurface);

      /*!
      \brief Copy Constructor

      Makes a deep copy of a Element

      */
      FluidPoroBoundary(const FluidPoroBoundary& old);

      /*!
      \brief Deep copy this instance of an element and return pointer to the copy

      The Clone() method is used from the virtual base class Element in cases
      where the type of the derived class is unknown and a copy-ctor is needed

      */
      DRT::Element* Clone() const override;

      /*!
      \brief Return unique ParObject id

      every class implementing ParObject needs a unique id defined at the
      top of the parobject.H file.
      */
      int UniqueParObjectId() const override
      {
        return FluidPoroBoundaryType::Instance().UniqueParObjectId();
      }

      //!@}

      //! @name Acess methods

      /*!
      \brief Print this element
      */
      void Print(std::ostream& os) const override;

      DRT::ElementType& ElementType() const override { return FluidPoroBoundaryType::Instance(); }

      //!@}

      //! @name Evaluation

      /*!
      \brief Evaluate element

      \param params (in/out): ParameterList for communication between control routine
                              and elements
      \param elemat1 (out)  : matrix to be filled by element. If nullptr on input,
                              the controling method does not epxect the element to fill
                              this matrix.
      \param elemat2 (out)  : matrix to be filled by element. If nullptr on input,
                              the controling method does not epxect the element to fill
                              this matrix.
      \param elevec1 (out)  : vector to be filled by element. If nullptr on input,
                              the controlling method does not epxect the element
                              to fill this vector
      \param elevec2 (out)  : vector to be filled by element. If nullptr on input,
                              the controlling method does not epxect the element
                              to fill this vector
      \param elevec3 (out)  : vector to be filled by element. If nullptr on input,
                              the controlling method does not epxect the element
                              to fill this vector
      \return 0 if successful, negative otherwise
      */
      int Evaluate(Teuchos::ParameterList& params, DRT::Discretization& discretization,
          std::vector<int>& lm, CORE::LINALG::SerialDenseMatrix& elemat1,
          CORE::LINALG::SerialDenseMatrix& elemat2, CORE::LINALG::SerialDenseVector& elevec1,
          CORE::LINALG::SerialDenseVector& elevec2,
          CORE::LINALG::SerialDenseVector& elevec3) override;

      //!@}

      //! @name Evaluate methods

      //!@}

      /*!
      \brief Return the location vector of this element

      The method computes degrees of freedom this element adresses.
      Degree of freedom ordering is as follows:<br>
      First all degrees of freedom of adjacent nodes are numbered in
      local nodal order, then the element internal degrees of freedom are
      given if present.<br>
      If a derived element has to use a different ordering scheme,
      it is welcome to overload this method as the assembly routines actually
      don't care as long as matrices and vectors evaluated by the element
      match the ordering, which is implicitly assumed.<br>
      Length of the output vector matches number of degrees of freedom
      exactly.<br>
      This version is intended to fill the LocationArray with the dofs
      the element will assemble into. In the standard case these dofs are
      the dofs of the element itself. For some special conditions (e.g.
      the weak dirichlet boundary condtion) a surface element will assemble
      into the dofs of a volume element.<br>

      \note The degrees of freedom returned are not necessarily only nodal dofs.
            Depending on the element implementation, output might also include
            element dofs.

      \param dis (in)        : the discretization this element belongs to
      \param la (out)        : location data for all dofsets of the discretization
      \param doDirichlet (in): whether to get the Dirichlet flags
      \param condstring (in) : Name of condition to be evaluated
      \param params (in)     : List of parameters for use at element level
      */
      void LocationVector(const Discretization& dis, LocationArray& la, bool doDirichlet,
          const std::string& condstring, Teuchos::ParameterList& params) const override;

     private:
      // don't want = operator
      FluidPoroBoundary& operator=(const FluidPoroBoundary& old);
    };

  }  // namespace ELEMENTS
}  // namespace DRT

BACI_NAMESPACE_CLOSE

#endif  // FLUID_ELE_PORO_H
