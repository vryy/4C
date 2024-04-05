/*-----------------------------------------------------------*/
/*! \file

\brief Implementation of enrichment-wall fluid elements.
In addition to that, it contains the interface between element call
and Gauss point loop (depending on the fluid implementation)
as well as some additional service routines (for the evaluation
of errors, turbulence statistics etc.)


\level 2

*/
/*-----------------------------------------------------------*/

#ifndef FOUR_C_FLUID_ELE_XWALL_HPP
#define FOUR_C_FLUID_ELE_XWALL_HPP

#include "baci_config.hpp"

#include "baci_fluid_ele.hpp"

#include <Teuchos_RCP.hpp>

BACI_NAMESPACE_OPEN


namespace DRT
{
  namespace ELEMENTS
  {
    class FluidXWallType : public FluidType
    {
     public:
      std::string Name() const override { return "FluidXWallType"; }

      static FluidXWallType& Instance();

      CORE::COMM::ParObject* Create(const std::vector<char>& data) override;

      Teuchos::RCP<DRT::Element> Create(const std::string eletype, const std::string eledistype,
          const int id, const int owner) override;

      Teuchos::RCP<DRT::Element> Create(const int id, const int owner) override;

      void NodalBlockInformation(Element* dwele, int& numdf, int& dimns, int& nv, int& np) override;

      CORE::LINALG::SerialDenseMatrix ComputeNullSpace(
          DRT::Node& node, const double* x0, const int numdof, const int dimnsp) override;

      void SetupElementDefinition(
          std::map<std::string, std::map<std::string, INPUT::LineDefinition>>& definitions)
          override;

     private:
      static FluidXWallType instance_;
    };


    /*!
    \brief A C++ wrapper for the fluid element
    */
    class FluidXWall : public Fluid
    {
     public:
      //! @name constructors and destructors and related methods

      /*!
      \brief standard constructor
      */
      FluidXWall(int id,  ///< A unique global id
          int owner       ///< ???
      );

      /*!
      \brief Copy Constructor

      Makes a deep copy of a Element

      */
      FluidXWall(const FluidXWall& old);

      /*!
      \brief Deep copy this instance of fluid and return pointer to the copy

      The Clone() method is used from the virtual base class Element in cases
      where the type of the derived class is unknown and a copy-ctor is needed

      */
      DRT::Element* Clone() const override;

      /*!
      \brief Get number of degrees of freedom of a certain node
             (implements pure virtual DRT::Element)

      The element decides how many degrees of freedom its nodes must have.
      As this may vary along a simulation, the element can redecide the
      number of degrees of freedom per node along the way for each of it's nodes
      separately.
      */
      int NumDofPerNode(const DRT::Node& node) const override
      {
        // number of Dof's is fluid-specific.
        const int nsd = CORE::FE::getDimension(distype_);
        if (nsd > 1)
          return 2 * nsd + 2;  // The enrichment dofs have to be a multiple of the usual dofs for
                               // the nullspace
        else
          dserror("1D Fluid elements are not supported");

        return 0;
      }

      DRT::ElementType& ElementType() const override { return FluidXWallType::Instance(); }

      /*!
      \brief Return value how expensive it is to evaluate this element

      \param double (out): cost to evaluate this element
      */
      // the standard element is 10.0
      double EvaluationCost() override { return 50.0; }

      int UniqueParObjectId() const override
      {
        return FluidXWallType::Instance().UniqueParObjectId();
      }

      /*!
      \brief Get vector of Teuchos::RCPs to the lines of this element
      */
      std::vector<Teuchos::RCP<DRT::Element>> Lines() override;

      /*!
      \brief Get vector of Teuchos::RCPs to the surfaces of this element
      */
      std::vector<Teuchos::RCP<DRT::Element>> Surfaces() override;

      //@}

      //! @name Acess methods

      /*!
      \brief Print this element
      */
      void Print(std::ostream& os) const override;



     protected:
      // don't want = operator
      FluidXWall& operator=(const FluidXWall& old);

    };  // class Fluid


    class FluidXWallBoundaryType : public FluidBoundaryType
    {
     public:
      std::string Name() const override { return "FluidXWallBoundaryType"; }

      Teuchos::RCP<DRT::Element> Create(const int id, const int owner) override;

      static FluidXWallBoundaryType& Instance();

     private:
      static FluidXWallBoundaryType instance_;
    };


    // class FluidXWallBoundary

    class FluidXWallBoundary : public FluidBoundary
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
      FluidXWallBoundary(int id, int owner, int nnode, const int* nodeids, DRT::Node** nodes,
          DRT::ELEMENTS::Fluid* parent, const int lsurface);

      /*!
      \brief Copy Constructor

      Makes a deep copy of a Element

      */
      FluidXWallBoundary(const FluidXWallBoundary& old);

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
        return FluidXWallBoundaryType::Instance().UniqueParObjectId();
      }

      //@}

      //! @name Acess methods

      /*!
      \brief Print this element
      */
      void Print(std::ostream& os) const override;

      DRT::ElementType& ElementType() const override { return FluidXWallBoundaryType::Instance(); }

      //@}

      //! @name Evaluation

      //@}

      //! @name Evaluate methods

      //@}

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
      FluidXWallBoundary& operator=(const FluidXWallBoundary& old);

    };  // class FluidXWallBoundary

  }  // namespace ELEMENTS
}  // namespace DRT



BACI_NAMESPACE_CLOSE

#endif  // FLUID_ELE_XWALL_H
