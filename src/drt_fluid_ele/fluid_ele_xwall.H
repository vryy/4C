/*-----------------------------------------------------------*/
/*!

\brief Implementation of enrichment-wall fluid elements.
In addition to that, it contains the interface between element call
and Gauss point loop (depending on the fluid implementation)
as well as some additional service routines (for the evaluation
of errors, turbulence statistics etc.)

\maintainer Martin Kronbichler

\level 2

*/
/*-----------------------------------------------------------*/

#ifndef DRT_FLUID_ELE_XWALL_H
#define DRT_FLUID_ELE_XWALL_H

#include "Teuchos_RCP.hpp"

#include "fluid_ele.H"


namespace DRT
{
  namespace ELEMENTS
  {
    class FluidXWallType : public FluidType
    {
     public:
      std::string Name() const { return "FluidXWallType"; }

      static FluidXWallType& Instance();

      virtual DRT::ParObject* Create(const std::vector<char>& data);

      virtual Teuchos::RCP<DRT::Element> Create(
          const std::string eletype, const std::string eledistype, const int id, const int owner);

      virtual Teuchos::RCP<DRT::Element> Create(const int id, const int owner);

      virtual void NodalBlockInformation(Element* dwele, int& numdf, int& dimns, int& nv, int& np);

      virtual void ComputeNullSpace(DRT::Discretization& dis, std::vector<double>& ns,
          const double* x0, int numdf, int dimns);

      virtual void SetupElementDefinition(
          std::map<std::string, std::map<std::string, DRT::INPUT::LineDefinition>>& definitions);

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
      virtual DRT::Element* Clone() const;

      /*!
      \brief Get number of degrees of freedom of a certain node
             (implements pure virtual DRT::Element)

      The element decides how many degrees of freedom its nodes must have.
      As this may vary along a simulation, the element can redecide the
      number of degrees of freedom per node along the way for each of it's nodes
      separately.
      */
      virtual int NumDofPerNode(const DRT::Node& node) const
      {
        // number of Dof's is fluid-specific.
        const int nsd = DRT::UTILS::getDimension(distype_);
        if (nsd > 1)
          return 2 * nsd + 2;  // The enrichment dofs have to be a multiple of the usual dofs for
                               // the nullspace
        else
          dserror("1D Fluid elements are not supported");

        return 0;
      }

      virtual DRT::ElementType& ElementType() const { return FluidXWallType::Instance(); }

      /*!
      \brief Return value how expensive it is to evaluate this element

      \param double (out): cost to evaluate this element
      */
      // the standard element is 10.0
      virtual double EvaluationCost() { return 50.0; }

      virtual int UniqueParObjectId() const
      {
        return FluidXWallType::Instance().UniqueParObjectId();
      }

      /*!
      \brief Get vector of Teuchos::RCPs to the lines of this element
      */
      virtual std::vector<Teuchos::RCP<DRT::Element>> Lines();

      /*!
      \brief Get vector of Teuchos::RCPs to the surfaces of this element
      */
      virtual std::vector<Teuchos::RCP<DRT::Element>> Surfaces();

      /*!
      \brief Get vector of Teuchos::RCPs to the volume of this element
      */
      virtual std::vector<Teuchos::RCP<DRT::Element>> Volumes();

      //@}

      //! @name Acess methods

      /*!
      \brief Print this element
      */
      void Print(std::ostream& os) const;

      /*!
      \brief Destructor

      */
      virtual ~FluidXWall();


     protected:
      // don't want = operator
      FluidXWall& operator=(const FluidXWall& old);

    };  // class Fluid


    class FluidXWallBoundaryType : public FluidBoundaryType
    {
     public:
      std::string Name() const { return "FluidXWallBoundaryType"; }

      Teuchos::RCP<DRT::Element> Create(const int id, const int owner);

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
      DRT::Element* Clone() const;

      /*!
      \brief Return unique ParObject id

      every class implementing ParObject needs a unique id defined at the
      top of the drt_parobject.H file.
      */
      virtual int UniqueParObjectId() const
      {
        return FluidXWallBoundaryType::Instance().UniqueParObjectId();
      }

      /*!
      \brief Destructor

      */
      virtual ~FluidXWallBoundary(){};

      //@}

      //! @name Acess methods

      /*!
      \brief Print this element
      */
      void Print(std::ostream& os) const;

      virtual DRT::ElementType& ElementType() const { return FluidXWallBoundaryType::Instance(); }

      //@}

      //! @name Evaluation

      /*!
      \brief Evaluate element

      \param params (in/out): ParameterList for communication between control routine
                              and elements
      \param elemat1 (out)  : matrix to be filled by element. If NULL on input,
                              the controling method does not epxect the element to fill
                              this matrix.
      \param elemat2 (out)  : matrix to be filled by element. If NULL on input,
                              the controling method does not epxect the element to fill
                              this matrix.
      \param elevec1 (out)  : vector to be filled by element. If NULL on input,
                              the controlling method does not epxect the element
                              to fill this vector
      \param elevec2 (out)  : vector to be filled by element. If NULL on input,
                              the controlling method does not epxect the element
                              to fill this vector
      \param elevec3 (out)  : vector to be filled by element. If NULL on input,
                              the controlling method does not epxect the element
                              to fill this vector
      \return 0 if successful, negative otherwise
      */
      //  int Evaluate(Teuchos::ParameterList&   params,
      //               DRT::Discretization&      discretization,
      //               std::vector<int>&         lm,
      //               Epetra_SerialDenseMatrix& elemat1,
      //               Epetra_SerialDenseMatrix& elemat2,
      //               Epetra_SerialDenseVector& elevec1,
      //               Epetra_SerialDenseVector& elevec2,
      //               Epetra_SerialDenseVector& elevec3);

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
      virtual void LocationVector(const Discretization& dis, LocationArray& la, bool doDirichlet,
          const std::string& condstring, Teuchos::ParameterList& params) const;

     private:
      // don't want = operator
      FluidXWallBoundary& operator=(const FluidXWallBoundary& old);

    };  // class FluidXWallBoundary

  }  // namespace ELEMENTS
}  // namespace DRT



#endif  // #ifndef FLUIDXWall_H
