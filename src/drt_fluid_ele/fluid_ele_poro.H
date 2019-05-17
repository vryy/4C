/*-----------------------------------------------------------*/
/*!

\brief A C++ wrapper for the poro fluid element

This file contains the element-specific service routines such as
Pack, Unpack, NumDofPerNode etc.

\maintainer Christoph Ager

\level 2

*/
/*-----------------------------------------------------------*/


#ifndef FLUID_ELE_PORO_H_
#define FLUID_ELE_PORO_H_

#include "../drt_inpar/inpar_structure.H"
#include "fluid_ele.H"

namespace DRT
{
  class Discretization;

  namespace ELEMENTS
  {
    class FluidPoroEleType : public FluidType
    {
     public:
      std::string Name() const { return "FluidPoroEleType"; }

      static FluidPoroEleType& Instance();

      virtual DRT::ParObject* Create(const std::vector<char>& data);

      virtual Teuchos::RCP<DRT::Element> Create(
          const std::string eletype, const std::string eledistype, const int id, const int owner);

      virtual Teuchos::RCP<DRT::Element> Create(const int id, const int owner);

      virtual void SetupElementDefinition(
          std::map<std::string, std::map<std::string, DRT::INPUT::LineDefinition>>& definitions);

      /// pre-evaluation
      virtual void PreEvaluate(DRT::Discretization& dis, Teuchos::ParameterList& p,
          Teuchos::RCP<LINALG::SparseOperator> systemmatrix1,
          Teuchos::RCP<LINALG::SparseOperator> systemmatrix2,
          Teuchos::RCP<Epetra_Vector> systemvector1, Teuchos::RCP<Epetra_Vector> systemvector2,
          Teuchos::RCP<Epetra_Vector> systemvector3);

     private:
      static FluidPoroEleType instance_;
    };

    /*!
    \brief A C++ wrapper for the fluid element
    */
    class FluidPoro : public virtual Fluid
    {
     public:
      //@}
      //! @name constructors and destructors and related methods

      /*!
      \brief standard constructor
      */
      FluidPoro(int id,  ///< A unique global id
          int owner      ///< number of owning processor
      );

      /*!
      \brief Copy Constructor

      Makes a deep copy of a Element

      */
      FluidPoro(const FluidPoro& old);

      /*!
      \brief Destructor

      */
      virtual ~FluidPoro(){};

      /*!
      \brief Deep copy this instance of fluid and return pointer to the copy

      The Clone() method is used from the virtual base class Element in cases
      where the type of the derived class is unknown and a copy-ctor is needed

      */
      DRT::Element* Clone() const;

      /*!
      \brief Return unique ParObject id

      every class implementing ParObject needs a unique id defined at the
      top of this file.
      */
      virtual int UniqueParObjectId() const
      {
        return FluidPoroEleType::Instance().UniqueParObjectId();
      }

      /*!
      \brief Pack this class so it can be communicated

      \ref Pack and \ref Unpack are used to communicate this element

      */
      virtual void Pack(DRT::PackBuffer& data) const;

      /*!
      \brief Unpack data from a char vector into this class

      \ref Pack and \ref Unpack are used to communicate this element
      */
      virtual void Unpack(const std::vector<char>& data);

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

      virtual DRT::ElementType& ElementType() const { return FluidPoroEleType::Instance(); }


      void SetKinematicType(INPAR::STR::KinemType kintype) { kintype_ = kintype; };

      INPAR::STR::KinemType KinematicType() const { return kintype_; };

      //@}

      //! @name Evaluation

      /*!
      \brief Evaluate an element, that is, call the element routines to evaluate fluid
      element matrices and vectors or evaluate errors, statistics or updates etc. directly.

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
      int Evaluate(Teuchos::ParameterList& params, DRT::Discretization& discretization,
          std::vector<int>& lm, Epetra_SerialDenseMatrix& elemat1,
          Epetra_SerialDenseMatrix& elemat2, Epetra_SerialDenseVector& elevec1,
          Epetra_SerialDenseVector& elevec2, Epetra_SerialDenseVector& elevec3);

      //@}

     private:
      //! don't want = operator
      FluidPoro& operator=(const FluidPoro& old);

      //! kinematic type
      INPAR::STR::KinemType kintype_;
    };

    //=======================================================================
    //=======================================================================
    //=======================================================================
    //=======================================================================


    /*!
    \brief An element representing a boundary element of a fluid element

    */

    class FluidPoroBoundaryType : public FluidBoundaryType
    {
     public:
      std::string Name() const { return "FluidPoroBoundaryType"; }

      Teuchos::RCP<DRT::Element> Create(const int id, const int owner);

      static FluidPoroBoundaryType& Instance();

     private:
      static FluidPoroBoundaryType instance_;
    };


    // class FluidPoroBoundary

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
      DRT::Element* Clone() const;

      /*!
      \brief Return unique ParObject id

      every class implementing ParObject needs a unique id defined at the
      top of the drt_parobject.H file.
      */
      virtual int UniqueParObjectId() const
      {
        return FluidPoroBoundaryType::Instance().UniqueParObjectId();
      }

      /*!
      \brief Destructor

      */
      virtual ~FluidPoroBoundary(){};

      //@}

      //! @name Acess methods

      /*!
      \brief Print this element
      */
      void Print(std::ostream& os) const;

      virtual DRT::ElementType& ElementType() const { return FluidPoroBoundaryType::Instance(); }

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
      int Evaluate(Teuchos::ParameterList& params, DRT::Discretization& discretization,
          std::vector<int>& lm, Epetra_SerialDenseMatrix& elemat1,
          Epetra_SerialDenseMatrix& elemat2, Epetra_SerialDenseVector& elevec1,
          Epetra_SerialDenseVector& elevec2, Epetra_SerialDenseVector& elevec3);

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
      FluidPoroBoundary& operator=(const FluidPoroBoundary& old);

    };  // class FluidBoundary

  }  // namespace ELEMENTS
}  // namespace DRT

#endif /* FLUID_ELE_PORO_H_ */
