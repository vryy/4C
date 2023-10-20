/*----------------------------------------------------------------------*/
/*! \file

\brief A C++ wrapper for the electromagnetic element

This file contains the element-specific service routines such as
Pack, Unpack, NumDofPerNode etc.


\level 2

*/
/*----------------------------------------------------------------------*/

#ifndef BACI_ELEMAG_ELE_H
#define BACI_ELEMAG_ELE_H

#include "baci_discretization_fem_general_utils_gausspoints.H"
#include "baci_discretization_fem_general_utils_integration.H"
#include "baci_discretization_fem_general_utils_local_connectivity_matrices.H"
#include "baci_discretization_fem_general_utils_polynomial.H"
#include "baci_lib_elementtype.H"
#include "baci_linalg_serialdensematrix.H"

#include <Teuchos_RCP.hpp>

const bool trac_with_pml_elemag = false;

namespace DRT
{
  // forward declarations
  class Discretization;

  namespace ELEMENTS
  {
    class ElemagIntFace;

    class ElemagType : public DRT::ElementType
    {
     public:
      std::string Name() const override { return "ElemagType"; }

      static ElemagType& Instance();

      DRT::ParObject* Create(const std::vector<char>& data) override;

      Teuchos::RCP<DRT::Element> Create(const std::string eletype, const std::string eledistype,
          const int id, const int owner) override;

      Teuchos::RCP<DRT::Element> Create(const int id, const int owner) override;

      void NodalBlockInformation(
          DRT::Element* dwele, int& numdf, int& dimns, int& nv, int& np) override;

      CORE::LINALG::SerialDenseMatrix ComputeNullSpace(
          DRT::Node& node, const double* x0, const int numdof, const int dimnsp) override;

      void SetupElementDefinition(
          std::map<std::string, std::map<std::string, DRT::INPUT::LineDefinition>>& definitions)
          override;

      /// pre-evaluation
      void PreEvaluate(DRT::Discretization& dis, Teuchos::ParameterList& p,
          Teuchos::RCP<CORE::LINALG::SparseOperator> systemmatrix1,
          Teuchos::RCP<CORE::LINALG::SparseOperator> systemmatrix2,
          Teuchos::RCP<Epetra_Vector> systemvector1, Teuchos::RCP<Epetra_Vector> systemvector2,
          Teuchos::RCP<Epetra_Vector> systemvector3) override;

     private:
      static ElemagType instance_;
    };

    /*!
    \brief A C++ wrapper for the electromagnetic element
    */
    class Elemag : public DRT::Element
    {
     public:
      //@}
      //! @name constructors and destructors and related methods

      /*!
      \brief standard constructor
      */
      Elemag(int id,  ///< A unique global id
          int owner   ///< ???
      );

      /*!
      \brief Copy Constructor

      Makes a deep copy of a Element

      */
      Elemag(const Elemag& old);

      /*!
      \brief Deep copy this instance and return pointer to the copy

      The Clone() method is used from the virtual base class Element in cases
      where the type of the derived class is unknown and a copy-ctor is needed

      */
      DRT::Element* Clone() const override;

      /*!
      \brief Destructor
      */
      ~Elemag() override;

      /*!
      \brief Get shape type of element
      */
      DiscretizationType Shape() const override { return distype_; };

      /*!
      \brief set discretization type of element
      */
      virtual void SetDisType(DRT::Element::DiscretizationType shape)
      {
        distype_ = shape;
        return;
      };

      /*!
      \brief Return number of lines of this element
      */
      int NumLine() const override { return CORE::DRT::UTILS::getNumberOfElementLines(distype_); }

      /*!
      \brief Return number of surfaces of this element
      */
      int NumSurface() const override
      {
        return CORE::DRT::UTILS::getNumberOfElementSurfaces(distype_);
      }

      /*!
      \brief Return number of volumes of this element (always 1)
      */
      int NumVolume() const override
      {
        return CORE::DRT::UTILS::getNumberOfElementVolumes(distype_);
      }

      /*!
      \brief Get vector of Teuchos::RCPs to the lines of this element
      */
      std::vector<Teuchos::RCP<DRT::Element>> Lines() override;

      /*!
      \brief Get vector of Teuchos::RCPs to the surfaces of this element
      */
      std::vector<Teuchos::RCP<DRT::Element>> Surfaces() override;

      /*!
      \brief Get vector of Teuchos::RCPs to the volume of this element
      */
      std::vector<Teuchos::RCP<DRT::Element>> Volumes() override;

      /*!
      \brief Get Teuchos::RCP to the internal face adjacent to this element as master element and
      the parent_slave element
      */
      Teuchos::RCP<DRT::Element> CreateFaceElement(
          DRT::Element* parent_slave,  //!< parent slave fluid3 element
          int nnode,                   //!< number of surface nodes
          const int* nodeids,          //!< node ids of surface element
          DRT::Node** nodes,           //!< nodes of surface element
          const int lsurface_master,   //!< local surface number w.r.t master parent element
          const int lsurface_slave,    //!< local surface number w.r.t slave parent element
          const std::vector<int>& localtrafomap  //! local trafo map
          ) override;

      /*!
      \brief Return unique ParObject id

      every class implementing ParObject needs a unique id defined at the
      top of this file.
      */
      int UniqueParObjectId() const override { return ElemagType::Instance().UniqueParObjectId(); }
      /*!
      \brief Pack this class so it can be communicated

      \ref Pack and \ref Unpack are used to communicate this element

      */
      void Pack(DRT::PackBuffer& data) const override;

      /*!
      \brief Unpack data from a char vector into this class

      \ref Pack and \ref Unpack are used to communicate this element
      */
      void Unpack(const std::vector<char>& data) override;

      //@}

      //! @name Geometry related methods

      //@}

      //! @name Acess methods

      /*!
      \brief Get number of degrees of freedom per node

      HDG element: No dofs are associated with nodes
      */
      int NumDofPerNode(const DRT::Node&) const override { return 0; }

      /*!
      \brief Get number of degrees of freedom per face
      This implementation is problem related. In the case of a fluid flow
      NumDofPerFace() reuturns the number of velocity components (one for each spatial dimension)
      times the number of nodes per face (given by NumDofPerComponent()).
      The number is nsd_ - 1 because the trace values are contained in the faces and
      therefore they have one component less than the spatial dimnsions of the problem.
      */
      int NumDofPerFace(const unsigned face) const override
      {
        return (CORE::DRT::UTILS::getDimension(distype_) - 1) *
               NumDofPerComponent(face);  // seems redundant, but this is different for
                                          // fluid_ele_hdg and Elemag_ele!
      }

      /*!
      \brief Get number of dofs per component per face
      It gives the number of nodes the unknowns have to be computed in.

      EXAMPLE: For a HEX8 NumDofPerComponent() will return 4 while for a TET4 will return 3.
      */
      int NumDofPerComponent(const unsigned face) const override
      {
        return CORE::DRT::UTILS::getBasisSize(CORE::DRT::UTILS::getEleFaceShapeType(distype_),
            (this->Faces()[face])->Degree(), completepol_);
      }

      /*!
      \brief Returns the number of unknowns to be computed in each element.

      EXAMPLE: For a HEX8 element and a 3D problem the return value will be 48.
      */
      virtual int NumDofPerElementAuxiliary() const
      {
        return (CORE::DRT::UTILS::getDimension(distype_) * 2) *
               CORE::DRT::UTILS::getBasisSize(distype_, degree_, completepol_);
      }

      /// Number of degrees of freedom per element
      int NumDofPerElement() const override { return 0; }

      /*!
       \brief Returns the degree of the element
       */
      int Degree() const override { return degree_; }

      /*!
       \brief Allows to set the element degree
       */
      void SetDegree(int degree)
      {
        degree_ = degree;
        return;
      }

      /*!
       \brief Returns the degree of the element
       */
      int UsesCompletePolynomialSpace() const { return completepol_; }

      /*!
      \brief Print this element
      */
      void Print(std::ostream& os) const override;

      DRT::ElementType& ElementType() const override { return ElemagType::Instance(); }

      //@}

      //! @name Input and Creation

      /*!
      \brief Read input for this element
      */
      bool ReadElement(const std::string& eletype, const std::string& distype,
          DRT::INPUT::LineDefinition* linedef) override;

      //@}

      //! @name Evaluation

      /*!
      \brief Evaluate an element, that is, call the element routines to evaluate Elemag
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


      /*!
      \brief Evaluate Neumann boundary condition

      \param params (in/out)    : ParameterList for communication between control routine
                                  and elements
      \param discretization (in): reference to the underlying discretization
      \param condition (in)     : condition to be evaluated
      \param lm (in)            : location vector of this element
      \param elevec1 (out)      : vector to be filled by element. If nullptr on input,

      \return 0 if successful, negative otherwise
      */
      int EvaluateNeumann(Teuchos::ParameterList& params, DRT::Discretization& discretization,
          DRT::Condition& condition, std::vector<int>& lm, CORE::LINALG::SerialDenseVector& elevec1,
          CORE::LINALG::SerialDenseMatrix* elemat1 = nullptr) override;

      //@}

      /// Trace values per element
      CORE::LINALG::SerialDenseVector elenodeTrace2d_;

      // state vectors which do not have to be stored as global distributed vector, since
      // we are in a discontinuous setting here
      //! @name interior stored states

      /// Electric field current step
      CORE::LINALG::SerialDenseVector eleinteriorElectric_;

      /// Electric field step n-1 for higher orders BDF
      CORE::LINALG::SerialDenseVector eleinteriorElectricnm1_;

      /// Electric field step n-2 for higher orders BDF
      CORE::LINALG::SerialDenseVector eleinteriorElectricnm2_;

      /// Electric field step n-3 for higher orders BDF
      CORE::LINALG::SerialDenseVector eleinteriorElectricnm3_;

      /// Post-processed solution
      CORE::LINALG::SerialDenseVector eleinteriorElectricPost_;

      /// Magnetic field current step
      CORE::LINALG::SerialDenseVector eleinteriorMagnetic_;

      /// Magnetic field step n-1 for higher orders BDF
      CORE::LINALG::SerialDenseVector eleinteriorMagneticnm1_;

      /// Magnetic field step n-2 for higher orders BDF
      CORE::LINALG::SerialDenseVector eleinteriorMagneticnm2_;

      /// Magnetic field step n-3 for higher orders BDF
      CORE::LINALG::SerialDenseVector eleinteriorMagneticnm3_;

      //@}

     protected:
      //! discretization type
      DRT::Element::DiscretizationType distype_;
      DRT::Element::DiscretizationType facedis_;

      /// stores the degree of the element
      unsigned char degree_;

      /// stores the polynomial type (tensor product or complete polynomial)
      bool completepol_;

      int fixedPower(const int number, const int pow) const
      {
        int result = 1;
        for (int i = 0; i < pow; ++i) result *= number;
        return result;
      }

     private:
      // don't want = operator
      Elemag& operator=(const Elemag& old);

    };  // class Elemag


    //=======================================================================
    //=======================================================================
    //=======================================================================
    //=======================================================================



    //=======================================================================
    //=======================================================================
    //=======================================================================
    //=======================================================================


    /*!
    \brief An element representing a boundary element of an Elemag element

    \note This is a pure Neumann boundary condition element. It's only
          purpose is to evaluate surface Neumann boundary conditions that might be
          adjacent to a parent Elemag element. It therefore does not implement
          the DRT::Element::Evaluate method and does not have its own ElementRegister class.

    */

    class ElemagBoundaryType : public DRT::ElementType
    {
     public:
      std::string Name() const override { return "ElemagBoundaryType"; }

      static ElemagBoundaryType& Instance();

      Teuchos::RCP<DRT::Element> Create(const int id, const int owner) override;

      void NodalBlockInformation(Element* dwele, int& numdf, int& dimns, int& nv, int& np) override
      {
      }

      CORE::LINALG::SerialDenseMatrix ComputeNullSpace(
          DRT::Node& node, const double* x0, const int numdof, const int dimnsp) override
      {
        CORE::LINALG::SerialDenseMatrix nullspace;
        dserror("method ComputeNullSpace not implemented");
        return nullspace;
      }

     private:
      static ElemagBoundaryType instance_;
    };

    /// class ElemagBoundary

    class ElemagBoundary : public DRT::FaceElement
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
      \param parent: The parent Elemag element of this surface
      \param lsurface: the local surface number of this surface w.r.t. the parent element
      */
      ElemagBoundary(int id, int owner, int nnode, const int* nodeids, DRT::Node** nodes,
          DRT::ELEMENTS::Elemag* parent, const int lsurface);

      /*!
      \brief Copy Constructor

      Makes a deep copy of a Element

      */
      ElemagBoundary(const ElemagBoundary& old);

      /*!
      \brief Deep copy this instance of an element and return pointer to the copy

      The Clone() method is used from the virtual base class Element in cases
      where the type of the derived class is unknown and a copy-ctor is needed

      */
      DRT::Element* Clone() const override;

      /*!
      \brief Get shape type of element
      */
      DiscretizationType Shape() const override;

      /*!
      \brief Return number of lines of this element
      */
      int NumLine() const override { return CORE::DRT::UTILS::getNumberOfElementLines(Shape()); }

      /*!
      \brief Return number of surfaces of this element
      */
      int NumSurface() const override
      {
        return CORE::DRT::UTILS::getNumberOfElementSurfaces(Shape());
      }

      /*!
      \brief Get vector of Teuchos::RCPs to the lines of this element

      */
      std::vector<Teuchos::RCP<DRT::Element>> Lines() override;

      /*!
      \brief Get vector of Teuchos::RCPs to the surfaces of this element
      */

      std::vector<Teuchos::RCP<DRT::Element>> Surfaces() override;

      /*!
      \brief Return unique ParObject id

      every class implementing ParObject needs a unique id defined at the
      top of the parobject.H file.
      */
      int UniqueParObjectId() const override
      {
        return ElemagBoundaryType::Instance().UniqueParObjectId();
      }

      /*!
      \brief Pack this class so it can be communicated

      \ref Pack and \ref Unpack are used to communicate this element

      */
      void Pack(DRT::PackBuffer& data) const override;

      /*!
      \brief Unpack data from a char vector into this class

      \ref Pack and \ref Unpack are used to communicate this element

      */
      void Unpack(const std::vector<char>& data) override;

      /*!
      \brief Destructor

      */
      ~ElemagBoundary() override;

      //@}

      //! @name Acess methods


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
        return ParentElement()->NumDofPerNode(node);
      }

      /// Number of degrees of freedom per element
      int NumDofPerElement() const override { return 0; }

      /*!
      \brief Print this element
      */
      void Print(std::ostream& os) const override;

      DRT::ElementType& ElementType() const override { return ElemagBoundaryType::Instance(); }

      //@}

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

      //@}

      //! @name Evaluate methods

      /*!
      \brief Evaluate Neumann boundary condition

      \param params (in/out)    : ParameterList for communication between control routine
                                  and elements
      \param discretization (in): reference to the underlying discretization
      \param condition (in)     : condition to be evaluated
      \param lm (in)            : location vector of this element
      \param elevec1 (out)      : vector to be filled by element. If nullptr on input,

      \return 0 if successful, negative otherwise
      */
      int EvaluateNeumann(Teuchos::ParameterList& params, DRT::Discretization& discretization,
          DRT::Condition& condition, std::vector<int>& lm, CORE::LINALG::SerialDenseVector& elevec1,
          CORE::LINALG::SerialDenseMatrix* elemat1 = nullptr) override;

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

      \note The degrees of freedom returned are not neccessarily only nodal dofs.
            Depending on the element implementation, output might also include
            element dofs.

      \param dis (in)      : the discretization this element belongs to
      \param la (out)      : location data for all dofsets of the discretization
      \param doDirichlet (in): whether to get the Dirichlet flags
      \param condstring (in): Name of condition to be evaluated
      \param condstring (in):  List of parameters for use at element level
      */
      void LocationVector(const Discretization& dis, LocationArray& la, bool doDirichlet,
          const std::string& condstring, Teuchos::ParameterList& params) const override;

      int degree_;

      /*!
     \brief Allows to set the element degree
     */
      void SetDegree(int degree)
      {
        degree_ = degree;
        return;
      }


     private:
      // don't want = operator
      ElemagBoundary& operator=(const ElemagBoundary& old);

    };  // class ElemagBoundary


    //=======================================================================
    //=======================================================================
    //=======================================================================
    //=======================================================================



    //=======================================================================
    //=======================================================================
    //=======================================================================
    //=======================================================================


    /*!
    \brief An element representing an internal face element between two Elemag elements
    */
    class ElemagIntFaceType : public DRT::ElementType
    {
     public:
      std::string Name() const override { return "ElemagIntFaceType"; }

      static ElemagIntFaceType& Instance();

      Teuchos::RCP<DRT::Element> Create(const int id, const int owner) override;

      void NodalBlockInformation(Element* dwele, int& numdf, int& dimns, int& nv, int& np) override
      {
      }

      CORE::LINALG::SerialDenseMatrix ComputeNullSpace(
          DRT::Node& node, const double* x0, const int numdof, const int dimnsp) override
      {
        CORE::LINALG::SerialDenseMatrix nullspace;
        dserror("method ComputeNullSpace not implemented");
        return nullspace;
      }

     private:
      static ElemagIntFaceType instance_;
    };


    // class ElemagIntFace

    class ElemagIntFace : public DRT::FaceElement
    {
     public:
      //! @name Constructors and destructors and related methods

      //! number of space dimensions
      /*!
      \brief Standard Constructor

      \param id: A unique global id
      \param owner: Processor owning this surface
      \param nnode: Number of nodes attached to this element
      \param nodeids: global ids of nodes attached to this element
      \param nodes: the discretizations map of nodes to build ptrs to nodes from
      \param master_parent: The master parent Elemag element of this surface
      \param slave_parent: The slave parent Elemag element of this surface
      \param lsurface_master: the local surface number of this surface w.r.t. the master parent
      element \param lsurface_slave: the local surface number of this surface w.r.t. the slave
      parent element \param localtrafomap: transformation map between the local coordinate systems
      of the face w.r.t the master parent element's face's coordinate system and the slave element's
      face's coordinate system
      */
      ElemagIntFace(int id, int owner, int nnode, const int* nodeids, DRT::Node** nodes,
          DRT::ELEMENTS::Elemag* parent_master, DRT::ELEMENTS::Elemag* parent_slave,
          const int lsurface_master, const int lsurface_slave,
          const std::vector<int> localtrafomap);

      /*!
      \brief Copy Constructor

      Makes a deep copy of a Element
      */
      ElemagIntFace(const ElemagIntFace& old);

      /*!
      \brief Deep copy this instance of an element and return pointer to the copy

      The Clone() method is used from the virtual base class Element in cases
      where the type of the derived class is unknown and a copy-ctor is needed

      */
      DRT::Element* Clone() const override;

      /*!
      \brief Get shape type of element
      */
      DiscretizationType Shape() const override;

      /*!
      \brief Return number of lines of this element
      */
      int NumLine() const override { return CORE::DRT::UTILS::getNumberOfElementLines(Shape()); }

      /*!
      \brief Return number of surfaces of this element
      */
      int NumSurface() const override
      {
        return CORE::DRT::UTILS::getNumberOfElementSurfaces(Shape());
      }

      /*!
      \brief Get vector of Teuchos::RCPs to the lines of this element
      */
      std::vector<Teuchos::RCP<DRT::Element>> Lines() override;

      /*!
      \brief Return unique ParObject id

      every class implementing ParObject needs a unique id defined at the
      top of the parobject.H file.
      */
      std::vector<Teuchos::RCP<DRT::Element>> Surfaces() override;

      /*!
      \brief Return unique ParObject id

      every class implementing ParObject needs a unique id defined at the
      top of the parobject.H file.
      */
      int UniqueParObjectId() const override
      {
        return ElemagIntFaceType::Instance().UniqueParObjectId();
      }

      /*!
      \brief Pack this class so it can be communicated

      \ref Pack and \ref Unpack are used to communicate this element
      */
      void Pack(DRT::PackBuffer& data) const override;

      /*!
      \brief Unpack data from a char vector into this class

      \ref Pack and \ref Unpack are used to communicate this element
      */
      void Unpack(const std::vector<char>& data) override;

      /*!
      \brief Destructor
      */
      ~ElemagIntFace() override;

      //@}

      //! @name Acess methods

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
        return std::max(
            ParentMasterElement()->NumDofPerNode(node), ParentSlaveElement()->NumDofPerNode(node));
      }

      /*!
      \brief Get number of degrees of freedom per element
             (implements pure virtual DRT::Element)

      The element decides how many element degrees of freedom it has.
      It can redecide along the way of a simulation.

      \note Element degrees of freedom mentioned here are dofs that are visible
            at the level of the total system of equations. Purely internal
            element dofs that are condensed internally should NOT be considered.
      */
      int NumDofPerElement() const override { return 0; }

      /*!
       \brief Returns the degree of the element
       */
      int Degree() const override { return degree_; }

      /*!
      \brief Allows to set the element degree
      */
      void SetDegree(int degree)
      {
        degree_ = degree;
        return;
      }

      /*!
      \brief create the location vector for patch of master and slave element

      \note All dofs shared by master and slave element are contained only once. Dofs from interface
      nodes are also included.
      */
      void PatchLocationVector(DRT::Discretization& discretization,  ///< discretization
          std::vector<int>& nds_master,        ///< nodal dofset w.r.t master parent element
          std::vector<int>& nds_slave,         ///< nodal dofset w.r.t slave parent element
          std::vector<int>& patchlm,           ///< local map for gdof ids for patch of elements
          std::vector<int>& master_lm,         ///< local map for gdof ids for master element
          std::vector<int>& slave_lm,          ///< local map for gdof ids for slave element
          std::vector<int>& face_lm,           ///< local map for gdof ids for face element
          std::vector<int>& lm_masterToPatch,  ///< local map between lm_master and lm_patch
          std::vector<int>& lm_slaveToPatch,   ///< local map between lm_slave and lm_patch
          std::vector<int>& lm_faceToPatch,    ///< local map between lm_face and lm_patch
          std::vector<int>&
              lm_masterNodeToPatch,  ///< local map between master nodes and nodes in patch
          std::vector<int>&
              lm_slaveNodeToPatch  ///< local map between slave nodes and nodes in patch
      );

      /*!
      \brief Print this element
      */
      void Print(std::ostream& os) const override;

      DRT::ElementType& ElementType() const override { return ElemagIntFaceType::Instance(); }

      //@}

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

      //@}

      //! @name Evaluate methods

      /*!
      \brief Evaluate Neumann boundary condition

      \param params (in/out)    : ParameterList for communication between control routine
                                  and elements
      \param discretization (in): reference to the underlying discretization
      \param condition (in)     : condition to be evaluated
      \param lm (in)            : location vector of this element
      \param elevec1 (out)      : vector to be filled by element. If nullptr on input,

      \return 0 if successful, negative otherwise
      */
      int EvaluateNeumann(Teuchos::ParameterList& params, DRT::Discretization& discretization,
          DRT::Condition& condition, std::vector<int>& lm, CORE::LINALG::SerialDenseVector& elevec1,
          CORE::LINALG::SerialDenseMatrix* elemat1 = nullptr) override;

      /*!
      \brief return the master parent Elemag element
      */
      DRT::ELEMENTS::Elemag* ParentMasterElement() const
      {
        DRT::Element* parent = this->DRT::FaceElement::ParentMasterElement();
        // make sure the static cast below is really valid
        dsassert(dynamic_cast<DRT::ELEMENTS::Elemag*>(parent) != nullptr,
            "Master element is no Elemag element");
        return static_cast<DRT::ELEMENTS::Elemag*>(parent);
      }

      /*!
      \brief return the slave parent Elemag element
      */
      DRT::ELEMENTS::Elemag* ParentSlaveElement() const
      {
        DRT::Element* parent = this->DRT::FaceElement::ParentSlaveElement();
        // make sure the static cast below is really valid
        dsassert(dynamic_cast<DRT::ELEMENTS::Elemag*>(parent) != nullptr,
            "Slave element is no Elemag element");
        return static_cast<DRT::ELEMENTS::Elemag*>(parent);
      }

      //@}

     private:
      // don't want = operator
      ElemagIntFace& operator=(const ElemagIntFace& old);

      // degree of this face element
      int degree_;

    };  // class ElemagIntFace



  }  // namespace ELEMENTS
}  // namespace DRT


#endif  // ELEMAG_ELE_H