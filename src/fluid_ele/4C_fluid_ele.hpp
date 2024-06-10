/*----------------------------------------------------------------------*/
/*! \file

\brief A C++ wrapper for the fluid element

This file contains the element-specific service routines such as
Pack, Unpack, NumDofPerNode etc.

In addition to that, it contains the interface between element call
and Gauss point loop (depending on the fluid implementation)
as well as some additional service routines (for the evaluation
of errors, turbulence statistics etc.).


\level 1

*/
/*----------------------------------------------------------------------*/

#ifndef FOUR_C_FLUID_ELE_HPP
#define FOUR_C_FLUID_ELE_HPP

#include "4C_config.hpp"

#include "4C_fem_general_elementtype.hpp"
#include "4C_fem_general_utils_local_connectivity_matrices.hpp"
#include "4C_linalg_serialdensematrix.hpp"

#include <Teuchos_RCP.hpp>

FOUR_C_NAMESPACE_OPEN

// forward declarations
namespace FLD
{
  class TDSEleData;
}

namespace Core::FE
{
  class Discretization;
}  // namespace Core::FE

namespace Discret
{
  namespace ELEMENTS
  {
    class FluidEleParameter;
    class FluidEleParameterStd;
    class FluidEleParameterPoro;
    class FluidEleParameterTimInt;
    class FluidBoundary;
    class FluidIntFace;

    class FluidType : public Core::Elements::ElementType
    {
     public:
      std::string Name() const override { return "FluidType"; }

      static FluidType& Instance();

      Core::Communication::ParObject* Create(const std::vector<char>& data) override;

      Teuchos::RCP<Core::Elements::Element> Create(const std::string eletype,
          const std::string eledistype, const int id, const int owner) override;

      Teuchos::RCP<Core::Elements::Element> Create(const int id, const int owner) override;

      void nodal_block_information(
          Core::Elements::Element* dwele, int& numdf, int& dimns, int& nv, int& np) override;

      Core::LinAlg::SerialDenseMatrix ComputeNullSpace(
          Core::Nodes::Node& node, const double* x0, const int numdof, const int dimnsp) override;

      void setup_element_definition(
          std::map<std::string, std::map<std::string, Input::LineDefinition>>& definitions)
          override;

      /// pre-evaluation
      void pre_evaluate(Core::FE::Discretization& dis, Teuchos::ParameterList& p,
          Teuchos::RCP<Core::LinAlg::SparseOperator> systemmatrix1,
          Teuchos::RCP<Core::LinAlg::SparseOperator> systemmatrix2,
          Teuchos::RCP<Epetra_Vector> systemvector1, Teuchos::RCP<Epetra_Vector> systemvector2,
          Teuchos::RCP<Epetra_Vector> systemvector3) override;

     private:
      static FluidType instance_;
    };

    /*!
    \brief A C++ wrapper for the fluid element
    */
    class Fluid : public Core::Elements::Element
    {
     public:
      /*!
       \brief Enrichment type
      */
      enum EnrichmentType
      {
        none,  // no enrichment
        xwall  // xwall: additional virtual nodes, every other is a virtual one
      };

      //! @name friends
      friend class FluidBoundary;
      friend class FluidIntFace;

      friend class FluidSystemEvaluator;

      //@}
      //! @name constructors and destructors and related methods

      /*!
      \brief standard constructor
      */
      Fluid(int id,  ///< A unique global id
          int owner  ///< ???
      );

      /*!
      \brief Copy Constructor

      Makes a deep copy of a Element

      */
      Fluid(const Fluid& old);

      /*!
      \brief Deep copy this instance of fluid and return pointer to the copy

      The Clone() method is used from the virtual base class Element in cases
      where the type of the derived class is unknown and a copy-ctor is needed

      */
      Core::Elements::Element* Clone() const override;

      /*!
      \brief Get shape type of element
      */
      Core::FE::CellType Shape() const override { return distype_; };

      /*!
      \brief set discretization type of element
      */
      virtual void SetDisType(Core::FE::CellType shape)
      {
        distype_ = shape;
        return;
      };

      /*!
      \brief Return number of lines of this element
      */
      int NumLine() const override { return Core::FE::getNumberOfElementLines(distype_); }

      /*!
      \brief Return number of surfaces of this element
      */
      int NumSurface() const override { return Core::FE::getNumberOfElementSurfaces(distype_); }

      /*!
      \brief Return number of volumes of this element (always 1)
      */
      int NumVolume() const override { return Core::FE::getNumberOfElementVolumes(distype_); }

      /*!
      \brief Get vector of Teuchos::RCPs to the lines of this element
      */
      std::vector<Teuchos::RCP<Core::Elements::Element>> Lines() override;

      /*!
      \brief Get vector of Teuchos::RCPs to the surfaces of this element
      */
      std::vector<Teuchos::RCP<Core::Elements::Element>> Surfaces() override;

      /*!
      \brief Get Teuchos::RCP to the internal face adjacent to this element as master element and
      the parent_slave element
      */
      Teuchos::RCP<Core::Elements::Element> CreateFaceElement(
          Core::Elements::Element* parent_slave,  //!< parent slave fluid3 element
          int nnode,                              //!< number of surface nodes
          const int* nodeids,                     //!< node ids of surface element
          Core::Nodes::Node** nodes,              //!< nodes of surface element
          const int lsurface_master,  //!< local surface number w.r.t master parent element
          const int lsurface_slave,   //!< local surface number w.r.t slave parent element
          const std::vector<int>& localtrafomap  //! local trafo map
          ) override;

      /*!
      \brief Return unique ParObject id

      every class implementing ParObject needs a unique id defined at the
      top of this file.
      */
      int UniqueParObjectId() const override { return FluidType::Instance().UniqueParObjectId(); }

      /*!
      \brief Pack this class so it can be communicated

      \ref Pack and \ref Unpack are used to communicate this element

      */
      void Pack(Core::Communication::PackBuffer& data) const override;

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
      \brief Get number of degrees of freedom of a certain node
             (implements pure virtual Core::Elements::Element)

      The element decides how many degrees of freedom its nodes must have.
      As this may vary along a simulation, the element can redecide the
      number of degrees of freedom per node along the way for each of it's nodes
      separately.
      */
      int NumDofPerNode(const Core::Nodes::Node& node) const override
      {
        // number of Dof's is fluid-specific.
        const int nsd = Core::FE::getDimension(distype_);
        if (nsd > 1)
          return nsd + 1;
        else
          FOUR_C_THROW("1D Fluid elements are not supported");

        return 0;
      }

      /*!
      \brief Get number of degrees of freedom per element
             (implements pure virtual Core::Elements::Element)

      The element decides how many element degrees of freedom it has.
      It can redecide along the way of a simulation.

      \note Element degrees of freedom mentioned here are dofs that are visible
            at the level of the total system of equations. Purely internal
            element dofs that are condensed internally should NOT be considered.
      */
      int num_dof_per_element() const override { return 0; }

      /*!
      \brief Print this element
      */
      void Print(std::ostream& os) const override;

      Core::Elements::ElementType& ElementType() const override { return FluidType::Instance(); }

      //@}

      //! @name Input and Creation

      /*!
      \brief Read input for this element
      */
      bool ReadElement(const std::string& eletype, const std::string& distype,
          Input::LineDefinition* linedef) override;

      //@}

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
      int Evaluate(Teuchos::ParameterList& params, Core::FE::Discretization& discretization,
          std::vector<int>& lm, Core::LinAlg::SerialDenseMatrix& elemat1,
          Core::LinAlg::SerialDenseMatrix& elemat2, Core::LinAlg::SerialDenseVector& elevec1,
          Core::LinAlg::SerialDenseVector& elevec2,
          Core::LinAlg::SerialDenseVector& elevec3) override;


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
      int evaluate_neumann(Teuchos::ParameterList& params, Core::FE::Discretization& discretization,
          Core::Conditions::Condition& condition, std::vector<int>& lm,
          Core::LinAlg::SerialDenseVector& elevec1,
          Core::LinAlg::SerialDenseMatrix* elemat1 = nullptr) override;

      //@}

      //! @name Time-dependent subgrid scales
      /*!
      \brief Memory allocation for subgrid-scale arrays
      */
      void ActivateTDS(int nquad, int nsd, double** saccn = nullptr, double** sveln = nullptr,
          double** svelnp = nullptr);

      /*!
      \brief Access to element-specific subgrid-scale arrays
      */
      Teuchos::RCP<FLD::TDSEleData>& TDS() { return tds_; }

      //@}


      //! @name Other

      /*!
      \brief Flag for ALE form of equations
      */
      bool IsAle() const { return is_ale_; }

      virtual void SetIsAle(bool is_ale)
      {
        is_ale_ = is_ale;

        return;
      };


     protected:
      //! discretization type
      Core::FE::CellType distype_;

      //! flag for euler/ale net algorithm
      bool is_ale_;

      //! time-dependent subgrid-scales (only allocated if needed)
      Teuchos::RCP<FLD::TDSEleData> tds_;

      // internal calculation methods

      // don't want = operator
      Fluid& operator=(const Fluid& old);

    };  // class Fluid

    //! Template Meta Programming version of switch over enrichment type
    template <Discret::ELEMENTS::Fluid::EnrichmentType ENRTYPE>
    struct MultipleNumNode
    {
    };
    template <>
    struct MultipleNumNode<Discret::ELEMENTS::Fluid::none>
    {
      static constexpr int multipleNode = 1;
    };
    template <>
    struct MultipleNumNode<Discret::ELEMENTS::Fluid::xwall>
    {
      static constexpr int multipleNode = 2;
    };

    //=======================================================================
    //=======================================================================
    //=======================================================================
    //=======================================================================



    //=======================================================================
    //=======================================================================
    //=======================================================================
    //=======================================================================


    /*!
    \brief An element representing a boundary element of a fluid element

    \note This is a pure Neumann boundary condition element. It's only
          purpose is to evaluate surface Neumann boundary conditions that might be
          adjacent to a parent fluid element. It therefore does not implement
          the Core::Elements::Element::Evaluate method and does not have its own ElementRegister
    class.

    */

    class FluidBoundaryType : public Core::Elements::ElementType
    {
     public:
      std::string Name() const override { return "FluidBoundaryType"; }

      static FluidBoundaryType& Instance();

      Core::Communication::ParObject* Create(const std::vector<char>& data) override;

      Teuchos::RCP<Core::Elements::Element> Create(const int id, const int owner) override;

      void nodal_block_information(
          Core::Elements::Element* dwele, int& numdf, int& dimns, int& nv, int& np) override
      {
      }

      Core::LinAlg::SerialDenseMatrix ComputeNullSpace(
          Core::Nodes::Node& node, const double* x0, const int numdof, const int dimnsp) override
      {
        Core::LinAlg::SerialDenseMatrix nullspace;
        FOUR_C_THROW("method ComputeNullSpace not implemented!");
        return nullspace;
      }

     private:
      static FluidBoundaryType instance_;
    };


    // class FluidBoundary

    class FluidBoundary : public Core::Elements::FaceElement
    {
     public:
      //! @name Friends
      friend class FluidBoundaryType;

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
      explicit FluidBoundary(int id, int owner, int nnode, const int* nodeids,
          Core::Nodes::Node** nodes, Discret::ELEMENTS::Fluid* parent, const int lsurface);

      /*!
      \brief Copy Constructor

      Makes a deep copy of a Element

      */
      FluidBoundary(const FluidBoundary& old);

      /*!
      \brief Deep copy this instance of an element and return pointer to the copy

      The Clone() method is used from the virtual base class Element in cases
      where the type of the derived class is unknown and a copy-ctor is needed

      */
      Core::Elements::Element* Clone() const override;

      /*!
      \brief Get shape type of element
      */
      Core::FE::CellType Shape() const override { return distype_; };

      /*!
      \brief Return number of lines of this element
      */
      int NumLine() const override { return Core::FE::getNumberOfElementLines(Shape()); }

      /*!
      \brief Return number of surfaces of this element
      */
      int NumSurface() const override { return Core::FE::getNumberOfElementSurfaces(Shape()); }

      /*!
      \brief Get vector of Teuchos::RCPs to the lines of this element

      */
      std::vector<Teuchos::RCP<Core::Elements::Element>> Lines() override;

      /*!
      \brief Return unique ParObject id

      every class implementing ParObject needs a unique id defined at the
      top of the parobject.H file.
      */

      std::vector<Teuchos::RCP<Core::Elements::Element>> Surfaces() override;

      /*!
      \brief Return unique ParObject id

      every class implementing ParObject needs a unique id defined at the
      top of the parobject.H file.
      */
      int UniqueParObjectId() const override
      {
        return FluidBoundaryType::Instance().UniqueParObjectId();
      }

      /*!
      \brief Pack this class so it can be communicated

      \ref Pack and \ref Unpack are used to communicate this element

      */
      void Pack(Core::Communication::PackBuffer& data) const override;

      /*!
      \brief Unpack data from a char vector into this class

      \ref Pack and \ref Unpack are used to communicate this element

      */
      void Unpack(const std::vector<char>& data) override;


      //@}

      //! @name Acess methods

      /*!
      \brief Get number of degrees of freedom of a certain node
             (implements pure virtual Core::Elements::Element)

      The element decides how many degrees of freedom its nodes must have.
      As this may vary along a simulation, the element can redecide the
      number of degrees of freedom per node along the way for each of it's nodes
      separately.
      */
      int NumDofPerNode(const Core::Nodes::Node& node) const override { return numdofpernode_; }

      /*!
      \brief Get number of degrees of freedom per element
             (implements pure virtual Core::Elements::Element)

      The element decides how many element degrees of freedom it has.
      It can redecide along the way of a simulation.

      \note Element degrees of freedom mentioned here are dofs that are visible
            at the level of the total system of equations. Purely internal
            element dofs that are condensed internally should NOT be considered.
      */
      int num_dof_per_element() const override { return 0; }

      /*!
      \brief Print this element
      */
      void Print(std::ostream& os) const override;

      Core::Elements::ElementType& ElementType() const override
      {
        return FluidBoundaryType::Instance();
      }

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
      int Evaluate(Teuchos::ParameterList& params, Core::FE::Discretization& discretization,
          std::vector<int>& lm, Core::LinAlg::SerialDenseMatrix& elemat1,
          Core::LinAlg::SerialDenseMatrix& elemat2, Core::LinAlg::SerialDenseVector& elevec1,
          Core::LinAlg::SerialDenseVector& elevec2,
          Core::LinAlg::SerialDenseVector& elevec3) override;

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
      int evaluate_neumann(Teuchos::ParameterList& params, Core::FE::Discretization& discretization,
          Core::Conditions::Condition& condition, std::vector<int>& lm,
          Core::LinAlg::SerialDenseVector& elevec1,
          Core::LinAlg::SerialDenseMatrix* elemat1 = nullptr) override;

      Discret::ELEMENTS::Fluid* parent_element() const
      {
        Core::Elements::Element* parent = this->Core::Elements::FaceElement::parent_element();
        // make sure the static cast below is really valid
        FOUR_C_ASSERT(dynamic_cast<Discret::ELEMENTS::Fluid*>(parent) != nullptr,
            "Master element is no fluid element");
        return static_cast<Discret::ELEMENTS::Fluid*>(parent);
      }

      int SurfaceNumber() { return FaceMasterNumber(); }

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
      void LocationVector(const Core::FE::Discretization& dis, LocationArray& la, bool doDirichlet,
          const std::string& condstring, Teuchos::ParameterList& params) const override;

     protected:
      //! discretization type
      Core::FE::CellType distype_;

      //! numdofpernode
      int numdofpernode_;

     private:
      explicit FluidBoundary(int id, int owner);

      // don't want = operator
      FluidBoundary& operator=(const FluidBoundary& old);

    };  // class FluidBoundary


    //=======================================================================
    //=======================================================================
    //=======================================================================
    //=======================================================================



    //=======================================================================
    //=======================================================================
    //=======================================================================
    //=======================================================================


    /*!
    \brief An element representing an internal face element between two fluid elements

    \note It's only purpose is to evaluate edge based stabilizations for XFEM.
    */
    class FluidIntFaceType : public Core::Elements::ElementType
    {
     public:
      std::string Name() const override { return "FluidIntFaceType"; }

      static FluidIntFaceType& Instance();

      Teuchos::RCP<Core::Elements::Element> Create(const int id, const int owner) override;

      void nodal_block_information(
          Core::Elements::Element* dwele, int& numdf, int& dimns, int& nv, int& np) override
      {
      }

      Core::LinAlg::SerialDenseMatrix ComputeNullSpace(
          Core::Nodes::Node& node, const double* x0, const int numdof, const int dimnsp) override
      {
        Core::LinAlg::SerialDenseMatrix nullspace;
        FOUR_C_THROW("method ComputeNullSpace not implemented!");
        return nullspace;
      }

      /// pre-evaluation
      void pre_evaluate(Core::FE::Discretization& dis, Teuchos::ParameterList& p,
          Teuchos::RCP<Core::LinAlg::SparseOperator> systemmatrix1,
          Teuchos::RCP<Core::LinAlg::SparseOperator> systemmatrix2,
          Teuchos::RCP<Epetra_Vector> systemvector1, Teuchos::RCP<Epetra_Vector> systemvector2,
          Teuchos::RCP<Epetra_Vector> systemvector3) override;

     private:
      static FluidIntFaceType instance_;
    };


    // class FluidIntFace

    class FluidIntFace : public Core::Elements::FaceElement
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
      \param master_parent: The master parent fluid element of this surface
      \param slave_parent: The slave parent fluid element of this surface
      \param lsurface_master: the local surface number of this surface w.r.t. the master parent
      element \param lsurface_slave: the local surface number of this surface w.r.t. the slave
      parent element \param localtrafomap: transformation map between the local coordinate systems
      of the face w.r.t the master parent element's face's coordinate system and the slave element's
      face's coordinate system
      */
      FluidIntFace(int id, int owner, int nnode, const int* nodeids, Core::Nodes::Node** nodes,
          Discret::ELEMENTS::Fluid* parent_master, Discret::ELEMENTS::Fluid* parent_slave,
          const int lsurface_master, const int lsurface_slave,
          const std::vector<int> localtrafomap);

      /*!
      \brief Copy Constructor

      Makes a deep copy of a Element
      */
      FluidIntFace(const FluidIntFace& old);

      /*!
      \brief Deep copy this instance of an element and return pointer to the copy

      The Clone() method is used from the virtual base class Element in cases
      where the type of the derived class is unknown and a copy-ctor is needed

      */
      Core::Elements::Element* Clone() const override;

      /*!
      \brief Get shape type of element
      */
      Core::FE::CellType Shape() const override;

      /*!
      \brief Return number of lines of this element
      */
      int NumLine() const override { return Core::FE::getNumberOfElementLines(Shape()); }

      /*!
      \brief Return number of surfaces of this element
      */
      int NumSurface() const override { return Core::FE::getNumberOfElementSurfaces(Shape()); }

      /*!
      \brief Get vector of Teuchos::RCPs to the lines of this element
      */
      std::vector<Teuchos::RCP<Core::Elements::Element>> Lines() override;

      /*!
      \brief Return unique ParObject id

      every class implementing ParObject needs a unique id defined at the
      top of the parobject.H file.
      */
      std::vector<Teuchos::RCP<Core::Elements::Element>> Surfaces() override;

      /*!
      \brief Return unique ParObject id

      every class implementing ParObject needs a unique id defined at the
      top of the parobject.H file.
      */
      int UniqueParObjectId() const override
      {
        return FluidIntFaceType::Instance().UniqueParObjectId();
      }

      /*!
      \brief Pack this class so it can be communicated

      \ref Pack and \ref Unpack are used to communicate this element
      */
      void Pack(Core::Communication::PackBuffer& data) const override;

      /*!
      \brief Unpack data from a char vector into this class

      \ref Pack and \ref Unpack are used to communicate this element
      */
      void Unpack(const std::vector<char>& data) override;



      //@}

      //! @name Acess methods


      /*!
      \brief Get number of degrees of freedom of a certain node
             (implements pure virtual Core::Elements::Element)

      The element decides how many degrees of freedom its nodes must have.
      As this may vary along a simulation, the element can redecide the
      number of degrees of freedom per node along the way for each of it's nodes
      separately.
      */
      int NumDofPerNode(const Core::Nodes::Node& node) const override
      {
        return std::max(
            ParentMasterElement()->NumDofPerNode(node), ParentSlaveElement()->NumDofPerNode(node));
      }

      /*!
      \brief Get number of degrees of freedom per element
             (implements pure virtual Core::Elements::Element)

      The element decides how many element degrees of freedom it has.
      It can redecide along the way of a simulation.

      \note Element degrees of freedom mentioned here are dofs that are visible
            at the level of the total system of equations. Purely internal
            element dofs that are condensed internally should NOT be considered.
      */
      int num_dof_per_element() const override { return 0; }

      /*!
      \brief create the location vector for patch of master and slave element

      \note All dofs shared by master and slave element are contained only once. Dofs from interface
      nodes are also included.
      */
      void PatchLocationVector(Core::FE::Discretization& discretization,  ///< discretization
          std::vector<int>& nds_master,        ///< nodal dofset w.r.t master parent element
          std::vector<int>& nds_slave,         ///< nodal dofset w.r.t slave parent element
          std::vector<int>& patchlm,           ///< local map for gdof ids for patch of elements
          std::vector<int>& lm_masterToPatch,  ///< local map between lm_master and lm_patch
          std::vector<int>& lm_slaveToPatch,   ///< local map between lm_slave and lm_patch
          std::vector<int>& lm_faceToPatch,    ///< local map between lm_face and lm_patch
          std::vector<int>&
              lm_masterNodeToPatch,  ///< local map between master nodes and nodes in patch
          std::vector<int>&
              lm_slaveNodeToPatch,  ///< local map between slave nodes and nodes in patch
          Teuchos::RCP<std::map<int, int>>
              pbcconnectivity  ///< connectivity between slave and PBC's master nodes
      );

      /*!
      \brief create the location vector for patch of master and slave element

      \note All dofs shared by master and slave element are contained only once. Dofs from interface
      nodes are also included.
      */
      void PatchLocationVector(Core::FE::Discretization& discretization,  ///< discretization
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
              lm_slaveNodeToPatch,  ///< local map between slave nodes and nodes in patch
          Teuchos::RCP<std::map<int, int>>
              pbcconnectivity  ///< connectivity between slave and PBC's master nodes
      );

      /*!
      \brief Print this element
      */
      void Print(std::ostream& os) const override;

      Core::Elements::ElementType& ElementType() const override
      {
        return FluidIntFaceType::Instance();
      }

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
      int Evaluate(Teuchos::ParameterList& params, Core::FE::Discretization& discretization,
          std::vector<int>& lm, Core::LinAlg::SerialDenseMatrix& elemat1,
          Core::LinAlg::SerialDenseMatrix& elemat2, Core::LinAlg::SerialDenseVector& elevec1,
          Core::LinAlg::SerialDenseVector& elevec2,
          Core::LinAlg::SerialDenseVector& elevec3) override;

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
      int evaluate_neumann(Teuchos::ParameterList& params, Core::FE::Discretization& discretization,
          Core::Conditions::Condition& condition, std::vector<int>& lm,
          Core::LinAlg::SerialDenseVector& elevec1,
          Core::LinAlg::SerialDenseMatrix* elemat1 = nullptr) override;

      /*!
      \brief return the master parent fluid element
      */
      Discret::ELEMENTS::Fluid* ParentMasterElement() const
      {
        Core::Elements::Element* parent = this->Core::Elements::FaceElement::ParentMasterElement();
        // make sure the static cast below is really valid
        FOUR_C_ASSERT(dynamic_cast<Discret::ELEMENTS::Fluid*>(parent) != nullptr,
            "Master element is no fluid element");
        return static_cast<Discret::ELEMENTS::Fluid*>(parent);
      }

      /*!
      \brief return the slave parent fluid element
      */
      Discret::ELEMENTS::Fluid* ParentSlaveElement() const
      {
        Core::Elements::Element* parent = this->Core::Elements::FaceElement::ParentSlaveElement();
        // make sure the static cast below is really valid
        FOUR_C_ASSERT(dynamic_cast<Discret::ELEMENTS::Fluid*>(parent) != nullptr,
            "Slave element is no fluid element");
        return static_cast<Discret::ELEMENTS::Fluid*>(parent);
      }

      //@}

     private:
      // don't want = operator
      FluidIntFace& operator=(const FluidIntFace& old);

    };  // class FluidIntFace



  }  // namespace ELEMENTS
}  // namespace Discret



FOUR_C_NAMESPACE_CLOSE

#endif
