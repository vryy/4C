/*----------------------------------------------------------------------------*/
/*! \file

\brief ALE element for 3D case


\level 1
*/
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
#ifndef ALE3_H
#define ALE3_H

/*----------------------------------------------------------------------------*/
/* header inclusions */
#include <Teuchos_RCP.hpp>
#include <Epetra_Vector.h>

#include "../drt_lib/drt_elementtype.H"
#include "../drt_fem_general/drt_utils_integration.H"
#include "../drt_fem_general/drt_utils_local_connectivity_matrices.H"
#include "../drt_lib/drt_element.H"
#include "../drt_fem_general/drt_utils_gausspoints.H"
#include "../drt_fem_general/drt_utils_integration.H"

const int NUMDIM_ALE3 = 3;  ///< number of dimensions
const int NODDOF_ALE3 = 3;  ///< number of dofs per node

/*----------------------------------------------------------------------------*/
/* forward declarations */
namespace DRT
{
  class Discretization;

  namespace ELEMENTS
  {
    // where should that be put?
    template <DRT::Element::DiscretizationType dtype>
    struct DisTypeToNumCornerNodes
    {
    };
    template <>
    struct DisTypeToNumCornerNodes<DRT::Element::tet4>
    {
      static const int numCornerNodes = 4;
    };
    template <>
    struct DisTypeToNumCornerNodes<DRT::Element::tet10>
    {
      static const int numCornerNodes = 4;
    };
    template <>
    struct DisTypeToNumCornerNodes<DRT::Element::hex8>
    {
      static const int numCornerNodes = 8;
    };
    template <>
    struct DisTypeToNumCornerNodes<DRT::Element::hex20>
    {
      static const int numCornerNodes = 8;
    };
    template <>
    struct DisTypeToNumCornerNodes<DRT::Element::hex27>
    {
      static const int numCornerNodes = 8;
    };
    template <>
    struct DisTypeToNumCornerNodes<DRT::Element::pyramid5>
    {
      static const int numCornerNodes = 5;
    };
    template <>
    struct DisTypeToNumCornerNodes<DRT::Element::wedge6>
    {
      static const int numCornerNodes = 6;
    };
    template <>
    struct DisTypeToNumCornerNodes<DRT::Element::wedge15>
    {
      static const int numCornerNodes = 6;
    };
    template <>
    struct DisTypeToNumCornerNodes<DRT::Element::nurbs8>
    {
      static const int numCornerNodes = 8;
    };
    template <>
    struct DisTypeToNumCornerNodes<DRT::Element::nurbs27>
    {
      static const int numCornerNodes = 8;
    };


    /*----------------------------------------------------------------------------*/
    /*----------------------------------------------------------------------------*/
    class Ale3Surface;
    class Ale3_Impl_Interface;
    template <DRT::Element::DiscretizationType distype>
    class Ale3_Impl;
    class Ale3Surface_Impl_Interface;
    template <DRT::Element::DiscretizationType distype>
    class Ale3Surface_Impl;


    class Ale3Type : public DRT::ElementType
    {
     public:
      std::string Name() const { return "Ale3Type"; }

      static Ale3Type& Instance();

      virtual DRT::ParObject* Create(const std::vector<char>& data);

      virtual Teuchos::RCP<DRT::Element> Create(
          const std::string eletype, const std::string eledistype, const int id, const int owner);

      virtual Teuchos::RCP<DRT::Element> Create(const int id, const int owner);

      virtual void NodalBlockInformation(
          DRT::Element* dwele, int& numdf, int& dimns, int& nv, int& np);

      virtual void ComputeNullSpace(DRT::Discretization& dis, std::vector<double>& ns,
          const double* x0, int numdf, int dimns);

      virtual void SetupElementDefinition(
          std::map<std::string, std::map<std::string, DRT::INPUT::LineDefinition>>& definitions);

     private:
      static Ale3Type instance_;
    };



    /*!
    \brief A C++ wrapper for the ale3 element
    */
    class Ale3 : public DRT::Element
    {
     public:
      //! @name Friends
      friend class Ale3Surface;
      // friend class Ale3_Impl_Interface;
      // friend class Ale3_Impl<

      //@}
      //! @name Constructors and destructors and related methods

      /*!
      \brief Standard Constructor

      \param id : A unique global id
      */
      Ale3(int id, int owner);

      /*!
      \brief Copy Constructor

      Makes a deep copy of a Element

      */
      Ale3(const Ale3& old);

      /*!
      \brief Deep copy this instance of Ale3 and return pointer to the copy

      The Clone() method is used from the virtual base class Element in cases
      where the type of the derived class is unknown and a copy-ctor is needed

      */
      DRT::Element* Clone() const;

      /*!
      \brief Get shape type of element
      */
      virtual DiscretizationType Shape() const;

      /*!
      \brief Return number of lines of this element
      */
      virtual int NumLine() const { return 0; }

      /*!
      \brief Return number of surfaces of this element
      */
      virtual int NumSurface() const
      {
        switch (NumNode())
        {
          case 8:
          case 20:
          case 27:
            return 6;  // hex
          case 4:
          case 10:
            return 4;  // tet
          case 6:
          case 15:
          case 5:
            return 5;  // wedge or pyramid
          default:
            dserror("Could not determine number of surfaces");
            return -1;
        }
      }

      /*!
      \brief Return number of volumes of this element (always 1)
      */
      virtual inline int NumVolume() const { return 1; }

      /*!
      \brief Get vector of Teuchos::RCPs to the surfaces of this element

      */
      virtual std::vector<Teuchos::RCP<DRT::Element>> Surfaces();

      /*!
      \brief Get vector of Teuchos::RCPs to the volume of this element

      */
      virtual std::vector<Teuchos::RCP<DRT::Element>> Volumes();

      /*!
      \brief Return unique ParObject id

      every class implementing ParObject needs a unique id defined at the
      top of this file.
      */
      virtual int UniqueParObjectId() const { return Ale3Type::Instance().UniqueParObjectId(); }

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
      \brief Destructor

      */
      virtual ~Ale3();

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
      virtual int NumDofPerNode(const DRT::Node& node) const { return 3; }

      /*!
      \brief Get number of degrees of freedom per element
             (implements pure virtual DRT::Element)

      The element decides how many element degrees of freedom it has.
      It can redecide along the way of a simulation.

      \note Element degrees of freedom mentioned here are dofs that are visible
            at the level of the total system of equations. Purely internal
            element dofs that are condensed internally should NOT be considered.
      */
      virtual int NumDofPerElement() const { return 0; }

      /*!
      \brief Print this element
      */
      void Print(std::ostream& os) const;

      virtual DRT::ElementType& ElementType() const { return Ale3Type::Instance(); }

      //@}

      //! @name Input and Creation

      /*!
      \brief Read input for this element
      */
      bool ReadElement(const std::string& eletype, const std::string& distype,
          DRT::INPUT::LineDefinition* linedef);

      //@}

      //! @name Evaluation

      /*!
      \brief Evaluate an element

      Evaluate ale3 element stiffness, mass, internal forces etc

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


      /*!
      \brief Evaluate a Neumann boundary condition

      this method evaluates a surfaces Neumann condition on the shell element

      \param params (in/out)    : ParameterList for communication between control routine
                                  and elements
      \param discretization (in): A reference to the underlying discretization
      \param condition (in)     : The condition to be evaluated
      \param lm (in)            : location vector of this element
      \param elevec1 (out)      : vector to be filled by element. If NULL on input,

      \return 0 if successful, negative otherwise
      */
      virtual int EvaluateNeumann(Teuchos::ParameterList& params,
          DRT::Discretization& discretization, DRT::Condition& condition, std::vector<int>& lm,
          Epetra_SerialDenseVector& elevec1, Epetra_SerialDenseMatrix* elemat1 = NULL);


      //@}

      //! @name Other

      //@}

      //! action parameters recognized by ale3
      enum ActionType
      {
        none,
        calc_ale_solid,         ///< compute stiffness based on fully nonlinear elastic solid with
                                ///< hyperelastic material law
        calc_ale_solid_linear,  ///< compute stiffness based on linear elastic solid with
                                ///< hyperelastic material law
        calc_ale_springs_material,  ///< compute stiffness based on springs algorithm in material
                                    ///< configuration
        calc_ale_springs_spatial,   ///< compute stiffness based on springs algorithm in spatial
                                    ///< configuration
        calc_ale_laplace_material,  ///< compute stiffness based on laplacian smoothing based on
                                    ///< material configuration
        calc_ale_laplace_spatial,   ///< compute stiffness based on laplacian smoothing based on
                                    ///< spatial configuration
        calc_ale_node_normal,       ///< Calculate boundary node normal
        ba_calc_ale_node_normal,    ///< Calculate boundary node normal TODO (@ HAHN unused action
                                    ///< please fix me!)
        setup_material,             ///< Setup material in case of ElastHyper Tool Box
        calc_det_jac  ///< calculate Jacobian determinant and mesh quality measure according to
                      ///< [Oddy et al. 1988]
      };

     private:
      //! don't want = operator
      Ale3& operator=(const Ale3& old);
    };

    //=======================================================================
    //=======================================================================
    //=======================================================================
    //=======================================================================

    class Ale3_Impl_Interface : SingletonDestruction
    {
     public:
      /// Empty constructor
      Ale3_Impl_Interface() {}

      virtual ~Ale3_Impl_Interface() {}

      /// Internal implementation class for fluid element
      static Ale3_Impl_Interface* Impl(DRT::ELEMENTS::Ale3* ele);

      virtual void static_ke_spring(Ale3* ele,       ///< pointer to element
          Epetra_SerialDenseMatrix& sys_mat,         ///< element stiffness matrix (to be filled)
          Epetra_SerialDenseVector& residual,        ///< element residual vector (to be filled)
          const std::vector<double>& displacements,  ///< nodal displacements
          const bool spatialconfiguration            ///< use spatial configuration (true), material
                                                     ///< configuration (false)
          ) = 0;

      virtual void static_ke_nonlinear(Ale3* ele,  ///< pointer to element
          DRT::Discretization& discretization,     ///< discretization
          std::vector<int>& lm,                    ///< node owner procs
          Epetra_SerialDenseMatrix& sys_mat,       ///< element stiffness matrix (to be filled)
          Epetra_SerialDenseVector& residual,      ///< element residual vector (to be filled)
          std::vector<double>& my_dispnp,          ///< nodal displacements
          Teuchos::ParameterList& params,          ///< parameter list
          const bool spatialconfiguration          ///< use spatial configuration (true), material
                                                   ///< configuration (false)
          ) = 0;

      virtual void static_ke_laplace(Ale3* ele,  ///< pointer to element
          DRT::Discretization& dis,              ///< discretization
          Epetra_SerialDenseMatrix& sys_mat,     ///< element stiffnes matrix (to be filled)
          Epetra_SerialDenseVector& residual,    ///< element residual vector (to be filled)
          std::vector<double>& my_dispnp,        ///< nodal displacements
          Teuchos::RCP<MAT::Material> material,  ///< material law
          const bool spatialconfiguration        ///< use spatial configuration (true), material
                                                 ///< configuration (false)
          ) = 0;

      virtual void ElementNodeNormal(
          Ale3* ele, Epetra_SerialDenseVector& elevec1, std::vector<double>& my_dispnp) = 0;
    };

    template <DRT::Element::DiscretizationType distype>
    class Ale3_Impl : public Ale3_Impl_Interface
    {
      Ale3_Impl() {}

     public:
      ~Ale3_Impl() {}

      /// Singleton access method
      static Ale3_Impl<distype>* Instance(bool create = true);

      /// called upon destruction
      virtual void Done();

      void static_ke_laplace(Ale3* ele,          ///< pointer to element
          DRT::Discretization& dis,              ///< discretization
          Epetra_SerialDenseMatrix& sys_mat,     ///< element stiffnes matrix (to be filled)
          Epetra_SerialDenseVector& residual,    ///< element residual vector (to be filled)
          std::vector<double>& my_dispnp,        ///< nodal displacements
          Teuchos::RCP<MAT::Material> material,  ///< material law
          const bool spatialconfiguration        ///< use spatial configuration (true), material
                                                 ///< configuration (false)
      );

      void static_ke_spring(Ale3* ele,                ///< pointer to element
          Epetra_SerialDenseMatrix& sys_mat_epetra,   ///< element stiffness matrix (to be filled)
          Epetra_SerialDenseVector& residual_epetra,  ///< element residual vector (to be filled)
          const std::vector<double>& displacements,   ///< nodal displacements
          const bool spatialconfiguration  ///< use spatial configuration (true), material
                                           ///< configuration (false)
      );

      void static_ke_nonlinear(Ale3* ele,       ///< pointer to element
          DRT::Discretization& discretization,  ///< discretization
          std::vector<int>& lm,                 ///< node owner procs
          Epetra_SerialDenseMatrix& sys_mat,    ///< element stiffness matrix (to be filled)
          Epetra_SerialDenseVector& residual,   ///< element residual vector (to be filled)
          std::vector<double>& my_dispnp,       ///< nodal displacements
          Teuchos::ParameterList& params,       ///< parameter list
          const bool spatialconfiguration       ///< use spatial configuration (true), material
                                                ///< configuration (false)
      );

      void ElementNodeNormal(Ale3* ele,       ///< pointer to element
          Epetra_SerialDenseVector& elevec1,  ///< normal vector (to be filled)
          std::vector<double>& my_dispnp      ///< nodal displacements
      );

      //! Calculate Jacobian matrix and its determinant
      void CalcJacobian(Ale3* ele,  ///< pointer to element
          double& detJ              ///< determinant of Jacobian matrix
      );

     private:
      static const int iel = DRT::UTILS::DisTypeToNumNodePerEle<distype>::numNodePerElement;
      static const int numcnd = DisTypeToNumCornerNodes<distype>::numCornerNodes;

      inline void ale3_edge_geometry(int i, int j, const LINALG::Matrix<3, iel>& xyze,
          double& length, double& dx, double& dy, double& dz);

      /*!
      \brief Prevents node s from penetrating face pqr.

      The calculated dynamic triangle sjq is perpendicular to edge pr and face
      pqr. According to Farhat et al.

      The interface of this function is rather wierd, I changed it to reduce
      re-calculation of vectors that are used multiple times.

      \param node_i (in)    : Determine position of triangle sjq in
      tetrahedra.
      \param sq             : edge vector
      \param len_sq         : length
      \param rp             : edge vector
      \param len_rp         : length
      \param qp             : edge vector
      \param local_x        : vector needed in calculation
      \param sysmat (in/out): The element's sys_mat
      */
      void ale3_add_tria_stiffness(int node_p, int node_q, int node_r, int node_s,
          const LINALG::Matrix<3, 1>& sq, const double len_sq, const LINALG::Matrix<3, 1>& rp,
          const double len_rp, const LINALG::Matrix<3, 1>& qp, const LINALG::Matrix<3, 1>& local_x,
          LINALG::Matrix<3 * iel, 3 * iel>& sys_mat);

      /*!
      \brief Prevents node-face-penetration for given nodes.

      Twelve-triangle configuration. According to Farhat et al.

      \param tet_i (in)    : Nodes. Arbitrary succession.
      \param sysmat (in/out): The element's sys_mat
      \param xyze (in)     : The actual element coordinates
      */
      void ale3_add_tetra_stiffness(int tet_0, int tet_1, int tet_2, int tet_3,
          LINALG::Matrix<3 * iel, 3 * iel>& sys_mat, const LINALG::Matrix<3, iel>& xyze);
      inline void ale3_tors_spring_tet4(
          LINALG::Matrix<3 * iel, 3 * iel>& sys_mat, const LINALG::Matrix<3, iel>& xyze);
      inline void ale3_tors_spring_pyramid5(
          LINALG::Matrix<3 * iel, 3 * iel>& sys_mat, const LINALG::Matrix<3, iel>& xyze);
      inline void ale3_tors_spring_wedge6(
          LINALG::Matrix<3 * iel, 3 * iel>& sys_mat, const LINALG::Matrix<3, iel>& xyze);
      inline void ale3_tors_spring_hex8(
          LINALG::Matrix<3 * iel, 3 * iel>& sys_mat, const LINALG::Matrix<3, iel>& xyze);
      inline void ale3_tors_spring_nurbs27(
          LINALG::Matrix<3 * iel, 3 * iel>& sys_mat, const LINALG::Matrix<3, iel>& xyze);


      inline DRT::UTILS::GaussRule3D getOptimalGaussrule();

      Ale3_Impl<distype> operator=(const Ale3_Impl<distype> other);
    };


    //=======================================================================
    //=======================================================================
    //=======================================================================
    //=======================================================================


    class Ale3SurfaceType : public DRT::ElementType
    {
     public:
      std::string Name() const { return "Ale3SurfaceType"; }

      static Ale3SurfaceType& Instance();

      virtual Teuchos::RCP<DRT::Element> Create(const int id, const int owner);

      virtual void NodalBlockInformation(
          DRT::Element* dwele, int& numdf, int& dimns, int& nv, int& np)
      {
      }

      virtual void ComputeNullSpace(
          DRT::Discretization& dis, std::vector<double>& ns, const double* x0, int numdf, int dimns)
      {
      }

     private:
      static Ale3SurfaceType instance_;
    };


    /*!
    \brief An element representing a surface of a ale3 element

    \note This is a pure Neumann boundary condition element. It's only
          purpose is to evaluate surface Neumann boundary conditions that might be
          adjacent to a parent ale3 element. It therefore does not implement
          the DRT::Element::Evaluate method and does not have its own ElementRegister class.

    */
    class Ale3Surface : public DRT::FaceElement
    {
     public:
      //! @name Constructors and destructors and related methods

      /*!
      \brief Standard Constructor

      \param id : A unique global id
      \param owner: Processor owning this surface
      \param nnode: Number of nodes attached to this element
      \param nodeids: global ids of nodes attached to this element
      \param nodes: the discretizations map of nodes to build ptrs to nodes from
      \param parent: The parent ale element of this surface
      \param lsurface: the local surface number of this surface w.r.t. the parent element
      */
      Ale3Surface(int id, int owner, int nnode, const int* nodeids, DRT::Node** nodes,
          DRT::ELEMENTS::Ale3* parent, const int lsurface);

      /*!
      \brief Copy Constructor

      Makes a deep copy of a Element

      */
      Ale3Surface(const Ale3Surface& old);

      /*!
      \brief Deep copy this instance of an element and return pointer to the copy

      The Clone() method is used from the virtual base class Element in cases
      where the type of the derived class is unknown and a copy-ctor is needed

      */
      DRT::Element* Clone() const;

      /*!
      \brief Get shape type of element
      */
      virtual DiscretizationType Shape() const;

      /*!
      \brief Return unique ParObject id

      every class implementing ParObject needs a unique id defined at the
      top of the drt_parobject.H file.
      */
      virtual int UniqueParObjectId() const
      {
        return Ale3SurfaceType::Instance().UniqueParObjectId();
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
      \brief Destructor

      */
      virtual ~Ale3Surface();

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
      virtual int NumDofPerNode(const DRT::Node& node) const { return 3; }

      /*!
      \brief Get number of degrees of freedom per element
             (implements pure virtual DRT::Element)

      The element decides how many element degrees of freedom it has.
      It can redecide along the way of a simulation.

      \note Element degrees of freedom mentioned here are dofs that are visible
            at the level of the total system of equations. Purely internal
            element dofs that are condensed internally should NOT be considered.
      */
      virtual int NumDofPerElement() const { return 0; }

      /*!
      \brief Print this element
      */
      void Print(std::ostream& os) const;

      virtual DRT::ElementType& ElementType() const { return Ale3SurfaceType::Instance(); }

      //@}

      //! @name Evaluate methods

      /*!
      \brief Evaluate the ale3 surface element

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
      virtual int Evaluate(Teuchos::ParameterList& params, DRT::Discretization& discretization,
          std::vector<int>& lm, Epetra_SerialDenseMatrix& elemat1,
          Epetra_SerialDenseMatrix& elemat2, Epetra_SerialDenseVector& elevec1,
          Epetra_SerialDenseVector& elevec2, Epetra_SerialDenseVector& elevec3);



      /*!
      \brief Evaluate a Neumann boundary condition

      this method evaluates a surface Neumann condition on the ale3 element

      \param params (in/out)    : ParameterList for communication between control routine
                                  and elements
      \param discretization (in): A reference to the underlying discretization
      \param condition (in)     : The condition to be evaluated
      \param lm (in)            : location vector of this element
      \param elevec1 (out)      : vector to be filled by element. If NULL on input,

      \return 0 if successful, negative otherwise
      */
      virtual int EvaluateNeumann(Teuchos::ParameterList& params,
          DRT::Discretization& discretization, DRT::Condition& condition, std::vector<int>& lm,
          Epetra_SerialDenseVector& elevec1, Epetra_SerialDenseMatrix* elemat1 = NULL);

      //@}

     private:
      // don't want = operator
      Ale3Surface& operator=(const Ale3Surface& old);

      //  compute kovariant metric tensor G for ale surface element
      //                                                  gammi 04/07
      void f3_metric_tensor_for_surface(const Epetra_SerialDenseMatrix xyze,
          const Epetra_SerialDenseMatrix deriv, Epetra_SerialDenseMatrix& metrictensor,
          double* drs);

    };  // class Ale3Surface

    class Ale3Surface_Impl_Interface : SingletonDestruction
    {
     public:
      /// Empty constructor
      Ale3Surface_Impl_Interface() {}

      virtual ~Ale3Surface_Impl_Interface() {}

      /// Internal implementation class for ale surface element
      static Ale3Surface_Impl_Interface* Impl(DRT::ELEMENTS::Ale3Surface* ele);

      virtual void ElementNodeNormal(Ale3Surface* ele, Teuchos::ParameterList& params,
          DRT::Discretization& discretization, std::vector<int>& lm,
          Epetra_SerialDenseVector& elevec1, std::vector<double>& mydispnp) = 0;
    };

    template <DRT::Element::DiscretizationType distype>
    class Ale3Surface_Impl : public Ale3Surface_Impl_Interface
    {
      Ale3Surface_Impl() {}

     public:
      ~Ale3Surface_Impl() {}

      /// Singleton access method
      static Ale3Surface_Impl<distype>* Instance(bool create = true);

      /// called upon destruction
      virtual void Done();

      //! number of element nodes
      static const int bdrynen_ = DRT::UTILS::DisTypeToNumNodePerEle<distype>::numNodePerElement;

      //! number of spatial dimensions of boundary element
      static const int bdrynsd_ = DRT::UTILS::DisTypeToDim<distype>::dim;

      //! number of spatial dimensions of parent element
      static const int nsd_ = bdrynsd_ + 1;

      //! number of degrees of freedom per node
      static const int numdofpernode_ = nsd_;

      void ElementNodeNormal(Ale3Surface* ele, Teuchos::ParameterList& params,
          DRT::Discretization& discretization, std::vector<int>& lm,
          Epetra_SerialDenseVector& elevec1, std::vector<double>& mydispnp);

     private:
      //  static const int iel = DRT::UTILS::DisTypeToNumNodePerEle<distype>::numNodePerElement;
      //  static const int numcnd = DisTypeToNumCornerNodes<distype>::numCornerNodes;

      //  inline void ale3_edge_geometry(int i, int j, const LINALG::Matrix<3, iel>& xyze,
      //                                 double& length,
      //                                 double& dx,
      //                                 double& dy,
      //                                 double& dz);

      //  /*!
      //  \brief Prevents node s from penetrating face pqr.
      //
      //  The calculated dynamic triangle sjq is perpendicular to edge pr and face
      //  pqr. According to Farhat et al.
      //
      //  The interface of this function is rather wierd, I changed it to reduce
      //  re-calculation of vectors that are used multiple times.
      //
      //  \param node_i (in)    : Determine position of triangle sjq in
      //  tetrahedra.
      //  \param sq             : edge vector
      //  \param len_sq         : length
      //  \param rp             : edge vector
      //  \param len_rp         : length
      //  \param qp             : edge vector
      //  \param local_x        : vector needed in calculation
      //  \param sysmat (in/out): The element's sys_mat
      //  */
      //  void ale3_add_tria_stiffness(
      //      int node_p, int node_q, int node_r, int node_s,
      //      const LINALG::Matrix<3, 1>& sq,
      //      const double len_sq,
      //      const LINALG::Matrix<3, 1>& rp,
      //      const double len_rp,
      //      const LINALG::Matrix<3, 1>& qp,
      //      const LINALG::Matrix<3, 1>& local_x,
      //      LINALG::Matrix<3*iel,3*iel>& sys_mat);
      //
      //  /*!
      //  \brief Prevents node-face-penetration for given nodes.
      //
      //  Twelve-triangle configuration. According to Farhat et al.
      //
      //  \param tet_i (in)    : Nodes. Arbitrary succession.
      //  \param sysmat (in/out): The element's sys_mat
      //  \param xyze (in)     : The actual element coordinates
      //  */
      //  void ale3_add_tetra_stiffness(int tet_0, int tet_1, int tet_2, int tet_3,
      //                                LINALG::Matrix<3*iel,3*iel>& sys_mat,
      //                                const LINALG::Matrix<3,iel>& xyze);
      //  inline void ale3_tors_spring_tet4(LINALG::Matrix<3*iel,3*iel>& sys_mat,
      //                                    const LINALG::Matrix<3,iel>& xyze);
      //  inline void ale3_tors_spring_pyramid5(LINALG::Matrix<3*iel,3*iel>& sys_mat,
      //                                        const LINALG::Matrix<3,iel>& xyze);
      //  inline void ale3_tors_spring_wedge6(LINALG::Matrix<3*iel,3*iel>& sys_mat,
      //                                      const LINALG::Matrix<3,iel>& xyze);
      //  inline void ale3_tors_spring_hex8(LINALG::Matrix<3*iel,3*iel>& sys_mat,
      //                                    const LINALG::Matrix<3,iel>& xyze);
      //  inline void ale3_tors_spring_nurbs27(LINALG::Matrix<3*iel,3*iel>& sys_mat,
      //               const LINALG::Matrix<3,iel>& xyze);
      //
      //
      //  inline DRT::UTILS::GaussRule3D getOptimalGaussrule();

      // Ale3Surface_Impl<distype> operator=(const Ale3Surface_Impl<distype> other);


     protected:
      //! array for shape functions for boundary element
      LINALG::Matrix<bdrynen_, 1> funct_;
      //! array for shape function derivatives for boundary element
      LINALG::Matrix<bdrynsd_, bdrynen_> deriv_;
      //! integration factor
      double fac_;
      //! normal vector pointing out of the domain
      LINALG::Matrix<nsd_, 1> unitnormal_;
      //! infinitesimal area element drs
      double drs_;
      //! coordinates of current integration point in reference coordinates
      LINALG::Matrix<bdrynsd_, 1> xsi_;
      //! node coordinates for boundary element
      LINALG::Matrix<nsd_, bdrynen_> xyze_;
    };



  }  // namespace ELEMENTS
}  // namespace DRT

#endif
