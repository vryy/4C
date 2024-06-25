/*----------------------------------------------------------------------*/
/*! \file

\brief dummy 3D boundary element without any physics


\level 2
*/
/*----------------------------------------------------------------------*/

#ifndef FOUR_C_BELE_BELE3_HPP
#define FOUR_C_BELE_BELE3_HPP


#include "4C_config.hpp"

#include "4C_comm_parobjectfactory.hpp"
#include "4C_fem_general_element.hpp"
#include "4C_fem_general_elementtype.hpp"
#include "4C_fem_general_node.hpp"
#include "4C_fem_general_utils_integration.hpp"
#include "4C_linalg_serialdensematrix.hpp"

#include <Epetra_Vector.h>
#include <Teuchos_RCP.hpp>

FOUR_C_NAMESPACE_OPEN

namespace Core::FE
{
  class Discretization;
}  // namespace Core::FE

namespace Discret
{
  namespace ELEMENTS
  {
    class Bele3Type : public Core::Elements::ElementType
    {
     public:
      std::string Name() const override { return "Bele3Type"; }

      static Bele3Type& Instance();

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

     private:
      static Bele3Type instance_;
    };

    /*!
     * A 3D boundary element with no physics attached
     *
     * This element is meant to have no physics. It can be used to have a boundary discretization
     * of surface/boundary elements. They can be of any 2d shape (quad4,quad9,tri3,...)
     *
     * The number of dof per node is set to 3 per default, so we can define displacement vectors by
     * using fill_complete on the boundary discretization. Furthermore numdofpernode can be adapted
     * if necessary.
     *
     */
    class Bele3 : public Core::Elements::Element
    {
      // friend class to fill number of dofs per node exclusively during creation
      friend class Bele3Type;

     public:
      //@}
      //! @name Constructors and destructors and related methods

      /*!
      \brief Standard Constructor
      */
      explicit Bele3(int id,  ///< A unique global id
          int owner           ///< proc num that owns this element
      );

      /*!
      \brief Copy Constructor

      Makes a deep copy of a Element

      */
      explicit Bele3(const Bele3& old);

      Core::Elements::Element* Clone() const override;
      Core::FE::CellType Shape() const override;
      int NumLine() const override
      {
        if (num_node() == 9 || num_node() == 8 || num_node() == 4)
          return 4;
        else if (num_node() == 3 || num_node() == 6)
          return 3;
        else
        {
          FOUR_C_THROW("Could not determine number of lines");
          return -1;
        }
      }
      int NumSurface() const override { return 1; }
      int NumVolume() const override { return -1; }
      std::vector<Teuchos::RCP<Core::Elements::Element>> Lines() override;
      std::vector<Teuchos::RCP<Core::Elements::Element>> Surfaces() override;
      int UniqueParObjectId() const override { return Bele3Type::Instance().UniqueParObjectId(); }
      void pack(Core::Communication::PackBuffer& data) const override;
      void unpack(const std::vector<char>& data) override;


      //@}

      //! @name Access methods

      int NumDofPerNode(const Core::Nodes::Node&) const override { return numdofpernode_; }
      int num_dof_per_element() const override { return 0; }
      void print(std::ostream& os) const override;
      Core::Elements::ElementType& ElementType() const override { return Bele3Type::Instance(); }

      //@}

      //! @name Evaluation

      int evaluate(Teuchos::ParameterList& params, Core::FE::Discretization& discretization,
          std::vector<int>& lm, Core::LinAlg::SerialDenseMatrix& elemat1,
          Core::LinAlg::SerialDenseMatrix& elemat2, Core::LinAlg::SerialDenseVector& elevec1,
          Core::LinAlg::SerialDenseVector& elevec2,
          Core::LinAlg::SerialDenseVector& elevec3) override;

      int evaluate_neumann(Teuchos::ParameterList& params, Core::FE::Discretization& discretization,
          Core::Conditions::Condition& condition, std::vector<int>& lm,
          Core::LinAlg::SerialDenseVector& elevec1,
          Core::LinAlg::SerialDenseMatrix* elemat1 = nullptr) override;

      /// Read input for this element
      bool ReadElement(const std::string& eletype, const std::string& distype,
          Input::LineDefinition* linedef) override;
      //@}

      //! @name Other
      //! does this element have non-zero displacements or not
      //  bool IsMoving() const { return is_moving_; }

      //@}


     private:
      /*!
        \brief Set number of dofs

        \param numdofpernode: number of degress of freedom for one node
       */
      virtual void set_num_dof_per_node(int numdofpernode) { numdofpernode_ = numdofpernode; }

      int numdofpernode_;  ///< number of degrees of freedom

      //! action parameters recognized by bele3
      enum ActionType
      {
        none,
        calc_struct_constrvol,
        calc_struct_volconstrstiff,
        calc_struct_stress
      };

      /*!
       * \brief check, whether higher order derivatives for shape functions (dxdx, dxdy, ...) are
       * necessary \return boolean indicating higher order status
       */
      bool is_higher_order_element(const Core::FE::CellType distype) const
      {
        bool hoel = true;
        switch (distype)
        {
          case Core::FE::CellType::quad4:
          case Core::FE::CellType::quad8:
          case Core::FE::CellType::quad9:
          case Core::FE::CellType::tri6:
            hoel = true;
            break;
          case Core::FE::CellType::tri3:
            hoel = false;
            break;
          default:
            FOUR_C_THROW("distype unknown!");
            break;
        }
        return hoel;
      };

      /*!
        \brief Create matrix with spatial configuration

        \param x     (out)  : nodal coords in spatial frame
        \param disp  (int)  : displacements
      */
      inline void spatial_configuration(
          Core::LinAlg::SerialDenseMatrix& x, const std::vector<double> disp) const
      {
        const int numnode = num_node();
        for (int i = 0; i < numnode; ++i)
        {
          x(i, 0) = Nodes()[i]->X()[0] + disp[i * 3 + 0];
          x(i, 1) = Nodes()[i]->X()[1] + disp[i * 3 + 1];
          x(i, 2) = Nodes()[i]->X()[2] + disp[i * 3 + 2];
        }
        return;
      }

      //! Submethod to compute the enclosed volume for volume constraint boundary condition
      double compute_constr_vols(
          const Core::LinAlg::SerialDenseMatrix& xc,  ///< current configuration
          const int numnode                           ///< num nodes
      );

      //! Submethod to compute constraint volume and its first and second derivatives w.r.t. the
      //! displacements
      void compute_vol_deriv(const Core::LinAlg::SerialDenseMatrix& x,  ///< spatial configuration
          const int numnode,                                            ///< number of nodes
          const int ndof,                                        ///< number of degrees of freedom
          double& V,                                             ///< volume
          Teuchos::RCP<Core::LinAlg::SerialDenseVector> Vdiff,   ///< first derivative
          Teuchos::RCP<Core::LinAlg::SerialDenseMatrix> Vdiff2,  ///< second derivative
          const int minind = 0,  ///< minimal index to compute enclosed volume with
          const int maxind = 2   ///< maximal index to compute enclosed volume with
      );

      //! vector with line elements
      //  std::vector<Teuchos::RCP<Core::Elements::Element> >                      lines_;

      //! flag for fixed or moving boundary
      //  const bool                                      is_moving_;

      //! don't want = operator
      Bele3& operator=(const Bele3& old);

      //! set number of gauss points to element shape default
      Core::FE::GaussRule2D get_optimal_gaussrule() const;

    };  // class Bele3



    //=======================================================================
    //=======================================================================
    //=======================================================================
    //=======================================================================

    class Bele3LineType : public Core::Elements::ElementType
    {
     public:
      std::string Name() const override { return "Bele3LineType"; }

      static Bele3LineType& Instance();

      Teuchos::RCP<Core::Elements::Element> Create(const int id, const int owner) override;

      void nodal_block_information(
          Core::Elements::Element* dwele, int& numdf, int& dimns, int& nv, int& np) override
      {
      }

      Core::LinAlg::SerialDenseMatrix ComputeNullSpace(
          Core::Nodes::Node& node, const double* x0, const int numdof, const int dimnsp) override
      {
        Core::LinAlg::SerialDenseMatrix nullspace;
        FOUR_C_THROW("method ComputeNullSpace not implemented for element type bele3!");
        return nullspace;
      }

     private:
      static Bele3LineType instance_;
    };


    /*!
    \brief An element representing a line of a bele3 element

    */
    class Bele3Line : public Core::Elements::FaceElement
    {
     public:
      //! @name Constructors and destructors and related methods

      /*!
      \brief Standard Constructor

      \param id : A unique global id
      \param owner: Processor owning this line
      \param nnode: Number of nodes attached to this element
      \param nodeids: global ids of nodes attached to this element
      \param nodes: the discretizations map of nodes to build ptrs to nodes from
      \param parent: The parent fluid element of this line
      \param lline: the local line number of this line w.r.t. the parent element
      */
      Bele3Line(int id, int owner, int nnode, const int* nodeids, Core::Nodes::Node** nodes,
          Discret::ELEMENTS::Bele3* parent, const int lline);

      /*!
      \brief Copy Constructor

      Makes a deep copy of a Element

      */
      Bele3Line(const Bele3Line& old);

      Core::Elements::Element* Clone() const override;
      Core::FE::CellType Shape() const override;
      int UniqueParObjectId() const override
      {
        return Bele3LineType::Instance().UniqueParObjectId();
      }
      void pack(Core::Communication::PackBuffer& data) const override;
      void unpack(const std::vector<char>& data) override;


      //@}

      //! @name Access methods


      /*!
      \brief Get number of degrees of freedom of a certain node
             (implements pure virtual Core::Elements::Element)

      For this 3D boundary element, we have 3 displacements, if needed
      */
      int NumDofPerNode(const Core::Nodes::Node&) const override { return numdofpernode_; }

      int num_dof_per_element() const override { return 0; }

      void print(std::ostream& os) const override;

      Core::Elements::ElementType& ElementType() const override
      {
        return Bele3LineType::Instance();
      }

      //@}

      //! @name Evaluation

      int evaluate(Teuchos::ParameterList& params, Core::FE::Discretization& discretization,
          std::vector<int>& lm, Core::LinAlg::SerialDenseMatrix& elemat1,
          Core::LinAlg::SerialDenseMatrix& elemat2, Core::LinAlg::SerialDenseVector& elevec1,
          Core::LinAlg::SerialDenseVector& elevec2,
          Core::LinAlg::SerialDenseVector& elevec3) override;

      //! @name Evaluate methods

      int evaluate_neumann(Teuchos::ParameterList& params, Core::FE::Discretization& discretization,
          Core::Conditions::Condition& condition, std::vector<int>& lm,
          Core::LinAlg::SerialDenseVector& elevec1,
          Core::LinAlg::SerialDenseMatrix* elemat1 = nullptr) override;

      //@}

     private:
      /*!
        \brief Set number of dofs

        \param numdofpernode: number of degress of freedom for one node
       */
      virtual void set_num_dof_per_node(int numdofpernode) { numdofpernode_ = numdofpernode; }

      int numdofpernode_;  ///< number of degrees of freedom

      //! action parameters recognized by Bele3Line
      enum ActionType
      {
        none,
        integrate_Shapefunction
      };

      //! don't want = operator
      Bele3Line& operator=(const Bele3Line& old);


      //! compute infintesimal line element dr for integration along the line
      double f2_substitution(const Core::LinAlg::SerialDenseMatrix xye,
          const Core::LinAlg::SerialDenseMatrix deriv, const int iel);

      //! Get Rule for Gaussintegration according to Discret::UTIL
      Core::FE::GaussRule1D get_optimal_gaussrule(const Core::FE::CellType& distype);

      //! integrate shape functions over a line
      void integrate_shape_function(Teuchos::ParameterList& params,
          Core::FE::Discretization& discretization, const std::vector<int>& lm,
          Core::LinAlg::SerialDenseVector& elevec1, const std::vector<double>& edispnp);


    };  // class Bele3Line



  }  // namespace ELEMENTS
}  // namespace Discret


FOUR_C_NAMESPACE_CLOSE

#endif
