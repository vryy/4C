/*----------------------------------------------------------------------*/
/*! \file
\brief 18-node hexahedral (bi-quadratic linear)
\level 1


*----------------------------------------------------------------------*/
#ifndef FOUR_C_SO3_HEX18_HPP
#define FOUR_C_SO3_HEX18_HPP

#include "4C_config.hpp"

#include "4C_fem_general_element.hpp"
#include "4C_fem_general_elementtype.hpp"
#include "4C_inpar_structure.hpp"
#include "4C_linalg_serialdensematrix.hpp"
#include "4C_mat_material_factory.hpp"
#include "4C_material_base.hpp"
#include "4C_so3_base.hpp"

FOUR_C_NAMESPACE_OPEN

// Several parameters which are fixed for Solid Hex18
const int NUMNOD_SOH18 = 18;         ///< number of nodes
const int NODDOF_SOH18 = 3;          ///< number of dofs per node
const int NUMDOF_SOH18 = 54;         ///< total dofs per element
const int NUMGPT_SOH18 = 2 * 3 * 3;  ///< total gauss points per element
const int NUMDIM_SOH18 = 3;          ///< number of dimensions

namespace Core::FE
{
  class Discretization;
}  // namespace Core::FE

namespace Discret
{
  class SoSh18Type;

  namespace ELEMENTS
  {
    class SoHex18Type : public Core::Elements::ElementType
    {
     public:
      std::string Name() const override { return "So_hex18Type"; }

      static SoHex18Type& Instance();

      Core::Communication::ParObject* Create(const std::vector<char>& data) override;

      Teuchos::RCP<Core::Elements::Element> Create(const std::string eletype,
          const std::string eledistype, const int id, const int owner) override;

      Teuchos::RCP<Core::Elements::Element> Create(const int id, const int owner) override;

      int Initialize(Core::FE::Discretization& dis) override;

      void nodal_block_information(
          Core::Elements::Element* dwele, int& numdf, int& dimns, int& nv, int& np) override;

      Core::LinAlg::SerialDenseMatrix ComputeNullSpace(
          Core::Nodes::Node& node, const double* x0, const int numdof, const int dimnsp) override;

      void setup_element_definition(
          std::map<std::string, std::map<std::string, Input::LineDefinition>>& definitions)
          override;

     private:
      static SoHex18Type instance_;

      std::string get_element_type_string() const { return "SOLIDH18"; }
    };  // class So_hex18Type

    class SoHex18 : public virtual SoBase
    {
     public:
      //! @name Friends
      friend class SoHex18Type;
      friend class SoSh18Type;

      //@}
      //! @name Constructors and destructors and related methods

      /*!
      \brief Standard Constructor

      \param id : A unique global id
      \param owner : elements owner
      */
      SoHex18(int id, int owner);

      /*!
      \brief Copy Constructor

      Makes a deep copy of a Element

      */
      SoHex18(const SoHex18& old);

      /*!
      \brief Deep copy this instance of Solid3 and return pointer to the copy

      The Clone() method is used from the virtual base class Element in cases
      where the type of the derived class is unknown and a copy-ctor is needed

      */
      Core::Elements::Element* Clone() const override;

      /*!
      \brief Get shape type of element
      */
      Core::FE::CellType Shape() const override { return Core::FE::CellType::hex18; };

      /*!
      \brief Return number of volumes of this element
      */
      int NumVolume() const override { return 1; }

      /*!
      \brief Return number of surfaces of this element
      */
      int NumSurface() const override { return 6; }

      /*!
      \brief Return number of lines of this element
      */
      int NumLine() const override { return 12; }

      /*!
      \brief Get vector of Teuchos::RCPs to the lines of this element

      */
      std::vector<Teuchos::RCP<Core::Elements::Element>> Lines() override;

      /*!
      \brief Get vector of Teuchos::RCPs to the surfaces of this element

      */
      std::vector<Teuchos::RCP<Core::Elements::Element>> Surfaces() override;

      /*!
      \brief Return unique ParObject id

      every class implementing ParObject needs a unique id defined at the
      top of this file.
      */
      int UniqueParObjectId() const override { return SoHex18Type::Instance().UniqueParObjectId(); }

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
      int NumDofPerNode(const Core::Nodes::Node& node) const override { return 3; }


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



      Core::Elements::ElementType& ElementType() const override { return SoHex18Type::Instance(); }

      //@}

      //! @name Input and Creation

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
        names.insert(std::pair<std::string,int>("Owner",1));
        // Name of data is 'StressesXYZ', dimension is 6 (sym. tensor value)
        names.insert(std::pair<std::string,int>("StressesXYZ",6));
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
      \brief Evaluate an element

      Evaluate so_hex8 element stiffness, mass, internal forces, etc.

      \param params (in/out): ParameterList for communication between control routine
                              and elements
      \param discretization : pointer to discretization for de-assembly
      \param lm (in)        : location matrix for de-assembly
      \param elemat1 (out)  : (stiffness-)matrix to be filled by element. If nullptr on input,
                              the controling method does not expect the element to fill
                              this matrix.
      \param elemat2 (out)  : (mass-)matrix to be filled by element. If nullptr on input,
                              the controling method does not expect the element to fill
                              this matrix.
      \param elevec1 (out)  : (internal force-)vector to be filled by element. If nullptr on input,
                              the controlling method does not expect the element
                              to fill this vector
      \param elevec2 (out)  : vector to be filled by element. If nullptr on input,
                              the controlling method does not expect the element
                              to fill this vector
      \param elevec3 (out)  : vector to be filled by element. If nullptr on input,
                              the controlling method does not expect the element
                              to fill this vector
      \return 0 if successful, negative otherwise
      */
      int evaluate(Teuchos::ParameterList& params, Core::FE::Discretization& discretization,
          std::vector<int>& lm, Core::LinAlg::SerialDenseMatrix& elemat1,
          Core::LinAlg::SerialDenseMatrix& elemat2, Core::LinAlg::SerialDenseVector& elevec1,
          Core::LinAlg::SerialDenseVector& elevec2,
          Core::LinAlg::SerialDenseVector& elevec3) override;

      /*!
      \brief Evaluate a Neumann boundary condition

      this method evaluates a surface Neumann condition on the solid3 element

      \param params (in/out)    : ParameterList for communication between control routine
                                  and elements
      \param discretization (in): A reference to the underlying discretization
      \param condition (in)     : The condition to be evaluated
      \param lm (in)            : location vector of this element
      \param elevec1 (out)      : vector to be filled by element. If nullptr on input,

      \return 0 if successful, negative otherwise
      */
      int evaluate_neumann(Teuchos::ParameterList& params, Core::FE::Discretization& discretization,
          Core::Conditions::Condition& condition, std::vector<int>& lm,
          Core::LinAlg::SerialDenseVector& elevec1,
          Core::LinAlg::SerialDenseMatrix* elemat1 = nullptr) override;

      /*!
      \brief Return value how expensive it is to evaluate this element

      \param double (out): cost to evaluate this element
      */
      double EvaluationCost() override
      {
        if (Material()->MaterialType() == Core::Materials::m_struct_multiscale)
          return 25000.0;
        else
          return 10.0;
      }
      //@}

     protected:
      //! action parameters recognized by so_hex8
      enum ActionType
      {
        none,
        calc_struct_linstiff,
        calc_struct_nlnstiff,
        calc_struct_internalforce,
        calc_struct_linstiffmass,
        calc_struct_nlnstiffmass,
        calc_struct_nlnstifflmass,  //!< internal force, its stiffness and lumped mass matrix
        calc_struct_stress,
        calc_struct_eleload,
        calc_struct_update_istep,
        calc_struct_reset_istep,  //!< reset elementwise internal variables
                                  //!< during iteration to last converged state
        calc_struct_reset_all,    //!< reset elementwise internal variables
                                  //!< to state in the beginning of the computation
        calc_struct_energy,       //!< compute internal energy
        calc_recover
      };

      //! vector of inverses of the jacobian in material frame
      std::vector<Core::LinAlg::Matrix<NUMDIM_SOH18, NUMDIM_SOH18>> invJ_;
      //! vector of inverses of the jacobian in material frame
      std::vector<double> detJ_;
      //! Gauss point coordinates
      std::vector<Core::LinAlg::Matrix<NUMDIM_SOH18, 1>> xsi_;
      //! Gauss point weights
      std::vector<double> wgt_;

      // internal calculation methods

      //! don't want = operator
      SoHex18& operator=(const SoHex18& old);


      //! init the inverse of the jacobian and its determinant in the material configuration
      virtual int init_jacobian_mapping();

      //! init the quadrature points
      virtual void init_gp();

      //! Calculate nonlinear stiffness and mass matrix
      virtual void nlnstiffmass(std::vector<int>& lm,  ///< location matrix
          std::vector<double>& disp,                   ///< current displacements
          std::vector<double>& residual,               ///< current residual displ
          Core::LinAlg::Matrix<NUMDOF_SOH18, NUMDOF_SOH18>*
              stiffmatrix,  ///< element stiffness matrix
          Core::LinAlg::Matrix<NUMDOF_SOH18, NUMDOF_SOH18>* massmatrix,  ///< element mass matrix
          Core::LinAlg::Matrix<NUMDOF_SOH18, 1>* force,  ///< element internal force vector
          Core::LinAlg::Matrix<NUMGPT_SOH18, Mat::NUM_STRESS_3D>* elestress,  ///< stresses at GP
          Core::LinAlg::Matrix<NUMGPT_SOH18, Mat::NUM_STRESS_3D>* elestrain,  ///< strains at GP
          Teuchos::ParameterList& params,         ///< algorithmic parameters e.g. time
          const Inpar::STR::StressType iostress,  ///< stress output option
          const Inpar::STR::StrainType iostrain   ///< strain output option
      );                                          ///< plastic strain output option

      //! Lump mass matrix (bborn 07/08)
      void lumpmass(Core::LinAlg::Matrix<NUMDOF_SOH18, NUMDOF_SOH18>* emass);

      void flip_t();

      // parameter space coords of one node
      Core::LinAlg::Matrix<3, 1> node_param_coord(const int node);
      // parameter space coords of all nodes
      Core::LinAlg::Matrix<18, 3> node_param_coord();

      // Update the element
      void update(){/* do nothing */};

      //@}


      /** recover elementwise stored stuff */
      virtual void recover(const std::vector<double>& residual) { return; }

     private:
      std::string get_element_type_string() const { return "SOLIDH18"; }
    };  // class So_hex18

  }  // namespace ELEMENTS
}  // namespace Discret



FOUR_C_NAMESPACE_CLOSE

#endif
