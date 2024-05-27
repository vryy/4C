/*----------------------------------------------------------------------*/
/*! \file

\brief implementation of the quadratic NURBS 27 element

\level 2


*----------------------------------------------------------------------*/
#ifndef FOUR_C_SO3_NURBS27_HPP
#define FOUR_C_SO3_NURBS27_HPP


#include "4C_config.hpp"

#include "4C_comm_parobjectfactory.hpp"
#include "4C_inpar_structure.hpp"
#include "4C_lib_element.hpp"
#include "4C_lib_elementtype.hpp"
#include "4C_lib_node.hpp"
#include "4C_linalg_fixedsizematrix.hpp"
#include "4C_linalg_serialdensematrix.hpp"
#include "4C_linalg_serialdensevector.hpp"
#include "4C_mat_material_factory.hpp"
#include "4C_material_base.hpp"
#include "4C_so3_base.hpp"

#include <Epetra_Vector.h>
#include <Teuchos_RCP.hpp>

FOUR_C_NAMESPACE_OPEN

/// Several parameters which are fixed for Solid nurbs27
const int NUMNOD_SONURBS27 = 27;  ///< number of nodes
const int NODDOF_SONURBS27 = 3;   ///< number of dofs per node
const int NUMDOF_SONURBS27 = 81;  ///< total dofs per element
const int NUMGPT_SONURBS27 = 27;  ///< total gauss points per element
const int NUMDIM_SONURBS27 = 3;   ///< number of dimensions


namespace DRT
{
  // forward declarations
  class Discretization;

  namespace ELEMENTS
  {
    namespace NURBS
    {
      class SoNurbs27Type : public DRT::ElementType
      {
       public:
        std::string Name() const override { return "So_nurbs27Type"; }

        static SoNurbs27Type& Instance();

        CORE::COMM::ParObject* Create(const std::vector<char>& data) override;

        Teuchos::RCP<DRT::Element> Create(const std::string eletype, const std::string eledistype,
            const int id, const int owner) override;

        Teuchos::RCP<DRT::Element> Create(const int id, const int owner) override;

        int Initialize(DRT::Discretization& dis) override;

        void nodal_block_information(
            DRT::Element* dwele, int& numdf, int& dimns, int& nv, int& np) override;

        CORE::LINALG::SerialDenseMatrix ComputeNullSpace(
            DRT::Node& node, const double* x0, const int numdof, const int dimnsp) override;

        void setup_element_definition(
            std::map<std::string, std::map<std::string, INPUT::LineDefinition>>& definitions)
            override;

       private:
        static SoNurbs27Type instance_;

        std::string get_element_type_string() const { return "SONURBS27"; }
      };


      /*!
      \brief A C++ version of a 27-control point nurbs solid element

      A structural 27-control points nurbs solid displacement element for large deformations.
      As its discretization is fixed many data structures are evaluated just once and kept
      for performance.

      \author gammi
      */
      class SoNurbs27 : public SoBase
      {
       public:
        //! @name Friends
        friend class SoNurbs27Type;

        //@}
        //! @name Constructors and destructors and related methods

        /*!
        \brief Standard Constructor

        \param id : A unique global id
        \param owner : elements owner
        */
        SoNurbs27(int id, int owner);

        /*!
        \brief Copy Constructor

        Makes a deep copy of a Element

        */
        SoNurbs27(const SoNurbs27& old);

        /*!
        \brief Deep copy this instance of Solid3 and return pointer to the copy

        The Clone() method is used from the virtual base class Element in cases
        where the type of the derived class is unknown and a copy-ctor is needed

        */
        DRT::Element* Clone() const override;

        /*!
        \brief Get shape type of element
        */
        CORE::FE::CellType Shape() const override;

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
        std::vector<Teuchos::RCP<DRT::Element>> Lines() override;

        /*!
        \brief Get vector of Teuchos::RCPs to the surfaces of this element

        */
        std::vector<Teuchos::RCP<DRT::Element>> Surfaces() override;

        /*!
        \brief Return unique ParObject id

        every class implementing ParObject needs a unique id defined at the
        top of this file.
        */
        int UniqueParObjectId() const override
        {
          return SoNurbs27Type::Instance().UniqueParObjectId();
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


        //@}

        //! @name Access methods


        /*!
        \brief Get number of degrees of freedom of a certain node
               (implements pure virtual DRT::Element)

        The element decides how many degrees of freedom its nodes must have.
        As this may vary along a simulation, the element can redecide the
        number of degrees of freedom per node along the way for each of it's nodes
        separately.
        */
        int NumDofPerNode(const DRT::Node& node) const override { return 3; }

        /*!
        \brief Get number of degrees of freedom per element
               (implements pure virtual DRT::Element)

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

        DRT::ElementType& ElementType() const override { return SoNurbs27Type::Instance(); }

        //! @name Input and Creation

        /*!
        \brief Read input for this element
        */
        bool ReadElement(const std::string& eletype, const std::string& distype,
            INPUT::LineDefinition* linedef) override;

        //@}

        //! @name Evaluation

        /*!
        \brief Evaluate an element

        Evaluate So_nurbs27 element stiffness, mass, internal forces, etc.

        If nullptr on input, the controling method does not expect the element
        to fill these matrices or vectors.

        \return 0 if successful, negative otherwise
        */
        int Evaluate(
            Teuchos::ParameterList&
                params,  ///< ParameterList for communication between control routine and elements
            DRT::Discretization& discretization,  ///< pointer to discretization for de-assembly
            std::vector<int>& lm,                 ///< location matrix for de-assembly
            CORE::LINALG::SerialDenseMatrix&
                elemat1,  ///< (stiffness-)matrix to be filled by element.
            CORE::LINALG::SerialDenseMatrix& elemat2,  ///< (mass-)matrix to be filled by element.
            CORE::LINALG::SerialDenseVector&
                elevec1,  ///< (internal force-)vector to be filled by element
            CORE::LINALG::SerialDenseVector& elevec2,  ///< vector to be filled by element
            CORE::LINALG::SerialDenseVector& elevec3   ///< vector to be filled by element
            ) override;


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
        int evaluate_neumann(Teuchos::ParameterList& params, DRT::Discretization& discretization,
            CORE::Conditions::Condition& condition, std::vector<int>& lm,
            CORE::LINALG::SerialDenseVector& elevec1,
            CORE::LINALG::SerialDenseMatrix* elemat1 = nullptr) override;

        //@}

        /*!
         \brief calculate a. scaled director matrix for thin shell-like structures
           or its inverse

         \return void
         */
        void do_calc_stc_matrix(CORE::LINALG::Matrix<81, 81>& elemat1,
            const INPAR::STR::StcScale stc_scaling, const int stc_layer, std::vector<int>& lm,
            DRT::Discretization& discretization, bool do_inverse);

       protected:
        //! action parameters recognized by So_nurbs27
        enum ActionType
        {
          none,
          calc_struct_linstiff,
          calc_struct_nlnstiff,
          calc_struct_internalforce,
          calc_struct_linstiffmass,
          calc_struct_nlnstiffmass,
          calc_struct_eleload,
          calc_struct_fsiload,
          calc_struct_update_istep,
          calc_stc_matrix,  //! calculate scaled director matrix for thin shell-like structures
          calc_stc_matrix_inverse,  //! calculate inverse scaled director matrix for thin shell-like
                                    //! structures
          calc_struct_reset_istep,
          calc_struct_energy,
          calc_struct_stifftemp,  //!< TSI specific: mechanical-thermal stiffness
          calc_struct_nlnstifflmass
        };

        //! vector of inverses of the jacobian in material frame
        std::vector<CORE::LINALG::Matrix<3, 3>> invJ_;
        //! determinant of Jacobian in material frame
        std::vector<double> detJ_;


        // internal calculation methods

        //! don't want = operator
        SoNurbs27& operator=(const SoNurbs27& old);

        //! init the inverse of the jacobian and its determinant in
        //! the material configuration
        virtual void init_jacobian_mapping(DRT::Discretization& dis);

        //! Calculate nonlinear stiffness and mass matrix
        virtual void sonurbs27_nlnstiffmass(std::vector<int>& lm,  ///< location matrix
            DRT::Discretization& discretization,        ///< discretisation to extract knot vector
            std::vector<double>& disp,                  ///< current displacements
            std::vector<double>& residual,              ///< current residual displ
            CORE::LINALG::Matrix<81, 81>* stiffmatrix,  ///< element stiffness matrix
            CORE::LINALG::Matrix<81, 81>* massmatrix,   ///< element mass matrix
            CORE::LINALG::Matrix<81, 1>* force,         ///< element internal force vector
            Teuchos::ParameterList& params);

        //! Evaluate Nurbs27 Shapefcts to keep them static
        std::vector<CORE::LINALG::Matrix<27, 1>> sonurbs27_shapefcts(
            const std::vector<CORE::LINALG::SerialDenseVector>& myknots,
            const CORE::LINALG::Matrix<27, 1>& weights);
        //! Evaluate Nurbs27 Derivs to keep them static
        std::vector<CORE::LINALG::Matrix<3, 27>> sonurbs27_derivs(
            const std::vector<CORE::LINALG::SerialDenseVector>& myknots,
            const CORE::LINALG::Matrix<27, 1>& weights);
        //! Evaluate Nurbs27 Weights to keep them static
        std::vector<double> sonurbs27_gpweights();

        //! calculate internal energy
        double calc_int_energy(DRT::Discretization& discretization, std::vector<double>& disp,
            Teuchos::ParameterList& params);

        //! Lump mass matrix (bborn 07/08)
        void lumpmass(CORE::LINALG::Matrix<NUMDOF_SONURBS27, NUMDOF_SONURBS27>* emass);

       private:
        std::string get_element_type_string() const { return "SONURBS27"; }
      };  // class So_nurbs27



    }  // namespace NURBS
  }    // namespace ELEMENTS
}  // namespace DRT

FOUR_C_NAMESPACE_CLOSE

#endif
