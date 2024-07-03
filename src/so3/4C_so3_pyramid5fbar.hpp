/*----------------------------------------------------------------------*/
/*! \file
\brief pyramid shaped solid element
\level 2


*----------------------------------------------------------------------*/
#ifndef FOUR_C_SO3_PYRAMID5FBAR_HPP
#define FOUR_C_SO3_PYRAMID5FBAR_HPP

#include "4C_config.hpp"

#include "4C_linalg_serialdensematrix.hpp"
#include "4C_so3_pyramid5.hpp"

FOUR_C_NAMESPACE_OPEN

namespace Core::FE
{
  class Discretization;
}  // namespace Core::FE

namespace Discret
{
  namespace ELEMENTS
  {
    // forward declarations
    class SoPyramid5fbarType : public Core::Elements::ElementType
    {
     public:
      std::string Name() const override { return "So_pyramid5fbarType"; }

      static SoPyramid5fbarType& Instance();

      Core::Communication::ParObject* Create(const std::vector<char>& data) override;

      Teuchos::RCP<Core::Elements::Element> Create(const std::string eletype,
          const std::string eledistype, const int id, const int owner) override;

      Teuchos::RCP<Core::Elements::Element> Create(const int id, const int owner) override;

      int initialize(Core::FE::Discretization& dis) override;

      void nodal_block_information(
          Core::Elements::Element* dwele, int& numdf, int& dimns, int& nv, int& np) override;

      Core::LinAlg::SerialDenseMatrix ComputeNullSpace(
          Core::Nodes::Node& node, const double* x0, const int numdof, const int dimnsp) override;

      void setup_element_definition(
          std::map<std::string, std::map<std::string, Input::LineDefinition>>& definitions)
          override;

     private:
      static SoPyramid5fbarType instance_;

      std::string get_element_type_string() const { return "SOLIDP5FBAR_DEPRECATED"; }
    };



    class SoPyramid5fbar : public SoPyramid5
    {
     public:
      //! @name Friends
      friend class SoPyramid5fbarType;

      //@}
      //! @name Constructors and destructors and related methods

      /*!
      \brief Standard Constructor

      \param id : A unique global id
      \param owner : elements owner
      */
      SoPyramid5fbar(int id, int owner);

      /*!
      \brief Copy Constructor

      Makes a deep copy of a Element

      */
      SoPyramid5fbar(const SoPyramid5fbar& old);

      /*!
      \brief Deep copy this instance of Solid3 and return pointer to the copy

      The Clone() method is used from the virtual base class Element in cases
      where the type of the derived class is unknown and a copy-ctor is needed

      */
      Core::Elements::Element* Clone() const override;

      /*!
      \brief Return unique ParObject id

      every class implementing ParObject needs a unique id defined at the
      top of this file.
      */
      int UniqueParObjectId() const override
      {
        return SoPyramid5fbarType::Instance().UniqueParObjectId();
      }

      /*!
      \brief Pack this class so it can be communicated

      \ref pack and \ref unpack are used to communicate this element

      */
      void pack(Core::Communication::PackBuffer& data) const override;

      /*!
      \brief Unpack data from a char vector into this class

      \ref pack and \ref unpack are used to communicate this element

      */
      void unpack(const std::vector<char>& data) override;


      //@}

      //! @name Access methods

      /*!
      \brief Print this element
      */
      void print(std::ostream& os) const override;

      Core::Elements::ElementType& ElementType() const override
      {
        return SoPyramid5fbarType::Instance();
      }

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

      Evaluate So_pyramid5fbar element stiffness, mass, internal forces, etc.

      If nullptr on input, the controlling method does not expect the element
      to fill these matrices or vectors.

      \return 0 if successful, negative otherwise
      */
      int evaluate(
          Teuchos::ParameterList&
              params,  ///< ParameterList for communication between control routine and elements
          Core::FE::Discretization& discretization,  ///< pointer to discretization for de-assembly
          std::vector<int>& lm,                      ///< location matrix for de-assembly
          Core::LinAlg::SerialDenseMatrix&
              elemat1,  ///< (stiffness-)matrix to be filled by element.
          Core::LinAlg::SerialDenseMatrix& elemat2,  ///< (mass-)matrix to be filled by element.
          Core::LinAlg::SerialDenseVector&
              elevec1,  ///< (internal force-)vector to be filled by element
          Core::LinAlg::SerialDenseVector& elevec2,  ///< vector to be filled by element
          Core::LinAlg::SerialDenseVector& elevec3   ///< vector to be filled by element
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
      int evaluate_neumann(Teuchos::ParameterList& params, Core::FE::Discretization& discretization,
          Core::Conditions::Condition& condition, std::vector<int>& lm,
          Core::LinAlg::SerialDenseVector& elevec1,
          Core::LinAlg::SerialDenseMatrix* elemat1 = nullptr) override;

      //@}

     protected:
      //! don't want = operator
      SoPyramid5fbar& operator=(const SoPyramid5fbar& old);

      // compute Jacobian mapping wrt to deformed configuration
      virtual void update_jacobian_mapping(
          const std::vector<double>& disp, Discret::ELEMENTS::PreStress& prestress);
      // compute defgrd in all gp for given disp
      virtual void def_gradient(const std::vector<double>& disp,
          Core::LinAlg::SerialDenseMatrix& gpdefgrd, Discret::ELEMENTS::PreStress& prestress);

      //! Calculate nonlinear stiffness and mass matrix
      virtual void nlnstiffmass(std::vector<int>& lm,  ///< location matrix
          std::vector<double>& disp,                   ///< current displacements
          std::vector<double>* vel,                    ///< current velocities
          std::vector<double>* acc,                    ///< current accelerations
          std::vector<double>& residual,               ///< current residual displ
          std::vector<double>& dispmat,                ///< current material displacements
          Core::LinAlg::Matrix<NUMDOF_SOP5, NUMDOF_SOP5>*
              stiffmatrix,                                             ///< element stiffness matrix
          Core::LinAlg::Matrix<NUMDOF_SOP5, NUMDOF_SOP5>* massmatrix,  ///< element mass matrix
          Core::LinAlg::Matrix<NUMDOF_SOP5, 1>* force,       ///< element internal force vector
          Core::LinAlg::Matrix<NUMDOF_SOP5, 1>* forceinert,  ///< element inertial force vector
          Core::LinAlg::Matrix<NUMDOF_SOP5, 1>* force_str,   ///< element structural force vector
          Core::LinAlg::Matrix<NUMGPT_SOP5, Mat::NUM_STRESS_3D>* elestress,  ///< stresses at GP
          Core::LinAlg::Matrix<NUMGPT_SOP5, Mat::NUM_STRESS_3D>* elestrain,  ///< strains at GP
          Core::LinAlg::Matrix<NUMGPT_SOP5, Mat::NUM_STRESS_3D>*
              eleplstrain,                             ///< plastic strains at GP
          Teuchos::ParameterList& params,              ///< algorithmic parameters e.g. time
          const Inpar::Solid::StressType iostress,     ///< stress output option
          const Inpar::Solid::StrainType iostrain,     ///< strain output option
          const Inpar::Solid::StrainType ioplstrain);  ///< plastic strain output option

      //! init the inverse of the jacobian and its determinant in the material configuration
      void init_jacobian_mapping() override;

     private:
      std::string get_element_type_string() const { return "SOLIDP5FBAR_DEPRECATED"; }
    };  // class So_pyramid5fbar

    //=======================================================================
    //=======================================================================
    //=======================================================================
    //=======================================================================

  }  // namespace ELEMENTS
}  // namespace Discret

FOUR_C_NAMESPACE_CLOSE

#endif
