/*----------------------------------------------------------------------*/
/*! \file
\brief 18 node solid shell with plasticity
\level 3
*/
/*----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 | definitions                                              seitz 11/14 |
 *----------------------------------------------------------------------*/
#ifndef FOUR_C_SO3_PLAST_SSN_SOSH18_HPP
#define FOUR_C_SO3_PLAST_SSN_SOSH18_HPP

#include "4C_config.hpp"

#include "4C_so3_plast_ssn.hpp"
#include "4C_so3_plast_ssn_eletypes.hpp"
#include "4C_so3_sh18.hpp"

FOUR_C_NAMESPACE_OPEN

namespace Core::FE
{
  class Discretization;
}  // namespace Core::FE

namespace Discret
{
  // forward declarations
  class So_sh18Plast;

  namespace ELEMENTS
  {
    class SoSh18PlastType : public SoSh18Type
    {
     public:
      std::string name() const override { return "SoSh18PlastType"; }

      static SoSh18PlastType& instance();

      Core::Communication::ParObject* create(const std::vector<char>& data) override;

      Teuchos::RCP<Core::Elements::Element> create(const std::string eletype,
          const std::string eledistype, const int id, const int owner) override;

      Teuchos::RCP<Core::Elements::Element> create(const int id, const int owner) override;

      int initialize(Core::FE::Discretization& dis) override;

      void setup_element_definition(
          std::map<std::string, std::map<std::string, Input::LineDefinition>>& definitions)
          override;

     private:
      static SoSh18PlastType instance_;

      std::string get_element_type_string() const { return "SOLIDSH18PLAST"; }
    };  // class SoSh18PlastType

    class SoSh18Plast : public virtual So3Plast<Core::FE::CellType::hex18>, public virtual SoSh18
    {
     public:
      //! @name Friends
      friend class SoSh18PlastType;


      //! Standard Constructor
      SoSh18Plast(int id,  //!< (i) this element's global id
          int owner        //!< elements owner
      );

      //! Copy Constructor
      //! Makes a deep copy of a Element
      SoSh18Plast(const SoSh18Plast& old);

      bool have_eas() const override { return (eastype_ != soh8p_easnone); };

      //! resolve "no unique final overrider"
      int num_volume() const override { return SoSh18::num_volume(); }
      Core::FE::CellType shape() const override { return Core::FE::CellType::hex18; };
      int num_surface() const override { return SoSh18::num_surface(); }
      int num_line() const override { return SoSh18::num_line(); }
      std::vector<Teuchos::RCP<Core::Elements::Element>> lines() override
      {
        return SoSh18::lines();
      }
      std::vector<Teuchos::RCP<Core::Elements::Element>> surfaces() override
      {
        return SoSh18::surfaces();
      }
      int num_dof_per_node(const Core::Nodes::Node& node) const override
      {
        return SoSh18::num_dof_per_node(node);
      }
      int num_dof_per_element() const override { return SoSh18::num_dof_per_element(); }
      void vis_names(std::map<std::string, int>& names) override
      {
        return SoSh18::vis_names(names);
      }
      bool vis_data(const std::string& name, std::vector<double>& data) override
      {
        return SoSh18::vis_data(name, data);
      }
      int evaluate_neumann(Teuchos::ParameterList& params, Core::FE::Discretization& discretization,
          Core::Conditions::Condition& condition, std::vector<int>& lm,
          Core::LinAlg::SerialDenseVector& elevec1,
          Core::LinAlg::SerialDenseMatrix* elemat1 = nullptr) override
      {
        return SoSh18::evaluate_neumann(params, discretization, condition, lm, elevec1, elemat1);
      }

      //! Deep copy this instance of Solid3 and return pointer to the copy
      //!
      //! The clone() method is used from the virtual base class Element in cases
      //! where the type of the derived class is unknown and a copy-ctor is needed
      Core::Elements::Element* clone() const override;


      //! Return unique ParObject id
      //!
      //! every class implementing ParObject needs a unique id defined at the top of
      //! this file.
      int unique_par_object_id() const override
      {
        return SoSh18PlastType::instance().unique_par_object_id();
      }

      //! Pack this class so it can be communicated
      //! Pack and \ref unpack are used to communicate this element
      void pack(Core::Communication::PackBuffer& data) const override;

      //! Unpack data from a char vector into this class
      //! Pack and \ref unpack are used to communicate this element
      void unpack(const std::vector<char>& data) override;

      //! Print this element
      void print(std::ostream& os) const override;

      //! return elementtype
      SoSh18PlastType& element_type() const override { return SoSh18PlastType::instance(); }

      //! read input for this element
      bool read_element(const std::string& eletype, const std::string& distype,
          const Core::IO::InputParameterContainer& container) override;

      //! synchronize the eas variables in the two base-classes
      void sync_eas();

      //! evaluate an element
      //! evaluate element stiffness, mass, internal forces, etc.
      //!
      //! if nullptr on input, the controlling method does not expect the element
      //!  to fill these matrices or vectors.
      //!
      //!  \return 0 if successful, negative otherwise
      int evaluate(
          Teuchos::ParameterList&
              params,  //!< ParameterList for communication between control routine and elements
          Core::FE::Discretization& discretization,  //!< pointer to discretization for de-assembly
          Core::Elements::Element::LocationArray& la,  //!< location array for de-assembly
          Core::LinAlg::SerialDenseMatrix&
              elemat1_epetra,  //!< (stiffness-)matrix to be filled by element.
          Core::LinAlg::SerialDenseMatrix&
              elemat2_epetra,  //!< (mass-)matrix to be filled by element.
          Core::LinAlg::SerialDenseVector&
              elevec1_epetra,  //!< (internal force-)vector to be filled by element
          Core::LinAlg::SerialDenseVector& elevec2_epetra,  //!< vector to be filled by element
          Core::LinAlg::SerialDenseVector& elevec3_epetra   //!< vector to be filled by element
          ) override
      {
        return Discret::ELEMENTS::So3Plast<Core::FE::CellType::hex18>::evaluate(params,
            discretization, la, elemat1_epetra, elemat2_epetra, elevec1_epetra, elevec2_epetra,
            elevec3_epetra);
      }


     private:
      // don't want = operator
      SoSh18Plast& operator=(const SoSh18Plast& old) = delete;

      std::string get_element_type_string() const { return "SOLIDSH18PLAST"; }

     protected:
      //! Calculate nonlinear stiffness and mass matrix with condensed plastic matrices
      void nln_stiffmass(std::vector<double>& disp,  // current displacements
          std::vector<double>& vel,                  // current velocities
          std::vector<double>& temp,                 // current temperatures
          Core::LinAlg::Matrix<numdofperelement_, numdofperelement_>*
              stiffmatrix,  // element stiffness matrix
          Core::LinAlg::Matrix<numdofperelement_, numdofperelement_>*
              massmatrix,                                         // element mass matrix
          Core::LinAlg::Matrix<numdofperelement_, 1>* force,      // element internal force vector
          Core::LinAlg::Matrix<numgpt_post, numstr_>* elestress,  // stresses at GP
          Core::LinAlg::Matrix<numgpt_post, numstr_>* elestrain,  // strains at GP
          Teuchos::ParameterList& params,           // algorithmic parameters e.g. time
          const Inpar::Solid::StressType iostress,  // stress output option
          const Inpar::Solid::StrainType iostrain   // strain output option
          ) override;

      //! don't want sosh18 nlnstiffmass
      void nlnstiffmass(std::vector<int>& lm,  ///< location matrix
          std::vector<double>& disp,           ///< current displacements
          std::vector<double>& residual,       ///< current residual displ
          Core::LinAlg::Matrix<NUMDOF_SOH18, NUMDOF_SOH18>*
              stiffmatrix,  ///< element stiffness matrix
          Core::LinAlg::Matrix<NUMDOF_SOH18, NUMDOF_SOH18>* massmatrix,  ///< element mass matrix
          Core::LinAlg::Matrix<NUMDOF_SOH18, 1>* force,  ///< element internal force vector
          Core::LinAlg::Matrix<NUMGPT_SOH18, Mat::NUM_STRESS_3D>* elestress,  ///< stresses at GP
          Core::LinAlg::Matrix<NUMGPT_SOH18, Mat::NUM_STRESS_3D>* elestrain,  ///< strains at GP
          Teuchos::ParameterList& params,           ///< algorithmic parameters e.g. time
          const Inpar::Solid::StressType iostress,  ///< stress output option
          const Inpar::Solid::StrainType iostrain   ///< strain output option
          ) override
      {
        FOUR_C_THROW("don't want this");
      }
    };

  }  // namespace ELEMENTS



}  // namespace Discret

FOUR_C_NAMESPACE_CLOSE

#endif
