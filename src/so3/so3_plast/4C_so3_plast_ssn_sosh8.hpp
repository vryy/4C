/*----------------------------------------------------------------------*/
/*! \file
\level 2
\brief elasto-plastic solid shell
*/


/*----------------------------------------------------------------------*
 | definitions                                              seitz 05/14 |
 *----------------------------------------------------------------------*/
#ifndef FOUR_C_SO3_PLAST_SSN_SOSH8_HPP
#define FOUR_C_SO3_PLAST_SSN_SOSH8_HPP

/*----------------------------------------------------------------------*
 | headers                                                  seitz 05/14 |
 *----------------------------------------------------------------------*/
#include "4C_config.hpp"

#include "4C_so3_plast_ssn.hpp"
#include "4C_so3_plast_ssn_eletypes.hpp"

FOUR_C_NAMESPACE_OPEN

namespace Core::FE
{
  class Discretization;
}  // namespace Core::FE

namespace Discret
{
  // forward declarations
  class So_sh8Plast;
  namespace ELEMENTS
  {
    class SoSh8PlastType : public SoHex8PlastType
    {
     public:
      std::string Name() const override { return "So_sh8PlastType"; }

      static SoSh8PlastType& Instance();

      Core::Communication::ParObject* Create(const std::vector<char>& data) override;

      Teuchos::RCP<Core::Elements::Element> Create(const std::string eletype,
          const std::string eledistype, const int id, const int owner) override;

      Teuchos::RCP<Core::Elements::Element> Create(const int id, const int owner) override;

      int initialize(Core::FE::Discretization& dis) override;

      void setup_element_definition(
          std::map<std::string, std::map<std::string, Input::LineDefinition>>& definitions)
          override;

     private:
      static SoSh8PlastType instance_;

      std::string get_element_type_string() const { return "SOLIDSH8PLAST"; }
    };  // class So_sh8PlastType

    class SoSh8Plast : public virtual So3Plast<Core::FE::CellType::hex8>
    {
     public:
      //! @name Friends
      friend class SoSh8PlastType;


      //! Standard Constructor
      SoSh8Plast(int id,  //!< (i) this element's global id
          int owner       //!< elements owner
      );

      //! Copy Constructor
      //! Makes a deep copy of a Element
      SoSh8Plast(const SoSh8Plast& old);

      //! Deep copy this instance of Solid3 and return pointer to the copy
      //!
      //! The Clone() method is used from the virtual base class Element in cases
      //! where the type of the derived class is unknown and a copy-ctor is needed
      Core::Elements::Element* Clone() const override;


      //! Return unique ParObject id
      //!
      //! every class implementing ParObject needs a unique id defined at the top of
      //! this file.
      int UniqueParObjectId() const override
      {
        return SoSh8PlastType::Instance().UniqueParObjectId();
      }

      //! Pack this class so it can be communicated
      //! Pack and \ref unpack are used to communicate this element
      void pack(Core::Communication::PackBuffer& data) const override;

      //! Unpack data from a char vector into this class
      //! Pack and \ref unpack are used to communicate this element
      void unpack(const std::vector<char>& data) override;

      //! Print this element
      void Print(std::ostream& os) const override;

      //! return elementtype
      SoSh8PlastType& ElementType() const override { return SoSh8PlastType::Instance(); }

      //! read input for this element
      bool ReadElement(const std::string& eletype, const std::string& distype,
          Input::LineDefinition* linedef) override;

      //! definition of shell-thickness direction
      enum ThicknessDirection
      {
        globx,      ///< global x
        globy,      ///< global y
        globz,      ///< global z
        autoj,      ///< find automatically by Jacobian
        autor,      ///< automatically set to x
        autos,      ///< automatically set to y
        autot,      ///< automatically set to z
        enfor,      ///< read-in r-direction is rearranged to t-dir
        enfos,      ///< read-in s-direction is rearranged to t-dir
        enfot,      ///< read-in t-direction stays t-dir
        undefined,  ///< no clear direction identified
        none        ///< no rearrangement
      };

      enum ANSType
      {
        anssosh8_p,
        ansnone_p
      };

     private:
      // don't want = operator
      SoSh8Plast& operator=(const SoSh8Plast& old) = delete;

      std::string get_element_type_string() const { return "SOLIDSH8PLAST"; }

     protected:
      static constexpr int num_sp = 8;  ///< number of ANS sampling points, here 8
      static constexpr int num_ans =
          3;  ///< number of modified ANS strains (E_rt,E_st,E_tt), here 3
      //! shell-thickness direction
      ThicknessDirection thickdir_;

      ANSType anstype_;

      //! in case of changed "thin" direction this is true
      bool nodes_rearranged_{};

      //! vector in thickness direction for compatibility with sosh8
      std::vector<double> thickvec_;

      static std::pair<bool, Core::LinAlg::Matrix<nsd_, nsd_>> jac_refe_;
      static std::pair<bool, Core::LinAlg::Matrix<nsd_, nsd_>> jac_curr_;
      static std::pair<bool, Core::LinAlg::Matrix<num_ans * num_sp, numdofperelement_>> B_ans_loc_;
      static std::pair<bool, Core::LinAlg::Matrix<numstr_, numstr_>> TinvT_;

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
          Teuchos::ParameterList& params,         // algorithmic parameters e.g. time
          const Inpar::STR::StressType iostress,  // stress output option
          const Inpar::STR::StrainType iostrain   // strain output option
          ) override;

      //! calculate nonlinear B-operator (potentially with ANS modification)
      void calculate_bop(Core::LinAlg::Matrix<numstr_, numdofperelement_>* bop,
          const Core::LinAlg::Matrix<nsd_, nsd_>* defgrd,
          const Core::LinAlg::Matrix<nsd_, nen_>* N_XYZ, const int gp) override;

      //! Evaluate all ANS related data at the ANS sampling points
      void anssetup(const Core::LinAlg::Matrix<nen_, nsd_>& xrefe,  ///< material element coords
          const Core::LinAlg::Matrix<nen_, nsd_>& xcurr,            ///< current element coords
          std::vector<Core::LinAlg::Matrix<nsd_, nen_>>**
              deriv_sp,  ///< derivs eval. at all sampling points
          std::vector<Core::LinAlg::Matrix<nsd_, nsd_>>& jac_sps,  ///< jac at all sampling points
          std::vector<Core::LinAlg::Matrix<nsd_, nsd_>>&
              jac_cur_sps,  ///< current jac at all sampling points
          Core::LinAlg::Matrix<num_ans * num_sp, numdofperelement_>& B_ans_loc);  ///< modified B

      //! Evaluate transformation matrix T (parameter->material) at gp
      void evaluate_t(const Core::LinAlg::Matrix<nsd_, nsd_>& jac,  ///< actual jacobian
          Core::LinAlg::Matrix<numstr_, numstr_>& TinvT);           ///< T^{-T}

      void ans_strains(const int gp,
          std::vector<Core::LinAlg::Matrix<nsd_, nsd_>>& jac_sps,  // jac at all sampling points
          std::vector<Core::LinAlg::Matrix<nsd_, nsd_>>&
              jac_cur_sps  // current jac at all sampling points
      );

      //! Find "thin"=thickness direction
      ThicknessDirection findthickdir();

      //! Find parametric co-ordinate which directs in enforced thickness direction
      ThicknessDirection enfthickdir(Core::LinAlg::Matrix<nsd_, 1>&
              thickdirglo  ///< global direction of enforced thickness direction
      );

      //! Re-initialize EAS data, needed for sosh8 morphing
      void re_init_eas(const Discret::ELEMENTS::So3PlastEasType EASType);

      void invalid_gp_data() override
      {
        So3Plast<Core::FE::CellType::hex8>::invalid_gp_data();
        jac_refe_.first = false;
        jac_curr_.first = false;
        TinvT_.first = false;
      }

      void invalid_ele_data()
      {
        So3Plast<Core::FE::CellType::hex8>::invalid_ele_data();
        B_ans_loc_.first = false;
      }

      const Core::LinAlg::Matrix<nsd_, nsd_>& jac_curr() const
      {
        if (jac_curr_.first == false) FOUR_C_THROW("jac_curr_ not valid");
        return jac_curr_.second;
      }
      Core::LinAlg::Matrix<nsd_, nsd_>& set_jac_curr()
      {
        jac_curr_.first = true;
        return jac_curr_.second;
      }

      const Core::LinAlg::Matrix<nsd_, nsd_>& jac_refe() const
      {
        if (jac_refe_.first == false) FOUR_C_THROW("jac_refe_ not valid");
        return jac_refe_.second;
      }
      Core::LinAlg::Matrix<nsd_, nsd_>& set_jac_refe()
      {
        jac_refe_.first = true;
        return jac_refe_.second;
      }

      const Core::LinAlg::Matrix<num_ans * num_sp, numdofperelement_>& b_ans_loc() const
      {
        if (B_ans_loc_.first == false) FOUR_C_THROW("B_ans_loc_ not valid");
        return B_ans_loc_.second;
      }
      Core::LinAlg::Matrix<num_ans * num_sp, numdofperelement_>& set_b_ans_loc()
      {
        B_ans_loc_.first = true;
        return B_ans_loc_.second;
      }

      const Core::LinAlg::Matrix<numstr_, numstr_>& tinv_t() const
      {
        if (TinvT_.first == false) FOUR_C_THROW("TinvT_ not valid");
        return TinvT_.second;
      }
      Core::LinAlg::Matrix<numstr_, numstr_>& set_tinv_t()
      {
        TinvT_.first = true;
        return TinvT_.second;
      }
    };  // class So_sh8Plast
  }     // namespace ELEMENTS
}  // namespace Discret



FOUR_C_NAMESPACE_CLOSE

#endif
