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

namespace DRT
{
  // forward declarations
  class So_sh8Plast;
  class Discretization;

  namespace ELEMENTS
  {
    class SoSh8PlastType : public SoHex8PlastType
    {
     public:
      std::string Name() const override { return "So_sh8PlastType"; }

      static SoSh8PlastType& Instance();

      CORE::COMM::ParObject* Create(const std::vector<char>& data) override;

      Teuchos::RCP<DRT::Element> Create(const std::string eletype, const std::string eledistype,
          const int id, const int owner) override;

      Teuchos::RCP<DRT::Element> Create(const int id, const int owner) override;

      int Initialize(DRT::Discretization& dis) override;

      void SetupElementDefinition(
          std::map<std::string, std::map<std::string, INPUT::LineDefinition>>& definitions)
          override;

     private:
      static SoSh8PlastType instance_;

      std::string GetElementTypeString() const { return "SOLIDSH8PLAST"; }
    };  // class So_sh8PlastType

    class SoSh8Plast : public virtual So3Plast<CORE::FE::CellType::hex8>
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
      DRT::Element* Clone() const override;


      //! Return unique ParObject id
      //!
      //! every class implementing ParObject needs a unique id defined at the top of
      //! this file.
      int UniqueParObjectId() const override
      {
        return SoSh8PlastType::Instance().UniqueParObjectId();
      }

      //! Pack this class so it can be communicated
      //! Pack and \ref Unpack are used to communicate this element
      void Pack(CORE::COMM::PackBuffer& data) const override;

      //! Unpack data from a char vector into this class
      //! Pack and \ref Unpack are used to communicate this element
      void Unpack(const std::vector<char>& data) override;

      //! Print this element
      void Print(std::ostream& os) const override;

      //! return elementtype
      SoSh8PlastType& ElementType() const override { return SoSh8PlastType::Instance(); }

      //! read input for this element
      bool ReadElement(const std::string& eletype, const std::string& distype,
          INPUT::LineDefinition* linedef) override;

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

      std::string GetElementTypeString() const { return "SOLIDSH8PLAST"; }

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

      static std::pair<bool, CORE::LINALG::Matrix<nsd_, nsd_>> jac_refe_;
      static std::pair<bool, CORE::LINALG::Matrix<nsd_, nsd_>> jac_curr_;
      static std::pair<bool, CORE::LINALG::Matrix<num_ans * num_sp, numdofperelement_>> B_ans_loc_;
      static std::pair<bool, CORE::LINALG::Matrix<numstr_, numstr_>> TinvT_;

      //! Calculate nonlinear stiffness and mass matrix with condensed plastic matrices
      void nln_stiffmass(std::vector<double>& disp,  // current displacements
          std::vector<double>& vel,                  // current velocities
          std::vector<double>& temp,                 // current temperatures
          CORE::LINALG::Matrix<numdofperelement_, numdofperelement_>*
              stiffmatrix,  // element stiffness matrix
          CORE::LINALG::Matrix<numdofperelement_, numdofperelement_>*
              massmatrix,                                         // element mass matrix
          CORE::LINALG::Matrix<numdofperelement_, 1>* force,      // element internal force vector
          CORE::LINALG::Matrix<numgpt_post, numstr_>* elestress,  // stresses at GP
          CORE::LINALG::Matrix<numgpt_post, numstr_>* elestrain,  // strains at GP
          Teuchos::ParameterList& params,         // algorithmic parameters e.g. time
          const INPAR::STR::StressType iostress,  // stress output option
          const INPAR::STR::StrainType iostrain   // strain output option
          ) override;

      //! calculate nonlinear B-operator (potentially with ANS modification)
      void CalculateBop(CORE::LINALG::Matrix<numstr_, numdofperelement_>* bop,
          const CORE::LINALG::Matrix<nsd_, nsd_>* defgrd,
          const CORE::LINALG::Matrix<nsd_, nen_>* N_XYZ, const int gp) override;

      //! Evaluate all ANS related data at the ANS sampling points
      void Anssetup(const CORE::LINALG::Matrix<nen_, nsd_>& xrefe,  ///< material element coords
          const CORE::LINALG::Matrix<nen_, nsd_>& xcurr,            ///< current element coords
          std::vector<CORE::LINALG::Matrix<nsd_, nen_>>**
              deriv_sp,  ///< derivs eval. at all sampling points
          std::vector<CORE::LINALG::Matrix<nsd_, nsd_>>& jac_sps,  ///< jac at all sampling points
          std::vector<CORE::LINALG::Matrix<nsd_, nsd_>>&
              jac_cur_sps,  ///< current jac at all sampling points
          CORE::LINALG::Matrix<num_ans * num_sp, numdofperelement_>& B_ans_loc);  ///< modified B

      //! Evaluate transformation matrix T (parameter->material) at gp
      void EvaluateT(const CORE::LINALG::Matrix<nsd_, nsd_>& jac,  ///< actual jacobian
          CORE::LINALG::Matrix<numstr_, numstr_>& TinvT);          ///< T^{-T}

      void AnsStrains(const int gp,
          std::vector<CORE::LINALG::Matrix<nsd_, nsd_>>& jac_sps,  // jac at all sampling points
          std::vector<CORE::LINALG::Matrix<nsd_, nsd_>>&
              jac_cur_sps  // current jac at all sampling points
      );

      //! Find "thin"=thickness direction
      ThicknessDirection findthickdir();

      //! Find parametric co-ordinate which directs in enforced thickness direction
      ThicknessDirection enfthickdir(CORE::LINALG::Matrix<nsd_, 1>&
              thickdirglo  ///< global direction of enforced thickness direction
      );

      //! Re-initialize EAS data, needed for sosh8 morphing
      void ReInitEas(const DRT::ELEMENTS::So3PlastEasType EASType);

      void InvalidGpData() override
      {
        So3Plast<CORE::FE::CellType::hex8>::InvalidGpData();
        jac_refe_.first = false;
        jac_curr_.first = false;
        TinvT_.first = false;
      }

      void InvalidEleData()
      {
        So3Plast<CORE::FE::CellType::hex8>::InvalidEleData();
        B_ans_loc_.first = false;
      }

      const CORE::LINALG::Matrix<nsd_, nsd_>& Jac_curr() const
      {
        if (jac_curr_.first == false) FOUR_C_THROW("jac_curr_ not valid");
        return jac_curr_.second;
      }
      CORE::LINALG::Matrix<nsd_, nsd_>& SetJac_curr()
      {
        jac_curr_.first = true;
        return jac_curr_.second;
      }

      const CORE::LINALG::Matrix<nsd_, nsd_>& Jac_refe() const
      {
        if (jac_refe_.first == false) FOUR_C_THROW("jac_refe_ not valid");
        return jac_refe_.second;
      }
      CORE::LINALG::Matrix<nsd_, nsd_>& SetJac_refe()
      {
        jac_refe_.first = true;
        return jac_refe_.second;
      }

      const CORE::LINALG::Matrix<num_ans * num_sp, numdofperelement_>& B_ans_loc() const
      {
        if (B_ans_loc_.first == false) FOUR_C_THROW("B_ans_loc_ not valid");
        return B_ans_loc_.second;
      }
      CORE::LINALG::Matrix<num_ans * num_sp, numdofperelement_>& SetB_ans_loc()
      {
        B_ans_loc_.first = true;
        return B_ans_loc_.second;
      }

      const CORE::LINALG::Matrix<numstr_, numstr_>& TinvT() const
      {
        if (TinvT_.first == false) FOUR_C_THROW("TinvT_ not valid");
        return TinvT_.second;
      }
      CORE::LINALG::Matrix<numstr_, numstr_>& SetTinvT()
      {
        TinvT_.first = true;
        return TinvT_.second;
      }
    };  // class So_sh8Plast
  }     // namespace ELEMENTS
}  // namespace DRT



FOUR_C_NAMESPACE_CLOSE

#endif
