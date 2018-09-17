/*----------------------------------------------------------------------*/
/*!
\file so3_ssn_plast_sosh8.H
\maintainer Matthias Mayr
\level 2
\brief elasto-plastic solid shell
*/


/*----------------------------------------------------------------------*
 | definitions                                              seitz 05/14 |
 *----------------------------------------------------------------------*/
#ifndef SO3_SSN_PLAST_SOSH8_H
#define SO3_SSN_PLAST_SOSH8_H

/*----------------------------------------------------------------------*
 | headers                                                  seitz 05/14 |
 *----------------------------------------------------------------------*/
#include "so3_ssn_plast.H"
#include "so3_ssn_plast_eletypes.H"

namespace DRT
{
  // forward declarations
  class So_sh8Plast;
  class Discretization;

  namespace ELEMENTS
  {
    class So_sh8PlastType : public So_hex8PlastType
    {
     public:
      std::string Name() const { return "So_sh8PlastType"; }

      static So_sh8PlastType& Instance();

      virtual DRT::ParObject* Create(const std::vector<char>& data);

      virtual Teuchos::RCP<DRT::Element> Create(
          const std::string eletype, const std::string eledistype, const int id, const int owner);

      virtual Teuchos::RCP<DRT::Element> Create(const int id, const int owner);

      virtual int Initialize(DRT::Discretization& dis);

      virtual void SetupElementDefinition(
          std::map<std::string, std::map<std::string, DRT::INPUT::LineDefinition>>& definitions);

     private:
      static So_sh8PlastType instance_;
    };  // class So_sh8PlastType

    class So_sh8Plast : public virtual So3_Plast<DRT::Element::hex8>
    {
     public:
      //! @name Friends
      friend class So_sh8PlastType;


      //! Standard Constructor
      So_sh8Plast(int id,  //!< (i) this element's global id
          int owner        //!< elements owner
      );

      //! Copy Constructor
      //! Makes a deep copy of a Element
      So_sh8Plast(const So_sh8Plast& old);

      //! Deep copy this instance of Solid3 and return pointer to the copy
      //!
      //! The Clone() method is used from the virtual base class Element in cases
      //! where the type of the derived class is unknown and a copy-ctor is needed
      DRT::Element* Clone() const;

      //! Destructor
      virtual ~So_sh8Plast();

      //! Return unique ParObject id
      //!
      //! every class implementing ParObject needs a unique id defined at the top of
      //! this file.
      virtual int UniqueParObjectId() const
      {
        return So_sh8PlastType::Instance().UniqueParObjectId();
      }

      //! Pack this class so it can be communicated
      //! Pack and \ref Unpack are used to communicate this element
      virtual void Pack(DRT::PackBuffer& data) const;

      //! Unpack data from a char vector into this class
      //! Pack and \ref Unpack are used to communicate this element
      virtual void Unpack(const std::vector<char>& data);

      //! Print this element
      void Print(std::ostream& os) const;

      //! return elementtype
      virtual So_sh8PlastType& ElementType() const { return So_sh8PlastType::Instance(); }

      //! read input for this element
      bool ReadElement(const std::string& eletype, const std::string& distype,
          DRT::INPUT::LineDefinition* linedef);

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
      So_sh8Plast& operator=(const So_sh8Plast& old);

     protected:
      static const int num_sp = 8;   ///< number of ANS sampling points, here 8
      static const int num_ans = 3;  ///< number of modified ANS strains (E_rt,E_st,E_tt), here 3
      //! shell-thickness direction
      ThicknessDirection thickdir_;

      ANSType anstype_;

      //! in case of changed "thin" direction this is true
      bool nodes_rearranged_;

      //! vector in thickness direction for compatibility with sosh8
      std::vector<double> thickvec_;

      static std::pair<bool, LINALG::Matrix<nsd_, nsd_>> jac_refe_;
      static std::pair<bool, LINALG::Matrix<nsd_, nsd_>> jac_curr_;
      static std::pair<bool, LINALG::Matrix<num_ans * num_sp, numdofperelement_>> B_ans_loc_;
      static std::pair<bool, LINALG::Matrix<numstr_, numstr_>> TinvT_;

      //! Calculate nonlinear stiffness and mass matrix with condensed plastic matrices
      virtual void nln_stiffmass(std::vector<double>& disp,  // current displacements
          std::vector<double>& vel,                          // current velocities
          std::vector<double>& temp,                         // current temperatures
          LINALG::Matrix<numdofperelement_, numdofperelement_>*
              stiffmatrix,  // element stiffness matrix
          LINALG::Matrix<numdofperelement_, numdofperelement_>* massmatrix,  // element mass matrix
          LINALG::Matrix<numdofperelement_, 1>* force,      // element internal force vector
          LINALG::Matrix<numgpt_post, numstr_>* elestress,  // stresses at GP
          LINALG::Matrix<numgpt_post, numstr_>* elestrain,  // strains at GP
          Teuchos::ParameterList& params,                   // algorithmic parameters e.g. time
          const INPAR::STR::StressType iostress,            // stress output option
          const INPAR::STR::StrainType iostrain             // strain output option
      );

      //! calculate nonlinear B-operator (potentially with ANS modification)
      virtual void CalculateBop(LINALG::Matrix<numstr_, numdofperelement_>* bop,
          const LINALG::Matrix<nsd_, nsd_>* defgrd, const LINALG::Matrix<nsd_, nen_>* N_XYZ,
          const int gp);

      //! Evaluate all ANS related data at the ANS sampling points
      void Anssetup(const LINALG::Matrix<nen_, nsd_>& xrefe,  ///< material element coords
          const LINALG::Matrix<nen_, nsd_>& xcurr,            ///< current element coords
          std::vector<LINALG::Matrix<nsd_, nen_>>**
              deriv_sp,  ///< derivs eval. at all sampling points
          std::vector<LINALG::Matrix<nsd_, nsd_>>& jac_sps,  ///< jac at all sampling points
          std::vector<LINALG::Matrix<nsd_, nsd_>>&
              jac_cur_sps,  ///< current jac at all sampling points
          LINALG::Matrix<num_ans * num_sp, numdofperelement_>& B_ans_loc);  ///< modified B

      //! Evaluate transformation matrix T (parameter->material) at gp
      void EvaluateT(const LINALG::Matrix<nsd_, nsd_>& jac,  ///< actual jacobian
          LINALG::Matrix<numstr_, numstr_>& TinvT);          ///< T^{-T}

      void AnsStrains(const int gp,
          std::vector<LINALG::Matrix<nsd_, nsd_>>& jac_sps,  // jac at all sampling points
          std::vector<LINALG::Matrix<nsd_, nsd_>>&
              jac_cur_sps  // current jac at all sampling points
      );

      //! Find "thin"=thickness direction
      ThicknessDirection findthickdir();

      //! Find parametric co-ordinate which directs in enforced thickness direction
      ThicknessDirection enfthickdir(LINALG::Matrix<nsd_, 1>&
              glothickdir  ///< global direction of enforced thickness direction
      );

      //! Re-initialize EAS data, needed for sosh8 morphing
      void ReInitEas(const DRT::ELEMENTS::So3Plast_EASType EASType);

      virtual void InvalidGpData()
      {
        So3_Plast<DRT::Element::hex8>::InvalidGpData();
        jac_refe_.first = false;
        jac_curr_.first = false;
        TinvT_.first = false;
      }

      void InvalidEleData()
      {
        So3_Plast<DRT::Element::hex8>::InvalidEleData();
        B_ans_loc_.first = false;
      }

      const LINALG::Matrix<nsd_, nsd_>& Jac_curr() const
      {
        if (jac_curr_.first == false) dserror("jac_curr_ not valid");
        return jac_curr_.second;
      }
      LINALG::Matrix<nsd_, nsd_>& SetJac_curr()
      {
        jac_curr_.first = true;
        return jac_curr_.second;
      }

      const LINALG::Matrix<nsd_, nsd_>& Jac_refe() const
      {
        if (jac_refe_.first == false) dserror("jac_refe_ not valid");
        return jac_refe_.second;
      }
      LINALG::Matrix<nsd_, nsd_>& SetJac_refe()
      {
        jac_refe_.first = true;
        return jac_refe_.second;
      }

      const LINALG::Matrix<num_ans * num_sp, numdofperelement_>& B_ans_loc() const
      {
        if (B_ans_loc_.first == false) dserror("B_ans_loc_ not valid");
        return B_ans_loc_.second;
      }
      LINALG::Matrix<num_ans * num_sp, numdofperelement_>& SetB_ans_loc()
      {
        B_ans_loc_.first = true;
        return B_ans_loc_.second;
      }

      const LINALG::Matrix<numstr_, numstr_>& TinvT() const
      {
        if (TinvT_.first == false) dserror("TinvT_ not valid");
        return TinvT_.second;
      }
      LINALG::Matrix<numstr_, numstr_>& SetTinvT()
      {
        TinvT_.first = true;
        return TinvT_.second;
      }
    };  // class So_sh8Plast
  }     // namespace ELEMENTS
}  // namespace DRT



#endif
