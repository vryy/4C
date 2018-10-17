/*!----------------------------------------------------------------------
\file so_sh18.H
\brief a in-plane bi-quadratic solid-like 7-parameter shell
\level 3

<pre>
\maintainer Christoph Meier
</pre>

*----------------------------------------------------------------------*/
#ifndef SH18_H
#define SH18_H

#include "../drt_lib/drt_element.H"
#include "../drt_lib/drt_elementtype.H"
#include "../drt_inpar/inpar_structure.H"
#include "../drt_mat/material.H"
#include "so_hex18.H"

const int num_eas = 9;

namespace DRT
{
  // forward declarations
  class Discretization;

  namespace ELEMENTS
  {
    class So_sh18Type : public So_hex18Type
    {
     public:
      //! @name Friends
      friend class So_sh18PlastType;

      std::string Name() const { return "So_sh18Type"; }

      static So_sh18Type& Instance();

      virtual DRT::ParObject* Create(const std::vector<char>& data);

      virtual Teuchos::RCP<DRT::Element> Create(
          const std::string eletype, const std::string eledistype, const int id, const int owner);

      virtual Teuchos::RCP<DRT::Element> Create(const int id, const int owner);

      virtual int Initialize(DRT::Discretization& dis);

      virtual void SetupElementDefinition(
          std::map<std::string, std::map<std::string, DRT::INPUT::LineDefinition>>& definitions);

     private:
      static So_sh18Type instance_;
    };

    /*!
    \brief A C++ version of the 18-node hex solid-shell element

    A structural 18-node hexahedral solid element for large deformations.
    Quadratic interpolation within a plane and linear interpolation in the third direction
    */
    class So_sh18 : public virtual So_hex18
    {
     public:
      //! @name Friends
      friend class So_sh18Type;
      friend class So_sh18PlastType;

      //@}
      //! @name Constructors and destructors and related methods

      /*!
      \brief Standard Constructor

      \param id : A unique global id
      \param owner : elements owner
      */
      So_sh18(int id, int owner);

      /*!
      \brief Copy Constructor

      Makes a deep copy of a Element

      */
      So_sh18(const So_sh18& old);

      /*!
      \brief Deep copy this instance of Solid3 and return pointer to the copy

      The Clone() method is used from the virtual base class Element in cases
      where the type of the derived class is unknown and a copy-ctor is needed

      */
      DRT::Element* Clone() const;


      /*!
      \brief Return unique ParObject id

      every class implementing ParObject needs a unique id defined at the
      top of this file.
      */
      virtual int UniqueParObjectId() const { return So_sh18Type::Instance().UniqueParObjectId(); }

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
      virtual ~So_sh18();

      //@}

      //! @name Acess methods

      /*!
       \brief Does this element use EAS?
       */
      virtual bool HaveEAS() const { return eas_; };

      /*!
      \brief Print this element
      */
      void Print(std::ostream& os) const;



      virtual DRT::ElementType& ElementType() const { return So_sh18Type::Instance(); }

      //@}

      //! @name Input and Creation

      //@}

      //! @name Input and Creation

      /*!
      \brief Read input for this element
      */
      virtual bool ReadElement(const std::string& eletype, const std::string& distype,
          DRT::INPUT::LineDefinition* linedef);

      //@}


     protected:
      //! use DSG for transverse shear locking
      bool dsg_shear_;
      //! use DSG for membrane locking
      bool dsg_membrane_;
      //! use DSG for curvature thickness locking
      bool dsg_ctl_;
      //! use EAS
      bool eas_;

      // EAS stuff *****************************************************
      LINALG::Matrix<num_eas, num_eas> KaaInv_;
      LINALG::Matrix<num_eas, 1> feas_;
      LINALG::Matrix<num_eas, NUMDOF_SOH18> Kad_;
      LINALG::Matrix<num_eas, 1> alpha_eas_;
      LINALG::Matrix<num_eas, 1> alpha_eas_last_timestep_;
      LINALG::Matrix<num_eas, 1> alpha_eas_delta_over_last_timestep_;
      LINALG::Matrix<num_eas, 1> alpha_eas_inc_;
      double old_step_length_;
      // EAS stuff *****************************************************

      // DSG factors ***************************************************
      // every DSG-modification is equivalent to a certain
      // linear combination of nodal coordinates / displacements
      std::vector<LINALG::Matrix<9, 9>> dsg_shear_r_;
      std::vector<LINALG::Matrix<9, 9>> dsg_shear_s_;
      std::vector<LINALG::Matrix<9, 9>> dsg_membrane_r_;
      std::vector<LINALG::Matrix<9, 9>> dsg_membrane_s_;
      std::vector<LINALG::Matrix<9, 9>> dsg_membrane_rs_;
      std::vector<LINALG::Matrix<9, 9>> dsg_transverse_t_;
      // DSG factors ***************************************************


      // internal calculation methods

      //! don't want = operator
      So_sh18& operator=(const So_sh18& old);


      //! init the inverse of the jacobian and its determinant in the material configuration
      //! return the number of negative eigenvalues
      virtual int InitJacobianMapping();

      //! Calculate nonlinear stiffness and mass matrix
      virtual void nlnstiffmass(std::vector<int>& lm,               ///< location matrix
          std::vector<double>& disp,                                ///< current displacements
          std::vector<double>& residual,                            ///< current residual displ
          LINALG::Matrix<NUMDOF_SOH18, NUMDOF_SOH18>* stiffmatrix,  ///< element stiffness matrix
          LINALG::Matrix<NUMDOF_SOH18, NUMDOF_SOH18>* massmatrix,   ///< element mass matrix
          LINALG::Matrix<NUMDOF_SOH18, 1>* force,  ///< element internal force vector
          LINALG::Matrix<NUMGPT_SOH18, MAT::NUM_STRESS_3D>* elestress,  ///< stresses at GP
          LINALG::Matrix<NUMGPT_SOH18, MAT::NUM_STRESS_3D>* elestrain,  ///< strains at GP
          Teuchos::ParameterList& params,         ///< algorithmic parameters e.g. time
          const INPAR::STR::StressType iostress,  ///< stress output option
          const INPAR::STR::StrainType iostrain   ///< strain output option
      );

      // parameter space coords of one node
      LINALG::Matrix<3, 1> NodeParamCoord(const int node);
      // parameter space coords of all nodes
      LINALG::Matrix<18, 3> NodeParamCoord();

      void FlipT();

      //! Lump mass matrix
      void soh18_lumpmass(LINALG::Matrix<NUMDOF_SOH18, NUMDOF_SOH18>* emass);

      void EvaluateT(const LINALG::Matrix<NUMDIM_SOH18, NUMDIM_SOH18>& jac,
          LINALG::Matrix<MAT::NUM_STRESS_3D, MAT::NUM_STRESS_3D>& TinvT);

      void EasSetup(std::vector<LINALG::Matrix<6, num_eas>>& M_gp,  // M-matrix evaluated at GPs
          LINALG::Matrix<3, 1>& G3_contra,  // contravariant basis vector G3 at element center
          const LINALG::Matrix<NUMNOD_SOH18, 3> xrefe);  // material element coords

      void SetupDSG();
      void Integrate_dsg_shear_r(const int gp, LINALG::Matrix<9, 9>& dsg_shear_r);
      void Integrate_dsg_shear_s(const int gp, LINALG::Matrix<9, 9>& dsg_shear_s);
      void Integrate_dsg_membrane_r(const int gp, LINALG::Matrix<9, 9>& dsg_membrane_r);
      void Integrate_dsg_membrane_s(const int gp, LINALG::Matrix<9, 9>& dsg_membrane_r);
      void Integrate_dsg_membrane_rs(const int gp, LINALG::Matrix<9, 9>& dsg_membrane_rs);
      void Integrate_dsg_transverse_t(const int gp, LINALG::Matrix<9, 9>& dsg_transverse_t);

      void CalcConsistentDefgrd(LINALG::Matrix<3, 3> defgrd_disp, LINALG::Matrix<6, 1> glstrain_mod,
          LINALG::Matrix<3, 3>& defgrd_mod);

      void CalculateBopLoc(const LINALG::Matrix<NUMNOD_SOH18, NUMDIM_SOH18>& xcurr,
          const LINALG::Matrix<NUMNOD_SOH18, NUMDIM_SOH18>& xrefe,
          const LINALG::Matrix<9, 1>& shape_q9, const LINALG::Matrix<2, 9>& deriv_q9, const int gp,
          LINALG::Matrix<MAT::NUM_STRESS_3D, NUMDOF_SOH18>& bop_loc);

      void CalculateLocStrain(const LINALG::Matrix<NUMNOD_SOH18, NUMDIM_SOH18>& xcurr,
          const LINALG::Matrix<NUMNOD_SOH18, NUMDIM_SOH18>& xrefe,
          const LINALG::Matrix<9, 1>& shape_q9, const LINALG::Matrix<2, 9>& deriv_q9, const int gp,
          LINALG::Matrix<MAT::NUM_STRESS_3D, 1>& lstrain);

      void CalculateGeoStiff(const LINALG::Matrix<9, 1>& shape_q9,
          const LINALG::Matrix<2, 9>& deriv_q9,
          LINALG::Matrix<MAT::NUM_STRESS_3D, MAT::NUM_STRESS_3D>& TinvT, const int gp,
          const double detJ_w, const LINALG::Matrix<MAT::NUM_STRESS_3D, 1>& stress,
          LINALG::Matrix<NUMDOF_SOH18, NUMDOF_SOH18>* stiffmatrix);

      void Update();


      /** recover elementwise stored stuff */
      virtual void Recover(const std::vector<double>& residual);

      //@}

    };  // class So_sh18



    //=======================================================================
    //=======================================================================
    //=======================================================================
    //=======================================================================



  }  // namespace ELEMENTS
}  // namespace DRT

#endif  // #ifndef SH8_H
