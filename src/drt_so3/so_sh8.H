/*----------------------------------------------------------------------*/
/*!
\brief solid shell8 element formulation

\maintainer Alexander Popp

\level 1
*/
/*----------------------------------------------------------------------*/
#ifndef SO_SH8_H
#define SO_SH8_H

#include "so_hex8.H"

// forward declarations
struct _SOH8_DATA;

namespace DRT
{
  // forward declarations
  class Discretization;

  namespace ELEMENTS
  {
    class So_sh8Type : public So_hex8Type
    {
     public:
      std::string Name() const { return "So_sh8Type"; }

      static So_sh8Type& Instance();

      virtual DRT::ParObject* Create(const std::vector<char>& data);

      virtual Teuchos::RCP<DRT::Element> Create(
          const std::string eletype, const std::string eledistype, const int id, const int owner);

      virtual Teuchos::RCP<DRT::Element> Create(const int id, const int owner);

      virtual int Initialize(DRT::Discretization& dis);

      virtual void NodalBlockInformation(
          DRT::Element* dwele, int& numdf, int& dimns, int& nv, int& np);

      virtual void ComputeNullSpace(DRT::Discretization& dis, std::vector<double>& ns,
          const double* x0, int numdf, int dimns);

      virtual void SetupElementDefinition(
          std::map<std::string, std::map<std::string, DRT::INPUT::LineDefinition>>& definitions);

     private:
      static So_sh8Type instance_;
    };

    /*!
    \brief A C++ 8-node Solid-Shell element inherited from so_hex8

    The Solid-Shell element technology is based on the work of
    (1) Vu-Quoc, Tan: "Optimal solid shells for non-linear analyses
                       of multilayer composites", CMAME 2003
    (2) Klinkel, Gruttmann, Wagner: "A robust non-linear solid shell element
                                     based on a mixed variational fromulation"

    Refer also to the Semesterarbeit of Alexander Popp, 2006

    \author frenzel (frenzel@lnm.mw.tum.de)
    */
    class So_sh8 : public So_hex8
    {
     public:
      //! @name Friends
      friend class So_sh8Type;
      friend class Soh8Surface;
      friend class Soh8Line;

      //@}
      //! @name Constructors and destructors and related methods

      /*!
      \brief Standard Constructor

      \param id : A unique global id
      \param owner : elements owning processor
      */
      So_sh8(int id, int owner);

      /*!
      \brief Copy Constructor

      Makes a deep copy of a Element

      */
      So_sh8(const So_sh8& old);

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
      virtual int UniqueParObjectId() const { return So_sh8Type::Instance().UniqueParObjectId(); }

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
      \brief Print this element
      */
      void Print(std::ostream& os) const;

      virtual So_sh8Type& ElementType() const { return So_sh8Type::Instance(); }

      //@}

      //! @name Input and Creation

      /*!
      \brief Read input for this element
      */
      bool ReadElement(const std::string& eletype, const std::string& distype,
          DRT::INPUT::LineDefinition* linedef);

      /*!
      \brief Destructor

      */
      virtual ~So_sh8();


      //@}

      //! @name Evaluation

      /*!
      \brief Evaluate an element

      Evaluate so_sh8 element stiffness, mass, internal forces, etc.

      \param params (in/out): ParameterList for communication between control routine
                              and elements
      \param discretization : pointer to discretization for de-assembly
      \param lm (in)        : location matrix for de-assembly
      \param elemat1 (out)  : (stiffness-)matrix to be filled by element. If NULL on input,
                              the controling method does not expect the element to fill
                              this matrix.
      \param elemat2 (out)  : (mass-)matrix to be filled by element. If NULL on input,
                              the controling method does not expect the element to fill
                              this matrix.
      \param elevec1 (out)  : (internal force-)vector to be filled by element. If NULL on input,
                              the controlling method does not expect the element
                              to fill this vector
      \param elevec2 (out)  : vector to be filled by element. If NULL on input,
                              the controlling method does not expect the element
                              to fill this vector
      \param elevec3 (out)  : vector to be filled by element. If NULL on input,
                              the controlling method does not expect the element
                              to fill this vector
      \return 0 if successful, negative otherwise
      */
      int Evaluate(Teuchos::ParameterList& params, DRT::Discretization& discretization,
          std::vector<int>& lm, Epetra_SerialDenseMatrix& elemat1,
          Epetra_SerialDenseMatrix& elemat2, Epetra_SerialDenseVector& elevec1,
          Epetra_SerialDenseVector& elevec2, Epetra_SerialDenseVector& elevec3);

      //! definition of shell-thickness direction
      enum ThicknessDirection
      {
        undefined,  ///< no clear direction identified
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
        none        ///< no rearrangement
      };

      enum ANSType
      {
        anssosh8,
        ansnone
      };

     private:
      // don't want = operator
      So_sh8& operator=(const So_sh8& old);

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

      //! Compute stiffness and mass matrix
      void sosh8_nlnstiffmass(std::vector<int>& lm,               ///< location matrix
          std::vector<double>& disp,                              ///< current displacements
          std::vector<double>& residual,                          ///< current residual displ
          LINALG::Matrix<NUMDOF_SOH8, NUMDOF_SOH8>* stiffmatrix,  ///< element stiffness matrix
          LINALG::Matrix<NUMDOF_SOH8, NUMDOF_SOH8>* massmatrix,   ///< element mass matrix
          LINALG::Matrix<NUMDOF_SOH8, 1>* force,                  ///< element internal force vector
          LINALG::Matrix<NUMDOF_SOH8, 1>* force_str,              // element structural force vector
          LINALG::Matrix<NUMGPT_SOH8, MAT::NUM_STRESS_3D>* elestress,  ///< stresses at GP
          LINALG::Matrix<NUMGPT_SOH8, MAT::NUM_STRESS_3D>* elestrain,  ///< strains at GP
          Teuchos::ParameterList& params,          ///< algorithmic parameters e.g. time
          const INPAR::STR::StressType iostress,   ///< stress output option
          const INPAR::STR::StrainType iostrain);  ///< strain output option

      //! Evaluate all ANS related data at the ANS sampling points
      void sosh8_anssetup(
          const LINALG::Matrix<NUMNOD_SOH8, NUMDIM_SOH8>& xrefe,  ///< material element coords
          const LINALG::Matrix<NUMNOD_SOH8, NUMDIM_SOH8>& xcurr,  ///< current element coords
          std::vector<LINALG::Matrix<NUMDIM_SOH8, NUMNOD_SOH8>>**
              deriv_sp,  ///< derivs eval. at all sampling points
          std::vector<LINALG::Matrix<NUMDIM_SOH8, NUMDIM_SOH8>>&
              jac_sps,  ///< jac at all sampling points
          std::vector<LINALG::Matrix<NUMDIM_SOH8, NUMDIM_SOH8>>&
              jac_cur_sps,  ///< current jac at all sampling points
          LINALG::Matrix<num_ans * num_sp, NUMDOF_SOH8>& B_ans_loc) const;  ///< modified B

      //! Evaluate transformation matrix T (parameter->material) at gp
      void sosh8_evaluateT(
          const LINALG::Matrix<NUMDIM_SOH8, NUMDIM_SOH8>& jac,             ///< actual jacobian
          LINALG::Matrix<MAT::NUM_STRESS_3D, MAT::NUM_STRESS_3D>& TinvT);  ///< T^{-T}

      //! Return true Cauchy-stress at gausspoint
      void sosh8_Cauchy(LINALG::Matrix<NUMGPT_SOH8, MAT::NUM_STRESS_3D>* elestress, const int gp,
          const LINALG::Matrix<NUMDIM_SOH8, NUMDIM_SOH8>& defgrd,
          const LINALG::Matrix<MAT::NUM_STRESS_3D, 1>& glstrain,
          const LINALG::Matrix<MAT::NUM_STRESS_3D, 1>& stress);

      //! Find "thin"=thickness direction
      ThicknessDirection sosh8_findthickdir();

      //! Find aspect ratio of the element
      double sosh8_calcaspectratio();

      //! Calculate the STC matrix
      virtual void CalcSTCMatrix(LINALG::Matrix<NUMDOF_SOH8, NUMDOF_SOH8>& elemat1,
          const INPAR::STR::STC_Scale stc_scaling, const int stc_layer, std::vector<int>& lm,
          DRT::Discretization& discretization, bool calcinverse);

      //! Find parametric co-ordinate which directs in enforced thickness direction
      ThicknessDirection sosh8_enfthickdir(LINALG::Matrix<NUMDIM_SOH8, 1>&
              glothickdir  ///< global direction of enforced thickness direction
      );

      //! return thickness direction
      const std::vector<double> GetThickvec() { return thickvec_; };

      //! Debug gmsh-plot to check thickness direction
      void sosh8_gmshplotlabeledelement(const int LabelIds[NUMNOD_SOH8]);

      const std::vector<LINALG::Matrix<NUMDIM_SOH8, NUMNOD_SOH8>> sosh8_derivs_sdc();

      /** \brief Evaluate the reference and current jacobian as well as the
       *  respective determinants */
      bool sosh8_evaluatejacobians(const unsigned gp,
          const std::vector<LINALG::Matrix<NUMDIM_SOH8, NUMNOD_SOH8>>& derivs,
          const LINALG::Matrix<NUMNOD_SOH8, NUMDIM_SOH8>& xrefe,
          const LINALG::Matrix<NUMNOD_SOH8, NUMDIM_SOH8>& xcurr,
          LINALG::Matrix<NUMDIM_SOH8, NUMDIM_SOH8>& jac_ref, double& detJ_ref,
          LINALG::Matrix<NUMDIM_SOH8, NUMDIM_SOH8>& jac_curr, double& detJ_curr) const;

      /** \brief evaluate the jacobian and the determinant for the given GP */
      void sosh8_evaluatejacobian(const unsigned gp,
          const std::vector<LINALG::Matrix<NUMDIM_SOH8, NUMNOD_SOH8>>& derivs,
          const LINALG::Matrix<NUMNOD_SOH8, NUMDIM_SOH8>& x,
          LINALG::Matrix<NUMDIM_SOH8, NUMDIM_SOH8>& jac, double& detJ) const;

      /** \brief Get the local B-operator */
      void sosh8_get_bop_loc(const unsigned gp,
          const std::vector<LINALG::Matrix<NUMDIM_SOH8, NUMNOD_SOH8>>& derivs,
          const LINALG::Matrix<NUMDIM_SOH8, NUMDIM_SOH8>& jac_curr, const double* r,
          const double* s, const LINALG::Matrix<num_ans * num_sp, NUMDOF_SOH8>& B_ans_loc,
          LINALG::Matrix<MAT::NUM_STRESS_3D, NUMDOF_SOH8>& bop_loc) const;

      /** \brief Get the local green lagrange strain */
      void sosh8_get_glstrain_loc(const unsigned gp,
          const LINALG::Matrix<NUMDIM_SOH8, NUMDIM_SOH8>& jac_curr,
          const LINALG::Matrix<NUMDIM_SOH8, NUMDIM_SOH8>& jac,
          const std::vector<LINALG::Matrix<NUMDIM_SOH8, NUMDIM_SOH8>>& jac_sps,
          const std::vector<LINALG::Matrix<NUMDIM_SOH8, NUMDIM_SOH8>>& jac_cur_sps, const double* r,
          const double* s, LINALG::Matrix<MAT::NUM_STRESS_3D, 1>& lstrain) const;

      void sosh8_get_deformationgradient(const unsigned gp,
          const std::vector<LINALG::Matrix<NUMDIM_SOH8, NUMNOD_SOH8>>& derivs,
          const LINALG::Matrix<NUMNOD_SOH8, NUMDIM_SOH8>& xcurr,
          const LINALG::Matrix<MAT::NUM_STRESS_3D, 1>& glstrain,
          LINALG::Matrix<NUMDIM_SOH8, NUMDIM_SOH8>& defgrd) const;

      double sosh8_calc_energy(const std::vector<double>& disp, Teuchos::ParameterList& params);

      double sosh8_third_invariant(const LINALG::Matrix<MAT::NUM_STRESS_3D, 1>& glstrain) const;
    };  // class So_sh8



  }  // namespace ELEMENTS
}  // namespace DRT


#endif  // #ifndef SO_SH8_H
