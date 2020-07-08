/*----------------------------------------------------------------------*/
/*! \file
\brief
\level 1


*----------------------------------------------------------------------*/

#include "so_weg6.H"

namespace DRT
{
  // forward declarations
  class Discretization;

  namespace ELEMENTS
  {
    class So_shw6Type : public DRT::ElementType
    {
     public:
      std::string Name() const { return "So_shw6Type"; }

      static So_shw6Type& Instance();

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
      static So_shw6Type instance_;
    };

    /*!
    \brief A C++ 6-node wedge Solid-Shell element inherited from so_weg6

    This element adopts the ANS approach to improve shear and curvature-thickness
    (trapezoidal) locking. However, being triangular inplane, the ANS evaluation
    points of the shear strains are choosen just on two of three surfaces (See
    also DRT::ELEMENTS::So_shw6::soshw6_anssetup).
    Thus, the element becomes mesh dependent (where in material space lies my
    parameter space), but this is in favour of less pronounced shear locking
    compared to a invariant formulation (locking free just for certain mesh configurations).

    See Diss. Frank Koschnick for details.

    */
    class So_shw6 : public So_weg6
    {
     public:
      //! @name Friends
      friend class So_shw6Type;
      friend class Soweg6Surface;
      friend class Soweg6Line;

      //@}
      //! @name Constructors and destructors and related methods

      /*!
      \brief Standard Constructor

      \param id : A unique global id
      \param owner : elements owning processor
      */
      So_shw6(int id, int owner);

      /*!
      \brief Copy Constructor

      Makes a deep copy of a Element

      */
      So_shw6(const So_shw6& old);

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
      virtual int UniqueParObjectId() const { return So_shw6Type::Instance().UniqueParObjectId(); }

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
      virtual ~So_shw6();

      //@}

      //! @name Acess methods

      /*!
       \brief Does this element use EAS?
       */
      virtual bool HaveEAS() const { return (eastype_ != soshw6_easnone); };


      /*!
      \brief Print this element
      */
      void Print(std::ostream& os) const;

      virtual DRT::ElementType& ElementType() const { return So_shw6Type::Instance(); }

      //@}

      //! @name Input and Creation

      /*!
      \brief Read input for this element
      */
      bool ReadElement(const std::string& eletype, const std::string& distype,
          DRT::INPUT::LineDefinition* linedef);

      //@}

      //! @name Evaluation



     private:
      /*! \brief EAS technology enhancement types of so_shw6
       * Solid-Shell Wedge6 has EAS enhancement of GL-strains to avoid locking
       * \param soshw6_poissonthick: 1 parameter for correct kinematic of solid-shell,
       *                             allows linear strain in thick-dir
       * \param soshw6_none:         no EAS
       */
      enum EASType
      {
        soshw6_easpoisthick = 1,
        soshw6_easnone = 0
      };

      //! type of EAS technology
      EASType eastype_;

      //! number of EAS parameters (alphas), defined by 'EASType"
      int neas_;
      // line search parameter (old step length)
      double old_step_length_;

      //! reorder nodes to optimally align material space with parameter space
      bool optimal_parameterspace_map_;

      //! in case of changed nodal order this is true
      bool nodes_rearranged_;

      //! don't want = operator
      So_shw6& operator=(const So_shw6& old);

      /*!
      \brief Evaluate an element

      Evaluate so_hex8 element stiffness, mass, internal forces, etc.

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
      virtual int Evaluate(Teuchos::ParameterList& params, DRT::Discretization& discretization,
          std::vector<int>& lm, Epetra_SerialDenseMatrix& elemat1,
          Epetra_SerialDenseMatrix& elemat2, Epetra_SerialDenseVector& elevec1,
          Epetra_SerialDenseVector& elevec2, Epetra_SerialDenseVector& elevec3);

      //! Compute stiffness and mass matrix
      virtual void soshw6_nlnstiffmass(std::vector<int>& lm,      // location matrix
          std::vector<double>& disp,                              // current displacements
          std::vector<double>& residual,                          // current residual displ
          LINALG::Matrix<NUMDOF_WEG6, NUMDOF_WEG6>* stiffmatrix,  // element stiffness matrix
          LINALG::Matrix<NUMDOF_WEG6, NUMDOF_WEG6>* massmatrix,   // element mass matrix
          LINALG::Matrix<NUMDOF_WEG6, 1>* force,                  // element internal force vector
          LINALG::Matrix<NUMDOF_WEG6, 1>* force_str,              // structure force
          LINALG::Matrix<NUMGPT_WEG6, MAT::NUM_STRESS_3D>* elestress,  // stresses at GP
          LINALG::Matrix<NUMGPT_WEG6, MAT::NUM_STRESS_3D>* elestrain,  // strains at GP
          Teuchos::ParameterList& params,          // algorithmic parameters e.g. time
          const INPAR::STR::StressType iostress,   // stress output option
          const INPAR::STR::StrainType iostrain);  // strain output option

      static const int num_sp = 5;   ///< number of ANS sampling points
      static const int num_ans = 3;  ///< number of modified ANS strains (E_rt,E_st,E_tt)

      //! Setup ANS interpolation (shear and trapezoidal)
      void soshw6_anssetup(
          const LINALG::Matrix<NUMNOD_WEG6, NUMDIM_WEG6>& xrefe,  ///< material element coords
          const LINALG::Matrix<NUMNOD_WEG6, NUMDIM_WEG6>& xcurr,  ///< current element coords
          std::vector<LINALG::Matrix<NUMDIM_WEG6, NUMNOD_WEG6>>**
              deriv_sp,  ///< derivs eval. at all sampling points
          std::vector<LINALG::Matrix<NUMDIM_WEG6, NUMDIM_WEG6>>&
              jac_sps,  ///< jac at all sampling points
          std::vector<LINALG::Matrix<NUMDIM_WEG6, NUMDIM_WEG6>>&
              jac_cur_sps,  ///< current jac at all sampling points
          LINALG::Matrix<num_ans * num_sp, NUMDOF_WEG6>& B_ans_loc);  ///< modified B

      //! Transformation matrix parameter->material space
      void soshw6_evaluateT(const LINALG::Matrix<NUMDIM_WEG6, NUMDIM_WEG6>& jac,
          LINALG::Matrix<MAT::NUM_STRESS_3D, MAT::NUM_STRESS_3D>& TinvT);

      //! EAS technology, init
      void soshw6_easinit();

      //! EAS technology, setup necessary data
      void soshw6_eassetup(
          std::vector<Epetra_SerialDenseMatrix>** M_GP,  // M-matrix evaluated at GPs
          double& detJ0,                                 // det of Jacobian at origin
          LINALG::Matrix<MAT::NUM_STRESS_3D, MAT::NUM_STRESS_3D>&
              T0invT,                                              // maps M(origin) local to global
          const LINALG::Matrix<NUMNOD_WEG6, NUMDIM_WEG6>& xrefe);  // material element coords

      //! Evaluate Cauchy stress (mapping with disp_based F per default)
      void soshw6_Cauchy(LINALG::Matrix<NUMGPT_WEG6, MAT::NUM_STRESS_3D>* elestress, const int gp,
          const LINALG::Matrix<NUMDIM_WEG6, NUMDIM_WEG6>& defgrd,
          const LINALG::Matrix<MAT::NUM_STRESS_3D, 1>& glstrain,
          const LINALG::Matrix<MAT::NUM_STRESS_3D, 1>& stress);

      //! Reorient wedge to optimally align material with parameter space
      int soshw6_findoptparmap();

      //! recover condensed EAS variables
      void soshw6_recover(const std::vector<double>& residual);

    };  // class So_shw6


  }  // namespace ELEMENTS
}  // namespace DRT
