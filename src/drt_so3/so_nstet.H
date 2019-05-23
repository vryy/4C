/*!----------------------------------------------------------------------

\brief NStet element

\level 3

\maintainer Christoph Meier
*----------------------------------------------------------------------*/
#ifndef SO_NSTET4_H
#define SO_NSTET4_H

#include "../drt_lib/drt_element.H"
#include "../drt_lib/drt_elementtype.H"
#include "../drt_lib/drt_parobjectfactory.H"
#include "../drt_inpar/inpar_structure.H"

#include "../drt_mat/material.H"

#include "Epetra_Time.h"

//#define PRINT_NSTET            ///< printout of setup
#define NSTET_NUMGREEDY 2    ///< number of sweeps of greedy patching
#define MIS_MIN_PATCHSIZE 3  ///< minimum MIS patch size required

// stabilization approach & parameters
//#define PUSOSOLBERG              ///< run the Puso&Solberg tet
#define ALPHA_NSTET 0.01  ///< stabilization parameter for vol-dev split stabilization & Puso-tet
#define BETA_NSTET 0.9    ///< stabilization parameter for MIS-pressure stabilization

// make parameters match for puso&solberg tet
#ifdef PUSOSOLBERG
#undef BETA_NSTET
#define BETA_NSTET 0.0  ///< beta has to be zero in puso&solberg tet
#undef MIS_MIN_PATCHSIZE
#define MIS_MIN_PATCHSIZE 1
#undef NSTET_NUMGREEDY
#define NSTET_NUMGREEDY 5
#endif

namespace DRT
{
  // forward declarations
  class Discretization;

  namespace ELEMENTS
  {
    class StructuralVolume;
    class NStet;


    //=======================================================================
    //=======================================================================

    class NStetType : public DRT::ElementType
    {
      //! allow NStet element to access the nodal data
      friend class DRT::ELEMENTS::NStet;

     public:
      std::string Name() const { return "NStetType"; }

      static NStetType& Instance();

      virtual DRT::ParObject* Create(const std::vector<char>& data);

      virtual Teuchos::RCP<DRT::Element> Create(
          const std::string eletype, const std::string eledistype, const int id, const int owner);

      virtual Teuchos::RCP<DRT::Element> Create(const int id, const int owner);

      virtual int Initialize(DRT::Discretization& dis);

      virtual void PreEvaluate(DRT::Discretization& dis, Teuchos::ParameterList& p,
          Teuchos::RCP<LINALG::SparseOperator> systemmatrix1,
          Teuchos::RCP<LINALG::SparseOperator> systemmatrix2,
          Teuchos::RCP<Epetra_Vector> systemvector1, Teuchos::RCP<Epetra_Vector> systemvector2,
          Teuchos::RCP<Epetra_Vector> systemvector3);

      virtual void NodalBlockInformation(
          DRT::Element* dwele, int& numdf, int& dimns, int& nv, int& np);

      virtual void ComputeNullSpace(DRT::Discretization& dis, std::vector<double>& ns,
          const double* x0, int numdf, int dimns);

      virtual void SetupElementDefinition(
          std::map<std::string, std::map<std::string, DRT::INPUT::LineDefinition>>& definitions);

     private:
      static NStetType instance_;

      // map of row nodes adjacent to NStet elements
      std::map<int, DRT::Node*> noderids_;

      // map of column NStet elements
      std::map<int, DRT::ELEMENTS::NStet*> elecids_;
      Teuchos::RCP<Epetra_Map> elecmap_;
      Teuchos::RCP<Epetra_Map> elermap_;

      // nodal stresses and strains for output
      Teuchos::RCP<Epetra_MultiVector> nstress_;
      Teuchos::RCP<Epetra_MultiVector> nstrain_;

      // vector of elements adjacent to each row node
      std::map<int, std::vector<DRT::ELEMENTS::NStet*>> adjele_;

      // map of nodes forming a patch around a node
      std::map<int, std::map<int, DRT::Node*>> adjnode_;

      // map of location vectors for patch around a node
      std::map<int, std::vector<int>> adjlm_;

      // vector of elements associated with an MIS node
      // key: mis node number
      std::map<int, std::vector<DRT::ELEMENTS::NStet*>> pstab_adjele_;

      // map indicating the weight of the mis node on an element
      // key: mis node number
      // values: weights for all elements on patch in pstab_adjele_
      std::map<int, std::vector<double>> pstab_adjele_weight_;

      // adjnode and adjlm for mis patches
      std::map<int, std::map<int, DRT::Node*>> pstab_adjnode_;
      std::map<int, std::vector<int>> pstab_adjlm_;

      // map indicating which element is associated with which MIS noderids
      // key is column element gid
      // value is MIS nodes this element is associated with
      std::map<int, std::vector<int>> pstab_cid_mis_;
      std::map<int, std::vector<double>> pstab_cid_mis_weight_;

      Teuchos::RCP<Epetra_Map> pstab_misstressout_;

      // nodal MIS stresses and strains for output
      Teuchos::RCP<Epetra_MultiVector> pstab_nstress_;
      Teuchos::RCP<Epetra_MultiVector> pstab_nstrain_;

      void InitElementsandMaps(std::map<int, DRT::ELEMENTS::NStet*>& elecids,
          std::map<int, DRT::Node*>& noderids, const int myrank, const int numproc,
          DRT::Discretization& dis);

      void InitAdjacency(std::map<int, DRT::ELEMENTS::NStet*>& elecids,
          std::map<int, DRT::Node*>& noderids,
          std::map<int, std::vector<DRT::ELEMENTS::NStet*>>& adjele,
          std::map<int, std::map<int, DRT::Node*>>& adjnode, std::map<int, std::vector<int>>& adjlm,
          DRT::Discretization& dis);

      void InitMISnode(std::map<int, int>& misnodesmap, std::map<int, DRT::Node*>& rnodes,
          const int myrank, const int numproc, DRT::Discretization& dis);

      void InitMISpatchesGreedyI(std::map<int, int>& misnodesmap,
          std::map<int, std::vector<DRT::ELEMENTS::NStet*>>& pstab_adjele,
          std::map<int, std::vector<int>>& pstab_cid_mis,
          std::map<int, std::vector<double>>& pstab_cid_mis_weight,
          std::map<int, DRT::ELEMENTS::NStet*>& elecids,
          std::map<int, DRT::ELEMENTS::NStet*>& elecids_full, std::map<int, DRT::Node*>& noderids,
          const int myrank, const int numproc, DRT::Discretization& dis);

      void InitMISpatchesGreedyII(std::map<int, int>& misnodesmap,
          std::map<int, std::vector<DRT::ELEMENTS::NStet*>>& pstab_adjele,
          std::map<int, std::vector<double>>& pstab_adjele_weight,
          std::map<int, std::vector<int>>& pstab_cid_mis,
          std::map<int, std::vector<double>>& pstab_cid_mis_weight,
          std::map<int, DRT::ELEMENTS::NStet*>& elecids,
          std::map<int, DRT::ELEMENTS::NStet*>& elecids_full, std::map<int, DRT::Node*>& noderids,
          const int myrank, const int numproc, DRT::Discretization& dis);

      void InitMISpatchesGreedyIII(std::map<int, int>& misnodesmap,
          std::map<int, std::vector<DRT::ELEMENTS::NStet*>>& pstab_adjele,
          std::map<int, std::vector<double>>& pstab_adjele_weight,
          std::map<int, std::vector<int>>& pstab_cid_mis,
          std::map<int, std::vector<double>>& pstab_cid_mis_weight,
          std::map<int, DRT::ELEMENTS::NStet*>& elecids,
          std::map<int, DRT::ELEMENTS::NStet*>& elecids_full, std::map<int, DRT::Node*>& noderids,
          const int myrank, const int numproc, DRT::Discretization& dis);


      void InitMISAdjacency(std::map<int, std::vector<DRT::ELEMENTS::NStet*>>& pstab_adjele,
          std::map<int, std::vector<DRT::ELEMENTS::NStet*>>& adjele,
          std::map<int, std::map<int, DRT::Node*>>& pstab_adjnode,
          std::map<int, std::map<int, DRT::Node*>>& adjnode,
          std::map<int, std::vector<int>>& pstab_adjlm, std::map<int, std::vector<int>>& adjlm,
          const int myrank, const int numproc, DRT::Discretization& dis);

      Teuchos::RCP<Epetra_Map> InitMISStressMap(
          std::map<int, std::vector<int>>& pstab_cid_mis, DRT::Discretization& dis);


      void ElementDeformationGradient(DRT::Discretization& dis);

      void NodalIntegration(Epetra_SerialDenseMatrix* stiff, Epetra_SerialDenseVector* force,
          std::map<int, DRT::Node*>& adjnode, std::vector<DRT::ELEMENTS::NStet*>& adjele,
          std::vector<int>& lm, const Epetra_Vector& disp, DRT::Discretization& dis,
          std::vector<double>* nodalstress, std::vector<double>* nodalstrain,
          const INPAR::STR::StressType iostress, const INPAR::STR::StrainType iostrain);


      void MISNodalIntegration(Epetra_SerialDenseMatrix* stiff, Epetra_SerialDenseVector* force,
          std::map<int, DRT::Node*>& adjnode, std::vector<DRT::ELEMENTS::NStet*>& adjele,
          std::vector<double>& weight, std::vector<int>& lm, const Epetra_Vector& disp,
          DRT::Discretization& dis, std::vector<double>* nodalstress,
          std::vector<double>* nodalstrain, const INPAR::STR::StressType iostress,
          const INPAR::STR::StrainType iostrain);

      void SelectMaterial(Teuchos::RCP<MAT::Material> mat, LINALG::Matrix<6, 1>& stress,
          LINALG::Matrix<6, 6>& cmat, double& density, LINALG::Matrix<6, 1>& glstrain,
          LINALG::Matrix<3, 3>& defgrd, int gp, const int eleGID);

      // compute deviatoric stresses and tangent
      static void DevStressTangent(LINALG::Matrix<6, 1>& Sdev, LINALG::Matrix<6, 6>& CCdev,
          LINALG::Matrix<6, 6>& CC, const LINALG::Matrix<6, 1>& S, const LINALG::Matrix<3, 3>& C);

      void StrainOutput(const INPAR::STR::StrainType iostrain, std::vector<double>& nodalstrain,
          LINALG::Matrix<3, 3>& F, const double& detF, const double volweight,
          const double devweight);

      void StrainOutput(const INPAR::STR::StrainType iostrain, std::vector<double>& nodalstrain,
          LINALG::Matrix<3, 3>& F, LINALG::Matrix<6, 1>& glstrain, const double weight);

      void StressOutput(const INPAR::STR::StressType iostress, std::vector<double>& nodalstress,
          LINALG::Matrix<6, 1>& stress, LINALG::Matrix<3, 3>& F, const double& detF);
    };

    /*!
    \brief A nodal-averaged strain 4-noded tet element

    */
    class NStet : public DRT::Element
    {
     public:
      friend class NStetType;

      /*!
      \brief Standard Constructor

      \param id : A unique global id
      */
      NStet(int id, int owner);

      /*!
      \brief Copy Constructor

      Makes a deep copy of a Element

      */
      NStet(const NStet& old);

      /*!
      \brief Deep copy this instance of Solid3 and return pointer to the copy

      The Clone() method is used from the virtual base class Element in cases
      where the type of the derived class is unknown and a copy-ctor is needed

      */
      inline DRT::Element* Clone() const { return new DRT::ELEMENTS::NStet(*this); }

      /*!
      \brief Get shape type of element
      */
      inline virtual DiscretizationType Shape() const { return tet4; }

      /*!
      \brief Return number of volumes of this element
      */
      virtual int NumVolume() const { return 1; }

      /*!
      \brief Return number of surfaces of this element
      */
      virtual int NumSurface() const { return 4; }

      /*!
      \brief Return number of lines of this element
      */
      virtual int NumLine() const { return 6; }

      /*!
      \brief Get vector of Teuchos::RCPs to the lines of this element

      */
      virtual std::vector<Teuchos::RCP<DRT::Element>> Lines();

      /*!
      \brief Get vector of Teuchos::RCPs to the surfaces of this element

      */
      virtual std::vector<Teuchos::RCP<DRT::Element>> Surfaces();

      /*!
      \brief Get vector of Teuchos::RCPs to the volumes of this element

      */
      virtual std::vector<Teuchos::RCP<DRT::Element>> Volumes();

      /*!
      \brief Return unique ParObject id

      every class implementing ParObject needs a unique id defined at the
      top of this file.
      */
      virtual int UniqueParObjectId() const { return NStetType::Instance().UniqueParObjectId(); }

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
      virtual ~NStet();

      /*!
      \brief Get number of degrees of freedom of a certain node
             (implements pure virtual DRT::Element)

      The element decides how many degrees of freedom its nodes must have.
      As this may vary along a simulation, the element can redecide the
      number of degrees of freedom per node along the way for each of it's nodes
      separately.
      */
      virtual inline int NumDofPerNode(const DRT::Node& node) const { return 3; }

      /*!
      \brief Get number of degrees of freedom per element
             (implements pure virtual DRT::Element)

      The element decides how many element degrees of freedom it has.
      It can redecide along the way of a simulation.

      \note Element degrees of freedom mentioned here are dofs that are visible
            at the level of the total system of equations. Purely internal
            element dofs that are condensed internally should NOT be considered.
      */
      virtual inline int NumDofPerElement() const { return 0; }

      /*!
      \brief Print this element
      */
      void Print(std::ostream& os) const;

      virtual NStetType& ElementType() const { return NStetType::Instance(); }

      /*!
      \brief Read input for this element
      */
      bool ReadElement(const std::string& eletype, const std::string& distype,
          DRT::INPUT::LineDefinition* linedef);

      /*!
      \brief Evaluate an element

      Evaluate so_tet4 element stiffness, mass, internal forces, etc.

      \param params (in/out): ParameterList for communication between control routine
                              and elements
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


      /*!
      \brief Evaluate a Neumann boundary condition

      this method evaluates a surface Neumann condition on the solid3 element

      \param params (in/out)    : ParameterList for communication between control routine
                                  and elements
      \param discretization (in): A reference to the underlying discretization
      \param condition (in)     : The condition to be evaluated
      \param lm (in)            : location vector of this element
      \param elevec1 (out)      : vector to be filled by element. If NULL on input,

      \return 0 if successful, negative otherwise
      */
      virtual int EvaluateNeumann(Teuchos::ParameterList& params,
          DRT::Discretization& discretization, DRT::Condition& condition, std::vector<int>& lm,
          Epetra_SerialDenseVector& elevec1, Epetra_SerialDenseMatrix* elemat1 = NULL);


     protected:
      //! action parameters recognized by this element
      enum ActionType
      {
        none,
        calc_struct_linstiff,
        calc_struct_nlnstiff,
        calc_struct_internalforce,
        calc_struct_linstiffmass,
        calc_struct_nlnstiffmass,
        calc_struct_nlnstifflmass,
        calc_struct_stress,
        postprocess_stress,
        calc_struct_eleload,
        calc_struct_fsiload,
        calc_struct_update_istep,
        calc_struct_reset_istep,  //!< reset elementwise internal variables
                                  //!< during iteration to last converged state
        calc_struct_B_and_F,
        multi_calc_dens,
        multi_readrestart
      };

      enum StressType
      {
        so_tet4_stress_none,
        so_tet4_stress_gpxyz,
        so_tet4_stress_gprst,
        so_tet4_stress_gp123,
        so_tet4_stress_ndxyz,
        so_tet4_stress_ndrst,
        so_tet4_stress_nd123
      };

      //! type of stress evaluation
      StressType stresstype_;

      //! number of the material law
      int material_;

      //! vector of history variables for each Gauss point
      std::vector<Teuchos::RCP<MAT::Material>> mat_;

      //! volume of element
      double V_;

      //! derivs of shape functions wrt to material coords
      LINALG::Matrix<4, 3> nxyz_;

      //! element's def-gradient
      LINALG::Matrix<3, 3> F_;

      // internal calculation methods

      // return volume of element
      inline double Vol() { return V_; }

     public:
      inline LINALG::Matrix<4, 3>& Nxyz() { return nxyz_; }

     protected:
      inline LINALG::Matrix<3, 3>& F() { return F_; }

      static LINALG::Matrix<3, 3> BuildF(LINALG::Matrix<4, 3>& xdisp, LINALG::Matrix<4, 3>& nxyz)
      {
        LINALG::Matrix<3, 3> F(false);
        F.MultiplyTN(xdisp, nxyz);
        F(0, 0) += 1.0;
        F(1, 1) += 1.0;
        F(2, 2) += 1.0;
        return F;
      }

      template <typename T>
      static LINALG::TMatrix<T, 3, 3> BuildF(
          LINALG::TMatrix<T, 4, 3>& xdisp, LINALG::Matrix<4, 3>& nxyz)
      {
        LINALG::TMatrix<T, 3, 3> F(true);
        for (int i = 0; i < 3; ++i)
        {
          for (int j = 0; j < 3; ++j)
            for (int k = 0; k < 4; ++k) F(i, j) += xdisp(k, i) * nxyz(k, j);
          F(i, i) += 1.0;
        }
        return F;
      }

      // don't want = operator
      NStet& operator=(const NStet& old);

      // init the inverse of the jacobian and its determinant
      // in the material configuration
      virtual void InitElement();

      //! Shape functions
      inline void ShapeFunction(LINALG::Matrix<4, 1>& funct, const double& e1, const double& e2,
          const double& e3, const double& e4)
      {
        // shape function is N_i = xsi_i
        funct(0) = e1;
        funct(1) = e2;
        funct(2) = e3;
        funct(3) = e4;
        return;
      }

      //! Shape function derivatives
      inline void ShapeFunctionDerivatives(LINALG::Matrix<4, 4>& deriv)
      {
        // Ni,j = 1.0 for i==j, otherwise 0.0
        deriv.Clear();
        deriv(0, 0) = 1.0;
        deriv(1, 1) = 1.0;
        deriv(2, 2) = 1.0;
        deriv(3, 3) = 1.0;
        return;
      }

      //! standards displ. tet4 calc routine
      virtual void nstetnlnstiffmass(std::vector<int>& lm, std::vector<double>& disp,
          std::vector<double>& residual, LINALG::Matrix<12, 12>* stiffmatrix,
          LINALG::Matrix<12, 12>* massmatrix, LINALG::Matrix<12, 1>* force,
          LINALG::Matrix<1, 6>* elestress, LINALG::Matrix<1, 6>* elestrain,
          const INPAR::STR::StressType iostress, const INPAR::STR::StrainType iostrain);


      //! lump mass matrix (bborn 07/08)
      void nstetlumpmass(LINALG::Matrix<12, 12>* emass);


      void so_nstet_expol(LINALG::Matrix<1, 6>& stresses, LINALG::Matrix<4, 6>& nodalstresses);



      void SelectMaterial(LINALG::Matrix<6, 1>& stress, LINALG::Matrix<6, 6>& cmat, double& density,
          LINALG::Matrix<6, 1>& glstrain, LINALG::Matrix<3, 3>& defgrd, int gp);


      //! @name Multi-scale related stuff

      /// Determine a homogenized material density for multi-scale analyses by averaging over the
      /// initial volume
      void nstet_homog(Teuchos::ParameterList& params);

      /// Read restart on the microscale
      void nstet_read_restart_multi();

      //@}

    };  // class NStet



  }  // namespace ELEMENTS
}  // namespace DRT

#endif  // #ifndef SO_NSTET_H_
