/*----------------------------------------------------------------------*/
/*! \file

\brief NStet5 element

\level 2

\maintainer Christoph Meier
*----------------------------------------------------------------------*/
#ifndef SO_NSTET5_H
#define SO_NSTET5_H

//#include "Sacado.hpp"

#include "../drt_lib/drt_element.H"
#include "../drt_lib/drt_elementtype.H"
#include "../drt_inpar/inpar_structure.H"

//#define PUSO_NSTET5              ///< run the Puso&Solberg style 5-node tet
#define ALPHA_NSTET5 0.1  ///< stabilization parameter for vol-dev split stabilization & Puso-tet



// forward declarations
namespace LINALG
{
  class SerialDenseVector;
  class SerialDenseMatrix;
}  // namespace LINALG
namespace DRT
{
  class Discretization;

  namespace ELEMENTS
  {
    class StructuralLine;
    class StructuralSurface;
    class StructuralVolume;
    class NStet5;
    class InvDesign;
    class PreStress;


    //=======================================================================
    //=======================================================================

    class NStet5Type : public DRT::ElementType
    {
      //! allow NStet5 element to access the nodal data
      friend class DRT::ELEMENTS::NStet5;

     public:
      std::string Name() const { return "NStet5Type"; }

      static NStet5Type& Instance();

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
      static NStet5Type instance_;

      //! map of row nodes adjacent to NStet5 elements
      std::map<int, DRT::Node*> noderids_;

      //! map of column NStet5 elements
      std::map<int, DRT::ELEMENTS::NStet5*> elecids_;

      //! nodal stresses and strains for output
      Teuchos::RCP<Epetra_MultiVector> nstress_;
      Teuchos::RCP<Epetra_MultiVector> nstrain_;

      //! map of nodes forming a patch around a node (includes center node)
      // std::map<centernodeid,std::map<nodeid,nodeptr> >
      std::map<int, std::map<int, DRT::Node*>> adjnode_;

      //! vector of elements adjacent to each row node
      // std::map<centernodeid,vector of ptrs to elements>
      std::map<int, std::vector<DRT::ELEMENTS::NStet5*>> adjele_;

      //! map of subelements of adjacent elements
      // std::map<centernodeid,std::map<eleid,vector of local subele numbers> >
      std::map<int, std::map<int, std::vector<int>>> adjsubele_;

      //! map of location vectors for patch around a node
      // std::map<centernodeid,vector of degrees of freedom on patch including element dofs>
      std::map<int, std::vector<int>> adjlm_;

      //! map of location vectors for patch around a node
      std::map<int, std::vector<std::vector<std::vector<int>>>> lmlm_;

      void InitElementsandMaps(std::map<int, DRT::ELEMENTS::NStet5*>& elecids,
          std::map<int, DRT::Node*>& noderids, const int myrank, const int numproc,
          DRT::Discretization& dis);

      void InitAdjacency(std::map<int, DRT::ELEMENTS::NStet5*>& elecids,
          std::map<int, DRT::Node*>& noderids,
          std::map<int, std::vector<DRT::ELEMENTS::NStet5*>>& adjele,
          std::map<int, std::map<int, DRT::Node*>>& adjnode, std::map<int, std::vector<int>>& adjlm,
          std::map<int, std::map<int, std::vector<int>>>& adjsubele,
          std::map<int, std::vector<std::vector<std::vector<int>>>>& adjlmlm,
          DRT::Discretization& dis);


      void ElementDeformationGradient(DRT::Discretization& dis);

      void NodalIntegration(Epetra_SerialDenseMatrix* stiff, Epetra_SerialDenseVector* force,
          std::map<int, DRT::Node*>& adjnode, std::vector<DRT::ELEMENTS::NStet5*>& adjele,
          std::map<int, std::vector<int>>& adjsubele, std::vector<int>& lm,
          std::vector<std::vector<std::vector<int>>>& lmlm, const Epetra_Vector& disp,
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


      //! build deformation gradient
      template <typename T>
      static LINALG::TMatrix<T, 3, 3> TBuildF(
          const LINALG::TMatrix<T, 4, 3>& xdisp, const LINALG::Matrix<4, 3>& nxyz)
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

      //! build and return determinant of a 3x3 matrix
      template <typename T>
      T Determinant3x3(LINALG::TMatrix<T, 3, 3>& A)
      {
        T b00 = A(0, 0);
        T b01 = A(0, 1);
        T b02 = A(0, 2);
        T b10 = A(1, 0);
        T b11 = A(1, 1);
        T b12 = A(1, 2);
        T b20 = A(2, 0);
        T b21 = A(2, 1);
        T b22 = A(2, 2);
        T a = b11 * b22 - b21 * b12;
        T b = -b10 * b22 + b20 * b12;
        T c = b10 * b21 - b20 * b11;
        T det = b00 * a + b01 * b + b02 * c;
        return det;
      }
    };

    //----------------------------------------------------------------------------
    /*!
    \brief A nodal-averaged strain 5-noded tet element

    */
    //----------------------------------------------------------------------------
    class NStet5 : public DRT::Element
    {
     public:
      friend class NStet5Type;
      // typedef Sacado::Fad::DFad<double> FAD;
      // typedef Sacado::Fad::DFad< Sacado::Fad::DFad<double> > FADFAD;

      /*!
      \brief Standard Constructor

      \param id : A unique global id
      */
      NStet5(int id, int owner);

      /*!
      \brief Copy Constructor

      Makes a deep copy of a Element

      */
      NStet5(const NStet5& old);

      /*!
      \brief Deep copy this instance of Solid3 and return pointer to the copy

      The Clone() method is used from the virtual base class Element in cases
      where the type of the derived class is unknown and a copy-ctor is needed

      */
      inline DRT::Element* Clone() const { return new DRT::ELEMENTS::NStet5(*this); }

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
      virtual int UniqueParObjectId() const { return NStet5Type::Instance().UniqueParObjectId(); }

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
      virtual ~NStet5();

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
      \brief The 3 degrees of freedom of the center node
      */
      virtual inline int NumDofPerElement() const { return 3; }

      /*!
      \brief Print this element
      */
      void Print(std::ostream& os) const;

      virtual NStet5Type& ElementType() const { return NStet5Type::Instance(); }

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
        prestress_update,
        inversedesign_update,
        inversedesign_switch,
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
      inline double Vol() { return V_; }

      ///----------------------------------------------- prestressing switch & time
      /// prestressing and invdesign assume 4 Gausspoints per element,
      /// that is, one per sub-element
      INPAR::STR::PreStress pstype_;
      double pstime_;
      double time_;
      /// Prestressing object
      Teuchos::RCP<DRT::ELEMENTS::PreStress> prestress_;
      /// compute Jacobian mapping wrt to deformed configuration
      void UpdateJacobianMapping(
          const std::vector<double>& disp, DRT::ELEMENTS::PreStress& prestress);
      /// compute defgrd in all gp for given disp
      void DefGradient(const std::vector<double>& disp, Epetra_SerialDenseMatrix& gpdefgrd,
          DRT::ELEMENTS::PreStress& prestress);

      friend class DRT::ELEMENTS::InvDesign;
      Teuchos::RCP<DRT::ELEMENTS::InvDesign> invdesign_;

      //---------------------------------------------- quantities related to subtets
      //! nodal connectivity  of subelements
      int sublm_[16];
      inline const int* SubLM(int i) const { return &(sublm_[i * 4]); }

      //! coordinates of middle node
      double midX_[3];
      inline const double* MidX() const { return midX_; }

      //! derivatives of shape functions for subtets
      LINALG::Matrix<4, 3> subnxyz_[4];
      inline const LINALG::Matrix<4, 3>& SubNxyz(int i) const { return subnxyz_[i]; }

      //! reference volume of subelements
      double subV_[4];
      inline const double& SubV(int i) const { return subV_[i]; }

      LINALG::Matrix<3, 3> subF_[4];
      inline LINALG::Matrix<3, 3>& SubF(int i) { return subF_[i]; }
      //---------------------------------------------------------------------------



      inline static LINALG::Matrix<3, 3> BuildF(
          const LINALG::Matrix<4, 3>& xdisp, const LINALG::Matrix<4, 3>& nxyz)
      {
        LINALG::Matrix<3, 3> F(false);
        F.MultiplyTN(xdisp, nxyz);
        F(0, 0) += 1.0;
        F(1, 1) += 1.0;
        F(2, 2) += 1.0;
        return F;
      }


      // don't want = operator
      NStet5& operator=(const NStet5& old);

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
      virtual void nstet5nlnstiffmass(std::vector<int>& lm, std::vector<double>& disp,
          LINALG::Matrix<15, 15>* stiffmatrix, LINALG::Matrix<15, 15>* massmatrix,
          LINALG::Matrix<15, 1>* force, LINALG::Matrix<1, 6>* elestress,
          LINALG::Matrix<1, 6>* elestrain, const INPAR::STR::StressType iostress,
          const INPAR::STR::StrainType iostrain);


      //! lump mass matrix (bborn 07/08)
      void nstet5lumpmass(LINALG::Matrix<15, 15>* emass);


      void so_nstet5_expol(LINALG::Matrix<1, 6>& stresses, LINALG::Matrix<4, 6>& nodalstresses);



      void SelectMaterial(LINALG::Matrix<6, 1>& stress, LINALG::Matrix<6, 6>& cmat, double& density,
          LINALG::Matrix<6, 1>& glstrain, LINALG::Matrix<3, 3>& defgrd, int gp);


      //! @name Multi-scale related stuff

      /// Determine a homogenized material density for multi-scale analyses by averaging over the
      /// initial volume
      void nstet5_homog(Teuchos::ParameterList& params);

      /// Read restart on the microscale
      void nstet5_read_restart_multi();

      //@}

    };  // class NStet5



  }  // namespace ELEMENTS
}  // namespace DRT

#endif  // #ifndef SO_NSTET5_H_
