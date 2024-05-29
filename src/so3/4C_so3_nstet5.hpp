/*----------------------------------------------------------------------*/
/*! \file

\brief NStet5 element

\level 2

*----------------------------------------------------------------------*/
#ifndef FOUR_C_SO3_NSTET5_HPP
#define FOUR_C_SO3_NSTET5_HPP


#include "4C_config.hpp"

#include "4C_discretization_fem_general_element.hpp"
#include "4C_discretization_fem_general_elementtype.hpp"
#include "4C_inpar_structure.hpp"

// #define PUSO_NSTET5              ///< run the Puso&Solberg style 5-node tet
#define ALPHA_NSTET5 0.1  ///< stabilization parameter for vol-dev split stabilization & Puso-tet

FOUR_C_NAMESPACE_OPEN

// forward declarations
namespace CORE::LINALG
{
  class SerialDenseVector;
  class SerialDenseMatrix;
}  // namespace CORE::LINALG
namespace DRT
{
  class Discretization;

  namespace ELEMENTS
  {
    class StructuralLine;
    class StructuralSurface;
    class StructuralVolume;
    class NStet5;
    class PreStress;


    //=======================================================================
    //=======================================================================

    class NStet5Type : public CORE::Elements::ElementType
    {
      //! allow NStet5 element to access the nodal data
      friend class DRT::ELEMENTS::NStet5;

     public:
      std::string Name() const override { return "NStet5Type"; }

      static NStet5Type& Instance();

      CORE::COMM::ParObject* Create(const std::vector<char>& data) override;

      Teuchos::RCP<CORE::Elements::Element> Create(const std::string eletype,
          const std::string eledistype, const int id, const int owner) override;

      Teuchos::RCP<CORE::Elements::Element> Create(const int id, const int owner) override;

      int Initialize(DRT::Discretization& dis) override;

      void pre_evaluate(DRT::Discretization& dis, Teuchos::ParameterList& p,
          Teuchos::RCP<CORE::LINALG::SparseOperator> systemmatrix1,
          Teuchos::RCP<CORE::LINALG::SparseOperator> systemmatrix2,
          Teuchos::RCP<Epetra_Vector> systemvector1, Teuchos::RCP<Epetra_Vector> systemvector2,
          Teuchos::RCP<Epetra_Vector> systemvector3) override;

      void nodal_block_information(
          CORE::Elements::Element* dwele, int& numdf, int& dimns, int& nv, int& np) override;

      CORE::LINALG::SerialDenseMatrix ComputeNullSpace(
          DRT::Node& node, const double* x0, const int numdof, const int dimnsp) override;

      void setup_element_definition(
          std::map<std::string, std::map<std::string, INPUT::LineDefinition>>& definitions)
          override;

     private:
      static NStet5Type instance_;

      std::string get_element_type_string() const { return "NSTET5"; }

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

      void init_elementsand_maps(std::map<int, DRT::ELEMENTS::NStet5*>& elecids,
          std::map<int, DRT::Node*>& noderids, const int myrank, const int numproc,
          DRT::Discretization& dis);

      void init_adjacency(std::map<int, DRT::ELEMENTS::NStet5*>& elecids,
          std::map<int, DRT::Node*>& noderids,
          std::map<int, std::vector<DRT::ELEMENTS::NStet5*>>& adjele,
          std::map<int, std::map<int, DRT::Node*>>& adjnode, std::map<int, std::vector<int>>& adjlm,
          std::map<int, std::map<int, std::vector<int>>>& adjsubele,
          std::map<int, std::vector<std::vector<std::vector<int>>>>& adjlmlm,
          DRT::Discretization& dis);


      void element_deformation_gradient(DRT::Discretization& dis);

      void nodal_integration(CORE::LINALG::SerialDenseMatrix* stiff,
          CORE::LINALG::SerialDenseVector* force, std::map<int, DRT::Node*>& adjnode,
          std::vector<DRT::ELEMENTS::NStet5*>& adjele, std::map<int, std::vector<int>>& adjsubele,
          std::vector<int>& lm, std::vector<std::vector<std::vector<int>>>& lmlm,
          const Epetra_Vector& disp, DRT::Discretization& dis, std::vector<double>* nodalstress,
          std::vector<double>* nodalstrain, const INPAR::STR::StressType iostress,
          const INPAR::STR::StrainType iostrain);


      void select_material(const Teuchos::RCP<CORE::MAT::Material>& mat,
          CORE::LINALG::Matrix<6, 1>& stress, CORE::LINALG::Matrix<6, 6>& cmat, double& density,
          CORE::LINALG::Matrix<6, 1>& glstrain, CORE::LINALG::Matrix<3, 3>& defgrd, int gp,
          const int eleGID);

      // compute deviatoric stresses and tangent
      static void dev_stress_tangent(CORE::LINALG::Matrix<6, 1>& Sdev,
          CORE::LINALG::Matrix<6, 6>& CCdev, CORE::LINALG::Matrix<6, 6>& CC,
          const CORE::LINALG::Matrix<6, 1>& S, const CORE::LINALG::Matrix<3, 3>& C);

      void strain_output(const INPAR::STR::StrainType iostrain, std::vector<double>& nodalstrain,
          CORE::LINALG::Matrix<3, 3>& F, const double& detF, const double volweight,
          const double devweight);

      void strain_output(const INPAR::STR::StrainType iostrain, std::vector<double>& nodalstrain,
          CORE::LINALG::Matrix<3, 3>& F, CORE::LINALG::Matrix<6, 1>& glstrain, const double weight);

      void stress_output(const INPAR::STR::StressType iostress, std::vector<double>& nodalstress,
          CORE::LINALG::Matrix<6, 1>& stress, CORE::LINALG::Matrix<3, 3>& F, const double& detF);


      //! build deformation gradient
      template <typename T>
      static CORE::LINALG::Matrix<3, 3, T> t_build_f(
          const CORE::LINALG::Matrix<4, 3, T>& xdisp, const CORE::LINALG::Matrix<4, 3>& nxyz)
      {
        CORE::LINALG::Matrix<3, 3, T> F(true);
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
      T determinant(CORE::LINALG::Matrix<3, 3, T>& A)
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
    class NStet5 : public CORE::Elements::Element
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
      inline CORE::Elements::Element* Clone() const override
      {
        return new DRT::ELEMENTS::NStet5(*this);
      }

      /*!
      \brief Get shape type of element
      */
      inline CORE::FE::CellType Shape() const override { return CORE::FE::CellType::tet4; }

      /*!
      \brief Return number of volumes of this element
      */
      int NumVolume() const override { return 1; }

      /*!
      \brief Return number of surfaces of this element
      */
      int NumSurface() const override { return 4; }

      /*!
      \brief Return number of lines of this element
      */
      int NumLine() const override { return 6; }

      /*!
      \brief Get vector of Teuchos::RCPs to the lines of this element

      */
      std::vector<Teuchos::RCP<CORE::Elements::Element>> Lines() override;

      /*!
      \brief Get vector of Teuchos::RCPs to the surfaces of this element

      */
      std::vector<Teuchos::RCP<CORE::Elements::Element>> Surfaces() override;

      /*!
      \brief Return unique ParObject id

      every class implementing ParObject needs a unique id defined at the
      top of this file.
      */
      int UniqueParObjectId() const override { return NStet5Type::Instance().UniqueParObjectId(); }

      /*!
      \brief Pack this class so it can be communicated

      \ref Pack and \ref Unpack are used to communicate this element

      */
      void Pack(CORE::COMM::PackBuffer& data) const override;

      /*!
      \brief Unpack data from a char vector into this class

      \ref Pack and \ref Unpack are used to communicate this element

      */
      void Unpack(const std::vector<char>& data) override;


      /*!
      \brief Get number of degrees of freedom of a certain node
             (implements pure virtual CORE::Elements::Element)

      The element decides how many degrees of freedom its nodes must have.
      As this may vary along a simulation, the element can redecide the
      number of degrees of freedom per node along the way for each of it's nodes
      separately.
      */
      inline int NumDofPerNode(const DRT::Node& node) const override { return 3; }

      /*!
      \brief The 3 degrees of freedom of the center node
      */
      inline int num_dof_per_element() const override { return 3; }

      /*!
      \brief Print this element
      */
      void Print(std::ostream& os) const override;

      NStet5Type& ElementType() const override { return NStet5Type::Instance(); }

      /*!
      \brief Read input for this element
      */
      bool ReadElement(const std::string& eletype, const std::string& distype,
          INPUT::LineDefinition* linedef) override;

      /*!
      \brief Evaluate an element

      Evaluate so_tet4 element stiffness, mass, internal forces, etc.

      \param params (in/out): ParameterList for communication between control routine
                              and elements
      \param elemat1 (out)  : (stiffness-)matrix to be filled by element. If nullptr on input,
                              the controling method does not expect the element to fill
                              this matrix.
      \param elemat2 (out)  : (mass-)matrix to be filled by element. If nullptr on input,
                              the controling method does not expect the element to fill
                              this matrix.
      \param elevec1 (out)  : (internal force-)vector to be filled by element. If nullptr on input,
                              the controlling method does not expect the element
                              to fill this vector
      \param elevec2 (out)  : vector to be filled by element. If nullptr on input,
                              the controlling method does not expect the element
                              to fill this vector
      \param elevec3 (out)  : vector to be filled by element. If nullptr on input,
                              the controlling method does not expect the element
                              to fill this vector
      \return 0 if successful, negative otherwise
      */
      int Evaluate(Teuchos::ParameterList& params, DRT::Discretization& discretization,
          std::vector<int>& lm, CORE::LINALG::SerialDenseMatrix& elemat1,
          CORE::LINALG::SerialDenseMatrix& elemat2, CORE::LINALG::SerialDenseVector& elevec1,
          CORE::LINALG::SerialDenseVector& elevec2,
          CORE::LINALG::SerialDenseVector& elevec3) override;


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
      int evaluate_neumann(Teuchos::ParameterList& params, DRT::Discretization& discretization,
          CORE::Conditions::Condition& condition, std::vector<int>& lm,
          CORE::LINALG::SerialDenseVector& elevec1,
          CORE::LINALG::SerialDenseMatrix* elemat1 = nullptr) override;


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
        calc_struct_eleload,
        calc_struct_fsiload,
        calc_struct_update_istep,
        calc_struct_reset_istep,  //!< reset elementwise internal variables
                                  //!< during iteration to last converged state
        calc_struct_B_and_F,
        multi_calc_dens,
        prestress_update,
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
      std::vector<Teuchos::RCP<CORE::MAT::Material>> mat_;

      //! volume of element
      double V_;
      inline double vol() { return V_; }

      ///----------------------------------------------- prestressing switch & time
      /// prestressing assume 4 Gausspoints per element,
      /// that is, one per sub-element
      INPAR::STR::PreStress pstype_;
      double pstime_;
      double time_;
      /// Prestressing object
      Teuchos::RCP<DRT::ELEMENTS::PreStress> prestress_;
      /// compute Jacobian mapping wrt to deformed configuration
      void update_jacobian_mapping(
          const std::vector<double>& disp, DRT::ELEMENTS::PreStress& prestress);
      /// compute defgrd in all gp for given disp
      void def_gradient(const std::vector<double>& disp, CORE::LINALG::SerialDenseMatrix& gpdefgrd,
          DRT::ELEMENTS::PreStress& prestress);

      //---------------------------------------------- quantities related to subtets
      //! nodal connectivity  of subelements
      int sublm_[16]{};
      inline const int* sub_lm(int i) const { return &(sublm_[i * 4]); }

      //! coordinates of middle node
      double midX_[3]{};
      inline const double* mid_x() const { return midX_; }

      //! derivatives of shape functions for subtets
      CORE::LINALG::Matrix<4, 3> subnxyz_[4];
      inline const CORE::LINALG::Matrix<4, 3>& sub_nxyz(int i) const { return subnxyz_[i]; }

      //! reference volume of subelements
      double subV_[4]{};
      inline const double& sub_v(int i) const { return subV_[i]; }

      CORE::LINALG::Matrix<3, 3> subF_[4];
      inline CORE::LINALG::Matrix<3, 3>& sub_f(int i) { return subF_[i]; }
      //---------------------------------------------------------------------------



      inline static CORE::LINALG::Matrix<3, 3> build_f(
          const CORE::LINALG::Matrix<4, 3>& xdisp, const CORE::LINALG::Matrix<4, 3>& nxyz)
      {
        CORE::LINALG::Matrix<3, 3> F(false);
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
      virtual void init_element();

      //! Shape functions
      inline void shape_function(CORE::LINALG::Matrix<4, 1>& funct, const double& e1,
          const double& e2, const double& e3, const double& e4)
      {
        // shape function is N_i = xsi_i
        funct(0) = e1;
        funct(1) = e2;
        funct(2) = e3;
        funct(3) = e4;
        return;
      }

      //! Shape function derivatives
      inline void shape_function_derivatives(CORE::LINALG::Matrix<4, 4>& deriv)
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
          CORE::LINALG::Matrix<15, 15>* stiffmatrix, CORE::LINALG::Matrix<15, 15>* massmatrix,
          CORE::LINALG::Matrix<15, 1>* force, CORE::LINALG::Matrix<1, 6>* elestress,
          CORE::LINALG::Matrix<1, 6>* elestrain, const INPAR::STR::StressType iostress,
          const INPAR::STR::StrainType iostrain);


      //! lump mass matrix (bborn 07/08)
      void nstet5lumpmass(CORE::LINALG::Matrix<15, 15>* emass);


      void so_nstet5_expol(
          CORE::LINALG::Matrix<1, 6>& stresses, CORE::LINALG::Matrix<4, 6>& nodalstresses);



      void select_material(CORE::LINALG::Matrix<6, 1>& stress, CORE::LINALG::Matrix<6, 6>& cmat,
          double& density, CORE::LINALG::Matrix<6, 1>& glstrain, CORE::LINALG::Matrix<3, 3>& defgrd,
          int gp);


      //! @name Multi-scale related stuff

      /// Determine a homogenized material density for multi-scale analyses by averaging over the
      /// initial volume
      void nstet5_homog(Teuchos::ParameterList& params);

      /// Read restart on the microscale
      void nstet5_read_restart_multi();

      //@}

     private:
      std::string get_element_type_string() const { return "NSTET5"; }
    };  // class NStet5



  }  // namespace ELEMENTS
}  // namespace DRT

FOUR_C_NAMESPACE_CLOSE

#endif
