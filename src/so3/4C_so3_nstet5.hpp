/*----------------------------------------------------------------------*/
/*! \file

\brief NStet5 element

\level 2

*----------------------------------------------------------------------*/
#ifndef FOUR_C_SO3_NSTET5_HPP
#define FOUR_C_SO3_NSTET5_HPP


#include "4C_config.hpp"

#include "4C_fem_general_element.hpp"
#include "4C_fem_general_elementtype.hpp"
#include "4C_inpar_structure.hpp"

// #define PUSO_NSTET5              ///< run the Puso&Solberg style 5-node tet
#define ALPHA_NSTET5 0.1  ///< stabilization parameter for vol-dev split stabilization & Puso-tet

FOUR_C_NAMESPACE_OPEN

// forward declarations
namespace Core::LinAlg
{
  class SerialDenseVector;
  class SerialDenseMatrix;
}  // namespace Core::LinAlg

namespace Core::FE
{
  class Discretization;
}  // namespace Core::FE

namespace Discret
{

  namespace ELEMENTS
  {
    class StructuralLine;
    class StructuralSurface;
    class StructuralVolume;
    class NStet5;
    class PreStress;


    //=======================================================================
    //=======================================================================

    class NStet5Type : public Core::Elements::ElementType
    {
      //! allow NStet5 element to access the nodal data
      friend class Discret::ELEMENTS::NStet5;

     public:
      std::string Name() const override { return "NStet5Type"; }

      static NStet5Type& Instance();

      Core::Communication::ParObject* Create(const std::vector<char>& data) override;

      Teuchos::RCP<Core::Elements::Element> Create(const std::string eletype,
          const std::string eledistype, const int id, const int owner) override;

      Teuchos::RCP<Core::Elements::Element> Create(const int id, const int owner) override;

      int initialize(Core::FE::Discretization& dis) override;

      void pre_evaluate(Core::FE::Discretization& dis, Teuchos::ParameterList& p,
          Teuchos::RCP<Core::LinAlg::SparseOperator> systemmatrix1,
          Teuchos::RCP<Core::LinAlg::SparseOperator> systemmatrix2,
          Teuchos::RCP<Epetra_Vector> systemvector1, Teuchos::RCP<Epetra_Vector> systemvector2,
          Teuchos::RCP<Epetra_Vector> systemvector3) override;

      void nodal_block_information(
          Core::Elements::Element* dwele, int& numdf, int& dimns, int& nv, int& np) override;

      Core::LinAlg::SerialDenseMatrix ComputeNullSpace(
          Core::Nodes::Node& node, const double* x0, const int numdof, const int dimnsp) override;

      void setup_element_definition(
          std::map<std::string, std::map<std::string, Input::LineDefinition>>& definitions)
          override;

     private:
      static NStet5Type instance_;

      std::string get_element_type_string() const { return "NSTET5"; }

      //! map of row nodes adjacent to NStet5 elements
      std::map<int, Core::Nodes::Node*> noderids_;

      //! map of column NStet5 elements
      std::map<int, Discret::ELEMENTS::NStet5*> elecids_;

      //! nodal stresses and strains for output
      Teuchos::RCP<Epetra_MultiVector> nstress_;
      Teuchos::RCP<Epetra_MultiVector> nstrain_;

      //! map of nodes forming a patch around a node (includes center node)
      // std::map<centernodeid,std::map<nodeid,nodeptr> >
      std::map<int, std::map<int, Core::Nodes::Node*>> adjnode_;

      //! vector of elements adjacent to each row node
      // std::map<centernodeid,vector of ptrs to elements>
      std::map<int, std::vector<Discret::ELEMENTS::NStet5*>> adjele_;

      //! map of subelements of adjacent elements
      // std::map<centernodeid,std::map<eleid,vector of local subele numbers> >
      std::map<int, std::map<int, std::vector<int>>> adjsubele_;

      //! map of location vectors for patch around a node
      // std::map<centernodeid,vector of degrees of freedom on patch including element dofs>
      std::map<int, std::vector<int>> adjlm_;

      //! map of location vectors for patch around a node
      std::map<int, std::vector<std::vector<std::vector<int>>>> lmlm_;

      void init_elementsand_maps(std::map<int, Discret::ELEMENTS::NStet5*>& elecids,
          std::map<int, Core::Nodes::Node*>& noderids, const int myrank, const int numproc,
          Core::FE::Discretization& dis);

      void init_adjacency(std::map<int, Discret::ELEMENTS::NStet5*>& elecids,
          std::map<int, Core::Nodes::Node*>& noderids,
          std::map<int, std::vector<Discret::ELEMENTS::NStet5*>>& adjele,
          std::map<int, std::map<int, Core::Nodes::Node*>>& adjnode,
          std::map<int, std::vector<int>>& adjlm,
          std::map<int, std::map<int, std::vector<int>>>& adjsubele,
          std::map<int, std::vector<std::vector<std::vector<int>>>>& adjlmlm,
          Core::FE::Discretization& dis);


      void element_deformation_gradient(Core::FE::Discretization& dis);

      void nodal_integration(Core::LinAlg::SerialDenseMatrix* stiff,
          Core::LinAlg::SerialDenseVector* force, std::map<int, Core::Nodes::Node*>& adjnode,
          std::vector<Discret::ELEMENTS::NStet5*>& adjele,
          std::map<int, std::vector<int>>& adjsubele, std::vector<int>& lm,
          std::vector<std::vector<std::vector<int>>>& lmlm, const Epetra_Vector& disp,
          Core::FE::Discretization& dis, std::vector<double>* nodalstress,
          std::vector<double>* nodalstrain, const Inpar::STR::StressType iostress,
          const Inpar::STR::StrainType iostrain);


      void select_material(const Teuchos::RCP<Core::Mat::Material>& mat,
          Core::LinAlg::Matrix<6, 1>& stress, Core::LinAlg::Matrix<6, 6>& cmat, double& density,
          Core::LinAlg::Matrix<6, 1>& glstrain, Core::LinAlg::Matrix<3, 3>& defgrd, int gp,
          const int eleGID);

      // compute deviatoric stresses and tangent
      static void dev_stress_tangent(Core::LinAlg::Matrix<6, 1>& Sdev,
          Core::LinAlg::Matrix<6, 6>& CCdev, Core::LinAlg::Matrix<6, 6>& CC,
          const Core::LinAlg::Matrix<6, 1>& S, const Core::LinAlg::Matrix<3, 3>& C);

      void strain_output(const Inpar::STR::StrainType iostrain, std::vector<double>& nodalstrain,
          Core::LinAlg::Matrix<3, 3>& F, const double& detF, const double volweight,
          const double devweight);

      void strain_output(const Inpar::STR::StrainType iostrain, std::vector<double>& nodalstrain,
          Core::LinAlg::Matrix<3, 3>& F, Core::LinAlg::Matrix<6, 1>& glstrain, const double weight);

      void stress_output(const Inpar::STR::StressType iostress, std::vector<double>& nodalstress,
          Core::LinAlg::Matrix<6, 1>& stress, Core::LinAlg::Matrix<3, 3>& F, const double& detF);


      //! build deformation gradient
      template <typename T>
      static Core::LinAlg::Matrix<3, 3, T> t_build_f(
          const Core::LinAlg::Matrix<4, 3, T>& xdisp, const Core::LinAlg::Matrix<4, 3>& nxyz)
      {
        Core::LinAlg::Matrix<3, 3, T> F(true);
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
      T determinant(Core::LinAlg::Matrix<3, 3, T>& A)
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
    class NStet5 : public Core::Elements::Element
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
      inline Core::Elements::Element* Clone() const override
      {
        return new Discret::ELEMENTS::NStet5(*this);
      }

      /*!
      \brief Get shape type of element
      */
      inline Core::FE::CellType Shape() const override { return Core::FE::CellType::tet4; }

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
      std::vector<Teuchos::RCP<Core::Elements::Element>> Lines() override;

      /*!
      \brief Get vector of Teuchos::RCPs to the surfaces of this element

      */
      std::vector<Teuchos::RCP<Core::Elements::Element>> Surfaces() override;

      /*!
      \brief Return unique ParObject id

      every class implementing ParObject needs a unique id defined at the
      top of this file.
      */
      int UniqueParObjectId() const override { return NStet5Type::Instance().UniqueParObjectId(); }

      /*!
      \brief Pack this class so it can be communicated

      \ref pack and \ref unpack are used to communicate this element

      */
      void pack(Core::Communication::PackBuffer& data) const override;

      /*!
      \brief Unpack data from a char vector into this class

      \ref pack and \ref unpack are used to communicate this element

      */
      void unpack(const std::vector<char>& data) override;


      /*!
      \brief Get number of degrees of freedom of a certain node
             (implements pure virtual Core::Elements::Element)

      The element decides how many degrees of freedom its nodes must have.
      As this may vary along a simulation, the element can redecide the
      number of degrees of freedom per node along the way for each of it's nodes
      separately.
      */
      inline int NumDofPerNode(const Core::Nodes::Node& node) const override { return 3; }

      /*!
      \brief The 3 degrees of freedom of the center node
      */
      inline int num_dof_per_element() const override { return 3; }

      /*!
      \brief Print this element
      */
      void print(std::ostream& os) const override;

      NStet5Type& ElementType() const override { return NStet5Type::Instance(); }

      /*!
      \brief Read input for this element
      */
      bool ReadElement(const std::string& eletype, const std::string& distype,
          Input::LineDefinition* linedef) override;

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
      int evaluate(Teuchos::ParameterList& params, Core::FE::Discretization& discretization,
          std::vector<int>& lm, Core::LinAlg::SerialDenseMatrix& elemat1,
          Core::LinAlg::SerialDenseMatrix& elemat2, Core::LinAlg::SerialDenseVector& elevec1,
          Core::LinAlg::SerialDenseVector& elevec2,
          Core::LinAlg::SerialDenseVector& elevec3) override;


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
      int evaluate_neumann(Teuchos::ParameterList& params, Core::FE::Discretization& discretization,
          Core::Conditions::Condition& condition, std::vector<int>& lm,
          Core::LinAlg::SerialDenseVector& elevec1,
          Core::LinAlg::SerialDenseMatrix* elemat1 = nullptr) override;


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
      std::vector<Teuchos::RCP<Core::Mat::Material>> mat_;

      //! volume of element
      double V_;
      inline double vol() { return V_; }

      ///----------------------------------------------- prestressing switch & time
      /// prestressing assume 4 Gausspoints per element,
      /// that is, one per sub-element
      Inpar::STR::PreStress pstype_;
      double pstime_;
      double time_;
      /// Prestressing object
      Teuchos::RCP<Discret::ELEMENTS::PreStress> prestress_;
      /// compute Jacobian mapping wrt to deformed configuration
      void update_jacobian_mapping(
          const std::vector<double>& disp, Discret::ELEMENTS::PreStress& prestress);
      /// compute defgrd in all gp for given disp
      void def_gradient(const std::vector<double>& disp, Core::LinAlg::SerialDenseMatrix& gpdefgrd,
          Discret::ELEMENTS::PreStress& prestress);

      //---------------------------------------------- quantities related to subtets
      //! nodal connectivity  of subelements
      int sublm_[16]{};
      inline const int* sub_lm(int i) const { return &(sublm_[i * 4]); }

      //! coordinates of middle node
      double midX_[3]{};
      inline const double* mid_x() const { return midX_; }

      //! derivatives of shape functions for subtets
      Core::LinAlg::Matrix<4, 3> subnxyz_[4];
      inline const Core::LinAlg::Matrix<4, 3>& sub_nxyz(int i) const { return subnxyz_[i]; }

      //! reference volume of subelements
      double subV_[4]{};
      inline const double& sub_v(int i) const { return subV_[i]; }

      Core::LinAlg::Matrix<3, 3> subF_[4];
      inline Core::LinAlg::Matrix<3, 3>& sub_f(int i) { return subF_[i]; }
      //---------------------------------------------------------------------------



      inline static Core::LinAlg::Matrix<3, 3> build_f(
          const Core::LinAlg::Matrix<4, 3>& xdisp, const Core::LinAlg::Matrix<4, 3>& nxyz)
      {
        Core::LinAlg::Matrix<3, 3> F(false);
        F.multiply_tn(xdisp, nxyz);
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
      inline void shape_function(Core::LinAlg::Matrix<4, 1>& funct, const double& e1,
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
      inline void shape_function_derivatives(Core::LinAlg::Matrix<4, 4>& deriv)
      {
        // Ni,j = 1.0 for i==j, otherwise 0.0
        deriv.clear();
        deriv(0, 0) = 1.0;
        deriv(1, 1) = 1.0;
        deriv(2, 2) = 1.0;
        deriv(3, 3) = 1.0;
        return;
      }

      //! standards displ. tet4 calc routine
      virtual void nstet5nlnstiffmass(std::vector<int>& lm, std::vector<double>& disp,
          Core::LinAlg::Matrix<15, 15>* stiffmatrix, Core::LinAlg::Matrix<15, 15>* massmatrix,
          Core::LinAlg::Matrix<15, 1>* force, Core::LinAlg::Matrix<1, 6>* elestress,
          Core::LinAlg::Matrix<1, 6>* elestrain, const Inpar::STR::StressType iostress,
          const Inpar::STR::StrainType iostrain);


      //! lump mass matrix (bborn 07/08)
      void nstet5lumpmass(Core::LinAlg::Matrix<15, 15>* emass);


      void so_nstet5_expol(
          Core::LinAlg::Matrix<1, 6>& stresses, Core::LinAlg::Matrix<4, 6>& nodalstresses);



      void select_material(Core::LinAlg::Matrix<6, 1>& stress, Core::LinAlg::Matrix<6, 6>& cmat,
          double& density, Core::LinAlg::Matrix<6, 1>& glstrain, Core::LinAlg::Matrix<3, 3>& defgrd,
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
}  // namespace Discret

FOUR_C_NAMESPACE_CLOSE

#endif
