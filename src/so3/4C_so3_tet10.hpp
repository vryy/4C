/*----------------------------------------------------------------------*/
/*! \file
\brief Solid Tet10 Element
\level 1
*----------------------------------------------------------------------*/
#ifndef FOUR_C_SO3_TET10_HPP
#define FOUR_C_SO3_TET10_HPP


#include "4C_config.hpp"

#include "4C_comm_parobjectfactory.hpp"
#include "4C_discretization_fem_general_element.hpp"
#include "4C_discretization_fem_general_elementtype.hpp"
#include "4C_discretization_fem_general_utils_integration.hpp"
#include "4C_inpar_structure.hpp"
#include "4C_linalg_serialdensematrix.hpp"
#include "4C_mat_material_factory.hpp"
#include "4C_material_base.hpp"
#include "4C_so3_base.hpp"

// gee: note that these are also defined in so_integrator.H
#define NUMNOD_SOTET10 10       ///< number of nodes
#define NODDOF_SOTET10 3        ///< number of dofs per node
#define NUMDOF_SOTET10 30       ///< total dofs per element
#define NUMGPT_SOTET10 4        ///< number gauss points per element for stiffness integration
#define NUMGPT_MASS_SOTET10 11  ///< number gauss points per element for mass integration
#define NUMDIM_SOTET10 3        ///< number of dimensions
#define NUMCOORD_SOTET10 4      ///< number of shape function coordinates (ksi1-ksi4)
#define NUMNOD_SOTET10_FACE 6   ///< number of nodes on a TET10 face (which is a TRI6)
#define NUMGPT_SOTET10_FACE 3   ///< number of GP    on a TET10 face (which is a TRI6)

FOUR_C_NAMESPACE_OPEN

namespace Discret
{
  // forward declarations
  class Discretization;

  namespace ELEMENTS
  {
    // forward declarations
    class PreStress;

    class SoTet10Type : public Core::Elements::ElementType
    {
     public:
      std::string Name() const override { return "So_tet10Type"; }

      static SoTet10Type& Instance();

      Core::Communication::ParObject* Create(const std::vector<char>& data) override;

      Teuchos::RCP<Core::Elements::Element> Create(const std::string eletype,
          const std::string eledistype, const int id, const int owner) override;

      Teuchos::RCP<Core::Elements::Element> Create(const int id, const int owner) override;

      int Initialize(Discret::Discretization& dis) override;

      void nodal_block_information(
          Core::Elements::Element* dwele, int& numdf, int& dimns, int& nv, int& np) override;

      Core::LinAlg::SerialDenseMatrix ComputeNullSpace(
          Core::Nodes::Node& node, const double* x0, const int numdof, const int dimnsp) override;

      void setup_element_definition(
          std::map<std::string, std::map<std::string, Input::LineDefinition>>& definitions)
          override;

     private:
      static SoTet10Type instance_;

      std::string get_element_type_string() const { return "SOLIDT10"; }
    };

    /*!
    \brief A C++ version of the 10-node tet solid element

    a structural 10-node tetrahedral solid displacement element
    \author biehler
    */

    class SoTet10 : public SoBase
    {
     public:
      //! @name Friends
      friend class SoTet10Type;

      //@}
      //! @name Constructors and destructors and related methods

      /*!
      \brief Standard Constructor

      \param id : A unique global id
      \param owner : elements owning processor
      */
      SoTet10(int id, int owner);

      /*!
      \brief Copy Constructor

      Makes a deep copy of a Element

      */
      SoTet10(const SoTet10& old);

      /*!
      \brief Deep copy this instance of Solid3 and return pointer to the copy

      The Clone() method is used from the virtual base class Element in cases
      where the type of the derived class is unknown and a copy-ctor is needed

      */
      Core::Elements::Element* Clone() const override;

      /*!
      \brief Get shape type of element
      */
      Core::FE::CellType Shape() const override;

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

      virtual std::vector<double> element_center_refe_coords();

      /*!
      \brief Return unique ParObject id

      every class implementing ParObject needs a unique id defined at the
      top of this file.
      */
      int UniqueParObjectId() const override { return SoTet10Type::Instance().UniqueParObjectId(); }

      /*!
      \brief Pack this class so it can be communicated

      \ref Pack and \ref Unpack are used to communicate this element

      */
      void Pack(Core::Communication::PackBuffer& data) const override;

      /*!
      \brief Unpack data from a char vector into this class

      \ref Pack and \ref Unpack are used to communicate this element

      */
      void Unpack(const std::vector<char>& data) override;

      //@}

      //! @name Acess methods


      /*!
      \brief Get number of degrees of freedom of a certain node
             (implements pure virtual Core::Elements::Element)

      The element decides how many degrees of freedom its nodes must have.
      As this may vary along a simulation, the element can redecide the
      number of degrees of freedom per node along the way for each of it's nodes
      separately.
      */
      int NumDofPerNode(const Core::Nodes::Node& node) const override { return 3; }

      /*!
      \brief Get number of degrees of freedom per element
             (implements pure virtual Core::Elements::Element)

      The element decides how many element degrees of freedom it has.
      It can redecide along the way of a simulation.

      \note Element degrees of freedom mentioned here are dofs that are visible
            at the level of the total system of equations. Purely internal
            element dofs that are condensed internally should NOT be considered.
      */
      int num_dof_per_element() const override { return 0; }

      /*!
      \brief Print this element
      */
      void Print(std::ostream& os) const override;

      Core::Elements::ElementType& ElementType() const override { return SoTet10Type::Instance(); }

      //@}

      //! @name Input and Creation

      /*!
      \brief Read input for this element
       */
      /*!
      \brief Query names of element data to be visualized using BINIO

      The element fills the provided map with key names of
      visualization data the element wants to visualize AT THE CENTER
      of the element geometry. The values is supposed to be dimension of the
      data to be visualized. It can either be 1 (scalar), 3 (vector), 6 (sym. tensor)
      or 9 (nonsym. tensor)

      Example:
      \code
        // Name of data is 'Owner', dimension is 1 (scalar value)
        names.insert(std::pair<std::string,int>("Owner",1));
        // Name of data is 'StressesXYZ', dimension is 6 (sym. tensor value)
        names.insert(std::pair<std::string,int>("StressesXYZ",6));
      \endcode*/

      void VisNames(std::map<std::string, int>& names)
          override;  ///< to be filled with key names of data to visualize and with int dimensions

      /*!
      \brief Query data to be visualized using BINIO of a given name

      The method is supposed to call this base method to visualize the owner of
      the element.
      If the derived method recognizes a supported data name, it shall fill it
      with corresponding data.
      If it does NOT recognizes the name, it shall do nothing.

      \warning The method must not change size of data

      */
      bool VisData(
          const std::string& name,  ///< Name of data that is currently processed for visualization
          std::vector<double>&
              data  ///< d ata to be filled by element if element recognizes the name
          ) override;
      //@}

      //! @name Input and Creation

      /*!
      \brief Read input for this element
      */
      bool ReadElement(const std::string& eletype, const std::string& distype,
          Input::LineDefinition* linedef) override;

      //@}

      //! @name Evaluation

      /*!
      \brief Evaluate an element

      Evaluate so_tet10 element stiffness, mass, internal forces, etc.

      \return 0 if successful, negative otherwise
      */
      int Evaluate(
          Teuchos::ParameterList&
              params,  ///< ParameterList for communication between control routine and elements
          Discret::Discretization& discretization,  ///< pointer to discretization for de-assembly
          std::vector<int>& lm,                     ///< location matrix for de-assembly
          Core::LinAlg::SerialDenseMatrix&
              elemat1,  ///< (stiffness-)matrix to be filled by element.
          Core::LinAlg::SerialDenseMatrix& elemat2,  ///< (mass-)matrix to be filled by element.
          Core::LinAlg::SerialDenseVector&
              elevec1,  ///< (internal force-)vector to be filled by element
          Core::LinAlg::SerialDenseVector& elevec2,  ///< vector to be filled by element
          Core::LinAlg::SerialDenseVector& elevec3   ///< vector to be filled by element
          ) override;


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
      int evaluate_neumann(Teuchos::ParameterList& params, Discret::Discretization& discretization,
          Core::Conditions::Condition& condition, std::vector<int>& lm,
          Core::LinAlg::SerialDenseVector& elevec1,
          Core::LinAlg::SerialDenseMatrix* elemat1 = nullptr) override;

      //@}

     protected:
      //! action parameters recognized by so_tet10
      enum ActionType
      {
        none,
        calc_struct_linstiff,
        calc_struct_nlnstiff,
        calc_struct_internalforce,
        calc_struct_linstiffmass,
        calc_struct_nlnstiffmass,
        calc_struct_nlnstifflmass,  //!< internal force, its stiffness and lumped mass matrix
        calc_struct_stress,
        calc_struct_eleload,
        calc_struct_fsiload,
        calc_struct_update_istep,
        calc_struct_reset_istep,  //!< reset elementwise internal variables
                                  //!< during iteration to last converged state
        calc_struct_reset_all,    //!< reset elementwise internal variables
                                  //!< to state in the beginning of the computation
        calc_struct_energy,       // compute energy
        prestress_update,
        calc_global_gpstresses_map,
        struct_init_gauss_point_data_output,
        struct_gauss_point_data_output
      };


      //! vector of inverses of the jacobian in material frame
      std::vector<Core::LinAlg::Matrix<NUMDIM_SOTET10, NUMDIM_SOTET10>> invJ_;
      // determinant of Jacobian in material frame
      std::vector<double> detJ_;
      // stuff consistent mass matrix
      //! vector of inverses of the jacobian in material frame
      std::vector<Core::LinAlg::Matrix<NUMDIM_SOTET10, NUMDIM_SOTET10>> invJ_mass_;
      //! determinant of Jacobian in material frame
      std::vector<double> detJ_mass_;

      //! vector of partial derivatives in material frame
      // vector<Core::LinAlg::Matrix<NUMNOD_SOTET4,NUMDIM_SOTET4> > nxyz_;
      Core::LinAlg::Matrix<NUMNOD_SOTET10, NUMDIM_SOTET10> nxyz_;

      /// prestressing switch & time
      Inpar::STR::PreStress pstype_;
      double pstime_;
      double time_;
      /// Prestressing object
      Teuchos::RCP<Discret::ELEMENTS::PreStress> prestress_;
      /// compute Jacobian mapping wrt to deformed configuration
      void update_jacobian_mapping(
          const std::vector<double>& disp, Discret::ELEMENTS::PreStress& prestress);

      //! Update history variables at the end of time step (fiber direction, inelastic deformation)
      //! (gebauer 07/19)
      void update_element(std::vector<double>& disp,      // current displacements
          Teuchos::ParameterList& params,                 // algorithmic parameters e.g. time
          const Teuchos::RCP<Core::Mat::Material>& mat);  // material

      /// compute defgrd in all gp for given disp
      void def_gradient(const std::vector<double>& disp, Core::LinAlg::SerialDenseMatrix& gpdefgrd,
          Discret::ELEMENTS::PreStress& prestress);


      // internal calculation methods

      // don't want = operator
      SoTet10& operator=(const SoTet10& old);

      //! init the inverse of the jacobian and its determinant in the material configuration
      virtual void init_jacobian_mapping();

      //! Calculate nonlinear stiffness and mass matrix
      virtual void so_tet10_nlnstiffmass(std::vector<int>& lm,  ///< location matrix
          std::vector<double>& disp,                            ///< current displacements
          std::vector<double>* vel,                             ///< current velocities
          std::vector<double>* acc,                             ///< current accelerations
          std::vector<double>& residual,                        ///< current residual displ
          std::vector<double>& dispmat,                         ///< current material displacements
          Core::LinAlg::Matrix<NUMDOF_SOTET10, NUMDOF_SOTET10>*
              stiffmatrix,  ///< element stiffness matrix
          Core::LinAlg::Matrix<NUMDOF_SOTET10, NUMDOF_SOTET10>*
              massmatrix,                                       ///< element mass matrix
          Core::LinAlg::Matrix<NUMDOF_SOTET10, 1>* force,       ///< element internal force vector
          Core::LinAlg::Matrix<NUMDOF_SOTET10, 1>* forceinert,  ///< element inertial force vector
          Core::LinAlg::Matrix<NUMDOF_SOTET10, 1>* force_str,   ///< element structural force vector
          Core::LinAlg::Matrix<NUMGPT_SOTET10, Mat::NUM_STRESS_3D>* elestress,  ///< stresses at GP
          Core::LinAlg::Matrix<NUMGPT_SOTET10, Mat::NUM_STRESS_3D>* elestrain,  ///< strains at GP
          Teuchos::ParameterList& params,          ///< algorithmic parameters e.g. time
          const Inpar::STR::StressType iostress,   ///< stress output option
          const Inpar::STR::StrainType iostrain);  ///< strain output option

      //! lump mass matrix (bborn 07/08)
      void so_tet10_lumpmass(Core::LinAlg::Matrix<NUMDOF_SOTET10, NUMDOF_SOTET10>* emass);

      // These functions are basically copied from the So_integrator,
      // along with the rather inconsitent matrix sizes.

      //! Evaluate Tet10 Shapefcts at 4 gausspoints to keep them static
      std::vector<Core::LinAlg::Matrix<NUMNOD_SOTET10, 1>> so_tet10_4gp_shapefcts();
      //! Evaluate Tet10 Derivs at 4 gausspoints to keep them static
      const std::vector<Core::LinAlg::Matrix<NUMDIM_SOTET10, NUMNOD_SOTET10>>&
      so_tet10_4gp_derivs();

      /*!
       * \brief Evaluate the first derivative of the shape functions at the Gauss point with the
       * 4 Gauss-point integration rule
       *
       * \param derivs First derivatives of the shape functions
       * \param gp Gauss point
       */
      template <Core::FE::GaussRule3D intrule>
      void so_tet10_derivs(
          Core::LinAlg::Matrix<NUMDIM_SOTET10, NUMNOD_SOTET10>& derivs, int gp) const;

      //! Evaluate Tet10 Weights at 4 gausspoints to keep them static
      const std::vector<double>& so_tet10_4gp_weights();

      //! Evaluate Tet10 Shapefcts at 10 gausspoints to keep them static
      const std::vector<Core::LinAlg::Matrix<NUMNOD_SOTET10, 1>>& so_tet10_11gp_shapefcts();
      //! Evaluate Tet10 Derivs at 10 gausspoints to keep them static
      const std::vector<Core::LinAlg::Matrix<NUMDIM_SOTET10, NUMNOD_SOTET10>>&
      so_tet10_11gp_derivs();
      //! Evaluate Tet10 Weights at 10 gausspoints to keep them static
      const std::vector<double>& so_tet10_11gp_weights();

      /*!
       * Executes the post setup call for all materials. This method will be called once per element
       * at the first Evaluate call.
       *
       * @param params ParameterList to be passed to the materials
       */
      void material_post_setup(Teuchos::ParameterList& params) override;

     private:
      std::string get_element_type_string() const { return "SOLIDT10"; }
    };  // class So_tet10



  }  // namespace ELEMENTS
}  // namespace Discret


FOUR_C_NAMESPACE_CLOSE

#endif
