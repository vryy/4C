/*----------------------------------------------------------------------*/
/*! \file

\brief Solid-scatra elements base class

\level 2


*----------------------------------------------------------------------*/
#ifndef FOUR_C_SO3_SCATRA_HPP
#define FOUR_C_SO3_SCATRA_HPP


#include "4C_config.hpp"

#include "4C_discretization_fem_general_utils_gausspoints.hpp"
#include "4C_discretization_fem_general_utils_integration.hpp"
#include "4C_inpar_scatra.hpp"
#include "4C_so3_scatra_eletypes.hpp"

FOUR_C_NAMESPACE_OPEN

namespace Discret
{
  // forward declarations
  class Discretization;

  namespace ELEMENTS
  {
    /*!
    \brief A C++ version of a 3 dimensional solid element with modifications for volume coupling
    with scatra

    A structural 3 dimensional solid displacement element for large deformations
    and (near)-incompressibility.

    */
    template <class so3_ele, Core::FE::CellType distype>
    class So3Scatra : public so3_ele
    {
      //! @name Friends
      friend class SoTet4ScatraType;
      friend class SoTet10ScatraType;
      friend class SoHex8ScatraType;
      friend class SoHex8fbarScatraType;
      friend class SoHex27ScatraType;
      friend class SoWeg6ScatraType;
      // friend class NStet5ScatraType;

     public:
      //@}
      //! @name Constructors and destructors and related methods

      /*!
      \brief Standard Constructor

      \param id : A unique global id
      \param owner : elements owner
      */
      So3Scatra(int id, int owner);

      /*!
      \brief Copy Constructor

      Makes a deep copy of a Element

      */
      So3Scatra(const So3Scatra& old);

      //! don't want = operator
      So3Scatra& operator=(const So3Scatra& old) = delete;
      //@}

      // static constexpr Core::FE::GaussIntegration intpoints_ =
      // Core::FE::GaussIntegration(distype);


      //! number of element nodes (
      static constexpr int numnod_ = Core::FE::num_nodes<distype>;
      //! number of space dimensions
      static constexpr int numdim_ = 3;
      //! number of dofs per node
      static constexpr int numdofpernode_ = 3;
      //! total dofs per element
      static constexpr int numdofperelement_ = numdofpernode_ * numnod_;
      //! number of strains/stresses
      static constexpr int numstr_ = 6;

      //! @name Acess methods

      /*!
      \brief Deep copy this instance of Solid3 and return pointer to the copy

      The Clone() method is used from the virtual base class Element in cases
      where the type of the derived class is unknown and a copy-ctor is needed

      */
      Core::Elements::Element* Clone() const override;

      /*!
      \brief Return unique ParObject id

      every class implementing ParObject needs a unique id defined at the
      top of this file.
      */
      int UniqueParObjectId() const override
      {
        switch (distype)
        {
          case Core::FE::CellType::tet4:
          {
            return SoTet4ScatraType::Instance().UniqueParObjectId();
          }
          case Core::FE::CellType::tet10:
            return SoTet10ScatraType::Instance().UniqueParObjectId();
          case Core::FE::CellType::hex8:
          {
            // cast the most specialised element
            // otherwise cast fails, because hex8fbar == hex8
            const auto* ele = dynamic_cast<const Discret::ELEMENTS::SoHex8fbar*>(this);

            if (ele != nullptr)
              return SoHex8fbarScatraType::Instance().UniqueParObjectId();
            else
              return SoHex8ScatraType::Instance().UniqueParObjectId();
          }
          case Core::FE::CellType::hex27:
            return SoHex27ScatraType::Instance().UniqueParObjectId();
          case Core::FE::CellType::wedge6:
            return SoWeg6ScatraType::Instance().UniqueParObjectId();
          default:
            FOUR_C_THROW("unknown element type!");
            break;
        }
        // Intel compiler needs a return so
        return -1;
      };

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

      //! @name Access methods

      /*!
      \brief Print this element
      */
      void Print(std::ostream& os) const override;

      //! returns a reference to the element type
      Core::Elements::ElementType& ElementType() const override
      {
        switch (distype)
        {
          case Core::FE::CellType::tet4:
          {
            const auto* ele = dynamic_cast<const Discret::ELEMENTS::SoTet4*>(this);
            if (ele != nullptr) return SoTet4ScatraType::Instance();
            break;
          }
          case Core::FE::CellType::tet10:
            return SoTet10ScatraType::Instance();
          case Core::FE::CellType::hex8:
          {
            // cast the most specialised element
            // otherwise cast fails, because hex8fbar == hex8
            const auto* ele = dynamic_cast<const Discret::ELEMENTS::SoHex8fbar*>(this);

            if (ele != nullptr)
              return SoHex8fbarScatraType::Instance();
            else
              return SoHex8ScatraType::Instance();
          }
          case Core::FE::CellType::hex27:
            return SoHex27ScatraType::Instance();
          case Core::FE::CellType::wedge6:
            return SoWeg6ScatraType::Instance();
          default:
            FOUR_C_THROW("unknown element type!");
            break;
        }
        // Intel compiler needs a return so
        return SoWeg6ScatraType::Instance();
      };


      //@}

      //! @name Input and Creation

      /*!
      \brief Read input for this element
      */
      bool ReadElement(const std::string& eletype, const std::string& eledistype,
          Input::LineDefinition* linedef) override;

      //! initialize the elements
      void InitElement();

      //@}

      //! @name Evaluation

      void pre_evaluate(
          Teuchos::ParameterList&
              params,  ///< ParameterList for communication between control routine and elements
          Discret::Discretization& discretization,    ///< pointer to discretization for de-assembly
          Core::Elements::Element::LocationArray& la  ///< location array for de-assembly
      );

      /*!
      \brief Evaluate an element

      Evaluate So3_poro element stiffness, mass, internal forces, etc.
      Templated evaluate routine of element matrixes

      If nullptr on input, the controlling method does not expect the element
      to fill these matrices or vectors.

      \return 0 if successful, negative otherwise
      */
      int Evaluate(
          Teuchos::ParameterList&
              params,  ///< ParameterList for communication between control routine and elements
          Discret::Discretization& discretization,  ///< pointer to discretization for de-assembly
          Core::Elements::Element::LocationArray& la,  ///< location array for de-assembly
          Core::LinAlg::SerialDenseMatrix&
              elemat1,  ///< (stiffness-)matrix to be filled by element.
          Core::LinAlg::SerialDenseMatrix& elemat2,  ///< (mass-)matrix to be filled by element.
          Core::LinAlg::SerialDenseVector&
              elevec1,  ///< (internal force-)vector to be filled by element
          Core::LinAlg::SerialDenseVector& elevec2,  ///< vector to be filled by element
          Core::LinAlg::SerialDenseVector& elevec3   ///< vector to be filled by element
          ) override;

      //! init the inverse of the jacobian and its determinant in the material configuration
      // virtual void init_jacobian_mapping();


      //@}

      /// Set element material
      /*!
        Material numbers are read from the input file. The element stores
        a corresponding material object. These material objects can be
        anything from very simple (just a little calculation) to highly
        sophisticated with history data. The material is packed and
        unpacked along with its element.

        \param matnum : material number from input file
       */
      void SetMaterial(int matnum, Teuchos::RCP<Core::Mat::Material> mat) override;

      /*!
       * @brief Evaluate Cauchy stress at given point in parameter space and calculate
       * linearizations
       *
       * @param[in] xi                   position in parameter space xi
       * @param[in] disp_nodal_values    vector containing nodal values of displacements
       * @param[in] scalar_nodal_values  vector containing nodal values of scalars
       * @param[in] n                    vector (\f[\mathbf{n}\f])
       * @param[in] dir                  direction vector (\f[\mathbf{v}\f]), can be either normal
       * or tangential vector
       * @param[out] cauchy_n_dir  cauchy stress tensor contracted using the vectors n and dir
       *                           (\f[ \mathbf{\sigma} \cdot \mathbf{n} \cdot \mathbf{v} \f])
       * @param[out] d_cauchyndir_dd  derivative of cauchy_n_dir w.r.t. displacements
       *                         (\f[ \frac{ \mathrm{d} \mathbf{\sigma} \cdot \mathbf{n} \cdot
       * \mathbf{v}} { \mathrm{d} \mathbf{d}} \f])
       * @param d_cauchyndir_ds  derivative of cauchy_n_dir w.r.t. vector of scalars s
       *                        (\f[ \frac{ \mathrm{d} \mathbf{\sigma} \cdot \mathbf{n} \cdot
       * \mathbf{v}} { \mathrm{d} \mathbf{s}} \f])
       * @param[out] d_cauchyndir_dn   derivative of cauchy_n_dir w.r.t. vector n
       *                        (\f[ \frac{ \mathrm{d} \mathbf{\sigma} \cdot \mathbf{n} \cdot
       * \mathbf{v}} { \mathrm{d} \mathbf{n}} \f])
       * @param[out] d_cauchyndir_ddir  derivative of cauchy_n_dir w.r.t. direction vector v
       *                        (\f[ \frac{ \mathrm{d} \mathbf{\sigma} \cdot \mathbf{n} \cdot
       * \mathbf{v}} { \mathrm{d} \mathbf{v}} \f])
       * @param[out] d_cauchyndir_dxi  derivative of cauchy_n_dir w.r.t. local parameter coord. xi
       *                        (\f[ \frac{ \mathrm{d} \mathbf{\sigma} \cdot \mathbf{n} \cdot
       * \mathbf{v}} { \mathrm{d} \mathbf{\xi}} \f])
       *
       * @note At the moment this method is only used for the nitsche contact formulation
       */
      void get_cauchy_n_dir_and_derivatives_at_xi(const Core::LinAlg::Matrix<3, 1>& xi,
          const std::vector<double>& disp_nodal_values,
          const std::vector<double>& scalar_nodal_values, const Core::LinAlg::Matrix<3, 1>& n,
          const Core::LinAlg::Matrix<3, 1>& dir, double& cauchy_n_dir,
          Core::LinAlg::SerialDenseMatrix* d_cauchyndir_dd,
          Core::LinAlg::SerialDenseMatrix* d_cauchyndir_ds,
          Core::LinAlg::Matrix<3, 1>* d_cauchyndir_dn,
          Core::LinAlg::Matrix<3, 1>* d_cauchyndir_ddir,
          Core::LinAlg::Matrix<3, 1>* d_cauchyndir_dxi);

      /// @name params
      /// return ScaTra::ImplType
      const Inpar::ScaTra::ImplType& ImplType() const { return impltype_; };

     private:
      //! scalar transport implementation type (physics)
      Inpar::ScaTra::ImplType impltype_;

      //! @name Evaluation
      //! Calculate mechanical scalar stiffness term needed for monolithic SSI K_dS
      void nln_kd_s_ssi(Core::Elements::Element::LocationArray& la,  //!< (i): location array
          std::vector<double>& disp,                                 //!< (i): current displacement
          Core::LinAlg::SerialDenseMatrix&
              stiffmatrix_kdS,            //!< (o): mechanical scalar stiffness term at current gp
          Teuchos::ParameterList& params  //!< parameter list
      );

      /*!
       * @brief calculate nonlinear B-operator (6x24)
       *
       * @param[out] bop    nonlinear B-operator
       * @param[in] defgrad  deformation gradient
       * @param[in] N_XYZ   (material) derivative of shape functions
       */
      void calculate_bop(Core::LinAlg::Matrix<numstr_, numdofperelement_>* bop,
          const Core::LinAlg::Matrix<numdim_, numdim_>* defgrad,
          const Core::LinAlg::Matrix<numdim_, numnod_>* N_XYZ) const;
      //@}

      // action type recognized by So3_Scatra elements
      enum ActionType
      {
        none,
        calc_struct_stiffscalar  //!< calculate coupling term k_dS for monolithic SSI
      };

      // integration points
      const Core::FE::IntegrationPoints3D intpoints_;
      //! total gauss points per element
      const int numgpt_;
      //! vector of coordinates of current integration point in reference coordinates
      std::vector<Core::LinAlg::Matrix<numdim_, 1>> xsi_;

      //! vector of inverses of the jacobian (J^{-1} = \frac{\mathrm{d} \vec{r}} {\mathrm{d} \vec{X}
      //! }) at each gauss point
      std::vector<Core::LinAlg::Matrix<numdim_, numdim_>> inv_j_;
      //! vector of determinants of the jacobian (\det[ \frac{\mathrm{d} \vec{X}} {\mathrm{d}
      //! \vec{r}} ]) at each gauss point
      std::vector<double> det_j_;

      Core::Nodes::Node** Nodes() override;

      Teuchos::RCP<Core::Mat::Material> material() const;

      int id() const;

    };  // class So3_Scatra


    //=======================================================================
    //=======================================================================
    //=======================================================================
    //=======================================================================

  }  // namespace ELEMENTS
}  // namespace Discret

FOUR_C_NAMESPACE_CLOSE

#endif
