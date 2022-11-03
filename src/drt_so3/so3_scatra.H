/*----------------------------------------------------------------------*/
/*! \file

\brief Solid-scatra elements base class

\level 2


*----------------------------------------------------------------------*/
#ifndef SO3_SCATRA_H
#define SO3_SCATRA_H


#include "so3_scatra_eletypes.H"

#include "../drt_fem_general/drt_utils_gausspoints.H"
#include "../drt_fem_general/drt_utils_integration.H"

#include "../drt_inpar/inpar_scatra.H"

namespace DRT
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
    template <class so3_ele, DRT::Element::DiscretizationType distype>
    class So3_Scatra : public so3_ele
    {
      //! @name Friends
      friend class So_tet4ScatraType;
      friend class So_tet10ScatraType;
      friend class So_hex8ScatraType;
      friend class So_hex8fbarScatraType;
      friend class So_hex27ScatraType;
      friend class So_weg6ScatraType;
      // friend class NStet5ScatraType;

     public:
      //@}
      //! @name Constructors and destructors and related methods

      /*!
      \brief Standard Constructor

      \param id : A unique global id
      \param owner : elements owner
      */
      So3_Scatra(int id, int owner);

      /*!
      \brief Copy Constructor

      Makes a deep copy of a Element

      */
      So3_Scatra(const So3_Scatra& old);

      //! don't want = operator
      So3_Scatra& operator=(const So3_Scatra& old) = delete;
      //@}

      // static const DRT::UTILS::GaussIntegration intpoints_ =
      // DRT::UTILS::GaussIntegration(distype);


      //! number of element nodes (
      static const int numnod_ = DRT::UTILS::DisTypeToNumNodePerEle<distype>::numNodePerElement;
      //! number of space dimensions
      static const int numdim_ = 3;
      //! number of dofs per node
      static const int numdofpernode_ = 3;
      //! total dofs per element
      static const int numdofperelement_ = numdofpernode_ * numnod_;
      //! number of strains/stresses
      static const int numstr_ = 6;

      //! @name Acess methods

      /*!
      \brief Deep copy this instance of Solid3 and return pointer to the copy

      The Clone() method is used from the virtual base class Element in cases
      where the type of the derived class is unknown and a copy-ctor is needed

      */
      DRT::Element* Clone() const override;

      /*!
      \brief Return unique ParObject id

      every class implementing ParObject needs a unique id defined at the
      top of this file.
      */
      int UniqueParObjectId() const override
      {
        switch (distype)
        {
          case DRT::Element::tet4:
          {
            return So_tet4ScatraType::Instance().UniqueParObjectId();
          }
          case DRT::Element::tet10:
            return So_tet10ScatraType::Instance().UniqueParObjectId();
          case DRT::Element::hex8:
          {
            // cast the most specialised element
            // otherwise cast fails, because hex8fbar == hex8
            const auto* ele = dynamic_cast<const DRT::ELEMENTS::So_hex8fbar*>(this);

            if (ele != nullptr)
              return So_hex8fbarScatraType::Instance().UniqueParObjectId();
            else
              return So_hex8ScatraType::Instance().UniqueParObjectId();
          }
          case DRT::Element::hex27:
            return So_hex27ScatraType::Instance().UniqueParObjectId();
          case DRT::Element::wedge6:
            return So_weg6ScatraType::Instance().UniqueParObjectId();
          default:
            dserror("unknown element type!");
            break;
        }
        // Intel compiler needs a return so
        return -1;
      };

      /*!
      \brief Pack this class so it can be communicated

      \ref Pack and \ref Unpack are used to communicate this element

      */
      void Pack(DRT::PackBuffer& data) const override;

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
      DRT::ElementType& ElementType() const override
      {
        switch (distype)
        {
          case DRT::Element::tet4:
          {
            const auto* ele = dynamic_cast<const DRT::ELEMENTS::So_tet4*>(this);
            if (ele != nullptr) return So_tet4ScatraType::Instance();
            break;
          }
          case DRT::Element::tet10:
            return So_tet10ScatraType::Instance();
          case DRT::Element::hex8:
          {
            // cast the most specialised element
            // otherwise cast fails, because hex8fbar == hex8
            const auto* ele = dynamic_cast<const DRT::ELEMENTS::So_hex8fbar*>(this);

            if (ele != nullptr)
              return So_hex8fbarScatraType::Instance();
            else
              return So_hex8ScatraType::Instance();
          }
          case DRT::Element::hex27:
            return So_hex27ScatraType::Instance();
          case DRT::Element::wedge6:
            return So_weg6ScatraType::Instance();
          default:
            dserror("unknown element type!");
            break;
        }
        // Intel compiler needs a return so
        return So_weg6ScatraType::Instance();
      };


      //@}

      //! @name Input and Creation

      /*!
      \brief Read input for this element
      */
      bool ReadElement(const std::string& eletype, const std::string& eledistype,
          DRT::INPUT::LineDefinition* linedef) override;

      //! initialize the elements
      void InitElement();

      //@}

      //! @name Evaluation

      void PreEvaluate(
          Teuchos::ParameterList&
              params,  ///< ParameterList for communication between control routine and elements
          DRT::Discretization& discretization,  ///< pointer to discretization for de-assembly
          DRT::Element::LocationArray& la       ///< location array for de-assembly
      );

      /*!
      \brief Evaluate an element

      Evaluate So3_poro element stiffness, mass, internal forces, etc.
      Templated evaluate routine of element matrixes

      If NULL on input, the controlling method does not expect the element
      to fill these matrices or vectors.

      \return 0 if successful, negative otherwise
      */
      int Evaluate(
          Teuchos::ParameterList&
              params,  ///< ParameterList for communication between control routine and elements
          DRT::Discretization& discretization,  ///< pointer to discretization for de-assembly
          DRT::Element::LocationArray& la,      ///< location array for de-assembly
          Epetra_SerialDenseMatrix& elemat1,    ///< (stiffness-)matrix to be filled by element.
          Epetra_SerialDenseMatrix& elemat2,    ///< (mass-)matrix to be filled by element.
          Epetra_SerialDenseVector& elevec1,    ///< (internal force-)vector to be filled by element
          Epetra_SerialDenseVector& elevec2,    ///< vector to be filled by element
          Epetra_SerialDenseVector& elevec3     ///< vector to be filled by element
          ) override;

      //! init the inverse of the jacobian and its determinant in the material configuration
      // virtual void InitJacobianMapping();


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
      void SetMaterial(int matnum) override;

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
      void GetCauchyNDirAndDerivativesAtXi(const LINALG::Matrix<3, 1>& xi,
          const std::vector<double>& disp_nodal_values,
          const std::vector<double>& scalar_nodal_values, const LINALG::Matrix<3, 1>& n,
          const LINALG::Matrix<3, 1>& dir, double& cauchy_n_dir,
          Epetra_SerialDenseMatrix* d_cauchyndir_dd, Epetra_SerialDenseMatrix* d_cauchyndir_ds,
          LINALG::Matrix<3, 1>* d_cauchyndir_dn, LINALG::Matrix<3, 1>* d_cauchyndir_ddir,
          LINALG::Matrix<3, 1>* d_cauchyndir_dxi);

      /// @name params
      /// return SCATRA::ImplType
      const INPAR::SCATRA::ImplType& ImplType() const { return impltype_; };

     private:
      //! scalar transport implementation type (physics)
      INPAR::SCATRA::ImplType impltype_;

      //! @name Evaluation
      //! Calculate mechanical scalar stiffness term needed for monolithic SSI K_dS
      void nln_kdS_ssi(DRT::Element::LocationArray& la,  //!< (i): location array
          std::vector<double>& disp,                     //!< (i): current displacement
          Epetra_SerialDenseMatrix&
              stiffmatrix_kdS,            //!< (o): mechanical scalar stiffness term at current gp
          Teuchos::ParameterList& params  //!< parameter list
      );

      /*!
       * @brief calculate nonlinear B-operator (6x24)
       *
       * @param[out] bop    nonlinear B-operator
       * @param[in] defgrd  deformation gradient
       * @param[in] N_XYZ   (material) derivative of shape functions
       */
      void CalculateBop(LINALG::Matrix<numstr_, numdofperelement_>* bop,
          const LINALG::Matrix<numdim_, numdim_>* defgrd,
          const LINALG::Matrix<numdim_, numnod_>* N_XYZ) const;
      //@}

      // action type recognized by So3_Scatra elements
      enum ActionType
      {
        none,
        calc_struct_stiffscalar  //!< calculate coupling term k_dS for monolithic SSI
      };

      // integration points
      const DRT::UTILS::IntegrationPoints3D intpoints_;
      //! total gauss points per element
      const int numgpt_;
      //! vector of coordinates of current integration point in reference coordinates
      std::vector<LINALG::Matrix<numdim_, 1>> xsi_;

      //! vector of inverses of the jacobian (J^{-1} = \frac{\mathrm{d} \vec{r}} {\mathrm{d} \vec{X}
      //! }) at each gauss point
      std::vector<LINALG::Matrix<numdim_, numdim_>> invJ_;
      //! vector of determinants of the jacobian (\det[ \frac{\mathrm{d} \vec{X}} {\mathrm{d}
      //! \vec{r}} ]) at each gauss point
      std::vector<double> detJ_;

      DRT::Node** Nodes() override;

      Teuchos::RCP<MAT::Material> Material() const;

      int Id() const;

    };  // class So3_Scatra


    //=======================================================================
    //=======================================================================
    //=======================================================================
    //=======================================================================

  }  // namespace ELEMENTS
}  // namespace DRT

#endif /* SO3_SCATRA_H */
