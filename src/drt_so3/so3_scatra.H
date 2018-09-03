/*!----------------------------------------------------------------------
\file so3_scatra.H

\brief Solid-scatra elements base class

\level 2

<pre>
   \maintainer Thon Moritz
               thon@mhpc.mw.tum.de
               089 - 289-10264
</pre>

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

      /*!
      \brief Destructor

      */
      virtual ~So3_Scatra(){};

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
      DRT::Element* Clone() const;

      /*!
      \brief Return unique ParObject id

      every class implementing ParObject needs a unique id defined at the
      top of this file.
      */
      // virtual int UniqueParObjectId() const = 0 ;
      int UniqueParObjectId() const
      {
        switch (distype)
        {
          case DRT::Element::tet4:
          {
            // const DRT::ELEMENTS::So_tet4* ele = dynamic_cast<const
            // DRT::ELEMENTS::So_tet4*>(this); if(ele != NULL)
            return So_tet4ScatraType::Instance().UniqueParObjectId();
            // else
            // return NStet5ScatraType::Instance().UniqueParObjectId();
            break;
          }
          case DRT::Element::tet10:
            return So_tet10ScatraType::Instance().UniqueParObjectId();
            break;
          case DRT::Element::hex8:
          {
            // cast the most specialised element
            // otherwise cast fails, because hex8fbar == hex8
            const DRT::ELEMENTS::So_hex8fbar* ele =
                dynamic_cast<const DRT::ELEMENTS::So_hex8fbar*>(this);

            if (ele != NULL)
              return So_hex8fbarScatraType::Instance().UniqueParObjectId();
            else
              return So_hex8ScatraType::Instance().UniqueParObjectId();

            break;
          }
          case DRT::Element::hex27:
            return So_hex27ScatraType::Instance().UniqueParObjectId();
            break;
          case DRT::Element::wedge6:
            return So_weg6ScatraType::Instance().UniqueParObjectId();
            break;
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
      void Pack(DRT::PackBuffer& data) const;

      /*!
      \brief Unpack data from a char vector into this class

      \ref Pack and \ref Unpack are used to communicate this element

      */
      void Unpack(const std::vector<char>& data);

      //@}

      //! @name Access methods

      /*!
      \brief Print this element
      */
      void Print(std::ostream& os) const;

      DRT::ElementType& ElementType() const
      {
        switch (distype)
        {
          case DRT::Element::tet4:
          {
            const DRT::ELEMENTS::So_tet4* ele = dynamic_cast<const DRT::ELEMENTS::So_tet4*>(this);
            if (ele != NULL) return So_tet4ScatraType::Instance();
            // else
            //  return NStet5ScatraType::Instance();
            break;
          }
          break;
          case DRT::Element::tet10:
            return So_tet10ScatraType::Instance();
            break;
          case DRT::Element::hex8:
          {
            // cast the most specialised element
            // otherwise cast fails, because hex8fbar == hex8
            const DRT::ELEMENTS::So_hex8fbar* ele =
                dynamic_cast<const DRT::ELEMENTS::So_hex8fbar*>(this);

            if (ele != NULL)
              return So_hex8fbarScatraType::Instance();
            else
              return So_hex8ScatraType::Instance();

            break;
          }
            return so3_ele::ElementType();
            break;
          case DRT::Element::hex27:
            return So_hex27ScatraType::Instance();
            break;
          case DRT::Element::wedge6:
            return So_weg6ScatraType::Instance();
            break;
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
          DRT::INPUT::LineDefinition* linedef);

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
      );

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
      void SetMaterial(int matnum);

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

      //! calculate nonlinear B-operator (6x24)
      void CalculateBop(
          LINALG::Matrix<numstr_, numdofperelement_>* bop,  //!< (o): nonlinear B-operator
          const LINALG::Matrix<numdim_, numdim_>* defgrd,   //!< (i): deformation gradient
          const LINALG::Matrix<numdim_, numnod_>*
              N_XYZ  //!< (i): (material) derivative of shape functions
          ) const;

      //@}

     protected:
      //! don't want = operator
      So3_Scatra& operator=(const So3_Scatra& old);

     private:
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

      //! container containing all kind of solid3 thisandthat
      // DRT::Container                                                        data_;

      DRT::Node** Nodes();

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
