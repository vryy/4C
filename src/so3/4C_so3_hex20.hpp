/*----------------------------------------------------------------------*/
/*! \file
\brief 3D quadratic serendipity element
\level 1

*----------------------------------------------------------------------*/
#ifndef FOUR_C_SO3_HEX20_HPP
#define FOUR_C_SO3_HEX20_HPP

#include "4C_config.hpp"

#include "4C_comm_parobjectfactory.hpp"
#include "4C_inpar_structure.hpp"
#include "4C_lib_element.hpp"
#include "4C_lib_elementtype.hpp"
#include "4C_linalg_serialdensematrix.hpp"
#include "4C_mat_material_factory.hpp"
#include "4C_material_base.hpp"
#include "4C_so3_base.hpp"

FOUR_C_NAMESPACE_OPEN

/// Several parameters which are fixed for Solid Hex20
const int NUMNOD_SOH20 = 20;  ///< number of nodes
const int NODDOF_SOH20 = 3;   ///< number of dofs per node
const int NUMDOF_SOH20 = 60;  ///< total dofs per element
const int NUMGPT_SOH20 = 27;  ///< total gauss points per element
const int NUMDIM_SOH20 = 3;   ///< number of dimensions


namespace DRT
{
  // forward declarations
  class Discretization;

  namespace ELEMENTS
  {
    // forward declarations
    class PreStress;

    class SoHex20Type : public DRT::ElementType
    {
     public:
      std::string Name() const override { return "So_hex20Type"; }

      static SoHex20Type& Instance();

      CORE::COMM::ParObject* Create(const std::vector<char>& data) override;

      Teuchos::RCP<DRT::Element> Create(const std::string eletype, const std::string eledistype,
          const int id, const int owner) override;

      Teuchos::RCP<DRT::Element> Create(const int id, const int owner) override;

      int Initialize(DRT::Discretization& dis) override;

      void NodalBlockInformation(
          DRT::Element* dwele, int& numdf, int& dimns, int& nv, int& np) override;

      CORE::LINALG::SerialDenseMatrix ComputeNullSpace(
          DRT::Node& node, const double* x0, const int numdof, const int dimnsp) override;

      void SetupElementDefinition(
          std::map<std::string, std::map<std::string, INPUT::LineDefinition>>& definitions)
          override;

     private:
      static SoHex20Type instance_;

      std::string GetElementTypeString() const { return "SOLIDH20"; }
    };

    /*!
    \brief A C++ version of a 20-node hex solid element

    A structural 20-node hexahedral solid displacement element for large deformations.
    As its discretization is fixed many data structures are evaluated just once and kept
    for performance. It heavily uses Epetra objects and methods and therefore relies
    on their performance.

    \author kloeppel
    */
    class SoHex20 : public SoBase
    {
     public:
      //! @name Friends
      friend class SoHex20Type;

      //@}
      //! @name Constructors and destructors and related methods

      /*!
      \brief Standard Constructor

      \param id : A unique global id
      \param owner : elements owner
      */
      SoHex20(int id, int owner);

      /*!
      \brief Copy Constructor

      Makes a deep copy of a Element

      */
      SoHex20(const SoHex20& old);

      /*!
      \brief Deep copy this instance of Solid3 and return pointer to the copy

      The Clone() method is used from the virtual base class Element in cases
      where the type of the derived class is unknown and a copy-ctor is needed

      */
      DRT::Element* Clone() const override;

      /*!
      \brief Get shape type of element
      */
      CORE::FE::CellType Shape() const override;

      /*!
      \brief Return number of volumes of this element
      */
      int NumVolume() const override { return 1; }

      /*!
      \brief Return number of surfaces of this element
      */
      int NumSurface() const override { return 6; }

      /*!
      \brief Return number of lines of this element
      */
      int NumLine() const override { return 12; }

      /*!
      \brief Get vector of Teuchos::RCPs to the lines of this element

      */
      std::vector<Teuchos::RCP<DRT::Element>> Lines() override;

      /*!
      \brief Get vector of Teuchos::RCPs to the surfaces of this element

      */
      std::vector<Teuchos::RCP<DRT::Element>> Surfaces() override;

      /*!
      \brief Return unique ParObject id

      every class implementing ParObject needs a unique id defined at the
      top of this file.
      */
      int UniqueParObjectId() const override { return SoHex20Type::Instance().UniqueParObjectId(); }

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


      //@}

      //! @name Acess methods


      /*!
      \brief Get number of degrees of freedom of a certain node
             (implements pure virtual DRT::Element)

      The element decides how many degrees of freedom its nodes must have.
      As this may vary along a simulation, the element can redecide the
      number of degrees of freedom per node along the way for each of it's nodes
      separately.
      */
      int NumDofPerNode(const DRT::Node& node) const override { return 3; }

      /*!
      \brief Get number of degrees of freedom per element
             (implements pure virtual DRT::Element)

      The element decides how many element degrees of freedom it has.
      It can redecide along the way of a simulation.

      \note Element degrees of freedom mentioned here are dofs that are visible
            at the level of the total system of equations. Purely internal
            element dofs that are condensed internally should NOT be considered.
      */
      int NumDofPerElement() const override { return 0; }

      /*!
      \brief Print this element
      */
      void Print(std::ostream& os) const override;

      DRT::ElementType& ElementType() const override { return SoHex20Type::Instance(); }

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
      \endcode

      */
      void VisNames(std::map<std::string, int>&
              names  ///< to be filled with key names of data to visualize and with int dimensions
          ) override;

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
          INPUT::LineDefinition* linedef) override;

      //@}

      //! @name Evaluation

      /*!
      \brief Evaluate an element

      Evaluate So_hex20 element stiffness, mass, internal forces, etc.

      If nullptr on input, the controling method does not expect the element
      to fill these matrices or vectors.

      \return 0 if successful, negative otherwise
      */
      int Evaluate(
          Teuchos::ParameterList&
              params,  ///< ParameterList for communication between control routine and elements
          DRT::Discretization& discretization,  ///< pointer to discretization for de-assembly
          std::vector<int>& lm,                 ///< location matrix for de-assembly
          CORE::LINALG::SerialDenseMatrix&
              elemat1,  ///< (stiffness-)matrix to be filled by element.
          CORE::LINALG::SerialDenseMatrix& elemat2,  ///< (mass-)matrix to be filled by element.
          CORE::LINALG::SerialDenseVector&
              elevec1,  ///< (internal force-)vector to be filled by element
          CORE::LINALG::SerialDenseVector& elevec2,  ///< vector to be filled by element
          CORE::LINALG::SerialDenseVector& elevec3   ///< vector to be filled by element
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
      int EvaluateNeumann(Teuchos::ParameterList& params, DRT::Discretization& discretization,
          CORE::Conditions::Condition& condition, std::vector<int>& lm,
          CORE::LINALG::SerialDenseVector& elevec1,
          CORE::LINALG::SerialDenseMatrix* elemat1 = nullptr) override;


      // const vector<double> GetFibervec(){return fiberdirection_;};
      std::vector<double> soh20_ElementCenterRefeCoords();

      //@}

     protected:
      //! action parameters recognized by So_hex20
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
        prestress_update,
        calc_struct_energy,  //!< compute internal energy
        multi_readrestart,   //!< multi-scale: read restart on microscale
        multi_calc_dens      //!< multi-scale: calculate homogenized density
      };

      //! vector of inverses of the jacobian in material frame
      std::vector<CORE::LINALG::Matrix<NUMDIM_SOH20, NUMDIM_SOH20>> invJ_;
      //! determinant of Jacobian in material frame
      std::vector<double> detJ_;

      /// prestressing switch & time
      INPAR::STR::PreStress pstype_;
      double pstime_;

      double time_;
      /// Prestressing object
      Teuchos::RCP<DRT::ELEMENTS::PreStress> prestress_;
      /// compute Jacobian mapping wrt to deformed configuration
      void UpdateJacobianMapping(
          const std::vector<double>& disp, DRT::ELEMENTS::PreStress& prestress);
      /// compute defgrd in all gp for given disp
      void DefGradient(const std::vector<double>& disp, CORE::LINALG::SerialDenseMatrix& gpdefgrd,
          DRT::ELEMENTS::PreStress& prestress);


      // internal calculation methods

      //! don't want = operator
      SoHex20& operator=(const SoHex20& old);


      //! init the inverse of the jacobian and its determinant in the material configuration
      virtual void InitJacobianMapping();

      //! Calculate linear stiffness and mass matrix
      virtual void soh20_linstiffmass(std::vector<int>& lm,  ///< location matrix
          std::vector<double>& disp,                         ///< current displacements
          std::vector<double>& residual,                     ///< current residual displ
          CORE::LINALG::Matrix<NUMDOF_SOH20, NUMDOF_SOH20>*
              stiffmatrix,  ///< element stiffness matrix
          CORE::LINALG::Matrix<NUMDOF_SOH20, NUMDOF_SOH20>* massmatrix,  ///< element mass matrix
          CORE::LINALG::Matrix<NUMDOF_SOH20, 1>* force,  ///< element internal force vector
          CORE::LINALG::Matrix<NUMGPT_SOH20, MAT::NUM_STRESS_3D>* elestress,  ///< stresses at GP
          CORE::LINALG::Matrix<NUMGPT_SOH20, MAT::NUM_STRESS_3D>* elestrain,  ///< strains at GP
          Teuchos::ParameterList& params,          ///< algorithmic parameters e.g. time
          const INPAR::STR::StressType iostress,   ///< stress output option
          const INPAR::STR::StrainType iostrain);  ///< strain output option

      //! Calculate nonlinear stiffness and mass matrix
      virtual void soh20_nlnstiffmass(std::vector<int>& lm,  ///< location matrix
          std::vector<double>& disp,                         ///< current displacements
          std::vector<double>* vel,                          ///< current velocities
          std::vector<double>* acc,                          ///< current accelerations
          std::vector<double>& residual,                     ///< current residual displ
          std::vector<double>& dispmat,                      ///< current material displacements
          CORE::LINALG::Matrix<NUMDOF_SOH20, NUMDOF_SOH20>*
              stiffmatrix,  ///< element stiffness matrix
          CORE::LINALG::Matrix<NUMDOF_SOH20, NUMDOF_SOH20>* massmatrix,  ///< element mass matrix
          CORE::LINALG::Matrix<NUMDOF_SOH20, 1>* force,       ///< element internal force vector
          CORE::LINALG::Matrix<NUMDOF_SOH20, 1>* forceinert,  ///< element inertial force vector
          CORE::LINALG::Matrix<NUMDOF_SOH20, 1>* force_str,   ///< element structural force vector
          CORE::LINALG::Matrix<NUMGPT_SOH20, MAT::NUM_STRESS_3D>* elestress,  ///< stresses at GP
          CORE::LINALG::Matrix<NUMGPT_SOH20, MAT::NUM_STRESS_3D>* elestrain,  ///< strains at GP
          Teuchos::ParameterList& params,          ///< algorithmic parameters e.g. time
          const INPAR::STR::StressType iostress,   ///< stress output option
          const INPAR::STR::StrainType iostrain);  ///< strain output option

      //! Lump mass matrix (bborn 07/08)
      void soh20_lumpmass(CORE::LINALG::Matrix<NUMDOF_SOH20, NUMDOF_SOH20>* emass);

      //! Evaluate Hex20 Shapefcts to keep them static
      std::vector<CORE::LINALG::Matrix<NUMNOD_SOH20, 1>> soh20_shapefcts();
      //! Evaluate Hex20 Derivs to keep them static
      std::vector<CORE::LINALG::Matrix<NUMDIM_SOH20, NUMNOD_SOH20>> soh20_derivs();
      //! Evaluate Hex20 Weights to keep them static
      std::vector<double> soh20_weights();

      //! Evaluate shapefunction, derivative and gaussweights
      void soh20_shapederiv(CORE::LINALG::Matrix<NUMNOD_SOH20, NUMGPT_SOH20>**
                                shapefct,  // pointer to pointer of shapefct
          CORE::LINALG::Matrix<NUMDOF_SOH20, NUMNOD_SOH20>** deriv,  // pointer to pointer of derivs
          CORE::LINALG::Matrix<NUMGPT_SOH20, 1>** weights);  // pointer to pointer of weights

      //   /*!
      //    * \brief Interface to material laws working with soh8.
      //    *
      //    * Here the interface to any material takes place. Stress and C-mat must be retrieved
      //    * whereas all necessary data is passed. Add whatever your material needs, but
      //    * make sure that exchange is not overdone performance-wise.
      //    * Right now EAS is based on GL-Strains which is ok for any hyperelastic material.
      //    * This means that e.g. the deformation gradient is not alleviated from locking
      //    * and therefore only to be used with disp-based soh8.
      //    *
      //    * \param *stress (out): Voigt-Vector of stresses at current gp
      //    * \param *cmat (out): Elasticity matrix at current gp
      //    * \param *density (out): density of material
      //    * \param *glstrain (in): Voigt-Vector of GL-strains at current gp
      //    * \param *defgrd (in): F at current gp, CAUTION! only for disp-based soh8
      //    * \param gp (in): current gp
      //    * \param params (in): parameterlist to access time, etc. in materials
      //    * */

      //   void soh20_mat_sel(CORE::LINALG::Matrix<6,1>* stress,
      //                     CORE::LINALG::Matrix<6,6>* cmat,
      //                     double* density,
      //                     CORE::LINALG::Matrix<6,1>* glstrain,
      //                     CORE::LINALG::Matrix<3,3>* defgrd,
      //                     const int gp,
      //                     Teuchos::ParameterList&  params);

      //! @name Multi-scale related stuff

      /*!
       * \brief Determine a homogenized material density for multi-scale
       * analyses by averaging over the initial volume
       * */
      void soh20_homog(Teuchos::ParameterList& params);

      /*!
       * \brief Read restart on the microscale
       * */
      void soh20_read_restart_multi();

      //@}

      /// temporary method for compatibility with solidshell, needs clarification
      std::vector<double> getthicknessvector() const
      {
        FOUR_C_THROW("not implemented");
        return std::vector<double>(3);
      };

     private:
      std::string GetElementTypeString() const { return "SOLIDH20"; }
    };  // class So_hex20



    //=======================================================================
    //=======================================================================
    //=======================================================================
    //=======================================================================



  }  // namespace ELEMENTS
}  // namespace DRT

FOUR_C_NAMESPACE_CLOSE

#endif
