/*----------------------------------------------------------------------*/
/*! \file

\brief Solid Hex8 element

\level 1


*/
/*----------------------------------------------------------------------*/
#ifndef FOUR_C_SO3_HEX8_HPP
#define FOUR_C_SO3_HEX8_HPP

#include "baci_config.hpp"

#include "baci_discretization_fem_general_utils_integration.hpp"
#include "baci_inpar_structure.hpp"
#include "baci_lib_discret.hpp"
#include "baci_lib_element.hpp"
#include "baci_lib_element_integration_select.hpp"
#include "baci_lib_elementtype.hpp"
#include "baci_linalg_serialdensematrix.hpp"
#include "baci_mat_material.hpp"
#include "baci_so3_base.hpp"

FOUR_C_NAMESPACE_OPEN

namespace STR
{
  namespace ELEMENTS
  {
    enum EvalErrorFlag : int;
  }  // namespace ELEMENTS
}  // namespace STR

// Several parameters which are fixed for Solid Hex8
const int NUMNOD_SOH8 = 8;   ///< number of nodes
const int NODDOF_SOH8 = 3;   ///< number of dofs per node
const int NUMDOF_SOH8 = 24;  ///< total dofs per element
const int NUMDIM_SOH8 = 3;   ///< number of dimensions

/// Gauss integration rule
struct GpRuleSoH8
{
  static constexpr enum CORE::FE::GaussRule3D rule =
      DRT::ELEMENTS::DisTypeToOptGaussRule<CORE::FE::CellType::hex8>::rule;
};
/// total gauss points per element
const unsigned NUMGPT_SOH8 = CORE::FE::GaussRule3DToNumGaussPoints<GpRuleSoH8::rule>::value;

namespace DRT
{
  // forward declarations
  class Discretization;

  namespace ELEMENTS
  {
    // forward declarations
    class PreStress;
    class SoSh8Type;


    class SoHex8Type : public DRT::ElementType
    {
     public:
      std::string Name() const override { return "So_hex8Type"; }

      static SoHex8Type& Instance();

      CORE::COMM::ParObject* Create(const std::vector<char>& data) override;

      Teuchos::RCP<DRT::Element> Create(
          std::string eletype, std::string eledistype, int id, int owner) override;

      Teuchos::RCP<DRT::Element> Create(int id, int owner) override;

      int Initialize(DRT::Discretization& dis) override;

      void NodalBlockInformation(
          DRT::Element* dwele, int& numdf, int& dimns, int& nv, int& np) override;

      CORE::LINALG::SerialDenseMatrix ComputeNullSpace(
          DRT::Node& node, const double* x0, const int numdof, const int dimnsp) override;

      void SetupElementDefinition(
          std::map<std::string, std::map<std::string, INPUT::LineDefinition>>& definitions)
          override;

     private:
      static SoHex8Type instance_;

      std::string GetElementTypeString() const { return "SOLIDH8"; }
    };

    /*!
    \brief A C++ version of the 8-node hex solid element

    A structural 8-node hexahedral solid element for large deformations. As its
    discretization is fixed many data structures are evaluated just once and kept
    for performance (e.g. shape functions, derivatives, etc.,
    see DRT::ELEMENTS::So_hex8::Integrator_So_hex8). It heavily uses
    Epetra objects and methods and therefore relies on their performance.

    There are 2 sets of EAS enhancements for GL-strains to alleviate locking
    (see DRT::ELEMENTS::So_hex8::EASType).

    */
    class SoHex8 : public SoBase
    {
     public:
      //! @name Friends
      friend class SoHex8Type;
      friend class SoSh8Type;

      //@}
      //! @name Constructors and destructors and related methods

      /*!
      \brief Standard Constructor

      \param id : A unique global id
      \param owner : elements owner
      */
      SoHex8(int id, int owner);

      /*!
      \brief Copy Constructor

      Makes a deep copy of a Element

      */
      SoHex8(const SoHex8& old);

      //! don't want = operator
      SoHex8& operator=(const SoHex8& old) = delete;

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


      virtual std::vector<double> ElementCenterRefeCoords();


      /*!
      \brief Return unique ParObject id

      every class implementing ParObject needs a unique id defined at the
      top of this file.
      */
      int UniqueParObjectId() const override { return SoHex8Type::Instance().UniqueParObjectId(); }

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
      \brief Does this element use EAS?
      */
      bool HaveEAS() const override { return (eastype_ != soh8_easnone); };

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


      DRT::ElementType& ElementType() const override { return SoHex8Type::Instance(); }

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

      \param names (out): On return, the derived class has filled names with
                          key names of data it wants to visualize and with int dimensions
                          of that data.
      */
      void VisNames(std::map<std::string, int>& names) override;

      /*!
      \brief Query data to be visualized using BINIO of a given name

      The method is supposed to call this base method to visualize the owner of
      the element.
      If the derived method recognizes a supported data name, it shall fill it
      with corresponding data.
      If it does NOT recognizes the name, it shall do nothing.

      \warning The method must not change size of data

      \param name (in):   Name of data that is currently processed for visualization
      \param data (out):  data to be filled by element if element recognizes the name
      */
      bool VisData(const std::string& name, std::vector<double>& data) override;

      //@}

      void MaterialPostSetup(Teuchos::ParameterList& params) override;

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

      Evaluate so_hex8 element stiffness, mass, internal forces, etc.

      \param params (in/out): ParameterList for communication between control routine
                              and elements
      \param discretization : pointer to discretization for de-assembly
      \param lm (in)        : location matrix for de-assembly
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
      int EvaluateNeumann(Teuchos::ParameterList& params, DRT::Discretization& discretization,
          DRT::Condition& condition, std::vector<int>& lm, CORE::LINALG::SerialDenseVector& elevec1,
          CORE::LINALG::SerialDenseMatrix* elemat1 = nullptr) override;


      // const vector<double> GetFibervec(){return fiberdirection_;};
      /// Evaluate center coordinates in reference system
      void soh8_ElementCenterRefeCoords(
          CORE::LINALG::Matrix<NUMDIM_SOH8, 1>&
              centercoord,  ///< center coordinates in reference system
          CORE::LINALG::Matrix<NUMNOD_SOH8, NUMDIM_SOH8> const& xrefe)
          const;  ///< material coord. of element

      /// Evaluate Gauss-Point coordinates in reference system
      void soh8_GaussPointRefeCoords(CORE::LINALG::Matrix<NUMDIM_SOH8, 1>&
                                         gpcoord,  ///< Gauss-Point coordinates in reference system
          CORE::LINALG::Matrix<NUMNOD_SOH8, NUMDIM_SOH8> const&
              xrefe,      ///< material coord. of element
          int gp) const;  ///< current Gauss-Point

      /*!
      \brief Return value how expensive it is to evaluate this element

      \param double (out): cost to evaluate this element
      */
      double EvaluationCost() override
      {
        if (Material()->MaterialType() == INPAR::MAT::m_struct_multiscale)
          return 25000.0;
        else
          return 10.0;
      }

      void GetCauchyNDirAndDerivativesAtXi(const CORE::LINALG::Matrix<3, 1>& xi,
          const std::vector<double>& disp, const CORE::LINALG::Matrix<3, 1>& n,
          const CORE::LINALG::Matrix<3, 1>& dir, double& cauchy_n_dir,
          CORE::LINALG::SerialDenseMatrix* d_cauchyndir_dd,
          CORE::LINALG::SerialDenseMatrix* d2_cauchyndir_dd2,
          CORE::LINALG::SerialDenseMatrix* d2_cauchyndir_dd_dn,
          CORE::LINALG::SerialDenseMatrix* d2_cauchyndir_dd_ddir,
          CORE::LINALG::SerialDenseMatrix* d2_cauchyndir_dd_dxi,
          CORE::LINALG::Matrix<3, 1>* d_cauchyndir_dn,
          CORE::LINALG::Matrix<3, 1>* d_cauchyndir_ddir,
          CORE::LINALG::Matrix<3, 1>* d_cauchyndir_dxi, const std::vector<double>* temp,
          CORE::LINALG::SerialDenseMatrix* d_cauchyndir_dT,
          CORE::LINALG::SerialDenseMatrix* d2_cauchyndir_dd_dT, const double* concentration,
          double* d_cauchyndir_dc) override;
      //@}

     protected:
      //! struct_ale
      bool structale_;

      //! action parameters recognized by so_hex8
      //  FixMe Deprecated: Should be replaced by the DRT::ELEMENTS::ActionType! hiermeier 04/16
      enum ActionType
      {
        none,
        calc_struct_linstiff,
        calc_struct_nlnstiff,
        calc_struct_internalforce,
        calc_struct_linstiffmass,
        calc_struct_nlnstiffmass,
        calc_struct_nlnstifflmass,  //!< internal force, its stiffness and lumped mass matrix
        calc_struct_nlnstiff_gemm,  //!< internal force, stiffness and mass for GEMM
        calc_struct_stress,
        calc_struct_eleload,
        calc_struct_fsiload,
        calc_struct_update_istep,
        calc_struct_reset_istep,    //!< reset elementwise internal variables
                                    //!< during iteration to last converged state
        calc_struct_store_istep,    //!< store internal information in history
        calc_struct_recover_istep,  //!< recover internal information from history
        calc_struct_reset_all,      //!< reset elementwise internal variables
                                    //!< to state in the beginning of the computation
        calc_struct_energy,         //!< compute internal energy
        prestress_update,
        multi_readrestart,           //!< multi-scale: read restart on microscale
        multi_eas_init,              //!< multi-scale: initialize EAS parameters on microscale
        multi_eas_set,               //!< multi-scale: set EAS parameters on microscale
        multi_calc_dens,             //!< multi-scale: calculate homogenized density
        calc_stc_matrix,             //! calculate scaled director matrix for thin shell structures
        calc_stc_matrix_inverse,     //! calculate inverse of scaled director matrix for thin shell
                                     //! structures
        calc_struct_stifftemp,       //!< TSI specific: mechanical-thermal stiffness
        calc_global_gpstresses_map,  //! basically calc_struct_stress but with assembly of global
                                     //! gpstresses map
        interpolate_velocity_to_point,  //! interpolate the structural velocity to a given point
        calc_struct_mass_volume,  //! calculate mass and volume for reference, material and spatial
                                  //! conf.
        calc_recover              //! recover condensed eas variables
      };

      /*!
       * \brief EAS technology enhancement types of so_hex8
       *
       * Solid Hex8 has EAS enhancement of GL-strains to avoid locking.
       */
      enum EASType  // with meaningfull value for matrix size info
      {
        soh8_easnone = 0,   //!< no EAS i.e. displacement based with tremendous locking
        soh8_eassosh8 = 7,  //!< related to Solid-Shell, 7 parameters to alleviate
                            //!< inplane (membrane) locking and main modes for Poisson-locking
        soh8_easmild = 9,   //!< 9 parameters consisting of modes to alleviate
                            //!< shear locking (bending) and main incompressibility modes
                            //!< (for Solid Hex8)
        soh8_easfull = 21,  //!< 21 parameters to prevent almost all locking modes.
                            //!< Equivalent to all 30 parameters to fully complete element
                            //!< with quadratic modes (see Andelfinger 1993 for details) and
                            //!< therefore also suitable for distorted elements. (for Solid Hex8)
        soh8_easa = 1000    //!< A-type EAS for sosh8p8
                            //!< NOTE: not the actual number of EAS parameters
      };

      //! type of EAS technology
      EASType eastype_;
      //! number of EAS parameters (alphas), defined by 'EASType'
      int neas_;

      struct EASData
      {
        CORE::LINALG::SerialDenseMatrix alpha = {};
        CORE::LINALG::SerialDenseMatrix alpha_backup = {};
        CORE::LINALG::SerialDenseMatrix alphao = {};
        CORE::LINALG::SerialDenseMatrix feas = {};
        CORE::LINALG::SerialDenseMatrix invKaa = {};
        CORE::LINALG::SerialDenseMatrix invKaao = {};
        CORE::LINALG::SerialDenseMatrix Kda = {};
        CORE::LINALG::SerialDenseMatrix Kdao = {};
        CORE::LINALG::SerialDenseMatrix eas_inc = {};
        CORE::LINALG::SerialDenseMatrix eas_inc_backup = {};
        CORE::LINALG::SerialDenseMatrix Kap = {};
      };
      //! EAS data
      EASData easdata_;

      //! vector of inverses of the jacobian in material frame
      std::vector<CORE::LINALG::Matrix<NUMDIM_SOH8, NUMDIM_SOH8>> invJ_;
      //! vector of inverses of the jacobian in material frame
      std::vector<CORE::LINALG::Matrix<NUMDIM_SOH8, NUMDIM_SOH8>> invJmat_;
      //! determinant of Jacobian in material frame
      std::vector<double> detJ_;

      //! evaluate the analytical constitutive matrix
      bool analyticalmaterialtangent_;

      /// prestressing switch & time
      INPAR::STR::PreStress pstype_;
      double pstime_;
      double time_;
      // line search parameter (old step length)
      double old_step_length_;
      /// Prestressing object
      Teuchos::RCP<DRT::ELEMENTS::PreStress> prestress_;
      // compute Jacobian mapping wrt to deformed configuration
      virtual void UpdateJacobianMapping(
          const std::vector<double>& disp, DRT::ELEMENTS::PreStress& prestress);

      //! Update history variables at the end of time step (fiber direction, inelastic deformation)
      //! (braeu 07/16)
      void Update_element(std::vector<double>& disp,  // current displacements
          Teuchos::ParameterList& params,             // algorithmic parameters e.g. time
          const Teuchos::RCP<MAT::Material>& mat);    // material

      // compute defgrd in all gp for given disp
      virtual void DefGradient(const std::vector<double>& disp,
          CORE::LINALG::SerialDenseMatrix& gpdefgrd, DRT::ELEMENTS::PreStress& prestress);


      // internal calculation methods

      /** \brief access one coordinate of all defined GPs
       *
       *  A pointer to the array of the GP-coordinate is returned. If dim==0
       *  the r-coordinate of all GPs will be returned. For dim==1 the s-coordinate
       *  and for dim==2 the t-coordinate.
       *
       *  \param(in) dim: Desired coordinate.
       *  \author hiermeier \date 12/17 */
      const double* soh8_get_coordinate_of_gausspoints(unsigned dim) const;

      /** \brief Evaluate the determinant of the current jacobian at each element
       *  corner and return the smallest value */
      double soh8_get_min_det_jac_at_corners(
          const CORE::LINALG::Matrix<NUMNOD_SOH8, NUMDIM_SOH8>& xcurr) const;

      void soh8_error_handling(const double& det_curr, Teuchos::ParameterList& params, int line_id,
          STR::ELEMENTS::EvalErrorFlag flag);

      //! init the inverse of the jacobian and its determinant in the material configuration
      virtual void InitJacobianMapping();

      //! init the inverse of the jacobian and its determinant in the material configuration
      virtual int InitJacobianMapping(std::vector<double>& dispmat);

      //! Calculate nonlinear stiffness and mass matrix
      virtual void nlnstiffmass(std::vector<int>& lm,  ///< location matrix
          std::vector<double>& disp,                   ///< current displacements
          std::vector<double>* vel,                    // current velocities
          std::vector<double>* acc,                    // current accelerations
          std::vector<double>& residual,               ///< current residual displ
          std::vector<double>& dispmat,                ///< current material displacements
          CORE::LINALG::Matrix<NUMDOF_SOH8, NUMDOF_SOH8>*
              stiffmatrix,                                             ///< element stiffness matrix
          CORE::LINALG::Matrix<NUMDOF_SOH8, NUMDOF_SOH8>* massmatrix,  ///< element mass matrix
          CORE::LINALG::Matrix<NUMDOF_SOH8, 1>* force,       ///< element internal force vector
          CORE::LINALG::Matrix<NUMDOF_SOH8, 1>* forceinert,  // element inertial force vector
          CORE::LINALG::Matrix<NUMDOF_SOH8, 1>*
              force_str,  // element structural force vector (no condensation; for NewtonLS)
          CORE::LINALG::Matrix<NUMGPT_SOH8, MAT::NUM_STRESS_3D>* elestress,  ///< stresses at GP
          CORE::LINALG::Matrix<NUMGPT_SOH8, MAT::NUM_STRESS_3D>* elestrain,  ///< strains at GP
          CORE::LINALG::Matrix<NUMGPT_SOH8, MAT::NUM_STRESS_3D>*
              eleplstrain,                     ///< plastic strains at GP
          Teuchos::ParameterList& params,      ///< algorithmic parameters e.g. time
          INPAR::STR::StressType iostress,     ///< stress output option
          INPAR::STR::StrainType iostrain,     ///< strain output option
          INPAR::STR::StrainType ioplstrain);  ///< plastic strain output option

      //! Calculate nonlinear stiffness and mass matrix (GEMM)
      virtual void soh8_nlnstiffmass_gemm(std::vector<int>& lm,  ///< location matrix
          std::vector<double>& dispo,                            ///< old displacements
          std::vector<double>& disp,                             ///< current displacements
          std::vector<double>& residual,                         ///< current residual displ
          CORE::LINALG::Matrix<NUMDOF_SOH8, NUMDOF_SOH8>*
              stiffmatrix,                                             ///< element stiffness matrix
          CORE::LINALG::Matrix<NUMDOF_SOH8, NUMDOF_SOH8>* massmatrix,  ///< element mass matrix
          CORE::LINALG::Matrix<NUMDOF_SOH8, 1>* force,  ///< element internal force vector
          CORE::LINALG::Matrix<NUMGPT_SOH8, MAT::NUM_STRESS_3D>* elestress,  ///< stresses at GP
          CORE::LINALG::Matrix<NUMGPT_SOH8, MAT::NUM_STRESS_3D>* elestrain,  ///< strains at GP
          CORE::LINALG::Matrix<NUMGPT_SOH8, MAT::NUM_STRESS_3D>*
              eleplstrain,                     ///< plastic strains at GP
          Teuchos::ParameterList& params,      ///< algorithmic parameters e.g. time
          INPAR::STR::StressType iostress,     ///< stress output option
          INPAR::STR::StrainType iostrain,     ///< strain output option
          INPAR::STR::StrainType ioplstrain);  ///< plastic strain output option

      //! Lump mass matrix (bborn 07/08)
      void soh8_lumpmass(CORE::LINALG::Matrix<NUMDOF_SOH8, NUMDOF_SOH8>* emass);

      //! Evaluate Hex8 Shapefcts to keep them static
      std::vector<CORE::LINALG::Matrix<NUMNOD_SOH8, 1>> soh8_shapefcts() const;
      //! Evaluate Hex8 Derivs to keep them static
      std::vector<CORE::LINALG::Matrix<NUMDIM_SOH8, NUMNOD_SOH8>> soh8_derivs() const;

      /*!
       * \brief Evaluate the first derivatives of the shape functions at the Gauss point gp
       *
       * \param derivs first derivatives of the shape functions
       * \param gp Gauss point
       */
      void soh8_derivs(CORE::LINALG::Matrix<NUMDIM_SOH8, NUMNOD_SOH8>& derivs, int gp) const;
      //! Evaluate Hex8 Weights to keep them static
      std::vector<double> soh8_weights() const;

      //! push forward of material stresses to the current, spatial configuration
      void PK2toCauchy(CORE::LINALG::Matrix<MAT::NUM_STRESS_3D, 1>* stress,
          CORE::LINALG::Matrix<NUMDIM_SOH8, NUMDIM_SOH8>* defgrd,
          CORE::LINALG::Matrix<NUMDIM_SOH8, NUMDIM_SOH8>* cauchystress);

      //! push forward of Green-Lagrange strain to Euler-Almansi strains
      void GLtoEA(CORE::LINALG::Matrix<MAT::NUM_STRESS_3D, 1>* glstrain,
          CORE::LINALG::Matrix<NUMDIM_SOH8, NUMDIM_SOH8>* defgrd,
          CORE::LINALG::Matrix<NUMDIM_SOH8, NUMDIM_SOH8>* euler_almansi);

      //@}

      //! EAS technology, init
      void soh8_easinit();

      //! Re-initialize EAS data, needed for sosh8 morphing
      void soh8_reiniteas(DRT::ELEMENTS::SoHex8::EASType EASType);

      //! EAS technology, setup necessary data
      void soh8_eassetup(
          std::vector<CORE::LINALG::SerialDenseMatrix>** M_GP,  // M-matrix evaluated at GPs
          double& detJ0,                                        // det of Jacobian at origin
          CORE::LINALG::Matrix<MAT::NUM_STRESS_3D, MAT::NUM_STRESS_3D>&
              T0invT,  // maps M(origin) local to global
          const CORE::LINALG::Matrix<NUMNOD_SOH8, NUMDIM_SOH8>& xrefe)
          const;  // material element coords

      //! EAS technology, update
      void soh8_easupdate();

      //! EAS technology, restore
      void soh8_easrestore();

      /*!
       * \brief Computes the deformation gradient at a Gauss point
       *
       * \param defgrd [out] : Deformation gradient
       * \param xdisp  [in] : Nodal displacements of the element
       * \param xcurr  [in] : Current nodal coordinates of the element
       * \param gp  [in] : Id of the Gauss point
       */
      void ComputeDeformationGradient(CORE::LINALG::Matrix<NUMDIM_SOH8, NUMDIM_SOH8>& defgrd,
          const CORE::LINALG::Matrix<NUMNOD_SOH8, NUMDIM_SOH8>& xdisp,
          const CORE::LINALG::Matrix<NUMNOD_SOH8, NUMDIM_SOH8>& xcurr, int gp);

      /*! \brief Calculate the deformation gradient that is consistent
       *         with modified (e.g. EAS) GL strain tensor.
       *         Expensive (two polar decomposition), but required, if
       *         the material evaluation is based on the deformation
       *         gradient rather than the GL strain tensor (e.g. plasticity).
       *
       * \param defgrd_disp  (in)  : displacement-based deformation gradient
       * \param glstrain_mod (in)  : modified GL strain tensor (strain-like Voigt notation)
       * \param defgrd_mod   (out) : consistent modified deformation gradient
       */
      void CalcConsistentDefgrd(const CORE::LINALG::Matrix<3, 3>& defgrd_disp,
          CORE::LINALG::Matrix<6, 1> glstrain_mod, CORE::LINALG::Matrix<3, 3>& defgrd_mod) const;


      //! @name Multi-scale related stuff

      /*!
       * \brief Determine a homogenized material density for multi-scale
       * analyses by averaging over the initial volume
       * */
      void soh8_homog(Teuchos::ParameterList& params);

      /*!
       * \brief Set EAS internal variables on the microscale
       *
       * Microscale internal EAS data has to be saved separately for every
       * macroscopic Gauss point and set before the determination of
       * microscale stiffness etc.
       * */
      void soh8_set_eas_multi(Teuchos::ParameterList& params);

      /*!
       * \brief Initialize EAS internal variables on the microscale
       * */
      void soh8_eas_init_multi(Teuchos::ParameterList& params);

      /*!
       * \brief Read restart on the microscale
       * */
      void soh8_read_restart_multi();

      //@}

      //! @name TSI related stuff
      //@{
      /*!
       * \brief Determine a homogenized material density for multi-scale
       * analyses by averaging over the initial volume
       * */
      void GetTemperatureForStructuralMaterial(
          const CORE::LINALG::Matrix<NUMNOD_SOH8, 1>&
              shapefcts,  ///< shape functions of current Gauss point
          Teuchos::ParameterList&
              params  ///< special material parameters e.g. scalar-valued temperature
      );

      //@}

      /// temporary method for compatibility with solidshell, needs clarification
      std::vector<double> getthicknessvector() const
      {
        FOUR_C_THROW("not implemented");
        return std::vector<double>(3);
      };

      /// FDCheck constitutive tensor and/or use the approximation as elastic stiffness matrix
      void EvaluateFiniteDifferenceMaterialTangent(
          CORE::LINALG::Matrix<NUMDOF_SOH8, NUMDOF_SOH8>* stiffmatrix,
          const CORE::LINALG::Matrix<MAT::NUM_STRESS_3D, 1>& stress,
          std::vector<double>& disp,  ///< current displacements
          double detJ_w,              ///< jacobian determinant times gauss weigth
          double detJ, double detJ0, double charelelength,
          const CORE::LINALG::Matrix<MAT::NUM_STRESS_3D, NUMDOF_SOH8>& bop,
          const CORE::LINALG::Matrix<6, NUMDOF_SOH8>& cb,
          const CORE::LINALG::Matrix<NUMDIM_SOH8, NUMNOD_SOH8>& N_XYZ,
          const CORE::LINALG::Matrix<MAT::NUM_STRESS_3D, MAT::NUM_STRESS_3D>& T0invT,
          const std::vector<CORE::LINALG::SerialDenseMatrix>* M_GP,
          const CORE::LINALG::SerialDenseMatrix* alpha, CORE::LINALG::SerialDenseMatrix& M, int gp,
          Teuchos::ParameterList& params);

     private:
      void PackEasData(CORE::COMM::PackBuffer& data) const
      {
        AddtoPack(data, easdata_.alpha);
        AddtoPack(data, easdata_.alpha_backup);
        AddtoPack(data, easdata_.alphao);
        AddtoPack(data, easdata_.feas);
        AddtoPack(data, easdata_.invKaa);
        AddtoPack(data, easdata_.invKaao);
        AddtoPack(data, easdata_.Kda);
        AddtoPack(data, easdata_.Kdao);
        AddtoPack(data, easdata_.eas_inc);
        AddtoPack(data, easdata_.eas_inc_backup);
        AddtoPack(data, easdata_.Kap);
      };

      void UnpackEasData(std::vector<char>::size_type& position, const std::vector<char>& data)
      {
        ExtractfromPack(position, data, easdata_.alpha);
        ExtractfromPack(position, data, easdata_.alpha_backup);
        ExtractfromPack(position, data, easdata_.alphao);
        ExtractfromPack(position, data, easdata_.feas);
        ExtractfromPack(position, data, easdata_.invKaa);
        ExtractfromPack(position, data, easdata_.invKaao);
        ExtractfromPack(position, data, easdata_.Kda);
        ExtractfromPack(position, data, easdata_.Kdao);
        ExtractfromPack(position, data, easdata_.eas_inc);
        ExtractfromPack(position, data, easdata_.eas_inc_backup);
        ExtractfromPack(position, data, easdata_.Kap);
      };

      /** recover elementwise stored stuff
       *
       * \author hiermeier
       * \date 04/16 */
      void soh8_recover(const std::vector<int>& lm, const std::vector<double>& residual);

      void soh8_computeEASInc(
          const std::vector<double>& residual, CORE::LINALG::SerialDenseMatrix* const eas_inc);

      /** \brief create a backup of the eas state consisting of enhanced strains and
       *  the enhanced strain increment
       *
       *  Since this function is optional, the backup data is only allocated if
       *  it is really needed.
       *
       *  \author hiermeier \date 12/17 */
      void soh8_create_eas_backup_state(const std::vector<double>& displ_incr);

      /** \brief recover a backup of the eas state consisting of enhanced strains
       *  and the enhanced strain increment
       *
       *  Since this function is optional, the backup data can only be accessed if
       *  soh8_create_eas_backup_state has been called previously.
       *
       *  \author hiermeier \date 12/17 */
      void soh8_recover_from_eas_backup_state();

     private:
      static const CORE::FE::IntPointsAndWeights<NUMDIM_SOH8> gp_rule_;
    };  // class So_hex8



    //=======================================================================
    //=======================================================================
    //=======================================================================
    //=======================================================================



  }  // namespace ELEMENTS
}  // namespace DRT

FOUR_C_NAMESPACE_CLOSE

#endif
