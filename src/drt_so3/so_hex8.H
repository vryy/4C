/*----------------------------------------------------------------------*/
/*!
\file so_hex8.H

\maintainer Michael Gee

*/
/*----------------------------------------------------------------------*/
#ifndef SOH8_H
#define SOH8_H

#include "../drt_lib/drt_element.H"
#include "../drt_lib/drt_elementtype.H"
#include "../drt_inpar/inpar_structure.H"
#include "../drt_mat/material.H"

#include "so_base.H"

#include "../drt_lib/drt_discret.H"

// Several parameters which are fixed for Solid Hex8
const int NUMNOD_SOH8 = 8;      ///< number of nodes
const int NODDOF_SOH8 = 3;      ///< number of dofs per node
const int NUMDOF_SOH8 = 24;     ///< total dofs per element
const int NUMGPT_SOH8 = 8;      ///< total gauss points per element
const int NUMDIM_SOH8 = 3;      ///< number of dimensions


namespace DRT
{
// forward declarations
class Discretization;

namespace ELEMENTS
{
// forward declarations
class InvDesign;
class PreStress;
class So_sh8Type;


class So_hex8Type : public DRT::ElementType
{
public:

  std::string Name() const { return "So_hex8Type"; }

  static So_hex8Type & Instance();

  virtual DRT::ParObject* Create( const std::vector<char> & data );

  virtual Teuchos::RCP<DRT::Element> Create( const std::string eletype,
                                             const std::string eledistype,
                                             const int id,
                                             const int owner );

  virtual Teuchos::RCP<DRT::Element> Create( const int id, const int owner );

  virtual int Initialize(DRT::Discretization& dis);

  virtual void NodalBlockInformation( DRT::Element * dwele, int & numdf, int & dimns, int & nv, int & np );

  virtual void ComputeNullSpace( DRT::Discretization & dis, std::vector<double> & ns, const double * x0, int numdf, int dimns );

  virtual void SetupElementDefinition( std::map<std::string,std::map<std::string,DRT::INPUT::LineDefinition> > & definitions );

private:

  static So_hex8Type instance_;
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

\author frenzel (frenzel@lnm.mw.tum.de)
*/
class So_hex8 : public So_base
{
public:

  //! @name Friends
  friend class So_hex8Type;
  friend class So_sh8Type;

  //@}
  //! @name Constructors and destructors and related methods

  /*!
  \brief Standard Constructor

  \param id : A unique global id
  \param owner : elements owner
  */
  So_hex8(int id, int owner);

  /*!
  \brief Copy Constructor

  Makes a deep copy of a Element

  */
  So_hex8(const So_hex8& old);

  /*!
  \brief Deep copy this instance of Solid3 and return pointer to the copy

  The Clone() method is used from the virtual base class Element in cases
  where the type of the derived class is unknown and a copy-ctor is needed

  */
  DRT::Element* Clone() const;

  /*!
  \brief Get shape type of element
  */
  virtual DiscretizationType Shape() const;

  /*!
  \brief Return number of volumes of this element
  */
  virtual int NumVolume() const {return 1;}

  /*!
  \brief Return number of surfaces of this element
  */
  virtual int NumSurface() const {return 6;}

  /*!
  \brief Return number of lines of this element
  */
  virtual int NumLine() const {return 12;}

  /*!
  \brief Get vector of Teuchos::RCPs to the lines of this element

  */
  virtual std::vector<Teuchos::RCP<DRT::Element> > Lines();

  /*!
  \brief Get vector of Teuchos::RCPs to the surfaces of this element

  */
  virtual std::vector<Teuchos::RCP<DRT::Element> > Surfaces();

  /*!
  \brief Get vector of Teuchos::RCPs to the volumes of this element

  */
  virtual std::vector<Teuchos::RCP<DRT::Element> > Volumes();


  virtual std::vector<double> ElementCenterRefeCoords();


  /*!
  \brief Return unique ParObject id

  every class implementing ParObject needs a unique id defined at the
  top of this file.
  */
  virtual int UniqueParObjectId() const { return So_hex8Type::Instance().UniqueParObjectId(); }

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
  virtual ~So_hex8();

  //@}

  //! @name Acess methods
  /*!
  \brief Does this element use EAS?
  */
  bool HaveEAS() const { return (eastype_!=soh8_easnone); };

  /*!
  \brief Get number of degrees of freedom of a certain node
         (implements pure virtual DRT::Element)

  The element decides how many degrees of freedom its nodes must have.
  As this may vary along a simulation, the element can redecide the
  number of degrees of freedom per node along the way for each of it's nodes
  separately.
  */
  virtual int NumDofPerNode(const DRT::Node& node) const {return 3;}


  /*!
  \brief Get number of degrees of freedom per element
         (implements pure virtual DRT::Element)

  The element decides how many element degrees of freedom it has.
  It can redecide along the way of a simulation.

  \note Element degrees of freedom mentioned here are dofs that are visible
        at the level of the total system of equations. Purely internal
        element dofs that are condensed internally should NOT be considered.
  */
  virtual int NumDofPerElement() const { return 0; }

  /*!
  \brief Print this element
  */
  void Print(std::ostream& os) const;



  virtual DRT::ElementType & ElementType() const
  { return So_hex8Type::Instance(); }

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
  virtual void VisNames(std::map<std::string,int>& names);

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
  virtual bool VisData(const std::string& name, std::vector<double>& data);

  //@}

  //! @name Input and Creation

  /*!
  \brief Read input for this element
  */
  virtual bool ReadElement(const std::string& eletype,
                           const std::string& distype,
                           DRT::INPUT::LineDefinition* linedef);

  //@}

  //! @name Evaluation

  /*!
  \brief Evaluate an element

  Evaluate so_hex8 element stiffness, mass, internal forces, etc.

  \param params (in/out): ParameterList for communication between control routine
                          and elements
  \param discretization : pointer to discretization for de-assembly
  \param lm (in)        : location matrix for de-assembly
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
  virtual int Evaluate(Teuchos::ParameterList& params,
                       DRT::Discretization&      discretization,
                       std::vector<int>&         lm,
                       Epetra_SerialDenseMatrix& elemat1,
                       Epetra_SerialDenseMatrix& elemat2,
                       Epetra_SerialDenseVector& elevec1,
                       Epetra_SerialDenseVector& elevec2,
                       Epetra_SerialDenseVector& elevec3);

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
                              DRT::Discretization&      discretization,
                              DRT::Condition&           condition,
                              std::vector<int>&         lm,
                              Epetra_SerialDenseVector& elevec1,
                              Epetra_SerialDenseMatrix* elemat1 = NULL);


  //const vector<double> GetFibervec(){return fiberdirection_;};
  const std::vector<double> soh8_ElementCenterRefeCoords();

  /*!
  \brief Return value how expensive it is to evaluate this element

  \param double (out): cost to evaluate this element
  */
  virtual double EvaluationCost()
  {
    if(Material()->MaterialType() == INPAR::MAT::m_struct_multiscale)
      return 25000.0;
    else
      return 10.0;
  }

  /// evaluate Cauchy stress at given point in parameter space
  virtual void GetCauchyAtXi(
      const LINALG::Matrix<NUMDIM_SOH8,1>& xi,
      const std::vector<double>& disp,
      LINALG::Matrix<NUMDIM_SOH8,NUMDIM_SOH8>& cauchy,
      Epetra_SerialDenseMatrix& dsdd);

  //! struct_ale
  bool structale_;
  //@}

protected:

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
    calc_struct_reset_istep,  //!< reset elementwise internal variables
                              //!< during iteration to last converged state
    calc_struct_store_istep, //!< store internal information in history
    calc_struct_recover_istep,  //!< recover internal information from history
    calc_struct_reset_all,  //!< reset elementwise internal variables
                            //!< to state in the beginning of the computation
    calc_struct_energy,       //!< compute internal energy
    calc_struct_errornorms,   //!< compute error norms (L2,H1,energy)
    postprocess_stress,       //!< postprocess stresses
    prestress_update,
    inversedesign_update,
    inversedesign_switch,
    multi_readrestart,        //!< multi-scale: read restart on microscale
    multi_eas_init,           //!< multi-scale: initialize EAS parameters on microscale
    multi_eas_set,            //!< multi-scale: set EAS parameters on microscale
    multi_calc_dens,          //!< multi-scale: calculate homogenized density
    calc_potential_stiff,     //!< volume interaction potential
    calc_stc_matrix,          //! calculate scaled director matrix for thin shell structures
    calc_stc_matrix_inverse,  //! calculate inverse of scaled director matrix for thin shell structures
    calc_struct_stifftemp,    //!< TSI specific: mechanical-thermal stiffness
    calc_global_gpstresses_map, //! basically calc_struct_stress but with assembly of global gpstresses map
    interpolate_velocity_to_point, //! interpolate the structural velocity to a given point
  };

  /*!
   * \brief EAS technology enhancement types of so_hex8
   *
   * Solid Hex8 has EAS enhancement of GL-strains to avoid locking.
   */
  enum EASType  // with meaningfull value for matrix size info
  {
    soh8_easnone = 0,  //!< no EAS i.e. displacement based with tremendous locking
    soh8_eassosh8 = 7,  //!< related to Solid-Shell, 7 parameters to alleviate
                        //!< inplane (membrane) locking and main modes for Poisson-locking
    soh8_easmild = 9,  //!< 9 parameters consisting of modes to alleviate
                       //!< shear locking (bending) and main incompressibility modes
                       //!< (for Solid Hex8)
    soh8_easfull = 21,  //!< 21 parameters to prevent almost all locking modes.
                        //!< Equivalent to all 30 parameters to fully complete element
                        //!< with quadratic modes (see Andelfinger 1993 for details) and therefore
                        //!< also suitable for distorted elements.
                        //!< (for Solid Hex8)
    soh8_easa = 1000  //!< A-type EAS for sosh8p8
                      //!< NOTE: not the actual number of EAS parameters
  };
  //! type of EAS technology
  EASType                                                               eastype_;
  //! number of EAS parameters (alphas), defined by 'EASType'
  int                                                                   neas_;
  //! container containing all kind of solid3 thisandthat
  DRT::Container                                                        data_;
  //! vector of inverses of the jacobian in material frame
  std::vector<LINALG::Matrix<NUMDIM_SOH8,NUMDIM_SOH8> >                 invJ_;
  //! vector of inverses of the jacobian in material frame
  std::vector<LINALG::Matrix<NUMDIM_SOH8,NUMDIM_SOH8> >                 invJmat_;
  //! determinant of Jacobian in material frame
  std::vector<double>                                                   detJ_;

  //! evaluate the analytical constitutive matrix
  bool                                                                  analyticalmaterialtangent_;

  /// prestressing switch & time
  INPAR::STR::PreStress                                                 pstype_;
  double                                                                pstime_;
  double                                                                time_;
  // line search parameter (old step length)
  double                                                                old_step_length_;
  /// Prestressing object
  Teuchos::RCP<DRT::ELEMENTS::PreStress>                                prestress_;
  // compute Jacobian mapping wrt to deformed configuration
  virtual void UpdateJacobianMapping(const std::vector<double>& disp,
                                     DRT::ELEMENTS::PreStress& prestress);
  // compute defgrd in all gp for given disp
  virtual void DefGradient(const std::vector<double>& disp,
                           Epetra_SerialDenseMatrix& gpdefgrd,
                           DRT::ELEMENTS::PreStress& prestress);

  // inverse design analysis objects
  friend class DRT::ELEMENTS::InvDesign;
  Teuchos::RCP<DRT::ELEMENTS::InvDesign>                                invdesign_;


  // internal calculation methods

  //! don't want = operator
  So_hex8& operator = (const So_hex8& old);


  //! init the inverse of the jacobian and its determinant in the material configuration
  virtual void InitJacobianMapping();

  //! init the inverse of the jacobian and its determinant in the material configuration
  virtual void InitJacobianMapping(std::vector<double>& dispmat);

  //! Calculate nonlinear stiffness and mass matrix
  virtual void nlnstiffmass(
      std::vector<int>& lm,  ///< location matrix
      std::vector<double>& disp,  ///< current displacements
      std::vector<double>*      vel,           // current velocities
      std::vector<double>*      acc,           // current accelerations
      std::vector<double>& residual,  ///< current residual displ
      std::vector<double>& dispmat,  ///< current material displacements
      LINALG::Matrix<NUMDOF_SOH8,NUMDOF_SOH8>* stiffmatrix,  ///< element stiffness matrix
      LINALG::Matrix<NUMDOF_SOH8,NUMDOF_SOH8>* massmatrix,  ///< element mass matrix
      LINALG::Matrix<NUMDOF_SOH8,1>* force,  ///< element internal force vector
      LINALG::Matrix<NUMDOF_SOH8,1>* forceinert,                 // element inertial force vector
      LINALG::Matrix<NUMDOF_SOH8,1>* force_str,                 // element structural force vector (no condensation; for NewtonLS)
      LINALG::Matrix<NUMGPT_SOH8,MAT::NUM_STRESS_3D>* elestress,  ///< stresses at GP
      LINALG::Matrix<NUMGPT_SOH8,MAT::NUM_STRESS_3D>* elestrain,  ///< strains at GP
      LINALG::Matrix<NUMGPT_SOH8,MAT::NUM_STRESS_3D>* eleplstrain,  ///< plastic strains at GP
      Teuchos::ParameterList& params,  ///< algorithmic parameters e.g. time
      const INPAR::STR::StressType iostress,  ///< stress output option
      const INPAR::STR::StrainType iostrain,  ///< strain output option
      const INPAR::STR::StrainType ioplstrain);  ///< plastic strain output option

  //! Calculate nonlinear stiffness and mass matrix (GEMM)
  virtual void soh8_nlnstiffmass_gemm(
      std::vector<int>& lm,  ///< location matrix
      std::vector<double>& dispo,  ///< old displacements
      std::vector<double>& disp,  ///< current displacements
      std::vector<double>& residual,  ///< current residual displ
      LINALG::Matrix<NUMDOF_SOH8,NUMDOF_SOH8>* stiffmatrix,  ///< element stiffness matrix
      LINALG::Matrix<NUMDOF_SOH8,NUMDOF_SOH8>* massmatrix,  ///< element mass matrix
      LINALG::Matrix<NUMDOF_SOH8,1>* force,  ///< element internal force vector
      LINALG::Matrix<NUMGPT_SOH8,MAT::NUM_STRESS_3D>* elestress,  ///< stresses at GP
      LINALG::Matrix<NUMGPT_SOH8,MAT::NUM_STRESS_3D>* elestrain,  ///< strains at GP
      LINALG::Matrix<NUMGPT_SOH8,MAT::NUM_STRESS_3D>* eleplstrain,  ///< plastic strains at GP
      Teuchos::ParameterList& params,  ///< algorithmic parameters e.g. time
      const INPAR::STR::StressType iostress,  ///< stress output option
      const INPAR::STR::StrainType iostrain,  ///< strain output option
      const INPAR::STR::StrainType ioplstrain);  ///< plastic strain output option

  //! Lump mass matrix (bborn 07/08)
  void soh8_lumpmass(LINALG::Matrix<NUMDOF_SOH8,NUMDOF_SOH8>* emass);

  /// Extrapolate Gauss-point values (e.g. stresses) to nodes and store results in elevectors
  void soh8_expol(
      LINALG::Matrix<NUMGPT_SOH8,MAT::NUM_STRESS_3D>& stresses,  ///< gp stresses
      Epetra_MultiVector& expolstresses);  ///< nodal stresses

  //! remodeling for fibers at the end of time step (st 01/10)
  void soh8_remodel(
        std::vector<int>&         lm,             // location matrix
        std::vector<double>&      disp,           // current displacements
        Teuchos::ParameterList&   params,         // algorithmic parameters e.g. time
        Teuchos::RCP<MAT::Material>        mat);            // material

  //! Evaluate Hex8 Shapefcts to keep them static
  const std::vector<LINALG::Matrix<NUMNOD_SOH8,1> > soh8_shapefcts();
  //! Evaluate Hex8 Derivs to keep them static
  const std::vector<LINALG::Matrix<NUMDIM_SOH8,NUMNOD_SOH8> > soh8_derivs();
  //! Evaluate Hex8 Weights to keep them static
  const std::vector<double> soh8_weights();

  //! push forward of material stresses to the current, spatial configuration
  void PK2toCauchy(
    LINALG::Matrix<MAT::NUM_STRESS_3D,1>* stress,
    LINALG::Matrix<NUMDIM_SOH8,NUMDIM_SOH8>* defgrd,
    LINALG::Matrix<NUMDIM_SOH8,NUMDIM_SOH8>* cauchystress
    );

  //! push forward of Green-Lagrange strain to Euler-Almansi strains
  void GLtoEA(
    LINALG::Matrix<MAT::NUM_STRESS_3D,1>* glstrain,
    LINALG::Matrix<NUMDIM_SOH8,NUMDIM_SOH8>* defgrd,
    LINALG::Matrix<NUMDIM_SOH8,NUMDIM_SOH8>* euler_almansi
    );

  //@}

  //! EAS technology, init
  void soh8_easinit();

  //! Re-initialize EAS data, needed for sosh8 morphing
  void soh8_reiniteas(const DRT::ELEMENTS::So_hex8::EASType EASType);

  //! EAS technology, setup necessary data
  void soh8_eassetup(
          std::vector<Epetra_SerialDenseMatrix>** M_GP,    // M-matrix evaluated at GPs
          double& detJ0,                      // det of Jacobian at origin
          LINALG::Matrix<MAT::NUM_STRESS_3D,MAT::NUM_STRESS_3D>& T0invT,   // maps M(origin) local to global
          const LINALG::Matrix<NUMNOD_SOH8,NUMDIM_SOH8>& xrefe);    // material element coords

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
  void CalcConsistentDefgrd(LINALG::Matrix<3,3> defgrd_disp,
                            LINALG::Matrix<6,1> glstrain_mod,
                            LINALG::Matrix<3,3>& defgrd_mod);


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
  void soh8_set_eas_multi(Teuchos::ParameterList&  params);

  /*!
   * \brief Initialize EAS internal variables on the microscale
   * */
  void soh8_eas_init_multi(Teuchos::ParameterList&  params);

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
    const LINALG::Matrix<NUMNOD_SOH8,1>& shapefcts,  ///< shape functions of current Gauss point
    Teuchos::ParameterList& params  ///< special material parameters e.g. scalar-valued temperature
    );

  //@}

  /// temporary method for compatibility with solidshell, needs clarification
  std::vector<double> getthicknessvector() const {dserror("not implemented"); return std::vector<double>(3);};

  /// FDCheck constitutive tensor and/or use the approximation as elastic stiffness matrix
  void EvaluateFiniteDifferenceMaterialTangent(
      LINALG::Matrix<NUMDOF_SOH8,NUMDOF_SOH8>* stiffmatrix,
      const LINALG::Matrix<MAT::NUM_STRESS_3D,1>& stress,
      std::vector<double>& disp, ///< current displacements
      const double detJ_w, ///< jacobian determinant times gauss weigth
      const double detJ,
      const double detJ0,
      const double charelelength,
      const LINALG::Matrix<MAT::NUM_STRESS_3D,NUMDOF_SOH8>& bop,
      const LINALG::Matrix<6,NUMDOF_SOH8>& cb,
      const LINALG::Matrix<NUMDIM_SOH8,NUMNOD_SOH8>& N_XYZ,
      const LINALG::Matrix<MAT::NUM_STRESS_3D,MAT::NUM_STRESS_3D>& T0invT,
      const std::vector<Epetra_SerialDenseMatrix>* M_GP,
      const Epetra_SerialDenseMatrix* alpha,
      LINALG::SerialDenseMatrix& M,
      const int gp,
      Teuchos::ParameterList& params
      );

private:
  /* recover elementwise stored stuff
   *
   * \author hiermeier
   * \date 04/16 */
  void soh8_recover(std::vector<int>&         lm,
                    std::vector<double>&      disp,
                    std::vector<double>&      residual);

}; // class So_hex8



//=======================================================================
//=======================================================================
//=======================================================================
//=======================================================================



} // namespace ELEMENTS
} // namespace DRT

#endif  // #ifndef SOH8_H
