/*----------------------------------------------------------------------*/
/*!
\file so3_ssn_plast.H
\maintainer Alexander Seitz

<pre>
   Maintainer: Alexander Seitz
               seitz@lnm.mw.tum.de
               http://www.lnm.mw.tum.de
               089 - 289-15271
</pre>
*/


/*----------------------------------------------------------------------*
 | definitions                                              seitz 07/13 |
 *----------------------------------------------------------------------*/
#ifndef SO3_SSN_PLAST_H
#define SO3_SSN_PLAST_H

/*----------------------------------------------------------------------*
 | headers                                                  seitz 07/13 |
 *----------------------------------------------------------------------*/
#include "../drt_plastic_ssn/plastic_ssn_manager.H"
#include "so3_ssn_plast_eletypes.H"
#include "../drt_thermo/thermo_ele_impl_utils.H"

// include header because number of Gauss points are determined
// dependent on distype
#include "../drt_fem_general/drt_utils_gausspoints.H"

/*----------------------------------------------------------------------*
 |                                                          seitz 07/13 |
 *----------------------------------------------------------------------*/
namespace DRT
{
class Discretization;

namespace UTILS
{
struct PlastSsnData;
}
namespace ELEMENTS
{

  //! A C++ version of a 3 dimensional solid element with modifications von Mises
  //! plasticity using a semi-smooth Newton method
  template<DRT::Element::DiscretizationType distype>
  class So3_Plast : public virtual So_base
  {
    //! @name Friends
    friend class So_hex8PlastType;
    friend class So_hex27PlastType;

  public:

    //@}
    //! @name Constructors and destructors and related methods


    //! Standard Constructor
    So3_Plast(
      int id,  //!< (i) this element's global id
      int owner  //!< elements owner
      );

    //! Copy Constructor
    //! Makes a deep copy of a Element
    So3_Plast(const So3_Plast& old);

    //! Deep copy this instance of Solid3 and return pointer to the copy
    //!
    //! The Clone() method is used from the virtual base class Element in cases
    //! where the type of the derived class is unknown and a copy-ctor is needed
    DRT::Element* Clone() const;

    //! Destructor
    virtual ~So3_Plast(){};

    //@}

    //! number of element nodes
    static const int nen_ = DRT::UTILS::DisTypeToNumNodePerEle<distype>::numNodePerElement;
    //! number of space dimensions
    static const int nsd_ = 3;
    //! number of dofs per node
    static const int numdofpernode_ = 3;
    //! total dofs per element
    static const int numdofperelement_ = numdofpernode_ * nen_;
    //! number of strains/stresses
    static const int numstr_ = 6;
    //! number of Gauss points per element (value is added in so3_thermo.cpp)
    int numgpt_;
    //! static const is required for fixedsizematrices
    //! TODO maybe more beauty is possible
    static const int numgpt_post = THR::DisTypeToSTRNumGaussPoints<distype>::nquad;


    //! @name Acess methods

    //! Return unique ParObject id
    //!
    //! every class implementing ParObject needs a unique id defined at the top of
    //! this file.
    virtual int UniqueParObjectId() const;


    //! Pack this class so it can be communicated
    //! Pack and \ref Unpack are used to communicate this element
    virtual void Pack(DRT::PackBuffer& data) const;

    //! Unpack data from a char vector into this class
    //! Pack and \ref Unpack are used to communicate this element
    virtual void Unpack(const std::vector<char>& data);

    //! Get number of degrees of freedom of a certain node in case of multiple
    //! dofsets (implements pure virtual DRT::Element)
    //!
    //! The element decides how many degrees of freedom its nodes must have.
    virtual int NumDofPerNode(const DRT::Node& node) const
    {
        return nsd_;
    };

    //! Get number of degrees of freedom of this element
    virtual int NumDofPerElement() const
    {
        return 0;
    };

    //@}

    //! @name Access methods

    //! Print this element
    void Print(std::ostream& os) const  ;

    //! return elementtype
     virtual DRT::ElementType& ElementType() const;

     //! return element shape
     virtual DiscretizationType Shape() const {return distype;};

     /*!
     \brief Return number of volumes of this element
     */
     virtual int NumVolume() const ;

     /*!
     \brief Return number of surfaces of this element
     */
     virtual int NumSurface() const ;

     /*!
     \brief Return number of lines of this element
     */
     virtual int NumLine() const ;

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
    //@}

    //! @name Input and Creation

    //! Query names of element data to be visualized using BINIO
    //!
    //! The element fills the provided map with key names of
    //! visualization data the element wants to visualize AT THE CENTER
    //! of the element geometry. The values is supposed to be dimension of the
    //! data to be visualized. It can either be 1 (scalar), 3 (vector), 6 (sym. tensor)
    //! or 9 (nonsym. tensor)
    //!
    //! Example:
    //! \code
    //!   // Name of data is 'Owner', dimension is 1 (scalar value)
    //!   names.insert(std::pair<string,int>("Owner",1));
    //!   // Name of data is 'StressesXYZ', dimension is 6 (sym. tensor value)
    //!   names.insert(std::pair<string,int>("StressesXYZ",6));
    //! \endcode
    //!
    //!  names (out): On return, the derived class has filled names with key
    //!               names of data it wants to visualize and with int dimensions
    //!               of that data.
    virtual void VisNames(std::map<std::string,int>& names);

    //!  Query data to be visualized using BINIO of a given name
    //!
    //! The method is supposed to call this base method to visualize the owner of
    //! the element.
    //! If the derived method recognizes a supported data name, it shall fill it
    //! with corresponding data.
    //! If it does NOT recognizes the name, it shall do nothing.
    //!
    //! warning: the method must not change size of data
    //!
    //!  name (in):   Name of data that is currently processed for visualization
    //! \param data (out):  data to be filled by element if element recognizes the name
    virtual bool VisData(const std::string& name, std::vector<double>& data);

    //! read input for this element
    bool ReadElement(
      const std::string& eletype,  //!< so3plast(fbar)
      const std::string& eledistype,  //!< hex8,tet4,...
      DRT::INPUT::LineDefinition* linedef  //!< what parameters have to be read
      );

    //@}

    //! @name Evaluation

    //! evaluate an element
    //! evaluate element stiffness, mass, internal forces, etc.
    //!
    //! if NULL on input, the controlling method does not expect the element
    //!  to fill these matrices or vectors.
    //!
    //!  \return 0 if successful, negative otherwise
    int Evaluate(
        Teuchos::ParameterList& params,//!< ParameterList for communication between control routine and elements
        DRT::Discretization& discretization, //!< pointer to discretization for de-assembly
        DRT::Element::LocationArray& la,  //!< location array for de-assembly
        Epetra_SerialDenseMatrix& elemat1_epetra, //!< (stiffness-)matrix to be filled by element.
        Epetra_SerialDenseMatrix& elemat2_epetra,//!< (mass-)matrix to be filled by element.
        Epetra_SerialDenseVector& elevec1_epetra,//!< (internal force-)vector to be filled by element
        Epetra_SerialDenseVector& elevec2_epetra, //!< vector to be filled by element
        Epetra_SerialDenseVector& elevec3_epetra//!< vector to be filled by element
      );

    //@}

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

    //! init the inverse of the jacobian and its determinant in the material
    //! configuration
    virtual void InitJacobianMapping();

    // get parameter list from ssn_plast_manager
    virtual void ReadParameterList(Teuchos::RCP<Teuchos::ParameterList> plparams);

    //@}

  protected:

    //! don't want = operator
    So3_Plast& operator = (const So3_Plast& old);

    //! action parameters recognized by so_hex8
    enum ActionType
    {
      none,
      calc_struct_linstiff,      //!< geometrical linear analysis: stiffness matrix
      calc_struct_nlnstiff,      //!< stiffness matrix
      calc_struct_internalforce, //!< internal force
      calc_struct_linstiffmass,  //!< geometrical linear analysis: internal force,
                                 //!< its stiffness and mass matrix
      calc_struct_nlnstiffmass,  //!< internal force, its stiffness and mass matrix
      calc_struct_nlnstifflmass, //!< internal force, its stiffness and lumped mass matrix
      calc_struct_stress,        //!< calculate stresses
      calc_struct_update_istep,  //!< update all at element level
      calc_struct_reset_istep,   //!< reset elementwise internal variables
                                 //!< during iteration to last converged state
                                 //!< needed for predictor TangDis
      postprocess_stress,        //!< postprocess stresses
      calc_struct_stifftemp,      //!< calculate the coupling block K_dT in monolithic TSI
      calc_struct_energy          //!< compute internal elastic energy
    };

    /*!
     * \brief EAS technology enhancement types of so_hex8
     *
     * Solid Hex8 has EAS enhancement of GL-strains to avoid locking.
     * New enum for so3_plast. Currently only for hex8
     */
    enum EASType  // with meaningful value for matrix size info
    {
      soh8p_easnone = 0,  //!< no EAS i.e. displacement based with tremendous locking
      soh8p_eassosh8 = 7, //!< related to Solid-Shell, 7 parameters to alleviate
                          //!< in-plane (membrane) locking and main modes for Poisson-locking
      soh8p_easmild = 9,  //!< 9 parameters consisting of modes to alleviate
                          //!< shear locking (bending) and main incompressibility modes
                          //!< (for Solid Hex8)
                          //!< The sosh18 also uses 9 eas parameters. Therefore, we can
                          //!< re-use this EASType here.
      soh8p_easfull = 21, //!< 21 parameters to prevent almost all locking modes.
                          //!< Equivalent to all 30 parameters to fully complete element
                          //!< with quadratic modes (see Andelfinger 1993 for details) and therefore
                          //!< also suitable for distorted elements.
                          //!< (for Solid Hex8)
    };

    //! number of plastic variables at each gauss point
    enum PlSpinType
    {
      zerospin=5, //!< 5 parameters for zero plastic spin (symmetric traceless tensor)
      plspin=8,   //!< 8 parameters for zero plastic spin (non-symmetric traceless tensor)
    };

    //! vector of inverses of the jacobian in material frame
    std::vector<LINALG::Matrix<nsd_,nsd_> > invJ_;
    //! determinant of Jacobian in material frame
    std::vector<double> detJ_;
    //! vector of coordinates of current integration point in reference coordinates
    std::vector<LINALG::Matrix<nsd_,1> > xsi_;
    //! Gauss point weights
    std::vector<double>                  wgt_;

    //! @name plasticity related stuff


    //! Calculate nonlinear stiffness and mass matrix with condensed plastic matrices
    virtual void nln_stiffmass(
        std::vector<int>&              lm,             // location matrix
        std::vector<double>&           disp,           // current displacements
        std::vector<double>&           residual,       // current residual displ
        std::vector<double>&           vel,            // current velocities
        std::vector<double>&           temp,           // current temperatures
        std::vector<double>&           temp_res,       // current residual temperatures
        LINALG::Matrix<numdofperelement_,numdofperelement_>* stiffmatrix, // element stiffness matrix
        LINALG::Matrix<numdofperelement_,numdofperelement_>* massmatrix,  // element mass matrix
        LINALG::Matrix<numdofperelement_,1>* force,                 // element internal force vector
        LINALG::Matrix<numdofperelement_,1>* force_str,          // structure force
        LINALG::Matrix<numgpt_post,numstr_>* elestress,   // stresses at GP
        LINALG::Matrix<numgpt_post,numstr_>* elestrain,   // strains at GP
        Teuchos::ParameterList&        params,         // algorithmic parameters e.g. time
        const INPAR::STR::StressType   iostress,  // stress output option
        const INPAR::STR::StrainType   iostrain,  // strain output option
        const int MyPID  // processor id
        );

    //! Calculate the coupling matrix K_dT for monolithic TSI
    virtual void nln_kdT_tsi(
        LINALG::Matrix<numdofperelement_,nen_>* k_dT,  //!< (o): mechanical thermal stiffness term at current gp
        Teuchos::ParameterList& params
        );

    // Add plastic increment of converged state to plastic history for nonlinear kinematics
    virtual void UpdatePlasticDeformation_nln(PlSpinType spintype);

    //@}

    //! calculate nonlinear B-operator
    virtual void CalculateBop(
      LINALG::Matrix<numstr_,numdofperelement_>* bop,
      LINALG::Matrix<nsd_,nsd_>* defgrd,
      LINALG::Matrix<nsd_,nen_>* N_XYZ
      );

    //! calculate linear B-operator
    virtual void CalculateBoplin(
      LINALG::Matrix<numstr_,numdofperelement_>* boplin,
      LINALG::Matrix<nsd_,nen_>* N_XYZ
      );

    /// Extrapolate Gauss-point values (e.g. stresses) to nodes and store results in elevectors
    /// todo: Unfortunately, there is no universal extrapolation function in all base elements,
    /// i.e. the hex8 elements would call soh8_expol, hex27 would call soh27_expol...
    /// For now, we copy the extrapolation function here and allow the use only for hex8 elements
    void soh8_expol(
        LINALG::Matrix<numgpt_post,numstr_>& stresses,  ///< gp stresses
        Epetra_MultiVector& expolstresses);  ///< nodal stresses

    // EAS element techonolgy *******************************************
    //! Initialize data for EAS (once)
    virtual void EasInit();

    //! setup EAS for each evaluation
    virtual void EasSetup(
        std::vector<Epetra_SerialDenseMatrix>** M_GP,    // M-matrix evaluated at GPs
        double& detJ0,                      // det of Jacobian at origin
        LINALG::Matrix<numstr_,numstr_>& T0invT,   // maps M(origin) local to global
        const LINALG::Matrix<nen_,nsd_>& xrefe)    // material element coords
    ;


    //! @name plasticity related stuff

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
    virtual void CalcConsistentDefgrd(LINALG::Matrix<3,3> defgrd_disp,
                              LINALG::Matrix<6,1> glstrain_mod,
                              LINALG::Matrix<3,3>& defgrd_mod);


    /*! \brief Evaluate the NCP function and the linearization
     *         and condense the additional degrees of freedom
     *         into the stiffness matrix block.
     *         Required data from the element evaluate is handed in,
     *         with potential EAS matrices or F-bar linearizations.
     *         The template argument "spinype" decides, if
     *         the additional evolution equation for the plastic
     *         spin is linearized and solved for.
     */
    template<int spintype>
    void CondensePlasticity(
        const LINALG::Matrix<nsd_,nsd_>& defgrd,
        const LINALG::Matrix<nsd_,nsd_>& deltaLp,
        const LINALG::Matrix<numstr_,numdofperelement_>& bop,
        const LINALG::Matrix<numstr_,1>* RCG,
        const int MyPID,
        const double detJ_w,
        const int gp,
        const double temp,
        Teuchos::ParameterList& params,
        LINALG::Matrix<numdofperelement_,1>* force,
        LINALG::Matrix<numdofperelement_,numdofperelement_>* stiffmatrix,
        int& num_active_gp,
        double& lp_res,
        Epetra_SerialDenseMatrix* M=NULL,
        Epetra_SerialDenseMatrix* Kda=NULL,
        std::vector<Epetra_SerialDenseVector>* dHda=NULL,
        const double* f_bar_factor=NULL,
        const LINALG::Matrix<numdofperelement_,1>* htensor=NULL
        );

    /*! \brief Recover the condensed additional degrees of freedom.
     *         Required data from the element evaluate is handed in,
     *         with potential EAS matrices.
     *         The template argument "spinype" decides, if
     *         the additional dofs for the plastic
     *         spin have to be recovered.
     */
    template<int spintype>
    void RecoverPlasticity(
        const LINALG::Matrix<numdofperelement_,1>& res_d,
        const INPAR::STR::PredEnum& pred,
        const int gp,
        const int MyPID,
        const double temp,
        const Teuchos::ParameterList& params,
        LINALG::Matrix<3,3>& deltaLp,
        double& lp_inc,
        bool recover
        );


    /*! \brief Return if plastic spin is solved for.
     */
    virtual bool HavePlasticSpin();


    //! calculate internal elastic energy
    double CalcIntEnergy(
        std::vector<double>&   disp,
        Teuchos::ParameterList&         params);

    //@}


    // algoirthmic parameters
    bool fbar_;

//     Data container of plastic_ssn_manager
    Teuchos::RCP<DRT::UTILS::PlastSsnData> data_;

    static const int VOIGT3X3SYM_[3][3];
    static const int VOIGT3X3NONSYM_[3][3];


    // plasticity ********************************************
    // Kbb^-1 at each Gauss point for recovery of inner variables
    std::vector<Epetra_SerialDenseMatrix> KbbInv_;

    // Kbd at each Gauss point for recovery of inner variables
    std::vector<Epetra_SerialDenseMatrix> Kbd_;

    // f_b at each Gauss point for recovery of inner variables
    std::vector<Epetra_SerialDenseVector> fbeta_;

    // plastic flow at each Gauss point at last Newton iteration
    std::vector<Epetra_SerialDenseVector> dDp_last_iter_;

    // increment of plastic flow over last Newton step (needed for Newton line search)
    std::vector<Epetra_SerialDenseVector> dDp_inc_;

    PlSpinType       plspintype_;

    // EAS element technology ************************************
    Teuchos::RCP<Epetra_SerialDenseMatrix> KaaInv_;
    Teuchos::RCP<Epetra_SerialDenseMatrix> Kad_;
    Teuchos::RCP<Epetra_SerialDenseMatrix> KaT_;
    Teuchos::RCP<LINALG::Matrix<numdofperelement_,nen_> > KdT_eas_;
    Teuchos::RCP<Epetra_SerialDenseVector> feas_;
    Teuchos::RCP<std::vector<Epetra_SerialDenseMatrix> > Kba_;
    Teuchos::RCP<Epetra_SerialDenseVector> alpha_eas_;
    Teuchos::RCP<Epetra_SerialDenseVector> alpha_eas_last_timestep_;
    Teuchos::RCP<Epetra_SerialDenseVector> alpha_eas_delta_over_last_timestep_;
    Teuchos::RCP<Epetra_SerialDenseVector> alpha_eas_inc_;
    EASType                  eastype_;
    int                      neas_;

    // TSI ******************************************************
    bool tsi_;

    /// derivative of the internal force vector w.r.t. temperature at this GP
    /// derivative w.r.t. gp temperature is sufficient as the gp temperature depends linearily
    /// on the nodal values.
    Teuchos::RCP<std::vector<LINALG::Matrix<numdofperelement_,1> > > dFintdT_;

    /// derivative of NCP w.r.t. temperatures at gp
    Teuchos::RCP<std::vector<Epetra_SerialDenseVector> > KbT_;

    /// Temperature at each gp in last Newton iteration
    /// this is needed for the recovery of the plastic flow using KbT_
    /// As the gp temperature depends linearly on the nodal temperature dofs
    /// a scalar value is sufficient
    Teuchos::RCP<std::vector<double> > temp_last_;

  };  // class So3_Plast

}  // namespace ELEMENTS
}  // namespace DRT

/*----------------------------------------------------------------------*/
#endif  // SO3_SSN_PLAST_H
