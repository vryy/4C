/*!----------------------------------------------------------------------
\file so_disp.H
\brief pure displacement based element for higher order shape functions (quadratic)

<pre>
Maintainer: Axel Gerstenberger
            gerstenberger@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15236
</pre>

*----------------------------------------------------------------------*/

#ifndef SO_DISP_H
#define SO_DISP_H

#include "../drt_fem_general/drt_utils_local_connectivity_matrices.H"
#include "../drt_lib/drt_elementtype.H"
#include "so_base.H"
#include "../drt_fem_general/drt_utils_integration.H"
#include "../drt_mat/material.H"
#include "../drt_inpar/inpar_structure.H"

// Several parameters which are fixed three dimensional solid elements
const int NUMDIM_DISP = 3;      ///< number of dimensions
const int NODDOF_DISP = 3;      ///< number of dofs per node

namespace DRT
{
// forward declarations
class Discretization;

namespace ELEMENTS
{
// forward declarations
class SoDispVolume;

class SoDispType : public DRT::ElementType
{
public:

  std::string Name() const { return "SoDispType"; }

  static SoDispType & Instance();

  virtual DRT::ParObject* Create( const std::vector<char> & data );

  virtual Teuchos::RCP<DRT::Element> Create( const std::string eletype,
                                             const std::string eledistype,
                                             const int id,
                                             const int owner );

  virtual Teuchos::RCP<DRT::Element> Create( const int id, const int owner );

  virtual void NodalBlockInformation( DRT::Element * dwele, int & numdf, int & dimns, int & nv, int & np );

  virtual void ComputeNullSpace( DRT::Discretization & dis, std::vector<double> & ns, const double * x0, int numdf, int dimns );

  virtual void SetupElementDefinition( std::map<std::string,std::map<std::string,DRT::INPUT::LineDefinition> > & definitions );

private:

  static SoDispType instance_;
};

/*!
\brief A general quadratic displacement based solid element

\author frenzel (frenzel@lnm.mw.tum.de)
*/
class SoDisp : public So_base
{
public:
  //! @name Friends
  //@}
  //! @name Constructors and destructors and related methods

  /*!
  \brief Standard Constructor

  \param id : A unique global id
  \param owner : elements owning processor
  */
  SoDisp(int id, int owner);

  /*!
  \brief Copy Constructor

  Makes a deep copy of a Element

  */
  SoDisp(const SoDisp& old);

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
  \brief Return number of lines of this element
  */
  virtual int NumLine() const {return DRT::UTILS::getNumberOfElementLines(Shape());}

  /*!
  \brief Return number of surfaces of this element
  */
  virtual int NumSurface() const {return DRT::UTILS::getNumberOfElementSurfaces(Shape());}

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

  /*!
  \brief Return unique ParObject id

  every class implementing ParObject needs a unique id defined at the
  top of this file.
  */
  virtual int UniqueParObjectId() const { return SoDispType::Instance().UniqueParObjectId(); }

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
  virtual ~SoDisp();

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
  { return SoDispType::Instance(); }

  //@}

  //! @name Input and Creation

  /*!
  \brief Read input for this element
  */
  bool ReadElement(const std::string& eletype,
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
  int Evaluate(Teuchos::ParameterList& params,
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


  //@}

protected:

  //! action parameters recognized by so_hex8
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
    calc_struct_reset_istep, //!< reset elementwise internal variables
                             //!< during iteration to last converged state
    calc_struct_errornorms,  //!< compute error norms (L2,H1,energy)
    calc_init_vol
  };

  //! stress evaluation types of so_disp
  enum StressType
  {
    sodisp_stress_none,
    sodisp_stress_gpxyz,
    sodisp_stress_gprst,
    sodisp_stress_gp123,
    sodisp_stress_ndxyz,
    sodisp_stress_ndrst,
    sodisp_stress_nd123
  };

  StressType                 stresstype_; ///< type of stress evaluation
  DRT::UTILS::GaussRule3D    gaussrule_;  ///< Gaussrule

  int numnod_disp_;      ///< number of nodes
  int numdof_disp_;      ///< total dofs per element
  int numgpt_disp_;      ///< total gauss points per element


  // internal calculation methods

  // don't want = operator
  SoDisp& operator = (const SoDisp& old);

  //! Calculate nonlinear stiffness and mass matrix
  virtual void sodisp_nlnstiffmass(std::vector<int>&       lm,
                                 std::vector<double>&      disp,
                                 std::vector<double>&      residual,
                                 Epetra_SerialDenseMatrix* stiffmatrix,
                                 Epetra_SerialDenseMatrix* massmatrix,
                                 Epetra_SerialDenseVector* force,
                                 Teuchos::ParameterList&   params);         // algorithmic parameters e.g. time

  //! evaluate shape functions and derivatives - non-static due to flexible shape
  void sodisp_shapederiv(
      std::vector<Epetra_SerialDenseVector>& shapefcts,
      std::vector<Epetra_SerialDenseMatrix>& derivs,
      std::vector<double>& weights);


}; // class SoDisp


} // namespace ELEMENTS
} // namespace DRT


#endif /*SO_DISP_H*/

