/*!----------------------------------------------------------------------
\file so_nstet.H

<pre>
Maintainer: Michael Gee
            gee@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15239
</pre>

*----------------------------------------------------------------------*/

#ifdef D_SOLID3
#ifdef CCADISCRET
#ifndef SO_NSTET4_H
#define SO_NSTET4_H


#include "Teuchos_RefCountPtr.hpp"
#include "Epetra_Vector.h"
#include "../linalg/linalg_serialdensevector.H"
#include "../linalg/linalg_serialdensematrix.H"
using namespace std;
using namespace Teuchos;

#include "../drt_lib/drt_element.H"
#include "../drt_lib/drt_elementtype.H"
#include "../drt_lib/drt_parobjectfactory.H"
#include "../drt_lib/drt_node.H"
#include "../drt_inpar/inpar_structure.H"

#include "../drt_mat/material.H"

#include "so_surface.H"
#include "so_line.H"

#include "../linalg/linalg_fixedsizematrix.H"

#define NUMNOD_NSTET       4    ///< number of nodes
#define NODDOF_NSTET       3    ///< number of dofs per node
#define NUMDOF_NSTET      12    ///< total dofs per element
#define NUMSTR_NSTET       6    ///< number of strains/stresses
#define NUMGPT_NSTET       1    ///< total gauss points per element
#define NUMDIM_NSTET       3    ///< number of dimensions
#define NUMCOORD_NSTET     4    ///< number of shape function cooordinates (ksi1-ksi4)

#define ALPHA_NSTET      0.2    ///< stabilization parameter
#define BETA_NSTET       0.2    ///< stabilization parameter

namespace DRT
{
// forward declarations
class Discretization;

namespace ELEMENTS
{
class StructuralVolume;
class NStet;

  enum NStet_StabType
  {
    so_nstet4_stab_none,
    so_nstet4_voldev,
    so_nstet4_puso,
    so_nstet4_vol,
    so_nstet4_dev
  };


//=======================================================================
//=======================================================================

class NStetType : public DRT::ElementType
{
  //! allow nstet element to access the nodal data
  friend class DRT::ELEMENTS::NStet;

public:

  std::string Name() const { return "NStetType"; }

  static NStetType & Instance() { return instance_; };

  virtual DRT::ParObject* Create( const vector<char> & data );

  virtual Teuchos::RCP<DRT::Element> Create( const string eletype,
                                             const string eledistype,
                                             const int id,
                                             const int owner );

  virtual int Initialize(DRT::Discretization& dis);

  virtual void PreEvaluate(DRT::Discretization& dis,
                           Teuchos::ParameterList& p,
                           RCP<LINALG::SparseOperator> systemmatrix1,
                           RCP<LINALG::SparseOperator> systemmatrix2,
                           RCP<Epetra_Vector>          systemvector1,
                           RCP<Epetra_Vector>          systemvector2,
                           RCP<Epetra_Vector>          systemvector3);

  virtual void NodalBlockInformation( DRT::Element * dwele, int & numdf, int & dimns, int & nv, int & np );

  virtual void ComputeNullSpace( DRT::Discretization & dis, std::vector<double> & ns, const double * x0, int numdf, int dimns );

private:

  static NStetType instance_;

  // type of stabilization used
  NStet_StabType                                      stabtype_;

  // map of row nodes adjacent to nstet elements
  map<int,DRT::Node*>                                 noderids_;
  // map of column nstet elements
  map<int,DRT::ELEMENTS::NStet*>                      elecids_;
  // map of nodal stresses
  map<int,vector<double> >                            nodestress_;
  // map of nodal stresses
  map<int,vector<double> >                            nodestrain_;

  // vector of elements adjacent to each row node
  map<int,vector<DRT::ELEMENTS::NStet*> >             adjele_;
  // map of nodes forming a patch around a row node
  map<int,map<int,DRT::Node*> >                       adjnode_;
  // map of location vectors for patch around a node
  map<int,vector<int> >                               adjlm_;

  // vector of nodal pressures in nodal col map
  RCP<Epetra_Vector>                                  pnodecol_;
  // vector of nodal J in nodal col map
  RCP<Epetra_Vector>                                  Jnodecol_;
  // vector of nodal F in nodal col map
  RCP<Epetra_MultiVector>                             Fnodecol_;

  void NodalIntegration(Epetra_SerialDenseMatrix*       stiff,
                        Epetra_SerialDenseVector*       force,
                        map<int,DRT::Node*>&            nodepatch,
                        vector<DRT::ELEMENTS::NStet*>&  adjele,
                        vector<double>*                 nodalstress,
                        vector<double>*                 nodalstrain,
                        const INPAR::STR::StressType    iostress,
                        const INPAR::STR::StrainType    iostrain,
                        double&                         J,
                        double&                         pressure,
                        LINALG::Matrix<9,1>&            Fvec);


  void SelectMaterial(RCP<MAT::Material> mat,
                      LINALG::Matrix<6,1>& stress,
                      LINALG::Matrix<6,6>& cmat,
                      double& density,
                      LINALG::Matrix<6,1>& glstrain,
                      LINALG::Matrix<3,3>& defgrd,
                      int gp);

   // compute deviatoric stresses and tangent
  static void DevStressTangent(LINALG::Matrix<NUMSTR_NSTET,1>& Sdev,
                               LINALG::Matrix<NUMSTR_NSTET,NUMSTR_NSTET>& CCdev,
                               LINALG::Matrix<NUMSTR_NSTET,NUMSTR_NSTET>& CC,
                               const LINALG::Matrix<NUMSTR_NSTET,1>& S,
                               const LINALG::Matrix<NUMDIM_NSTET,NUMDIM_NSTET>& C,
                               double& pressure);

   // compute deviatoric and volumetric stresses and tangent
  static void DevVolStressTangent(
                               LINALG::Matrix<NUMSTR_NSTET,1>& Sdev,
                               LINALG::Matrix<NUMSTR_NSTET,1>& Svol,
                               LINALG::Matrix<NUMSTR_NSTET,NUMSTR_NSTET>& CCdev,
                               LINALG::Matrix<NUMSTR_NSTET,NUMSTR_NSTET>& CCvol,
                               LINALG::Matrix<NUMSTR_NSTET,NUMSTR_NSTET>& CC,
                               const LINALG::Matrix<NUMSTR_NSTET,1>& S,
                               const LINALG::Matrix<NUMDIM_NSTET,NUMDIM_NSTET>& C,
                               double& elepressure);

  /// volumetric and isochoric stabilization methods
  void VolDevStab(LINALG::Matrix<NUMDIM_NSTET,NUMDIM_NSTET>& cauchygreen,
                  LINALG::Matrix<NUMSTR_NSTET,1>& stress,
                  LINALG::Matrix<NUMSTR_NSTET,NUMSTR_NSTET>& cmat,
                  double& pressure);

  /// isochoric stabilization methods
  void DevStab(LINALG::Matrix<NUMDIM_NSTET,NUMDIM_NSTET>& cauchygreen,
               LINALG::Matrix<NUMSTR_NSTET,1>& stress,
               LINALG::Matrix<NUMSTR_NSTET,NUMSTR_NSTET>& cmat);


  /// little helper routine
  inline void MatrixtoVector(const LINALG::Matrix<3,3>& F, LINALG::Matrix<9,1>& Fvec)
  {
    Fvec(0) = F(0,0);
    Fvec(1) = F(0,1);
    Fvec(2) = F(0,2);
    Fvec(3) = F(1,0);
    Fvec(4) = F(1,1);
    Fvec(5) = F(1,2);
    Fvec(6) = F(2,0);
    Fvec(7) = F(2,1);
    Fvec(8) = F(2,2);
    return;
  }
  inline void VectortoMatrix(const LINALG::Matrix<9,1>& Fvec, LINALG::Matrix<3,3>& F)
  {
    F(0,0) = Fvec(0);
    F(0,1) = Fvec(1);
    F(0,2) = Fvec(2);
    F(1,0) = Fvec(3);
    F(1,1) = Fvec(4);
    F(1,2) = Fvec(5);
    F(2,0) = Fvec(6);
    F(2,1) = Fvec(7);
    F(2,2) = Fvec(8);
    return;
  }

};


//=======================================================================
//=======================================================================


/*!
\brief A nodal-averaged strain 4-noded tet element according to Mike Puso

*/

class NStet : public DRT::Element
{
public:

  friend class NStetType;


  /*!
  \brief Standard Constructor

  \param id : A unique global id
  */
  NStet(int id, int owner);

  /*!
  \brief Copy Constructor

  Makes a deep copy of a Element

  */
  NStet(const NStet& old);

  /*!
  \brief Deep copy this instance of Solid3 and return pointer to the copy

  The Clone() method is used from the virtual base class Element in cases
  where the type of the derived class is unknown and a copy-ctor is needed

  */
  inline DRT::Element* Clone() const
  { return new DRT::ELEMENTS::NStet(*this); }

  /*!
  \brief Get shape type of element
  */
  inline virtual DiscretizationType Shape() const { return tet4; }

  /*!
  \brief Return number of volumes of this element
  */
  virtual int NumVolume() const {return 1;}

  /*!
  \brief Return number of surfaces of this element
  */
  virtual int NumSurface() const {return 4;}

  /*!
  \brief Return number of lines of this element
  */
  virtual int NumLine() const {return 6;}

  /*!
  \brief Get vector of RefCountPtrs to the lines of this element

  */
  virtual vector<RCP<DRT::Element> > Lines();

  /*!
  \brief Get vector of RefCountPtrs to the surfaces of this element

  */
  virtual vector<RCP<DRT::Element> > Surfaces();

  /*!
  \brief Get vector of RefCountPtrs to the volumes of this element

  */
  virtual vector<RCP<DRT::Element> > Volumes();

  /*!
  \brief Return unique ParObject id

  every class implementing ParObject needs a unique id defined at the
  top of this file.
  */
  virtual int UniqueParObjectId() const { return NStetType::Instance().UniqueParObjectId(); }

  /*!
  \brief Pack this class so it can be communicated

  \ref Pack and \ref Unpack are used to communicate this element

  */
  virtual void Pack(vector<char>& data) const;

  /*!
  \brief Unpack data from a char vector into this class

  \ref Pack and \ref Unpack are used to communicate this element

  */
  virtual void Unpack(const vector<char>& data);

  /*!
  \brief Destructor

  */
  virtual ~NStet();

  /*!
  \brief Get number of degrees of freedom of a certain node
         (implements pure virtual DRT::Element)

  The element decides how many degrees of freedom its nodes must have.
  As this may vary along a simulation, the element can redecide the
  number of degrees of freedom per node along the way for each of it's nodes
  separately.
  */
  virtual inline int NumDofPerNode(const DRT::Node& node) const {return 3;}

  /*!
  \brief Get number of degrees of freedom per element
         (implements pure virtual DRT::Element)

  The element decides how many element degrees of freedom it has.
  It can redecide along the way of a simulation.

  \note Element degrees of freedom mentioned here are dofs that are visible
        at the level of the total system of equations. Purely internal
        element dofs that are condensed internally should NOT be considered.
  */
  virtual inline int NumDofPerElement() const { return 0; }

  /*!
  \brief Print this element
  */
  void Print(ostream& os) const;

  virtual NStetType & ElementType() const
  { return NStetType::Instance(); }

  /*!
  \brief Read input for this element
  */
  //bool ReadElement();

  /*!
  \brief Read input for this element
  */
  bool ReadElement(const std::string& eletype,
                   const std::string& distype,
                   DRT::INPUT::LineDefinition* linedef);

  /*!
  \brief Evaluate an element

  Evaluate so_tet4 element stiffness, mass, internal forces, etc.

  \param params (in/out): ParameterList for communication between control routine
                          and elements
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
  int Evaluate(ParameterList& params,
               DRT::Discretization&      discretization,
               vector<int>&              lm,
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
  virtual int EvaluateNeumann(ParameterList& params,
                              DRT::Discretization&      discretization,
                              DRT::Condition&           condition,
                              vector<int>&              lm,
                              Epetra_SerialDenseVector& elevec1,
                              Epetra_SerialDenseMatrix* elemat1 = NULL);


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
    postprocess_stress,
    calc_struct_eleload,
    calc_struct_fsiload,
    calc_struct_update_istep,
    calc_struct_update_imrlike,
    calc_struct_reset_istep,  //!< reset elementwise internal variables
                              //!< during iteration to last converged state
    calc_struct_B_and_F
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
  StressType                                               stresstype_;
  //! type of stabilization
  NStet_StabType                                           stabtype_;
  //! number of the material law
  int                                                      material_;
  //! vector of history variables for each Gauss point
  vector<RCP<MAT::Material> >                              mat_;
  //! material volume of element
  double                                                   V_;
  //! current spatial volume of element
  //bool                                                     visNew_;
  //double                                                   v_;
  //! derivs of shape functions wrt to material coords
  LINALG::Matrix<NUMNOD_NSTET,NUMDIM_NSTET>                nxyz_;
  //! element's def-gradient and flag indicating age of it
  bool                                                     FisNew_;
  LINALG::Matrix<NUMDIM_NSTET,NUMDIM_NSTET>                F_;

  // internal calculation methods

  /// volumetric and isochoric stabilization methods
  void VolDevStab(LINALG::Matrix<NUMDIM_NSTET,NUMDIM_NSTET>& defgrd,
                  LINALG::Matrix<NUMDIM_NSTET,NUMDIM_NSTET>& cauchygreen,
                  LINALG::Matrix<NUMSTR_NSTET,1>& stress,
                  LINALG::Matrix<NUMSTR_NSTET,NUMSTR_NSTET>& cmat,
                  double& density);
  void VolDevStabLinear(LINALG::Matrix<NUMDIM_NSTET,NUMDIM_NSTET>& defgrd,
                        LINALG::Matrix<NUMDOF_NSTET,1>* force);

  /// isochoric stabilization methods
  void DevStab(LINALG::Matrix<NUMDIM_NSTET,NUMDIM_NSTET>& defgrd,
               LINALG::Matrix<NUMDIM_NSTET,NUMDIM_NSTET>& cauchygreen,
               LINALG::Matrix<NUMSTR_NSTET,1>& stress,
               LINALG::Matrix<NUMSTR_NSTET,NUMSTR_NSTET>& cmat,
               double& density);

  // don't want = operator
  NStet& operator = (const NStet& old);

  // init the inverse of the jacobian and its determinant
  // in the material configuration
  virtual void InitElement();

  // return volume of element
  virtual inline double Volume() { return V_; }

  //! Shape functions
  inline void ShapeFunction(LINALG::Matrix<NUMNOD_NSTET,1>& funct,
                            const double& e1,
                            const double& e2,
                            const double& e3,
                            const double& e4)
  {
    // shape function is N_i = xsi_i
    funct(0) = e1;
    funct(1) = e2;
    funct(2) = e3;
    funct(3) = e4;
    return;
  }

  //! Shape function derivatives
  inline void ShapeFunctionDerivatives(LINALG::Matrix<NUMNOD_NSTET,NUMCOORD_NSTET>& deriv)
  {
    // Ni,j = 1.0 for i==j, otherwise 0.0
    deriv.Clear();
    deriv(0,0) = 1.0;
    deriv(1,1) = 1.0;
    deriv(2,2) = 1.0;
    deriv(3,3) = 1.0;
    return;
  }

  //! volumetric stabilization
  virtual void VolStabilization(map<int,DRT::ELEMENTS::NStet*>&          adjele,
                                map<int,DRT::Node*>&                     adjnode,
                                vector<int>&                             adjlm,
                                vector<int>&                             adjlmowner,
                                vector<map<int,DRT::ELEMENTS::NStet*> >& nodaladjele,
                                vector<map<int,DRT::Node*> >&            nodaladjnode,
                                LINALG::SerialDenseMatrix&               stiff,
                                LINALG::SerialDenseVector&               force);


  //! almost standard displ. tet4 calc routine
  virtual void nstetnlnstiffmass(vector<int>&              lm,
                                vector<double>&           disp,
                                vector<double>&           residual,
                                LINALG::Matrix<NUMDOF_NSTET,NUMDOF_NSTET>* stiffmatrix,
                                LINALG::Matrix<NUMDOF_NSTET,NUMDOF_NSTET>* massmatrix,
                                LINALG::Matrix<NUMDOF_NSTET,          1>* force);

  //! almost standard displ. tet4 calc routine
  virtual void calcstress(vector<int>&              lm,
                          vector<double>&           disp,
                          LINALG::Matrix<NUMGPT_NSTET,NUMSTR_NSTET>* stresses,
                          LINALG::Matrix<NUMGPT_NSTET,NUMSTR_NSTET>* elestrain,
                          const INPAR::STR::StressType     iostress,
                          const INPAR::STR::StrainType     iostrain);

  /// Extrapolate Gauss-point values (e.g. stresses) to nodes and store them in elevecs
  void so_tet4_expol
  (
      LINALG::Matrix<NUMGPT_NSTET,NUMSTR_NSTET>& stresses, ///< gp stresses
      LINALG::Matrix<NUMDOF_NSTET,1>& elevec1, ///< nodal normal stresses
      LINALG::Matrix<NUMDOF_NSTET,1>& elevec2  ///< nodal shear stresses
  );

  //! lump mass matrix (bborn 07/08)
  void nstetlumpmass(LINALG::Matrix<NUMDOF_NSTET,NUMDOF_NSTET>* emass);

  void DeformationGradient(vector<double>& disp);

  //void DeformedVolume(vector<double>& disp);

  void so_nstet_expol(LINALG::Matrix<NUMGPT_NSTET,NUMSTR_NSTET>& stresses,
                     LINALG::Matrix<NUMNOD_NSTET,NUMSTR_NSTET>& nodalstresses);



  void SelectMaterial(LINALG::Matrix<6,1>& stress,
                      LINALG::Matrix<6,6>& cmat,
                      double& density,
                      LINALG::Matrix<6,1>& glstrain,
                      LINALG::Matrix<3,3>& defgrd,
                      int gp);

}; // class NStet






} // namespace ELEMENTS
} // namespace DRT

#endif  // #ifndef SO_NSTET4_H
#endif  // #ifdef  CCADISCRET
#endif  // #ifdef  D_SOLID3


