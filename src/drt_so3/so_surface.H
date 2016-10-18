/*!----------------------------------------------------------------------
\file so_surface.H

\brief class for evaluation of equations on the structural surface

\maintainer Michael Gee

\level 1
*----------------------------------------------------------------------*/
#ifndef SO_SURFACE_H
#define SO_SURFACE_H

#include "../drt_lib/drt_element.H"
#include "../drt_lib/drt_elementtype.H"
#include "../drt_fem_general/drt_utils_integration.H"
#include "../drt_lib/drt_node.H"


namespace DRT
{
namespace ELEMENTS
{

class StructuralSurfaceType : public DRT::ElementType
{
public:

  std::string Name() const { return "StructuralSurfaceType"; }

  static StructuralSurfaceType & Instance();

  virtual Teuchos::RCP<DRT::Element> Create( const int id, const int owner );

  virtual void NodalBlockInformation( DRT::Element * dwele, int & numdf, int & dimns, int & nv, int & np ) {}

  virtual void ComputeNullSpace( DRT::Discretization & dis, std::vector<double> & ns, const double * x0, int numdf, int dimns ) {}

private:

  static StructuralSurfaceType instance_;
};

/*!
\brief A general structural surface element (tri3/tri6/quad4/quad8/quad9/nurbs9)

This element implements all common shapes. It also implements all
common types of loads currently demanded in baci

\author gee (gee@lnm.mw.tum.de)
*/
class StructuralSurface : public DRT::FaceElement
{
public:

  /*!
  \brief Standard Constructor

  \param id : A unique global id
  \param owner: Processor owning this surface
  \param nnode: Number of nodes attached to this element
  \param nodeids: global ids of nodes attached to this element
  \param nodes: the discretizations map of nodes to build ptrs to nodes from
  \param parent: The parent fluid element of this surface
  \param lsurface: the local surface number of this surface w.r.t. the parent element
  */
  explicit StructuralSurface(int id, int owner, int nnode, const int* nodeids,
                             DRT::Node** nodes, DRT::Element* parent, const int lsurface);

  /*!
  \brief Copy Constructor

  Makes a deep copy of a Element

  */
  explicit StructuralSurface(const StructuralSurface& old);

  /*!
  \brief Deep copy this instance of an element and return pointer to the copy

  The Clone() method is used from the virtual base class Element in cases
  where the type of the derived class is unknown and a copy-ctor is needed

  */
  DRT::Element* Clone() const;

  /*!
  \brief Get shape type of element
  */
  virtual DiscretizationType Shape() const;

  /*!
  \brief Return unique ParObject id

  every class implementing ParObject needs a unique id defined at the
  top of the drt_parobject.H file.
  */
  virtual inline int UniqueParObjectId() const
  { return StructuralSurfaceType::Instance().UniqueParObjectId(); }

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
  virtual ~StructuralSurface() {}

  /*!
  \brief Get number of degrees of freedom of a certain node
         (implements pure virtual DRT::Element)

  The element decides how many degrees of freedom its nodes must have.
  As this may vary along a simulation, the element can redecide the
  number of degrees of freedom per node along the way for each of it's nodes
  separately.
  */
  inline virtual int NumDofPerNode(const DRT::Node& node) const
  {return ParentMasterElement()->NumDofPerNode(node);}

  /*!
  \brief Get number of degrees of freedom per element
         (implements pure virtual DRT::Element)

  The element decides how many element degrees of freedom it has.
  It can redecide along the way of a simulation.

  \note Element degrees of freedom mentioned here are dofs that are visible
        at the level of the total system of equations. Purely internal
        element dofs that are condensed internally should NOT be considered.
  */
  inline virtual int NumDofPerElement() const { return 0; }

  /*!
  \brief Print this element
  */
  void Print(std::ostream& os) const;

  virtual DRT::ElementType & ElementType() const
  { return StructuralSurfaceType::Instance(); }

  virtual std::vector<Teuchos::RCP<DRT::Element> > Lines();

  virtual int NumLine() const;

  /*!
  \brief Evaluate a Neumann boundary condition

  this method evaluates a surface Neumann condition on the fluid3 element

  \param params (in/out)    : ParameterList for communication between control routine
                              and elements
  \param discretization (in): A reference to the underlying discretization
  \param condition (in)     : The condition to be evaluated
  \param lm (in)            : location vector of this element
  \param elevec1 (out)      : vector to be filled by element. If NULL on input,

  \return 0 if successful, negative otherwise
  */
  virtual int EvaluateNeumann(Teuchos::ParameterList&   params,
                              DRT::Discretization&      discretization,
                              DRT::Condition&           condition,
                              std::vector<int>&         lm,
                              Epetra_SerialDenseVector& elevec1,
                              Epetra_SerialDenseMatrix* elemat1 = NULL);

  //! Evaluate method for StructuralSurface-Elements
  virtual int Evaluate(Teuchos::ParameterList& params,
                       DRT::Discretization&      discretization,
                       std::vector<int>&         lm,
                       Epetra_SerialDenseMatrix& elematrix1,
                       Epetra_SerialDenseMatrix& elematrix2,
                       Epetra_SerialDenseVector& elevector1,
                       Epetra_SerialDenseVector& elevector2,
                       Epetra_SerialDenseVector& elevector3);

  //! Evaluate method for StructuralSurface-Elements
  virtual int Evaluate(Teuchos::ParameterList& params,
                       DRT::Discretization&      discretization,
                       LocationArray&            la,
                       Epetra_SerialDenseMatrix& elematrix1,
                       Epetra_SerialDenseMatrix& elematrix2,
                       Epetra_SerialDenseVector& elevector1,
                       Epetra_SerialDenseVector& elevector2,
                       Epetra_SerialDenseVector& elevector3);

  /*!
  \brief Evaluate J-integral

  this method evaluate and return components of J-integral used to define strength of singularity at crack tip

  \param discret (in)      : underlying structural discretization
  \param dispvec (in)      : column vector holding the displacement
  \param supp_func (in)    : holds the value of support function at each node
  \param normalcrack (in)  : normal to crack in crack-coordinate-system
  \param tangentcrack (in) : tangent to crack in crack-coordinate-system
   */
  std::vector<double> Evaluate_J_integral_crack( Teuchos::RCP<DRT::Discretization>& discret,
                                                 Teuchos::RCP<Epetra_Vector>& dispvec,
                                                 std::map<int,double> supp_func,
                                                 LINALG::Matrix<3,1>& normalcrack,
                                                 LINALG::Matrix<3,1>& tangentcrack );

  //! Return a pointer to the parent element of this boundary element
  virtual DRT::Element* ParentElement() const
  {return ParentMasterElement();}

  //! Return local surface number
  int LSurfNumber() const {return FaceMasterNumber();}

public: // changed to public rauch 05/2014

  //! action for surface evaluation
  enum ActionType
  {
    none,
    calc_struct_constrvol,
    calc_struct_volconstrstiff,
    calc_struct_constrarea,
    calc_init_vol,
    calc_surfstress_stiff,
    calc_struct_areaconstrstiff,
    calc_struct_monitarea,
    prestress_update,
    prestress_writerestart,
    prestress_readrestart,
    calc_brownian_motion,
    calc_brownian_motion_damping,
    calc_struct_centerdisp,
    calc_struct_rotation,
    calc_undo_struct_rotation,
    calc_struct_area,
    calc_ref_nodal_normals,
    calc_cur_nodal_normals,
    calc_struct_area_poro,
    calc_fluid_traction,
    calc_ecm_traction,
    calc_cur_normal_at_point,
    mark_immersed_elements,
    calc_struct_robinforcestiff,
    calc_cell_growth
  };

protected:

  // don't want = operator
  StructuralSurface& operator = (const StructuralSurface& old);

  // variables

  //! gaussian integration to be used
  DRT::UTILS::GaussRule2D     gaussrule_;

  /*!
  \brief Evaluate a Neumann boundary condition

  this method evaluates normal and detA at gaussian point

  \param detA   (out) : are at gaussian point
  \param normal (out) : The normal at gaussian point, length is detA!
  \param x      (in)  : nodal coords in either material or spatial frame
  \param deriv  (in)  : derivatives of shape functions
  */
  void SurfaceIntegration(double& detA,
                          std::vector<double>& normal,
                          const Epetra_SerialDenseMatrix& x,
                          const Epetra_SerialDenseMatrix& deriv);

  /*!
  \brief Evaluate a Neumann boundary condition

  this method evaluates the normal at gaussian point

  \param normal (out) : The normal at gaussian point, length is detA!
  \param x      (in)  : nodal coords in either material or spatial frame
  \param deriv  (in)  : derivatives of shape functions
  */
  void SurfaceIntegration(std::vector<double>& normal,
                          const Epetra_SerialDenseMatrix& x,
                          const Epetra_SerialDenseMatrix& deriv);

  /*!
  \brief Linearize a Neumann boundary condition analytically

  This method computes the analytical surface derivative if necessary (for follower loads).
  */
  void analytical_DSurfaceIntegration(Epetra_SerialDenseMatrix& d_normal,
                                      const Epetra_SerialDenseMatrix& x,
                                      const Epetra_SerialDenseMatrix& deriv);

  /*!
  \brief Linearize a Neumann boundary condition with automatic differentiation

  This method computes the automatic (Sacado) surface derivative if necessary (for follower loads).
  This method is for development purposes only, and currently not used.
  */
  void automatic_DSurfaceIntegration(Epetra_SerialDenseMatrix& d_normal,
                                     const Epetra_SerialDenseMatrix& x,
                                     const Epetra_SerialDenseMatrix& deriv);

  /*!
  \brief Create matrix with material configuration

  \param x  (out)  : nodal coords in material frame
  */
  inline void MaterialConfiguration(LINALG::SerialDenseMatrix& x) const
  {
    const int numnode = NumNode();
    for (int i=0; i<numnode; ++i)
    {
      x(i,0) = Nodes()[i]->X()[0];
      x(i,1) = Nodes()[i]->X()[1];
      x(i,2) = Nodes()[i]->X()[2];
    }
    return;
  }

  /*!
  \brief Create matrix with spatial configuration

  \param x     (out)  : nodal coords in spatial frame
  \param disp  (int)  : displacements
  */
  inline void SpatialConfiguration(LINALG::SerialDenseMatrix& x,
                                   const std::vector<double> disp) const
  {
    const int numnode = NumNode();
    for (int i=0; i<numnode; ++i)
    {
      x(i,0) = Nodes()[i]->X()[0] + disp[i*3+0];
      x(i,1) = Nodes()[i]->X()[1] + disp[i*3+1];
      x(i,2) = Nodes()[i]->X()[2] + disp[i*3+2];
    }
    return;
  }

  /*!
  \brief Create matrix with spatial configuration

  \param x     (out)  : nodal coords in spatial frame
  \param xrefe (out)  : nodal coords in material frame
  \param disp  (int)  : displacements
  */
  inline void SpatialConfiguration(LINALG::SerialDenseMatrix& x,
                                   LINALG::SerialDenseMatrix& xrefe,
                                   const std::vector<double> disp) const
  {
    const int numnode = NumNode();
    for (int i=0; i<numnode; ++i)
    {
      x(i,0) = xrefe(i,0) + disp[i*3+0];
      x(i,1) = xrefe(i,1) + disp[i*3+1];
      x(i,2) = xrefe(i,2) + disp[i*3+2];
    }
    return;
  }

  //! Submethod to compute the enclosed volume for volume constraint boundary condition
  double ComputeConstrVols
  (
      const LINALG::SerialDenseMatrix& xscurr, ///< current configuration
      const int numnode ///< num nodes
  );


  /*!
  \brief Submethod to compute interfacial area and its first and second
         derivatives w.r.t. the displacements (this is needed for surface energy problems)
  */
  void ComputeAreaDeriv
  (
    const LINALG::SerialDenseMatrix& x, ///< spatial configuration
    const int numnode, ///< number of nodes
    const int ndof,  ///< number of degrees of freedom
    double& A,  ///< area
    Teuchos::RCP<Epetra_SerialDenseVector> Adiff, ///< first derivative
    Teuchos::RCP<Epetra_SerialDenseMatrix> Adiff2 ///< second derivative
  );

  //! Submethod to compute constraint volume and its first and second derivatives w.r.t. the displacements
  void ComputeVolDeriv
  (
      const LINALG::SerialDenseMatrix& x, ///< spatial configuration
      const int numnode, ///< number of nodes
      const int ndof,  ///< number of degrees of freedom
      double& V,  ///< volume
      Teuchos::RCP<Epetra_SerialDenseVector> Vdiff, ///< first derivative
      Teuchos::RCP<Epetra_SerialDenseMatrix> Vdiff2, ///< second derivative
      const int minind = 0,  ///<minimal index to compute enclosed volume with
      const int maxind = 2   ///<maximal index to compute enclosed volume with
  );

  /// Submethod to compute normal vectors at nodes (for reference configuration) which can be assembled
  /// (i.e., added) to obtain averaged nodal normal vectors
  void BuildNormalsAtNodes(Epetra_SerialDenseVector& nodenormals,
                           std::vector<double> mydisp,
                           bool refconfig);

  //! Submethod to compute surface porosity
  void CalculateSurfacePorosity(Teuchos::ParameterList& params,
                                DRT::Discretization&      discretization,
                                LocationArray&            la);
}; // class StructuralSurface : public DRT::Element


} // namespace ELEMENTS
} // namespace DRT


#endif // #ifndef SO_SURFACE_H
