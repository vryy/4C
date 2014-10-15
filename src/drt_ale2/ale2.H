/*----------------------------------------------------------------------------*/
/*!
\file ale2.H

<pre>
Maintainer: Matthias Mayr
            mayr@mhpc.mw.tum.de
            089 - 289 15362
</pre>
*/
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
#ifndef ALE2_H
#define ALE2_H

/*----------------------------------------------------------------------------*/
/* header inclusions */
#include <Teuchos_RCP.hpp>
#include <Epetra_Vector.h>

#include "../drt_lib/drt_elementtype.H"
#include "../drt_fem_general/drt_utils_integration.H"

/*----------------------------------------------------------------------------*/
/* forward declarations */
namespace DRT
{

class Discretization;

namespace ELEMENTS
{

class Ale2Line;

/*----------------------------------------------------------------------------*/
/* definition of classes */
class Ale2Type : public DRT::ElementType
{
public:

  std::string Name() const { return "Ale2Type"; }

  static Ale2Type & Instance() { return instance_; };

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

  static Ale2Type instance_;
};

/*----------------------------------------------------------------------------*/
/*!
\brief A C++ wrapper for the ale2 element
*/
class Ale2 : public DRT::Element
{
public:

  //! @name Friends
  friend class Ale2Surface;

  //@}
  //! @name Constructors and destructors and related methods

  /*!
  \brief Standard Constructor

  \param id : A unique global id
  */
  Ale2(int id, int owner);

  /*!
  \brief Copy Constructor

  Makes a deep copy of a Element

  */
  Ale2(const Ale2& old);

  /*!
  \brief Deep copy this instance of Ale2 and return pointer to the copy

  The Clone() method is used from the virtual base class Element in cases
  where the type of the derived class is unknown and a copy-ctor is needed

  */
  DRT::Element* Clone() const;

  /*!
  \brief Get shape type of element
  */
  virtual DiscretizationType Shape() const;

  /*!
  \brief Return number of lines of this element
  */
  virtual int NumLine() const
  {if (NumNode()==9 || NumNode()==8 || NumNode()==4) return 4;
    else if (NumNode()==3 || NumNode()==6) return 3;
    else {
      dserror("Could not determine number of lines");
      return -1;
    }
  }

  /*!
  \brief Return number of surfaces of this element
  */
  virtual int NumSurface() const {return 1;}

  /*!
  \brief Return number of volumes of this element (always 1)
  */
  virtual int NumVolume() const {return -1;}

  /*!
  \brief Get vector of Teuchos::RCPs to the lines of this element

  */
  virtual std::vector<Teuchos::RCP<DRT::Element> > Lines();

  /*!
  \brief Get vector of Teuchos::RCPs to the surfaces of this element

  */
  virtual std::vector<Teuchos::RCP<DRT::Element> > Surfaces();

  /*!
  \brief Return unique ParObject id

  every class implementing ParObject needs a unique id defined at the
  top of this file.
  */
  virtual int UniqueParObjectId() const { return Ale2Type::Instance().UniqueParObjectId(); }

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
  virtual ~Ale2();

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
  virtual int NumDofPerNode(const DRT::Node& node) const {return 2;}

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
  { return Ale2Type::Instance(); }

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

  Evaluate ale2 element stiffness, mass, internal forces etc

  \param params (in/out): ParameterList for communication between control routine
                          and elements
  \param elemat1 (out)  : matrix to be filled by element. If NULL on input,
                          the controling method does not epxect the element to fill
                          this matrix.
  \param elemat2 (out)  : matrix to be filled by element. If NULL on input,
                          the controling method does not epxect the element to fill
                          this matrix.
  \param elevec1 (out)  : vector to be filled by element. If NULL on input,
                          the controlling method does not epxect the element
                          to fill this vector
  \param elevec2 (out)  : vector to be filled by element. If NULL on input,
                          the controlling method does not epxect the element
                          to fill this vector
  \param elevec3 (out)  : vector to be filled by element. If NULL on input,
                          the controlling method does not epxect the element
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

  this method evaluates a surfaces Neumann condition on the shell element

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

  //! @name Other

  /// set number of gauss points to element shape default
  DRT::UTILS::GaussRule2D getOptimalGaussrule(const DiscretizationType& distype);

  //@}


private:

  void static_ke(
      DRT::Discretization& dis,
      std::vector<int>& lm,
      Epetra_SerialDenseMatrix* sys_mat, Epetra_SerialDenseVector& residual,
      bool incremental, std::vector<double>& my_dispnp,
      Teuchos::ParameterList& params);

  void static_ke_laplace(
      DRT::Discretization& dis,
      std::vector<int>& lm,
      Epetra_SerialDenseMatrix* sys_mat, Epetra_SerialDenseVector& residual,
      std::vector<double>& displacements,
      Teuchos::ParameterList& params);

  void static_ke_spring(Epetra_SerialDenseMatrix* sys_mat,
      Epetra_SerialDenseVector& residual,
      bool incremental,
      std::vector<double>& displacements);

  void static_ke_nonlinear( const std::vector<int> &lm,
      Epetra_SerialDenseMatrix *sys_mat, Epetra_SerialDenseVector& residual,
      const std::vector<double>& my_dispnp,
      Teuchos::ParameterList &params);

  void edge_geometry(int i, int j, const Epetra_SerialDenseMatrix& xyze,
                     double* length, double* sin_alpha, double* cos_alpha);
  double ale2_area_tria(const Epetra_SerialDenseMatrix& xyze, int i, int j, int k);
  void ale2_tors_spring_tri3(Epetra_SerialDenseMatrix* sys_mat, const Epetra_SerialDenseMatrix& xyze);
  void ale2_tors_spring_quad4(Epetra_SerialDenseMatrix* sys_mat, const Epetra_SerialDenseMatrix& xyze);
  void ale2_torsional(int i, int j, int k, const Epetra_SerialDenseMatrix& xyze,
                      Epetra_SerialDenseMatrix* k_torsion);


  void ale2_jacobianmatrix(const Epetra_SerialDenseMatrix& xrefe,
      const Epetra_SerialDenseMatrix& deriv, Epetra_SerialDenseMatrix& xjm,
      double* det, const int numnode);

  void ale2_jacobianmatrix_cure(const Epetra_SerialDenseMatrix& xcure,
      const Epetra_SerialDenseMatrix& deriv, Epetra_SerialDenseMatrix& xjm_cure,
      double* det_cure, const int numnode);

  void ale2_bop(Epetra_SerialDenseMatrix& bop, Epetra_SerialDenseMatrix& deriv,
      Epetra_SerialDenseMatrix& xjm, Epetra_SerialDenseMatrix& xji,
      const int numnode);

  void ale2_greenlag(Epetra_SerialDenseMatrix& defgrad,
      Epetra_SerialDenseMatrix& greenlag, Epetra_SerialDenseMatrix& deriv,
      Epetra_SerialDenseMatrix& xji, const std::vector<double>& my_dispnp);


  //! action parameters recognized by ale2
  enum ActionType
  {
    none,
    calc_ale_solid, ///< compute stiffness based on fully nonlinear elastic solid with hyperelastic material law
    calc_ale_springs, ///< compute stiffness based on springs algorithm on last time step configuration
    calc_ale_laplace,///< compute stiffness based on laplacian smoothing based on initial configuration
    setup_material, ///< Setup material in case of ElastHyper Tool Box

    // old types will be removed in the near future
    calc_ale_lin, ///< compute stiffness based on linear elastic solid at current configuration
    calc_ale_lin_fixed_ref, ///< compute stiffness based on linear elastic solid in initial configuration
    calc_ale_springs_fixed_ref ///< compute stiffness based on springs algorithm on initial configuration
  };

  //! don't want = operator
  Ale2& operator = (const Ale2& old);


};


//=======================================================================
//=======================================================================
//=======================================================================
//=======================================================================


class Ale2LineType : public DRT::ElementType
{
public:

  std::string Name() const { return "Ale2LineType"; }

  static Ale2LineType & Instance() { return instance_; };

  virtual Teuchos::RCP<DRT::Element> Create( const int id, const int owner );

  virtual void NodalBlockInformation( DRT::Element * dwele, int & numdf, int & dimns, int & nv, int & np ) {}

  virtual void ComputeNullSpace( DRT::Discretization & dis, std::vector<double> & ns, const double * x0, int numdf, int dimns ) {}

private:

  static Ale2LineType instance_;
};


/*!
\brief An element representing a line of a ale2 element
*/
class Ale2Line : public DRT::Element
{
public:

  //! @name Constructors and destructors and related methods

  /*!
  \brief Standard Constructor

  \param id : A unique global id
  \param owner: Processor owning this line
  \param nnode: Number of nodes attached to this element
  \param nodeids: global ids of nodes attached to this element
  \param nodes: the discretizations map of nodes to build ptrs to nodes from
  \param parent: The parent ale element of this line
  \param lline: the local line number of this line w.r.t. the parent element
  */
  Ale2Line(int id, int owner, int nnode, const int* nodeids,
           DRT::Node** nodes, DRT::ELEMENTS::Ale2* parent, const int lline);

  /*!
  \brief Copy Constructor

  Makes a deep copy of a Element

  */
  Ale2Line(const Ale2Line& old);

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
  virtual int UniqueParObjectId() const { return Ale2LineType::Instance().UniqueParObjectId(); }

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
  virtual ~Ale2Line();

  //@}

  //! @name Access methods


  /*!
  \brief Get number of degrees of freedom of a certain node
         (implements pure virtual DRT::Element)

  The element decides how many degrees of freedom its nodes must have.
  As this may vary along a simulation, the element can redecide the
  number of degrees of freedom per node along the way for each of it's nodes
  separately.
  */
  virtual int NumDofPerNode(const DRT::Node& node) const
  {return 2;}

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
  { return Ale2LineType::Instance(); }

  //@}

  //! @name Evaluate methods

  /*!
  \brief Evaluate a Neumann boundary condition

  this method evaluates a line Neumann condition on the ale2 element

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

private:

  // don't want = operator
  Ale2Line& operator = (const Ale2Line& old);

}; // class Ale2Line




} // namespace ELEMENTS
} // namespace DRT

#endif
