/*----------------------------------------------------------------------*/
/*!
\file red_airway.H

\brief Implements an airway element

\level 3

\maintainer Lena Yoshihara
*/
/*----------------------------------------------------------------------*/
#ifndef RED_AIRWAY_H_
#define RED_AIRWAY_H_

#include "Teuchos_RCP.hpp"
#include "Epetra_Vector.h"


#include "../drt_lib/drt_element.H"
#include "../drt_lib/drt_parobjectfactory.H"
#include "../drt_lib/drt_elementtype.H"
#include "../drt_lib/drt_node.H"
#include "../drt_fem_general/drt_utils_integration.H"


// forward declarations
namespace DRT
{
// forward declarations
class Discretization;

namespace ELEMENTS
{
// forward declarations
class RedAirwayImplInterface;
template <DRT::Element::DiscretizationType distype> class AirwayImpl;

class RedAirwayType : public DRT::ElementType
{
public:

  std::string Name() const { return "RedAirwayType"; }

  static RedAirwayType& Instance();

  virtual DRT::ParObject* Create( const std::vector<char> & data );

  virtual Teuchos::RCP<DRT::Element> Create( const std::string eletype,
                                             const std::string eledistype,
                                             const int id,
                                             const int owner );

  virtual Teuchos::RCP<DRT::Element> Create( const int id, const int owner );

  virtual void NodalBlockInformation( DRT::Element * dwele, int & numdf, int & dimns, int & nv, int & np ) {}

  virtual void ComputeNullSpace( DRT::Discretization & dis, std::vector<double> & ns, const double * x0, int numdf, int dimns ) {}

  virtual void SetupElementDefinition( std::map<std::string,std::map<std::string,DRT::INPUT::LineDefinition> > & definitions );

private:

  static RedAirwayType instance_;
};

/*!
\brief A C++ wrapper for the airway element

\author ismail (ismail@lnm.mw.tum.de)
*/
class RedAirway : public DRT::Element
{
public:

  //! @name Friends
  friend class RedAirwayImplInterface;
  friend class AirwayImpl<DRT::Element::line2>;

  //@}
  //! @name Constructors and destructors and related methods

  /*!
  \brief Standard Constructor

  \param id : A unique global id
  */
  RedAirway(int id, int owner);

  /*!
  \brief Copy Constructor

  Makes a deep copy of a Element

  */
  RedAirway(const RedAirway& old);

  /*!
  \brief Deep copy this instance of RedAirway and return pointer to the copy

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
  {if (NumNode()==2) return 1;
      else {
          dserror("Could not determine number of lines");
          return -1;
      }
  }

  /*!
  \brief Return number of surfaces of this element (always 1)
  */
  virtual int NumSurface() const {return -1;}

  /*!
  \brief Return number of volumes of this element (always 1)
  */
  virtual int NumVolume() const {return -1;}

  /*!
  \brief Get vector of Teuchos::RCPs to the lines of this element
  */
  virtual std::vector<Teuchos::RCP<DRT::Element> > Lines();


  /*!
  \brief Return unique ParObject id

  every class implementing ParObject needs a unique id defined at the
  top of this file.
  */
  virtual int UniqueParObjectId() const { return RedAirwayType::Instance().UniqueParObjectId(); }

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
  virtual ~RedAirway();

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
  virtual int NumDofPerNode(const DRT::Node& node) const {return 1;}

  /*!
  \brief Get number of degrees of freedom per element
         (implements pure virtual DRT::Element)

  The element decides how many element degrees of freedom it has.
  It can redecide along the way of a simulation.

  \note Element degrees of freedom mentioned here are dofs that are visible
        at the level of the total system of equations. Purely internal
        element dofs that are condensed internally should NOT be considered.
  */
  virtual int NumDofPerElement() const
  {
    return 0;
  }

  /*!
  \brief Print this element
  */
  void Print(std::ostream& os) const;

  virtual RedAirwayType & ElementType() const
  { return RedAirwayType::Instance(); }

  //@}

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

  Evaluate airway element stiffness, mass, internal forces etc

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

  An element derived from this class uses the EvaluateNeumann method to receive commands
  and parameters from some control routine in params and evaluates a Neumann boundary condition
  given in condition

  \note This class implements a dummy of this method that prints a warning and
        returns false.

  \param params (in/out)    : ParameterList for communication between control routine
                              and elements
  \param discretization (in): A reference to the underlying discretization
  \param condition (in)     : The condition to be evaluated
  \param lm (in)            : location vector of this element
  \param elevec1 (out)      : Force vector to be filled by element

  \return 0 if successful, negative otherwise
  */
  virtual int EvaluateNeumann(Teuchos::ParameterList& params,
                              DRT::Discretization&      discretization,
                              DRT::Condition&           condition,
                              std::vector<int>&         lm,
                              Epetra_SerialDenseVector& elevec1,
                              Epetra_SerialDenseMatrix* elemat1 = NULL);

  /*!
  \brief Evaluate a Neumann boundary condition

  this method evaluates a line Neumann condition on the airway element

  \return 0 if successful, negative otherwise
  */
  virtual int EvaluateDirichlet(Teuchos::ParameterList&   params,
                                DRT::Discretization&      discretization,
                                DRT::Condition&           condition,
                                std::vector<int>&         lm,
                                Epetra_SerialDenseVector& elevec1);


  //@}

  //! @name Other
  std::string Type(){return elemType_;}

  //! @name Other
  std::string Resistance(){return resistance_;}

  //! @name Other
  std::string ElemSolvingType(){return elemsolvingType_;}

  //@}

  //! get cross-sectional area
  //  double         getA(){ return A_;}

  //! get youngs modulus of the wall
  //  double         getEw(){ return Ew_;}

  //! get youngs modulus of the air
  //  double         getEa(){ return Ea_;}

  //! get wall thickness
  //  double         gettw(){ return tw_;}

  //! set qin
  //  void           setqin (double Qin ){ qin_ = Qin;}

  //! set qout
  //  void           setqout(double Qout){ qout_ = Qout;}

  //! get element parameters
  void           getParams(std::string name, double & var);

  //! get element parameters
  void           getParams(std::string name, int    & var);

private:

  //! action parameters recognized by airway
  enum ActionType
  {
    none,
    calc_sys_matrix_rhs,
    calc_sys_matrix_rhs_iad,
    get_initial_state,
    set_bc,
    calc_flow_rates,
    get_coupled_values,
    get_junction_volume_mix,
    solve_scatra,
    solve_junction_scatra,
    calc_cfl,
    solve_blood_air_transport,
    update_scatra,
    update_elem12_scatra,
    eval_nodal_ess_vals,
    eval_PO2_from_concentration,
    calc_elem_volumes
  };



  //! Element Type
  std::string                                     elemType_;

  //! Resistance Type
  std::string                                     resistance_;

  //! Solver Type
  std::string                                     elemsolvingType_;

  //! container containing all kind of airway thisandthat
  DRT::Container                                  data_;


  // data
  std::map<std::string, double>                   elemParams_;


  // element tree generation
  int                                             generation_;

  // internal calculation methods

  // don't want = operator
  RedAirway& operator = (const RedAirway& old);


  /// set number of gauss points to element shape default
  DRT::UTILS::GaussRule1D getOptimalGaussrule(const DiscretizationType& distype);

  /*!
   * \brief check, whether higher order derivatives for shape functions (dxdx, dxdy, ...) are necessary
   * \return boolean indicating higher order status
   */
  bool isHigherOrderElement(
    const DRT::Element::DiscretizationType  distype  ///< discretization type
    ) const;


}; // class RedAirway


class RedAcinusType : public DRT::ElementType
{
public:

  std::string Name() const { return "RedAcinusType"; }

  static RedAcinusType & Instance();

  virtual DRT::ParObject* Create( const std::vector<char> & data );

  virtual Teuchos::RCP<DRT::Element> Create( const std::string eletype,
                                             const std::string eledistype,
                                             const int id,
                                             const int owner );

  virtual Teuchos::RCP<DRT::Element> Create( const int id, const int owner );

  virtual void NodalBlockInformation( DRT::Element * dwele, int & numdf, int & dimns, int & nv, int & np ) {}

  virtual void ComputeNullSpace( DRT::Discretization & dis, std::vector<double> & ns, const double * x0, int numdf, int dimns ) {}

  virtual void SetupElementDefinition( std::map<std::string,std::map<std::string,DRT::INPUT::LineDefinition> > & definitions );

private:

  static RedAcinusType instance_;
};

/*!
\brief A C++ wrapper for the acinus element

\author ismail (ismail@lnm.mw.tum.de)
*/
class RedAcinus : public DRT::Element
{
public:

  //! @name Friends
  friend class RedAirwayImplInterface;
  friend class AirwayImpl<DRT::Element::line2>;

  //@}
  //! @name Constructors and destructors and related methods

  /*!
  \brief Standard Constructor

  \param id : A unique global id
  */
  RedAcinus(int id, int owner);

  /*!
  \brief Copy Constructor

  Makes a deep copy of a Element

  */
  RedAcinus(const RedAcinus& old);

  /*!
  \brief Deep copy this instance of Redacinus and return pointer to the copy

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
  {if (NumNode()==2) return 1;
      else {
          dserror("Could not determine number of lines");
          return -1;
      }
  }

  /*!
  \brief Return number of surfaces of this element (always 1)
  */
  virtual int NumSurface() const {return -1;}

  /*!
  \brief Return number of volumes of this element (always 1)
  */
  virtual int NumVolume() const {return -1;}

  /*!
  \brief Get vector of Teuchos::RCPs to the lines of this element
  */
  virtual std::vector<Teuchos::RCP<DRT::Element> > Lines();

  /*!
  \brief Return center coordinates of element
  */
  virtual std::vector<double> ElementCenterRefeCoords();

  /*!
  \brief Return unique ParObject id

  every class implementing ParObject needs a unique id defined at the
  top of this file.
  */
  virtual int UniqueParObjectId() const { return RedAcinusType::Instance().UniqueParObjectId(); }

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
  virtual ~RedAcinus();

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
  virtual int NumDofPerNode(const DRT::Node& node) const {return 1;}

  /*!
  \brief Get number of degrees of freedom per element
         (implements pure virtual DRT::Element)

  The element decides how many element degrees of freedom it has.
  It can redecide along the way of a simulation.

  \note Element degrees of freedom mentioned here are dofs that are visible
        at the level of the total system of equations. Purely internal
        element dofs that are condensed internally should NOT be considered.
  */
  virtual int NumDofPerElement() const
    {
      return 0;
    }

  /*!
  \brief Print this element
  */
  void Print(std::ostream& os) const;

  virtual RedAcinusType & ElementType() const
  { return RedAcinusType::Instance(); }

  //@}

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

  Evaluate acinus element stiffness, mass, internal forces etc

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

  An element derived from this class uses the EvaluateNeumann method to receive commands
  and parameters from some control routine in params and evaluates a Neumann boundary condition
  given in condition

  \note This class implements a dummy of this method that prints a warning and
        returns false.

  \param params (in/out)    : ParameterList for communication between control routine
                              and elements
  \param discretization (in): A reference to the underlying discretization
  \param condition (in)     : The condition to be evaluated
  \param lm (in)            : location vector of this element
  \param elevec1 (out)      : Force vector to be filled by element

  \return 0 if successful, negative otherwise
  */
  virtual int EvaluateNeumann(Teuchos::ParameterList& params,
                              DRT::Discretization&      discretization,
                              DRT::Condition&           condition,
                              std::vector<int>&         lm,
                              Epetra_SerialDenseVector& elevec1,
                              Epetra_SerialDenseMatrix* elemat1 = NULL);

  /*!
  \brief Evaluate a Neumann boundary condition

  this method evaluates a line Neumann condition on the acinus element

  \return 0 if successful, negative otherwise
  */
  virtual int EvaluateDirichlet(Teuchos::ParameterList&   params,
                                DRT::Discretization&      discretization,
                                DRT::Condition&           condition,
                                std::vector<int>&         lm,
                                Epetra_SerialDenseVector& elevec1);


  //@}


  //! @name Other
  std::string Type(){return elemType_;}

  //! @name Other
  std::string Resistance(){return resistance_;}

  //@}

  //! get cross-sectional area
  //  double         getA(){ return A_;}

  //! get youngs modulus of the wall
  //  double         getEw(){ return Ew_;}

  //! get youngs modulus of the air
  //  double         getEa(){ return Ea_;}

  //! get wall thickness
  //  double         gettw(){ return tw_;}

  //! set qin
  //  void           setqin (double Qin ){ qin_ = Qin;}

  //! set qout
  //  void           setqout(double Qout){ qout_ = Qout;}

  //! set element parameters
  void           setParams(std::string name, double & var);

  //! get element parameters
  void           getParams(std::string name, double & var);

  //! get element parameters
  void           getParams(std::string name, int    & var);

  const std::map<std::string, double> getAllParams(){return elemParams_;}

private:

  //! action parameters recognized by acinus
  enum ActionType
  {
    none,
    calc_sys_matrix_rhs,
    calc_sys_matrix_rhs_iad,
    get_initial_state,
    set_bc,
    calc_flow_rates,
    get_coupled_values,
    get_junction_volume_mix,
    solve_scatra,
    solve_junction_scatra,
    calc_cfl,
    solve_blood_air_transport,
    update_scatra,
    update_elem12_scatra,
    eval_nodal_ess_vals,
    eval_PO2_from_concentration,
    calc_elem_volumes
  };


  //! Element Type
  std::string                                     elemType_;

  //! Resistance Type
  std::string                                     resistance_;

  //! container containing all kind of airway thisandthat
  DRT::Container                                  data_;


  // data
  std::map<std::string, double>                   elemParams_;


  // element tree generation
  int                                             generation_;

  // internal calculation methods

  // don't want = operator
  RedAcinus& operator = (const RedAcinus& old);


  /// set number of gauss points to element shape default
  DRT::UTILS::GaussRule1D getOptimalGaussrule(const DiscretizationType& distype);

  /*!
   * \brief check, whether higher order derivatives for shape functions (dxdx, dxdy, ...) are necessary
   * \return boolean indicating higher order status
   */
  bool isHigherOrderElement(
    const DRT::Element::DiscretizationType  distype  ///< discretization type
    ) const;


}; // class RedAcinus

/*!
\brief A C++ wrapper for the inter acinar dependency element

\author ismail (ismail@lnm.mw.tum.de)
*/

class RedInterAcinarDepType : public DRT::ElementType
{
public:

  std::string Name() const { return "RedInterAcinarDepType"; }

  static RedInterAcinarDepType & Instance();

  virtual DRT::ParObject* Create( const std::vector<char> & data );

  virtual Teuchos::RCP<DRT::Element> Create( const std::string eletype,
                                             const std::string eledistype,
                                             const int id,
                                             const int owner );

  virtual Teuchos::RCP<DRT::Element> Create( const int id, const int owner );

  virtual void NodalBlockInformation( DRT::Element * dwele, int & numdf, int & dimns, int & nv, int & np ) {}

  virtual void ComputeNullSpace( DRT::Discretization & dis, std::vector<double> & ns, const double * x0, int numdf, int dimns ) {}

  virtual void SetupElementDefinition( std::map<std::string,std::map<std::string,DRT::INPUT::LineDefinition> > & definitions );

private:

  static RedInterAcinarDepType instance_;
};

/*!
\brief A C++ wrapper for the inter acinar dependency element

\author ismail (ismail@lnm.mw.tum.de)
*/
class RedInterAcinarDep : public DRT::Element
{
public:

  //! @name Friends
  friend class RedAirwayImplInterface;
  friend class AirwayImpl<DRT::Element::line2>;

  //@}
  //! @name Constructors and destructors and related methods

  /*!
  \brief Standard Constructor

  \param id : A unique global id
  */
  RedInterAcinarDep(int id, int owner);

  /*!
  \brief Copy Constructor

  Makes a deep copy of a Element

  */
  RedInterAcinarDep(const RedInterAcinarDep& old);

  /*!
  \brief Deep copy this instance of Redacinus and return pointer to the copy

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
  {if (NumNode()==2) return 1;
      else {
          dserror("Could not determine number of lines");
          return -1;
      }
  }

  /*!
  \brief Return number of surfaces of this element (always 1)
  */
  virtual int NumSurface() const {return -1;}

  /*!
  \brief Return number of volumes of this element (always 1)
  */
  virtual int NumVolume() const {return -1;}

  /*!
  \brief Get vector of Teuchos::RCPs to the lines of this element
  */
  virtual std::vector<Teuchos::RCP<DRT::Element> > Lines();

  /*!
  \brief Return unique ParObject id

  every class implementing ParObject needs a unique id defined at the
  top of this file.
  */
  virtual int UniqueParObjectId() const { return RedInterAcinarDepType::Instance().UniqueParObjectId(); }

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
  virtual ~RedInterAcinarDep();

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
  virtual int NumDofPerNode(const DRT::Node& node) const {return 1;}

  /*!
  \brief Get number of degrees of freedom per element
         (implements pure virtual DRT::Element)

  The element decides how many element degrees of freedom it has.
  It can redecide along the way of a simulation.

  \note Element degrees of freedom mentioned here are dofs that are visible
        at the level of the total system of equations. Purely internal
        element dofs that are condensed internally should NOT be considered.
  */
  virtual int NumDofPerElement() const
    {

      return 0;
    }

  /*!
  \brief Print this element
  */
  void Print(std::ostream& os) const;

  virtual RedInterAcinarDepType & ElementType() const
  { return RedInterAcinarDepType::Instance(); }

  //@}

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

  Evaluate inter acinar dependency element stiffness, mass, internal forces etc

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

  An element derived from this class uses the EvaluateNeumann method to receive commands
  and parameters from some control routine in params and evaluates a Neumann boundary condition
  given in condition

  \note This class implements a dummy of this method that prints a warning and
        returns false.

  \param params (in/out)    : ParameterList for communication between control routine
                              and elements
  \param discretization (in): A reference to the underlying discretization
  \param condition (in)     : The condition to be evaluated
  \param lm (in)            : location vector of this element
  \param elevec1 (out)      : Force vector to be filled by element

  \return 0 if successful, negative otherwise
  */
  virtual int EvaluateNeumann(Teuchos::ParameterList& params,
                              DRT::Discretization&      discretization,
                              DRT::Condition&           condition,
                              std::vector<int>&         lm,
                              Epetra_SerialDenseVector& elevec1,
                              Epetra_SerialDenseMatrix* elemat1 = NULL);

  /*!
  \brief Evaluate a Neumann boundary condition

  this method evaluates a line Neumann condition on the acinus element

  \return 0 if successful, negative otherwise
  */
  virtual int EvaluateDirichlet(Teuchos::ParameterList&   params,
                                DRT::Discretization&      discretization,
                                DRT::Condition&           condition,
                                std::vector<int>&         lm,
                                Epetra_SerialDenseVector& elevec1);


  //@}


  //! get cross-sectional area
  //  double         getA(){ return A_;}

  //! get youngs modulus of the wall
  //  double         getEw(){ return Ew_;}

  //! get youngs modulus of the air
  //  double         getEa(){ return Ea_;}

  //! get wall thickness
  //  double         gettw(){ return tw_;}

  //! set qin
  //  void           setqin (double Qin ){ qin_ = Qin;}

  //! set qout
  //  void           setqout(double Qout){ qout_ = Qout;}

  //! get element parameters
  void           getParams(std::string name, double & var);

  //! get element parameters
  void           getParams(std::string name, int    & var);

private:

  //! action parameters recognized by inter acinar dependency
  enum ActionType
  {
    none,
    calc_sys_matrix_rhs,
    calc_sys_matrix_rhs_iad,
    get_initial_state,
    set_bc,
    calc_flow_rates,
    get_coupled_values,
    get_junction_volume_mix,
    solve_scatra,
    calc_cfl,
    solve_blood_air_transport,
    update_scatra,
    update_elem12_scatra,
    eval_nodal_ess_vals,
    eval_PO2_from_concentration,
    calc_elem_volumes
  };

  //! container containing all kind of airway thisandthat
  DRT::Container                                  data_;

  // data
  std::map<std::string, double>                        elemParams_;

  // element tree generation
  int                                             generation_;

  // internal calculation methods

  // don't want = operator
  RedInterAcinarDep& operator = (const RedInterAcinarDep& old);


  /// set number of gauss points to element shape default
  DRT::UTILS::GaussRule1D getOptimalGaussrule(const DiscretizationType& distype);

  /*!
   * \brief check, whether higher order derivatives for shape functions (dxdx, dxdy, ...) are necessary
   * \return boolean indicating higher order status
   */
  bool isHigherOrderElement(
    const DRT::Element::DiscretizationType  distype  ///< discretization type
    ) const;


}; // class RedInterAcinarDep

/*!
\brief A C++ wrapper for the Air/Blood diffusion element

\author ismail (ismail@lnm.mw.tum.de)
*/

class RedAirBloodScatraType : public DRT::ElementType
{
public:

  std::string Name() const { return "RedAirBloodScatraType"; }

  static RedAirBloodScatraType & Instance();

  virtual DRT::ParObject* Create( const std::vector<char> & data );

  virtual Teuchos::RCP<DRT::Element> Create( const std::string eletype,
                                             const std::string eledistype,
                                             const int id,
                                             const int owner );

  virtual Teuchos::RCP<DRT::Element> Create( const int id, const int owner );

  virtual void NodalBlockInformation( DRT::Element * dwele, int & numdf, int & dimns, int & nv, int & np ) {}

  virtual void ComputeNullSpace( DRT::Discretization & dis, std::vector<double> & ns, const double * x0, int numdf, int dimns ) {}

  virtual void SetupElementDefinition( std::map<std::string,std::map<std::string,DRT::INPUT::LineDefinition> > & definitions );

private:

  static RedAirBloodScatraType instance_;
}; // class RedAirBloodScatraType

class RedAirBloodScatra : public DRT::Element
{
public:

  //! @name Friends
  friend class RedAirwayImplInterface;
  friend class AirwayImpl<DRT::Element::line2>;

  //@}
  //! @name Constructors and destructors and related methods

  /*!
  \brief Standard Constructor

  \param id : A unique global id
  */
  RedAirBloodScatra(int id, int owner);

  /*!
  \brief Copy Constructor

  Makes a deep copy of a Element

  */
  RedAirBloodScatra(const RedAirBloodScatra& old);

  /*!
  \brief Deep copy this instance of Redacinus and return pointer to the copy

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
  {if (NumNode()==2) return 1;
      else {
          dserror("Could not determine number of lines");
          return -1;
      }
  }

  /*!
  \brief Return number of surfaces of this element (always 1)
  */
  virtual int NumSurface() const {return -1;}

  /*!
  \brief Return number of volumes of this element (always 1)
  */
  virtual int NumVolume() const {return -1;}

  /*!
  \brief Get vector of Teuchos::RCPs to the lines of this element
  */
  virtual std::vector<Teuchos::RCP<DRT::Element> > Lines();

  /*!
  \brief Return unique ParObject id

  every class implementing ParObject needs a unique id defined at the
  top of this file.
  */
  virtual int UniqueParObjectId() const { return RedAirBloodScatraType::Instance().UniqueParObjectId(); }

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
  virtual ~RedAirBloodScatra();

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
  virtual int NumDofPerNode(const DRT::Node& node) const {return 1;}

  /*!
  \brief Get number of degrees of freedom per element
         (implements pure virtual DRT::Element)

  The element decides how many element degrees of freedom it has.
  It can redecide along the way of a simulation.

  \note Element degrees of freedom mentioned here are dofs that are visible
        at the level of the total system of equations. Purely internal
        element dofs that are condensed internally should NOT be considered.
  */
  virtual int NumDofPerElement() const
    {

      return 0;
    }

  /*!
  \brief Print this element
  */
  void Print(std::ostream& os) const;

  virtual RedAirBloodScatraType & ElementType() const
  { return RedAirBloodScatraType::Instance(); }

  //@}

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
    names.insert(std::pair<string,int>("Owner",1));
    // Name of data is 'StressesXYZ', dimension is 6 (sym. tensor value)
    names.insert(std::pair<string,int>("StressesXYZ",6));
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

  Evaluate inter acinar dependency element stiffness, mass, internal forces etc

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

  An element derived from this class uses the EvaluateNeumann method to receive commands
  and parameters from some control routine in params and evaluates a Neumann boundary condition
  given in condition

  \note This class implements a dummy of this method that prints a warning and
        returns false.

  \param params (in/out)    : ParameterList for communication between control routine
                              and elements
  \param discretization (in): A reference to the underlying discretization
  \param condition (in)     : The condition to be evaluated
  \param lm (in)            : location vector of this element
  \param elevec1 (out)      : Force vector to be filled by element

  \return 0 if successful, negative otherwise
  */
  virtual int EvaluateNeumann(Teuchos::ParameterList& params,
                              DRT::Discretization&      discretization,
                              DRT::Condition&           condition,
                              std::vector<int>&         lm,
                              Epetra_SerialDenseVector& elevec1,
                              Epetra_SerialDenseMatrix* elemat1 = NULL);

  /*!
  \brief Evaluate a Neumann boundary condition

  this method evaluates a line Neumann condition on the acinus element

  \return 0 if successful, negative otherwise
  */
  virtual int EvaluateDirichlet(Teuchos::ParameterList&   params,
                                DRT::Discretization&      discretization,
                                DRT::Condition&           condition,
                                std::vector<int>&         lm,
                                Epetra_SerialDenseVector& elevec1);


  //@}


  //! get cross-sectional area
  //  double         getA(){ return A_;}

  //! get youngs modulus of the wall
  //  double         getEw(){ return Ew_;}

  //! get youngs modulus of the air
  //  double         getEa(){ return Ea_;}

  //! get wall thickness
  //  double         gettw(){ return tw_;}

  //! set qin
  //  void           setqin (double Qin ){ qin_ = Qin;}

  //! set qout
  //  void           setqout(double Qout){ qout_ = Qout;}

  //! get element parameters
  void           getParams(std::string name, double & var);

  //! get element parameters
  void           getParams(std::string name, int    & var);

private:

  //! action parameters recognized by inter acinar dependency
  enum ActionType
  {
    none,
    calc_sys_matrix_rhs,
    calc_sys_matrix_rhs_iad,
    get_initial_state,
    set_bc,
    calc_flow_rates,
    get_coupled_values,
    get_junction_volume_mix,
    solve_scatra,
    calc_cfl,
    solve_junction_scatra,
    solve_blood_air_transport,
    update_scatra,
    update_elem12_scatra,
    eval_nodal_ess_vals,
    eval_PO2_from_concentration,
    calc_elem_volumes
  };

  //! container containing all kind of airway thisandthat
  DRT::Container                                  data_;

  // data
  std::map<std::string, double>                        elemParams_;

  // element tree generation
  int                                             generation_;

  // internal calculation methods

  // don't want = operator
  RedAirBloodScatra  & operator = (const RedAirBloodScatra& old);


  /// set number of gauss points to element shape default
  DRT::UTILS::GaussRule1D getOptimalGaussrule(const DiscretizationType& distype);

  /*!
   * \brief check, whether higher order derivatives for shape functions (dxdx, dxdy, ...) are necessary
   * \return boolean indicating higher order status
   */
  bool isHigherOrderElement(
    const DRT::Element::DiscretizationType  distype  ///< discretization type
    ) const;

}; // class RedAirBloodScatra


/*!
\brief A C++ wrapper for the Air/Blood diffusion element

\author ismail (ismail@lnm.mw.tum.de)
*/
class RedAirBloodScatraLine3Type : public DRT::ElementType
{
public:

  std::string Name() const { return "RedAirBloodScatraLine3Type"; }

  static RedAirBloodScatraLine3Type & Instance();

  virtual DRT::ParObject* Create( const std::vector<char> & data );

  virtual Teuchos::RCP<DRT::Element> Create( const std::string eletype,
                                             const std::string eledistype,
                                             const int id,
                                             const int owner );

  virtual Teuchos::RCP<DRT::Element> Create( const int id, const int owner );

  virtual void NodalBlockInformation( DRT::Element * dwele, int & numdf, int & dimns, int & nv, int & np ) {}

  virtual void ComputeNullSpace( DRT::Discretization & dis, std::vector<double> & ns, const double * x0, int numdf, int dimns ) {}

  virtual void SetupElementDefinition( std::map<std::string,std::map<std::string,DRT::INPUT::LineDefinition> > & definitions );

private:

  static RedAirBloodScatraLine3Type instance_;
}; // class RedAirBloodScatraLine3Type


/*!
\brief A C++ wrapper for the inter acinar dependency element

\author ismail (ismail@lnm.mw.tum.de)
*/
class RedAirBloodScatraLine3 : public DRT::Element
{
public:

  //! @name Friends
  friend class RedAirwayImplInterface;
  friend class AirwayImpl<DRT::Element::line3>;

  //@}
  //! @name Constructors and destructors and related methods

  /*!
  \brief Standard Constructor

  \param id : A unique global id
  */
  RedAirBloodScatraLine3(int id, int owner);

  /*!
  \brief Copy Constructor

  Makes a deep copy of a Element

  */
  RedAirBloodScatraLine3(const RedAirBloodScatraLine3& old);

  /*!
  \brief Deep copy this instance of Redacinus and return pointer to the copy

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
  {if (NumNode()==3)
    {
      return 2;
    }
    else
    {
      dserror("Could not determine number of lines");
      return -1;
    }
  }

  /*!
  \brief Return number of surfaces of this element (always 1)
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
  \brief Return unique ParObject id

  every class implementing ParObject needs a unique id defined at the
  top of this file.
  */
  virtual int UniqueParObjectId() const { return RedAirBloodScatraLine3Type::Instance().UniqueParObjectId(); }

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
  virtual ~RedAirBloodScatraLine3();

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
  virtual int NumDofPerNode(const DRT::Node& node) const {return 1;}

  /*!
  \brief Get number of degrees of freedom per element
         (implements pure virtual DRT::Element)

  The element decides how many element degrees of freedom it has.
  It can redecide along the way of a simulation.

  \note Element degrees of freedom mentioned here are dofs that are visible
        at the level of the total system of equations. Purely internal
        element dofs that are condensed internally should NOT be considered.
  */
  virtual int NumDofPerElement() const
    {

      return 0;
    }

  /*!
  \brief Print this element
  */
  void Print(std::ostream& os) const;

  virtual RedAirBloodScatraLine3Type & ElementType() const
  { return RedAirBloodScatraLine3Type::Instance(); }

  //@}

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
    names.insert(std::pair<string,int>("Owner",1));
    // Name of data is 'StressesXYZ', dimension is 6 (sym. tensor value)
    names.insert(std::pair<string,int>("StressesXYZ",6));
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

  Evaluate inter acinar dependency element stiffness, mass, internal forces etc

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

  An element derived from this class uses the EvaluateNeumann method to receive commands
  and parameters from some control routine in params and evaluates a Neumann boundary condition
  given in condition

  \note This class implements a dummy of this method that prints a warning and
        returns false.

  \param params (in/out)    : ParameterList for communication between control routine
                              and elements
  \param discretization (in): A reference to the underlying discretization
  \param condition (in)     : The condition to be evaluated
  \param lm (in)            : location vector of this element
  \param elevec1 (out)      : Force vector to be filled by element

  \return 0 if successful, negative otherwise
  */
  virtual int EvaluateNeumann(Teuchos::ParameterList& params,
                              DRT::Discretization&      discretization,
                              DRT::Condition&           condition,
                              std::vector<int>&         lm,
                              Epetra_SerialDenseVector& elevec1,
                              Epetra_SerialDenseMatrix* elemat1 = NULL);

  /*!
  \brief Evaluate a Neumann boundary condition

  this method evaluates a line Neumann condition on the acinus element

  \return 0 if successful, negative otherwise
  */
  virtual int EvaluateDirichlet(Teuchos::ParameterList&   params,
                                DRT::Discretization&      discretization,
                                DRT::Condition&           condition,
                                std::vector<int>&         lm,
                                Epetra_SerialDenseVector& elevec1);


  //@}


  //! get cross-sectional area
  //  double         getA(){ return A_;}

  //! get youngs modulus of the wall
  //  double         getEw(){ return Ew_;}

  //! get youngs modulus of the air
  //  double         getEa(){ return Ea_;}

  //! get wall thickness
  //  double         gettw(){ return tw_;}

  //! set qin
  //  void           setqin (double Qin ){ qin_ = Qin;}

  //! set qout
  //  void           setqout(double Qout){ qout_ = Qout;}

  //! get element parameters
  void           getParams(std::string name, double & var);

  //! get element parameters
  void           getParams(std::string name, int    & var);

private:

  //! action parameters recognized by inter acinar dependency
  enum ActionType
  {
    none,
    calc_sys_matrix_rhs,
    calc_sys_matrix_rhs_iad,
    get_initial_state,
    set_bc,
    calc_flow_rates,
    get_coupled_values,
    get_junction_volume_mix,
    solve_scatra,
    calc_cfl,
    solve_junction_scatra,
    solve_blood_air_transport,
    update_scatra,
    update_elem12_scatra,
    eval_nodal_ess_vals,
    eval_PO2_from_concentration,
    calc_elem_volumes
  };

  //! container containing all kind of airway thisandthat
  DRT::Container                                  data_;

  // data
  std::map<std::string, double>                        elemParams_;

  // element tree generation
  int                                             generation_;

  // internal calculation methods

  // don't want = operator
  RedAirBloodScatraLine3  & operator = (const RedAirBloodScatraLine3& old);


  /// set number of gauss points to element shape default
  DRT::UTILS::GaussRule1D getOptimalGaussrule(const DiscretizationType& distype);

  /*!
   * \brief check, whether higher order derivatives for shape functions (dxdx, dxdy, ...) are necessary
   * \return boolean indicating higher order status
   */
  bool isHigherOrderElement(
    const DRT::Element::DiscretizationType  distype  ///< discretization type
    ) const;

}; // class RedAirBloodScatraLine3


} // namespace ELEMENTS
} // namespace DRT


#endif  // #ifndef RED_AIRWAY_H_
