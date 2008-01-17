/*!----------------------------------------------------------------------
\file drt_celement.H

\class CONTACT::CElement

\brief A contact element

<pre>
Maintainer: Michael Gee
            gee@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15239
</pre>

*----------------------------------------------------------------------*/
#ifdef CCADISCRET
#ifndef DRT_CELEMENT_H
#define DRT_CELEMENT_H


#include "../drt_lib/drt_element.H"
#include "../drt_lib/drt_elementregister.H"
#include "drt_cnode.H"
#include "../drt_lib/linalg_serialdensematrix.H"

using namespace std;
using namespace Teuchos;

// forward declarations

/*!
\brief CONTACT: namespace of the baci contact algorithms

*/
namespace CONTACT
{

// forward declarations

/*!
\class CElement

\brief A contact element


\author gee (gee@lnm.mw.tum.de), popp (popp@lnm.mw.tum.de)
*/
class CElement : public DRT::Element
{
public:

  //! @name Enums and Friends
  //@}

  //! @name Constructors and destructors and related methods

  /*!
  \brief Standard Constructor

  \param id    (in): A globally unique element id
  \param etype (in): Type of element
  \param owner (in): owner processor of the element
  \param shape (in): shape of this element
  \param numnode (in): Number of nodes to this element
  \param nodeids (in): ids of nodes adjacent to this element
  \param isslave (in): flag indicating whether element is slave or master side
  */
  explicit CElement(int id, ElementType etype, int owner, 
                    const DRT::Element::DiscretizationType& shape, 
                    const int numnode,
                    const int* nodeids, 
                    const bool isslave);

  /*!
  \brief Copy Constructor

  Makes a deep copy of this class

  */
  CElement(const CONTACT::CElement& old);

  /*!
  \brief Destructor

  */
  virtual ~CElement() { return; }

  /*!
  \brief Deep copy the derived class and return pointer to it

  */
  virtual CONTACT::CElement* Clone() const;


  /*!
  \brief Return unique ParObject id

  Every class implementing ParObject needs a unique id defined at the
  top of drt_parobject.H
  */
  virtual int UniqueParObjectId() const { return ParObject_CElement; }

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
  \brief Allocate and return a ElementRegister class

   Assuming the contact elements do not need initialization, we return a 
   dummy base class here.

  */
  virtual RefCountPtr<DRT::ElementRegister> ElementRegister() const
  { return rcp(new DRT::ElementRegister(Type())); }

  //@}

  //! @name Query methods


  /*!
  \brief Get shape type of element
  */
  virtual DiscretizationType Shape() const { return shape_;}

  /*!
  \brief Return number of lines to this element
  */
  virtual int NumLine() const {return 0;}

  /*!
  \brief Return number of surfaces to this element
  */
  virtual int NumSurface() const {return 0;}


  /*!
  \brief Get vector of ptrs to the lines of this element
  */
  virtual Element** Lines() { return NULL; }

  /*!
  \brief Get vector of ptrs to the surfaces of this element
  */
  virtual Element** Surfaces() { return NULL; }


  /*!
  \brief Get number of degrees of freedom of a certain node
  
  This CElement is picky: It cooperates only with CNodes, not with
  standard Node objects!
  
  */
  virtual int NumDofPerNode(const DRT::Node& node) const
  {
    const CONTACT::CNode* cnode = dynamic_cast<const CONTACT::CNode*>(&node);
    if (!cnode) dserror("Node is not a CNode");
    return cnode->NumDof();
  }

  /*!
  \brief Get number of degrees of freedom per element

  For now contact elements do not have degrees of freedom independent of
  the nodes
  */
  virtual int NumDofPerElement() const { return 0; }

  /*!
  \brief Print this element

  */
  virtual void Print(ostream& os) const;
  
  /*!
  \brief Return reference area

  */
  virtual inline double& RefArea() { return refarea_; }
  
  /*!
  \brief Return current area

  */
  virtual inline double& Area() { return area_; }
  
  /*!
  \brief Return slave (true) or master status

  */
  virtual inline bool IsSlave() { return isslave_; }
  
  /*!
  \brief Return number of potentially contacting elements
  
  */
  virtual inline int NumSearchElements() const { return (int)searchelements_.size(); }
  
  /*!
  \brief Return global ids of potentially contacting elements
  
  */
  virtual inline vector<int>& SearchElements() {return searchelements_; }
  
  //@}

  //! @name Construction


  //@}

  //! @name Conditions

  //@}

  //! @name Evaluation methods



  /*!
  \brief Evaluate an element

  An element derived from this class uses the Evaluate method to receive commands
  and parameters from some control routine in params and evaluates element matrices and
  vectors accoring to the command in params.

  \note This class implements a dummy of this method that prints a warning and
        returns false.

  \param params (in/out)    : ParameterList for communication between control routine
                              and elements
  \param discretization (in): A reference to the underlying discretization
  \param lm (in)            : location vector of this element
  \param elemat1 (out)      : matrix to be filled by element depending on commands
                              given in params
  \param elemat2 (out)      : matrix to be filled by element depending on commands
                              given in params
  \param elevec1 (out)      : vector to be filled by element depending on commands
                              given in params
  \param elevec2 (out)      : vector to be filled by element depending on commands
                              given in params
  \param elevec3 (out)      : vector to be filled by element depending on commands
                              given in params
  \return 0 if successful, negative otherwise
  */
  virtual int Evaluate(ParameterList&            params,
                       DRT::Discretization&      discretization,
                       vector<int>&              lm,
                       Epetra_SerialDenseMatrix& elemat1,
                       Epetra_SerialDenseMatrix& elemat2,
                       Epetra_SerialDenseVector& elevec1,
                       Epetra_SerialDenseVector& elevec2,
                       Epetra_SerialDenseVector& elevec3);

  /*!
  \brief Get local coordinates for local node id
  */
  virtual bool LocalCoordinatesOfNode(int lid, double* xi);
    
  /*!
  \brief Get local numbering for global node id
  */
  virtual int GetLocalNodeId(int nid);
  
  /*!
  \brief Build element normal at node passed in
  */
  virtual void BuildNormalAtNode(int nid, vector<double>& n);
  
  /*!
  \brief Compute element normal at local coordinate xi
  */
  virtual void ComputeNormalAtXi(double* xi, vector<double>& n);
  
  /*!
  \brief Compute element normal at local coordinate xi
  */
  virtual LINALG::SerialDenseMatrix GetNodalCoords();
  	
  /*!
  \brief Evaluate Jacobian det for parameter space integration - LINEAR / QUAD 1D
  */
  virtual double Jacobian_1D(const vector<double>& val,
  													 const vector<double>& deriv,
  													 const LINALG::SerialDenseMatrix& coord);
  
  /*!
  \brief Compute length (in 3D area) of the element
  */
  virtual double ComputeArea();
  
  /*!
  \brief Evaluate shape functions and derivatives - LINEAR / QUAD 1D 
  */
  virtual bool EvaluateShape_1D(const double* xi, vector<double>& val,
  															vector<double>& deriv, const int valdim);
  
  /*!
  \brief Evaluate dual shape functions and derivatives - LINEAR / QUAD 1D 
  */
  virtual bool EvaluateShape_Dual1D(const double* xi, vector<double>& val,
  																 	vector<double>& deriv, const int valdim);
  
  /*!
  \brief Interpolate global coordinates for given local element coordinates
  
  This method interpolates global coordinates for a given local element
  coordinate variable using the element node coordinates. For interpolation
  one can choose between shape functions or shape function derivatives!
  
  \param xi (in)				: local element coordinates
  \param inttype (in) 	: set to true for shape function usage,
  												set to false for derivative usage
  \param globccord (out): interpolated global coordinates
  */
  virtual bool LocalToGlobal(const double* xi, double* globcoord, bool inttype);

  /*!
  \brief Add CElements to this CElement's potentially contact partners
  */
  virtual bool AddSearchElements(const vector<int>& gids);
  
  //@}


private:

  DRT::Element::DiscretizationType shape_;   // shape of this element
  
  bool               isslave_; 							 // indicating slave or master side
  double						 refarea_; 							 // element length or area in reference configuration
  double						 area_;	   							 // element length or area in current configuration
  
  vector<int>		  	 searchelements_; 			 // global ids of potentially contacting elements 

}; // class CElement
} // namespace CONTACT


// << operator
ostream& operator << (ostream& os, const CONTACT::CElement& ele);









#endif  // #ifndef DRT_CELEMENT_H
#endif  // #ifdef CCADISCRET
