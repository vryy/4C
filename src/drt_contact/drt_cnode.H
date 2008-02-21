/*!----------------------------------------------------------------------
\file drt_cnode.H

\class CONTACT::CNode

\brief A class for a contact node

This class represents a finite element node capable of contact.

<pre>
Maintainer: Michael Gee
            gee@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15239
</pre>

*----------------------------------------------------------------------*/
#ifdef CCADISCRET
#ifndef CNODE_H
#define CNODE_H


#include "../drt_lib/drt_node.H"
#include "../drt_lib/drt_discret.H"

using namespace std;
using namespace Teuchos;


/*!
\brief CONTACT: namespace of the baci contact algorithms

*/
namespace CONTACT
{

/*!
\class CNode

\brief A class for a contact node derived from DRT::Node

\author gee (gee@lnm.mw.tum.de), popp (popp@lnm.mw.tum.de)
*/
class CNode : public DRT::Node
{
public:

  //! @name Enums and Friends

  /*!
  \brief The Discretization is a friend of CNode
  */
  friend class DRT::Discretization;


  //@}

  //! @name Constructors and destructors and related methods

  /*!
  \brief Standard Constructor

  \param id     (in): A globally unique node id
  \param coords (in): vector of nodal coordinates, length 3
  \param owner  (in): Owner of this node.
  \param numdof (in): Number of degrees of freedom to this node
  \param dofs   (in): list of global degrees of freedom
  \param isslave(in): flag indicating whether node is slave or master side
  */
  CNode(int id, const double* coords, const int owner, 
        const int numdof, const vector<int>& dofs, const bool isslave);

  /*!
  \brief Copy Constructor

  Makes a deep copy of a CNode

  */
  CNode(const CONTACT::CNode& old);

  /*!
  \brief Deep copy the derived class and return pointer to it

  */
  virtual CONTACT::CNode* Clone() const;

  /*!
  \brief Destructor

  */
  virtual ~CNode() { return; }

  /*!
  \brief Return unique ParObject id

  every class implementing ParObject needs a unique id defined at the
  top of drt_lib/drt_parobject.H.
  */
  virtual int UniqueParObjectId() const { return ParObject_CNode;}

  /*!
  \brief Pack this class so it can be communicated

  \ref Pack and \ref Unpack are used to communicate this node

  */
  virtual void Pack(vector<char>& data) const;

  /*!
  \brief Unpack data from a char vector into this class

  \ref Pack and \ref Unpack are used to communicate this node

  */
  virtual void Unpack(const vector<char>& data);

  //@}

  //! @name Access methods

  /*!
  \brief Print this cnode
  */
  virtual void Print(ostream& os) const ;

  /*!
  \brief Is Node on slave or master side of contact interface
  */
  inline virtual bool IsSlave() const { return isslave_; }

  /*!
  \brief Get no. of degrees of freedom
  */
  inline virtual int NumDof() const { return numdof_; }

  /*!
  \brief Get predefined degrees of freedom
  */
  inline virtual const int* Dofs() const { return &dofs_[0]; }

  /*!
  \brief Return current configuration (length 3)
  */
  virtual inline double* xspatial() { return xspatial_; }
  
  /*!
  \brief Return current displacement (length 3)
  */
  virtual inline double* u() { return u_; }
  
  /*!
  \brief Return current nodal normal (only for slave side!) (length 3)
  */
  virtual inline double* n() { return n_; }
  
  /*!
  \brief Return current closest node's global id (by reference)
  */
  virtual inline int& ClosestNode() { return closestnode_; }
  
  /*!
  \brief Return projection status of this node (slave nodes only)
  */
  virtual inline bool& HasProj() { return hasproj_; }
    
  /*!
  \brief Return contact status of this node (active=true)
  */
  virtual inline bool& Active() { return active_; }
  
  /*!
  brief Add a value to the 'D' map of this node
  
  The 'D' map is later assembled to the D matrix.
  Note that drows_ here is a vector.
    
  \param row : local dof row id to add to (rowwise)
  \param val : value to be added
  \param col : global dof column id of the value added
  */
  void AddDValue(int row, int col, double val);
  
  /*!
  brief Add a value to the 'M' map of this node
    
  The 'M' map is later assembled to the M matrix.
  Note that mrows_ here is a vector.
      
  \param row : local dof row id to add to (rowwise)
  \param val : value to be added
  \param col : global dof column id of the value added
  */
  void AddMValue(int row, int col, double val);
    
  /*!
  brief Add a value to the 'Mmod' map of this node
     
  The 'Mmod' map is later assembled to the M matrix.
  Note that mmodrows_ here is a vector.
       
  \param row : local dof row id to add to (rowwise)
  \param val : value to be added
  \param col : global dof column id of the value added
  */
  void AddMmodValue(int row, int col, double val);  

  /*!
  brief Add a value to the weighted gap of this node
    
  This value is later assembled to the weighted gap vec.
  Note that grow_ here is a scalar.
     
  \param val : value to be added
  */
  void AddgValue(double val);
  
  /*!
  \brief Return the 'D' map (vector) of this node
  */
  virtual inline vector<map<int,double> >& GetD() { return drows_; }
  
  /*!
  \brief Return the 'M' map (vector) of this node
  */
  virtual inline vector<map<int,double> >& GetM() { return mrows_; }
    
  /*!
  \brief Return the 'Mmod' map (vector) of this node
  */
  virtual inline vector<map<int,double> >& GetMmod() { return mmodrows_; }
  
  /*!
  \brief Return the weighted gap (scalar) of this node
  */
  virtual inline double& Getg() { return grow_; }
  
  //@}

  //! @name Construction methods

  /*!
  \brief Build nodal normal
  */
  virtual void BuildAveragedNormal();
  
  /*!
  \brief Find closest node from given node set and return global ID
    
  This method will compute the distance of the active CNode to all
  nodes of the given Epetra_Map on the given DRT::Discretization
    
  \param intdis (in):         Master Node to project
  \param nodesearchmap (in) : Slave Celement to project on
  \param xi (out) : Local coordinates of projection on Celement
  */
  virtual CONTACT::CNode* FindClosestNode(const RCP<DRT::Discretization> intdis,
  														            const RCP<Epetra_Map> nodesearchmap,
  														            double& mindist);
    
  //@}

protected:

  bool            isslave_;        // true if cnode is on slave side of contact problem
  int             numdof_;         // number of degrees of freedom (length of dofs_)
  vector<int>     dofs_;           // degrees of freedom of this node
  double          n_[3];				   // nodal normal for contact methods
  double					xspatial_[3];	   // current configuration
  double					u_[3];				   // current displacement
  int							closestnode_;	   // closest node id on the other side of interface
  bool						hasproj_;				 // true if feasible projection within search set (slave nodes only)
  bool					  active_;				 // true if cnode is in contact (active set strategy)
  
  vector<map<int,double> >     drows_;   		  // nodal rows of D matrix
  vector<map<int,double> > 		 mrows_;   			// nodal rows of M matrix
  vector<map<int,double> >		 mmodrows_;   	// nodal rows of Mmod matrix
  double											 grow_;				  // nodal entry of g vector
  

}; // class CNode
} // namespace CONTACT


// << operator
ostream& operator << (ostream& os, const CONTACT::CNode& cnode);









#endif  // #ifndef CNODE_H
#endif  // #ifdef CCADISCRET
