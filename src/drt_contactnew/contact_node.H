/*!----------------------------------------------------------------------
\file contact_node.H

<pre>
-------------------------------------------------------------------------
                        BACI Contact library
            Copyright (2008) Technical University of Munich

Under terms of contract T004.008.000 there is a non-exclusive license for use
of this work by or on behalf of Rolls-Royce Ltd & Co KG, Germany.

This library is proprietary software. It must not be published, distributed,
copied or altered in any form or any media without written permission
of the copyright holder. It may be used under terms and conditions of the
above mentioned license by or on behalf of Rolls-Royce Ltd & Co KG, Germany.

This library contains and makes use of software copyrighted by Sandia Corporation
and distributed under LGPL licence. Licensing does not apply to this or any
other third party software used here.

Questions? Contact Dr. Michael W. Gee (gee@lnm.mw.tum.de)
                   or
                   Prof. Dr. Wolfgang A. Wall (wall@lnm.mw.tum.de)

http://www.lnm.mw.tum.de

-------------------------------------------------------------------------
</pre>

<pre>
Maintainer: Alexander Popp
            popp@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15264
</pre>

*----------------------------------------------------------------------*/
#ifdef CCADISCRET
#ifndef CONTACT_NODE_H
#define CONTACT_NODE_H

#include "../drt_mortar/mortar_node.H"

using namespace std;
using namespace Teuchos;


/*!
\brief CONTACT: namespace of the baci contact and meshtying algorithms

*/
namespace CONTACT
{

/*!
\brief A class for a contact node derived from MORTAR::MortarNode

This class represents a finite element node capable of contact.

\author popp (popp@lnm.mw.tum.de)
*/
class CoNode : public MORTAR::MortarNode
{
public:

  //! @name Enums and Friends

  /*!
  \brief The Discretization is a friend of CoNode
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
  \param isslave(in): flag indicating whether node is slave or master
  \param initactive (in): flag indicating whether initially set to active

  */
  CoNode(int id, const double* coords, const int owner,
         const int numdof, const vector<int>& dofs, const bool isslave,
         const bool initactive);

  /*!
  \brief Copy Constructor

  Makes a deep copy of a CoNode

  */
  CoNode(const CONTACT::CoNode& old);

  /*!
  \brief Deep copy the derived class and return pointer to it

  */
  virtual CONTACT::CoNode* Clone() const;

  /*!
  \brief Destructor

  */
  virtual ~CoNode() { return; }

  /*!
  \brief Return unique ParObject id

  every class implementing ParObject needs a unique id defined at the
  top of drt_lib/drt_parobject.H.

  */
  virtual int UniqueParObjectId() const { return ParObject_CoNode;}

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
  \brief Is Node initialized as active node (only slave nodes)
  */
  virtual bool IsInitActive() const
  { if (!IsSlave()) dserror("ERROR: InitActive requested for Master node");
    return initactive_; }

  /*!
  \brief Modify initial active status of slave node

  This belated modification is necessary to be able to use
  the binary search tree for contact initialization in the
  load-controlled quasistatic case (instead of input file
  information Active/Inactive)

  */
  virtual bool& SetInitActive()
  { if (!IsSlave()) dserror("ERROR: InitActive requested for Master node");
    return initactive_; }

  /*!
  \brief Return current nodal tangent t_xi (only for slave side!) (length 3)
  */
  virtual double* txi() { return txi_; }

  /*!
  \brief Return current nodal tangent t_eta (only for slave side!) (length 3)
  */
  virtual double* teta() { return teta_; }

  /*!
  \brief Return contact status of this node (active=true)
  */
  virtual bool& Active() { return active_; }

  /*!
  \brief Return scaling factor for weighted gap (penalty strategy)
  */
  virtual double& Kappa() { return kappa_; }

  /*!
  \brief Add a value to the weighted gap of this node

  This value is later assembled to the weighted gap vec.
  Note that grow_ here is a scalar.

  \param val : value to be added

  */
  void AddgValue(double& val);
 
  /*!
  \brief Add a value to the map of LM derivatives of this node
     
  The 'DerivZ' map is later assembled to the global DerivZ matrix.
  Note that derivz_ here is a vector.
  
  Note: This is only calculated when performing a penalty strategy
     
  \param row : local dof row id to add to (rowwise)
  \param val : value to be added
  \param col : global dof column id of the value added
  
  */
  void AddDerivZValue(int& row, const int& col, double val);

  /*!
  \brief Return the weighted gap (scalar) of this node
  */
  virtual double& Getg() { return grow_; }

  /*!
  \brief Return the 'DerivN' map (vector) of this node

  These maps contain the directional derivatives of the node's
  averaged normal vector with respect to the slave displacements.
  A vector is used because the normal itself is a vector (2 or 3 components).

  */
  virtual vector<map<int,double> >& GetDerivN() { return derivn_; }

  /*!
  \brief Return the 'DerivTxi' map (vector) of this node

  These maps contain the directional derivatives of the node's
  unit tangent vector t_xi with respect to the slave displacements.
  A vector is used because the tangent itself is a vector (2 or 3 components).

  */
  virtual vector<map<int,double> >& GetDerivTxi() { return derivtxi_; }

  /*!
  \brief Return the 'DerivTeta' map (vector) of this node

  These maps contain the directional derivatives of the node's
  unit tangent vector t_eta with respect to the slave displacements.
  A vector is used because the tangent itself is a vector (2 or 3 components).

  */
  virtual vector<map<int,double> >& GetDerivTeta() { return derivteta_; }

  /*!
  \brief Return the 'DerivD' map of this node

  D-matrix entries with respect to the slave/master displacements.
  It is a "map of maps", with the outer map containing all master
  node indices k adjacent to this node i and the inner map containing
  all directional derivatives l existing for D_ik.
  
  Note: In earlier versions this was just a simple map accounting for the
  diagonality of D when using dual shape functions. Due to the extension to
  support arbitrary types of shape functions, this is not possible anymore.

  */
  virtual map<int,map<int,double> >& GetDerivD() { return derivd_; }

  /*!
  \brief Return the 'DerivM' map of maps of this node

  This map contains the directional derivatives of the node's
  M-matrix entries with respect to the slave/master displacements.
  It is a "map of maps", with the outer map containing all master
  node indices k adjacent to this node i and the inner map containing
  all directional derivatives l existing for M_ik.

  */
  virtual map<int,map<int,double> >& GetDerivM() { return derivm_; }

  /*!
  \brief Return one specific 'DerivD' map of this node

  This method returns the map of directional derivatives of one
  specific D-matrix D_ik entry of this node i.

  */
  virtual map<int,double>& GetDerivD(int& k)
  {
    typedef map<int,map<int,double> >::const_iterator CI;
    CI p =derivd_.find(k);
    if (p==derivd_.end())
     dserror("ERROR: GetDerivD: No map entry existing for given index");
    return derivd_[k];
  }
  
  /*!
  \brief Return one specific 'DerivM' map of this node

  This method returns the map of directional derivatives of one
  specific M-matrix M_ik entry of this node i.

  */
  virtual map<int,double>& GetDerivM(int& k)
  {
    typedef map<int,map<int,double> >::const_iterator CI;
    CI p =derivm_.find(k);
    if (p==derivm_.end())
     dserror("ERROR: GetDerivM: No map entry existing for given index");
    return derivm_[k];
  }

  /*!
  \brief Return the 'DerivG' map of this node

  This map contains the directional derivatives of the node's
  weighted gap entry g~ with respect to the slave/master displacements.

  */
  virtual map<int,double>& GetDerivG() { return derivg_; }

  /*!
  \brief Return the 'DerivZ map of this node
  
  This map contains the directional derivatives of the node's
  lagrange multiplier entries with respect to the slave/master displacements.

  Note: This is only calculated when performing a penalty strategy

   */
  virtual vector<map<int,double> >& GetDerivZ() { return derivz_; }
    
  //@}

  //! @name Construction methods

  /*!
  \brief Build nodal normal
  */
  void BuildAveragedNormal();

  /*!
  \brief Build directional derivative of nodal normal + tangents

  This method will be called after having finished the method
  BuildAveragedNormal() for the computation of nodal normals. The result
  (directional derivatives) will be stored in the nodal maps derivn_ and
  will later be assembled to the global system of equations.
  Please note that we also compute the directional derivative of the
  nodal tangent here,and the results will be stored analogously in
  the nodal maps derivtxi_ (2D and 3D) and derivteta_ (only 3D).

  \param elens (in):  Matrix containing normals of adjacent elements
  \param length (in): Length of the nodal averaged normal
  \param ltxi (in):   Length of the nodal tangent txi

  */
  void DerivAveragedNormal(Epetra_SerialDenseMatrix& elens,
                           double length, double ltxi);
  
  //@}
  
protected:

  bool            active_;         // true if cnode is in contact (active set strategy)
  bool            initactive_;     // true if node is initialized as active node

  double          txi_[3];         // nodal tangent t_xi for contact methods
  double          teta_[3];        // nodal tangent t_eta for contact methods

  double          grow_;           // nodal entry of g vector

  vector<map<int,double> >      derivn_;     // directional derivative of nodal normal
  vector<map<int,double> >      derivtxi_;   // directional derivative of nodal tangent t_xi
  vector<map<int,double> >      derivteta_;  // directional derivative of nodal tangent t_eta
  map<int,map<int,double> >     derivd_;     // directional derivative of nodal D-matrix value
  map<int,map<int,double> >     derivm_;     // directional derivative of nodal M-matrix values
  map<int,double>               derivg_;     // directional derivative of nodal weighted gap value

  // Penalty-related quantities
  vector<map<int,double> >      derivz_;     // direction derivative of nodal z-matrix value

}; // class CoNode
} // namespace CONTACT


// << operator
ostream& operator << (ostream& os, const CONTACT::CoNode& cnode);

#endif  // #ifndef CNODE_H
#endif  // #ifdef CCADISCRET
