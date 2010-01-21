/*!----------------------------------------------------------------------
\file drt_cnodedatacontainer.H

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
Maintainer: Markus Gitterle
            gee@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15251
</pre>

*----------------------------------------------------------------------*/
#ifdef CCADISCRET
#ifndef CNODEDATACONTAINER_H
#define CNODEDATACONTAINER_H


using namespace std;
using namespace Teuchos;


/*!
\brief CONTACT: namespace of the baci contact algorithms

*/
namespace CONTACT
{

/*!
\brief A class containing additional data from contact nodes

This class contains additional information from contact nodes which are
are not needed for contact search and therefore are only available on the
node's processor (ColMap).

\author gitterle (gitterle@lnm.mw.tum.de)
*/
class CNodeDataContainer
{
public:

	//@}

  //! @name Constructors and destructors and related methods

  /*!
  \brief Standard Constructor

  */
  CNodeDataContainer();

  /*!
  \brief Copy Constructor

  Makes a deep copy of a CNodeDataContainer

  */
  CNodeDataContainer(const CONTACT::CNodeDataContainer& old);

  /*!
  \brief Deep copy the derived class and return pointer to it

  */
  virtual CONTACT::CNodeDataContainer* Clone() const;

  /*!
  \brief Destructor

  */
  virtual ~CNodeDataContainer() { return; }

  /*!
  \brief Pack this class so that it can be communicated

  This function packs the datacontainer. This is only called
  when the class has been initialized and the pointer to this
  class exists.

  */
  virtual void Pack(vector<char>& data) const;

  /*!
  \brief Unpack data from a vector into this class

  This function unpacks the datacontainer. This is only called
  when the class has been initialized and the pointer to this
  class exists.

  */
  virtual void Unpack(int& position, const vector<char>& data);

  //@}

  //! @name Access methods

  /*!
  \brief Return current nodal normal (only for slave side!) (length 3)
  */
  virtual inline double* n() { return n_; }

  /*!
  \brief Return current nodal tangent t_xi (only for slave side!) (length 3)
  */
  virtual inline double* txi() { return txi_; }

  /*!
  \brief Return current nodal tangent t_eta (only for slave side!) (length 3)
  */
  virtual inline double* teta() { return teta_; }

  /*!
  \brief Return current Lagrange mutliplier in step n+1 (only for slave side!) (length 3)
  */
  virtual inline double* lm() { return lm_; }

  /*!
  \brief Return old Lagrange mutliplier from step n (only for slave side!) (length 3)
  */
  virtual inline double* lmold() { return lmold_; }

  /*!
  \brief Return Lagrange mutliplier from last Uzawa step (only for slave side!) (length 3)
  */
  virtual inline double* lmuzawa() { return lmuzawa_; }

  /*!
  \brief Return jump per time step (only for slave side!) (length 3)
  */
  virtual inline double* jump() { return jump_; }

  /*!
  \brief Return contact status of last converged state n (active=true)
  */
  virtual inline bool& ActiveOld() { return activeold_; }

  /*!
  \brief Return friction status of this node (slip=true)
  */
  virtual inline bool& Slip() { return slip_; }

  /*!
  \brief Return the 'D' map (vector) of this node
  */
  virtual inline vector<map<int,double> >& GetD() { return drows_; }

  /*!
  \brief Return the 'M' map (vector) of this node
  */
  virtual inline vector<map<int,double> >& GetM() { return mrows_; }

  /*!
  \brief Return the old 'D' map (vector) of this node (last converged state)
  */
  virtual inline vector<map<int,double> >& GetDOld() { return drowsold_; }

  /*!
  \brief Return the old 'M' map (vector) of this node (last converged state)
  */
  virtual inline vector<map<int,double> >& GetMOld() { return mrowsold_; }

  /*!
  \brief Return the 'D' map (vector) of this node (Petrov-Galerkin-approach)
  */
  virtual inline vector<map<int,double> >& GetDPG() { return drowsPG_; }

  /*!
  \brief Return the 'M' map (vector) of this node (Petrov-Galerkin-approach)
  */
  virtual inline vector<map<int,double> >& GetMPG() { return mrowsPG_; }

  /*!
  \brief Return the old 'D' map (vector) of this node (last converged state) (Petrov-Galerkin-approach)
  */
  virtual inline vector<map<int,double> >& GetDOldPG() { return drowsoldPG_; }

  /*!
  \brief Return the old 'M' map (vector) of this node (last converged state) (Petrov-Galerkin-approach)
  */
  virtual inline vector<map<int,double> >& GetMOldPG() { return mrowsoldPG_; }

  /*!
  \brief Return the 'Mmod' map (vector) of this node
  */
  virtual inline vector<map<int,double> >& GetMmod() { return mmodrows_; }

  /*!
  \brief Return the map with according master nodes
  */
  virtual inline set<int>& GetSNodes() { return snodes_; }

  /*!
  \brief Return the map with according master nodes
  */
  virtual inline set<int>& GetMNodes() { return mnodes_; }

  /*!
  \brief Return the old map with according master nodes
  */
  virtual inline set<int>& GetMNodesOld() { return mnodesold_; }

  /*!
  \brief Return the weighted gap (scalar) of this node
  */
  virtual inline double& Getg() { return grow_; }

  /*!
  \brief Return the 'DerivN' map (vector) of this node

  These maps contain the directional derivatives of the node's
  averaged normal vector with respect to the slave displacements.
  A vector is used because the normal itself is a vector (2 or 3 components).

  */
  virtual inline vector<map<int,double> >& GetDerivN() { return derivn_; }

  /*!
  \brief Return the 'DerivTxi' map (vector) of this node

  These maps contain the directional derivatives of the node's
  unit tangent vector t_xi with respect to the slave displacements.
  A vector is used because the tangent itself is a vector (2 or 3 components).

  */
  virtual inline vector<map<int,double> >& GetDerivTxi() { return derivtxi_; }

  /*!
  \brief Return the 'DerivTeta' map (vector) of this node

  These maps contain the directional derivatives of the node's
  unit tangent vector t_eta with respect to the slave displacements.
  A vector is used because the tangent itself is a vector (2 or 3 components).

  */
  virtual inline vector<map<int,double> >& GetDerivTeta() { return derivteta_; }

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
  virtual inline map<int,map<int,double> >& GetDerivD() { return derivd_; }

  /*!
  \brief Return one specific 'DerivD' map of this node

  This method returns the map of directional derivatives of one
  specific D-matrix D_ik entry of this node i.

  */
  virtual inline map<int,double>& GetDerivD(int& k)
  {
    typedef map<int,map<int,double> >::const_iterator CI;
    CI p =derivd_.find(k);
    if (p==derivd_.end())
     dserror("ERROR: GetDerivD: No map entry existing for given index");
    return derivd_[k];
  }

  /*!
  \brief Return the 'DerivM' map of maps of this node

  This map contains the directional derivatives of the node's
  M-matrix entries with respect to the slave/master displacements.
  It is a "map of maps", with the outer map containing all master
  node indices k adjacent to this node i and the inner map containing
  all directional derivatives l existing for M_ik.

  */
  virtual inline map<int,map<int,double> >& GetDerivM() { return derivm_; }

  /*!
  \brief Return one specific 'DerivM' map of this node

  This method returns the map of directional derivatives of one
  specific M-matrix M_ik entry of this node i.

  */
  virtual inline map<int,double>& GetDerivM(int& k)
  {
    typedef map<int,map<int,double> >::const_iterator CI;
    CI p =derivm_.find(k);
    if (p==derivm_.end())
     dserror("ERROR: GetDerivM: No map entry existing for given index");
    return derivm_[k];
  }

  /*!
  \brief Return the 'DerivD' map of this node (Petrov-Galerkin-approach)

  D-matrix entries with respect to the slave/master displacements.
  It is a "map of maps", with the outer map containing all master
  node indices k adjacent to this node i and the inner map containing
  all directional derivatives l existing for D_ik.

  Note: In earlier versions this was just a simple map accounting for the
  diagonality of D when using dual shape functions. Due to the extension to
  support arbitrary types of shape functions, this is not possible anymore.

  */
  virtual inline map<int,map<int,double> >& GetDerivDPG() { return derivdPG_; }

  /*!
  \brief Return one specific 'DerivD' map of this node (Petrov-Galerkin-approach)

  This method returns the map of directional derivatives of one
  specific D-matrix D_ik entry of this node i.

  */
  virtual inline map<int,double>& GetDerivDPG(int& k)
  {
    typedef map<int,map<int,double> >::const_iterator CI;
    CI p =derivdPG_.find(k);
    if (p==derivdPG_.end())
     dserror("ERROR: GetDerivDPG: No map entry existing for given index");
    return derivdPG_[k];
  }

  /*!
  \brief Return the 'DerivM' map of maps of this node (Petrov-Galerkin-approach)

  This map contains the directional derivatives of the node's
  M-matrix entries with respect to the slave/master displacements.
  It is a "map of maps", with the outer map containing all master
  node indices k adjacent to this node i and the inner map containing
  all directional derivatives l existing for M_ik.

  */
  virtual inline map<int,map<int,double> >& GetDerivMPG() { return derivmPG_; }

  /*!
  \brief Return one specific 'DerivM' map of this node (Petrov-Galerkin-approach)

  This method returns the map of directional derivatives of one
  specific M-matrix M_ik entry of this node i.

  */
  virtual inline map<int,double>& GetDerivMPG(int& k)
  {
    typedef map<int,map<int,double> >::const_iterator CI;
    CI p =derivmPG_.find(k);
    if (p==derivmPG_.end())
     dserror("ERROR: GetDerivM: No map entry existing for given index");
    return derivmPG_[k];
  }

  /*!
  \brief Return the 'DerivG' map of this node

  This map contains the directional derivatives of the node's
  weighted gap entry g~ with respect to the slave/master displacements.

  */
  virtual inline map<int,double>& GetDerivG() { return derivg_; }

  /*!
  \brief Return the 'DerivJump' map (vector) of this node

  These maps contain the directional derivatives of the node's
  relative movement (jump).
  A vector is used because the jump itself is a vector (2 or 3 components).

  */
  virtual inline vector<map<int,double> >& GetDerivJump() { return derivjump_; }

  /*!
  \brief Return current penalty traction (length 3)
  */
  virtual inline double* traction() { return traction_; }

  /*!
  \brief Return old penalty traction (length 3)
  */
  virtual inline double* tractionold() { return tractionold_; }

  /*!
  \brief Return scaling factor for weighted gap (penalty strategy)
  */
  virtual inline double& Kappa() { return kappa_; }

  /*!
  \brief Return the 'DerivZ map of this node

  This map contains the directional derivatives of the node's
  lagrange multiplier entries with respect to the slave/master displacements.

  Note: This is only calculated when performing a penalty strategy

   */
  virtual inline vector<map<int,double> >& GetDerivZ() { return derivz_; }

  //@}

  //! @name Construction methods

  //@}

protected:

  double          n_[3];           // nodal normal for contact methods
  double          txi_[3];         // nodal tangent t_xi for contact methods
  double          teta_[3];        // nodal tangent t_eta for contact methods
  double          lm_[3];          // current Lagrange multiplier value (n+1)
  double          lmold_[3];       // old Lagrange multiplier value (last converged state n)
	double          lmuzawa_[3];     // Uzawa Lagrange multiplier value (last Uzawa step k)
  double          jump_[3];        // jump per time step

  bool            activeold_;      // true if cnode was in contact (last converged state n)
  bool            slip_;           // true if cnode slips (contact force reaches friction bound)

  vector<map<int,double> >      drows_;          // nodal rows of D matrix
	vector<map<int,double> >      mrows_;          // nodal rows of M matrix
  vector<map<int,double> >      drowsold_;       // nodal rows of old D matrix
  vector<map<int,double> >      mrowsold_;       // nodal rows of old M matrix
  vector<map<int,double> >      mmodrows_;       // nodal rows of Mmod matrix
  vector<map<int,double> >      drowsPG_;        // nodal rows of D matrix with petrov galerkin approach
  vector<map<int,double> >      mrowsPG_;        // nodal rows of M matrix with petrov galerkin approach
  vector<map<int,double> >      drowsoldPG_;     // nodal rows of old D matrix with petrov galerkin approach
  vector<map<int,double> >      mrowsoldPG_;     // nodal rows of old M matrix with petrov galerkin approach
  set<int>                      snodes_;         // nodal set of according slave nodes
  set<int>                      mnodes_;         // nodal set of according master nodes
  set<int>                      mnodesold_;      // nodal set of old according masternodes
  double                        grow_;           // nodal entry of g vector

  vector<map<int,double> >      derivn_;     // directional derivative of nodal normal
  vector<map<int,double> >      derivtxi_;   // directional derivative of nodal tangent t_xi
  vector<map<int,double> >      derivteta_;  // directional derivative of nodal tangent t_eta
  map<int,map<int,double> >     derivd_;     // directional derivative of nodal D-matrix value
  map<int,map<int,double> >     derivm_;     // directional derivative of nodal M-matrix values
  map<int,map<int,double> >     derivdPG_; // directional derivative of nodal D-matrix value
  map<int,map<int,double> >     derivmPG_; // directional derivative of nodal M-matrix values
  map<int,double>               derivg_;     // directional derivative of nodal weighted gap value
  vector<map<int,double> >      derivjump_;  // directional derivative of nodal weighted jump vector

  // Penalty-related quantities
  double                        traction_[3];     // traction vector of current time step
  double                        tractionold_[3];  // traction vector of previous time step
  double                        kappa_;           // gap-scaling factor
  vector<map<int,double> >      derivz_;          // direction derivative of nodal z-matrix value

}; // class CNodeDataContainer
} // namespace CONTACT



#endif  // #ifndef CNODEDATACONTAINER_H
#endif  // #ifdef CCADISCRET
