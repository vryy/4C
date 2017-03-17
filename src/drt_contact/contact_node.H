/*----------------------------------------------------------------------*/
/*!
\file contact_node.H

\brief A class for a contact node

\level 2

\maintainer Philipp Farah, Alexander Seitz

*/
/*----------------------------------------------------------------------*/
#ifndef CONTACT_NODE_H
#define CONTACT_NODE_H

#include "../drt_mortar/mortar_node.H"

/*!
 \brief CONTACT: namespace of the BACI contact and meshtying algorithms

 */
namespace CONTACT
{

// foward declaration
class CoNode;

class CoNodeType: public DRT::ParObjectType
{
public:

  std::string Name() const { return "CoNodeType"; }

  static CoNodeType & Instance() { return instance_; }
  ;

  virtual DRT::ParObject* Create(const std::vector<char> & data);

private:

  static CoNodeType instance_;
};

/*!
 \brief A class containing additional data from contact nodes

 This class contains additional information from contact nodes which are
 are not needed for contact search and therefore are only available on the
 node's processor (ColMap). The class CoNodeDataContainer must be declared
 before the CoNode itself.

 \author gitterle (gitterle@lnm.mw.tum.de)
 */
class CoNodeDataContainer
{
public:

  //! @name Constructors and destructors and related methods

  /*!
   \brief Standard Constructor

   */
  CoNodeDataContainer();

  /*!
   \brief Destructor

   */
  virtual ~CoNodeDataContainer() { return; }

  /*!
   \brief Pack this class so that it can be communicated

   This function packs the datacontainer. This is only called
   when the class has been initialized and the pointer to this
   class exists.

   */
  virtual void Pack(DRT::PackBuffer& data) const;

  /*!
   \brief Unpack data from a vector into this class

   This function unpacks the datacontainer. This is only called
   when the class has been initialized and the pointer to this
   class exists.

   */
  virtual void Unpack(std::vector<char>::size_type& position,
      const std::vector<char>& data);

  //@}

  //! @name Access methods

  /*!
   \brief Return current nodal tangent t_xi (only for slave side!) (length 3)
   */
  virtual double* txi() { return txi_; }

  /*!
   \brief Return current nodal tangent t_eta (only for slave side!) (length 3)
   */
  virtual double* teta() { return teta_; }

  /*!
   \brief Return the weighted gap (scalar) of this node
   */
  virtual double& Getg() { return grow_; }

  /*!
   \brief Return the weighted gap (scalar) of this node for NTS
   */
  virtual double& Getgnts() { return gnts_; }

  /*!
   \brief Return the weighted gap (scalar) of this node for LTS
   */
  virtual double& Getglts() { return glts_; }

  /*!
  \brief Return the weighted gap (scalar) of this node for LTL
  */
 virtual double& Getgltl() { return gltl_; }

  //! @name Extended contact
  //! @{

  //! Return the constraint residuum of this node
  virtual double& GetWcLm() {return Wc_lm_;}

  //! Return linearization of this node (slave part)
  virtual GEN::pairedvector<int, double>& GetWcSuLm() {return Wc_su_lm_;}

  //! Return linearization of this node (master part)
  virtual std::map<int, double>& GetWcMuLm() {return Wc_mu_lm_;}

  //! Return linearization of this node (slave part)
  virtual std::map<int, std::map<int, double> >& GetWcSuU() {return Wc_su_u_;}

  //! Return linearization of this node (master part)
  virtual std::map<int, std::map<int, double> >& GetWcMuU() {return Wc_mu_u_;}


  //! @}

  /*!
  \brief Return contact status of last converged state n (active=true)

  */
  virtual inline bool& ActiveOld() { return activeold_; }

  /*!
   \brief Return the 'DerivN' map (vector) of this node

   These maps contain the directional derivatives of the node's
   averaged normal vector with respect to the slave displacements.
   A vector is used because the normal itself is a vector (2 or 3 components).

   */
  virtual std::vector<GEN::pairedvector<int, double> >& GetDerivN() { return derivn_; }

  /*!
   \brief Return the 'DerivTxi' map (vector) of this node

   These maps contain the directional derivatives of the node's
   unit tangent vector t_xi with respect to the slave displacements.
   A vector is used because the tangent itself is a vector (2 or 3 components).

   */
  virtual std::vector<GEN::pairedvector<int, double> >& GetDerivTxi() { return derivtxi_; }

  /*!
   \brief Return the 'DerivTeta' map (vector) of this node

   These maps contain the directional derivatives of the node's
   unit tangent vector t_eta with respect to the slave displacements.
   A vector is used because the tangent itself is a vector (2 or 3 components).

   */
  virtual std::vector<GEN::pairedvector<int, double> >& GetDerivTeta() { return derivteta_; }

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
  virtual std::map<int, std::map<int, double> >& GetDerivD() { return derivd_; }
  virtual std::map<int, std::map<int, double> >& GetDerivDlts() { return derivdlts_; }
  virtual std::map<int, std::map<int, double> >& GetDerivDltl() { return derivdltl_; }

  /*!
   \brief Return the 'DerivM' map of maps of this node

   This map contains the directional derivatives of the node's
   M-matrix entries with respect to the slave/master displacements.
   It is a "map of maps", with the outer map containing all master
   node indices k adjacent to this node i and the inner map containing
   all directional derivatives l existing for M_ik.

   */
  virtual std::map<int, std::map<int, double> >& GetDerivM() { return derivm_; }
  virtual std::map<int, std::map<int, double> >& GetDerivMnts() { return derivmnts_; }
  virtual std::map<int, std::map<int, double> >& GetDerivMlts() { return derivmlts_; }
  virtual std::map<int, std::map<int, double> >& GetDerivMltl() { return derivmltl_; }

  /*!
   \brief Return one specific 'DerivD' map of this node

   This method returns the map of directional derivatives of one
   specific D-matrix D_ik entry of this node i.

   */
  virtual std::map<int, double>& GetDerivD(int& k)
  {
    typedef std::map<int, std::map<int, double> >::const_iterator CI;
    CI p = derivd_.find(k);
    if (p == derivd_.end())
      dserror("ERROR: GetDerivD: No map entry existing for given index");
    return derivd_[k];
  }

  /*!
   \brief Return one specific 'DerivM' map of this node

   This method returns the map of directional derivatives of one
   specific M-matrix M_ik entry of this node i.

   */
  virtual std::map<int, double>& GetDerivM(int& k)
  {
    typedef std::map<int, std::map<int, double> >::const_iterator CI;
    CI p = derivm_.find(k);
    if (p == derivm_.end())
      dserror("ERROR: GetDerivM: No map entry existing for given index");
    return derivm_[k];
  }

  /*!
   \brief Return the 'DerivG' map of this node

   This map contains the directional derivatives of the node's
   weighted gap entry g~ with respect to the slave/master displacements.

   */
  virtual std::map<int, double>& GetDerivG() { return derivg_; }
  virtual std::map<int, double>& GetDerivGnts() { return derivgnts_;}
  virtual std::map<int, double>& GetDerivGlts() { return derivglts_;}
  virtual std::map<int, double>& GetDerivGltl() { return derivgltl_;}

  /*!
   \brief Return the 'DerivGlm' map of this node

   */
  virtual std::map<int, double>& GetDerivGW() { return derivgw_; }

  /*!
   \brief Return the 'DerivW' map of this node

   This map contains the directional derivatives of the node's
   weighted wear increment entry w~ with respect to the slave/master displacements.

   */
  virtual std::map<int, double>& GetDerivW() { return derivw_; }

  /*!
   \brief Return the 'DerivW' map of this node

   This map contains the directional derivatives of the node's
   weighted wear increment entry w~ with respect to the lagr. mult.

   */
  virtual std::map<int, double>& GetDerivWlm() { return derivw_lm_; }

  /*!
   \brief Return scaling factor for weighted gap

   Note: This is only calculated when performing a penalty strategy

   */
  virtual double& Kappa() { return kappa_; }

  /*!
   \brief Return the 'DerivZ map of this node

   This map contains the directional derivatives of the node's
   lagrange multiplier entries with respect to the slave/master displacements.

   Note: This is only calculated when performing a penalty strategy

   */
  virtual std::vector<std::map<int, double> >& GetDerivZ() { return derivz_; }

  /*!
   \brief Return the 'Derivscalefac_' map of this node

   This map contains the directional derivatives of the node's
   scaling factor witch is needed for better conditioning

   */
  virtual std::map<int, double>& GetDerivScale() { return derivscalefac_; }

  virtual GEN::pairedvector<int, double>& GetAlpha() {return alpha_;};
  virtual double& GetAlphaN() {return nalpha_;};

  //@}

protected:

  // don't want = operator and cctor
  CoNodeDataContainer operator =(const CoNodeDataContainer& old);
  CoNodeDataContainer(const CoNodeDataContainer& old);

  double txi_[3];  // nodal tangent t_xi for contact methods at newton step n
  double teta_[3]; // nodal tangent t_eta for contact methods at newton step n
  double grow_;    // nodal entry of g vector
  double gnts_;    // nodal entry of g vector for nts
  double glts_;    // nodal entry of g vector for lts
  double gltl_;    // nodal entry of g vector for lts

  double nalpha_;  // scaling factor for nts-mortar hybrid formulation
  bool activeold_; // true if cnode was in contact (last converged state n)

  std::vector<GEN::pairedvector<int, double> > derivn_;    // directional derivative of nodal normal
  std::vector<GEN::pairedvector<int, double> > derivtxi_;  // directional derivative of nodal tangent t_xi
  std::vector<GEN::pairedvector<int, double> > derivteta_; // directional derivative of nodal tangent t_eta
  std::map<int, std::map<int, double> > derivd_;           // directional derivative of nodal D-matrix value
  std::map<int, std::map<int, double> > derivdlts_;           // directional derivative of nodal D-matrix value lts
  std::map<int, std::map<int, double> > derivdltl_;           // directional derivative of nodal D-matrix value ltl
  std::map<int, std::map<int, double> > derivm_;           // directional derivative of nodal M-matrix values
  std::map<int, std::map<int, double> > derivmnts_;        // directional derivative of nodal M-matrix values nts
  std::map<int, std::map<int, double> > derivmlts_;        // directional derivative of nodal M-matrix values lts
  std::map<int, std::map<int, double> > derivmltl_;        // directional derivative of nodal M-matrix values ltl
  std::map<int, double> derivg_;                           // directional derivative of nodal weighted gap value
  std::map<int, double> derivgnts_;                        // directional derivative of nts gap value
  std::map<int, double> derivglts_;                        // directional derivative of lts gap value
  std::map<int, double> derivgltl_;                        // directional derivative of ltl gap value
  std::map<int, double> derivgw_;                          // directional derivative of nodal weighted gap value
  std::map<int, double> derivw_;                           // directional derivative of nodal weighted wear value
  std::map<int, double> derivw_lm_;                        // directional derivative of nodal weighted wear value
  std::map<int, double> derivscalefac_;                    // directional derivative of nodal scale factor value

  GEN::pairedvector<int, double> alpha_;    // lin. of scaling factor for hybrid formulation

  // *** Penalty-related quantities **************************
  double kappa_;                               // gap-scaling factor kappa
  std::vector<std::map<int, double> > derivz_; // direction derivative of nodal z-matrix value


  //! @name Extended contact formulation
  //! @{

  // TODO: Create own data container (or avoid storing in nodes) for XContact formulation

  //! Contact potential varied after Lagrange multiplier (constraint residuum vector)
  double Wc_lm_;

  //! Contact potential varied after slave displacements and linearized after Lagrange multipliers
  GEN::pairedvector<int, double> Wc_su_lm_;

  //! Contact potential varied after master displacements and linearized after Lagrange multipliers
  std::map<int, double> Wc_mu_lm_;

  //! Contact potential varied after slave displacements and linearized after all displacements
  std::map<int, std::map<int, double> > Wc_su_u_;

  //! Contact potential varied after master displacements and linearized after all displacements
  std::map<int, std::map<int, double> > Wc_mu_u_;

  //! @}
};

namespace AUG {
  class NodeDataContainer
  {
  public:
    /// @name General methods
    /// @{

    /// constructor
    NodeDataContainer( CoNode& parentNode );

    /// constructor
    NodeDataContainer( CoNode& parentNode, int maxNumMasterElements );

    virtual ~NodeDataContainer() {};

    void Setup();

    /*! \brief Pack this class so that it can be communicated
     *
     *  This function packs the datacontainer. This is only called
     *  when the class has been initialized and the pointer to this
     *  class exists. */
    void Pack( DRT::PackBuffer& data ) const;

    /*! \brief Unpack data from a vector into this class
     *
     *  This function unpacks the data container. This is only called
     *  when the class has been initialized and the pointer to this
     *  class exists. */
    void Unpack( std::vector<char>::size_type& position,
                 const std::vector<char>& data );

    /// @}

    /// @name accessors
    /// @{

    /*! \brief Return the weighted gap (scalar) of this node */
    inline double& GetWGap() { return wGap_; }
    inline double GetWGap() const { return wGap_; }

    /*! \brief Return the scaling factor kappa for this node */
    inline double& GetKappa() { return kappa_; }
    inline double GetKappa() const { return kappa_; }

    /*!
    \brief Return the linearization of the scaling factor kappa

    */
    std::map<int,double>& GetKappaLin() { return kappaLin_; }
    const std::map<int,double>& GetKappaLin() const { return kappaLin_; }

    /*! \brief Return the scaling factor augA
     *
     *  Integration over the whole slave interface, without consideration of
     *  the segments */
    inline double& GetAugA() { return augA_; }
    inline double GetAugA() const { return augA_; }

    /*! \brief Return the linearization of the scaling factor augA
     *
     *  Integration over the whole slave interface, without consideration of
     *  the segments */
    inline GEN::pairedvector<int,double>& GetAugALin() { return augALin_; }
    inline const GEN::pairedvector<int,double>& GetAugALin() const { return augALin_; }

    /*! \brief Return the linearization of the averaged weighted gap */
    std::map<int,double>& GetAWGapLin() {return aWGapLin_; }
    const std::map<int,double>& GetAWGapLin() const {return aWGapLin_; }

    /*! \brief Return the linearization of the weighted gap */
    std::map<int,double>& GetWGapLin() {return wGapLin_; }
    const std::map<int,double>& GetWGapLin() const {return wGapLin_; }

    /*! \brief Return the variation of the weighted gap (slave part) */
    inline GEN::pairedvector<int,std::pair<int,double> >& GetVarWGapSl() {return varWGapSl_; }
    inline const GEN::pairedvector<int,std::pair<int,double> >& GetVarWGapSl() const {return varWGapSl_; }

    /*! \brief Return the variation of the weighted gap (master part) */
    std::map<int,std::pair<int,double> >& GetVarWGapMa() {return varWGapMa_; }
    const std::map<int,std::pair<int,double> >& GetVarWGapMa() const {return varWGapMa_; }

    /*! \brief Return the linearization of the variation of the weighted gap
     *  (slave part) */
    std::map<int,std::map<int,double> >& GetVarWGapLinSl() {return varWGapLinSl_; }
    const std::map<int,std::map<int,double> >& GetVarWGapLinSl() const {return varWGapLinSl_; }

    /*! \brief Return the linearization of the variation of the weighted gap
     *  (master part) */
    std::map<int,std::map<int,double> >& GetVarWGapLinMa() {return varWGapLinMa_; }
    const std::map<int,std::map<int,double> >& GetVarWGapLinMa() const {return varWGapLinMa_; }

    /// @}

  private:
    // *** Augmented Lagrangian formulation ********************

    /// one slave element projects approximately maximal in this number of master elements
    int maxNumMasterElements_;

    // nodal scalar values
    double kappa_;        ///< gap-scaling factor kappa
    double wGap_;         ///< nodal entry of weighted gap vector
    double augA_;         ///< nodal scaling factor

    // Variation -- RHS
    GEN::pairedvector<int,std::pair<int,double> > varWGapSl_; ///< variation of the weighted gap (slave)
    std::map<int,std::pair<int,double> > varWGapMa_;          ///< variation of the weighted gap (master)

    // Linearization
    std::map<int,std::map<int,double> > varWGapLinSl_;  ///< linearization of the weighted gap variation (slave)
    std::map<int,std::map<int,double> > varWGapLinMa_;  ///< linearization of the weighted gap variation (master)
    std::map<int,double> kappaLin_;     ///< linearization of the scaling factor (segment based integration)
    GEN::pairedvector<int,double> augALin_;             ///< linearization of the scaling factor (element based integration)
    std::map<int,double> aWGapLin_;     ///< linearization of the averaged weighted gap
    std::map<int,double> wGapLin_;      ///< linearization of the weighted gap

    const CoNode& parentNode_;
  };  // class NodeDataContainer
} // namespace AUG

// class CoNodeDataContainer

/*!
\brief A class containing additional data from poro contact nodes

This class contains additional information from poro contact nodes which are
are not needed for normal contact nodes.
CoNodePoroDataContainer must be declared before the CoNode itself.

\author ager (ager@lnm.mw.tum.de)
*/
class CoNodePoroDataContainer
{
public:

  //! @name Constructors and destructors and related methods

  /*!
  \brief Standard Constructor

  */
  CoNodePoroDataContainer();

  /*!
  \brief Destructor

  */
  virtual ~CoNodePoroDataContainer() { return; }

  /*!
  \brief Pack this class so that it can be communicated

  This function packs the datacontainer. This is only called
  when the class has been initialized and the pointer to this
  class exists.

  */
  virtual void Pack(DRT::PackBuffer& data) const;

  /*!
  \brief Unpack data from a vector into this class

  This function unpacks the datacontainer. This is only called
  when the class has been initialized and the pointer to this
  class exists.

  */
  virtual void Unpack(std::vector<char>::size_type& position, const std::vector<char>& data);

  //! @name Access methods

  /*!
  \brief Return the normal coupling condition (scalar) of this node
  */
  virtual double& GetnCoup() { return ncouprow_; }

  /*!
  \brief Return the 'DerivnCoup' map of this node

  This map contains the directional derivatives of the node's
  normal coupling condition with respect to the slave/master displacements.

  */
  virtual std::map<int,double>& GetDerivnCoup() { return derivncoup_; }

  /*!
  \brief Return the 'VelDerivnCoup' map of this node

  This map contains the derivatives of the node's
  normal coupling condition with respect to the velocities. (for one sided contact just slave!)

  */
  virtual std::map<int,double>& GetVelDerivnCoup() { return velderivncoup_; }

  /*!
  \brief Return the 'PresDerivnCoup' map of this node

  This map contains the derivatives of the node's
  normal coupling condition with respect to the pressures. (for one sided contact just slave!)
  //h.Willmann
  */
  virtual std::map<int,double>& GetPresDerivnCoup() { return presderivncoup_; }

  /*!
  \brief Return current nodal fluid pressure (porous media!) (length 3)
  */
  virtual double* fpres() { return &fpres_; }

  /*!
  \brief Return current nodal fluid velocity (porous media!) (length 3)
  */
  virtual double* fvel() { return fvel_; }

  /*!
  \brief Return current nodal structural velocity (porous media!) (length 3)
  */
  virtual double* svel() { return svel_; }

  /*!
  \brief Return current nodal lagrangean multiplier (length 3)
  */
  virtual double* poroLM() { return porolm_; }

protected:

  // don't want = operator and cctor
  CoNodePoroDataContainer operator = (const CoNodePoroDataContainer& old);
  CoNodePoroDataContainer(const CoNodePoroDataContainer& old);


  double          ncouprow_;           // nodal entry of n-coupling vector
  std::map<int,double>                    derivncoup_;    // directional derivative of nodal weighted n-coupling vector
  std::map<int,double>                    velderivncoup_; // velocity derivative of nodal weighted n-coupling vector
  std::map<int,double>                    presderivncoup_; // pressure derivative of nodal weighted n-coupling vector //h.Willmann

  double           fvel_[3];        //fluid velocity for porous problem
  double           fpres_;       //fluid pressure for porous problem
  double           svel_[3];        //structural velocity for porous problem
  double           porolm_[3];      //lagrange multiplier from poro no penetration condition!
}; //class CoNodePoroDataContainer

/*!
\brief A class containing additional data from tsi contact nodes

This class contains additional information from tsi contact nodes which are
are not needed for normal contact nodes.

\author seitz (seitz@lnm.mw.tum.de)
*/
class CoNodeTSIDataContainer
{
public:

  //! @name Constructors and destructors and related methods

  /*!
  \brief Standard Constructor

  */
  CoNodeTSIDataContainer(double t_ref, double t_dam);

  /*!
  \brief empty Constructor: unpack data later

  */
  CoNodeTSIDataContainer(){};

  /*!
  \brief Destructor

  */
  virtual ~CoNodeTSIDataContainer() { return; }

  /*!
  \brief Pack this class so that it can be communicated

  This function packs the datacontainer. This is only called
  when the class has been initialized and the pointer to this
  class exists.

  */
  virtual void Pack(DRT::PackBuffer& data) const;

  /*!
  \brief Unpack data from a vector into this class

  This function unpacks the datacontainer. This is only called
  when the class has been initialized and the pointer to this
  class exists.

  */
  virtual void Unpack(std::vector<char>::size_type& position, const std::vector<char>& data);


  //! @name Access methods

  /*!
  \brief Return max (Temp_slave , Temp_master)
  */
  double& TempMaster() { return temp_master_; }

  /*!
  \brief Return temperature
  */
  double& Temp() { return temp_; }

  /*!
  \brief Return reference temperature
  */
  double& Temp_Ref() { return t_ref_; }

  /*!
  \brief Return temperature
  */
  double& Temp_Dam() { return t_dam_; }

  /*!
  \brief Return thermo Lagrange multiplier
  */
  double& ThermoLM() { return thermo_lm_; }

  std::map<int,double>& DerivTempMasterDisp() {return derivTempMasterDisp_;}
  std::map<int,double>& DerivTempMasterTemp() {return derivTempMasterTemp_;}

  void Clear();


protected:

  // don't want = operator and cctor
  CoNodeTSIDataContainer operator = (const CoNodeTSIDataContainer& old);
  CoNodeTSIDataContainer(const CoNodeTSIDataContainer& old);

  double temp_;
  double t_ref_;
  double t_dam_;
  double thermo_lm_;

  double temp_master_;
  std::map<int,double>derivTempMasterDisp_;
  std::map<int,double>derivTempMasterTemp_;

}; //class CoNodeTSIDataContainer

/*!
 \brief A class for a contact node derived from MORTAR::MortarNode

 This class represents a finite element node capable of contact.

 \author popp (popp@lnm.mw.tum.de)
 */
class CoNode: public MORTAR::MortarNode
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
  CoNode(int id, const double* coords, const int owner, const int numdof,
      const std::vector<int>& dofs, const bool isslave, const bool initactive);

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
  virtual ~CoNode()
  {
    return;
  }

  /*!
   \brief Return unique ParObject id

   every class implementing ParObject needs a unique id defined at the
   top of drt_lib/drt_parobject.H.

   */
  virtual int UniqueParObjectId() const
  {
    return CoNodeType::Instance().UniqueParObjectId();
  }

  /*!
   \brief Pack this class so it can be communicated

   \ref Pack and \ref Unpack are used to communicate this node

   */
  virtual void Pack(DRT::PackBuffer& data) const;

  /*!
   \brief Unpack data from a char vector into this class

   \ref Pack and \ref Unpack are used to communicate this node

   */
  virtual void Unpack(const std::vector<char>& data);

  //@}

  //! @name Access methods

  /*!
   \brief Print this cnode
   */
  virtual void Print(std::ostream& os) const;

  /*!
   \brief Is Node initialized as active node (only slave nodes)
   */
  virtual bool IsInitActive() const
  {
    if (!IsSlave())
      dserror("ERROR: InitActive requested for Master node");
    return initactive_;
  }

  /*!
   \brief Modify initial active status of slave node

   This belated modification is necessary to be able to use
   the binary search tree for contact initialization in the
   load-controlled quasistatic case (instead of input file
   information Active/Inactive)

   */
  virtual bool& SetInitActive()
  {
    if (!IsSlave())
      dserror("ERROR: InitActive requested for Master node");
    return initactive_;
  }

  /*!
   \brief Return contact status of this node (active=true)
   */
  virtual bool& Active()
  {
    return active_;
  }

  virtual bool& InvolvedM() { return involvedm_; }

  /*!
   \brief Return data container of this node

   This method returns the data container of this node where additional
   contact specific quantities/information are stored.

   */
  inline CONTACT::CoNodeDataContainer& CoData() { return *codata_; }

  /*! \brief Return data container of this augmented contact node
   *
   *  This method returns the data container of this node where additional
   *  augmented contact specific quantities/information are stored.
   *
   *  \author hiermeier \date 03/17 */
  inline CONTACT::AUG::NodeDataContainer& AugData()
  {
    if ( augdata_.is_null() )
      dserror("There are no augmented contact node data available!" );
    return *augdata_;
  };
  inline const CONTACT::AUG::NodeDataContainer& AugData() const
  {
    if ( augdata_.is_null() )
      dserror("There are no augmented contact node data available!" );
    return *augdata_;
  };

  inline CONTACT::CoNodePoroDataContainer& CoPoroData() { return *coporodata_; }
  inline CONTACT::CoNodeTSIDataContainer& CoTSIData() { return *cTSIdata_; }

  //@}

  //! @name Evaluation methods

  /*!
   \brief Add a value to the weighted gap of this node

   This value is later assembled to the weighted gap vec.
   Note that grow_ here is a scalar.

   \param val : value to be added

   */
  void AddgValue(double& val);

  /*!
   \brief Add a value to the point-wise gap of this node (NTS)

   \param val : value to be added

   */
  void AddntsGapValue(double& val);

  /*!
   \brief Add a value to the line-weighted gap of this node (LTS)

   \param val : value to be added

   */
  void AddltsGapValue(double& val);

  /*!
   \brief Add a value to the point-wise gap of this node (LTL)

   \param val : value to be added

   */
  void AddltlGapValue(double& val);


  /*!
  \brief Add a value to the weighted gap of this node (augmented Lagrange)

  This value is later assembled to the averaged weighted gap vec.
  Note that grow_ here is a scalar.

  \param val : value to be added

  */
  void AddWGapValue(double& val);

  /*!
  \brief Add a value to the scaling factor kappa of this node (augmented Lagrange)

  This value is later assembled needed during the assembling of the averaged
  weighted gap. Note that kappa_ is a scalar.

  \param val : value to be added
  */
  void AddKappaValue(double& val);

  /*!
  \brief Add a value to the variation of weighted gap (slave)

  \param cal : slave dof id
  \param gid : slave gid
  \param val : value to be added
  */
  void AddVarWGapSl(int col, int gid, double val);

  /*!
  \brief Add a value to the variation of weighted gap (master)

  \param cal : master dof id
  \param gid : master gid
  \param val : value to be added
  */
  void AddVarWGapMa(int col, int gid, double val);

  /*!
  \brief TODO
  */
  void AddWcLm( double val );

  /*!
  \brief TODO
  */
  void AddWcSuLm( int col, double val );

  /*!
  \brief TODO
  */
  void AddWcMuLm(int col, double val);

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
  \brief Add a value to the NCoup of this node

  \param val : value to be added

  */
  void AddNcoupValue(double& val);

  /*!
   \brief Build nodal normal
   */
  void BuildAveragedNormal();

  /*!
   \brief Initializes the data container of the node

   With this function, the container with contact specific quantities/information
   is initialized.

   */
  virtual void InitializeDataContainer();

  /*!
   \brief Initializes the data container of the augmented contact node

   With this function, the container with augmented contact specific
   quantities/information is initialized.

   */
  void InitializeAugDataContainer( int maxNumMasterEles );

  /*!
  \brief Initializes the poro data container of the node

  With this function, the container with contact specific quantities/information
  is initialized. --- Used to initialize PoroDataContainer for master nodes!

  */
  virtual void InitializePoroDataContainer();

  /*!
  \brief Initializes the TSI data container of the node
  */
  virtual void InitializeTSIDataContainer(double t_ref, double t_dam);

  /*!
   \brief Resets the data container of the node

   With this function, the container with contact specific quantities/information
   is deleted / reset to Teuchos::null pointer

   */
  virtual void ResetDataContainer();

  /*!
   \brief Get number of linearization entries

   */
  virtual int& GetLinsize() { return linsize_; };
  int GetLinsize() const { return linsize_; };
  //@}

  //! @name Empty functions (friction only)

  // All these functions only have functionality for friction nodes, thus they are
  // defined as empty here in the general mortar node. They can be called whenever you like.

  virtual void AddSNode(int node) {}
  virtual void AddMNode(int node) {}

  virtual bool HasCoTSIData(){return (cTSIdata_!=Teuchos::null);}
  //@}

private:
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
  void DerivAveragedNormal(Epetra_SerialDenseMatrix& elens, double length,
      double ltxi);

  void DerivAveragedNormalHermit(double length, double ltxi);

protected:

  bool active_;         // true if cnode is in contact (active set strategy)
  bool initactive_;     // true if node is initialized as active node
  bool involvedm_;      // cnode is an master node for which the mortar integration is performed (both-sided wear)
  int  linsize_;     // number of lin entries per direction vector

  Teuchos::RCP<CONTACT::CoNodeDataContainer> codata_; // additional information of proc's contact nodes

  Teuchos::RCP<CONTACT::AUG::NodeDataContainer> augdata_; // additional information of proc's augmented contact nodes

  Teuchos::RCP<CONTACT::CoNodePoroDataContainer> coporodata_;   // additional information of proc's poro contact nodes

  Teuchos::RCP<CONTACT::CoNodeTSIDataContainer> cTSIdata_;   // additional information of TSI contact nodes

};
// class CoNode
}// namespace CONTACT

// << operator
std::ostream& operator <<(std::ostream& os, const CONTACT::CoNode& cnode);

#endif  // #ifndef CNODE_H
