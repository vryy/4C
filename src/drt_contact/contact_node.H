/*----------------------------------------------------------------------*/
/*! \file
\brief A class for a contact node

\level 2


*/
/*----------------------------------------------------------------------*/
#ifndef CONTACT_NODE_H
#define CONTACT_NODE_H

#include "../drt_mortar/mortar_node.H"
#include "../drt_contact_aug/contact_aug_utils.H"
#include <unordered_map>

namespace CONTACT
{
  // forward declaration
  class CoNode;

  class CoNodeType : public DRT::ParObjectType
  {
   public:
    std::string Name() const final { return "CoNodeType"; }

    static CoNodeType& Instance() { return instance_; };

    virtual DRT::ParObject* Create(const std::vector<char>& data);

   private:
    static CoNodeType instance_;
  };

  /*!
   \brief A class containing additional data from contact nodes

   This class contains additional information from contact nodes which are
   are not needed for contact search and therefore are only available on the
   node's processor (ColMap). The class CoNodeDataContainer must be declared
   before the CoNode itself.

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

     This function packs the data container. This is only called
     when the class has been initialized and the pointer to this
     class exists.

     */
    virtual void Pack(DRT::PackBuffer& data) const;

    /*!
     \brief Unpack data from a vector into this class

     This function unpacks the data container. This is only called
     when the class has been initialized and the pointer to this
     class exists.

     */
    virtual void Unpack(std::vector<char>::size_type& position, const std::vector<char>& data);

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
    \brief Return the gap (vector) of this node for LTL
    */
    virtual double* Getgltl() { return gltl_; }

    /*!
    \brief Return the gap (vector) of this node for LTL
    */
    virtual double* Getjumpltl() { return jumpltl_; }


    //! @name Extended contact
    //! @{

    //! Return the constraint residuum of this node
    virtual double& GetWcLm() { return Wc_lm_; }

    //! Return linearization of this node (slave part)
    virtual GEN::pairedvector<int, double>& GetWcSuLm() { return Wc_su_lm_; }

    //! Return linearization of this node (master part)
    virtual std::map<int, double>& GetWcMuLm() { return Wc_mu_lm_; }

    //! Return linearization of this node (slave part)
    virtual std::map<int, std::map<int, double>>& GetWcSuU() { return Wc_su_u_; }

    //! Return linearization of this node (master part)
    virtual std::map<int, std::map<int, double>>& GetWcMuU() { return Wc_mu_u_; }


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
    virtual std::vector<GEN::pairedvector<int, double>>& GetDerivN() { return derivn_; }
    virtual std::vector<GEN::pairedvector<int, double>>& GetDerivTangent() { return derivEdge_; }

    /*!
     \brief Return the 'DerivTxi' map (vector) of this node

     These maps contain the directional derivatives of the node's
     unit tangent vector t_xi with respect to the slave displacements.
     A vector is used because the tangent itself is a vector (2 or 3 components).

     */
    virtual std::vector<GEN::pairedvector<int, double>>& GetDerivTxi() { return derivtxi_; }

    /*!
     \brief Return the 'DerivTeta' map (vector) of this node

     These maps contain the directional derivatives of the node's
     unit tangent vector t_eta with respect to the slave displacements.
     A vector is used because the tangent itself is a vector (2 or 3 components).

     */
    virtual std::vector<GEN::pairedvector<int, double>>& GetDerivTeta() { return derivteta_; }

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
    virtual std::map<int, std::map<int, double>>& GetDerivD() { return derivd_; }
    virtual std::map<int, std::map<int, double>>& GetDerivDlts() { return derivdlts_; }
    virtual std::map<int, std::map<int, double>>& GetDerivDltl() { return derivdltl_; }

    /*!
     \brief Return the 'DerivM' map of maps of this node

     This map contains the directional derivatives of the node's
     M-matrix entries with respect to the slave/master displacements.
     It is a "map of maps", with the outer map containing all master
     node indices k adjacent to this node i and the inner map containing
     all directional derivatives l existing for M_ik.

     */
    virtual std::map<int, std::map<int, double>>& GetDerivM() { return derivm_; }
    virtual std::map<int, std::map<int, double>>& GetDerivMnts() { return derivmnts_; }
    virtual std::map<int, std::map<int, double>>& GetDerivMlts() { return derivmlts_; }
    virtual std::map<int, std::map<int, double>>& GetDerivMltl() { return derivmltl_; }

    /*!
     \brief Return one specific 'DerivD' map of this node

     This method returns the map of directional derivatives of one
     specific D-matrix D_ik entry of this node i.

     */
    virtual std::map<int, double>& GetDerivD(int& k)
    {
      typedef std::map<int, std::map<int, double>>::const_iterator CI;
      CI p = derivd_.find(k);
      if (p == derivd_.end()) dserror("ERROR: GetDerivD: No map entry existing for given index");
      return derivd_[k];
    }

    /*!
     \brief Return one specific 'DerivM' map of this node

     This method returns the map of directional derivatives of one
     specific M-matrix M_ik entry of this node i.

     */
    virtual std::map<int, double>& GetDerivM(int& k)
    {
      typedef std::map<int, std::map<int, double>>::const_iterator CI;
      CI p = derivm_.find(k);
      if (p == derivm_.end()) dserror("ERROR: GetDerivM: No map entry existing for given index");
      return derivm_[k];
    }

    /*!
     \brief Return the 'DerivG' map of this node

     This map contains the directional derivatives of the node's
     weighted gap entry g~ with respect to the slave/master displacements.

     */
    virtual std::map<int, double>& GetDerivG() { return derivg_; }
    virtual std::map<int, double>& GetDerivGnts() { return derivgnts_; }
    virtual std::map<int, double>& GetDerivGlts() { return derivglts_; }
    virtual std::vector<std::map<int, double>>& GetDerivGltl() { return derivgltl_; }
    virtual std::vector<std::map<int, double>>& GetDerivJumpltl() { return derivjumpltl_; }

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
    virtual std::vector<std::map<int, double>>& GetDerivZ() { return derivz_; }

    virtual GEN::pairedvector<int, double>& GetAlpha() { return alpha_; };
    virtual double& GetAlphaN() { return nalpha_; };

    /*!
     \brief Return old normal vector (length 3)
     */
    virtual inline double* Normal_old() { return n_old_; }

    //! @}

   protected:
    //! nodal tangent t_xi for contact methods at newton step n
    double txi_[3];

    //! nodal tangent t_eta for contact methods at newton step n
    double teta_[3];

    //! nodal entry of g vector
    double grow_;

    //! nodal entry of g vector for node-to-segment
    double gnts_;

    //! nodal entry of g vector for line-to-segment
    double glts_;

    //! nodal entry of g vector for line-to-segment
    double gltl_[3];

    //! nodal entry of jump vector for line-to-segment
    double jumpltl_[3];

    //! scaling factor for node-to-segment-mortar hybrid formulation
    double nalpha_;

    //! true if contact node was in contact (last converged state n)
    bool activeold_;

    //! directional derivative of nodal normal
    std::vector<GEN::pairedvector<int, double>> derivn_;

    //! directional derivative of nodal tangent
    std::vector<GEN::pairedvector<int, double>> derivEdge_;

    //! directional derivative of nodal tangent t_xi
    std::vector<GEN::pairedvector<int, double>> derivtxi_;

    //! directional derivative of nodal tangent t_eta
    std::vector<GEN::pairedvector<int, double>> derivteta_;

    //! directional derivative of nodal D-matrix value
    std::map<int, std::map<int, double>> derivd_;

    //! directional derivative of nodal D-matrix value line-to-segment
    std::map<int, std::map<int, double>> derivdlts_;

    //! directional derivative of nodal D-matrix value line-to-line
    std::map<int, std::map<int, double>> derivdltl_;

    //! directional derivative of nodal M-matrix values
    std::map<int, std::map<int, double>> derivm_;

    //! directional derivative of nodal M-matrix values node-to-segment
    std::map<int, std::map<int, double>> derivmnts_;

    //! directional derivative of nodal M-matrix values line-to-segment
    std::map<int, std::map<int, double>> derivmlts_;

    //! directional derivative of nodal M-matrix values line-to-line
    std::map<int, std::map<int, double>> derivmltl_;

    //! directional derivative of nodal weighted gap value
    std::map<int, double> derivg_;

    //! directional derivative of node-to-segment gap value
    std::map<int, double> derivgnts_;

    //! directional derivative of line-to-segment gap value
    std::map<int, double> derivglts_;

    //! directional derivative of line-to-line gap value
    std::vector<std::map<int, double>> derivgltl_;

    //! directional derivative of line-to-line gap value
    std::vector<std::map<int, double>> derivjumpltl_;

    //! directional derivative of nodal weighted gap value
    std::map<int, double> derivgw_;

    //! directional derivative of nodal weighted wear value
    std::map<int, double> derivw_;

    //! directional derivative of nodal weighted wear value
    std::map<int, double> derivw_lm_;

    //! lin. of scaling factor for hybrid formulation
    GEN::pairedvector<int, double> alpha_;

    //! @name Penalty-related quantities
    //! @{

    //! gap-scaling factor kappa
    double kappa_;

    //! @}

    //! direction derivative of nodal z-matrix value
    std::vector<std::map<int, double>> derivz_;


    //! @name Extended contact formulation
    //! @{

    // TODO: Create own data container (or avoid storing in nodes) for XContact formulation

    //! Contact potential varied after Lagrange multiplier (constraint residuum vector)
    double Wc_lm_;

    //! Contact potential varied after slave displacements and linearized after Lagrange multipliers
    GEN::pairedvector<int, double> Wc_su_lm_;

    //! Contact potential varied after master displacements and linearized after Lagrange
    //! multipliers
    std::map<int, double> Wc_mu_lm_;

    //! Contact potential varied after slave displacements and linearized after all displacements
    std::map<int, std::map<int, double>> Wc_su_u_;

    //! Contact potential varied after master displacements and linearized after all displacements
    std::map<int, std::map<int, double>> Wc_mu_u_;

    //! @}

    double n_old_[3];  // normal of previous time step

   private:
    //! Operator= (disabled)
    CoNodeDataContainer operator=(const CoNodeDataContainer& old) = delete;

    //! Copy constructor (disabled)
    CoNodeDataContainer(const CoNodeDataContainer& old) = delete;
  };

  namespace AUG
  {
    /** \brief Nodal data container for all contact contributions related to the
     *  augmented framework
     *
     *  \author hiermeier \date 02/18 */
    class NodeDataContainer
    {
     public:
      enum SideType
      {
        slave,
        master
      };

      /// @name General methods
      /// @{

      /// constructor
      NodeDataContainer(CoNode& parentNode);

      /// constructor
      NodeDataContainer(
          CoNode& parentNode, const int slMaElementAreaRatio, const bool isTriangleOnMaster);

      virtual ~NodeDataContainer(){};

      void Setup();

      /*! \brief Pack this class so that it can be communicated
       *
       *  This function packs the datacontainer. This is only called
       *  when the class has been initialized and the pointer to this
       *  class exists. */
      void Pack(DRT::PackBuffer& data) const;

      /*! \brief Unpack data from a vector into this class
       *
       *  This function unpacks the data container. This is only called
       *  when the class has been initialized and the pointer to this
       *  class exists. */
      void Unpack(std::vector<char>::size_type& position, const std::vector<char>& data);

      void Complete();
      /// @}

      /// @name accessors
      /// @{

      /*! \brief Return the weighted gap (scalar) of this node */
      inline double& GetWGap() { return wGap_; }
      inline double GetWGap() const { return wGap_; }

      /*! \brief Return the scaling factor kappa for this node
       *
       *  This factor is equivalent to the tributary area of the (potentially)
       *  active nodes. The term is evaluated for all found slave and master
       *  element pairs. Therefore, contributions of non-projectable GP are
       *  missing. This stands in contrast to the augA quantities!*/
      inline double& GetKappa() { return kappa_; }
      inline double GetKappa() const { return kappa_; }

      /*! \brief Return the linearization of the scaling factor kappa
       *
       *  Only potentially active parts are considered. See note for
       *  GetKappa(). */
      inline Deriv1stMap& GetDeriv1st_Kappa() { return d_kappa_; }
      inline const Deriv1stMap& GetDeriv1st_Kappa() const { return d_kappa_; }

      /*! \brief Return the 2-nd order derivative of the tributary area
       *
       *  Only potentially active parts are considered. See note for
       *  GetKappa(). */
      inline Deriv2ndMap& GetDeriv2nd_Kappa() { return dd_kappa_; }
      inline const Deriv2ndMap& GetDeriv2nd_Kappa() const { return dd_kappa_; }

      /*! \brief Return the scaling factor augA
       *
       *  Integration over the whole slave interface, without consideration of
       *  the segments or projections. */
      inline double& GetAugA() { return augA_; }
      inline double GetAugA() const { return augA_; }

      /*! \brief Return the first order derivative of the scaling factor augA
       *
       *  Integration over the whole slave interface, without consideration of
       *  the segments or projections. */
      inline Deriv1stMap& GetDeriv1st_A() { return d_augA_; }
      inline const Deriv1stMap& GetDeriv1st_A() const { return d_augA_; }

      /*! \brief Return the second order derivative of the scaling factor augA
       *
       *  Integration over the whole slave interface, without consideration of
       *  the segments or projections. */
      inline Deriv2ndMap& GetDeriv2nd_A() { return dd_augA_; }
      inline const Deriv2ndMap& GetDeriv2nd_A() const { return dd_augA_; }

      /*! \brief Return the 1-st order derivative of smooth averaged unit normal */
      inline Deriv1stVecMap& GetDeriv1st_N() { return d_avg_unit_normal_; }
      inline const Deriv1stVecMap& GetDeriv1st_N() const { return d_avg_unit_normal_; }

      /*! \brief Return the 2-nd order derivative of smooth averaged unit normal */
      inline Deriv2ndVecMap& GetDeriv2nd_N() { return dd_avg_unit_normal_; }
      inline const Deriv2ndVecMap& GetDeriv2nd_N() const { return dd_avg_unit_normal_; }

      /// @}

      inline Deriv1stMap& GetDeriv1st_WGapSl() { return d_wgap_sl_; }
      inline const Deriv1stMap& GetDeriv1st_WGapSl() const { return d_wgap_sl_; }

      inline Deriv1stMap& GetDeriv1st_WGapMa() { return d_wgap_ma_; }
      inline const Deriv1stMap& GetDeriv1st_WGapMa() const { return d_wgap_ma_; }

      inline Teuchos::RCP<Deriv1stMap>& GetDeriv1st_WGapSl_Complete_Ptr()
      {
        return d_wgap_sl_complete_;
      }

      inline Deriv1stMap& GetDeriv1st_WGapSl_Complete()
      {
        if (d_wgap_sl_complete_.is_null())
          dserror(
              "The complete 1-st order derivative of the weighted gap "
              "(slave) is not initialized! [NULL pointer]");
        return *d_wgap_sl_complete_;
      }

      inline const Deriv1stMap& GetDeriv1st_WGapSl_Complete() const
      {
        if (d_wgap_sl_complete_.is_null())
          dserror(
              "The complete 1-st order derivative of the weighted gap "
              "(slave) is not initialized! [NULL pointer]");
        return *d_wgap_sl_complete_;
      }

      inline Teuchos::RCP<Deriv1stMap>& GetDeriv1st_WGapMa_Complete_Ptr()
      {
        return d_wgap_ma_complete_;
      }

      inline Deriv1stMap& GetDeriv1st_WGapMa_Complete()
      {
        if (d_wgap_ma_complete_.is_null())
          dserror(
              "The complete 1-st order derivative of the weighted gap "
              "(master) is not initialized! [NULL pointer]");
        return *d_wgap_ma_complete_;
      }

      inline const Deriv1stMap& GetDeriv1st_WGapMa_Complete() const
      {
        if (d_wgap_ma_complete_.is_null())
          dserror(
              "The complete 1-st order derivative of the weighted gap "
              "(master) is not initialized! [NULL pointer]");
        return *d_wgap_ma_complete_;
      }

      inline Deriv2ndMap& GetDeriv2nd_WGapSl() { return dd_wgap_sl_; }
      inline const Deriv2ndMap& GetDeriv2nd_WGapSl() const { return dd_wgap_sl_; }

      inline Deriv2ndMap& GetDeriv2nd_WGapMa() { return dd_wgap_ma_; }
      inline const Deriv2ndMap& GetDeriv2nd_WGapMa() const { return dd_wgap_ma_; }

      /// @name Debug nodal values
      /// @{

      inline std::pair<int, double>& Get_Debug() { return debug_.x_; }
      inline const std::pair<int, double>& Get_Debug() const { return debug_.x_; }

      inline std::vector<std::pair<int, double>>& Get_DebugVec() { return debug_.x_vec_; }
      inline const std::vector<std::pair<int, double>>& Get_DebugVec() const
      {
        return debug_.x_vec_;
      }

      inline Deriv1stMap& GetDeriv1st_Debug() { return debug_.d_; }
      inline const Deriv1stMap& GetDeriv1st_Debug() const { return debug_.d_; }

      inline Deriv2ndMap& GetDeriv2nd_Debug() { return debug_.dd_; }
      inline const Deriv2ndMap& GetDeriv2nd_Debug() const { return debug_.dd_; }

      inline Deriv1stVecMap& GetDeriv1st_DebugVec() { return debug_.d_vec_; }
      inline const Deriv1stVecMap& GetDeriv1st_DebugVec() const { return debug_.d_vec_; }

      inline Deriv2ndVecMap& GetDeriv2nd_DebugVec() { return debug_.dd_vec_; }
      inline const Deriv2ndVecMap& GetDeriv2nd_DebugVec() const { return debug_.dd_vec_; }

      /// @}

      inline int GetNumMEntries() const { return mentries_; }

      inline void SetNumMEntries(int mentries) { mentries_ = mentries; }

     private:
      int ApproximateMEntries(const int slMaElementAreaRatio, const bool isTriangleOnMaster) const;

     private:
      // *** Augmented Lagrangian formulation ********************
      int mentries_;

      // nodal scalar values
      double kappa_;  ///< gap-scaling factor kappa
      double wGap_;   ///< nodal entry of weighted gap vector
      double augA_;   ///< nodal scaling factor

      /// @name Linearization
      /// @{

      /// 1-st order derivative of the tributary area (active part)
      Deriv1stMap d_kappa_;

      /// 2-nd order derivative of the tributary area (active part)
      Deriv2ndMap dd_kappa_;

      /// 1-st order derivative of the tributary area (inactive part)
      Deriv1stMap d_augA_;

      /// 2-nd order derivative of the tributary area (inactive part)
      Deriv2ndMap dd_augA_;

      // 1-st derivative of the smooth averaged nodal unit normal
      Deriv1stVecMap d_avg_unit_normal_;

      // 2-nd derivative of the smooth averaged nodal unit normal
      Deriv2ndVecMap dd_avg_unit_normal_;

      /// @}

      /** In correspondence with the chosen integration strategy this variable
       *  contains the complete or the incomplete contributions of the slave side
       *  to the weighted gap */
      Deriv1stMap d_wgap_sl_;
      Deriv1stMap d_wgap_ma_;

      /** This variable contains always the complete contributions of the slave
       *  side to the weighted gap */
      Teuchos::RCP<Deriv1stMap> d_wgap_sl_complete_;
      Teuchos::RCP<Deriv1stMap> d_wgap_ma_complete_;

      Deriv2ndMap dd_wgap_sl_;
      Deriv2ndMap dd_wgap_ma_;

      /// internal debug data structure, holding debug information only
      struct Debug
      {
        Debug() : d_(0), dd_(0), d_vec_(0), dd_vec_(0){};

        ~Debug(){};

        void Complete();

        std::pair<int, double> x_;
        std::vector<std::pair<int, double>> x_vec_;

        Deriv1stMap d_;
        Deriv2ndMap dd_;

        Deriv1stVecMap d_vec_;
        Deriv2ndVecMap dd_vec_;
      };

      Debug debug_;

      /// constant reference to the parent node
      const CoNode& parentNode_;

    };  // class NodeDataContainer
  }     // namespace AUG

  // class CoNodeDataContainer

  /*!
  \brief A class containing additional data from poro contact nodes

  This class contains additional information from poro contact nodes which are
  are not needed for normal contact nodes.
  CoNodePoroDataContainer must be declared before the CoNode itself.

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
    virtual std::map<int, double>& GetDerivnCoup() { return derivncoup_; }

    /*!
    \brief Return the 'VelDerivnCoup' map of this node

    This map contains the derivatives of the node's
    normal coupling condition with respect to the velocities. (for one sided contact just slave!)

    */
    virtual std::map<int, double>& GetVelDerivnCoup() { return velderivncoup_; }

    /*!
    \brief Return the 'PresDerivnCoup' map of this node

    This map contains the derivatives of the node's
    normal coupling condition with respect to the pressures. (for one sided contact just slave!)
    //h.Willmann
    */
    virtual std::map<int, double>& GetPresDerivnCoup() { return presderivncoup_; }

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
    CoNodePoroDataContainer operator=(const CoNodePoroDataContainer& old) = delete;
    CoNodePoroDataContainer(const CoNodePoroDataContainer& old) = delete;


    double ncouprow_;  // nodal entry of n-coupling vector
    std::map<int, double>
        derivncoup_;  // directional derivative of nodal weighted n-coupling vector
    std::map<int, double>
        velderivncoup_;  // velocity derivative of nodal weighted n-coupling vector
    std::map<int, double>
        presderivncoup_;  // pressure derivative of nodal weighted n-coupling vector //h.Willmann

    double fvel_[3];    // fluid velocity for porous problem
    double fpres_;      // fluid pressure for porous problem
    double svel_[3];    // structural velocity for porous problem
    double porolm_[3];  // lagrange multiplier from poro no penetration condition!
  };                    // class CoNodePoroDataContainer

  /*!
  \brief A class containing additional data from tsi contact nodes

  This class contains additional information from tsi contact nodes which are
  are not needed for normal contact nodes.

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

    std::map<int, double>& DerivTempMasterDisp() { return derivTempMasterDisp_; }
    std::map<int, double>& DerivTempMasterTemp() { return derivTempMasterTemp_; }

    void Clear();


   protected:
    // don't want = operator and cctor
    CoNodeTSIDataContainer operator=(const CoNodeTSIDataContainer& old) = delete;
    CoNodeTSIDataContainer(const CoNodeTSIDataContainer& old) = delete;

    double temp_;
    double t_ref_;
    double t_dam_;
    double thermo_lm_;

    double temp_master_;
    std::map<int, double> derivTempMasterDisp_;
    std::map<int, double> derivTempMasterTemp_;

  };  // class CoNodeTSIDataContainer

  /*!
  \brief A class containing additional data from EHL contact nodes

  This class contains additional information from EHL contact nodes which are
  are not needed for normal contact nodes.

  */
  class CoNodeEhlDataContainer
  {
   public:
    //! empty c-tor
    CoNodeEhlDataContainer() {}

    //! d-tor
    virtual ~CoNodeEhlDataContainer(){};

    //! pack for parallel communication
    virtual void Pack(DRT::PackBuffer& data) const {
        /* no need to pack, since terms are re-evaluated after parallel communication */};

    //! unpack and re-init after parallel comunication
    virtual void Unpack(std::vector<char>::size_type& position, const std::vector<char>& data){
        /* no need to pack, since terms are re-evaluated after parallel communication */};

    //! clear all stored data
    void Clear()
    {
      weighted_relTangVel_.Clear();
      deriv_weighted_relTangVel_.clear();
      weighted_avTangVel_.Clear();
      deriv_weighted_avTangVel_.clear();
      tang_grad_.clear();
      tang_grad_deriv_.clear();
      weighted_relTangVel_.Clear();
      weighted_avTangVel_.Clear();
    }

    LINALG::Matrix<3, 1>& GetWeightedRelTangVel() { return weighted_relTangVel_; }
    std::unordered_map<int, LINALG::Matrix<3, 1>>& GetWeightedRelTangVelDeriv()
    {
      return deriv_weighted_relTangVel_;
    }

    LINALG::Matrix<3, 1>& GetWeightedAvTangVel() { return weighted_avTangVel_; }
    std::unordered_map<int, LINALG::Matrix<3, 1>>& GetWeightedAvTangVelDeriv()
    {
      return deriv_weighted_avTangVel_;
    }

    std::unordered_map<int, LINALG::Matrix<3, 1>>& GetSurfGrad() { return tang_grad_; }
    std::unordered_map<int, std::unordered_map<int, LINALG::Matrix<3, 1>>>& GetSurfGradDeriv()
    {
      return tang_grad_deriv_;
    }

   protected:
    // don't want = operator and cctor
    CoNodeEhlDataContainer operator=(const CoNodeEhlDataContainer& old) = delete;
    CoNodeEhlDataContainer(const CoNodeEhlDataContainer& old) = delete;

    // actual data
    LINALG::Matrix<3, 1> weighted_relTangVel_;
    std::unordered_map<int, LINALG::Matrix<3, 1>> deriv_weighted_relTangVel_;

    LINALG::Matrix<3, 1> weighted_avTangVel_;
    std::unordered_map<int, LINALG::Matrix<3, 1>> deriv_weighted_avTangVel_;

    std::unordered_map<int, LINALG::Matrix<3, 1>> tang_grad_;
    std::unordered_map<int, std::unordered_map<int, LINALG::Matrix<3, 1>>> tang_grad_deriv_;
  };

  /*!
   \brief A class for a contact node derived from MORTAR::MortarNode

   This class represents a finite element node capable of contact.

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
    virtual ~CoNode() = default;

    /*!
     \brief Return unique ParObject id

     every class implementing ParObject needs a unique id defined at the
     top of drt_lib/drt_parobject.H.

     */
    virtual int UniqueParObjectId() const { return CoNodeType::Instance().UniqueParObjectId(); }

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
     \brief Print this contact node
     */
    virtual void Print(std::ostream& os) const;

    /*!
     \brief Is Node initialized as active node (only slave nodes)
     */
    virtual bool IsInitActive() const
    {
      if (!IsSlave()) dserror("ERROR: InitActive requested for Master node");
      return initactive_;
    }

    /*!
     \brief Modify initial active status of slave node

     This belated modification is necessary to be able to use
     the binary search tree for contact initialization in the
     load-controlled quasi-static case (instead of input file
     information Active/Inactive)

     */
    virtual bool& SetInitActive()
    {
      if (!IsSlave()) dserror("ERROR: InitActive requested for Master node");
      return initactive_;
    }

    /*!
     \brief Return contact status of this node (active=true)
     */
    virtual bool& Active() { return active_; }

    virtual bool& InvolvedM() { return involvedm_; }

    /*!
     \brief Return data container of this node

     This method returns the data container of this node where additional
     contact specific quantities/information are stored.

     */
    inline CONTACT::CoNodeDataContainer& CoData() { return *codata_; }
    inline CONTACT::CoNodeDataContainer& CoData() const { return *codata_; }

    /*! \brief Return data container of this augmented contact node
     *
     *  This method returns the data container of this node where additional
     *  augmented contact specific quantities/information are stored.
     *
     *  \author hiermeier \date 03/17 */
    inline CONTACT::AUG::NodeDataContainer& AugData()
    {
      if (augdata_.is_null()) dserror("There are no augmented contact node data available!");
      return *augdata_;
    };
    inline const CONTACT::AUG::NodeDataContainer& AugData() const
    {
      if (augdata_.is_null()) dserror("There are no augmented contact node data available!");
      return *augdata_;
    };

    inline CONTACT::CoNodePoroDataContainer& CoPoroData() { return *coporodata_; }
    inline CONTACT::CoNodeTSIDataContainer& CoTSIData() { return *cTSIdata_; }
    inline CONTACT::CoNodeEhlDataContainer& CoEhlData() { return *cEHLdata_; }

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
    void AddltlGapValue(double* val);

    /*!
     \brief Add a value to the point-wise jump of this node (LTL)

     \param val : value to be added

     */
    void AddltlJumpValue(double* val);

    /*!
    \brief Add a value to the weighted gap of this node (augmented Lagrange)

    This value is later assembled to the averaged weighted gap vector.
    Note that grow_ here is a scalar.

    \param val : value to be added

    */
    void AddWGapValue(const double val);

    /*!
    \brief Add a value to the scaling factor kappa of this node (augmented Lagrange)

    This value is later assembled needed during the assembling of the averaged
    weighted gap. Note that kappa_ is a scalar.

    \param val : value to be added
    */
    void AddKappaValue(double& val);

    /*!
    \brief TODO
    */
    void AddWcLm(double val);

    /*!
    \brief TODO
    */
    void AddWcSuLm(int col, double val);

    /*!
    \brief TODO
    */
    void AddWcMuLm(int col, double val);

    /*!
    \brief Add a value to the map of LM derivatives of this node

    The 'DerivZ' map is later assembled to the global DerivZ matrix.
    Note that derivz_ here is a vector.

    Note: This is only calculated when performing a penalty strategy

    \param row : local dof row id to add to (row-wise)
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
     \brief Build nodal edge tangent
     */
    void BuildAveragedEdgeTangent();

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
    void InitializeAugDataContainer(const int slMaElementAreaRatio, const bool isTriangleOnMaster);

    /*!
    \brief Initializes the poro data container of the node

    With this function, the container with contact specific quantities/information
    is initialized. --- Used to initialize PoroDataContainer for master nodes!

    */
    virtual void InitializePoroDataContainer();

    /*!
    \brief Initializes the ehl data container of the node

    With this function, the container with contact specific quantities/information
    is initialized.

    */
    virtual void InitializeEhlDataContainer();

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

    //! @}

    /*! @name Empty functions (friction only)
     *
     * All these functions only have functionality for friction nodes, thus they are
     * defined as empty here in the general mortar node. They can be called whenever you like.
     */

    virtual void AddSNode(int node) {}
    virtual void AddMNode(int node) {}

    //! @}

    virtual bool HasCoTSIData() { return (cTSIdata_ != Teuchos::null); }

    /*!
     \brief Write nodal normals to old ones
     */
    void StoreOldNormal();

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
    void DerivAveragedNormal(Epetra_SerialDenseMatrix& elens, double length, double ltxi);

   protected:
    //! true if contact node is in contact (active set strategy)
    bool active_;

    //! true if node is initialized as active node
    bool initactive_;

    //! cnode is an master node for which the mortar integration is performed (both-sided wear)
    bool involvedm_;

    //! number of lin entries per direction vector
    int linsize_;

    //! Additional information of proc's contact nodes
    Teuchos::RCP<CONTACT::CoNodeDataContainer> codata_;

    //! Additional information of proc's augmented contact nodes
    Teuchos::RCP<CONTACT::AUG::NodeDataContainer> augdata_;

    //! Additional information of proc's poro contact nodes
    Teuchos::RCP<CONTACT::CoNodePoroDataContainer> coporodata_;

    //! Additional information of TSI contact nodes
    Teuchos::RCP<CONTACT::CoNodeTSIDataContainer> cTSIdata_;

    //! Additional information of EHL contact nodes
    Teuchos::RCP<CONTACT::CoNodeEhlDataContainer> cEHLdata_;
  };
  // class CoNode
}  // namespace CONTACT

// << operator
std::ostream& operator<<(std::ostream& os, const CONTACT::CoNode& cnode);

#endif  // #ifndef CNODE_H
