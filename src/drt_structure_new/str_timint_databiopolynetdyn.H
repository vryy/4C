/*---------------------------------------------------------------------*/
/*!
\file str_timint_databiopolynetdyn.H

\brief Statmech dynamics data container for the structural (time)
       integration

\maintainer Jonas Eichinger

\date Jun 22, 2016

\level 3

*/
/*---------------------------------------------------------------------*/

#ifndef STR_TIMINT_DATABIOPOLYNETDYN_
#define STR_TIMINT_DATABIOPOLYNETDYN_

#include "../drt_lib/drt_dserror.H"
#include <NOX_Abstract_Vector.H>
#include "../drt_lib/drt_discret.H"
#include "../drt_inpar/inpar_statmech.H"

// forward declaration
namespace STR {
namespace TIMINT {
  class BaseDataGlobalState;
  class BaseDataSDyn;

  /*! StatMech data container for the statistical mechanic model evaluation procedure.
   * this data container holds all variables from the statmech section of the input file.
   * It is used by the modelevaluator, statmech manager as well as on element level. On
   * element level only these functions can be called that are implemented in the
   * STATMECH ParamsInterface class an overloaded here.
   * \author Jonas Eichinger
   * \date 06/16 */
  class DataSMDyn
  {
  public:
    //! constructor
    DataSMDyn();

    //! destructor
    virtual ~DataSMDyn() {};
    //---------------------------------------------------------------------------
    // @name set routines which are used to set the parameters of the data container
    //
    //  These functions are NOT allowed to be called by the elements
    //---------------------------------------------------------------------------
    //! @{
    //! initialize the stuff coming from outside
    void Init(const Teuchos::RCP<const STR::TIMINT::BaseDataSDyn>& str_data_ptr);

    //! setup member variables
    void Setup();

    //! returns the isinit_ flag
    inline const bool& IsInit() const
    { return isinit_; };

    //! returns the issetup_ flag
    inline const bool& IsSetup() const
    { return issetup_; };

    //! Checks the init and setup status
    inline void CheckInitSetup() const
    {
      if (!IsInit() or !IsSetup())
        dserror("Call Init() and Setup() first!");
    }

    //! Checks the init status
    inline void CheckInit() const
    {
      if (!IsInit())
        dserror("Init() has not been called, yet!");
    }

    //! get basedatasdyn
    const Teuchos::RCP<const BaseDataSDyn>& DataSDyn() const
    {
      CheckInitSetup();
      return data_sdyn_ptr_;
    }

    /// return flag for statmech model
    const bool& StatMechProb() const
    {
      CheckInitSetup();
      return statmechprob_;
    };

    /// return type of predefined neumann boundary condition
    const enum INPAR::STATMECH::NBCType& NbcType() const
    {
      CheckInitSetup();
      return nbctype_;
    };

      /// return type of network
    const enum INPAR::STATMECH::NetworkType& NetworkType() const
    {
      CheckInitSetup();
      return networktype_;
    };

    /// return resolution for search routine
    const int& SearchRes() const
    {
      CheckInitSetup();
      return searchres_;
    };

    /// return resolution for search routine
    const Teuchos::RCP<std::vector<int> >& SearchResDir() const
    {
      CheckInitSetup();
      return searchresdir_;
    };

    /// return flag to check orientation between two filaments before linking them
    const bool& CheckOrient() const
    {
      CheckInitSetup();
      return checkorient_;
    };

    /// binding site search method: volpart, binning, octree
    const enum INPAR::STATMECH::BSSearchType& BindingSiteSearch() const
    {
      CheckInitSetup();
      return bindingsitesearch_;
    };

    /// return flag to check orientation between two filaments before linking them
    const double& MaxRandRorce() const
    {
      CheckInitSetup();
      return maxrandforce_;
    };

    /// flag for gmsh output
    const bool& GmshOutput() const
    {
      CheckInitSetup();
      return gmshoutput_;
    };

    /// number of interpolation points for higher order plotting
    const int& GmshNumIntpt() const
    {
      CheckInitSetup();
      return gmshnumintpt_;
    };

    /// flag for visualization of structure detection
    const bool& GmshNetStruct() const
    {
      CheckInitSetup();
      return gmshnetstruct_;
    };

    /// factor for thickness
    const double& PlotFactorThick() const
    {
      CheckInitSetup();
      return plotfactorthick_;
    };

    /// flag for special output
    const enum INPAR::STATMECH::StatOutput& SpecialOutput() const
    {
      CheckInitSetup();
      return specialoutput_;
    };

    /// output interval for gmsh output
    const int& GmshOutInterval() const
    {
      CheckInitSetup();
      return gmshoutinterval_;
    };

    /// output intervall for special output
    const int& OutputInterval() const
    {
      CheckInitSetup();
      return outputinterval_;
    };

    /// start time for writing output
    const double& StartTimeOut() const
    {
      CheckInitSetup();
      return starttimeout_;
    };

    /// start actions equilibration time, kt switch time,...,starttimeact
    const Teuchos::RCP<std::vector<double> >& ActionTime() const
    {
      CheckInitSetup();
      return actiontime_;
    };

    /// time step sizes corresponding to actiontime
    const Teuchos::RCP<std::vector<double> >& ActionDt() const
    {
      CheckInitSetup();
      return actiondt_;
    };

    /// controls entry in actiontime at which dbcs are applied
    const int& BcTimeIndex() const
    {
      CheckInitSetup();
      return bctimeindex_;
    };

    /// time when action is supposed to happen
    const double& StartTimeAction() const
    {
      CheckInitSetup();
      return actiontime_->at(bctimeindex_);
    };

    /// controls entry in actiontime at which dbcs are applied
    const double& TimeIntConstRandNumb_() const
    {
      CheckInitSetup();
      return timeintconstrandnumb_;
    };

    /// number of bins for specialoutput densitydensitycorr
    const int& HistogramBins() const
    {
      CheckInitSetup();
      return histogrambins_;
    };

    /// creates numrasterpoints^3 rasterpoints, preferably numrasterpoints=2*n+1, n=1,2,3,...
    const int& NumrasterPoints() const
    {
      CheckInitSetup();
      return numrasterpoints_;
    };

    /// thermal energy
    const double& KT() const
    {
      CheckInitSetup();
      return kt_;
    };

    /// thermal energy
    void UpdateKT()
    {
      CheckInitSetup();
      kt_ = ktact_;
      return;
    };

    /// thermal energy after starttimeact_
    const double& KTAct() const
    {
      CheckInitSetup();
      return ktact_;
    };

    /// chemical association-rate before ktswitchtime
    const double& KOnStart() const
    {
      CheckInitSetup();
      return k_on_start_;
    };

    /// chemical dissociation-rate after ktswitchtime // 100.0 1000.0
    const double& KOnEnd() const
    {
      CheckInitSetup();
      return k_on_end_;
    };

    /// chemical association-rate in case of linker with both ends on the same filament
    const double& KOnSelf() const
    {
      CheckInitSetup();
      return k_on_self_;
    };

    /// in analogy to k_on
    const double& KOffStart() const
    {
      CheckInitSetup();
      return k_off_start_;
    };

    /// in analogy to k_on
    const double& KOffEnd() const
    {
      CheckInitSetup();
      return k_off_end_;
    };

    /// applied time curve for nbc after starttimeact_;
    const int& NbcCurveNumber() const
    {
      CheckInitSetup();
      return nbccurvenumber_;
    };

    /// absolute value of applied constant creep force
    const double& NbcForceAmp() const
    {
      CheckInitSetup();
      return nbcforceamp_;
    };

    /// number of nbc nodes to which point forces are applied in sequence
    const int& NumNbcNodes() const
    {
      CheckInitSetup();
      return numnbcnodes_;
    };

    /// filament model: std, helical
    const enum INPAR::STATMECH::FilamentModel& FilamentModel() const
    {
      CheckInitSetup();
      return filamentmodel_;
    };

    /// toggles filament polarity
    const bool& FilamentPolarity() const
    {
      CheckInitSetup();
      return filamentpolarity_;
    };

    /// number of initial beam elements
    const int& NumEvalElements() const
    {
      CheckInitSetup();
      return num_eval_elements_;
    };

    /// rise per monomer in the actin helix
    const Teuchos::RCP<std::vector<double> >& RisePerBinSpot() const
    {
      CheckInitSetup();
      return riseperbinspot_;
    };

    ///  rotation per monomer in the actin helix
    const double& RotPerBinSpot() const
    {
      CheckInitSetup();
      return rotperbinspot_;
    };

    /// angle betw. binding spot orientation of dist. vector bindinspot->crosslinke
    const double& PhiBinSpot() const
    {
      CheckInitSetup();
      return phibinspot_;
    };

    /// every n-th binding spot is available for binding
    const int& BinSpotInterval() const
    {
      CheckInitSetup();
      return binspotinterval_;
    };

    /// get linker model: none, std, bellseq, active, stdintpol
    const enum INPAR::STATMECH::LinkerModel& LinkerModel() const
    {
      CheckInitSetup();
      return linkermodel_;
    };

    /// constrain brownian motion of linkers to 2d
    const bool& PlaneLinkerMotion() const
    {
      CheckInitSetup();
      return planelinkermotion_;
    };

    /// constrain brownian motion of linkers to 2d
    const double& ActiveLinkerFraction() const
    {
      CheckInitSetup();
      return activelinkerfraction_;
    };

    /// number of crosslinkers in the simulated volume
    const int& NumCrosslink() const
    {
      CheckInitSetup();
      return numcrosslink_;
    };

    /// number of crosslinker elements set up before the first time step
    const int& InitOccupiedBinSpots() const
    {
      CheckInitSetup();
      return initoccupiedbinspots_;
    };

    /// number of substrate filaments
    const int& NumSubstrateFil() const
    {
      CheckInitSetup();
      return numsubstratefil_;
    };

    /// distance between the two binding domains of a linker
    const double& RLink() const
    {
      CheckInitSetup();
      return r_link_;
    };

    /// deviation from r_link
    const double& DeltaRLink() const
    {
      CheckInitSetup();
      return deltar_link_;
    };

    /// assuming circular cross section of crosslinker with diameter according to phd tharmann
    const double& ALink() const
    {
      CheckInitSetup();
      return alink_;
    };

    /// must set ilink=0.0 ->truss!
    const double& ILink() const
    {
      CheckInitSetup();
      return ilink_;
    };

    /// 2xilink
    const double& IpLink() const
    {
      CheckInitSetup();
      return iplink_;
    };

//    /// torsion stiffness from nishizaka 2000 and ilink from phd tharmann
//    const double& ELink() const
//    {
//      CheckInitSetup();
//      return elink_;
//    };

    /// length compression factor for active crosslinkers
    const double& LinkerScaleFactor() const
    {
     CheckInitSetup();
     return linkerscalefactor_;
    };

    /// stroke distance of an active linker
    const double& StrokeDistance() const
    {
     CheckInitSetup();
     return strokedistance_;
    };

    /// deltad in bell's equation, change in reaction variable
    const double& DeltaBellSeq() const
    {
     CheckInitSetup();
     return deltabellseq_;
    };

    /// stiffness of the potential which determines probability of a certain crosslinker orientation
    const double& COrient() const
    {
     CheckInitSetup();
     return corient_;
    };

    /// preferred binding angle betw. two filaments
    const double& PhiZero() const
    {
     CheckInitSetup();
     return phizero_;
    };

    /// deviation from phi0
    const double& PhiZeroDev() const
    {
     CheckInitSetup();
     return phizerodev_;
    };

    /// active linker time cycle (myosin)
    const double& ActiveLinkerCycle() const
    {
     CheckInitSetup();
     return activelinkercycle_;
    };

    /// active linker recovery time fraction (myosin)
    const double& ActiveRecoveryFraction() const
    {
     CheckInitSetup();
     return activerecoveryfraction_;
    };

    /// rate from long to short
    const double& KActiveShortStart() const
    {
     CheckInitSetup();
     return k_active_short_start_;
    };

    /// in analogy to K_ACT_SHORT_start
    const double& KActiveShortEnd() const
    {
     CheckInitSetup();
     return k_active_short_end_;
    };

    /// rate from short to long
    const double& KActiveLongStart() const
    {
     CheckInitSetup();
     return k_active_long_start_;
    };

    /// in analogy to K_ACT_LONG_start
    const double& KActiveLongEnd() const
    {
     CheckInitSetup();
     return k_active_long_end_;
    };

    /// reduce the linker number by this number at starttimeact_
    const int& ReduceCrossLinksBy() const
    {
     CheckInitSetup();
     return reducecrosslinksby_;
    };

    /// toggles rotation of motors in addition to contraction
    const bool& CrossBridgeModel() const
    {
     CheckInitSetup();
     return crossbridgemodel_;
    };
    //! @}

    /// return type of predefined dirichlet boundary condition
    const enum INPAR::STATMECH::DBCType& DbcType() const
    {
      CheckInitSetup();
      return dbctype_;
    };

    /// edge lengths of the cuboid simulation volume
    const Teuchos::RCP<std::vector<double> >& PeriodLength() const
    {
      CheckInitSetup();
      return periodlength_;
    };

    /// amplitude of shear in micrometers
    const double& ShearAmplitude() const
    {
      CheckInitSetup();
      return shearamplitude_;
    };

    /// ~ 1e-3 / 2.27 according to cyron2011 eq 52 ff, viscosity of surrounding fluid
    const double& Eta() const
    {
      CheckInitSetup();
      return eta_;
    };

    /// direction of dbc displacement after starttimeact_; (x==1, y==2, z==3)
    const int& DbcDispDir() const
    {
      CheckInitSetup();
      return dbcdispdir_;
    };

    /// applied time curve for dbc after starttimeact_;
    const int& CurveNumber() const
    {
      CheckInitSetup();
      return curvenumber_;
    };

  private:

    bool isinit_;

    bool issetup_;

    //---------------------------------------------------------------------------
    /// parameter list of statmech dynamics (read only)
    //---------------------------------------------------------------------------
    Teuchos::RCP<const Teuchos::ParameterList> smdynparams_;
    Teuchos::RCP<const BaseDataSDyn> data_sdyn_ptr_;
    //---------------------------------------------------------------------------
    // setup statmech model
    //---------------------------------------------------------------------------
    bool statmechprob_;
    //---------------------------------------------------------------------------
    // general statmech simulation flags
    //---------------------------------------------------------------------------
    // none shearfixed sheartrans pinnodes,  choose dbc type
    enum INPAR::STATMECH::DBCType dbctype_;
    // constcreep, choose nbc type
    enum INPAR::STATMECH::NBCType nbctype_;
    // specify network / simulation type
    enum INPAR::STATMECH::NetworkType networktype_;
    // resolution for search routine linkers x binding spots
    int searchres_;
    // resolution for search routine linkers x binding spots in each direction
    Teuchos::RCP<std::vector<int> > searchresdir_;
    // check orientation between two filaments before linking them
    bool checkorient_;
    // binding site search method: volpart, binning, octree
    enum INPAR::STATMECH::BSSearchType bindingsitesearch_;
    // any random force beyond MAXRANDFORCE*(standdev) will be omitted and redrawn. -1.0 means no bounds.
    double maxrandforce_;
    //---------------------------------------------------------------------------
    // output set up
    //---------------------------------------------------------------------------
    // visualization via gmsh
    bool gmshoutput_;
    // number of interpolation points for higher order plotting (in case of beam3eb)
    int gmshnumintpt_;
    // visualize structure detection
    bool gmshnetstruct_;
    // 0 = gmsh line output, >0 = gmsh wedge output (c * (real thickness))
    double plotfactorthick_;
    // structanaly viscoelasticity
    enum INPAR::STATMECH::StatOutput specialoutput_;
    // gmsh output interval
    int gmshoutinterval_;
    // special output interval
    int outputinterval_;
    //---------------------------------------------------------------------------
    // event triggers
    //---------------------------------------------------------------------------
    // start timt for writing output
    double starttimeout_;
    // start actions equilibration time, kt switch time,...,starttimeact
    Teuchos::RCP<std::vector<double> > actiontime_;
    // time step sizes corresponding to actiontime
    Teuchos::RCP<std::vector<double> > actiondt_;
    // controls entry in actiontime at which dbcs are applied
    int bctimeindex_;
    // within this time interval the random numbers remain constant. -1.0 means no prescribed time interval
    double timeintconstrandnumb_;
    //---------------------------------------------------------------------------
    // simulation specs
    //---------------------------------------------------------------------------
    // edge lengths of the cuboid volume
    Teuchos::RCP<std::vector<double> > periodlength_;
    // number of bins for specialoutput densitydensitycorr
    int histogrambins_;
    // creates numrasterpoints^3 rasterpoints, preferably numrasterpoints=2*n+1, n=1,2,3,...
    int numrasterpoints_;
    // amplitude of shear in micrometers
    double shearamplitude_;
    // ~ 1e-3 / 2.27 according to cyron2011 eq 52 ff, viscosity of surrounding fluid
    double eta_;
    // thermal energy
    double kt_;
    // thermal energy after starttimeact_;
    double ktact_;
    // chemical association-rate before ktswitchtime
    double k_on_start_;
    // chemical dissociation-rate after ktswitchtime // 100.0 1000.0
    double k_on_end_;
    // chemical association-rate in case of linker with both ends on the same filament
    double k_on_self_;
    // in analogy to k_on
    double k_off_start_;
    // in analogy to k_on
    double k_off_end_;
    // direction of dbc displacement after starttimeact_; (x==1, y==2, z==3)
    int dbcdispdir_;
    // applied time curve for dbc after starttimeact_;
    int curvenumber_;
    // applied time curve for nbc after starttimeact_;
    int nbccurvenumber_;
    // absolute value of applied constant creep force
    double nbcforceamp_;
    // number of nbc nodes to which point forces are applied in sequence
    int numnbcnodes_;
    //---------------------------------------------------------------------------
    // filament model
    //---------------------------------------------------------------------------
    // select filament model: std, helical
    enum INPAR::STATMECH::FilamentModel filamentmodel_;
    // toggles filament polarity
    bool filamentpolarity_;
    // number of (initial) beam elements
    int num_eval_elements_;
    // rise per monomer in the actin helix
    Teuchos::RCP<std::vector<double> > riseperbinspot_;
    // rotation per monomer in the actin helix
    double rotperbinspot_;
    // angle betw. binding spot orientation of dist. vector bindinspot->crosslinker
    double phibinspot_;
    // every n-th binding spot is available for binding
    int binspotinterval_;
    //---------------------------------------------------------------------------
    // crosslinker model
    //---------------------------------------------------------------------------
    // set linker model: none, std, bellseq, active, stdintpol
    enum INPAR::STATMECH::LinkerModel linkermodel_;
    // constrain brownian motion of linkers to 2d
    bool planelinkermotion_;
    // fraction of crosslinkers which is active
    double activelinkerfraction_;
    // number of crosslinkers in the simulated volume
    int numcrosslink_;
    // number of crosslinker elements set up before the first time step
    int initoccupiedbinspots_;
    // number of substrate filaments
    int numsubstratefil_;
    // distance between the two binding domains of a linker
    double r_link_;
    // deviation from r_link
    double deltar_link_;
    // assuming circular cross section of crosslinker with diameter according to phd tharmann
    double alink_;
    // must set ilink=0.0 ->truss!
    double ilink_;
    // 2xilink
    double iplink_;
    // torsion stiffness from nishizaka 2000 and ilink from phd tharmann
    double elink_;
    // length compression factor for active crosslinkers
    double linkerscalefactor_;
    // stroke distance of an active linker
    double strokedistance_;
    // deltad in bell's equation, change in reaction variable
    double deltabellseq_;
    // stiffness of the potential which determines probability of a certain crosslinker orientation
    double corient_;
    // preferred binding angle betw. two filaments
    double phizero_;
    // deviation from phi0
    double phizerodev_;
    // active linker time cycle (myosin)
    double activelinkercycle_;
    // active linker recovery time fraction (myosin)
    double activerecoveryfraction_;
    // rate from long to short
    double k_active_short_start_;
    // in analogy to K_ACT_SHORT_start
    double k_active_short_end_;
    // rate from short to long
    double k_active_long_start_;
    // in analogy to K_ACT_LONG_start
    double k_active_long_end_;
    // reduce the linker number by this number at starttimeact_
    int reducecrosslinksby_;
    // toggles rotation of motors in addition to contraction
    bool crossbridgemodel_;
  }; // class StatMechData

} // namespace MODELEVALUATOR
} // namespace STR


#endif /* STR_TIMINT_DATABIOPOLYNETDYN_ */
