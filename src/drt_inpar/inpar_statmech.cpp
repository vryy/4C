/*-----------------------------------------------------------*/
/*!
\file inpar_statmech.cpp

\brief input parameter for statistical mechanic problem

\maintainer Jonas Eichinger

\level 2

*/
/*-----------------------------------------------------------*/


#include "drt_validparameters.H"
#include "inpar_statmech.H"


void INPAR::STATMECH::SetValidParameters(Teuchos::RCP<Teuchos::ParameterList> list)
{
  using namespace DRT::INPUT;
  using Teuchos::tuple;
  using Teuchos::setStringToIntegralParameter;

  Teuchos::Array<std::string> yesnotuple = tuple<std::string>("Yes","No","yes","no","YES","NO");
  Teuchos::Array<int> yesnovalue = tuple<int>(true,false,true,false,true,false);

  Teuchos::ParameterList& statmech = list->sublist("STATISTICAL MECHANICS",false,"");

  setStringToIntegralParameter<int>("STATMECHPROB","No",
                                 "Statistical mechanics problem",
                                 yesnotuple,yesnovalue,&statmech);

  setStringToIntegralParameter<int>("CROSSLINKER","No",
                                 "Crosslinker in problem",
                                 yesnotuple,yesnovalue,&statmech);

  //Reading which kind of special output should be written to files
  setStringToIntegralParameter<int>("SPECIALOUTPUT","None","kind of special statistical output data written into files",
                                 //listing possible std::strings in input file in category SPECIALOUTPUT
                                 tuple<std::string>("None","none",
                                                    "endtoend_log",
                                                    "anisotropic",
                                                    "orientationcorrelation",
                                                    "endtoend_const",
//                                                    "viscoelasticity",
//                                                    "networkcreep",
//                                                    "networkrelax",
                                                    "networkdispfield",
                                                    "structanaly",
                                                    "octree",
                                                    "loom",
                                                    "loomelnrg",
                                                    "motassay",
                                                    "linkerlength",
                                                    "deltatheta",
                                                    "vesceqshapes"
                                                    ),
                                 //translating input std::strings into BACI input parameters
                                 tuple<int>(statout_none,statout_none,
                                            statout_endtoendlog,
                                            statout_anisotropic,
                                            statout_orientationcorrelation,
                                            statout_endtoendconst,
//                                            statout_viscoelasticity,
//                                            statout_networkcreep,
//                                            statout_networkrelax,
                                            statout_networkdispfield,
                                            statout_structanaly,
                                            statout_octree,
                                            statout_loom,
                                            statout_loomelnrg,
                                            statout_motassay,
                                            statout_linkerlength,
                                            statout_deltatheta,
                                            statout_vesceqshapes),
                                 &statmech);

  //Reading which kind of Dirichlet boundary condition should be applied
  setStringToIntegralParameter<int>("DBCTYPE","std","Dirichlet BC type applied",
                                 //listing possible std::strings in input file in category DBCTYPE
                                 tuple<std::string>("none",
                                                    "std",
                                                    "shearfixed",
                                                    "shearfixeddel",
                                                    "sheartrans",
                                                    "pinnodes" ,
                                                    "affineshear",
                                                    "affinesheardel",
                                                    "movablesupport1d"),
                                 //translating input std::strings into BACI input parameters
                                 tuple<int>(dbctype_none,
                                            dbctype_std,
                                            dbctype_shearfixed,
                                            dbctype_shearfixeddel,
                                            dbctype_sheartrans,
                                            dbctype_pinnodes,
                                            dbctype_affineshear,
                                            dbctype_affinesheardel,
                                            dbctype_movablesupport1d),
                                            &statmech);

  //Reading which kind of Dirichlet boundary condition should be applied
  setStringToIntegralParameter<int>("NBCTYPE","std","Neumann BC type applied",
                                 //listing possible std::strings in input file in category DBCTYPE
                                 tuple<std::string>("std",
                                                    "constcreep",
                                                    "randompointforce"),
                                 //translating input std::strings into BACI input parameters
                                 tuple<int>(nbctype_std,
                                            nbctype_constcreep,
                                            nbctype_randompointforce),
                                            &statmech);

  //Reading which kind of biopolymer network will be simulated
  setStringToIntegralParameter<int>("NETWORKTYPE","std","Network type simulated",
                                 //listing possible std::strings in input file in category FILAMENTMODEL
                                 tuple<std::string>("std",
                                                    "loom",
                                                    "motassay"),
                                 //translating input std::strings into BACI input parameters
                                 tuple<int>(networktype_std,
                                            networktype_casimir,
                                            networktype_motassay),
                                            &statmech);

  //Reading which kind of linker model should be applied
  setStringToIntegralParameter<int>("LINKERMODEL","none","Linker model applied in Statmech simulations",
                                 //listing possible std::strings in input file in category LINKERMODEL
                                 tuple<std::string>("none",
                                                    "std",
                                                    "stdintpol",
                                                    "bellseq",
                                                    "bellseqintpol",
                                                    "active",
                                                    "activeintpol",
                                                    "myosinthick"),
                                 //translating input std::strings into BACI input parameters
                                 tuple<int>(linkermodel_none,
                                            linkermodel_std,
                                            linkermodel_stdintpol,
                                            linkermodel_bellseq,
                                            linkermodel_bellseqintpol,
                                            linkermodel_active,
                                            linkermodel_activeintpol,
                                            linkermodel_myosinthick),
                                            &statmech);

  // linker are only allowed to move in 2d, third coordinate is not
  setStringToIntegralParameter<int>("PLANELINKERMOTION","No",
                                 "Plane Brownian Motion of linkers",
                                 yesnotuple,yesnovalue,&statmech);

  // toggles rotation of motors in addition to contraction
  setStringToIntegralParameter<int>("CROSSBRIDGEMODEL","No",
                                 "swinging cross bridge model for active linkers",
                                 yesnotuple,yesnovalue,&statmech);

  //Reading which kind of filament model should be applied
  setStringToIntegralParameter<int>("FILAMENTMODEL","std","Filament model applied in Statmech simulations",
                                 //listing possible std::strings in input file in category FILAMENTMODEL
                                 tuple<std::string>("std",
                                                    "helical"),
                                 //translating input std::strings into BACI input parameters
                                 tuple<int>(filamentmodel_std,
                                            filamentmodel_helical),
                                            &statmech);

  //Reading which kind of search routine is used
  setStringToIntegralParameter<int>("BINDINGSITESEARCH","volpart","binding site search method",
                                 //listing possible std::strings in input file in category DBCTYPE
                                 tuple<std::string>("volpart",
                                                    "binning",
                                                    "octree"),
                                 //translating input std::strings into BACI input parameters
                                 tuple<int>(bsstype_volpart,
                                            bsstype_binning,
                                            bsstype_octree),
                                            &statmech);

  //time after which writing of statistical output is started
  DoubleParameter("STARTTIMEOUT",0.0,"Time after which writing of statistical output is started",&statmech);
  // Time values at which certain actions are carried out
  setNumericStringParameter("ACTIONTIME","-1.0","Points in time (corresponding to ACTIONDT values), where certain actions are carried out. Order: [t_equilib; t_ktswitch; ...; t_act]",&statmech);
  // time step sizes corresponding to ACTIONTIME
  setNumericStringParameter("ACTIONDT","-1.0","Time step sizes corresponding to ACTIONTIME values.",&statmech);
  // index controlling the start of BC application (see ACTIONTIME)
  IntParameter("BCTIMEINDEX",-1,"Integer refers to the n-th entry of ACTIONTIME. States beginning of BC application. Starts counting at '1' !",&statmech);
  //Reading whether fixed seed for random numbers should be applied
  setStringToIntegralParameter<int>("FILAMENTPOLARITY","No","toggles filament polarity",yesnotuple,yesnovalue,&statmech);
  //Rise per monomer in the actin double helix according to Howard, p. 125
  setNumericStringParameter("RISEPERBINSPOT","0.00277","rise per monomer in the actin one-start helix",&statmech);
  //Rotation per monomer in the actin double helix according to Howard, p. 125
  DoubleParameter("ROTPERBINSPOT",-2.8999,"rotation per monomer in the actin double-helix",&statmech);
  //angular offset of the binding spot orientation (constant for each filament)
  DoubleParameter("BINSPOTOFFSET",0.0,"angular offset of the binding spot orientation (constant for each filament)",&statmech);
  //angle between binding spot orientation and the surface of the cone-shaped binding spot reactive volume
  DoubleParameter("PHIBINSPOT",0.524,"angle between binding spot orientation and the surface of the cone-shaped binding spot reactive volume",&statmech);
  //Reading double parameter for shear flow field
  DoubleParameter("SHEARAMPLITUDE",0.0,"Shear amplitude of flow in z-direction; note: not amplitude of displacement, but of shear strain!",&statmech);
  //Reading double parameter for viscosity of background fluid
  DoubleParameter("ETA",0.0,"viscosity",&statmech);
  //Reading double parameter for thermal energy in background fluid (temperature * Boltzmann constant)
  DoubleParameter("KT",0.0,"thermal energy",&statmech);
  //Reading double parameter for thermal energy in background fluid (temperature * Boltzmann constant)
  DoubleParameter("KTACT",0.0,"thermal energy for t>=STARTTIMEACT",&statmech);
  //Reading double parameter for crosslinker on-rate at the beginning
  DoubleParameter("K_ON_start",0.0,"crosslinker on-rate at the end",&statmech);
  //Reading double parameter for crosslinker on-rate at the end
  DoubleParameter("K_ON_end",0.0,"crosslinker on-rate at the end",&statmech);
  //Reading double parameter for crosslinker off-rate at the beginning
  DoubleParameter("K_OFF_start",0.0,"crosslinker off-rate at the beginning",&statmech);
  //Reading double parameter for crosslinker off-rate at the end
  DoubleParameter("K_OFF_end",0.0,"crosslinker off-rate at the end",&statmech);
  //Reading double parameter for crosslinker off-rate at the end
  DoubleParameter("K_ON_SELF",0.0,"crosslinker on-rate for crosslinkers with both bonds on same filament",&statmech);
  // chemical rate for contractile conformation change
  DoubleParameter("K_ACT_SHORT_start",0.0,"rate",&statmech);
  // chemical rate for contractile conformation change
  DoubleParameter("K_ACT_SHORT_end",0.0,"rate",&statmech);
  // chemical rate for extensional conformation change
  DoubleParameter("K_ACT_LONG_start",0.0,"rate",&statmech);
  // chemical rate for extensional conformation change
  DoubleParameter("K_ACT_LONG_end",0.0,"rate",&statmech);
  // active linker stroke distance
  DoubleParameter("STROKEDISTANCE",0.005,"active linker stroke distance",&statmech);
  // scaling factor for linker length changes
  DoubleParameter("LINKERSCALEFACTOR",0.0,"Scaling factor for active linker length changes. No effect on BEAM3CL elements!",&statmech);
  //displacement in the reaction coordinate used in Bell's euqations
  DoubleParameter("DELTABELLSEQ",0.0,"displacement in the reaction coordinate used in Bell's eqation (<0.0 -> catch bond, >0 -> std bond, 0 == no Bell's equation ",&statmech);
  // active linker fraction
  DoubleParameter("ACTIVELINKERFRACTION",0.0,"Fraction of linkers that show active behavior", &statmech);
  // cycle time
  DoubleParameter("ACTIVELINKERCYCLE",0.04,"duration of a work cycle of an active linker",&statmech);
  // time fraction during which no bonding is possible due to the linker being in its recovery state
  DoubleParameter("ACTIVERECOVERYFRACTION",0.95,"fraction of ACTIVELINKERCYCLE during which the linker recovers (i.e. is unbound)",&statmech);
  //number of overall crosslink molecules in the boundary volume
  IntParameter("NUMCROSSLINK",0,"number of crosslinkers for switching on- and off-rates; if molecule diffusion model is used: number of crosslink molecules",&statmech);
  //number of overall crosslink molecules in the boundary volume
  IntParameter("INITOCCUPIEDBINSPOTS",0,"binding spots occupied by (singly-bound) crosslinkers before the first time step",&statmech);
  //number of filaments used as substrate for motility assay setups
  IntParameter("NUMSUBSTRATEFIL",0,"Number of filaments used as substrate filaments",&statmech);
  //number by which the number of crosslinkers is reduced.
  IntParameter("REDUCECROSSLINKSBY",0,"number of crosslinker elements by which the overall number of crosslinker is reduced.",&statmech);
  IntParameter("BINSPOTINTERVAL",1,"determines every n-th binding spot available for crosslinking",&statmech);
  //Reading double parameter for crosslinker protein mean length
  DoubleParameter("R_LINK",0.0,"Mean distance between two nodes connected by a crosslinker",&statmech);
  //Absolute value of difference between maximal/minimal and mean cross linker length
  DoubleParameter("DeltaR_LINK",0.0,"Absolute value of difference between maximal/minimal and mean cross linker length",&statmech);
  // Three values representing the size of the periodic box in each spatial direction
  setNumericStringParameter("PERIODLENGTH","0.0 0.0 0.0", "Values representing the size of the periodic box in each spatial direction",&statmech);
  //angle between filament axes at crosslinked points with zero potential energy
  DoubleParameter("PHIZERO",0.0,"equilibrium angle between crosslinker axis and filament at each binding site",&statmech);
  //only angles in the range PHIZERO +/- PHIODEV are admitted at all for the angle PHI between filament axes at crosslinked points;
  //the default value for this parameter is 2*pi so that by default any value is admitted
  DoubleParameter("PHIZERODEV",6.28,"only angles in the range PHIZERO +/- PHIODEV",&statmech);
  //stiffness of orientation potential of crosslinkers
  DoubleParameter("CORIENT",0.0,"stiffness of orientation potential of crosslinkers",&statmech);
 //Young's modulus of crosslinkers
  DoubleParameter("ELINK",0.0,"Moment of inertia of area of crosslinkers",&statmech);
  //Moment of inertia of area of crosslinkers
  DoubleParameter("ILINK",0.0,"Moment of inertia of area of crosslinkers",&statmech);
  //Polar moment of inertia of area of crosslinkers
  DoubleParameter("IPLINK",0.0,"Polar moment of inertia of area of crosslinkers",&statmech);
  //Cross section of crosslinkers
  DoubleParameter("ALINK",0.0,"Cross section of crosslinkers",&statmech);
  //Torsional stiffness of crosslinkers
  DoubleParameter("KTOR1_LINK",0.0,"Torsional stiffness at binding spots between Truss linkers and filaments",&statmech);
  //Torsional stiffness of crosslinkers
  DoubleParameter("KTOR2_LINK",0.0,"Torsional stiffness between tangents of filaments",&statmech);
  //Makes filaments and crosslinkers be plotted by that factor thicker than they are acutally
  DoubleParameter("PlotFactorThick",0.0,"Makes filaments and crosslinkers be plotted by that factor thicker than they are acutally",&statmech);
  //Reading whether fixed seed for random numbers should be applied
  setStringToIntegralParameter<int>("CHECKORIENT","No","If chosen crosslinkers are set only after check of orientation of linked filaments", yesnotuple, yesnovalue, &statmech);
  //Gmsh Output switch
  setStringToIntegralParameter<int>("GMSHOUTPUT","No","If chosen gmsh output is generated.", yesnotuple,yesnovalue,&statmech);
  // toggling Gmsh Output for structure detection
  setStringToIntegralParameter<int>("GMSHNETSTRUCT","No","If chosen, special gmsh visualization for network structure types is generated.", yesnotuple,yesnovalue,&statmech);
  //Number of time steps between two special outputs written
  IntParameter("OUTPUTINTERVAL",1,"Number of time steps between two special outputs written",&statmech);
  //Number of interpolation points for higher order plotting of filaments, applicable only in case of Kirchchoff
  //type of beam elements reconstructed with hermite polynomials
  IntParameter("GMSHNUMINTPT",10,"Number of interpolation points for higher order plotting of filaments",&statmech);
  //Number of time steps between two gmsh outputs written
  IntParameter("GMSHOUTINTERVAL",100,"Number of time steps between two gmsh outputs written",&statmech);
  //Reading direction of oscillatory motion that DBC nodes are subjected to (we need this when using periodic BCs)
  IntParameter("DBCDISPDIR",0,"Global spatial direction of oscillatory motion by Dirichlet BCs",&statmech);
  //Reading time curve number for Dirichlet boundary conditions
  IntParameter("CURVENUMBER",0,"Specifies Time Curve number of imposed Dirichlet BCs",&statmech);
  //Reading time curve number for Neumann boundary conditions
  IntParameter("NBCCURVENUMBER",0,"Specifies Time Curve number of Neumann BCs",&statmech);
  // amplitude of Neumann boundary force
  DoubleParameter("NBCFORCEAMP",0.0,"constant creep force in NBCs",&statmech);
  // amplitude of Neumann boundary force
  IntParameter("NUMNBCNODES",0,"Number of nodes to which Neumann point forces are applied sequentially.",&statmech);
  //Reading number of elements that are taken into account when applying Dirichlet Conditions (useful to avoid redundant evaluation)
  // when Crosslink elements are added or the bead-spring-model is used
  IntParameter("NUM_EVAL_ELEMENTS",-1,"number of elements that are taken into account when applying Dirichlet Conditions",&statmech);
  // number of partitions along the edge length of the volume determining the resolution of the search grid
  IntParameter("SEARCHRES",1,"leads to the indexing of SEARCHRES^3 cubic volume partitions",&statmech);
  // number of histogram bins for post-analysis
  IntParameter("HISTOGRAMBINS",1,"number of bins for histograms showing the density-density-correlation-function",&statmech);
  // number of raster point along for ddcorr output boundary box shift -> n³ points in volume
  IntParameter("NUMRASTERPOINTS",3,"number of bins for histograms showing the density-density-correlation-function",&statmech);
  //time interval in which random numbers are constant
  DoubleParameter("TIMEINTCONSTRANDNUMB",-1.0,"Within this time interval the random numbers remain constant. -1.0 means no prescribed time interval.'",&statmech);
  //cutoff for random forces, which determines the maximal value
  DoubleParameter("MAXRANDFORCE",-1.0,"Any random force beyond MAXRANDFORCE*(standard dev.) will be omitted and redrawn. -1.0 means no bounds.'",&statmech);
}
