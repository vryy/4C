/*!----------------------------------------------------------------------
\file statmech_manager.H

<pre>
Maintainer: Christian Cyron
            cyron@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15234
</pre>

*----------------------------------------------------------------------*/
#ifdef CCADISCRET
#ifndef STATMECH_MANAGER_H
#define STATMECH_MANAGER_H

#include <Teuchos_StandardParameterEntryValidators.hpp>
#include "../drt_geometry/searchtree.H"
#include "../drt_geometry/intersection_service.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_beamcontact/beam3contact_manager.H"
#include "../drt_io/io.H"
#include "../drt_io/io_control.H"

//including random number library of blitz for statistical forces
#include <random/uniform.h>
#include <random/normal.h>

//forward declaration of class Beam3 used in CheckOrientation(...)
namespace DRT
{
  namespace ELEMENTS
  {
    class Beam3;
  }
}


/*!
\brief management and auxiliary functions for statistical mechanics
\author cyron (cyron@lnm.mw.tum.de)
*/
using namespace CONTACT;

class StatMechManager
{
public:

  //! \brief Standard Constructor
  explicit StatMechManager(ParameterList& params, DRT::Discretization& discret);

  //! \brief Destructor
  virtual ~StatMechManager() {};

  //! \brief write special output for statistical mechanics into files
  void Output(ParameterList& params, const int ndim, const double& time,const int& istep, const double& dt, const Epetra_Vector& dis, const Epetra_Vector& fint, RCP<CONTACT::Beam3cmanager> beamcmanager=Teuchos::null);

  //! \brief saving all necessary class variables for restart
  void WriteRestart(IO::DiscretizationWriter& output);

  //! \brief write restart information for fully overlapping Epetra_Multivector with name "name"
  void WriteRestartRedundantMultivector(IO::DiscretizationWriter& output,const string name, RCP<Epetra_MultiVector> multivector);

  //! \brief read restart information for fully overlapping Epetra_Multivector with name "name"
  void ReadRestartRedundantMultivector(IO::DiscretizationReader& reader,const string name, RCP<Epetra_MultiVector> multivector);

  //! \brief reading all necessary information after a restart
  void ReadRestart(IO::DiscretizationReader& reader);

  //! \brief saves all relevant variables *_ as *conv_ to allow for returning to the beginning of a time step
  void WriteConv();

  //! \brief resets all relevant variables *_ to *conv_ in order to return to the state at the beginning of the time step
  void RestoreConv(RCP<LINALG::SparseOperator>& stiff);

  /*! \brief special update for statistical mechanics including setting and deleting crosslinkers in the discretization
   *  \param dt				(in)	current time step size
   *  \param disrow		(in)	current displacement vector
   *  \param stiff		(in)
   *  \param ndim			(in)
   */
  void Update(const int& istep, const double dt, Epetra_Vector& disrow, RCP<LINALG::SparseOperator>& stiff, int ndim);

  //! \brief initializing special output for statistical mechanics by looking for a suitable name of the outputfile and setting up an empty file with this name
  void InitOutput(const int ndim, const double& dt);

  //! \brief seed all random generators of this object properly
  void SeedRandomGenerators(const int seedparameter);

	//! \brief write output for Gmsh postprocessing
	void GmshOutput(const Epetra_Vector& disrow, const std::ostringstream& filename, const int& step, RCP<CONTACT::Beam3cmanager> beamcmanager=Teuchos::null);

	/*! \brief gmsh output in case of periodic boundary conditions
		 *  \param coord 						(in)
		 *  \param color 						(in)
		 *  \param gmshfilecontent  (in)
		 *  \param eleid 						(in)
		 *  \param crosslinkdiff		(in) determines the visualization mode (actual elements vs. crosslink molecule representations)
		 */
	void GmshOutputPeriodicBoundary(const LINALG::SerialDenseMatrix& coord, const double& color, std::stringstream& gmshfilecontent, int eleid, bool ignoreeleid=false);


  //! \brief Searches and saves in variable  crosslinkerneighbours neighbours for each col node
  void SearchAndSetCrosslinkers(const int& istep,
																const double& dt,
																const Epetra_Map& noderowmap,
																const Epetra_Map& nodecolmap,
																const std::map<int,LINALG::Matrix<3,1> >& currentpositions,
																const std::map<int,LINALG::Matrix<3,1> >& currentrotations);

  //! \brief generate gaussian randomnumbers with mean "meanvalue" and standarddeviation "standarddeviation" for parallel use
  void GenerateGaussianRandomNumbers(RCP<Epetra_MultiVector> randomnumbers,const double meanvalue, const double standarddeviation);

  //! \brief Initializes beam3 elements in case of perioddic boundary conditions
  void PeriodicBoundaryBeam3Init(DRT::Element* element);

  //! \brief Initializes beam3ii elements in case of perioddic boundary conditions
  void PeriodicBoundaryBeam3iiInit(DRT::Element* element);

  //! \brief Initializes truss3 elements in case of perioddic boundary conditions
  void PeriodicBoundaryTruss3Init(DRT::Element* element);

  //! \brief Initializes trusslm elements in case of perioddic boundary conditions
  void PeriodicBoundaryTrussLmInit(DRT::Element* element);

  //! \brief Shifts current position of nodes so that they comply with periodic boundary conditions
  void PeriodicBoundaryShift(Epetra_Vector& disrow, int ndim, const double &dt);

  //! \brief checks orientation of crosslinker relative to linked filaments
  bool CheckOrientation(const LINALG::Matrix<3,1> direction, const Epetra_MultiVector& nodaltriadscol, const LINALG::Matrix<2,1>& LID, RCP<double> phifil = Teuchos::null);

	//! \brief counts the number of possible crosslinking spots and does an Orientation Check
	void OrientationCorrelation(const Epetra_Vector& disrow, const int &istep);

  /*! \brief adjusted DBC evaluation for periodic case
   *  \param params					 a parameter list*/
  void EvaluateDirichletPeriodic(ParameterList& params,
																RCP<Epetra_Vector> disn,
																RCP<Epetra_Vector> dirichtoggle,
																RCP<Epetra_Vector> invtoggle);

  /*! \brief Determine Dirichlet condition at given time and apply its values to a system vector
   *  \param nodeids				 vector containing DBC Node IDs
   *  \param onoff					 vector containing onoff values for DOFs*/
  void DoDirichletConditionPeriodic(vector<int>*								nodeids,
																		vector<int>*								onoff,
																		RCP<Epetra_Vector> disn,
																		RCP<Epetra_Vector> dirichtoggle,
																		RCP<Epetra_Vector> invtoggle,
																		Epetra_Vector& deltadbc);

	//! brief computes the internal energy of the network
	void ComputeInternalEnergy(const RCP<Epetra_Vector> dis, double& energy, const double& dt, const std::ostringstream& filename);

  //! \brief containing statistical mechanics input parameters
  Teuchos::ParameterList statmechparams_;

  //! \brief current simulation time
  double time_;

  //! \brief how often a time step had to be restarted with new random numbers due to unconverged solution of nonlinear system of equations since simulation start
  int unconvergedsteps_;

protected:

  // don't want = operator and cctor
  StatMechManager operator = (const StatMechManager& old);
  StatMechManager(const StatMechManager& old);

  /*!\brief check whether element is "broken" due to periodic boundary conditions
    \param coord  (in)  set of nodal coordinates
		\param cut    (in/out)  in case of a broken element, cut holds component and information about the location of the nodes of the element
		\param broken (in/out)  continuous==false, broken==true*/
  void CheckForBrokenElement(LINALG::SerialDenseMatrix& coord, LINALG::SerialDenseMatrix& cut, bool *broken);

  /*!\brief return the current coordinates and DOF LIDs of all nodes of given element
		 \param element (in) 	   an element
		 \param dis		  (in)	   current displacement vector
		 \param coord	  (in/out) current node coordinates of element
		 \param lids		(in/out) nodal DOF LIDs of the element*/
  void GetElementNodeCoords(DRT::Element* element, RCP<Epetra_Vector> dis, LINALG::SerialDenseMatrix& coord, vector<int>* lids=NULL);
  /*!\brief update locations of force sensors
     \param sensornodes   (in) vector holding the GIDs of the sensor nodes
     \param oscdir        (in) direction (DOF) of oscillation*/
  void UpdateForceSensors(vector<int>& sensornodes, int oscdir);

  //! \brief delete crosslinkers according to delcrosslinkercol (if entry for a node equals 1 its crosslinker is deleted)
  void SearchAndDeleteCrosslinkers(const double& dt, const Epetra_Map& noderowmap, const Epetra_Map& nodecolmap,
																	 const std::map<int,LINALG::Matrix<3,1> >& currentpositions);

  //! \brief force based unbinding of crosslinkers
  void ForcedCrosslinkerUnbinding(const double& dt, Epetra_Vector* unbindbyforce, Epetra_Vector* delcrosselement);

  /*! \brief reduce number of currently existing crosslinker elements by numtoreduce
   *  \param (in) numtoreduce number of crosslinkers by which overall number of crosslinkers is reduced*/
  void ReduceNumOfCrosslinkersBy(const int numtoreduce);

  /*! \brief Search neighbour filament nodes to a crosslink molecule (when crosslinker diffusion is simulated)
   *  \param currentpositions 		(in) column map node positions
   *  \param nodeinpartition  		(in) vector indicating the partition of the nodes
   *  \param numbond							(in) vector holding the number of established bonds for each crosslink molecule
   *  \param crosslinkerpositions	(in) vector holding the positions of the crosslink molecules
   *  \param crosslinkpartitions	(in) vector holding the information to which partial volume the crosslink molecules belong
   *  \param neighboursLID		(in/out) vector holding the LIDs of the nodes close to crosslink molecules*/
  void DetectNeighbourNodes(const std::map<int,LINALG::Matrix<3,1> >& currentpositions,
													  std::vector<std::vector<std::vector<int> > >* nodeinpartition,
													  Epetra_Vector& numbond,
													  Epetra_MultiVector& crosslinkerpositions,
													  Epetra_MultiVector& crosslinkpartitions,
														RCP<Epetra_MultiVector>& neighboursLID);

  //! \brief initializes vector with crosslinker positions
  void CrosslinkerMoleculeInit();

  //! \brief shifts crosslinker reprentations (spatial points) if they drift out of the boundary box in case of periodic BCs
	void CrosslinkerPeriodicBoundaryShift(Epetra_MultiVector& crosslinkerpositions);

  //! \brief generates a vector with a random permutation of the numbers between 0 and N - 1
  std::vector<int> Permutation(const int& N);

  /*! \brief writes element orientation and coordinates of fiber elements to an output file
			\param discol     (in) displacement in column map format
			\param filename  (in) filename for the first output file (coordinates, element GIDs and filament numbers)
	 */
	void FilamentOrientations(const Epetra_Vector& discol, std::vector<LINALG::Matrix<3,1> >* normedvectors, const std::ostringstream& filename, bool fileoutput);

	//! \brief output for density-density-correlation output
	void DDCorrOutput(const Epetra_Vector& disrow, const std::ostringstream& filename, const int& istep, const double& dt);

	//! \brief obtain Density-Density-Correlation function (radial distribution function) for given crosslinkers crosslinksperbinrow
	void DDCorrFunction(Epetra_MultiVector& crosslinksperbinrow, Epetra_MultiVector& crosslinksperbinrotrow, LINALG::Matrix<3,1>* centershift);

	//! \brief obtain distribution for spherical coordinates (all vectors in and out)
	void SphericalCoordsDistribution(const Epetra_Vector& disrow, Epetra_Vector& phibinsrow, Epetra_Vector& thetabinsrow, Epetra_Vector& costhetabinsrow);

	//! \brief Selects raster point with the smallest average distance to all crosslinker elements and makes it the new center of the boundary box
	void DDCorrShift(LINALG::Matrix<3,1>* boxcenter, LINALG::Matrix<3,1>* centershift, std::vector<int>* crosslinkerentries);

	//! \brief obtain radial crosslinker density distribution
	void RadialDensityDistribution(Epetra_Vector& radialdistancesrow, LINALG::Matrix<3,1>& centershift);

	/*! \brief determine the current network structure
	 * 	\param disrow  (in) displacement in row format
	 * 	\param centershift (in) shift vector between old and new boundary box center
	 * 	\param crosslinkerentries (in) crosslinkermap_ lids of crosslinker elements
	 * 	\param istep (in) i-th time step
	 * 	\param filename (in) output file name
	 * 	\param filorientoutput (in) determine whether or not to output filament element directional vectors
	 */
	void DDCorrCurrentStructure(const Epetra_Vector& disrow,
															LINALG::Matrix<3,1>* cog,
															LINALG::Matrix<3,1>* centershift,
															std::vector<int>* crosslinkerentries,
															const int& istep,
															const std::ostringstream& filename,
															bool filorientoutput=true);

	//! \brief iterate vector in order to obtain close approximation to actual (structure) vector
	void DDCorrIterateVector(const Epetra_Vector& discol, LINALG::Matrix<3,1>* vectorj, const int& maxiterations);

	//! \brief simply returns the number of free, one-bonded, and two-bonded crosslink molecules
	void CrosslinkCount(const std::ostringstream& filename);

	//! \brief returns polymer length per unit volume with respect to radial distance to COG of the network
	void ComputeLocalMeshSize(const Epetra_Vector& disrow, LINALG::Matrix<3,1>& centershift, const int &istep);

	//! \brief check whether a crosslink connects two filaments or occupies two binding spots of the same element and therefore is to be visualized differently
	bool CheckForKinkedVisual(int eleid);

	//! \brief gmsh output of (periodic boundary) box around "boxcenter" with edgelength "length"
	void GmshOutputBox(double color, LINALG::Matrix<3,1>* boxcenter, double length, const std::ostringstream *filename);

	//! \brief gmsh output for crosslink molecule diffusion
	void GmshOutputCrosslinkDiffusion(double color, const std::ostringstream *filename, const Epetra_Vector& disrow);

	/*! \brief prepare visualization vectors
	 * 	\param dis        (in) discplacement
	 */
	void GmshPrepareVisualization(const Epetra_Vector& dis);

	//! \brief GmshOutput for crosslinkers occupying two binding spots on the same filament
	void GmshKinkedVisual(const LINALG::SerialDenseMatrix& coord, const double& color, int eleid, std::stringstream& gmshfilecontent);

	//! \brief Gmsh wedge output for filament and crosslinkers as solid lines using prisms
	void GmshWedge(const int& n, const Epetra_SerialDenseMatrix& coord, DRT::Element* thisele, std::stringstream& gmshfilecontent, const double color, bool ignoreeleid=false, bool drawsphere=true);

	/*! \brief Gmsh output of volume boundaries of detected structure type (detection in DDCorrCurrStructure())
	 *  \param n               (in) number of edges for approximation of the cluster sphere or the bundle cylinder
	 *  \param gmshfilecontent (in) where output is written to
	 *  \param color           (in) output color
	 */
	void GmshNetworkStructVolume(const int& n, std::stringstream& gmshfilecontent, const double color);

	/*! \brief gmsh output for test volumes in presence of periodic boundary conditions
	 *  \param coord          (in) coordinates of the end points of the given line
	 *  \param numsections    (in) number of sections by which the given line is divided
	 *  \param gmshfilcontent (in) where output is written to
	 *  \param color          (in) output color
	 */
	void GmshNetworkStructVolumePeriodic(const Epetra_SerialDenseMatrix& coord, const int numsections, std::stringstream& gmshfilecontent,const double color);


	/*! \brief Handles Crosslinker diffusion
	 *  \param dis             (in) latest displacement
	 *  \param disi						 (in) latest displacement increment in row map format
	 *  \param mean						 (in) mean of the normal distribution
	 *  \param standarddev     (in) standard deviation of the normal distribution
	 *  \param init						 (in) switch between methods
	 */
	void CrosslinkerDiffusion(const Epetra_Vector& dis,
														double mean,
														double standarddev,
														const double &dt);

	/*! \brief update positions of crosslink molecules if they constitute a crosslinker element
	 *  \param currentpositions			(in) column map node positions
	 *  \param LID									(in) column map node LIDs
	 *  \param crosslinkernumber		(in) position in crosslink molecule vector*/
	void CrosslinkerIntermediateUpdate(const std::map<int,LINALG::Matrix<3,1> >& currentpositions, const LINALG::SerialDenseMatrix& LID,
																		 int crosslinkernumber, bool coupledmovement=true);

	/*! \brief map crosslink molecules and filament nodes to the according volume partition and detect nodes in proximity of the molecules
	 * 	\param currentpositions			(in) 			column map node positions
	 * 	\param neighbourslid				(in/out)	vector holding the node lids within search proximity of crosslink molecules*/
	void PartitioningAndSearch(const std::map<int,LINALG::Matrix<3,1> >& currentpositions, RCP<Epetra_MultiVector>& neighbourslid);

  //! \brief point of time at which writing special statistical output is started
  double starttimeoutput_;

  //! \brief end to end distance of a single filament at the beginning of statistical output
  double endtoendref_;

  //! \brief time step number at the beginning of statistical output
  int istart_;

  //! \brief position of first node in last time step
  LINALG::Matrix<3,1> beginold_;

  //! \brief position of last nodein last time step
  LINALG::Matrix<3,1> endold_;

  //! \brief number of nodes which are part of the discretization already before starting adding crosslinkers
  int basisnodes_;

  //! \brief number of elements which are part of the discretization already before starting adding crosslinkers
  int basiselements_;

  //! \brief number of file into which output is written (only important for single filament dynamics)
  int outputfilenumber_;

  //! \brief generates Gaussian random numbers with zero mean and unit standard deviation
  ranlib::Normal<double> normalgen_;

  //! \brief generates uniformly distributed random numbers in interval [0;1]
  ranlib::UniformClosed<double> uniformclosedgen_;

  //! \brief generates uniformly distributed random numbers in interval [0;1[
  ranlib::UniformClosedOpen<double> uniformclosedopengen_;

  //! \brief reference to current discretization
  DRT::Discretization& discret_;

  //! \brief storing for each node column map LID number of filament to which it belongs
  RCP<Epetra_Vector> filamentnumber_;

  //! \brief storing for degree of freedom whether force has to be measured
  RCP<Epetra_Vector> forcesensor_;

  //! \store accumulated disp
  double sumsquareincpar_;
  double sumsquareincort_;
  double sumrotmiddle_;
  double sumsquareincmid_;
  double sumsquareincrot_;

  //store displacement of midpoint compared to initial state as a 3-vector
  LINALG::Matrix<3,1>sumdispmiddle_;

  //! brief SearchTree which allows accelerated search for neighbouring nodes in case of dynamic crosslinkers
  Teuchos::RCP<GEO::SearchTree>  octTree_;

  //! brief rotation matrix for transformation from global to local coordinates
  RCP<LINALG::SerialDenseMatrix> trafo_;

  // \brief column Epetra_Map handling crosslink molecules
  RCP<Epetra_Map> crosslinkermap_;
  // \brief row Epetra_Map handling crosslink molecules
  RCP<Epetra_Map> transfermap_;
  // \brief an Epetra_Map handling processor-specific data of the density-density-correlation-function
  RCP<Epetra_Map> ddcorrrowmap_;
  // \brief an Epetra_Map handling cumulated data of the density-density-correlation-function from all processors
  RCP<Epetra_Map> ddcorrcolmap_;
	// \brief a vector memorizing the bonding status of the crosslink molecules
	RCP<Epetra_MultiVector>		crosslinkerbond_;
	// \brief a vector containing the current postion of the crosslink molecules
	RCP<Epetra_MultiVector>		crosslinkerpositions_;
	// \brief a vector indicating the number of attached crosslinkers to a node
	RCP<Epetra_Vector>	numcrossnodes_;
	// \brief a vector indicating the number of attached nodes to a crosslinker
	RCP<Epetra_Vector>	numbond_;
	// \brief a vector indicating crosslinker elements that cover two binding spots on the same filament
	RCP<Epetra_Vector>	crosslinkonsamefilament_;
	// \brief a vector mapping crosslink molecule IDs to crosslinker element IDs
	RCP<Epetra_Vector> crosslink2element_;
	// \brief a vector just for the sake of visual quality
	RCP<Epetra_MultiVector> visualizepositions_;
	// \brief a vector holding information on whether or not a crosslink molecule searches for neighbours (switch for active/passive crosslinker)
	RCP<Epetra_Vector> searchforneighbours_;
	// \brief a vector holding indices used in DensityDensityCorrOutput()
	RCP<std::vector<double> > startindex_;

	//the following variables are necessary only for repeating a time step with different random numbers after its having failed one time

	// \brief elements deleted at the beginning when this time step was conducted the last time are packed into variable deletedelements_
	std::vector<std::vector<char> > deletedelements_;

	// \brief GIDs of elements added at the beginning when this time step was conducted the last time
	std::vector<int> addedelements_;

	//variables to save certain class variables *_ in *conv_ in order to allow for repeating a time step with completely new random numbers and actions
	RCP<Epetra_MultiVector> crosslinkerbondconv_;
  RCP<Epetra_MultiVector> crosslinkerpositionsconv_;
  RCP<Epetra_Vector> numcrossnodesconv_;
  RCP<Epetra_Vector> numbondconv_;
  RCP<Epetra_Vector> crosslinkonsamefilamentconv_;
	RCP<Epetra_Vector> crosslink2elementconv_;
  RCP<Epetra_Vector> searchforneighboursconv_;

  // \brief structure type
  int structuretype_;
  // \brief store test volume visualization data for Gmsh Output
  std::vector<LINALG::Matrix<3,1> > testvolumepos_;
  // \brief characteristic length (radii, thickness) of the detected structure
  double characlength_;
  // \brief center of gravity of the detected structure
  LINALG::Matrix<3,1> cog_;

  // glue initial set of broken elements at the beginning of oscillation to the upper and lower plate
  std::vector<int> oscillnodes_;
  std::vector<int> fixednodes_;
  std::vector<int> freenodes_;
  bool initialset_;

}; // class StatMechManager

#endif /*STATMECH_MANAGER_H*/
#endif  // #ifdef CCADISCRET
