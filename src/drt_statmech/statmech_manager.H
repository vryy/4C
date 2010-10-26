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

class StatMechManager
{
public:

  //! \brief Standard Constructor
  explicit StatMechManager(ParameterList& params, DRT::Discretization& discret);

  //! \brief Destructor
  virtual ~StatMechManager() {};

  //! \brief write special output for statistical mechanics into files
  void StatMechOutput(ParameterList& params, const int ndim, const double& time,const int& istep, const double& dt, const Epetra_Vector& dis, const Epetra_Vector& fint);

  //! \brief saving all necessary class variables for restart
  void StatMechWriteRestart(IO::DiscretizationWriter& output);

  //! \brief reading all necessary information after a restart
  void StatMechReadRestart(IO::DiscretizationReader& reader);

  /*! \brief special update for statistical mechanics including setting and deleting crosslinkers in the discretization
   *  \param dt				(in)	current time step size
   *  \param disrow		(in)	current displacement vector
   *  \param stiff		(in)
   *  \param ndim			(in)
   */
  void StatMechUpdate(const int& istep, const double dt, Epetra_Vector& disrow, RCP<LINALG::SparseOperator>& stiff, int ndim);

  //! \brief initializing special output for statistical mechanics by looking for a suitable name of the outputfile and setting up an empty file with this name
  void StatMechInitOutput(const int ndim, const double& dt);

	//! \brief write output for Gmsh postprocessing
	void GmshOutput(const Epetra_Vector& disrow, const std::ostringstream& filename, const int& step);

	/*! \brief gmsh output in case of periodic boundary conditions
		 *  \param coord 						(in)
		 *  \param color 						(in)
		 *  \param gmshfilecontent  (in)
		 *  \param eleid 						(in)
		 *  \param crosslinkdiff		(in) determines the visualization mode (actual elements vs. crosslink molecule representations)
		 */
		void GmshOutputPeriodicBoundary(const LINALG::SerialDenseMatrix& coord, const double& color, std::stringstream& gmshfilecontent, int eleid, bool ignoreeleid=false);

  //! \brief set crosslinker if two nodes are close enough and probability check is passed
  void SetCrosslinkers(const double& dt, const Epetra_Map& noderowmap, const Epetra_Map& nodecolmap,
                       const std::map<int,LINALG::Matrix<3,1> >& currentpositions,const std::map<int,
                       LINALG::Matrix<3,1> >& currentrotations, Epetra_MultiVector& crosslinkerneighbours);

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

  /*!\brief check whether element is "broken" due to periodic boundary conditions
    \param coord  (in)  set of nodal coordinates
		\param cut    (in/out)  in case of a broken element, cut holds component and information about the location of the nodes of the element
		\param broken (in/out)  continuous==false, broken==true
		\date 03/10
  */
  void CheckForBrokenElement(LINALG::SerialDenseMatrix& coord, LINALG::SerialDenseMatrix& cut, bool *broken);

  /*!\brief return the current coordinates and DOF LIDs of all nodes of given element
		 \param element (in) 	   an element
		 \param dis		  (in)	   current displacement vector
		 \param coord	  (in/out) current node coordinates of element
		 \param lids		(in/out) nodal DOF LIDs of the element
  */
  void GetElementNodeCoords(DRT::Element* element, RCP<Epetra_Vector> dis, LINALG::SerialDenseMatrix& coord, vector<int>* lids=NULL);

  /*!\brief update locations of force sensors
     \param sensornodes   (in) vector holding the GIDs of the sensor nodes
     \param oscdir        (in) direction (DOF) of oscillation
   */
  void UpdateForceSensors(vector<int>& sensornodes, int oscdir);

  //! \brief checks orientation of crosslinker relative to linked filaments
  bool CheckOrientation(const LINALG::Matrix<3,1> direction, const Epetra_MultiVector& nodaltriadscol, const LINALG::Matrix<2,1>& LID);

  //! \brief containing statistical mechanics input parameters
  Teuchos::ParameterList statmechparams_;

  //! \brief return initial number of elements
  virtual int BasisElements(){return basiselements_;}
  //! \brief return initial number of nodes
  virtual int BasisNodes(){return basisnodes_;}
  //! \brief current simulation time
  double time_;

protected:

  // don't want = operator and cctor
  StatMechManager operator = (const StatMechManager& old);
  StatMechManager(const StatMechManager& old);

  //! \brief delete crosslinkers according to delcrosslinkercol (if entry for a node equals 1 its crosslinker is deleted)
  void DelCrosslinkers(const double& dt, const Epetra_Map& noderowmap, const Epetra_Map& nodecolmap);

  //! \brief delete crosslinkers according to delcrosslinkercol (if entry for a node equals 1 its crosslinker is deleted)
  void SearchAndDeleteCrosslinkers(const double& dt, const Epetra_Map& noderowmap, const Epetra_Map& nodecolmap,
																	 const std::map<int,LINALG::Matrix<3,1> >& currentpositions);

  //! \brief Searches and saves in variable  crosslinkerneighbours neighbours for each col node
  void SearchNeighbours(const std::map<int,LINALG::Matrix<3,1> > currentpositions, RCP<Epetra_MultiVector>& crosslinkerneighbours);

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

  //! \brief Shifts current position of nodes so that they comply with periodic boundary conditions
  void PeriodicBoundaryShift(Epetra_Vector& disrow, int ndim, const double &dt);

  //! \brief initializes vector with crosslinker positions
  void CrosslinkerMoleculeInit();

  //! \brief shifts crosslinker reprentations (spatial points) if they drift out of the boundary box in case of periodic BCs
	void CrosslinkerPeriodicBoundaryShift(Epetra_MultiVector& crosslinkerpositions);

  //! \brief generates a vector with a random permutation of the numbers between 0 and N - 1
  std::vector<int> Permutation(const int& N);

  /*! \brief writes element orientation and coordinates of fiber elements to an output file
			\param disrow     (in) displacement in row format
			\param filename1  (in) filename for the first output file (coordinates, element GIDs and filament numbers)
			\param filename2  (in) filename for the second output file (matrix containing crosslinker numbers between filaments)
	 */
	void StructPolymorphOutput(const Epetra_Vector& disrow, const std::ostringstream& filename1, const std::ostringstream& filename2);

	void DensityDensityCorrOutput(const std::ostringstream& filename);

	//! \brief check whether a crosslink connects two filaments or occupies two binding spots of the same element and therefore is to be visualized differently
	bool CheckForKinkedVisual(int eleid);

	//! \brief gmsh output of periodic boundary box
	void GmshOutputBoundaryBox(double color, const std::ostringstream *filename);

	//! \brief gmsh output for crosslink molecule diffusion
	void GmshOutputCrosslinkDiffusion(double color, const std::ostringstream *filename, const Epetra_Vector& disrow);

	/*! \brief prepare visualization vectors
	 * 	\param dis        (in) discplacement
	 */
	void GmshPrepareVisualization(const Epetra_Vector& dis);

	//! \brief GmshOutput for crosslinkers occupying two binding spots on the same filament
	void GmshKinkedVisual(const LINALG::SerialDenseMatrix& coord, const double& color, int eleid, std::stringstream& gmshfilecontent);

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

	//! \brief initialization
	bool isinit_;

  //! \brief switch between different on and off-rates
  bool konswitch_;

  //! \brief how often SearchNeighbours() has already been called
  int nsearch_;

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

  //! \brief number of elements which are part of the discretization in the beginning of this time step
  int currentelements_;

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

  //! \brief storing for each col map node the GIDs of nodes to which a crosslinker has been established
  RCP<Epetra_MultiVector>  crosslinkerpartner_;


  //! \brief storing for each col map node how many crosslinkers it has currently
  RCP<Epetra_Vector>  numcrosslinkerpartner_;

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

  //! SearchTree which allows accelerated search for neighbouring nodes in case of dynamic crosslinkers
  Teuchos::RCP<GEO::SearchTree>  octTree_;

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
	// \brief a vector just for the sake of visual quality
	RCP<Epetra_MultiVector> visualizepositions_;
	// \brief a vector holding information on whether or not a crosslink molecule searches for neighbours
	RCP<Epetra_Vector> searchforneighbours_;
	// \brief a vector holding indices used in DensityDensityCorrOutput()
	std::vector<std::vector<int> > startindex_;
}; // class StatMechManager

#endif /*STATMECH_MANAGER_H*/
#endif  // #ifdef CCADISCRET
