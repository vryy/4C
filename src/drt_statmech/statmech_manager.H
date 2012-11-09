/*!----------------------------------------------------------------------
\file statmech_manager.H

<pre>
Maintainer: Kei Müller
            mueller@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15276
</pre>

*----------------------------------------------------------------------*/
#ifndef STATMECH_MANAGER_H
#define STATMECH_MANAGER_H

#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterListAcceptorDefaultBase.hpp>
#include <Teuchos_StandardParameterEntryValidators.hpp>
#include "Teuchos_RefCountPtr.hpp"
#include "Epetra_SerialDenseMatrix.h"
#include "Epetra_SerialDenseVector.h"
#include "Epetra_Vector.h"
#include "Epetra_MultiVector.h"
#include "Epetra_Map.h"

// include random number libraries of BOOST for statistical forces
#include <boost/random.hpp>
#include <boost/random/mersenne_twister.hpp>

#include "../linalg/linalg_fixedsizematrix.H"

/*!
\brief management and auxiliary functions for statistical mechanics
\author cyron (cyron@lnm.mw.tum.de)
*/

// forward declarations
namespace DRT
{
  class Condition;
  class Element;
  class Discretization;
}
namespace CONTACT
{
  class Beam3cmanager;
}

namespace IO
{
  class InputControl;
  class OutputControl;
  class DiscretizationWriter;
  class DiscretizationReader;
}

namespace LINALG
{
  class SparseOperator;
  class SparseMatrix;
  class SerialDenseMatrix;
  class MapExtractor;
}

class Beam3ContactOctree;

// random generator type
typedef boost::mt19937 randnumgen;

namespace STATMECH
{

class StatMechManager
{
public:
  //! \brief Standard Constructor
  explicit StatMechManager(Teuchos::RCP<DRT::Discretization> discret);

  //! \brief Destructor
  virtual ~StatMechManager() {};

  /*! \brief write special output for statistical mechanics into files
   *  \param params         (in) input parameter list
   *  \param ndim           (in) spatial dimensions
   *  \param time           (in) total physical time
   *  \param istep          (in) time step
   *  \param dt             (in) stime step size
   *  \param dis            (in) displacement vector
   *  \param fint           (in) internal force vector
   *  \param beamcmanager   (in) beam contact manager
   *  \param printscreen    (in) toggle output to console / log file
   */
  void Output(const int                            ndim,
              const double&                        time,
              const int&                           istep,
              const double&                        dt,
              const Epetra_Vector&                 dis,
              const Epetra_Vector&                 fint,
              Teuchos::RCP<CONTACT::Beam3cmanager> beamcmanager=Teuchos::null,
              bool                                 printscreen=true);

  //! \brief saving all necessary class variables for restart
  void WriteRestart(Teuchos::RCP<IO::DiscretizationWriter> output, double& dt);

  //! \brief write restart information for fully overlapping Epetra_Multivector with name "name"
  void WriteRestartRedundantMultivector(Teuchos::RCP<IO::DiscretizationWriter> output,const std::string name, Teuchos::RCP<Epetra_MultiVector> multivector);

  //! \brief read restart information for fully overlapping Epetra_Multivector with name "name"
  void ReadRestartRedundantMultivector(IO::DiscretizationReader& reader,const std::string name, Teuchos::RCP<Epetra_MultiVector> multivector);

  //! \brief reading all necessary information after a restart
  void ReadRestart(IO::DiscretizationReader& reader, double& dt);

  //! \brief saves all relevant variables *_ as *conv_ to allow for returning to the beginning of a time step
  void WriteConv();

  //! \brief resets all relevant variables *_ to *conv_ in order to return to the state at the beginning of the time step
  void RestoreConv(Teuchos::RCP<LINALG::SparseOperator>& stiff, Teuchos::RCP<CONTACT::Beam3cmanager> beamcmanager);

  /*! \brief special update for statistical mechanics including setting and deleting crosslinkers in the discretization
   *  \param istep    (in)  i-th time step
   *  \param timen    (in)  target time
   *  \param dt       (in)	current time step size
   *  \param disrow   (in)	current displacement vector
   *  \param stiff    (in)
   *  \param ndim     (in)
   */
  void Update(const int& istep,
              const double& timen,
              const double& dt,
              Epetra_Vector& disrow,
              Teuchos::RCP<LINALG::SparseOperator>& stiff,
              int& ndim,
              Teuchos::RCP<CONTACT::Beam3cmanager> beamcmanager=Teuchos::null,
              bool rebuildoctree=false,
              bool printscreen=true);

  /*! \brief update time step size in time integration scheme
   *  \param dt            (in) time step size
   *  \param timeconverged (in) time of last converged timestep   */
  void UpdateTimeAndStepSize(double& dt, double& timeconverged);

  //! \brief initializing special output for statistical mechanics by looking for a suitable name of the outputfile and setting up an empty file with this name
  void InitOutput(const int& ndim, const double& dt);

  //! \brief Get current binding spot positions
  void GetNodePositions(Epetra_Vector& discol, std::map<int,LINALG::Matrix<3,1> >& currentpositions, std::map<int, LINALG::Matrix<3, 1> >& currentrotations, bool positionsonly=false);

  //! \brief sets crosslinkers  whereever possible before the first time step
  void SetInitialCrosslinkers(Teuchos::RCP<CONTACT::Beam3cmanager> beamcmanager=Teuchos::null);

  //! \brief seed all random generators of this object properly
  void SeedRandomGenerators(const int seedparameter);

  //! \brief write output for Gmsh postprocessing
  void GmshOutput(const Epetra_Vector& disrow, const std::ostringstream& filename, const int& step, Teuchos::RCP<CONTACT::Beam3cmanager> beamcmanager=Teuchos::null);

  /*! \brief gmsh output in case of periodic boundary conditions
     *  \param coord            (in) coordinates of element's nodes
     *  \param color            (in) as is
     *  \param gmshfilecontent  (in) stream to be written in output file
     *  \param eleid            (in) element gid
     *  \param ignoreeleid      (in) signals whether or not we have an element to visualize
     */
  void GmshOutputPeriodicBoundary(const LINALG::SerialDenseMatrix& coord, const double& color, std::stringstream& gmshfilecontent, int eleid, bool ignoreeleid=false);

  //! \brief update nodal triads
  void GetBindingSpotTriads(Epetra_MultiVector* nodaltriadscol);

  //! \brief generate gaussian randomnumbers with mean "meanvalue" and standarddeviation "standarddeviation" for parallel use
  void GenerateGaussianRandomNumbers(Teuchos::RCP<Epetra_MultiVector> randomnumbers,const double meanvalue, const double standarddeviation);

  //! \brief Initializes beam3 elements in case of periodic boundary conditions
  void PeriodicBoundaryBeam3Init(DRT::Element* element);

  //! \brief Initializes beam3ii elements in case of periodic boundary conditions
  void PeriodicBoundaryBeam3iiInit(DRT::Element* element);

  //! \brief Initializes truss3 elements in case of periodic boundary conditions
  void PeriodicBoundaryTruss3Init(DRT::Element* element);

  /*! \brief Shifts current position of nodes so that they comply with periodic boundary conditions
   * \param disrow  (in)  row map displacement vector
   * \param ndim    (in)  number of spatial dimensions
   * \param timen   (in)  target time
   * \param dt      (in)  time step size
   */
  void PeriodicBoundaryShift(Epetra_Vector& disrow, int ndim, const double &timen, const double &dt);

  //! \brief checks orientation of crosslinker relative to linked filaments
  bool CheckOrientation(const LINALG::Matrix<3,1> direction, const Epetra_MultiVector& nodaltriadscol, const Epetra_SerialDenseMatrix& LID, Teuchos::RCP<double> phifil = Teuchos::null);

	//! \brief counts the number of possible crosslinking spots and does an Orientation Check
	void OrientationCorrelation(const Epetra_Vector& disrow, const int &istep);

  //! \brief create a new crosslinker element and add it to the discretization
  void AddNewCrosslinkerElement(const int& crossgid,
                                int* globalnodeids,
                                const std::vector<double>& xrefe,
                                const std::vector<double>& rotrefe,
                                DRT::Discretization& mydiscret,
                                bool addinitlinks=false);

  //! \brief remove crosslinker element from a chosen discretization
  void RemoveCrosslinkerElements(DRT::Discretization& mydiscret,
                                 Epetra_Vector& delcrosselement,
                                 std::vector<std::vector<char> >* deletedelements);

  /*! \brief adjusted DBC evaluation for StatMech
   *  \param params          (in)     a parameter list
   *  \param disn            (in)     displacement vector
   *  \param dbcmapextractor (in/out) map extractor for Dirichlet DOFs (new method)*/
  void EvaluateDirichletStatMech(Teuchos::ParameterList&                     params,
                                 Teuchos::RCP<Epetra_Vector>        dis,
                                 Teuchos::RCP<LINALG::MapExtractor> dbcmapextractor=Teuchos::null);

  /*! \brief adjusted DBC evaluation for periodic case
   *  \param params          (in)     a parameter list
   *  \param disn            (in)     displacement vector
   *  \param dbcmapextractor (in/out) map extractor for Dirichlet DOFs (new method)*/
  void EvaluateDirichletPeriodic(Teuchos::ParameterList&            params,
                                 Teuchos::RCP<Epetra_Vector>        dis,
                                 Teuchos::RCP<LINALG::MapExtractor> dbcmapextractor=Teuchos::null);

  /* \brief signals whether or not to start DBC application according to given time
   * \param params  (in)  Paramter list*/
  bool DBCStart(Teuchos::ParameterList& params);

  //! brief computes the internal energy of the network
  void ComputeInternalEnergy(const Teuchos::RCP<Epetra_Vector> dis, double& energy, const double& dt, const std::ostringstream& filename, bool fillzeros=true, bool writefile=true);

  //! \brief return number of initial elements
  const int& NumBasisElements() const {return basiselements_;}

  //! \brief return vector holding the dimensions of the periodic boundary box
  const Teuchos::RCP<std::vector<double> > GetPeriodLength() const {return periodlength_;}

  //! \brief return statistical mechanics parameter list
  const Teuchos::ParameterList& GetStatMechParams() const {return statmechparams_;}

  //! \brief return number of unconverged steps since beginning of the simulation
  const int& NumberOfUnconvergedSteps() const {return unconvergedsteps_;}

  //! \brief return number of unconverged steps since beginning of the simulation
  void UpdateNumberOfUnconvergedSteps(bool flag = true);

  //! \brief pass statistical mechanics parameters to parameter list
  void AddStatMechParamsTo(Teuchos::ParameterList& params, Teuchos::RCP<Epetra_MultiVector> randomnumbers = Teuchos::null);

protected:

  // don't want = operator and cctor
  StatMechManager operator = (const StatMechManager& old);
  StatMechManager(const StatMechManager& old);

  //! \brief Initialize periodlength_ vector and search resolution
  void InitializeStatMechValues();

  /*! \brief Searches and saves in variable  crosslinkerneighbours neighbours for each col node
   *  \param istep            (in) time step
   *  \param timen            (in) target time
   *  \param dt               (in) step size
   *  \param noderowmap       (in) Row Map
   *  \param nodecolmap       (in) Column Map
   *  \param currentpositions (in) positions of the nodes in column map format (fully overlapping)
   *  \param currentrotations (in) rotations of the nodes in column map format (fully overlapping)
   *  \param beamcmanager     (in) beam contact manager
   */
  void SearchAndSetCrosslinkers(const int&                                istep,
                                const double&                             timen,
                                const double&                             dt,
                                const Epetra_Map&                         noderowmap,
                                const Epetra_Map&                         nodecolmap,
                                const std::map<int,LINALG::Matrix<3,1> >& currentpositions,
                                const std::map<int,LINALG::Matrix<3,1> >& currentrotations,
                                Teuchos::RCP<CONTACT::Beam3cmanager>      beamcmanager,
                                bool                                      printscreen=true);

  /*! \brief delete crosslinkers according to delcrosslinkercol (if entry for a node equals 1 its crosslinker is deleted)
   *  \param timen            (in) target time
   *  \param dt               (in) step size
   *  \param noderowmap       (in) Row Map
   *  \param nodecolmap       (in) Column Map
   *  \param currentpositions (in) positions of the nodes in column map format (fully overlapping)
   *  \param discol           (in) node discplacement in column map format (fully overlapping)
   *  \param beamcmanager     (in) beam contact manager
   */
  void SearchAndDeleteCrosslinkers(const double&                             timen,
                                   const double&                             dt,
                                   const Epetra_Map&                         noderowmap,
                                   const Epetra_Map&                         nodecolmap,
                                   const std::map<int,LINALG::Matrix<3,1> >& currentpositions,
                                   Epetra_Vector&                            discol,
                                   Teuchos::RCP<CONTACT::Beam3cmanager>      beamcmanager,
                                   bool                                      printscreen=true);

  /*! \brief determines whether or not a crosslinker element can be set in case of a loom-type network (one horizontal filament/ n vertical filaments)
   *  \param LID              (in) local Ids of nodes
   *  \param currentpositions (in) node positions in column map format (fully overlapping)
   *  \param bspottriadscol   (in) nodal triads as quaternions in column map format (fully overlapping)
   */
  bool SetCrosslinkerLoom(Epetra_SerialDenseMatrix& LID, const std::map<int, LINALG::Matrix<3, 1> >& currentpositions, Epetra_MultiVector& bspottriadscol);

  /*!\brief check whether element is "broken" due to periodic boundary conditions
    \param coord  (in)  set of nodal coordinates
    \param cut    (in/out)  in case of a broken element, cut holds component and information about the location of the nodes of the element*/
  bool CheckForBrokenElement(LINALG::SerialDenseMatrix& coord, LINALG::SerialDenseMatrix& cut);

  /*!\brief return the current coordinates and DOF LIDs of all nodes of given element
     \param element (in)      an element
     \param dis     (in)      current displacement vector (row map)
     \param coord   (in/out)  current node coordinates of element
     \param lids    (in/out)  nodal DOF LIDs of the element*/
  void GetElementNodeCoords(DRT::Element* element, Teuchos::RCP<Epetra_Vector> dis, LINALG::SerialDenseMatrix& coord, std::vector<int>* lids=NULL);
  /*!\brief update locations of force sensors
     \param sensornodes   (in) vector holding the GIDs of the sensor nodes
     \param oscdir        (in) direction (DOF) of oscillation*/
  void UpdateForceSensors(std::vector<int>& sensornodes, int oscdir);

  //! \brief Update positions and rotations of binding spots
  void UpdateBindingSpots(const Epetra_Vector& discol, std::map<int,LINALG::Matrix<3,1> >& currentpositions);

  //! \brief adjust off-rate by considering forces in the crosslinker
  void ForceDependentOffRate(const double& dt, const double& koff0, Epetra_Vector* punlink, Epetra_Vector& discol);

  /*! \brief reduce number of currently existing crosslinker elements by numtoreduce
   *  \param (in) numtoreduce number of crosslinkers by which overall number of crosslinkers is reduced*/
  void ReduceNumOfCrosslinkersBy(const int numtoreduce);

  /*! \brief Search neighbour filament nodes to a crosslink molecule (when crosslinker diffusion is simulated)
   *  \param currentpositions     (in) column map node positions (fully overlapping)
   *  \param nodeinpartition      (in) vector indicating the partition of the nodes
   *  \param numbond              (in) vector holding the number of established bonds for each crosslink molecule
   *  \param crosslinkerpositions (in) vector holding the positions of the crosslink molecules
   *  \param crosslinkpartitions  (in) vector holding the information to which partial volume the crosslink molecules belong
   *  \param bspottriadscol       (in) nodal triads as quaternions in column map format (fully overlapping)
   *  \param neighboursLID        (in/out) vector holding the LIDs of the nodes close to crosslink molecules
   *  */
  void DetectNeighbourNodes(const std::map<int,LINALG::Matrix<3,1> >&     currentpositions,
                            std::vector<std::vector<std::vector<int> > >* nodeinpartition,
                            Epetra_Vector&                                numbond,
                            Epetra_MultiVector&                           crosslinkerpositions,
                            Epetra_MultiVector&                           crosslinkpartitions,
                            Epetra_MultiVector&                           bspottriadscol,
                            Teuchos::RCP<Epetra_MultiVector>&             neighboursLID);

  //! \brief rotation of a vector around a given fixed axis by angle phirot
  void RotationAroundFixedAxis(LINALG::Matrix<3,1>& axis, LINALG::Matrix<3,1>& vector, const double& phirot);

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
  void DDCorrCurrentStructure(const Epetra_Vector&      disrow,
                              LINALG::Matrix<3,1>*      cog,
                              LINALG::Matrix<3,1>*      centershift,
                              std::vector<int>*         crosslinkerentries,
                              const int&                istep,
                              const std::ostringstream& filename,
                              bool                      filorientoutput=false);

  //! \brief iterate vector in order to obtain close approximation to actual (structure) vector
  void DDCorrIterateVector(const Epetra_Vector& discol, LINALG::Matrix<3,1>* vectorj, const int& maxiterations);

  //! \brief output distances between doubly bound crosslinkers in case of a loom setup
  void LoomOutput(const Epetra_Vector& disrow, const std::ostringstream& filename);

  //! \brief output the coverage of crosslinker binding sites (nodes) and the distribution of bound crosslinkers
  void CrosslinkCoverageOutput(const Epetra_Vector& disrow, const std::ostringstream& filename, bool coverageonly=true);

  //! \brief measure the force between two loom filaments and output it
  void LoomOutputAttraction(const Epetra_Vector& disrow, const std::ostringstream& filename, const int& step);

  //! \brief measure the filaments' internal elastic energy and separation between the hoop and the fixed end
  void LoomOutputElasticEnergy(const Epetra_Vector& disrow, const double& dt, const std::ostringstream& filename);

  //! \brief simply returns the number of free, one-bonded, and two-bonded crosslink molecules
  void CrosslinkCount(const std::ostringstream& filename);

  //! \brief returns polymer length per unit volume with respect to radial distance to COG of the network
  void ComputeLocalMeshSize(const Epetra_Vector& disrow, LINALG::Matrix<3,1>& centershift, const int &istep);

  //! \brief Viscoelastic output
  void ViscoelasticityOutput(const double& time, const Epetra_Vector& dis, const Epetra_Vector& fint, std::ostringstream& filename);

  //! \brief check whether a crosslink connects two filaments or occupies two binding spots of the same element and therefore is to be visualized differently
  bool CheckForKinkedVisual(int eleid);

  //! \brief gmsh output of (periodic boundary) box around "boxcenter" with edgelength "length"
  void GmshOutputBox(double color, LINALG::Matrix<3,1>* boxcenter, std::vector<double>& dimension, const std::ostringstream *filename, bool barrier=true);

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
   *  \param dis              (in) latest displacement
   *  \param disi             (in) latest displacement increment in row map format
   *  \param mean             (in) mean of the normal distribution
   *  \param standarddev      (in) standard deviation of the normal distribution
   *  \param init             (in) switch between methods
   */
  void CrosslinkerDiffusion(const Epetra_Vector& dis,
                            double               mean,
                            double               standarddev,
                            const double&        dt);

  /*! \brief update positions of crosslink molecules if they constitute a crosslinker element
   *  \param currentpositions     (in) column map node positions (fully overlapping)
   *  \param LID                  (in) column map node LIDs
   *  \param crosslinkernumber    (in) position in crosslink molecule vector*/
  void CrosslinkerIntermediateUpdate(const std::map<int,LINALG::Matrix<3,1> >& currentpositions,
                                     const LINALG::SerialDenseMatrix&          LID,
                                     const int&                                crosslinkernumber,
                                     Epetra_MultiVector&                       bspottriadscol,
                                     bool                                      coupledmovement=true);

  /*! \brief map crosslink molecules and filament nodes to the according volume partition and detect nodes in proximity of the molecules
   *  \param currentpositions     (in)      column map node positions (fully overlapping)
   *  \param bspottriadscol      (in)      nodal triads as quaternions in column map format (fully overlapping)
   *  \param neighbourslid        (in/out)  vector holding the node lids within search proximity of crosslink molecules*/
  void PartitioningAndSearch(const std::map<int,LINALG::Matrix<3,1> >& currentpositions, Epetra_MultiVector& bspottriadscol, Teuchos::RCP<Epetra_MultiVector>& neighbourslid);

  /*! \brief comunicate Vector to all participating processors
   *  \param InVec    (in) Source/Input vector
   *  \param OutVec   (in) Target/Output vector
   *  \param doexport (in) export flag
   *  \param doimport (in) import flag
   */
  void CommunicateVector(Epetra_Vector& InVec,
                         Epetra_Vector& OutVec,
                         bool doexport=true,
                         bool doimport=true,
                         bool zerofy=true);

  /*! \brief comunicate MultiVector to all participating processors
   *  \param InVec    (in) Source/Input vector
   *  \param OutVec   (in) Target/Output vector
   *  \param doexport (in) export flag
   *  \param doimport (in) import flag
   */
  void CommunicateMultiVector(Epetra_MultiVector& InVec,
                              Epetra_MultiVector& OutVec,
                              bool doexport=true,
                              bool doimport=true,
                              bool zerofy=true);

  /* \brief set Dirichlet values according to given Dirichlet map)
   *  \param disn            (in)     displacement vector
   *  \param  dbcgids        (in/out)  set holding dbc DOF GIDs (new method)*/
  void DBCSetValues(Teuchos::RCP<Epetra_Vector>        dis,
                    Teuchos::RCP<Epetra_Vector>        deltadbc,
                    Teuchos::RCP<std::set<int> >       dbcgids);

  /* \brief evaluate Dirichlet conditions*/
  void DoDirichletConditionPredefined(DRT::Condition&             cond,
                                      DRT::Discretization&        dis,
                                      const bool                  usetime,
                                      const double                time,
                                      Teuchos::RCP<Epetra_Vector> systemvector,
                                      Teuchos::RCP<Epetra_Vector> systemvectord,
                                      Teuchos::RCP<Epetra_Vector> systemvectordd,
                                      Teuchos::RCP<Epetra_Vector> toggle,
                                      Teuchos::RCP<std::set<int> > dbcgids);

  /*! \brief Determine Dirichlet condition at given time and apply its values to a system vector
   *  \param nodeids      (in)      vector containing DBC Node IDs
   *  \param onoff        (in)      vector containing onoff values for DOFs
   *  \param disn         (in)      displacement vector
   *  \param dirichtoggle (in/out)  vector marking the Dirichlet DOFs (old method)
   *  \param invtoggle    (in/out)  the inverse vector of dirichtoggle (old method)
   *  \param deltadbc     (in)      vector holding the increment in Dirichlet values
   *  \param dbcgids      (in/out)  set holding dbc DOF GIDs (new method)*/
  void DoDirichletConditionPeriodic(std::vector<int>*            nodeids,
                                    std::vector<int>*            onoff,
                                    Teuchos::RCP<Epetra_Vector>  dis,
                                    Teuchos::RCP<Epetra_Vector>  deltadbc,
                                    Teuchos::RCP<std::set<int> > dbcgids);

  /*! \brief Dirichlet conditions needed for viscoelastic measurements*/
  void DBCOscillatoryMotion(Teuchos::ParameterList&     params,
                            Teuchos::RCP<Epetra_Vector> dis,
                            Teuchos::RCP<Epetra_Vector> deltadbc);

  /*! \brief Get DBCs defined in Input file and ad the DBCs DOFs to dbcgids
   *  \param params   (in)  Parameter List
   *  \param dis      (in)  displacement vector
   *  \param dbcgids  (in/out)  set of DOF GIDs on which Dirichlet values are imposed*/
  void DBCGetPredefinedConditions(Teuchos::ParameterList& params,
                            Teuchos::RCP<Epetra_Vector>   dis,
                            Teuchos::RCP<std::set<int> >  dbcgids);

  /*! \brief pin down and release nodes at their location by Dirichlet Values*/
  void DBCPinNodes();

  /*! \brief pin down a node at its location by Dirichlet Values
   *  \param dbcgids               (in/out)  set holding dbc DOF GIDs (new method)
   *  \param dbcmapextractor       (in)      vector containing onoff values for free DOFs */
  void DBCCreateMap(Teuchos::RCP<std::set<int> > dbcgids,
                    Teuchos::RCP<LINALG::MapExtractor> dbcmapextractor);

  //! \brief containing statistical mechanics input parameters
  Teuchos::ParameterList statmechparams_;

  //! \brief how often a time step had to be restarted with new random numbers due to unconverged solution of nonlinear system of equations since simulation start
  int unconvergedsteps_;

  //! \brief time at which output starts
  double starttimeoutput_;

  //! \brief time step size vector (entries are sorted as: [t_eq; t_ktswith; t_i; t_i+1; t_act]
  Teuchos::RCP<std::vector<double> > timestepsizes_;

  //! \brief points in time associated with certain actions
  Teuchos::RCP<std::vector<double> > actiontime_;

  //! \brief holds current position in ACTIONTIME / ACTIONDT
  int timeintervalstep_;

  //! \brief end to end distance of a single filament at the beginning of statistical output
  double endtoendref_;

  //! \brief time step number at the beginning of statistical output
  int istart_;

  //! \brief position of first node in last time step
  LINALG::Matrix<3,1> beginold_;

  //! \brief position of last nodein last time step
  LINALG::Matrix<3,1> endold_;

  //! \brief vector holding the dimensions of the periodic box
  Teuchos::RCP<std::vector<double> > periodlength_;

  //! \brief vector holding the spatial resolution of the search algorithm binding spots x crosslinkers
  Teuchos::RCP<std::vector<int> > searchres_;

  //! \brief octree
  Teuchos::RCP<Beam3ContactOctree> octree_;

  //! \brief number of nodes which are part of the discretization already before starting adding crosslinkers
  int basisnodes_;

  //! \brief number of elements which are part of the discretization already before starting adding crosslinkers
  int basiselements_;

  //! \brief number of file into which output is written (only important for single filament dynamics)
  int outputfilenumber_;

  //! \brief a random number generator object (boost)
  randnumgen randomnumbergen_;

#if (BOOST_MAJOR_VERSION == 1) && (BOOST_MINOR_VERSION >= 47)
  //! \brief generates uniformly distributed random numbers in interval [0;1[
  Teuchos::RCP<boost::uniform_01<randnumgen&> > uniformgen_;
#else
  //! \brief generates uniformly distributed random numbers in interval [0;1[
  Teuchos::RCP<boost::variate_generator<randnumgen&,boost::uniform_01<> > > uniformgen_;
#endif

  //! \brief generates Gaussian random numbers with zero mean and unit standard deviation
  Teuchos::RCP<boost::variate_generator<randnumgen&,boost::normal_distribution<> > > normalgen_;

  //! \brief reference to current discretization
  Teuchos::RCP<DRT::Discretization> discret_;

  //! \brief storing for each node column map LID number of filament to which it belongs
  Teuchos::RCP<Epetra_Vector> filamentnumber_;

  //! \brief storing for degree of freedom whether force has to be measured
  Teuchos::RCP<Epetra_Vector> forcesensor_;

  //! \store accumulated disp
  double sumsquareincpar_;
  double sumsquareincort_;
  double sumrotmiddle_;
  double sumsquareincmid_;
  double sumsquareincrot_;

  //store displacement of midpoint compared to initial state as a 3-vector
  LINALG::Matrix<3,1>sumdispmiddle_;

  //! \brief rotation matrix for transformation from global to local coordinates
  Teuchos::RCP<LINALG::SerialDenseMatrix> trafo_;

  //! \brief column Epetra_Map handling crosslink molecules
  Teuchos::RCP<Epetra_Map> crosslinkermap_;

  //! \brief row Epetra_Map handling crosslink molecules
  Teuchos::RCP<Epetra_Map> transfermap_;

  //! \brief an Epetra_Map handling processor-specific data of the density-density-correlation-function
  Teuchos::RCP<Epetra_Map> ddcorrrowmap_;

  //! \brief an Epetra_Map handling cumulated data of the density-density-correlation-function from all processors
  Teuchos::RCP<Epetra_Map> ddcorrcolmap_;

  // note: the following two maps will be of importance later on. For now, we use noderowmap_ and nodecolmap_
  //! \brief an Epetra_Map handling processor-specific data of the binding spot geometry
  Teuchos::RCP<Epetra_Map> bspotrowmap_;

  //! \brief an Epetra_Map handling cumulated data of the binding spot geometry from all processors
  Teuchos::RCP<Epetra_Map> bspotcolmap_;

  //! \brief a vector memorizing the bonding status of the crosslink molecules
  Teuchos::RCP<Epetra_MultiVector>		crosslinkerbond_;

  //! \brief a vector containing the current postion of the crosslink molecules
  Teuchos::RCP<Epetra_MultiVector>		crosslinkerpositions_;

  //! \brief a vector indicating the number of attached nodes to a crosslinker
  Teuchos::RCP<Epetra_Vector>	numbond_;

  //! \brief a vector indicating crosslinker elements that cover two binding spots on the same filament
  Teuchos::RCP<Epetra_Vector> crosslinkonsamefilament_;

  //! \brief a vector mapping crosslink molecule IDs to crosslinker element IDs
  Teuchos::RCP<Epetra_Vector> crosslink2element_;

  //! \brief a vector mapping crosslinker element IDs to crosslink molecule IDs
  Teuchos::RCP<Epetra_Vector> element2crosslink_;

  //! \brief a vector indicating whether a binding spot is occupied or not
  Teuchos::RCP<Epetra_Vector> bspotstatus_;

  //! \brief a vector managing the binding spot geometry (Column Map)
  Teuchos::RCP<Epetra_Vector> bspotorientations_;

  //! \brief a vector storing the curve parameters xi of the respective binding spots (Column Map)
  Teuchos::RCP<Epetra_Vector> bspotxi_;

  //! \brief a vector managing binding spots to their respective element (Row Map)
  Teuchos::RCP<Epetra_Vector> bspot2element_;

  //! \brief a vector just for the sake of visual quality
  Teuchos::RCP<Epetra_MultiVector> visualizepositions_;

  //! \brief a vector holding information on whether or not a crosslink molecule searches for neighbours (switch for active/passive crosslinker)
  Teuchos::RCP<Epetra_Vector> searchforneighbours_;

  //! \brief a vector holding indices used in DensityDensityCorrOutput()
  Teuchos::RCP<std::vector<double> > startindex_;

  //the following variables are necessary only for repeating a time step with different random numbers after its having failed one time

  //! \brief elements deleted at the beginning when this time step was conducted the last time are packed into variable deletedelements_
  std::vector<std::vector<char> > deletedelements_;

  //! \brief contact elements deleted at the beginning when this time step was conducted the last time are packed into variable deletedelements_
  std::vector<std::vector<char> > deletedcelements_;

  //! \brief GIDs of elements added at the beginning when this time step was conducted the last time
  std::vector<int> addedelements_;

  //! \brief GIDs of elements added at the beginning when this time step was conducted the last time
  std::vector<int> addedcelements_;

  //variables to save certain class variables *_ in *conv_ in order to allow for repeating a time step with completely new random numbers and actions
  //! \brief vector holding converged state of crosslinkerbond_// \brief vector holding converged state of crosslinkerbond_
  Teuchos::RCP<Epetra_MultiVector> crosslinkerbondconv_;
  //! \brief vector holding converged state of crosslinkerpositions_
  Teuchos::RCP<Epetra_MultiVector> crosslinkerpositionsconv_;
  //! \brief vector holding converged state of bspotstatus_
  Teuchos::RCP<Epetra_Vector> bspotstatusconv_;
  //! \brief vector holding converged state of numbond_
  Teuchos::RCP<Epetra_Vector> numbondconv_;
  //! \brief vector holding converged state of crosslinkonsamefilament_
  Teuchos::RCP<Epetra_Vector> crosslinkonsamefilamentconv_;
  //! \brief vector holding converged state of crosslink2element_
  Teuchos::RCP<Epetra_Vector> crosslink2elementconv_;
  //! \brief vector holding converged state of searchforneighbour_
  Teuchos::RCP<Epetra_Vector> searchforneighboursconv_;

  //! \brief column map displacement of the previous time step (needed to calculate off-rate in case of force dependency)
  Teuchos::RCP<Epetra_Vector> disprev_;

  //! \brief structure type
  int structuretype_;
  //! \brief store test volume visualization data for Gmsh Output
  std::vector<LINALG::Matrix<3,1> > testvolumepos_;
  //! \brief characteristic length (radii, thickness) of the detected structure
  double characlength_;
  //! \brief center of gravity of the detected structure
  LINALG::Matrix<3,1> cog_;

  // glue initial set of broken elements at the beginning of oscillation to the upper and lower plate
  //! \brief vector with oscillating Dirichlet nodes
  std::vector<int> oscillnodes_;
  //! \brief vector with fixed Dirichlet nodes
  std::vector<int> fixednodes_;
  //! \brief vector with free nodes
  std::vector<int> freenodes_;
  //! \brief vector of node sets for DBC application
  std::vector<std::vector<int> > dbcnodesets_;
  //! \brief flag signaling the use of the initially found DBC node set throughout the simulation
  bool useinitdbcset_;
}; // class StatMechManager
}

#endif /*STATMECH_MANAGER_H*/
