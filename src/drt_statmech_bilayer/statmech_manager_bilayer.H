/*!----------------------------------------------------------------------
\file statmech_manager_bilayer.H

<pre>
Maintainer: Dhrubajoti Mukherjee
            mukherjee@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15276
</pre>

*----------------------------------------------------------------------*/
#ifndef STATMECH_MANAGER_BILAYER_H
#define STATMECH_MANAGER_BILAYER_H

#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterListAcceptorDefaultBase.hpp>
#include <Teuchos_StandardParameterEntryValidators.hpp>
#include "Epetra_SerialDenseMatrix.h"
#include "Epetra_SerialDenseVector.h"
#include "Epetra_Vector.h"
#include "Epetra_MultiVector.h"
#include "Epetra_Map.h"

// include random number libraries of BOOST for statistical forces
#include <boost/random.hpp>
#include <boost/random/mersenne_twister.hpp>

#include "../linalg/linalg_fixedsizematrix.H"
#include "../drt_lib/drt_node.H"

#include <iostream>
#include <iomanip>
#include <cstdio>
#include <sys/stat.h>
#include <math.h>


// forward declarations
namespace DRT
{
  class Condition;
  class Element;
  class Discretization;
  class DiscretizationFaces;
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


// random generator type
typedef boost::mt19937 randnumgen;

namespace STATMECH
{

class StatMechManagerBilayer
{
public:
  //! \brief Standard Constructor
  explicit StatMechManagerBilayer(Teuchos::RCP<DRT::Discretization> discret);

  //! \brief Destructor
  virtual ~StatMechManagerBilayer() {};

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
              bool                                 printscreen=true);


  // Print current area of closed vesicle
  void PrintTotalArea(const Epetra_Vector&  dis);

  // Print current volume of closed vesicle
  void PrintTotalVol(const Epetra_Vector&  dis);

  // Return Gmsh and StatMech Output root path
  const std::string& StatMechRootPath() const{return outputrootpath_;}

  /*! \brief special update for statistical mechanics including setting and deleting crosslinkers in the discretization
   *  \param istep    (in)  i-th time step
   *  \param timen    (in)  target time
   *  \param dt       (in)  current time step size
   *  \param disrow   (in)  current displacement vector
   *  \param stiff    (in)
   *  \param ndim     (in)
   */
  void Update(const int&                                      istep,
              const double&                                   timen,
              const double&                                   dt,
              Epetra_Vector&                                  disrow,
              Teuchos::RCP<LINALG::SparseOperator>&           stiff,
              bool                                            rebuildoctree= true,
              bool                                            printscreen= true);

  /*! \brief update time step size in time integration scheme
   *  \param dt            (in) time step size
   *  \param timeconverged (in) time of last converged timestep
   *  \param initialaset   (in) different behavior depending on call*/
  void UpdateTimeAndStepSize(double& dt,
                             double  timeconverged,
                             bool initialset=false);

  //! \brief initializing special output for statistical mechanics by looking for a suitable name of the outputfile and setting up an empty file with this name
  void InitOutput(const int& ndim,
                  const Epetra_Vector& dis,
                  const int& istep,
                  const double& dt);

  void BuildStatMechRootPath();


  //! \brief seed all random generators of this object properly
  void SeedRandomGenerators(const int seedparameter, const int seedparameter2);

  //! \brief generate gaussian randomnumbers with mean "meanvalue" and standarddeviation "standarddeviation" for parallel use
  void GenerateGaussianRandomNumbers(Teuchos::RCP<Epetra_MultiVector> randomnumbers,
                                     const double                     meanvalue,
                                     const double                     standarddeviation);

  /*! \brief adjusted DBC evaluation for StatMech
   *  \param params          (in)     a parameter list
   *  \param dis             (in)     displacement vector
   *  \param vel             (in)     velocity vector
   *  \param dbcmapextractor (in/out) map extractor for Dirichlet DOFs (new method)*/
  void EvaluateDirichletStatMech(Teuchos::ParameterList&            params,
                                 Teuchos::RCP<Epetra_Vector>        dis,
                                 Teuchos::RCP<Epetra_Vector>        vel,
                                 Teuchos::RCP<LINALG::MapExtractor> dbcmapextractor=Teuchos::null);

  /*! \brief adjusted NBC evaluation for StatMech
   *  \param params       (in) a parameter list
   *  \param disn         (in) new displacement state
   *  \param systemvector (in) force vector
   *  \param systemmatrix (in) stiffness matrix
   * */
  void EvaluateNeumannStatMech(Teuchos::ParameterList&              params,
                               Teuchos::RCP<Epetra_Vector>          disn,
                               Teuchos::RCP<Epetra_Vector>          systemvector,
                               Teuchos::RCP<LINALG::SparseOperator> systemmatrix);

  //! \brief return vector holding the dimensions of the periodic boundary box
//  const Teuchos::RCP<std::vector<double> > GetPeriodLength() const {return periodlength_;}

  //! \brief return statistical mechanics parameter list
  const Teuchos::ParameterList& GetStatMechBilayerParams() const {return statmechBilayerparams_;}

  //! \brief return number of unconverged steps since beginning of the simulation
  const int& NumberOfUnconvergedSteps() const {return unconvergedsteps_;}

  //! \brief return number of unconverged steps since beginning of the simulation
  void UpdateNumberOfUnconvergedSteps();

  //! \brief pass statistical mechanics parameters to parameter list
  void AddStatMechParamsTo(Teuchos::ParameterList&          params,
                           Teuchos::RCP<Epetra_MultiVector> randomnumbers = Teuchos::null);

  void GmshOutput(const Epetra_Vector&          disrow,
                  const std::ostringstream&     filename,
                  const int&                    step,
                  const double&                 time);

protected:

  // don't want = operator and cctor
  StatMechManagerBilayer operator = (const StatMechManagerBilayer& old);

  StatMechManagerBilayer(const StatMechManagerBilayer& old);

  //! \brief starting from the input file, build the path trunk to the StatMech subdirectories
//  void BuildStatMechRootPath();

  //! \brief Initialize periodlength_ vector and search resolution
  void InitializeStatMechValues();

  //! \brief creation of a fully overlapping node column map during statmech manager construction
  void CreateFullyOverlappingNodeMap();

  bool IfBondFlip(const int&    RandEdgeGID,
                  Epetra_Vector&     discol);

  bool IfBondFlipParallel(const int&    RandEdgeGID,
                          Epetra_Vector&     discol);

  bool CheckBoundingBoxQuality(const int&    RandEdgeGID,
                               Epetra_Vector&     discol);

  bool CheckBoundingBoxQualityParallel(const int&    RandEdgeGID,
                                       Epetra_Vector&     discol);

  void DoBondFlip(const int&    RandEdgeGID);

  void DoBondFlipParallel(DRT::Discretization& mydiscret,
                          DRT::DiscretizationFaces& myfacediscret,
                          const int&    RandEdgeGID);

  void ModifyElementConnectivity(DRT::Discretization& mydiscret,
                                 std::map<int,int>& delEle,
                                 std::map<int,int>& oldnew );

  void ModifyElementConnectivityNew(DRT::Discretization& mydiscret,
                                 std::map<int,int>& delEle,
                                 std::vector< std::map<int,int> >& oldnew);

  //! \brief initialize filament numbers
  void InitializeFilamentNumbers();

  //! \brief initialize force sensors
  void InitializeForceSensors();

  //! \brief set start step of the simulation (important after restart)
  void SetStartStep(Teuchos::ParameterList& parameters);

  /*!\brief return the inclusive angle between two given vectors */
  double GetTheta(LINALG::Matrix<3,1> &T1,
                  LINALG::Matrix<3,1> &T2);

  //! \brief generates a vector with a random permutation of the numbers between 0 and N - 1
  std::vector<int> Permutation(const int& N);

  /*! \brief writes element orientation and coordinates of fiber elements to an output file
      \param discol     (in) displacement in column map format
      \param filename  (in) filename for the first output file (coordinates, element GIDs and filament numbers)
   */
  void FilamentOrientations(const Epetra_Vector&               discol,
                            std::vector<LINALG::Matrix<3,1> >* normedvectors,
                            const std::ostringstream&          filename,
                            bool                               fileoutput);

  LINALG::Matrix <1,3> GetSpatialPosition(DRT::Node * node, Epetra_Vector&  discol);

  // Activate bond-flip mechanism
  void EvaluateBondFlips(const int&                             istep,
                        const double&                          timen,
                        const double&                          dt,
                        Epetra_Vector&                         discol,
                        bool                                   printscreen= true);

  // Activate bond-flip mechanism
  void EvaluateBondFlipsParallel(const int&                     istep,
                                 const double&                  timen,
                                 const double&                  dt,
                                 Epetra_Vector&                 discol,
                                 bool                           printscreen= true);

  //! \brief obtain distribution for spherical coordinates (all vectors in and out)
  void SphericalCoordsDistribution(const Epetra_Vector& disrow,
                                   Epetra_Vector&       phibinsrow,
                                   Epetra_Vector&       thetabinsrow,
                                   Epetra_Vector&       costhetabinsrow);

  /*! \brief comunicate MultiVector to all participating processors
   *  \param InVec    (in) Source/Input vector
   *  \param OutVec   (in) Target/Output vector
   *  \param doexport (in) export flag
   *  \param doimport (in) import flag
   *  \param zerofy   (in) zero out OutVec flag
   */
  void CommunicateMultiVector(Teuchos::RCP<Epetra_MultiVector> InVec,
                              Teuchos::RCP<Epetra_MultiVector> OutVec,
                              bool                             doexport=true,
                              bool                             doimport=true,
                              bool                             zerofy=true,
                              bool                             exportinsert=false);


  /******************Begin Methods for GMSH visualization*****************/

//  void GmshPrepareVisualization(const Epetra_Vector& dis);

  void PrintGmshTriangleToStream(std::stringstream& gmshfilecontent,
                                 LINALG::SerialDenseMatrix coord);

  void PrintGmshNewBondsToStream(std::stringstream& gmshfilecontent,
                                 LINALG::Matrix <1,3> & MasterNodeSpatial,
                                 LINALG::Matrix <1,3> & SlaveNodeSpatial,
                                 const double color);

  /******************End Methods for GMSH visualization*****************/

  //! \brief containing statistical mechanics input parameters
  Teuchos::ParameterList statmechBilayerparams_;

  //! \brief output root path for Gmsh and StatMech output
  std::string outputrootpath_;

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

  //! \brief holds position in ACTIONTIME from where on BCs are applied
  int bctimeindex_;

  //! \brief time step number at the beginning of statistical output
  int istart_;

  //! \brief position of first node in last time step
  LINALG::Matrix<3,1> beginold_;

  //! \brief position of last nodein last time step
  LINALG::Matrix<3,1> endold_;

  //! \brief number of time steps where output has already been written
  int numoutputstep_;

  //! \brief vector holding the dimensions of the periodic box
  Teuchos::RCP<std::vector<double> > periodlength_;

  //! \brief vector holding the NodeIDs of new faces
//  Teuchos::RCP<std::vector<LINALG::Matrix<1,2> > > new_faces_NodeID_;

  //! \brief vector holding the NodeIDs of new faces
  std::map<int,std::vector<int> > new_faces_NodeID_;

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

  /// face discretization (only initialized for edge-based stabilization)
  Teuchos::RCP<DRT::DiscretizationFaces> facediscret_;

}; // class StatMechManagerBilayer
}

#endif /*STATMECH_MANAGER_BILAYER_H*/
