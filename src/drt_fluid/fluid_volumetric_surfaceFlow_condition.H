/*!----------------------------------------------------------------------
\file fluidwomersleycondition.H
\brief Method to deal with womersley flow profiles

<pre>
Maintainer: Mahmoud Ismail
            ismail@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15268
</pre>

*----------------------------------------------------------------------*/

#ifdef CCADISCRET
#ifndef FLUIDWOMERSLEYCONDITION_H
#define FLUIDWOMERSLEYCONDITION_H


#ifdef PARALLEL
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif /* PARALLEL */

#include <Epetra_Map.h>

#include "fluid_utils_mapextractor.H"
#include "fluid_utils.H"
#include "Teuchos_RefCountPtr.hpp"
#include "../drt_lib/drt_discret.H"
#include "../linalg/linalg_utils.H"
#include "../drt_io/io.H"

using namespace std;
using namespace Teuchos;

namespace FLD
{
namespace UTILS
{

//--------------------------------------------------------------------
// Wrapper class (to be called from outside) for volumetric surface
// flow
//--------------------------------------------------------------------

/*!
\brief Womersley boundary condition wrapper
this class is meant to do some organisation stuff

\author ismail (ismail@lnm.mw.tum.de)
*/
class FluidVolumetricSurfaceFlowWrapper
{
  friend class FluidImplicitTimeInt;


public:
  /*!
  \brief Standard Constructor
  */
  FluidVolumetricSurfaceFlowWrapper( RCP<DRT::Discretization>  actdis,
                                     IO::DiscretizationWriter& output,
                                     double dta);

  /*!
  \brief Destructor
  */
  virtual ~FluidVolumetricSurfaceFlowWrapper();

  /*!
  \brief Wrapper for FluidVolumetricSurfaceBc::FindBoundaryNodes
  */
  void FindBoundaryNodes();

  /*!
  \brief Wrapper for FluidVolumetricSurfaceFlowBc::CenterOfMassCalculation
  */
  void CenterOfMassCalculation();

  /*!
  \brief Wrapper for FluidVolumetricSurfaceFlowBc::RadiiCalculation
  */
  void RadiiCalculation();

  /*!
  \brief Wrapper for FluidVolumetricSurfaceFlowBc::EvaluateVelocities
  */
  void EvaluateVelocities(RCP<Epetra_Vector> bcdof,
                          double             time);

  /*!
  \brief Wrapper for FluidVolumetricSurfaceFlowBc::EvaluateVelocities
  */
  void EvaluateCondMap(RCP<Epetra_Map> &   bcmap );

  void InsertCondVector(Epetra_Vector & vec1, Epetra_Vector & vec2);
  /*!
  \brief Returns the map extractor of the conditions
  */
  void EvaluateMapExtractor(RCP<FLD::UTILS::MapExtractor>&  mapextractor );

  /*!
  \brief Wrapper for FluidVolumetricSurfaceFlowBc::UpdateResidual
  */
  void UpdateResidual(RCP<Epetra_Vector>  residual );


  /*!
  \brief Wrapper for FluidVolumetricSurfaceFlowBc::Output
  */
  void Output(IO::DiscretizationWriter&  output);

  /*!
  \brief Wrapper for FluidVolumetricSurfaceFlowBc::ReadRestart
  */
  void ReadRestart(IO::DiscretizationReader& reader);


private:
  /*!
  \brief all single fluid volumetric surface flow conditions
  */
  map<const int, RCP<class FluidVolumetricSurfaceFlowBc> > fvsf_map_;

  //! fluid discretization
  RCP<DRT::Discretization>      discret_;

  //! output
  IO::DiscretizationWriter&     output_;

  //! 
  RCP<Epetra_Vector>            womersleyvbc_;

  //!
  RCP<FLD::UTILS::MapExtractor> womersley_mp_extractor_;

}; // class FluidWomersleyWrapper



//--------------------------------------------------------------------
// Actual Womersley bc calculation stuff
//--------------------------------------------------------------------
/*!
\brief Womersley boundary condition

\author ismail (ismail@lnm.mw.tum.de)
*/
class FluidVolumetricSurfaceFlowBc
{
  friend class FluidVolumetricSurfaceFlowWrapper;

public:

  /*!
  \brief Standard Constructor
  */
  FluidVolumetricSurfaceFlowBc( RCP<DRT::Discretization>  actdis,
                                IO::DiscretizationWriter& output,
                                double                    dta,
                                int                       condid,
                                int                       surf_numcond,
                                int                       line_numcond);

  /*!
  \brief Empty Constructor
  */
  FluidVolumetricSurfaceFlowBc();

  /*!
  \brief Destructor
  */
  virtual ~FluidVolumetricSurfaceFlowBc() {}

  /*!
  \brief finds the boundary nodes
  */
  void FindBoundaryNodes();

  /*!
  \brief calculates the center of mass
  */
  void CenterOfMassCalculation(RCP<std::vector<double> > coords,
                               RCP<std::vector<double> > normal);


  /*!
  \brief calculates the local radii of all nodes
  */
  void EvalLocalNormalizedRadii();

  /*!
  \brief get the node row map of the womersley condition
  */
  void BuildConditionNodeRowMap(RCP<DRT::Discretization>  dis,
                                string                    condname,
                                int                       condid,
                                int                       condnum,
                                RCP<Epetra_Map> &         cond_noderowmap);

  /*!
  \brief get the dof row map of the womersley condition
  */
  void BuildConditionDofRowMap(RCP<DRT::Discretization>      dis,
                               const string                  condname,
                               int                           condid,
                               int                           condnum  ,
                               RCP<Epetra_Map> &             cond_dofrowmap);

#if 0
  /*!
  \brief Broadcast node information over all processors
  */
  void BroadcastNodeInformation(RCP<DRT::Discretization>  dis,
                                string                    noderowmap,
                                int                       condid,
                                int                       condnum,
                                RCP<Epetra_Map>           cond_noderowmap);
#endif

  /*!
  \brief Evaluate velocities
  */
  void EvaluateVelocities(RCP<Epetra_Vector> bcdof,
                          double             time);

  /*!
  \brief Evaluate velocities
  */
  void Velocities(  RCP<DRT::Discretization>  disc,
                    RCP<Epetra_Vector>        bcdof,
                    RCP<Epetra_Map>           cond_noderowmap,
                    RCP<Epetra_Vector>        local_radii,
                    RCP<Epetra_Vector>        border_radii,
                    RCP<std::vector<double> > normal,
                    RCP<ParameterList>        params);

  /*!
  \brief Polynomail shaped velocity profile
  */
  double PolynomailVelocity(double r, int order);

  /*!
  \brief Womersley shaped velocity profile
  */
  double WomersleyVelocity(double r,
                           double R,
                           double Bn,
                           //complex<double> Bn,
                           int    n,
                           double t);

  /*!
  \brief Corrects the flow rate
  */
  void CorrectFlowRate(RCP<Epetra_Vector> bcdof,
                       double             time);

  /*!
  \brief Calculate the Flowrate on a boundary
  */
  double FlowRateCalculation(double time, int condid);


  /*!
  \brief Calculate the Flowrate on a boundary
  */
  void   DFT(RCP<vector<double> >   f,
             RCP<vector<complex<double> > > & F);



protected:


private:
  /*!
  \brief calculate area at outflow boundary
  */
  double Area( double& density, double& viscosity, int condid);

  /*!
  \brief output
  */
  void Output(IO::DiscretizationWriter&  output, int condnum );

  /*!
  \brief Read restart
  */
  void ReadRestart(IO::DiscretizationReader& reader, int condnum);

  /*!
  \brief Bessel function of orders 0 and 1
  */
  std::complex<double> BesselJ01(std::complex<double> z, bool order);

  /*!
  \brief Interpolation function
  */
  void Interpolate(RCP<std::vector<double> > V1,
                   RCP<std::vector<double> > V2,
                   int index1,
                   int & index2,
                   double time);

protected:

  RCP<Epetra_Vector> womersleybc_;

private:

  //! ID of present condition
  int condid_;

  //! Number of present surface condition
  int condnum_s_;

  //! Number of present line condition
  int condnum_l_;

  //! time period of present cyclic problem
  double period_;

  //! fluid viscosity
  double viscosity_;

  //! fluid density
  double density_;

  //! time step size
  double dta_;

  //! the processor ID from the communicator
  int                           myrank_;

  //! fluid discretization
  RCP<DRT::Discretization>      discret_;

  //! output
  IO::DiscretizationWriter&     output_;

  //! Flowrate array for Womersley conditions
  RCP<std::vector<double> >     flowrates_;

  //! Position at which the next element should be replaced
  //! initialised to zero as the first element will be replaced first
  int flowratespos_;

  //! Number of time steps which fit (exactly) into one period
  int cyclesteps_;

  //! center of mass coordinates
  RCP<std::vector<double> >     cmass_;

  //! avarage normal of the surface
  RCP<std::vector<double> >     normal_;

  //! direction normal of the velocity
  RCP<std::vector<double> >     vnormal_;

  //! a Node row map of the nodes that belong to the current condition
  RCP<Epetra_Map> cond_surfnoderowmap_;

  //! a Node row map of the nodes that belong to border of the current condition
  RCP<Epetra_Map> cond_linenoderowmap_;

  //! a Dof row map of the degrees of freedom that belong to the current condition
  RCP<Epetra_Map> cond_dofrowmap_;

  //! A map of the local radii
  RCP<Epetra_Vector> local_radii_;

  //! A map of corresponding border radii
  RCP<Epetra_Vector> border_radii_;
  
  //! A correction factor that will scale the velocities to the correct volumetric-flow-rate
  double correction_factor_;

  //! initial area of the codition surface
  double area_;

  //! Number of modes
  int n_harmonics_;

  //! order of a polynomial velocity profile
  int order_;

  //! Type of the flow profile
  string flowprofile_type_;

  //! is +1 if inflow, else -1
  double flow_dir_;

  //! flag to compute the avarage surface normal
  bool computeNormal_;

  //! flag to compute the surface center of mass
  bool computeCenterMass_;

  //! flag to correct the flowprofile
  bool correct_flow_;

};  // FluidVolumetricSurfaceFlowBc

}
}

#endif /* FLUIDWOMERSLEYCONDITION_H */

#endif /* CCADISCRET */
