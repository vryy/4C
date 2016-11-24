/*!-----------------------------------------------------------------------------------------------*
\file turbulence_statistics_oracles.H

\brief statistical data processing for ORACLES premixed combustion test case

<pre>
\level 2
\maintainer Benjamin Krank
            krank@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15252
</pre>
 *------------------------------------------------------------------------------------------------*/


#ifndef TURBULENCE_STATISTICS_ORACLES_H_
#define TURBULENCE_STATISTICS_ORACLES_H_

#ifdef PARALLEL
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "../drt_lib/drt_discret.H"
#include "../linalg/linalg_utils.H"
#include "../drt_inpar/inpar_fluid.H"

namespace DRT
{
class DofSet;
}

namespace COMBUST
{
  /*!

    \brief statistics manager for ORACLES problem

    \detailed

    computes profiles of time averaged velocities

    \author henke
    \date 06/08
   */
class TurbulenceStatisticsORACLES
{

public:

    /*!
  \brief constructor

  */
  TurbulenceStatisticsORACLES(Teuchos::RCP<DRT::Discretization> discret,
                              Teuchos::ParameterList&           params,
                              const std::string&                statistics_outfilename,
                              const std::string&                geotype,
                              const bool                        withscatra);

  /*!
  \brief destructor

  */
  virtual ~TurbulenceStatisticsORACLES();


  //! @name functions for averaging

  void DoTimeSample(
      Teuchos::RCP<Epetra_Vector>     velnp,
      Teuchos::RCP<Epetra_Vector>     force,
      Teuchos::RCP<Epetra_Vector>     phi,
      Teuchos::RCP<const DRT::DofSet> stddofset
  );

  void ExtractProfiles();

  void ExtractSetOfProfiles(
    const std::vector<double>&  x1locations,
    const std::vector<double>&  x2locations,
    LINALG::SerialDenseMatrix& profilesu,
    LINALG::SerialDenseMatrix& profilesv,
    LINALG::SerialDenseMatrix& profilesw,
    LINALG::SerialDenseMatrix& profilesp,
    LINALG::SerialDenseMatrix& profilesuu,
    LINALG::SerialDenseMatrix& profilesvv,
    LINALG::SerialDenseMatrix& profilesww,
    LINALG::SerialDenseMatrix& profilespp,
    LINALG::SerialDenseMatrix& profilesuv,
    LINALG::SerialDenseMatrix& profilesuw,
    LINALG::SerialDenseMatrix& profilesvw,
    LINALG::SerialDenseMatrix& profilesg
  );

  void TimeAverageStatistics();

  void OutputStatistics(int step);

  void WriteStatisticsFile(
      const int               step,
      const std::string&           suffix,
      const int               x1pos,
      const std::vector<double>&   x2locations,
      double* profileu,
      double* profilev,
      double* profilew,
      double* profilep,
      double* profileuu,
      double* profilevv,
      double* profileww,
      double* profilepp,
      double* profileuv,
      double* profileuw,
      double* profilevw,
      double* profileg
  );

  void ClearStatistics();

private:

  /*!
  \brief sort criterium for double values up to a tolerance of 10-9

  This is used to create sets of doubles (e.g. coordinates)

  */
  class LineSortCriterion
  {
  public:
    bool operator() (
        const double& p1,
        const double& p2) const
    {
      return (p1 < p2 - 1E-9);
    }

  protected:
  private:
  };

  //! export location from this proc to all procs
  void ExportLocation(std::set<double,LineSortCriterion>& location);

  //! pointer to discretisation
  Teuchos::RCP<DRT::Discretization> discret_;

  //! dofset containing fluid standard dofs (no XFEM dofs)
  Teuchos::RCP<const DRT::DofSet> stddofset_;

  //! pointer to ScaTra discretisation
  Teuchos::RCP<DRT::Discretization> scatradis_;

  //! parameter list of time integration
  Teuchos::ParameterList& params_;

  //! name of statistics output file, despite the ending
  const std::string statistics_outfilename_;

  //! pointer to fluid solution vel/pres and residual (force)
  Teuchos::RCP<Epetra_Vector> vel_;
  Teuchos::RCP<Epetra_Vector> force_;

  //! pointer to scalar transport solution phi (G-function) on fluid RowNodeMap
  Teuchos::RCP<Epetra_Vector> phi_;

  //! boolean indicating scalar (G-function) transport statistics
  const bool withscatra_;

  //! step height
  const double h_; // 29.9mm

  //! bounds of ORACLES domain in x1-direction
  const double x1min_; // -149.5mm
  const double x1max_; //  478.4mm
  //! bounds of ORACLES domain in x2-direction
  const double x2min_; //   0.0mm
  const double x2max_; // 130.6mm
  //! bounds of ORACLES domain in x2-direction
  const double x2inflowmin_; //  29.9mm
  const double x2inflowmax_; // 100.7mm
  //!
  const double x2inflowchannelmin_;
  const double x2inflowchannelmax_;
  //! bounds of ORACLES domain in x3-direction
  const double x3min_; // -75.25mm
  const double x3max_; //  75.25mm

  const double utau_;

  const double midupchan_;   // x2-coordinate mid line upper channel
  const double midlowchan_;  // x2-coordinate mid line lower channel
  const double midchamber_;  // x2-coordinate mid line chamber
  const double x2first_;     // x2-coordinate first node away from wall
  const double x2max_first_; // x2-coordinate first node away from top wall of chamber
  const double x2min_first_; // x2-coordinate first node away from bottom wall of chamber
  const double x2inflowchannelmax_first_; // x2-coordinate first node away from top wall of upper inflow channel
  const double x2inflowchannelmin_first_; // x2-coordinate first node away from bottom wall of upper inflow channel

  // coordinates of nodes in x2-direction
  Teuchos::RCP<std::vector<double> > x2_inflow_; // x2-coordinates of nodes at x =-5h (both inflow channels)
  Teuchos::RCP<std::vector<double> > x2_0h_;     // x2-coordinates of nodes at x = 0h (step/expansion)
  Teuchos::RCP<std::vector<double> > x2_1h_;     // x2-coordinates of nodes at x = 1h
  Teuchos::RCP<std::vector<double> > x2_2h_;     // x2-coordinates of nodes at x = 2h
  Teuchos::RCP<std::vector<double> > x2_3h_;     // x2-coordinates of nodes at x > 3h

  // coordinates of nodes in x1-direction
  Teuchos::RCP<std::vector<double> > x1_midupperchannel_;   // x1-coordinates of nodes at y = 85.5mm = 2.86h (mid line upper channel)
  Teuchos::RCP<std::vector<double> > x1_midlowerchannel_;   // x1-coordinates of nodes at y = 45.1mm = 1.51h (mid line lower channel)
  Teuchos::RCP<std::vector<double> > x1_midchamber_;        // x1-coordinates of nodes at y = 65.3mm = 2.18h (mid line chamber)
  Teuchos::RCP<std::vector<double> > x1_wallchamber_;       // x1-coordinates of nodes at wall of chamber
  Teuchos::RCP<std::vector<double> > x1_wallinflowchannel_; // x1-coordinates of nodes at wall of inflow channels
  // number of nodes in x1-direction
  size_t numx1_midinflowchannel_;
  size_t numx1_midchamber_;
  size_t numx1_wallchamber_;
  size_t numx1_wallinflowchannel_;


  //! coordinates of locations in x1-, x2 and x3-direction for statistical evaluation
  LINALG::Matrix<21,1> x1positions_;
  const double x3position_;

  //! sum over nodal forces at wall
  Teuchos::RCP<LINALG::SerialDenseMatrix> wallforceu_;
  Teuchos::RCP<LINALG::SerialDenseMatrix> wallforcev_;
  Teuchos::RCP<LINALG::SerialDenseMatrix> wallforcew_;

  Teuchos::RCP<LINALG::SerialDenseMatrix> wallvelu_;
  Teuchos::RCP<LINALG::SerialDenseMatrix> wallvelv_;
  Teuchos::RCP<LINALG::SerialDenseMatrix> wallvelw_;
  Teuchos::RCP<LINALG::SerialDenseMatrix> wallp_;

  Teuchos::RCP<LINALG::SerialDenseMatrix> wallinflowchannelforceu_;
  Teuchos::RCP<LINALG::SerialDenseMatrix> wallinflowchannelforcev_;
  Teuchos::RCP<LINALG::SerialDenseMatrix> wallinflowchannelforcew_;

  Teuchos::RCP<LINALG::SerialDenseMatrix> wallinflowchannelvelu_;
  Teuchos::RCP<LINALG::SerialDenseMatrix> wallinflowchannelvelv_;
  Teuchos::RCP<LINALG::SerialDenseMatrix> wallinflowchannelvelw_;
  Teuchos::RCP<LINALG::SerialDenseMatrix> wallinflowchannelp_;

  //! first order momentum
  Teuchos::RCP<LINALG::SerialDenseMatrix> vertinflowu_;
  Teuchos::RCP<LINALG::SerialDenseMatrix> vertinflowv_;
  Teuchos::RCP<LINALG::SerialDenseMatrix> vertinfloww_;
  Teuchos::RCP<LINALG::SerialDenseMatrix> vertinflowp_;
  Teuchos::RCP<LINALG::SerialDenseMatrix> vertinflowg_;

  Teuchos::RCP<LINALG::SerialDenseMatrix> vertmixingu_;
  Teuchos::RCP<LINALG::SerialDenseMatrix> vertmixingv_;
  Teuchos::RCP<LINALG::SerialDenseMatrix> vertmixingw_;
  Teuchos::RCP<LINALG::SerialDenseMatrix> vertmixingp_;
  Teuchos::RCP<LINALG::SerialDenseMatrix> vertmixingg_;

  Teuchos::RCP<LINALG::SerialDenseMatrix> vert1hu_;
  Teuchos::RCP<LINALG::SerialDenseMatrix> vert1hv_;
  Teuchos::RCP<LINALG::SerialDenseMatrix> vert1hw_;
  Teuchos::RCP<LINALG::SerialDenseMatrix> vert1hp_;
  Teuchos::RCP<LINALG::SerialDenseMatrix> vert1hg_;

  Teuchos::RCP<LINALG::SerialDenseMatrix> vert2hu_;
  Teuchos::RCP<LINALG::SerialDenseMatrix> vert2hv_;
  Teuchos::RCP<LINALG::SerialDenseMatrix> vert2hw_;
  Teuchos::RCP<LINALG::SerialDenseMatrix> vert2hp_;
  Teuchos::RCP<LINALG::SerialDenseMatrix> vert2hg_;

  Teuchos::RCP<LINALG::SerialDenseMatrix> vertchamberu_;
  Teuchos::RCP<LINALG::SerialDenseMatrix> vertchamberv_;
  Teuchos::RCP<LINALG::SerialDenseMatrix> vertchamberw_;
  Teuchos::RCP<LINALG::SerialDenseMatrix> vertchamberp_;
  Teuchos::RCP<LINALG::SerialDenseMatrix> vertchamberg_;

  //! second order momentum
  Teuchos::RCP<LINALG::SerialDenseMatrix> vertinflowuu_;
  Teuchos::RCP<LINALG::SerialDenseMatrix> vertinflowvv_;
  Teuchos::RCP<LINALG::SerialDenseMatrix> vertinflowww_;
  Teuchos::RCP<LINALG::SerialDenseMatrix> vertinflowpp_;

  Teuchos::RCP<LINALG::SerialDenseMatrix> vertmixinguu_;
  Teuchos::RCP<LINALG::SerialDenseMatrix> vertmixingvv_;
  Teuchos::RCP<LINALG::SerialDenseMatrix> vertmixingww_;
  Teuchos::RCP<LINALG::SerialDenseMatrix> vertmixingpp_;

  Teuchos::RCP<LINALG::SerialDenseMatrix> vert1huu_;
  Teuchos::RCP<LINALG::SerialDenseMatrix> vert1hvv_;
  Teuchos::RCP<LINALG::SerialDenseMatrix> vert1hww_;
  Teuchos::RCP<LINALG::SerialDenseMatrix> vert1hpp_;

  Teuchos::RCP<LINALG::SerialDenseMatrix> vert2huu_;
  Teuchos::RCP<LINALG::SerialDenseMatrix> vert2hvv_;
  Teuchos::RCP<LINALG::SerialDenseMatrix> vert2hww_;
  Teuchos::RCP<LINALG::SerialDenseMatrix> vert2hpp_;

  Teuchos::RCP<LINALG::SerialDenseMatrix> vertchamberuu_;
  Teuchos::RCP<LINALG::SerialDenseMatrix> vertchambervv_;
  Teuchos::RCP<LINALG::SerialDenseMatrix> vertchamberww_;
  Teuchos::RCP<LINALG::SerialDenseMatrix> vertchamberpp_;

  //! second order momentum
  Teuchos::RCP<LINALG::SerialDenseMatrix> vertinflowuv_;
  Teuchos::RCP<LINALG::SerialDenseMatrix> vertinflowuw_;
  Teuchos::RCP<LINALG::SerialDenseMatrix> vertinflowvw_;

  Teuchos::RCP<LINALG::SerialDenseMatrix> vertmixinguv_;
  Teuchos::RCP<LINALG::SerialDenseMatrix> vertmixinguw_;
  Teuchos::RCP<LINALG::SerialDenseMatrix> vertmixingvw_;

  Teuchos::RCP<LINALG::SerialDenseMatrix> vert1huv_;
  Teuchos::RCP<LINALG::SerialDenseMatrix> vert1huw_;
  Teuchos::RCP<LINALG::SerialDenseMatrix> vert1hvw_;

  Teuchos::RCP<LINALG::SerialDenseMatrix> vert2huv_;
  Teuchos::RCP<LINALG::SerialDenseMatrix> vert2huw_;
  Teuchos::RCP<LINALG::SerialDenseMatrix> vert2hvw_;

  Teuchos::RCP<LINALG::SerialDenseMatrix> vertchamberuv_;
  Teuchos::RCP<LINALG::SerialDenseMatrix> vertchamberuw_;
  Teuchos::RCP<LINALG::SerialDenseMatrix> vertchambervw_;

  //! viscosity
  double visc_;
  //!density
  double dens_;

  //! number of samples taken in this record (corresponds to number of averaged time steps)
  size_t numsamp_;

  //! number of records written
  size_t countrecord_;

};

} // namespace COMBUST


#endif
