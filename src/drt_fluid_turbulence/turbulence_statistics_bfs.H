/*----------------------------------------------------------------------*/
/*!

\brief Write (time and space) averaged values to file for
turbulent flow over a backward-facing step

o Create sets for various evaluation lines in domain
  (Construction based on a round robin communication pattern):
  - 21 lines in x2-direction
  - lines along upper and lower wall

o loop nodes closest to centerlines

  - generate 4 toggle vectors (u,v,w,p), for example

                            /  1  u dof in homogeneous plane
                 toggleu_  |
                            \  0  elsewhere

  - pointwise multiplication velnp.*velnp for second order
    moments

o values on lines are averaged in time over all steps between two
  outputs

Required parameters are the number of velocity degrees of freedom (3)
and the basename of the statistics outfile. These parameters are
expected to be contained in the fluid time integration parameter list
given on input.

This method is intended to be called every upres_ steps during fluid
output.

\maintainer Martin Kronbichler

\level 2

*/
/*----------------------------------------------------------------------*/

#ifndef TURBULENCE_STATISTICS_BFS_H
#define TURBULENCE_STATISTICS_BFS_H

#ifdef PARALLEL
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif /* PARALLEL */

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "../drt_lib/drt_discret.H"
#include "../linalg/linalg_utils.H"
#include "../drt_inpar/inpar_fluid.H"


namespace FLD
{
  class TurbulenceStatisticsBfs
  {
   public:
    /*!
    \brief Standard Constructor (public)

        o Create sets for lines

    o Allocate distributed vector for squares

    */
    TurbulenceStatisticsBfs(Teuchos::RCP<DRT::Discretization> actdis,
        Teuchos::ParameterList& params, const std::string& statistics_outfilename,
        const std::string& geotype);

    /*!
    \brief Destructor

    */
    virtual ~TurbulenceStatisticsBfs();


    //! @name functions for averaging

    /*!
    \brief The values of velocity and its squared values are added to
    global vectors. This method allows to do the time average of the
    nodal values after a certain amount of timesteps.
    */
    void DoTimeSample(Teuchos::RCP<Epetra_Vector> velnp, Teuchos::RCP<Epetra_Vector> stresses);

    /*!
    \brief The values of velocity, pressure, temperature and its squared
    values are added to global vectors. This method allows to do the time
    average of the nodal values after a certain amount of timesteps.
    */
    void DoLomaTimeSample(
        Teuchos::RCP<Epetra_Vector> velnp, Teuchos::RCP<Epetra_Vector> scanp, const double eosfac);

    /*!
    \brief The values of velocity, pressure, phi and its squared
    values are added to global vectors. This method allows to do the time
    average of the nodal values after a certain amount of timesteps.
    */
    void DoScatraTimeSample(Teuchos::RCP<Epetra_Vector> velnp, Teuchos::RCP<Epetra_Vector> scanp);

    /*!
    \brief Dump the result to file.

    step on input is used to print the timesteps which belong to the
    statistic to the file
    */

    void DumpStatistics(int step);

    /*!
    \brief Dump the result to file for low-Mach-number flow.

    step on input is used to print the timesteps which belong to the
    statistic to the file
    */

    void DumpLomaStatistics(int step);

    /*!
    \brief Dump the result to file for turbulent flow with passive scalar.

    step on input is used to print the timesteps which belong to the
    statistic to the file
    */

    void DumpScatraStatistics(int step);


   protected:
    /*!
    \brief sort criterium for double values up to a tolerance of 10-9

    This is used to create sets of doubles (e.g. coordinates)

    */
    class LineSortCriterion
    {
     public:
      bool operator()(const double& p1, const double& p2) const { return (p1 < p2 - 1E-9); }

     protected:
     private:
    };

   private:
    //! geometry of DNS of incompressible flow over bfs by Le, Moin and Kim or geometry of Avancha
    //! and Pletcher of LES of flow over bfs with heating
    enum GeoType
    {
      none,
      geometry_DNS_incomp_flow,
      geometry_LES_flow_with_heating,
      geometry_EXP_vogel_eaton
    };

    //! number of samples taken
    int numsamp_;

    //! number of coordinates in x1- and x2-direction
    int numx1coor_;
    int numx2coor_;

    //! number of locations in x1- and x2-direction for statistical evaluation
    int numx1statlocations_;
    int numx2statlocations_;
    int numx1supplocations_;

    //! bounds for extension of backward-facing step in x2-direction
    double x2min_;
    double x2max_;

    //! bounds for extension of backward-facing step in x3-direction
    double x3min_;
    double x3max_;

    //! The discretisation (required for nodes, dofs etc;)
    Teuchos::RCP<DRT::Discretization> discret_;

    //! parameter list
    Teuchos::ParameterList& params_;

    //! geometry of DNS of incompressible flow over bfs by Le, Moin and Kim or geometry of Avancha
    //! and Pletcher of LES of flow over bfs with heating
    FLD::TurbulenceStatisticsBfs::GeoType geotype_;

    //! boolean indicating turbulent inflow channel discretization
    const bool inflowchannel_;
    //! x-coordinate of outflow of inflow channel
    const double inflowmax_;

    //! name of statistics output file, despite the ending
    const std::string statistics_outfilename_;

    //! pointer to vel/pres^2 field (space allocated in constructor)
    Teuchos::RCP<Epetra_Vector> squaredvelnp_;
    //! pointer to T^2 field (space allocated in constructor)
    Teuchos::RCP<Epetra_Vector> squaredscanp_;
    //! pointer to 1/T field (space allocated in constructor)
    Teuchos::RCP<Epetra_Vector> invscanp_;
    //! pointer to (1/T)^2 field (space allocated in constructor)
    Teuchos::RCP<Epetra_Vector> squaredinvscanp_;

    //! toogle vectors: sums are computed by scalarproducts
    Teuchos::RCP<Epetra_Vector> toggleu_;
    Teuchos::RCP<Epetra_Vector> togglev_;
    Teuchos::RCP<Epetra_Vector> togglew_;
    Teuchos::RCP<Epetra_Vector> togglep_;

    //! available x1- and x2-coordinates
    Teuchos::RCP<std::vector<double>> x1coordinates_;
    Teuchos::RCP<std::vector<double>> x2coordinates_;

    //! coordinates of locations in x1- and x2-direction for statistical evaluation
    LINALG::Matrix<21, 1> x1statlocations_;
    LINALG::Matrix<2, 1> x2statlocations_;

    //! coordinates of supplementary locations in x2-direction for velocity derivative
    LINALG::Matrix<2, 1> x2supplocations_;

    //! coordinates of supplementary locations in x1-direction for statistical evaluation
    //! (check of inflow profile)
    LINALG::Matrix<10, 1> x1supplocations_;

    //! matrices containing values
    Teuchos::RCP<Epetra_SerialDenseMatrix> x1sumu_;
    Teuchos::RCP<Epetra_SerialDenseMatrix> x1sump_;

    Teuchos::RCP<Epetra_SerialDenseMatrix> x1sumrho_;
    Teuchos::RCP<Epetra_SerialDenseMatrix> x1sumT_;

    Teuchos::RCP<Epetra_SerialDenseMatrix> x1sumtauw_;

    Teuchos::RCP<Epetra_SerialDenseMatrix> x2sumu_;
    Teuchos::RCP<Epetra_SerialDenseMatrix> x2sumv_;
    Teuchos::RCP<Epetra_SerialDenseMatrix> x2sumw_;
    Teuchos::RCP<Epetra_SerialDenseMatrix> x2sump_;

    Teuchos::RCP<Epetra_SerialDenseMatrix> x2sumrho_;
    Teuchos::RCP<Epetra_SerialDenseMatrix> x2sumT_;

    Teuchos::RCP<Epetra_SerialDenseMatrix> x2sumsqu_;
    Teuchos::RCP<Epetra_SerialDenseMatrix> x2sumsqv_;
    Teuchos::RCP<Epetra_SerialDenseMatrix> x2sumsqw_;
    Teuchos::RCP<Epetra_SerialDenseMatrix> x2sumsqp_;

    Teuchos::RCP<Epetra_SerialDenseMatrix> x2sumsqrho_;
    Teuchos::RCP<Epetra_SerialDenseMatrix> x2sumsqT_;

    Teuchos::RCP<Epetra_SerialDenseMatrix> x2sumuv_;
    Teuchos::RCP<Epetra_SerialDenseMatrix> x2sumuw_;
    Teuchos::RCP<Epetra_SerialDenseMatrix> x2sumvw_;

    Teuchos::RCP<Epetra_SerialDenseMatrix> x2sumrhou_;
    Teuchos::RCP<Epetra_SerialDenseMatrix> x2sumuT_;
    Teuchos::RCP<Epetra_SerialDenseMatrix> x2sumrhov_;
    Teuchos::RCP<Epetra_SerialDenseMatrix> x2sumvT_;

    void convertStringToGeoType(const std::string& geotype);
  };

}  // namespace FLD

#endif  // not TURBULENCE_STATISTICS_BFS_H
