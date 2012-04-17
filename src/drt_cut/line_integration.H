#ifndef LINE_INTEGRATION_H
#define LINE_INTEGRATION_H

#include<vector>
#include "../linalg/linalg_fixedsizematrix.H"

#include "../../src/drt_fem_general/drt_utils_gausspoints.H"

/*!
\brief Performs integration of  base functions along the line using the standard Gaussian rule
*/
class LineIntegration
{
  public:
    LineIntegration(LINALG::Matrix<2,2> endPts, int inte_num, std::vector<double> alpha, bool bcellInt)
      :end_pts_(endPts),
       inte_num_(inte_num),
       alpha_(alpha),
       bcellInt_(bcellInt)
      {
      }

    ~LineIntegration(){};

    /*!
    \brief Integration of a function along the line using standard Gaussian rule
    */
    double integrate_line();

    /*!
    \brief Choose the base function to be integrated
    */
    void set_integ_type(std::string inttype){intType_ = inttype;}

    /*!
    \brief Transform the Gauss integration point available in the limit (-1,1) to the actual line coordinates
     */
    template <DRT::Element::DiscretizationType celldistype>
    void Transform( const LINALG::Matrix<2,2> &xyze,
                    const double                    & eta,
                    LINALG::Matrix<2,1>       & x_gp_lin,
                    LINALG::Matrix<2,1>       & normal,
                    double                    & drs);

  private:
    /*!
    \brief Compute the normal vector of the line
    */
    LINALG::Matrix<2,1> compute_normal();

    //first index decides the x or y coordinate, second index decides the start point or end point
    LINALG::Matrix<2,2> end_pts_;
    int inte_num_;
    std::vector<double> alpha_;
    bool bcellInt_;
    std::string intType_;
                
};

#endif
