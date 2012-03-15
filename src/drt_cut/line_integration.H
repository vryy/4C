#ifndef LINE_INTEGRATION_H
#define LINE_INTEGRATION_H

#include<vector>
#include "../linalg/linalg_fixedsizematrix.H"

/*!
\brief Performs integration of  base functions along the line using the standard Gaussian rule
*/
class LineIntegration
{
  public:
    LineIntegration(std::vector<double> coord1, std::vector<double> coord2, int inte_num, std::vector<double> alpha, bool bcellInt)
      :point_begin_(coord1),
       point_end_(coord2),
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

  private:
    /*!
    \brief Compute the normal vector of the line
    */
    LINALG::Matrix<2,1> compute_normal();

    /*!
    \brief Get the weights of the standard Gaussian integration rule
    */
    std::vector<double> get_Gauss_weights();

    /*!
    \brief Get the location of the standard Gaussian integration rule
    */
    std::vector<double> get_Gauss_line_pts();

    /*!
    \brief Get the integration rule for the considered line from Gaussian quadrature in the interval (-1,1)
    */
    LINALG::Matrix<2,8> find_line_integration_pts();

    /*!
    \brief The Jacobian of transformation of weights from standard quadrature defined in (-1,1) to the considered interval
    */
    double half_length();

    std::vector<double> point_begin_;
    std::vector<double> point_end_;
    int inte_num_;
    std::vector<double> alpha_;
    bool bcellInt_;
    std::string intType_;
                
};

#endif
