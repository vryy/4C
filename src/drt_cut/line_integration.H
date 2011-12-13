#ifndef LINE_INTEGRATION_H
#define LINE_INTEGRATION_H

#include<vector>
#include "../linalg/linalg_fixedsizematrix.H"

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
    double integrate_line();
    void set_integ_type(std::string inttype){intType_ = inttype;}

  private:
    std::vector<double> point_begin_;
    std::vector<double> point_end_;
    int inte_num_;
    std::vector<double> alpha_;
    bool bcellInt_;
    std::string intType_;
    LINALG::Matrix<2,1> compute_normal();
    std::vector<double> get_Gauss_weights();
    std::vector<double> get_Gauss_line_pts();
    LINALG::Matrix<2,8> find_line_integration_pts();
    double half_length();
                
};

#endif
