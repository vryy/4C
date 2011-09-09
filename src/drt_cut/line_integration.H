#ifndef LINE_INTEGRATION_H
#define LINE_INTEGRATION_H

#include<vector>
namespace GEO
{
        namespace CUT
        {
        class LineIntegration
        {
                public:
                        LineIntegration(std::vector<double> coord1, std::vector<double> coord2, int inte_num, std::vector<double> alpha)
                                :point_begin_(coord1),
                                 point_end_(coord2),
                                 inte_num_(inte_num),
                                 alpha_(alpha)
                        {
                        }
                        ~LineIntegration(){};
                        double integrate_line();

                private:
                        std::vector<double> point_begin_;
                        std::vector<double> point_end_;
                        int inte_num_;
                        std::vector<double> alpha_;
                        std::vector<double> compute_normal();
                        std::vector<double> get_Gauss_weights();
                        std::vector<double> get_Gauss_line_pts();
                        std::vector<std::vector<double> > find_line_integration_pts();
                        double half_length();
                
        };
        }
}

#endif
