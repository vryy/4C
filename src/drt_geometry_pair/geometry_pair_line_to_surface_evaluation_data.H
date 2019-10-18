/*----------------------------------------------------------------------*/
/*! \file

\brief Container for parameters for line to surface pairs, as well as global evaluation data.

\level 1
\maintainer Ivo Steinbrecher
*/
// End doxygen header.


#ifndef GEOMETRY_PAIR_LINE_TO_SURFACE_EVALUATION_DATA_GLOBAL_H
#define GEOMETRY_PAIR_LINE_TO_SURFACE_EVALUATION_DATA_GLOBAL_H


#include "geometry_pair_line_to_3D_evaluation_data.H"

#include <Teuchos_RCP.hpp>
#include <unordered_map>


namespace GEOMETRYPAIR
{
  /**
   * \brief Class to manage input parameters and evaluation data for line to surface interactions.
   */
  template <typename scalar_type>
  class LineToSurfaceEvaluationData : public LineTo3DEvaluationData
  {
   public:
    /**
     * \brief Constructor (derived).
     */
    LineToSurfaceEvaluationData(const Teuchos::ParameterList& input_parameter_list)
        : LineTo3DEvaluationData(input_parameter_list), node_normal_map_(){};

    /**
     * \brief Destructor.
     */
    ~LineToSurfaceEvaluationData() {}

    /**
     * \brief Reset the evaluation data (derived).
     */
    void Reset() override;

    /**
     * \brief Return a mutable reference to the node map.
     */
    inline std::unordered_map<int, std::vector<Teuchos::RCP<NodalNormal<scalar_type>>>>&
    GetNodeNormalMapMutable()
    {
      return node_normal_map_;
    }

   private:
    //! The keys to this map are the global node IDs and the values are the (averaged) normals of
    //! all connected elements at that node.
    std::unordered_map<int, std::vector<Teuchos::RCP<NodalNormal<scalar_type>>>> node_normal_map_;
  };
}  // namespace GEOMETRYPAIR

#endif
