/*----------------------------------------------------------------------*/
/*! \file

\brief base class for geometry pair evaluation data subcontainers.

\level 1
\maintainer Ivo Steinbrecher
*/
// End doxygen header.


#ifndef GEOMETRY_PAIR_EVALUATION_DATA_BASE_H
#define GEOMETRY_PAIR_EVALUATION_DATA_BASE_H


#include <Teuchos_ParameterList.hpp>


namespace GEOMETRYPAIR
{
  /**
   * \brief A base class that all geometry pair evaluation data container have to inherit from.
   */
  class GeometryEvaluationDataBase
  {
   public:
    /**
     * \brief Constructor.
     *
     * \param input_parameter_list (in) Parameter list with the geometry evaluation parameters.
     */
    GeometryEvaluationDataBase(const Teuchos::ParameterList& input_parameter_list){};

    /**
     * \brief Destructor.
     */
    virtual ~GeometryEvaluationDataBase(){};

    /**
     * \brief Clear data in this evaluation container.
     */
    virtual void Clear(){};
  };
}  // namespace GEOMETRYPAIR

#endif
