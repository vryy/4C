/*----------------------------------------------------------------------*/
/*! \file

\brief base class for geometry pair evaluation data subcontainers.

\level 1
\maintainer Ivo Steinbrecher
*/
// End doxygen header.


#ifndef GEOMETRY_PAIR_EVALUATION_DATA_BASE_H
#define GEOMETRY_PAIR_EVALUATION_DATA_BASE_H


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
     */
    GeometryEvaluationDataBase(){};

    /**
     * \brief Destructor.
     */
    virtual ~GeometryEvaluationDataBase(){};

    /**
     * \brief Reset data in this evaluation container.
     */
    virtual void Reset(){};
  };
}  // namespace GEOMETRYPAIR

#endif
