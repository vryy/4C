/*!
\file geometry_pair_evaluation_data_base.H

\brief Data container holding pointers to all subcontainers that in turn hold
       all evaluation data specific to their problem type,
       that can not be stored pairwise.

<pre>
\level 3
\maintainer Ivo Steinbrecher
            ivo.steinbrecher@unibw.de
            +49 89 6004-4403
</pre>
*/


#ifndef GEOMETRY_PAIR_EVALUATION_DATA_GLOBAL_H
#define GEOMETRY_PAIR_EVALUATION_DATA_GLOBAL_H


#include <Teuchos_RCP.hpp>


// Forward declaration.
namespace GEOMETRYPAIR
{
  class GeometryEvaluationDataBase;
  class LineToVolumeEvaluationData;
}  // namespace GEOMETRYPAIR


namespace GEOMETRYPAIR
{
  /**
   * \brief Class to store all geometry evaluation data subcontainers.
   */
  class GeometryEvaluationDataGlobal
  {
   public:
    /**
     * \brief Empty constructor.
     */
    GeometryEvaluationDataGlobal();

    /**
     * \brief Destructor.
     */
    virtual ~GeometryEvaluationDataGlobal(){};

    /**
     * \brief Reset the subcontainers.
     */
    void ResetGeometryEvaluationDataGlobal();

    /**
     * \brief Create the line to volume sub container.
     */
    void BuildLineToVolumeEvaluationData();

    /**
     * \brief Return the pointer to the LineToVolumeEvaluationData object
     * @return Line to volume evaluation data container.
     */
    inline Teuchos::RCP<GEOMETRYPAIR::LineToVolumeEvaluationData> LineToVolumeEvaluationData() const
    {
      return line_to_volume_evaluation_data_;
    }

   private:
    //! Pointer to the evaluation data container class of line to volume interaction.
    Teuchos::RCP<GEOMETRYPAIR::LineToVolumeEvaluationData> line_to_volume_evaluation_data_;

    //! Vector with all evaluation data subcontainers.
    std::vector<Teuchos::RCP<GEOMETRYPAIR::GeometryEvaluationDataBase>> sub_container_list_;
  };
}  // namespace GEOMETRYPAIR

#endif
