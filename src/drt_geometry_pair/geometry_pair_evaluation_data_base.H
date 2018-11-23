/*!
\file geometry_pair_evaluation_data_base.H

\brief base class for geometry pair evaluation data subcontainers.

<pre>
\level 3
\maintainer Ivo Steinbrecher
            ivo.steinbrecher@unibw.de
            +49 89 6004-4403
</pre>
*/


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
    GeometryEvaluationDataBase();

    /**
     * \brief Destructor.
     */
    virtual ~GeometryEvaluationDataBase(){};

    /**
     * \brief Initialize the container.
     */
    virtual void Init();

    /**
     * \brief Setup the container.
     */
    virtual void Setup();


   protected:
    /**
     * \brief Check if Init was already called.
     */
    void CheckInit() const;

    /**
     * \brief Check if Init and Setup were already called.
     */
    void CheckInitSetup() const;


   protected:
    //! Indicates if the Init() function has been called.
    bool isinit_;

    //! Indicates if the Setup() function has been called.
    bool issetup_;
  };
}  // namespace GEOMETRYPAIR

#endif
