/*!
\file geometry_pair.H

\brief base class for two or more elements that interact geometrically, e.g. contact or meshtying
problems.

\level 1
\maintainer Ivo Steinbrecher
*/
// End doxygen header.


#ifndef GEOMETRY_PAIR_H
#define GEOMETRY_PAIR_H


#include "../drt_lib/drt_dserror.H"

#include <Teuchos_RCP.hpp>


// Forward declarations.
namespace DRT
{
  class Element;
}
namespace GEOMETRYPAIR
{
  class GeometryEvaluationDataGlobal;
}


namespace GEOMETRYPAIR
{
  /**
   * \brief Base class for all geometry pairs.
   */
  class GeometryPair
  {
   public:
    /**
     * \brief Constructor.
     */
    GeometryPair();

    /**
     * \brief Destructor.
     */
    virtual ~GeometryPair(){};

    /**
     * \brief Initialize the pair.
     * @param evaluation_data_ptr Pointer to the evaluation data.
     * @param element1 Pointer to the first element.
     * @param element2 Pointer to the second element.
     */
    virtual void Init(Teuchos::RCP<GEOMETRYPAIR::GeometryEvaluationDataGlobal> evaluation_data_ptr,
        const DRT::Element* element1, const DRT::Element* element2);

    /**
     * \brief Setup the pair.
     */
    virtual void Setup();

    /*!
    \brief Get the evaluation data container.
    */
    inline Teuchos::RCP<GEOMETRYPAIR::GeometryEvaluationDataGlobal> EvaluationData() const
    {
      return evaluation_data_;
    };


   protected:
    /**
     * \brief Check if Init was already called.
     */
    inline void CheckInit() const
    {
      if (!isinit_) dserror("Init() has not been called yet!");
    };

    /**
     * \brief Check if Init and Setup were already called.
     */
    inline void CheckInitSetup() const
    {
      if (!isinit_ || !issetup_) dserror("Init() and Setup() have not been called yet!");
    };

    /**
     * \brief  Pointer to the first element.
     */
    inline const DRT::Element* Element1() const { return element1_; };

    /**
     * \brief  Pointer to the second element.
     */
    inline const DRT::Element* Element2() const { return element2_; };


   protected:
    //! Indicates if the Init() function has been called.
    bool isinit_;

    //! Indicates if the Setup() function has been called.
    bool issetup_;


   private:
    //! Link to the geometry evaluation container.
    Teuchos::RCP<GEOMETRYPAIR::GeometryEvaluationDataGlobal> evaluation_data_;

    //! Link to the element containing the first geometry.
    const DRT::Element* element1_;

    //! Link to the element containing the second geometry.
    const DRT::Element* element2_;
  };
}  // namespace GEOMETRYPAIR

#endif
