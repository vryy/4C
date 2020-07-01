/*-----------------------------------------------------------------------*/
/*! \file
\brief Contact utility functions related to parallel runs


\level 1

*/
/*----------------------------------------------------------------------------*/

#ifndef SRC_DRT_CONTACT_CONTACT_UTILS_PARALLEL_H_
#define SRC_DRT_CONTACT_CONTACT_UTILS_PARALLEL_H_

namespace Teuchos
{
  class ParameterList;
}

namespace CONTACT
{
  namespace UTILS
  {
    /*!
    \brief Decide whether to use the new code path that performs ghosting in a safe way or not

    The new code path performing redistribution and ghosting in a safe way, i.e. such that ghosting
    is extended often and far enough, is not working for all contact scenarios, yet.
    Use this function to check, whether the scenario given in the input file can use the new path or
    has to stick to the old path (with bugs in the extension of the interface ghosting).

    @param[in] contactParams Parameter list with all contact-relevant input parameters
    @return True, if new path is chosen. False otherwise.
    */
    bool UseSafeRedistributeAndGhosting(const Teuchos::ParameterList& contactParams);

  }  // namespace UTILS
}  // namespace CONTACT

#endif /* SRC_DRT_CONTACT_CONTACT_UTILS_PARALLEL_H_ */
