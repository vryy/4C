/*----------------------------------------------------------------------------*/
/*!
\file nln_utils.H

\brief Collection of utilities for nonlinear solver objects

<pre>
Maintainer: Matthias Mayr
            mayr@mhpc.mw.tum.de
            089 - 289-10362
</pre>
*/

/*----------------------------------------------------------------------------*/

#ifndef NLN_UTILS_H
#define NLN_UTILS_H

/*----------------------------------------------------------------------------*/
/* headers */

// Teuchos
#include <Teuchos_RCP.hpp>

/*----------------------------------------------------------------------------*/
/* forward declarations */

namespace Teuchos
{
  class ParameterList;
} // namespace Teuchos

/*----------------------------------------------------------------------------*/
/* definition of methods */
namespace NLNSOL
{
//! NLNSOL::UTILS: Utilities for nonlinear solver NLNSOL
namespace UTILS
{

  /*! \brief Create Teuchos::ParameterList from XML-File
   *
   *  The parameter list is read from a file that is specified in the input file
   *  in the section "--NONLINEAR SOLVER" in the parameter "XML_FILE"
   *
   *  \author mayr.mt \date 06/2014
   */
  Teuchos::RCP<Teuchos::ParameterList> CreateParamListFromXML();

  /*! \brief Create Teuchos::ParameterList from XML-File
   *
   *  The parameter list is read from a file whoose filename is given as input
   *  to this routine.
   *
   *  \author mayr.mt \date 06/2014
   */
  Teuchos::RCP<Teuchos::ParameterList> CreateParamListFromXML(
      const std::string filename ///< filename of xml file to be read
      );

  /*! \brief Create Teuchos::ParameterList from XML-File
   *
   *  The parameter list is read into parameters list 'params' from a file
   *  whoose filename is given as input to this routine.
   *
   *  We allow that sublists are read from separate xml-files. So, we need to
   *  loop over the read parameter list and check, whether we need to read
   *  additional sublists from another xml-file. This is done recursively.
   *
   *  \Note: We introduce the convention, that sublists, which shall be read from
   *  a file, a marked by an std::string parameter entry
   *  "sublistfile: mypath/myfilename.xml". This helps to distinguish between
   *  standard std::string parameters and those, that replace a sublist.
   *
   *  \author mayr.mt \date 01/2015
   */
  void CreateParamListFromXML(
      const std::string filename, ///< filename of xml file to be read
      Teuchos::ParameterList& params ///< parameter list to be filled
      );

} // namespace UTILS
} // namespace NLNSOL

#endif /* NLN_UTILS_H */
