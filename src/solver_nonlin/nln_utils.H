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
  /*! \brief Detect stagnation during iterations
   *
   *  We monitor the reduction in the residual norm to heuristically detect
   *  stagnation during the nonlinear iterations. An iterative process is
   *  considered in stagnation if the residual reduction in #stagitermax_
   *  subsequent iterations is below #stagthreshold_.
   *
   *  \author mayr.mt \date 03/2015
   */
  class StagnationDetection
  {
  public:

    //! @name Construction / Destruction
    //@{

    //! constructor
    StagnationDetection();

    //! destructor
    virtual ~StagnationDetection(){};

    //@}

    //! @name Setup
    //@{

    //! Initialization
    virtual void Init(
        const double norminitial ///< initial residual norm
        );

    //@}

    //! @name Check for stagnation
    //@{

    /*! \brief Check for stagnation
     *
     *  Compute ration of two subsequent residual norms and decide whether this
     *  decrease is considered as acceptable or as stagnation.
     *
     *  \return Bool that indicates stagnation (true) or sufficient residual
     *  norm decrease (false)
     *
     *  \author mayr.mt \date 03/2015
     */
    virtual const bool Check(const double norm);

    /*! \brief Get status of stagnation test
     *
     *  \return Bool that indicates stagnation (true) or sufficient residual
     *  norm decrease (false)
     *
     *  \author mayr.mt \date 03/2015
     */
    virtual const bool Status() const;

    //@}

  protected:

  private:

    //! Access to members
    //@{

    //! Has Init() already been called?
    virtual const bool IsInit() const { return isinit_; }

    //@}

    //! Has Init() been called?
    bool isinit_;

    //! number of subsequent iterations with stagnation
    int stagiter_;

    //! maximally allowed number of subsequent iterations with stagnation
    int stagitermax_;

    //! threshold to declare stagnation
    double stagthreshold_;

    //! residual norm of previous iteration
    double normprev_;


  }; // class StagnationDetection

  // ---------------------------------------------------------------------------

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

  /*! Status of ApplyInverse()
   *
   *  'opstatus_converged' on position 0 to be consistent with error code of
   *  Epetra::Operator::ApplyInverse()
   *
   *  \author mayr.mt \date 03/2015
   */
  enum OperatorStatus
  {
    opstatus_converged, ///< converged successfully
    opstatus_undefined, ///< undefined status
    opstatus_unconverged, ///< did not converge in given number of iterations
    opstatus_stagnation, ///< reached stagnation during iterations
    opstatus_failed ///< failure due to non-specified error
  };

} // namespace UTILS
} // namespace NLNSOL

#endif /* NLN_UTILS_H */
