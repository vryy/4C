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

#ifdef HAVE_MueLu
// MueLu
#include <MueLu_VerbosityLevel.hpp>
#endif

// Teuchos
#include <Teuchos_RCP.hpp>
#include <Teuchos_VerboseObject.hpp>

// baci
#include "../drt_lib/drt_dserror.H"

/*----------------------------------------------------------------------------*/
/* forward declarations */

// Teuchos
namespace Teuchos
{
  class ParameterList;
}  // namespace Teuchos

/*----------------------------------------------------------------------------*/
/* definition of methods */
namespace NLNSOL
{
  //! NLNSOL::UTILS: Utilities for nonlinear solver NLNSOL
  namespace UTILS
  {
    /*! \brief Store and manage access to configuration of NLNSOL parameter list
     *
     *  Input to the NLNSOL framework is done via an xml-file. The path to this
     *  file is specified in the dat-file parameter --NONLINEAR SOLVER/XML_FILE.
     *
     *  The xml-file has just a one-level hierarchy, i.e. the global list contains
     *  only one layer of sublists. Each sublist describes a single algorithmic
     *  part/object. In case, that further (sublist) specifications are required
     *  for nested algorithms/objects, we set at std::string parameter with the
     *  name of the sublist that contains these additional parameters.
     *
     *  \author mayr.mt \date 10/2015
     */
    class NlnConfig : public Teuchos::VerboseObject<NlnConfig>
    {
     public:
      //! @name Construction / Destruction
      //@{

      //! constructor
      NlnConfig();

      //! destructor
      virtual ~NlnConfig(){};

      //@}

      //! @name Setup
      //@{

      /*! \brief Setup the configuration object
       *
       *  Read a user-given xml-file into the Teuchos::ParameterList #params_.
       *
       *  ToDo (mayr) Add some sense of parameter validation and input checking.
       *
       *  \author mayr.mt \date 10/2015
       */
      virtual bool Setup(const std::string filename  ///< file name of xml-file
      );

      //@}

      //! @name Access to parameter lists
      //@{

      /*! \brief Access to parameter (sub-)list \c listname
       *
       *  This is the recommended way of accessing a sublist, since it protects
       *  the user-given configuration from changes during runtime.
       *
       *  \author mayr.mt \date 10/2015
       */
      virtual const Teuchos::ParameterList& GetSubList(
          const std::string listname  ///< name of sublist
          ) const;

      /*! \brief Access to parameter (sub-)list \c listname (non-const)
       *
       *  \note Returning a non-const list is necessary for certain NOX packages.
       *  If these issues are resolved, this routine should go away.
       *
       *  \warning Do not use this routine since it allows for modifications of
       *  the requested list during runtime which might be dangerous.
       *
       *  \author mayr.mt \date 10/2015
       */
      virtual Teuchos::ParameterList& GetSubListNonConst(
          const std::string listname  ///< name of sublist
          ) const;

      /*! \brief Access to parameter \c param in sublist \c list
       *
       *  First, we check whether the requested parameter \c param in the sublist
       *  \c list is available. Otherwise, we throw an error.
       *
       *  If available, the requested parameter is finally returned.
       *
       *  ToDo (mayr) Return default value, if parameter does not exist.
       *
       *  \author mayr.mt \date 10/2015
       */
      template <typename ParamType>
      ParamType GetParameter(const std::string list,  ///< name of sublist
          const std::string param                     ///< name of parameter
          ) const
      {
        if (not params_->isSublist(list)) dserror("Sublist '%s' does not exist.", list.c_str());

        const Teuchos::ParameterList& params = params_->sublist(list);

        return GetParameter<ParamType>(params, param);
      };

      /*! \brief Access to parameter \c param in top level parameter list
       *
       *  First, we check, whether the requested parameter \c param in the top
       *  level parameter list #params_ is available. Otherwise, we throw an error.
       *
       *  If available, the requested parameter is finally returned.
       *
       *  ToDo (mayr) Return default value, if parameter does not exist.
       *
       *  \author mayr.mt \date 10/2015
       */
      template <typename ParamType>
      ParamType GetParameter(const std::string param  ///< name of parameter
          ) const
      {
        return GetParameter<ParamType>(*params_, param);
      };

      /*! \brief Access to parameter \c param in given parameter list \c list
       *
       *  First, we check, whether the requested parameter \c param in the
       *  parameter list \c list is available. Otherwise, we throw an error.
       *
       *  If available, the requested parameter is finally returned.
       *
       *  ToDo (mayr) Return default value, if parameter does not exist.
       *
       *  \author mayr.mt \date 10/2015
       */
      template <typename ParamType>
      ParamType GetParameter(const Teuchos::ParameterList& list,  ///< sublist
          const std::string param                                 ///< name of parameter
          ) const
      {
        if (not list.isParameter(param.c_str()))
          dserror("Parameter '%s' does not exist in parameter list '%s'.", param.c_str(),
              list.name().c_str());

        return list.get<ParamType>(param);
      };

      /*! \brief Access to global parameter list
       *
       *  \note Returning a non-const list is necessary for certain NOX packages.
       *  If these issues are resolved, this routine should go away.
       *
       *  \warning Do not use this routine since it allows for modifications of
       *  the requested list during runtime which might be dangerous.
       *
       *  \author mayr.mt \date 10/2015
       */
      virtual Teuchos::RCP<Teuchos::ParameterList> GetAllNonConstRcp() const;

      //@}

     private:
      //! global parameter list
      Teuchos::RCP<Teuchos::ParameterList> params_;
    };

    // ---------------------------------------------------------------------------

    /*! \brief Detect stagnation during iterations
     *
     *  We monitor the reduction in the residual norm to heuristically detect
     *  stagnation during the nonlinear iterations. An iterative process is
     *  considered in stagnation if the residual reduction in #stagitermax_
     *  subsequent iterations is below #stagthreshold_.
     *
     *  \author mayr.mt \date 03/2015
     */
    class StagnationDetection : public Teuchos::VerboseObject<StagnationDetection>
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
      virtual void Init(Teuchos::RCP<const NLNSOL::UTILS::NlnConfig>
                            config,    ///< configuration of nonlinear solver
          const std::string listname,  ///< name of sublist with configuration of this operator
          const double norminitial     ///< initial residual norm
      );

      //@}

      //! @name Check for stagnation
      //@{

      /*! \brief Check for stagnation
       *
       *  Compute ratio of two subsequent residual norms and decide whether this
       *  decrease is considered as acceptable or as stagnation.
       *
       *  Update #normprev_ with current residual norm to make it available for
       *  the stagnation test in the next iteration.
       *
       *  \return Bool that indicates stagnation (true) or sufficient residual
       *  norm decrease (false)
       *
       *  \author mayr.mt \date 03/2015
       */
      virtual bool Check(const double norm);

      /*! \brief Get status of stagnation test
       *
       *  An iterative process is considered to be in stagnation if
       *  <ul>
       *  <li> the ratio of #stagitermax_ subsequent residuals is above a given
       *       threshold #stagthreshold_ </li>
       *  <li> the residual increases </li>
       *  </ul>
       *
       *  \return Bool that indicates stagnation (true) or sufficient residual
       *  norm decrease (false)
       *
       *  \author mayr.mt \date 03/2015
       */
      virtual bool Status(Teuchos::RCP<Teuchos::ParameterList> oparams =
                              Teuchos::null  ///< parameter list to provide detailed output
          ) const;

      /*! \brief Get parameter list with status of stagnation test
       *
       *  Create new Teuchos::RCP<const Teuchos::ParameterList> and pass it to
       *  Status() to be filled. Then return.
       *
       *  \return Parameter list describing the status of the stagnation test
       *
       *  \author mayr.mt \date 03/2015
       */
      virtual Teuchos::RCP<const Teuchos::ParameterList> StatusParams() const;

      //@}

      //! @name Access routines
      //@{

      //! Has stagnation detection been activated?
      virtual bool IsActive() const { return active_; }

      //@}

     protected:
     private:
      //! Access to members
      //@{

      //! Has Init() already been called?
      virtual bool IsInit() const { return isinit_; }

      //@}

      //! Has Init() been called?
      bool isinit_;

      //! Teuchos::ParameterList for configuration
      Teuchos::RCP<const NLNSOL::UTILS::NlnConfig> config_;

      //! number of subsequent iterations with stagnation
      int stagiter_;

      //! maximally allowed number of subsequent iterations with stagnation
      int stagitermax_;

      //! ratio of residual norms of subsequent iterations
      double stagratio_;

      //! threshold to declare stagnation
      double stagthreshold_;

      //! residual norm of previous iteration
      double normprev_;

      //! Has stagation detection been activated by the user?
      bool active_;

    };  // class StagnationDetection

    // ---------------------------------------------------------------------------

    //! @name Input from XML-File
    //@{

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
     *  The parameter list is read from a file whose filename is given as input
     *  to this routine.
     *
     *  \author mayr.mt \date 06/2014
     */
    Teuchos::RCP<Teuchos::ParameterList> CreateParamListFromXML(
        const std::string filename  ///< filename of xml file to be read
    );

    /*! \brief Create Teuchos::ParameterList from XML-File
     *
     *  The parameter list is read into parameters list 'params' from a file
     *  whose filename is given as input to this routine.
     *
     *  We follow the 'advanced XML style' in MueLu.
     *
     *  \author mayr.mt \date 10/2015
     */
    void CreateParamListFromXML(const std::string filename,  ///< filename of xml file to be read
        Teuchos::ParameterList& params                       ///< parameter list to be filled
    );

    /*! \brief Create Teuchos::ParameterList from XML-File
     *
     *  The parameter list is read into parameters list 'params' from a file
     *  whose filename is given as input to this routine.
     *
     *  We allow that sublists are read from separate xml-files. So, we need to
     *  loop over the read parameter list and check, whether we need to read
     *  additional sublists from another xml-file. This is done recursively.
     *
     *  \note: We introduce the convention, that sublists, which shall be read from
     *  a file, a marked by an std::string parameter entry
     *  "sublistfile: mypath/myfilename.xml". This helps to distinguish between
     *  standard std::string parameters and those, that replace a sublist.
     *
     *  \author mayr.mt \date 01/2015
     */
    void CreateParamListFromXMLOld(const std::string filename,  ///< filename of xml file to be read
        Teuchos::ParameterList& params                          ///< parameter list to be filled
    );

    //@}

    /*! \brief Translate verbosity level from std::string to
     *  Teuchos::EVerbosityLevel
     *
     *  Verbosity level is given as std::string in xml-input. We need to translate
     *  this to Teuchos::EVerbosityLevel to feed it to the Teuchos::VerboseObject.
     *
     *  Basically, we allow for all verbosity levels that are known in
     *  Teuchos::VerboseObject, though we use the std::string notation used in
     *  MueLu.
     *
     *  \author mayr.mt \date 03/2015
     */
    Teuchos::EVerbosityLevel TranslateVerbosityLevelToTeuchos(
        const std::string verblevel = "default"  ///< std::string describing the verosity level
    );

#ifdef HAVE_MueLu
    /*! \brief Translate verbosity level from std::string to
     *  MueLu::MsgType
     *
     *  Verbosity level is given as std::string in xml-input. We need to translate
     *  this to MueLu::MsgType to feed it to MueLu.
     *
     *  \author mayr.mt \date 10/2015
     */
    MueLu::MsgType TranslateVerbosityLevelToMueLu(
        const std::string verblevel = "default"  ///< std::string describing the verosity level
    );
#endif

    /*! \brief Status of ApplyInverse()
     *
     *  'opstatus_converged' on position 0 to be consistent with error code of
     *  Epetra::Operator::ApplyInverse()
     *
     *  \author mayr.mt \date 03/2015
     */
    enum OperatorStatus
    {
      opstatus_converged,    ///< converged successfully
      opstatus_undefined,    ///< undefined status
      opstatus_unconverged,  ///< did not converge in given number of iterations
      opstatus_stagnation,   ///< reached stagnation during iterations
      opstatus_failed        ///< failure due to non-specified error
    };

  }  // namespace UTILS
}  // namespace NLNSOL

#endif /* NLN_UTILS_H */
