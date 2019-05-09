/*----------------------------------------------------------------------*/
/*!
\level 2

\brief Input parameters for multi-level monte carlo

\maintainer Jonas Biehler
*/

/*----------------------------------------------------------------------*/



#include "drt_validparameters.H"
#include "inpar_mlmc.H"



namespace INPAR
{
  namespace MLMC
  {
    /*----------------------------------------------------------------------*/
    /*----------------------------------------------------------------------*/
    void SetValidRandomFieldParameters(Teuchos::ParameterList& list)
    {
      using Teuchos::setStringToIntegralParameter;
      using Teuchos::tuple;
      using namespace DRT::INPUT;

      // define some tuples that are often used to account for different writing of certain key
      // words
      Teuchos::Array<std::string> yesnotuple =
          tuple<std::string>("Yes", "No", "yes", "no", "YES", "NO");
      Teuchos::Array<int> yesnovalue = tuple<int>(true, false, true, false, true, false);

      // safety measure to prevent using random fields that were not properly setup or meant to be
      // used at all
      setStringToIntegralParameter<int>(
          "ACTIVE", "No", "Do we want to use this field ?  ", yesnotuple, yesnovalue, &list);

      // Parameters to simulate random fields
      IntParameter("RANDOM_FIELD_DIMENSION", 3, "Dimension of Random Field 2 or 3", &list);
      IntParameter("SIZE_PER_DIM", 512, "Number of points per dimension", &list);
      DoubleParameter("CORRLENGTH", 30, "Correlation length of Random Field", &list);
      IntParameter("NUM_COS_TERMS", 64, "Number of terms in geometric row ", &list);
      setStringToIntegralParameter<int>("SPECTRAL_MATCHING", "Yes",
          "Perform spectral matching with psd 1 yes ", yesnotuple, yesnovalue, &list);
      DoubleParameter("SIGMA", 1.0, "sigma of random field", &list);
      DoubleParameter("MEAN", 0.0, "Mean value of random field", &list);
      setStringToIntegralParameter<int>("CORRSTRUCT", "gaussian",
          "Correlation structure of random field", tuple<std::string>("Gaussian", "gaussian"),
          tuple<int>(INPAR::MLMC::corr_gaussian, INPAR::MLMC::corr_gaussian), &list);
      setStringToIntegralParameter<int>("MARGINALPDF", "gaussian",
          "Target marginal probability distribution function",
          tuple<std::string>("Gaussian", "gaussian", "Beta", "beta", "lognormal", "Lognormal"),
          tuple<int>(INPAR::MLMC::pdf_gaussian, INPAR::MLMC::pdf_gaussian, INPAR::MLMC::pdf_beta,
              INPAR::MLMC::pdf_beta, INPAR::MLMC::pdf_lognormal, INPAR::MLMC::pdf_lognormal),
          &list);
      DoubleParameter("NONGAUSSPARAM1", 0, "First parameter for non-gaussian pdf", &list);
      DoubleParameter("NONGAUSSPARAM2", 0, "Second parameter for non-gaussian pdf", &list);
      DoubleParameter("KAPPA_U", 6.283185307, "CUTOFF WAVE NUMBER FOR PSD", &list);
      setStringToIntegralParameter<int>("CALC_METHOD", "FFT",
          "Calculation method for the random field", tuple<std::string>("FFT", "COS", "FOURIER"),
          tuple<int>(INPAR::MLMC::calc_m_fft, INPAR::MLMC::calc_m_cos, INPAR::MLMC::calc_m_fourier),
          &list);

      DoubleParameter("PERIODICITY_FOURIER", 1.0,
          "Periodic Length of Random Field in case of Fourier Series Expansion", &list);
      // For testing
      setStringToIntegralParameter<int>("USEDETVALUE", "No",
          "Instead of doing proper MC simulation use DETVALUE for the stochastic parameter",
          yesnotuple, yesnovalue, &list);
      DoubleParameter("CONTBLENDVALUE", 1.5, "Use this values for parameter continuation", &list);

      IntParameter("FOURIER_TRUNCATION_THRESHOLD", 130,
          "Truncation threshold for fourier series expansion", &list);

      // define cutoff values
      setStringToIntegralParameter<int>("BOUNDED", "No",
          "Cutoff Randomfield to prevent unrealistically low or high values", yesnotuple,
          yesnovalue, &list);
      DoubleParameter("LOWERBOUND", 1.5, "Lower cutoff value", &list);
      DoubleParameter("UPPERBOUND", 1.5, "Uower cutoff value", &list);

      // do have a locally varying mean/median value for the field
      setStringToIntegralParameter<int>("SPATIAL_VARIABLE_MEDIAN", "No",
          "Field has spatially varying median value", yesnotuple, yesnovalue, &list);
    }


    /*----------------------------------------------------------------------*/
    /*----------------------------------------------------------------------*/
    void SetValidRandomVariableParameters(Teuchos::ParameterList& list)
    {
      using Teuchos::setStringToIntegralParameter;
      using Teuchos::tuple;
      using namespace DRT::INPUT;

      // define some tuples that are often used to account for different writing of certain key
      // words
      Teuchos::Array<std::string> yesnotuple =
          tuple<std::string>("Yes", "No", "yes", "no", "YES", "NO");
      Teuchos::Array<int> yesnovalue = tuple<int>(true, false, true, false, true, false);

      // safety measure to prevent using random variables that were not properly
      // setup or meant to be used at all
      setStringToIntegralParameter<int>(
          "ACTIVE", "No", "Do we want to use this variable ?  ", yesnotuple, yesnovalue, &list);

      // Parameters of random variable
      DoubleParameter(
          "PARAM_1", 0.0, "First parameter of variable pdf (e.g. mean if gaussian )", &list);
      DoubleParameter(
          "PARAM_2", 1.0, "Second parameter of variable pdf (e.g. sigma if gaussian )", &list);
      setStringToIntegralParameter<int>("PDF", "gaussian",
          "Target probability distribution function",
          tuple<std::string>("Gaussian", "gaussian", "Beta", "beta", "lognormal", "Lognormal"),
          tuple<int>(INPAR::MLMC::pdf_gaussian, INPAR::MLMC::pdf_gaussian, INPAR::MLMC::pdf_beta,
              INPAR::MLMC::pdf_beta, INPAR::MLMC::pdf_lognormal, INPAR::MLMC::pdf_lognormal),
          &list);
      DoubleParameter("CONTBLENDVALUE", 1.5, "Use this values for parameter continuation", &list);

      // define cutoff values
      setStringToIntegralParameter<int>("BOUNDED", "No",
          "Cutoff random variable to prevent unrealistically low or high values", yesnotuple,
          yesnovalue, &list);
      DoubleParameter("LOWERBOUND", 1.5, "Lower cutoff value", &list);
      DoubleParameter("UPPERBOUND", 1.5, "Upper cutoff value", &list);
    }


    void SetValidParameters(Teuchos::RCP<Teuchos::ParameterList> list)
    {
      using namespace DRT::INPUT;
      using Teuchos::setStringToIntegralParameter;
      using Teuchos::tuple;

      Teuchos::Array<std::string> yesnotuple =
          tuple<std::string>("Yes", "No", "yes", "no", "YES", "NO");
      Teuchos::Array<int> yesnovalue = tuple<int>(true, false, true, false, true, false);

      /* parameters for multi-level monte carlo */
      Teuchos::ParameterList& mlmcp = list->sublist("MULTI LEVEL MONTE CARLO", false, "");

      IntParameter("NUMRUNS", 200, "Number of Monte Carlo runs", &mlmcp);

      IntParameter(
          "START_RUN", 0, "Run to start calculating the difference to lower level", &mlmcp);

      IntParameter("INITRANDOMSEED", 1000, "Random seed for first Monte Carlo run", &mlmcp);


      setStringToIntegralParameter<int>("REDUCED_OUTPUT", "NO",
          "Write reduced Coarse Level Output, i.e. no mesh stresses, just disp", yesnotuple,
          yesnovalue, &mlmcp);

      IntParameter(
          "CONTNUMMAXTRIALS", 8, "Half stepsize CONTNUMMAXTRIALS times before giving up", &mlmcp);

      setStringToIntegralParameter<int>("FWDPROBLEM", "structure",
          "WHICH KIND OF FORWARD DO WE WANT", tuple<std::string>("structure", "red_airways"),
          tuple<int>(structure, red_airways), &mlmcp);

      setStringToIntegralParameter<int>("UQSTRATEGY", "MC_PLAIN", "WHICH UQ STRATEGY WILL BE USED",
          tuple<std::string>("MC_PLAIN", "mc_plain", "MC_PARAMETERCONTINUATION",
              "mc_parametercontinuation", "MC_SCALEDTHICKNESS", "mc_scaledthickness"),
          tuple<int>(mc_plain, mc_plain, mc_paracont, mc_paracont, mc_scaledthick, mc_scaledthick),
          &mlmcp);

      IntParameter("NUMCONTSTEPS", 2, "Number of continuation steps", &mlmcp);

      // list of materials with stochastic constitutive parameters and
      // the stochastic parameters for the respective material
      StringParameter("PARAMLIST_R_FIELD", "none",
          "list of std::string of parameters that are to be modelled as random fields, 1 YOUNG "
          "BETA",
          &mlmcp);

      // list of materials with stochastic constitutive parameters and
      // the stochastic parameters for the respective material
      StringParameter("PARAMLIST_R_VAR", "none",
          "list of std::string of parameters that are to be modelled as random fields, 1 YOUNG "
          "BETA",
          &mlmcp);

      setNumericStringParameter("OUTPUT_ELEMENT_IDS", "-1",
          "Set ID's of Output Elements, default is -1 which is none", &mlmcp);

      setStringToIntegralParameter<int>("WRITE_RV_TO_FILE", "NO",
          "Write random variables used for random field generation to file", yesnotuple, yesnovalue,
          &mlmcp);

      // For variable geometry/wall thickness
      setStringToIntegralParameter<int>("RANDOMGEOMETRY", "No",
          "Do consider random geometry/wall thickness", yesnotuple, yesnovalue, &mlmcp);

      IntParameter("NUMALESTEPS", 1,
          "How many ALE steps do we want to use to compute uncertain geometry", &mlmcp);

      DoubleParameter("INITIALTHICKNESS", 10., "wall thickness in input file", &mlmcp);

      DoubleParameter("Z_POS_AAA_START_RF", -1078.,
          "Location of bifurcation of AAA including an offset ", &mlmcp);

      DoubleParameter(
          "TRANSITION_WIDTH", 15., "Transition domain to blend in random wall thickness", &mlmcp);
      // For variable geometry/wall thickness
      setStringToIntegralParameter<int>("START_RF_ABOVE_BIFURCATION", "No",
          "Start random geometry above bifurcation", yesnotuple, yesnovalue, &mlmcp);

      /*----------------------------------------------------------------------*/
      // set valid parameters for random fields

      // Note: the maximum number of random fields is hardwired here. If you change this,
      // don't forget to edit the corresponding parts in globalproblems.cpp, too.
      for (int i = 1; i < 4; i++)
      {
        std::stringstream ss;
        ss << "RANDOM FIELD " << i;
        std::stringstream ss_description;
        ss_description << "random field parameters for uncertainty quantification " << i;
        Teuchos::ParameterList& randomfieldlist =
            list->sublist(ss.str(), false, ss_description.str());
        SetValidRandomFieldParameters(randomfieldlist);
      }

      /*----------------------------------------------------------------------*/
      // set valid parameters for random fields

      // Note: the maximum number of random fields is hardwired here. If you change this,
      // don't forget to edit the corresponding parts in globalproblems.cpp, too.
      for (int i = 1; i < 4; i++)
      {
        std::stringstream ss;
        ss << "RANDOM VARIABLE " << i;
        std::stringstream ss_description;
        ss_description << "random variable parameters for uncertainty quantification " << i;
        Teuchos::ParameterList& randomvariablelist =
            list->sublist(ss.str(), false, ss_description.str());
        SetValidRandomVariableParameters(randomvariablelist);
      }
    }
  }  // namespace MLMC
}  // namespace INPAR
