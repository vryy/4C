/*----------------------------------------------------------------------*/
/*!
\file drt_function.cpp

\brief Managing and evaluating of spatial functions

<pre>
\level 1

\maintainer Martin Kronbichler
            http://www.lnm.mw.tum.de
            089 - 289-15235
</pre>
*/
/*----------------------------------------------------------------------*/


#include <Sacado.hpp>
#include "drt_function.H"
#include "drt_functionvariables.H"
#include "drt_globalproblem.H"
#include "standardtypes_cpp.H"
#include "drt_linedefinition.H"
#include "../drt_combust/combust_functions.H"
#include "../drt_fluid/fluid_functions.H"
#include "../drt_fluid_xfluid/xfluid_functions.H"
#include "../drt_io/io.H"

/*----------------------------------------------------------------------*/
//! Print function
/*----------------------------------------------------------------------*/
void PrintFunctionDatHeader()
{
  DRT::UTILS::FunctionManager functionmanager;
  Teuchos::RCP<DRT::INPUT::Lines> lines = functionmanager.ValidFunctionLines();

  lines->Print(std::cout);
}


Teuchos::RCP<DRT::INPUT::Lines> DRT::UTILS::FunctionManager::ValidFunctionLines()
{

  DRT::INPUT::LineDefinition onecomponentexpr;
  onecomponentexpr
  .AddNamedString("FUNCTION")
  ;

  DRT::INPUT::LineDefinition componentexpr;
  componentexpr
  .AddNamedInt("COMPONENT")
  .AddNamedString("FUNCTION")
  ;

  DRT::INPUT::LineDefinition variableexpr;
  variableexpr
  .AddNamedInt("VARIABLE")
  .AddNamedString("NAME")
  .AddNamedString("TYPE")
  .AddOptionalNamedString("DESCRIPTION")
  .AddOptionalNamedInt("NUMPOINTS")
  .AddOptionalNamedString("BYNUM")
  .AddOptionalNamedDoubleVector("TIMERANGE",2)
  .AddOptionalNamedDoubleVector("TIMES","NUMPOINTS")
  .AddOptionalNamedDoubleVector("VALUES","NUMPOINTS")
  .AddOptionalNamedString("PERIODIC")
  .AddOptionalNamedDouble("T1")
  .AddOptionalNamedDouble("T2")
  ;

  DRT::INPUT::LineDefinition variableexprmulti;
  variableexprmulti
  .AddNamedInt("VARIABLE")
  .AddNamedString("NAME")
  .AddNamedString("TYPE")
  .AddOptionalNamedInt("NUMPOINTS")
  .AddOptionalNamedString("BYNUM")
  .AddOptionalNamedDoubleVector("TIMERANGE",2)
  .AddOptionalNamedDoubleVector("TIMES","NUMPOINTS")
  .AddOptionalNamedDoubleVector("VALUES","NUMPOINTS")
  .AddOptionalNamedStringVector("DESCRIPTION","NUMPOINTS") // only NUMPOINTS-1 are taken
  .AddOptionalNamedString("PERIODIC")
  .AddOptionalNamedDouble("T1")
  .AddOptionalNamedDouble("T2")
  ;

  Teuchos::RCP<DRT::INPUT::Lines> lines = Teuchos::rcp(new DRT::INPUT::Lines("FUNCT"));
  lines->Add(onecomponentexpr);
  lines->Add(componentexpr);
  lines->Add(variableexpr);
  lines->Add(variableexprmulti);

  // Old function +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  DRT::INPUT::LineDefinition componentvarexpr;
  componentvarexpr
    .AddNamedInt("COMPONENT")
    .AddNamedDoubleVector("VAREXPR",3)
    .AddNamedString("FUNCTION")
    .AddOptionalNamedInt("NUMCONSTANTS")
    .AddOptionalNamedPairOfStringAndDoubleVector("CONSTANTS","NUMCONSTANTS")
    ;

  DRT::INPUT::LineDefinition varexpr;
  varexpr
    .AddNamedDoubleVector("VAREXPR",3)
    .AddNamedString("FUNCTION")
    .AddOptionalNamedInt("NUMCONSTANTS")
    .AddOptionalNamedPairOfStringAndDoubleVector("CONSTANTS","NUMCONSTANTS")
    ;

  DRT::INPUT::LineDefinition linelin;
  linelin
  .AddNamedDoubleVector("LINE_LIN",8)
  ;

  DRT::INPUT::LineDefinition radiuslin;
  radiuslin
  .AddNamedDoubleVector("RADIUS_LIN",8)
  ;

  DRT::INPUT::LineDefinition radiusquad;
  radiusquad
  .AddNamedDoubleVector("RADIUS_QUAD",6)
  ;

  DRT::INPUT::LineDefinition beltrami;
  beltrami
  .AddTag("BELTRAMI")
  .AddNamedDouble("c1")
  ;

  DRT::INPUT::LineDefinition kimmoin;
  kimmoin
  .AddTag("KIM-MOIN")
  ;

  DRT::INPUT::LineDefinition bochevup;
  bochevup
  .AddTag("BOCHEV-UP")
  ;

  DRT::INPUT::LineDefinition bochevrhs;
  bochevrhs
  .AddTag("BOCHEV-RHS")
  ;

  DRT::INPUT::LineDefinition beltramiup;
  beltramiup
  .AddTag("BELTRAMI-UP")
  .AddNamedInt("MAT")
  .AddNamedInt("ISSTAT")
  ;

  DRT::INPUT::LineDefinition beltramigradu;
  beltramigradu
  .AddTag("BELTRAMI-GRADU")
  .AddNamedInt("MAT")
  .AddNamedInt("ISSTAT")
  ;

  DRT::INPUT::LineDefinition beltramirhs;
  beltramirhs
  .AddTag("BELTRAMI-RHS")
  .AddNamedInt("MAT")
  .AddNamedInt("ISSTAT")
  .AddNamedInt("ISSTOKES")
  ;

  DRT::INPUT::LineDefinition kimmoinup;
  kimmoinup
  .AddTag("KIMMOIN-UP")
  .AddNamedInt("MAT")
  .AddNamedInt("ISSTAT")
  ;

  DRT::INPUT::LineDefinition kimmoingradu;
  kimmoingradu
  .AddTag("KIMMOIN-GRADU")
  .AddNamedInt("MAT")
  .AddNamedInt("ISSTAT")
  ;

  DRT::INPUT::LineDefinition kimmoinrhs;
  kimmoinrhs
  .AddTag("KIMMOIN-RHS")
  .AddNamedInt("MAT")
  .AddNamedInt("ISSTAT")
  .AddNamedInt("ISSTOKES")
  ;

  DRT::INPUT::LineDefinition turbboulayer;
  turbboulayer
  .AddTag("TURBBOULAYER")
  ;

  DRT::INPUT::LineDefinition turbboulayerbfs;
  turbboulayerbfs
  .AddTag("TURBBOULAYER-BFS")
  ;

  DRT::INPUT::LineDefinition turbboulayeroracles;
  turbboulayeroracles
  .AddTag("TURBBOULAYERORACLES")
  ;

  DRT::INPUT::LineDefinition jefferyhamel;
  jefferyhamel
  .AddTag("JEFFERY-HAMEL")
  ;

  DRT::INPUT::LineDefinition womersley;
  womersley
  .AddTag("WOMERSLEY")
  .AddNamedInt("Local")
  .AddNamedInt("MAT")
  .AddNamedInt("CURVE")
  .AddNamedString("FSI")
  ;

  DRT::INPUT::LineDefinition localwomersley;
  localwomersley
  .AddTag("WOMERSLEY")
  .AddNamedInt("Local")
  .AddNamedDouble("Radius")
  .AddNamedInt("MAT")
  .AddNamedInt("CURVE")
  ;

  DRT::INPUT::LineDefinition cylinder3d;
  cylinder3d
  .AddNamedDouble("CYLINDER_3D")
  ;

  DRT::INPUT::LineDefinition controlledrotation;
  controlledrotation
  .AddTag("CONTROLLEDROTATION")
  .AddNamedString("FILE")
  .AddNamedString("TYPE")
  .AddNamedDoubleVector("ORIGIN",3)
  ;

  DRT::INPUT::LineDefinition accelerationprofile;
  accelerationprofile
  .AddTag("ACCELERATIONPROFILE")
  .AddNamedString("FILE")
  ;

  DRT::INPUT::LineDefinition ramptovalue;
  ramptovalue
  .AddTag("RAMPTOVALUE")
  .AddNamedDouble("VALUE")
  .AddNamedDouble("STARTTIME")
  .AddNamedDouble("DURATION")
  .AddNamedString("TYPE")
  ;

  DRT::INPUT::LineDefinition nodenormal;
  nodenormal
  .AddTag("NODENORMAL")
  .AddNamedString("GEOMETRY")
  .AddNamedDoubleVector("ORIGIN",3)
  .AddNamedDouble("RADIUS")
  .AddNamedDouble("CYLINDERHEIGHT")
  .AddNamedDoubleVector("ORIENTATION",3)
  .AddNamedDouble("CASSINIA")
  ;

  lines->Add(componentvarexpr);
  lines->Add(varexpr);
  lines->Add(linelin);
  lines->Add(radiuslin);
  lines->Add(radiusquad);
  lines->Add(beltrami);
  lines->Add(kimmoin);
  lines->Add(bochevup);
  lines->Add(bochevrhs);
  lines->Add(beltramiup);
  lines->Add(beltramigradu);
  lines->Add(beltramirhs);
  lines->Add(kimmoinup);
  lines->Add(kimmoingradu);
  lines->Add(kimmoinrhs);
  lines->Add(turbboulayer);
  lines->Add(turbboulayerbfs);
  lines->Add(turbboulayeroracles);
  lines->Add(jefferyhamel);
  lines->Add(womersley);
  lines->Add(localwomersley);
  lines->Add(cylinder3d);
  lines->Add(controlledrotation);
  lines->Add(accelerationprofile);
  lines->Add(ramptovalue);
  lines->Add(nodenormal);
  // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  COMBUST::CombustValidFunctionLines(lines);
  DRT::UTILS::XfluidValidFunctionLines(lines);

  return lines;
}


void DRT::UTILS::FunctionManager::ReadInput(DRT::INPUT::DatFileReader& reader)
{

  functions_.clear();

  Teuchos::RCP<DRT::INPUT::Lines> lines = ValidFunctionLines();

  // test for as many functions as there are
  for (int i = 1;; ++i)
  {
    std::vector<Teuchos::RCP<DRT::INPUT::LineDefinition> > functions = lines->Read(reader, i);

    if (functions.size()==0)
      break;

    else
    {
      Teuchos::RCP<DRT::INPUT::LineDefinition> function = functions[0];

      // Old function +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      if (function->HaveNamed("BELTRAMI"))
      {
        double c1;
        function->ExtractDouble("c1",c1);

        functions_.push_back(Teuchos::rcp(new FLD::BeltramiFunction(c1)));
      }
      else if (function->HaveNamed("KIM-MOIN"))
      {
        functions_.push_back(Teuchos::rcp(new FLD::KimMoinFunction()));
      }
      else if (function->HaveNamed("BOCHEV-UP"))
      {
        functions_.push_back(Teuchos::rcp(new FLD::BochevUPFunction()));
      }
      else if (function->HaveNamed("BOCHEV-RHS"))
      {
        functions_.push_back(Teuchos::rcp(new FLD::BochevRHSFunction()));
      }
      else if (function->HaveNamed("BELTRAMI-UP"))
      {
        // read material
        int mat_id = -1;

        function->ExtractInt("MAT",mat_id);

        if(mat_id<=0) dserror("Please give a (reasonable) 'MAT'/material in BELTRAMI-UP");

        functions_.push_back(Teuchos::rcp(new FLD::BeltramiUP(mat_id)));
      }
      else if (function->HaveNamed("BELTRAMI-GRADU"))
      {
        // read material
        int mat_id = -1;

        function->ExtractInt("MAT",mat_id);

        if(mat_id<=0) dserror("Please give a (reasonable) 'MAT'/material in BELTRAMI-GRADU");

        functions_.push_back(Teuchos::rcp(new FLD::BeltramiGradU(mat_id)));
      }
      else if (function->HaveNamed("BELTRAMI-RHS"))
      {
        // read material
        int mat_id = -1;
        int is_stokes     = 0;

        function->ExtractInt("MAT",mat_id);
        function->ExtractInt("ISSTOKES",is_stokes);

        if(mat_id<=0) dserror("Please give a (reasonable) 'MAT'/material in BELTRAMI-GRADU");

        functions_.push_back(Teuchos::rcp(new FLD::BeltramiRHS(mat_id, (bool)is_stokes)));
      }
      else if (function->HaveNamed("KIMMOIN-UP"))
      {
        // read material
        int mat_id = -1;
        int is_stationary = 0;

        function->ExtractInt("MAT",mat_id);
        function->ExtractInt("ISSTAT",is_stationary);

        if(mat_id<=0) dserror("Please give a (reasonable) 'MAT'/material in KIMMOIN-UP");

        functions_.push_back(Teuchos::rcp(new FLD::KimMoinUP(mat_id, (bool)is_stationary)));
      }
      else if (function->HaveNamed("KIMMOIN-GRADU"))
      {
        // read material
        int mat_id = -1;
        int is_stationary = 0;

        function->ExtractInt("MAT",mat_id);
        function->ExtractInt("ISSTAT",is_stationary);

        if(mat_id<=0) dserror("Please give a (reasonable) 'MAT'/material in KIMMOIN-GRADU");

        functions_.push_back(Teuchos::rcp(new FLD::KimMoinGradU(mat_id, (bool)is_stationary)));
      }
      else if (function->HaveNamed("KIMMOIN-RHS"))
      {
        // read material
        int mat_id = -1;
        int is_stationary = 0;
        int is_stokes     = 0;

        function->ExtractInt("MAT",mat_id);
        function->ExtractInt("ISSTAT",is_stationary);
        function->ExtractInt("ISSTOKES",is_stokes);

        if(mat_id<=0) dserror("Please give a (reasonable) 'MAT'/material in KIMMOIN-GRADU");

        functions_.push_back(Teuchos::rcp(new FLD::KimMoinRHS(mat_id, (bool)is_stationary, (bool)is_stokes)));
      }
      else if (function->HaveNamed("TURBBOULAYER"))
      {
        functions_.push_back(Teuchos::rcp(new FLD::TurbBouLayerFunction()));
      }
      else if (function->HaveNamed("TURBBOULAYER-BFS"))
      {
        functions_.push_back(Teuchos::rcp(new FLD::TurbBouLayerFunctionBFS()));
      }
      else if (function->HaveNamed("TURBBOULAYERORACLES"))
      {
        functions_.push_back(Teuchos::rcp(new FLD::TurbBouLayerFunctionORACLES()));
      }
      else if (function->HaveNamed("JEFFERY-HAMEL"))
      {
        functions_.push_back(Teuchos::rcp(new FLD::JefferyHamelFlowFunction()));
      }
      else if (function->HaveNamed("ZALESAKSDISK"))
      {
        functions_.push_back(Teuchos::rcp(new COMBUST::ZalesaksDiskFunction()));
      }
      else if (function->HaveNamed("CIRCULARFLAME2"))
      {
        functions_.push_back(Teuchos::rcp(new COMBUST::CircularFlame2Function()));
      }
      else if (function->HaveNamed("CIRCULARFLAME3"))
      {
        functions_.push_back(Teuchos::rcp(new COMBUST::CircularFlame3Function()));
      }
      else if (function->HaveNamed("CIRCULARFLAME4"))
      {
        functions_.push_back(Teuchos::rcp(new COMBUST::CircularFlame4Function()));
      }
      else if (function->HaveNamed("DAMBREAKOBSTACLE"))
      {
        functions_.push_back(Teuchos::rcp(new COMBUST::DamBreakObstacle()));
      }
      else if (function->HaveNamed("COLLAPSINGWATERCOLUMN"))
      {
        functions_.push_back(Teuchos::rcp(new COMBUST::CollapsingWaterColumnFunction()));
      }
      else if (function->HaveNamed("COLLAPSINGWATERCOLUMNCOARSE"))
      {
        functions_.push_back(Teuchos::rcp(new COMBUST::CollapsingWaterColumnFunctionCoarse()));
      }
      else if (function->HaveNamed("IMPACTDROP"))
      {
        functions_.push_back(Teuchos::rcp(new COMBUST::ImpactFunction()));
      }
      else if (function->HaveNamed("BUBBLES"))
      {
        functions_.push_back(Teuchos::rcp(new COMBUST::BubbleFunction()));
      }
      else if (function->HaveNamed("ORACLESGFUNC"))
      {
        functions_.push_back(Teuchos::rcp(new COMBUST::ORACLESGFunction()));
      }
      else if (function->HaveNamed("ROTATINGCONE"))
      {
        functions_.push_back(Teuchos::rcp(new COMBUST::RotatingConeFunction()));
      }
      else if (function->HaveNamed("LEVELSETCUTTEST"))
      {
        functions_.push_back(Teuchos::rcp(new COMBUST::LevelSetCutTestFunction()));
      }
      else if (function->HaveNamed("FORWARDFACINGSTEP"))
      {
        functions_.push_back(Teuchos::rcp(new GerstenbergerForwardfacingStep()));
      }
      else if (function->HaveNamed("SLIPLENGTHFUNCTION"))
      {
        functions_.push_back(Teuchos::rcp(new SlipLengthLevelSetManipulator()));
      }
      else if (function->HaveNamed("MOVINGLEVELSETCYLINDER"))
      {
        std::vector<double> origin;
        function->ExtractDoubleVector("ORIGIN",origin);

        double radius;
        function->ExtractDouble("RADIUS",radius);

        std::vector<double> direction;
        function->ExtractDoubleVector("DIRECTION",direction);

        double distance;
        function->ExtractDouble("DISTANCE",distance);

        double maxspeed;
        function->ExtractDouble("MAXSPEED",maxspeed);

        functions_.push_back(Teuchos::rcp(new MovingLevelSetCylinder( &origin, radius, &direction, distance, maxspeed)));
      }
      else if (function->HaveNamed("TAYLORCOUETTEFLOW"))
      {
        double radius_i;
        function->ExtractDouble("RADIUS_I",radius_i);
        double radius_o;
        function->ExtractDouble("RADIUS_O",radius_o);

        double vel_theta_i;
        function->ExtractDouble("VEL_THETA_I",vel_theta_i);
        double vel_theta_o;
        function->ExtractDouble("VEL_THETA_O",vel_theta_o);

        double sliplength_i;
        function->ExtractDouble("SLIPLENGTH_I",sliplength_i);
        double sliplength_o;
        function->ExtractDouble("SLIPLENGTH_O",sliplength_o);

        double traction_theta_i;
        function->ExtractDouble("TRACTION_THETA_I",traction_theta_i);
        double traction_theta_o;
        function->ExtractDouble("TRACTION_THETA_O",traction_theta_o);

        double viscosity;
        function->ExtractDouble("VISCOSITY",viscosity);

        functions_.push_back(Teuchos::rcp(new TaylorCouetteFlow( radius_i,
            radius_o,
            vel_theta_i,
            vel_theta_o,
            sliplength_i,
            sliplength_o,
            traction_theta_i,
            traction_theta_o,
            viscosity)));

      }
      else if (function->HaveNamed("CONTROLLEDROTATION"))
      {
        std::string fileName;
        function->ExtractString("FILE", fileName);

        std::string type;
        function->ExtractString("TYPE", type);

        std::vector<double> origin;
        function->ExtractDoubleVector("ORIGIN",origin);

        functions_.push_back(Teuchos::rcp(new ControlledRotationFunction(fileName, type, origin[0], origin[1], origin[2])));
      }
      else if (function->HaveNamed("ACCELERATIONPROFILE"))
      {
        std::string fileName;
        function->ExtractString("FILE", fileName);

        functions_.push_back(Teuchos::rcp(new AccelerationProfileFunction(fileName)));
      }
      else if (function->HaveNamed("VAREXPR"))
      {
        Teuchos::RCP<VariableExprFunction> vecfunc = Teuchos::rcp(new VariableExprFunction());

        std::vector<double> origin;
        function->ExtractDoubleVector("VAREXPR",origin);
        std::string component;
        function->ExtractString("FUNCTION",component);

        std::vector<std::pair<std::string,double> > constants;
        if (function->HaveNamed("CONSTANTS"))
          function->ExtractPairOfStringAndDoubleVector("CONSTANTS",constants);

        vecfunc->AddExpr(component,origin[0],origin[1],origin[2],constants);
        functions_.push_back(vecfunc);
      }
      else if (COMBUST::CombustFunctionHaveNamed(function, &functions_))
      {
      }
      else if (DRT::UTILS::XfluidFunctionHaveNamed(function, &functions_))
      {
      }
      else
      {
        // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

        // define a new vector of functions
        Teuchos::RCP<ExprFunction> vecfunc = Teuchos::rcp(new ExprFunction());

        // evaluate the maximum component and the number of variables
        int maxcomp = 0;
        int maxvar = -1;
        for (unsigned int n=0; n<functions.size(); ++n)
        {
          Teuchos::RCP<DRT::INPUT::LineDefinition> functnumcompvar = functions[n];
          functnumcompvar->ExtractInt("COMPONENT",maxcomp);
          functnumcompvar->ExtractInt("VARIABLE",maxvar);
        }

        // evaluate the number of rows used for the definition of the variables
        int numrowsvar;
        numrowsvar = functions.size() - maxcomp - 1;

        // define a vector of strings
        std::vector<std::string> functstring(maxcomp+1);

        // read each row where the components of the i-th function are defined
        for (int n = 0; n <= maxcomp; ++n)
        {
          // update the current row
          Teuchos::RCP<DRT::INPUT::LineDefinition> functcomp = functions[n];

          // check the validity of the n-th component
          int compid = 0;
          functcomp->ExtractInt("COMPONENT",compid);
          if (compid!=n) dserror("expected COMPONENT %d but got COMPONENT %d", n, compid);


          // read the expression of the n-th component of the i-th function
          functcomp->ExtractString("FUNCTION",functstring[n]);
        }

        // define the structure functvarvector
        std::vector<std::vector<Teuchos::RCP<FunctionVariable> > > functvarvector;

        // define the structure functvar
        std::vector<Teuchos::RCP<FunctionVariable> > functvar;

        // define the structure vardef
        Teuchos::RCP<FunctionVariable> vardef;

        int vardefinition = 1;
        int varidold = -1;

        // read each row where the variables of the i-th function are defined
        for (int j = 1;j <= numrowsvar; ++j)
        {
          // update the current row
          Teuchos::RCP<DRT::INPUT::LineDefinition> timevar = functions[maxcomp+j];

          // read the number of the variable
          int varid;
          timevar->ExtractInt("VARIABLE",varid);

          // evaluate the number of the definition for the variable
          if (varid == varidold)
          {
            ++vardefinition;
          }
          else
          {
            vardefinition = 1;
          }

          // update the old varid
          varidold = varid;

          // read the name of the variable
          std::string varname;
          timevar->ExtractString("NAME",varname);

          // read the type of the variable
          std::string vartype;
          timevar->ExtractString("TYPE",vartype);

          // read periodicity data
          struct periodicstruct periodicdata;
          periodicdata.periodic = timevar->FindString("PERIODIC");
          if (periodicdata.periodic)
          {
            timevar->ExtractDouble("T1",periodicdata.t1);
            timevar->ExtractDouble("T2",periodicdata.t2);
          }
          else
          {
            periodicdata.t1 = 0;
            periodicdata.t2 = 0;
          }

          // distinguish the type of the variable
          if (vartype=="expression")
          {
            std::string description;
            timevar->ExtractString("DESCRIPTION",description);
            vardef = Teuchos::rcp(new ParsedFunctionVariable(varname,description));
          }
          else if(vartype=="linearinterpolation")
          {
            // read the number of points
            int numpoints;
            timevar->ExtractInt("NUMPOINTS",numpoints);

            // read times
            std::vector<double> times;
            bool bynum = timevar->FindString("BYNUM");

            if (bynum) // times defined by number of points
            {
              // read the time range
              std::vector<double> timerange;
              timevar->ExtractDoubleVector("TIMERANGE",timerange);

              // get initial and final time
              double t_initial = timerange[0];
              double t_final = timerange[1];

              // build the vector of times
              times.push_back(t_initial);
              int n = 0;
              double dt = (t_final - t_initial) / (numpoints - 1);
              while(times[n] + dt <= t_final+1.0e-14)
              {
                if (times[n] + 2*dt <= t_final+1.0e-14)
                {
                  times.push_back(times[n] + dt);
                }
                else
                {
                  times.push_back(t_final);
                }
                ++n;
              }
            }
            else // times defined by vector
            {
              timevar->ExtractDoubleVector("TIMES",times);
            }

            // check if the times are in ascending order
            for (unsigned int k=1; k<times.size(); ++k)
            {
              if (times[k] <= times[k-1])
                dserror("the TIMES must be in ascending order");
            }

            // read values
            std::vector<double> values;
            timevar->ExtractDoubleVector("VALUES",values);

            vardef = Teuchos::rcp(new LinearInterpolationVariable(varname,times,values,periodicdata));
          }
          else if(vartype=="multifunction")
          {
            // read the number of points
            int numpoints;
            timevar->ExtractInt("NUMPOINTS",numpoints);

            // read times
            std::vector<double> times;
            bool bynum = timevar->FindString("BYNUM");

            if (bynum) // times defined by number of points
            {
              // read the time range
              std::vector<double> timerange;
              timevar->ExtractDoubleVector("TIMERANGE",timerange);

              // get initial and final time
              double t_initial = timerange[0];
              double t_final = timerange[1];

              // build the vector of times
              times.push_back(t_initial);
              int n = 0;
              double dt = (t_final - t_initial) / (numpoints - 1);
              while(times[n] + dt <= t_final+1.0e-14)
              {
                if (times[n] + 2*dt <= t_final+1.0e-14)
                {
                  times.push_back(times[n] + dt);
                }
                else
                {
                  times.push_back(t_final);
                }
                ++n;
              }
            }
            else // times defined by vector
            {
              timevar->ExtractDoubleVector("TIMES",times);
            }

            // check if the times are in ascending order
            for (unsigned int k=1; k<times.size(); ++k)
            {
              if (times[k] <= times[k-1])
                dserror("the TIMES must be in ascending order");
            }

            // read descriptions (strings separated with spaces)
            std::vector<std::string> description_vec;
            timevar->ExtractStringVector("DESCRIPTION",description_vec);

            // check if the number of times = number of descriptions + 1
            int numtimes = times.size();
            int numdescriptions = description_vec.size();
            if (numtimes!=numdescriptions+1) dserror("the number of TIMES and the number of DESCRIPTIONs must be consistent");

            vardef = Teuchos::rcp(new MultiFunctionVariable(varname,times,description_vec,periodicdata));
          }
          else if(vartype=="fourierinterpolation")
          {
            // read the number of points
            int numpoints;
            timevar->ExtractInt("NUMPOINTS",numpoints);

            // read times
            std::vector<double> times;
            bool bynum = timevar->FindString("BYNUM");

            if (bynum) // times defined by number of points
            {
              // read the time range
              std::vector<double> timerange;
              timevar->ExtractDoubleVector("TIMERANGE",timerange);

              // get initial and final time
              double t_initial = timerange[0];
              double t_final = timerange[1];

              // build the vector of times
              times.push_back(t_initial);
              int n = 0;
              double dt = (t_final - t_initial) / (numpoints - 1);
              while(times[n] + dt <= t_final+1.0e-14)
              {
                if (times[n] + 2*dt <= t_final+1.0e-14)
                {
                  times.push_back(times[n] + dt);
                }
                else
                {
                  times.push_back(t_final);
                }
                ++n;
              }
            }
            else // times defined by vector
            {
              timevar->ExtractDoubleVector("TIMES",times);
            }

            // check if the times are in ascending order
            for (unsigned int k=1; k<times.size(); ++k)
            {
              if (times[k] <= times[k-1])
                dserror("the TIMES must be in ascending order");
            }

            // read values
            std::vector<double> values;
            timevar->ExtractDoubleVector("VALUES",values);

            vardef = Teuchos::rcp(new FourierInterpolationVariable(varname,times,values,periodicdata));
          }
          else
          {
            dserror("unknown variable type");
          }

          // insert the variable in the vector of the variables of the function
          if (vardefinition == 1)
          {
            if (varid != 0)
            {
              functvarvector.push_back(functvar);
              functvar.clear();
            }
          }
          functvar.push_back(vardef);

          if (j == numrowsvar)
          {
            functvarvector.push_back(functvar);
          }
        }

        // add the expressions to the function vector
        for (int n = 0;n <= maxcomp; ++n)
        {
          vecfunc->AddExpr(functstring[n],functvarvector);
        }

        functions_.push_back(vecfunc);
      }
    }
  }
}


DRT::UTILS::Function& DRT::UTILS::FunctionManager::Funct(int num)
{
  // ensure that desired function is available (prevents segmentation fault)
  if (functions_.size()< (unsigned int)(num+1) || num<0)
    dserror("function %d not available",num+1);

  return *(functions_[num]);
}
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
DRT::UTILS::ExprFunction::ExprFunction()
{
  dim_ = DRT::Problem::Instance()->NDim();
  expr_.clear();
  exprdd_.clear();
  variables_.clear();
  isparsed_ = false;
}


DRT::UTILS::ExprFunction::~ExprFunction()
{
}


void DRT::UTILS::ExprFunction::AddExpr(std::string buf,std::vector<std::vector<Teuchos::RCP<FunctionVariable> > > variables)
{
  variables_ = variables;

  Teuchos::RCP< DRT::PARSER::Parser<double> > parser = Teuchos::rcp(new DRT::PARSER::Parser<double>(buf));
  parser->AddVariable("x",0);
  parser->AddVariable("y",0);
  parser->AddVariable("z",0);
  parser->AddVariable("t",0);

  Teuchos::RCP< DRT::PARSER::Parser<Sacado::Fad::DFad<Sacado::Fad::DFad<double> > > > parserdd = Teuchos::rcp(new DRT::PARSER::Parser<Sacado::Fad::DFad<Sacado::Fad::DFad<double> > >(buf));
  parserdd->AddVariable("x",0);
  parserdd->AddVariable("y",0);
  parserdd->AddVariable("z",0);
  parserdd->AddVariable("t",0);

  // add variables from the defined VARIABLES
  for(unsigned int i=0; i<variables.size(); ++i)
  {
    parser->AddVariable(variables_[i][0]->Name(),0);
    parserdd->AddVariable(variables_[i][0]->Name(),0);
  }

  // save parsers
  expr_.push_back(parser);
  exprdd_.push_back(parserdd);

  isparsed_ = false;

  return;
}


double DRT::UTILS::ExprFunction::Evaluate(const int index, const double* x, double t)
{
  // we consider a function of the type F = F ( x, y, z, t, v1(t), ..., vn(t) )
  if(not isparsed_)
    ParseExpressions();

  double index_mod = index;

  if (expr_.size() == 1)
  {
    index_mod = 0;
  }

  // set spatial variables
  if(dim_>0)
    expr_[index_mod]->SetValue("x",x[0]);
  if(dim_>1)
    expr_[index_mod]->SetValue("y",x[1]);
  if(dim_>2)
    expr_[index_mod]->SetValue("z",x[2]);

  // set temporal variable
  expr_[index_mod]->SetValue("t",t);

  // set the values of the variables at time t
  for (unsigned int i = 0; i<variables_.size(); ++i)
  {
    // find the right definition of the variable according to the hierarchy
    unsigned int n = 0;
    bool containtime = false;
    while (containtime == false)
    {
      if (n == variables_[i].size())
      {
        dserror("the VARIABLE %d is not defined at time %f", i,t);
      }
      else
      {
        containtime = variables_[i][n]->ContainTime(t);
      }
      if (containtime == false)
      {
        ++n;
      }
    }

    expr_[index_mod]->SetValue(variables_[i][0]->Name(),variables_[i][n]->Value(t));
  }

  // evaluation of F = F ( x, y, z, t, v1, ..., vn )
  return expr_[index_mod]->Evaluate();
}


std::vector<double> DRT::UTILS::ExprFunction::EvaluateSpatialDerivative(
    const int index,
    const double* x,
    const double t
    )
{

  // parse expression if not already parsed
  if(not isparsed_)
    ParseExpressions();

  double index_mod = index;

  if (expr_.size() == 1)
  {
    index_mod = 0;
  }

  // define Fad object for evaluation
  typedef Sacado::Fad::DFad<Sacado::Fad::DFad<double> >  FAD;

  // define FAD variables
  // we consider a function of the type F = F ( x, y, z, t, v1(t), ..., vn(t) )
  FAD xfad(4+variables_.size(), 0, x[0]);
  FAD yfad(4+variables_.size(), 1, x[1]);
  FAD zfad(4+variables_.size(), 2, x[2]);
  FAD tfad(4+variables_.size(), 3, t);

  xfad.val() = Sacado::Fad::DFad<double>(4+variables_.size(), 0, x[0]);
  yfad.val() = Sacado::Fad::DFad<double>(4+variables_.size(), 1, x[1]);
  zfad.val() = Sacado::Fad::DFad<double>(4+variables_.size(), 2, x[2]);
  tfad.val() = Sacado::Fad::DFad<double>(4+variables_.size(), 3, t);

  std::vector<FAD> fadvectvars(variables_.size());
  for(unsigned int i=0; i<variables_.size(); ++i)
  {
    // find the right definition of the variable according to the hierarchy
    unsigned int n = 0;
    bool containtime = false;
    while (containtime == false)
    {
      if (n == variables_[i].size())
      {
        dserror("the variable %d is not defined in the time considered", i);
      }
      else
      {
        containtime = variables_[i][n]->ContainTime(t);
      }
      if (containtime == false)
      {
        ++n;
      }
    }
    fadvectvars[i] = FAD(4+variables_.size(), 4+i, variables_[i][n]->Value(t));
    fadvectvars[i].val() = Sacado::Fad::DFad<double>(4+variables_.size(), 4+i, variables_[i][n]->Value(t));
  }
  FAD fdfad;

  // set spatial variables
  switch(dim_)
  {
  case 3:
  {
    exprdd_[index_mod]->SetValue("x",xfad);
    exprdd_[index_mod]->SetValue("y",yfad);
    exprdd_[index_mod]->SetValue("z",zfad);
    break;
  }
  case 2:
  {
    exprdd_[index_mod]->SetValue("x",xfad);
    exprdd_[index_mod]->SetValue("y",yfad);
    exprdd_[index_mod]->SetValue("z",0);
    break;
  }
  case 1:
  {
    exprdd_[index_mod]->SetValue("x",xfad);
    exprdd_[index_mod]->SetValue("y",0);
    exprdd_[index_mod]->SetValue("z",0);
    break;
  }
  default: dserror("Problem dimension has to be 1, 2, or 3."); break;
  }

  // set temporal variable
  exprdd_[index_mod]->SetValue("t",tfad);

  // set the values of the variables at time t
  for(unsigned int i=0; i<variables_.size(); ++i)
  {
    exprdd_[index_mod]->SetValue(variables_[i][0]->Name(),fadvectvars[i]);
  }

  // evaluation of derivatives
  fdfad = exprdd_[index_mod]->Evaluate();

  // result vector
  std::vector<double> res(3,0.0);
  for(unsigned int d=0; d<3; ++d)
  {
    res[d] = fdfad.dx(d).val();
  }

  // return derivatives
  return res;
}


std::vector<double> DRT::UTILS::ExprFunction::EvaluateTimeDerivative(
    const int index,
    const double* x,
    const double t,
    const unsigned deg
)
{
  // resulting vector holding
  std::vector<double> res(deg+1);

  // add the value at time t
  res[0] = Evaluate(index,x,t);

  // parse expression if not already parsed
  if(not isparsed_)
    ParseExpressions();

  double index_mod = index;
  if (expr_.size() == 1)
  {
    index_mod = 0;
  }

  // define Fad object for evaluation
  typedef Sacado::Fad::DFad<Sacado::Fad::DFad<double> >  FAD;

  // define FAD variables
  // we consider a function of the type F = F ( x, y, z, t, v1(t), ..., vn(t) )
  FAD xfad(4+variables_.size(), 0, x[0]);
  FAD yfad(4+variables_.size(), 1, x[1]);
  FAD zfad(4+variables_.size(), 2, x[2]);
  FAD tfad(4+variables_.size(), 3, t);

  xfad.val() = Sacado::Fad::DFad<double>(4+variables_.size(), 0, x[0]);
  yfad.val() = Sacado::Fad::DFad<double>(4+variables_.size(), 1, x[1]);
  zfad.val() = Sacado::Fad::DFad<double>(4+variables_.size(), 2, x[2]);
  tfad.val() = Sacado::Fad::DFad<double>(4+variables_.size(), 3, t);

  std::vector<FAD> fadvectvars(variables_.size());
  for(unsigned int i=0; i<variables_.size(); ++i)
  {
    // find the right definition of the variable according to the hierarchy
    unsigned int n = 0;
    bool containtime = false;
    while (containtime == false)
    {
      if (n == variables_[i].size())
      {
        dserror("the variable %d is not defined in the time considered", i);
      }
      else
      {
        containtime = variables_[i][n]->ContainTime(t);
      }
      if (containtime == false)
      {
        ++n;
      }
    }
    fadvectvars[i] = FAD(4+variables_.size(), 4+i, variables_[i][n]->Value(t));
    fadvectvars[i].val() = Sacado::Fad::DFad<double>(4+variables_.size(), 4+i, variables_[i][n]->Value(t));
  }
  FAD fdfad;

  // set spatial variables
  switch(dim_)
  {
  case 3:
  {
    exprdd_[index_mod]->SetValue("x",xfad);
    exprdd_[index_mod]->SetValue("y",yfad);
    exprdd_[index_mod]->SetValue("z",zfad);
    break;
  }
  case 2:
  {
    exprdd_[index_mod]->SetValue("x",xfad);
    exprdd_[index_mod]->SetValue("y",yfad);
    exprdd_[index_mod]->SetValue("z",0);
    break;
  }
  case 1:
  {
    exprdd_[index_mod]->SetValue("x",xfad);
    exprdd_[index_mod]->SetValue("y",0);
    exprdd_[index_mod]->SetValue("z",0);
    break;
  }
  default: dserror("Problem dimension has to be 1, 2, or 3."); break;
  }

  // set temporal variable
  exprdd_[index_mod]->SetValue("t",tfad);

  // set the values of the variables at time t
  for(unsigned int i=0; i<variables_.size(); ++i)
  {
    exprdd_[index_mod]->SetValue(variables_[i][0]->Name(),fadvectvars[i]);
  }

  // add the 1st time derivative at time t
  if (deg >= 1)
  {
    // evaluation of derivatives
    fdfad = exprdd_[index_mod]->Evaluate();

    // evaluation of dF/dt applying the chain rule:
    // dF/dt = dF*/dt + sum_i(dF/dvi*dvi/dt)
    double fdfad_dt = fdfad.dx(3).val();
    for(unsigned int i=0; i<variables_.size(); ++i)
    {
      // find the right definition of the variable according to the hierarchy
      unsigned int n = 0;
      bool containtime = false;
      while (containtime == false)
      {
        containtime = variables_[i][n]->ContainTime(t);
        if (n == variables_[i].size())
        {
          dserror("the variable %d is not defined in the time considered", i);
        }
        if (containtime == false)
        {
          ++n;
        }
      }
      fdfad_dt += fdfad.dx(4+i).val() * variables_[i][n]->TimeDerivativeValue(t);
    }

    res[1]=fdfad_dt;
  }

  // add the 2nd time derivative at time t
  if (deg >= 2)
  {
    // evaluation of d^2F/dt^2 applying the chain rule:
    //d^2F/dt^2 = d(dF*/dt)/dt + sum_i{
    //                                [d(dF*/dt)/dvi + d(dF/dvi)/dt + sum_j[d(dF/dvi)/dvj * dvj/dt]] * dvi/dt +
    //                                 + dF/dvi * d^2vi/dt^2
    //                                }
    double fdfad_dt2 = fdfad.dx(3).dx(3);
    std::vector<double> fdfad_dt2_term(variables_.size());

    for(unsigned int i=0; i<variables_.size(); ++i)
    {
      fdfad_dt2_term[i] = 0;
      // find the right definition of the variable according to the hierarchy
      unsigned int n = 0;
      bool containtime = false;
      while (containtime == false)
      {
        containtime = variables_[i][n]->ContainTime(t);
        if (n == variables_[i].size())
        {
          dserror("the variable %d is not defined in the time considered", i);
        }
        if (containtime == false)
        {
          ++n;
        }
      }
      fdfad_dt2_term[i] += fdfad.dx(3).dx(4+i);
      fdfad_dt2_term[i] += fdfad.dx(4+i).dx(3);
      for(unsigned int j=0; j<variables_.size(); ++j)
      {
        // find the right definition of the variable according to the hierarchy
        unsigned int m = 0;
        bool containtime = false;
        while (containtime == false)
        {
          containtime = variables_[j][m]->ContainTime(t);
          if (m == variables_[j].size())
          {
            dserror("the variable %d is not defined in the time considered", j);
          }
          if (containtime == false)
          {
            ++m;
          }
        }
        fdfad_dt2_term[i] += fdfad.dx(4+i).dx(4+j) * variables_[j][m]->TimeDerivativeValue(t);
      }
      fdfad_dt2_term[i] *= variables_[i][n]->TimeDerivativeValue(t);
      fdfad_dt2_term[i] += fdfad.dx(4+i).val() * variables_[i][n]->TimeDerivativeValue(t,2);
      fdfad_dt2 += fdfad_dt2_term[i];
    }

    res[2]=fdfad_dt2;
  }

  // error for higher derivatives
  if (deg >= 3)
  {
    dserror("Higher time derivatives than second not supported!");
  }

  return res;
}


bool DRT::UTILS::ExprFunction::IsVariable(const int index,const std::string& varname) const
{
  if(index>(int)expr_.size()-1 || index<0)
    dserror("Tried to add a variable to a function in a not available dimension.");

  return expr_[index]->IsVariable(varname);
}


void DRT::UTILS::ExprFunction::AddVariable(const int index,const std::string& varname,double varvalue)
{
  if(index>(int)expr_.size()-1 || index<0)
    dserror("Tried to add a variable to a function in a not available dimension.");
  if(isparsed_)
    dserror("Function has already been parsed! Variables can no longer be added!");

  expr_[index]->AddVariable(varname,varvalue);
  exprdd_[index]->AddVariable(varname,varvalue);
}


void DRT::UTILS::ExprFunction::ParseExpressions()
{
  // define iterators
  std::vector<Teuchos::RCP<DRT::PARSER::Parser<double> > >::iterator it;
  std::vector<Teuchos::RCP<DRT::PARSER::Parser<Sacado::Fad::DFad<Sacado::Fad::DFad<double> > > > >::iterator itdd;

  // loop over expressions and parse them
  for( it=expr_.begin() ; it!=expr_.end() ; it++)
    (*it)->ParseFunction();

  // loop over expressions for derivatives and parse them
  for( itdd=exprdd_.begin() ; itdd!=exprdd_.end() ; itdd++)
    (*itdd)->ParseFunction();

  isparsed_ = true;
}
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
DRT::UTILS::VariableExprFunction::VariableExprFunction():
    isparsed_(false)
{
}


DRT::UTILS::VariableExprFunction::~VariableExprFunction()
{
}


void DRT::UTILS::VariableExprFunction::AddExpr(std::string buf,
                                       double x,
                                       double y,
                                       double z,
                                       std::vector<std::pair<std::string,double> > constants
  )
{
  // do the almost same as the expression function (base class) but do not yet parse!

  // build the parser for the function evaluation
  Teuchos::RCP< DRT::PARSER::Parser<double> > parser = Teuchos::rcp(new DRT::PARSER::Parser<double>(buf));

  // build the parser for the function derivative evaluation
  Teuchos::RCP< DRT::PARSER::Parser<Sacado::Fad::DFad<double> > > parserd =
      Teuchos::rcp(new DRT::PARSER::Parser<Sacado::Fad::DFad<double> >(buf));

  // add constants
  for(std::vector<std::pair<std::string,double> >::iterator it=constants.begin();it!=constants.end();it++)
    parser->AddVariable(it->first,it->second);
  for(std::vector<std::pair<std::string,double> >::iterator it=constants.begin();it!=constants.end();it++)
    parserd->AddVariable(it->first,it->second);

  // save the parsers
  expr_.push_back(parser);
  exprd_.push_back(parserd);

  isparsed_=false;

  return;
}


bool DRT::UTILS::VariableExprFunction::IsVariable(
    int index,
    const std::string& varname) const
{
  if(index>(int)expr_.size()-1 || index<0)
    dserror("Tried to add a variable to a function in a not available dimension.");

  return expr_[index]->IsVariable(varname);
}


void DRT::UTILS::VariableExprFunction::AddVariable(
    int index,
    const std::string& varname,
    double varvalue)
{
  if(index>(int)expr_.size()-1 || index<0)
    dserror("Tried to add a variable to a function in a not available dimension.");
  if(isparsed_)
    dserror("Function has already been parsed! Variables can no longer be added!");

  expr_[index]->AddVariable(varname,varvalue);
  exprd_[index]->AddVariable(varname,varvalue);

}


void DRT::UTILS::VariableExprFunction::ParseExpressions()
{
  // define iterators
  std::vector<Teuchos::RCP<DRT::PARSER::Parser<double> > >::iterator it;
  std::vector<Teuchos::RCP<DRT::PARSER::Parser<Sacado::Fad::DFad<double> > > >::iterator itd;

  // loop over expressions and parse them
  for( it=expr_.begin() ; it!=expr_.end() ; it++)
    (*it)->ParseFunction();
  // loop over expressions for derivatives and parse them
  for( itd=exprd_.begin() ; itd!=exprd_.end() ; itd++)
    (*itd)->ParseFunction();

  isparsed_=true;
}


double DRT::UTILS::VariableExprFunction::Evaluate(
    const int index,
    const std::vector<std::pair<std::string,double> >& variables)
{
  if(not isparsed_)
    ParseExpressions();
  // set the values of the variables
  std::vector<std::pair<std::string,double> >::const_iterator it;
  for( it=variables.begin() ; it!=variables.end() ; it++)
    expr_[index]->SetValue(it->first,it->second);

  // evaluate the function and return the result
  return expr_[index]->Evaluate();
}


double DRT::UTILS::VariableExprFunction::Evaluate(
    const int index,
    const std::vector<std::pair<std::string,double> >& variables,
    const std::vector<std::pair<std::string,double> >& constants)
{
  if(not isparsed_)
    ParseExpressions();
  // set the values of the variables
  std::vector<std::pair<std::string,double> >::const_iterator it;
  for( it=variables.begin() ; it!=variables.end() ; it++)
    expr_[index]->SetValue(it->first,it->second);
  // set the values of the constants
  for( it=constants.begin() ; it!=constants.end() ; it++)
  {
    if( expr_[index]->IsVariable(it->first) )
      expr_[index]->SetValue(it->first,it->second);
  }

  // evaluate the function and return the result
  return expr_[index]->Evaluate();
}


std::vector<double> DRT::UTILS::VariableExprFunction::EvaluateDerivative(
    const int index,
    const std::vector<std::pair<std::string,double> >& variables)
{
  if(not isparsed_)
    ParseExpressions();

  // Fad object for evaluation
  // sacado data type replaces "double"
  typedef Sacado::Fad::DFad<double>  FAD;

  // number of variables
  int numvariables = variables.size();

  // counter for variable numbering
  int counter = 0;

  // set the values of the variables
  std::vector<std::pair<std::string,double> >::const_iterator it;
  for( it=variables.begin() ; it!=variables.end() ; it++)
  {
    // for 1st order derivatives
    FAD varfad(numvariables, counter, it->second);
    // set the value in expression
    exprd_[index]->SetValue(it->first,varfad);
    //update counter
    counter++;
  }

  // evaluate the expression
  FAD fdfad = exprd_[index]->Evaluate();

  // resulting vector
  std::vector<double> res(numvariables);

  // fill the result vector
  for(int i=0;i<numvariables;i++)
    res[i]=fdfad.dx(i);

  return res;
}


std::vector<double> DRT::UTILS::VariableExprFunction::EvaluateDerivative(
    int index,
    const std::vector<std::pair<std::string,double> >& variables,
    const std::vector<std::pair<std::string,double> >& constants)
{
  if(not isparsed_)
    ParseExpressions();

  // Fad object for evaluation
  // sacado data type replaces "double"
  typedef Sacado::Fad::DFad<double>  FAD;

  // number of variables
  int numvariables = variables.size();

  // counter for variable numbering
  int counter = 0;

  // set the values of the variables
  std::vector<std::pair<std::string,double> >::const_iterator it;
  for( it=variables.begin() ; it!=variables.end() ; it++)
  {
    // for 1st order derivatives
    FAD varfad(numvariables, counter, it->second);
    // set the value in expression
    exprd_[index]->SetValue(it->first,varfad);
    //update counter
    counter++;
  }

  // set the values of the constants
  for( it=constants.begin() ; it!=constants.end() ; it++)
  {
    if( exprd_[index]->IsVariable(it->first) )
      // set the value in expression
      exprd_[index]->SetValue(it->first,it->second);
  }

  // evaluate the expression
  FAD fdfad = exprd_[index]->Evaluate();

  // resulting vector
  std::vector<double> res(numvariables);

  // fill the result vector
  for(int i=0;i<numvariables;i++)
    res[i]=fdfad.dx(i);

  return res;
}


double DRT::UTILS::VariableExprFunction::Evaluate(
    const int index,
    const double* x,
    const double t
)
{
  std::vector<std::pair<std::string,double> > variables;
  variables.reserve(dim_);

  switch(dim_)
  {
  case 3:
  {
    variables.push_back(std::pair<std::string,double>("x",x[0])); //-x_[index]
    variables.push_back(std::pair<std::string,double>("y",x[1])); //-y_[index]
    variables.push_back(std::pair<std::string,double>("z",x[2])); //-z_[index]
  }
  case 2:
  {
    variables.push_back(std::pair<std::string,double>("x",x[0])); //-x_[index]
    variables.push_back(std::pair<std::string,double>("y",x[1])); //-y_[index]
  }
  case 1:
  {
    variables.push_back(std::pair<std::string,double>("x",x[0])); //-x_[index]
  }
  }

  variables.push_back(std::pair<std::string,double>("t",t));

  return Evaluate(index,variables);
}


std::vector<double> DRT::UTILS::VariableExprFunction::EvaluateSpatialDerivative(
    int                  index,
    const double*        x,
    const double         t
    )
{
  std::vector<std::pair<std::string,double> > variables(4);

  variables[0] = std::pair<std::string,double>("x",x[0]);
  variables[1] = std::pair<std::string,double>("y",x[1]);
  variables[2] = std::pair<std::string,double>("z",x[2]);
  variables[3] = std::pair<std::string,double>("t",t);

  return EvaluateDerivative(index,variables);
}


/*----------------------------------------------------------------------*
 | Constructor of ControlledRotation                         hahn 04/13 |
 *----------------------------------------------------------------------*/
DRT::UTILS::ControlledRotationFunction::ControlledRotationFunction(
    std::string fileName, std::string type, double origin_x,
    double origin_y, double origin_z) :
Function(), NUMMANEUVERCELLS_(4)
{
    // Initialize variables
    // *****************************************************************************

    // Initialize condition type
    if (type == "STRUCTURE") {    // structure
        type_ = 1;
    } else if (type == "FLUID") { // fluid
        type_ = 2;
    } else {
        dserror("When using the function CONTROLLEDROTATION, the type must be either 'STRUCTURE' or 'FLUID'");
    }

    // Initialize origin, about which the rotation is performed
    origin_(0,0) = origin_x;
    origin_(1,0) = origin_y;
    origin_(2,0) = origin_z;

    // Initialize time of previous time step (at t-deltaT)
    timeOld_ = 0.0;

    // Initialize previous angular acceleration (at t-deltaT)
    omegaDotOld_B_(0,0) = 0.0;
    omegaDotOld_B_(1,0) = 0.0;
    omegaDotOld_B_(2,0) = 0.0;

    // Current angular rate (at t)
    omega_B_(0,0) = 0.0;
    omega_B_(1,0) = 0.0;
    omega_B_(2,0) = 0.0;

    // Initialize satellite's current attitude trafo matrix from B- to I-system (at t)
    satAtt_dcm_IB_(0,0) = 1.0; satAtt_dcm_IB_(0,1) = 0.0; satAtt_dcm_IB_(0,2) = 0.0;
    satAtt_dcm_IB_(1,0) = 0.0; satAtt_dcm_IB_(1,1) = 1.0; satAtt_dcm_IB_(1,2) = 0.0;
    satAtt_dcm_IB_(2,0) = 0.0; satAtt_dcm_IB_(2,1) = 0.0; satAtt_dcm_IB_(2,2) = 1.0;

    // Initialize satellite's current attitude quaternion from B- to I-system (at t)
    satAtt_q_IB_(0,0) = 0.0;
    satAtt_q_IB_(1,0) = 0.0;
    satAtt_q_IB_(2,0) = 0.0;
    satAtt_q_IB_(3,0) = 1.0;

    // Initialize maneuvers
    maneuvers_.clear();

    // Read maneuver file and fill maneuvers variable
    // *****************************************************************************

    std::string line;
    std::stringstream lineStream;
    std::string cell;

    // Open file
    std::ifstream file (fileName.c_str());
    if (!file.is_open()) {
        dserror("Unable to open file: %s", fileName.c_str());
    }

    // Loop through all lines
    while (getline (file,line)) {
        if (!line.empty()) {
            // Clear local variables
            lineStream.clear();
            cell.clear();

            // Obtain all numManeuverCells=4 values from current line (t, omegaDot_x_B, omegaDot_y_B, omegaDot_z_B)
            lineStream << line;
            for (int i=0; i<NUMMANEUVERCELLS_; i++) {
                // Obtain i-th cell from current line
                getline(lineStream, cell, ' ');

                // If empty cell, than either empty line or one cell in the line
                // missing, anyhow an error.
                if (cell.empty()) {
                    dserror("Error during reading of file: %s", fileName.c_str());
                }

                // Convert input cell from string to double
                double cellDouble = (double)strtod(cell.c_str(), NULL);

                // Add cell to maneuvers vector
                maneuvers_.push_back(cellDouble);
            }
        }
    }

    // Close file
    file.close();

    // Determine number of maneuvers
    numManeuvers_ = (int)(maneuvers_.size()) / (int)NUMMANEUVERCELLS_;

    // Output maneuver list
    printf("\n=================================================================================\n");
    printf("ControlledRotation - %s\n", type.c_str());
    printf("---------------------------------------------------------------------------------\n");
    printf("The following %d maneuvers have been loaded from file %s:\n", numManeuvers_, fileName.c_str());
    for (int i=0; i<numManeuvers_; i++) {
        printf("Time: %e  OmegaDot_B: %e, %e, %e \n",
                maneuvers_[i*NUMMANEUVERCELLS_+0], maneuvers_[i*NUMMANEUVERCELLS_+1],
                maneuvers_[i*NUMMANEUVERCELLS_+2], maneuvers_[i*NUMMANEUVERCELLS_+3]);
    }
    printf("=================================================================================\n\n");

}

/*----------------------------------------------------------------------*
 | Evaluate ControlledRotation and return for structures the current    |
 | displacement and for fluids the current velocity of the respective   |
 | node for the given index.                                 hahn 04/13 |
 *----------------------------------------------------------------------*/
double DRT::UTILS::ControlledRotationFunction::Evaluate(const int index,
    const double* xp, double t)
{
    // Check, if a restart has been performed
    // *****************************************************************************
    const int step = DRT::Problem::Instance()->Restart();
    if ((step > 0) && (timeOld_ == 0.0)) {
      dserror("When using the function CONTROLLEDROTATION, the restart functionality cannot be used!");
    }

    // If new time step, apply angular acceleration (if desired) and determine
    // new attitude satAtt_dcm_IB_
    // *****************************************************************************
    // Determine time difference
    double deltaT = t - timeOld_;

    if (deltaT > 1e-12) { // new time step

        // Determine current angular acceleration (at t)
        // -----------------------------------------------------------------------------
        LINALG::Matrix<3,1> omegaDot_B;
        omegaDot_B(0,0) = 0.0;
        omegaDot_B(1,0) = 0.0;
        omegaDot_B(2,0) = 0.0;

        for (int i=0; i<numManeuvers_; i++) {
            if (t >= maneuvers_[i*NUMMANEUVERCELLS_+0]) {
                omegaDot_B(0,0) = maneuvers_[i*NUMMANEUVERCELLS_+1];
                omegaDot_B(1,0) = maneuvers_[i*NUMMANEUVERCELLS_+2];
                omegaDot_B(2,0) = maneuvers_[i*NUMMANEUVERCELLS_+3];
            }
        }

        // Calculate current angular rate (at t) by integration
        // of angular acceleration (trapezoidal rule):
        // omega_(t) = deltaT * (omegaDotOld + omegaDot) / 2 + omega_(t-deltaT)
        // -----------------------------------------------------------------------------
        LINALG::Matrix<3,1> deltaOmega;
        deltaOmega.Update(omegaDotOld_B_, omegaDot_B); // 1) deltaOmega <- omegaDotOld_ + omegaDot
        deltaOmega.Scale(deltaT/2.0);                  // 2) deltaOmega <- deltaOmega * deltaT / 2.0
        omega_B_ += deltaOmega;                        // 3) omega_ <- omega_ + deltaOmega

        /* // Debugging output
           cout << "omegaDot: "; omegaDot_B.Print(cout); // Print omegaDot_B
           cout << "omega:    "; omega_B_.Print(cout);   // Print omega_B_
        */

        omegaDotOld_B_ = omegaDot_B; // Set omegaDotOld_B_ for next time step

        // Calculate new attitude quaternion satAtt_q_IB_ [Wertz, p. 511f]
        // -----------------------------------------------------------------------------
        LINALG::Matrix<4,4> mOmega; // Skew-symmetric matrix containing angular velocity components
        mOmega(0,0) =  0.0;           mOmega(0,1) =  omega_B_(2,0); mOmega(0,2) = -omega_B_(1,0); mOmega(0,3) = omega_B_(0,0);
        mOmega(1,0) = -omega_B_(2,0); mOmega(1,1) =  0.0;           mOmega(1,2) =  omega_B_(0,0); mOmega(1,3) = omega_B_(1,0);
        mOmega(2,0) =  omega_B_(1,0); mOmega(2,1) = -omega_B_(0,0); mOmega(2,2) =  0.0;           mOmega(2,3) = omega_B_(2,0);
        mOmega(3,0) = -omega_B_(0,0); mOmega(3,1) = -omega_B_(1,0); mOmega(3,2) = -omega_B_(2,0); mOmega(3,3) = 0.0;

        mOmega.Scale(deltaT/2.0);
        mOmega(0,0) = 1.0;
        mOmega(1,1) = 1.0;
        mOmega(2,2) = 1.0;
        mOmega(3,3) = 1.0;

        LINALG::Matrix<4,1> satAtt_q_IB_TMP(true);
        satAtt_q_IB_TMP.Multiply(mOmega, satAtt_q_IB_);
        satAtt_q_IB_ = satAtt_q_IB_TMP;

        satAtt_q_IB_.Scale(1/satAtt_q_IB_.Norm2()); // Normalize attitude quaternion

        // Create transformation matrix satAtt_dcm_IB_ [Wertz, (E-8)]
        // -----------------------------------------------------------------------------
        const double q1 = satAtt_q_IB_(0,0);
        const double q2 = satAtt_q_IB_(1,0);
        const double q3 = satAtt_q_IB_(2,0);
        const double q4 = satAtt_q_IB_(3,0);

        satAtt_dcm_IB_(0,0) = q1*q1-q2*q2-q3*q3+q4*q4; satAtt_dcm_IB_(0,1) = 2.0*(q1*q2-q3*q4);        satAtt_dcm_IB_(0,2) = 2.0*(q1*q3+q2*q4);
        satAtt_dcm_IB_(1,0) = 2.0*(q1*q2+q3*q4);       satAtt_dcm_IB_(1,1) = -q1*q1+q2*q2-q3*q3+q4*q4; satAtt_dcm_IB_(1,2) = 2.0*(q2*q3-q1*q4);
        satAtt_dcm_IB_(2,0) = 2.0*(q1*q3-q2*q4);       satAtt_dcm_IB_(2,1) = 2.0*(q2*q3+q1*q4);        satAtt_dcm_IB_(2,2) = -q1*q1-q2*q2+q3*q3+q4*q4;

        // Update time of last time step
        // -----------------------------------------------------------------------------
        timeOld_ = t;
    }

    // Obtain the current node position in the inertial system
    // *****************************************************************************
    // NOTE: 1) Here it is assumed that the difference between the mesh system and the
    //          body system is just given by a constant displacement named origin_.
    //          Hence the variable origin_ specifies the point, given in the mesh
    //          system, about which is being rotated.
    //       2) The inertial system used here has position and attitude of the body
    //          system at initial time. Thus it is deviating from the ECI system J2000
    //          at most by a constant displacement that is due to the satellite position
    //          at initial time and a constant rotation that is due to the satellite's
    //          attitude at initial time. Due to Galilei-invariance this shouldn't be
    //          a problem and it can easily be converted to the J2000 system.

    // Node reference position, given in the body system
    LINALG::Matrix<3,1> nodeReferencePos_B;
    nodeReferencePos_B(0,0) = xp[0] - origin_(0,0);
    nodeReferencePos_B(1,0) = xp[1] - origin_(1,0);
    nodeReferencePos_B(2,0) = xp[2] - origin_(2,0);

    // Node position, given in the inertial system
    LINALG::Matrix<3,1> nodePos_I;
    nodePos_I.Multiply(satAtt_dcm_IB_, nodeReferencePos_B);

    // Calculate and return the displacement/velocity of the node for the given index
    // *****************************************************************************

    if (t<=0) {
        // Return zero displacement/translational velocity of the node for the given index
        return 0.0;
    }

    // Displacement of the node for the given index
    const double dispNodePos_I = nodePos_I(index,0) - nodeReferencePos_B(index,0);

    if (type_ == 1) { // structure

        // Return the displacement of the node for the given index
        return dispNodePos_I;

    } else if (type_ == 2) { // fluid

        // Node velocity, given in the inertial system: v = omega x r
        double nodeVel_I[3];
        nodeVel_I[0] = omega_B_(1,0) * nodePos_I(2,0) - omega_B_(2,0) * nodePos_I(1,0);
        nodeVel_I[1] = omega_B_(2,0) * nodePos_I(0,0) - omega_B_(0,0) * nodePos_I(2,0);
        nodeVel_I[2] = omega_B_(0,0) * nodePos_I(1,0) - omega_B_(1,0) * nodePos_I(0,0);

        // Return the translational velocity of the node for the given index
        return (nodeVel_I[index]);

    } else {
        dserror("When using the function CONTROLLEDROTATION, the type must be either 'STRUCTURE' or 'FLUID'");
        return 0.0;
    }
}

/*----------------------------------------------------------------------*
 | Constructor of AccelerationProfile                        hahn 09/13 |
 *----------------------------------------------------------------------*/
DRT::UTILS::AccelerationProfileFunction::AccelerationProfileFunction(std::string fileName) :
Function(), NUMACCELERATIONCELLS_(4)
{
    // Initialize variables
    // *****************************************************************************

    // Initialize time of previous time step (at t-deltaT)
    timeOld_ = 0.0;

    // Initialize accelerations
    accelerations_.clear();

    // Initialize current acceleration (at t)
    acc_B_(0,0) = 0.0;
    acc_B_(1,0) = 0.0;
    acc_B_(2,0) = 0.0;

    // Read acceleration profile file and fill acceleration variable
    // *****************************************************************************

    std::string line;
    std::stringstream lineStream;
    std::string cell;

    // Open file
    std::ifstream file (fileName.c_str());
    if (!file.is_open()) {
        dserror("Unable to open file: %s", fileName.c_str());
    }

    // Loop through all lines
    while (getline (file,line)) {
        if (!line.empty()) {
            // Clear local variables
            lineStream.clear();
            cell.clear();

            // Obtain all numAccelerationCells=4 values from current line (t, acc_x_B, acc_y_B, acc_z_B)
            lineStream << line;
            for (int i=0; i<NUMACCELERATIONCELLS_; i++) {
                // Obtain i-th cell from current line
                getline(lineStream, cell, ' ');

                // If empty cell, than either empty line or one cell in the line
                // missing, anyhow an error.
                if (cell.empty()) {
                    dserror("Error during reading of file: %s", fileName.c_str());
                }

                // Convert input cell from string to double
                double cellDouble = (double)strtod(cell.c_str(), NULL);

                // Add cell to acceleration vector
                accelerations_.push_back(cellDouble);
            }
        }
    }

    // Close file
    file.close();

    // Determine number of accelerations
    numAccelerations_ = (int)(accelerations_.size()) / (int)NUMACCELERATIONCELLS_;

    // Output acceleration list
    printf("\n=================================================================================\n");
    printf("AccelerationProfile\n");
    printf("---------------------------------------\n");
    printf("The following %d acceleration rows have been loaded from file %s:\n", numAccelerations_, fileName.c_str());
    for (int i=0; i<numAccelerations_; i++) {
        printf("Time: %e  acc_B: %e, %e, %e \n",
                accelerations_[i*NUMACCELERATIONCELLS_+0], accelerations_[i*NUMACCELERATIONCELLS_+1],
                accelerations_[i*NUMACCELERATIONCELLS_+2], accelerations_[i*NUMACCELERATIONCELLS_+3]);
    }
    printf("=================================================================================\n\n");

}

/*---------------------------------------------------------------------*
 | Evaluate AccelerationProfile and return the respective acceleration |
 | of the node for the given index.                         hahn 09/13 |
 *---------------------------------------------------------------------*/
double DRT::UTILS::AccelerationProfileFunction::Evaluate(const int index,
    const double* xp, double t)
{
    // Determine time difference
    double deltaT = t - timeOld_;

    // If new time step, determine current acceleration
    // *****************************************************************************
    if (deltaT > 1e-9) { // new time step

        // Determine current acceleration (at t)
        // -------------------------------------------------------------------------
        acc_B_(0,0) = 0.0;
        acc_B_(1,0) = 0.0;
        acc_B_(2,0) = 0.0;

        for (int i=0; i<numAccelerations_; i++) {
            if (t >= accelerations_[i*NUMACCELERATIONCELLS_+0]) {
                acc_B_(0,0) = accelerations_[i*NUMACCELERATIONCELLS_+1];
                acc_B_(1,0) = accelerations_[i*NUMACCELERATIONCELLS_+2];
                acc_B_(2,0) = accelerations_[i*NUMACCELERATIONCELLS_+3];
            }
        }

        // Update time of last time step
        // -------------------------------------------------------------------------
        timeOld_ = t;
    }

    // Return the acceleration of the node for the given index
    // *****************************************************************************
    return acc_B_(index,0);
}

