/*----------------------------------------------------------------------*/
/*!
\file inpar_s2i.cpp

\brief input quantities and globally accessible enumerations for scatra-scatra interface coupling

\level 2

<pre>
\maintainer Rui Fang
            fang@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/
            089-289-15251
</pre>
*/
/*----------------------------------------------------------------------*/
#include "inpar_s2i.H"

#include "drt_validparameters.H"

#include "../drt_lib/drt_conditiondefinition.H"

/*------------------------------------------------------------------------*
 | set valid parameters for scatra-scatra interface coupling   fang 01/16 |
 *------------------------------------------------------------------------*/
void INPAR::S2I::SetValidParameters(Teuchos::RCP<Teuchos::ParameterList> list)
{
  using namespace DRT::INPUT;
  using Teuchos::tuple;
  using Teuchos::setStringToIntegralParameter;

  Teuchos::ParameterList& s2icoupling = list->sublist("SCALAR TRANSPORT DYNAMIC",true).sublist(
      "S2I COUPLING",
      false,
      "control parameters for scatra-scatra interface coupling"
      );

  // type of global system matrix in global system of equations
  setStringToIntegralParameter<int>(
      "MATRIXTYPE",
      "sparse",
      "type of global system matrix in global system of equations",
      tuple<std::string>(
          "sparse",
          "block_geometry",
          "block_condition",
          "block_condition_dof"
          ),
      tuple<int>(
          matrix_sparse,
          matrix_block_geometry,
          matrix_block_condition,
          matrix_block_condition_dof
          ),
      &s2icoupling
      );

  // type of mortar meshtying
  setStringToIntegralParameter<int>(
      "COUPLINGTYPE",
      "Undefined",
      "type of mortar meshtying",
      tuple<std::string>(
          "Undefined",
          "MatchingNodes",
          "StandardMortar",
          "SaddlePointMortar_Petrov",
          "SaddlePointMortar_Bubnov",
          "CondensedMortar_Petrov",
          "CondensedMortar_Bubnov",
          "StandardNodeToSegment"
          ),
      tuple<int>(
          coupling_undefined,
          coupling_matching_nodes,
          coupling_mortar_standard,
          coupling_mortar_saddlepoint_petrov,
          coupling_mortar_saddlepoint_bubnov,
          coupling_mortar_condensed_petrov,
          coupling_mortar_condensed_bubnov,
          coupling_nts_standard
          ),
      &s2icoupling
      );

  // flag for equilibration of global system of equations
  setStringToIntegralParameter<int>(
      "EQUILIBRATION",
      "none",
      "flag for equilibration of global system of equations",
      tuple<std::string>(
          "none",
          "rows",
          "columns",
          "full"
          ),
      tuple<int>(
          equilibration_none,
          equilibration_rows,
          equilibration_columns,
          equilibration_full
          ),
      &s2icoupling
      );

  // flag for interface side underlying Lagrange multiplier definition
  setStringToIntegralParameter<int>(
      "LMSIDE",
      "slave",
      "flag for interface side underlying Lagrange multiplier definition",
      tuple<std::string>(
          "slave",
          "master"
          ),
      tuple<int>(
          side_slave,
          side_master
          ),
      &s2icoupling
      );

  // flag for evaluation of interface linearizations and residuals on slave side only
  BoolParameter("SLAVEONLY","No","flag for evaluation of interface linearizations and residuals on slave side only",&s2icoupling);

  // node-to-segment projection tolerance
  DoubleParameter("NTSPROJTOL",0.0,"node-to-segment projection tolerance",&s2icoupling);
}


/*------------------------------------------------------------------------*
 | set valid conditions for scatra-scatra interface coupling   fang 01/16 |
 *------------------------------------------------------------------------*/
void INPAR::S2I::SetValidConditions(std::vector<Teuchos::RCP<DRT::INPUT::ConditionDefinition> >& condlist)
{
  using namespace DRT::INPUT;

  /*--------------------------------------------------------------------*/
  // scatra-scatra interface coupling
  {
    // definition of scatra-scatra interface coupling line condition
    Teuchos::RCP<ConditionDefinition> s2iline =
        Teuchos::rcp(new ConditionDefinition("DESIGN S2I COUPLING LINE CONDITIONS",
                                             "S2ICoupling",
                                             "Scatra-scatra line interface coupling",
                                             DRT::Condition::S2ICoupling,
                                             true,
                                             DRT::Condition::Line));

    // definition of scatra-scatra interface coupling surface condition
    Teuchos::RCP<ConditionDefinition> s2isurf =
        Teuchos::rcp(new ConditionDefinition("DESIGN S2I COUPLING SURF CONDITIONS",
                                             "S2ICoupling",
                                             "Scatra-scatra surface interface coupling",
                                             DRT::Condition::S2ICoupling,
                                             true,
                                             DRT::Condition::Surface));

    // equip condition definitions with input file line components
    std::vector<Teuchos::RCP<ConditionComponent> > s2icomponents;

    {
      // interface ID
      s2icomponents.push_back(Teuchos::rcp(new IntConditionComponent("ConditionID")));

      // interface sides for scatra-scatra interface coupling
      std::vector<Teuchos::RCP<CondCompBundle> > interfacesides;
      {
        {
          // undefined side
          std::vector<Teuchos::RCP<ConditionComponent> > undefinedside;

          // insert undefined-side condition components into vector of interface sides
          interfacesides.push_back(Teuchos::rcp(new CondCompBundle("Undefined",undefinedside,INPAR::S2I::side_undefined)));
        }

        {
          // slave side
          std::vector<Teuchos::RCP<ConditionComponent> > slaveside;

          // kinetic models for scatra-scatra interface coupling
          std::vector<Teuchos::RCP<CondCompBundle> > kineticmodels;
          {
            {
              // constant permeability
              std::vector<Teuchos::RCP<ConditionComponent> > constperm;
              constperm.push_back(Teuchos::rcp(new SeparatorConditionComponent("numscal")));                // total number of existing scalars
              std::vector<Teuchos::RCP<SeparatorConditionComponent> > intsepcomp;                           // empty vector --> no separators for integer vectors needed
              std::vector<Teuchos::RCP<IntVectorConditionComponent> > intvectcomp;                          // empty vector --> no integer vectors needed
              std::vector<Teuchos::RCP<SeparatorConditionComponent> > realsepcomp;
              realsepcomp.push_back(Teuchos::rcp(new SeparatorConditionComponent("permeabilities")));       // string separator in front of real permeability vector in input file line
              std::vector<Teuchos::RCP<RealVectorConditionComponent> > realvectcomp;
              realvectcomp.push_back(Teuchos::rcp(new RealVectorConditionComponent("permeabilities",0)));   // real vector of constant permeabilities
              constperm.push_back(Teuchos::rcp(new IntRealBundle(
                  "permeabilities",
                  Teuchos::rcp(new IntConditionComponent("numscal")),
                  intsepcomp,
                  intvectcomp,
                  realsepcomp,
                  realvectcomp
              )));

              kineticmodels.push_back(Teuchos::rcp(new CondCompBundle("ConstantPermeability",constperm,INPAR::S2I::kinetics_constperm)));
            }

            {
              // Butler-Volmer
              std::vector<Teuchos::RCP<ConditionComponent> > butlervolmer;
              butlervolmer.push_back(Teuchos::rcp(new SeparatorConditionComponent("numscal")));            // total number of existing scalars
              std::vector<Teuchos::RCP<SeparatorConditionComponent> > intsepcomp;
              intsepcomp.push_back(Teuchos::rcp(new SeparatorConditionComponent("stoichiometries")));
              std::vector<Teuchos::RCP<IntVectorConditionComponent> > intvectcomp;                         // string separator in front of integer stoichiometry vector in input file line
              intvectcomp.push_back(Teuchos::rcp(new IntVectorConditionComponent("stoichiometries",0)));   // integer vector of stoichiometric coefficients
              std::vector<Teuchos::RCP<SeparatorConditionComponent> > realsepcomp;                         // empty vector --> no separators for real vectors needed
              std::vector<Teuchos::RCP<RealVectorConditionComponent> > realvectcomp;                       // empty vector --> no real vectors needed
              butlervolmer.push_back(Teuchos::rcp(new IntRealBundle(
                  "stoichiometries",
                  Teuchos::rcp(new IntConditionComponent("numscal")),
                  intsepcomp,
                  intvectcomp,
                  realsepcomp,
                  realvectcomp
              )));
              butlervolmer.push_back(Teuchos::rcp(new SeparatorConditionComponent("e-")));
              butlervolmer.push_back(Teuchos::rcp(new IntConditionComponent("e-")));
              butlervolmer.push_back(Teuchos::rcp(new SeparatorConditionComponent("k_r")));
              butlervolmer.push_back(Teuchos::rcp(new RealConditionComponent("k_r")));
              butlervolmer.push_back(Teuchos::rcp(new SeparatorConditionComponent("alpha_a")));
              butlervolmer.push_back(Teuchos::rcp(new RealConditionComponent("alpha_a")));
              butlervolmer.push_back(Teuchos::rcp(new SeparatorConditionComponent("alpha_c")));
              butlervolmer.push_back(Teuchos::rcp(new RealConditionComponent("alpha_c")));

              kineticmodels.push_back(Teuchos::rcp(new CondCompBundle("Butler-Volmer",butlervolmer,INPAR::S2I::kinetics_butlervolmer)));
            }

            {
              // Butler-Volmer-Peltier
              std::vector<Teuchos::RCP<ConditionComponent> > butlervolmerpeltier;
              butlervolmerpeltier.push_back(Teuchos::rcp(new SeparatorConditionComponent("numscal")));            // total number of existing scalars
              std::vector<Teuchos::RCP<SeparatorConditionComponent> > intsepcomp;
              intsepcomp.push_back(Teuchos::rcp(new SeparatorConditionComponent("stoichiometries")));
              std::vector<Teuchos::RCP<IntVectorConditionComponent> > intvectcomp;                         // string separator in front of integer stoichiometry vector in input file line
              intvectcomp.push_back(Teuchos::rcp(new IntVectorConditionComponent("stoichiometries",0)));   // integer vector of stoichiometric coefficients
              std::vector<Teuchos::RCP<SeparatorConditionComponent> > realsepcomp;                         // empty vector --> no separators for real vectors needed
              std::vector<Teuchos::RCP<RealVectorConditionComponent> > realvectcomp;                       // empty vector --> no real vectors needed
              butlervolmerpeltier.push_back(Teuchos::rcp(new IntRealBundle(
                  "stoichiometries",
                  Teuchos::rcp(new IntConditionComponent("numscal")),
                  intsepcomp,
                  intvectcomp,
                  realsepcomp,
                  realvectcomp
              )));
              butlervolmerpeltier.push_back(Teuchos::rcp(new SeparatorConditionComponent("e-")));
              butlervolmerpeltier.push_back(Teuchos::rcp(new IntConditionComponent("e-")));
              butlervolmerpeltier.push_back(Teuchos::rcp(new SeparatorConditionComponent("k_r")));
              butlervolmerpeltier.push_back(Teuchos::rcp(new RealConditionComponent("k_r")));
              butlervolmerpeltier.push_back(Teuchos::rcp(new SeparatorConditionComponent("alpha_a")));
              butlervolmerpeltier.push_back(Teuchos::rcp(new RealConditionComponent("alpha_a")));
              butlervolmerpeltier.push_back(Teuchos::rcp(new SeparatorConditionComponent("alpha_c")));
              butlervolmerpeltier.push_back(Teuchos::rcp(new RealConditionComponent("alpha_c")));
              butlervolmerpeltier.push_back(Teuchos::rcp(new SeparatorConditionComponent("peltier")));
              butlervolmerpeltier.push_back(Teuchos::rcp(new RealConditionComponent("peltier")));

              kineticmodels.push_back(Teuchos::rcp(new CondCompBundle("Butler-Volmer-Peltier",butlervolmerpeltier,INPAR::S2I::kinetics_butlervolmerpeltier)));
            }
          } // kinetic models for scatra-scatra interface coupling

          // insert kinetic models into vector with slave-side condition components
          slaveside.push_back(Teuchos::rcp(new SeparatorConditionComponent("KineticModel")));
          slaveside.push_back(Teuchos::rcp(new CondCompBundleSelector(
              "kinetic models for scatra-scatra interface coupling",
              Teuchos::rcp(new StringConditionComponent(
                  "kinetic model",
                  "ConstantPermeability",
                  Teuchos::tuple<std::string>("ConstantPermeability","Butler-Volmer","Butler-Volmer-Peltier"),
                  Teuchos::tuple<int>(INPAR::S2I::kinetics_constperm,INPAR::S2I::kinetics_butlervolmer,INPAR::S2I::kinetics_butlervolmerpeltier))),
              kineticmodels)));

          // insert slave-side condition components into vector of interface sides
          interfacesides.push_back(Teuchos::rcp(new CondCompBundle("Slave",slaveside,INPAR::S2I::side_slave)));
        }

        {
          // master side
          std::vector<Teuchos::RCP<ConditionComponent> > masterside;

          // insert master-side condition components into vector of interface sides
          interfacesides.push_back(Teuchos::rcp(new CondCompBundle("Master",masterside,INPAR::S2I::side_master)));
        }
      } // interface sides for scatra-scatra interface coupling

      // insert interface sides into vector with input file line components
      s2icomponents.push_back(Teuchos::rcp(new CondCompBundleSelector(
          "interface sides for scatra-scatra interface coupling",
          Teuchos::rcp(new StringConditionComponent(
              "interface side",
              "Undefined",
              Teuchos::tuple<std::string>("Undefined","Slave","Master"),
              Teuchos::tuple<int>(INPAR::S2I::side_undefined,INPAR::S2I::side_slave,INPAR::S2I::side_master))),
          interfacesides)));
    }

    // insert input file line components into condition definitions
    for (unsigned i=0; i<s2icomponents.size(); ++i)
    {
      s2iline->AddComponent(s2icomponents[i]);
      s2isurf->AddComponent(s2icomponents[i]);
    }

    // insert condition definitions into global list of valid condition definitions
    condlist.push_back(s2iline);
    condlist.push_back(s2isurf);
  }


  /*--------------------------------------------------------------------*/
  // scatra-scatra interface coupling (domain partitioning for block preconditioning of global system matrix)
  {
    // partitioning of 2D domain into 2D subdomains
    Teuchos::RCP<ConditionDefinition> s2ilinepartitioning =
        Teuchos::rcp(new ConditionDefinition("DESIGN S2I COUPLING SURF CONDITIONS / PARTITIONING",
                                             "S2ICouplingPartitioning",
                                             "Scatra-scatra line interface coupling (domain partitioning)",
                                             DRT::Condition::S2ICouplingPartitioning,
                                             false,
                                             DRT::Condition::Surface));

    // partitioning of 3D domain into 3D subdomains
    Teuchos::RCP<ConditionDefinition> s2isurfpartitioning =
        Teuchos::rcp(new ConditionDefinition("DESIGN S2I COUPLING VOL CONDITIONS / PARTITIONING",
                                             "S2ICouplingPartitioning",
                                             "Scatra-scatra surface interface coupling (domain partitioning)",
                                             DRT::Condition::S2ICouplingPartitioning,
                                             false,
                                             DRT::Condition::Volume));

    // insert condition definitions into global list of valid condition definitions
    condlist.push_back(s2ilinepartitioning);
    condlist.push_back(s2isurfpartitioning);
  }

  return;
}
