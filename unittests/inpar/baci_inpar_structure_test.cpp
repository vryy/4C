/*----------------------------------------------------------------------*/
/*! \file
\brief Unit tests for input parameters of structure field.

\level 0
*/
/*----------------------------------------------------------------------*/
#include <gtest/gtest.h>

#include "baci_inpar_structure.H"

namespace
{
  TEST(InparStructureTest, String2ModelType)
  {
    using namespace BACI;

    EXPECT_EQ(INPAR::STR::String2ModelType("Structure"), INPAR::STR::model_structure);
    EXPECT_EQ(INPAR::STR::String2ModelType("Contact"), INPAR::STR::model_contact);
    EXPECT_EQ(INPAR::STR::String2ModelType("Meshtying"), INPAR::STR::model_meshtying);
    EXPECT_EQ(INPAR::STR::String2ModelType("SpringDashpot"), INPAR::STR::model_springdashpot);
    EXPECT_EQ(
        INPAR::STR::String2ModelType("Lag-Pen-Constraint"), INPAR::STR::model_lag_pen_constraint);
    EXPECT_EQ(INPAR::STR::String2ModelType("Basic-Coupling"), INPAR::STR::model_basic_coupling);
    EXPECT_EQ(
        INPAR::STR::String2ModelType("Monolithic-Coupling"), INPAR::STR::model_monolithic_coupling);
    EXPECT_EQ(INPAR::STR::String2ModelType("Partitioned-Coupling"),
        INPAR::STR::model_partitioned_coupling);
    EXPECT_EQ(
        INPAR::STR::String2ModelType("BeamInteractionOld"), INPAR::STR::model_beam_interaction_old);
    EXPECT_EQ(INPAR::STR::String2ModelType("BrownianDyn"), INPAR::STR::model_browniandyn);
    EXPECT_EQ(INPAR::STR::String2ModelType("BeamInteraction"), INPAR::STR::model_beaminteraction);
  }

  TEST(InparStructureTest, ModelTypeString)
  {
    using namespace BACI;

    EXPECT_EQ(INPAR::STR::ModelTypeString(INPAR::STR::model_structure), "Structure");
    EXPECT_EQ(INPAR::STR::ModelTypeString(INPAR::STR::model_contact), "Contact");
    EXPECT_EQ(INPAR::STR::ModelTypeString(INPAR::STR::model_meshtying), "Meshtying");
    EXPECT_EQ(INPAR::STR::ModelTypeString(INPAR::STR::model_springdashpot), "SpringDashpot");
    EXPECT_EQ(
        INPAR::STR::ModelTypeString(INPAR::STR::model_lag_pen_constraint), "Lag-Pen-Constraint");
    EXPECT_EQ(INPAR::STR::ModelTypeString(INPAR::STR::model_basic_coupling), "Basic-Coupling");
    EXPECT_EQ(
        INPAR::STR::ModelTypeString(INPAR::STR::model_monolithic_coupling), "Monolithic-Coupling");
    EXPECT_EQ(INPAR::STR::ModelTypeString(INPAR::STR::model_partitioned_coupling),
        "Partitioned-Coupling");
    EXPECT_EQ(
        INPAR::STR::ModelTypeString(INPAR::STR::model_beam_interaction_old), "BeamInteractionOld");
    EXPECT_EQ(INPAR::STR::ModelTypeString(INPAR::STR::model_browniandyn), "BrownianDyn");
    EXPECT_EQ(INPAR::STR::ModelTypeString(INPAR::STR::model_beaminteraction), "BeamInteraction");
  }
}  // namespace
