/*----------------------------------------------------------------------*/
/*! \file
\brief Unit tests for input parameters of structure field.

\level 0
*/
/*----------------------------------------------------------------------*/
#include <gtest/gtest.h>

#include "4C_inpar_structure.hpp"

namespace
{
  TEST(InparStructureTest, String2ModelType)
  {
    using namespace FourC;

    EXPECT_EQ(Inpar::STR::String2ModelType("Structure"), Inpar::STR::model_structure);
    EXPECT_EQ(Inpar::STR::String2ModelType("Contact"), Inpar::STR::model_contact);
    EXPECT_EQ(Inpar::STR::String2ModelType("Meshtying"), Inpar::STR::model_meshtying);
    EXPECT_EQ(Inpar::STR::String2ModelType("SpringDashpot"), Inpar::STR::model_springdashpot);
    EXPECT_EQ(
        Inpar::STR::String2ModelType("Lag-Pen-Constraint"), Inpar::STR::model_lag_pen_constraint);
    EXPECT_EQ(Inpar::STR::String2ModelType("Basic-Coupling"), Inpar::STR::model_basic_coupling);
    EXPECT_EQ(
        Inpar::STR::String2ModelType("Monolithic-Coupling"), Inpar::STR::model_monolithic_coupling);
    EXPECT_EQ(Inpar::STR::String2ModelType("Partitioned-Coupling"),
        Inpar::STR::model_partitioned_coupling);
    EXPECT_EQ(
        Inpar::STR::String2ModelType("BeamInteractionOld"), Inpar::STR::model_beam_interaction_old);
    EXPECT_EQ(Inpar::STR::String2ModelType("BrownianDyn"), Inpar::STR::model_browniandyn);
    EXPECT_EQ(Inpar::STR::String2ModelType("BeamInteraction"), Inpar::STR::model_beaminteraction);
  }

  TEST(InparStructureTest, ModelTypeString)
  {
    using namespace FourC;

    EXPECT_EQ(Inpar::STR::ModelTypeString(Inpar::STR::model_structure), "Structure");
    EXPECT_EQ(Inpar::STR::ModelTypeString(Inpar::STR::model_contact), "Contact");
    EXPECT_EQ(Inpar::STR::ModelTypeString(Inpar::STR::model_meshtying), "Meshtying");
    EXPECT_EQ(Inpar::STR::ModelTypeString(Inpar::STR::model_springdashpot), "SpringDashpot");
    EXPECT_EQ(
        Inpar::STR::ModelTypeString(Inpar::STR::model_lag_pen_constraint), "Lag-Pen-Constraint");
    EXPECT_EQ(Inpar::STR::ModelTypeString(Inpar::STR::model_basic_coupling), "Basic-Coupling");
    EXPECT_EQ(
        Inpar::STR::ModelTypeString(Inpar::STR::model_monolithic_coupling), "Monolithic-Coupling");
    EXPECT_EQ(Inpar::STR::ModelTypeString(Inpar::STR::model_partitioned_coupling),
        "Partitioned-Coupling");
    EXPECT_EQ(
        Inpar::STR::ModelTypeString(Inpar::STR::model_beam_interaction_old), "BeamInteractionOld");
    EXPECT_EQ(Inpar::STR::ModelTypeString(Inpar::STR::model_browniandyn), "BrownianDyn");
    EXPECT_EQ(Inpar::STR::ModelTypeString(Inpar::STR::model_beaminteraction), "BeamInteraction");
  }
}  // namespace
