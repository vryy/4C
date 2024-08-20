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

    EXPECT_EQ(Inpar::Solid::string_to_model_type("Structure"), Inpar::Solid::model_structure);
    EXPECT_EQ(Inpar::Solid::string_to_model_type("Contact"), Inpar::Solid::model_contact);
    EXPECT_EQ(Inpar::Solid::string_to_model_type("Meshtying"), Inpar::Solid::model_meshtying);
    EXPECT_EQ(
        Inpar::Solid::string_to_model_type("SpringDashpot"), Inpar::Solid::model_springdashpot);
    EXPECT_EQ(Inpar::Solid::string_to_model_type("Lag-Pen-Constraint"),
        Inpar::Solid::model_lag_pen_constraint);
    EXPECT_EQ(
        Inpar::Solid::string_to_model_type("Basic-Coupling"), Inpar::Solid::model_basic_coupling);
    EXPECT_EQ(Inpar::Solid::string_to_model_type("Monolithic-Coupling"),
        Inpar::Solid::model_monolithic_coupling);
    EXPECT_EQ(Inpar::Solid::string_to_model_type("Partitioned-Coupling"),
        Inpar::Solid::model_partitioned_coupling);
    EXPECT_EQ(Inpar::Solid::string_to_model_type("BeamInteractionOld"),
        Inpar::Solid::model_beam_interaction_old);
    EXPECT_EQ(Inpar::Solid::string_to_model_type("BrownianDyn"), Inpar::Solid::model_browniandyn);
    EXPECT_EQ(
        Inpar::Solid::string_to_model_type("BeamInteraction"), Inpar::Solid::model_beaminteraction);
  }

  TEST(InparStructureTest, ModelTypeString)
  {
    using namespace FourC;

    EXPECT_EQ(Inpar::Solid::model_type_string(Inpar::Solid::model_structure), "Structure");
    EXPECT_EQ(Inpar::Solid::model_type_string(Inpar::Solid::model_contact), "Contact");
    EXPECT_EQ(Inpar::Solid::model_type_string(Inpar::Solid::model_meshtying), "Meshtying");
    EXPECT_EQ(Inpar::Solid::model_type_string(Inpar::Solid::model_springdashpot), "SpringDashpot");
    EXPECT_EQ(Inpar::Solid::model_type_string(Inpar::Solid::model_lag_pen_constraint),
        "Lag-Pen-Constraint");
    EXPECT_EQ(
        Inpar::Solid::model_type_string(Inpar::Solid::model_basic_coupling), "Basic-Coupling");
    EXPECT_EQ(Inpar::Solid::model_type_string(Inpar::Solid::model_monolithic_coupling),
        "Monolithic-Coupling");
    EXPECT_EQ(Inpar::Solid::model_type_string(Inpar::Solid::model_partitioned_coupling),
        "Partitioned-Coupling");
    EXPECT_EQ(Inpar::Solid::model_type_string(Inpar::Solid::model_beam_interaction_old),
        "BeamInteractionOld");
    EXPECT_EQ(Inpar::Solid::model_type_string(Inpar::Solid::model_browniandyn), "BrownianDyn");
    EXPECT_EQ(
        Inpar::Solid::model_type_string(Inpar::Solid::model_beaminteraction), "BeamInteraction");
  }
}  // namespace
