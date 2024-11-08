// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include <gtest/gtest.h>

#include "4C_comm_pack_helpers.hpp"

using namespace FourC;

TEST(PackUnpackTest, Strings)
{
  std::string data = "Hello, World!";
  Core::Communication::PackBuffer pack_buffer;
  Core::Communication::add_to_pack(pack_buffer, data);

  std::string unpacked_data;
  Core::Communication::UnpackBuffer unpack_buffer(pack_buffer());
  Core::Communication::extract_from_pack(unpack_buffer, unpacked_data);

  EXPECT_EQ(data, unpacked_data);
}


TEST(PackUnpackTest, Map)
{
  using MapType = std::map<std::string, bool>;
  const MapType data{{"1", true}, {"key", false}};

  Core::Communication::PackBuffer pack_buffer;
  Core::Communication::add_to_pack(pack_buffer, data);

  MapType unpacked_data;
  Core::Communication::UnpackBuffer unpack_buffer(pack_buffer());
  Core::Communication::extract_from_pack(unpack_buffer, unpacked_data);

  EXPECT_EQ(data, unpacked_data);
}
