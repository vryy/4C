# This file is part of 4C multiphysics licensed under the
# GNU Lesser General Public License v3.0 or later.
#
# See the LICENSE.md file in the top-level for license information.
#
# SPDX-License-Identifier: LGPL-3.0-or-later

FROM ghcr.io/4c-multiphysics/4c-dependencies-ubuntu24.04:latest
LABEL org.opencontainers.image.description="Image containing the built 4C code"
LABEL org.4c-multiphysics.project=4C

RUN git clone https://github.com/4C-multiphysics/4C.git

RUN mkdir 4C/build && cd 4C/build && \
  cmake ../ --fresh --preset=docker && \
  cmake --build . --target full -- -j `nproc`
