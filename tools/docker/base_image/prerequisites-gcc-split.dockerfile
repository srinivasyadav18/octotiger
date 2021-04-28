# Copyright (c) 2021 Parsa Amini
#
# Distributed under the Boost Software License, Version 1.0. (See accompanying
# file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
#
# SYNOPSIS
#
# docker build [-t tag] [-build-arg GCC_RELEASE=8|9|10]
#   [-build-arg BUILD_TYPE=Debug,Release,RelWithDebInfo]
#   [-build-arg CMAKE_VERSION=version] [...]
#
# DESCRIPTION
# This is a Docker file that is used for building Octo-Tiger on CircleCI. It can
# be configured to use GCC 8, 9, or 10 a desired CMake version, with desired
# build types for Boost, Vc, and HPX

ARG GCC_RELEASE=8

FROM gcc:${GCC_RELEASE} AS hwloc-stage
ARG HWLOC_VERSION=2.4.1

RUN curl -JL https://download.open-mpi.org/release/hwloc/v${HWLOC_VERSION%.*}/hwloc-${HWLOC_VERSION}.tar.gz \
        | tar xz \
    && ( \
        cd hwloc-${HWLOC_VERSION} \
        && ./configure --prefix=/local/hwloc --enable-silent-rules \
        && make -j && make install \
    ) \
    && rm -rf hwloc-hwloc-${HWLOC_VERSION}



ARG BUILD_TYPE=Release
FROM gcc:${GCC_RELEASE} AS boost-stage
ARG BOOST_VERSION=1.73.0

RUN curl -JL http://downloads.sourceforge.net/project/boost/boost/${BOOST_VERSION}/boost_$(echo $BOOST_VERSION | tr . _).tar.gz \
        | tar xz \
    && ( \
        cd boost_$(echo $BOOST_VERSION | tr . _) \
        && ./bootstrap.sh --prefix=/local/boost \
        && ./b2 -j22 -d0 --with-iostreams  --with-atomic --with-filesystem --with-program_options \
            --with-regex --with-system --with-chrono --with-date_time \
            --with-thread $(echo ${BUILD_TYPE/%WithDebInfo/ease} | tr '[:upper:]' '[:lower:]') install \
    ) \
    && rm -rf boost_$(echo $BOOST_VERSION | tr . _)



FROM gcc:${GCC_RELEASE} AS silo-stage

ARG CMAKE_VERSION=3.17.0

RUN curl -JLO https://github.com/Kitware/CMake/releases/download/v${CMAKE_VERSION}/cmake-${CMAKE_VERSION}-Linux-x86_64.sh \
    && bash cmake-${CMAKE_VERSION}-Linux-x86_64.sh --prefix=/usr --exclude-subdir \
    && rm -f cmake-${CMAKE_VERSION}-Linux-x86_64.sh

RUN git clone https://github.com/HDFGroup/hdf5.git --depth=1 --branch=hdf5-1_10_4 \
    && cmake -Shdf5 -Bhdf5/build \
        -DBUILD_TESTING=OFF \
        -DHDF5_BUILD_CPP_LIB=OFF \
        -DCMAKE_INSTALL_PREFIX=/local/hdf5 \
        -DCMAKE_BUILD_TYPE=Release \
        -Wno-dev \
    && cmake --build hdf5/build --parallel \
    && cmake --install hdf5/build \
    && rm -rf hdf5

RUN curl -JL https://wci.llnl.gov/sites/wci/files/2021-01/silo-4.10.2.tgz \
        | tar xz \
    && ( \
        cd silo-4.10.2 \
        && ./configure --disable-fortran --prefix=/local/silo --enable-silent-rules \
            --with-hdf5=/local/hdf5/include,/local/hdf5/lib --enable-optimization \
        && make -j install \
    ) \
    && rm -rf silo-4.10.2


FROM gcc:${GCC_RELEASE}

COPY --from=hwloc-stage /local /local
COPY --from=boost-stage /local /local
COPY --from=silo-stage /local /local

ARG CMAKE_VERSION=3.17.0
ARG HPX_BRANCH=1.6.0
ARG BUILD_TYPE=Release

ARG DEBIAN_FRONTEND=noninteractive
RUN apt-get -qq update \
    && apt-get -qq install --yes --no-install-recommends \
        libgoogle-perftools-dev \
        ninja-build \
        doxygen \
        python3-pip \
        locales \
    && rm -rf /var/lib/apt/lists/* \
    && locale-gen en_US.UTF-8

RUN curl -JLO https://github.com/Kitware/CMake/releases/download/v${CMAKE_VERSION}/cmake-${CMAKE_VERSION}-Linux-x86_64.sh \
    && bash cmake-${CMAKE_VERSION}-Linux-x86_64.sh --prefix=/usr --exclude-subdir \
    && rm -f cmake-${CMAKE_VERSION}-Linux-x86_64.sh

RUN git clone https://github.com/VcDevel/Vc.git --depth=1 --branch=1.4.1 \
    && cmake -SVc -BVc/build \
        -DBUILD_TESTING=OFF \
        -DCMAKE_INSTALL_PREFIX=/local/vc \
        -DCMAKE_BUILD_TYPE=${BUILD_TYPE} \
        -GNinja \
    && cmake --build Vc/build --target install \
    && rm -rf Vc

RUN test -d /local/boost \
    && git clone https://github.com/STEllAR-GROUP/hpx.git --depth=1 --branch=${HPX_BRANCH} \
    && cmake -Shpx -Bhpx/build \
        -DBOOST_ROOT=/local/boost \
        -DHPX_WITH_EXAMPLES=OFF \
        -DHWLOC_ROOT=/local/hwloc \
        -DCMAKE_INSTALL_PREFIX=/local/hpx \
        -DCMAKE_BUILD_TYPE=${BUILD_TYPE} \
        -GNinja \
    && cmake --build hpx/build --target install \
    && rm -rf hpx

ENV PATH=/local/silo/bin:/local/hdf5/bin:/local/hpx/bin:$PATH \
    LD_LIBRARY_PATH=/local/silo/lib:/local/hdf5/lib:/local/boost/lib:/local/vc/lib:/local/hpx/lib
