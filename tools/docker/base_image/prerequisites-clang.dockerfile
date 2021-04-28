# Copyright (c) 2018-2021 Parsa Amini
#
# Distributed under the Boost Software License, Version 1.0. (See accompanying
# file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
#
# SYNOPSIS
#
# docker build [-t tag] [-build-arg LLVM_RELEASE=11|12|13]
#   [-build-arg BUILD_TYPE=Debug,Release,RelWithDebInfo]
#   [-build-arg CMAKE_VERSION=version] [...]
#
# DESCRIPTION
# This is a Docker file that can be used for building Octo-Tiger on CircleCI.
# It be configured to use Clang 11, 12, or 13, a desired CMake version, with
# desired build types for Boost, Vc, and HPX

ARG DEBIAN_RELEASE=buster
FROM debian:$DEBIAN_RELEASE

ARG LLVM_RELEASE=13
ARG BUILD_TYPE=Release
ARG HPX_BRANCH=1.6
ARG BOOST_VERSION=1.73.0
ARG CMAKE_VERSION=3.17.0
ARG HWLOC_VERSION=2.4.1

RUN apt-get -qq update \
    && apt-get -q install --yes \
        libgoogle-perftools-dev \
        ninja-build \
        ca-certificates \
        curl \
        wget \
        lsb-release \
        software-properties-common \
        git \
        gnupg2 \
    && curl -JL https://apt.llvm.org/llvm-snapshot.gpg.key | apt-key add - \
    && curl -JL https://apt.llvm.org/llvm.sh | bash -s -- ${LLVM_RELEASE} \
    && apt-get -qq install -y clang-format-${LLVM_RELEASE} clang-tidy-${LLVM_RELEASE} libc++-${LLVM_RELEASE}-dev libc++abi-${LLVM_RELEASE}-dev \
    && rm -rf /var/lib/apt/lists/*

ENV PATH=/local/silo/bin:/local/hdf5/bin:/local/hpx/bin:/llvm/bin:$PATH \
    LD_LIBRARY_PATH=/local/silo/lib:/local/hdf5/lib:/local/boost/lib:/local/vc/lib:/local/hpx/lib:/llvm/lib

RUN curl -JLO https://github.com/Kitware/CMake/releases/download/v${CMAKE_VERSION}/cmake-${CMAKE_VERSION}-Linux-x86_64.sh \
    && bash cmake-${CMAKE_VERSION}-Linux-x86_64.sh --prefix=/usr --exclude-subdir \
    && rm -f cmake-${CMAKE_VERSION}-Linux-x86_64.sh

RUN curl -JL https://download.open-mpi.org/release/hwloc/v${HWLOC_VERSION%.*}/hwloc-${HWLOC_VERSION}.tar.gz \
        | tar xz \
    && ( \
        cd hwloc-${HWLOC_VERSION} \
        && ./configure --prefix=/local/hwloc \
        && make -j22 && make install \
    ) \
    && rm -rf hwloc-${HWLOC_VERSION}

# --slave   /usr/bin/clang++ clang++ /usr/bin/clang++-${LLVM_RELEASE}  \
RUN update-alternatives \
        --install /usr/bin/llvm-config       llvm-config      /usr/bin/llvm-config-${LLVM_RELEASE} 100 \
        --slave   /usr/bin/llvm-ar           llvm-ar          /usr/bin/llvm-ar-${LLVM_RELEASE} \
        --slave   /usr/bin/llvm-as           llvm-as          /usr/bin/llvm-as-${LLVM_RELEASE} \
        --slave   /usr/bin/llvm-bcanalyzer   llvm-bcanalyzer  /usr/bin/llvm-bcanalyzer-${LLVM_RELEASE} \
        --slave   /usr/bin/llvm-cov          llvm-cov         /usr/bin/llvm-cov-${LLVM_RELEASE} \
        --slave   /usr/bin/llvm-diff         llvm-diff        /usr/bin/llvm-diff-${LLVM_RELEASE} \
        --slave   /usr/bin/llvm-dis          llvm-dis         /usr/bin/llvm-dis-${LLVM_RELEASE} \
        --slave   /usr/bin/llvm-dwarfdump    llvm-dwarfdump   /usr/bin/llvm-dwarfdump-${LLVM_RELEASE} \
        --slave   /usr/bin/llvm-extract      llvm-extract     /usr/bin/llvm-extract-${LLVM_RELEASE} \
        --slave   /usr/bin/llvm-link         llvm-link        /usr/bin/llvm-link-${LLVM_RELEASE} \
        --slave   /usr/bin/llvm-mc           llvm-mc          /usr/bin/llvm-mc-${LLVM_RELEASE} \
        --slave   /usr/bin/llvm-nm           llvm-nm          /usr/bin/llvm-nm-${LLVM_RELEASE} \
        --slave   /usr/bin/llvm-objdump      llvm-objdump     /usr/bin/llvm-objdump-${LLVM_RELEASE} \
        --slave   /usr/bin/llvm-ranlib       llvm-ranlib      /usr/bin/llvm-ranlib-${LLVM_RELEASE} \
        --slave   /usr/bin/llvm-readobj      llvm-readobj     /usr/bin/llvm-readobj-${LLVM_RELEASE} \
        --slave   /usr/bin/llvm-rtdyld       llvm-rtdyld      /usr/bin/llvm-rtdyld-${LLVM_RELEASE} \
        --slave   /usr/bin/llvm-size         llvm-size        /usr/bin/llvm-size-${LLVM_RELEASE} \
        --slave   /usr/bin/llvm-stress       llvm-stress      /usr/bin/llvm-stress-${LLVM_RELEASE} \
        --slave   /usr/bin/llvm-symbolizer   llvm-symbolizer  /usr/bin/llvm-symbolizer-${LLVM_RELEASE} \
        --slave   /usr/bin/llvm-tblgen       llvm-tblgen      /usr/bin/llvm-tblgen-${LLVM_RELEASE} \
    && update-alternatives \
        --install /usr/bin/clang                 clang                 /usr/bin/clang-${LLVM_RELEASE} 100 \
        --slave   /usr/bin/asan_symbolize        asan_symbolize        /usr/bin/asan_symbolize-${LLVM_RELEASE} \
        --slave   /usr/bin/c-index-test          c-index-test          /usr/bin/c-index-test-${LLVM_RELEASE} \
        --slave   /usr/bin/clang-check           clang-check           /usr/bin/clang-check-${LLVM_RELEASE} \
        --slave   /usr/bin/clang-cl              clang-cl              /usr/bin/clang-cl-${LLVM_RELEASE} \
        --slave   /usr/bin/clang-cpp             clang-cpp             /usr/bin/clang-cpp-${LLVM_RELEASE} \
        --slave   /usr/bin/clang-format          clang-format          /usr/bin/clang-format-${LLVM_RELEASE} \
        --slave   /usr/bin/clang-format-diff     clang-format-diff     /usr/bin/clang-format-diff-${LLVM_RELEASE} \
        --slave   /usr/bin/clang-include-fixer   clang-include-fixer   /usr/bin/clang-include-fixer-${LLVM_RELEASE} \
        --slave   /usr/bin/clang-offload-bundler clang-offload-bundler /usr/bin/clang-offload-bundler-${LLVM_RELEASE} \
        --slave   /usr/bin/clang-query           clang-query           /usr/bin/clang-query-${LLVM_RELEASE} \
        --slave   /usr/bin/clang-rename          clang-rename          /usr/bin/clang-rename-${LLVM_RELEASE} \
        --slave   /usr/bin/clang-reorder-fields  clang-reorder-fields  /usr/bin/clang-reorder-fields-${LLVM_RELEASE} \
        --slave   /usr/bin/clang-tidy            clang-tidy            /usr/bin/clang-tidy-${LLVM_RELEASE} \
        --slave   /usr/bin/lldb                  lldb                  /usr/bin/lldb-${LLVM_RELEASE} \
        --slave   /usr/bin/lldb-server           lldb-server           /usr/bin/lldb-server-${LLVM_RELEASE} \
    && update-alternatives --install /usr/bin/ld ld /usr/bin/ld.lld-${LLVM_RELEASE} 100 \
    && rm -f /usr/bin/clang++ \
    && printf "#!/bin/sh\nexec /usr/bin/clang++-${LLVM_RELEASE} -stdlib=libc++ \${@}" >/usr/bin/clang++ \
    && chmod +x /usr/bin/clang++

RUN git clone https://github.com/HDFGroup/hdf5.git --depth=1 --branch=hdf5-1_10_4 \
    && cmake -Shdf5 -Bhdf5/build \
        -DBUILD_TESTING=OFF \
        -DHDF5_BUILD_CPP_LIB=OFF \
        -DCMAKE_INSTALL_PREFIX=/local/hdf5 \
        -DCMAKE_BUILD_TYPE=Release \
        -GNinja \
        -Wno-dev \
    && cmake --build hdf5/build --target install \
    && rm -rf hdf5

RUN curl -JL https://wci.llnl.gov/sites/wci/files/2021-01/silo-4.10.2.tgz \
        | tar xz \
    && ( \
        cd silo-4.10.2 \
        && ./configure --disable-fortran --prefix=/local/silo \
            --with-hdf5=/local/hdf5/include,/local/hdf5/lib --enable-optimization \
        && make -j22 install \
    ) \
    && rm -rf silo-4.10.2

#    update-alternatives --install /usr/bin/cc  cc  /usr/bin/clang 100 \
# && update-alternatives --install /usr/bin/c++ c++ /usr/bin/clang++ 100

RUN curl -JL http://downloads.sourceforge.net/project/boost/boost/${BOOST_VERSION}/boost_$(echo $BOOST_VERSION | tr . _).tar.gz \
        | tar xz \
    && ( \
        cd boost_$(echo $BOOST_VERSION | tr . _) \
        && ./bootstrap.sh --prefix=/local/boost --with-toolset=clang \
        && ./b2 -j22 --with-iostreams --with-atomic --with-filesystem --with-program_options \
            --with-regex --with-system --with-chrono --with-date_time \
            --with-thread $(echo ${BUILD_TYPE/%WithDebInfo/ease} | tr '[:upper:]' '[:lower:]') \
            toolset=clang \
            install \
    ) \
    && rm -rf boost_$(echo $BOOST_VERSION | tr . _)

RUN git clone https://github.com/VcDevel/Vc.git --depth=1 --branch=1.4.1 \
    && cmake -SVc -BVc/build \
        -DBUILD_TESTING=OFF \
        -DCMAKE_INSTALL_PREFIX=/local/vc \
        -DCMAKE_BUILD_TYPE=${BUILD_TYPE} \
        -GNinja \
    && cmake --build Vc/build --target install \
    && rm -rf Vc

RUN git clone https://github.com/STEllAR-GROUP/hpx.git --depth=1 --branch=${HPX_BRANCH} \
    && cmake -Shpx -Bhpx/build \
        -DBOOST_ROOT=/local/boost \
        -DHPX_WITH_EXAMPLES=OFF \
        -DHWLOC_ROOT=/local/hwloc \
        -DCMAKE_INSTALL_PREFIX=/local/hpx \
        -DCMAKE_BUILD_TYPE=${BUILD_TYPE} \
        -GNinja \
    && cmake --build hpx/build --target install \
    && rm -rf hpx
