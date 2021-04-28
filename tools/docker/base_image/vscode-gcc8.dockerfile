ARG GCC_RELEASE=8
FROM stellargroup/octotiger:prerequisites-gcc${GCC_RELEASE}

RUN git clone https://github.com/stellar-group/octotiger.git --depth=1 --branch=master /octotiger

RUN mkdir -p /root/.local/share/CMakeTools/ \
    && echo '[{"name":"GCC","compilers":{"C":"/usr/local/bin/gcc","CXX":"/usr/local/bin/c++"},"keep":true}]' >/root/.local/share/CMakeTools/cmake-tools-kits.json \
    && mkdir -p /octotiger/.vscode \
    && echo '{"configurations":[{"name":"Linux","includePath":["${workspaceFolder}/**","/local/hpx/include/","/local/boost/include/","/local/hdf5/include","/local/hwloc/include","/local/silo/include/","/local/vc/include/"],"defines":[],"compilerPath":"/usr/local/bin/gcc","cStandard":"gnu17","cppStandard":"gnu++14","intelliSenseMode":"linux-gcc-x64","configurationProvider":"ms-vscode.cmake-tools"}],"version":4}' >/octotiger/.vscode/c_cpp_properties.json \
    && echo '{"cmake.configureSettings":{"HPX_DIR":"/local/hpx/lib/cmake/HPX","Vc_DIR":"/local/vc/lib/cmake/Vc","Silo_DIR":"/local/silo","HDF5_ROOT":"/local/hdf5","BOOST_ROOT":"/local/boost","OCTOTIGER_WITH_DOCU":"ON"}}' >/octotiger/.vscode/settings.json
