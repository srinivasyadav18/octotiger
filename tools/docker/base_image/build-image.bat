:: Copyright (c) 2018-2021 Parsa Amini
::
:: Distributed under the Boost Software License, Version 1.0. (See accompanying
:: file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

:: Description
:: This Batch script is used to build the Docker images used to build and test
:: Octo-Tiger on CircleCI

:: Useful docker build diagnostics argument(s):
::     --no-cache    Re-run every stage

::docker build -t stellargroup/octotiger:prerequisites-gcc8 -f prerequisites-gcc.dockerfile --build-arg GCC_RELEASE=8 .
docker build -t stellargroup/octotiger:prerequisites-gcc8 -f prerequisites-gcc-split.dockerfile --build-arg GCC_RELEASE=8 .
::docker build -t stellargroup/octotiger:prerequisites-gcc10 -f prerequisites-gcc.dockerfile --build-arg GCC_RELEASE=10 .
::docker build -t stellargroup/octotiger:prerequisites-gcc10-debug -f prerequisites-gcc.dockerfile --build-arg GCC_RELEASE=10 --build-arg BUILD_TYPE=Debug .
::docker build -t stellargroup/octotiger:prerequisites-gcc10-relwithdebinfo -f prerequisites-gcc.dockerfile --build-arg GCC_RELEASE=10 --build-arg BUILD_TYPE=RelWithDebInfo .

::docker build -t stellargroup/octotiger:prerequisites-clang13-debug -f prerequisites-clang.dockerfile --build-arg DEBIAN_RELEASE=buster --build-arg LLVM_RELEASE=13 --build-arg BUILD_TYPE=Debug .
::docker build -t stellargroup/octotiger:prerequisites-clang13 -f prerequisites-clang.dockerfile --build-arg DEBIAN_RELEASE=buster --build-arg LLVM_RELEASE=13 --build-arg BUILD_TYPE=Release .

docker build -t stellargroup/octotiger:vscode-gcc8 -f vscode-gcc8.dockerfile --build-arg GCC_RELEASE=8 .

::docker push stellargroup/octotiger:prerequisites-gcc8
::docker push stellargroup/octotiger:vscode-gcc8
::docker push stellargroup/octotiger:prerequisites-clang13
pause