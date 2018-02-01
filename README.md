# RAJA Proxy Applications

This repo contains a selection of proxy applications written using the RAJA
performance portability layer. These applications are examples of how RAJA can
be used in real codes, as well as providing a convenient vehicle for testing
features and analyzing performance.

For more details about the concepts and design behind RAJA, please take a look
at the GitHub [repository](https://github.com/LLNL/RAJA).

## Quick Start

This repository is hosted on [GitHub](https://github.com/LLNL/RAJAProxies).  To
clone the repo into your local working space, use the command:

    $ git clone --recursive https://github.com/LLNL/RAJAProxies.git 

The `--recursive` argument is used to download the repository's submodules, RAJA
and [BLT](https://github.com/LLNL/blt).  Once this command has run, the `master`
branch will be cloned into the `RAJA-Proxies` directory.  RAJA uses CMake and
BLT to handle builds. From the root of the repository, configuration looks like:

    $ mkdir build && cd build
    $ cmake ../

CMake will provide output about the compiler that has been detected, and what
features are discovered. Some features, like OpenMP, will be enabled if they are
discovered.  Once CMake has completed, the applications can be compiled using
Make:

    $ make

For more details about enabling other features, please see the user
documentation.

Each application will be built in the `bin` directory. The binary names will
include the version and features of each application. To run the application,
simply run the relevant executable. For example:

    $ ./bin/lulesh-v1.0-RAJA-OMP.exe

would run the RAJA version of LULESH v1.0 with the OpenMP backend.

For more detail on the input parameters supported by each application, take a
look at the user documentation.

## User Documentation

The [documentation]() is hosted on ReadTheDocs, and is the best place to learn
more about these applications.

## Questions?

If you have any questions about this repo, please contact raja-dev@llnl.gov.

## Authors

This repository was developed by:

* David Beckingsale (david@llnl.gov)
* Rich Hornung (hornung1@llnl.gov)

Additional application developers include:

* Jeff Keasler (keasler1@llnl.gov)

## Release

This repository is part of the RAJA project, released as:

```
Copyright (c) 2017-2018, Lawrence Livermore National Security, LLC.

Produced at the Lawrence Livermore National Laboratory.

All rights reserved.

LLNL-CODE-689114  OCEC-16-063

Unlimited Open Source - BSD Distribution
```

The include applications have their own licenses and release numbers, included
in the appropriate subdirectory.
