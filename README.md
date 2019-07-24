# RAJA Proxy Applications

This project contains a collection of proxy applications written using the 
[RAJA Rerformance Portability Layer](https://github.com/LLNL/RAJA). 
These applications are examples of how RAJA is used in real codes and 
provideng a convenient vehicle for testing features and analyzing performance.

## Quick Start

This repository is hosted on [GitHub](https://github.com/LLNL/RAJAProxies).
To clone the repo into your local working space, use the command:

    $ git clone --recursive https://github.com/LLNL/RAJAProxies.git 

The `--recursive` argument is used to download the repository's submodules, 
RAJA and the build system [BLT](https://github.com/LLNL/blt). After you 
execute this command, you will see the `master` branch in the `raja-proxies` 
directory. 

Then, you can build RAJA and the proxy applications like any other CMake 
project, provided you have a C++ compiler that supports the C++11 standard. 
The simplest way to build the code is to do the following in the top-level 
`raja-proxies` directory (in-source builds are not allowed!):

    $ mkdir build
    $ cd build
    $ cmake ../
    $ make

More details about RAJA configuration options are located in the 
[**RAJA User Guide and Tutorial**](http://raja.readthedocs.io/en/master/).

The executable for Each application will be located in the `bin` directory
of your build space. The executable names will include the name of the proxy
app, its version and parallel programming model it is using. To run an
application, simply run the desired executable. For example, to run 
run the RAJA version of LULESH v1.0 with the OpenMP backend, execute the
following command:

    $ ./lulesh-v1.0-RAJA-omp.exe

## Proxy Application Information

Information about each available proxy application is available here
[RAJA_Proxy_Apps.md](RAJA_Proxy_Apps.md).

## Questions?

If you have any questions about this repo, please send email to 
raja-dev@llnl.gov or contact one of the individuals listed below.

## Authors

This repository is maintained by:

* David Beckingsale (david@llnl.gov)
* Rich Hornung (hornung1@llnl.gov)

## Release Information

Please see [RELEASE.md](./RELEASE.md) for release information for each proxy
application.
