
blimpy Maintenance & Regression Testing
=======================================


### Introduction

The purpose of the regression testing suite is to exercise and validate results from blimpy functional modules.  This is important in order to minimize potential inadvertent breakage when new development has occured. It is always best to catch bugs as soon as possible after they are introduced.
<br><br>
The primary method of launching regression testing is through the use of the `pytest` executable.  This is invoked in the following ways:
* Manually by a developer, on the command line in a terminal window.  This would follow downloading blimpy and setting up the development/testing environment (discussed later). 
* Automatically as part of a Github Pull Request (PR) after finalizing a fork of blimpy.
* Automatically as part of a Github Merge after a PR is approved.
<br>

### Development/Test Preparation

* The development of an amendment to `blimpy` begins with taking a fork from a github site, normally from `https://github.com/UCBerkeleySETI/blimpy`.
* Also, from the same site, `blimpy` is downloaded to a local computer.  The download operations can be performed in a few different ways but the simplest might be to download the zip file by clicking on the `Code` button and selecting `Download ZIP`.  Once the the zip file is in a local directory, unzip it and move the blimpy directory tree to wherever is appropriate for testing.  The zip file can now be discarded.
* Change directory into the `tests` directory (where this file is located) and execute `bash download_data.sh` which will perform all required regression testing initialization.
* When the previous step has completed, change directory up one level to the top of the `blimpy` directory tree.
* Execute: ```python3 setup.py install```
* Then, install `pytest` and `pyslalib` from pypi.org: `python3 -m pip install pytest pyslalib`.

### Regression Test Operations

* Running the full suite of regression tests is invoked by executing `pytest` with no parameters specified.  It is possible to run a single regression test file by specifying it as an argument to `pytest`.  For example, if one wishes to only run the plotting tests, the following is the command line to use: `pytest tests/test_plotting.py`.
* It is **highly encouraged** for developers to perform regression testing frequently in order to avoid surprises later on.
* Once, development activity on the local machine is complete and the last regression test has run verifying the absence of negative side effects, then the new and/or modified blimpy files can be uploaded to the developer's fork github site.
* At the fork github site, the developer can request a pull clicking on the `Pull request` button.  This automatically starts the PR process mentioned in the introduction section.
