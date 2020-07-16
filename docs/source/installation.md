# Installation

## Prerequisites

PyProp is dependent upon the AirfoilDatabase package, also distributed by the USU AeroLab. The source code for AirfoilDatabase can be found on our [Github page](https://github.com/usuaero/AirfoilDatabase). Please make sure you have the latest version of AirfoilDatabase downloaded and installed before installing PyProp.

## Getting Python

If you do not have Python installed on your machine, it can be downloaded from a number of locations. We use [https://www.anaconda.com/distribution/](https://www.anaconda.com/distribution/). Please be sure you have Python 3.6 or later.

## Getting the Source Code

You can either download the source as a ZIP file and extract the contents, or clone the PyProp repository using Git. If your system does not already have a version of Git installed, you will not be able to use this second option unless you first download and install Git. If you are unsure, you can check by typing `git --version` into a command prompt.

### Downloading source as a ZIP file (Not recommended)

1. Open a web browser and navigate to [https://github.com/usuaero/PyProp](https://github.com/usuaero/PyProp)
2. Make sure the Branch is set to `Master`
3. Click the `Clone or download` button
4. Select `Download ZIP`
5. Extract the downloaded ZIP file to a local directory on your machine

### Cloning the Github repository (Recommended)

1. From the command prompt navigate to the directory where PyProp will be installed. Note: git will automatically create a folder within this directory called PyProp. Keep this in mind if you do not want multiple nested folders called PyProp.
2. Execute

    $ git clone https://github.com/usuaero/PyProp

We recommend cloning from the repository, as this will allow you to most easily download and install periodic updates (because we will always be updating PyProp!). This can be done using the following command

    $ git pull

Please reinstall after pulling from the repository.

## Installing

Once you have the source code downloaded, navigate to the root (PyProp/) directory and execute

    $ pip install .

Please note that any time you update the source code (e.g. after executing a git pull), PyProp will need to be reinstalled by executing the above command.

## Testing the Installation

Once the installation is complete, run

    $ py.test

to verify PyProp is working properly on your machine. Some warnings and depreciation notices are normal.