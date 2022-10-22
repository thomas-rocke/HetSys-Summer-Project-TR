# Installing QUIP & quippy

quippy (Python wrapped QUIP) is the amin means of using GAP as an ASE calculator. To get the most recent version of QUIP, the program needs to be installed from source.


# Instructions (assuming no previous QUIP install)
A more general version of install instructions can be found at: https://github.com/libAtoms/QUIP/blob/public/README.md


Fetch the QUIP repo using:

git clone --recursive https://github.com/libAtoms/QUIP.git

cd into the generated dir (cd QUIP)

Set the following environment variable to tell QUIP the desired architecture to install for (others are available, look in the QUIP/arch for options):

export QUIP_ARCH=linux_x86_64_gfortran

Then run the following command to configure the QUIP install:

make config

Most of the defaults should be reasonable for the install we need. The main thing to change are the questions on the python and pip paths - setting these to "python3" and "pip3" respectively should ensure a correct install.

Once QUIP is correctly configured, we then need to build and install it. This can be done by the following command:

make && make install-quippy

This should install quippy as a python module, which should then give access to the generic quippy.potential.Potential object we use to test GAP models.
