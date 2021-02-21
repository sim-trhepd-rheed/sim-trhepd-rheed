# sim-trhepd-rheed

## Overview

The software 'sim-trhepd-rheed' is a simulator for total-reflection high-energy positron diffraction (TRHEPD) 
and reflection high-energy electron diffraction (RHEED).
Please cite Ref. [1], if you use the software in a research paper.

The software was written originally by Takashi Hanada and the present package on github is maintained by Takeo Hoshi.
The software is based on the theory written in Ref.[2].

Refernce list

[1] T. Hanada, H. Daimon, and S. Ino， Rocking-curve analysis of reflection high-energy electron diffraction from the Si(111)-(√3 × √3 )R30°-Al, -Ga, and -In surfaces,  Phys. Rev. B 51, 13320–13325 (1995).

[2] A. Ichimiya, Many-beam calculation of reflection high energy electron diffraction (RHEED) intensities by the multi-slice method, Jpn. J. Appl. Phys. 22, 176-180 (1983).

## Preparation 

### Download the source code

Download the source code with git command in a  working directory

$ git clone http://github.com/sim-trhepd-rheed/sim-trhepd-rheed

The following directory will be generated

sim-trhepd-rheed

### Building the executable files 

In the sim-trhepd-rheed directory, execute the command below in order.

$ cd src

$ make

Then the following files will be generated 

bulk.exe

surf.exe

As an option, when you build atom coordinate visioning tool, execute the command below. 

$ make xyz

Then the following files will be generated 

xyzb.exe


