# MLOC

![AUR](https://img.shields.io/badge/version-v10.4.5-brightgreen.svg)
![AUR](https://img.shields.io/badge/release-10%2F15%2F2018-orange.svg)
[![AUR](https://img.shields.io/badge/license-GPLv3-blue.svg)](https://www.gnu.org/licenses/quick-guide-gplv3.en.html)
[![AUR](https://img.shields.io/badge/Docs-latest-48C9B0.svg)](https://github.com/shipengcheng1230/MLOC/tree/master/Docs)


## Introduction

**MLOC** is a multiple-events hypercentroid decomposition relocation program developed and maintained by [Eric Bergman](https://www.researchgate.net/profile/Eric_Bergman2) since 1990 whose underlying idea is based on [Jordan and Sverdrup (1981)](https://pubs.geoscienceworld.org/ssa/bssa/article-abstract/71/4/1105/102070/teleseismic-location-techniques-and-their). This tutorial is a skeleton iteration based on the **MLOC** Training Workshop organized by Eric and [Ezgi](https://inside.mines.edu/~ekarasoz/) that was held in Albuquerque, NM from 10/15/18 to 10/19/18. Further update of this software will be reflected upon Eric's notification.

## Distribution

Eric has never personally distributed this software community-widely but he is fine with redistributing it by those who attended his **MLOC** workshop. Therefore, I keep the track of this software under my personal private repository and impose *GPLv3* license here. Those who share this repository shall respect it. 

## Setup

Before start, some prerequisites, including [GMT5](http://gmt.soest.hawaii.edu/projects/gmt) and a Fortran compiler, must be met, for instance:

```console
➜  MLOC git:(master) ✗ gmt --version
5.4.4

➜  MLOC git:(master) ✗ gfortran --version
GNU Fortran (Homebrew GCC 8.2.0) 8.2.0
Copyright (C) 2018 Free Software Foundation, Inc.
This is free software; see the source for copying conditions.  There is NO
warranty; not even for MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
```

First clone this repository to your local:

```Shell
git clone https://github.com/shipengcheng1230/MLOC.git
```

Then download the tables file [here](https://drive.google.com/drive/folders/15Vr0Gi_0WSK73DNHGmBA49EgFkn84YZl?usp=sharing) and put it under `MLOC/mloc_working/tables`. Your directory should look like this:

```console
├── Docs
│   ├── Focal Depth
│   ├── HD calibration cartoons
│   ├── Hi-Rez topographic datasets.rtfd
│   └── Papers
├── MNF utilities
│   ├── mnf search
│   └── to mnf
├── clusters
│   └── Tunisia
├── mloc utilities
│   ├── lres
│   ├── null_test
│   ├── rstat
│   └── xdat
├── mloc_gfortran
├── mloc_intel
├── mloc_src
└── mloc_working
    ├── Utilities
    ├── jsa4
    └── tables
```

A brief summary of this program is as below:
- `Docs` contains all the related references and manuals. 
- `MNF utilities` contains functionalities for transforming ISC bulletin to `.mnf` format. 
- `mloc utilities` contains functionalities for **MLOC** processing procedures. 
- `mloc_gfortran` and `mloc_intel` contain `makefile` for building `mloc` program .
- `mloc_src` contains the source code.
- `mloc_working` will be your primary working directory, in which `table` has all the geological models used by this program.

Now compile source code and any other utilities, for example using `gfortran`:

```Shell
cd mloc_gfortran && make && cp ./mloc_g ../mloc_working
```

Modify `mloc.conf` such that `WORKING_DIR` points to your `mloc_working` directory and `AUTHOR` to your name not exceeding 8 characters (for more configuration see docs):

```
WORKING_DIR: /Users/spc/Softwares/mloc/mloc_working
AUTHOR: SPC
```

Now it should not have problem running the program:

```console
➜  mloc_working git:(master) ✗ ./mloc_g 

mloc v10.4.5, release date 9/28/2018                                            

Current program limits: 
  nevmax =   200
  nqmax  =  4000
  ntmax1 = 35000

Enter a basename for this run: 
```

If you cannot come to this step, please pull up an issue.

## Quick Start

### 1. Download Data

We will use [ISC bulletin](http://www.isc.ac.uk/iscbulletin/search/bulletin/) data which contain all the possible phases identified by the agency. It is also desired to use own waveform data to pick phases manually to proliferate our data coverage but we will not cover this here for now.

![ISC Search](https://github.com/shipengcheng1230/MLOC/blob/master/RDfigures/ISC_Search.png)

The database contains two type: *Reviewed ISC Bulletin* and *ISC Bulletin* where the first one are those reviewed for correction. For output format we will choose *ISF Bulletin*. *QuakeML* cannot work with **MLOC** for now. Within the search region box, it falls on your own choice to draw the region you prefer. After clicking *Search bulletin*, you will be redirected towards the results page. Copy and paste all the contents and save to your local in which the file looks like:

```
International Seismological Centre
ISC: On-Line Bulletin
Any use of data from the ISC should be cited. The correct format for citations may be found on our citation page.

The ISC Bulletin has been rebuilt for the period 1964-1979. All ISC searches will now return the upgraded set of data for this period. The work on the Rebuild project continues for the period 1980-2010 and will result in further gradual bulletin updates in due time.

Once the search has completed, a compressed KML file will be available to view the results in Google Earth.

Make an event map

Search summary:
Database: ISC Bulletin
Search type: Circular search
Central latitude: 34.1
Central longtitude: 9.9
Radius: 200 km
Start date: 1960-01-01 00:00:00
End date: 2018-10-19 00:00:00
Events found: 215
DATA_TYPE BULLETIN IMS1.0:short
ISC Bulletin
Event   876000 Tunisia
   Date       Time        Err   RMS Latitude Longitude  Smaj  Smin  Az Depth   Err Ndef Nsta Gap  mdist  Mdist Qual   Author      OrigID
1961/01/21 03:45:25                  35.2500   10.5000                                                             uk BCIS       1901692
```

For simplicity, you will find a piece of sample data named `test.dat` in `MLOC/clusters/Tunisia/Test/Data` that centered somewhere in Tunisia with radius of 200 km dated from 1960 to 2018. Then, use the utility named `isc_ims2mnf` to obtain a `.mnf` style data:

```console

```

### 2. Convert to `.mnf` format

### 3. Fulfill command file

### 4. Run the program

### 5. Check the results

### 6. Trim the data

### 7. Second run

### 8. Modify crustal model


## Additional Commands

### 1. Stations
It is also worth mentioning the [International Registry of Seismograph Stations (IR)](http://www.isc.ac.uk/registries/). Currently the update-to-date list is stored in `MLOC/mloc_working/tables/stn/master_stn.dat` which looks like:

```
# Deployment codes
# INSN: Iran National Seismograph Network
# ITSN: Iran Telemetered Seismograph Network
# Kharkheh: Kharkheh Dam Network, Iran
# TA: Transportable Array <http://ds.iris.edu/ds/nodes/dmc/earthscope/usarray/#transportable_array>
#
034A   27.0647   -98.6833    155      IR       IRIS .TA      .                   Hebbronville, Texas, U.S.A.
035A   26.9379   -98.1023     29      IR       IRIS .TA      .                   Encino, Texas, U.S.A.
035Z   26.46300  -98.06831    19      IR       IRIS .TA      .                   Hargill, Texas, U.S.A.
```

To supplement your own stations, put them in 

### 2. Indirect calibration

### 3. Phase identification

### 4. More

Use the `help` command within `mloc` program to view more:

```Shell
The commands are:
  anno
  bdps bias bloc bptc
  cal  ccat cfil clim comm corr cptf ctyp cvff cvou cvtt
  datf dbug dcal dem1 dem2 dep  diff
  ellp epap eplt even
  fdhp flag fmap fplt frec freh
  help hlim
  inpu
  kill
  lat  lgtt lmod long lonr lres
  mare mdou mech memb
  nsmd
  oldr
  pert phid phyp plot pltt ppri pttt puke
  radf rdpp rels revi rfil rhdf run 
  secv shcl skip splt sstn star stat step stop subc
  taup tfil tikh time tomo tptt tt5e tt5s ttou
  vect vlog vscr
  weig wind
  xsec

  For more information, follow the "help" command with or without a command name
```