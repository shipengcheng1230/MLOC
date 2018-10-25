# MLOC

[![AUR](https://img.shields.io/badge/version-v10.4.5-brightgreen.svg)](https://github.com/shipengcheng1230/MLOC/blob/master/mloc_src/mloc_version_history.txt)
![AUR](https://img.shields.io/badge/release-10%2F15%2F2018-orange.svg)
[![AUR](https://img.shields.io/badge/license-GPLv3-blue.svg)](https://www.gnu.org/licenses/quick-guide-gplv3.en.html)
[![AUR](https://img.shields.io/badge/Docs-latest-48C9B0.svg)](https://github.com/shipengcheng1230/MLOC/tree/master/Docs)


## Introduction

**MLOC** is a multiple-events hypercentroid decomposition relocation program developed and maintained by [Eric Bergman](https://www.researchgate.net/profile/Eric_Bergman2) since 1990 whose underlying idea is based on [Jordan and Sverdrup (1981)](https://pubs.geoscienceworld.org/ssa/bssa/article-abstract/71/4/1105/102070/teleseismic-location-techniques-and-their). This tutorial is a skeleton iteration based on the **MLOC** Training Workshop organized by Eric and [Ezgi](https://inside.mines.edu/~ekarasoz/) that was held in Albuquerque, NM from 10/15/18 to 10/19/18. Further update of this software will be reflected upon Eric's notification.

## Distribution

Eric has never personally distributed this software community-widely but he is fine with redistributing it by those who attended his **MLOC** workshop. Therefore, I keep the track of this software under my personal private repository and impose *GPLv3* license here.

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

For simplicity, you will find a piece of sample data named `test.dat` in `MLOC/clusters/Tunisia/Test/Data` that centered somewhere in Tunisia with radius of 200 km dated from 1960 to 2018.

### 2. Convert to `.mnf` format

Then, use the utility named `isc_ims2mnf` to obtain a `.mnf` style data:

```console
➜  Data git:(master) ✗ ./isc_ims2mnf 
Release date April 6, 2018, writing MNF v1.3.3 

Enter input filename: 
test.dat
Bulletin output?: 
y
Bulletin comment: 
Tunisia test data
EOF reached after    215 events
```

This will generate a single `.mnf` bulletin if you choose `y` for *Bulletin output?*, which looks like:

```
B   Tunisia test data                                                                                                    
F   MNF v1.3.3                                                                                                           
E   Tunisia                                                                                                              
I   ISC        876000                                                                                                    
H = 1961/01/21 03 45 25.00         35.2500   10.5000                   0.0                    BCIS                1901692
STOP                                                                                                                     
E   Tunisia                                                                                                              
I   ISC        853630                                                                                                    
H = 1965/09/05 22 06 55.58         34.1967    8.6501                  10.0                    ISC                00876034
M   4.3  mb    ISC                                                                                               00876034
P   ISO     10.05 353  Pn       1965  9  5 22  9 20.8   -1   0.8 Pn            .        .ISO  .  .    ISC        28078860
P   IFR     11.48 271  Pn       1965  9  5 22  9 40.6   -1   0.9 Pn            .        .IFR  .  .    ISC        28078861
P   WRM     15.81 352  P        1965  9  5 22 10 48.6   -1   6.6 P             .        .WRM  .  .    ISC        28078862
P   DOU     16.17 351  P        1965  9  5 22 10 48.3   -1   2.3 P             .        .DOU  .  .    ISC        28078863
P   EKA     22.67 342  P        1965  9  5 22 11 57.0   -1   0.1 P             .        .EKA  .  .    ISC        28078864
P   UPP     26.34  10  P        1965  9  5 22 12 30.8   -1  -0.6 P             .        .UPP  .  .    ISC        28078865
P   UME     30.50  10  P        1965  9  5 22 13  6.6   -1  -1.9 P             .        .UME  .  .    ISC        28078866
P   KIR     34.36   8  P        1965  9  5 22 13 35.3   -1  -7.0 P             .        .KIR  .  .    ISC        28078867
P   SHI     37.29  85  P        1965  9  5 22 14  9.0   -1   0.9 P             .        .SHI  .  .    ISC        28078868
P   SCH     54.77 317  P        1965  9  5 22 16 23.9   -1  -1.4 P             .        .SCH  .  .    ISC        28078869
P   CPO     74.07 302  P        1965  9  5 22 18 31.0   -1  -1.0 P             .        .CPO  .  .    ISC        28078870
P   WMO     83.41 308  P        1965  9  5 22 19 23.0   -1  -0.4 P             .        .WMO  .  .    ISC        28078871
P   UBO     86.51 318  P        1965  9  5 22 19 40.0   -1   1.0 P             .        .UBO  .  .    ISC        28078872
P   BMO     87.30 325  P        1965  9  5 22 19 41.0   -1  -1.7 P             .        .BMO  .  .    ISC        28078873
P   TFO     91.65 314  P        1965  9  5 22 20  5.0   -1   1.6 P             .        .TFO  .  .    ISC        28078874
STOP                          
```

We then use `mnf_search` to refine our search, for instance removing those without much phases read:

```console
➜  Data git:(master) ✗ ./mnf_search
Release date October 7, 2018

Enter input filename: 
test.dat.mnf
Create a new bulletin (1) or individual event files (2)?
2
Create mloc command file? 
y      
Enter command file basename: 
test1.1
Use lat-lon limits?
n
Use focal depth limit?
n
Use nearest station distance?
n
   215 events read
   205 events pass the search criteria
  80 events  that pass the search criteria and have  10 or more phase readings
  51 events  that pass the search criteria and have  20 or more phase readings
  41 events  that pass the search criteria and have  30 or more phase readings
  32 events  that pass the search criteria and have  40 or more phase readings
  30 events  that pass the search criteria and have  50 or more phase readings
  24 events  that pass the search criteria and have  60 or more phase readings
  21 events  that pass the search criteria and have  70 or more phase readings
  18 events  that pass the search criteria and have  80 or more phase readings
  15 events  that pass the search criteria and have  90 or more phase readings
  13 events  that pass the search criteria and have 100 or more phase readings
  11 events  that pass the search criteria and have 110 or more phase readings
  11 events  that pass the search criteria and have 120 or more phase readings
  11 events  that pass the search criteria and have 130 or more phase readings
  11 events  that pass the search criteria and have 140 or more phase readings
  10 events  that pass the search criteria and have 150 or more phase readings
   8 events  that pass the search criteria and have 160 or more phase readings
   7 events  that pass the search criteria and have 170 or more phase readings
   6 events  that pass the search criteria and have 180 or more phase readings
   6 events  that pass the search criteria and have 190 or more phase readings
   6 events  that pass the search criteria and have 200 or more phase readings
Minimum number of phase arrivals: 
50
Event number selection: beginning and end numbers: 
1 215
EOF reached after    215 events
  30 events selected
```

Now we will have a list of `.mnf` file each of which contains one event with at lease 50 phases identified. Those will be the events we are going to relocate.

### 3. Fulfill command file

You should now see a command file (all command file are ended with `.cfil`) named `test1.1.cfil` as you previously specified in your `mnf_search`, which looks like a bunch of repetition of:

```
memb
even 19670126.1611.42
inpu 19670126.1611.42.mnf
memb
even 19680423.2230.24
inpu 19680423.2230.24.mnf
memb
```

Here, `memb`, `even` and `input` together specify each of the input events with their unique names and input `.mnf` files, which will be used by `mloc` for obtaining phase information as well as travel time.

### 4. Run the program

First of all, let's make a new folder named `test` in our primary working directory, again which is `MLOC/mloc_working`, and copy paste all the `.mnf` file as well as the command file `.cfil` there from previous data folder. Let's take a look at `cfil header.txt` which is in the working folder:

```
pltt 1 2 4 5 6 7
dem1 globe
bdps tables/stn/bdps.dat
ppri pP
ppri sP
dcal on
phyp off
hlim 0. 1.0
clim 0. 180.
wind 3 4
frec 1 1 0 1
freh 1 1 0 1
depc 16
```

This file serves as a list of common commands that will be used almost everywhere in your analysis with **MLOC**. A brief summary is as below:
- `pltt`: Control which kind of plot will be produced
- `dem1`: Determine which global velocity model will be used, here is *AK135*. For regions in oceans, `etopo1` is more preferential.
- `bdps`: Provide a list of stations that are known to produce bogus phase information.
- `ppri`: Prevent phase re-identification. **MLOC** has a number of procedures that will re-identify phases that are not reasonable w.r.t. velocity models we use.
- `dcal`: Direction calibration. Since we have local data with numerous Pn/Pg, it is possible to do direct calibration as oppose somewhere in the middle of ocean where only tele-seismic phase are available.
- `phyp`: P phases only for hypercentroid. For some places where S phases are prone to large reading errors, it is better to be turned on.
- `hlim`: Distance limit (in degrees) in data used for hypercentroid
- `clim`: Distance limit (in degrees) in data used for cluster vector (all data are used here)
- `wind`: Determine the weight coefficient w.r.t. sigma (confidential interval)
- `frec`: Free parameters for clusters: latitude longitude depth origin_time. It is often recommended not to do free depth relocation at first before removing outliers.
- `freh`: Free parameters for hypercentroid, the same as above.
- `depc`: Assign a uniform depth to all the events otherwise depths information in `.mnf` file will be used instead.

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
