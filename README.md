# MLOC

![AUR](https://img.shields.io/badge/version-v10.4.5-brightgreen.svg)
![AUR](https://img.shields.io/badge/release-10%2F15%2F2018-orange.svg)
[![AUR](https://img.shields.io/badge/license-GPLv3-blue.svg)](https://www.gnu.org/licenses/quick-guide-gplv3.en.html)
[![AUR](https://img.shields.io/badge/Docs-latest-7f00ff.svg)](https://github.com/shipengcheng1230/MLOC/tree/master/Docs)


## Introduction

**MLOC** is a multiple-events relocation program developed and maintained by [Eric Bergman](https://www.researchgate.net/profile/Eric_Bergman2) since 1990 whose underlying idea is based on [Jordan and Sverdrup (1981)](https://pubs.geoscienceworld.org/ssa/bssa/article-abstract/71/4/1105/102070/teleseismic-location-techniques-and-their). This tutorial is a skeleton iteration based on the **MLOC** Training Workshop organized by Eric and [Ezgi](https://inside.mines.edu/~ekarasoz/) that was held in Albuquerque, NM from 10/15/18 to 10/19/18. Further update of this software will be reflected upon Eric's notification.

## Distribution

Eric has never personally distributed this software community-widely but he is fine with redistributing it by those who attended his **MLOC** workshop. Therefore, I keep the track of this software under my personal private repository and impose *GPLv3* license here. Those who shared this repository shall respect it. 

## Setup

First clone this repository to your local:
```
git clone https://github.com/shipengcheng1230/MLOC.git
```

Then download the tables file [here](https://drive.google.com/drive/folders/15Vr0Gi_0WSK73DNHGmBA49EgFkn84YZl?usp=sharing) and put it under `MLOC/mloc_working/tables`. Your directory should looks like:

```
├── Docs
│   ├── Focal Depth
│   ├── HD calibration cartoons
│   │   ├── Direct
│   │   └── Indirect
│   │       └── old ones
│   ├── Hi-Rez topographic datasets.rtfd
│   └── Papers
├── MNF utilities
│   ├── mnf search
│   └── to\ mnf
│       ├── any2mnf
│       ├── isc_ims2mnf
│       └── seisan2mnf
├── clusters
│   └── Tunisia
│       └── Test
│           └── Data
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
    │   └── jsa4.1_example_gmt_scripts
    └── tables
        ├── crust
        ├── ellipticity
        ├── faults
        │   └── iran_faults_ghods
        ├── gmt
        │   ├── cpt
        │   │   └── Examples
        │   │       └── Qeshm
        │   └── dem
        │       ├── ETOPO
        │       ├── GEBCO
        │       │   └── gebco_docs
        │       ├── GINA
        │       ├── GLOBE
        │       └── custom
        ├── kml
        ├── spread
        ├── stn
        │   └── archives
        │       ├── convert_geog
        │       └── ere_stn_geog
        └── tau-p
```

