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

```shell
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

`Docs` contains all the related references and manuals. `MNF utilities` contains functionalities for transforming ISC bulletin to `.mnf` format. `mloc utilities` contains functionalities for **MLOC** processing procedures. `mloc_gfortran` and `mloc_intel` contain `makefile` for building `mloc` program from source code listed under `mloc_src`. `mloc_working` will be your primary working directory, in which the `table` has all the geological models used by this program.

Compile source code and any other utilities, for example using `gfortran`:

```
cd mloc_gfortran && make && cp ./mloc_g ../mloc_working
```

Modify `mloc.conf` file such that `WORKING_DIR` points to your `mloc_working` directory and `AUTHOR` to your name not exceeding 8 characters (for more configuration see docs):

```
WORKING_DIR: /Users/spc/Softwares/mloc/mloc_working
AUTHOR: SPC
```

Now you shouldn't have problem running the program:

```shell
➜  mloc_working git:(master) ✗ ./mloc_g 

mloc v10.4.5, release date 9/28/2018                                            

Current program limits: 
  nevmax =   200
  nqmax  =  4000
  ntmax1 = 35000

Enter a basename for this run: 
```