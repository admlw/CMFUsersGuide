# CMF Users' Guide
## Thomas Carroll, Harry Hausner, Adam Lister, Brian Rebel and Jenny Thomas


Table of Contents
=================

   * [Introduction](#introduction)
   * [Event Lists](#event-lists)
         * [cmf::EventContainer](#cmfeventcontainer)
         * [cmf::MetaData](#cmfmetadata)
         * [cmf::SpillSummary](#cmfspillsummary)
   * [CMF Overview](#cmf-overview)
      * [Structure](#structure)
      * [core](#core)
         * [Event](#event)
         * [ShifterAndWeighter](#shifterandweighter)
         * [VarVals](#varvals)
         * [Spectrum Tools](#spectrum-tools)
      * [data](#data)
      * [dataProducts](#dataproducts)
      * [fhicl](#fhicl)
      * [macros](#macros)
      * [modules](#modules)
      * [scripts](#scripts)
      * [utilities](#utilities)
   * [Instructions](#instructions)
      * [Generating Event Lists](#generating-event-lists)
         * [Producing EventListTree Files](#producing-eventlisttree-files)
            * [art Module: CMFCAFToEventLists](#art-module-cmfcaftoeventlists)
            * [script: makeProductionXEventLists.sh](#script-makeproductionxeventlistssh)
            * [script: startEventLists.sh](#script-starteventlistssh)
         * [Producing CAF Text Lists](#producing-caf-text-lists)
            * [cafe macro: get_eventlist.C](#cafe-macro-get_eventlistc)
            * [macro: compareEvents.C](#macro-compareeventsc)
            * [script: compareEventLists.sh](#script-compareeventlistssh)
         * [Looking for Duplicate Events Among Selections](#looking-for-duplicate-events-among-selections)
            * [macro: compareSelectedSets.C](#macro-compareselectedsetsc)
            * [script: compareSelections.sh](#script-compareselectionssh)
      * [Generating Covariance Matrices](#generating-covariance-matrices)
         * [Locally](#locally)
         * [On The Grid](#on-the-grid)
      * [Random Universes](#random-universes)
         * [Random Universe Generation](#random-universe-generation)
      * [Sensitivity Contour Production](#sensitivity-contour-production)

# Introduction



This Users' guide is intended to be a partner document to the CMF technote, also contained within this docdb entry. It documents the the Covariance Matrix Fit (CMF)  source code as it stands for the production 5 analysis, and provides an example for how the analysis is run start-to-end. 

Development of the CMF code base is primarily done in the nux development branch (`nux-dev-br`).

CMF is developed using the CMake build of NOvASoft, instructions for which can be found [here](https://cdcvs.fnal.gov/redmine/projects/novaart/wiki/Editing_Code_with_CMake_and_buildtool). Help is nearby on slack at #cmakebuild and #cmf.



# Event Lists

`cmf::EventList`s are the event record container used by the CMF framework. 

Inside an EventList ROOT file, there is a set of trees, the first tree is the `metadata` tree, which contains two branches, `metadata` and `spillsummary`.

The `metadata` tree contains information related to the `cmf::MetaData` class (detector, filetype, selection, interaction type, and period of the selected events), while the `spillsummary` tree contains information related to the `cmf::SpillSummary`  class (POT, livetime, and number of spills for each selection).

In addition, there are a set of trees which correspond to a `cmf::EventContainer` for each selection type, for example `MCNearEpoch5FHCBeamNuMuSelQ4NuEBarCC`, which holds information related to the `cmf::EventList` for the selection for the Near Detector/Epoch 5/FHC/NuMuSel Quartile 4/ NuEBarrCC interactions. These trees contain a branch for the event ID (run, subrun, event information), a branch for reconstructed information (`dataVars`), a branch for truth information (`truthVars`),  a branch for different event weights (`weightVars`), and a branch for each of the GENIE weights. Each GENIE weight is stored separately due to the need to store "-2", "-1", "+1" and "+2" sigma values. 

This section goes into more detail about the `cmf::EventContainer`, `cmf::MetaData` and `cmf::SpillSummary` classes.

### cmf::EventContainer

The event container contains a vector of `cmf::Event` objects, which each contain:
```
EventID     : the event ID
DataVarVals : the reconstructed variables for the event
MCVarVals   : the truth-level variables for the event 
              (truth, wieght, and genie variables)
```

**N.B.** The names of `DataVarVals` and `MCVarVals` are somewhat misnomers. `DataVarVals`is more akin to reconstructed variables, while the `MCVarVals` is more akin to truth information, though it also includes information related to the weights.

### cmf::MetaData

The `cmf::MetaData` class holds information which fully categorises 
```
novadaq::cnv::DetID     : detector ID (only ND or FD are really used for CMF):
                          kUNKNOWN_DET
                          kNEARDET 
                          kFARDET     
                          kNDOS      
                          kNDSBTEST
                          kTESTBEAM
                          kNDetector
                          kFCCDAQ

cmf::Filetype_t         : type of file:
                          kBeam
                          kSwap
                          kTauSwap
                          kDataFile
                          kCosmicBackgroundFile
                          kRockFluxSwap
                          kRockNonSwap
                          kNumFileTypes
                          kBadFileType
                          kUnknownFileType
                          
cmf::SelectionType_t    : choice of selection:
                          kNuMuSelection
                          kNuESelection
                          kNuESelectionLowPID
                          kNuESelectionMidPID
                          kNuESelectionHighPID
                          kNuESelectionPeripheral
                          kNuMuSelectionQ1
                          kNuMuSelectionQ2
                          kNuMuSelectionQ3
                          kNuMuSelectionQ4
                          kNCSelection
                          kUnknownSelection
                          
cmf::InteractionType_t  : type of interaction (from Genie):
                          kUnknownInteraction
                          kNuMuCC,
                          kNuMuBarCC,
                          kNuECC,
                          kNuEBarCC,
                          kNuTauCC,
                          kNuTauBarCC,
                          kNC,
                          kCosmicMuon,
                          kRockMuon


std::string (epoch)     : the epoch of the data
```

### cmf::SpillSummary
The spill summary holds information about the beam spill:
```
totalPOT        : total POT for this data
goodPOT         : total good POT (passes good run selection) for this data
liveTime        : total live time recorded
totalNumSpills  : total number of spills
numGoodSpills   : total number of good spills (passes good run selection)
```



# CMF Overview

## Structure

The CMF code lives within `novasoft/CovarianceMatrixFit`. Within this directory there are several sub-directories:

```
core        : the core classes that CMF is built on: EventLists, VarVals, ShifterAndWeighter
data        : contains root files for calibration systematic uncertainties
dataProducts: data products and structs commonly used in CMF analysis
fhicl       : contains all fhicl files
macros 	    : useful .C files
modules     : art modules and plugins which are run with fhicl files
scripts     : scripts for, e.g. submitting to the grid
utilities   : helper utility classes and methods
```



## core

### Event

The event record container, `cmf::Event` is defined in `Core/Event.h` and `Core/Event.cxx`. As previously discussed, the `cmf::Event` contains information related to the event ID as well as the dataVarVals and mcVarVals.

The Event defines some useful typedefs:

- `typedef std::unique_ptr<Event> EventPtr;`
- `typedef std::vector<EventPtr>  EventContainer;`
- `typedef std::map<cmf::MetaData, cmf::EventList> EventListMap;`

The Event .cxx file also contains two methods which write EventLists into ROOT files, 

- `SerializeEventListMap`: this helps to create the relevant trees for each selection type
- `SerializeEvents`: helper function to setup the different branches for the above trees

### ShifterAndWeighter

The ShifterAndWeighter does the work of evaluating how the relative weight of each simulated event changes based on the various systematic uncertainty parameters or the oscillation weights and determining the simulated event energy if a shift in either the hadronic or leptonic portion of the energy is desired.

### VarVals

VarVals defines how event information is stored in the ROOT files and retrieved at analysis time.  There are several structs and classes used, 

- `DataVars`		: a struct which contains a small set reconstructed quantities for each event. It also holds a "fake weight" which holds a weight for the events based on exposure or in the case of fake data the total weight we want to assign to the event;
- `DataVarVals`  : a class containing a DataVars and a method to determine the reconstructed neutrino energy based on what the ShifterAndWeighter says should be shifted;
- `EventID`          : a class that uniquely identifies the event based on run, subrun, event, slice and in the case of simulation, MC cycle number;
- `TruthVars`  : a struct containing relevant truth information such as the neutrino flavor, energy, and parent particle;
- `WeightVars`: a struct containing weighting information that is not contained within the portion of MCVars 
- `MCVarVals`  : a class containing TruthVars, WeightVars, and a vector of floats corresponding to  systematic uncertainties

### Spectrum Tools

The SpectrumTools class has one purpose: it takes in an event list and a given point in parameter space and fills a spectrum, applying oscillation and systematic weights based on the point in parameter space.

## data

This directory houses the files which store the CMF file systematic uncertainties, for example, calibration systematic uncertainties

## dataProducts

This directory contains all of the cmf classes, structs and other definitions used by CMF modules. 

- `Constants.h` contains constant variables which are used through the rest of CMF.
- `StaticFuncs.h` contains helper functions.
- `Structs.h` contains structs, primarily the cmf::MetaData and cmf::SpillSummary.
- `InputPoint.cxx/h` contains the InputPoint class, which defines a given point in parameter space including both oscillation parameters and systematic parameters.
- `PointResult.cxx/h` contains the PointResult class, which contains information related to the fit result at a given point in parameter space, for example the &chi;<sup>2</sup> value, and predicted spectrum. This is used in the case where a fit is performed using the CovarianceMatrixFitter plugin.
- `GridPointResult.cxx/h` contains the LibraryPoint and PredictionPoint classes, as well as classes they depend on (GridPoint, HiddenParameters, HiddenParamPoint, LPUniverseFitResult). These are used in the case where sensitivities are produced using a prediction library.
- `FakeUniverse.cxx/h` contains the FakeUniverse class, which wraps up an InputPoint with the asimov and poisson-fluctuated fake data.



## fhicl

The fhicl directory contains fhicl files for the CMF analysis. The format of the fhicl files is as follows:

- Lower case fhicl files with the format  `cmf_*job.fcl` are fhicl files which you may run.
- Upper case fhicl files with the format `CMF_*.fcl` are fhicl files which contain configurations.

In general, there is a configuration-level fhicl file for each job-level fhicl file.



## macros

The macros directory contains handy macros for producing plots. Note that some of the macros must be compiled with the rest of novasoft as they hook into some of the CMF code. Each of the macros which must be compiled have a corresponding `art_make_exec` listing in the CMakeLists.txt.

## modules

This directory contains the CMF modules and plugins which perform the bulk of the CMF analysis. A non-exhaustive list of important modules/plugins is as follows:

- `CovarianceMatrixMaker_plugin.cc` handles production of the covariance matrices.
- `CMFRandomUniverses_plugin.cc` handles generation of poisson-fluctuated and systematically shifted random universes.
- `CMFPredictionLibrary_plugin.cc` generates libraries of predictions at a given point in parameter space with a configurable number of predictions for each hidden parameter.  This is generally used for sensitivities.
- `CovarianceMatrixFitter_plugin.cc` calls MINUIT to perform a fit for a given set of oscillation parameters
- `CMFLibraryContourMaker_plugin.cc` uses the output of the prediction library or the fitter in order to construct a &Delta;&chi;<sup>2</sup> space, and contours.
- `EventListManipulator.cxx/h` is used for interpreting the Event List data format 

## scripts

The scripts directory holds a set of bash scripts to make job running and submission easier. There is a README which provides a description of what each of the scripts does.



## utilities

This directory holds a set of utility classes which are generally accessed by multiple CMF modules. 



# Instructions

## Generating Event Lists

This section details how to generate a CMF-style tree of selected events using the CAF selections and how to compare the generated tree to the CAF selection.

### Producing EventListTree Files

The CovarianceMatrixFit package provides an art module to read the decaf or concatenated CAF files from a SAM dataset and produce output files containing the necessary TTrees to use as event lists in CMF.  It also has a set of scripts that can be used to run that process.

#### art Module: CMFCAFToEventLists
The module is located at `CovarianceMatrixFit/modules/CMFCAFToEventLists_module.cc`.  The module is an EDAnalyzer, but does not actually do anything in the analyze method because the CAFs are not art files.  We use a module here to be able to configure the conversion job and because art provides a nice state machine that processes jobs in a particular order.  All the work of the module is done from methods called by the module's `endJob` method.  The `CMFCAFToEventLists::FillVariables()` method calls other methods to fill the `cmf::DataVarVals`, `cmf::EventId` and `cmf::MCVarVals` objects.  It also calls methods to determine if the event passes the event selection determined through the configuration.  The code is the best resource for understanding how these steps are accomplished and details will not be given here.  Once the events have been selected they are written into ROOT TTrees based on the metadata information for each event.  

The module must undergo rewrites for each new production due to changing selections and changing lists of systematic uncertainties. 

#### script: makeProductionXEventLists.sh

The script that runs the above module are located at
```
CovarianceMatrixFit/scripts/makeProuction{4,5}EventLists.sh
```
It has six functions; one function prints the necessary command line inputs to the terminal when the script is called without inputs.  The output from that function is
```
 Usage                : makeEventLists.sh 
                        <outdir> <detector> <file type> 
                        <selection> <horn current> <files per iteration> 
                        <tag> <systematic name>
 <outdir>             : directory where the output root files will go
 <detector>           : Near or Far
 <file type>          : data, nonswap, fluxswap, tauswap, 
                        cosmicbackground, rockfluxswap, rocknonswap
 <selection>          : NuESel, NCSel, or NuMuSel
 <horn current>       : fhc or rhc
 <files per iteration>: number of files for each subjob
 <tag>                : analysis tag, eg 2018, 2019
 <systematic name>    : optional name for systematically shifted set
```
The last argument is optional and is only used when making event lists for systematically shifted datasets.

The remaining functions handle converting files from SAM datasets containing the concatenated CAF files. The most important function for converting the concatenated CAF files is `makeXROOTDFile` which creates a text file containing the xrootd location of the concatenated CAF files.  This method checks which combination of detector, file type, selection, and horn current is being requested and then creates the corresponding SAM dataset name from those inputs. 

**N.B** the SAM dataset name is not dynamically determined based on the production run desired.  The base of the name must be changed by hand for each production run.

The script configures the art job to produce a text file containing the event identification information, ie run, subrun, event, slice, and cycle numbers, the event energy, and the CVN value for the event.  This text file has a name of the form `selectedEvents_<detector>_<selection>_<horn current>_<file type>_<systematic>_<analysis tag>.txt`. 

#### script: startEventLists.sh

The script that calls the `makeEventListsFromDecaf.sh` script is located at `CovarianceMatrixFit/scripts/startEventLists.sh`.  If the script is called without arguments it prints the necessary command line inputs to the terminal,
```
  Usage                 : startEventLists.sh 
                          <production> <detector> <file type> <files per iteration> 
                          <tag> <systematic_name>
  <production>          : Production number of the analysis, i.e. 4,5. 
  						  This configures which makeProductionXEventLists.sh script
  						  to run.
  <detector>            : Near or Far
  <file type>           : data, nonswap, fluxswap, 
                          tauswap, cosmicbackground, 
                          rockfluxswap, rocknonswap
  <files per iteration> : number of files for each subjob
  <tag>                 : analysis tag, eg 2018, 2019
  <outdir>              : location where the output root files should go
  <systematic name>     : optional name for systematically shifted set
```
The script will start a series of jobs that run locally and in the background (using `nohup`) to convert each of the possible event selections and horn currents available in the `makeEventListsFromDecaf.sh` script.

This script must be run multiple times to create a complete eventlist. For example, for the `calib-shift-(fd/nd)-xyview-neg-offset` systematic uncertainty, the script must be run in the following configurations:

- near
  - nonswap
- far
  - nonswap
  - fluxswap
  - tauswap

This script runs the `CMFCAFToEventLists` module locally. Check each dataset is complete with `htop` to avoid overwhelming the GPVMs. As a general guide, 10 files per iteration works reasonably well.

### Producing CAF Text Lists

#### cafe macro: get_eventlist.C

The `CovarianceMatrixFit/macros/cafe/get_eventlist.C` macro creates a text file list of events from a given SAM CAF concatenated dataset.  The macro must be run with the `cafe` executable.  The macro has to be provided with the desired SAM dataset name, the name of the output text file, the selection to use and the detector. Those values are given nonsense default values which can be replaced easily using a `sed` command as is done in the `CovarianceMatrixFit/scripts/compareEventLists.sh` script discussed next.  If one wants to change the values at run time they can be given as arguments to the macro at run time as well. The output file name has the form `cafEvents_<detector>_<selection>_<horn current>_<file type>_<systematic>_<analysis tag>.txt}`.

This macro must be run from a terminal that is set up to run novasoft in an SRT environment.  One cannot run the SRT and ups/CMake environments at the same time in the same terminal.

#### macro: compareEvents.C

The `CovarianceMattrixFit/macros/compareEvents.C` macro compares the output text files of selected events produced by the `CovarianceMatrixFit/macros/cafe/get_eventlist.C` macro and the `CovarianceMatrixFit/modules/CMFCAFToEventLists_module.cc` module. The file names for the CAF and CMF text files are assumed to follow the forms indicated above.

This macro makes histograms comparing the selected number of events by CAF and CMF, including any events that are in one sample but not the other.  It also produces histograms showing the differences between the reconstructed energy, product of PPFX and cross section central value weights, and the CVN score.  The macro produces comparisons of the weighted and unweighted spectra as well.

#### script: compareEventLists.sh

The `CovarianceMatrixFit/scripts/compareEventLists.sh` script will run the `CovarianceMatrixFit/macros/cafe/get_eventlist.C` macro and compare the output from that macro to the output text file from the `CovarianceMatrixFit/modules/CMFCAFToEventLists_module.cc`.  This script also runs the `CovarianceMattrixFit/macros/compareEvents.C` macro described above. If the script is called without arguments it prints the necessary command line arguments to the screen,
```
  Usage             : compareEventLists.sh 
                      <selection> <filetype> <tag> <systematic name>
  <selection>       : nue, numu, nus
  <filetype>        : data, nonswap, fluxswap, tauswap, cosmicbackground, rockfluxswap, rocknonswap
  <tag>             : analysis tag, eg 2018, 2019
  <systematic name> : optional name for systematically shifted set
```
Both the SRT build of novasoft and samweb must be set up or the script will exit with a corresponding message. The script will compare all combinations of horn current and detector for the given set of selection and file  type chosen at run time.

### Looking for Duplicate Events Among Selections

#### macro: compareSelectedSets.C

The `CovarianceMatrixFit/macros/compareSelectSets.C` macro compares the event lists in the text files produced by the `CovarianceMatrixFit/modules/CMFCAFToEventListsd_module.cc` module to see if any events in one PID selection, ie $/nu_{/mu}$, $/nu_{e}$ or NC, overlaps with any other selection.  It requires the detector, the selections to be compared, the file type, horn current and analysis tag be specified.  The macro produces a text file containing the duplicated event identification information.  The produced file name has the form `duplicate_list_<detector>_<selection 1>_<selection 2>_<horn current>_<file type>_<analysis tag>.txt`.  If one wanted to compare the CAF selected events, the macro would have to be edited such that the input lists to compare have file names starting with `cafEvents` rather than `selectedEvents`.

#### script: compareSelections.sh

The `CovarianceMatrixFit/scripts/compareSelections.sh` script calls the `CovarianceMatrixFit/macros/compareSelectSets.C` macro for all combinations of detector, horn current and file type for the analysis tag and set of selections to compare. If the script is called without arguments it prints the necessary command line arguments to the screen,
```
  Usage         : compareSelections.sh <tag> <selection 1> <selection 2>
  <tag>         : analysis tag, eg 2018, 2019
  <selection X> : nue, numu, nus
``````


## Generating Covariance Matrices

This section contains information on how to generate covariance matrices for a given systematic uncertainty both locally.

### Locally

In order to generate a covariance matrix locally, a single fhicl file can be run 

```
cmf_covariancematrixmakerjob.fcl
```

**N.B.** Ensure you're using the CMF version of this file. Another version, `covariancematrixmakerjob.fcl` exists, but is related to the FNEX framework and is deprecated.

Inside this fhicl file, there are four options that a user should configure:

```
TREEFILE  : Path to EventList file. 
            For now these are located in /nova/ana/users/brebel/skimmed
SYSTPAR   : Systematic parameter to vary
NUMITER   : Number of iterations of the systematic to run 
            (i.e. number of universes)
DBS_SET   : The selections to use for the universes (e.g. dbs_nd_fd_fhc_numuconcat_nc)
```

The options for which systematics you can choose can be found in `CMF_SystematicParameters.fcl`. Note that by default the EventList is assumed to be uncapped, i.e. the `TreeDirectories` label is set to `list/full`.

Once these substitutions have been made, a covariance matrix can be generated with the following command

```
art -c cmf_covariancematrixmakerjob.fcl
```

### On The Grid

Running on the grid is made easy by the existence of a bash script, located in 

```
CovarianceMatrixFit/scripts/CMF_Run_Covariance_Grid.sh
```

which has a usage:

```
usage            : CMF_Run_Covariance_Grid.sh <version> <qualifiers> <systematics> <eventlist file> <local products> <output dir>"
<version>        : specify version of novasoft to use"
<qualifiers>     : specify qualifiers for the version of novasoft"
<systematics>    : specify one of calib, genie, mec, nue, norm, reco, xsec1, xsec2, xsec3"
<selections>     : specify selections to use, e.g. dbs_nd_fd_fhc_numuconcat_nc"
<eventlist file> : full path to event list tree file in pnfs"
<local products> : name of the local products directory"
<output dir>     : top level directory for output matrix root files"
```

where the `<local products>` is the name of a tar file (without the `.tar.bz2`) which should be stored in your resilient area (`/pnfs/nova/resilient/users/${USER}`).

This script by default will run a set of several uncertainties defined by the `<systematics>` tag. Each systematic uncertainty will have by default 2000 systematic universes across 200 grid jobs. This can be modified by changing the variables in the `setVariables()` function in the bash script.



## Random Universes

We want to be able to test how an analysis responds to the conditions of having different values for systematic uncertainties as well as statistical fluctuations.  Each combination of systematic uncertainties and oscillation parameters is referred to as a "random universe".  

### Random Universe Generation

The `CovarianceMatrixFit/modules/CMFRandomUniverses_plugin.cc` handles the creation of random universes.  It can be configured to either produce universes with a fixed set of input oscillation parameters, as one would want when testing Asimov and median sensitivities, or universes with any of the oscillation parameters varied as one would want when doing a Feldman-Cousins style analysis.  In either case, the systematic parameters can be varied within their allowed ranges for each simulated universe.  A single set of values is used for each universe.  When generating combinations of oscillation parameters, typically a few parameters would be held fixed while others would vary.  For example, you might fix  &Delta;m<sup>2</sup><sub>32</sub> and &theta;<sub>23</sub> but vary &delta;<sub>CP</sub> and &theta;<sub>13</sub> from universe to universe.  

The plugin can be run interactively to produce a small (< 100) number of random universes or on the grid to produce a small number of universes for each instance of the grid submission. In either case, the `CovarianceMatrixFit/fhicl/cmf_randomuniversesjob.fcl` is used to control the job.  That configuration file cannot be run as it is in the repository as it has some placeholder values that are made to be easily changed by sed when submitting to the grid. The individual configurations for either Feldman-Cousins style or Asimov/Median style running are located in `CovarianceMatrixFit/fhicl/CMF_RandomUniverses.fcl`.

Note that generating random universes with systematic fluctuations takes significantly longer than a Poisson fluctuated universe (approx. 24 hours locally for a single systematicly fluctuated universe). Bare this in mind when submitting grid jobs.

## Sensitivity Contour Production

There are several steps to producing sensitivity contours from the random universes and prediction libraries.

The first step is to run the `cmf_gridpointsbestfitjob.fcl` over the random universe file and prediction libraries as `art -c cmf_gridpointsbestfitjob.fcl <Random Universe File> <Prediction Library>`. This is best done on the grid with the `CMF_Submit_GridPoints.sh` script. This script configures the fhicl configuration to replace `MATRIXFILE` with the name of the matrix file to use, `PAR1PAR2` with the name of the parameter configuration to use, `STARTUNI` with the universe to start with, `FIRSTRUN` with a run number which will not overlap with another job's output, `NUMUNI` with the number of universes per job to run, and `<X/Y><MIN/MAX>VAL` with the range of the parameters to use. Additionally, the job can be configured to sample only a fraction of the prediction library. It is recommended to set the sampling fraction to 0.10 for internal studies as it can significantly speed up the job. For the proper analysis however the sampling fraction should be set to 1.00. The script has the usage:
```
Usage: CMF_Submit_GridPoints.sh <grid type> <novasoft version> <qualifiers> <analysis tag> <local products directory name> <output directory> <random universe file> <covariance matrix file> <parameters> <min for param 1> <max for param 1> <min for param 2> <max for param 2> <number of universes> <universes per job> <scanning fraction> <scanning seed> <scaled ND expsoure?> [OPTIONAL] <filling in missing?>
<grid type>                     : FNAL or UW
<novasoft version>              : the version of novasoft to use
<qualifiers>                    : the qualifiers for the novasoft version to use, ie eXX:grid:prof
<analysis tag>                  : tag for the analysis run, eg 2018
<local products directory name> : the base name of the directory containing your build of the novasoft ups product
<output directory>              : the top level directory where the job output should go, ie /nova/app/<your username> 
<random universe file>          : full path to the random universe file in pnfs resilient area
<covariance matrix file>        : full path to the covariance matrix file in pnfs resilient area
<parameters>                    : which parameters form your grid, e.g. Th24Dmsq41 or Th24Dmsq41_statonly
<min for param 1>               : minimum value of the first parameter to be use
<max for param 1>               : maximum value of the first parameter to be use
<min for param 2>               : minimum value of the second parameter to be use
<max for param 2>               : maximum value of the second parameter to be use
<number of universes>           : how many universes to look at total
<universes per job>             : how many universes to look at per job
<sampling fraction>             : what fraction of the prediction libraries to sample
<sampling seed>                 : what seed to use for the prediction library sample
<scaled ND expsoure?>           : Scale the ND exposure? true or false
The following arguments are required for filling in missing outputs
<filling in missing?>           : just say true to redo missing jobs
```
The jobs will return outputs containing the raw &chi;<sup>2</sup> values for a section of XY-grid space for some of the random universes. These outputs are then supplied as the imput for the `cmf_contourfromgridjob.fcl` job, which is to be run locally. It is recommended that the locations for the grid point outputs are listed in a file (using the xrootd path, if nesicarry) which can then be fed into the job as `art -c cmf_contourfromgridjob.fcl -S <source list>`. The fhicl need only be configured to correct grid parameters by changing `PAR1PAR2` to the correct configuration from `CMF_ContourFromGrid.fcl`. As this job takes several hours it is recommended to run the job under `nohup` or your perferred method for ensure your jobs are not killed should you be disconnected from the GPVM.
  
This process should give you a root file with the sensitivity contours for your random universes.
