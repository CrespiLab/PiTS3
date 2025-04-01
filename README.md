# TS_pipeline
**March 31st, 2025**

 (Semi)automatic pipeline for transition state (TS) search in molecular photoswitches.

## Installation Instructions
### Local
```
python -m venv TS_pipeline.venv
git clone https://github.com/CrespiLab/TS_pipeline/
source TS_pipeline.venv/bin/activate
sudo apt install libopenbabel-dev swig
pip install openbabel --global-option=build_ext --global-option="-I/usr/include/openbabel3" \
                      --global-option="-L/usr/lib/x86_64-linux-gnu/openbabel/3.1.1/"
cd TS_pipeline
pip install -e .
```
Adjust according to your folder structure:
--global-option="-I/usr/include/openbabel3" and --global-option="-L/usr/lib/x86_64-linux-gnu/openbabel/3.1.1/"

### Global environment
Note that xTB and CREST depend on OMP_NUM_THREADS global variable, which controls the number of cores available to those utilities. Adjust it accordingly.

By default, the tool will check the accessibility of crest, xtb, orca and pysis in your environment (in PATH). If you have multiple installations/versions of these utilities, make sure the right one is in your PATH or change ts_pipe/config.py accordingly (see instructions in the file).

### Programmes required
These programmes need to be callable f***REMOVED*** the command-line (i.e. should be available in your PATH).
If not: adjust config.py and add paths to programmes.
- OpenBabel
- xTB
-  CREST
-  ORCA



### Clu***REMOVED***r



# *TBD*

##TEST
```
../../../src/fragments_combiner.py nbd_*

../../../main.py C12H9NS.xyz
../../../main.py C12H9NS.xyz -m "TSModeFragment"
../../../main.py C12H9NS.xyz -m 'C1C=CCC=C1'

```

## Instructions
### Fragments Combiner Tool

### TS Mode Selection
The tool 

## RUN
### Local
main.py *.xyz | tee output.log


