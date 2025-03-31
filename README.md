# TS_pipeline
**March 31st, 2025**

 (Semi)automatic pipeline for transition state (TS) search in molecular photoswitches.

## Installation Instructions
### Local
```
python -m venv TS_pipeline.venv
git clone git@github.com:CrespiLab/TS_pipeline
source TS_pipeline.venv/bin/activate
sudo apt install libopenbabel-dev swig
pip install openbabel --global-option=build_ext --global-option="-I/usr/include/openbabel3" --global-option="-L/usr/lib/x86_64-linux-gnu/openbabel/3.1.1/"
cd TS_pipeline
pip install .
```
Adjust according to your folder structure:
--global-option="-I/usr/include/openbabel3" and --global-option="-L/usr/lib/x86_64-linux-gnu/openbabel/3.1.1/"

### Programmes required
These programmes need to be callable f***REMOVED*** the command-line.
If not: adjust config.py and add paths to programmes.
- OpenBabel
If you need to compile it, see here:

- xTB

-  CREST

-  ORCA



### Clu***REMOVED***r




## TEST
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


