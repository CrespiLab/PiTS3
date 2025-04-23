# PiTS<sup>3</sup>
**March 31st, 2025**

 (Semi)automatic pipeline for transition state (TS) search in molecular photoswitches.

## Installation Instructions
### Local
```
python -m venv PiTS3.venv
git clone https://github.com/CrespiLab/PiTS3/
source PiTS3.venv/bin/activate
sudo apt install libopenbabel-dev swig
pip install openbabel --global-option=build_ext --global-option="-I/usr/include/openbabel3" \
                      --global-option="-L/usr/lib/x86_64-linux-gnu/openbabel/3.1.1/"
cd PiTS3
pip install -e .
```
Adjust according to your folder structure:
--global-option="-I/usr/include/openbabel3" and --global-option="-L/usr/lib/x86_64-linux-gnu/openbabel/3.1.1/"

### Global environment
Note that xTB and CREST depend on **OMP_NUM_THREADS** global variable, which controls the number of cores available to those utilities. Adjust it accordingly.

**By default**, the tool will check the accessibility of **crest**, **xtb**, **orca** and **pysis** in your environment (in PATH). If you have multiple installations/versions of these utilities, make sure the right one is in your PATH or change ts_pipe/config.py accordingly (see instructions in the file). It is especially relevant for working in a supercomputer cluster environment, where you might have a module system for utility loading.

### Programmes required
These programmes need to be callable from the command-line (i.e. should be available in your PATH).
If not: adjust config.py and add paths to programmes.
- OpenBabel
- xTB
-  CREST
-  ORCA

### Cluster

## TEST
We provide a set of test examples. These can be run in the test folder using the provided bash script. 
```
(PiTS3.venv) /PiTS3$ cd test
(PiTS3.venv) /PiTS3/test$ ./test.sh
```

# Instructions

## Tools
- tspipe
- ts_data_collector
- fragments_combiner
- cdxml_to_xyz

### Fragments Combiner Tool (fragments_combiner)

## TS Mode Selection
The tool

## Run
###

### Local (needs to be updated)
main.py *.xyz | tee output.log

## TEST
```
../../../src/fragments_combiner.py nbd_*

../../../main.py C12H9NS.xyz
../../../main.py C12H9NS.xyz -m "TSModeFragment"
../../../main.py C12H9NS.xyz -m 'C1C=CCC=C1'

```




 
