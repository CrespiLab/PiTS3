# PiTS<sup>3</sup>
**April 23<sup>rd</sup>, 2025**

 (Semi)automatic pipeline for transition state (TS) search in molecular photoswitches.

## Installation Instructions
### Local
```
python -m venv PiTS3.venv
git clone https://github.com/CrespiLab/PiTS3/
source PiTS3.venv/bin/activate
(PiTS3.venv) sudo apt install libopenbabel-dev swig
(PiTS3.venv) pip install openbabel --global-option=build_ext --global-option="-I/usr/include/openbabel3" \
                      --global-option="-L/usr/lib/x86_64-linux-gnu/openbabel/3.1.1/"
(PiTS3.venv) cd PiTS3
(PiTS3.venv) pip install -e .
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

## Test run
We provide a set of test examples. These can be run in the test folder using the provided bash script. 
```
(PiTS3.venv) /PiTS3/test/test1_cdxml-to-xyz$ ./test1_cdxml-to-xyz.sh
(PiTS3.venv) /PiTS3/test/test2_fragments-combiner_tspipe_noorca$ ./test2_fragments-combiner_tspipe_noorca.sh
(PiTS3.venv) /PiTS3/test/test3_fragments-combiner_tspipe_orca$ ./test3_fragments-combiner_tspipe_orca.sh
```

# Instructions
## Full run (locally)
Example for stilbene, starting from `fragments_1.cdxml` (1 fragment) and `fragments_2.cdxml` (2 fragments).
The fragment combiner tool yields Cartesian coordinate files of two molecules, named after their molecular formulas: 
`C15H11N.xyz` and `C15H14O.xyz`.
```
fragments_combiner stilbene_fragments_1.cdxml stilbene_fragments_2.cdxml
tspipe -m 'C-C=C-C' C15H11N.xyz C15H14O.xyz
ts_data_collector -m 'C-C=C-C' C15H11N C15H14O
```

## List of Tools
- fragments_combiner
- tspipe
- ts_data_collector
- cdxml_to_xyz

### Fragments Combiner Tool (fragments_combiner)
```
fragments_combiner frag1.cdxml frag2.cdxml
```

### tspipe
```
tspipe -m 'C-C=C-C'    stilbene.xyz       | tee output_stilbene.log
tspipe -m 'C-C=N-C'    diarylimine.xyz    | tee output_diarylimine.log
tspipe -m 'C1=CCNC=C1' ADAE.xyz           | tee output_ADAE.log
tspipe -m 'C1C=CCC=C1' NBD.xyz            | tee output_NBD.log
```

### ts_data_collector
```
ts_data_collector -m 'C-C=C-C'    stilbene
ts_data_collector -m 'C-C=N-C'    diarylimine 
ts_data_collector -m 'C1=CCNC=C1' ADAE
ts_data_collector -m 'C1C=CCC=C1' NBD
```
