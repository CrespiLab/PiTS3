# To-Do List

### To Do
- [ ] Merge everything
- [ ] Improve file and directory handling in Python
  - [ ] Switch from os.system to subprocess.run
  - [ ] Better file names
- [ ] Implement restart option
  - [ ] Write checkpoint (.json) file during run
  - [ ] "Passive" restart for cancelled jobs
  - [ ] Restart from a certain step with different conditions - is it
  possible without input files?
  - [ ] Hardcode step numbers
  - [ ] ORCA stage statistics: if reached; paths to start ORCA jobs; check
  ORCA termination
- [ ] Clean in ORCA directory (keep .err, .inp, .out, .xyz, .hess, .engrad,
.log, .gbw)
  - [ ] Rethink scratch usage when input files are implemented
- [ ] Scrape data immediately after finishing of ORCA compound job: put relevant data into one output file (e.g., energies)
  - [ ] Transfer part of the data_collector code to core to get ΔE and
  barriers when the user follows the standard pipe
- [ ] Check ORCA in the config, but do not stop if it is not found
- [ ] Scrape generated .xyz files and output into .json results (actual
geometries)

- [ ] Switch over to using ORCA for all steps
  - [ ] Investigate possibility of using NEB
- [ ] Implement three-point TS search to allow search for higher-energy mechanisms
- [ ] Add option to perform conformational analysis before TS search (or
triple-crest)
- [ ] Switch over to using input files for all steps

- [ ] Investigate possibility to first find TS and then apply substitutions 

### In Progress :)
- [x] Include data collector tool

### Completed ✓
- [x] Clean up pysisyphus/qm_calcs folder after completion
- [x] Fix fragment combiner (RDKit)

