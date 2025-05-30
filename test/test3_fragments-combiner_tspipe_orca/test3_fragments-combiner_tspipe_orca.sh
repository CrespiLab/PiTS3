#!/usr/bin/bash

# 1) fragments combiner
# 2) tspipe with orca
# 3) ts_data_collector: collect data from ORCA runs

fragments_combiner frag1.cdxml frag2.cdxml
tspipe -m 'C-C=C-C'    C14H12.xyz        --yes -o orca_three_points_BS_SP_last.inp | tee test_C14H12.log
tspipe -m 'C-C=N-C'    C13H11N.xyz       --yes -o orca_three_points_BS_SP_last.inp | tee test_C13H11N.log
tspipe -m 'C1=CCNC=C1' C19H12F6N2S2.xyz  --yes -o orca_three_points_BS_SP_last.inp | tee test_C19H12F6N2S2.log
tspipe -m 'C1C=CCC=C1' C14H11N.xyz       --yes -o orca_three_points_BS_SP_last.inp | tee test_C14H11N.log

ts_data_collector -m 'C-C=C-C'    C14H12
ts_data_collector -m 'C-C=N-C'    C13H11N      
ts_data_collector -m 'C1=CCNC=C1' C19H12F6N2S2
ts_data_collector -m 'C1C=CCC=C1' C14H11N   
