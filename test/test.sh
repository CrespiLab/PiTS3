#!/usr/bin/bash

cdxml_to_xyz test1.cdxml
fragments_combiner test2_frag1.cdxml test2_frag2.cdxml
tspipe -m 'C-C=C-C'    C13H11N.xyz          --postpone-orca -o orca_three_points_BS_SP_last.inp | tee test_C13H11N.xyz
tspipe -m 'C-C=N-C'    C14H12.xyz           --postpone-orca -o orca_three_points_BS_SP_last.inp | tee test_C14H12.xyz
tspipe -m 'C1C=CCC=C1' C14H11N.xyz          --postpone-orca -o orca_three_points_BS_SP_last.inp | tee test_C14H11N.xyz
tspipe -m 'C1=CCNC=C1' C19H12F6N2S2.xyz     --postpone-orca -o orca_three_points_BS_SP_last.inp | tee test_C19H12F6N2S2.xyz
                    
ts_data_collector -m 'C-C=C-C'    C13H11N
ts_data_collector -m 'C-C=N-C'    C14H12      
ts_data_collector -m 'C1C=CCC=C1' C14H11N     
ts_data_collector -m 'C1=CCNC=C1' C19H12F6N2S2





