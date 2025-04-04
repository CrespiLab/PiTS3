#!/usr/bin/bash

cdxml_to_xyz test1.cdxml
fragments_combiner test2_frag1.cdxml test2_frag2.cdxml
tspipe -m 'C-C=N-C' C8H9N.xyz -o /proj/scgrp/users/x_***REMOVED***pe/TS_pipeline/ts_pipeline/templates/orca_three_points_BS_SP_last.inp | tee test.out
data_collector -m 'C-C=N-C' --individual 