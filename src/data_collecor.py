#!/usr/bin/env python3

import os, sys, re, argparse, pdb

parser = argparse.ArgumentParser(
                    prog='TS_pipeline',
                    description='Semi-automatically finds TS for molecular switches',)

parser.add_argument('filename', nargs='+', help='XYZ file to process')
args = parser.parse_args()



if __name__ == '__main__':
    for file in args.filename:
        with open(file, 'r') as f:
            text = f.read()
        text = re.split(r'COMPOUND JOB  ?[0-9]', text)
        G_TS               = float(*re.findall(r'Final Gibbs free.*(-\d+\.\d+) Eh',      text[3]))
#        G_minus_E_TS       = float(*re.findall(r'G-E\(el\).*(-?\d+\.\d+) Eh',              text[3]))
#        E_high_TS          = float(*re.findall(r'FINAL SINGLE POINT ENERGY +(-\d+\.\d+)', text[4]))
        G_react_1          = float(*re.findall(r'Final Gibbs free.*(-\d+\.\d+) Eh',      text[5]))
#        G_minus_E_react_1  = float(*re.findall(r'G-E\(el\).*(-?\d+\.\d+) Eh',              text[5]))
#        E_high_react_1     = float(*re.findall(r'FINAL SINGLE POINT ENERGY +(-\d+\.\d+)', text[6]))
        G_react_2          = float(*re.findall(r'Final Gibbs free.*(-\d+\.\d+) Eh',      text[7]))
#        G_minus_E_react_2  = float(*re.findall(r'G-E\(el\).*(-?\d+\.\d+) Eh',              text[7]))
#        E_high_react_2     = float(*re.findall(r'FINAL SINGLE POINT ENERGY +(-\d+\.\d+)', text[8]))
#        print('G          at TS:', G_TS                              )        
#        print('G-E(el)    at TS:', G_minus_E_TS                      )        
#        print('E high lev at TS:', E_high_TS                         )        
#        print('G          at R :', G_react_1                         )        
#        print('G-E(el)    at R :', G_minus_E_react_1                 )        
#        print('E high lev at R :', E_high_react_1                    )        
#        print('G          at P :', G_react_2                         )        
#        print('G-E(el)    at P :', G_minus_E_react_2                 )        
#        print('E high lev at P :', E_high_react_2                    )        
#        print('=====================================================')
#        print('ΔG(P-R),  kcal/mol:', (G_react_1 - G_react_2) * 627.509)
#        print('ΔG(TS-P), kcal/mol:', (G_TS - G_react_1)      * 627.509)
#        print('ΔG(TS-R), kcal/mol:', (G_TS - G_react_2)      * 627.509)
#        print('=====================================================')
        print('G(abs)TS,hartree:',  G_TS                           )
        print('ΔG(P-R),  kJ/mol:', (G_react_1 - G_react_2) * 2625.5)
        print('ΔG(TS-P), kJ/mol:', (G_TS - G_react_1)      * 2625.5)
        print('ΔG(TS-R), kJ/mol:', (G_TS - G_react_2)      * 2625.5)
#        print('G          at R,  kcal/mol:', (G_minus_E_react_1 + E_high_react_1)*627.509)
#        print('G          at P,  kcal/mol:', (G_minus_E_react_2 + E_high_react_2)*627.509)
        
#        breakpoint()

