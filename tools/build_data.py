"""Deterministic teaching data; no measurements or student records are used."""
import csv
import math
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
DATA=ROOT/'data'
DATA.mkdir(exist_ok=True)
def write(name,columns,rows):
    with (DATA/name).open('w',newline='') as f:
        writer=csv.writer(f); writer.writerow(columns); writer.writerows(rows)

write('amounts.csv',['reagent','role','mass_g','molar_mass_g_mol','stoich'],[
    ['salicylic_acid','reactant',1.3812,138.12,1],
    ['acetic_anhydride','reactant',1.53135,102.09,1]])
write('isolation.csv',['sample_id','product','isolated_mass_g','product_molar_mass_g_mol','product_stoich'],[
    ['P1','aspirin',1.20,180.16,1]])
write('conformers.csv',['conformer','delta_G_J_mol'],[['anti',0],['gauche_combined',3500]])
write('enantiomers.csv',['sample_id','R_percent','S_percent'],[['E1',70,30],['E2',50,50],['E3',95,5]])
write('tlc.csv',['sample_id','spot_distance_cm','front_distance_cm'],[
    ['starting',1.2,5],['product',3.2,5],['crude',1.2,5],['crude',3.2,5]])
write('pka.csv',['compound','pKa','solvent','status'],[
    ['ethanoic_acid',4.8,'water','approximate_teaching_value'],
    ['methylammonium',10.6,'water','approximate_teaching_value']])

# Schematic transmittance traces, generated from Gaussian dips. These are not
# reference spectra and must not be used to identify a real chemical sample.
peaks={
    'ethanol':[(3350,160,52),(2970,55,16),(1050,32,48)],
    'ethanal':[(1725,24,67),(2820,23,15),(2720,22,12),(2970,50,18)],
    'ethanoic_acid':[(3000,280,47),(1710,24,60),(1260,35,36)],
    'ethyl_ethanoate':[(1740,24,68),(1240,40,46),(1050,28,28),(2980,55,19)]}
rows=[]
for sample,bands in peaks.items():
    for wn in range(600,4001,4):
        trans=96-sum(depth*math.exp(-0.5*((wn-centre)/width)**2) for centre,width,depth in bands)
        rows.append([sample,wn,round(max(0,trans),4)])
write('ir.csv',['sample_id','wavenumber_cm1','transmittance_percent'],rows)
write('nmr_peaks.csv',['sample_id','shift_ppm','integration','multiplicity','assignment'],[
    ['unknown_U1',1.25,3,'triplet','ethyl_terminal_CH3'],
    ['unknown_U1',2.05,3,'singlet','acetyl_CH3'],
    ['unknown_U1',4.12,2,'quartet','O_CH2']])
print('Wrote seven small teaching tables and one simulated IR dataset.')
