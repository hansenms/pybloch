import matplotlib
matplotlib.use('Agg')

import bloch
import numpy as np

T1 = 1000.0e-3
T2 = 50.0e-3

seq = bloch.load_seq('ssfp_store.p')

TEs = bloch.find_echo_times(seq)

sim_full = bloch.ZSimulator(seq, T1, T2)
sim_echos = bloch.ZSimulator(seq, T1, T2, sample_times=TEs)

with bloch.Timer('Single Z full'):
    Mxy, Mz = sim_full(0.0)

bloch.plot_simulation(seq,Mxy,filename='single_z_full.pdf')

with bloch.Timer('Single Z echos'):
    Mxy, Mz = sim_echos(0.0)

offset = 400000 #We will plot from 400ms and on
bloch.plot_simulation(seq[offset:,:],Mxy,sample_times=TEs-offset, filename='single_z_echo_times.pdf')
