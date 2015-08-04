import matplotlib
matplotlib.use('Agg')

import bloch
import numpy as np

T1 = 1000.0e-3
T2 = 50.0e-3

seq = bloch.load_seq('ssfp_store.p')

TEs = bloch.find_echo_times(seq)
sim_echos = bloch.ZSimulator(seq, T1, T2, sample_times=TEs)

with bloch.Timer('Single Z echos'):
    Mxy, Mz = sim_echos(0.0)

bloch.plot_simulation(seq,Mxy,sample_times=TEs, filename='single_z_echo_times.pdf')
