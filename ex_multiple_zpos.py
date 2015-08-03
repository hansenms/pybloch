import matplotlib
matplotlib.use('Agg')

import bloch
import numpy as np
import multiprocessing

T1 = 1000.0e-3
T2 = 50.0e-3

seq = bloch.load_seq('ssfp_ss_startup.p')
TEs = bloch.find_echo_times(seq)

zpos = np.linspace(-10.0e-3,10e-3,100)
sim = bloch.ZSimulator(seq, T1, T2, sample_times=TEs)

number_of_processes = multiprocessing.cpu_count()

with bloch.Timer('Multiprocessing using map'):
    pool = multiprocessing.Pool(number_of_processes)
    tmp = pool.map(sim,zpos)
    Mxy_reduced, Mz_reduced = reduce(lambda x,y: (x[0]+y[0], x[1]+y[1]),tmp)
    Mxy_reduced = Mxy_reduced/zpos.shape[0]
    Mz_reduced = Mz_reduced/zpos.shape[0]

bloch.plot_simulation(seq,Mxy_reduced,sample_times=TEs, filename='multiple_z_ssfp.pdf')
