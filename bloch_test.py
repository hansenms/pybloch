import matplotlib
import matplotlib.pyplot as plt
import bloch
import poet_sim
import numpy as np
import cmath
import multiprocessing

%matplotlib inline

reload(bloch)
reload(poet_sim)

seq = poet_sim.read_poet_simulation('/Users/hansenms/temp/pulse_seqs/ssfp_ss_startup.txt')
plt.plot(np.abs(seq[:,1]))

TEs = poet_sim.find_echo_times(seq)

seq_short = seq[:,:]
zpos = np.linspace(-10.0e-3,10e-3,100)
sim = bloch.ZSimulator(seq_short, 1000.0e-3, 50.0e-3, sample_times=TEs)
#sim = bloch.ZSimulator(seq_short, 1000.0e-3, 50.0e-3)

with bloch.Timer('Single z'):
    Mxy, Mz = sim(0.0)

plt.plot(np.abs(Mxy))

number_of_processes = multiprocessing.cpu_count()

with bloch.Timer('Multiprocessing map'):
    pool = multiprocessing.Pool(number_of_processes)
    tmp = pool.map(sim,zpos)
    Mxy_reduced, Mz_reduced = reduce(lambda x,y: (x[0]+y[0], x[1]+y[1]),tmp)
    Mxy_reduced = Mxy_reduced/zpos.shape[0]
    Mz_reduced = Mz_reduced/zpos.shape[0]

plt.plot(np.abs(Mxy_reduced))
