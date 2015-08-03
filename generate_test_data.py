import bloch
import poet_sim
import matplotlib
import matplotlib.pyplot as plt
import scipy.io as sio

seq = poet_sim.read_poet_simulation('/Users/hansenms/temp/pulse_seqs/ssfp_ss_startup.txt')
seq[:,0] = 6.146890*seq[:,0]
bloch.save_seq(seq,'ssfp_new.p')

plt.plot(np.abs(seq[92000:96000,0]))

70.0/(360.0*(np.sum(np.abs(seq[92000:96000,0]),0)*1.0e-6*42.576e6)/(np.pi*2))
