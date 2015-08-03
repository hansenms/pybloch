import bloch
import poet_sim
import matplotlib
import matplotlib.pyplot as plt
import scipy.io as sio

seq = poet_sim.read_poet_simulation('/Users/hansenms/temp/pulse_seqs/ssfp_ss_startup_store_mag.txt')
bloch.save_seq(seq,'ssfp_store.p')

seq = poet_sim.read_poet_simulation('/Users/hansenms/temp/pulse_seqs/ssfp_ss_startup.txt')
bloch.save_seq(seq,'ssfp_ss_startup.p')

seq = poet_sim.read_poet_simulation('/Users/hansenms/temp/pulse_seqs/ssfp_ss_tr10.txt')
bloch.save_seq(seq,'ssfp_ss_tr10.p')

seq = poet_sim.read_poet_simulation('/Users/hansenms/temp/pulse_seqs/flash.txt')
bloch.save_seq(seq,'flash.p')
