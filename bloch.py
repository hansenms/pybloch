import numpy as np
import time

def zrot(phi):
    R = np.matrix([[np.cos(phi), -np.sin(phi), 0],[np.sin(phi), np.cos(phi), 0],[0, 0, 1]],dtype=np.float32)
    return R

def yrot(phi):
    R = np.matrix([[np.cos(phi), 0, np.sin(phi)], [0, 1, 0],[-np.sin(phi), 0, np.cos(phi)]],dtype=np.float32)
    return R

def xrot(phi):
    R = np.matrix([[1, 0, 0], [0, np.cos(phi), -np.sin(phi)],[0, np.sin(phi), np.cos(phi)]],dtype=np.float32)
    return R

def throt(phi,theta):
    Rz = zrot(-theta)
    Rx = xrot(phi)
    R = np.linalg.inv(Rz)*Rx*Rz
    return R

#Longitudinal relaxation
#
#    Mz(t) = M0 + (Mz(0)-M0)*exp(t/T1)
#    Mz(t) = M0 + Mz(0)*exp(t/T1)-M0*exp(t/T1)
#    Mz(t) = Mz(0)*exp(t/T1) + (M0-M0*exp(t/T1))

def freeprecess(T,T1,T2,df,M0=1.0):
    phi = 2*np.pi*df*T;
    E1 = np.exp(-T/T1);
    E2 = np.exp(-T/T2);

    Afp = np.matrix([[E2, 0, 0], [0, E2, 0], [0, 0, E1]], dtype=np.float32)*zrot(phi)
    Bfp = np.array([[0],[0],[M0-M0*E1]], dtype=np.float64)
    return (Afp, Bfp)

# M1 is the magnetization before flip
# M2 is magnetization after flip
# M3 is magnetization at echo time
# Then
#   M2 = Ry * M1
#   M3 = Afp_te* M2 + Bfp_te
#   M1 = Afp_tr*M3 + Bfp_tr
#   M2 = Ry*(Afp_tr*M3 + Bfp_tr)
#   M3 = Afp_te*(Ry*(Afp_tr*M3 + Bfp_tr)) + Bfp_te
#   M3 = Afp_te*Ry*Afp_tr*M3 + Afp_te*Ry*Bfp_tr + Bfp_te
#   (eye(1)-Afp_te*Ry*Afp_tr)*M3 = Afp_te*Ry*Bfp_tr + Bfp_te
#   Mss = np.linalg.inv((np.eye(3)-Afp_te*Ry*Afp_tr))*(Afp_te*Ry*Bfp_tr+Bfp_te)
def sssignal(flip,T1,T2,TE,TR,dfreq):
    Ry = yrot(flip)

    (Afp_te, Bfp_te) = freeprecess(TE,T1,T2,df)
    (Afp_tr, Bfp_tr) = freeprecess(TR-TE,T1,T2,df)

    Mss = np.linalg.inv((np.eye(3)-Afp_te*Ry*Afp_tr))*(Afp_te*Ry*Bfp_tr+Bfp_te)

    return Mss

def srsignal(flip,T1,T2,TE,TR,df):
    Ry = yrot(flip)

    (Afp_te, Bfp_te) = freeprecess(TE,T1,T2,df)
    (Afp_tr, Bfp_tr) = freeprecess(TR-TE,T1,T2,df)
    Afp_tr = np.matrix([[0, 0, 0], [0, 0, 0], [0, 0, 1]], dtype=np.float32)*Afp_tr
    Mss = np.linalg.inv((np.eye(3)-Afp_te*Ry*Afp_tr))*(Afp_te*Ry*Bfp_tr+Bfp_te)

    return Mss

def sesignal(T1,T2,TE,TR,df = 0.0):
    Rflip = yrot(np.pi/2)
    Rrefoc = xrot(np.pi)

    (Atr,Btr) = freeprecess(TR-TE,T1,T2,df)
    (Ate2,Bte2) = freeprecess(TE/2,T1,T2,df)
    Atr = np.matrix([[0, 0, 0], [0, 0, 0], [0, 0, 1]], dtype=np.float32)*Atr

    Mss = np.linalg.inv(np.eye(3)-Ate2*Rrefoc*Ate2*Rflip*Atr) * (Bte2+Ate2*Rrefoc*(Bte2+Ate2*Rflip*Btr))

    return Mss

def ssfpsignal(flip,T1,T2,TE,TR,df,phi=np.pi):
    Ry = yrot(flip)

    (Afp_te, Bfp_te) = freeprecess(TE,T1,T2,df)
    (Afp_tr, Bfp_tr) = freeprecess(TR-TE,T1,T2,df)
    Afp_tr = zrot(phi)*Afp_tr
    Mss = np.linalg.inv((np.eye(3)-Afp_te*Ry*Afp_tr))*(Afp_te*Ry*Bfp_tr+Bfp_te)

    return Mss

def find_relevant_events(seq,z = 0.0, y = 0.0, x = 0.0, sample_times = None):

    if sample_times is None:
        return np.ones((seq.shape[0],))
    else:
        sample_times = np.extract(sample_times < seq.shape[0],sample_times)

    mask = np.reshape(np.array([1.0,z,y,x,0.0]),(1,5))
    sample_events = np.zeros((seq.shape[0],))
    sample_events[sample_times] = 1

    seq_events = (np.sum(np.abs(seq * mask),1) > 0.0) + (sample_events > 0)
    return seq_events

def bloch_sim(seq, T1, T2,  M = 1.0, df = 0.0, z = 0.0, y = 0.0, x = 0.0, dt = 1.0e-6, gamma = 42.576e6, sample_times = None):
    [Afp, Bfp] = freeprecess(dt,T1,T2,df)
    time_points = seq.shape[0]

    #Figure out when we actually have something relevant going on
    seq_events = find_relevant_events(seq,z,y,x,sample_times)

    if sample_times == None:
        sample_times = range(time_points)
    else:
        sample_times = np.extract(sample_times < seq.shape[0],sample_times)

    output_length = len(sample_times)
    Mz = np.zeros((output_length,), dtype=np.float32)
    Mx = np.zeros((output_length,), dtype=np.float32)
    My = np.zeros((output_length,), dtype=np.float32)

    next_sample_idx = 0

    M_current = np.array([[0.0], [0.0], [M]])

    total_flip = 0.0;

    t = 0
    while t < time_points:
        #Any RF excitation
        if (np.abs(seq[t,0]) > 0.0):
            flip = np.abs(seq[t,0])*gamma*dt*np.pi*2
            Rrf = throt(flip,np.angle(seq[t,0]))
            M_current = Rrf*M_current

        #Z-gradient
        if (np.abs(z) > 0.0) and (np.abs(seq[t,1]) > 0.0):
            #print "Apply z rotation"
            Rz = zrot(z*seq[t,1]*gamma*dt*np.pi*2)
            M_current = Rz * M_current

        #Y-gradient
        if (np.abs(y) > 0.0) and (np.abs(seq[t,2]) > 0.0):
            #print "Apply z rotation"
            Ry = yrot(y*seq[t,2]*gamma*dt*np.pi*2)
            M_current = Ry * M_current

        #X-gradient
        if (np.abs(x) > 0.0) and (np.abs(seq[t,3]) > 0.0):
            #print "Apply z rotation"
            Rx = xrot(x*seq[t,3]*gamma*dt*np.pi*2)
            M_current = Rx * M_current

        #Relaxation
        M_current = Afp*M_current + Bfp

        if t == sample_times[next_sample_idx]:
            Mz[next_sample_idx] = M_current[2]
            Mx[next_sample_idx] = M_current[0]
            My[next_sample_idx] = M_current[1]

            next_sample_idx = next_sample_idx + 1

            #We have all the samples we need. No need to simutale anymore
            if (next_sample_idx >= output_length):
                break;

        #Now let's try to look ahead to see if we can just freeprecess for a
        #few samples, if there are no relevant sequence events and no samples

        t2 = 1;
        while not (t+t2 <= time_points) and seq_events[t+t2]:
            t2 = t2 + 1

        #We can just free precess here
        if (t2 > 1):
            [Afp_noevents, Bfp_noevents] = freeprecess(dt*(t2-1),T1,T2,df)
            M_current = Afp_noevents*M_current + Bfp_noevents
            t = t + (t2-1)

        t = t + 1

    return (Mx + 1j*My, Mz)

class ZSimulator(object):
    def __init__(self, seq, T1, T2,  M = 1.0, df = 0.0, y = 0.0, x = 0.0, dt = 1.0e-6, gamma = 42.576e6, sample_times = None):
        self.seq = seq
        self.T1 = T1
        self.T2 = T2
        self.M = M
        self.df = df
        self.y = y
        self.x = x
        self.dt = dt
        self.gamma = gamma
        self.sample_times = sample_times
    def __call__(self, zp):
        return bloch_sim(self.seq, self.T1, self.T2,  M = self.M, df = self.df, z = zp, y = self.y, x = self.x, dt = self.dt, gamma = self.gamma, sample_times = self.sample_times)


class Timer(object):
    def __init__(self, name=None):
        self.name = name

    def __enter__(self):
        self.tstart = time.time()

    def __exit__(self, type, value, traceback):
        if self.name:
            print '[%s]' % self.name,
        print 'Elapsed: %s' % (time.time() - self.tstart)
