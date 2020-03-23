#!/usr/bin/env python3

import numpy as np
import multiprocessing as mp
import petar
import time
import sys
import getopt
import os

class NBSingle(petar.SimpleParticle):

    def __init__(self, _dat=None, _offset=int(0), _append=False, **kwargs):
        keys = [['id',1]]
        petar.DictNpArrayMix.__init__(self, keys, _dat, _offset, _append, **kwargs)
        petar.SimpleParticle.__init__(self, _dat, _offset+self.ncols, True, **kwargs)
        keys = [['pot',1]]
        petar.DictNpArrayMix.__init__(self, keys, _dat, _offset+self.ncols, True, **kwargs)

    def scaleToNB(self, **kwargs):
        if 'scale' in kwargs.keys():
            rscale=kwargs['scale'][0]
            mscale=kwargs['scale'][1]
            vscale=kwargs['scale'][2]
            self.mass /= mscale
            self.pos  /= rscale
            self.vel  /= vscale

    def calcEtot(self):
        if (not 'etot' in self.__dict__.keys()): 
            self.ncols += 1
            self.keys.append(['etot',1])
        self.etot = self.ekin + self.mass*self.pot

class NBBinary(petar.DictNpArrayMix):
    def __init__(self, _dat=None, _offset=int(0), _append=False, **kwargs):
        keys = [['id1',1],['id2',1],['id',1],['m1',1],['m2',1],['pos',3], ['vel',3], ['dpos_in',3], ['dvel_in',3], ['pot',1], ['semi',1], ['ecc',1], ['period',1]]
        petar.DictNpArrayMix.__init__(self, keys, _dat, _offset, _append, **kwargs)

    def convertToSingle(self):
        bm = self.m1+self.m2
        bmr1 = ( self.m2/bm).reshape((bm.size,1))
        bmr2 = (-self.m1/bm).reshape((bm.size,1))
        bmr1a = bmr1/206265.0
        bmr2a = bmr2/206265.0

        single=NBSingle()
        single.id = np.append(self.id1,self.id2)
        single.size = single.id.size
        single.mass = np.append(self.m1,self.m2)
        single.pos = np.append(self.pos+self.dpos_in*bmr1a, self.pos+self.dpos_in*bmr2a, axis=0)
        invr = 1.0/np.sqrt(petar.vecDot(self.dpos_in,self.dpos_in))
        single.vel = np.append(self.vel+self.dvel_in*bmr1,  self.vel+self.dvel_in*bmr2, axis=0)
        single.pot = np.append(self.pot-self.m2*invr,  self.pot-self.m1*invr)

        return single

class NBMerger(petar.DictNpArrayMix):
    def __init__(self, _dat=None, _offset=int(0), _append=False, **kwargs):
        keys = [['id1',1],['id2',1],['id3',1],['id',1],['m1',1],['m2',1],['m3',1],['pos',3], ['vel',3], ['dpos_out',3], ['dvel_out',3], ['dpos_in',3], ['dvel_in',3], ['pot',1], ['semi_out',1], ['ecc_out',1], ['period_out',1], ['semi_in',1], ['ecc_in',1], ['period_in',1]]
        petar.DictNpArrayMix.__init__(self, keys, _dat, _offset, _append, **kwargs)

    def convertToSingle(self):
        micm = self.m1 + self.m2
        mm = micm + self.m3
        micmr = (self.m3/mm).reshape((micm.size,1))
        m3r   = (-micm/mm).reshape((micm.size,1))
        micmra = micmr/206265.0
        m3ra =   m3r/206265.0
        mxicm = self.pos + micmra*self.dpos_out
        mx3   = self.pos + m3ra*self.dpos_out
        mvicm = self.vel + micmr*self.dvel_out
        mv3   = self.vel + m3r*self.dvel_out
        mx1r  = (self.m2/micm).reshape((micm.size,1))
        mx2r  = (-self.m1/micm).reshape((micm.size,1))
        mx1ra =  mx1r/206265.0
        mx2ra =  mx2r/206265.0
        mx1   = mxicm + mx1ra*self.dpos_in
        mx2   = mxicm + mx2ra*self.dpos_in
        mv1   = mvicm + mx1r*self.dvel_in
        mv2   = mvicm + mx2r*self.dvel_in
        
        single=NBSingle()
        single.id = np.concatenate((self.id1,self.id2,self.id3))
        single.size = single.id.size
        single.mass = np.concatenate((self.m1,self.m2,self.m3))
        single.pos  = np.concatenate((mx1, mx2, mx3))
        single.vel  = np.concatenate((mv1, mv2, mv3))
        dx13 = mx1-mx3
        dx23 = mx2-mx3
        invr13 = 1.0/np.sqrt(petar.vecDot(dx13,dx13))
        invr23 = 1.0/np.sqrt(petar.vecDot(dx23,dx23))
        invrin = 1.0/np.sqrt(petar.vecDot(self.dpos_in,self.dpos_in))
        single.pot  = np.concatenate((self.pot-self.m2*invrin-self.m3*invr13, self.pot-self.m1*invrin-self.m3*invr23, self.pot-self.m1*invr13-self.m2*invr23))

        return single

def dataProcessOne(file_path, lagr, core, esc, time_profile, **kwargs): 
    m_frac = lagr.initargs['mass_fraction']
    G=1.0
    r_bin=0.1
    average_mode='sphere'

    if ('G' in kwargs.keys()): G=kwargs['G']
    if ('r_max_binary' in kwargs.keys()): r_bin=kwargs['r_max_binary']
    if ('average_mode' in kwargs.keys()): average_mode=kwargs['average_mode']

    start_time = time.time()
    #print('Loadfile')
    path,t=file_path.split()
    particle = NBSingle()
    particle.loadtxt(path+'single.40_'+t)
    if (os.path.isfile(path+'binary.40_'+t)):
        binary=NBBinary()
        binary.loadtxt(path+'binary.40_'+t)
        single=binary.convertToSingle()
        particle.append(single)
    if (os.path.isfile(path+'merger.40_'+t)):
        merger=NBMerger()
        merger.loadtxt(path+'merger.40_'+t)
        particle.append(merger.convertToSingle())
    particle.scaleToNB(**kwargs)
    read_time = time.time()

    # find binary
    #print('Find pair')
    kdtree,single,binary=petar.findPair(particle,G,r_bin,True)
    find_pair_time = time.time()
    
    # get cm, density
    #print('Get density')
    cm_pos, cm_vel=core.calcDensityAndCenter(particle,kdtree)
    #print('cm pos:',cm_pos,' vel:',cm_vel)
    get_density_time = time.time()

    #print('Correct center')
    particle.correctCenter(cm_pos, cm_vel)

    # r2
    particle.calcR2()
    # rc

    #print('Core radius')
    rc = core.calcCoreRadius(particle)
    #print('rc: ',rc)

    core.addTime(float(t))
    core.size+=1

    n_frac=m_frac.size+1
    cm_vel=np.array([0,0,0]) # avoid kinetic energy jump 
    single.correctCenter(cm_pos, cm_vel)
    binary.correctCenter(cm_pos, cm_vel)
    center_and_r2_time = time.time()

    single.savetxt('data.'+t+'.single')
    binary.savetxt('data.'+t+'.binary')

    #print('Lagrangian radius')
    lagr.calcOneSnapshot(float(t), single, binary, rc, average_mode)
    lagr_time = time.time()

    rhindex=np.where(m_frac==0.5)[0]
    esc.calcRCutIsolate(lagr.all.r[-1,rhindex])
    esc.findEscaper(float(t),single,binary,G)

    time_profile['read'] += read_time-start_time
    time_profile['find_pair'] += find_pair_time-read_time
    time_profile['density'] += get_density_time-find_pair_time
    time_profile['center_core'] += center_and_r2_time-get_density_time
    time_profile['lagr'] += lagr_time-center_and_r2_time

    return time_profile

def dataProcessList(file_list, **kwargs):
    """ process lagragian calculation for a list of file snapshots
    file_list: file path list
    """
    lagr=petar.LagrangianMultiple(**kwargs)
    time_profile=dict()
    time_profile['read'] = 0.0
    time_profile['find_pair'] = 0.0
    time_profile['density'] = 0.0   
    time_profile['center_core'] = 0.0   
    time_profile['lagr'] = 0.0
    core=petar.Core()
    esc=petar.Escaper(NBSingle)
    for path in file_list:
        #print(' data:',path)
        dataProcessOne(path, lagr, core, esc, time_profile, **kwargs)
    for key, item in time_profile.items():
        item /= len(file_list)
    return lagr, core, esc, time_profile


def parallelDataProcessList(file_list, n_cpu=int(0), **kwargs):
    """ parellel process lagragian calculation for a list of file snapshots
    file_list: file path list
    """
    if (n_cpu==int(0)):
        n_cpu = mp.cpu_count()
        #print('n_cpu:',n_cpu)
    pool = mp.Pool(n_cpu)

    n_files=len(file_list)
    n_pieces = np.ones(n_cpu)*int(n_files/n_cpu)
    n_left = n_files%n_cpu
    n_pieces[:n_left]+=1
    n_offset=np.append([0],n_pieces.cumsum()).astype(int)
    #print('work_pieces',n_pieces)

    file_part = [file_list[n_offset[i]:n_offset[i+1]] for i in range(n_cpu)]

    result=[None]*n_cpu
    for rank in range(n_cpu):
        result[rank] = pool.apply_async(dataProcessList, (file_part[rank],), kwargs)

    # Step 3: Don't forget to close
    pool.close()

    lagri=[]
    corei=[]
    esci=[]
    time_profilei=[]
    for i in range(n_cpu):
        lagri.append(result[i].get()[0])
        corei.append(result[i].get()[1])
        esci.append(result[i].get()[2])
        time_profilei.append(result[i].get()[3])
    lagr = petar.join(*lagri)
    core = petar.join(*corei)
    esc  = petar.joinEscaper(*esci)
    time_profile=time_profilei[0]
    for key in time_profile.keys():
        for i in range(1,n_cpu):
            time_profile[key] += time_profilei[i][key]/n_cpu
    return lagr, core, esc, time_profile

if __name__ == '__main__':

    filename='dat.lst'
    lagr_filename='lagr.dat'
    core_filename='core.dat'
    esc_prefix='esc'
    average_mode='sphere'
    n_cpu=0

    def usage():
        print("A tool for processing a list of snapshot data to detect binaries, calculate Langragian radii and properties, get the density center and core radius")
        print("Usage: petar.data.process [options] data_filename")
        print("data_filename: A list of snapshot data path and suffix(time[NB]), each line for one snapshot")
        print("option:")
        print("  -h(--help): help")
        print("  -l(--lagr-filename): Lagrangian radii data filename (lagr.dat)")
        print("  -m(--mass-fraction): Lagrangian radii mass fraction (0.1,0.3,0.5,0.7,0.9)")
        print("  -G(--gravitational-constant): Gravitational constant (1.0)")
        print("  -r(--r-max-binary): maximum sepration for detecting binaries (0.1)")
        print("  -a(--average-mode): Lagrangian properity average mode: sphere: average from center to Lagragian radii; shell: average between two neighbor radii (sphere)")
        print("  -c(--core-filename): core data (time, density center, core radius) filename (core.dat)")
        print("  -e(--esc-filename): esc data filename prefix (esc)")
        print("  -s(--scale): scale factor from NB to Astro: r[pc],m[Msun],v[kms/s],t[Myr] (1.0,1.0,1.0,1.0)")
        print("  -n(--n-cpu): number of CPU threads for parallel processing (all threads)")

    try:
        shortargs = 'l:m:G:r:a:c:e:s:n:h'
        longargs = ['lagr-filename=','mass-fraction=','gravitational-constant=','r-max-binary=','average-mode=', 'core-filename=','esc-filename=','scale=','n-cpu=','help']
        opts,remainder= getopt.getopt( sys.argv[1:], shortargs, longargs)

        kwargs=dict()
        for opt,arg in opts:
            if opt in ('-h','--help'):
                usage()
                sys.exit(1)
            elif opt in ('-l','--lagr-filename'):
                lagr_filename = arg
            elif opt in ('-m','--mass-fraction'):
                kwargs['mass_fraction'] = np.array([float(x) for x in arg.split(',')])
            elif opt in ('-G','--gravitational-constant'):
                kwargs['G'] = float(arg)
            elif opt in ('-r','--r-max-binary'):
                kwargs['r_max_binary'] = float(arg)
            elif opt in ('-a','--average-mode'):
                kwargs['average_mode'] = arg
            elif opt in ('-c','--core-filename'):
                core_filename = arg
            elif opt in ('-e','--esc-filename'):
                esc_filename = arg
            elif opt in ('-s','--scale'):
                kwargs['scale'] = np.array([float(x) for x in arg.split(',')])
            elif opt in ('-n','--n-cpu'):
                n_cpu = int(arg)
            else:
                assert False, "unhandeld option"

    except getopt.GetoptError:
        print('getopt error!')
        usage()
        sys.exit(1)

    filename = remainder[0]

    for key, item in kwargs.items(): print(key,':',item)

    fl = open(filename,'r')
    file_list = fl.read()
    path_list = file_list.splitlines()
     
    lagr,core,esc,time_profile = parallelDataProcessList(path_list,n_cpu=n_cpu,**kwargs)
     
    lagr.savetxt(lagr_filename)
    core.savetxt(core_filename)
    esc.single.savetxt(esc_prefix+'.single')
    esc.binary.savetxt(esc_prefix+'.binary')
     
    print ('Time profile:')
    for key, item in time_profile.items():
        print (key,item,)
     
    print ("Lagr data is saved in file: ",lagr_filename, ' core data is saved in file: ', core_filename, ' esc data are saved in file: ', esc_prefix,'.[single/binary]')

