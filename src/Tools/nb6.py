import numpy as np
import collections

def dict_tree(names,data):
    if (len(names)>1):
        nk = len(names[0])
        if (len(data.shape) < 2): raise Exception("Input data is not two-dimensional array.")
        if (data.shape[0]%nk > 0): 
            raise Exception("Input data columns size %i and names size %i unmatched" %(data.shape[0],nk))
        ncols = data.shape[0]/nk
        dtemp = collections.OrderedDict() 
        for i in range(nk): dtemp[names[0][i]]=dict_tree(names[1:],data[i*ncols:(i+1)*ncols])
        return dtemp
    else:
        nk = len(names[0])
        if (nk!=data.shape[0]): raise Exception("The inner most name size %i and data columns %i unmatched" %(nk,data[0].size))
        dtemp = collections.OrderedDict() 
        for i in range(nk): dtemp[names[0][i]]=data[i] 
        return dtemp

def status(path,seflag):
    data = np.transpose(np.loadtxt(path))
    dtemp = collections.OrderedDict()

    ifoff = 16
    dtemp['t']    = dict_tree((('NB','Myr'),),data[0:2].astype('f'))    # Time in NB unit and Myr
    dtemp['tcr']  = data[2].astype('f')                                 # crossing time in Myr
    dtemp['trh']  = data[3].astype('f')                                 # half-mass relaxation time in Myr
    dtemp['mass'] = dict_tree((('T','S','B'),),data[4:7].astype('f'))   # total, single, binary mass in M_sun
    dtemp['q']    = data[7].astype('f')                                 # Virial ration
    dtemp['rh']   = data[8].astype('f')                                 # half-mass radius in pc
    dtemp['rt']   = data[9].astype('f')                                 # tidal radius in pc
    dtemp['rden'] = dict_tree((('x','y','z'),),data[10:13].astype('f')) # Density center position
    dtemp['rhod'] = data[13].astype('f')                                # Density weighted average density ΣRHO2/ΣRHO 
    dtemp['rhom'] = data[14].astype('f')                                # Maximum mass density / half mass mean value
    dtemp['mmax'] = data[15].astype('f')                                # Maxium stellar mass

    ieoff = 12 + ifoff
    dtemp['energy'] = dict_tree( (('Etot','Ekin','Epot','Ebin','Etid','Em','Ecol','Ece','Ekick','Eesc','Ebesc','Emesc'),), data[ifoff:ieoff].astype('f') )                                                  # Energy parameters
    
    inoff = 5 + ieoff
    dtemp['n']    = dict_tree((('T','S','B','M','Tu'),), data[ieoff:inoff].astype('i'))     # Number of stars Total(resolved),Single,Binary,Triple,Total(unresolved)

    iloff = 486 + inoff
    ffrac=('0.1%','1%','10%','30%','50%','70%','90%','100%','Rc')
    ffracn=('0.1%','1%','10%','30%','50%','70%','90%','100%')
    lpars=('r','n','m','v','vx','vy','vz','vr','vt','vrot','s','sx','sy','sz','sr','st','srot','e')
    dtemp['lagr'] = dict_tree( (('T','S','B'),lpars,ffrac), data[inoff:iloff].astype('f') )      # Lagrangian radii (total, single, binary)    
    
    ibfoff = 16 + iloff
    dtemp['binary'] = dict_tree( (('m',),ffracn), data[iloff:iloff+8].astype('f') )              # binary mass in global lagr
    dtemp['binary'].update(dict_tree( (('n',),ffracn), data[iloff+8:ibfoff].astype('i')) )       # binary number in global lagr    

    ipbfoff = 18 + ibfoff
    dtemp['pbinary'] = dict_tree( (('m',),ffrac), data[ibfoff:ibfoff+9].astype('f') )           # primordinary binary mass in global lagr
    dtemp['pbinary'].update(dict_tree( (('n',),ffrac), data[ibfoff+9:ipbfoff].astype('i')) )    # primordinary binary number in global lagr

    iebfoff = 17 + ipbfoff
    dtemp['binary'].update(dict_tree( (('Ebin',),ffrac), data[ipbfoff:ipbfoff+9].astype('f') ))   # binary binding energy in global lagr (negative)
    dtemp['binary'].update(dict_tree( (('Ebinb',),ffracn), data[ipbfoff+9:iebfoff].astype('f') )) # binary binding energy in binary lagr (negative)

    iepbfoff = 17 + iebfoff
    dtemp['pbinary'].update(dict_tree( (('Ebin',),ffrac), data[iebfoff:iebfoff+9].astype('f') ))    # binary binding energy in global lagr (negative)
    dtemp['pbinary'].update(dict_tree( (('Ebinb',),ffracn), data[iebfoff+9:iepbfoff].astype('f') )) # binary binding energy in binary lagr (negative)
                                     
    iaoff = 27 + iepbfoff
    dtemp['A'] = dict_tree( (('x','y','z'),ffrac), data[iepbfoff:iaoff].astype('f') )   # Angular momentum in global lagr

    if (seflag):
        isfoff = 58 + iaoff
        dtemp['mass'].update(dict_tree( (('SE',),('dM',)), data[iaoff:iaoff+1].astype('f') ))           # stellar mass loss
        # New stellar mass (red giant, helium star, red supergiant, naked helium star, white dwarf, neutron star)
        dtemp['mass']['SE'].update(dict_tree( (('New',),('RG','He','RSG','NHe','WD','SN')), data[iaoff+1:iaoff+7].astype('f') ))  
        # Current stellar mass from KZ type -1 to 15
        dtemp['mass']['SE'].update(dict_tree( (('Current',),('PMS','LMS','HMS','HG','RG','CHB','FAGB','SAGB','HeMS','HeHG','HeGB','HeWD','COWD','ONWD','NS','BH','SNR'),('S','B','BB')), data[iaoff+7:isfoff].astype('f') ))

        isnoff = 72 + isfoff
        # event counts
        #  1. NDISS: Tidal dissipation at pericenter (#27 >0)
        #  2. NTIDE: Tidal captures (#27 >0)
        #  3. NSYNC: Synchronous binaries (#27 >0)
        #  4. NCOLL: stellar collision
        #  5. NCOAL: stellar coalescence
        #  6. NCIRC: circularized binaries (#27 >0)
        #  7. NROCHE: Roche stage triggered times
        #  8. NRO: Roche binary events
        #  9. NCE: Common envelope binaries
        # 10. NHYP: Hyperbolic collision
        # 11. NHYPC: Hyperbolic common envelope binaries
        # 12. NKICK: WD/NS/BH kick
        # 13. NMDOT: Stellar mass loss event
        dtemp['n'].update(dict_tree( (('SE',),('Events',),('Diss','Tid','Syn','Coll','Coal','Circ','Roche','RoBin','CE','Hcoll','HCE','Kick','dM')), data[isfoff:isfoff+13].astype('i') ))
        # New stellar types
        # 14. NRG: New red giants
        # 15. NHE: New helium stars
        # 16. NRS: New red supergiants
        # 17. NNH: New naked helium stars
        # 18. NWD: New white dwarfs
        # 19. NSN: New neutron stars
        # 20. NBH: New black holes
        # 21. NBS: New blue stragglers
        dtemp['n']['SE'].update(dict_tree( (('New',),('RG','He','RSG','NHe','WD','NS','BH','BS')), data[isfoff+13:isfoff+21].astype('i') ))
        # Current stellar number from KZ type -1 to 15
        dtemp['n']['SE'].update(dict_tree( (('Current',),('PMS','LMS','HMS','HG','RG','CHB','FAGB','SAGB','HeMS','HeHG','HeGB','HeWD','COWD','ONWD','NS','BH','SNR'),('S','B','BB')), data[isfoff+21:isnoff].astype('i') ))

        isloff = 1782 + isnoff
        # lagraigan radii parameters for stellar types
        # 1.  Low mass main sequence (M < 0.7) (0)
        # 2.  High mass main sequence  (1)
        # 3.  Hertzsprung gap (HG). (2)
        # 4.  Red giant. (3)
        # 5.  Core Helium burning. (HB) (4)
        # 6.  AGB (5-6)
        # 7.  Helium types (7-9)
        # 8.  White dwarf (10-12)
        # 9.  Neutron star (13)
        # 10. Black hole (14)
        # 11. Pre main sequence (-1)
        dtemp['lagr'].update(dict_tree( (('SE',),('LMS','HMS','HG','RG','CHB','AGB','He','WD','NS','BH','PMS'),lpars,ffrac), data[isnoff:isloff].astype('f') ))
        
        return dtemp
