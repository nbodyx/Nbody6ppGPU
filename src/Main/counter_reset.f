      subroutine counter_reset
*
*     Reset counters to zero
*     ---------------------------
*
      include 'common6.h'
*
      NSTEPI = 0 
      NSTEPR = 0 
      NSTEPU = 0 
      NSTEPT = 0 
      NSTEPQ = 0 
      NSTEPC = 0 
      NBLOCK = 0 
      NBLCKR = 0
      NNPRED = 0 
C      NBPRED = 0 
      NIRRF  = 0
      NBCORR = 0 
      NBFLUX = 0
      NBFULL = 0 
      NBVOID = 0 
      NICONV = 0 
      NLSMIN = 0 
      NBSMIN = 0 
      NBDIS  = 0 
      NBDIS2 = 0 
      NCMDER = 0 
C      NBDER  = 0 
      NFAST  = 0 
      NBFAST = 0
C      NBSTAT = 0 
      NKSTRY = 0 
      NKSREG = 0 
C      NEWKS  = 0 
      NKSHYP = 0 
      NKSPER = 0 
C      NPRECT = 0 
C      NKSREF = 0 
      NKSMOD = 0 
      NTTRY  = 0 
      NTRIP  = 0 
      NQUAD  = 0 
      NCHAIN = 0 
      NMERG  = 0 
      NEWHI  = 0 
      NDISS  = 0 
      NTIDE  = 0 
      NSYNC  = 0 
      NCOLL  = 0 
      NCOAL  = 0
      NDD    = 0
      NCIRC  = 0
      NROCHE = 0
      NRO    = 0
      NCE    = 0
      NHYP   = 0
      NHYPC  = 0
      NSESC  = 0 
      NBESC  = 0 
      NMESC  = 0 
      NMDOT  = 0 
      NRG    = 0 
      NHE    = 0 
      NRS    = 0 
      NNH    = 0 
      NWD    = 0 
      NSN    = 0 
      NBH    = 0 
      NBS    = 0 
      ZMRG   = 0
      ZMHE   = 0
      ZMRS   = 0
      ZMNH   = 0
      ZMWD   = 0
      ZMSN   = 0
      ZMDOT  = 0

      RETURN

      END
