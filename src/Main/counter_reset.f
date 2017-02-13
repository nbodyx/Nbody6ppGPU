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
C      NDISS  = 0 
C      NTIDE  = 0 
C      NSYNC  = 0 
C      NCOLL  = 0 
C      NCOAL  = 0
C      NDD    = 0
C      NCIRC  = 0
C      NROCHE = 0
C      NRO    = 0
C      NCE    = 0
C      NHYP   = 0
C      NHYPC  = 0
C      NSESC  = 0 
C      NBESC  = 0 
C      NMESC  = 0 
C      NMDOT  = 0 
C      NRG    = 0 
C      NHE    = 0 
C      NRS    = 0 
C      NNH    = 0 
C      NWD    = 0 
C      NSN    = 0 
C      NBH    = 0 
C      NBS    = 0 
C      ZMRG   = 0
C      ZMHE   = 0
C      ZMRS   = 0
C      ZMNH   = 0
C      ZMWD   = 0
C      ZMSN   = 0
C      ZMDOT  = 0

      RETURN

      END
