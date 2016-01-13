*     TIMING header
*
*
*     COMMON TLIST for irregular block time step list (Need include params.h before include tlist.h)
*
      parameter(NDELAYMAX=20*LMAX)
      COMMON/TLIST/ NXTLST(NMAX),NXTLEN,NDTK(64),NDTMIN,NDTMAX,
     &              NXTLIMIT,NXTLEVEL,NLSTDELAY(NDELAYMAX),NGHOSTS
      INTEGER NXTLST,NXTLEN,NDTK,NDTMIN,NDTMAX,NXTLIMIT,NXTLEVEL
      INTEGER NLSTDELAY,NGHOSTS

*     NXTLST:   Step sorted index list
*     NXTLEN:   Length of nxtlst of next block
C*     NXTK(K):  step level of nxtlst(K)
*     NDTK(K):  End position of step level DTK(K) in nxtlst
*               (Be carefule here larger NDTK means smaller step), 
*               NDTK+1 will be next step level
*     NDTMIN:   Current Minimum step level 
*     NDTMAX:   Current Maximum step level
*     NXTLIMIT: Current total number of particles (sorting limit of nxtlst)
*     NXTLEVEL: Next integrating block step level
*     NLSTDELAY: Delay list for adding particles into NXTLST, this is for safe of NXTLST order
*     NGHOSTS:  Ghost particle numbers
