*     KS MPI communication flags parameters.
*     ------------------
*     
*     Index the variables as below, (n) means need n space, first n-1 is index for variable array, last is the value

**   INTEGER:
***   single variables:
****   1      2      3     4      5      6      7
****   DELAY  JCLOSE JCMAX IQCOLL IPHASE KSPAIR KICK
***   1 dimension array:
****   100
****   KSTAR

**   REAL*8:
***   single variables:
****   1      2     3     4     5     6     7     8     9     10   
****   EMERGE BE(3) ETIDE ECOLL EGRAV E(10) EBCH0 ECDOT EKICK RMAX
***   1 dimension array:
****   100    101    102   103          104
****   TEV    RADIUS SPIN  STEPS(CLUMP) TMDIS(BINARY)
***   2 dimension array:      
****   200    201        202
****   X0DOT  CM(BINARY) CM(MODES)

**     argument meanings:
***     Operator: 1: save 2: communicate 3: reset
***     par_type: 1: integer 2: real*8

*     Operator flags
      parameter (K_store=1,K_comm=2,K_reset=3)
*     Parameter types flags
      parameter (K_int=1,K_real8=2)
      
*     Splitter for different dimensions
      parameter (IISPLIT=100,IRSPLIT=100,IRSPLIT2=200)

*     Integer flag
      parameter (K_DELAY=1,K_JCLOSE=2,K_JCMAX=3,K_IQCOLL=4,K_IPHASE=5,
     &     K_KSPAIR=6,K_KICK=7,K_KSTAR=100,
*     REAL flag
     &     K_EMERGE=1,K_BE3=2,K_ETIDE=3,K_ECOLL=4,K_EGRAV=5,K_E10=6,
     &     K_EBCH0=7,K_ECDOT=8,K_EKICK=9,K_RMAX=10,
     &     K_TEV=100,K_RADIUS=101,K_SPIN=102,K_STEPS=103,K_TMDIS=104,
     &     K_X0DOT=200,K_CM_B=201,K_CM_M=202)
