      SUBROUTINE custom_output(vec_idx, vec_len, vec_offset)
*
*
*     Custom output (HDF5/CSV) for NBODY6/NBODY6++.
*     -----------------------------------------------
*     Developed by Maxwell Xu CAI (NAOC/KIAA)
*     For bug report/feedback, please email maxwellemail@gmail.com
*
*     Version 2.0 -- Bug fixes & finetuning, June 12, 2014
*     Version 1.4 -- June 11, 2014
*     Version 1.0 -- Initial release, Nov 2013
*

#ifdef H5OUTPUT
      USE HDF5
      INCLUDE 'common6.h'
      COMMON/STEVOL/ XLUM(NMAX),XRSTAR(NMAX),XTEFF(NMAX),
     &               XZMET(NMAX),XCMASS(NMAX),XRC(NMAX)
      COMMON/POTDEN/  RHO(NMAX),XNDBL(NMAX),PHIDBL(NMAX)
      INTEGER :: vec_len ! The length of the current vector passed thought the ar g     
      INTEGER :: vec_offset
      INTEGER :: vec_idx(vec_len+vec_offset-1)
      REAL*8 :: h5_current_time
      REAL*8 :: current_time
      REAL*8 :: time_step
      REAL*8 :: XLUM, XRSTAR, XTEFF, XZMET, XCMASS, XRC
      REAL*8 :: RHO, XNDBL, PHIDBL
      logical first
      data first/.true./
      SAVE h5_current_time,first

*       Calculate the current time from the NBODY6++ global variable
      current_time = TIME

*     Initialization of h5_current_time (L. Wang)
      if(first) then
         IF(KZ(47).GT.0) THEN
            h5_current_time = REAL(INT(TIME+TOFF) +
     &           INT((TIME-INT(TIME))/DTK(KZ(47)+1))*DTK(KZ(47)+1))
         ELSE
            h5_current_time = TIME + TOFF
         END IF
         first = .false.
      end if

      IF(current_time+toff .GE. h5_current_time) THEN
          IF(KZ(47).GT.0) THEN
              time_step = DTK(KZ(47)+1)
          ELSE
              time_step = 0
          ENDIF
*     Call the write back subroutine
      IF(KZ(46).EQ.1 .OR. KZ(46).EQ.3) THEN 
          CALL HDF5_output(vec_idx, vec_len, vec_offset,X0,X0dot,X,XDOT,
     &         f,fdot,T0,BODY,KZ,NAME,current_time, time_step,toff, N, 
     &         kstar, rho, phidbl, xlum, xrstar, xteff,xcmass,xrc,
     &         .TRUE.)
      ELSE IF(KZ(46).EQ.2 .OR. KZ(46).EQ.4) THEN
          CALL CSV_output(vec_idx, vec_len, vec_offset,X0,X0dot,x,xdot,
     &         f,fdot,T0,BODY,KZ,NAME,current_time, time_step,toff, N, 
     &         kstar, rho, phidbl, xlum, xrstar, xteff,xcmass,xrc)
      ENDIF
*===============WRITE BACK DONE========================
*       Leap the h5_current_time forward, if KZ(47) != 0
          IF(KZ(47) .GT. 0) THEN
             h5_current_time = DTK(KZ(47)+1) + h5_current_time
          ELSE ! KZ(47) == 0, output ALL
             h5_current_time = current_time + TOFF
          ENDIF
      ELSE

      ENDIF
#endif
      END

*-----------------------------------------------------------------------*
      SUBROUTINE CSV_output(vec_idx,vec_len,vec_offset,X0,X0dot,X,XDOT,
     &     f,fdot,T0,MASS,KZ,id_name, current_time, time_step,
     &     toff, total_n,
     &     kstar, rho, phi,xlum,xrstar,xteff,xcmass,xrc)
*
*      
      include 'params.h'
  
      REAL*8 :: X0(3,NMAX),X0DOT(3,NMAX), X(3,NMAX), XDOT(3,NMAX)
      REAL*8 :: f(3,NMAX),fdot(3,NMAX)
      REAL*8 :: T0(NMAX)
      REAL*8 :: MASS(NMAX)
      INTEGER :: KZ(50)
      INTEGER :: vec_len ! The length of the current vector passed thought the arguemnt
      INTEGER :: vec_offset
      INTEGER :: vec_idx(vec_len)
      INTEGER :: vec_id_list(NMAX) ! the generated list of updated IDs
      INTEGER :: vec_name_list(NMAX) ! the generated list of updated names
      INTEGER :: id_name(NMAX)  ! the corresponding name of the ID
      REAL*8 :: current_time
      REAL*8 :: toff
      REAL*8 :: time_step
      INTEGER :: total_n ! total number of particles
      INTEGER :: I, J ! Loop variables
*     Stellar evolution
      INTEGER :: kstar(NMAX) ! Stellar type indicator, see define.f
      REAL*8 :: phi(NMAX)    ! Local potential
      REAL*8 :: rho(NMAX)    ! Local density
      REAL*8 :: xlum(NMAX)   ! Luminosity
      REAL*8 :: xrstar(NMAX) ! Stellar radius 
      REAL*8 :: xteff(NMAX)  ! Effective temperature
      REAL*8 :: xcmass(NMAX) ! Core mass
      REAL*8 :: xrc(NMAX)    ! Core radius


*       Generate a list of updated particle IDs
          J = 1
          IF (KZ(46).EQ.4) THEN !output all particles
*       Need predict all particles, use xbpreadall (L.Wang)
            CALL xbpredall
*            
            DO 81 I = vec_offset, total_n ! traverse the whole particle set
                 vec_id_list(J) = I
                 vec_name_list(J) = id_name(I)
                 J = J + 1
  81       CONTINUE
          ELSE 
            DO 96 K = 1, vec_len
              I = vec_idx(K)
*      Suppress the if statement to select particle since we need to save
*      all active particles in the nxtlst (L.Wnag)
*              IF ((time_step.GT.0 .AND. T0(I).LT.current_time
*     &            .AND.T0(I).GE.(current_time-time_step)) .OR. 
*     &            (time_step.EQ.0 .AND. T0(I).EQ.current_time))  THEN
                vec_id_list(J) = I
                vec_name_list(J) = id_name(I)
                J = J + 1
*              ENDIF
  96       CONTINUE
          ENDIF

          IF(KZ(46).EQ.2) THEN
              DO 101 I = 1, J-1
                 WRITE (40,666) vec_name_list(I),current_time+toff,
     &             (X0(K,vec_id_list(I)),K=1,3),
     &             (X0DOT(K,vec_id_list(I)),K=1,3),
     &             MASS(vec_id_list(I)),
     &             (F(K,vec_id_list(I))*2,K=1,3),
     &             (FDOT(K,vec_id_list(I))*6,K=1,3),
     &             KSTAR(vec_id_list(I)),
     &             RHO(vec_id_list(I)),
     &             PHI(vec_id_list(I)),
     &             XLUM(vec_id_list(I)),
     &             XRSTAR(vec_id_list(I)),
     &             XTEFF(vec_id_list(I)),
     &             XCMASS(vec_id_list(I)),
     &             XRC(vec_id_list(I))
  101         CONTINUE
          ELSE
              DO 106 I = 1, J-1
                 WRITE (40,666) vec_name_list(I),current_time+toff,
     &             (X(K,vec_id_list(I)),K=1,3),
     &             (XDOT(K,vec_id_list(I)),K=1,3),
     &             MASS(vec_id_list(I)),
     &             (F(K,vec_id_list(I))*2,K=1,3),
     &             (FDOT(K,vec_id_list(I))*6,K=1,3),
     &             KSTAR(vec_id_list(I)),
     &             RHO(vec_id_list(I)),
     &             PHI(vec_id_list(I)),
     &             XLUM(vec_id_list(I)),
     &             XRSTAR(vec_id_list(I)),
     &             XTEFF(vec_id_list(I)),
     &             XCMASS(vec_id_list(I)),
     &             XRC(vec_id_list(I))
  106         CONTINUE

          ENDIF
  666    FORMAT (I9,",",ES14.7,",",
     &     ES14.7,",",ES14.7,",",ES14.7,",",
     &     ES14.7,",",ES14.7,",",ES14.7,",",
     &     ES14.7,",",
     &     ES14.7,",",ES14.7,",",ES14.7,",",
     &     ES14.7,",",ES14.7,",",ES14.7,",",
     &     I2,",",
     &     ES14.7,",",
     &     ES14.7,",",
     &     ES14.7,",",
     &     ES14.7,",",
     &     ES14.7,",",
     &     ES14.7,",",
     &     ES14.7)
*       Now J-1 is the actually size of vector to be written
      END



      SUBROUTINE HDF5_INIT
#ifdef H5OUTPUT
      USE HDF5
      INCLUDE 'params.h'
*       SAVE Block
      INTEGER :: h5_file_id, h5_step, h5_group_id, h5_dset_ids(32)
      INTEGER :: h5_vec_len
      REAL*8 :: h5_current_time
      LOGICAL h5_particle_updated(NMAX)
      COMMON/h5part/ h5_file_id, h5_step, h5_group_id, h5_dset_ids,
     &     h5_vec_len, h5_particle_updated
 
      INTEGER ERROR
*     Initialize FORTRAN interface.
      CALL h5open_f(ERROR) 

*     Create a new file using default properties.
      CALL h5fcreate_f('data.h5part', H5F_ACC_TRUNC_F, h5_file_id, 
     &      ERROR)

      h5_step = 0
      print *, 'h5_file_id_init', h5_file_id
      h5_group_id = 0
      h5_prev_group_id = 0
      h5_current_time = 0
*       set the update_particle marker array
      DO 350 I = 1, NMAX
         h5_particle_updated(I) = .FALSE.
  350 CONTINUE


#endif
      END

      
      SUBROUTINE HDF5_output(vec_idx,vec_len,vec_offset,X0,X0dot,X,XDOT,
     &     f, fdot,T0,MASS,KZ,id_name, current_time, time_step,
     &     toff, total_n,
     &     kstar, rho, phi,xlum,xrstar,xteff,xcmass,xrc,
     &     finalize)
#ifdef H5OUTPUT
      USE HDF5
      include 'params.h'
  
      REAL*8 :: X0(3,NMAX),X0DOT(3,NMAX), X(3,NMAX), XDOT(3,NMAX)
      REAL*8 :: f(3,NMAX),fdot(3,NMAX)
      REAL*8 :: T0(NMAX)
      REAL*8 :: MASS(NMAX)
      INTEGER :: KZ(50)
      INTEGER :: vec_len ! The length of the current vector passed thought the arguemnt
      INTEGER :: vec_offset
      INTEGER :: original_vec_len
      INTEGER :: vec_idx(vec_len)
      INTEGER :: vec_id_list(NMAX) ! the generated list of updated IDs
      INTEGER :: vec_name_list(NMAX) ! the generated list of updated names
      INTEGER :: id_name(NMAX)  ! the corresponding name of the ID
      LOGICAL :: finalize ! if TRUE, finalize the datasets
      INTEGER :: error
      REAL*8  :: vec_x(NMAX)
      REAL*8  :: vec_y(NMAX)
      REAL*8  :: vec_z(NMAX)
      INTEGER  :: vec_integer(NMAX)
      REAL*8 :: current_time
      REAL*8 :: toff
      REAL*8 :: time_step
      CHARACTER(LEN=16) :: h5_step_name
      CHARACTER(LEN=20) :: h5_step_group_name
      INTEGER :: total_n ! total number of particles
      INTEGER :: I, J ! Loop variables
*     Stellar evolution
      INTEGER :: kstar(NMAX) ! Stellar type indicator, see define.f
      REAL*8 :: rho(NMAX)    ! Local density
      REAL*8 :: phi(NMAX)    ! Local potential
      REAL*8 :: xlum(NMAX)   ! Luminosity
      REAL*8 :: xrstar(NMAX) ! Stellar radius 
      REAL*8 :: xteff(NMAX)  ! Effective temperature
      REAL*8 :: xcmass(NMAX) ! Core mass
      REAL*8 :: xrc(NMAX)    ! Core radius
*     SAVE Block
      INTEGER :: h5_file_id, h5_step, h5_group_id, h5_dset_ids(32)
      INTEGER :: h5_vec_len
      LOGICAL :: h5_particle_updated(NMAX)
      LOGICAL :: h5_file_inited
      COMMON/h5part/ h5_file_id, h5_step, h5_group_id, h5_dset_ids,
     &     h5_vec_len, h5_particle_updated
      DATA h5_file_inited /.false./

*       CALL INIT if not yet done.
      IF (.not.h5_file_inited) then
          CALL HDF5_INIT
          h5_file_inited = .true.
      END IF 
*       IF file cannot be initilized, quit the subroutine
      IF (h5_file_id .EQ. 0) RETURN


*      IF(current_time+toff .GT. h5_current_time) THEN
          DO 161 I = 1, 32 ! Close any previously opened dset
             IF(h5_dset_ids(I) .GT. 0) THEN
*                 CALL h5dclose_f(h5_dset_ids(I), error)
                 h5_dset_ids(I) = 0
             ENDIF
  161     CONTINUE
*     Prepare new group
          WRITE (h5_step_name, *), h5_step
          h5_step_group_name = 'Step#' // ADJUSTL(h5_step_name)
          CALL h5gcreate_f(h5_file_id, h5_step_group_name,
     &       h5_group_id, error) ! Create the group
          h5_step = h5_step + 1

          CALL HDF5_write_attribute_scalar_real(h5_group_id, 'Time',
     &        current_time+toff) ! Write attribute to the group

          CALL HDF5_write_attribute_scalar_integer(h5_group_id,'TotalN',
     &        total_n) ! Write attribute to the group

*================New Group Prep. done===================
*       Generate a list of updated particle IDs
          J = 1
          IF (KZ(46).EQ.3) THEN !output all particles
*       Need predict all particles, use xbpreadall (L.Wang)
             CALL xbpredall
*
            DO 181 I = vec_offset, total_n ! traverse the whole particle set
                 vec_id_list(J) = I
                 vec_name_list(J) = id_name(I)
                 J = J + 1
  181       CONTINUE
          ELSE 
            DO 196 K = 1, vec_len
              I = vec_idx(K)
*      Suppress the if statement to select particle since we need to save
*      all active particles in the nxtlst (L.Wnag)
*              IF ((time_step.GT.0 .AND. T0(I).LT.current_time
*     &            .AND.T0(I).GE.(current_time-time_step)) .OR. 
*     &            (time_step.EQ.0 .AND. T0(I).EQ.current_time))  THEN
                vec_id_list(J) = I
                vec_name_list(J) = id_name(I)
                J = J + 1
*              ENDIF
  196       CONTINUE
          ENDIF
*       Now J-1 is the actually size of vector to be written
*          vec_len = J-1
          original_vec_len = 0
          CALL HDF5_write_integer_vector_as_dset(h5_group_id, 'ID',
     &        vec_name_list, J-1, 1, h5_dset_ids(1),
     &        original_vec_len,finalize)

*       Write X, Y, Z according to the generated ID list
          DO 221 I = 1, J-1
             IF (KZ(46).EQ.1) THEN
               vec_x(i) = X0(1, vec_id_list(I))
               vec_y(i) = X0(2, vec_id_list(I))
               vec_z(i) = X0(3, vec_id_list(I))
             ELSE 
               vec_x(i) = X(1, vec_id_list(I))
               vec_y(i) = X(2, vec_id_list(I))
               vec_z(i) = X(3, vec_id_list(I))
             ENDIF
  221     CONTINUE
          CALL HDF5_write_real_vector_as_dset(h5_group_id, 'X',
     &        vec_x, J-1, 1, h5_dset_ids(2),
     &        original_vec_len,finalize)
          CALL HDF5_write_real_vector_as_dset(h5_group_id, 'Y',
     &        vec_y, J-1, 1, h5_dset_ids(3),
     &        original_vec_len,finalize)
          CALL HDF5_write_real_vector_as_dset(h5_group_id, 'Z',
     &        vec_z, J-1, 1, h5_dset_ids(4),
     &        original_vec_len,finalize)
*       Write VX, VY, VZ according to the generated ID list
          DO 231 I = 1, J-1
             IF (KZ(46).EQ.1) THEN
               vec_x(i) = X0DOT(1, vec_id_list(I))
               vec_y(i) = X0DOT(2, vec_id_list(I))
               vec_z(i) = X0DOT(3, vec_id_list(I))
             ELSE 
               vec_x(i) = XDOT(1, vec_id_list(I))
               vec_y(i) = XDOT(2, vec_id_list(I))
               vec_z(i) = XDOT(3, vec_id_list(I))
             ENDIF
  231     CONTINUE
          CALL HDF5_write_real_vector_as_dset(h5_group_id, 'VX',
     &        vec_x, J-1, 1, h5_dset_ids(5),
     &        original_vec_len,finalize)
          CALL HDF5_write_real_vector_as_dset(h5_group_id, 'VY',
     &        vec_y, J-1, 1, h5_dset_ids(6),
     &        original_vec_len,finalize)
          CALL HDF5_write_real_vector_as_dset(h5_group_id, 'VZ',
     &        vec_z, J-1, 1, h5_dset_ids(7),
     &        original_vec_len,finalize)
*       Write individual mass for each particle (needed when stellar evolition is on)
         DO 241 I = 1, J-1
             vec_x(i) = MASS(vec_id_list(I))
  241    CONTINUE
          CALL HDF5_write_real_vector_as_dset(h5_group_id, 'Mass',
     &        vec_x, J-1, 1, h5_dset_ids(8),
     &        original_vec_len,finalize)
*       Write AX, AY, AZ according to the generated ID list
          DO 251 I = 1, J-1
             vec_x(i) = F(1, vec_id_list(I))*2
             vec_y(i) = F(2, vec_id_list(I))*2
             vec_z(i) = F(3, vec_id_list(I))*2
  251     CONTINUE
          CALL HDF5_write_real_vector_as_dset(h5_group_id, 'AX',
     &        vec_x, J-1, 1, h5_dset_ids(9),
     &        original_vec_len,finalize)
          CALL HDF5_write_real_vector_as_dset(h5_group_id, 'AY',
     &        vec_y, J-1, 1, h5_dset_ids(10),
     &        original_vec_len,finalize)
          CALL HDF5_write_real_vector_as_dset(h5_group_id, 'AZ',
     &        vec_z, J-1, 1, h5_dset_ids(11),
     &        original_vec_len,finalize)
*       Write JX, JY, JZ according to the generated ID list
          DO 261 I = 1, J-1
             vec_x(i) = FDOT(1, vec_id_list(I))*6
             vec_y(i) = FDOT(2, vec_id_list(I))*6
             vec_z(i) = FDOT(3, vec_id_list(I))*6
  261     CONTINUE
          CALL HDF5_write_real_vector_as_dset(h5_group_id, 'JX',
     &        vec_x, J-1, 1, h5_dset_ids(12),
     &        original_vec_len,finalize)
          CALL HDF5_write_real_vector_as_dset(h5_group_id, 'JY',
     &        vec_y, J-1, 1, h5_dset_ids(13),
     &        original_vec_len,finalize)
          CALL HDF5_write_real_vector_as_dset(h5_group_id, 'JZ',
     &        vec_z, J-1, 1, h5_dset_ids(14),
     &        original_vec_len,finalize)
*     Write stellar evolution data if KZ(12)=1 and KZ(19)=1 or 3
          IF(KZ(12).GE.1 .AND. (KZ(19).EQ.1 .OR. KZ(19).GE.3)) THEN
*     Write kstar
            DO 308 I = 1, J-1
               vec_integer(i) = kstar(vec_id_list(I))
  308       CONTINUE
            CALL HDF5_write_integer_vector_as_dset(h5_group_id, 'KSTAR',
     &          vec_integer, J-1, 1, h5_dset_ids(15),
     &          original_vec_len,finalize)
*     Write rho
            DO 318 I = 1, J-1
               vec_x(i) = rho(vec_id_list(I))
  318       CONTINUE
            CALL HDF5_write_real_vector_as_dset(h5_group_id, 'RHO',
     &          vec_x, J-1, 1, h5_dset_ids(16),
     &          original_vec_len,finalize)
*     Write phi
            DO 328 I = 1, J-1
               vec_x(i) = phi(vec_id_list(I))
  328       CONTINUE
            CALL HDF5_write_real_vector_as_dset(h5_group_id, 'PHI',
     &          vec_x, J-1, 1, h5_dset_ids(17),
     &          original_vec_len,finalize)
*     Write xlum
            DO 338 I = 1, J-1
               vec_x(i) = xlum(vec_id_list(I))
  338       CONTINUE
            CALL HDF5_write_real_vector_as_dset(h5_group_id, 'LUM',
     &          vec_x, J-1, 1, h5_dset_ids(18),
     &          original_vec_len,finalize)
*     Write xrstar
            DO 348 I = 1, J-1
               vec_x(i) = xrstar(vec_id_list(I))
  348       CONTINUE
            CALL HDF5_write_real_vector_as_dset(h5_group_id, 'RSTAR',
     &          vec_x, J-1, 1, h5_dset_ids(19),
     &          original_vec_len,finalize)
*     Write xteff
            DO 358 I = 1, J-1
               vec_x(i) = xteff(vec_id_list(I))
  358       CONTINUE
            CALL HDF5_write_real_vector_as_dset(h5_group_id, 'TEFF',
     &          vec_x, J-1, 1, h5_dset_ids(20),
     &          original_vec_len,finalize)
*     Write xzmet
*     This value comes simply from the input file, same everywhere, no need write
*            DO 368 I = 1, J-1
*               vec_x(i) = xzmet(vec_id_list(I))
*  368       CONTINUE
*            CALL HDF5_write_real_vector_as_dset(h5_group_id, 'ZMET',
*     &          vec_x, J-1, 1, h5_dset_ids(21),
*     &          original_vec_len,finalize)
*     Write xcmass
            DO 378 I = 1, J-1
               vec_x(i) = xcmass(vec_id_list(I))
  378       CONTINUE
            CALL HDF5_write_real_vector_as_dset(h5_group_id, 'CMASS',
     &          vec_x, J-1, 1, h5_dset_ids(22),
     &          original_vec_len,finalize)
*     Write xrc
            DO 388 I = 1, J-1
               vec_x(i) = xrc(vec_id_list(I))
  388       CONTINUE
            CALL HDF5_write_real_vector_as_dset(h5_group_id, 'CRADIUS',
     &          vec_x, J-1, 1, h5_dset_ids(23),
     &          original_vec_len,finalize)


            GOTO 555
  555     CONTINUE

            ENDIF !if KZ(12)=1 and KZ(19)=1 or 3 

          CALL h5fflush_f(h5_group_id, H5F_SCOPE_LOCAL_F, error)
          CALL h5gclose_f(h5_group_id, error) 
*==================Write back done=====================
*     ! ELSE (current_time <= h5_current_time) THEN

*      ENDIF ! current_time .GT. h5_current_time
#endif
      END




      SUBROUTINE HDF5_write_real_vector_as_dset(group_id, dset_name,
     &     vec, vec_len, offset, dset_id, original_dset_len, finalize)
#ifdef H5OUTPUT
      USE HDF5
      IMPLICIT NONE
      INTEGER :: group_id ! the group ID to be written upon
      CHARACTER(LEN=*), INTENT(IN) :: dset_name ! The name of the dset
      INTEGER :: vec_len ! the data array
      REAL*8  :: vec(vec_len)     ! the data array
      REAL  :: vec_written(vec_len)   ! the data actually written
      INTEGER :: offset  ! offset of vec (by default, start from 1)
      INTEGER :: dset_id ! IF not 0, write to that dset; otherwise create new
      INTEGER :: original_dset_len ! The original dset length, if expanding
      LOGICAL :: finalize ! Finalize the dataset, no more data can be added to it
      INTEGER :: I ! Loop variable
      INTEGER :: error
      INTEGER :: dspace_id
      INTEGER :: memspace_id
      INTEGER :: crp_list ! Dataset creation property identifier
      INTEGER(8), DIMENSION(2) :: data_dims
      INTEGER(8), DIMENSION(2) :: data_maxdims
      INTEGER(8), DIMENSION(2) :: data_chunkdims
      INTEGER(8), DIMENSION(2) :: data_start
      INTEGER(8), DIMENSION(2) :: data_count

      data_dims(1) = vec_len
      data_dims(2) = 1
      data_maxdims(1) = H5S_UNLIMITED_F
      data_maxdims(2) = 1
      data_chunkdims(1) = 10
      data_chunkdims(2) = 1


      IF(dset_id .EQ. 0) THEN
*     Create new dataset
         IF(finalize .EQV. .TRUE.) THEN 
            CALL h5screate_simple_f(1, data_dims, dspace_id,error)
            CALL h5dcreate_f(group_id, dset_name, H5T_NATIVE_REAL,
     &         dspace_id, dset_id, error)
         ELSE ! Finalize = .FALSE.
*       Create simple dataspace with 1D extensible dimension
            CALL h5screate_simple_f(1, data_dims, dspace_id,
     &         error, data_maxdims)
*       Modify dataset creation properties, i.e. enable chunking
            CALL h5pcreate_f(H5P_DATASET_CREATE_F, crp_list, error)
            CALL h5pset_chunk_f(crp_list, 1, data_chunkdims, error)
            CALL h5dcreate_f(group_id, dset_name, H5T_NATIVE_REAL,
     &         dspace_id, dset_id, error, crp_list)
         ENDIF ! finalize == .TRUE.

         IF(offset .GT. 1) THEN
*     Shift the data leftward & reduce precision
            DO 10 I = 1, vec_len
               vec_written(I) = REAL(vec(I + offset - 1))
  10        CONTINUE
            CALL h5dwrite_f(dset_id, H5T_NATIVE_REAL, vec_written,
     &         data_dims, error)
         ELSE ! offset = 1
            DO 15 I = 1, vec_len ! only reduce precision
               vec_written(I) = REAL(vec(I))
  15        CONTINUE
            CALL h5dwrite_f(dset_id, H5T_NATIVE_REAL, vec_written,
     &         data_dims, error)
         ENDIF ! offset > 1
         IF(finalize .EQV. .TRUE.) CALL h5dclose_f(dset_id, error)
      ELSE ! dset_id > 0
*     Expand the existing dataset
          data_start(1) = original_dset_len ! Orignal dset length (offset)
          data_start(2) = 1
          data_count(1) = vec_len ! length of the new data
          data_count(2) = 1
          data_dims(1) = vec_len + original_dset_len
          data_dims(2) = 1

          CALL h5dset_extent_f(dset_id, data_dims, error)
*       Create memspace to indicate the size of the buffer, i.e. data_count
          CALL h5screate_simple_f(1, data_count, memspace_id, error)
          CALL h5dget_space_f(dset_id, dspace_id, error)
          CALL h5sselect_hyperslab_f(dspace_id, H5S_SELECT_SET_F,
     &       data_start, data_count, error)

         IF(offset .GT. 1) THEN
*     Shift the data leftward & reduce precision
          DO 20 I = 1, vec_len
               vec_written(I) = REAL(vec(I + offset - 1))
  20      CONTINUE
          CALL h5dwrite_f(dset_id, H5T_NATIVE_REAL, vec_written,
     &        data_dims, error, memspace_id, dspace_id)
         ELSE ! offset = 1
          DO 25 I = 1, vec_len ! only reduce precision
               vec_written(I) = REAL(vec(I))
  25      CONTINUE
           CALL h5dwrite_f(dset_id, H5T_NATIVE_REAL, vec_written,
     &         data_dims, error, memspace_id, dspace_id)
         ENDIF ! offset > 1
*         CALL h5sclose_f(dspace_id, error)
         IF(finalize .EQV. .TRUE.) CALL h5dclose_f(dset_id, error)
      ENDIF ! dset_id == 0
#endif
      END


      SUBROUTINE HDF5_write_integer_vector_as_dset(group_id, dset_name,
     &     vec, vec_len, offset, dset_id, original_dset_len,finalize)
#ifdef H5OUTPUT
      USE HDF5
      INTEGER :: group_id ! the group ID to be written upon
      CHARACTER(LEN=*), INTENT(IN) :: dset_name ! The name of the dset
      INTEGER :: vec_len ! the data array
      INTEGER :: offset  ! offset of vec (by default, start from 1)
      INTEGER :: vec(vec_len+offset-1)     ! the data array
      INTEGER :: dset_id ! IF not 0, write to that dset; otherwise create new
      INTEGER :: original_dset_len ! The original dset length, before expanding
      LOGICAL :: finalize ! Finalize the dataset, no more data can be added to it
      INTEGER  :: vec_written(vec_len+offset-1)   !if vec starts not from 1, use this
      INTEGER :: I ! Loop variable
      INTEGER :: error
      INTEGER :: dspace_id
      INTEGER :: memspace_id
      INTEGER :: crp_list ! Dataset creation property identifier
      INTEGER(8), DIMENSION(2) :: data_dims
      INTEGER(8), DIMENSION(2) :: data_maxdims
      INTEGER(8), DIMENSION(2) :: data_chunkdims
      INTEGER(8), DIMENSION(2) :: data_start
      INTEGER(8), DIMENSION(2) :: data_count

      data_dims(1) = vec_len
      data_dims(2) = 1
      data_maxdims(1) = H5S_UNLIMITED_F
      data_maxdims(2) = 1
      data_chunkdims(1) = 100
      data_chunkdims(2) = 1


      IF(dset_id .EQ. 0) THEN
*     Create new dataset
         IF(finalize .EQV. .TRUE.) THEN 
            CALL h5screate_simple_f(1, data_dims, dspace_id,error)
            CALL h5dcreate_f(group_id, dset_name, H5T_NATIVE_INTEGER,
     &         dspace_id, dset_id, error)
         ELSE 
*       Create simple dataspace with 1D extensible dimension
            CALL h5screate_simple_f(1, data_dims, dspace_id,
     &         error, data_maxdims)
*       Modify dataset creation properties, i.e. enable chunking
            CALL h5pcreate_f(H5P_DATASET_CREATE_F, crp_list, error)
            CALL h5pset_chunk_f(crp_list, 1, data_chunkdims, error)
            CALL h5dcreate_f(group_id, dset_name, H5T_NATIVE_INTEGER,
     &         dspace_id, dset_id, error, crp_list)
         ENDIF


         IF(offset .GT. 1) THEN
*     Shift the data leftward
            DO 30 I = 1, vec_len
               vec_written(I) = vec(I + offset - 1)
  30        CONTINUE
            CALL h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, vec_written,
     &         data_dims, error)
         ELSE
            CALL h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, vec,
     &         data_dims, error)
         ENDIF
         IF(finalize .EQV. .TRUE.) CALL h5dclose_f(dset_id, error)
      ELSE  ! dset_id > 0
*     Expand the existing dataset
          data_start(1) = original_dset_len ! Orignal dset length (offset)
          data_start(2) = 1
          data_count(1) = vec_len ! length of the new data
          data_count(2) = 1
          data_dims(1) = vec_len + original_dset_len
          data_dims(2) = 1

          CALL h5dset_extent_f(dset_id, data_dims, error)
*       Create memspace to indicate the size of the buffer, i.e. data_count
          CALL h5screate_simple_f(2, data_count, memspace_id, error)
          CALL h5dget_space_f(dset_id, dspace_id, error)
          CALL h5sselect_hyperslab_f(dspace_id, H5S_SELECT_SET_F,
     &       data_start, data_count, error)

         IF(offset .GT. 1) THEN
*     Shift the data leftward
          DO 40 I = 1, vec_len
               vec_written(I) = vec(I + offset - 1)
  40      CONTINUE
          CALL h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, vec_written,
     &        data_dims, error, memspace_id, dspace_id)
         ELSE ! offset = 1
           CALL h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, vec,
     &         data_dims, error, memspace_id, dspace_id)
         ENDIF ! offset > 1
         CALL h5sclose_f(dspace_id, error)
         IF(finalize .EQV. .TRUE.) CALL h5dclose_f(dset_id, error)
      ENDIF
#endif
      END

      SUBROUTINE HDF5_write_attribute_scalar_real(loc_id, att_name, val)
#ifdef H5OUTPUT
      USE HDF5
      IMPLICIT NONE
      CHARACTER(LEN=*) :: att_name
      REAL*8 :: val
      INTEGER :: loc_id
      INTEGER :: attrib_space_id
      INTEGER :: attrib_id
      INTEGER :: error

      INTEGER(8), DIMENSION(2) :: data_dims
      data_dims(1) = 1
      data_dims(2) = 1

*       Write attributes to the group
      CALL H5Screate_f(H5S_SCALAR_F, attrib_space_id, error)
      CALL H5Acreate_f(loc_id, TRIM(att_name), H5T_NATIVE_DOUBLE,
     &     attrib_space_id, attrib_id, error, H5P_DEFAULT_F,
     &     H5P_DEFAULT_F)
      CALL H5Awrite_f(attrib_id, H5T_NATIVE_DOUBLE, val, data_dims,
     &      error)
      CALL H5Aclose_f(attrib_id, error)
      CALL H5Sclose_f(attrib_space_id, error)

#endif
      END


      SUBROUTINE HDF5_write_attribute_scalar_integer(loc_id, 
     &                att_name, val)
#ifdef H5OUTPUT
      USE HDF5
      IMPLICIT NONE
      CHARACTER(LEN=*) :: att_name
      INTEGER :: val
      INTEGER :: loc_id
      INTEGER :: attrib_space_id
      INTEGER :: attrib_id
      INTEGER :: error

      INTEGER(8), DIMENSION(2) :: data_dims
      data_dims(1) = 1
      data_dims(2) = 1

*       Write attributes to the group
      CALL H5Screate_f(H5S_SCALAR_F, attrib_space_id, error)
      CALL H5Acreate_f(loc_id, TRIM(att_name), H5T_NATIVE_INTEGER,
     &     attrib_space_id, attrib_id, error, H5P_DEFAULT_F,
     &     H5P_DEFAULT_F)
      CALL H5Awrite_f(attrib_id, H5T_NATIVE_INTEGER, val, data_dims,
     &      error)
      CALL H5Aclose_f(attrib_id, error)
      CALL H5Sclose_f(attrib_space_id, error)

#endif
      END
