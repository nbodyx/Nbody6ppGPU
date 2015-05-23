*     Developed by Maxwell Xu CAI (NAOC/KIAA)
*     For bug report/feedback, please email maxwellemail@gmail.com
*     The binary format is developed by Long Wang (May 11, 2015)
*
*     Version 2.0 -- Bug fixes & finetuning, June 12, 2014
*     Version 1.4 -- June 11, 2014
*     Version 1.0 -- Initial release, Nov 2013
*

      SUBROUTINE HDF5_INIT(h5_fname)
*
*
*       Initialize a new HDF5 file given by its file name.
*       --------------------------
*
*
#ifdef H5OUTPUT
      USE HDF5
      INCLUDE 'params.h'
*       SAVE Block
      INTEGER :: h5_file_id, h5_step, h5_group_id, h5_dset_ids(256)
      INTEGER :: h5_vec_len
      REAL*8 :: h5_current_time
      CHARACTER(LEN=64) :: h5_file_name
      CHARACTER(LEN=64) :: h5_fname
      COMMON/h5part/ h5_file_id, h5_step, h5_group_id, h5_dset_ids,
     &     h5_vec_len, h5_file_name, h5_current_time
      INTEGER ERROR
      
*       Close any previously opened file (if any)
      IF (h5_file_id .GT. 0) CALL HDF5_CLOSE

*     Initialize FORTRAN interface.
      CALL h5open_f(ERROR) 

*     Create a new file using default properties.
      h5_file_name = h5_fname
      CALL h5fcreate_f(TRIM(h5_file_name), H5F_ACC_TRUNC_F, h5_file_id, 
     &      ERROR)

      h5_step = -1
      print*, 'HDF5 Init ID:',h5_file_id,' Snapshot file: ',h5_file_name
      h5_group_id = 0
      h5_prev_group_id = 0
      h5_current_time = -999.0


#endif
      END

      SUBROUTINE HDF5_CLOSE
*
*
*       Close the currently opened HDF5 handle.
*       No argument needed. This subroutine share the common block with
*       HDF5_INIT(). If the current HDF5 indicated by the h5_file_id is
*       closed, this subroutine will do nothing.
*       --------------------------
*
*
#ifdef H5OUTPUT
      USE HDF5
      INCLUDE 'params.h'
*       SAVE Block
      INTEGER :: h5_file_id, h5_step, h5_group_id, h5_dset_ids(256)
      INTEGER :: h5_vec_len
      REAL*8 :: h5_current_time
      CHARACTER(LEN=64) :: h5_file_name
      COMMON/h5part/ h5_file_id, h5_step, h5_group_id, h5_dset_ids,
     &     h5_vec_len, h5_file_name, h5_current_time
      INTEGER ERROR

      IF (h5_file_id .GT. 0) THEN
        CALL h5fclose_f(h5_file_id, ERROR)
        h5_group_id = 0
        DO 5 I = 1, 256
            h5_dset_ids(I) = 0
  5    CONTINUE
        h5_step = -1
        h5_group_id = 0
        h5_prev_group_id = 0
        h5_current_time = 0
      ENDIF
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
      REAL*4  :: vec(vec_len)     ! the data array
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
               vec_written(I) = vec(I + offset - 1)
  10        CONTINUE
            CALL h5dwrite_f(dset_id, H5T_NATIVE_REAL, vec_written,
     &         data_dims, error)
         ELSE ! offset = 1
c$$$            DO 15 I = 1, vec_len ! only reduce precision
c$$$               vec_written(I) = REAL(vec(I))
c$$$  15        CONTINUE
            CALL h5dwrite_f(dset_id, H5T_NATIVE_REAL, vec,
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
               vec_written(I) = vec(I + offset - 1)
  20      CONTINUE
          CALL h5dwrite_f(dset_id, H5T_NATIVE_REAL, vec_written,
     &        data_dims, error, memspace_id, dspace_id)
         ELSE ! offset = 1
c$$$          DO 25 I = 1, vec_len ! only reduce precision
c$$$               vec_written(I) = REAL(vec(I))
c$$$  25      CONTINUE
           CALL h5dwrite_f(dset_id, H5T_NATIVE_REAL, vec,
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


      SUBROUTINE custom_update_file(TTOT,DELTAT)
#ifdef H5OUTPUT
      REAL*8 :: TTOT, DELTAT
      CHARACTER(LEN=20) :: TCHAR
      CHARACTER(LEN=64) :: h5_file_name
      COMMON/h5part/ h5_file_id, h5_step, h5_group_id, h5_dset_ids,
     &     h5_vec_len, h5_file_name, h5_current_time
      INTEGER :: h5_file_id, h5_step, h5_group_id, h5_dset_ids(256)
      INTEGER :: h5_vec_len
      REAL*8 :: h5_current_time

*     Close any previously opened file (if any)
      IF (h5_file_id .GT. 0) CALL HDF5_CLOSE
      h5_file_id = 0

      call string_left(TCHAR,TTOT,DELTAT)
      h5_file_name='snap.40_'//trim(TCHAR)//'.h5part'
#else
      COMMON/BINARYOUT/ DTOUT
      REAL*8 DELTAT,DTOUT
      DTOUT = DELTAT
#endif
      
      END


      SUBROUTINE output_single(TTOT,N_SINGLE,KSEV,KMODE)
#ifdef H5OUTPUT
      USE HDF5
      INCLUDE 'params.h'
      INTEGER :: h5_file_id, h5_step, h5_group_id, h5_dset_ids(256)
      INTEGER :: h5_vec_len
      REAL*8 :: h5_current_time
      REAL*8 :: TTOT
      INTEGER :: N_SINGLE,KSEV,KMODE
      CHARACTER(LEN=64) :: h5_file_name
      CHARACTER(LEN=16) :: h5_step_name
      CHARACTER(LEN=20) :: h5_step_group_name
      COMMON/h5part/ h5_file_id, h5_step, h5_group_id, h5_dset_ids,
     &     h5_vec_len, h5_file_name, h5_current_time
      INTEGER ERROR
      INTEGER :: original_vec_len
      LOGICAL :: finalize

*     single
      COMMON/OUTSINGLE/ S_M(NMAX),S_X1(NMAX), S_X2(NMAX), S_X3(NMAX),
     &                  S_V1(NMAX), S_V2(NMAX), S_V3(NMAX),
C     &                  S_F1(NMAX), S_F2(NMAX), S_F3(NMAX),
     &                  S_RS(NMAX), S_L(NMAX), S_TE(NMAX),
     &                  S_RC(NMAX), S_MC(NMAX),
     &                  NS_KW(NMAX),NS_NAM(NMAX)
      REAL*4 S_M, S_X1, S_X2, S_X3, S_V1, S_V2, S_V3
C      REAL*4 S_F1, S_F2, S_F3
      REAL*4 S_RS, S_L, S_TE, S_RC, S_MC
      INTEGER NS_KW,NS_NAM

*       CALL INIT if not yet done.
      IF (h5_file_id.eq.0) then
          CALL HDF5_INIT(h5_file_name)
      END IF 
*       IF file cannot be initilized, quit the subroutine
      IF (h5_file_id .EQ. 0) RETURN

      
      IF (TTOT.GT.h5_current_time) THEN
*       Needs to create new group
          h5_step = h5_step + 1
          WRITE (h5_step_name, *), h5_step
          h5_step_group_name = 'Step#' // ADJUSTL(h5_step_name)
          CALL h5gcreate_f(h5_file_id, h5_step_group_name,
     &       h5_group_id, error) ! Create the group
*       Close any previously opened datasets
          DO 161 I = 1, 256 ! Close any previously opened dset
            IF(h5_dset_ids(I) .GT. 0) THEN
*               CALL h5dclose_f(h5_dset_ids(I), error)
                h5_dset_ids(I) = 0
            ENDIF
  161     CONTINUE

*       Write group attributes (more to follow?)
          CALL HDF5_write_attribute_scalar_real(h5_group_id, 'Time',
     &        TTOT) ! Write attribute to the group

          CALL HDF5_write_attribute_scalar_integer(h5_group_id,
     &    'N_SINGLE', N_SINGLE) ! Write attribute to the group
          
          h5_current_time = TTOT
      ELSE
*       Open existing group
          WRITE (h5_step_name, *), h5_step
          h5_step_group_name = '/Step#' // ADJUSTL(h5_step_name)
          CALL h5gopen_f(h5_file_id, h5_step_group_name, h5_group_id,
     &        ERROR)
      ENDIF 
      original_vec_len = 0
      finalize = .TRUE.
      IF (N_SINGLE.LE.0) RETURN

*     Write datasets
      CALL HDF5_write_real_vector_as_dset(h5_group_id, 'X1',
     &     S_X1, N_SINGLE, 1, h5_dset_ids(2),
     &     original_vec_len,finalize)
      CALL HDF5_write_real_vector_as_dset(h5_group_id, 'X2',
     &     S_X2, N_SINGLE, 1, h5_dset_ids(3),
     &     original_vec_len,finalize)
      CALL HDF5_write_real_vector_as_dset(h5_group_id, 'X3',
     &     S_X3, N_SINGLE, 1, h5_dset_ids(4),
     &     original_vec_len,finalize)

      CALL HDF5_write_real_vector_as_dset(h5_group_id, 'V1',
     &     S_V1, N_SINGLE, 1, h5_dset_ids(5),
     &     original_vec_len,finalize)
      CALL HDF5_write_real_vector_as_dset(h5_group_id, 'V2',
     &     S_V2, N_SINGLE, 1, h5_dset_ids(6),
     &     original_vec_len,finalize)
      CALL HDF5_write_real_vector_as_dset(h5_group_id, 'V3',
     &     S_V3, N_SINGLE, 1, h5_dset_ids(7),
     &     original_vec_len,finalize)

      CALL HDF5_write_real_vector_as_dset(h5_group_id, 'M',
     &     S_M, N_SINGLE, 1, h5_dset_ids(8),
     &     original_vec_len,finalize)

*     stellar evolution data
      IF(KSEV.GT.0) THEN
         CALL HDF5_write_real_vector_as_dset(h5_group_id, 'RS',
     &        S_RS, N_SINGLE, 1, h5_dset_ids(9),
     &        original_vec_len,finalize)
         CALL HDF5_write_real_vector_as_dset(h5_group_id, 'L',
     &        S_L, N_SINGLE, 1, h5_dset_ids(10),
     &        original_vec_len,finalize)
         CALL HDF5_write_real_vector_as_dset(h5_group_id, 'TE',
     &        S_TE, N_SINGLE, 1, h5_dset_ids(11),
     &        original_vec_len,finalize)
         CALL HDF5_write_real_vector_as_dset(h5_group_id, 'RC',
     &        S_RC, N_SINGLE, 1, h5_dset_ids(12),
     &        original_vec_len,finalize)
         CALL HDF5_write_real_vector_as_dset(h5_group_id, 'MC',
     &        S_MC, N_SINGLE, 1, h5_dset_ids(13),
     &        original_vec_len,finalize)
         CALL HDF5_write_integer_vector_as_dset(h5_group_id, 'KW',
     &        NS_KW, N_SINGLE, 1, h5_dset_ids(14),
     &        original_vec_len,finalize)
      END IF

      CALL HDF5_write_integer_vector_as_dset(h5_group_id, 'NAM',
     &     NS_NAM, N_SINGLE, 1, h5_dset_ids(15),
     &     original_vec_len,finalize)

      CALL h5fflush_f(h5_group_id, H5F_SCOPE_LOCAL_F, error)
      CALL h5gclose_f(h5_group_id, error) 
#else
      INCLUDE 'params.h'
*     single
      COMMON/OUTSINGLE/ S_M(NMAX),S_X1(NMAX), S_X2(NMAX), S_X3(NMAX),
     &                  S_V1(NMAX), S_V2(NMAX), S_V3(NMAX),
C     &                  S_F1(NMAX), S_F2(NMAX), S_F3(NMAX),
     &                  S_RS(NMAX), S_L(NMAX), S_TE(NMAX),
     &                  S_RC(NMAX), S_MC(NMAX),
     &                  NS_KW(NMAX),NS_NAM(NMAX)
      COMMON/BINARYOUT/ DTOUT
      REAL*8 DTOUT,TTOT
      REAL*4 S_M, S_X1, S_X2, S_X3, S_V1, S_V2, S_V3
C      REAL*4 S_F1, S_F2, S_F3
      REAL*4 S_RS, S_L, S_TE, S_RC, S_MC
      INTEGER NS_KW,NS_NAM,N_SINGLE,KSEV,KMODE
      CHARACTER(LEN=20) :: TCHAR
      CHARACTER(LEN=64) :: filename
      call string_left(TCHAR,TTOT,DTOUT)
      filename='single.40_'//trim(TCHAR)
      IF(KMODE.EQ.1.OR.KMODE.EQ.3) THEN
       OPEN (UNIT=40,STATUS='UNKNOWN',FORM='UNFORMATTED',FILE=filename)
       WRITE (40) N_SINGLE, NS_NAM(1:N_SINGLE), S_M(1:N_SINGLE),
     &      S_X1(1:N_SINGLE), S_X2(1:N_SINGLE), S_X3(1:N_SINGLE),
     &      S_V1(1:N_SINGLE), S_V2(1:N_SINGLE), S_V3(1:N_SINGLE)
       IF(KSEV.GT.0) THEN
         WRITE (40) S_RS(1:N_SINGLE), S_L(1:N_SINGLE), S_TE(1:N_SINGLE),
     &         S_RC(1:N_SINGLE), S_MC(1:N_SINGLE), NS_KW(1:N_SINGLE)
       END IF
      ELSE IF(KMODE.EQ.2.OR.KMODE.EQ.4) THEN
       OPEN (UNIT=40,STATUS='UNKNOWN',FORM='FORMATTED',FILE=filename)
       WRITE(40,*) '## N_SINGLE ',N_SINGLE,' Time[NB]',TTOT
       IF(KSEV.GT.0) THEN
          DO K=1,N_SINGLE
             WRITE (40,*) NS_NAM(K), S_M(K),
     &            S_X1(K), S_X2(K), S_X3(K),
     &            S_V1(K), S_V2(K), S_V3(K),
     &            S_RS(K), S_L(K), S_TE(K),
     &            S_MC(K), S_RC(K), NS_KW(K)
          END DO
       ELSE
          DO K=1,N_SINGLE
             WRITE (40,*) NS_NAM(K), S_M(K),
     &            S_X1(K), S_X2(K), S_X3(K),
     &            S_V1(K), S_V2(K), S_V3(K)
          END DO
       END IF
      END IF
      CLOSE(40)
#endif
      END





      SUBROUTINE output_binary(TTOT,N_BINARY,KSEV,KMODE)
#ifdef H5OUTPUT
      USE HDF5
      INCLUDE 'params.h'
      INTEGER :: h5_file_id, h5_step, h5_group_id, h5_dset_ids(256)
      INTEGER :: h5_subgroup_id ! ID for the Step#i groups
      INTEGER :: h5_vec_len
      REAL*8 :: h5_current_time,TTOT
      INTEGER :: N_BINARY,KSEV,KMODE
      CHARACTER(LEN=64) :: h5_file_name
      CHARACTER(LEN=64) :: h5_step_name
      CHARACTER(LEN=20) :: h5_step_group_name
      COMMON/h5part/ h5_file_id, h5_step, h5_group_id, h5_dset_ids,
     &     h5_vec_len, h5_file_name, h5_current_time
      INTEGER ERROR
      INTEGER :: original_vec_len
      LOGICAL :: finalize ! if TRUE, finalize the datasets

*     Binary
      COMMON/OUTBINARY/ B_M1(KMAX), B_M2(KMAX),
     &                  B_XC1(KMAX), B_XC2(KMAX), B_XC3(KMAX),
     &                  B_VC1(KMAX), B_VC2(KMAX), B_VC3(KMAX), 
     &                  B_XR1(KMAX), B_XR2(KMAX), B_XR3(KMAX), 
     &                  B_VR1(KMAX), B_VR2(KMAX), B_VR3(KMAX), 
C     &                  B_FC1(KMAX), B_FC2(KMAX), B_FC3(KMAX), 
     &                  B_RS1(KMAX), B_L1(KMAX), B_TE1(KMAX),
     &                  B_RS2(KMAX), B_L2(KMAX), B_TE2(KMAX),
     &                  B_RC1(KMAX), B_MC1(KMAX), B_RC2(KMAX), 
     &                  B_MC2(KMAX), B_A(KMAX), B_ECC(KMAX), B_P(KMAX),
     &                  NB_KW1(KMAX), NB_NAM1(KMAX), NB_KW2(KMAX), 
     &                  NB_NAM2(KMAX), NB_KWC(KMAX), NB_NAMC(KMAX)
      REAL*4 B_M1, B_M2, B_XC1, B_XC2, B_XC3, B_VC1, B_VC2, B_VC3
      REAL*4 B_XR1, B_XR2, B_XR3, B_VR1, B_VR2, B_VR3
C      REAL*4 B_FC1, B_FC2, B_FC3
      REAL*4 B_RS1, B_L1, B_TE1, B_RS2, B_L2, B_TE2
      REAL*4 B_RC1, B_MC1, B_RC2, B_MC2,B_A, B_ECC, B_P
      INTEGER NB_KW1, NB_NAM1, NB_KW2, NB_NAM2, NB_KWC, NB_NAMC

*       CALL INIT if not yet done.
      IF (h5_file_id.eq.0) then
          CALL HDF5_INIT(h5_file_name)
      END IF 
*       IF file cannot be initilized, quit the subroutine
      IF (h5_file_id .EQ. 0) RETURN

      IF (TTOT.GT.h5_current_time) THEN
*       Needs to create new group
          h5_step = h5_step + 1
          WRITE (h5_step_name, *), h5_step
          h5_step_group_name = 'Step#' // ADJUSTL(h5_step_name)
          CALL h5gcreate_f(h5_file_id, h5_step_group_name,
     &       h5_group_id, error) ! Create the group

*       Create subgroup for binaries
          CALL h5gcreate_f(h5_group_id, 'Binaries',
     &       h5_subgroup_id, error) ! Create the subgroup

*       Close any previously opened datasets
*     Close any previously opened dset
          DO 361 I = 1, 256     
             IF(h5_dset_ids(I) .GT. 0) THEN
*               CALL h5dclose_f(h5_dset_ids(I), error)
                h5_dset_ids(I) = 0
             ENDIF
  361     CONTINUE
          original_vec_len = 0
          finalize = .TRUE.
*       Write group attributes (more to follow?)
          CALL HDF5_write_attribute_scalar_real(h5_group_id, 'Time',
     &        TTOT) ! Write attribute to the group

          CALL HDF5_write_attribute_scalar_integer(h5_subgroup_id,
     &    'N_BINARY', N_BINARY) ! Write attribute to the subgroup

          h5_current_time = TTOT
      ELSE
*       Open existing group
          WRITE (h5_step_name, *), h5_step
          h5_step_group_name = '/Step#' // ADJUSTL(h5_step_name)
          CALL h5gopen_f(h5_file_id, h5_step_group_name, h5_group_id,
     &        ERROR)

*       Create subgroup for binaries in the opened group
          CALL h5gcreate_f(h5_group_id, 'Binaries',
     &       h5_subgroup_id, error) ! Create the subgroup
          CALL HDF5_write_attribute_scalar_integer(h5_subgroup_id,
     &    'N_BINARY', N_BINARY) ! Write attribute to the group
          original_vec_len = 0
          finalize = .TRUE.
      ENDIF

      IF (N_BINARY.LE.0) RETURN
      
*       Write datasets
      CALL HDF5_write_real_vector_as_dset(h5_subgroup_id, 'XC1',
     &     B_XC1, N_BINARY, 1, h5_dset_ids(32),
     &     original_vec_len,finalize)
      CALL HDF5_write_real_vector_as_dset(h5_subgroup_id, 'XC2',
     &     B_XC2, N_BINARY, 1, h5_dset_ids(33),
     &     original_vec_len,finalize)
      CALL HDF5_write_real_vector_as_dset(h5_subgroup_id, 'XC3',
     &     B_XC3, N_BINARY, 1, h5_dset_ids(34),
     &     original_vec_len,finalize)

      CALL HDF5_write_real_vector_as_dset(h5_subgroup_id, 'VC1',
     &     B_VC1, N_BINARY, 1, h5_dset_ids(35),
     &     original_vec_len,finalize)
      CALL HDF5_write_real_vector_as_dset(h5_subgroup_id, 'VC2',
     &     B_VC2, N_BINARY, 1, h5_dset_ids(36),
     &     original_vec_len,finalize)
      CALL HDF5_write_real_vector_as_dset(h5_subgroup_id, 'VC3',
     &     B_VC3, N_BINARY, 1, h5_dset_ids(37),
     &     original_vec_len,finalize)

      CALL HDF5_write_real_vector_as_dset(h5_subgroup_id, 'M1',
     &     B_M1, N_BINARY, 1, h5_dset_ids(38),
     &     original_vec_len,finalize)
      CALL HDF5_write_real_vector_as_dset(h5_subgroup_id, 'M2',
     &     B_M2, N_BINARY, 1, h5_dset_ids(39),
     &     original_vec_len,finalize)

      CALL HDF5_write_real_vector_as_dset(h5_subgroup_id, 'XR1',
     &     B_XR1, N_BINARY, 1, h5_dset_ids(40),
     &     original_vec_len,finalize)
      CALL HDF5_write_real_vector_as_dset(h5_subgroup_id, 'XR2',
     &     B_XR2, N_BINARY, 1, h5_dset_ids(41),
     &     original_vec_len,finalize)
      CALL HDF5_write_real_vector_as_dset(h5_subgroup_id, 'XR3',
     &     B_XR3, N_BINARY, 1, h5_dset_ids(42),
     &     original_vec_len,finalize)

      CALL HDF5_write_real_vector_as_dset(h5_subgroup_id, 'VR1',
     &     B_VR1, N_BINARY, 1, h5_dset_ids(43),
     &     original_vec_len,finalize)
      CALL HDF5_write_real_vector_as_dset(h5_subgroup_id, 'VR2',
     &     B_VR2, N_BINARY, 1, h5_dset_ids(44),
     &     original_vec_len,finalize)
      CALL HDF5_write_real_vector_as_dset(h5_subgroup_id, 'VR3',
     &     B_VR3, N_BINARY, 1, h5_dset_ids(45),
     &     original_vec_len,finalize)

*     stellar evolution data
      IF(KSEV.GT.0) THEN
         CALL HDF5_write_real_vector_as_dset(h5_subgroup_id, 'RS1',
     &        B_RS1, N_BINARY, 1, h5_dset_ids(46),
     &        original_vec_len,finalize)
         CALL HDF5_write_real_vector_as_dset(h5_subgroup_id, 'L1',
     &        B_L1, N_BINARY, 1, h5_dset_ids(47),
     &        original_vec_len,finalize)
         CALL HDF5_write_real_vector_as_dset(h5_subgroup_id, 'TE1',
     &        B_TE1, N_BINARY, 1, h5_dset_ids(48),
     &        original_vec_len,finalize)
         CALL HDF5_write_real_vector_as_dset(h5_subgroup_id, 'RS2',
     &        B_RS2, N_BINARY, 1, h5_dset_ids(49),
     &        original_vec_len,finalize)
         CALL HDF5_write_real_vector_as_dset(h5_subgroup_id, 'L2',
     &        B_L2, N_BINARY, 1, h5_dset_ids(50),
     &        original_vec_len,finalize)
         CALL HDF5_write_real_vector_as_dset(h5_subgroup_id, 'TE2',
     &        B_TE2, N_BINARY, 1, h5_dset_ids(51),
     &        original_vec_len,finalize)

         CALL HDF5_write_real_vector_as_dset(h5_subgroup_id, 'RC1',
     &        B_RC1, N_BINARY, 1, h5_dset_ids(52),
     &        original_vec_len,finalize)
         CALL HDF5_write_real_vector_as_dset(h5_subgroup_id, 'MC1',
     &        B_MC1, N_BINARY, 1, h5_dset_ids(53),
     &        original_vec_len,finalize)
         CALL HDF5_write_real_vector_as_dset(h5_subgroup_id, 'RC2',
     &        B_RC2, N_BINARY, 1, h5_dset_ids(54),
     &        original_vec_len,finalize)
         CALL HDF5_write_real_vector_as_dset(h5_subgroup_id, 'MC2',
     &        B_MC2, N_BINARY, 1, h5_dset_ids(55),
     &        original_vec_len,finalize)
      END IF

      CALL HDF5_write_real_vector_as_dset(h5_subgroup_id, 'A',
     &     B_A, N_BINARY, 1, h5_dset_ids(56),
     &     original_vec_len,finalize)
      CALL HDF5_write_real_vector_as_dset(h5_subgroup_id, 'ECC',
     &     B_ECC, N_BINARY, 1, h5_dset_ids(57),
     &     original_vec_len,finalize)
      CALL HDF5_write_real_vector_as_dset(h5_subgroup_id, 'P',
     &     B_P, N_BINARY, 1, h5_dset_ids(58),
     &     original_vec_len,finalize)
      
*     stellar evolution data
      IF(KSEV.GT.0) THEN
         CALL HDF5_write_integer_vector_as_dset(h5_subgroup_id, 'KW1',
     &        NB_KW1, N_BINARY, 1, h5_dset_ids(59),
     &        original_vec_len,finalize)
         CALL HDF5_write_integer_vector_as_dset(h5_subgroup_id, 'KW2',
     &        NB_KW2, N_BINARY, 1, h5_dset_ids(61),
     &        original_vec_len,finalize)
      END IF

      CALL HDF5_write_integer_vector_as_dset(h5_subgroup_id, 'KWC',
     &     NB_KWC, N_BINARY, 1, h5_dset_ids(63),
     &     original_vec_len,finalize)

      CALL HDF5_write_integer_vector_as_dset(h5_subgroup_id, 'NAM1',
     &     NB_NAM1, N_BINARY, 1, h5_dset_ids(60),
     &     original_vec_len,finalize)
      CALL HDF5_write_integer_vector_as_dset(h5_subgroup_id, 'NAM2',
     &     NB_NAM2, N_BINARY, 1, h5_dset_ids(62),
     &     original_vec_len,finalize)
      CALL HDF5_write_integer_vector_as_dset(h5_subgroup_id, 'NAMC',
     &     NB_NAMC, N_BINARY, 1, h5_dset_ids(64),
     &     original_vec_len,finalize)

      CALL h5fflush_f(h5_subgroup_id, H5F_SCOPE_LOCAL_F, error)
      CALL h5fflush_f(h5_group_id, H5F_SCOPE_LOCAL_F, error)
      CALL h5gclose_f(h5_subgroup_id, error) 
      CALL h5gclose_f(h5_group_id, error) 
#else
      INCLUDE 'params.h'
*     Binary
      COMMON/OUTBINARY/ B_M1(KMAX), B_M2(KMAX),
     &                  B_XC1(KMAX), B_XC2(KMAX), B_XC3(KMAX),
     &                  B_VC1(KMAX), B_VC2(KMAX), B_VC3(KMAX), 
     &                  B_XR1(KMAX), B_XR2(KMAX), B_XR3(KMAX), 
     &                  B_VR1(KMAX), B_VR2(KMAX), B_VR3(KMAX), 
C     &                  B_FC1(KMAX), B_FC2(KMAX), B_FC3(KMAX), 
     &                  B_RS1(KMAX), B_L1(KMAX), B_TE1(KMAX),
     &                  B_RS2(KMAX), B_L2(KMAX), B_TE2(KMAX),
     &                  B_RC1(KMAX), B_MC1(KMAX), B_RC2(KMAX), 
     &                  B_MC2(KMAX), B_A(KMAX), B_ECC(KMAX), B_P(KMAX),
     &                  NB_KW1(KMAX), NB_NAM1(KMAX), NB_KW2(KMAX), 
     &                  NB_NAM2(KMAX), NB_KWC(KMAX), NB_NAMC(KMAX)
      COMMON/BINARYOUT/ DTOUT
      REAL*8 DTOUT,TTOT

      REAL*4 B_M1, B_M2, B_XC1, B_XC2, B_XC3, B_VC1, B_VC2, B_VC3
      REAL*4 B_XR1, B_XR2, B_XR3, B_VR1, B_VR2, B_VR3
C      REAL*4 B_FC1, B_FC2, B_FC3
      REAL*4 B_RS1, B_L1, B_TE1, B_RS2, B_L2, B_TE2
      REAL*4 B_RC1, B_MC1, B_RC2, B_MC2,B_A, B_ECC, B_P
      INTEGER NB_KW1, NB_NAM1, NB_KW2, NB_NAM2, NB_KWC, NB_NAMC
      INTEGER N_BINARY,KSEV,KMODE
      CHARACTER(LEN=64) :: filename
      CHARACTER(LEN=20) :: TCHAR

      IF (N_BINARY.EQ.0) RETURN
      
      call string_left(TCHAR,TTOT,DTOUT)
      filename='binary.40_'//trim(TCHAR)
      IF(KMODE.EQ.1.OR.KMODE.EQ.3) THEN
       OPEN (UNIT=40,STATUS='UNKNOWN',FORM='UNFORMATTED',FILE=filename)
       WRITE (40) N_BINARY, NB_NAM1(1:N_BINARY), NB_NAM2(1:N_BINARY),
     &     NB_NAMC(1:N_BINARY), B_M1(1:N_BINARY), B_M2(1:N_BINARY),
     &     B_XC1(1:N_BINARY), B_XC2(1:N_BINARY), B_XC3(1:N_BINARY),
     &     B_VC1(1:N_BINARY), B_VC2(1:N_BINARY), B_VC3(1:N_BINARY), 
     &     B_XR1(1:N_BINARY), B_XR2(1:N_BINARY), B_XR3(1:N_BINARY), 
     &     B_VR1(1:N_BINARY), B_VR2(1:N_BINARY), B_VR3(1:N_BINARY),
     &     B_A(1:N_BINARY), B_ECC(1:N_BINARY), B_P(1:N_BINARY)
       IF(KSEV.GT.0) THEN
         WRITE (40) B_RS1(1:N_BINARY), B_RS2(1:N_BINARY),
     &        B_L1(1:N_BINARY), B_L2(1:N_BINARY),
     &        B_TE1(1:N_BINARY), B_TE2(1:N_BINARY),
     &        B_MC1(1:N_BINARY), B_MC2(1:N_BINARY),
     &        B_RC1(1:N_BINARY), B_RC2(1:N_BINARY), 
     &        NB_KW1(1:N_BINARY), NB_KW2(1:N_BINARY), NB_KWC(1:N_BINARY)
       END IF
      ELSE IF(KMODE.EQ.2.OR.KMODE.EQ.4) THEN
       OPEN (UNIT=40,STATUS='UNKNOWN',FORM='FORMATTED',FILE=filename)
       WRITE (40,*) '## N_BINARY ',N_BINARY, ' Time[NB]', TTOT
       IF(KSEV.GT.0) THEN
          DO K=1,N_BINARY
             WRITE (40,*) NB_NAM1(K), NB_NAM2(K),
     &         NB_NAMC(K), B_M1(K), B_M2(K),
     &         B_XC1(K), B_XC2(K), B_XC3(K),
     &         B_VC1(K), B_VC2(K), B_VC3(K), 
     &         B_XR1(K), B_XR2(K), B_XR3(K), 
     &         B_VR1(K), B_VR2(K), B_VR3(K),
     &         B_A(K), B_ECC(K), B_P(K),
     &         B_RS1(K), B_RS2(K),
     &         B_L1(K), B_L2(K),
     &         B_TE1(K), B_TE2(K),
     &         B_MC1(K), B_MC2(K),
     &         B_RC1(K), B_RC2(K), 
     &         NB_KW1(K), NB_KW2(K), NB_KWC(K)
          END DO
        ELSE
           DO K=1,N_BINARY
              WRITE (40,*) NB_NAM1(K), NB_NAM2(K),
     &          NB_NAMC(K), B_M1(K), B_M2(K),
     &          B_XC1(K), B_XC2(K), B_XC3(K),
     &          B_VC1(K), B_VC2(K), B_VC3(K), 
     &          B_XR1(K), B_XR2(K), B_XR3(K), 
     &          B_VR1(K), B_VR2(K), B_VR3(K),
     &          B_A(K), B_ECC(K), B_P(K)
           END DO
        END IF
       END IF

         
      CLOSE(40)
      
#endif
      END





      SUBROUTINE output_merger(TTOT,N_MERGER,KSEV,KMODE)
#ifdef H5OUTPUT
      USE HDF5
      INCLUDE 'params.h'
      INTEGER :: h5_file_id, h5_step, h5_group_id, h5_dset_ids(256)
      INTEGER :: h5_subgroup_id
      INTEGER :: h5_vec_len
      REAL*8 :: h5_current_time,TTOT
      INTEGER :: N_MERGER,KSEV,KMODE
      CHARACTER(LEN=64) :: h5_file_name
      CHARACTER(LEN=64) :: h5_step_name
      CHARACTER(LEN=20) :: h5_step_group_name
      COMMON/h5part/ h5_file_id, h5_step, h5_group_id, h5_dset_ids,
     &     h5_vec_len, h5_file_name, h5_current_time
      INTEGER ERROR
      INTEGER :: original_vec_len
      LOGICAL :: finalize

*     Merger
      COMMON/OUTMERGER/ M_M1(MMAX), M_M2(MMAX), M_M3(MMAX),
     &                  M_XC1(MMAX), M_XC2(MMAX), M_XC3(MMAX),
     &                  M_VC1(MMAX), M_VC2(MMAX), M_VC3(MMAX),
     &                  M_XR01(MMAX), M_XR02(MMAX), M_XR03(MMAX), 
     &                  M_VR01(MMAX), M_VR02(MMAX), M_VR03(MMAX), 
     &                  M_XR11(MMAX), M_XR12(MMAX), M_XR13(MMAX), 
     &                  M_VR11(MMAX), M_VR12(MMAX), M_VR13(MMAX), 
C     &                  M_FC1(MMAX), M_FC2(MMAX), M_FC3(MMAX),
     &                  M_RS1(MMAX), M_L1(MMAX), M_TE1(MMAX), 
     &                  M_RS2(MMAX), M_L2(MMAX), M_TE2(MMAX), 
     &                  M_RS3(MMAX), M_L3(MMAX), M_TE3(MMAX), 
     &                  M_RC1(MMAX), M_MC1(MMAX), M_RC2(MMAX), 
     &                  M_MC2(MMAX), M_RC3(MMAX), M_MC3(MMAX),
     &                  M_A0(MMAX), M_ECC0(MMAX), M_P0(MMAX),
     &                  M_A1(MMAX), M_ECC1(MMAX), M_P1(MMAX),
     &                  NM_KW1(MMAX), NM_NAM1(MMAX), NM_KW2(MMAX),
     &                  NM_NAM2(MMAX), NM_KW3(MMAX), NM_NAM3(MMAX),
     &                  NM_KWC(MMAX), NM_NAMC(MMAX)
      REAL*4 M_M1, M_M2, M_M3, M_XC1, M_XC2, M_XC3, M_VC1, M_VC2, M_VC3
      REAL*4 M_XR01, M_XR02, M_XR03, M_VR01, M_VR02, M_VR03
      REAL*4 M_XR11, M_XR12, M_XR13, M_VR11, M_VR12, M_VR13
C      REAL*4 M_FC1, M_FC2, M_FC3
      REAL*4 M_RS1, M_L1, M_TE1, M_RS2, M_L2, M_TE2, M_RS3, M_L3, M_TE3
      REAL*4 M_RC1, M_MC1, M_RC2, M_MC2, M_RC3, M_MC3
      REAL*4 M_A0, M_ECC0, M_P0, M_A1, M_ECC1, M_P1
      INTEGER NM_KW1, NM_NAM1, NM_KW2, NM_NAM2
      INTEGER NM_KW3, NM_NAM3, NM_KWC, NM_NAMC

*       CALL INIT if not yet done.
      IF (h5_file_id.eq.0) then
          CALL HDF5_INIT(h5_file_name)
      END IF 
*       IF file cannot be initilized, quit the subroutine
      IF (h5_file_id .EQ. 0) RETURN

      
      IF (TTOT.GT.h5_current_time) THEN
*       Needs to create new group
           h5_step = h5_step + 1 
          WRITE (h5_step_name, *), h5_step
          h5_step_group_name = 'Step#' // ADJUSTL(h5_step_name)
          CALL h5gcreate_f(h5_file_id, h5_step_group_name,
     &       h5_group_id, error) ! Create the group

*       Create subgroup for mergers
          CALL h5gcreate_f(h5_group_id, 'Mergers',
     &       h5_subgroup_id, error) ! Create the subgroup

*       Close any previously opened datasets
          DO 361 I = 1, 256     ! Close any previously opened dset
             IF(h5_dset_ids(I) .GT. 0) THEN
*     CALL h5dclose_f(h5_dset_ids(I), error)
                h5_dset_ids(I) = 0
             ENDIF
  361    CONTINUE
         original_vec_len = 0
         finalize = .TRUE.
*       Write group attributes (more to follow?)
          CALL HDF5_write_attribute_scalar_real(h5_group_id, 'Time',
     &        TTOT) ! Write attribute to the group

          CALL HDF5_write_attribute_scalar_integer(h5_subgroup_id,
     &    'N_MERGER', N_MERGER) ! Write attribute to the group

          h5_current_time = TTOT

      ELSE
*       Open existing group
          WRITE (h5_step_name, *), h5_step
          h5_step_group_name = '/Step#' // ADJUSTL(h5_step_name)
          CALL h5gopen_f(h5_file_id, h5_step_group_name, h5_group_id,
     &        ERROR)

*       Create subgroup for binaries in the opened group
          CALL h5gcreate_f(h5_group_id, 'Mergers',
     &       h5_subgroup_id, error)  ! Create the subgroup
          CALL HDF5_write_attribute_scalar_integer(h5_subgroup_id,
     &    'N_MERGER', N_MERGER) ! Write attribute to the group
          original_vec_len = 0
          finalize = .TRUE.
      ENDIF 

      IF (N_MERGER.LE.0) RETURN

*       Write datasets
      CALL HDF5_write_real_vector_as_dset(h5_subgroup_id, 'XC1',
     &     M_XC1, N_MERGER, 1, h5_dset_ids(82),
     &     original_vec_len,finalize)
      CALL HDF5_write_real_vector_as_dset(h5_subgroup_id, 'XC2',
     &     M_XC2, N_MERGER, 1, h5_dset_ids(83),
     &     original_vec_len,finalize)
      CALL HDF5_write_real_vector_as_dset(h5_subgroup_id, 'XC3',
     &     M_XC3, N_MERGER, 1, h5_dset_ids(84),
     &     original_vec_len,finalize)

      CALL HDF5_write_real_vector_as_dset(h5_subgroup_id, 'VC1',
     &     M_VC1, N_MERGER, 1, h5_dset_ids(85),
     &     original_vec_len,finalize)
      CALL HDF5_write_real_vector_as_dset(h5_subgroup_id, 'VC2',
     &     M_VC2, N_MERGER, 1, h5_dset_ids(86),
     &     original_vec_len,finalize)
      CALL HDF5_write_real_vector_as_dset(h5_subgroup_id, 'VC3',
     &     M_VC3, N_MERGER, 1, h5_dset_ids(87),
     &     original_vec_len,finalize)

      CALL HDF5_write_real_vector_as_dset(h5_subgroup_id, 'M1',
     &     M_M1, N_MERGER, 1, h5_dset_ids(88),
     &     original_vec_len,finalize)
      CALL HDF5_write_real_vector_as_dset(h5_subgroup_id, 'M2',
     &     M_M2, N_MERGER, 1, h5_dset_ids(89),
     &     original_vec_len,finalize)
      CALL HDF5_write_real_vector_as_dset(h5_subgroup_id, 'M3',
     &     M_M3, N_MERGER, 1, h5_dset_ids(90),
     &     original_vec_len,finalize)

      CALL HDF5_write_real_vector_as_dset(h5_subgroup_id, 'XR01',
     &     M_XR01, N_MERGER, 1, h5_dset_ids(91),
     &     original_vec_len,finalize)
      CALL HDF5_write_real_vector_as_dset(h5_subgroup_id, 'XR02',
     &     M_XR02, N_MERGER, 1, h5_dset_ids(92),
     &     original_vec_len,finalize)
      CALL HDF5_write_real_vector_as_dset(h5_subgroup_id, 'XR03',
     &     M_XR03, N_MERGER, 1, h5_dset_ids(93),
     &     original_vec_len,finalize)

      CALL HDF5_write_real_vector_as_dset(h5_subgroup_id, 'VR01',
     &     M_VR01, N_MERGER, 1, h5_dset_ids(94),
     &     original_vec_len,finalize)
      CALL HDF5_write_real_vector_as_dset(h5_subgroup_id, 'VR02',
     &     M_VR02, N_MERGER, 1, h5_dset_ids(95),
     &     original_vec_len,finalize)
      CALL HDF5_write_real_vector_as_dset(h5_subgroup_id, 'VR03',
     &     M_VR03, N_MERGER, 1, h5_dset_ids(96),
     &     original_vec_len,finalize)

      CALL HDF5_write_real_vector_as_dset(h5_subgroup_id, 'XR11',
     &     M_XR11, N_MERGER, 1, h5_dset_ids(97),
     &     original_vec_len,finalize)
      CALL HDF5_write_real_vector_as_dset(h5_subgroup_id, 'XR12',
     &     M_XR12, N_MERGER, 1, h5_dset_ids(98),
     &     original_vec_len,finalize)
      CALL HDF5_write_real_vector_as_dset(h5_subgroup_id, 'XR13',
     &     M_XR13, N_MERGER, 1, h5_dset_ids(99),
     &     original_vec_len,finalize)

      CALL HDF5_write_real_vector_as_dset(h5_subgroup_id, 'VR11',
     &     M_VR11, N_MERGER, 1, h5_dset_ids(100),
     &     original_vec_len,finalize)
      CALL HDF5_write_real_vector_as_dset(h5_subgroup_id, 'VR12',
     &     M_VR12, N_MERGER, 1, h5_dset_ids(101),
     &     original_vec_len,finalize)
      CALL HDF5_write_real_vector_as_dset(h5_subgroup_id, 'VR13',
     &     M_VR13, N_MERGER, 1, h5_dset_ids(102),
     &     original_vec_len,finalize)

*     stellar evolution data
      IF(KSEV.GT.0) THEN
         CALL HDF5_write_real_vector_as_dset(h5_subgroup_id, 'RS1',
     &        M_RS1, N_MERGER, 1, h5_dset_ids(103),
     &        original_vec_len,finalize)
         CALL HDF5_write_real_vector_as_dset(h5_subgroup_id, 'L1',
     &        M_L1, N_MERGER, 1, h5_dset_ids(104),
     &        original_vec_len,finalize)
         CALL HDF5_write_real_vector_as_dset(h5_subgroup_id, 'TE1',
     &        M_TE1, N_MERGER, 1, h5_dset_ids(105),
     &        original_vec_len,finalize)
         CALL HDF5_write_real_vector_as_dset(h5_subgroup_id, 'RS2',
     &        M_RS2, N_MERGER, 1, h5_dset_ids(106),
     &        original_vec_len,finalize)
         CALL HDF5_write_real_vector_as_dset(h5_subgroup_id, 'L2',
     &        M_L2, N_MERGER, 1, h5_dset_ids(107),
     &        original_vec_len,finalize)
         CALL HDF5_write_real_vector_as_dset(h5_subgroup_id, 'TE2',
     &        M_TE2, N_MERGER, 1, h5_dset_ids(108),
     &        original_vec_len,finalize)
         CALL HDF5_write_real_vector_as_dset(h5_subgroup_id, 'RS3',
     &        M_RS3, N_MERGER, 1, h5_dset_ids(109),
     &        original_vec_len,finalize)
         CALL HDF5_write_real_vector_as_dset(h5_subgroup_id, 'L3',
     &        M_L3, N_MERGER, 1, h5_dset_ids(110),
     &        original_vec_len,finalize)
         CALL HDF5_write_real_vector_as_dset(h5_subgroup_id, 'TE3',
     &        M_TE3, N_MERGER, 1, h5_dset_ids(111),
     &        original_vec_len,finalize)

         CALL HDF5_write_real_vector_as_dset(h5_subgroup_id, 'RC1',
     &        M_RC1, N_MERGER, 1, h5_dset_ids(112),
     &        original_vec_len,finalize)
         CALL HDF5_write_real_vector_as_dset(h5_subgroup_id, 'MC1',
     &        M_MC1, N_MERGER, 1, h5_dset_ids(113),
     &        original_vec_len,finalize)
         CALL HDF5_write_real_vector_as_dset(h5_subgroup_id, 'RC2',
     &        M_RC2, N_MERGER, 1, h5_dset_ids(114),
     &        original_vec_len,finalize)
         CALL HDF5_write_real_vector_as_dset(h5_subgroup_id, 'MC2',
     &        M_MC2, N_MERGER, 1, h5_dset_ids(115),
     &        original_vec_len,finalize)
         CALL HDF5_write_real_vector_as_dset(h5_subgroup_id, 'RC3',
     &        M_RC3, N_MERGER, 1, h5_dset_ids(116),
     &        original_vec_len,finalize)
         CALL HDF5_write_real_vector_as_dset(h5_subgroup_id, 'MC3',
     &        M_MC3, N_MERGER, 1, h5_dset_ids(117),
     &        original_vec_len,finalize)
      END IF

      CALL HDF5_write_real_vector_as_dset(h5_subgroup_id, 'A0',
     &     M_A0, N_MERGER, 1, h5_dset_ids(118),
     &     original_vec_len,finalize)
      CALL HDF5_write_real_vector_as_dset(h5_subgroup_id, 'ECC0',
     &     M_ECC0, N_MERGER, 1, h5_dset_ids(119),
     &     original_vec_len,finalize)
      CALL HDF5_write_real_vector_as_dset(h5_subgroup_id, 'P0',
     &     M_P0, N_MERGER, 1, h5_dset_ids(120),
     &     original_vec_len,finalize)
      CALL HDF5_write_real_vector_as_dset(h5_subgroup_id, 'A1',
     &     M_A1, N_MERGER, 1, h5_dset_ids(121),
     &     original_vec_len,finalize)
      CALL HDF5_write_real_vector_as_dset(h5_subgroup_id, 'ECC1',
     &     M_ECC1, N_MERGER, 1, h5_dset_ids(122),
     &     original_vec_len,finalize)
      CALL HDF5_write_real_vector_as_dset(h5_subgroup_id, 'P1',
     &     M_P1, N_MERGER, 1, h5_dset_ids(123),
     &     original_vec_len,finalize)

*     stellar evolution data
      IF(KSEV.GT.0) THEN
         CALL HDF5_write_integer_vector_as_dset(h5_subgroup_id, 'KW1',
     &        NM_KW1, N_MERGER, 1, h5_dset_ids(124),
     &        original_vec_len,finalize)
         CALL HDF5_write_integer_vector_as_dset(h5_subgroup_id, 'KW2',
     &        NM_KW2, N_MERGER, 1, h5_dset_ids(126),
     &        original_vec_len,finalize)
         CALL HDF5_write_integer_vector_as_dset(h5_subgroup_id, 'KW3',
     &        NM_KW3, N_MERGER, 1, h5_dset_ids(128),
     &        original_vec_len,finalize)
      END IF

      CALL HDF5_write_integer_vector_as_dset(h5_subgroup_id, 'KWC',
     &     NM_KWC, N_MERGER, 1, h5_dset_ids(130),
     &     original_vec_len,finalize)
      CALL HDF5_write_integer_vector_as_dset(h5_subgroup_id, 'NAM1',
     &     NM_NAM1, N_MERGER, 1, h5_dset_ids(125),
     &     original_vec_len,finalize)
      CALL HDF5_write_integer_vector_as_dset(h5_subgroup_id, 'NAM2',
     &     NM_NAM2, N_MERGER, 1, h5_dset_ids(127),
     &     original_vec_len,finalize)
      CALL HDF5_write_integer_vector_as_dset(h5_subgroup_id, 'NAM3',
     &     NM_NAM3, N_MERGER, 1, h5_dset_ids(129),
     &     original_vec_len,finalize)
      CALL HDF5_write_integer_vector_as_dset(h5_subgroup_id, 'NAMC',
     &     NM_NAMC, N_MERGER, 1, h5_dset_ids(131),
     &     original_vec_len,finalize)

      CALL h5fflush_f(h5_subgroup_id, H5F_SCOPE_LOCAL_F, error)
      CALL h5fflush_f(h5_group_id, H5F_SCOPE_LOCAL_F, error)
      CALL h5gclose_f(h5_subgroup_id, error) 
      CALL h5gclose_f(h5_group_id, error) 

#else
      INCLUDE 'params.h'
*     Merger
      COMMON/OUTMERGER/ M_M1(MMAX), M_M2(MMAX), M_M3(MMAX),
     &                  M_XC1(MMAX), M_XC2(MMAX), M_XC3(MMAX),
     &                  M_VC1(MMAX), M_VC2(MMAX), M_VC3(MMAX),
     &                  M_XR01(MMAX), M_XR02(MMAX), M_XR03(MMAX), 
     &                  M_VR01(MMAX), M_VR02(MMAX), M_VR03(MMAX), 
     &                  M_XR11(MMAX), M_XR12(MMAX), M_XR13(MMAX), 
     &                  M_VR11(MMAX), M_VR12(MMAX), M_VR13(MMAX), 
C     &                  M_FC1(MMAX), M_FC2(MMAX), M_FC3(MMAX),
     &                  M_RS1(MMAX), M_L1(MMAX), M_TE1(MMAX), 
     &                  M_RS2(MMAX), M_L2(MMAX), M_TE2(MMAX), 
     &                  M_RS3(MMAX), M_L3(MMAX), M_TE3(MMAX), 
     &                  M_RC1(MMAX), M_MC1(MMAX), M_RC2(MMAX), 
     &                  M_MC2(MMAX), M_RC3(MMAX), M_MC3(MMAX),
     &                  M_A0(MMAX), M_ECC0(MMAX), M_P0(MMAX),
     &                  M_A1(MMAX), M_ECC1(MMAX), M_P1(MMAX),
     &                  NM_KW1(MMAX), NM_NAM1(MMAX), NM_KW2(MMAX),
     &                  NM_NAM2(MMAX), NM_KW3(MMAX), NM_NAM3(MMAX),
     &                  NM_KWC(MMAX), NM_NAMC(MMAX)
      COMMON/BINARYOUT/ DTOUT
      REAL*8 DTOUT,TTOT

      REAL*4 M_M1, M_M2, M_M3, M_XC1, M_XC2, M_XC3, M_VC1, M_VC2, M_VC3
      REAL*4 M_XR01, M_XR02, M_XR03, M_VR01, M_VR02, M_VR03
      REAL*4 M_XR11, M_XR12, M_XR13, M_VR11, M_VR12, M_VR13
C      REAL*4 M_FC1, M_FC2, M_FC3
      REAL*4 M_RS1, M_L1, M_TE1, M_RS2, M_L2, M_TE2, M_RS3, M_L3, M_TE3
      REAL*4 M_RC1, M_MC1, M_RC2, M_MC2, M_RC3, M_MC3
      REAL*4 M_A0, M_ECC0, M_P0, M_A1, M_ECC1, M_P1
      INTEGER NM_KW1, NM_NAM1, NM_KW2, NM_NAM2
      INTEGER NM_KW3, NM_NAM3, NM_KWC, NM_NAMC
      INTEGER N_MERGER,KSEV,KMODE
      CHARACTER(LEN=20) :: TCHAR
      CHARACTER(LEN=64) :: filename

      IF (N_MERGER.EQ.0) RETURN
      
      call string_left(TCHAR,TTOT,DTOUT)
      filename='merger.40_'//trim(TCHAR)
      IF(KMODE.EQ.1.OR.KMODE.EQ.3) THEN
       OPEN (UNIT=40,STATUS='UNKNOWN',FORM='UNFORMATTED',FILE=filename)
       WRITE (40) NM_NAM1(1:N_MERGER), NM_NAM2(1:N_MERGER),
     &     NM_NAM3(1:N_MERGER), NM_NAMC(1:N_MERGER),
     &     M_M1(1:N_MERGER), M_M2(1:N_MERGER), M_M3(1:N_MERGER),
     &     M_XC1(1:N_MERGER), M_XC2(1:N_MERGER), M_XC3(1:N_MERGER),
     &     M_VC1(1:N_MERGER), M_VC2(1:N_MERGER), M_VC3(1:N_MERGER),
     &     M_XR01(1:N_MERGER), M_XR02(1:N_MERGER), M_XR03(1:N_MERGER), 
     &     M_VR01(1:N_MERGER), M_VR02(1:N_MERGER), M_VR03(1:N_MERGER), 
     &     M_XR11(1:N_MERGER), M_XR12(1:N_MERGER), M_XR13(1:N_MERGER), 
     &     M_VR11(1:N_MERGER), M_VR12(1:N_MERGER), M_VR13(1:N_MERGER)
       IF (KSEV.GT.0) THEN
      WRITE(40) M_RS1(1:N_MERGER), M_RS2(1:N_MERGER), M_RS3(1:N_MERGER),
     &        M_L1(1:N_MERGER), M_L2(1:N_MERGER), M_L3(1:N_MERGER), 
     &        M_TE1(1:N_MERGER), M_TE2(1:N_MERGER), M_TE3(1:N_MERGER), 
     &        M_RC1(1:N_MERGER), M_RC2(1:N_MERGER), M_RC3(1:N_MERGER),
     &        M_MC1(1:N_MERGER), M_MC2(1:N_MERGER), M_MC3(1:N_MERGER),
     &        M_A0(1:N_MERGER), M_ECC0(1:N_MERGER), M_P0(1:N_MERGER),
     &        M_A1(1:N_MERGER), M_ECC1(1:N_MERGER), M_P1(1:N_MERGER),
     &        NM_KW1(1:N_MERGER), NM_KW2(1:N_MERGER),
     &        NM_KW3(1:N_MERGER), NM_KWC(1:N_MERGER)
       END IF
      ELSE IF(KMODE.EQ.2.OR.KMODE.EQ.4) THEN
       OPEN (UNIT=40,STATUS='UNKNOWN',FORM='FORMATTED',FILE=filename)
       WRITE (40,*) '## N_MERGER ',N_MERGER, ' Time[NB]', TTOT
       IF(KSEV.GT.0) THEN
          DO K=1,N_MERGER
             WRITE (40) NM_NAM1(K), NM_NAM2(K),
     &            NM_NAM3(K), NM_NAMC(K),
     &            M_M1(K), M_M2(K), M_M3(K),
     &            M_XC1(K), M_XC2(K), M_XC3(K),
     &            M_VC1(K), M_VC2(K), M_VC3(K),
     &            M_XR01(K), M_XR02(K), M_XR03(K), 
     &            M_VR01(K), M_VR02(K), M_VR03(K), 
     &            M_XR11(K), M_XR12(K), M_XR13(K), 
     &            M_VR11(K), M_VR12(K), M_VR13(K),
     &            M_A0(K), M_ECC0(K), M_P0(K),
     &            M_A1(K), M_ECC1(K), M_P1(K),
     &            NM_KW1(K), NM_KW2(K),
     &            NM_KW3(K), NM_KWC(K),
     &            M_RS1(K), M_RS2(K), M_RS3(K),
     &            M_L1(K), M_L2(K), M_L3(K), 
     &            M_TE1(K), M_TE2(K), M_TE3(K), 
     &            M_MC1(K), M_MC2(K), M_MC3(K),
     &            M_RC1(K), M_RC2(K), M_RC3(K)
          END DO
       ELSE
          DO K=1,N_MERGER
             WRITE (40) NM_NAM1(K), NM_NAM2(K),
     &            NM_NAM3(K), NM_NAMC(K),
     &            M_M1(K), M_M2(K), M_M3(K),
     &            M_XC1(K), M_XC2(K), M_XC3(K),
     &            M_VC1(K), M_VC2(K), M_VC3(K),
     &            M_XR01(K), M_XR02(K), M_XR03(K), 
     &            M_VR01(K), M_VR02(K), M_VR03(K), 
     &            M_XR11(K), M_XR12(K), M_XR13(K), 
     &            M_VR11(K), M_VR12(K), M_VR13(K),
     &            M_A0(K), M_ECC0(K), M_P0(K),
     &            M_A1(K), M_ECC1(K), M_P1(K)
          END DO
       END IF
      END IF
      CLOSE(40)
      
#endif
      END
