c ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c RCS version control information:
c $Header: ErrPack.f,v 1.2 97/03/18 17:06:50 wiscombe Exp $
c ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      SUBROUTINE  ErrMsg( MESSAG, FATAL )

c        Print out a warning or error message;  abort if error

      LOGICAL       FATAL, MsgLim
      CHARACTER*(*) MESSAG
      INTEGER       MaxMsg, NumMsg
      SAVE          MaxMsg, NumMsg, MsgLim
      DATA NumMsg / 0 /,  MaxMsg / 100 /,  MsgLim / .FALSE. /


      IF ( FATAL )  THEN
         WRITE ( *, '(/,2A,/)' )  ' ******* ERROR >>>>>>  ', MESSAG
         STOP
      END IF

      NumMsg = NumMsg + 1
      IF( MsgLim )  RETURN

      IF ( NumMsg.LE.MaxMsg )  THEN
         WRITE ( *, '(/,2A,/)' )  ' ******* WARNING >>>>>>  ', MESSAG
      ELSE
         WRITE ( *,99 )
         MsgLim = .True.
      ENDIF

      RETURN

   99 FORMAT( //,' >>>>>>  TOO MANY WARNING MESSAGES --  ',
     &   'They will no longer be printed  <<<<<<<', // )
      END

      LOGICAL FUNCTION  WrtBad ( VarNam )

c          Write names of erroneous variables and return 'TRUE'

c      INPUT :   VarNam = Name of erroneous variable to be written
c                         ( CHARACTER, any length )

      CHARACTER*(*)  VarNam
      INTEGER        MaxMsg, NumMsg
      SAVE  NumMsg, MaxMsg
      DATA  NumMsg / 0 /,  MaxMsg / 50 /


      WrtBad = .TRUE.
      NumMsg = NumMsg + 1
      WRITE ( *, '(3A)' )  ' ****  Input variable  ', VarNam,
     &                     '  in error  ****'
      IF ( NumMsg.EQ.MaxMsg )
     &   CALL  ErrMsg ( 'Too many input errors.  Aborting...', .TRUE. )

      RETURN
      END

      LOGICAL FUNCTION  WrtDim ( DimNam, MinVal )

c          Write name of too-small symbolic dimension and
c          the value it should be increased to;  return 'TRUE'

c      INPUT :  DimNam = Name of symbolic dimension which is too small
c                        ( CHARACTER, any length )
c               Minval = Value to which that dimension should be
c                        increased (at least)

      CHARACTER*(*)  DimNam
      INTEGER        MinVal


      WRITE ( *, '(/,3A,I7)' )  ' ****  Symbolic dimension  ', DimNam,
     &                     '  should be increased to at least ', MinVal
      WrtDim = .TRUE.

      RETURN
      END

      LOGICAL FUNCTION  TstBad( VarNam, RelErr )

c       Write name (VarNam) of variable failing self-test and its
c       percent error from the correct value;  return  'FALSE'.

      CHARACTER*(*)  VarNam
      REAL           RelErr


      TstBad = .FALSE.
      WRITE( *, '(/,3A,1P,E11.2,A)' )
     &       ' Output variable ', VarNam,' differed by ', 100.*RelErr,
     &       ' per cent from correct value.  Self-test failed.'

      RETURN
      END

c ---------------------------------------------------------------------
c  Fortran-90 versions of machine-constant routines R1MACH, D1MACH, I1MACH
c
c  {R,D,I}1MACH revisited: no more uncommenting DATA statements
c  
c  Presented at the IFIP WG 2.5 International Workshop on 
c  "Current Directions in Numerical Software and High Performance 
c  Computing", 19 - 20 October 1995, Kyoto, Japan. 
c  
c  The widely-used original routines were modified to use Fortran-90 
c  intrinsic functions.  This was not completely possible with I1MACH, 
c  which returns some parameters (logical unit numbers of standard
c  input, standard output, and standard error) that may require
c  user customization. 
c  
c  David Gay (dmg@bell-labs.com)
c  Eric Grosse (ehg@bell-labs.com)
c  Bell Laboratories
c  700 Mountain Avenue
c  Murray Hill, New Jersey 07974-0636
c  USA 
c  
c  References:
c  
c  David Gay and Eric Grosse, Comment on Algorithm 528, Bell Labs, Murray 
c    Hill, NJ. submitted to ACM Transactions on Mathematical Software,
c    August 1996.
c
c http://www.nsc.liu.se/~boein/ifip/kyoto/workshop-info/proceedings/einarsson
c    /d1mach.html  (THIS WEB SITE WORKED AS OF APR 2000)
c -------------------------------------------------------------------------


      REAL FUNCTION R1MACH (I)
c
c   R1MACH can be used to obtain machine-dependent parameters for
c   single precision numbers.  The results for various values of I are:
c
c   R1MACH(1) = B**(EMIN-1), the smallest positive magnitude.
c   R1MACH(2) = B**EMAX*(1 - B**(-T)), the largest magnitude.
c   R1MACH(3) = B**(-T), the smallest relative spacing.
c   R1MACH(4) = B**(1-T), the largest relative spacing.
c   R1MACH(5) = LOG10(B)
c
c   Assume single precision numbers are represented in the T-digit,
c   base-B form
c
c              sign (B**E)*( (X(1)/B) + ... + (X(T)/B**T) )
c
c   where 0 <= X(I) < B for I=1,...,T; 0 < X(1); and EMIN <= E <= EMAX.
c
c   The values of B, T, EMIN and EMAX are provided in I1MACH as follows:
c   I1MACH(10) = B, the base.
c   I1MACH(11) = T, the number of base-B digits.
c   I1MACH(12) = EMIN, the smallest exponent E.
c   I1MACH(13) = EMAX, the largest exponent E.
c
c***REFERENCES  
c
c  P. Fox, A. Hall and N. Schryer, Framework for a portable library,
c     ACM Transactions on Mathematical Software 4, 177-188 (1978).
c
c  David Gay and Eric Grosse, Comment on Algorithm 528, Bell Labs, Murray 
c     Hill, NJ. submitted to ACM Transactions on Mathematical Software,
c     August 1996. 
c
c***REVISION HISTORY  (YYMMDD)
c   790101  DATE WRITTEN
c   960329  Modified for Fortran 90 (BE after suggestions by Eric Grosse)      
c --------------------------------------------------------------------

      IMPLICIT NONE
      INTEGER :: I
      REAL :: B, X = 1.0

      B = RADIX(X)

      SELECT CASE (I)
        CASE (1)
          R1MACH = TINY(X)            ! smallest positive magnitude.
        CASE (2)
          R1MACH = HUGE(X)            ! largest magnitude.
        CASE (3)
          R1MACH = B**(-DIGITS(X))    ! smallest relative spacing.
        CASE (4)
          R1MACH = B**(1-DIGITS(X))   ! largest relative spacing.
        CASE (5)
          R1MACH = LOG10(B)
        CASE DEFAULT
          STOP 'R1MACH -- input argument out of bounds'
      END SELECT

      RETURN
      END FUNCTION R1MACH


      DOUBLE PRECISION FUNCTION D1MACH (I)
c
c   D1MACH can be used to obtain machine-dependent parameters for
c   double precision numbers.  The results for various values of I are:
c
c   D1MACH(1) = B**(EMIN-1), the smallest positive magnitude.
c   D1MACH(2) = B**EMAX*(1 - B**(-T)), the largest magnitude.
c   D1MACH(3) = B**(-T), the smallest relative spacing.
c   D1MACH(4) = B**(1-T), the largest relative spacing.
c   D1MACH(5) = LOG10(B)
c
c   Assume double precision numbers are represented in the T-digit,
c   base-B form
c
c        sign (B**E)*( (X(1)/B) + ... + (X(T)/B**T) )
c
c   where 0 <= X(I) < B for I=1,...,T; 0 < X(1); and EMIN <= E <= EMAX.
c
c   The values of B, T, EMIN and EMAX are provided in I1MACH as follows:
c   I1MACH(10) = B, the base.
c   I1MACH(11) = T, the number of base-B digits.
c   I1MACH(12) = EMIN, the smallest exponent E.
c   I1MACH(13) = EMAX, the largest exponent E.
c
c***REFERENCES  
c
c  P. Fox, A. Hall and N. Schryer, Framework for a portable library,
c     ACM Transactions on Mathematical Software 4, 177-188 (1978).
c
c  David Gay and Eric Grosse, Comment on Algorithm 528, Bell Labs, Murray 
c    Hill, NJ. submitted to ACM Transactions on Mathematical Software,
c    August 1996. 
c
c***REVISION HISTORY  (YYMMDD)
c   790101  DATE WRITTEN
c   960329  Modified for Fortran 90 (BE after suggestions by Eric Grosse)      
c --------------------------------------------------------------------

      IMPLICIT NONE
      INTEGER :: I
      DOUBLE PRECISION :: B, X = 1.D0

      B = RADIX(X)

      SELECT CASE (I)
        CASE (1)
          D1MACH = TINY(X)            ! smallest positive magnitude.
        CASE (2)
          D1MACH = HUGE(X)            ! largest magnitude.
        CASE (3)
          D1MACH = B**(-DIGITS(X))    ! smallest relative spacing.
        CASE (4)
          D1MACH = B**(1-DIGITS(X))   ! largest relative spacing.
        CASE (5)
          D1MACH = LOG10(B)
        CASE DEFAULT
          STOP 'D1MACH -- input arg out of bounds'
      END SELECT

      RETURN
      END FUNCTION D1MACH


      INTEGER FUNCTION I1MACH (I)
c
c   I1MACH can be used to obtain machine-dependent parameters for the
c   local machine environment.  The results for various values of I are:
c
c   I/O unit numbers (**MAY REQUIRE USER CUSTOMIZATION**):
c     I1MACH( 1) = the standard input unit.
c     I1MACH( 2) = the standard output unit.
c     I1MACH( 3) = the standard punch unit (obsolete, will cause error)
c     I1MACH( 4) = the standard error message unit.
c                  (the error message unit is usually 0 in UNIX systems)
c
c   Words:
c     I1MACH( 5) = the number of bits per integer storage unit.
c     I1MACH( 6) = the number of characters per integer storage unit.
c                  (obsolete, will cause an error)
c
c   Integers:
c     assume integers are represented in the S-digit, base-A form
c
c          sign ( X(S-1)*A**(S-1) + ... + X(1)*A + X(0) )
c
c     where 0 <= X(I) < A for I=0,...,S-1.
c
c     I1MACH( 7) = A, the base.
c     I1MACH( 8) = S, the number of base-A digits.
c     I1MACH( 9) = A**S - 1, the largest magnitude.
c
c   Floating-Point Numbers:
c     Assume floating-point numbers are represented in the T-digit,
c     base-B form
c                sign (B**E)*( (X(1)/B) + ... + (X(T)/B**T) )
c
c     where 0 <= X(I) .LT. B for I=1,...,T; 0 < X(1); and EMIN <= E <= EMAX.
c
c     I1MACH(10) = B, the base.
c
c   Single-Precision:
c     I1MACH(11) = T, the number of base-B digits.
c     I1MACH(12) = EMIN, the smallest exponent E.
c     I1MACH(13) = EMAX, the largest exponent E.
c
c   Double-Precision:
c     I1MACH(14) = T, the number of base-B digits.
c     I1MACH(15) = EMIN, the smallest exponent E.
c     I1MACH(16) = EMAX, the largest exponent E.
c
c***REFERENCES  
c
c  P. Fox, A. Hall and N. Schryer, Framework for a portable library,
c     ACM Transactions on Mathematical Software 4, 177-188 (1978).
c
c  David Gay and Eric Grosse, Comment on Algorithm 528, Bell Labs, Murray 
c    Hill, NJ. submitted to ACM Transactions on Mathematical Software,
c    August 1996. 
c
c***REVISION HISTORY  (YYMMDD)
c   750101  DATE WRITTEN
c   960411  Modified for Fortran 90 (BE after suggestions by Eric Grosse)    
c --------------------------------------------------------------------

      IMPLICIT NONE
      INTEGER :: I
      REAL :: X_single  = 1.0
      DOUBLE PRECISION :: X_double = 1.D0

      SELECT CASE (I)
        CASE (1)
          I1MACH = 5 ! Input unit
        CASE (2)
          I1MACH = 6 ! Output unit
        CASE (3)
          STOP 'I1MACH: input arg = 3 is obsolete'
        CASE (4)
          I1MACH = 0 ! Error message unit
        CASE (5)
          I1MACH = BIT_SIZE(I)
        CASE (6)
          STOP 'I1MACH: input arg = 6 is obsolete'
        CASE (7)
          I1MACH = RADIX(1)
        CASE (8)
          I1MACH = BIT_SIZE(I) - 1
        CASE (9)
          I1MACH = HUGE(1)
        CASE (10)
          I1MACH = RADIX(X_single)
        CASE (11)
          I1MACH = DIGITS(X_single)
        CASE (12)
          I1MACH = MINEXPONENT(X_single)
        CASE (13)
          I1MACH = MAXEXPONENT(X_single)
        CASE (14)
          I1MACH = DIGITS(X_double)
        CASE (15)
          I1MACH = MINEXPONENT(X_double)
        CASE (16)
          I1MACH = MAXEXPONENT(X_double) 
        CASE DEFAULT
          STOP 'I1MACH: input argument out of bounds'
      END SELECT

      RETURN
      END FUNCTION I1MACH

c ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c
c     LINEAR-EQUATION-SOLVING LAPACK (v3.0) ROUTINES FOR DISORT v1.3
c     (current as of Nov 99)
c
c NOTES:  
c
c (1) If possible, use locally optimized versions of these routines
c     instead of these Fortran versions.  A significant portion 
c     of DISORT computer time is spent in these routines.
c
c (2) Prebuilt LAPACK libraries are available for a variety of
c     computers.  These will run much faster.  See:
c
c        http://www.netlib.org/lapack/archives/
c
c (3) Upgrades to LAPACK are available from 
c
c        http://www.netlib.org/lapack/
c        http://netlib.bell-labs.com/netlib/master/readme.html
c
c ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c
c CHANGES MADE TO ORIGINAL LAPACK:
c
c (1) Major simplifications in routine ILAENV.
c
c (2) Removal of the following comment block from the beginning of
c     each routine:
c  -- LAPACK routine (version 2.0) --
c     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
c     Courant Institute, Argonne National Lab, and Rice University
c
c (3) Other cosmetic edits.  
c
c ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c
c Call tree (omitting calls to error msg routine XERBLA):
c
c    SGBSV
c        SGBTRF
c            ILAENV
c            ISAMAX [BLAS-1]
c            SCOPY [BLAS-1]
c            SGBTF2
c                ISAMAX [BLAS-1]
c                SGER [BLAS-2]
c                SSCAL [BLAS-1]
c                SSWAP [BLAS-1]
c            SGEMM [BLAS-3]
c                LSAME
c            SGER [BLAS-2]
c            SLASWP
c                SSWAP [BLAS-1]
c            SSCAL [BLAS-1]
c            SSWAP [BLAS-1]
c            STRSM [BLAS-3]
c                LSAME
c        SGBTRS
c            LSAME [BLAS-2]
c            SGEMV [BLAS-2]
c                LSAME
c            SGER [BLAS-2]
c            SSWAP [BLAS-1]
c            STBSV [BLAS-2]
c                LSAME
c    SGESV
c        SGETRF
c            SGEMM [BLAS-3]
c                LSAME
c            SGETF2
c                ISAMAX [BLAS-1]
c                SGER [BLAS-2]
c                SSCAL [BLAS-1]
c                SSWAP [BLAS-1]
c            SLASWP
c                SSWAP [BLAS-1]
c            STRSM [BLAS-3]
c            ILAENV
c        SGETRS
c            LSAME [BLAS-2]
c            SLASWP
c                SSWAP [BLAS-1]
c            STRSM [BLAS-3]
c                LSAME

c  The BLAS-1,2,3 routines are in a separate file.  See if you have
c  optimized local versions of BLAS routines before using that file.
c ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


      SUBROUTINE SGBSV( N, KL, KU, NRHS, AB, LDAB, IPIV, B, LDB, INFO )
c
c  -- LAPACK driver routine (version 2.0) --
c     March 31, 1993 
c
c     .. Scalar Arguments ..
      INTEGER            INFO, KL, KU, LDAB, LDB, N, NRHS
c     ..
c     .. Array Arguments ..
      INTEGER            IPIV( * )
      REAL               AB( LDAB, * ), B( LDB, * )
c     ..
c
c  Purpose
c  =======
c
c  SGBSV computes the solution to a real system of linear equations
c  A * X = B, where A is a band matrix of order N with KL subdiagonals
c  and KU superdiagonals, and X and B are N-by-NRHS matrices.
c
c  The LU decomposition with partial pivoting and row interchanges is
c  used to factor A as A = L * U, where L is a product of permutation
c  and unit lower triangular matrices with KL subdiagonals, and U is
c  upper triangular with KL+KU superdiagonals.  The factored form of A
c  is then used to solve the system of equations A * X = B.
c
c  Arguments
c  =========
c
c  N       (input) INTEGER
c          The number of linear equations, i.e., the order of the
c          matrix A.  N >= 0.
c
c  KL      (input) INTEGER
c          The number of subdiagonals within the band of A.  KL >= 0.
c
c  KU      (input) INTEGER
c          The number of superdiagonals within the band of A.  KU >= 0.
c
c  NRHS    (input) INTEGER
c          The number of right hand sides, i.e., the number of columns
c          of the matrix B.  NRHS >= 0.
c
c  AB      (input/output) REAL array, dimension (LDAB,N)
c          On entry, the matrix A in band storage, in rows KL+1 to
c          2*KL+KU+1; rows 1 to KL of the array need not be set.
c          The j-th column of A is stored in the j-th column of the
c          array AB as follows:
c          AB(KL+KU+1+i-j,j) = A(i,j) for max(1,j-KU)<=i<=min(N,j+KL)
c          On exit, details of the factorization: U is stored as an
c          upper triangular band matrix with KL+KU superdiagonals in
c          rows 1 to KL+KU+1, and the multipliers used during the
c          factorization are stored in rows KL+KU+2 to 2*KL+KU+1.
c          See below for further details.
c
c  LDAB    (input) INTEGER
c          The leading dimension of the array AB.  LDAB >= 2*KL+KU+1.
c
c  IPIV    (output) INTEGER array, dimension (N)
c          The pivot indices that define the permutation matrix P;
c          row i of the matrix was interchanged with row IPIV(i).
c
c  B       (input/output) REAL array, dimension (LDB,NRHS)
c          On entry, the N-by-NRHS right hand side matrix B.
c          On exit, if INFO = 0, the N-by-NRHS solution matrix X.
c
c  LDB     (input) INTEGER
c          The leading dimension of the array B.  LDB >= max(1,N).
c
c  INFO    (output) INTEGER
c          = 0:  successful exit
c          < 0:  if INFO = -i, the i-th argument had an illegal value
c          > 0:  if INFO = i, U(i,i) is exactly zero.  The factorization
c                has been completed, but the factor U is exactly
c                singular, and the solution has not been computed.
c
c  Further Details
c  ===============
c
c  The band storage scheme is illustrated by the following example, when
c  M = N = 6, KL = 2, KU = 1:
c
c  On entry:                       On exit:
c
c      *    *    *    +    +    +       *    *    *   u14  u25  u36
c      *    *    +    +    +    +       *    *   u13  u24  u35  u46
c      *   a12  a23  a34  a45  a56      *   u12  u23  u34  u45  u56
c     a11  a22  a33  a44  a55  a66     u11  u22  u33  u44  u55  u66
c     a21  a32  a43  a54  a65   *      m21  m32  m43  m54  m65   *
c     a31  a42  a53  a64   *    *      m31  m42  m53  m64   *    *
c
c  Array elements marked * are not used by the routine; elements marked
c  + need not be set on entry, but are required by the routine to store
c  elements of U because of fill-in resulting from the row interchanges.
c
c  =====================================================================
c
c     .. External Subroutines ..
      EXTERNAL           SGBTRF, SGBTRS, XERBLA
c     ..
c     .. Intrinsic Functions ..
      INTRINSIC          MAX
c     ..

c     Test the input parameters.

      INFO = 0
      IF( N.LT.0 ) THEN
         INFO = -1
      ELSE IF( KL.LT.0 ) THEN
         INFO = -2
      ELSE IF( KU.LT.0 ) THEN
         INFO = -3
      ELSE IF( NRHS.LT.0 ) THEN
         INFO = -4
      ELSE IF( LDAB.LT.2*KL+KU+1 ) THEN
         INFO = -6
      ELSE IF( LDB.LT.MAX( N, 1 ) ) THEN
         INFO = -9
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'SGBSV ', -INFO )
         RETURN
      END IF

c     Compute the LU factorization of the band matrix A.

      CALL SGBTRF( N, N, KL, KU, AB, LDAB, IPIV, INFO )
      IF( INFO.EQ.0 ) THEN

c        Solve the system A*X = B, overwriting B with X.

         CALL SGBTRS( 'No transpose', N, KL, KU, NRHS, AB, LDAB, IPIV,
     &                B, LDB, INFO )
      END IF
      RETURN

c          End of SGBSV
      END

      SUBROUTINE SGBTRF( M, N, KL, KU, AB, LDAB, IPIV, INFO )
c
c     February 29, 1992
c
c     .. Scalar Arguments ..
      INTEGER            INFO, KL, KU, LDAB, M, N
c     ..
c     .. Array Arguments ..
      INTEGER            IPIV( * )
      REAL               AB( LDAB, * )
c     ..
c
c  Purpose
c  =======
c
c  SGBTRF computes an LU factorization of a real m-by-n band matrix A
c  using partial pivoting with row interchanges.
c
c  This is the blocked version of the algorithm, calling Level 3 BLAS.
c
c  Arguments
c  =========
c
c  M       (input) INTEGER
c          The number of rows of the matrix A.  M >= 0.
c
c  N       (input) INTEGER
c          The number of columns of the matrix A.  N >= 0.
c
c  KL      (input) INTEGER
c          The number of subdiagonals within the band of A.  KL >= 0.
c
c  KU      (input) INTEGER
c          The number of superdiagonals within the band of A.  KU >= 0.
c
c  AB      (input/output) REAL array, dimension (LDAB,N)
c          On entry, the matrix A in band storage, in rows KL+1 to
c          2*KL+KU+1; rows 1 to KL of the array need not be set.
c          The j-th column of A is stored in the j-th column of the
c          array AB as follows:
c          AB(kl+ku+1+i-j,j) = A(i,j) for max(1,j-ku)<=i<=min(m,j+kl)
c
c          On exit, details of the factorization: U is stored as an
c          upper triangular band matrix with KL+KU superdiagonals in
c          rows 1 to KL+KU+1, and the multipliers used during the
c          factorization are stored in rows KL+KU+2 to 2*KL+KU+1.
c          See below for further details.
c
c  LDAB    (input) INTEGER
c          The leading dimension of the array AB.  LDAB >= 2*KL+KU+1.
c
c  IPIV    (output) INTEGER array, dimension (min(M,N))
c          The pivot indices; for 1 <= i <= min(M,N), row i of the
c          matrix was interchanged with row IPIV(i).
c
c  INFO    (output) INTEGER
c          = 0: successful exit
c          < 0: if INFO = -i, the i-th argument had an illegal value
c          > 0: if INFO = +i, U(i,i) is exactly zero. The factorization
c               has been completed, but the factor U is exactly
c               singular, and division by zero will occur if it is used
c               to solve a system of equations.
c
c  Further Details
c  ===============
c
c  The band storage scheme is illustrated by the following example, when
c  M = N = 6, KL = 2, KU = 1:
c
c  On entry:                       On exit:
c
c      *    *    *    +    +    +       *    *    *   u14  u25  u36
c      *    *    +    +    +    +       *    *   u13  u24  u35  u46
c      *   a12  a23  a34  a45  a56      *   u12  u23  u34  u45  u56
c     a11  a22  a33  a44  a55  a66     u11  u22  u33  u44  u55  u66
c     a21  a32  a43  a54  a65   *      m21  m32  m43  m54  m65   *
c     a31  a42  a53  a64   *    *      m31  m42  m53  m64   *    *
c
c  Array elements marked * are not used by the routine; elements marked
c  + need not be set on entry, but are required by the routine to store
c  elements of U because of fill-in resulting from the row interchanges.
c
c  =====================================================================
c
c     .. Parameters ..
      REAL               ONE, ZERO
      PARAMETER          ( ONE = 1.0E+0, ZERO = 0.0E+0 )
      INTEGER            NBMAX, LDWORK
      PARAMETER          ( NBMAX = 64, LDWORK = NBMAX+1 )
c     ..
c     .. Local Scalars ..
      INTEGER            I, I2, I3, II, IP, J, J2, J3, JB, JJ, JM, JP,
     &                   JU, K2, KM, KV, NB, NW
      REAL               TEMP
c     ..
c     .. Local Arrays ..
      REAL               WORK13( LDWORK, NBMAX ),
     &                   WORK31( LDWORK, NBMAX )
c     ..
c     .. External Functions ..
      INTEGER            ILAENV, ISAMAX
      EXTERNAL           ILAENV, ISAMAX
c     ..
c     .. External Subroutines ..
      EXTERNAL           SCOPY, SGBTF2, SGEMM, SGER, SLASWP, SSCAL,
     &                   SSWAP, STRSM, XERBLA
c     ..
c     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
c     ..


c     KV is the number of superdiagonals in the factor U, allowing for
c     fill-in

      KV = KU + KL

c     Test the input parameters.

      INFO = 0
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( KL.LT.0 ) THEN
         INFO = -3
      ELSE IF( KU.LT.0 ) THEN
         INFO = -4
      ELSE IF( LDAB.LT.KL+KV+1 ) THEN
         INFO = -6
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'SGBTRF', -INFO )
         RETURN
      END IF

c     Quick return if possible

      IF( M.EQ.0 .OR. N.EQ.0 )  RETURN

c     Determine the block size for this environment

      NB = ILAENV( 1, 'SGBTRF', ' ', M, N, KL, KU )

c     The block size must not exceed the limit set by the size of the
c     local arrays WORK13 and WORK31.

      NB = MIN( NB, NBMAX )

      IF( NB.LE.1 .OR. NB.GT.KL ) THEN

c        Use unblocked code

         CALL SGBTF2( M, N, KL, KU, AB, LDAB, IPIV, INFO )
      ELSE

c        Use blocked code

c        Zero the superdiagonal elements of the work array WORK13

         DO 20 J = 1, NB
            DO 10 I = 1, J - 1
               WORK13( I, J ) = ZERO
   10       CONTINUE
   20    CONTINUE

c        Zero the subdiagonal elements of the work array WORK31

         DO 40 J = 1, NB
            DO 30 I = J + 1, NB
               WORK31( I, J ) = ZERO
   30       CONTINUE
   40    CONTINUE

c        Gaussian elimination with partial pivoting

c        Set fill-in elements in columns KU+2 to KV to zero

         DO 60 J = KU + 2, MIN( KV, N )
            DO 50 I = KV - J + 2, KL
               AB( I, J ) = ZERO
   50       CONTINUE
   60    CONTINUE

c        JU is the index of the last column affected by the current
c        stage of the factorization

         JU = 1

         DO 180 J = 1, MIN( M, N ), NB
            JB = MIN( NB, MIN( M, N )-J+1 )

c           The active part of the matrix is partitioned
c
c              A11   A12   A13
c              A21   A22   A23
c              A31   A32   A33
c
c           Here A11, A21 and A31 denote the current block of JB columns
c           which is about to be factorized. The number of rows in the
c           partitioning are JB, I2, I3 respectively, and the numbers
c           of columns are JB, J2, J3. The superdiagonal elements of A13
c           and the subdiagonal elements of A31 lie outside the band.

            I2 = MIN( KL-JB, M-J-JB+1 )
            I3 = MIN( JB, M-J-KL+1 )

c           J2 and J3 are computed after JU has been updated.

c           Factorize the current block of JB columns

            DO 80 JJ = J, J + JB - 1

c              Set fill-in elements in column JJ+KV to zero

               IF( JJ+KV.LE.N ) THEN
                  DO 70 I = 1, KL
                     AB( I, JJ+KV ) = ZERO
   70             CONTINUE
               END IF

c              Find pivot and test for singularity. KM is the number of
c              subdiagonal elements in the current column.

               KM = MIN( KL, M-JJ )
               JP = ISAMAX( KM+1, AB( KV+1, JJ ), 1 )
               IPIV( JJ ) = JP + JJ - J
               IF( AB( KV+JP, JJ ).NE.ZERO ) THEN
                  JU = MAX( JU, MIN( JJ+KU+JP-1, N ) )
                  IF( JP.NE.1 ) THEN

c                    Apply interchange to columns J to J+JB-1

                     IF( JP+JJ-1.LT.J+KL ) THEN

                        CALL SSWAP( JB, AB( KV+1+JJ-J, J ), LDAB-1,
     &                              AB( KV+JP+JJ-J, J ), LDAB-1 )
                     ELSE

c                       The interchange affects columns J to JJ-1 of A31
c                       which are stored in the work array WORK31

                        CALL SSWAP( JJ-J, AB( KV+1+JJ-J, J ), LDAB-1,
     &                              WORK31( JP+JJ-J-KL, 1 ), LDWORK )
                        CALL SSWAP( J+JB-JJ, AB( KV+1, JJ ), LDAB-1,
     &                              AB( KV+JP, JJ ), LDAB-1 )
                     END IF
                  END IF

c                 Compute multipliers

                  CALL SSCAL( KM, ONE / AB( KV+1, JJ ), AB( KV+2, JJ ),
     &                        1 )

c                 Update trailing submatrix within the band and within
c                 the current block. JM is the index of the last column
c                 which needs to be updated.

                  JM = MIN( JU, J+JB-1 )
                  IF( JM.GT.JJ )
     &               CALL SGER( KM, JM-JJ, -ONE, AB( KV+2, JJ ), 1,
     &                          AB( KV, JJ+1 ), LDAB-1,
     &                          AB( KV+1, JJ+1 ), LDAB-1 )
               ELSE

c                 If pivot is zero, set INFO to the index of the pivot
c                 unless a zero pivot has already been found.

                  IF( INFO.EQ.0 )
     &               INFO = JJ
               END IF

c              Copy current column of A31 into the work array WORK31

               NW = MIN( JJ-J+1, I3 )
               IF( NW.GT.0 )
     &            CALL SCOPY( NW, AB( KV+KL+1-JJ+J, JJ ), 1,
     &                        WORK31( 1, JJ-J+1 ), 1 )
   80       CONTINUE
            IF( J+JB.LE.N ) THEN

c              Apply the row interchanges to the other blocks.

               J2 = MIN( JU-J+1, KV ) - JB
               J3 = MAX( 0, JU-J-KV+1 )

c              Use SLASWP to apply the row interchanges to A12, A22, and
c              A32.

               CALL SLASWP( J2, AB( KV+1-JB, J+JB ), LDAB-1, 1, JB,
     &                      IPIV( J ), 1 )

c              Adjust the pivot indices.

               DO 90 I = J, J + JB - 1
                  IPIV( I ) = IPIV( I ) + J - 1
   90          CONTINUE

c              Apply the row interchanges to A13, A23, and A33
c              columnwise.

               K2 = J - 1 + JB + J2
               DO 110 I = 1, J3
                  JJ = K2 + I
                  DO 100 II = J + I - 1, J + JB - 1
                     IP = IPIV( II )
                     IF( IP.NE.II ) THEN
                        TEMP = AB( KV+1+II-JJ, JJ )
                        AB( KV+1+II-JJ, JJ ) = AB( KV+1+IP-JJ, JJ )
                        AB( KV+1+IP-JJ, JJ ) = TEMP
                     END IF
  100             CONTINUE
  110          CONTINUE

c              Update the relevant part of the trailing submatrix

               IF( J2.GT.0 ) THEN

c                 Update A12

                  CALL STRSM( 'Left', 'Lower', 'No transpose', 'Unit',
     &                        JB, J2, ONE, AB( KV+1, J ), LDAB-1,
     &                        AB( KV+1-JB, J+JB ), LDAB-1 )

                  IF( I2.GT.0 ) THEN

c                    Update A22

                     CALL SGEMM( 'No transpose', 'No transpose', I2, J2,
     &                           JB, -ONE, AB( KV+1+JB, J ), LDAB-1,
     &                           AB( KV+1-JB, J+JB ), LDAB-1, ONE,
     &                           AB( KV+1, J+JB ), LDAB-1 )
                  END IF

                  IF( I3.GT.0 ) THEN

c                    Update A32

                     CALL SGEMM( 'No transpose', 'No transpose', I3, J2,
     &                           JB, -ONE, WORK31, LDWORK,
     &                           AB( KV+1-JB, J+JB ), LDAB-1, ONE,
     &                           AB( KV+KL+1-JB, J+JB ), LDAB-1 )
                  END IF
               END IF

               IF( J3.GT.0 ) THEN

c                 Copy the lower triangle of A13 into the work array
c                 WORK13

                  DO 130 JJ = 1, J3
                     DO 120 II = JJ, JB
                        WORK13( II, JJ ) = AB( II-JJ+1, JJ+J+KV-1 )
  120                CONTINUE
  130             CONTINUE

c                 Update A13 in the work array

                  CALL STRSM( 'Left', 'Lower', 'No transpose', 'Unit',
     &                        JB, J3, ONE, AB( KV+1, J ), LDAB-1,
     &                        WORK13, LDWORK )

                  IF( I2.GT.0 ) THEN

c                    Update A23

                     CALL SGEMM( 'No transpose', 'No transpose', I2, J3,
     &                           JB, -ONE, AB( KV+1+JB, J ), LDAB-1,
     &                           WORK13, LDWORK, ONE, AB( 1+JB, J+KV ),
     &                           LDAB-1 )
                  END IF

                  IF( I3.GT.0 ) THEN

c                    Update A33

                     CALL SGEMM( 'No transpose', 'No transpose', I3, J3,
     &                           JB, -ONE, WORK31, LDWORK, WORK13,
     &                           LDWORK, ONE, AB( 1+KL, J+KV ), LDAB-1 )
                  END IF

c                 Copy the lower triangle of A13 back into place

                  DO 150 JJ = 1, J3
                     DO 140 II = JJ, JB
                        AB( II-JJ+1, JJ+J+KV-1 ) = WORK13( II, JJ )
  140                CONTINUE
  150             CONTINUE
               END IF
            ELSE

c              Adjust the pivot indices.

               DO 160 I = J, J + JB - 1
                  IPIV( I ) = IPIV( I ) + J - 1
  160          CONTINUE
            END IF

c           Partially undo the interchanges in the current block to
c           restore the upper triangular form of A31 and copy the upper
c           triangle of A31 back into place

            DO 170 JJ = J + JB - 1, J, -1
               JP = IPIV( JJ ) - JJ + 1
               IF( JP.NE.1 ) THEN

c                 Apply interchange to columns J to JJ-1

                  IF( JP+JJ-1.LT.J+KL ) THEN

c                    The interchange does not affect A31

                     CALL SSWAP( JJ-J, AB( KV+1+JJ-J, J ), LDAB-1,
     &                           AB( KV+JP+JJ-J, J ), LDAB-1 )
                  ELSE

c                    The interchange does affect A31

                     CALL SSWAP( JJ-J, AB( KV+1+JJ-J, J ), LDAB-1,
     &                           WORK31( JP+JJ-J-KL, 1 ), LDWORK )
                  END IF
               END IF

c              Copy the current column of A31 back into place

               NW = MIN( I3, JJ-J+1 )
               IF( NW.GT.0 )
     &            CALL SCOPY( NW, WORK31( 1, JJ-J+1 ), 1,
     &                        AB( KV+KL+1-JJ+J, JJ ), 1 )
  170       CONTINUE
  180    CONTINUE
      END IF

      RETURN

c            End of SGBTRF
      END


      SUBROUTINE SGBTRS( TRANS, N, KL, KU, NRHS, AB, LDAB, IPIV, B, LDB,
     &                   INFO )
c
c     March 31, 1993 
c
c     .. Scalar Arguments ..
      CHARACTER          TRANS
      INTEGER            INFO, KL, KU, LDAB, LDB, N, NRHS
c     ..
c     .. Array Arguments ..
      INTEGER            IPIV( * )
      REAL               AB( LDAB, * ), B( LDB, * )
c     ..
c
c  Purpose
c  =======
c
c  SGBTRS solves a system of linear equations
c     A * X = B  or  A' * X = B
c  with a general band matrix A using the LU factorization computed
c  by SGBTRF.
c
c  Arguments
c  =========
c
c  TRANS   (input) CHARACTER*1
c          Specifies the form of the system of equations.
c          = 'N':  A * X = B  (No transpose)
c          = 'T':  A'* X = B  (Transpose)
c          = 'C':  A'* X = B  (Conjugate transpose = Transpose)
c
c  N       (input) INTEGER
c          The order of the matrix A.  N >= 0.
c
c  KL      (input) INTEGER
c          The number of subdiagonals within the band of A.  KL >= 0.
c
c  KU      (input) INTEGER
c          The number of superdiagonals within the band of A.  KU >= 0.
c
c  NRHS    (input) INTEGER
c          The number of right hand sides, i.e., the number of columns
c          of the matrix B.  NRHS >= 0.
c
c  AB      (input) REAL array, dimension (LDAB,N)
c          Details of the LU factorization of the band matrix A, as
c          computed by SGBTRF.  U is stored as an upper triangular band
c          matrix with KL+KU superdiagonals in rows 1 to KL+KU+1, and
c          the multipliers used during the factorization are stored in
c          rows KL+KU+2 to 2*KL+KU+1.
c
c  LDAB    (input) INTEGER
c          The leading dimension of the array AB.  LDAB >= 2*KL+KU+1.
c
c  IPIV    (input) INTEGER array, dimension (N)
c          The pivot indices; for 1 <= i <= N, row i of the matrix was
c          interchanged with row IPIV(i).
c
c  B       (input/output) REAL array, dimension (LDB,NRHS)
c          On entry, the right hand side matrix B.
c          On exit, the solution matrix X.
c
c  LDB     (input) INTEGER
c          The leading dimension of the array B.  LDB >= max(1,N).
c
c  INFO    (output) INTEGER
c          = 0:  successful exit
c          < 0: if INFO = -i, the i-th argument had an illegal value
c
c  =====================================================================
c
c     .. Parameters ..
      REAL               ONE
      PARAMETER          ( ONE = 1.0E+0 )
c     ..
c     .. Local Scalars ..
      LOGICAL            LNOTI, NOTRAN
      INTEGER            I, J, KD, L, LM
c     ..
c     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
c     ..
c     .. External Subroutines ..
      EXTERNAL           SGEMV, SGER, SSWAP, STBSV, XERBLA
c     ..
c     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
c     ..


c     Test the input parameters.

      INFO = 0
      NOTRAN = LSAME( TRANS, 'N' )
      IF( .NOT.NOTRAN .AND. .NOT.LSAME( TRANS, 'T' ) .AND. .NOT.
     &    LSAME( TRANS, 'C' ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( KL.LT.0 ) THEN
         INFO = -3
      ELSE IF( KU.LT.0 ) THEN
         INFO = -4
      ELSE IF( NRHS.LT.0 ) THEN
         INFO = -5
      ELSE IF( LDAB.LT.( 2*KL+KU+1 ) ) THEN
         INFO = -7
      ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
         INFO = -10
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'SGBTRS', -INFO )
         RETURN
      END IF

c     Quick return if possible

      IF( N.EQ.0 .OR. NRHS.EQ.0 )  RETURN

      KD = KU + KL + 1
      LNOTI = KL.GT.0

      IF( NOTRAN ) THEN

c        Solve  A*X = B.
c
c        Solve L*X = B, overwriting B with X.
c
c        L is represented as a product of permutations and unit lower
c        triangular matrices L = P(1) * L(1) * ... * P(n-1) * L(n-1),
c        where each transformation L(i) is a rank-one modification of
c        the identity matrix.

         IF( LNOTI ) THEN
            DO 10 J = 1, N - 1
               LM = MIN( KL, N-J )
               L = IPIV( J )
               IF( L.NE.J )
     &            CALL SSWAP( NRHS, B( L, 1 ), LDB, B( J, 1 ), LDB )
               CALL SGER( LM, NRHS, -ONE, AB( KD+1, J ), 1, B( J, 1 ),
     &                    LDB, B( J+1, 1 ), LDB )
   10       CONTINUE
         END IF

         DO 20 I = 1, NRHS

c           Solve U*X = B, overwriting B with X.

            CALL STBSV( 'Upper', 'No transpose', 'Non-unit', N, KL+KU,
     &                  AB, LDAB, B( 1, I ), 1 )
   20    CONTINUE

      ELSE

c        Solve A'*X = B.

         DO 30 I = 1, NRHS

c           Solve U'*X = B, overwriting B with X.

            CALL STBSV( 'Upper', 'Transpose', 'Non-unit', N, KL+KU, AB,
     &                  LDAB, B( 1, I ), 1 )
   30    CONTINUE

c        Solve L'*X = B, overwriting B with X.

         IF( LNOTI ) THEN
            DO 40 J = N - 1, 1, -1
               LM = MIN( KL, N-J )
               CALL SGEMV( 'Transpose', LM, NRHS, -ONE, B( J+1, 1 ),
     &                     LDB, AB( KD+1, J ), 1, ONE, B( J, 1 ), LDB )
               L = IPIV( J )
               IF( L.NE.J )
     &            CALL SSWAP( NRHS, B( L, 1 ), LDB, B( J, 1 ), LDB )
   40       CONTINUE
         END IF
      END IF
      RETURN

c           End of SGBTRS
      END


      SUBROUTINE SGBTF2( M, N, KL, KU, AB, LDAB, IPIV, INFO )
c
c     February 29, 1992
c
c     .. Scalar Arguments ..
      INTEGER            INFO, KL, KU, LDAB, M, N
c     ..
c     .. Array Arguments ..
      INTEGER            IPIV( * )
      REAL               AB( LDAB, * )
c     ..
c
c  Purpose
c  =======
c
c  SGBTF2 computes an LU factorization of a real m-by-n band matrix A
c  using partial pivoting with row interchanges.
c
c  This is the unblocked version of the algorithm, calling Level 2 BLAS.
c
c  Arguments
c  =========
c
c  M       (input) INTEGER
c          The number of rows of the matrix A.  M >= 0.
c
c  N       (input) INTEGER
c          The number of columns of the matrix A.  N >= 0.
c
c  KL      (input) INTEGER
c          The number of subdiagonals within the band of A.  KL >= 0.
c
c  KU      (input) INTEGER
c          The number of superdiagonals within the band of A.  KU >= 0.
c
c  AB      (input/output) REAL array, dimension (LDAB,N)
c          On entry, the matrix A in band storage, in rows KL+1 to
c          2*KL+KU+1; rows 1 to KL of the array need not be set.
c          The j-th column of A is stored in the j-th column of the
c          array AB as follows:
c          AB(kl+ku+1+i-j,j) = A(i,j) for max(1,j-ku)<=i<=min(m,j+kl)
c
c          On exit, details of the factorization: U is stored as an
c          upper triangular band matrix with KL+KU superdiagonals in
c          rows 1 to KL+KU+1, and the multipliers used during the
c          factorization are stored in rows KL+KU+2 to 2*KL+KU+1.
c          See below for further details.
c
c  LDAB    (input) INTEGER
c          The leading dimension of the array AB.  LDAB >= 2*KL+KU+1.
c
c  IPIV    (output) INTEGER array, dimension (min(M,N))
c          The pivot indices; for 1 <= i <= min(M,N), row i of the
c          matrix was interchanged with row IPIV(i).
c
c  INFO    (output) INTEGER
c          = 0: successful exit
c          < 0: if INFO = -i, the i-th argument had an illegal value
c          > 0: if INFO = +i, U(i,i) is exactly zero. The factorization
c               has been completed, but the factor U is exactly
c               singular, and division by zero will occur if it is used
c               to solve a system of equations.
c
c  Further Details
c  ===============
c
c  The band storage scheme is illustrated by the following example, when
c  M = N = 6, KL = 2, KU = 1:
c
c  On entry:                       On exit:
c
c      *    *    *    +    +    +       *    *    *   u14  u25  u36
c      *    *    +    +    +    +       *    *   u13  u24  u35  u46
c      *   a12  a23  a34  a45  a56      *   u12  u23  u34  u45  u56
c     a11  a22  a33  a44  a55  a66     u11  u22  u33  u44  u55  u66
c     a21  a32  a43  a54  a65   *      m21  m32  m43  m54  m65   *
c     a31  a42  a53  a64   *    *      m31  m42  m53  m64   *    *
c
c  Array elements marked * are not used by the routine; elements marked
c  + need not be set on entry, but are required by the routine to store
c  elements of U, because of fill-in resulting from the row
c  interchanges.
c
c  =====================================================================
c
c     .. Parameters ..
      REAL               ONE, ZERO
      PARAMETER          ( ONE = 1.0E+0, ZERO = 0.0E+0 )
c     ..
c     .. Local Scalars ..
      INTEGER            I, J, JP, JU, KM, KV
c     ..
c     .. External Functions ..
      INTEGER            ISAMAX
      EXTERNAL           ISAMAX
c     ..
c     .. External Subroutines ..
      EXTERNAL           SGER, SSCAL, SSWAP, XERBLA
c     ..
c     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
c     ..


c     KV is the number of superdiagonals in the factor U, allowing for
c     fill-in.

      KV = KU + KL

c     Test the input parameters.

      INFO = 0
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( KL.LT.0 ) THEN
         INFO = -3
      ELSE IF( KU.LT.0 ) THEN
         INFO = -4
      ELSE IF( LDAB.LT.KL+KV+1 ) THEN
         INFO = -6
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'SGBTF2', -INFO )
         RETURN
      END IF

c     Quick return if possible

      IF( M.EQ.0 .OR. N.EQ.0 )  RETURN

c     Gaussian elimination with partial pivoting

c     Set fill-in elements in columns KU+2 to KV to zero.

      DO 20 J = KU + 2, MIN( KV, N )
         DO 10 I = KV - J + 2, KL
            AB( I, J ) = ZERO
   10    CONTINUE
   20 CONTINUE

c     JU is the index of the last column affected by the current stage
c     of the factorization.

      JU = 1

      DO 40 J = 1, MIN( M, N )

c        Set fill-in elements in column J+KV to zero.

         IF( J+KV.LE.N ) THEN
            DO 30 I = 1, KL
               AB( I, J+KV ) = ZERO
   30       CONTINUE
         END IF

c        Find pivot and test for singularity. KM is the number of
c        subdiagonal elements in the current column.

         KM = MIN( KL, M-J )
         JP = ISAMAX( KM+1, AB( KV+1, J ), 1 )
         IPIV( J ) = JP + J - 1
         IF( AB( KV+JP, J ).NE.ZERO ) THEN
            JU = MAX( JU, MIN( J+KU+JP-1, N ) )

c           Apply interchange to columns J to JU.

            IF( JP.NE.1 )
     &         CALL SSWAP( JU-J+1, AB( KV+JP, J ), LDAB-1,
     &                     AB( KV+1, J ), LDAB-1 )

            IF( KM.GT.0 ) THEN

c              Compute multipliers.

               CALL SSCAL( KM, ONE / AB( KV+1, J ), AB( KV+2, J ), 1 )

c              Update trailing submatrix within the band.

               IF( JU.GT.J )
     &            CALL SGER( KM, JU-J, -ONE, AB( KV+2, J ), 1,
     &                       AB( KV, J+1 ), LDAB-1, AB( KV+1, J+1 ),
     &                       LDAB-1 )
            END IF
         ELSE

c           If pivot is zero, set INFO to the index of the pivot
c           unless a zero pivot has already been found.

            IF( INFO.EQ.0 )
     &         INFO = J
         END IF
   40 CONTINUE
      RETURN

c           End of SGBTF2
      END


      SUBROUTINE SGESV( N, NRHS, A, LDA, IPIV, B, LDB, INFO )
c
c  -- LAPACK driver routine (version 2.0) --
c     March 31, 1993 
c
c     .. Scalar Arguments ..
      INTEGER            INFO, LDA, LDB, N, NRHS
c     ..
c     .. Array Arguments ..
      INTEGER            IPIV( * )
      REAL               A( LDA, * ), B( LDB, * )
c     ..
c
c  Purpose
c  =======
c
c  SGESV computes the solution to a real system of linear equations
c     A * X = B,
c  where A is an N-by-N matrix and X and B are N-by-NRHS matrices.
c
c  The LU decomposition with partial pivoting and row interchanges is
c  used to factor A as
c     A = P * L * U,
c  where P is a permutation matrix, L is unit lower triangular, and U is
c  upper triangular.  The factored form of A is then used to solve the
c  system of equations A * X = B.
c
c  Arguments
c  =========
c
c  N       (input) INTEGER
c          The number of linear equations, i.e., the order of the
c          matrix A.  N >= 0.
c
c  NRHS    (input) INTEGER
c          The number of right hand sides, i.e., the number of columns
c          of the matrix B.  NRHS >= 0.
c
c  A       (input/output) REAL array, dimension (LDA,N)
c          On entry, the N-by-N coefficient matrix A.
c          On exit, the factors L and U from the factorization
c          A = P*L*U; the unit diagonal elements of L are not stored.
c
c  LDA     (input) INTEGER
c          The leading dimension of the array A.  LDA >= max(1,N).
c
c  IPIV    (output) INTEGER array, dimension (N)
c          The pivot indices that define the permutation matrix P;
c          row i of the matrix was interchanged with row IPIV(i).
c
c  B       (input/output) REAL array, dimension (LDB,NRHS)
c          On entry, the N-by-NRHS matrix of right hand side matrix B.
c          On exit, if INFO = 0, the N-by-NRHS solution matrix X.
c
c  LDB     (input) INTEGER
c          The leading dimension of the array B.  LDB >= max(1,N).
c
c  INFO    (output) INTEGER
c          = 0:  successful exit
c          < 0:  if INFO = -i, the i-th argument had an illegal value
c          > 0:  if INFO = i, U(i,i) is exactly zero.  The factorization
c                has been completed, but the factor U is exactly
c                singular, so the solution could not be computed.
c
c  =====================================================================
c
c     .. External Subroutines ..
      EXTERNAL           SGETRF, SGETRS, XERBLA
c     ..
c     .. Intrinsic Functions ..
      INTRINSIC          MAX
c     ..


c     Test the input parameters.

      INFO = 0
      IF( N.LT.0 ) THEN
         INFO = -1
      ELSE IF( NRHS.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -4
      ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
         INFO = -7
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'SGESV ', -INFO )
         RETURN
      END IF

c     Compute the LU factorization of A.

      CALL SGETRF( N, N, A, LDA, IPIV, INFO )
      IF( INFO.EQ.0 ) THEN

c        Solve the system A*X = B, overwriting B with X.

         CALL SGETRS( 'No transpose', N, NRHS, A, LDA, IPIV, B, LDB,
     &                INFO )
      END IF
      RETURN

c           End of SGESV
      END


      SUBROUTINE SGETRF( M, N, A, LDA, IPIV, INFO )
c
c     March 31, 1993 
c
c     .. Scalar Arguments ..
      INTEGER            INFO, LDA, M, N
c     ..
c     .. Array Arguments ..
      INTEGER            IPIV( * )
      REAL               A( LDA, * )
c     ..
c
c  Purpose
c  =======
c
c  SGETRF computes an LU factorization of a general M-by-N matrix A
c  using partial pivoting with row interchanges.
c
c  The factorization has the form
c     A = P * L * U
c  where P is a permutation matrix, L is lower triangular with unit
c  diagonal elements (lower trapezoidal if m > n), and U is upper
c  triangular (upper trapezoidal if m < n).
c
c  This is the right-looking Level 3 BLAS version of the algorithm.
c
c  Arguments
c  =========
c
c  M       (input) INTEGER
c          The number of rows of the matrix A.  M >= 0.
c
c  N       (input) INTEGER
c          The number of columns of the matrix A.  N >= 0.
c
c  A       (input/output) REAL array, dimension (LDA,N)
c          On entry, the M-by-N matrix to be factored.
c          On exit, the factors L and U from the factorization
c          A = P*L*U; the unit diagonal elements of L are not stored.
c
c  LDA     (input) INTEGER
c          The leading dimension of the array A.  LDA >= max(1,M).
c
c  IPIV    (output) INTEGER array, dimension (min(M,N))
c          The pivot indices; for 1 <= i <= min(M,N), row i of the
c          matrix was interchanged with row IPIV(i).
c
c  INFO    (output) INTEGER
c          = 0:  successful exit
c          < 0:  if INFO = -i, the i-th argument had an illegal value
c          > 0:  if INFO = i, U(i,i) is exactly zero. The factorization
c                has been completed, but the factor U is exactly
c                singular, and division by zero will occur if it is used
c                to solve a system of equations.
c
c  =====================================================================
c
c     .. Parameters ..
      REAL               ONE
      PARAMETER          ( ONE = 1.0E+0 )
c     ..
c     .. Local Scalars ..
      INTEGER            I, IINFO, J, JB, NB
c     ..
c     .. External Subroutines ..
      EXTERNAL           SGEMM, SGETF2, SLASWP, STRSM, XERBLA
c     ..
c     .. External Functions ..
      INTEGER            ILAENV
      EXTERNAL           ILAENV
c     ..
c     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
c     ..


c     Test the input parameters.

      INFO = 0
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -4
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'SGETRF', -INFO )
         RETURN
      END IF

c     Quick return if possible

      IF( M.EQ.0 .OR. N.EQ.0 )  RETURN

c     Determine the block size for this environment.

      NB = ILAENV( 1, 'SGETRF', ' ', M, N, -1, -1 )
      IF( NB.LE.1 .OR. NB.GE.MIN( M, N ) ) THEN

c        Use unblocked code.

         CALL SGETF2( M, N, A, LDA, IPIV, INFO )
      ELSE

c        Use blocked code.

         DO 20 J = 1, MIN( M, N ), NB
            JB = MIN( MIN( M, N )-J+1, NB )

c           Factor diagonal and subdiagonal blocks and test for exact
c           singularity.

            CALL SGETF2( M-J+1, JB, A( J, J ), LDA, IPIV( J ), IINFO )

c           Adjust INFO and the pivot indices.

            IF( INFO.EQ.0 .AND. IINFO.GT.0 )
     &         INFO = IINFO + J - 1
            DO 10 I = J, MIN( M, J+JB-1 )
               IPIV( I ) = J - 1 + IPIV( I )
   10       CONTINUE

c           Apply interchanges to columns 1:J-1.

            CALL SLASWP( J-1, A, LDA, J, J+JB-1, IPIV, 1 )

            IF( J+JB.LE.N ) THEN

c              Apply interchanges to columns J+JB:N.

               CALL SLASWP( N-J-JB+1, A( 1, J+JB ), LDA, J, J+JB-1,
     &                      IPIV, 1 )

c              Compute block row of U.

               CALL STRSM( 'Left', 'Lower', 'No transpose', 'Unit', JB,
     &                     N-J-JB+1, ONE, A( J, J ), LDA, A( J, J+JB ),
     &                     LDA )
               IF( J+JB.LE.M ) THEN

c                 Update trailing submatrix.

                  CALL SGEMM( 'No transpose', 'No transpose', M-J-JB+1,
     &                        N-J-JB+1, JB, -ONE, A( J+JB, J ), LDA,
     &                        A( J, J+JB ), LDA, ONE, A( J+JB, J+JB ),
     &                        LDA )
               END IF
            END IF
   20    CONTINUE
      END IF
      RETURN

c           End of SGETRF
      END


      SUBROUTINE SGETRS( TRANS, N, NRHS, A, LDA, IPIV, B, LDB, INFO )
c
c     March 31, 1993 
c
c     .. Scalar Arguments ..
      CHARACTER          TRANS
      INTEGER            INFO, LDA, LDB, N, NRHS
c     ..
c     .. Array Arguments ..
      INTEGER            IPIV( * )
      REAL               A( LDA, * ), B( LDB, * )
c     ..
c
c  Purpose
c  =======
c
c  SGETRS solves a system of linear equations
c     A * X = B  or  A' * X = B
c  with a general N-by-N matrix A using the LU factorization computed
c  by SGETRF.
c
c  Arguments
c  =========
c
c  TRANS   (input) CHARACTER*1
c          Specifies the form of the system of equations:
c          = 'N':  A * X = B  (No transpose)
c          = 'T':  A'* X = B  (Transpose)
c          = 'C':  A'* X = B  (Conjugate transpose = Transpose)
c
c  N       (input) INTEGER
c          The order of the matrix A.  N >= 0.
c
c  NRHS    (input) INTEGER
c          The number of right hand sides, i.e., the number of columns
c          of the matrix B.  NRHS >= 0.
c
c  A       (input) REAL array, dimension (LDA,N)
c          The factors L and U from the factorization A = P*L*U
c          as computed by SGETRF.
c
c  LDA     (input) INTEGER
c          The leading dimension of the array A.  LDA >= max(1,N).
c
c  IPIV    (input) INTEGER array, dimension (N)
c          The pivot indices from SGETRF; for 1<=i<=N, row i of the
c          matrix was interchanged with row IPIV(i).
c
c  B       (input/output) REAL array, dimension (LDB,NRHS)
c          On entry, the right hand side matrix B.
c          On exit, the solution matrix X.
c
c  LDB     (input) INTEGER
c          The leading dimension of the array B.  LDB >= max(1,N).
c
c  INFO    (output) INTEGER
c          = 0:  successful exit
c          < 0:  if INFO = -i, the i-th argument had an illegal value
c
c  =====================================================================
c
c     .. Parameters ..
      REAL               ONE
      PARAMETER          ( ONE = 1.0E+0 )
c     ..
c     .. Local Scalars ..
      LOGICAL            NOTRAN
c     ..
c     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
c     ..
c     .. External Subroutines ..
      EXTERNAL           SLASWP, STRSM, XERBLA
c     ..
c     .. Intrinsic Functions ..
      INTRINSIC          MAX
c     ..


c     Test the input parameters.

      INFO = 0
      NOTRAN = LSAME( TRANS, 'N' )
      IF( .NOT.NOTRAN .AND. .NOT.LSAME( TRANS, 'T' ) .AND. .NOT.
     &    LSAME( TRANS, 'C' ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( NRHS.LT.0 ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -5
      ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
         INFO = -8
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'SGETRS', -INFO )
         RETURN
      END IF

c     Quick return if possible

      IF( N.EQ.0 .OR. NRHS.EQ.0 )  RETURN

      IF( NOTRAN ) THEN

c        Solve A * X = B.
c
c        Apply row interchanges to the right hand sides.

         CALL SLASWP( NRHS, B, LDB, 1, N, IPIV, 1 )

c        Solve L*X = B, overwriting B with X.

         CALL STRSM( 'Left', 'Lower', 'No transpose', 'Unit', N, NRHS,
     &               ONE, A, LDA, B, LDB )

c        Solve U*X = B, overwriting B with X.

         CALL STRSM( 'Left', 'Upper', 'No transpose', 'Non-unit', N,
     &               NRHS, ONE, A, LDA, B, LDB )
      ELSE

c        Solve A' * X = B.
c
c        Solve U'*X = B, overwriting B with X.

         CALL STRSM( 'Left', 'Upper', 'Transpose', 'Non-unit', N, NRHS,
     &               ONE, A, LDA, B, LDB )

c        Solve L'*X = B, overwriting B with X.

         CALL STRSM( 'Left', 'Lower', 'Transpose', 'Unit', N, NRHS, ONE,
     &               A, LDA, B, LDB )

c        Apply row interchanges to the solution vectors.

         CALL SLASWP( NRHS, B, LDB, 1, N, IPIV, -1 )
      END IF

      RETURN

c           End of SGETRS
      END


      SUBROUTINE SGETF2( M, N, A, LDA, IPIV, INFO )
c
c     June 30, 1992
c
c     .. Scalar Arguments ..
      INTEGER            INFO, LDA, M, N
c     ..
c     .. Array Arguments ..
      INTEGER            IPIV( * )
      REAL               A( LDA, * )
c     ..
c
c  Purpose
c  =======
c
c  SGETF2 computes an LU factorization of a general m-by-n matrix A
c  using partial pivoting with row interchanges.
c
c  The factorization has the form
c     A = P * L * U
c  where P is a permutation matrix, L is lower triangular with unit
c  diagonal elements (lower trapezoidal if m > n), and U is upper
c  triangular (upper trapezoidal if m < n).
c
c  This is the right-looking Level 2 BLAS version of the algorithm.
c
c  Arguments
c  =========
c
c  M       (input) INTEGER
c          The number of rows of the matrix A.  M >= 0.
c
c  N       (input) INTEGER
c          The number of columns of the matrix A.  N >= 0.
c
c  A       (input/output) REAL array, dimension (LDA,N)
c          On entry, the m by n matrix to be factored.
c          On exit, the factors L and U from the factorization
c          A = P*L*U; the unit diagonal elements of L are not stored.
c
c  LDA     (input) INTEGER
c          The leading dimension of the array A.  LDA >= max(1,M).
c
c  IPIV    (output) INTEGER array, dimension (min(M,N))
c          The pivot indices; for 1 <= i <= min(M,N), row i of the
c          matrix was interchanged with row IPIV(i).
c
c  INFO    (output) INTEGER
c          = 0: successful exit
c          < 0: if INFO = -k, the k-th argument had an illegal value
c          > 0: if INFO = k, U(k,k) is exactly zero. The factorization
c               has been completed, but the factor U is exactly
c               singular, and division by zero will occur if it is used
c               to solve a system of equations.
c
c  =====================================================================
c
c     .. Parameters ..
      REAL               ONE, ZERO
      PARAMETER          ( ONE = 1.0E+0, ZERO = 0.0E+0 )
c     ..
c     .. Local Scalars ..
      INTEGER            J, JP
c     ..
c     .. External Functions ..
      INTEGER            ISAMAX
      EXTERNAL           ISAMAX
c     ..
c     .. External Subroutines ..
      EXTERNAL           SGER, SSCAL, SSWAP, XERBLA
c     ..
c     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
c     ..


c     Test the input parameters.

      INFO = 0
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -4
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'SGETF2', -INFO )
         RETURN
      END IF

c     Quick return if possible

      IF( M.EQ.0 .OR. N.EQ.0 )  RETURN

      DO 10 J = 1, MIN( M, N )

c        Find pivot and test for singularity.

         JP = J - 1 + ISAMAX( M-J+1, A( J, J ), 1 )
         IPIV( J ) = JP
         IF( A( JP, J ).NE.ZERO ) THEN

c           Apply the interchange to columns 1:N.

            IF( JP.NE.J )
     &         CALL SSWAP( N, A( J, 1 ), LDA, A( JP, 1 ), LDA )

c           Compute elements J+1:M of J-th column.

            IF( J.LT.M )
     &         CALL SSCAL( M-J, ONE / A( J, J ), A( J+1, J ), 1 )

         ELSE IF( INFO.EQ.0 ) THEN

            INFO = J
         END IF

         IF( J.LT.MIN( M, N ) ) THEN

c           Update trailing submatrix.

            CALL SGER( M-J, N-J, -ONE, A( J+1, J ), 1, A( J, J+1 ), LDA,
     &                 A( J+1, J+1 ), LDA )
         END IF
   10 CONTINUE
      RETURN

c           End of SGETF2
      END


      SUBROUTINE SLASWP( N, A, LDA, K1, K2, IPIV, INCX )
c
c     October 31, 1992
c
c     .. Scalar Arguments ..
      INTEGER            INCX, K1, K2, LDA, N
c     ..
c     .. Array Arguments ..
      INTEGER            IPIV( * )
      REAL               A( LDA, * )
c     ..
c
c  Purpose
c  =======
c
c  SLASWP performs a series of row interchanges on the matrix A.
c  One row interchange is initiated for each of rows K1 through K2 of A.
c
c  Arguments
c  =========
c
c  N       (input) INTEGER
c          The number of columns of the matrix A.
c
c  A       (input/output) REAL array, dimension (LDA,N)
c          On entry, the matrix of column dimension N to which the row
c          interchanges will be applied.
c          On exit, the permuted matrix.
c
c  LDA     (input) INTEGER
c          The leading dimension of the array A.
c
c  K1      (input) INTEGER
c          The first element of IPIV for which a row interchange will
c          be done.
c
c  K2      (input) INTEGER
c          The last element of IPIV for which a row interchange will
c          be done.
c
c  IPIV    (input) INTEGER array, dimension (M*abs(INCX))
c          The vector of pivot indices.  Only the elements in positions
c          K1 through K2 of IPIV are accessed.
c          IPIV(K) = L implies rows K and L are to be interchanged.
c
c  INCX    (input) INTEGER
c          The increment between successive values of IPIV.  If IPIV
c          is negative, the pivots are applied in reverse order.
c
c =====================================================================
c
c     .. Local Scalars ..
      INTEGER            I, IP, IX
c     ..
c     .. External Subroutines ..
      EXTERNAL           SSWAP
c     ..


c     Interchange row I with row IPIV(I) for each of rows K1 through K2.

      IF( INCX.EQ.0 )
     &   RETURN
      IF( INCX.GT.0 ) THEN
         IX = K1
      ELSE
         IX = 1 + ( 1-K2 )*INCX
      END IF

      IF( INCX.EQ.1 ) THEN

         DO 10 I = K1, K2
            IP = IPIV( I )
            IF( IP.NE.I )
     &         CALL SSWAP( N, A( I, 1 ), LDA, A( IP, 1 ), LDA )
   10    CONTINUE

      ELSE IF( INCX.GT.1 ) THEN
         DO 20 I = K1, K2
            IP = IPIV( IX )
            IF( IP.NE.I )
     &         CALL SSWAP( N, A( I, 1 ), LDA, A( IP, 1 ), LDA )
            IX = IX + INCX
   20    CONTINUE

      ELSE IF( INCX.LT.0 ) THEN
         DO 30 I = K2, K1, -1
            IP = IPIV( IX )
            IF( IP.NE.I )
     &         CALL SSWAP( N, A( I, 1 ), LDA, A( IP, 1 ), LDA )
            IX = IX + INCX
   30    CONTINUE
      END IF

      RETURN

c           End of SLASWP
      END


      INTEGER FUNCTION ILAENV( ISPEC, NAME, OPTS, N1, N2, N3, N4 )
c
C       ** SPECIALIZED TO JUST THOSE CASES NEEDED BY DISORT V1.3
C       ** (ISPEC=1 and NAME= SGBTRF or SGETRF)
c       ** Only possible values are 1, 32, and 64 in this case.
c
c     September 30, 1994
c
c     .. Scalar Arguments ..
      CHARACTER*( * )    NAME, OPTS
      INTEGER            ISPEC, N1, N2, N3, N4
c     ..
c
c  Purpose
c  =======
c
c  ILAENV is called from the LAPACK routines to choose problem-dependent
c  parameters for the local environment.  See ISPEC for a description of
c  the parameters.
c
c  This version provides a set of parameters which should give good,
c  but not optimal, performance on many of the currently available
c  computers.  Users are encouraged to modify this subroutine to set
c  the tuning parameters for their particular machine using the option
c  and problem size information in the arguments.
c
c  This routine will not function correctly if it is converted to all
c  lower case.  Converting it to all upper case is allowed.
c
c  Arguments
c  =========
c
c  ISPEC   (input) INTEGER
c          Specifies the parameter to be returned as the value of
c          ILAENV.
c          = 1: the optimal blocksize; if this value is 1, an unblocked
c               algorithm will give the best performance.
c
c  NAME    (input) CHARACTER*(*)
c          The name of the calling subroutine, in either upper case or
c          lower case.
c
c  OPTS    (input) CHARACTER*(*)
c          The character options to the subroutine NAME, concatenated
c          into a single character string.  For example, UPLO = 'U',
c          TRANS = 'T', and DIAG = 'N' for a triangular routine would
c          be specified as OPTS = 'UTN'.
c
c  N1      (input) INTEGER
c  N2      (input) INTEGER
c  N3      (input) INTEGER
c  N4      (input) INTEGER
c          Problem dimensions for the subroutine NAME; these may not all
c          be required.
c
c (ILAENV) (output) INTEGER
c          >= 0: the value of the parameter specified by ISPEC
c          < 0:  ISPEC had an illegal value.
c
c  Further Details
c  ===============
c
c  The following conventions have been used when calling ILAENV from the
c  LAPACK routines:
c  1)  OPTS is a concatenation of all of the character options to
c      subroutine NAME, in the same order that they appear in the
c      argument list for NAME, even if they are not used in determining
c      the value of the parameter specified by ISPEC.
c  2)  The problem dimensions N1, N2, N3, N4 are specified in the order
c      that they appear in the argument list for NAME.  N1 is used
c      first, N2 second, and so on, and unused problem dimensions are
c      passed a value of -1.
c  3)  The parameter value returned by ILAENV is checked for validity in
c      the calling subroutine.  For example, ILAENV is used to retrieve
c      the optimal blocksize for STRTRI as follows:
c
c      NB = ILAENV( 1, 'STRTRI', UPLO // DIAG, N, -1, -1, -1 )
c      IF( NB.LE.1 ) NB = MAX( 1, N )
c
c  =====================================================================
c
c     .. Local Scalars ..
      CHARACTER*2        C2
      CHARACTER*3        C3
      CHARACTER*6        SUBNAM
      INTEGER            I, IC, IZ, NB
c     ..
c     .. Intrinsic Functions ..
      INTRINSIC          CHAR, ICHAR
c     ..

      IF( ISPEC.NE.1 ) THEN
c                              Invalid value for ISPEC
         ILAENV = -1
         RETURN
      END IF

c     Convert NAME to upper case if the first character is lower case.

      ILAENV = 1
      SUBNAM = NAME
      IC = ICHAR( SUBNAM( 1:1 ) )
      IZ = ICHAR( 'Z' )
      IF( IZ.EQ.90 .OR. IZ.EQ.122 ) THEN

c                          ASCII character set

         IF( IC.GE.97 .AND. IC.LE.122 ) THEN
            SUBNAM( 1:1 ) = CHAR( IC-32 )
            DO 10 I = 2, 6
               IC = ICHAR( SUBNAM( I:I ) )
               IF( IC.GE.97 .AND. IC.LE.122 )
     &            SUBNAM( I:I ) = CHAR( IC-32 )
   10       CONTINUE
         END IF

      ELSE IF( IZ.EQ.233 .OR. IZ.EQ.169 ) THEN

c                          EBCDIC character set

         IF( ( IC.GE.129 .AND. IC.LE.137 ) .OR.
     &       ( IC.GE.145 .AND. IC.LE.153 ) .OR.
     &       ( IC.GE.162 .AND. IC.LE.169 ) ) THEN
            SUBNAM( 1:1 ) = CHAR( IC+64 )
            DO 20 I = 2, 6
               IC = ICHAR( SUBNAM( I:I ) )
               IF( ( IC.GE.129 .AND. IC.LE.137 ) .OR.
     &             ( IC.GE.145 .AND. IC.LE.153 ) .OR.
     &             ( IC.GE.162 .AND. IC.LE.169 ) )
     &            SUBNAM( I:I ) = CHAR( IC+64 )
   20       CONTINUE
         END IF

      ELSE IF( IZ.EQ.218 .OR. IZ.EQ.250 ) THEN

c        Prime machines:  ASCII+128

         IF( IC.GE.225 .AND. IC.LE.250 ) THEN
            SUBNAM( 1:1 ) = CHAR( IC-32 )
            DO 30 I = 2, 6
               IC = ICHAR( SUBNAM( I:I ) )
               IF( IC.GE.225 .AND. IC.LE.250 )
     &            SUBNAM( I:I ) = CHAR( IC-32 )
   30       CONTINUE
         END IF
      END IF

      C2 = SUBNAM( 2:3 )
      C3 = SUBNAM( 4:6 )

      NB = 1

      IF( C2.EQ.'GE' ) THEN
      
         IF( C3.EQ.'TRF' ) THEN
               NB = 64
         END IF
         
      ELSE IF( C2.EQ.'GB' ) THEN
      
         IF( C3.EQ.'TRF' ) THEN
               IF( N4.LE.64 ) THEN
                  NB = 1
               ELSE
                  NB = 32
               END IF
         END IF
         
      END IF
      
      ILAENV = NB

      RETURN

c          End of ILAENV
      END

c ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c
c     BLAS routines needed by LAPACK routines used in DISORT v1.3
c     (current as of Nov 99)
c     (level 3 routines first, then level 2, then level 1)
c
c NOTES:  
c
c (1) Try hard to use or implement local libraries of these BLAS
c     routines rather than these Fortran versions; Fortran versions
c     will slow LAPACK down considerably; for information on 
c     available optimized BLAS libraries for a variety of
c     computers, refer to:
c
c        http://www.netlib.org/blas/faq.html
c
c (2) The BLAS routines are available from 
c
c        http://www.netlib.org/blas
c        http://netlib.bell-labs.com/netlib/master/readme.html
c
c ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c
c The following authorship attribution was removed from the beginning
c of each BLAS-2,3 routine:
c     Jack Dongarra, Argonne National Laboratory.
c     Iain Duff, AERE Harwell.
c     Jeremy Du Croz, Numerical Algorithms Group Ltd.
c     Sven Hammarling, Numerical Algorithms Group Ltd.
c     Richard Hanson, Sandia National Labs.
c
c Other cosmetic edits were made, so these will not 'diff' perfectly
c with the original BLAS routines.
c
c ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c
c Routine list + call tree (omitting calls to error msg routine XERBLA):
c
c    SGEMM
c        LSAME
c    STRSM
c        LSAME
c    SGEMV
c        LSAME
c    STBSV
c        LSAME
c    SGER
c    ISAMAX
c    SCOPY
c    SSCAL
c    SSWAP
c    XERBLA
c ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


      SUBROUTINE SGEMM( TRANSA, TRANSB, M, N, K, ALPHA, A, LDA, B, LDB,
     &                   BETA, C, LDC )

c  Level 3 Blas routine.
c
c  -- Written on 8-February-1989.

c     .. Scalar Arguments ..
      CHARACTER*1        TRANSA, TRANSB
      INTEGER            M, N, K, LDA, LDB, LDC
      REAL               ALPHA, BETA
c     .. Array Arguments ..
      REAL               A( LDA, * ), B( LDB, * ), C( LDC, * )
c     ..
c
c  Purpose
c  =======
c
c  SGEMM  performs one of the matrix-matrix operations
c
c     C := alpha*op( A )*op( B ) + beta*C,
c
c  where  op( X ) is one of
c
c     op( X ) = X   or   op( X ) = X',
c
c  alpha and beta are scalars, and A, B and C are matrices, with op( A )
c  an m by k matrix,  op( B )  a  k by n matrix and  C an m by n matrix.
c
c  Parameters
c  ==========
c
c  TRANSA - CHARACTER*1.
c           On entry, TRANSA specifies the form of op( A ) to be used in
c           the matrix multiplication as follows:
c
c              TRANSA = 'N' or 'n',  op( A ) = A.
c
c              TRANSA = 'T' or 't',  op( A ) = A'.
c
c              TRANSA = 'C' or 'c',  op( A ) = A'.
c
c           Unchanged on exit.
c
c  TRANSB - CHARACTER*1.
c           On entry, TRANSB specifies the form of op( B ) to be used in
c           the matrix multiplication as follows:
c
c              TRANSB = 'N' or 'n',  op( B ) = B.
c
c              TRANSB = 'T' or 't',  op( B ) = B'.
c
c              TRANSB = 'C' or 'c',  op( B ) = B'.
c
c           Unchanged on exit.
c
c  M      - INTEGER.
c           On entry,  M  specifies  the number  of rows  of the  matrix
c           op( A )  and of the  matrix  C.  M  must  be at least  zero.
c           Unchanged on exit.
c
c  N      - INTEGER.
c           On entry,  N  specifies the number  of columns of the matrix
c           op( B ) and the number of columns of the matrix C. N must be
c           at least zero.
c           Unchanged on exit.
c
c  K      - INTEGER.
c           On entry,  K  specifies  the number of columns of the matrix
c           op( A ) and the number of rows of the matrix op( B ). K must
c           be at least  zero.
c           Unchanged on exit.
c
c  ALPHA  - REAL            .
c           On entry, ALPHA specifies the scalar alpha.
c           Unchanged on exit.
c
c  A      - REAL             array of DIMENSION ( LDA, ka ), where ka is
c           k  when  TRANSA = 'N' or 'n',  and is  m  otherwise.
c           Before entry with  TRANSA = 'N' or 'n',  the leading  m by k
c           part of the array  A  must contain the matrix  A,  otherwise
c           the leading  k by m  part of the array  A  must contain  the
c           matrix A.
c           Unchanged on exit.
c
c  LDA    - INTEGER.
c           On entry, LDA specifies the first dimension of A as declared
c           in the calling (sub) program. When  TRANSA = 'N' or 'n' then
c           LDA must be at least  max( 1, m ), otherwise  LDA must be at
c           least  max( 1, k ).
c           Unchanged on exit.
c
c  B      - REAL             array of DIMENSION ( LDB, kb ), where kb is
c           n  when  TRANSB = 'N' or 'n',  and is  k  otherwise.
c           Before entry with  TRANSB = 'N' or 'n',  the leading  k by n
c           part of the array  B  must contain the matrix  B,  otherwise
c           the leading  n by k  part of the array  B  must contain  the
c           matrix B.
c           Unchanged on exit.
c
c  LDB    - INTEGER.
c           On entry, LDB specifies the first dimension of B as declared
c           in the calling (sub) program. When  TRANSB = 'N' or 'n' then
c           LDB must be at least  max( 1, k ), otherwise  LDB must be at
c           least  max( 1, n ).
c           Unchanged on exit.
c
c  BETA   - REAL            .
c           On entry,  BETA  specifies the scalar  beta.  When  BETA  is
c           supplied as zero then C need not be set on input.
c           Unchanged on exit.
c
c  C      - REAL             array of DIMENSION ( LDC, n ).
c           Before entry, the leading  m by n  part of the array  C must
c           contain the matrix  C,  except when  beta  is zero, in which
c           case C need not be set on entry.
c           On exit, the array  C  is overwritten by the  m by n  matrix
c           ( alpha*op( A )*op( B ) + beta*C ).
c
c  LDC    - INTEGER.
c           On entry, LDC specifies the first dimension of C as declared
c           in  the  calling  (sub)  program.   LDC  must  be  at  least
c           max( 1, m ).
c           Unchanged on exit.
c
c
c     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
c     .. External Subroutines ..
      EXTERNAL           XERBLA
c     .. Intrinsic Functions ..
      INTRINSIC          MAX
c     .. Local Scalars ..
      LOGICAL            NOTA, NOTB
      INTEGER            I, INFO, J, L, NCOLA, NROWA, NROWB
      REAL               TEMP
c     .. Parameters ..
      REAL               ONE         , ZERO
      PARAMETER        ( ONE = 1.0E+0, ZERO = 0.0E+0 )
c     ..
c     .. Executable Statements ..
c
c     Set  NOTA  and  NOTB  as  true if  A  and  B  respectively are not
c     transposed and set  NROWA, NCOLA and  NROWB  as the number of rows
c     and  columns of  A  and the  number of  rows  of  B  respectively.
c
      NOTA  = LSAME( TRANSA, 'N' )
      NOTB  = LSAME( TRANSB, 'N' )
      IF( NOTA )THEN
         NROWA = M
         NCOLA = K
      ELSE
         NROWA = K
         NCOLA = M
      END IF
      IF( NOTB )THEN
         NROWB = K
      ELSE
         NROWB = N
      END IF
c
c     Test the input parameters.
c
      INFO = 0
      IF(      ( .NOT.NOTA                 ).AND.
     &         ( .NOT.LSAME( TRANSA, 'C' ) ).AND.
     &         ( .NOT.LSAME( TRANSA, 'T' ) )      )THEN
         INFO = 1
      ELSE IF( ( .NOT.NOTB                 ).AND.
     &         ( .NOT.LSAME( TRANSB, 'C' ) ).AND.
     &         ( .NOT.LSAME( TRANSB, 'T' ) )      )THEN
         INFO = 2
      ELSE IF( M  .LT.0               )THEN
         INFO = 3
      ELSE IF( N  .LT.0               )THEN
         INFO = 4
      ELSE IF( K  .LT.0               )THEN
         INFO = 5
      ELSE IF( LDA.LT.MAX( 1, NROWA ) )THEN
         INFO = 8
      ELSE IF( LDB.LT.MAX( 1, NROWB ) )THEN
         INFO = 10
      ELSE IF( LDC.LT.MAX( 1, M     ) )THEN
         INFO = 13
      END IF
      IF( INFO.NE.0 )THEN
         CALL XERBLA( 'SGEMM ', INFO )
         RETURN
      END IF
c
c     Quick return if possible.
c
      IF( ( M.EQ.0 ).OR.( N.EQ.0 ).OR.
     &    ( ( ( ALPHA.EQ.ZERO ).OR.( K.EQ.0 ) ).AND.( BETA.EQ.ONE ) ) )
     &   RETURN
c
c     And if  alpha.eq.zero.
c
      IF( ALPHA.EQ.ZERO )THEN
         IF( BETA.EQ.ZERO )THEN
            DO 20, J = 1, N
               DO 10, I = 1, M
                  C( I, J ) = ZERO
   10          CONTINUE
   20       CONTINUE
         ELSE
            DO 40, J = 1, N
               DO 30, I = 1, M
                  C( I, J ) = BETA*C( I, J )
   30          CONTINUE
   40       CONTINUE
         END IF
         RETURN
      END IF
c
c     Start the operations.
c
      IF( NOTB )THEN

         IF( NOTA )THEN
c
c           Form  C := alpha*A*B + beta*C.
c
            DO 90, J = 1, N
               IF( BETA.EQ.ZERO )THEN
                  DO 50, I = 1, M
                     C( I, J ) = ZERO
   50             CONTINUE
               ELSE IF( BETA.NE.ONE )THEN
                  DO 60, I = 1, M
                     C( I, J ) = BETA*C( I, J )
   60             CONTINUE
               END IF
               DO 80, L = 1, K
                  IF( B( L, J ).NE.ZERO )THEN
                     TEMP = ALPHA*B( L, J )
                     DO 70, I = 1, M
                        C( I, J ) = C( I, J ) + TEMP*A( I, L )
   70                CONTINUE
                  END IF
   80          CONTINUE
   90       CONTINUE

         ELSE
c
c           Form  C := alpha*A'*B + beta*C
c
            DO 120, J = 1, N
               DO 110, I = 1, M
                  TEMP = ZERO
                  DO 100, L = 1, K
                     TEMP = TEMP + A( L, I )*B( L, J )
  100             CONTINUE
                  IF( BETA.EQ.ZERO )THEN
                     C( I, J ) = ALPHA*TEMP
                  ELSE
                     C( I, J ) = ALPHA*TEMP + BETA*C( I, J )
                  END IF
  110          CONTINUE
  120       CONTINUE
         END IF

      ELSE

         IF( NOTA )THEN
c
c           Form  C := alpha*A*B' + beta*C
c
            DO 170, J = 1, N
               IF( BETA.EQ.ZERO )THEN
                  DO 130, I = 1, M
                     C( I, J ) = ZERO
  130             CONTINUE
               ELSE IF( BETA.NE.ONE )THEN
                  DO 140, I = 1, M
                     C( I, J ) = BETA*C( I, J )
  140             CONTINUE
               END IF
               DO 160, L = 1, K
                  IF( B( J, L ).NE.ZERO )THEN
                     TEMP = ALPHA*B( J, L )
                     DO 150, I = 1, M
                        C( I, J ) = C( I, J ) + TEMP*A( I, L )
  150                CONTINUE
                  END IF
  160          CONTINUE
  170       CONTINUE

         ELSE
c
c           Form  C := alpha*A'*B' + beta*C
c
            DO 200, J = 1, N
               DO 190, I = 1, M
                  TEMP = ZERO
                  DO 180, L = 1, K
                     TEMP = TEMP + A( L, I )*B( J, L )
  180             CONTINUE
                  IF( BETA.EQ.ZERO )THEN
                     C( I, J ) = ALPHA*TEMP
                  ELSE
                     C( I, J ) = ALPHA*TEMP + BETA*C( I, J )
                  END IF
  190          CONTINUE
  200       CONTINUE
         END IF
      END IF
c
      RETURN
c
c          End of SGEMM .
      END


      SUBROUTINE STRSM( SIDE, UPLO, TRANSA, DIAG, M, N, ALPHA, A, LDA,
     &                   B, LDB )

c  Level 3 Blas routine.
c
c  -- Written on 8-February-1989.

c     .. Scalar Arguments ..
      CHARACTER*1        SIDE, UPLO, TRANSA, DIAG
      INTEGER            M, N, LDA, LDB
      REAL               ALPHA
c     .. Array Arguments ..
      REAL               A( LDA, * ), B( LDB, * )
c     ..
c
c  Purpose
c  =======
c
c  STRSM  solves one of the matrix equations
c
c     op( A )*X = alpha*B,   or   X*op( A ) = alpha*B,
c
c  where alpha is a scalar, X and B are m by n matrices, A is a unit, or
c  non-unit,  upper or lower triangular matrix  and  op( A )  is one  of
c
c     op( A ) = A   or   op( A ) = A'.
c
c  The matrix X is overwritten on B.
c
c  Parameters
c  ==========
c
c  SIDE   - CHARACTER*1.
c           On entry, SIDE specifies whether op( A ) appears on the left
c           or right of X as follows:
c
c              SIDE = 'L' or 'l'   op( A )*X = alpha*B.
c
c              SIDE = 'R' or 'r'   X*op( A ) = alpha*B.
c
c           Unchanged on exit.
c
c  UPLO   - CHARACTER*1.
c           On entry, UPLO specifies whether the matrix A is an upper or
c           lower triangular matrix as follows:
c
c              UPLO = 'U' or 'u'   A is an upper triangular matrix.
c
c              UPLO = 'L' or 'l'   A is a lower triangular matrix.
c
c           Unchanged on exit.
c
c  TRANSA - CHARACTER*1.
c           On entry, TRANSA specifies the form of op( A ) to be used in
c           the matrix multiplication as follows:
c
c              TRANSA = 'N' or 'n'   op( A ) = A.
c
c              TRANSA = 'T' or 't'   op( A ) = A'.
c
c              TRANSA = 'C' or 'c'   op( A ) = A'.
c
c           Unchanged on exit.
c
c  DIAG   - CHARACTER*1.
c           On entry, DIAG specifies whether or not A is unit triangular
c           as follows:
c
c              DIAG = 'U' or 'u'   A is assumed to be unit triangular.
c
c              DIAG = 'N' or 'n'   A is not assumed to be unit
c                                  triangular.
c
c           Unchanged on exit.
c
c  M      - INTEGER.
c           On entry, M specifies the number of rows of B. M must be at
c           least zero.
c           Unchanged on exit.
c
c  N      - INTEGER.
c           On entry, N specifies the number of columns of B.  N must be
c           at least zero.
c           Unchanged on exit.
c
c  ALPHA  - REAL            .
c           On entry,  ALPHA specifies the scalar  alpha. When  alpha is
c           zero then  A is not referenced and  B need not be set before
c           entry.
c           Unchanged on exit.
c
c  A      - REAL             array of DIMENSION ( LDA, k ), where k is m
c           when  SIDE = 'L' or 'l'  and is  n  when  SIDE = 'R' or 'r'.
c           Before entry  with  UPLO = 'U' or 'u',  the  leading  k by k
c           upper triangular part of the array  A must contain the upper
c           triangular matrix  and the strictly lower triangular part of
c           A is not referenced.
c           Before entry  with  UPLO = 'L' or 'l',  the  leading  k by k
c           lower triangular part of the array  A must contain the lower
c           triangular matrix  and the strictly upper triangular part of
c           A is not referenced.
c           Note that when  DIAG = 'U' or 'u',  the diagonal elements of
c           A  are not referenced either,  but are assumed to be  unity.
c           Unchanged on exit.
c
c  LDA    - INTEGER.
c           On entry, LDA specifies the first dimension of A as declared
c           in the calling (sub) program.  When  SIDE = 'L' or 'l'  then
c           LDA  must be at least  max( 1, m ),  when  SIDE = 'R' or 'r'
c           then LDA must be at least max( 1, n ).
c           Unchanged on exit.
c
c  B      - REAL             array of DIMENSION ( LDB, n ).
c           Before entry,  the leading  m by n part of the array  B must
c           contain  the  right-hand  side  matrix  B,  and  on exit  is
c           overwritten by the solution matrix  X.
c
c  LDB    - INTEGER.
c           On entry, LDB specifies the first dimension of B as declared
c           in  the  calling  (sub)  program.   LDB  must  be  at  least
c           max( 1, m ).
c           Unchanged on exit.
c
c
c     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
c     .. External Subroutines ..
      EXTERNAL           XERBLA
c     .. Intrinsic Functions ..
      INTRINSIC          MAX
c     .. Local Scalars ..
      LOGICAL            LSIDE, NOUNIT, UPPER
      INTEGER            I, INFO, J, K, NROWA
      REAL               TEMP
c     .. Parameters ..
      REAL               ONE         , ZERO
      PARAMETER        ( ONE = 1.0E+0, ZERO = 0.0E+0 )
c     ..

c
c     Test the input parameters.
c
      LSIDE  = LSAME( SIDE  , 'L' )
      IF( LSIDE )THEN
         NROWA = M
      ELSE
         NROWA = N
      END IF
      NOUNIT = LSAME( DIAG  , 'N' )
      UPPER  = LSAME( UPLO  , 'U' )
c
      INFO   = 0
      IF(      ( .NOT.LSIDE                ).AND.
     &         ( .NOT.LSAME( SIDE  , 'R' ) )      )THEN
         INFO = 1
      ELSE IF( ( .NOT.UPPER                ).AND.
     &         ( .NOT.LSAME( UPLO  , 'L' ) )      )THEN
         INFO = 2
      ELSE IF( ( .NOT.LSAME( TRANSA, 'N' ) ).AND.
     &         ( .NOT.LSAME( TRANSA, 'T' ) ).AND.
     &         ( .NOT.LSAME( TRANSA, 'C' ) )      )THEN
         INFO = 3
      ELSE IF( ( .NOT.LSAME( DIAG  , 'U' ) ).AND.
     &         ( .NOT.LSAME( DIAG  , 'N' ) )      )THEN
         INFO = 4
      ELSE IF( M  .LT.0               )THEN
         INFO = 5
      ELSE IF( N  .LT.0               )THEN
         INFO = 6
      ELSE IF( LDA.LT.MAX( 1, NROWA ) )THEN
         INFO = 9
      ELSE IF( LDB.LT.MAX( 1, M     ) )THEN
         INFO = 11
      END IF
      IF( INFO.NE.0 )THEN
         CALL XERBLA( 'STRSM ', INFO )
         RETURN
      END IF
c
c     Quick return if possible.
c
      IF( N.EQ.0 )  RETURN
c
c     And when  alpha.eq.zero.
c
      IF( ALPHA.EQ.ZERO )THEN
         DO 20, J = 1, N
            DO 10, I = 1, M
               B( I, J ) = ZERO
   10       CONTINUE
   20    CONTINUE
         RETURN
      END IF
c
c     Start the operations.
c
      IF( LSIDE )THEN

         IF( LSAME( TRANSA, 'N' ) )THEN
c
c           Form  B := alpha*inv( A )*B.
c
            IF( UPPER )THEN
               DO 60, J = 1, N
                  IF( ALPHA.NE.ONE )THEN
                     DO 30, I = 1, M
                        B( I, J ) = ALPHA*B( I, J )
   30                CONTINUE
                  END IF
                  DO 50, K = M, 1, -1
                     IF( B( K, J ).NE.ZERO )THEN
                        IF( NOUNIT )
     &                     B( K, J ) = B( K, J )/A( K, K )
                        DO 40, I = 1, K - 1
                           B( I, J ) = B( I, J ) - B( K, J )*A( I, K )
   40                   CONTINUE
                     END IF
   50             CONTINUE
   60          CONTINUE

            ELSE

               DO 100, J = 1, N
                  IF( ALPHA.NE.ONE )THEN
                     DO 70, I = 1, M
                        B( I, J ) = ALPHA*B( I, J )
   70                CONTINUE
                  END IF
                  DO 90 K = 1, M
                     IF( B( K, J ).NE.ZERO )THEN
                        IF( NOUNIT )
     &                     B( K, J ) = B( K, J )/A( K, K )
                        DO 80, I = K + 1, M
                           B( I, J ) = B( I, J ) - B( K, J )*A( I, K )
   80                   CONTINUE
                     END IF
   90             CONTINUE
  100          CONTINUE
            END IF

         ELSE
c
c           Form  B := alpha*inv( A' )*B.
c
            IF( UPPER )THEN
               DO 130, J = 1, N
                  DO 120, I = 1, M
                     TEMP = ALPHA*B( I, J )
                     DO 110, K = 1, I - 1
                        TEMP = TEMP - A( K, I )*B( K, J )
  110                CONTINUE
                     IF( NOUNIT )
     &                  TEMP = TEMP/A( I, I )
                     B( I, J ) = TEMP
  120             CONTINUE
  130          CONTINUE
            ELSE
               DO 160, J = 1, N
                  DO 150, I = M, 1, -1
                     TEMP = ALPHA*B( I, J )
                     DO 140, K = I + 1, M
                        TEMP = TEMP - A( K, I )*B( K, J )
  140                CONTINUE
                     IF( NOUNIT )
     &                  TEMP = TEMP/A( I, I )
                     B( I, J ) = TEMP
  150             CONTINUE
  160          CONTINUE
            END IF
         END IF

      ELSE

         IF( LSAME( TRANSA, 'N' ) )THEN
c
c           Form  B := alpha*B*inv( A ).
c
            IF( UPPER )THEN
               DO 210, J = 1, N
                  IF( ALPHA.NE.ONE )THEN
                     DO 170, I = 1, M
                        B( I, J ) = ALPHA*B( I, J )
  170                CONTINUE
                  END IF
                  DO 190, K = 1, J - 1
                     IF( A( K, J ).NE.ZERO )THEN
                        DO 180, I = 1, M
                           B( I, J ) = B( I, J ) - A( K, J )*B( I, K )
  180                   CONTINUE
                     END IF
  190             CONTINUE
                  IF( NOUNIT )THEN
                     TEMP = ONE/A( J, J )
                     DO 200, I = 1, M
                        B( I, J ) = TEMP*B( I, J )
  200                CONTINUE
                  END IF
  210          CONTINUE

            ELSE

               DO 260, J = N, 1, -1
                  IF( ALPHA.NE.ONE )THEN
                     DO 220, I = 1, M
                        B( I, J ) = ALPHA*B( I, J )
  220                CONTINUE
                  END IF
                  DO 240, K = J + 1, N
                     IF( A( K, J ).NE.ZERO )THEN
                        DO 230, I = 1, M
                           B( I, J ) = B( I, J ) - A( K, J )*B( I, K )
  230                   CONTINUE
                     END IF
  240             CONTINUE
                  IF( NOUNIT )THEN
                     TEMP = ONE/A( J, J )
                     DO 250, I = 1, M
                       B( I, J ) = TEMP*B( I, J )
  250                CONTINUE
                  END IF
  260          CONTINUE
            END IF

         ELSE
c
c           Form  B := alpha*B*inv( A' ).
c
            IF( UPPER )THEN
               DO 310, K = N, 1, -1
                  IF( NOUNIT )THEN
                     TEMP = ONE/A( K, K )
                     DO 270, I = 1, M
                        B( I, K ) = TEMP*B( I, K )
  270                CONTINUE
                  END IF
                  DO 290, J = 1, K - 1
                     IF( A( J, K ).NE.ZERO )THEN
                        TEMP = A( J, K )
                        DO 280, I = 1, M
                           B( I, J ) = B( I, J ) - TEMP*B( I, K )
  280                   CONTINUE
                     END IF
  290             CONTINUE
                  IF( ALPHA.NE.ONE )THEN
                     DO 300, I = 1, M
                        B( I, K ) = ALPHA*B( I, K )
  300                CONTINUE
                  END IF
  310          CONTINUE

            ELSE

               DO 360, K = 1, N
                  IF( NOUNIT )THEN
                     TEMP = ONE/A( K, K )
                     DO 320, I = 1, M
                        B( I, K ) = TEMP*B( I, K )
  320                CONTINUE
                  END IF
                  DO 340, J = K + 1, N
                     IF( A( J, K ).NE.ZERO )THEN
                        TEMP = A( J, K )
                        DO 330, I = 1, M
                           B( I, J ) = B( I, J ) - TEMP*B( I, K )
  330                   CONTINUE
                     END IF
  340             CONTINUE
                  IF( ALPHA.NE.ONE )THEN
                     DO 350, I = 1, M
                        B( I, K ) = ALPHA*B( I, K )
  350                CONTINUE
                  END IF
  360          CONTINUE
            END IF
         END IF
      END IF
c
      RETURN
c
c          End of STRSM .
      END


      SUBROUTINE STBSV( UPLO, TRANS, DIAG, N, K, A, LDA, X, INCX )

c  Level 2 Blas routine.
c
c  -- Written on 22-October-1986.
c
c     .. Scalar Arguments ..
      INTEGER            INCX, K, LDA, N
      CHARACTER*1        DIAG, TRANS, UPLO
c     .. Array Arguments ..
      REAL               A( LDA, * ), X( * )
c     ..
c
c  Purpose
c  =======
c
c  STBSV  solves one of the systems of equations
c
c     A*x = b,   or   A'*x = b,
c
c  where b and x are n element vectors and A is an n by n unit, or
c  non-unit, upper or lower triangular band matrix, with ( k + 1 )
c  diagonals.
c
c  No test for singularity or near-singularity is included in this
c  routine. Such tests must be performed before calling this routine.
c
c  Parameters
c  ==========
c
c  UPLO   - CHARACTER*1.
c           On entry, UPLO specifies whether the matrix is an upper or
c           lower triangular matrix as follows:
c
c              UPLO = 'U' or 'u'   A is an upper triangular matrix.
c
c              UPLO = 'L' or 'l'   A is a lower triangular matrix.
c
c           Unchanged on exit.
c
c  TRANS  - CHARACTER*1.
c           On entry, TRANS specifies the equations to be solved as
c           follows:
c
c              TRANS = 'N' or 'n'   A*x = b.
c
c              TRANS = 'T' or 't'   A'*x = b.
c
c              TRANS = 'C' or 'c'   A'*x = b.
c
c           Unchanged on exit.
c
c  DIAG   - CHARACTER*1.
c           On entry, DIAG specifies whether or not A is unit
c           triangular as follows:
c
c              DIAG = 'U' or 'u'   A is assumed to be unit triangular.
c
c              DIAG = 'N' or 'n'   A is not assumed to be unit
c                                  triangular.
c
c           Unchanged on exit.
c
c  N      - INTEGER.
c           On entry, N specifies the order of the matrix A.
c           N must be at least zero.
c           Unchanged on exit.
c
c  K      - INTEGER.
c           On entry with UPLO = 'U' or 'u', K specifies the number of
c           super-diagonals of the matrix A.
c           On entry with UPLO = 'L' or 'l', K specifies the number of
c           sub-diagonals of the matrix A.
c           K must satisfy  0 .le. K.
c           Unchanged on exit.
c
c  A      - REAL             array of DIMENSION ( LDA, n ).
c           Before entry with UPLO = 'U' or 'u', the leading ( k + 1 )
c           by n part of the array A must contain the upper triangular
c           band part of the matrix of coefficients, supplied column by
c           column, with the leading diagonal of the matrix in row
c           ( k + 1 ) of the array, the first super-diagonal starting at
c           position 2 in row k, and so on. The top left k by k triangle
c           of the array A is not referenced.
c           The following program segment will transfer an upper
c           triangular band matrix from conventional full matrix storage
c           to band storage:
c
c                 DO 20, J = 1, N
c                    M = K + 1 - J
c                    DO 10, I = MAX( 1, J - K ), J
c                       A( M + I, J ) = matrix( I, J )
c              10    CONTINUE
c              20 CONTINUE
c
c           Before entry with UPLO = 'L' or 'l', the leading ( k + 1 )
c           by n part of the array A must contain the lower triangular
c           band part of the matrix of coefficients, supplied column by
c           column, with the leading diagonal of the matrix in row 1 of
c           the array, the first sub-diagonal starting at position 1 in
c           row 2, and so on. The bottom right k by k triangle of the
c           array A is not referenced.
c           The following program segment will transfer a lower
c           triangular band matrix from conventional full matrix storage
c           to band storage:
c
c                 DO 20, J = 1, N
c                    M = 1 - J
c                    DO 10, I = J, MIN( N, J + K )
c                       A( M + I, J ) = matrix( I, J )
c              10    CONTINUE
c              20 CONTINUE
c
c           Note that when DIAG = 'U' or 'u' the elements of the array A
c           corresponding to the diagonal elements of the matrix are not
c           referenced, but are assumed to be unity.
c           Unchanged on exit.
c
c  LDA    - INTEGER.
c           On entry, LDA specifies the first dimension of A as declared
c           in the calling (sub) program. LDA must be at least
c           ( k + 1 ).
c           Unchanged on exit.
c
c  X      - REAL             array of dimension at least
c           ( 1 + ( n - 1 )*abs( INCX ) ).
c           Before entry, the incremented array X must contain the n
c           element right-hand side vector b. On exit, X is overwritten
c           with the solution vector x.
c
c  INCX   - INTEGER.
c           On entry, INCX specifies the increment for the elements of
c           X. INCX must not be zero.
c           Unchanged on exit.
c
c
c     .. Parameters ..
      REAL               ZERO
      PARAMETER        ( ZERO = 0.0E+0 )
c     .. Local Scalars ..
      REAL               TEMP
      INTEGER            I, INFO, IX, J, JX, KPLUS1, KX, L
      LOGICAL            NOUNIT
c     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
c     .. External Subroutines ..
      EXTERNAL           XERBLA
c     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
c     ..
c     .. Executable Statements ..
c
c     Test the input parameters.
c
      INFO = 0
      IF     ( .NOT.LSAME( UPLO , 'U' ).AND.
     &         .NOT.LSAME( UPLO , 'L' )      )THEN
         INFO = 1
      ELSE IF( .NOT.LSAME( TRANS, 'N' ).AND.
     &         .NOT.LSAME( TRANS, 'T' ).AND.
     &         .NOT.LSAME( TRANS, 'C' )      )THEN
         INFO = 2
      ELSE IF( .NOT.LSAME( DIAG , 'U' ).AND.
     &         .NOT.LSAME( DIAG , 'N' )      )THEN
         INFO = 3
      ELSE IF( N.LT.0 )THEN
         INFO = 4
      ELSE IF( K.LT.0 )THEN
         INFO = 5
      ELSE IF( LDA.LT.( K + 1 ) )THEN
         INFO = 7
      ELSE IF( INCX.EQ.0 )THEN
         INFO = 9
      END IF
      IF( INFO.NE.0 )THEN
         CALL XERBLA( 'STBSV ', INFO )
         RETURN
      END IF
c
c     Quick return if possible.
c
      IF( N.EQ.0 )  RETURN
c
      NOUNIT = LSAME( DIAG, 'N' )
c
c     Set up the start point in X if the increment is not unity. This
c     will be  ( N - 1 )*INCX  too small for descending loops.
c
      IF( INCX.LE.0 )THEN
         KX = 1 - ( N - 1 )*INCX
      ELSE IF( INCX.NE.1 )THEN
         KX = 1
      END IF
c
c     Start the operations. In this version the elements of A are
c     accessed by sequentially with one pass through A.
c
      IF( LSAME( TRANS, 'N' ) )THEN
c
c        Form  x := inv( A )*x.
c
         IF( LSAME( UPLO, 'U' ) )THEN
            KPLUS1 = K + 1
            IF( INCX.EQ.1 )THEN
               DO 20, J = N, 1, -1
                  IF( X( J ).NE.ZERO )THEN
                     L = KPLUS1 - J
                     IF( NOUNIT )
     &                  X( J ) = X( J )/A( KPLUS1, J )
                     TEMP = X( J )
                     DO 10, I = J - 1, MAX( 1, J - K ), -1
                        X( I ) = X( I ) - TEMP*A( L + I, J )
   10                CONTINUE
                  END IF
   20          CONTINUE
            ELSE
               KX = KX + ( N - 1 )*INCX
               JX = KX
               DO 40, J = N, 1, -1
                  KX = KX - INCX
                  IF( X( JX ).NE.ZERO )THEN
                     IX = KX
                     L  = KPLUS1 - J
                     IF( NOUNIT )
     &                  X( JX ) = X( JX )/A( KPLUS1, J )
                     TEMP = X( JX )
                     DO 30, I = J - 1, MAX( 1, J - K ), -1
                        X( IX ) = X( IX ) - TEMP*A( L + I, J )
                        IX      = IX      - INCX
   30                CONTINUE
                  END IF
                  JX = JX - INCX
   40          CONTINUE
            END IF
         ELSE
            IF( INCX.EQ.1 )THEN
               DO 60, J = 1, N
                  IF( X( J ).NE.ZERO )THEN
                     L = 1 - J
                     IF( NOUNIT )
     &                  X( J ) = X( J )/A( 1, J )
                     TEMP = X( J )
                     DO 50, I = J + 1, MIN( N, J + K )
                        X( I ) = X( I ) - TEMP*A( L + I, J )
   50                CONTINUE
                  END IF
   60          CONTINUE
            ELSE
               JX = KX
               DO 80, J = 1, N
                  KX = KX + INCX
                  IF( X( JX ).NE.ZERO )THEN
                     IX = KX
                     L  = 1  - J
                     IF( NOUNIT )
     &                  X( JX ) = X( JX )/A( 1, J )
                     TEMP = X( JX )
                     DO 70, I = J + 1, MIN( N, J + K )
                        X( IX ) = X( IX ) - TEMP*A( L + I, J )
                        IX      = IX      + INCX
   70                CONTINUE
                  END IF
                  JX = JX + INCX
   80          CONTINUE
            END IF
         END IF
      ELSE
c
c        Form  x := inv( A')*x.
c
         IF( LSAME( UPLO, 'U' ) )THEN
            KPLUS1 = K + 1
            IF( INCX.EQ.1 )THEN
               DO 100, J = 1, N
                  TEMP = X( J )
                  L    = KPLUS1 - J
                  DO 90, I = MAX( 1, J - K ), J - 1
                     TEMP = TEMP - A( L + I, J )*X( I )
   90             CONTINUE
                  IF( NOUNIT )
     &               TEMP = TEMP/A( KPLUS1, J )
                  X( J ) = TEMP
  100          CONTINUE
            ELSE
               JX = KX
               DO 120, J = 1, N
                  TEMP = X( JX )
                  IX   = KX
                  L    = KPLUS1  - J
                  DO 110, I = MAX( 1, J - K ), J - 1
                     TEMP = TEMP - A( L + I, J )*X( IX )
                     IX   = IX   + INCX
  110             CONTINUE
                  IF( NOUNIT )
     &               TEMP = TEMP/A( KPLUS1, J )
                  X( JX ) = TEMP
                  JX      = JX   + INCX
                  IF( J.GT.K )
     &               KX = KX + INCX
  120          CONTINUE
            END IF
         ELSE
            IF( INCX.EQ.1 )THEN
               DO 140, J = N, 1, -1
                  TEMP = X( J )
                  L    = 1      - J
                  DO 130, I = MIN( N, J + K ), J + 1, -1
                     TEMP = TEMP - A( L + I, J )*X( I )
  130             CONTINUE
                  IF( NOUNIT )
     &               TEMP = TEMP/A( 1, J )
                  X( J ) = TEMP
  140          CONTINUE
            ELSE
               KX = KX + ( N - 1 )*INCX
               JX = KX
               DO 160, J = N, 1, -1
                  TEMP = X( JX )
                  IX   = KX
                  L    = 1       - J
                  DO 150, I = MIN( N, J + K ), J + 1, -1
                     TEMP = TEMP - A( L + I, J )*X( IX )
                     IX   = IX   - INCX
  150             CONTINUE
                  IF( NOUNIT )
     &               TEMP = TEMP/A( 1, J )
                  X( JX ) = TEMP
                  JX      = JX   - INCX
                  IF( ( N - J ).GE.K )
     &               KX = KX - INCX
  160          CONTINUE
            END IF
         END IF
      END IF
c
      RETURN
c
c          End of STBSV .
      END


      SUBROUTINE SGEMV( TRANS, M, N, ALPHA, A, LDA, X, INCX,
     &                   BETA, Y, INCY )

c  Level 2 Blas routine.
c
c  -- Written on 22-October-1986.
c
c     .. Scalar Arguments ..
      REAL               ALPHA, BETA
      INTEGER            INCX, INCY, LDA, M, N
      CHARACTER*1        TRANS
c     .. Array Arguments ..
      REAL               A( LDA, * ), X( * ), Y( * )
c     ..
c
c  Purpose
c  =======
c
c  SGEMV  performs one of the matrix-vector operations
c
c     y := alpha*A*x + beta*y,   or   y := alpha*A'*x + beta*y,
c
c  where alpha and beta are scalars, x and y are vectors and A is an
c  m by n matrix.
c
c  Parameters
c  ==========
c
c  TRANS  - CHARACTER*1.
c           On entry, TRANS specifies the operation to be performed as
c           follows:
c
c              TRANS = 'N' or 'n'   y := alpha*A*x + beta*y.
c
c              TRANS = 'T' or 't'   y := alpha*A'*x + beta*y.
c
c              TRANS = 'C' or 'c'   y := alpha*A'*x + beta*y.
c
c           Unchanged on exit.
c
c  M      - INTEGER.
c           On entry, M specifies the number of rows of the matrix A.
c           M must be at least zero.
c           Unchanged on exit.
c
c  N      - INTEGER.
c           On entry, N specifies the number of columns of the matrix A.
c           N must be at least zero.
c           Unchanged on exit.
c
c  ALPHA  - REAL            .
c           On entry, ALPHA specifies the scalar alpha.
c           Unchanged on exit.
c
c  A      - REAL             array of DIMENSION ( LDA, n ).
c           Before entry, the leading m by n part of the array A must
c           contain the matrix of coefficients.
c           Unchanged on exit.
c
c  LDA    - INTEGER.
c           On entry, LDA specifies the first dimension of A as declared
c           in the calling (sub) program. LDA must be at least
c           max( 1, m ).
c           Unchanged on exit.
c
c  X      - REAL             array of DIMENSION at least
c           ( 1 + ( n - 1 )*abs( INCX ) ) when TRANS = 'N' or 'n'
c           and at least
c           ( 1 + ( m - 1 )*abs( INCX ) ) otherwise.
c           Before entry, the incremented array X must contain the
c           vector x.
c           Unchanged on exit.
c
c  INCX   - INTEGER.
c           On entry, INCX specifies the increment for the elements of
c           X. INCX must not be zero.
c           Unchanged on exit.
c
c  BETA   - REAL            .
c           On entry, BETA specifies the scalar beta. When BETA is
c           supplied as zero then Y need not be set on input.
c           Unchanged on exit.
c
c  Y      - REAL             array of DIMENSION at least
c           ( 1 + ( m - 1 )*abs( INCY ) ) when TRANS = 'N' or 'n'
c           and at least
c           ( 1 + ( n - 1 )*abs( INCY ) ) otherwise.
c           Before entry with BETA non-zero, the incremented array Y
c           must contain the vector y. On exit, Y is overwritten by the
c           updated vector y.
c
c  INCY   - INTEGER.
c           On entry, INCY specifies the increment for the elements of
c           Y. INCY must not be zero.
c           Unchanged on exit.
c
c
c     .. Parameters ..
      REAL               ONE         , ZERO
      PARAMETER        ( ONE = 1.0E+0, ZERO = 0.0E+0 )
c     .. Local Scalars ..
      REAL               TEMP
      INTEGER            I, INFO, IX, IY, J, JX, JY, KX, KY, LENX, LENY
c     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
c     .. External Subroutines ..
      EXTERNAL           XERBLA
c     .. Intrinsic Functions ..
      INTRINSIC          MAX
c     ..
c     .. Executable Statements ..
c
c     Test the input parameters.
c
      INFO = 0
      IF     ( .NOT.LSAME( TRANS, 'N' ).AND.
     &         .NOT.LSAME( TRANS, 'T' ).AND.
     &         .NOT.LSAME( TRANS, 'C' )      )THEN
         INFO = 1
      ELSE IF( M.LT.0 )THEN
         INFO = 2
      ELSE IF( N.LT.0 )THEN
         INFO = 3
      ELSE IF( LDA.LT.MAX( 1, M ) )THEN
         INFO = 6
      ELSE IF( INCX.EQ.0 )THEN
         INFO = 8
      ELSE IF( INCY.EQ.0 )THEN
         INFO = 11
      END IF
      IF( INFO.NE.0 )THEN
         CALL XERBLA( 'SGEMV ', INFO )
         RETURN
      END IF
c
c     Quick return if possible.
c
      IF( ( M.EQ.0 ).OR.( N.EQ.0 ).OR.
     &    ( ( ALPHA.EQ.ZERO ).AND.( BETA.EQ.ONE ) ) )
     &   RETURN
c
c     Set  LENX  and  LENY, the lengths of the vectors x and y, and set
c     up the start points in  X  and  Y.
c
      IF( LSAME( TRANS, 'N' ) )THEN
         LENX = N
         LENY = M
      ELSE
         LENX = M
         LENY = N
      END IF
      IF( INCX.GT.0 )THEN
         KX = 1
      ELSE
         KX = 1 - ( LENX - 1 )*INCX
      END IF
      IF( INCY.GT.0 )THEN
         KY = 1
      ELSE
         KY = 1 - ( LENY - 1 )*INCY
      END IF
c
c     Start the operations. In this version the elements of A are
c     accessed sequentially with one pass through A.
c
c     First form  y := beta*y.
c
      IF( BETA.NE.ONE )THEN
         IF( INCY.EQ.1 )THEN
            IF( BETA.EQ.ZERO )THEN
               DO 10, I = 1, LENY
                  Y( I ) = ZERO
   10          CONTINUE
            ELSE
               DO 20, I = 1, LENY
                  Y( I ) = BETA*Y( I )
   20          CONTINUE
            END IF
         ELSE
            IY = KY
            IF( BETA.EQ.ZERO )THEN
               DO 30, I = 1, LENY
                  Y( IY ) = ZERO
                  IY      = IY   + INCY
   30          CONTINUE
            ELSE
               DO 40, I = 1, LENY
                  Y( IY ) = BETA*Y( IY )
                  IY      = IY           + INCY
   40          CONTINUE
            END IF
         END IF
      END IF

      IF( ALPHA.EQ.ZERO )  RETURN

      IF( LSAME( TRANS, 'N' ) )THEN
c
c        Form  y := alpha*A*x + y.
c
         JX = KX
         IF( INCY.EQ.1 )THEN
            DO 60, J = 1, N
               IF( X( JX ).NE.ZERO )THEN
                  TEMP = ALPHA*X( JX )
                  DO 50, I = 1, M
                     Y( I ) = Y( I ) + TEMP*A( I, J )
   50             CONTINUE
               END IF
               JX = JX + INCX
   60       CONTINUE
         ELSE
            DO 80, J = 1, N
               IF( X( JX ).NE.ZERO )THEN
                  TEMP = ALPHA*X( JX )
                  IY   = KY
                  DO 70, I = 1, M
                     Y( IY ) = Y( IY ) + TEMP*A( I, J )
                     IY      = IY      + INCY
   70             CONTINUE
               END IF
               JX = JX + INCX
   80       CONTINUE
         END IF
      ELSE
c
c        Form  y := alpha*A'*x + y.
c
         JY = KY
         IF( INCX.EQ.1 )THEN
            DO 100, J = 1, N
               TEMP = ZERO
               DO 90, I = 1, M
                  TEMP = TEMP + A( I, J )*X( I )
   90          CONTINUE
               Y( JY ) = Y( JY ) + ALPHA*TEMP
               JY      = JY      + INCY
  100       CONTINUE
         ELSE
            DO 120, J = 1, N
               TEMP = ZERO
               IX   = KX
               DO 110, I = 1, M
                  TEMP = TEMP + A( I, J )*X( IX )
                  IX   = IX   + INCX
  110          CONTINUE
               Y( JY ) = Y( JY ) + ALPHA*TEMP
               JY      = JY      + INCY
  120       CONTINUE
         END IF
      END IF
c
      RETURN
c
c          End of SGEMV .
      END


      SUBROUTINE SGER( M, N, ALPHA, X, INCX, Y, INCY, A, LDA )

c  Level 2 Blas routine.
c
c  -- Written on 22-October-1986.
c
c     .. Scalar Arguments ..
      REAL               ALPHA
      INTEGER            INCX, INCY, LDA, M, N
c     .. Array Arguments ..
      REAL               A( LDA, * ), X( * ), Y( * )
c     ..
c
c  Purpose
c  =======
c
c  SGER   performs the rank 1 operation
c
c     A := alpha*x*y' + A,
c
c  where alpha is a scalar, x is an m element vector, y is an n element
c  vector and A is an m by n matrix.
c
c  Parameters
c  ==========
c
c  M      - INTEGER.
c           On entry, M specifies the number of rows of the matrix A.
c           M must be at least zero.
c           Unchanged on exit.
c
c  N      - INTEGER.
c           On entry, N specifies the number of columns of the matrix A.
c           N must be at least zero.
c           Unchanged on exit.
c
c  ALPHA  - REAL            .
c           On entry, ALPHA specifies the scalar alpha.
c           Unchanged on exit.
c
c  X      - REAL             array of dimension at least
c           ( 1 + ( m - 1 )*abs( INCX ) ).
c           Before entry, the incremented array X must contain the m
c           element vector x.
c           Unchanged on exit.
c
c  INCX   - INTEGER.
c           On entry, INCX specifies the increment for the elements of
c           X. INCX must not be zero.
c           Unchanged on exit.
c
c  Y      - REAL             array of dimension at least
c           ( 1 + ( n - 1 )*abs( INCY ) ).
c           Before entry, the incremented array Y must contain the n
c           element vector y.
c           Unchanged on exit.
c
c  INCY   - INTEGER.
c           On entry, INCY specifies the increment for the elements of
c           Y. INCY must not be zero.
c           Unchanged on exit.
c
c  A      - REAL             array of DIMENSION ( LDA, n ).
c           Before entry, the leading m by n part of the array A must
c           contain the matrix of coefficients. On exit, A is
c           overwritten by the updated matrix.
c
c  LDA    - INTEGER.
c           On entry, LDA specifies the first dimension of A as declared
c           in the calling (sub) program. LDA must be at least
c           max( 1, m ).
c           Unchanged on exit.
c
c
c     .. Parameters ..
      REAL               ZERO
      PARAMETER        ( ZERO = 0.0E+0 )
c     .. Local Scalars ..
      REAL               TEMP
      INTEGER            I, INFO, IX, J, JY, KX
c     .. External Subroutines ..
      EXTERNAL           XERBLA
c     .. Intrinsic Functions ..
      INTRINSIC          MAX
c     ..
c     .. Executable Statements ..
c
c     Test the input parameters.
c
      INFO = 0
      IF     ( M.LT.0 )THEN
         INFO = 1
      ELSE IF( N.LT.0 )THEN
         INFO = 2
      ELSE IF( INCX.EQ.0 )THEN
         INFO = 5
      ELSE IF( INCY.EQ.0 )THEN
         INFO = 7
      ELSE IF( LDA.LT.MAX( 1, M ) )THEN
         INFO = 9
      END IF
      IF( INFO.NE.0 )THEN
         CALL XERBLA( 'SGER  ', INFO )
         RETURN
      END IF
c
c     Quick return if possible.
c
      IF( ( M.EQ.0 ).OR.( N.EQ.0 ).OR.( ALPHA.EQ.ZERO ) )  RETURN
c
c     Start the operations. In this version the elements of A are
c     accessed sequentially with one pass through A.
c
      IF( INCY.GT.0 )THEN
         JY = 1
      ELSE
         JY = 1 - ( N - 1 )*INCY
      END IF
      IF( INCX.EQ.1 )THEN
         DO 20, J = 1, N
            IF( Y( JY ).NE.ZERO )THEN
               TEMP = ALPHA*Y( JY )
               DO 10, I = 1, M
                  A( I, J ) = A( I, J ) + X( I )*TEMP
   10          CONTINUE
            END IF
            JY = JY + INCY
   20    CONTINUE
      ELSE
         IF( INCX.GT.0 )THEN
            KX = 1
         ELSE
            KX = 1 - ( M - 1 )*INCX
         END IF
         DO 40, J = 1, N
            IF( Y( JY ).NE.ZERO )THEN
               TEMP = ALPHA*Y( JY )
               IX   = KX
               DO 30, I = 1, M
                  A( I, J ) = A( I, J ) + X( IX )*TEMP
                  IX        = IX        + INCX
   30          CONTINUE
            END IF
            JY = JY + INCY
   40    CONTINUE
      END IF
c
      RETURN
c
c           End of SGER  .
      END


      INTEGER FUNCTION ISAMAX( N, SX, INCX )

c  Level 1 Blas routine.

c     Finds the index of element having max. absolute value.
c     Jack Dongarra, LINPACK, 3/11/78.

c     .. Scalar Arguments ..

      INTEGER   INCX, N
c     ..
c     .. Array Arguments ..

      REAL      SX( 1 )
c     ..
c     .. Local Scalars ..

      INTEGER   I, IX
      REAL      SMAX
c     ..
c     .. Intrinsic Functions ..

      INTRINSIC ABS
c     ..

      ISAMAX = 0
      IF( N.LT.1 ) RETURN

      ISAMAX = 1
      IF( N.EQ.1 ) RETURN

      IF( INCX.EQ.1 ) GO TO 30

c        code for increment not equal to 1

      IX    = 1
      SMAX  = ABS( SX( 1 ) )
      IX    = IX + INCX

      DO 20 I = 2, N

         IF( ABS( SX( IX ) ).LE.SMAX ) GO TO 10

         ISAMAX = I
         SMAX   = ABS( SX( IX ) )

   10    CONTINUE
         IX  = IX + INCX

   20 CONTINUE

      RETURN

c        code for increment equal to 1

   30 CONTINUE
      SMAX  = ABS( SX( 1 ) )

      DO 40 I = 2, N

         IF( ABS( SX(I) ).LE.SMAX ) GO TO 40

         ISAMAX = I
         SMAX   = ABS( SX(I) )

   40 CONTINUE

      RETURN

      END


      SUBROUTINE SCOPY( N, SX, INCX, SY, INCY )

c  Level 1 BLAS routine.

c     Copies a vector, x, to a vector, y.
c     Uses unrolled loops for increments equal to 1.
c     Jack Dongarra, LINPACK, 3/11/78.

c     .. Scalar Arguments ..

      INTEGER   INCX, INCY, N
c     ..
c     .. Array Arguments ..

      REAL      SX( 1 ), SY( 1 )
c     ..
c     .. Local Scalars ..

      INTEGER   I, IX, IY, M, MP1
c     ..
c     .. Intrinsic Functions ..

      INTRINSIC MOD
c     ..

      IF( N.LE.0 ) RETURN

      IF( INCX.EQ.1 .AND. INCY.EQ.1 ) GO TO 20

c        code for unequal increments or equal increments
c        not equal to 1

      IX  = 1
      IY  = 1

      IF( INCX.LT.0 ) IX  = ( - N + 1 )*INCX + 1

      IF( INCY.LT.0 ) IY  = ( - N + 1 )*INCY + 1

      DO 10 I = 1, N

         SY( IY ) = SX( IX )
         IX  = IX + INCX
         IY  = IY + INCY

   10 CONTINUE

      RETURN

c        code for both increments equal to 1

c        clean-up loop

   20 CONTINUE
      M  = MOD( N, 7 )

      IF( M.EQ.0 ) GO TO 40

      DO 30 I = 1, M
         SY( I ) = SX( I )
   30 CONTINUE

      IF( N.LT.7 ) RETURN

   40 CONTINUE
      MP1  = M + 1

      DO 50 I = MP1, N, 7

         SY( I )     = SX( I )
         SY( I + 1 ) = SX( I + 1 )
         SY( I + 2 ) = SX( I + 2 )
         SY( I + 3 ) = SX( I + 3 )
         SY( I + 4 ) = SX( I + 4 )
         SY( I + 5 ) = SX( I + 5 )
         SY( I + 6 ) = SX( I + 6 )

   50 CONTINUE

      RETURN

      END


      SUBROUTINE SSCAL( N, SA, SX, INCX )

c  Level 1 Blas routine.

c     Scales a vector by a constant.
c     Uses unrolled loops for increment equal to 1.
c     Jack Dongarra, LINPACK, 3/11/78.

c     .. Scalar Arguments ..

      INTEGER   INCX, N
      REAL      SA
c     ..
c     .. Array Arguments ..

      REAL      SX( 1 )
c     ..
c     .. Local Scalars ..

      INTEGER   I, M, MP1, NINCX
c     ..
c     .. Intrinsic Functions ..

      INTRINSIC MOD
c     ..

      IF( N.LE.0 ) RETURN

      IF( INCX.EQ.1 ) GO TO 20

c        code for increment not equal to 1


      NINCX  = N*INCX

      DO 10 I = 1, NINCX, INCX
         SX( I ) = SA * SX( I )
   10 CONTINUE

      RETURN

c        code for increment equal to 1

c        clean-up loop

   20 CONTINUE
      M  = MOD( N, 5 )

      IF( M.EQ.0 ) GO TO 40

      DO 30 I = 1, M
         SX( I ) = SA * SX( I )
   30 CONTINUE

      IF( N.LT.5 ) RETURN

   40 CONTINUE
      MP1  = M + 1

      DO 50 I = MP1, N, 5

         SX( I )     = SA*SX( I )
         SX( I + 1 ) = SA*SX( I + 1 )
         SX( I + 2 ) = SA*SX( I + 2 )
         SX( I + 3 ) = SA*SX( I + 3 )
         SX( I + 4 ) = SA*SX( I + 4 )

   50 CONTINUE

      RETURN

      END


      SUBROUTINE SSWAP( N, SX, INCX, SY, INCY )

c  Level 1 Blas routine.

c     Interchanges two vectors.
c     Uses unrolled loops for increments equal to 1.
c     Jack Dongarra, LINPACK, 3/11/78.

c     .. Scalar Arguments ..

      INTEGER   INCX, INCY, N
c     ..
c     .. Array Arguments ..

      REAL      SX( 1 ), SY( 1 )
c     ..
c     .. Local Scalars ..

      INTEGER   I, IX, IY, M, MP1
      REAL      STEMP
c     ..
c     .. Intrinsic Functions ..

      INTRINSIC MOD
c     ..

      IF( N.LE.0 ) RETURN

      IF( INCX.EQ.1 .AND. INCY.EQ.1 ) GO TO 20

c         code for unequal increments or equal increments not equal to 1

      IX  = 1
      IY  = 1

      IF( INCX.LT.0 ) IX  = ( - N + 1 )*INCX + 1

      IF( INCY.LT.0 ) IY  = ( - N + 1 )*INCY + 1

      DO 10 I = 1, N

         STEMP    = SX( IX )
         SX( IX ) = SY( IY )
         SY( IY ) = STEMP
         IX  = IX + INCX
         IY  = IY + INCY

   10 CONTINUE

      RETURN

c       code for both increments equal to 1

c       clean-up loop

   20 CONTINUE
      M  = MOD( N, 3 )

      IF( M.EQ.0 ) GO TO 40

      DO 30 I = 1, M

         STEMP   = SX( I )
         SX( I ) = SY( I )
         SY( I ) = STEMP

   30 CONTINUE

      IF( N.LT.3 ) RETURN

   40 CONTINUE
      MP1  = M + 1

      DO 50 I = MP1, N, 3

         STEMP   = SX( I )
         SX( I ) = SY( I )
         SY( I ) = STEMP

         STEMP       = SX( I + 1 )
         SX( I + 1 ) = SY( I + 1 )
         SY( I + 1 ) = STEMP

         STEMP       = SX( I + 2 )
         SX( I + 2 ) = SY( I + 2 )
         SY( I + 2 ) = STEMP

   50 CONTINUE

      RETURN

      END


      LOGICAL FUNCTION LSAME( CA, CB )
c
c     September 30, 1994
c     Auxiliary routine for Level 2 Blas.
c
c     .. Scalar Arguments ..
      CHARACTER          CA, CB
c     ..
c
c  Purpose
c  =======
c
c  LSAME returns .TRUE. if CA is the same letter as CB regardless of
c  case.
c
c  Arguments
c  =========
c
c  CA      (input) CHARACTER*1
c  CB      (input) CHARACTER*1
c          CA and CB specify the single characters to be compared.
c
c =====================================================================

c     .. Intrinsic Functions ..
      INTRINSIC          ICHAR
c     ..
c     .. Local Scalars ..
      INTEGER            INTA, INTB, ZCODE
c     ..

c     Test if the characters are equal

      LSAME = CA.EQ.CB
      IF( LSAME )  RETURN

c     Now test for equivalence if both characters are alphabetic.

      ZCODE = ICHAR( 'Z' )

c     Use 'Z' rather than 'A' so that ASCII can be detected on Prime
c     machines, on which ICHAR returns a value with bit 8 set.
c     ICHAR('A') on Prime machines returns 193 which is the same as
c     ICHAR('A') on an EBCDIC machine.

      INTA = ICHAR( CA )
      INTB = ICHAR( CB )

      IF( ZCODE.EQ.90 .OR. ZCODE.EQ.122 ) THEN

c        ASCII is assumed - ZCODE is the ASCII code of either lower or
c        upper case 'Z'.

         IF( INTA.GE.97 .AND. INTA.LE.122 ) INTA = INTA - 32
         IF( INTB.GE.97 .AND. INTB.LE.122 ) INTB = INTB - 32

      ELSE IF( ZCODE.EQ.233 .OR. ZCODE.EQ.169 ) THEN

c        EBCDIC is assumed - ZCODE is the EBCDIC code of either lower or
c        upper case 'Z'.

         IF( INTA.GE.129 .AND. INTA.LE.137 .OR.
     &       INTA.GE.145 .AND. INTA.LE.153 .OR.
     &       INTA.GE.162 .AND. INTA.LE.169 ) INTA = INTA + 64
         IF( INTB.GE.129 .AND. INTB.LE.137 .OR.
     &       INTB.GE.145 .AND. INTB.LE.153 .OR.
     &       INTB.GE.162 .AND. INTB.LE.169 ) INTB = INTB + 64

      ELSE IF( ZCODE.EQ.218 .OR. ZCODE.EQ.250 ) THEN

c        ASCII is assumed, on Prime machines - ZCODE is the ASCII code
c        plus 128 of either lower or upper case 'Z'.

         IF( INTA.GE.225 .AND. INTA.LE.250 ) INTA = INTA - 32
         IF( INTB.GE.225 .AND. INTB.LE.250 ) INTB = INTB - 32
      END IF

      LSAME = INTA.EQ.INTB

      RETURN

c          End of LSAME
      END


      SUBROUTINE XERBLA( SRNAME, INFO )
c
c     September 30, 1994
c     Auxiliary routine for Level 2 Blas.
c
c     .. Scalar Arguments ..
      CHARACTER*(*)      SRNAME
      INTEGER            INFO
c     ..

c  Arguments
c  =========
c
c  SRNAME  (input) CHARACTER
c          The name of the routine which called XERBLA.
c
c  INFO    (input) INTEGER
c          The position of the invalid parameter in the parameter list
c          of the calling routine.
c
c =====================================================================


      WRITE( *, '(3A,I2,A)' ) ' ** On entry to ', SRNAME, 
     &   ' argument number ', INFO, ' had an illegal value'

      STOP

c         End of XERBLA
      END


