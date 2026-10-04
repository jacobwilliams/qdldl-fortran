!*****************************************************************************************
!>
!  Real and integer kinds of the QDLDL library, chosen by preprocessor macros:
!
!  * `REAL32`, `REAL64` (the default), or `REAL128` selects the real kind `wp`.
!  * `INT64` selects 64-bit integers for indices and counts (`ip`), needed when the
!    nonzeros of \(L\) can exceed \(2^{31}-1\). The default is 32-bit.
!
!  For example: `fpm test --flag "-DREAL32 -DINT64"`.

    module qdldl_kinds

    use, intrinsic :: iso_fortran_env, only: real32, real64, real128, int32, int64

    implicit none

    private

#ifdef REAL32
    integer,parameter,public :: wp = real32   !! Real working precision [4 bytes]
#elif REAL64
    integer,parameter,public :: wp = real64   !! Real working precision [8 bytes]
#elif REAL128
    integer,parameter,public :: wp = real128  !! Real working precision [16 bytes]
#else
    integer,parameter,public :: wp = real64   !! Real working precision if not specified [8 bytes]
#endif

#ifdef INT64
    integer,parameter,public :: ip = int64    !! Integer kind of indices and counts [8 bytes]
#else
    integer,parameter,public :: ip = int32    !! Integer kind of indices and counts if not specified [4 bytes]
#endif

    end module qdldl_kinds
!*****************************************************************************************
