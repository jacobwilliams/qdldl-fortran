!*****************************************************************************************
!>
!  Approximate minimum degree (AMD) ordering: a modernized version of the Fortran 77
!  `AMD` routine (`AMD/Source/amd.f`) of SuiteSparse.
!
!  The algorithm and its variable names are unchanged, so it can be compared with the
!  original line by line. Changes:
!
!  * Free-form source in a module, with `intent` on every argument, and the kind
!    `ip` of [[qdldl_kinds]] for every integer (so it runs with 64-bit indices).
!  * Structured loops (`do while`, `exit`, `cycle`) instead of `GOTO`.
!  * The tests for integer overflow no longer rely on wrap-around, which is undefined
!    in Fortran: `wflg + n <= wflg` becomes `wflg > huge(wflg) - n`, and the hash
!    of a supervariable is reduced modulo `hmod` as it is summed (the same value as
!    the original's `mod(hash, hmod)`, without overflow for large `n`).
!
!### License
!
!  AMD, Copyright (c) 1996-2022, Timothy A. Davis, Patrick R. Amestoy, and
!  Iain S. Duff. All Rights Reserved. SPDX-License-Identifier: BSD-3-clause.
!  Authors of the original Fortran version, and Copyright (C) 1995 by:
!  Timothy A. Davis, Patrick Amestoy, Iain S. Duff, & John K. Reid.
!  Changed from the original: translated to modern Fortran, as described above.
!
!  Redistribution and use in source and binary forms, with or without modification,
!  are permitted provided that the following conditions are met: (1) redistributions
!  of source code must retain the above copyright notice, this list of conditions and
!  the following disclaimer; (2) redistributions in binary form must reproduce the
!  above copyright notice, this list of conditions and the following disclaimer in the
!  documentation and/or other materials provided with the distribution; (3) neither
!  the name of the organizations to which the authors are affiliated, nor the names
!  of its contributors may be used to endorse or promote products derived from this
!  software without specific prior written permission. THIS SOFTWARE IS PROVIDED BY
!  THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED
!  WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
!  MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT
!  SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
!  INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
!  LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
!  PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
!  WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
!  ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
!  POSSIBILITY OF SUCH DAMAGE.

    module qdldl_amd

    use qdldl_kinds, only: ip

    implicit none

    private

    public :: amd

    contains
!*****************************************************************************************

!*****************************************************************************************
!>
!  Given a representation of the nonzero pattern of a symmetric matrix, \(A\)
!  (excluding the diagonal), performs an approximate minimum (UMFPACK/MA38-style)
!  degree ordering to compute a pivot order such that the introduction of nonzeros
!  (fill-in) in the Cholesky factors \(A = LL^T\) are kept low. At each step, the
!  pivot selected is the one with the minimum UMFPACK/MA38-style upper-bound on the
!  external degree. Aggressive absorption is used to tighten the bound on the degree.
!
!  **The arguments are not checked for errors on input.**
!
!  The input pattern is held in `iw(1:pfree-1)`: row `i` (`len(i)` column indices,
!  excluding the diagonal, in any order, without duplicates) is
!  `iw(pe(i):pe(i)+len(i)-1)`. Both triangles must be present. There may be empty
!  space between the rows. `iwlen >= pfree + n` is recommended (better yet,
!  `iwlen > 1.2*pfree`), or excessive compressions take place; the algorithm does
!  not run at all if `iwlen < pfree - 1`.
!
!  On output, `last` holds the permutation (row `last(k)` of \(A\) is the `k`-th
!  row of \(PAP^T\)) and `elen` its inverse.
!
!### References
!
!  * T. A. Davis and I. S. Duff, "An unsymmetric-pattern multifrontal method for
!    sparse LU factorization", SIAM J. Matrix Anal. Appl. 18(1), 140-158.
!  * P. R. Amestoy, T. A. Davis, and I. S. Duff, "An approximate minimum degree
!    ordering algorithm", SIAM J. Matrix Anal. Appl. 17(4), 886-905, 1996.
!  * A. George and J. Liu, "The evolution of the minimum degree ordering algorithm",
!    SIAM Review 31(1), 1-19, 1989.
!
!### Data structures
!
!  During execution `pe` holds, for a principal supervariable `i`, the index in `iw`
!  of its description; for a non-principal supervariable absorbed into `j`, `-j`; for
!  an unabsorbed element `e`, the index in `iw` of its description; for an element
!  absorbed into `e2`, `-e2` (or `0` for a null element). On output it holds the
!  assembly tree. The list of supervariable `i` is: its elements
!  `iw(pe(i):pe(i)+elen(i)-1)`, then its supervariables up to `pe(i)+len(i)-1`. The
!  list of element `e` is its supervariables `iw(pe(e):pe(e)+len(e)-1)`.
!  `abs(nv(i))` is the number of rows represented by principal supervariable `i` (`0`
!  if non-principal; negative while `i` is in the pattern `Lme` of the current pivot
!  element `me`). `degree(i)` is the approximate external degree of supervariable
!  `i`, or `|Le|` for an element. `head`, `next`, and `last` hold the degree lists
!  and the hash buckets. `w` flags elements and variables: for an element `e`,
!  `w(e) = 0` if absorbed, `w(e) - wflg = |Le \ Lme|` if `w(e) >= wflg`.

    subroutine amd(n, pe, iw, len, iwlen, pfree, nv, next, last, head, elen, degree, ncmpa, w)

    integer(ip),intent(in)    :: n             !! the matrix order (`1 <= n < huge/2 - 2`)
    integer(ip),intent(inout) :: pe(n)         !! on input, the start of each row in `iw` (or 0 if empty);
                                               !! on output, the assembly tree
    integer(ip),intent(in)    :: iwlen         !! the length of `iw`
    integer(ip),intent(inout) :: iw(iwlen)     !! on input, the pattern; undefined on output
    integer(ip),intent(inout) :: len(n)        !! on input, the number of entries of each row; undefined on output
    integer(ip),intent(inout) :: pfree         !! on input, `iw(pfree:iwlen)` is empty; on output, the size of `iw`
                                               !! that would have needed no compressions
    integer(ip),intent(out)   :: nv(n)         !! on output, `nv(e)` is the degree of element `e` when it was created
    integer(ip),intent(out)   :: next(n)       !! work array: degree lists and hash buckets
    integer(ip),intent(out)   :: last(n)       !! on output, the permutation
    integer(ip),intent(out)   :: head(n)       !! work array: heads of the degree lists and hash buckets
    integer(ip),intent(out)   :: elen(n)       !! on output, the inverse permutation
    integer(ip),intent(out)   :: degree(n)     !! work array: approximate degrees
    integer(ip),intent(out)   :: ncmpa         !! the number of times `iw` was compressed
    integer(ip),intent(out)   :: w(n)          !! work array: flags

    ! local integers (see the original for their full description)
    integer(ip) :: deg      ! the degree of a variable or element
    integer(ip) :: degme    ! size, |Lme|, of the current element, me (= degree (me))
    integer(ip) :: dext     ! external degree, |Le \ Lme|, of some element e
    integer(ip) :: dmax     ! largest |Le| seen so far
    integer(ip) :: e        ! an element
    integer(ip) :: elenme   ! the length, elen (me), of element list of pivotal var.
    integer(ip) :: eln      ! the length, elen (...), of an element list
    integer(ip) :: hash     ! the computed value of the hash function
    integer(ip) :: hmod     ! the hash function is computed modulo hmod = max (1,n-1)
    integer(ip) :: i        ! a supervariable
    integer(ip) :: ilast    ! the entry in a link list preceding i
    integer(ip) :: inext    ! the entry in a link list following i
    integer(ip) :: j        ! a supervariable
    integer(ip) :: jlast    ! the entry in a link list preceding j
    integer(ip) :: jnext    ! the entry in a link list, or path, following j
    integer(ip) :: k        ! the pivot order of an element or variable
    integer(ip) :: knt1     ! loop counter used during element construction
    integer(ip) :: knt2     ! loop counter used during element construction
    integer(ip) :: knt3     ! loop counter used during compression
    integer(ip) :: lenj     ! len (j)
    integer(ip) :: ln       ! length of a supervariable list
    integer(ip) :: maxmem   ! amount of memory needed for no compressions
    integer(ip) :: me       ! current supervariable being eliminated, and the current element
    integer(ip) :: mem      ! memory in use assuming no compressions have occurred
    integer(ip) :: mindeg   ! current minimum degree
    integer(ip) :: nel      ! number of pivots selected so far
    integer(ip) :: newmem   ! amount of new memory needed for current pivot element
    integer(ip) :: nleft    ! n - nel, the number of nonpivotal rows/columns remaining
    integer(ip) :: nvi      ! the number of variables in a supervariable i (= nv (i))
    integer(ip) :: nvj      ! the number of variables in a supervariable j (= nv (j))
    integer(ip) :: nvpiv    ! number of pivots in current element
    integer(ip) :: slenme   ! number of variables in variable list of pivotal variable
    integer(ip) :: we       ! w (e)
    integer(ip) :: wflg     ! used for flagging the w array
    integer(ip) :: wnvi     ! wflg - nv (i)
    logical     :: absorb   ! j can be absorbed into i (supervariable detection)

    ! local pointers (indices into iw)
    integer(ip) :: p        ! pointer into lots of things
    integer(ip) :: p1       ! pe (i) for some variable i (start of element list)
    integer(ip) :: p2       ! pe (i) + elen (i) -  1 for some var. i (end of el. list)
    integer(ip) :: p3       ! index of first supervariable in clean list
    integer(ip) :: pdst     ! destination pointer, for compression
    integer(ip) :: pend     ! end of memory to compress
    integer(ip) :: pj       ! pointer into an element or variable
    integer(ip) :: pme      ! pointer into the current element (pme1...pme2)
    integer(ip) :: pme1     ! the current element, me, is stored in iw (pme1...pme2)
    integer(ip) :: pme2     ! the end of the current element
    integer(ip) :: pn       ! pointer into a "clean" variable, also used to compress
    integer(ip) :: psrc     ! source pointer, for compression

    !=======================================================================
    !  INITIALIZATIONS
    !=======================================================================

    wflg = 2
    mindeg = 1
    ncmpa = 0
    nel = 0
    hmod = max(1_ip, n-1)
    dmax = 0
    mem = pfree - 1
    maxmem = mem
    me = 0

    do i = 1, n
        last(i) = 0
        head(i) = 0
        nv(i) = 1
        w(i) = 1
        elen(i) = 0
        degree(i) = len(i)
    end do

    ! ----------------------------------------------------------------
    ! initialize degree lists and eliminate rows with no off-diag. nz.
    ! ----------------------------------------------------------------

    do i = 1, n
        deg = degree(i)
        if (deg > 0) then
            ! place i in the degree list corresponding to its degree
            inext = head(deg)
            if (inext /= 0) last(inext) = i
            next(i) = inext
            head(deg) = i
        else
            ! we have a variable that can be eliminated at once because
            ! there is no off-diagonal non-zero in its row.
            nel = nel + 1
            elen(i) = -nel
            pe(i) = 0
            w(i) = 0
        end if
    end do

    !=======================================================================
    !  WHILE (selecting pivots) DO
    !=======================================================================

    do while (nel < n)

        !=======================================================================
        !  GET PIVOT OF MINIMUM DEGREE
        !=======================================================================

        ! find next supervariable for elimination
        do deg = mindeg, n
            me = head(deg)
            if (me > 0) exit
        end do
        mindeg = deg

        ! remove chosen variable from link list
        inext = next(me)
        if (inext /= 0) last(inext) = 0
        head(deg) = inext

        ! me represents the elimination of pivots nel+1 to nel+nv(me).
        ! place me itself as the first in this set.  It will be moved
        ! to the nel+nv(me) position when the permutation vectors are
        ! computed.
        elenme = elen(me)
        elen(me) = -(nel + 1)
        nvpiv = nv(me)
        nel = nel + nvpiv

        !=======================================================================
        !  CONSTRUCT NEW ELEMENT
        !=======================================================================

        ! At this point, me is the pivotal supervariable.  It will be
        ! converted into the current element.  Scan list of the
        ! pivotal supervariable, me, setting tree pointers and
        ! constructing new list of supervariables for the new element,
        ! me.  p is a pointer to the current position in the old list.

        ! flag the variable "me" as being in Lme by negating nv (me)
        nv(me) = -nvpiv
        degme = 0

        if (elenme == 0) then

            ! construct the new element in place
            pme1 = pe(me)
            pme2 = pme1 - 1

            do p = pme1, pme1 + len(me) - 1
                i = iw(p)
                nvi = nv(i)
                if (nvi > 0) then
                    ! i is a principal variable not yet placed in Lme.
                    ! store i in new list
                    degme = degme + nvi
                    ! flag i as being in Lme by negating nv (i)
                    nv(i) = -nvi
                    pme2 = pme2 + 1
                    iw(pme2) = i
                    ! remove variable i from degree list.
                    ilast = last(i)
                    inext = next(i)
                    if (inext /= 0) last(inext) = ilast
                    if (ilast /= 0) then
                        next(ilast) = inext
                    else
                        ! i is at the head of the degree list
                        head(degree(i)) = inext
                    end if
                end if
            end do
            ! this element takes no new memory in iw:
            newmem = 0

        else

            ! construct the new element in empty space, iw (pfree ...)
            p = pe(me)
            pme1 = pfree
            slenme = len(me) - elenme

            do knt1 = 1, elenme + 1

                if (knt1 > elenme) then
                    ! search the supervariables in me.
                    e = me
                    pj = p
                    ln = slenme
                else
                    ! search the elements in me.
                    e = iw(p)
                    p = p + 1
                    pj = pe(e)
                    ln = len(e)
                end if

                ! search for different supervariables and add them to the
                ! new list, compressing when necessary. this loop is
                ! executed once for each element in the list and once for
                ! all the supervariables in the list.
                do knt2 = 1, ln
                    i = iw(pj)
                    pj = pj + 1
                    nvi = nv(i)
                    if (nvi > 0) then

                        ! compress iw, if necessary
                        if (pfree > iwlen) then
                            ! prepare for compressing iw by adjusting
                            ! pointers and lengths so that the lists being
                            ! searched in the inner and outer loops contain
                            ! only the remaining entries.
                            pe(me) = p
                            len(me) = len(me) - knt1
                            ! nothing left of supervariable me
                            if (len(me) == 0) pe(me) = 0
                            pe(e) = pj
                            len(e) = ln - knt2
                            ! nothing left of element e
                            if (len(e) == 0) pe(e) = 0

                            ncmpa = ncmpa + 1
                            ! store first item in pe
                            ! set first entry to -item
                            do j = 1, n
                                pn = pe(j)
                                if (pn > 0) then
                                    pe(j) = iw(pn)
                                    iw(pn) = -j
                                end if
                            end do

                            ! psrc/pdst point to source/destination
                            pdst = 1
                            psrc = 1
                            pend = pme1 - 1

                            do while (psrc <= pend)
                                ! search for next negative entry
                                j = -iw(psrc)
                                psrc = psrc + 1
                                if (j > 0) then
                                    iw(pdst) = pe(j)
                                    pe(j) = pdst
                                    pdst = pdst + 1
                                    ! copy from source to destination
                                    lenj = len(j)
                                    do knt3 = 0, lenj - 2
                                        iw(pdst + knt3) = iw(psrc + knt3)
                                    end do
                                    pdst = pdst + lenj - 1
                                    psrc = psrc + lenj - 1
                                end if
                            end do

                            ! move the new partially-constructed element
                            p1 = pdst
                            do psrc = pme1, pfree - 1
                                iw(pdst) = iw(psrc)
                                pdst = pdst + 1
                            end do
                            pme1 = p1
                            pfree = pdst
                            pj = pe(e)
                            p = pe(me)
                        end if

                        ! i is a principal variable not yet placed in Lme
                        ! store i in new list
                        degme = degme + nvi
                        ! flag i as being in Lme by negating nv (i)
                        nv(i) = -nvi
                        iw(pfree) = i
                        pfree = pfree + 1

                        ! remove variable i from degree link list
                        ilast = last(i)
                        inext = next(i)
                        if (inext /= 0) last(inext) = ilast
                        if (ilast /= 0) then
                            next(ilast) = inext
                        else
                            ! i is at the head of the degree list
                            head(degree(i)) = inext
                        end if

                    end if
                end do

                if (e /= me) then
                    ! set tree pointer and flag to indicate element e is
                    ! absorbed into new element me (the parent of e is me)
                    pe(e) = -me
                    w(e) = 0
                end if
            end do

            pme2 = pfree - 1
            ! this element takes newmem new memory in iw (possibly zero)
            newmem = pfree - pme1
            mem = mem + newmem
            maxmem = max(maxmem, mem)
        end if

        ! me has now been converted into an element in iw (pme1..pme2)

        ! degme holds the external degree of new element
        degree(me) = degme
        pe(me) = pme1
        len(me) = pme2 - pme1 + 1

        ! make sure that wflg is not too large.  With the current
        ! value of wflg, wflg+n must not cause integer overflow
        if (wflg > huge(wflg) - n) call reset_w()

        !=======================================================================
        !  COMPUTE (w (e) - wflg) = |Le\Lme| FOR ALL ELEMENTS
        !=======================================================================

        ! Scan 1:  compute the external degrees of previous elements
        ! with respect to the current element.  That is:
        !      (w (e) - wflg) = |Le \ Lme|
        ! for each element e that appears in any supervariable in Lme.
        ! If (w (e) - wflg) becomes zero, then the element e will be
        ! absorbed in scan 2.

        do pme = pme1, pme2
            i = iw(pme)
            eln = elen(i)
            if (eln > 0) then
                ! note that nv (i) has been negated to denote i in Lme:
                nvi = -nv(i)
                wnvi = wflg - nvi
                do p = pe(i), pe(i) + eln - 1
                    e = iw(p)
                    we = w(e)
                    if (we >= wflg) then
                        ! unabsorbed element e has been seen in this loop
                        we = we - nvi
                    else if (we /= 0) then
                        ! e is an unabsorbed element
                        ! this is the first we have seen e in all of Scan 1
                        we = degree(e) + wnvi
                    end if
                    w(e) = we
                end do
            end if
        end do

        !=======================================================================
        !  DEGREE UPDATE AND ELEMENT ABSORPTION
        !=======================================================================

        ! Scan 2:  for each i in Lme, sum up the degree of Lme (which
        ! is degme), plus the sum of the external degrees of each Le
        ! for the elements e appearing within i, plus the
        ! supervariables in i.  Place i in hash list.

        do pme = pme1, pme2
            i = iw(pme)
            p1 = pe(i)
            p2 = p1 + elen(i) - 1
            pn = p1
            hash = 0
            deg = 0

            ! scan the element list associated with supervariable i
            do p = p1, p2
                e = iw(p)
                ! dext = | Le \ Lme |
                dext = w(e) - wflg
                if (dext > 0) then
                    deg = deg + dext
                    iw(pn) = e
                    pn = pn + 1
                    hash = mod(hash + e, hmod)
                else if (dext == 0) then
                    ! aggressive absorption: e is not adjacent to me, but
                    ! the |Le \ Lme| is 0, so absorb it into me
                    pe(e) = -me
                    w(e) = 0
                end if
                ! (else: element e has already been absorbed, due to
                ! regular absorption, in the element construction. Ignore it.)
            end do

            ! count the number of elements in i (including me):
            elen(i) = pn - p1 + 1

            ! scan the supervariables in the list associated with i
            p3 = pn
            do p = p2 + 1, p1 + len(i) - 1
                j = iw(p)
                nvj = nv(j)
                if (nvj > 0) then
                    ! j is unabsorbed, and not in Lme.
                    ! add to degree and add to new list
                    deg = deg + nvj
                    iw(pn) = j
                    pn = pn + 1
                    hash = mod(hash + j, hmod)
                end if
            end do

            ! update the degree and check for mass elimination
            if (deg == 0) then

                ! mass elimination: there is nothing left of this node except
                ! for an edge to the current pivot element.  elen (i) is 1,
                ! and there are no variables adjacent to node i.
                ! Absorb i into the current pivot element, me.
                pe(i) = -me
                nvi = -nv(i)
                degme = degme - nvi
                nvpiv = nvpiv + nvi
                nel = nel + nvi
                nv(i) = 0
                elen(i) = 0

            else

                ! update the upper-bound degree of i
                ! the following degree does not yet include the size
                ! of the current element, which is added later:
                degree(i) = min(degree(i), deg)

                ! add me to the list for i
                ! move first supervariable to end of list
                iw(pn) = iw(p3)
                ! move first element to end of element part of list
                iw(p3) = iw(p1)
                ! add new element to front of list.
                iw(p1) = me
                ! store the new length of the list in len (i)
                len(i) = pn - p1 + 1

                ! place in hash bucket.  Save hash key of i in last (i).
                hash = hash + 1
                j = head(hash)
                if (j <= 0) then
                    ! the degree list is empty, hash head is -j
                    next(i) = -j
                    head(hash) = -i
                else
                    ! degree list is not empty
                    ! use last (head (hash)) as hash head
                    next(i) = last(j)
                    last(j) = i
                end if
                last(i) = hash
            end if
        end do

        degree(me) = degme

        ! Clear the counter array, w (...), by incrementing wflg.
        dmax = max(dmax, degme)
        wflg = wflg + dmax

        ! make sure that wflg+n does not cause integer overflow
        if (wflg > huge(wflg) - n) call reset_w()
        ! at this point, w (1..n) .lt. wflg holds

        !=======================================================================
        !  SUPERVARIABLE DETECTION
        !=======================================================================

        do pme = pme1, pme2
            i = iw(pme)
            if (nv(i) < 0) then
                ! i is a principal variable in Lme

                ! examine all hash buckets with 2 or more variables.  We
                ! do this by examing all unique hash keys for super-
                ! variables in the pattern Lme of the current element, me
                hash = last(i)
                ! let i = head of hash bucket, and empty the hash bucket
                j = head(hash)
                if (j == 0) cycle
                if (j < 0) then
                    ! degree list is empty
                    i = -j
                    head(hash) = 0
                else
                    ! degree list is not empty, restore last () of head
                    i = last(j)
                    last(j) = 0
                end if
                if (i == 0) cycle

                do while (next(i) /= 0)

                    ! this bucket has one or more variables following i.
                    ! scan all of them to see if i can absorb any entries
                    ! that follow i in hash bucket.  Scatter i into w.
                    ln = len(i)
                    eln = elen(i)
                    ! do not flag the first element in the list (me)
                    do p = pe(i) + 1, pe(i) + ln - 1
                        w(iw(p)) = wflg
                    end do

                    ! scan every other entry j following i in bucket
                    jlast = i
                    j = next(i)

                    do while (j /= 0)

                        ! check if j and i have identical nonzero pattern
                        ! (the same size data structure, the same number of
                        ! adjacent elements, and the same entries)
                        absorb = len(j) == ln .and. elen(j) == eln
                        if (absorb) then
                            ! do not flag the first element in the list (me)
                            do p = pe(j) + 1, pe(j) + ln - 1
                                if (w(iw(p)) /= wflg) then
                                    ! an entry (iw(p)) is in j but not in i
                                    absorb = .false.
                                    exit
                                end if
                            end do
                        end if

                        if (absorb) then
                            ! found it!  j can be absorbed into i
                            pe(j) = -i
                            ! both nv (i) and nv (j) are negated since they
                            ! are in Lme, and the absolute values of each
                            ! are the number of variables in i and j:
                            nv(i) = nv(i) + nv(j)
                            nv(j) = 0
                            elen(j) = 0
                            ! delete j from hash bucket
                            j = next(j)
                            next(jlast) = j
                        else
                            ! j cannot be absorbed into i
                            jlast = j
                            j = next(j)
                        end if
                    end do

                    ! no more variables can be absorbed into i
                    ! go to next i in bucket and clear flag array
                    wflg = wflg + 1
                    i = next(i)
                    if (i == 0) exit
                end do
            end if
        end do

        !=======================================================================
        !  RESTORE DEGREE LISTS AND REMOVE NONPRINCIPAL SUPERVAR. FROM ELEMENT
        !=======================================================================

        p = pme1
        nleft = n - nel
        do pme = pme1, pme2
            i = iw(pme)
            nvi = -nv(i)
            if (nvi > 0) then
                ! i is a principal variable in Lme
                ! restore nv (i) to signify that i is principal
                nv(i) = nvi

                ! compute the external degree (add size of current elem)
                deg = min(degree(i) + degme - nvi, nleft - nvi)

                ! place the supervariable at the head of the degree list
                inext = head(deg)
                if (inext /= 0) last(inext) = i
                next(i) = inext
                last(i) = 0
                head(deg) = i

                ! save the new degree, and find the minimum degree
                mindeg = min(mindeg, deg)
                degree(i) = deg

                ! place the supervariable in the element pattern
                iw(p) = i
                p = p + 1
            end if
        end do

        !=======================================================================
        !  FINALIZE THE NEW ELEMENT
        !=======================================================================

        nv(me) = nvpiv + degme
        ! nv (me) is now the degree of pivot (including diagonal part)
        ! save the length of the list for the new element me
        len(me) = p - pme1
        if (len(me) == 0) then
            ! there is nothing left of the current pivot element
            pe(me) = 0
            w(me) = 0
        end if
        if (newmem /= 0) then
            ! element was not constructed in place: deallocate part
            ! of it (final size is less than or equal to newmem,
            ! since newly nonprincipal variables have been removed).
            pfree = p
            mem = mem - newmem + len(me)
        end if

    end do   ! END WHILE (selecting pivots)

    !=======================================================================
    !  COMPUTE THE PERMUTATION VECTORS
    !=======================================================================

    ! The time taken by the following code is O(n).  At this
    ! point, elen (e) = -k has been done for all elements e,
    ! and elen (i) = 0 has been done for all nonprincipal
    ! variables i.  At this point, there are no principal
    ! supervariables left, and all elements are absorbed.

    ! compute the ordering of unordered nonprincipal variables
    do i = 1, n
        if (elen(i) == 0) then

            ! i is an un-ordered row.  Traverse the tree from i until
            ! reaching an element, e.  The element, e, was the
            ! principal supervariable of i and all nodes in the path
            ! from i to when e was selected as pivot.
            j = -pe(i)
            ! while (j is a variable) do:
            do while (elen(j) >= 0)
                j = -pe(j)
            end do
            e = j

            ! get the current pivot ordering of e
            k = -elen(e)

            ! traverse the path again from i to e, and compress the
            ! path (all nodes point to e).  Path compression allows
            ! this code to compute in O(n) time.  Order the unordered
            ! nodes in the path, and place the element e at the end.
            j = i
            ! while (j is a variable) do:
            do while (elen(j) >= 0)
                jnext = -pe(j)
                pe(j) = -e
                if (elen(j) == 0) then
                    ! j is an unordered row
                    elen(j) = k
                    k = k + 1
                end if
                j = jnext
            end do
            ! leave elen (e) negative, so we know it is an element
            elen(e) = -k
        end if
    end do

    ! reset the inverse permutation (elen (1..n)) to be positive,
    ! and compute the permutation (last (1..n)).
    do i = 1, n
        k = abs(elen(i))
        last(k) = i
        elen(i) = k
    end do

    !=======================================================================
    !  RETURN THE MEMORY USAGE IN IW
    !=======================================================================

    ! If maxmem is less than or equal to iwlen, then no compressions
    ! occurred, and iw (maxmem+1 ... iwlen) was unused.  Otherwise
    ! compressions did occur, and iwlen would have had to have been
    ! greater than or equal to maxmem for no compressions to occur.
    ! Return the value of maxmem in the pfree argument.
    pfree = maxmem

    contains

        subroutine reset_w()
        !! Resets the flags of the unabsorbed elements and variables, and `wflg`.
        integer(ip) :: x
        do x = 1, n
            if (w(x) /= 0) w(x) = 1
        end do
        wflg = 2
        end subroutine reset_w

    end subroutine amd
!*****************************************************************************************

!*****************************************************************************************
    end module qdldl_amd
!*****************************************************************************************
