c
c     Maintainer: Steffen Genkinger 
c                 genkinger@statik.uni-stuttgart.de 
c                 http://www.uni-stuttgart.de/ibs/members/genkinger/ 
c                 0711 - 685-6127 
c
c     ---------------------------------------------------------------  
	subroutine qat2v2(pcell, wcell, kcell, wnode, knode, wface, 
     &                    mtri, mptri, medge, mpedge, mface, mpface,
     &                    ptri, ctri, itri, ktri, pptri, kptri,
     &                    pedge, cedge, kedge, ppedge, kpedge,
     &                    pface, kface, ppface, kpface)
c
c	Make PolyTri and PolyLine strips from random quad and tri cell
c	     data.  The result is suitable to use with Visual2.
c
	integer pcell(4,kcell)
	integer wnode(knode), wface(6,mface), wcell(kcell)
	integer ptri(mtri), ctri(mtri), itri(mtri), pptri(mptri)
	integer pedge(medge), cedge(medge), ppedge(mpedge)
	integer pface(mface), ppface(mpface)
c
c       all the following arguments are either work arrays or required
c          input to this routine:
c
c	where:	pcell  - cell pointers.  Each entry contains the node
c			   numbers that make up the cell.  For triangular
c			   cells pcell(1,n) must equal pcell(4,n).  For
c			   quad cells the pointers must be filled either
c			   clockwise or counter-clockwise (i.e. the 
c			   diagonals must be pcell(1,n)-pcell(3,n) and
c			   pcell(2,n)-pcell(4,n).  This array must be 
c			   atleast kcell in length. 
c		wcell  - work array (atleast kcell in length).
c		kcell  - the total number of cells.
c		wnode  - work array (atleast knode in length).
c		knode  - the total number of nodes used.
c		wface  - work array (atleast 6*mface in length).
c		mtri   - dimension of the PolyTri vectors.
c		mptri  - dimension of the pointer PolyTri vector (pptri).
c		medge  - dimension of the PolyLine Edge vectors.
c		mpedge - dimension of ppedge.
c		mface  - dimension of the PolyLine Face vectors.
c		mpface - dimension of ppface.
c
c	all the following arguments are returned, filled with Visual2
c	   structures:
c
c	where:	ptri   - PolyTri strip array.
c		ctri   - PolyTri connection array.
c		itri   - PolyTri cell number array.
c		ktri   - number of elements used in the PolyTri arrays.
c		pptri  - individual strip pointers into the PolyTri
c			   arrays.
c		kptri  - the total number of elements used in pptri.
c			   i.e. the number of strips.
c		pedge  - PolyLine Edge strip array.
c		cedge  - PolyLine Edge connection array.  Note: this is
c			   not a Visual2 requirement, but may be useful
c			   in patching together regions.
c		kedge  - number of elements used in the edge arrays.
c		ppedge - individual strip pointers into the PolyLine Edge
c			   arrays.
c		kpedge - the total number of elements used in ppedge.
c		pface  - PolyLine Face strip array.
c		kface  - number of elements used in pface.
c		ppface - individual strip pointers into the PolyLine Face
c			   arrays.
c		kpface - the total number of elements used in ppface.
c
c       Copyright 1990, Massachusetts Institute of Technology.
c
	common /blkvis/ kc, kfill
	integer indx(3,3), scell(100), kmax(2)
	logical rev
c
	data indx/1, 2, 3,  2, 3, 1,  3, 1, 2/
c
c	zero counters
c
	ktri   = 0
	kptri  = 0
	kedge  = 0
	kpedge = 0
	kface  = 0
	kpface = 0
c
c	zero hash table
c
	do 1 kn = 1, knode
	  wnode(kn) = 0
 1	continue
c
c	fill work face array:
c	  1	  2		3	4		5	6
c	node#1  node#2 	     cell#   cell/edge#      next node entrys
c
	kfill = 0
	do 2 kc = 1, kcell
	  k1 = pcell(1,kc)
	  k2 = pcell(2,kc)
	  k3 = pcell(3,kc)
	  k4 = pcell(4,kc)
	  call setfac(k1, k2, wnode, wface, mface)
	  if(k1 .eq. k4) then
            call setfac(k1, k3, wnode, wface, mface)
            call setfac(k2 ,k3, wnode ,wface, mface)
	  else
	    call setfac(k1, k4, wnode, wface, mface)
	    call setfac(k2, k3, wnode, wface, mface)
	    call setfac(k3, k4, wnode, wface, mface)
	  endif
 2	continue
c
c	fill in edge structures
c
	do 3 kf = 1, kfill
	  if(wface(4,kf) .eq. 0) then
	    k1 = wface(1,kf)
	    k2 = wface(2,kf)
c
c	    trace back to beginning of edge
	    kl = kf
	    kb = k1
	    kw = wnode(k1)
 31	    ks = kw
	    n = 5
	    if(wface(1,ks) .ne. kb) n = 6
	    kw = wface(n,ks)
	    if(ks.ne.kl .and. ks.ne.kf .and. wface(4,ks).eq.0) then
	      kl = ks
	      kb = wface(7-n,ks)
	      kw = wnode(kb)
	      go to 31
	    endif
	    if(ks .eq. kf .and. kb .eq. k2) go to 32
	    if(kw .ne. 0) go to 31
	    k1 = kb
c
c	    go foward and fill up PolyLine
 32	    kedge = kedge + 1
	    if(kedge .gt. medge) then
	      write(*,*) ' ERROR - MEDGE too small for fill'
	      stop
	    endif
	    pedge(kedge) = k1
	    cedge(kedge) = -ks
	    wface(4,ks) = -kedge
	    if(k1 .eq. wface(1,ks)) then
	      k1 = wface(2,ks)
	    else
	      k1 = wface(1,ks)
	    endif
	    kw = wnode(k1)
 33	    ks = kw
	    n = 5
	    if(wface(1,ks) .ne. k1) n = 6
	    kw = wface(n,ks)
	    if(wface(4,ks) .eq. 0) then
	      kedge = kedge + 1
	      if(kedge .gt. medge) then
	        write(*,*) ' ERROR - MEDGE too small for fill'
	        stop
	      endif
	      pedge(kedge) = k1
	      cedge(kedge) = -ks
	      wface(4,ks) = -kedge
	      k1 = wface(7-n,ks)
	      kw = wnode(k1)
	    endif
	    if(kw .ne. 0) go to 33
	    kedge = kedge + 1
	    if(kedge .gt. medge) then
	      write(*,*) ' ERROR - MEDGE too small for fill'
	      stop
	    endif
	    pedge(kedge) = k1
	    cedge(kedge) = 0
	    kpedge = kpedge + 1
	    if(kpedge .gt. mpedge) then
	      write(*,*) ' ERROR - MPEDGE too small for fill'
	      stop
	    endif
	    ppedge(kpedge) = kedge
	  endif
 3	continue
c
c	fill in PolyTri structures
c
	do 4 kc = 1, kcell
	  wcell(kc) = 0
 4	continue
c
c	start at edges
	do 5 ke = 1, kedge
	  kf = -cedge(ke)
	  if(kf .le. 0) go to 5
	  kc = wface(3,kf)
	  if(wcell(kc) .ne. 0) go to 5
	  k1 = wface(2,kf)
	  k2 = wface(1,kf)
	  ksav = ktri
	  itry = 0
 51	  itry = itry + 1
	  if(itry .eq. 4) go to 56
	  if(itry .eq. 3) then
	    if(icell .ge. 100) go to 56
	    if(ktri .ge. kmaxx) go to 56
	    kf = -cedge(ke)
	    k1 = wface(2,kf)
	    k2 = wface(1,kf)
	  endif
	  if(itry .eq. 2) then
	    kmaxx = ktri
	    if(icell .ge. 100) go to 56
	    kf = -cedge(ke)
	    k1 = wface(1,kf)
	    k2 = wface(2,kf)
	  endif
	  if(itry .eq. 2 .or. itry .eq. 3) then
	    kc = wface(3,kf)
	    do 52 i = 1, icell
	      ks = scell(i)
	      wcell(ks) = 0
 52	    continue
	    ktri = ksav
	  endif
c
	  ktri = ktri + 1
	  wcell(kc) = ktri + 1
	  icell = 1
	  scell(icell) = kc
          if(ktri .gt. mtri) then
            write(*,*) ' ERROR - MTRI too small for fill'
            stop
          endif
	  itri(ktri) = 0
	  ctri(ktri) = 0
	  ptri(ktri) = k1
          ktri = ktri + 1
          if(ktri .gt. mtri) then
            write(*,*) ' ERROR - MTRI too small for fill'
            stop
          endif
          itri(ktri) = kc
          ptri(ktri) = k2
	  ctri(ktri) = 0
 53	  if(pcell(1,kc) .eq. pcell(4,kc)) then
	    kn = pcell(1,kc)
	    if(kn .eq. k2 .or. kn .eq. k1) kn = pcell(2,kc)
	    if(kn .eq. k2 .or. kn .eq. k1) kn = pcell(3,kc)
	  else
	    kk = 1
	    if(pcell(4,kc) .eq. k2) kk = 2
	    if(pcell(1,kc) .eq. k2) kk = 3
	    if(pcell(2,kc) .eq. k2) kk = 4
	    k3 = pcell(kk,kc)
            ktri = ktri + 1
            if(ktri .gt. mtri) then
              write(*,*) ' ERROR - MTRI too small for fill'
              stop
            endif
            itri(ktri) = kc
            ptri(ktri) = k3
	    ctri(ktri) = 0
	    kn = pcell(1,kc)
	    if(kn.eq.k1 .or. kn.eq.k2 .or. kn.eq.k3) kn = pcell(2,kc)
	    if(kn.eq.k1 .or. kn.eq.k2 .or. kn.eq.k3) kn = pcell(3,kc)
	    if(kn.eq.k1 .or. kn.eq.k2 .or. kn.eq.k3) kn = pcell(4,kc)
	    k1 = k2
	    k2 = k3
	  endif
	  kw = wnode(k2)
 54	  n = 5
	  if(wface(1,kw) .ne. k2) n = 6
	  ks = kw
	  kw = wface(n,kw)
	  if(wface(7-n,ks) .eq. kn) then
	    if(wface(3,ks) .eq. kc) then
	      kc = wface(4,ks)
	    else
	      kc = wface(3,ks)
	    endif
	    if(kc .lt. 0) go to 55
	    if(wcell(kc) .ne. 0) go to 55
	    ktri = ktri + 1
	    wcell(kc) = ktri
	    icell = min0(icell+1,100)
	    scell(icell) = kc
            if(ktri .gt. mtri) then
              write(*,*) ' ERROR - MTRI too small for fill'
              stop
            endif
            itri(ktri) = kc
            ptri(ktri) = kn
	    ctri(ktri) = 0
	    k1 = k2
	    k2 = kn
	    go to 53
	  endif
	  if(kw .ne. 0) go to 54
 55       ktri = ktri + 1
          if(ktri .gt. mtri) then
            write(*,*) ' ERROR - MTRI too small for fill'
            stop
          endif
          itri(ktri) = 0
          ptri(ktri) = kn
          ctri(ktri) = 0
	  go to 51
c
 56	  kptri = kptri + 1
          if(kptri .gt. mptri) then
            write(*,*) ' ERROR - MPTRI too small for fill'
            stop
          endif
	  pptri(kptri) = ktri
 5      continue
c
c	start at triangular untouched cells
	do 6 ke = 1, kcell
	  if(wcell(ke) .ne. 0) go to 6
	  if(pcell(1,ke) .ne. pcell(4,ke)) go to 6
	  ksav = ktri
	  itry = 0
 61	  itry = itry + 1
	  if(itry .eq. 5) go to 68
	  if(itry .le. 3) then
	    indx1 = indx(1,itry)
	    indx2 = indx(2,itry)
	    indx3 = indx(3,itry)
	    if(itry .ne. 1) then
	      if(icell .ge. 100) go to 68
	      kmax(itry-1) = ktri
	      do 62 i = 1, icell
	        kc = scell(i)
	        wcell(kc) = 0
 62	      continue
	      ktri = ksav
	    endif
	    icell = 0
	  endif
	  if(itry .eq. 4) then
	    if(icell .ge. 100) go to 68
	    if(ktri .ge. kmax(1) .and. ktri .ge. kmax(2)) go to 68
	    it = 1
	    if(kmax(1) .le. kmax(2)) it = 2
	    indx1 = indx(1,it)
	    indx2 = indx(2,it)
	    indx3 = indx(3,it)
	    do 63 i = 1, icell
	      kc = scell(i)
	      wcell(kc) = 0
 63	    continue
	    ktri = ksav
	    icell = 0
	  endif
c	      
	  rev = .false.
	  k1 = pcell(indx1,ke)
	  k2 = pcell(indx2,ke)
	  kc = ke
	  ktri = ktri + 1
	  kstri = ktri
	  wcell(ke) = ktri + 1
	  icell = min0(icell+1,100)
	  scell(icell) = ke
          if(ktri .gt. mtri) then
            write(*,*) ' ERROR - MTRI too small for fill'
            stop
          endif
	  itri(ktri) = 0
	  ctri(ktri) = 0
	  ptri(ktri) = k1
          ktri = ktri + 1
          if(ktri .gt. mtri) then
            write(*,*) ' ERROR - MTRI too small for fill'
            stop
          endif
          itri(ktri) = ke
          ptri(ktri) = k2
	  ctri(ktri) = 0
 64	  if(pcell(1,kc) .eq. pcell(4,kc)) then
	    kn = pcell(1,kc)
	    if(kn .eq. k2 .or. kn .eq. k1) kn = pcell(2,kc)
	    if(kn .eq. k2 .or. kn .eq. k1) kn = pcell(3,kc)
	  else
	    kk = 1
	    if(pcell(4,kc) .eq. k2) kk = 2
	    if(pcell(1,kc) .eq. k2) kk = 3
	    if(pcell(2,kc) .eq. k2) kk = 4
	    k3 = pcell(kk,kc)
	    ktri = ktri + 1
	    if(ktri .gt. mtri) then
	      write(*,*) ' ERROR - MTRI too small for fill'
	      stop
	     endif
	    itri(ktri) = kc
	    ptri(ktri) = k3
	    ctri(ktri) = 0
	    kn = pcell(1,kc)
	    if(kn.eq.k1 .or. kn.eq.k2 .or. kn.eq.k3) kn = pcell(2,kc)
	    if(kn.eq.k1 .or. kn.eq.k2 .or. kn.eq.k3) kn = pcell(3,kc)
	    if(kn.eq.k1 .or. kn.eq.k2 .or. kn.eq.k3) kn = pcell(4,kc)
	    k1 = k2
	    k2 = k3
	  endif
	  kw = wnode(k2)
 65	  n = 5
	  if(wface(1,kw) .ne. k2) n = 6
	  ks = kw
	  kw = wface(n,kw)
	  if(wface(7-n,ks) .eq. kn) then
	    if(wface(3,ks) .eq. kc) then
	      kc = wface(4,ks)
	    else
	      kc = wface(3,ks)
	    endif
	    if(kc .lt. 0) go to 66
	    if(wcell(kc) .ne. 0) go to 66
	    ktri = ktri + 1
	    wcell(kc) = ktri
	    icell = min0(icell+1,100)
	    scell(icell) = kc
            if(ktri .gt. mtri) then
              write(*,*) ' ERROR - MTRI too small for fill'
              stop
            endif
            itri(ktri) = kc
            ptri(ktri) = kn
	    ctri(ktri) = 0
	    k1 = k2
	    k2 = kn
	    go to 64
	  endif
	  if(kw .ne. 0) go to 65
 66       ktri = ktri + 1
          if(ktri .gt. mtri) then
            write(*,*) ' ERROR - MTRI too small for fill'
            stop
          endif
          itri(ktri) = 0
          ptri(ktri) = kn
          ctri(ktri) = 0
	  if(.not. rev) then
	    rev = .true.
            kh = (ktri - kstri)/2
c $ dir no_recurrence
            do 67 k = kstri, kstri+kh
              kx = ktri + kstri - k
              isave = ptri(k)
              ptri(k) = ptri(kx)
              ptri(kx) = isave
              isave = itri(k)
              itri(k) = itri(kx)
              itri(kx) = isave
	      if(isave .ne. 0) wcell(isave) = kx
	      isave = itri(k)
	      if(isave .ne. 0) wcell(isave) = k
 67         continue
	    ktri = ktri - 1
	    kc = ke
	    kn = pcell(indx1,kc)
	    k2 = pcell(indx2,kc)
	    k1 = pcell(indx3,kc)
	    kw = wnode(k2)
            go to 65
	  endif
	  go to 61
c
 68	  kptri = kptri + 1
          if(kptri .gt. mptri) then
            write(*,*) ' ERROR - MPTRI too small for fill'
            stop
          endif
	  pptri(kptri) = ktri
 6      continue
c
c	start at quad untouched cells
	do 7 ke = 1, kcell
	  if(wcell(ke) .ne. 0) go to 7
	  ksav = ktri
	  itry = 0
 71	  itry = itry + 1
	  if(itry .eq. 4) go to 78
	  if(itry .le. 2) then
	    indx1 = itry
	    indx2 = itry + 1
	    if(itry .ne. 1) then
	      if(icell .ge. 100) go to 78
	      kmaxx = ktri
	      do 72 i = 1, icell
	        kc = scell(i)
	        wcell(kc) = 0
 72	      continue
	      ktri = ksav
	    endif
	    icell = 0
	  endif
	  if(itry .eq. 3) then
	    if(icell .ge. 100) go to 78
	    if(ktri .ge. kmaxx) go to 78
	    indx1 = 1
	    indx2 = 2
	    do 73 i = 1, icell
	      kc = scell(i)
	      wcell(kc) = 0
 73	    continue
	    ktri = ksav
	    icell = 0
	  endif
c	      
	  rev = .false.
	  k1 = pcell(indx1,ke)
	  k2 = pcell(indx2,ke)
	  kc = ke
	  ktri = ktri + 1
	  kstri = ktri
	  wcell(ke) = ktri + 1
	  icell = min0(icell+1,100)
	  scell(icell) = ke
          if(ktri .gt. mtri) then
            write(*,*) ' ERROR - MTRI too small for fill'
            stop
          endif
	  itri(ktri) = 0
	  ctri(ktri) = 0
	  ptri(ktri) = k1
          ktri = ktri + 1
          if(ktri .gt. mtri) then
            write(*,*) ' ERROR - MTRI too small for fill'
            stop
          endif
          itri(ktri) = ke
          ptri(ktri) = k2
	  ctri(ktri) = 0
 74	  if(pcell(1,kc) .eq. pcell(4,kc)) then
	    kn = pcell(1,kc)
	    if(kn .eq. k2 .or. kn .eq. k1) kn = pcell(2,kc)
	    if(kn .eq. k2 .or. kn .eq. k1) kn = pcell(3,kc)
	  else
	    kk = 1
	    if(pcell(4,kc) .eq. k2) kk = 2
	    if(pcell(1,kc) .eq. k2) kk = 3
	    if(pcell(2,kc) .eq. k2) kk = 4
	    k3 = pcell(kk,kc)
	    ktri = ktri + 1
	    if(ktri .gt. mtri) then
	      write(*,*) ' ERROR - MTRI too small for fill'
	      stop
	     endif
	    itri(ktri) = kc
	    ptri(ktri) = k3
	    ctri(ktri) = 0
	    kn = pcell(1,kc)
	    if(kn.eq.k1 .or. kn.eq.k2 .or. kn.eq.k3) kn = pcell(2,kc)
	    if(kn.eq.k1 .or. kn.eq.k2 .or. kn.eq.k3) kn = pcell(3,kc)
	    if(kn.eq.k1 .or. kn.eq.k2 .or. kn.eq.k3) kn = pcell(4,kc)
	    k1 = k2
	    k2 = k3
	  endif
	  kw = wnode(k2)
 75	  n = 5
	  if(wface(1,kw) .ne. k2) n = 6
	  ks = kw
	  kw = wface(n,kw)
	  if(wface(7-n,ks) .eq. kn) then
	    if(wface(3,ks) .eq. kc) then
	      kc = wface(4,ks)
	    else
	      kc = wface(3,ks)
	    endif
	    if(kc .lt. 0) go to 76
	    if(wcell(kc) .ne. 0) go to 76
	    ktri = ktri + 1
	    wcell(kc) = ktri
	    icell = min0(icell+1,100)
	    scell(icell) = kc
            if(ktri .gt. mtri) then
              write(*,*) ' ERROR - MTRI too small for fill'
              stop
            endif
            itri(ktri) = kc
            ptri(ktri) = kn
	    ctri(ktri) = 0
	    k1 = k2
	    k2 = kn
	    go to 74
	  endif
	  if(kw .ne. 0) go to 75
 76       ktri = ktri + 1
          if(ktri .gt. mtri) then
            write(*,*) ' ERROR - MTRI too small for fill'
            stop
          endif
          itri(ktri) = 0
          ptri(ktri) = kn
          ctri(ktri) = 0
	  if(.not. rev) then
	    rev = .true.
            kh = (ktri - kstri)/2
c $ dir no_recurrence
            do 77 k = kstri, kstri+kh
              kx = ktri + kstri - k
              isave = ptri(k)
              ptri(k) = ptri(kx)
              ptri(kx) = isave
              isave = itri(k)
              itri(k) = itri(kx)
              itri(kx) = isave
	      if(isave .ne. 0) wcell(isave) = kx
	      isave = itri(k)
	      if(isave .ne. 0) wcell(isave) = k
 77         continue
	    ktri = ktri - 1
	    kc = ke
	    kn = ptri(ktri+1)
	    k2 = ptri(ktri)
	    k1 = ptri(ktri-1)
	    kw = wnode(k2)
            go to 75
	  endif
	  go to 71
c
 78	  kptri = kptri + 1
          if(kptri .gt. mptri) then
            write(*,*) ' ERROR - MPTRI too small for fill'
            stop
          endif
	  pptri(kptri) = ktri
 7      continue
c
c	patch up ctri and cedge
	do 8 kpp = 1, kptri
	  kpl = 1
	  if(kpp .ne. 1) kpl = pptri(kpp-1) + 1
	  do 81 kp = kpl, pptri(kpp)
	    if(ctri(kp) .ne. 0) go to 81
	    if(kp .eq. kpl) then
	      k1 = ptri(kp)
	      k2 = ptri(kp+1)
	      kc = itri(kp+1)
	    else if(kp .eq. pptri(kpp)) then
	      k1 = ptri(kp-1)
	      k2 = ptri(kp)
	      kc = itri(kp-1)
	    else
	      k1 = ptri(kp-1)
	      k2 = ptri(kp+1)
	      kc = itri(kp)
	    endif
	    kw = wnode(k2)
 811	    n = 5
	    if(wface(1,kw) .ne. k2) n = 6
	    ks = kw
	    kw = wface(n,kw)
	    if(wface(7-n,ks) .eq. k1) then
	      if(wface(3,ks) .eq. kc) then
	        kc = wface(4,ks)
	      else
	        kc = wface(3,ks)
	      endif
	      if(kc .lt. 0) then
	        ctri(kp) = kc
	        cedge(-kc) = kp
	      else
	        ks = wcell(kc) - 1
	        if(pcell(1,kc) .ne. pcell(4,kc)) then
	          if(ptri(ks).ne.k1 .and. ptri(ks).ne.k2) ks = ks + 1
	          kn = ptri(ks)
	          if(kn .eq. k2 .or. kn .eq. k1) kn = ptri(ks+1)
	          if(kn .eq. k2 .or. kn .eq. k1) kn = ptri(ks+2)
	        else
	          kn = pcell(1,kc)
	          if(kn .eq. k2 .or. kn .eq. k1) kn = pcell(2,kc)
	          if(kn .eq. k2 .or. kn .eq. k1) kn = pcell(3,kc)
	        endif
	        if((itri(ks).eq.0 .or. itri(ks+2).eq.0) .and.
     &             ptri(ks+1) .ne. kn) then
	          if(ptri(ks) .eq. kn) ks = ks + 2
	        else
	          if(ptri(ks) .ne. kn) ks = ks + 1
	          if(ptri(ks) .ne. kn) ks = ks + 1
	        endif
	        ctri(kp) = ks
	        ctri(ks) = kp
	      endif
	      go to 81
	    endif
	    if(kw .ne. 0) go to 811
	    write(*,*) 'Error Connection Not Found!',k1,k2,kc
	    stop
 81	  continue
 8	continue
c
c	make face structures
c
c	first from strips
	do 82 kpt = 1, kptri
	  kpl = 1
	  if(kpt .ne. 1) kpl = pptri(kpt-1) + 1
	  if(pptri(kpt)-kpl .lt. 20) go to 82
	  kface = kface + 1
          if(kface .gt. mface) then
            write(*,*) ' ERROR - MFACE too small for fill'
            stop
          endif
	  pface(kface) = ptri(kpl)
	  do 822 kp = kpl, pptri(kpt)-2, 2
	    k1 = ptri(kp)
	    k2 = ptri(kp+2)
	    kw = wnode(k2)
 821	    n = 5
	    if(wface(1,kw) .ne. k2) n = 6
	    ks = kw
	    kw = wface(n,kw)
	    if(wface(7-n,ks) .eq. k1) then
	      kface = kface + 1
              if(kface .gt. mface) then
                write(*,*) ' ERROR - MFACE too small for fill'
                stop
              endif
	      pface(kface) = k2
	      wface(4,ks) = -kface
	      go to 822
	    endif
	    if(kw .ne. 0) go to 821
	    write(*,*) 'Error - No Face!',k1,k2,kpt,kp
	    stop
 822	  continue
	  kpface = kpface + 1
          if(kpface .gt. mpface) then
            write(*,*) ' ERROR - MPFACE too small for fill'
            stop
          endif
	  ppface(kpface) = kface
c
	  kface = kface + 1
          if(kface .gt. mface) then
            write(*,*) ' ERROR - MFACE too small for fill'
            stop
          endif
	  pface(kface) = ptri(kpl+1)
	  do 824 kp = kpl+1, pptri(kpt)-2, 2
	    k1 = ptri(kp)
	    k2 = ptri(kp+2)
	    kw = wnode(k2)
 823	    n = 5
	    if(wface(1,kw) .ne. k2) n = 6
	    ks = kw
	    kw = wface(n,kw)
	    if(wface(7-n,ks) .eq. k1) then
	      kface = kface + 1
              if(kface .gt. mface) then
                write(*,*) ' ERROR - MFACE too small for fill'
                stop
              endif
	      pface(kface) = k2
	      wface(4,ks) = -kface
	      go to 824
	    endif
	    if(kw .ne. 0) go to 823
	    write(*,*) 'Error - No Face2!',k1,k2,kpt,kp
	    stop
 824	  continue
	  kpface = kpface + 1
          if(kpface .gt. mpface) then
            write(*,*) ' ERROR - MPFACE too small for fill'
            stop
          endif
	  ppface(kpface) = kface
 82	continue
c
c	the rest (cross - strips)
	do 9 kf = 1, kfill
	  if(wface(4,kf) .lt. 0) go to 9
	  rev = .false.
	  kface = kface + 1
	  ksface = kface
          if(kface .gt. mface) then
            write(*,*) ' ERROR - MFACE too small for fill'
            stop
          endif
	  pface(kface) = wface(1,kf)
	  k2 = wface(2,kf)
	  ks = kf
 91	  kface = kface + 1
          if(kface .gt. mface) then
            write(*,*) ' ERROR - MFACE too small for fill'
            stop
          endif
	  pface(kface) = k2
	  wface(4,ks) = -kface
	  kw = wnode(k2)
	  kl = 0
 92	  ks = kw
          n = 5
          if(wface(1,ks) .ne. k2) n = 6
	  kw = wface(n,ks)
	  if(wface(4,ks) .gt. 0) then
	    kl = wface(7-n,ks)
	    ksl = ks
	  endif
	  if(kw .ne. 0) go to 92
	  if(kl .ne. 0) then
	    k2 = kl
	    ks = ksl
	    go to 91
	  endif
	  if(.not. rev) then
	    rev = .true.
	    k2 = wface(1,kf)
	    kw = wnode(k2)
	    kh = (kface - ksface)/2
c $ dir no_recurrence
	    do 93 k = ksface, ksface+kh
	      ke = kface + ksface - k
	      isave = pface(k)
	      pface(k) = pface(ke)
	      pface(ke) = isave
 93	    continue
	    kl = 0
	    go to 92
	  endif
	  kpface = kpface + 1
          if(kpface .gt. mpface) then
            write(*,*) ' ERROR - MPFACE too small for fill'
            stop
          endif
	  ppface(kpface) = kface
 9	continue
	return
	end

	subroutine setfac(k1, k2, wnode, wface, mface)
c
c	sets the face into wface structure
c
	common /blkvis/ kc, kfill
	integer wnode(1), wface(6,mface)
c
	if(wnode(k1) .eq. 0) then
c
c	  first time node has been hit
	  kfill = kfill + 1
	  if(kfill .gt. mface) then
	    write(*,*) ' ERROR - MFACE too small for work storage'
	    stop
	  endif
	  wnode(k1) = kfill
	  wface(1,kfill) = k1
	  wface(2,kfill) = k2
	  wface(3,kfill) = kc
	  wface(4,kfill) = 0
	  wface(5,kfill) = 0
	  wface(6,kfill) = 0
	  if(wnode(k2) .eq. 0) then
	    wnode(k2) = kfill
	  else
	    kw = wnode(k2)
 1	    n = 5
	    if(wface(1,kw) .ne. k2) n = 6
	    kf = kw
	    kw = wface(n,kw)
	    if(kw .ne. 0) go to 1
	    wface(n,kf) = kfill
	  endif
	else
c
c	  old node
          kw = wnode(k1)
 2	  n = 1
          if(wface(1,kw) .ne. k1) n = 2
	  if(wface(3-n,kw) .eq. k2) then
	    if(wface(4,kw) .ne. 0) then
	      write(*,*) ' ERROR - More than two cells share a face',
     &                   wface(3,kw),wface(4,kw),kc
	      stop
	    endif
	    wface(4,kw) = kc
	  else
	    n = n + 4
            kf = kw
            kw = wface(n,kw)
            if(kw .ne. 0) go to 2
	    kfill = kfill + 1
            if(kfill .gt. mface) then
              write(*,*) ' ERROR - MFACE too small for work storage'
              stop
            endif
	    wface(n,kf) = kfill
            wface(1,kfill) = k1
            wface(2,kfill) = k2
            wface(3,kfill) = kc
            wface(4,kfill) = 0
            wface(5,kfill) = 0
            wface(6,kfill) = 0
            if(wnode(k2) .eq. 0) then
              wnode(k2) = kfill
            else
              kw = wnode(k2)
 3	      n = 5
              if(wface(1,kw) .ne. k2) n = 6
              kf = kw
              kw = wface(n,kw)
              if(kw .ne. 0) go to 3
              wface(n,kf) = kfill
            endif
	  endif
	endif
c
	return
	end
C=======================================================================
        SUBROUTINE v2call(IOPT,CMNCOL,CMUNIT,
     &	                   XYPIX,XYMIN,XYMAX,
     &  		   NKEYS,IKEYS,FKEYS,FLIMS,
     &			   MNODE,MPTRI,MPPTRI,
     &			   MFACE,MPFACE,MEDGE,MPEDGE,BGCOLOUR)
C       ****************************************************************
C       *                                                              *
C       *  SUBROUTINE TO CALL V2_INIT                                  *
C       *  THIS IS BETTER DONE WITH A FORTRAN ROUTINE, SINCE IT IS     *
C       *  A PAIN IN THE ASS TO PASS CHAR VARIABLES FROM C TO FORTRAN  *
C       *                                                              *
C       ****************************************************************
        
	INTEGER IOPT,CMNCOL,CMUNIT,BGCOLOUR
	INTEGER XYPIX(2),NKEYS,IKEYS(14),FKEYS(14)
	INTEGER MNODE,MPTRI,MPPTRI,MFACE,MPFACE,MEDGE,MPEDGE
	REAL    XYMIN(2), XYMAX(2),FLIMS(2,14)
        CHARACTER CMFILE*53
	CHARACTER*16 TKEYS(14)
	CHARACTER*80 TITL
	
	DATA TKEYS / 'PRESSURE        ',
     &	             'STREAMFUNCTION  ',
     &		     'VORTICITY       ',
     &		     'FLOW VECTORS    ',
     &		     'VELOCITY Ux     ',
     &		     'VELOCITY Uy     ', 
     &		     'ABS. VEL. |U|   ',
     &		     'STRL./STAT.STRL.',
     &		     'GRID-VEL. VECT. ',
     &		     'STOPPING TIME   ',
     &		     'MOVIE CREATION  ',
     &                 'TUR. KIN. ENERGY',
     &                 'TUR. DISS. RATE',
     &                 'TUR. VISCOSITY' /

       TITL='HOT PICS!!!!'
C UNIX Stuttgart
       IF(BGCOLOUR.EQ.0) THEN
          CMFILE='/bau/stat16/users/statik/ccarat/viscol/spec_black.col'
       ELSE IF (BGCOLOUR.EQ.1) THEN
          CMFILE='/bau/stat16/users/statik/ccarat/viscol/spec_white.col'
       ELSE IF (BGCOLOUR.EQ.2) THEN
          CMFILE='/bau/stat16/users/statik/ccarat/viscol/spec_grey_.col'
C LINUX Stuttgart
       ELSE IF (BGCOLOUR.EQ.10) THEN
          CMFILE='/home/ccarat/viscol/spec_black.col'
       ELSE IF (BGCOLOUR.EQ.11) THEN
          CMFILE='/home/ccarat/viscol/spec_white.col'
       ELSE
          WRITE(*,*) 'cannot find colour file - STOPPING!'
          STOP
       ENDIF
C 
C      WRITE(*,*) 'call of V2_INIT not compiled in!!!'
C      WRITE(*,*) 'change in  src/fortran/vis2_qat2v2.f'
C      STOP
C 
	CALL V2_INIT(TITL,IOPT,CMNCOL,CMFILE,CMUNIT,
     &  	    XYPIX,XYMIN,XYMAX,
     &  	   NKEYS,IKEYS,TKEYS,FKEYS,FLIMS,
     &  	   MNODE,MPTRI,MPPTRI,
     &  	   MFACE,MPFACE,MEDGE,MPEDGE)
	
	RETURN
	END

C=======================================================================
        SUBROUTINE v3call(IOPT,WIN3TMP,CMUNIT,
     &	                  NKEYS,IKEYS,FKEYS,FLIMS,
     &			  MIRROR,KNODE,KCEL1,KCEL2,KCEL3,KCEL4,
     &			  KSURF,BGCOLOUR)
C       ****************************************************************
C       *                                                              *
C       *  SUBROUTINE TO CALL V3_INIT                                  *
C       *  THIS IS BETTER DONE WITH A FORTRAN ROUTINE, SINCE IT IS     *
C       *  A PAIN IN THE ASS TO PASS CHAR VARIABLES FROM C TO FORTRAN  *
C       *                                                              *
C       ****************************************************************
        
	INTEGER IOPT,CMUNIT,BGCOLOUR,WIN3TMP
	INTEGER NKEYS,IKEYS(7),FKEYS(7)
	INTEGER MIRROR,KNODE,KCEL1,KCEL2,KCEL3,KCEL4
	INTEGER KEQUIV,KNPTET,KPTET,KNBLOCK,BLOCKS,KSURF
	REAL   FLIMS(2,7)
        LOGICAL WIN3D
	CHARACTER CMFILE*53
	CHARACTER*16 TKEYS(7)
	CHARACTER*80 TITL
	
	DATA TKEYS / 'VELOCITY Ux     ',
     &		     'VELOCITY Uy     ',
     &		     'VELOCITY Uz     ', 
     &               'PRESSURE        ',     
     &		     'FLOW VECTORS    ',
     &		     'ABS. VEL. |U|   ',
     &		     'MOVIE CREATION  '/

       TITL='HOT PICS!!!!'
C UNIX Stuttgart
       IF(BGCOLOUR.EQ.0) THEN
          CMFILE='/bau/stat16/users/statik/ccarat/viscol/spec_black.col'
       ELSE IF (BGCOLOUR.EQ.1) THEN
          CMFILE='/bau/stat16/users/statik/ccarat/viscol/spec_white.col'
       ELSE IF (BGCOLOUR.EQ.2) THEN
          CMFILE='/bau/stat16/users/statik/ccarat/viscol/spec_grey_.col'
C LINUX Stuttgart
       ELSE IF (BGCOLOUR.EQ.10) THEN
          CMFILE='/home/ccarat/viscol/spec_black.col'
       ELSE IF (BGCOLOUR.EQ.11) THEN
          CMFILE='/home/ccarat/viscol/spec_white.col'
       ELSE
          WRITE(*,*) 'cannot find colour file - STOPPING!'
          STOP
       ENDIF
       WIN3D=.TRUE.
       IF (WIN3TMP.EQ.1) THEN
          WIN3D=.FALSE.
       ENDIF
       KEQUIV=0
       KNPTET=0
       KPTET=0
       KNBLOCK=0
       BLOCKS=0
       KNSURF=0
C 
C      WRITE(*,*) 'call of V3_INIT not compiled in!!!'
C      WRITE(*,*) 'change in  src/fortran/vis2_qat2v2.f'
C      STOP
C 
	CALL V3_INIT(TITL,IOPT,CMFILE,CMUNIT,WIN3D,
     &  	     NKEYS,IKEYS,TKEYS,FKEYS,FLIMS,MIRROR,
     &  	     KNODE,KEQUIV,KCEL1,KCEL2,KCEL3,KCEL4,
     &  	     KNPTET,KPTET,KNBLOCK,BLOCKS,KSURF,KNSURF)
	
	RETURN
	END


C=======================================================================
        SUBROUTINE v3call_struct(IOPT,WIN3TMP,CMUNIT,
     &	                  NKEYS,IKEYS,FKEYS,FLIMS,
     &			  MIRROR,KNODE,KCEL1,KCEL2,KCEL3,KCEL4,
     &			  KSURF,BGCOLOUR)
C       ****************************************************************
C       *                                                              *
C       *  SUBROUTINE TO CALL V3_INIT                                  *
C       *  THIS IS BETTER DONE WITH A FORTRAN ROUTINE, SINCE IT IS     *
C       *  A PAIN IN THE ASS TO PASS CHAR VARIABLES FROM C TO FORTRAN  *
C       *                                                              *
C       ****************************************************************
        
	INTEGER IOPT,CMUNIT,BGCOLOUR,WIN3TMP
	INTEGER NKEYS,IKEYS(5),FKEYS(5)
	INTEGER MIRROR,KNODE,KCEL1,KCEL2,KCEL3,KCEL4
	INTEGER KEQUIV,KNPTET,KPTET,KNBLOCK,BLOCKS,KSURF
	REAL   FLIMS(2,5)
        LOGICAL WIN3D
	CHARACTER CMFILE*53
	CHARACTER*16 TKEYS(5)
	CHARACTER*80 TITL
	
	DATA TKEYS / 'DISPL. Dx       ',
     &		     'DISPL. Dy       ',
     &		     'DISPL. Dz       ', 
     &		     'ABS. DISP. |D|  ',
     &		     'MOVIE CREATION  '/

       TITL='HOT PICS!!!!'
C UNIX Stuttgart
       IF(BGCOLOUR.EQ.0) THEN
          CMFILE='/bau/stat16/users/statik/ccarat/viscol/spec_black.col'
       ELSE IF (BGCOLOUR.EQ.1) THEN
          CMFILE='/bau/stat16/users/statik/ccarat/viscol/spec_white.col'
       ELSE IF (BGCOLOUR.EQ.2) THEN
          CMFILE='/bau/stat16/users/statik/ccarat/viscol/spec_grey_.col'
C LINUX Stuttgart
       ELSE IF (BGCOLOUR.EQ.10) THEN
          CMFILE='/home/ccarat/viscol/spec_black.col'
       ELSE IF (BGCOLOUR.EQ.11) THEN
          CMFILE='/home/ccarat/viscol/spec_white.col'
       ELSE
          WRITE(*,*) 'cannot find colour file - STOPPING!'
          STOP
       ENDIF
       WIN3D=.TRUE.
       IF (WIN3TMP.EQ.1) THEN
          WIN3D=.FALSE.
       ENDIF
       KEQUIV=0
       KNPTET=0
       KPTET=0
       KNBLOCK=0
       BLOCKS=0
       KNSURF=0
C 
C      WRITE(*,*) 'call of V3_INIT not compiled in!!!'
C      WRITE(*,*) 'change in  src/fortran/vis2_qat2v2.f'
C      STOP
C 
	CALL V3_INIT(TITL,IOPT,CMFILE,CMUNIT,WIN3D,
     &  	     NKEYS,IKEYS,TKEYS,FKEYS,FLIMS,MIRROR,
     &  	     KNODE,KEQUIV,KCEL1,KCEL2,KCEL3,KCEL4,
     &  	     KNPTET,KPTET,KNBLOCK,BLOCKS,KSURF,KNSURF)
	
	RETURN
	END
