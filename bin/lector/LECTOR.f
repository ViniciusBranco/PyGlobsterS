	PROGRAM LECTOR
C Version 01-11-2008
C------------------------------------------------------------------
C Copyright: Alexandre Vazdekis
C------------------------------------------------------------------
C Instituto de Astrofisica de Canarias
C La Laguna 38200, Tenerife, Spain
C E-mail: vazdekis@ll.iac.es
C
C------------------------------------------------------------------
C This is free software; you can redistribute or modify it as you wish.
C------------------------------------------------------------------
C
C DESCRIPTION:
C 
C This program computes line-strength indices in 1D ascii spectra (logarithmic
C wavelength scale is also accepted) . In particular it measures the Lick
C indices (Worthey et al.1994, ApJS,94,687; Worthey & Ottaviani 1997,
C ApJS,111,377), Hgamma age indicators of Vazdekis & Arimoto 1999 (ApJ,525,144)
C and Vazdekis et al. 2001 (ApJ, 549,274), the indices of Rose 1994 (AJ,107,206)
C and Jones & Worthey 1995 (ApJ,446,L31) and the CaII triplet indices
C of Cenarro et al 2001 (MNRAS,326,959) and the Hbeta_o index of
C Cervantes & Vazdekis 2008 (MNRAS, in press, arXiv:0810.3240)
C
C------------------------------------------------------------------
C INPUT FILES:
C
C "BANDS"
C
C which include the definitions of the Lick indices or your own Lick style index
C definition. These indices are composed of two pseudocontinua at each side of
C the feature and the feature itself. Please have a look at the file BANDS. In
C column 7  you should type "A" when the index is expressed in Angstrom and "M"
C for magnitude.
C
C A file (named as you wish) including a column with the names of the spectra.
C Optionally, you can include the corresponding redshift to each spectrum in a
C second column (it can be given in km/s or z). If a second column with the
C redshift values is not provided the program will ask you for the value to be
C redshifted. Alternatively, you can run the program for a single spectrum.
C
C OUTPUT FILES:
C
C There are three sets of output files. The first two sets are meant to be used
C meant for plotting purposes:
C
C 1.- Index measurements can be found in:
C
C "GIVENLIST_LINE" 
C
C which include the Hgamma age indicators of 
C   Vazdekis et al. 2001 and Vazdekis & Arimoto 1999
C + Lick indices of Worthey etal 1994 and Worthey & Ottaviani 1997
C + Cenarro et al 2001 CaII triplet indices
C GIVENLIST_ROSE (Rose 1994 + Jones & Worthey 1995 indices)
C
C 2.- Index error estimates can be found in:
C
C "GIVENLIST_LINE_ERR"
C "GIVENLIST_ROSE_ERR"
C
C 3.- Finally, it is provided an easy reading version of these 
C tables:
C
C "GIVENLIST_INDICES"
C
C If any index falls outside the spectral range its value is
C 99.999. If an index falls outside the spectral range it is
C ignored.
C

C------------------------------------------------------------------
C COMPILATION:
C
C f77 -o LECTOR LECTOR.f
C or
C g77 -o LECTOR LECTOR.f 
C 
C------------------------------------------------------------------
C
C EXECUTING LECTOR:
C
C You just need to type:
C
C LECTOR
C
C Do not forget to include the file "BANDS" in the same 
C directory as the LECTOR program.
C
C Then the program asks you whether you have a single spectrum or a list
C of spectra. If it is a list, the program will ask  you whether a second
C column with the redshift values is given in the list. Then the program
C asks whether the redshift is given in km/s or z. Next asks whether you
C wish to calculate the index errors. If you agree then it asks for the
C readout noise and electrons/ADU factor.
C
C------------------------------------------------------------------
C UPDATES:
C 
C 1-November-2008: Line-strengths for atomic indices from 
C redshifted spectra are now correct. Thanks goes to Patricia
C Sanchez-Blazquez & Alexander Hansson.
C  
C 1-June-2003: a more friendly version for the use of this 
C program. The program can also measure indices in spectra 
C with logarithmic wavelength scale.
C
C 1-July-2002: added "indexf" subroutine written by
C N. Cardiel and J. Gorgas for measuring the CaII triplet
C and Paschen features defined in Cenarro et al. 2001
C (MNRAS,326,959).
C
C 25-July-2001: error estimate corrections on the basis
C of Jim Rose's suggestions to include the readout noise
C of the total number of pixels in the selected aperture.
C------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION b(50,6),bnom(100),mag(50),bsin(50,6)
      DIMENSION rose(20),erose(20),rosen(20),ncat(5),cat(5),ecat(5)
      DIMENSION xb(99000,2),xi(99000,2),xr(99000,2),ewidx(50)
      DIMENSION eewidx(50),tonte(99000,2),s11(99000),s22(99000)
      CHARACTER*80 SPC,SPC0,tso,ewtso,lick,jrose,licker,roser,expla,chau
      character*80 SPCt
      CHARACTER*50 bnom
      CHARACTER*30 rosen,ncat
      LOGICAL sl
      character*1 cv(80),ansmin,ansmax
      character*1 sll,s2,mag,EWY,eory,zkms
      COMMON/cmnmx/s(99000,2),nnr
      COMMON/croz/s1(99000),sf(99000),reds
      COMMON/eerr/RN,EADU
C
 	nirose=15
 	rosen(1)='Hdelta/4045                   '
 	rosen(2)='Hdelta/4063                   '
 	rosen(3)='SrII/4045                     '
 	rosen(4)='SrII/4063                     '
 	rosen(5)='Hgamma/Gband                  '
 	rosen(6)='Hgamma/4325                   '
 	rosen(7)='4289/4271                     '
 	rosen(8)='4384/4352                     '
 	rosen(9)='p[Fe/H]                       '
 	rosen(10)='CaII                          '
 	rosen(11)='3888/3859                     '
 	rosen(12)='p4220/4208                    '
 	rosen(13)='FeI_HR                        '
 	rosen(14)='CaI_HR                        '
 	rosen(15)='HgammaHR                      '
C
 	nicat=3
 	ncat(1)='CaT*                          '
 	ncat(2)='PaT                           '
 	ncat(3)='CaT                           '
C
 	vluz=2.997925e18
 	Cst=-2.5d0*dlog10(2.718281828d0)
 	write(*,*)'Single spectrum or list of spectra(s/l)?'
 	read(*,*)sll
        if ((sll(1:1).eq.'l').or.(sll(1:1).eq.'L')) then
           sl=.true.
        else
           sl=.false.
        end if
 	IF(sl)THEN
 	  write(*,*)'Name of the list of spectra?'
 	  read(*,'(A)')tso
 	  write(*,*)'Is any redshift written in column 2 of',tso
 	  write(*,*)'(y/n)?'
 	  read(*,'(A1)')s2
	  if(s2.eq.'y'.or.s2.eq.'Y') then
	    s2='y'
	  else
	    s2='n'
	  endif
 	ELSE
 	  write(*,*)'Name of the spectrum?'
 	  read(*,'(A)')tso
 	  s2='n'
 	ENDIF
 	if(s2.ne.'y')then
 	  s2='n'
 	  write(*,*)'Redshift to be applied ?'
 	  read(*,*)vrcg
 	else
 	  s2='y'
 	endif
 	if(vrcg.ne.0.0d0.or.s2.eq.'y')then
 	 write(*,*)'Is the redshift given in z or velocity (z/v)?'
 	 read(*,*)zkms
 	endif
C
C Formatted ascii file with index names & values
	EWY='y'
C
 	write(*,*)'Estimate errors (y/n)?'
 	read(*,'(A1)')eory
 	if(eory.eq.'y'.or.eory.eq.'Y')then
 	 eory='y'
 	 write(*,*)'Please provide the following Readout Noise:'
 	 write(*,*)' RN*SQRT(total number of pixels) '
 	 write(*,*)'in the aperture taking into account the slit width'
 	 write(*,*)'Readout Noise in e?'
 	 read(*,*)RN
 	 write(*,*)'e/ADU?'
 	 read(*,*)EADU
 	else
 	 eory='n'
	 EADU=1.0d0
	 RN=1.0d0
 	endif
C
 	ewtso=tso
 	write(ewtso(lnblnk(tso)+1:lnblnk(tso)+1+8),'(A8)')'_INDICES'
 	lick=tso
 	write(lick(lnblnk(tso)+1:lnblnk(tso)+1+5),'(A5)')'_LINE'
 	jrose=tso
 	write(jrose(lnblnk(tso)+1:lnblnk(tso)+1+5),'(A5)')'_ROSE'
 	licker=tso
 	write(licker(lnblnk(tso)+1:lnblnk(tso)+1+9),'(A9)')'_LINE_ERR'
 	roser=tso
 	write(roser(lnblnk(tso)+1:lnblnk(tso)+1+9),'(A9)')'_ROSE_ERR'
C
	if(EWY.eq.'y'.or.EWY.eq.'Y')then
	  EWY='y'
	  OPEN(79,FILE=ewtso,STATUS='OLD',ERR=8956)
	  CLOSE(79,STATUS='DELETE')
8956	  OPEN(79,IOSTAT=IOS,FILE=ewtso,STATUS='NEW')
       	else
	  EWY='n'
	endif
	OPEN(97,FILE=lick,STATUS='OLD',ERR=8955)
	CLOSE(97,STATUS='DELETE')
8955	OPEN(97,IOSTAT=IOS,FILE=lick,STATUS='NEW')
	OPEN(95,FILE=jrose,STATUS='OLD',ERR=8954)
	CLOSE(95,STATUS='DELETE')
8954	OPEN(95,IOSTAT=IOS,FILE=jrose,STATUS='NEW')
	if(eory.eq.'y')then
	  OPEN(89,FILE=licker,STATUS='OLD',ERR=8953)
	  CLOSE(89,STATUS='DELETE')
8953	  OPEN(89,IOSTAT=IOS,FILE=licker,STATUS='NEW')
	  OPEN(63,FILE=roser,STATUS='OLD',ERR=8952)
	  CLOSE(63,STATUS='DELETE')
8952	  OPEN(63,IOSTAT=IOS,FILE=roser,STATUS='NEW')
	endif
C
 	open(99,file='BANDS',status='old')
 	do ibas=1,4
 	   read(99,'(A)')expla
 	enddo
 	numbd=0
 	do i=1,999
 	  read(99,*,end=19)(bsin(i,j),j=1,6),mag(i),bnom(i)
 	  numbd=numbd+1
 	enddo
19     close(99)
C
32	format(A30,2x,f8.3,1x,f8.3,1x,f7.2,1x,f7.2)
C
       IF(sl)THEN
 	 open(99,file=tso,status='old')
 	 jlist=99999
       ELSE
 	 jlist=1
       ENDIF
C Bucle general:
 	DO n=1,jlist
C
       IF(sl)THEN
 	 read(99,'(A80)',end=300)SPC0
	 SPCt=SPC0
	 do j=1,len(SPC0)-2
	   if(index(SPCt,' ').gt.j)then
 	     if(s2.eq.'y')then
	       read(SPC0(j:index(SPCt,' ')),'(A)')SPC
	       read(SPC0(index(SPCt,' '):len(SPC0)),*)vrcg
	       goto 233
	     else
	       read(SPC0(j:index(SPCt,' ')),'(A)')SPC
	       goto 233
	     endif
	   else
             write(SPCt(j:j),'(A1)')'c'
	   endif
	 enddo
       ELSE
 	  SPC=tso
       ENDIF
233    vrec=vrcg
C
 	if(EWY.eq.'y')then
      write(79,*)'________________________________________________',SPC
      write(79,*)' '
 	endif
      write(*,*)'________________________________________________',SPC
      write(*,*)' '
      open(90,file=SPC,status='old')
      nnr=0
      do nni=1,90000
 	        read(90,'(A)',end=191,ERR=8976)chau
		read(chau,'(80(A1))')(cv(lnni),lnni=1,80)
		ansmax=cv(1)
		ansmin=cv(1)
		do lnni=2,80
		 if(lgt(cv(lnni),ansmax)) ansmax=cv(lnni)
		 if(llt(cv(lnni),ansmax)) ansmin=cv(lnni)
		enddo
		if(ansmin.eq.ansmax.and.ansmin.eq.' ') goto 8976		 
	        read(chau,*,ERR=8976)xaton1,xaton2
 		nnr=nnr+1
 		s(nnr,1)=xaton1
 		s(nnr,2)=xaton2
8976		continue
      enddo
191   close(90)
c        
      disps=abs(s(int(1+nnr*0.5),1)-s(int(nnr*0.5),1))
      IF(disps.lt.0.001)THEN
       do l=1,nnr
 	  s(l,1)=10**s(l,1)
 	  s11(l)=s(l,1)
	  s22(l)=s(l,2)
       enddo
       disps=dble(abs(s(nnr,1)-s(1,1))/(nnr-1))
       dlamb=s(1,1)-disps
c       call LOCT(s1,nnr,slog(l),k1,k2)
       do l=1,nnr
       	  dlamb=dlamb+disps
	  s1(l)=dlamb
c 	call LOCT(s11,nnr,s1(l),k1,k2)
 	  call hunt(s11,nnr,s1(l),k1)
	  if(l.eq.1)then
	   call polint(s11(k1),s22(k1),2,s1(l),sf(l),dy)
	  elseif(l.eq.2)then
	   call polint(s11(k1),s22(k1),3,s1(l),sf(l),dy)
	  elseif(l.eq.nnr-1)then
	   call polint(s11(k1-1),s22(k1-1),3,s1(l),sf(l),dy)
	  elseif(l.eq.nnr)then
	   call polint(s11(k1-1),s22(k1-1),2,s1(l),sf(l),dy)
	  else
	   call polint(s11(k1-1),s22(k1-1),4,s1(l),sf(l),dy)
	  endif
       enddo
       do ll=1,nnr
	 s(ll,1)=s1(ll)
	 s(ll,2)=sf(ll)
       enddo
      ELSE
        do ll=1,nnr
 	  s1(ll)=s(ll,1)
	  sf(ll)=s(ll,2)
        enddo        
      ENDIF
C
      vcc=0.0d0
      if(zkms.eq.'z'.or.zkms.eq.'Z')then
 	   reds=1.0d0+vrec
      else
           vcc=vrec*1.0e13/vluz
 	   reds=(1.0d0+vcc)/sqrt(1.0d0-vcc*vcc)
      endif
C------------Indexf-----------
      do irse=1,nicat
 	  cat(irse)=99.999
 	  ecat(irse)=99.999
      enddo
      IF((8461.*reds).gt.s(1,1).and.(8792.*reds).lt.s(nnr,1))THEN
 	call NICO(cat(1),cat(2),cat(3))
      ENDIF
C---------Rose indices--------
       do irs=1,nirose
 	 rose(irs)=99.999
 	 erose(irs)=99.999
       enddo
       IF((4041.*reds).gt.s(1,1).and.(4108.*reds).lt.s(nnr,1))THEN
 	call Hd4045(rose(1),erose(1))
       ENDIF
       IF((4059.*reds).gt.s(1,1).and.(4108.*reds).lt.s(nnr,1))THEN
 	call Hd4063(rose(2),erose(2))
       ENDIF
       IF((4041.*reds).gt.s(1,1).and.(4083.*reds).lt.s(nnr,1))THEN
 	call SrII4045(rose(3),erose(3))
       ENDIF
       IF((4059.*reds).gt.s(1,1).and.(4083.*reds).lt.s(nnr,1))THEN
 	call SrII4063(rose(4),erose(4))
       ENDIF
       IF((4303.*reds).gt.s(1,1).and.(4347.*reds).lt.s(nnr,1))THEN
 	call HgG(rose(5),erose(5))
       ENDIF
       IF((4321.*reds).gt.s(1,1).and.(4347.*reds).lt.s(nnr,1))THEN
 	call Hg4325(rose(6),erose(6))
       ENDIF
       IF((4267.*reds).gt.s(1,1).and.(4295.*reds).lt.s(nnr,1))THEN
 	call s42894271(rose(7),erose(7))
       ENDIF
       IF((4347.*reds).gt.s(1,1).and.(4389.*reds).lt.s(nnr,1))THEN
 	call s43844352(rose(8),erose(8))
       ENDIF
       IF((4033.*reds).gt.s(1,1).and.(4094.*reds).lt.s(nnr,1))THEN
 	call pFeH(rose(9),erose(9))
       ENDIF
       IF((3928.*reds).gt.s(1,1).and.(3973.*reds).lt.s(nnr,1))THEN
 	call CaII(rose(10),erose(10))
       ENDIF
       IF((3856.*reds).gt.s(1,1).and.(3893.*reds).lt.s(nnr,1))THEN
 	call s38883859(rose(11),erose(11))
       ENDIF
       IF((4203.*reds).gt.s(1,1).and.(4225.*reds).lt.s(nnr,1))THEN
 	call p42204208(rose(12),erose(12))
       ENDIF
       IF((4033.*reds).gt.s(1,1).and.(4055.*reds).lt.s(nnr,1))THEN
 	call FeIHR(rose(13))
       ENDIF
       IF((4215.*reds).gt.s(1,1).and.(4236.*reds).lt.s(nnr,1))THEN
 	call CaIHR(rose(14))
       ENDIF
       IF((4328.*reds).gt.s(1,1).and.(4353.*reds).lt.s(nnr,1))THEN
 	call HgHR(rose(15))
       ENDIF
       do ke=1,nirose
 	 if(erose(ke).gt.99..or.eory.ne.'y') erose(ke)=99.999
       enddo
 	 write(95,'(A50,1x,40(f7.3))')SPC,(rose(iii),iii=1,nirose)
       if(eory.eq.'y')then
 	 write(63,'(A50,1x,40(f7.3))')SPC,(erose(iii),iii=1,nirose)
       endif
       IF(EWY.eq.'y')THEN
 	 do ir=1,nirose
      if(rose(ir).lt.99.)write(79,32)rosen(ir),rose(ir),erose(ir)
 	 enddo
 	 write(79,*)' '
       ENDIF
       do ir=1,nirose
      if(rose(ir).lt.99.)write(*,32)rosen(ir),rose(ir),erose(ir)
       enddo
 	 write(*,*)' '
C---------Lick style indices---------
 	 DO i=1,numbd
 	   ewidx(i)=99.999
 	   eewidx(i)=99.999
 	   ntonto=0
 	   do kaku=1,6
 		b(i,kaku)=bsin(i,kaku)*reds
 	   enddo
 	 IF(b(i,1).ge.s(1,1).and.b(i,6).le.s(nnr,1))THEN
 	  ntonto=ntonto+1
 	  lb=0
 	  li=0
 	  lr=0
 	  do j=1,nnr
c 	    if(j.eq.nnr)then
c 	     d=abs(s(j,1)-s(j-1,1))
c 	    else
c 	     d=abs(s(j+1,1)-s(j,1))
c 	    endif
 	    call LOCT(s1,nnr,b(i,1),k1,k2)
	    d1=abs(s(k2,1)-s(k1,1))	    
 	    call LOCT(s1,nnr,b(i,2),k1,k2)
	    d2=abs(s(k2,1)-s(k1,1))	    
 	    IF(s(j,1).ge.b(i,1)-d1.and.s(j,1).le.b(i,2)+d2)THEN
 	      lb=lb+1
 	      xb(lb,1)=s(j,1)
 	      xb(lb,2)=s(j,2)
 	    ENDIF
 	    call LOCT(s1,nnr,b(i,3),k1,k2)
	    d3=abs(s(k2,1)-s(k1,1))	    
 	    call LOCT(s1,nnr,b(i,4),k1,k2)
	    d4=abs(s(k2,1)-s(k1,1))	    
 	    IF(s(j,1).ge.b(i,3)-d3.and.s(j,1).le.b(i,4)+d4)THEN
 	      li=li+1
 	      xi(li,1)=s(j,1)
 	      xi(li,2)=s(j,2)
 	    ENDIF
 	    call LOCT(s1,nnr,b(i,5),k1,k2)
	    d5=abs(s(k2,1)-s(k1,1))	    
 	    call LOCT(s1,nnr,b(i,6),k1,k2)
	    d6=abs(s(k2,1)-s(k1,1))	    
 	    IF(s(j,1).ge.b(i,5)-d5.and.s(j,1).le.b(i,6)+d6)THEN
 	      lr=lr+1
 	      xr(lr,1)=s(j,1)
 	      xr(lr,2)=s(j,2)
 	    ENDIF
 	  enddo
 	  wb=(b(i,1)+b(i,2))*0.5d0
 	  wr=(b(i,5)+b(i,6))*0.5d0
 	  wrb=wr-wb
 	  ew=0.0d0
 	  sb=0.0d0
 	  sr=0.0d0
 	  CALL ITCT(xb,lb,b(i,1),b(i,2),sb)
 	  CALL ITCT(xr,lr,b(i,5),b(i,6),sr)
 	  sb=sb/abs(b(i,2)-b(i,1))
 	  sr=sr/abs(b(i,6)-b(i,5))
C -----------------errors-----------------
 	if(eory.eq.'y')then
C  Cardiel etal 1998 + Vazdekis & Arimoto 1999
 	Uwb=(b(i,1)+b(i,2))*0.5d0
 	Uwr=(b(i,5)+b(i,6))*0.5d0
 	Uwc=(b(i,3)+b(i,4))*0.5d0
 	Urca=(Uwr-Uwc)/(Uwr-Uwb)
 	Urcb=(Uwc-Uwb)/(Uwr-Uwb)
 	Udwc=abs(b(i,4)-b(i,3))
 	Udwb=abs(b(i,2)-b(i,1))
 	Udwr=abs(b(i,6)-b(i,5))
 	c2=sqrt((1./Udwc)+(Urca*Urca/Udwb)+(Urcb*Urcb/Udwr))
 	c1=Udwc*c2
 	SNA=0.
 	   do j=1,li
 		tonte(j,1)=xi(j,1)
 		tonte(j,2)=(xi(j,2)/EADU)+(RN/EADU)*(RN/EADU)
 		tonte(j,2)=xi(j,2)/sqrt(tonte(j,2))
 	   enddo
 	   call ITCT(tonte,li,b(i,3),b(i,4),snerr)
 	   SNA=snerr+SNA
 	   do j=1,lb
 		tonte(j,1)=xb(j,1)
 		tonte(j,2)=(xb(j,2)/EADU)+(RN/EADU)*(RN/EADU)
 		tonte(j,2)=xb(j,2)/sqrt(tonte(j,2))
 	   enddo
 	   call ITCT(tonte,lb,b(i,1),b(i,2),snerr)
 	   SNA=snerr+SNA
 	   do j=1,lr
 		tonte(j,1)=xr(j,1)
 		tonte(j,2)=(xr(j,2)/EADU)+(RN/EADU)**2
 		tonte(j,2)=xr(j,2)/sqrt(tonte(j,2))
 	   enddo
 	   call ITCT(tonte,lr,b(i,5),b(i,6),snerr)
 	   SNA=snerr+SNA
	   d=(d1+d6)*0.5d0
           SNA=SNA/(d*sqrt(d)*((Udwc+Udwb+Udwr)/d))
 	endif
C -------fin error---------------------------
 	  if(mag(i).eq.'A'.or.mag(i).eq.'a')then
 	   do j=1,li
 	    continuo=(1.0d0/wrb)*(sb*(wr-xi(j,1))+sr*(xi(j,1)-wb))
 	    xi(j,2)=(1.0d0-(xi(j,2)/continuo))
 	   enddo
 	  elseif(mag(i).eq.'M'.or.mag(i).eq.'m')then
 	   do j=1,li
 	    continuo=(1.0d0/wrb)*(sb*(wr-xi(j,1))+sr*(xi(j,1)-wb))
 	    xi(j,2)=xi(j,2)/continuo
 	   enddo
 	  else
      write(*,*)'Write A or M in column 7 of BANDS for',bnom(i)
 	    goto 300
 	  endif
 	  CALL ITCT(xi,li,b(i,3),b(i,4),ew)
 	  if(mag(i).eq.'M'.or.mag(i).eq.'m')then
 	    ew=-2.50d0*dlog10(ew/(b(i,4)-b(i,3)))
 	    if(eory.eq.'y')then
 		c3=2.50d0*c2*dlog10(2.718281828d0)
 		eew=c3/SNA
 	    endif
 	  else
	    ew=ew/reds
 	    if(eory.eq.'y')then
 		eew=(c1-c2*ew)/SNA
 	    else
 		eew=99.999
 	    endif
 	  endif
 	  ewidx(i)=ew
 	  eewidx(i)=eew
 	  if(eewidx(i).gt.99.) eewidx(i)=99.999
 	  write(*,32)bnom(i),ewidx(i),eewidx(i)
 	  if(EWY.eq.'y')then
 	   if(ewidx(i).lt.99.) then
 	      write(79,32)bnom(i),ewidx(i),eewidx(i)
 	   endif
 	  endif
 	 ENDIF
 	ENDDO

       IF((8461.*reds).gt.s(1,1).and.(8792.*reds).lt.s(nnr,1))THEN
 	do irse=1,nicat
 	 write(*,32)ncat(irse),cat(irse),ecat(irse)
 	enddo
 	if(EWY.eq.'y')then
 	 do irse=1,nicat
 	  write(79,32)ncat(irse),cat(irse),ecat(irse)
 	 enddo
 	endif
       ENDIF
       write(97,'(A50,1x,44(f7.3))')SPC,(ewidx(iii),iii=1,numbd)
     & ,(cat(iii),iii=1,nicat)
 	if(eory.eq.'y')then
       write(89,'(A50,1x,44(f7.3))')SPC,(eewidx(iii),iii=1,numbd)
     & ,(ecat(iii),iii=1,nicat)
 	endif
 	ENDDO
300	close(99)
 	close(97)
 	close(95)
 	if(eory.eq.'y')then
 		close(89)
 		close(63)
 	endif
 	if(EWY.eq.'y')then
 		close(79)
 	endif
 	write(*,*)' '
 	write(*,*)'1.- Index measurements can be found in:'
 	write(*,*)lick,jrose
 	if(eory.eq.'y')then
 	 write(*,*)'2.- Index error estimates can be found in:'
 	 write(*,*)licker,roser
 	endif
 	if(EWY.eq.'y')then
 	 write(*,*)'3.- File containing index names and values:'
 	 write(*,*)ewtso
 	endif
 	END

 	SUBROUTINE ITCT(x,nn,cont1,cont2,sum)
 	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
 	DIMENSION x(99000,2)
 	sum=0.
 	f1=0
 	f2=0.
c Flux x(i,2) is assumed to be centered on lambda x(i,1)
c Therefore this flux applies to x(i,1)-d*0.5, x(i,1)+d*0.5
 	do j=3,nn-2
		dlan=abs(0.5d0*(x(j+1,1)-x(j-1,1)))
 		sum=sum+x(j,2)*dlan
 	enddo
	dlan=abs(x(2,1)-x(1,1))
	if(cont1.ge.x(2,1)-dlan*0.5d0) then
	 sum=sum+x(2,2)*abs((x(2,1)+dlan*0.5d0)-cont1)
	else
	 sum=sum+x(2,2)*dlan
	 sum=sum+x(1,2)*abs((x(2,1)-dlan*0.5d0)-cont1)
	endif
	dlan=abs(x(nn,1)-x(nn-1,1))
	if(cont2.lt.x(nn-1,1)+dlan*0.5d0)then
	 sum=sum+x(nn-1,2)*abs(cont2-(x(nn-1,1)-dlan*0.5d0))
	else
	 sum=sum+x(nn-1,2)*dlan
	 sum=sum+x(nn,2)*abs(cont2-(x(nn-1,1)+dlan*0.5d0))
	endif
 	return
 	end

      SUBROUTINE polint(xa,ya,n,x,y,dy)
C    Numerical Recipes subroutine
      INTEGER n,NMAX
      DOUBLE PRECISION dy,x,y,xa(n),ya(n)
      PARAMETER (NMAX=10)
      INTEGER i,m,ns
      DOUBLE PRECISION den,dif,dift,ho,hp,w,c(NMAX),d(NMAX)
      ns=1
      dif=abs(x-xa(1))
      do 11 i=1,n
        dift=abs(x-xa(i))
        if (dift.lt.dif) then
          ns=i
          dif=dift
        endif
        c(i)=ya(i)
        d(i)=ya(i)
11    continue
      y=ya(ns)
      ns=ns-1
      do 13 m=1,n-1
        do 12 i=1,n-m
          ho=xa(i)-x
          hp=xa(i+m)-x
          w=c(i+1)-d(i)
          den=ho-hp
          if(den.eq.0.)stop 'failure in polint'
          den=w/den
          d(i)=hp*den
          c(i)=ho*den
12      continue
        if (2*ns.lt.n-m)then
          dy=c(ns+1)
        else
          dy=d(ns)
          ns=ns-1
        endif
        y=y+dy
13    continue
      return
      END

      SUBROUTINE hunt(xx,n,x,jlo)
C subroutine from the Numerical Recipes
      INTEGER jlo,n
      double precision x,xx(n)
      INTEGER inc,jhi,jm
      LOGICAL ascnd
      ascnd=xx(n).ge.xx(1)
      if(jlo.le.0.or.jlo.gt.n)then
        jlo=0
        jhi=n+1
        goto 3
      endif
      inc=1
      if(x.ge.xx(jlo).eqv.ascnd)then
1       jhi=jlo+inc
        if(jhi.gt.n)then
          jhi=n+1
        else if(x.ge.xx(jhi).eqv.ascnd)then
          jlo=jhi
          inc=inc+inc
          goto 1
        endif
      else
        jhi=jlo
2       jlo=jhi-inc
        if(jlo.lt.1)then
          jlo=0
        else if(x.lt.xx(jlo).eqv.ascnd)then
          jhi=jlo
          inc=inc+inc
          goto 2
        endif
      endif
3     if(jhi-jlo.eq.1)then
        if(x.eq.xx(n))jlo=n-1
        if(x.eq.xx(1))jlo=1
        return
      endif
      jm=(jhi+jlo)/2
      if(x.ge.xx(jm).eqv.ascnd)then
        jlo=jm
      else
        jhi=jm
      endif
      goto 3
      END

 	SUBROUTINE LOCT(xx,n,x,ii1,ii2)
C     A slightly modified version of the 
C     subroutine "locate" of the Numerical Recipes
 	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
 	DIMENSION xx(n)
 	jl=0
 	k=0
 	ju=n+1
6111	if(ju-jl.gt.1)then
 		jm=(ju+jl)/2.
 		if((xx(n).gt.xx(1)).eqv.(x.gt.xx(jm)))then
 			jl=jm
 		else
 			ju=jm
 		endif
 		goto 6111
 	endif
 	j=jl
 	ii1=j
 	ii2=j+1
 	return
 	end

 	SUBROUTINE INTERP(x1,x2,y1,y2,x,y)
 	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
 	y=y1+((y2-y1)/(x2-x1))*(x-x1)
 	return
 	end

 	SUBROUTINE MINMAX(k,ixx,npix,ans,xlambda)
 	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
 	COMMON/cmnmx/s(99000,2),nnr
	ans=s(k-npix,2)
	xlambda=s(k-npix,1)
 	if(ixx.eq.0)then
 	  do j=k-npix+1,k+npix
	   ans=min(s(j,2),ans)
	   if(ans.eq.s(j,2)) xlambda=s(j,1)
 	  enddo
 	else
 	  do j=k-npix+1,k+npix
	   ans=max(s(j,2),ans)
	   if(ans.eq.s(j,2)) xlambda=s(j,1)
 	  enddo
 	endif
 	return
 	end

 	subroutine Hd4045(roz,eroz)
C...calculo de Hdelta/4045
 	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
 	COMMON/cmnmx/s(99000,2),nnr
 	COMMON/croz/s1(99000),sf(99000),reds
 	COMMON/eerr/RN,EADU
 	IF((4041.*reds).gt.s(1,1).and.(4108.*reds).lt.s(nnr,1))THEN
 	  tito=4102.875*reds
 	  call LOCT(s1,nnr,tito,k1,k2)
 	  nnnp=int(6.5/abs(s1(k2)-s1(k1)))/2
 	  CALL MINMAX(k1,0,nnnp,aaa,xlambda)
 	  tito=4045.8*reds
 	  call LOCT(s1,nnr,tito,k1,k2)
 	  nnnp=int(6.5/abs(s1(k2)-s1(k1)))/2
 	  CALL MINMAX(k1,0,nnnp,bbb,xlambda)
 	  roz=aaa/bbb
 	saaa=(aaa/EADU)+(RN/EADU)*(RN/EADU)
 	sbbb=(bbb/EADU)+(RN/EADU)*(RN/EADU)
 	era=sqrt(saaa)
 	erb=sqrt(sbbb)
 	eroz=(1./(bbb*bbb))*sqrt(bbb*bbb*era*era+aaa*aaa*erb*erb)
 	ELSE
 	  roz=99.999
 	  eroz=99.999
 	ENDIF
 	return
 	end

 	subroutine Hd4063(roz,eroz)
C...calculo de Hdelta/4063
 	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
 	COMMON/cmnmx/s(99000,2),nnr
 	COMMON/croz/s1(99000),sf(99000),reds
 	COMMON/eerr/RN,EADU
 	IF((4059.*reds).gt.s(1,1).and.(4108.*reds).lt.s(nnr,1))THEN
 	  tito=4102.875*reds
 	  call LOCT(s1,nnr,tito,k1,k2)
 	  nnnp=int(6.5/abs(s1(k2)-s1(k1)))/2
 	  CALL MINMAX(k1,0,nnnp,aaa,xlambda)
 	  tito=4063.6*reds
 	  call LOCT(s1,nnr,tito,k1,k2)
 	  nnnp=int(6.5/abs(s1(k2)-s1(k1)))/2
 	  CALL MINMAX(k1,0,nnnp,bbb,xlambda)
 	  roz=aaa/bbb
 	saaa=(aaa/EADU)+(RN/EADU)*(RN/EADU)
 	sbbb=(bbb/EADU)+(RN/EADU)*(RN/EADU)
 	era=sqrt(saaa)
 	erb=sqrt(sbbb)
 	eroz=(1./(bbb*bbb))*sqrt(bbb*bbb*era*era+aaa*aaa*erb*erb)
 	ELSE
 	  roz=99.999
 	  eroz=99.999
 	ENDIF
 	return
 	end

 	subroutine SrII4045(roz,eroz)
C...calculo de SrII/4045
 	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
 	COMMON/cmnmx/s(99000,2),nnr
 	COMMON/croz/s1(99000),sf(99000),reds
 	COMMON/eerr/RN,EADU
 	IF((4041.*reds).gt.s(1,1).and.(4083.*reds).lt.s(nnr,1))THEN
 	  tito=4077.7*reds
 	  call LOCT(s1,nnr,tito,k1,k2)
 	  nnnp=int(6.5/abs(s1(k2)-s1(k1)))/2
 	  CALL MINMAX(k1,0,nnnp,aaa,xlambda)
 	  tito=4045.8*reds
 	  call LOCT(s1,nnr,tito,k1,k2)
 	  nnnp=int(6.5/abs(s1(k2)-s1(k1)))/2
 	  CALL MINMAX(k1,0,nnnp,bbb,xlambda)
 	  roz=aaa/bbb
 	saaa=(aaa/EADU)+(RN/EADU)*(RN/EADU)
 	sbbb=(bbb/EADU)+(RN/EADU)*(RN/EADU)
 	era=sqrt(saaa)
 	erb=sqrt(sbbb)
 	eroz=(1./(bbb*bbb))*sqrt(bbb*bbb*era*era+aaa*aaa*erb*erb)
 	ELSE
 	  roz=99.999
 	  eroz=99.999
 	ENDIF
 	return
 	end

 	subroutine SrII4063(roz,eroz)
C...calculo de SrII/4063
 	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
 	COMMON/cmnmx/s(99000,2),nnr
 	COMMON/croz/s1(99000),sf(99000),reds
 	COMMON/eerr/RN,EADU
 	IF((4059.*reds).gt.s(1,1).and.(4083.*reds).lt.s(nnr,1))THEN
 	  tito=4077.7*reds
 	  call LOCT(s1,nnr,tito,k1,k2)
 	  nnnp=int(6.5/abs(s1(k2)-s1(k1)))/2
 	  CALL MINMAX(k1,0,nnnp,aaa,xlambda)
 	  tito=4063.6*reds
 	  call LOCT(s1,nnr,tito,k1,k2)
 	  nnnp=int(6.5/abs(s1(k2)-s1(k1)))/2
 	  CALL MINMAX(k1,0,nnnp,bbb,xlambda)
 	  roz=aaa/bbb
 	saaa=(aaa/EADU)+(RN/EADU)*(RN/EADU)
 	sbbb=(bbb/EADU)+(RN/EADU)*(RN/EADU)
 	era=sqrt(saaa)
 	erb=sqrt(sbbb)
 	eroz=(1./(bbb*bbb))*sqrt(bbb*bbb*era*era+aaa*aaa*erb*erb)
 	ELSE
 	  roz=99.999
 	  eroz=99.999
 	ENDIF
 	return
 	end

 	subroutine HgG(roz,eroz)
C...calculo de Hgamma/Gband
 	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
 	COMMON/cmnmx/s(99000,2),nnr
 	COMMON/croz/s1(99000),sf(99000),reds
 	COMMON/eerr/RN,EADU
 	IF((4303.*reds).gt.s(1,1).and.(4347.*reds).lt.s(nnr,1))THEN
 	  tito=4341.625*reds
 	  call LOCT(s1,nnr,tito,k1,k2)
 	  nnnp=int(6.5/abs(s1(k2)-s1(k1)))/2
 	  CALL MINMAX(k1,0,nnnp,aaa,xlambda)
 	  tito=4308.*reds
 	  call LOCT(s1,nnr,tito,k1,k2)
 	  nnnp=int(6.5/abs(s1(k2)-s1(k1)))/2
 	  CALL MINMAX(k1,0,nnnp,bbb,xlambda)
 	  roz=aaa/bbb
 	saaa=(aaa/EADU)+(RN/EADU)*(RN/EADU)
 	sbbb=(bbb/EADU)+(RN/EADU)*(RN/EADU)
 	era=sqrt(saaa)
 	erb=sqrt(sbbb)
 	eroz=(1./(bbb*bbb))*sqrt(bbb*bbb*era*era+aaa*aaa*erb*erb)
 	ELSE
 	  roz=99.999
 	  eroz=99.999
 	ENDIF
 	return
 	end

 	subroutine Hg4325(roz,eroz)
C...calculo de Hgamma/4325
 	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
 	COMMON/cmnmx/s(99000,2),nnr
 	COMMON/croz/s1(99000),sf(99000),reds
 	COMMON/eerr/RN,EADU
 	IF((4321.*reds).gt.s(1,1).and.(4347.*reds).lt.s(nnr,1))THEN
 	  tito=4341.625*reds
 	  call LOCT(s1,nnr,tito,k1,k2)
 	  nnnp=int(6.5/abs(s1(k2)-s1(k1)))/2
 	  CALL MINMAX(k1,0,nnnp,aaa,xlambda)
 	  tito=4325.8*reds
 	  call LOCT(s1,nnr,tito,k1,k2)
 	  nnnp=int(6.5/abs(s1(k2)-s1(k1)))/2
 	  CALL MINMAX(k1,0,nnnp,bbb,xlambda)
 	  roz=aaa/bbb
 	saaa=(aaa/EADU)+(RN/EADU)*(RN/EADU)
 	sbbb=(bbb/EADU)+(RN/EADU)*(RN/EADU)
 	era=sqrt(saaa)
 	erb=sqrt(sbbb)
 	eroz=(1./(bbb*bbb))*sqrt(bbb*bbb*era*era+aaa*aaa*erb*erb)
 	ELSE
 	  roz=99.999
 	  eroz=99.999
 	ENDIF
 	return
 	end

 	subroutine s42894271(roz,eroz)
C...calculo de 4289/4271
 	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
 	COMMON/cmnmx/s(99000,2),nnr
 	COMMON/croz/s1(99000),sf(99000),reds
 	COMMON/eerr/RN,EADU
 	IF((4267.*reds).gt.s(1,1).and.(4295.*reds).lt.s(nnr,1))THEN
 	  tito=4289.7*reds
 	  call LOCT(s1,nnr,tito,k1,k2)
 	  nnnp=int(6.5/abs(s1(k2)-s1(k1)))/2
 	  CALL MINMAX(k1,0,nnnp,aaa,xlambda)
 	  tito=4271.8*reds
 	  call LOCT(s1,nnr,tito,k1,k2)
 	  nnnp=int(6.5/abs(s1(k2)-s1(k1)))/2
 	  CALL MINMAX(k1,0,nnnp,bbb,xlambda)
 	  roz=aaa/bbb
 	saaa=(aaa/EADU)+(RN/EADU)*(RN/EADU)
 	sbbb=(bbb/EADU)+(RN/EADU)*(RN/EADU)
 	era=sqrt(saaa)
 	erb=sqrt(sbbb)
 	eroz=(1./(bbb*bbb))*sqrt(bbb*bbb*era*era+aaa*aaa*erb*erb)
 	ELSE
 	  roz=99.999
 	  eroz=99.999
 	ENDIF
 	return
 	end

 	subroutine s43844352(roz,eroz)
C...calculo de 4384/4352
 	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
 	COMMON/cmnmx/s(99000,2),nnr
 	COMMON/croz/s1(99000),sf(99000),reds
 	COMMON/eerr/RN,EADU
 	IF((4347.*reds).gt.s(1,1).and.(4389.*reds).lt.s(nnr,1))THEN
 	  tito=4383.6*reds
 	  call LOCT(s1,nnr,tito,k1,k2)
 	  nnnp=int(6.5/abs(s1(k2)-s1(k1)))/2
 	  CALL MINMAX(k1,0,nnnp,aaa,xlambda)
 	  tito=4351.9*reds
 	  call LOCT(s1,nnr,tito,k1,k2)
 	  nnnp=int(6.5/abs(s1(k2)-s1(k1)))/2
 	  CALL MINMAX(k1,0,nnnp,bbb,xlambda)
 	  roz=aaa/bbb
 	saaa=(aaa/EADU)+(RN/EADU)*(RN/EADU)
 	sbbb=(bbb/EADU)+(RN/EADU)*(RN/EADU)
 	era=sqrt(saaa)
 	erb=sqrt(sbbb)
 	eroz=(1./(bbb*bbb))*sqrt(bbb*bbb*era*era+aaa*aaa*erb*erb)
 	ELSE
 	  roz=99.999
 	  eroz=99.999
 	ENDIF
 	return
 	end

 	subroutine pFeH(roz,eroz)
C...calculo de p[Fe/H]
 	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
 	COMMON/cmnmx/s(99000,2),nnr
 	COMMON/croz/s1(99000),sf(99000),reds
 	COMMON/eerr/RN,EADU
 	IF((4033.*reds).gt.s(1,1).and.(4094.*reds).lt.s(nnr,1))THEN
 	  tito=4038.*reds
 	  call LOCT(s1,nnr,tito,k1,k2)
 	  nnnp=int(6.5/abs(s1(k2)-s1(k1)))/2
 	  CALL MINMAX(k1,1,nnnp,aaa,xlambda)
 	  tito=4060.*reds
 	  call LOCT(s1,nnr,tito,k1,k2)
 	  nnnp=int(6.5/abs(s1(k2)-s1(k1)))/2
 	  CALL MINMAX(k1,1,nnnp,bbb,xlambda)
 	  tito=4089.*reds
 	  call LOCT(s1,nnr,tito,k1,k2)
 	  nnnp=int(6.5/abs(s1(k2)-s1(k1)))/2
 	  CALL MINMAX(k1,1,nnnp,aaa2,xlambda)
 	  roz=(aaa+aaa2)*0.5/bbb
 	saaa=(aaa/EADU)+(RN/EADU)*(RN/EADU)
 	saaa2=(aaa2/EADU)+(RN/EADU)*(RN/EADU)
 	sbbb=(bbb/EADU)+(RN/EADU)*(RN/EADU)
 	era=sqrt(saaa)
 	era2=sqrt(saaa2)
 	erb=sqrt(sbbb)
 	eroz=(0.5/bbb)*sqrt(era*era+era2*era2+
     &   (erb*erb*(aaa+aaa2)**2.)/(bbb*bbb))
 	ELSE
 	  roz=99.999
 	  eroz=99.999
 	ENDIF
 	return
 	end

 	subroutine CaII(roz,eroz)
C...calculo de CaII
 	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
 	COMMON/cmnmx/s(99000,2),nnr
 	COMMON/croz/s1(99000),sf(99000),reds
 	COMMON/eerr/RN,EADU
 	IF((3928.*reds).gt.s(1,1).and.(3973.*reds).lt.s(nnr,1))THEN
 	  tito=3968.*reds
 	  call LOCT(s1,nnr,tito,k1,k2)
 	  nnnp=int(6.5/abs(s1(k2)-s1(k1)))/2
 	  CALL MINMAX(k1,0,nnnp,aaa,xlambda)
 	  tito=3933.*reds
 	  call LOCT(s1,nnr,tito,k1,k2)
 	  nnnp=int(6.5/abs(s1(k2)-s1(k1)))/2
 	  CALL MINMAX(k1,0,nnnp,bbb,xlambda)
 	  roz=aaa/bbb
 	saaa=(aaa/EADU)+(RN/EADU)*(RN/EADU)
 	sbbb=(bbb/EADU)+(RN/EADU)*(RN/EADU)
 	era=sqrt(saaa)
 	erb=sqrt(sbbb)
 	eroz=(1./(bbb*bbb))*sqrt(bbb*bbb*era*era+aaa*aaa*erb*erb)
 	ELSE
 	  roz=99.999
 	  eroz=99.999
 	ENDIF
 	return
 	end

 	subroutine s38883859(roz,eroz)
C...calculo de 3888/3859
 	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
 	COMMON/cmnmx/s(99000,2),nnr
 	COMMON/croz/s1(99000),sf(99000),reds
 	COMMON/eerr/RN,EADU
 	IF((3856.*reds).gt.s(1,1).and.(3893.*reds).lt.s(nnr,1))THEN
 	  tito=3888.*reds
 	  call LOCT(s1,nnr,tito,k1,k2)
 	  nnnp=int(6.5/abs(s1(k2)-s1(k1)))/2
 	  CALL MINMAX(k1,0,nnnp,aaa,xlambda)
 	  tito=3859.9*reds
 	  call LOCT(s1,nnr,tito,k1,k2)
 	  nnnp=int(6.5/abs(s1(k2)-s1(k1)))/2
 	  CALL MINMAX(k1,0,nnnp,bbb,xlambda)
 	  roz=aaa/bbb
 	saaa=(aaa/EADU)+(RN/EADU)*(RN/EADU)
 	sbbb=(bbb/EADU)+(RN/EADU)*(RN/EADU)
 	era=sqrt(saaa)
 	erb=sqrt(sbbb)
 	eroz=(1./(bbb*bbb))*sqrt(bbb*bbb*era*era+aaa*aaa*erb*erb)
 	ELSE
 	  roz=99.999
 	  eroz=99.999
 	ENDIF
 	return
 	end

 	subroutine p42204208(roz,eroz)
C...calculo de p4220/4208
 	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
 	COMMON/cmnmx/s(99000,2),nnr
 	COMMON/croz/s1(99000),sf(99000),reds
 	COMMON/eerr/RN,EADU
 	IF((4203.*reds).gt.s(1,1).and.(4225.*reds).lt.s(nnr,1))THEN
 	  tito=4220.*reds
 	  call LOCT(s1,nnr,tito,k1,k2)
 	  nnnp=int(6.5/abs(s1(k2)-s1(k1)))/2
 	  CALL MINMAX(k1,1,nnnp,aaa,xlambda)
 	  tito=4208.*reds
 	  call LOCT(s1,nnr,tito,k1,k2)
 	  nnnp=int(6.5/abs(s1(k2)-s1(k1)))/2
 	  CALL MINMAX(k1,1,nnnp,bbb,xlambda)
 	  roz=aaa/bbb
 	saaa=(aaa/EADU)+(RN/EADU)*(RN/EADU)
 	sbbb=(bbb/EADU)+(RN/EADU)*(RN/EADU)
 	era=sqrt(saaa)
 	erb=sqrt(sbbb)
 	eroz=(1./(bbb*bbb))*sqrt(bbb*bbb*era*era+aaa*aaa*erb*erb)
 	ELSE
 	  roz=99.999
 	  eroz=99.999
 	ENDIF
 	return
 	end

 	subroutine FeIHR(roz)
C...calculo de FeI_HR
 	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
 	DIMENSION h(99000,2)
 	COMMON/cmnmx/s(99000,2),nnr
 	COMMON/croz/s1(99000),sf(99000),reds
 	IF((4033.*reds).gt.s(1,1).and.(4055.*reds).lt.s(nnr,1))THEN
 	  tito=4038.*reds
 	  call LOCT(s1,nnr,tito,ii1,kk1)
 	  nnnp=int(6.5/abs(s1(kk1)-s1(ii1)))/2
 	  CALL MINMAX(ii1,1,nnnp,sb,sbl)
 	  tito=4050.*reds
 	  call LOCT(s1,nnr,tito,ii2,kk2)
 	  nnnp=int(6.5/abs(s1(kk2)-s1(ii2)))/2
 	  CALL MINMAX(ii2,1,nnnp,sr,srl)
 	  tito=4045.825
 	  tito1=tito-1.87
 	  tito2=tito+1.87
 	  lh=0
 	  do j=1,nnr
 	    if(j.eq.nnr)then
 	     d=abs(s(j,1)-s(j-1,1))
 	    else
 	     d=abs(s(j+1,1)-s(j,1))
 	    endif
 	    IF(s(j,1).ge.tito1-d.and.s(j,1).le.tito2+d)THEN
 	      lh=lh+1
 	      h(lh,1)=s(j,1)
 	      h(lh,2)=s(j,2)
 	    ENDIF
 	  enddo
 	  roz=0.
 	  do j=1,lh
 	    CALL INTERP(sbl,srl,sb,sr,h(lh,1),conthr)
 	    h(j,2)=(1.-(h(j,2))/conthr)
 	  enddo
 	  CALL ITCT(h,lh,tito1,tito2,roz)
 	ELSE
 	  roz=99.999
 	ENDIF
 	return
 	end

 	subroutine CaIHR(roz)
C...calculo de CaI_HR
 	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
 	DIMENSION h(99000,2)
 	COMMON/cmnmx/s(99000,2),nnr
 	COMMON/croz/s1(99000),sf(99000),reds
 	IF((4215.*reds).gt.s(1,1).and.(4236.*reds).lt.s(nnr,1))THEN
 	  tito=4220.*reds
 	  call LOCT(s1,nnr,tito,ii1,kk1)
 	  nnnp=int(6.5/abs(s1(kk1)-s1(ii1)))/2
 	  CALL MINMAX(ii1,1,nnnp,sb,sbl)
 	  tito=4231.*reds
 	  call LOCT(s1,nnr,tito,ii2,kk2)
 	  nnnp=int(6.5/abs(s1(kk2)-s1(ii2)))/2
 	  CALL MINMAX(ii2,1,nnnp,sr,srl)
 	  tito=4226.740
 	  tito1=tito-1.87
 	  tito2=tito+1.87
 	  lh=0
 	  do j=1,nnr
 	    if(j.eq.nnr)then
 	     d=abs(s(j,1)-s(j-1,1))
 	    else
 	     d=abs(s(j+1,1)-s(j,1))
 	    endif
 	    IF(s(j,1).ge.tito1-d.and.s(j,1).le.tito2+d)THEN
 	      lh=lh+1
 	      h(lh,1)=s(j,1)
 	      h(lh,2)=s(j,2)
 	    ENDIF
 	  enddo
 	  roz=0.
 	  do j=1,lh
 	    CALL INTERP(sbl,srl,sb,sr,h(lh,1),conthr)
 	    h(j,2)=(1.-(h(j,2))/conthr)
 	  enddo
 	  CALL ITCT(h,lh,tito1,tito2,roz)
 	ELSE
 	  roz=99.999
 	ENDIF
 	return
 	end

 	subroutine HgHR(roz)
C...calculo de HgammaHR
 	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
 	DIMENSION h(99000,2)
 	COMMON/cmnmx/s(99000,2),nnr
 	COMMON/croz/s1(99000),sf(99000),reds
 	IF((4328.*reds).gt.s(1,1).and.(4353.*reds).lt.s(nnr,1))THEN
 	  tito=4333.*reds
 	  call LOCT(s1,nnr,tito,ii1,kk1)
 	  nnnp=int(6.5/abs(s1(kk1)-s1(ii1)))/2
 	  CALL MINMAX(ii1,1,nnnp,sb,sbl)
 	  tito=4348.*reds
 	  call LOCT(s1,nnr,tito,ii2,kk2)
 	  nnnp=int(6.5/abs(s1(kk2)-s1(ii2)))/2
 	  CALL MINMAX(ii2,1,nnnp,sr,srl)
 	  tito=4340.477
 	  tito1=tito-1.87
 	  tito2=tito+1.87
 	  lh=0
 	  do j=1,nnr
 	    if(j.eq.nnr)then
 	     d=abs(s(j,1)-s(j-1,1))
 	    else
 	     d=abs(s(j+1,1)-s(j,1))
 	    endif
 	    IF(s(j,1).ge.tito1-d.and.s(j,1).le.tito2+d)THEN
 	      lh=lh+1
 	      h(lh,1)=s(j,1)
 	      h(lh,2)=s(j,2)
 	    ENDIF
 	  enddo
 	  roz=0.
 	  do j=1,lh
 	    CALL INTERP(sbl,srl,sb,sr,h(lh,1),conthr)
 	    h(j,2)=(1.-(h(j,2))/conthr)
 	  enddo
 	  CALL ITCT(h,lh,tito1,tito2,roz)
 	ELSE
 	  roz=99.999
 	ENDIF
 	return
 	end

C******************************************************************************
C******************************************************************************
 	SUBROUTINE NICO(cats,pat,cat)
 	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
 	PARAMETER (NXMAX=20420)
 	COMMON/cmnmx/s(99000,2),nnr
        COMMON/croz/s1(99000),sf(99000),reds
      IF((8460.0d0*reds).gt.s(1,1).and.(8793.0d0*reds).lt.s(nnr,1))THEN
 	  DISP=abs(s(INT(nnr*0.5),1)-s(INT(nnr*0.5-1.),1))
 	  call INDICE(40,nnr,sf,s(1,1),DISP,reds,cats)
 	  call INDICE(39,nnr,sf,s(1,1),DISP,reds,pat)
 	  call INDICE(35,nnr,sf,s(1,1),DISP,reds,cat)
 	ELSE
 	  cats=99.999
 	  pat=99.999
 	  cats=99.999
 	ENDIF
 	return
 	end
C******************************************************************************
 	SUBROUTINE INDICE(NINDEX,NPIXELS,S,STWV,DISP,RVEL,FINDEX)
 	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C Version 19-January-2000
C------------------------------------------------------------------------------
C Copyright N. Cardiel & J. Gorgas, Departamento de Astrofisica
C Universidad Complutense de Madrid, 28040-Madrid, Spain
C E-mail: ncl@astrax.fis.ucm.es or fjg@astrax.fis.ucm.es
C------------------------------------------------------------------------------
C This program is free software; you can redistribute it and/or modify it
C under the terms of the GNU General Public License as published by the Free
C Software Foundation; either version 2 of the License, or (at your option) any
C later version. See the file gnu-public-license.txt for details.
C******************************************************************************
C This subroutine computes a line-strength index in a single spectrum.
C INPUT:
C	 NINDEX: number of index to be computed (see list below)
C	 S(NPIXELS): spectrum to be measured
C	 STWV: central wavelength of first pixel
C	 DISP: dispersion (Angs./pixel)   Note: the scale must be linear
C	 RVEL: spectrum radial velocity (km/sec)
C OUTPUT:
C	 FINDEX: measured index
C------------------------------------------------------------------------------
C List of defined line-strength indices (see also subroutine INDEXDEF)
C (01) D4000...: 4000-Angs. break (Bruzual 1983)
C (02) B4000...: modified 4000-Angs. break (Gorgas et al. 1998)
C (03) HdA.....: Hdelta A (Worthey & Ottaviani 1997)
C (04) HdF.....: Hdelta F (Worthey & Ottaviani 1997)
C (05) CN1.....: Lick
C (06) CN2.....: Lick
C (07) Ca4227..: Lick
C (08) G4300...: Lick
C (09) HgA.....: Hgamma A (Worthey & Ottaviani 1997)
C (10) HgF.....: Hgamma F (Worthey & Ottaviani 1997)
C (11) Fe4383..: Lick
C (12) Ca4455..: Lick
C (13) Fe4531..: Lick
C (14) Fe4668..: Lick
C (15) Hbeta...: Lick
C (16) Hbeta_p.: Hbeta plus from Gonzalez thesis (p116)
C (17) OIII_1..: OIII_1 from Gonzales thesis (p116)
C (18) OIII_2..: OIII_2 from Gonzales thesis (p116)
C (19) Fe5015..: Lick
C (20) Mg1.....: Lick
C (21) Mg2.....: Lick
C (22) Mgb5177.: Lick
C (23) Fe5270..: Lick
C (24) Fe5335..: Lick
C (25) Fe5406..: Lick
C (26) Fe5709..: Lick
C (27) Fe5782..: Lick
C (28) Na5895..: Lick
C (29) TiO1....: Lick
C (30) TiO2....: Lick
C (31) Ca1(DTT): Diaz, Terlevich & Terlevich (1989)
C (32) Ca2(DTT): Diaz, Terlevich & Terlevich (1989)
C (33) Ca3(DTT): Diaz, Terlevich & Terlevich (1989)
C (34) MgI(DTT): Diaz, Terlevich & Terlevich (1989)
C (35) CaT.....: CGCVP et al. (2000)
C (36) CaT1....: CGCVP et al. (2000)
C (37) CaT2....: CGCVP et al. (2000)
C (38) CaT3....: CGCVP et al. (2000)
C (39) PaT.....: CGCVP et al. (2000)
C (40) CaT*....: CGCVP et al. (2000)
C------------------------------------------------------------------------------
C parametros dimensionales
 	INTEGER NXMAX		   !maximum dimension in the spectral direction
 	PARAMETER (NXMAX=20420)
C
 	INTEGER NBDMAX        !maximum no. of bands allowed for a generic index
 	PARAMETER (NBDMAX=198)  		     !note that NBDMAX=NWVMAX/2
C
 	INTEGER NWVMAX      !2*maximum no. of bands allowed for a generic index
 	PARAMETER (NWVMAX=396)
C constantes
C	 REAL C 					  !speed of light (km/s)
 	PARAMETER (C=2.9979246E+5)
C------------------------------------------------------------------------------
C parametros de la subrutina
 	INTEGER NINDEX
 	INTEGER NPIXELS
 	DIMENSION S(NXMAX)
C------------------------------------------------------------------------------
C variables locales:
 	INTEGER ITI						!tipo de indice
 	INTEGER J
 	INTEGER NBAND,NB
 	INTEGER J1(NBDMAX),J2(NBDMAX)
 	INTEGER J1MIN,J2MAX	 !calcula los limites reales del indice a medir
 	INTEGER NCEFF
 	INTEGER NCONTI,NABSOR
 	DIMENSION WV(NWVMAX)				       !l.d.o. de las bandas
 	DIMENSION  FWV(NWVMAX/4)   !constantes para multiplicar se\~{n}al de asbsorc.
 	DIMENSION  SS(NXMAX)
 	DIMENSION D1(NBDMAX),D2(NBDMAX),RL(NBDMAX),RG(NBDMAX)
 	DIMENSION FX(NBDMAX)					!flujo en las bandas
 	DIMENSION  SC(NXMAX)
 	DOUBLE PRECISION MWB,MWR
 	DIMENSION  WL(NXMAX)					 !pesos para el D4000
 	CHARACTER*8 INDEXNAME
 	DOUBLE PRECISION SMEAN
 	LOGICAL IFCHAN(NXMAX)	!indica los canales usados para medir el indice
C
C------------------------------------------------------------------------------
C------------------------------------------------------------------------------
C proteccion
 	IF(NPIXELS.GT.NXMAX)THEN
 	  WRITE(*,101) 'FATAL ERROR: NPIXELS.GT.NXMAX'
 	  WRITE(*,101) '=> Enlarge NXMAX and recompile the program.'
 	  STOP
 	END IF
C------------------------------------------------------------------------------
 	FINDEX=0.			   !hasta que se demuestre lo contrario
C Definimos el indice que queremos medir
 	CALL INDEXDEF(NINDEX,ITI,INDEXNAME,WV,FWV)
C Numero de bandas, limites y anchura teniendo en cuenta la velocidad radial
 	NBAND=3 	      !indices atomicos y moleculares (lo mas probable)
 	IF(ITI.EQ.3) NBAND=2			!el D4000 solo tiene dos bandas
 	IF(ITI.EQ.4) NBAND=2			!el B4000 solo tiene dos bandas
 	IF((ITI.GE.101).AND.(ITI.LE.9999))THEN     !para los indices genericos:
 	  NCONTI=(ITI/100)			  !numero de bandas de continuo
 	  NABSOR=ITI-NCONTI*100 	       !numero de bandas de absorciones
 	  NBAND=NCONTI+NABSOR
 	END IF
C..............................................................................
 	WLMIN=STWV-DISP*0.5d0
c Me traigo la correccion relativista calculada en el 
c programa principal (LECTOR):
	RCVEL1=RVEL
c y anulamos el calculo original de Cardiel:
c 	RCVEL=RVEL/C
c 	RCVEL1=1.+RCVEL
c 	RCVEL1=RCVEL1/SQRT(1.-RCVEL*RCVEL)		!correccion relativista
 	DO NB=1,NBAND					       !para cada banda
 	  CA=WV(2*NB-1)*RCVEL1  			 !redshifted wavelength
 	  CB=WV(2*NB)*RCVEL1				 !redshifted wavelength
 	  C3=(CA-WLMIN)/DISP+1.d0			    !band limit (channel)
 	  C4=(CB-WLMIN)/DISP				  !band limit (channel)
 	  IF((C3.LT.1.d0).OR.(C4.GT.dble(NPIXELS)))THEN       !index out of range
 	    WRITE(*,101) 'ERROR: index out of range in subroutine '//
     >       'INDICES.'
 	    WRITE(*,100) 'Press <CR> to continue...'
 	    READ(*,*)
 	    RETURN
 	  END IF
 	  J1(NB)=INT(C3)			  !band limit: integer(channel)
 	  J2(NB)=INT(C4)			  !band limit: integer(channel)
 	  D1(NB)=C3-dble(J1(NB))	    !fraction (excess) of first channel
 	  D2(NB)=C4-dble(J2(NB))	     !fraction (defect) of last channel
 	  RL(NB)=CB-CA  		     !redshifted band width (angstroms)
 	  RG(NB)=C4-C3  		      !redshifted band width (channels)
 	END DO
C..............................................................................
 	J1MIN=J1(1)	!calculamos limites por si las bandas no estan en orden
 	J2MAX=J2(1)
 	DO NB=2,NBAND
 	  IF(J1MIN.GT.J1(NB)) J1MIN=J1(NB)
 	  IF(J2MAX.LT.J2(NB)) J2MAX=J2(NB)
 	END DO
C------------------------------------------------------------------------------
C Fijamos los canales a usar para medir el indice (utilizando la variable
C logica evitamos el problema de la posible superposicion de las bandas)
 	DO J=1,NPIXELS
 	  IFCHAN(J)=.FALSE.		 !inicializamos: ningun canal utilizado
 	END DO
C
 	DO NB=1,NBAND					 !recorremos las bandas
 	  DO J=J1(NB),J2(NB)+1
 	    IFCHAN(J)=.TRUE.
 	  END DO
 	END DO
C
 	NCEFF=0 	    !contamos el numero de canales effectivo a utilizar
 	DO J=1,NPIXELS
 	  IF(IFCHAN(J)) NCEFF=NCEFF+1
 	END DO
C------------------------------------------------------------------------------
C Normalizamos espectro a medir, usando la se\~{n}al solo en la region
C del indice
 	SMEAN=0.D0
 	DO J=1,NPIXELS
 	  IF(IFCHAN(J)) SMEAN=SMEAN+DBLE(S(J))      !solo canales en las bandas
 	END DO
 	SMEAN=SMEAN/DBLE(NCEFF) 	     !valor promedio (DOUBLE PRECISION)
 	FSMEAN=dble(SMEAN)				 !valor promedio (REAL)
C
 	DO J=1,NPIXELS
 	  SS(J)=S(J)/FSMEAN				 !normalizamos espectro
 	END DO
C------------------------------------------------------------------------------
C Calculamos pseudo continuo en indices moleculares y atomicos (ITI=1,2)
C (formulas en Tesis de JJGG, pag. 35)
C
C NOTA: al transformar las integrales en sumatorios, habria que multiplicar
C cada valor de la funcion a integrar por el incremento en longitud de onda
C (que en el sumatorio coincide con DISP). Sin embargo, este valor es factor
C comun y puede salir fuera del sumatorio, por lo que solo hace falta
C introducirlo al final. Asimismo, como al calcular los limites de las bandas,
C J1() y J2(), hemos tenido presente la vel. radial, la anchura de las bandas,
C RL(), es mayor en un factor (1+z) que la anchura cuando el objeto no presenta
C velocidad radial. Es decir, el sumatorio se extiende sobre una region (en
C longitud de onda) algo mayor. Sin embargo el efecto queda anulado al dividir
C por RL() que se encuentra ensanchado en el mismo factor.
C
 	IF((ITI.EQ.1).OR.(ITI.EQ.2))THEN
C..............................................................................
 	  SB=0.d0			       !cuentas promedio en la banda azul
 	  DO J=J1(1),J2(1)+1			      !recorremos la banda azul
 	    IF(J.EQ.J1(1))THEN
 	      F=1.-D1(1)
 	    ELSEIF(J.EQ.J2(1)+1)THEN
 	      F=D2(1)
 	    ELSE
 	      F=1.d0
 	    END IF
 	    SB=SB+F*SS(J)
 	  END DO
 	  SB=SB*DISP					 !completamos sumatorio
 	  SB=SB/RL(1)			!dividimos por anchura de la banda azul
C..............................................................................
 	  SR=0.d0			       !cuentas promedio en la banda roja
 	  DO J=J1(3),J2(3)+1			      !recorremos la banda roja
 	    IF(J.EQ.J1(3))THEN
 	      F=1.-D1(3)
 	    ELSEIF(J.EQ.J2(3)+1)THEN
 	      F=D2(3)
 	    ELSE
 	      F=1.d0
 	    END IF
 	    SR=SR+F*SS(J)
 	  END DO
 	  SR=SR*DISP					 !completamos sumatorio
 	  SR=SR/RL(3)			!dividimos por anchura de la banda roja
C..............................................................................
C Trabajamos en la escala en l.d.o. sin corregir de Vrad (es decir, con las
C bandas de los indices desplazadas a Vrad correspondiente). Se obtiene lo
C mismo si mantenemos las bandas de los indices a Vrad=0 pero al calcular WLA
C para cada canal J dividimos por RCVEL1 para obtener la escala en l.d.o.
C corregida de Vrad.
 	  MWB=(WV(1)+WV(2))/2.d0	     !mean wavelength blue band at Vrad=0
 	  MWB=MWB*RCVEL1			       !idem a Vrad considerado
 	  MWR=(WV(5)+WV(6))/2.d0	      !mean wavelength red band at Vrad=0
 	  MWR=MWR*RCVEL1			       !idem a Vrad considerado
C Calculamos el valor del pseudo continuo desde el borde mas azul de todas las
C bandas al borde mas rojo de todas las bandas (no importa que las banda
C "azul" no sea la mas azul, etc.)
 	  DO J=J1MIN,J2MAX+1
 	    WLA=dble(J-1)*DISP+STWV			    !l.d.o. del canal J
 	    SC(J)=(SB*(MWR-WLA)+SR*(WLA-MWB))/(MWR-MWB)
 	  END DO
C..............................................................................
 	END IF
C------------------------------------------------------------------------------
C Pesos para el D4000: debido a que la integral es  F(ldo)*d(nu), hay que
C multiplicar el flujo por el cuadrado de la longitud de onda, y de este
C modo la integral se transforma en F(ldo)*d(ldo).
C
C NOTA: estamos corrigiendo WLA de Vrad, aunque luego medimos el indice con
C los limites de las bandas multiplicados por (1+z). Esto no es importante
C porque solo hay un factor Cte=(1+z) en los pesos que no influye en el
C indice. Sin embargo, vamos a trabajar con los pesos asi porque de esta
C manera la normalizacion a 4000 siempre nos proporciona pesos cercanos a uno
C independientemente del valor de Vrad.
 	IF(ITI.EQ.3)THEN
 	  DO J=1,NPIXELS
 	    WLA=dble(J-1)*DISP+STWV			    !l.d.o. del canal J
 	    WLA=WLA/RCVEL1		 !l.d.o. del canal J corrigiendo a Vrad
 	    WLA=WLA/4000.d0	 !normalizamos a 4000 para tener todo proximo a 1
 	    WL(J)=WLA*WLA
 	  END DO
 	END IF
C------------------------------------------------------------------------------
C Para el B4000 usaremos codigo comun con el D4000, donde los pesos son iguales
C a uno.
 	IF(ITI.EQ.4)THEN
 	  DO J=1,NPIXELS
 	    WL(J)=1.0d0
 	  END DO
 	END IF
C------------------------------------------------------------------------------
C Para los indices genericos calculamos el pseudocontinuo ajustando por minimos
C cuadrados (pesando con errores si procede) a todos los pixels de las bandas
C de continuo.
 	IF((ITI.GE.101).AND.(ITI.LE.9999))THEN
C..............................................................................
C calculamos la recta del continuo mediante minimos cuadrados (la recta es
C de la forma y= amc * x + bmc
C (NOTA: para la variable "x" utilizamos el valor del numero de pixel en lugar
C de la longitud de onda porque, en principio, son numeros mas pequenhos)
 	  SUM0=0.d0
 	  SUMX=0.d0
 	  SUMY=0.d0
 	  SUMXY=0.d0
 	  SUMXX=0.d0
 	  DO NB=1,NCONTI	       !recorremos todas las bandas de continuo
 	    DO J=J1(NB),J2(NB)+1  !recorremos todos los pixels de dichas bandas
 	      IF(J.EQ.J1(NB))THEN	 !comprobamos efecto de borde izquierdo
 		F=1.-D1(NB)
 	      ELSEIF(J.EQ.J2(NB)+1)THEN 	       !efecto de borde derecho
 		F=D2(NB)
 	      ELSE				    !pixels sin efecto de borde
 		F=1.d0
 	      END IF
 	      SUM0=SUM0+F
 	      SUMX=SUMX+F*dble(J)
 	      SUMY=SUMY+F*SS(J)
 	      SUMXY=SUMXY+F*dble(J)*SS(J)
 	      SUMXX=SUMXX+F*dble(J)*dble(J)
 	    END DO
 	  END DO
 	  DETER=SUM0*SUMXX-SUMX*SUMX
 	  AMC=(SUM0*SUMXY-SUMX*SUMY)/DETER
 	  BMC=(SUMXX*SUMY-SUMX*SUMXY)/DETER
C..............................................................................
C calculamos el pseudocontinuo desde el borde mas azul de todas las bandas
C hasta el borde mas rojo
 	  DO J=J1MIN,J2MAX+1
 	    SC(J)=AMC*dble(J)+BMC
 	  END DO
C..............................................................................
 	END IF
C------------------------------------------------------------------------------
C------------------------------------------------------------------------------
C medimos indices
 	IF((ITI.EQ.1).OR.(ITI.EQ.2))THEN	!indices moleculares y atomicos
 	  TC=0.d0
 	  DO J=J1(2),J2(2)+1			   !recorremos la banda central
 	    IF(J.EQ.J1(2))THEN
 	      F=1.-D1(2)
 	    ELSEIF(J.EQ.J2(2)+1)THEN
 	      F=D2(2)
 	    ELSE
 	      F=1.d0
 	    END IF
 	    TC=TC+F*SS(J)/SC(J)
 	  END DO
 	  TC=TC*DISP					 !completamos sumatorio
 	  IF(ITI.EQ.1)THEN				      !indice molecular
C	    FINDEX=-2.5d0*ALOG10(TC/RL(2))
 	    FINDEX=-2.5d0*DLOG10(TC/RL(2))
 	  ELSEIF(ITI.EQ.2)THEN  				!indice atomico
 	    FINDEX=RL(2)-TC
 	    FINDEX=FINDEX/RCVEL1		       !correccion por redshift
 	  END IF
C------------------------------------------------------------------------------
 	ELSEIF((ITI.EQ.3).OR.(ITI.EQ.4))THEN				 !D4000
C NOTA: el sumatorio TC no incluye el factor DISP, que corresponderia
C con el incremento (diferencial en la integral), dado que al ser un factor
C constante tampoco altera el resultado a la hora de computar el D4000 (o el
C B4000).
 	  DO NB=1,NBAND 				       !bucle en bandas
 	    TC=0.d0
 	    DO J=J1(NB),J2(NB)+1			!recorremos la banda NB
 	      IF(J.EQ.J1(NB))THEN
 		F=1.-D1(NB)
 	      ELSEIF(J.EQ.J2(NB)+1)THEN
 		F=D2(NB)
 	      ELSE
 		F=1
 	      END IF
 	      TC=TC+F*SS(J)*WL(J)
 	    END DO
 	    FX(NB)=TC
 	  END DO
 	  FINDEX=FX(2)/FX(1)
C------------------------------------------------------------------------------
 	ELSEIF((ITI.GE.101).AND.(ITI.LE.9999))THEN
C recorremos las bandas con absorciones
 	  TC=0.d0
 	  SUMRL=0.d0
 	  DO NB=1,NABSOR
 	    DO J=J1(NCONTI+NB),J2(NCONTI+NB)+1     !recorremos todos los pixels
 	      IF(J.EQ.J1(NCONTI+NB))THEN !comprobamos efecto de borde izquierdo
 		F=1.-D1(NCONTI+NB)
 	      ELSEIF(J.EQ.J2(NCONTI+NB)+1)THEN         !efecto de borde derecho
 		F=D2(NCONTI+NB)
 	      ELSE				    !pixels sin efecto de borde
 		F=1.d0
 	      END IF
 	      TC=TC+F*FWV(NB)*SS(J)/SC(J)   !multiplicamos absorcion por factor
 	    END DO
 	    SUMRL=SUMRL+FWV(NB)*RL(NCONTI+NB)
 	  END DO
 	  TC=TC*DISP				      !completamos el sumatorio
 	  FINDEX=SUMRL-TC	  !medimos el indice generico como los atomicos
 	  FINDEX=FINDEX/RCVEL1  		       !correccion por redshift
C------------------------------------------------------------------------------
 	END IF
C------------------------------------------------------------------------------
C Fin de subrutina
100	FORMAT(A,$)
101	FORMAT(A)
 	END
C
C******************************************************************************
C Retorna las "bandpasses" del indice numero NINDEX, a la vez que el tipo de
C indice:
C ITI=3 o 4: D4000 o B4000
C ITI=2: indice atomico
C ITI=1: indice molecular
C------------------------------------------------------------------------------
 	SUBROUTINE INDEXDEF(NINDEX,ITI,NAME,WV,FWV)
 	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C	IMPLICIT NONE
C parametros generales
 	INTEGER NWVMAX      !2*maximum no. of bands allowed for a generic index
 	PARAMETER (NWVMAX=396)  		     !see subroutine selindex.f
C parametros de la subrutina
 	INTEGER NINDEX
 	INTEGER ITI
 	CHARACTER*8 NAME
 	DIMENSION WV(NWVMAX)				       !l.d.o. de las bandas
 	DIMENSION FWV(NWVMAX/4)   !constantes para multiplicar se\~{n}al de asbsorc.
C------------------------------------------------------------------------------
 	IF(NINDEX.EQ.1)THEN
 	  NAME='D4000'
 	  ITI=3
 	  WV(1)=3750.000d0
 	  WV(2)=3950.000d0
 	  WV(3)=4050.000d0
 	  WV(4)=4250.000d0
 	ELSEIF(NINDEX.EQ.2)THEN
 	  NAME='B4000'
 	  ITI=4
 	  WV(1)=3750.000d0
 	  WV(2)=3950.000d0
 	  WV(3)=4050.000d0
 	  WV(4)=4250.000d0
 	ELSEIF(NINDEX.EQ.3)THEN
 	  NAME='HdA'
 	  ITI=2
 	  WV(1)=4041.600d0
 	  WV(2)=4079.750d0
 	  WV(3)=4083.500d0
 	  WV(4)=4122.250d0
 	  WV(5)=4128.500d0
 	  WV(6)=4161.000d0
 	ELSEIF(NINDEX.EQ.4)THEN
 	  NAME='HdF'
 	  ITI=2
 	  WV(1)=4057.250d0
 	  WV(2)=4088.500d0
 	  WV(3)=4091.000d0
 	  WV(4)=4112.250d0
 	  WV(5)=4114.750d0
 	  WV(6)=4137.250d0
 	ELSEIF(NINDEX.EQ.5)THEN
 	  NAME='CN1'
 	  ITI=1
 	  WV(1)=4080.125d0
 	  WV(2)=4117.625d0
 	  WV(3)=4142.125d0
 	  WV(4)=4177.125d0
 	  WV(5)=4244.125d0
 	  WV(6)=4284.125d0
 	ELSEIF(NINDEX.EQ.6)THEN
 	  NAME='CN2'
 	  ITI=1
 	  WV(1)=4083.875d0
 	  WV(2)=4096.375d0
 	  WV(3)=4142.125d0
 	  WV(4)=4177.125d0
 	  WV(5)=4244.125d0
 	  WV(6)=4284.125d0
 	ELSEIF(NINDEX.EQ.7)THEN
 	  NAME='Ca4227'
 	  ITI=2
 	  WV(1)=4211.000d0
 	  WV(2)=4219.750d0
 	  WV(3)=4222.250d0
 	  WV(4)=4234.750d0
 	  WV(5)=4241.000d0
 	  WV(6)=4251.000d0
 	ELSEIF(NINDEX.EQ.8)THEN
 	  NAME='G4300'
 	  ITI=2
 	  WV(1)=4266.375d0
 	  WV(2)=4282.625d0
 	  WV(3)=4281.375d0
 	  WV(4)=4316.375d0
 	  WV(5)=4318.875d0
 	  WV(6)=4335.125d0
 	ELSEIF(NINDEX.EQ.9)THEN
 	  NAME='HgA'
 	  ITI=2
 	  WV(1)=4283.500d0
 	  WV(2)=4319.750d0
 	  WV(3)=4319.750d0
 	  WV(4)=4363.500d0
 	  WV(5)=4367.250d0
 	  WV(6)=4419.750d0
 	ELSEIF(NINDEX.EQ.10)THEN
 	  NAME='HgF'
 	  ITI=2
 	  WV(1)=4283.500d0
 	  WV(2)=4319.750d0
 	  WV(3)=4331.250d0
 	  WV(4)=4352.250d0
 	  WV(5)=4354.750d0
 	  WV(6)=4384.750d0
 	ELSEIF(NINDEX.EQ.11)THEN
 	  NAME='Fe4383'
 	  ITI=2
 	  WV(1)=4359.125d0
 	  WV(2)=4370.375d0
 	  WV(3)=4369.125d0
 	  WV(4)=4420.375d0
 	  WV(5)=4442.875d0
 	  WV(6)=4455.375d0
 	ELSEIF(NINDEX.EQ.12)THEN
 	  NAME='Ca4455'
 	  ITI=2
 	  WV(1)=4445.875d0
 	  WV(2)=4454.625d0
 	  WV(3)=4452.125d0
 	  WV(4)=4474.625d0
 	  WV(5)=4477.125d0
 	  WV(6)=4492.125d0
 	ELSEIF(NINDEX.EQ.13)THEN
 	  NAME='Fe4531'
 	  ITI=2
 	  WV(1)=4504.250d0
 	  WV(2)=4514.250d0
 	  WV(3)=4514.250d0
 	  WV(4)=4559.250d0
 	  WV(5)=4560.500d0
 	  WV(6)=4579.250d0
 	ELSEIF(NINDEX.EQ.14)THEN
 	  NAME='Fe4668'
 	  ITI=2
 	  WV(1)=4611.500d0
 	  WV(2)=4630.250d0
 	  WV(3)=4634.000d0
 	  WV(4)=4720.250d0
 	  WV(5)=4742.750d0
 	  WV(6)=4756.500d0
 	ELSEIF(NINDEX.EQ.15)THEN
 	  NAME='Hbeta'
 	  ITI=2
 	  WV(1)=4827.875d0
 	  WV(2)=4847.875d0
 	  WV(3)=4847.875d0
 	  WV(4)=4876.625d0
 	  WV(5)=4876.625d0
 	  WV(6)=4891.625d0
 	ELSEIF(NINDEX.EQ.16)THEN
 	  NAME='Hbeta_p'
 	  ITI=2
 	  WV(1)=4815.000d0
 	  WV(2)=4845.000d0
 	  WV(3)=4851.320d0
 	  WV(4)=4871.320d0
 	  WV(5)=4880.000d0
 	  WV(6)=4930.000d0
 	ELSEIF(NINDEX.EQ.17)THEN
 	  NAME='OIII_1'
 	  ITI=2
 	  WV(1)=4885.000d0
 	  WV(2)=4935.000d0
 	  WV(3)=4948.920d0
 	  WV(4)=4978.920d0
 	  WV(5)=5030.000d0
 	  WV(6)=5070.000d0
 	ELSEIF(NINDEX.EQ.18)THEN
 	  NAME='OIII_2'
 	  ITI=2
 	  WV(1)=4885.000d0
 	  WV(2)=4935.000d0
 	  WV(3)=4996.850d0
 	  WV(4)=5016.850d0
 	  WV(5)=5030.000d0
 	  WV(6)=5070.000d0
 	ELSEIF(NINDEX.EQ.19)THEN
 	  NAME='Fe5015'
 	  ITI=2
 	  WV(1)=4946.500d0
 	  WV(2)=4977.750d0
 	  WV(3)=4977.750d0
 	  WV(4)=5054.000d0
 	  WV(5)=5054.000d0
 	  WV(6)=5065.250d0
 	ELSEIF(NINDEX.EQ.20)THEN
 	  NAME='Mg1'
 	  ITI=1
 	  WV(1)=4895.125d0
 	  WV(2)=4957.625d0
 	  WV(3)=5069.125d0
 	  WV(4)=5134.125d0
 	  WV(5)=5301.125d0
 	  WV(6)=5366.125d0
 	ELSEIF(NINDEX.EQ.21)THEN
 	  NAME='Mg2'
 	  ITI=1
 	  WV(1)=4895.125d0
 	  WV(2)=4957.625d0
 	  WV(3)=5154.125d0
 	  WV(4)=5196.625d0
 	  WV(5)=5301.125d0
 	  WV(6)=5366.125d0
 	ELSEIF(NINDEX.EQ.22)THEN
 	  NAME='Mgb5177'
 	  ITI=2
 	  WV(1)=5142.625d0
 	  WV(2)=5161.375d0
 	  WV(3)=5160.125d0
 	  WV(4)=5192.625d0
 	  WV(5)=5191.375d0
 	  WV(6)=5206.375d0
 	ELSEIF(NINDEX.EQ.23)THEN
 	  NAME='Fe5270'
 	  ITI=2
 	  WV(1)=5233.150d0
 	  WV(2)=5248.150d0
 	  WV(3)=5245.650d0
 	  WV(4)=5285.650d0
 	  WV(5)=5285.650d0
 	  WV(6)=5318.150d0
 	ELSEIF(NINDEX.EQ.24)THEN
 	  NAME='Fe5335'
 	  ITI=2
 	  WV(1)=5304.625d0
 	  WV(2)=5315.875d0
 	  WV(3)=5312.125d0
 	  WV(4)=5352.125d0
 	  WV(5)=5353.375d0
 	  WV(6)=5363.375d0
 	ELSEIF(NINDEX.EQ.25)THEN
 	  NAME='Fe5406'
 	  ITI=2
 	  WV(1)=5376.250d0
 	  WV(2)=5387.500d0
 	  WV(3)=5387.500d0
 	  WV(4)=5415.000d0
 	  WV(5)=5415.000d0
 	  WV(6)=5425.000d0
 	ELSEIF(NINDEX.EQ.26)THEN
 	  NAME='Fe5709'
 	  ITI=2
 	  WV(1)=5672.875d0
 	  WV(2)=5696.625d0
 	  WV(3)=5696.625d0
 	  WV(4)=5720.375d0
 	  WV(5)=5722.875d0
 	  WV(6)=5736.625d0
 	ELSEIF(NINDEX.EQ.27)THEN
 	  NAME='Fe5782'
 	  ITI=2
 	  WV(1)=5765.375d0
 	  WV(2)=5775.375d0
 	  WV(3)=5776.625d0
 	  WV(4)=5796.625d0
 	  WV(5)=5797.875d0
 	  WV(6)=5811.625d0
 	ELSEIF(NINDEX.EQ.28)THEN
 	  NAME='Na5895'
 	  ITI=2
 	  WV(1)=5860.625d0
 	  WV(2)=5875.625d0
 	  WV(3)=5876.875d0
 	  WV(4)=5909.375d0
 	  WV(5)=5922.125d0
 	  WV(6)=5948.125d0
 	ELSEIF(NINDEX.EQ.29)THEN
 	  NAME='TiO1'
 	  ITI=1
 	  WV(1)=5816.625d0
 	  WV(2)=5849.125d0
 	  WV(3)=5936.625d0
 	  WV(4)=5994.125d0
 	  WV(5)=6038.625d0
 	  WV(6)=6103.625d0
 	ELSEIF(NINDEX.EQ.30)THEN
 	  NAME='TiO2'
 	  ITI=1
 	  WV(1)=6066.625d0
 	  WV(2)=6141.625d0
 	  WV(3)=6189.625d0
 	  WV(4)=6272.125d0
 	  WV(5)=6372.625d0
 	  WV(6)=6415.125d0
 	ELSEIF(NINDEX.EQ.31)THEN
 	  NAME='Ca1(DTT)'
 	  ITI=2
 	  WV(1)=8447.500d0
 	  WV(2)=8462.500d0
 	  WV(3)=8483.000d0
 	  WV(4)=8513.000d0
 	  WV(5)=8842.500d0
 	  WV(6)=8857.500d0
 	ELSEIF(NINDEX.EQ.32)THEN
 	  NAME='Ca2(DTT)'
 	  ITI=2
 	  WV(1)=8447.500d0
 	  WV(2)=8462.500d0
 	  WV(3)=8527.000d0
 	  WV(4)=8557.000d0
 	  WV(5)=8842.500d0
 	  WV(6)=8857.500d0
 	ELSEIF(NINDEX.EQ.33)THEN
 	  NAME='Ca3(DTT)'
 	  ITI=2
 	  WV(1)=8447.500d0
 	  WV(2)=8462.500d0
 	  WV(3)=8647.000d0
 	  WV(4)=8677.000d0
 	  WV(5)=8842.500d0
 	  WV(6)=8857.500d0
 	ELSEIF(NINDEX.EQ.34)THEN
 	  NAME='MgI(DTT)'
 	  ITI=2
 	  WV(1)=8775.000d0
 	  WV(2)=8787.000d0
 	  WV(3)=8799.500d0
 	  WV(4)=8814.500d0
 	  WV(5)=8845.000d0
 	  WV(6)=8855.000d0
 	ELSEIF(NINDEX.EQ.35)THEN
 	  NAME='CaT'
 	  ITI=503
 	  WV(1)=8474.000d0
 	  WV(2)=8484.000d0
 	  WV(3)=8563.000d0
 	  WV(4)=8577.000d0
 	  WV(5)=8619.000d0
 	  WV(6)=8642.000d0
 	  WV(7)=8700.000d0
 	  WV(8)=8725.000d0
 	  WV(9)=8776.000d0
 	  WV(10)=8792.00d0
 	  WV(11)=8484.00d0
 	  WV(12)=8513.00d0
 	  FWV(1)=1.0d0
 	  WV(13)=8522.00d0
 	  WV(14)=8562.00d0
 	  FWV(2)=1.0d0
 	  WV(15)=8642.00d0
 	  WV(16)=8682.00d0
 	  FWV(3)=1.0d0
 	ELSEIF(NINDEX.EQ.36)THEN
 	  NAME='CaT1'
 	  ITI=501
 	  WV(1)=8474.000d0
 	  WV(2)=8484.000d0
 	  WV(3)=8563.000d0
 	  WV(4)=8577.000d0
 	  WV(5)=8619.000d0
 	  WV(6)=8642.000d0
 	  WV(7)=8700.000d0
 	  WV(8)=8725.000d0
 	  WV(9)=8776.000d0
 	  WV(10)=8792.00d0
 	  WV(11)=8484.00d0
 	  WV(12)=8513.00d0
 	  FWV(1)=1.0d0
 	ELSEIF(NINDEX.EQ.37)THEN
 	  NAME='CaT2'
 	  ITI=501
 	  WV(1)=8474.000d0
 	  WV(2)=8484.000d0
 	  WV(3)=8563.000d0
 	  WV(4)=8577.000d0
 	  WV(5)=8619.000d0
 	  WV(6)=8642.000d0
 	  WV(7)=8700.000d0
 	  WV(8)=8725.000d0
 	  WV(9)=8776.000d0
 	  WV(10)=8792.00d0
 	  WV(11)=8522.00d0
 	  WV(12)=8562.00d0
 	  FWV(1)=1.0d0
 	ELSEIF(NINDEX.EQ.38)THEN
 	  NAME='CaT3'
 	  ITI=501
 	  WV(1)=8474.000d0
 	  WV(2)=8484.000d0
 	  WV(3)=8563.000d0
 	  WV(4)=8577.000d0
 	  WV(5)=8619.000d0
 	  WV(6)=8642.000d0
 	  WV(7)=8700.000d0
 	  WV(8)=8725.000d0
 	  WV(9)=8776.000d0
 	  WV(10)=8792.00d0
 	  WV(11)=8642.00d0
 	  WV(12)=8682.00d0
 	  FWV(1)=1.0d0
 	ELSEIF(NINDEX.EQ.39)THEN
 	  NAME='PaT'
 	  ITI=503
 	  WV(1)=8474.000d0
 	  WV(2)=8484.000d0
 	  WV(3)=8563.000d0
 	  WV(4)=8577.000d0
 	  WV(5)=8619.000d0
 	  WV(6)=8642.000d0
 	  WV(7)=8700.000d0
 	  WV(8)=8725.000d0
 	  WV(9)=8776.000d0
 	  WV(10)=8792.00d0
 	  WV(11)=8461.00d0
 	  WV(12)=8474.00d0
 	  FWV(1)=1.0d0
 	  WV(13)=8577.00d0
 	  WV(14)=8619.00d0
 	  FWV(2)=1.0d0
 	  WV(15)=8730.00d0
 	  WV(16)=8772.00d0
 	  FWV(3)=1.0d0
 	ELSEIF(NINDEX.EQ.40)THEN
 	  NAME='CaT*'
 	  ITI=506
 	  WV(1)=8474.000d0
 	  WV(2)=8484.000d0
 	  WV(3)=8563.000d0
 	  WV(4)=8577.000d0
 	  WV(5)=8619.000d0
 	  WV(6)=8642.000d0
 	  WV(7)=8700.000d0
 	  WV(8)=8725.000d0
 	  WV(9)=8776.000d0
 	  WV(10)=8792.00d0
 	  WV(11)=8484.00d0
 	  WV(12)=8513.00d0
 	  FWV(1)=1.0d0
 	  WV(13)=8522.00d0
 	  WV(14)=8562.00d0
 	  FWV(2)=1.0d0
 	  WV(15)=8642.00d0
 	  WV(16)=8682.00d0
 	  FWV(3)=1.0d0
 	  WV(17)=8461.00d0
 	  WV(18)=8474.00d0
 	  FWV(4)=-0.93d0
 	  WV(19)=8577.00d0
 	  WV(20)=8619.00d0
 	  FWV(5)=-0.93d0
 	  WV(21)=8730.00d0
 	  WV(22)=8772.00d0
 	  FWV(6)=-0.93d0
 	ELSEIF(NINDEX.EQ.41)THEN
 	  NAME='AZ1'
 	  ITI=2
 	  WV(1)=8474.000d0
 	  WV(2)=8489.000d0
 	  WV(3)=8490.000d0
 	  WV(4)=8506.000d0
 	  WV(5)=8521.000d0
 	  WV(6)=8531.000d0
 	ELSEIF(NINDEX.EQ.42)THEN
 	  NAME='AZ2'
 	  ITI=2
 	  WV(1)=8521.000d0
 	  WV(2)=8531.000d0
 	  WV(3)=8532.000d0
 	  WV(4)=8552.000d0
 	  WV(5)=8555.000d0
 	  WV(6)=8595.000d0
 	ELSEIF(NINDEX.EQ.43)THEN
 	  NAME='AZ3'
 	  ITI=2
 	  WV(1)=8626.000d0
 	  WV(2)=8650.000d0
 	  WV(3)=8653.000d0
 	  WV(4)=8671.000d0
 	  WV(5)=8695.000d0
 	  WV(6)=8725.000d0
 	ELSE
 	  WRITE(*,100) 'NINDEX='
 	  WRITE(*,*) NINDEX
 	  WRITE(*,101) 'FATAL ERROR: index number out of range.'
 	  WRITE(*,101) '(subroutine INDEXDEF)'
 	  STOP
 	END IF
C
100	FORMAT(A,$)
101	FORMAT(A)
 	END
