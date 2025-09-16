       program ro
C=============================================================================C
C                                                                             C
C      program for the course "photogrammetry exercise":                      C
C              the topics: determination of relative orientation parameters.  C
C                                                                             C
C      program written by Jaan-Rong Tsay on 12-14 December 2001 in NCKU       C
C                                                                             C
C=============================================================================C

       character*80 name
       real*8       b,pl,kl,wr,pr,kr,f,xl,yl,xr,yr,hd(6),e,s0,vol
       real*4       a(21)
       real*8       spl,skl,swr,spr,skr
       integer*4    no,ir,ib,an,id1,nl(6),iop,lf,hl(6),ic,itr

       write(*,1)
 1     format(/' C',77('='),'C'/' C',77x,'C'/' C',5x,
     +   'program for the course "photogrammetry exercise":',23x,'C'/
     +   ' C',13x,'the topics: determination of relative orientation',
     +   ' parameters.   C'/' C',77x,'C'/' C',5x,'program written by ',
     +   'Jaan-Rong Tsay on 12 December 2001 in NCKU           C'/,
     +   ' C',' program adapted by Yi-Hsien,Tsai in October 2024',
     +   ' based on the original work C'/,
     +   ' C',77x,'C'/' C',77('='),'C'//
     +   ' Data format of the input file :'/
     +         6x,'Focal length (mm)'/
     +         6x,'Point number, xl(mm),yl(mm), xr(mm),yr(mm)'/
     +         6x,'     ...'/
     +         ' =================================='/
     +         ' > Input file name = ? (e.g. ro.dat)')
       read(*,'(a80)') name
       open(1,file=name,err=99)

       write(*,*) ' > Output file name = ? (e.g. ro.out)'
       read(*,'(a80)') name
       open(2,file=name,err=99)

       write(2,3)
 3     format('C',77('='),'C'/'C',77x,'C'/'C    PHGR_EX3 ',
     +  ': DETERMINATION OF RELATIVE ORIENTATION PARAMETERS          ',
     +  '    C'/'C',77x,'C'/'C             TO COMPUTE THE ',
     +  'R.O. PARAMETERS BY SWING-SWING METHOD            C'/
     +  'C',77x,'C'/'C',24x,'SCHEME-II:VOLUME DISCREPANCY',25x,'C'/
     +  'C',77x,'C'/'C',77('.'),'C'/'C',77x,'C'/'C       PROGRAM ',
     +  'WRITTEN BY Jaan-Rong Tsay on 12 December 2001 in NCKU',9x,'C'/
     +  'C',7x,'ADAPTED BY Yi-Hsien Tsai in October 2024 based on the',
     +  ' original work','   C'/
     +  'C',77x,'C'/'C',77('='),'C'//
     +  ' Relative Orientation Parameters:'//
     +  '         phie_L : phie  angle of the left  photo'/
     +  '         kapa_L : kapa  angle of the left  photo'/
     +  '         omega_R: omega angle of the right photo'/
     +  '         phie_R : phie  angle of the right photo'/
     +  '         kapa_R : kapa  angle of the right photo'/)

       write(*,*) ' > initial value of phie_L  = ? deg. (e.g. 0.0)'
       read(*,*) pl
       write(*,*) ' > initial value of kapa_L  = ? deg. (e.g. 0.0)'
       read(*,*) kl
       write(*,*) ' > initial value of omega_R = ? deg. (e.g. 0.0)'
       read(*,*) wr
       write(*,*) ' > initial value of phie_R  = ? deg. (e.g. 0.0)'
       read(*,*) pr
       write(*,*) ' > initial value of kapa_R  = ? deg. (e.g. 0.0)'
       read(*,*) kr
       write(2,16) pl,kl,wr,pr,kr
 16    format(' Initial values of 5 R.O. parameters (units :degree):'//
     +   6x,'Left  photo:                  phie=',f10.6,' kapa=',f10.6/
     +   6x,'Right photo: omega=',f10.6,' phie=',f10.6,' kapa=',f10.6)

C..................... change the angle unit to radian(s)
       call degrad(pl)
       call degrad(kl)
       call degrad(wr)
       call degrad(pr)
       call degrad(kr)

       write(*,26)
 26    format(/' > Threshold value adopted in convergence condition =',
     +         ' ? (e.g. 0.00000001)')
       read(*,*) e
       write(2,17) e
 17    format(//' Computation converges if max. |X| is less than',
     +          f20.17/)

       write(2,18)
       write(*,18)
 18    format('  ITER    d_phie(L)    d_kapa(L)    d_omega(R)   ',
     +        ' d_phie(R)    d_kapa(R)'/
     +        '  ==== ============ ============ ============= ',
     +        '============ ============')

C..... itr = number of iterations
       itr=0
 10    continue
       itr=itr+1
       do 4 i=1,21
 4     a(i)=0.e0

C............... f = focal length ; a(21) = augmented matrix [N|U]
C...............                    N = At A ; U = At L for V = A X - L
       read(1,*) f
 5     read(1,*,end=6,err=6) no,xl,yl,xr,yr
       call norsum(a,f,pl,kl,wr,pr,kr,xl,yl,xr,yr)
       goto 5
 6     continue

       ir=0
       ib=5
       an=21
       id1=6
       nl(1)=0
       nl(2)=1
       nl(3)=3
       nl(4)=6
       nl(5)=10
       nl(6)=15   
C............ iop = 1 : solve the system
C............ iop = 2 : solve the system and determine its inverse
       iop=1
C............ lf = 0 : no output of detailed messages
C............ lf = 2 : output detailed messages
       lf=0
CRT       write(*,7)
CRT 7     format(/' Figures concerned with solving the normal equations:')

       call spalg4(a,nl,hd,hl,id1,ir,ib,an,lf,1.d-15,iop,6)

       if(lf.lt.0) then
         write(*,*) ' *** N-matrix is singular ! ***'
         write(2,*) ' *** N-matrix is singular ! ***'
         stop
       endif

       pl=pl+a(16)
       kl=kl+a(17)
       wr=wr+a(18)
       pr=pr+a(19)
       kr=kr+a(20)

       write(2,19) itr,(a(15+i),i=1,5)
       write(*,19) itr,(a(15+i),i=1,5)
 19    format(i6,2f13.9,f14.9,2f13.9)

       j=0
       do 8 i=16,20
         if(abs(a(i)).gt.e) then
           j=1
           goto 9
         endif
 8     continue
 9     if(j.eq.1) then
         rewind(1)
         goto 10
       endif

       iop=2
       lf=0
CRT       write(*,11)
CRT 11    format(/' Figures concerned with computing the inverse:')

       call spalg4(a,nl,hd,hl,id1,ir,ib,an,lf,1.d-15,iop,6)

CRT       write(*,12)
CRT 12    format(/' Solution and inversion operations are completed.')

       write(2,20)
 20    format(//' vol = pseudo-observation of a volume defined by the',
     + ' 3 vectors'/7x,'of baseline vector and two homolog ray-ones'/
     + ' F   = vl*wr - vr*wl where (ul,vl,wl) and (ur,vr,wr) are ',
     + 'defined in the lecture'/' vol = F * B where B = baseline length'
     + //7x,'NO     xl(mm)     yl(mm)     xr(mm)     yr(mm)',11x,
     + 'F(mm^2)'/' ======== ========== ========== ========== ==========',
     + ' =================')

       ic=0
       s0=0.d0
       rewind(1)
       read(1,*) f
 13    read(1,*,end=14,err=14) no,xl,yl,xr,yr
       call volum(vol,f,pl,kl,wr,pr,kr,xl,yl,xr,yr)
       write(2,21) no,xl,yl,xr,yr,vol
 21    format(i9,4f11.6,f18.12)
       s0=s0+vol*vol
       ic=ic+1
       goto 13

 14    write(*,*) ' > A reference baseline B = ? mm. (e.g. 40)'
       read(*,*) b
       if(ic.ge.1) then
         write(2,15) b*dsqrt(s0/ic)
 15      format(/' RMS value of all pseudo-observations of "volumes" =',
     +           f15.6,' mm^3')
       endif
       if(ic.gt.5) then
         write(2,22) dsqrt(s0/(ic-5))
 22      format(/' standard deviation of unit weight =',f16.8,' mm^2')
       endif
       if(ic.ge.5) then
         write(2,23) ic-5
 23      format(/' Degree of freedom =',i4/)
       endif

C.......... change the angle unit to degree
       call radeg(pl)
       call radeg(kl)
       call radeg(wr)
       call radeg(pr)
       call radeg(kr)

       if(ic.gt.5) then
         s0=dsqrt(s0/(ic-5))
         spl=s0*sqrt(a(1))
         skl=s0*sqrt(a(3))
         swr=s0*sqrt(a(6))
         spr=s0*sqrt(a(10))
         skr=s0*sqrt(a(15))
         call radeg(spl)
         call radeg(skl)
         call radeg(swr)
         call radeg(spr)
         call radeg(skr)
         write(2,24) pl,spl,kl,skl,wr,swr,pr,spr,kr,skr
 24      format(' Phie_L  =',f12.6,' +/-',f10.6,' deg.'/
     +          ' Kapa_L  =',f12.6,' +/-',f10.6,' deg.'/
     +          ' Omega_R =',f12.6,' +/-',f10.6,' deg.'/
     +          ' Phie_R  =',f12.6,' +/-',f10.6,' deg.'/
     +          ' Kapa_R  =',f12.6,' +/-',f10.6,' deg.')
       elseif(ic.eq.5) then
         write(2,25) pl,kl,wr,pr,kr
 25      format(' Phie_L  =',f12.6,' deg.'/
     +          ' Kapa_L  =',f12.6,' deg.'/
     +          ' Omega_R =',f12.6,' deg.'/
     +          ' Phie_R  =',f12.6,' deg.'/
     +          ' Kapa_R  =',f12.6,' deg.')
       endif

C................. compute the model coordinates
       rewind(1)
       read(1,*) f
       write(2,29)
 29    format(//' (U,V,W) = model coordinates'//9x,'NO',
     +          '        U(mm)        V(mm)        W(mm)'/
     +          ' ========== ============ ============ ============')

C...................... change the angle unit to radian(s)
       call degrad(pl)
       call degrad(kl)
       call degrad(wr)
       call degrad(pr)
       call degrad(kr)

 27    read(1,*,end=28,err=28) no,xl,yl,xr,yr
       call mc(no,f,pl,kl,wr,pr,kr,xl,yl,xr,yr,b)
       goto 27
 28    continue

       goto 1000
 99    write(*,*) ' *** invalid file-opening !!!'
       stop

 1000  write(2,2)
 2     format(//'C',77('='),'C'/'C',77x,'C'/'C',35x,
     +        'THE END'35x,'C'/'C',77x,'C'/'C',77('='),'C')
       end

C................ mc : compute the model coordinates
       subroutine mc(no,f,pl,kl,wr,pr,kr,xl,yl,xr,yr,b)
       real*8     f,wl,pl,kl,wr,pr,kr,xl,yl,xr,yr,b,s1,s2,u,v,w
       real*8     cwl,swl,cpl,spl,ckl,skl,cwr,swr,cpr,spr,ckr,skr
       real*8     u1,v1,w1,u2,v2,w2,d
       integer*4  no

       wl=0.d0
       cwl=cos(wl)
       swl=sin(wl)
       cpl=cos(pl)
       spl=sin(pl)
       ckl=cos(kl)
       skl=sin(kl)
       u1=xl*(cpl*ckl-spl*swl*skl)+yl*(-cpl*skl-spl*swl*ckl)+f*spl*cwl
       v1=xl*cwl*skl+yl*cwl*ckl+f*swl
       w1=xl*(spl*ckl+cpl*swl*skl)+yl*(-spl*skl+cpl*swl*ckl)-f*cpl*cwl

       cwr=cos(wr)
       swr=sin(wr)
       cpr=cos(pr)
       spr=sin(pr)
       ckr=cos(kr)
       skr=sin(kr)
       u2=xr*(cpr*ckr-spr*swr*skr)+yr*(-cpr*skr-spr*swr*ckr)+f*spr*cwr
       v2=xr*cwr*skr+yr*cwr*ckr+f*swr
       w2=xr*(spr*ckr+cpr*swr*skr)+yr*(-spr*skr+cpr*swr*ckr)-f*cpr*cwr

C...... compute the scale estimate for the right photo by averaging since
C......        there still exists residuals on photo coordinates
       d=u1*v2-v1*u2
       if(d.eq.0.d0) then
         write(*,*) ' ***** d1 = 0 ! [subroutine: mc]'
         stop
       endif
       s2=v1*b/d
       d=u1*w2-u2*w1
       if(d.eq.0.d0) then
         write(*,*) ' ***** d2 = 0 ! [subroutine: mc]'
         stop
       endif
       s2=(s2+w1*b/d)/2.d0

C...... output the model coordinates which are average values since there
C......        still exists residuals on photo coordinates
       u=b+u2*s2
       v=v2*s2
       w=w2*s2

       if(u1.eq.0.d0) then
         write(*,*) ' ***** u1 = 0 ! [subroutine: mc]'
         return
       endif
       s1=(b+u2*s2)/u1
       write(2,1) no,(u+s1*u1)/2.d0,(v+s1*v1)/2.d0,(w+s1*w1)/2.d0
 1     format(i11,3f13.3)

       return
       end

C............... volum : compute the pseudo-observation of a F-value
C...............         F = vol/B = vl*wr - vr*wl
      subroutine volum(vol,f,pl,kl,wr,pr,kr,xl,yl,xr,yr)
      real*8     vol,f,wl,pl,kl,wr,pr,kr,xl,yl,xr,yr

      wl=0.d0
      cwl=cos(wl)
      swl=sin(wl)
      cpl=cos(pl)
      spl=sin(pl)
      ckl=cos(kl)
      skl=sin(kl)
      u1=xl*(cpl*ckl-spl*swl*skl)+yl*(-cpl*skl-spl*swl*ckl)+f*spl*cwl
      v1=xl*cwl*skl+yl*cwl*ckl+f*swl
      w1=xl*(spl*ckl+cpl*swl*skl)+yl*(-spl*skl+cpl*swl*ckl)-f*cpl*cwl

      cwr=cos(wr)
      swr=sin(wr)
      cpr=cos(pr)
      spr=sin(pr)
      ckr=cos(kr)
      skr=sin(kr)
      u2=xr*(cpr*ckr-spr*swr*skr)+yr*(-cpr*skr-spr*swr*ckr)+f*spr*cwr
      v2=xr*cwr*skr+yr*cwr*ckr+f*swr
      w2=xr*(spr*ckr+cpr*swr*skr)+yr*(-spr*skr+cpr*swr*ckr)-f*cpr*cwr

      vol=v1*w2-v2*w1

      return
      end

C............... norsum : build the normal system by summing all elements
      subroutine norsum(a,f,pl,kl,wr,pr,kr,xl,yl,xr,yr)
      real*8 f,wl,pl,kl,wr,pr,kr,xl,yl,xr,yr
      real*8 u1,v1,w1,u2,v2,w2,b(6)
      real*8 cpl,spl,cwl,swl,ckl,skl,cpr,spr,cwr,swr,ckr,skr
      real*4 a(21)
      integer*4 i,j,ij

      wl=0.d0
      cwl=cos(wl)
      swl=sin(wl)
      cpl=cos(pl)
      spl=sin(pl)
      ckl=cos(kl)
      skl=sin(kl)
      u1=xl*(cpl*ckl-spl*swl*skl)+yl*(-cpl*skl-spl*swl*ckl)+f*spl*cwl
      v1=xl*cwl*skl+yl*cwl*ckl+f*swl
      w1=xl*(spl*ckl+cpl*swl*skl)+yl*(-spl*skl+cpl*swl*ckl)-f*cpl*cwl

      cwr=cos(wr)
      swr=sin(wr)
      cpr=cos(pr)
      spr=sin(pr)
      ckr=cos(kr)
      skr=sin(kr)
      u2=xr*(cpr*ckr-spr*swr*skr)+yr*(-cpr*skr-spr*swr*ckr)+f*spr*cwr
      v2=xr*cwr*skr+yr*cwr*ckr+f*swr
      w2=xr*(spr*ckr+cpr*swr*skr)+yr*(-spr*skr+cpr*swr*ckr)-f*cpr*cwr

C....   1:phie_L   2:kapa_L   3:omega_R   4:phie_R   5:kapa_R       l
C....     b(1)        b(2)       b(3)        b(4)      b(5)        b(6)
C....      V = B X - L  (observation equations)
      b(1)=u1*v2
      b(2)=-w2*(u1*cpl*cwl+w1*spl*cwl)+v2*(u1*swl-v1*spl*cwl)
      b(3)=-v1*v2*cpr+w1*(u2*spr-w2*cpr)
      b(4)=-u2*v1
      b(5)=-v1*(u2*swr-v2*spr*cwr)+w1*(u2*cpr*cwr+w2*spr*cwr)
      b(6)=v1*w2-v2*w1

C....     a(21) = [ N | U ] with N=At*A & U=At*L (normal equations)
      do 1 i=1,5
        a(15+i)=a(15+i)+b(i)*b(6)
        do 2 j=1,i
          ij=(i-1)*i/2+j
          a(ij)=a(ij)+b(i)*b(j)
 2      continue
 1    continue

      return
      end

C............... degrad = transform a degree unit to a radian unit
      subroutine degrad(a)
      real*8 a

      a=a*3.141592654/180.d0

      return
      end

C............... radeg = transform a radian unit to a degree unit
      subroutine radeg(a)
      real*8 a

      a=a*180.d0/3.141592654

      return
      end

      SUBROUTINE SPALG4 (A,N,HD,L,ID1,IR,IP,NX,LF,E,IOP,KOUT)
C- jobfv4.f:call spalg4 (vekngl,nl,hd,hl,id1,ir,ib,an,lf,epsi,iop,outpv)
C---- IMPLICIT  LOGICAL(A-Z)
      INTEGER*4 NX,ID1
      REAL*4    A(NX)
      REAL*8    HD(ID1),AII,R,HDI,DIF,AAA,HDJ
      INTEGER*4 N(ID1),L(ID1),IR,IP,LF,IOP
      REAL*8    E,E0,EPS
      integer*4 KOUT
      INTEGER*4 IA,IE,ID,I1,IS,NULLV,NULLN,NN11,NNUL
      INTEGER*4 I,J,K,I0,K0,II,IJ,KJ,MJ,IND,JND,KND,KKK,III
      INTEGER*4 IED,NA,NU,IED0,IEDM1,JA,JE,JJJ,KI,JND1,KK

C .......INVERSION EINER HERMITESCHEN PROFILMATRIX NG ..................
C  A     NG-MATRIX ALS VEKTOR                                          .
C        SPEICHERUNG DES UNTEREN (oberen) DREIECKS IM PROFIL,          .
C        EINSCHLIESSLICH RECHTER SEITEN UND PLL/PVV                    .
C  N     INDEXFELD FUER 0-TE SPALTE (Zeile)                            .
C  HD,L  HILFSVEKTOREN                                                 .
C  ID1   ANZAHL DER UNBEKANNTEN +1 (au+1)                              .
C  IR    BREITE DES RANDES (ir=az+da*(ba-1))                           .
C  IP    BREITE DES BANDES OHNE RAND (ib=s*t+3=l+3, ohne G-Krue.)      .
C  NX    LAENGE DES NG-VEKTORS (an)                                    .
C  LF    STEUERUNG DES AUSDRUCKS                                       .
C        0  KEIN AUSDRUCK                                              .
C        2  AUSFUEHRLICHER AUSDRUCK                                    .
C        NACH VERLASSEN DES UNTERPROGRAMMS                             .
C        1  O.K.                                                       .
C        -1 NG-MATRIX SINGULAER                                        .
C  E     ABBRUCHSCHRANKE                                               .
C  IOP   STEUERUNG DER BERECHNUNG                                      .
C        0  LOESUNG UND INVERSION                                      .
C        1  LOESUNG UND REDUKTION                                      .
C        2  INVERSION (VORR. REDUZIERTE MATRIX)                        .
C ......................................................................

C------------------------------------------------- NULLELEMENTE VOR REDUKTION
       NULLV=0
C------------------------------------------------ NULLELEMENTE NACH REDUKTION
       NULLN=0
C----------------------------------------------------- ANZAHL DER UNBEKANNTEN
       ID=ID1-1
C-------------------------------------- 1 + ANZAHL DER UNBEKANNTEN MINUS RAND
       IS=ID-IR+1
C------------------------------------------------------------------ GRENZWERT
       EPS=1.D-8
C----------------------------------------------------------------------- NULL
       E0=0.E0
C----------------------------------------------------- LAENGE DES NGL VEKTORS
       NN11=N(ID1)+ID1

C------------------------------------------------- NULLELEMENTE VOR REDUKTION
       DO 4 I=1,NN11
         IF( ABS( A(I) ).NE.E0) NULLV=NULLV+1
 4     CONTINUE

C------------------------------------------------ (IOP.EQ.2) => NUR INVERSION
       IF(IOP.EQ.2) GOTO 51

       DO 10 I=1,ID
C---------------------------------- NULLINDEX + ZEILENINDEX => HAUPTDIAGONALE
C******************************* volle NG-Matrix der indirekten M.(ohne N_GG)
C*        II=N(I)+I
C*   10   HD(I)=A(II)
         IF(I.GT.IP.AND.I.LT.IS) THEN
           II=N(I)+IP
         ELSE
           II=N(I)+I
         ENDIF
 10    HD(I)=A(II)

c     write(*,'(a)') 'Hauptdiagonale ...'
c     write(*,'(6f12.6)') (HD(I),I=1,ID)
c     pause

      DO 50 I=1,ID
C----------------------------------------------- NULLTER INDEX AKTUELLE ZEILE
        I0=N(I)
C---------------------------------------------------- NULLINDEX + ZEILENINDEX
C********************      II=I0+I
        IF(I.GT.IP.AND.I.LT.IS) THEN
          II=I0+IP
        ELSE
          II=I0+I
        ENDIF
C------------------------------------------------------- HAUPTDIAGONALELEMENT
        HDI=HD(I)
C------------------------------------------- REDUZIERTES HAUPTDIAGONALELEMENT
        AII=A(II)
C-------------------------------------------------------------------- AUSGABE
        DIF=HDI-AII
        IF(LF.EQ.0) GOTO 12
C--------- LF.EQ.0 : KEIN AUSDRUCK
        WRITE(KOUT,6000) I,HDI,AII,DIF
 6000   FORMAT(1X,I5,6E15.8)
        IF(ABS(HDI).LT.EPS) GOTO 12
        DIF=AII*100./HDI
        WRITE(KOUT,6010) DIF
 6010   FORMAT(1X,50X,2F10.2)
C-------------------------------------------------- FEHLERAUSGANG NULLABFRAGE
 12     IF(ABS(AII).LE.ABS(HDI)*E.OR.ABS(AII).EQ.0.) GOTO 99
        AII=1./AII
        A(II)=AII
        IND=0
C------------------------------------------------------------- NAECHSTE ZEILE
        I1=I+1
C--------------------------------------------------- NAECHSTE ZEILE IM RAND ?
        IF(I1.GE.IS) GOTO 13
C----------------------------------------------- INDEX I1 LIEGT NICHT IM RAND
C********************************* IE=I+IP oder I1+IP
        IE=I+IP-1
C--------------------------------------- ENDINDEX = I1 + BANDBREITE IM RAND ?
        IF(IE.GE.IS) IE=IS-1
C...........................................................................
C           ANFANG BIS ENDE IM BAND
C.......... NULLINDEX
        NU=N(I1-1)
        II=0
        DO 14 J=I1,IE
C-------------------------------------------------------------------- ZAEHLER
          II=II+1
C--------------------------------------------------------------------- ANFANG
          NA=N(J)
C------------------------------------------------------------- ANFANG + ZEILE
C************************ IJ=NA+I
          IF(J.GT.IP.AND.J.LT.IS) THEN
            IJ=NA+I+IP-J
          ELSE
            IJ=NA+I
          ENDIF
C-------------------------------------------------------------- ANFANG - NULL
          MJ=NA-NU
          NU=NA
          IF(MJ.LT.II) GOTO 14
          IF( ABS( A(IJ) ).EQ.E0) GOTO 14
          IND=IND+1
          L(IND)=J
 14     CONTINUE
C...........................................................................
C       INDEX I1 LIEGT IM RAND

C----------------------------------------------------------- RAND GLEICH NULL ?
 13     IF(IR.EQ.0) GOTO 11
        IA=IS
        IF(I1.GE.IS) IA=I1
C----------------------------------------------------------------- AUSSERHALB
        IF(IA.GT.ID) GOTO 11
C .....................................................................
C       RAND

        DO 9 J=IA,ID
          IND=IND+1
 9      L(IND)=J

C .....................................................................
C       RAND = NULL  ODER ABGEARBEITET

 11     IND=IND+1
        L(IND)=ID1
        IED=IND
        IED0=IED

C ............................................................................

        DO 30 JND=1,IED
          J=L(JND)
C****************************** IJ=N(J)+I
          IF(J.GT.IP.AND.J.LT.IS) THEN
            IJ=N(J)+I+IP-J
          ELSE
            IJ=N(J)+I
          ENDIF
          R=A(IJ)*AII
          IED0=IED0+1

          DO 25 KKK=JND,IED
            KND=IED0-KKK
            K=L(KND)
            K0=N(K)
C********************************* KJ=K0+J
            IF(K.GT.IP.AND.K.LT.IS) THEN
              KJ=K0+J+IP-K
              KI=K0+I+IP-K
            ELSE
              KJ=K0+J
              KI=K0+I
            ENDIF
C********************************* 25 A(KJ)=A(KJ)-R*A(K0+I)
 25       A(KJ)=A(KJ)-R*A(KI)

 30     A(IJ)=R

C ............................................................................

 50   CONTINUE

C ............................................................................
CRT--------------------------------------- Inversion: Vorr. reduzierte Matrix
 51   CONTINUE
      I0=N(ID1)

C ............................................................................
      DO 100 III=2,ID
        I=ID1-III
        II=I0+I
        AII=A(II)
        I1=I+1
        IND=0
        IF(I1.GE.IS) GOTO 63
        IE=I+IP-1
        IF(IE.GE.IS) IE=IS-1
        II=0
C .............................................................................
        NU=N(I1-1)
        DO 64 J=I1,IE
          II=II+1
          NA=N(J)
C************************************** IJ=NA+I
          IF(J.GT.IP.AND.J.LT.IS) THEN
            IJ=NA+I+IP-J
          ELSE
            IJ=NA+I
          ENDIF
          AAA=A(IJ)
          HD(J)=AAA
          MJ=NA-NU
          NU=NA
          IF(MJ.LT.II) GOTO 64
          IF( ABS( AAA ).EQ.E0) GOTO 64
          IND=IND+1
          L(IND)=J
          IF(IOP.NE.1) A(IJ)=E0
 64     CONTINUE
C ............................................................................
 63     IF(IR.EQ.0) GOTO 61
        IA=IS
        IF(I1.GE.IS) IA=I1
        IF(IA.GT.ID) GOTO 61
C ............................................................................
        DO 62 J=IA,ID
          IND=IND+1
          IJ=N(J)+I
          HD(J)=A(IJ)
          IF(IOP.NE.1) A(IJ)=E0
 62     L(IND)=J
C...........................................................................
 61     IND=IND+1
        L(IND)=ID1
        NULLN=NULLN+IND
        IJ=N(ID1)+I
        HD(ID1)=A(IJ)
        IED=IND
        IEDM1=IED-1
        IF(IEDM1.LT.1) GOTO 100
C...........................................................................

        DO 80 JJJ=1,IEDM1
          JND=IED-JJJ
          J=L(JND)
          HDJ=HD(J)
          JND1=JND-1
          K0=N(J)
C------------------------------------------------------ INVERSION ODER JND>1
          IF(JND1.LE.0.OR.IOP.EQ.1) GOTO 72
C
          DO 70 KND=1,JND1
            K=L(KND)
C****************************************** KI=N(K)+I
            IF(K.GT.IP.AND.K.LT.IS) THEN
              KI=N(K)+I+IP-K
            ELSE
              KI=N(K)+I
            ENDIF
            IF(J.GT.IP.AND.J.LT.IS) THEN
              KK=K0+K+IP-J
            ELSE
              KK=K0+K
            ENDIF
C****************************************** 70 A(KI)=A(KI)-HDJ*A(K0+K)
 70       A(KI)=A(KI)-HDJ*A(KK)
C
 72       JA=JND
          JE=IED
          IF(IOP.EQ.1) JA=IED
          IF(IOP.EQ.2) JE=IEDM1
          IF(JA.GT.JE) GOTO 80
C
          DO 75 KND=JA,JE
            K=L(KND)
            K0=N(K)
C******************************************* KI=K0+I
            IF(K.GT.IP.AND.K.LT.IS) THEN
              KI=K0+I+IP-K
              KJ=K0+J+IP-K
            ELSE
              KI=K0+I
              KJ=K0+J
            ENDIF
C******************************************* 75 A(KI)=A(KI)-HDJ*A(K0+J)
 75       A(KI)=A(KI)-HDJ*A(KJ)
C
 80     CONTINUE
C ............................................................................

        IF(IOP.EQ.1) GOTO 100
C********************************************* II=N(I)+I
        IF(I.GT.IP.AND.I.LT.IS) THEN
          II=N(I)+IP
        ELSE
          II=N(I)+I
        ENDIF

C ............................................................................

        DO 90 JJJ=1,IEDM1
          JND=IED-JJJ
          J=L(JND)
C*********************************************** IJ=N(J)+I
          IF(J.GT.IP.AND.J.LT.IS) THEN
            IJ=N(J)+I+IP-J
          ELSE
            IJ=N(J)+I
          ENDIF
 90     A(II)=A(II)-HD(J)*A(IJ)

C .............................................................................

 100  CONTINUE

C ............................................................................
      NNUL=ID1*(ID1+1)/2
      NULLN=NULLN+ID1

CRT--------------------------------------------------- Ausfuehrliche Ausdruck
      IF(LF .EQ. 2) WRITE(KOUT,6990) NNUL,NN11,NULLV,NULLN
 6990 FORMAT(11X,'Total number ',T44,I10/
     +       11X,'Profile'      ,T44,I10/
     +       11X,'non-zeros before reduction', T44,I10/
     +       11X,'non-zeros after reduction',T44,I10//)

C-------------------------------------------------------------------- AUSGANG
      LF=1
      RETURN

C------------------------------------------------------------- FEHLERAUSGANG
 99   LF=-1
      RETURN

      END
