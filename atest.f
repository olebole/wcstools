	Real*8 CD(4), RCON
	Integer*4 I, J

	RCON = 1.74532925199433d-2

c  ngc1275-hst.fits
c	CD(1) = -5.354878e-05
c	CD(2) =  1.3992766e-05
c	CD(3) = -1.4005416e-05
c	CD(4) = -5.35004e-05
c  hsttest2.fits
c	CD(1) = -3.726666e-05
c	CD(2) =  4.052544e-05
c	CD(3) = -4.052544e-05
c	CD(4) = -3.726666e-05
c  m51hst.fits
	CD(1) = 6.37990400000001E-6
	CD(2) = 2.36777600000000E-5
	CD(3) = 2.36777600000000E-5
	CD(4) = -6.3799040000000E-6
	Write(*,*) CD(1)*CD(1) - CD(2)*CD(2)
	Write(*,*) CD(1)*CD(4) - CD(2)*CD(3)
	Do I = 1, 4
	    Do J = 1, 4
		If (I .ne. J) Then
		    Write(*,*) I, J, Datan2 (CD(I), CD(J)) / RCON
		    Endif
		Enddo
	    Enddo
	Do I = 1, 4
	    Do J = 1, 4
		If (I .ne. J) Then
		    Write(*,*) -I, J, Datan2 (-CD(I), CD(J)) / RCON
		    Endif
		Enddo
	    Enddo
	Do I = 1, 4
	    Do J = 1, 4
		If (I .ne. J) Then
		    Write(*,*) I, -J, Datan2 (CD(I), -CD(J)) / RCON
		    Endif
		Enddo
	    Enddo
	Stop
	End
