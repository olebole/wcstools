CFLAGS= -g

OBJS =	imsetwcs.o findreg.o findcoor.o findstar.o wcscon.o getwcs.o \
	wcs.o gscread.o fitsio.o irafio.o hget.o hput.o worldpos.o platepos.o

libwcs.a:	$(OBJS)
	ar rv $@ $?
	ranlib $@

findcoor.o:	wcs.h
findreg.o:	wcs.h
findstar.o:	fitshead.h
fitsio.o:	fitshead.h
getwcs.o:	fitshead.h wcs.h
gscread.o:	fitshead.h
hget.o:		fitshead.h
hput.o:		fitshead.h
imsetwcs.o:	fitshead.h wcs.h
irafio.o:	fitshead.h
platepos.o:	wcs.h
wcs.o:		fitshead.h wcs.h
worldpos.o:	wcs.h
