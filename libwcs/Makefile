CFLAGS= -g
$CC= cc

OBJS =	imsetwcs.o matchstar.o findstar.o wcscon.o \
	fitswcs.o wcs.o gscread.o ujcread.o tabread.o imrotate.o \
	fitsio.o irafio.o hget.o hput.o imio.o worldpos.o platepos.o sortstar.o

libwcs.a:	$(OBJS)
	ar rv $@ $?
	ranlib $@

findstar.o:	fitshead.h
fitsio.o:	fitshead.h
fitswcs.o:	fitshead.h wcs.h
gscread.o:	fitshead.h
hget.o:		fitshead.h
hput.o:		fitshead.h
imrotate.o:	fitshead.h
imsetwcs.o:	fitshead.h wcs.h
irafio.o:	fitshead.h
matchstar.o:	wcs.h
platepos.o:	wcs.h
tabread.o:	fitshead.h
wcs.o:		fitshead.h wcs.h
wcscon.o:	wcs.h
worldpos.o:	wcs.h
