CFLAGS= -g
$CC= cc

OBJS =	imsetwcs.o imgetwcs.o matchstar.o findstar.o daoread.o wcscon.o \
	slasubs.o fitswcs.o wcs.o gscread.o ujcread.o uacread.o tabread.o \
	imrotate.o fitsio.o imhio.o hget.o hput.o imio.o worldpos.o \
	platepos.o sortstar.o cel.o proj.o sph.o wcstrig.o

libwcs.a:	$(OBJS)
	ar rv $@ $?
	ranlib $@

cel.o:		cel.h proj.h wcstrig.h
findstar.o:	fitshead.h lwcs.h
fitsio.o:	fitshead.h
fitswcs.o:	fitshead.h wcs.h
gscread.o:	fitshead.h wcs.h
hget.o:		fitshead.h
hput.o:		fitshead.h
imrotate.o:	fitshead.h
imgetwcs.o:	fitshead.h wcs.h lwcs.h
imsetwcs.o:	fitshead.h wcs.h lwcs.h
imhio.o:	fitshead.h
matchstar.o:	wcs.h lwcs.h
platepos.o:	wcs.h
proj.o:		proj.h wcstrig.h
sph.o:		wcstrig.h
tabread.o:	fitshead.h
uacread.o:	wcs.h
ujcread.o:	wcs.h
wcs.o:		fitshead.h wcs.h cel.h proj.h
wcscon.o:	wcs.h
wcstrig.o:	wcstrig.h
worldpos.o:	wcs.h
