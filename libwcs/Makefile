CFLAGS= -g
$CC= cc

OBJS =	imsetwcs.o imgetwcs.o matchstar.o findstar.o daoread.o wcscon.o \
	slasubs.o fitswcs.o wcs.o gscread.o ujcread.o uacread.o tabread.o \
	imrotate.o fitsio.o imhio.o hget.o hput.o imio.o worldpos.o tnxpos.o \
	dsspos.o platepos.o sortstar.o platefit.o iget.o \
	wcslib.o lin.o cel.o proj.o sph.o wcstrig.o 

libwcs.a:	$(OBJS)
	ar rv $@ $?
	ranlib $@

cel.o:		cel.h proj.h wcstrig.h
dsspos.o:	wcs.h
findstar.o:	fitshead.h lwcs.h
fitsio.o:	fitshead.h
fitswcs.o:	fitshead.h wcs.h
gscread.o:	fitshead.h wcs.h
hget.o:		fitshead.h
hput.o:		fitshead.h
iget.o:		fitshead.h
imgetwcs.o:	fitshead.h wcs.h lwcs.h
imhio.o:	fitshead.h
imio.o:		fitshead.h
imsetwcs.o:	fitshead.h wcs.h lwcs.h
imrotate.o:	fitshead.h
lin.o:		lin.h
matchstar.o:	fitshead.h wcs.h lwcs.h
platefit.o:	wcs.h
platepos.o:	wcs.h
proj.o:		proj.h wcstrig.h
sph.o:		wcstrig.h
tabread.o:	fitshead.h wcs.h
tnxpos.o:	wcs.h
uacread.o:	fitshead.h wcs.h
ujcread.o:	fitshead.h wcs.h
wcs.o:		fitshead.h wcs.h cel.h proj.h
wcscon.o:	wcs.h
wcslib.o:	wcslib.h
wcstrig.o:	wcstrig.h
worldpos.o:	wcs.h
wcs.h:		wcslib.h
wcslib.h:	cel.h lin.h
cel.h:		proj.h
