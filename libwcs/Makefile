CFLAGS= -g
CC= cc

OBJS =	imsetwcs.o imgetwcs.o matchstar.o findstar.o daoread.o wcscon.o \
	slasubs.o fitswcs.o wcsinit.o wcs.o ty2read.o webread.o tmcread.o \
	gscread.o gsc2read.o ujcread.o uacread.o tabread.o binread.o ctgread.o \
	actread.o catutil.o imrotate.o fitsfile.o imhfile.o hget.o hput.o \
	imio.o dateutil.o worldpos.o tnxpos.o dsspos.o platepos.o \
	sortstar.o platefit.o tnxfit.o iget.o fileutil.o polfit.o \
	wcslib.o lin.o cel.o proj.o sph.o wcstrig.o 

libwcs.a:	$(OBJS)
	ar rv $@ $?
	ranlib $@

actread.o:	fitsfile.h wcscat.h wcs.h
binread.o:	wcscat.h wcs.h
ctgread.o:	wcscat.h wcs.h
catutil.o:	wcscat.h wcs.h
cel.o:		wcslib.h
dateutil.o:	fitsfile.h wcscat.h
dsspos.o:	wcs.h
fileutil.o:	fitshead.h
findstar.o:	fitshead.h lwcs.h imio.h
fitsfile.o:	fitsfile.h
fitswcs.o:	fitsfile.h wcs.h
gscread.o:	fitsfile.h wcs.h wcscat.h
gsc2read.o:	fitsfile.h wcs.h wcscat.h
hget.o:		fitshead.h
hput.o:		fitshead.h
iget.o:		fitshead.h
imgetwcs.o:	fitshead.h wcs.h lwcs.h
imhfile.o:	fitsfile.h
imio.o:		fitshead.h imio.h
imsetwcs.o:	fitshead.h wcs.h lwcs.h wcscat.h
imrotate.o:	fitshead.h imio.h
lin.o:		wcslib.h
matchstar.o:	fitshead.h wcs.h lwcs.h
platefit.o:	wcs.h
platepos.o:	wcs.h
polfit.o:	
proj.o:		wcslib.h
sortstar.o:	wcscat.h
sph.o:		wcslib.h
tabread.o:	fitshead.h wcs.h
tmcread.o:	fitsfile.h wcs.h wcscat.h
tnxpos.o:	wcs.h
ty2read.o:	fitsfile.h wcs.h wcscat.h
uacread.o:	fitshead.h wcs.h wcscat.h
ujcread.o:	fitshead.h wcs.h wcscat.h
wcs.o:		fitshead.h wcs.h
wcsinit.o:	fitshead.h wcs.h
wcscon.o:	wcs.h
wcslib.o:	wcslib.h
wcstrig.o:	wcslib.h
webread.o:	wcscat.h
worldpos.o:	wcs.h
fitsfile.h:	fitshead.h imio.h
wcs.h:		wcslib.h fitshead.h
