CFLAGS= -g
CC=cc
LIBWCS = libwcs/libwcs.a
LIBS = $(LIBWCS) -lm
#CATLIBS = $(LIBS) -lnsl -lsocket
CATLIBS = $(LIBS)
BIN = bin
.PRECIOUS: ${LIBWCS}
.c.o:
	$(CC) -c $(CFLAGS) $(DEFS) $<

all:	cphead delwcs delhead edhead fixpix gethead i2f imcat imhead immatch \
	imrot imsize imstar imwcs scat sethead addpix getpix setpix sky2xy \
	keyhead skycoor subpix xy2sky wcshead conpix gettab newfits \
	imstack imextract sumpix remap getcol getdate

addpix: addpix.c $(LIBWCS) libwcs/fitsfile.h
	$(CC) $(CFLAGS) -o $(BIN)/addpix addpix.c $(LIBS)

conpix: conpix.c $(LIBWCS) libwcs/fitsfile.h
	$(CC) $(CFLAGS) -o $(BIN)/conpix conpix.c $(LIBS)

cphead: cphead.c $(LIBWCS) libwcs/fitsfile.h
	$(CC) $(CFLAGS) -o $(BIN)/cphead cphead.c $(LIBS)

delwcs: delwcs.c $(LIBWCS) libwcs/fitsfile.h
	$(CC) $(CFLAGS) -o $(BIN)/delwcs delwcs.c $(LIBS)

delhead: delhead.c $(LIBWCS) libwcs/fitsfile.h
	$(CC) $(CFLAGS) -o $(BIN)/delhead delhead.c $(LIBS)

edhead: edhead.c $(LIBWCS) libwcs/fitsfile.h
	$(CC) $(CFLAGS) -o $(BIN)/edhead edhead.c $(LIBS)

fixpix: fixpix.c $(LIBWCS) libwcs/fitsfile.h
	$(CC) $(CFLAGS) -o $(BIN)/fixpix fixpix.c $(LIBS)

getcol: getcol.c $(LIBWCS) libwcs/wcscat.h
	$(CC) $(CFLAGS) -o $(BIN)/getcol getcol.c $(CATLIBS)

getdate: getdate.c $(LIBWCS) libwcs/wcscat.h
	$(CC) $(CFLAGS) -o $(BIN)/getdate getdate.c $(CATLIBS)

gethead: gethead.c $(LIBWCS) libwcs/fitsfile.h libwcs/wcscat.h
	$(CC) $(CFLAGS) -o $(BIN)/gethead gethead.c $(CATLIBS)

getpix: getpix.c $(LIBWCS) libwcs/fitsfile.h libwcs/wcscat.h
	$(CC) $(CFLAGS) -o $(BIN)/getpix getpix.c $(CATLIBS)

gettab: gettab.c $(LIBWCS) libwcs/wcscat.h
	$(CC) $(CFLAGS) -o $(BIN)/gettab gettab.c $(CATLIBS)

i2f: i2f.c $(LIBWCS) libwcs/fitsfile.h
	$(CC) $(CFLAGS) -o $(BIN)/i2f i2f.c $(LIBS)

imcat: imcat.c $(LIBWCS) libwcs/fitsfile.h libwcs/wcs.h
	$(CC) $(CFLAGS) -o $(BIN)/imcat imcat.c $(CATLIBS)

imhead: imhead.c $(LIBWCS) libwcs/fitsfile.h libwcs/wcs.h
	$(CC) $(CFLAGS) -o $(BIN)/imhead imhead.c $(LIBS)

imrot: imrot.c $(LIBWCS) libwcs/fitsfile.h
	$(CC) $(CFLAGS) -o $(BIN)/imrot imrot.c $(LIBS)

imsize: imsize.c $(LIBWCS) libwcs/fitsfile.h libwcs/wcs.h
	$(CC) $(CFLAGS) -o $(BIN)/imsize imsize.c $(LIBS)

imstack: imstack.c $(LIBWCS) libwcs/fitsfile.h
	$(CC) $(CFLAGS) -o $(BIN)/imstack imstack.c $(LIBS)

imextract: imextract.c $(LIBWCS) libwcs/fitsfile.h libwcs/wcscat.h
	$(CC) $(CFLAGS) -o $(BIN)/imextract imextract.c $(CATLIBS)

imstar: imstar.c $(LIBWCS) libwcs/fitsfile.h libwcs/wcs.h libwcs/lwcs.h libwcs/wcscat.h
	$(CC) $(CFLAGS) -o $(BIN)/imstar imstar.c $(CATLIBS)

imwcs: imwcs.c $(LIBWCS) libwcs/fitsfile.h libwcs/lwcs.h
	$(CC) $(CFLAGS) -o $(BIN)/imwcs imwcs.c $(CATLIBS)

immatch: immatch.c $(LIBWCS) libwcs/fitsfile.h libwcs/lwcs.h
	$(CC) $(CFLAGS) -o $(BIN)/immatch immatch.c $(CATLIBS)

immwcs: immwcs.c $(LIBWCS) libwcs/fitsfile.h libwcs/lwcs.h
	$(CC) $(CFLAGS) -o $(BIN)/immwcs immwcs.c $(CATLIBS)

keyhead: keyhead.c $(LIBWCS) libwcs/fitsfile.h libwcs/wcs.h
	$(CC) $(CFLAGS) -o $(BIN)/keyhead keyhead.c $(LIBS)

newfits: newfits.c $(LIBWCS) libwcs/fitshead.h libwcs/wcs.h
	$(CC) $(CFLAGS) -o $(BIN)/newfits newfits.c $(LIBS)

remap: remap.c $(LIBWCS) libwcs/fitsfile.h libwcs/wcs.h
	$(CC) $(CFLAGS) -o $(BIN)/remap remap.c $(LIBS)

scat: scat.c $(LIBWCS) libwcs/wcscat.h libwcs/wcs.h
	$(CC) $(CFLAGS) -o $(BIN)/scat scat.c $(CATLIBS)

sethead: sethead.c $(LIBWCS) libwcs/fitsfile.h
	$(CC) $(CFLAGS) -o $(BIN)/sethead sethead.c $(LIBS)

setpix: setpix.c $(LIBWCS) libwcs/fitsfile.h libwcs/wcscat.h
	$(CC) $(CFLAGS) -o $(BIN)/setpix setpix.c $(CATLIBS)

sky2xy: sky2xy.c $(LIBWCS) libwcs/wcs.h libwcs/fitsfile.h
	$(CC) $(CFLAGS) -o $(BIN)/sky2xy sky2xy.c $(LIBS)

skycoor: skycoor.c $(LIBWCS) libwcs/wcs.h libwcs/wcscat.h
	$(CC) $(CFLAGS) -o $(BIN)/skycoor skycoor.c $(CATLIBS)

subpix: subpix.c $(LIBWCS) libwcs/fitsfile.h
	$(CC) $(CFLAGS) -o $(BIN)/subpix subpix.c $(LIBS)

sumpix: sumpix.c $(LIBWCS) libwcs/fitsfile.h libwcs/wcscat.h
	$(CC) $(CFLAGS) -o $(BIN)/sumpix sumpix.c $(CATLIBS)

wcshead: wcshead.c $(LIBWCS) libwcs/fitsfile.h
	$(CC) $(CFLAGS) -o $(BIN)/wcshead wcshead.c $(LIBS)

xy2sky: xy2sky.c $(LIBWCS) libwcs/wcs.h
	$(CC) $(CFLAGS) -o $(BIN)/xy2sky xy2sky.c $(LIBS)

$(LIBWCS): libwcs/*.c libwcs/fitshead.h libwcs/wcs.h
	cd libwcs; make
