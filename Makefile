CFLAGS= -g
CC=cc
LIBWCS = libwcs/libwcs.a
LIBS = $(LIBWCS) -lm
BIN = bin
.PRECIOUS: ${LIBWCS}

all:	delwcs delhead edhead fixpix gethead i2f imcat imhead immatch imrot \
	imsize imstar imwcs scat sethead addpix getpix setpix sky2xy \
	keyhead skycoor subpix xy2sky wcshead conpix

addpix: addpix.c $(LIBWCS) libwcs/fitshead.h
	$(CC) $(CFLAGS) -o $(BIN)/addpix addpix.c $(LIBS)

conpix: conpix.c $(LIBWCS) libwcs/fitshead.h
	$(CC) $(CFLAGS) -o $(BIN)/conpix conpix.c $(LIBS)

delwcs: delwcs.c $(LIBWCS) libwcs/fitshead.h
	$(CC) $(CFLAGS) -o $(BIN)/delwcs delwcs.c $(LIBS)

delhead: delhead.c $(LIBWCS) libwcs/fitshead.h
	$(CC) $(CFLAGS) -o $(BIN)/delhead delhead.c $(LIBS)

edhead: edhead.c $(LIBWCS) libwcs/fitshead.h
	$(CC) $(CFLAGS) -o $(BIN)/edhead edhead.c $(LIBS)

fixpix: fixpix.c $(LIBWCS) libwcs/fitshead.h
	$(CC) $(CFLAGS) -o $(BIN)/fixpix fixpix.c $(LIBS)

gethead: gethead.c $(LIBWCS) libwcs/fitshead.h
	$(CC) $(CFLAGS) -o $(BIN)/gethead gethead.c $(LIBS)

getpix: getpix.c $(LIBWCS) libwcs/fitshead.h
	$(CC) $(CFLAGS) -o $(BIN)/getpix getpix.c $(LIBS)

i2f: i2f.c $(LIBWCS) libwcs/fitshead.h
	$(CC) $(CFLAGS) -o $(BIN)/i2f i2f.c $(LIBS)

imcat: imcat.c $(LIBWCS) libwcs/fitshead.h libwcs/wcs.h
	$(CC) $(CFLAGS) -o $(BIN)/imcat imcat.c $(LIBS)

imhead: imhead.c $(LIBWCS) libwcs/fitshead.h libwcs/wcs.h
	$(CC) $(CFLAGS) -o $(BIN)/imhead imhead.c $(LIBS)

imrot: imrot.c $(LIBWCS) libwcs/fitshead.h
	$(CC) $(CFLAGS) -o $(BIN)/imrot imrot.c $(LIBS)

imsize: imsize.c $(LIBWCS) libwcs/fitshead.h libwcs/wcs.h
	$(CC) $(CFLAGS) -o $(BIN)/imsize imsize.c $(LIBS)

imstack: imstack.c $(LIBWCS) libwcs/fitshead.h
	$(CC) $(CFLAGS) -o $(BIN)/imstack imstack.c $(LIBS)

imstar: imstar.c $(LIBWCS) libwcs/fitshead.h libwcs/wcs.h libwcs/lwcs.h
	$(CC) $(CFLAGS) -o $(BIN)/imstar imstar.c $(LIBS)

imwcs: imwcs.c $(LIBWCS) libwcs/fitshead.h libwcs/lwcs.h
	$(CC) $(CFLAGS) -o $(BIN)/imwcs imwcs.c $(LIBS)

immatch: immatch.c $(LIBWCS) libwcs/fitshead.h libwcs/lwcs.h
	$(CC) $(CFLAGS) -o $(BIN)/immatch immatch.c $(LIBS)

keyhead: keyhead.c $(LIBWCS) libwcs/fitshead.h libwcs/wcs.h
	$(CC) $(CFLAGS) -o $(BIN)/keyhead keyhead.c $(LIBS)

scat: scat.c $(LIBWCS) libwcs/fitshead.h libwcs/wcs.h
	$(CC) $(CFLAGS) -o $(BIN)/scat scat.c $(LIBS)

sethead: sethead.c $(LIBWCS) libwcs/fitshead.h
	$(CC) $(CFLAGS) -o $(BIN)/sethead sethead.c $(LIBS)

setpix: setpix.c $(LIBWCS) libwcs/fitshead.h
	$(CC) $(CFLAGS) -o $(BIN)/setpix setpix.c $(LIBS)

sky2xy: sky2xy.c $(LIBWCS) libwcs/wcs.h
	$(CC) $(CFLAGS) -o $(BIN)/sky2xy sky2xy.c $(LIBS)

skycoor: skycoor.c $(LIBWCS) libwcs/wcs.h
	$(CC) $(CFLAGS) -o $(BIN)/skycoor skycoor.c $(LIBS)

subpix: subpix.c $(LIBWCS) libwcs/fitshead.h
	$(CC) $(CFLAGS) -o $(BIN)/subpix subpix.c $(LIBS)

wcshead: wcshead.c $(LIBWCS) libwcs/fitshead.h
	$(CC) $(CFLAGS) -o $(BIN)/wcshead wcshead.c $(LIBS)

xy2sky: xy2sky.c $(LIBWCS) libwcs/wcs.h
	$(CC) $(CFLAGS) -o $(BIN)/xy2sky xy2sky.c $(LIBS)

$(LIBWCS): libwcs/*.c libwcs/fitshead.h libwcs/wcs.h
	cd libwcs; make
