CFLAGS= -g
CC=cc
LIBWCS = libwcs/libwcs.a
LIBS = $(LIBWCS) -lm
.PRECIOUS: ${LIBWCS}

all:	delwcs edhead fixpix gethead i2f imcat imhead immatch imrot \
	imsize imstar imwcs scat sethead addpix getpix setpix sky2xy \
	keyhead skycoor subpix xy2sky

addpix: addpix.c $(LIBWCS) libwcs/fitshead.h
	$(CC) $(CFLAGS) -o addpix addpix.c $(LIBS)

delwcs: delwcs.c $(LIBWCS) libwcs/fitshead.h
	$(CC) $(CFLAGS) -o delwcs delwcs.c $(LIBS)

edhead: edhead.c $(LIBWCS) libwcs/fitshead.h
	$(CC) $(CFLAGS) -o edhead edhead.c $(LIBS)

fixpix: fixpix.c $(LIBWCS) libwcs/fitshead.h
	$(CC) $(CFLAGS) -o fixpix fixpix.c $(LIBS)

gethead: gethead.c $(LIBWCS) libwcs/fitshead.h
	$(CC) $(CFLAGS) -o gethead gethead.c $(LIBS)

getpix: getpix.c $(LIBWCS) libwcs/fitshead.h
	$(CC) $(CFLAGS) -o getpix getpix.c $(LIBS)

i2f: i2f.c $(LIBWCS) libwcs/fitshead.h
	$(CC) $(CFLAGS) -o i2f i2f.c $(LIBS)

imcat: imcat.c $(LIBWCS) libwcs/fitshead.h libwcs/wcs.h
	$(CC) $(CFLAGS) -o imcat imcat.c $(LIBS)

imhead: imhead.c $(LIBWCS) libwcs/fitshead.h libwcs/wcs.h
	$(CC) $(CFLAGS) -o imhead imhead.c $(LIBS)

imrot: imrot.c $(LIBWCS) libwcs/fitshead.h
	$(CC) $(CFLAGS) -o imrot imrot.c $(LIBS)

imsize: imsize.c $(LIBWCS) libwcs/fitshead.h libwcs/wcs.h
	$(CC) $(CFLAGS) -o imsize imsize.c $(LIBS)

imstack: imstack.c $(LIBWCS) libwcs/fitshead.h
	$(CC) $(CFLAGS) -o imstack imstack.c $(LIBS)

imstar: imstar.c $(LIBWCS) libwcs/fitshead.h libwcs/wcs.h libwcs/lwcs.h
	$(CC) $(CFLAGS) -o imstar imstar.c $(LIBS)

imwcs: imwcs.c $(LIBWCS) libwcs/fitshead.h libwcs/lwcs.h
	$(CC) $(CFLAGS) -o imwcs imwcs.c $(LIBS)

immatch: immatch.c $(LIBWCS) libwcs/fitshead.h libwcs/lwcs.h
	$(CC) $(CFLAGS) -o immatch immatch.c $(LIBS)

keyhead: keyhead.c $(LIBWCS) libwcs/fitshead.h libwcs/wcs.h
	$(CC) $(CFLAGS) -o keyhead keyhead.c $(LIBS)

scat: scat.c $(LIBWCS) libwcs/fitshead.h libwcs/wcs.h
	$(CC) $(CFLAGS) -o scat scat.c $(LIBS)

sethead: sethead.c $(LIBWCS) libwcs/fitshead.h
	$(CC) $(CFLAGS) -o sethead sethead.c $(LIBS)

setpix: setpix.c $(LIBWCS) libwcs/fitshead.h
	$(CC) $(CFLAGS) -o setpix setpix.c $(LIBS)

sky2xy: sky2xy.c $(LIBWCS) libwcs/wcs.h
	$(CC) $(CFLAGS) -o sky2xy sky2xy.c $(LIBS)

skycoor: skycoor.c $(LIBWCS) libwcs/wcs.h
	$(CC) $(CFLAGS) -o skycoor skycoor.c $(LIBS)

subpix: subpix.c $(LIBWCS) libwcs/fitshead.h
	$(CC) $(CFLAGS) -o subpix subpix.c $(LIBS)

wcshead: wcshead.c $(LIBWCS) libwcs/fitshead.h
	$(CC) $(CFLAGS) -o wcshead wcshead.c $(LIBS)

xy2sky: xy2sky.c $(LIBWCS) libwcs/wcs.h
	$(CC) $(CFLAGS) -o xy2sky xy2sky.c $(LIBS)

$(LIBWCS): libwcs/*.c libwcs/fitshead.h libwcs/wcs.h
	cd libwcs; make
