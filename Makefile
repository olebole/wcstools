CFLAGS= -g
CC=cc
LIBWCS = libwcs/libwcs.a
LIBS = $(LIBWCS) -lm
.PRECIOUS: ${LIBWCS}

all:	delwcs edhead gethead i2f imcat imgsc imhead imrot imsize imstar \
	imuac imwcs scat sethead addpix getpix setpix sgsc sky2xy \
	skycoor suac subpix xy2sky

addpix: addpix.c $(LIBWCS) libwcs/fitshead.h
	$(CC) $(CFLAGS) -o addpix addpix.c $(LIBS)

delwcs: delwcs.c $(LIBWCS) libwcs/fitshead.h
	$(CC) $(CFLAGS) -o delwcs delwcs.c $(LIBS)

edhead: edhead.c $(LIBWCS) libwcs/fitshead.h
	$(CC) $(CFLAGS) -o edhead edhead.c $(LIBS)

gethead: gethead.c $(LIBWCS) libwcs/fitshead.h
	$(CC) $(CFLAGS) -o gethead gethead.c $(LIBS)

getpix: getpix.c $(LIBWCS) libwcs/fitshead.h
	$(CC) $(CFLAGS) -o getpix getpix.c $(LIBS)

i2f: i2f.c $(LIBWCS) libwcs/fitshead.h
	$(CC) $(CFLAGS) -o i2f i2f.c $(LIBS)

imcat: imcat.c $(LIBWCS) libwcs/fitshead.h libwcs/wcs.h
	$(CC) $(CFLAGS) -o imcat imcat.c $(LIBS)

imgsc: imgsc.c $(LIBWCS) libwcs/fitshead.h libwcs/wcs.h
	$(CC) $(CFLAGS) -o imgsc imgsc.c $(LIBS)

imhead: imhead.c $(LIBWCS) libwcs/fitshead.h libwcs/wcs.h
	$(CC) $(CFLAGS) -o imhead imhead.c $(LIBS)

imrot: imrot.c $(LIBWCS) libwcs/fitshead.h
	$(CC) $(CFLAGS) -o imrot imrot.c $(LIBS)

imsize: imsize.c $(LIBWCS) libwcs/fitshead.h libwcs/wcs.h
	$(CC) $(CFLAGS) -o imsize imsize.c $(LIBS)

imstar: imstar.c $(LIBWCS) libwcs/fitshead.h libwcs/wcs.h libwcs/lwcs.h
	$(CC) $(CFLAGS) -o imstar imstar.c $(LIBS)

imuac: imuac.c $(LIBWCS) libwcs/fitshead.h libwcs/wcs.h
	$(CC) $(CFLAGS) -o imuac imuac.c $(LIBS)

imwcs: imwcs.c $(LIBWCS) libwcs/fitshead.h libwcs/lwcs.h
	$(CC) $(CFLAGS) -o imwcs imwcs.c $(LIBS)

scat: scat.c $(LIBWCS) libwcs/fitshead.h libwcs/wcs.h
	$(CC) $(CFLAGS) -o scat scat.c $(LIBS)

sethead: sethead.c $(LIBWCS) libwcs/fitshead.h
	$(CC) $(CFLAGS) -o sethead sethead.c $(LIBS)

setpix: setpix.c $(LIBWCS) libwcs/fitshead.h
	$(CC) $(CFLAGS) -o setpix setpix.c $(LIBS)

sgsc: sgsc.c $(LIBWCS) libwcs/fitshead.h libwcs/wcs.h
	$(CC) $(CFLAGS) -o sgsc sgsc.c $(LIBS)

sky2xy: sky2xy.c $(LIBWCS) libwcs/wcs.h
	$(CC) $(CFLAGS) -o sky2xy sky2xy.c $(LIBS)

skycoor: skycoor.c $(LIBWCS) libwcs/wcs.h
	$(CC) $(CFLAGS) -o skycoor skycoor.c $(LIBS)

suac: suac.c $(LIBWCS) libwcs/fitshead.h libwcs/wcs.h
	$(CC) $(CFLAGS) -o suac suac.c $(LIBS)

subpix: subpix.c $(LIBWCS) libwcs/fitshead.h
	$(CC) $(CFLAGS) -o subpix subpix.c $(LIBS)

xy2sky: xy2sky.c $(LIBWCS) libwcs/wcs.h
	$(CC) $(CFLAGS) -o xy2sky xy2sky.c $(LIBS)

$(LIBWCS): libwcs/*.c libwcs/fitshead.h libwcs/wcs.h
	cd libwcs; make

