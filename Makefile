CFLAGS= -g
CC=cc
LIBWCS = libwcs/libwcs.a
LIBS = $(LIBWCS) -lm
.PRECIOUS: ${LIBWCS}

all:	imwcs delwcs xy2sky sky2xy imstar imrot i2f imgsc imujc imsize imhead imtab

delwcs: delwcs.c $(LIBWCS) libwcs/fitshead.h
	$(CC) $(CFLAGS) -o $@ delwcs.c $(LIBS)

edhead: edhead.c $(LIBWCS) libwcs/fitshead.h
	$(CC) $(CFLAGS) -o $@ edhead.c $(LIBS)

gethead: gethead.c $(LIBWCS) libwcs/fitshead.h
	$(CC) $(CFLAGS) -o $@ gethead.c $(LIBS)

i2f: i2f.c $(LIBWCS) libwcs/fitshead.h
	$(CC) $(CFLAGS) -o $@ i2f.c $(LIBS)

imcat: imcat.c $(LIBWCS) libwcs/fitshead.h libwcs/wcs.h
	$(CC) $(CFLAGS) -o $@ imcat.c $(LIBS)

imgsc: imgsc.c $(LIBWCS) libwcs/fitshead.h libwcs/wcs.h
	$(CC) $(CFLAGS) -o $@ imgsc.c $(LIBS)

imhead: imhead.c $(LIBWCS) libwcs/fitshead.h libwcs/wcs.h
	$(CC) $(CFLAGS) -o $@ imhead.c $(LIBS)

imrot: imrot.c $(LIBWCS) libwcs/fitshead.h
	$(CC) $(CFLAGS) -o $@ imrot.c $(LIBS)

imsize: imsize.c $(LIBWCS) libwcs/fitshead.h libwcs/wcs.h
	$(CC) $(CFLAGS) -o $@ imsize.c $(LIBS)

imstar: imstar.c $(LIBWCS) libwcs/fitshead.h libwcs/wcs.h libwcs/lwcs.h
	$(CC) $(CFLAGS) -o $@ imstar.c $(LIBS)

imtab: imtab.c $(LIBWCS) libwcs/fitshead.h libwcs/wcs.h
	$(CC) $(CFLAGS) -o $@ imtab.c $(LIBS)

imujc: imujc.c $(LIBWCS) libwcs/fitshead.h libwcs/wcs.h
	$(CC) $(CFLAGS) -o $@ imujc.c $(LIBS)

imwcs: imwcs.c $(LIBWCS) libwcs/fitshead.h libwcs/lwcs.h
	$(CC) $(CFLAGS) -o $@ imwcs.c $(LIBS)

sky2xy: sky2xy.c $(LIBWCS) libwcs/wcs.h
	$(CC) $(CFLAGS) -o $@ sky2xy.c $(LIBS)

xy2sky: xy2sky.c $(LIBWCS) libwcs/wcs.h
	$(CC) $(CFLAGS) -o $@ xy2sky.c $(LIBS)

$(LIBWCS): libwcs/*.c libwcs/fitshead.h libwcs/wcs.h
	cd libwcs; make

