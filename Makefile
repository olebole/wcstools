CLDFLAGS= -g
CFLAGS= $(CLDFLAGS) -Ilibwcs
LIBWCS = libwcs/libwcs.a
LIBS = $(LIBWCS) -lm

setwcs: setwcs.o $(LIBWCS) libwcs/fitshead.h
	cc $(CLDFLAGS) -o $@ setwcs.o $(LIBS)

delwcs: delwcs.o $(LIBWCS) libwcs/fitshead.h
	cc $(CLDFLAGS) -o $@ delwcs.o $(LIBS)

xy2sky: xy2sky.o $(LIBWCS)
	cc $(CLDFLAGS) -o $@ xy2sky.o $(LIBS)

sky2xy: sky2xy.o $(LIBWCS)
	cc $(CLDFLAGS) -o $@ sky2xy.o $(LIBS)

imstar: imstar.o $(LIBWCS) libwcs/fitshead.h libwcs/wcs.h
	cc $(CLDFLAGS) -o $@ imstar.o $(LIBS)

libwcs/libwcs.a: libwcs/*.c
	cd libwcs; make
