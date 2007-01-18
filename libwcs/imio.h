/*** imio.h  memory access subroutines
 *** January 5, 2007
 *** By Doug Mink, dmink@cfa.harvard.edu
 *** Harvard-Smithsonian Center for Astrophysics
 *** Copyright (C) 1996-2007
 *** Smithsonian Astrophysical Observatory, Cambridge, MA, USA

    This library is free software; you can redistribute it and/or
    modify it under the terms of the GNU Lesser General Public
    License as published by the Free Software Foundation; either
    version 2 of the License, or (at your option) any later version.

    This library is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    Lesser General Public License for more details.
    
    You should have received a copy of the GNU Lesser General Public
    License along with this library; if not, write to the Free Software
    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

    Correspondence concerning WCSTools should be addressed as follows:
           Internet email: dmink@cfa.harvard.edu
           Postal address: Doug Mink
                           Smithsonian Astrophysical Observatory
                           60 Garden St.
                           Cambridge, MA 02138 USA
 */

#ifndef imio_h_
#define imio_h_

/* Image pixel access subroutines in imio.c */

#ifdef __cplusplus /* C++ prototypes */
extern "C" {
#endif

#ifdef __STDC__   /* Full ANSI prototypes */

/* Image pixel access subroutines in imio.c */

    double getpix(	/* Read one pixel from any data type 2-D array (0,0)*/
	char *image,	/* Image array as 1-D vector */
	int bitpix,	/* FITS bits per pixel
			 *  16 = short, -16 = unsigned short, 32 = int
			 * -32 = float, -64 = double */
	int w,		/* Image width in pixels */
	int h,		/* Image height in pixels */
	double bzero,	/* Zero point for pixel scaling */
	double bscale,	/* Scale factor for pixel scaling */
	int x,		/* Zero-based horizontal pixel number */
	int y);		/* Zero-based vertical pixel number */
    double getpix1(	/* Read one pixel from any data type 2-D array (1,1)*/
	char *image,	/* Image array as 1-D vector */
	int bitpix,	/* FITS bits per pixel */
	int w,		/* Image width in pixels */
	int h,		/* Image height in pixels */
	double bzero,	/* Zero point for pixel scaling */
	double bscale,	/* Scale factor for pixel scaling */
	int x,		/* One-based horizontal pixel number */
	int y);		/* One-based vertical pixel number */
    double maxvec(	/* Get maximum value in vector from a image */
	char *image,	/* Image array from which to extract vector */
	int bitpix,	/* Number of bits per pixel in image */
	double bzero,	/* Zero point for pixel scaling */
	double bscale,	/* Scale factor for pixel scaling */
	int pix1,	/* Offset of first pixel to extract */
	int npix);	/* Number of pixels to extract */
    void putpix(	/* Write one pixel to any data type 2-D array (0,0)*/
	char *image,	/* Image array as 1-D vector */
	int bitpix,	/* FITS bits per pixel */
	int w,		/* Image width in pixels */
	int h,		/* Image height in pixels */
	double bzero,	/* Zero point for pixel scaling */
	double bscale,	/* Scale factor for pixel scaling */
	int x,		/* Zero-based horizontal pixel number */
	int y,		/* Zero-based vertical pixel number */
	double dpix);	/* Value to put into image pixel */
    void putpix1(	/* Write one pixel to any data type 2-D array (1,1) */
	char *image,	/* Image array as 1-D vector */
	int bitpix,	/* FITS bits per pixel */
	int w,		/* Image width in pixels */
	int h,		/* Image height in pixels */
	double bzero,	/* Zero point for pixel scaling */
	double bscale,	/* Scale factor for pixel scaling */
	int x,		/* One-based horizontal pixel number */
	int y,		/* One-based vertical pixel number */
	double dpix);	/* Value to put into image pixel */
    void addpix(	/* Add to one pixel in any data type 2-D array (0,0)*/
	char *image,	/* Image array as 1-D vector */
	int bitpix,	/* FITS bits per pixel */
	int w,		/* Image width in pixels */
	int h,		/* Image height in pixels */
	double bzero,	/* Zero point for pixel scaling */
	double bscale,	/* Scale factor for pixel scaling */
	int x,		/* Zero-based horizontal pixel number */
	int y,		/* Zero-based vertical pixel number */
	double dpix);	/* Value to add to image pixel */
    void addpix1(	/* Add to one pixel in any data type 2-D array (1,1)*/
	char *image,	/* Image array as 1-D vector */
	int bitpix,	/* FITS bits per pixel */
	int w,		/* Image width in pixels */
	int h,		/* Image height in pixels */
	double bzero,	/* Zero point for pixel scaling */
	double bscale,	/* Scale factor for pixel scaling */
	int x,		/* One-based horizontal pixel number */
	int y,		/* One-based vertical pixel number */
	double dpix);	/* Value to add to image pixel */
    void movepix(	/* Move one pixel value between two 2-D arrays (0,0) */
	char *image1,	/* Pointer to first pixel in input image */
	int bitpix1,	/* Bits per input pixel (FITS codes) */
	int w1,		/* Number of horizontal pixels in input image */
	int x1,		/* Zero-based row for input pixel */
	int y1,		/* Zero-based column for input pixel */
	char *image2,	/* Pointer to first pixel in output image */
	int bitpix2,	/* Bits per output pixel (FITS codes) */
	int w2,		/* Number of horizontal pixels in output image */
	int x2,		/* Zero-based row for output pixel */
	int y2);	/* Zero-based column for output pixel */
    void movepix1(	/* Move one pixel value between two 2-D arrays (1,1) */
	char *image1,	/* Pointer to first pixel in input image */
	int bitpix1,	/* Bits per input pixel (FITS codes) */
	int w1,		/* Number of horizontal pixels in input image */
	int x1,		/* One-based row for input pixel */
	int y1,		/* One-based column for input pixel */
	char *image2,	/* Pointer to first pixel in output image */
	int bitpix2,	/* Bits per output pixel (FITS codes) */
	int w2,		/* Number of horizontal pixels in output image */
	int x2,		/* One-based row for output pixel */
	int y2);	/* One-based column for output pixel */

/* Image vector processing subroutines in imio.c */

    void addvec(	/* Add constant to vector from 2-D array */
	char *image,	/* Image array as 1-D vector */
	int bitpix,	/* FITS bits per pixel */
	double bzero,	/* Zero point for pixel scaling */
	double bscale,	/* Scale factor for pixel scaling */
	int pix1,	/* Offset of first pixel to which to add */
	int npix,	/* Number of pixels to which to add */
	double dpix);	/* Value to add to pixels */
    void multvec(	/* Multiply vector from 2-D array by a constant */
	char *image,	/* Image array as 1-D vector */
	int bitpix,	/* FITS bits per pixel */
	double bzero,	/* Zero point for pixel scaling */
	double bscale,	/* Scale factor for pixel scaling */
	int pix1,	/* Offset of first pixel to multiply */
	int npix,	/* Number of pixels to multiply */
	double dpix);	/* Value to add to pixels */
    void getvec(	/* Read vector from 2-D array */
	char *image,	/* Image array as 1-D vector */
	int bitpix,	/* FITS bits per pixel */
	double bzero,	/* Zero point for pixel scaling */
	double bscale,	/* Scale factor for pixel scaling */
	int pix1,	/* Offset of first pixel to extract */
	int npix,	/* Number of pixels to extract */
	double *dvec0);	/* Vector of pixels (returned) */
    void putvec(	/* Write vector into 2-D array */
	char *image,	/* Image array as 1-D vector */
	int bitpix,	/* FITS bits per pixel */
	double bzero,	/* Zero point for pixel scaling */
	double bscale,	/* Scale factor for pixel scaling */
	int pix1,	/* Offset of first pixel to insert */
	int npix,	/* Number of pixels to insert */
	double *dvec0);	/* Vector of pixels to insert */
    void fillvec(	/* Write constant into a vector */
	char *image,	/* Image array as 1-D vector */
	int bitpix,	/* FITS bits per pixel */
	double bzero,	/* Zero point for pixel scaling */
	double bscale,	/* Scale factor for pixel scaling */
	int pix1,	/* Zero-based offset of first pixel to multiply */
	int npix,	/* Number of pixels to multiply */
	double dpix);	/* Value to which to set pixels */
    void fillvec1(	/* Write constant into a vector */
	char *image,	/* Image array as 1-D vector */
	int bitpix,	/* FITS bits per pixel */
	double bzero,	/* Zero point for pixel scaling */
	double bscale,	/* Scale factor for pixel scaling */
	int pix1,	/* One-based offset of first pixel to multiply */
	int npix,	/* Number of pixels to multiply */
	double dpix);	/* Value to which to set pixels */

/* Image pixel byte-swapping subroutines in imio.c */

    void imswap(	/* Swap alternating bytes in a vector */
	int bitpix,	/* Number of bits per pixel */
	char *string,	/* Address of starting point of bytes to swap */
	int nbytes);	/* Number of bytes to swap */
    void imswap2(	/* Swap bytes in a vector of 2-byte (short) integers */
	char *string,	/* Address of starting point of bytes to swap */
	int nbytes);	/* Number of bytes to swap */
    void imswap4(	/* Reverse bytes in a vector of 4-byte numbers */
	char *string,	/* Address of starting point of bytes to swap */
	int nbytes);	/* Number of bytes to swap */
    void imswap8(	/* Reverse bytes in a vector of 8-byte numbers */
	char *string,	/* Address of starting point of bytes to swap */
	int nbytes);	/* Number of bytes to swap */
    int imswapped();	/* Return 1 if machine byte order is not FITS order */

#else /* __STDC__ */ /* K&R declarations */

extern double getpix(); /* Read one pixel from any data type 2-D array (0,0)*/
extern double getpix1(); /* Read one pixel from any data type 2-D array (1,1)*/
extern void putpix();   /* Write one pixel to any data type 2-D array (0,0)*/
extern void putpix1();  /* Write one pixel to any data type 2-D array (1,1) */
extern void addpix();   /* Add to one pixel in any data type 2-D array (0,0)*/
extern void addpix1();  /* Add to one pixel in any data type 2-D array (1,1)*/
extern void movepix();  /* Move one pixel value between two 2-D arrays (0,0) */
extern void movepix1(); /* Move one pixel value between two 2-D arrays (1,1) */
extern void getvec();   /* Read vector from a 2-D array */
extern void putvec();   /* Write vector into a 2-D array */
extern void fillvec();   /* Write constant into a vector */
extern void fillvec1();   /* Write constant into a vector */
extern void imswap();   /* Swap alternating bytes in a vector */
extern void imswap2();  /* Swap bytes in a vector of 2-byte (short) integers */
extern void imswap4();  /* Reverse bytes in a vector of 4-byte numbers */
extern void imswap8();  /* Reverse bytes in a vector of 8-byte numbers */
extern int imswapped(); /* Return 1 if machine byte order is not FITS order */

#endif  /* __STDC__ */

#ifdef __cplusplus
}
#endif  /* __cplusplus */

#endif	/* imio_h_ */

/* May 31 1996	Use stream I/O for reading as well as writing
 * Jun 12 1996	Add byte-swapping subroutines
 * Aug  6 1996	Add MOVEPIX, HDEL and HCHANGE declarations
 *
 * May 27 1998	Split off imio subroutines to imio.h

 * Sep 27 1999	Add Fortran-indexed (1,1), not (0,0) image access *1()
 * Sep 28 1999	Add addpix()
 *
 * Feb 27 2004	Add fillvec()
 *
 * Jan  5 2007	Add prototype declarations for ANSI C and C++
 */
