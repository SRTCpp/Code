#include <strstream>
#include <iostream>
#include <valarray>

#include <string.h>
#include "Jcrap.h"
#include <fitsio.h>


void Jcrap::printfitserror( int status)
{
    /*****************************************************/
    /* Print out cfitsio error messages and exit program */
    /*****************************************************/


    if (status)
    {
       fits_report_error(stderr, status); /* print error report */
//       exit( status );    /* terminate the program, returning error status */
    }
    return;
}



void cube::FITShdrread(char const d[])
{
//cfitsio way
// open the file
	fitsfile *fptr(0);
	int anynull, status(0);
	cout << "Opening FITS file\n"; cout.flush();
	if ( fits_open_file(&fptr, d, READONLY, &status) )
		printfitserror( status );
	  
	int keysexist, morekeys, hdus, hdutype;
	fits_get_num_hdus(fptr, &hdus, &status);
	for (int hdu(1);hdu<=hdus;hdu++) {
		fits_get_hdrspace(fptr, &keysexist, &morekeys, &status);
		cout << "Got header info HDU " << hdu << ":  " << keysexist << " keys exist\n";
		for (int i(1);i<=keysexist;i++) {
			char keyname[8000], keyvalue[8000], comment[8000];
			fits_read_keyn(fptr, i, keyname, keyvalue, comment, &status);
			if (keyword(keyname).size()>0)
				strcat(keyname, (string("_HDU")+int2str(hdu)).c_str());
			(*this).keyword(keyname, keyvalue);
			cout << i << " Setting " << keyname << " = " << keyvalue << "\n";
		}
		
		if (hdus > hdu) 
			fits_movrel_hdu(fptr, +1, &hdutype, &status);
		
	}
	
//	cout << "Read this info so far:  \n";
//	cout << infile << "\n";
	
	nx=ny=nz=1;
	if (str2int(keyword("NAXIS")) > 0) 
		nx = str2int(keyword("NAXIS1"));
	if (str2int(keyword("NAXIS")) > 1) 
		ny = str2int(keyword("NAXIS2"));
	if (str2int(keyword("NAXIS")) > 2) 			
		ny = str2int(keyword("NAXIS3"));

	fits_close_file(fptr, &status);
}

void cube::FITSreaddata(char const d[])
{
	int pixels=N(X)*N(Y)*N(Z);
//	cout << "before FITS instantiation\n";
	
//cfitsio way
	fitsfile *fptr(0);
	int fitsstatus(0), anynull;
	if ( fits_open_file(&fptr, d, READONLY, &fitsstatus) )
		printfitserror( fitsstatus );
	
	valarray<float> valdata(pixels+1);

	cout << "before data read\n";
	float nullval;   // dunno what this does. seems cfitsio wants it tho
	if ( fits_read_img(fptr, TFLOAT, 1, pixels, &nullval,
        &valdata[0], &anynull, &fitsstatus) )	// cfitsio
		printfitserror( fitsstatus );
	
	cout << "DAta read in.  First values = " << valdata[3] << "\n";
	
	if (osuppress<1) printpercent(100,100);
	if (osuppress<1) cout << "\n";
	int i=0;
	for (int z=0;z<N(Z);z++)
		for (int y=0;y<N(Y);y++)
			for (int x=0;x<N(X);x++,i++)
				(*this)(x,y,z) = valdata[i];
	
	fits_close_file(fptr, &fitsstatus);
}

void cube::FITSwrite(char const d[]) const
{/*
	long naxis=0;
	long naxes[3];	
	
	naxis = long(N(X)>1) + long(N(Y)>1) + long(N(Z)>1);
	
	int i=0;
	if (N(X)>1) { naxes[i] = N(X); i++; }
	if (N(Y)>1) { naxes[i] = N(Y); i++; }
	if (N(Z)>1) { naxes[i] = N(Z); i++; }
	
	FITS jfits(d, FLOAT_IMG, naxis, naxes);
	
	valarray<float> datavalarr(N(X)*N(Y)*N(Z));
	i=0;
	for (int z=0;z<N(Z);z++)
		for (int y=0;y<N(Y);y++)
			for (int x=0;x<N(X);x++,i++)
				datavalarr[i] = (*this)(x,y,z);
	
	for (map<string, string, less<string> >::iterator i=hdr.begin();i!=hdr.end();i++) {
		if (i->first == string("BITPIX")) { hdr.erase(i); continue; }
		if (i->first == string("EXTEND")) { hdr.erase(i); continue; }
		if (i->first == string("SIMPLE")) { hdr.erase(i); continue; }
		if (i->first == string("BSCALE")) { hdr.erase(i); continue; }
		if (i->first == string("BZERO")) { hdr.erase(i); continue; }
		int hdrvaltype=0;		// kludgy.  Key:  0 = string  1 = float  2 = int  3 = bool
		if (i->second[0] == '\'' ) {
			hdrvaltype=0;
			string outstr(i->second.substr(1,1000));
			for (int j=0;j < outstr.size();j++) 
				if (outstr[j] == '\'' ) outstr[j] = '\0';
			jfits.pHDU().addKey(i->first, outstr, "");
		}
		else {
			hdrvaltype = 2;
			for (int j=0;j < i->second.size();j++) {
			if (i->second[j] == 'T') {
				jfits.pHDU().addKey(i->first, bool(1), "");
				hdrvaltype = 3;
				break;
			}
			if (i->second[j] == 'F') {
				jfits.pHDU().addKey(i->first, bool(0), "");
				hdrvaltype = 3;
				break;
			}
			if (!isdigit(i->second[j]) && i->second[j]!=' ') hdrvaltype = 1;
			}
		}
		if (hdrvaltype == 1)
			jfits.pHDU().addKey(i->first, str2float(i->second), "");
		if (hdrvaltype == 2)
			jfits.pHDU().addKey(i->first, str2int(i->second), "");
//		cout << i->first << " = " << i->second << " -- type("<<hdrvaltype <<")\n";
	}

//	cout << jfits.pHDU() << "\n";
	jfits.pHDU().write(1, N(X)*N(Y)*N(Z), datavalarr);

	*/
}
