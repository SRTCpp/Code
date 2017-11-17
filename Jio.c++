#include<string>
#include<cstdio>
#include<map>
#include<iostream>
#include<math.h>
#include<strstream>
#include<sstream>
#include<fstream>
#include<utility>
#include<algorithm>
#include<vector>
#include<string.h>

#include<gdal_priv.h>

#include"Jcube.h"

extern int graphno;

using namespace Jcrap;

ostream& operator<<(ostream& out, ::filetype& f)
{
	string filestr;
	if (f == vicar)  filestr = "vicar";
	if (f == isis)   filestr = "isis";
	if (f == isis3)   filestr = "isis3";
	if (f == text)   filestr = "text";
	if (f == Jcube1) filestr = "Jcube1";
	if (f == jpg)    filestr = "jpg";
	if (f == gif)    filestr = "gif";
	if (f == tif)    filestr = "tif";
	if (f == fits)   filestr = "fits";
	if (f == Jcube1a)filestr = "Jcube1a";
	if (f == Jcube1b)filestr = "Jcube1b";
	if (f == bsqraw) filestr = "BSQraw";
	if (f == graphout) filestr = "graphout";
	
	return out << filestr;
}

ostream& operator<<(ostream& out, const cube& c)
{
	out << "OUTPUTTING ENTIRE CUBE:  DIMS (" << c.Nx() << ", " << c.Ny()
			<< ", " << c.Nz() << ")\n";
	for (int z=0;z<c.Nz();z++){
		for (int y=0;y<c.Ny();y++){
			for (int x=0;x<c.Nx();x++){
				out << "{" << c.Axis(X,x) << "," << c.Axis(Y,y) << "," << c.Axis(Z,z) << "}" << c(x,y,z) << "  ";
			}
			out << "\n";
		}
		out << "---------------------------------------\n";
	}
	return (out <<
			"===========================================================\n");
}

ostream& operator<<(ostream& out, axis a)
{
	char achar;
	if (a==X) achar = 'X';
	if (a==Y) achar = 'Y';
	if (a==Z) achar = 'Z';
	if (a==UNK) achar='?';
	return (out << achar);
}

ostream& operator<<(ostream& out, datatype d)
{
	string dname;
	if (d==NONE) dname="NONE";
	if (d==DN) dname="DN";
	if (d==FLUX) dname="W / m\\sp2\\ep\\*mm";
	if (d==PHOTONS) dname="Photons";
	if (d==FLUXCGS) dname="erg / s*cm\\sp2\\ep*cm";
	if (d==FLUXCGSWAVENUM) dname="erg / s*cm\\sp2\\ep*cm\\sp-1\\ep";
	if (d==FLUXCGSHZ) dname="erg / s*cm\\sp2\\ep*Hz\\sp-1\\ep";
	if (d==ARB) dname="Arbitrary Units";
	if (d==KM)			dname="Diameter (km)";
	if (d==PHOTONFLUX) dname="Photons / s / \\ep\\*mm";
	return (out << dname);
}

ostream& operator<<(ostream& out, axistype t)
{
	string tname;
	if (t==ARRAYVAL) 	tname="Array Values";
	if (t==MICRONS) 	tname="Wavelength in Microns (\\*mm)";
	if (t==WAVENUM) 	tname="Wavelength in Wavenumbers (cm\\sp-1\\ep)";
	if (t==CM) 			tname="Wavelength in Centimeters (cm)";
	if (t==M)			tname="Wavelength in Meters (m)";
	if (t==KMPERS)    tname="\\*Dv (km/s)";
	if (t==TIME)		tname="Time (hours)";
	if (t==HZ)			tname="Frequency (1/T)";
	return (out << tname);
}

void cube::write(string outfile, float minval, float maxval, float (*s)(float)) 
{
	write(outfile.c_str(), minval, maxval, s);
}

void cube::write(const char outfile[1000], float minval, float maxval, float (*s)(float)) 
{
	cout <<"inside write\n";
	if (cubetype==jpg || cubetype==gif || cubetype==tif || cubetype==tif32) {
		cout << "calling IMWrite\n";cout.flush();
		IMwrite(outfile,minval,maxval,s);
	}
	else if (cubetype==fits){
		 FITSwrite(outfile);
		 cout<<"In the fits file write."<<endl;
	} 
	else {
		ofstream fout(outfile);
		cout << outfile << " -- ";
      hdr["Nx"]=int2str(Nx());
      hdr["Ny"]=int2str(Ny());
      hdr["Nz"]=int2str(Nz());
		putheader(&fout);
		writedata(&fout);
	}
}

void cube::graph(const char* arg, otype dest)
{
	static vector<string> argrecord;
	string argstr(arg);		
	if (keyword("graph_linecolor").size()) {	
		cout << "Graphing in a special color:  \"" << keyword("graph_linecolor") << "\"\n";

		argstr += " --pen-color \"";
		argstr += int2str((graphno%5)+1);
		argstr += "=";
		argstr += keyword("graph_linecolor");
		argstr += "\" -m ";
		argstr += int2str((graphno%5)+1);
		argstr += " ";
	}
	argrecord.push_back(argstr);
	
	bool combine=0;
	if (arg[0]<=int('9') && arg[0]>int('0')) {
		combine=1;
		argstr = argstr.substr(argstr.find_first_of(" "));
	}
	
	FILE *gfile;
	char gname[12];
	string command("graph ");
	
	if (!combine) {
		sprintf(gname, "graph%d.out", graphno);
		command += argrecord[graphno];
//		command += " ";
//		command += gname;
		command += " ";
		graphno++;
		gfile=fopen(gname, "w");
		if (gfile==0) {
			cout << "FAILED TO OPEN GRAPH FILE -- PERMISSION ERROR?\n";
			return; // don't proceed into a segfault
		}
	}
	else {
		int prevno;
		sscanf(arg, "%d", &prevno);
		cout << "Combining previous " << prevno << " graphs.\n";
		
		for (int i=prevno;i>=1;i--) {
			sprintf(gname, "graph%d.out", graphno-i);
			command += argrecord[graphno-i];
			command += " ";
			command += gname;
			command += " ";
		}
		graphno++;
	}

	
	// special titan graph mode
	if (keyword("titangraph") == string("yes")) {
		string command("graph -T fig --blankout 0.0 -C -E x -E y -X \"Wavelength (microns)\" -Y \"I / F\" ");
		command += string("-f 0.11 -W 0.008 ");
		
		cout << "titangraphing\n"; cout.flush();
		
	// 1.6 um
		command += string("--reposition +0.64 +0.4 0.4 -w 0.55 -h 1.2 -x 1.52 1.62 0.03 -y 0. 0.28");
		int prevno;
		if (combine)
			sscanf(arg, "%d", &prevno);
		else prevno = 1;
		
	
		for (int i=prevno+1;i>1;i--) {
			sprintf(gname, "graph%d.out", graphno-i);
			command += argrecord[graphno-i];
			command += " ";
			command += gname;
			command += " ";
		}
		command += " -X \"\" -Y \"\" -E y ";
		
	// 1.26 um
		command += "--reposition +0.42 +0.4 0.4 -w 0.55 -h 1.2 -x 1.21 1.31 0.03 -y 0. 0.28 ";
		for (int i=prevno+1;i>1;i--) {
			sprintf(gname, "graph%d.out", graphno-i);
			command += argrecord[graphno-i];
			command += " ";
			command += gname;
			command += " ";
		}
		
	// 1.05 um
		command += "--reposition +0.20 +0.4 0.4 -w 0.55 -h 1.2 -x 1.01 1.10 0.03 -y 0. 0.28 ";
		for (int i=prevno+1;i>1;i--) {
			sprintf(gname, "graph%d.out", graphno-i);
			command += argrecord[graphno-i];
			command += " ";
			command += gname;
			command += " ";
		}
		command += " -Y \"I / F\" ";
		
	// 0.92 um
		command += " --reposition -0.02 +0.4 0.4 -w 0.55 -h 1.2 -x 0.87 0.96 0.03 -y 0. 0.28 ";
		for (int i=prevno+1;i>1;i--) {
			sprintf(gname, "graph%d.out", graphno-i);
			command += argrecord[graphno-i];
			command += " ";
			command += gname;
			command += " ";
		}
		command += " -E x -E y -X \"Wavelength (microns)\" -Y \"I / F\" -f 0.0525 ";
		
	// 5 um
		command += " --reposition +0.42 +0. 0.4 -w 1.1 -h 1.0 -x 4.8 5.1 0.05 -y 0 0.12 ";
		for (int i=prevno+1;i>1;i--) {
			sprintf(gname, "graph%d.out", graphno-i);
			command += argrecord[graphno-i];
			command += " ";
			command += gname;
			command += " ";
		}
		command += " -X "" -Y "" -E y -f 0.11 ";
		
	// 2.7 um
		command += " --reposition +0.20 +0. 0.4 -w 0.55 -h 1.0 -x 2.6 2.82 0.05 -y 0 0.12 ";
		for (int i=prevno+1;i>1;i--) {
			sprintf(gname, "graph%d.out", graphno-i);
			command += argrecord[graphno-i];
			command += " ";
			command += gname;
			command += " ";
		}
		command += "-Y \"I / F\"";
		
	// 2.0 um
		command += " --reposition -0.02 +0. 0.4 -w 0.55 -h 1.0 -x 1.9 2.1 0.05 -y 0 0.19 ";
		for (int i=prevno+1;i>1;i--) {
			sprintf(gname, "graph%d.out", graphno-i);
			command += argrecord[graphno-i];
			command += " ";
			command += gname;
			command += " ";
		}
		command += " > titangraph";
		command += int2str(graphno);
		command += ".fig";
		
		
//		graphno++;
		cout << command << "\n\n\n";
		system(command.c_str());

//		ofstream commandfile("graph.commands", ios::out | ios::app);
//		commandfile << command << "\n";

		return;		
	}
	
	char labels[200];
	for (int i=0;i<200;i++) labels[i]='\0';
	strstream labelstr(labels, 199, ios::out); 
	labelstr << "-L \"" << hdr["title"] << "\" ";
	if (ny==1 && nz==1){
//		cout << "Graphing using only x\n";
		labelstr << "-X \"";
		if (hdr["xlabel"].size()) labelstr << hdr["xlabel"];
		else labelstr << Xtype();
		labelstr << "\" -Y \"";
		if (hdr["ylabel"].size()) labelstr << hdr["ylabel"];
		else labelstr << Dtype();
		labelstr << "\" \0"; 
		if (!combine) for (int i=0;i<nx;i++) fprintf(gfile, "%le %e\n", Xaxis(i), (*this)(i,0,0));
	}
	else if (nx==1 && nz==1){
		labelstr << "-X \"" << Ytype() << "\" -Y \"" << Dtype() << "\" \0";
		if (!combine) for (int i=0;i<ny;i++) fprintf(gfile, "%le %e\n", yaxis[i], (*this)(0,i,0));
	}
	else if (nx==1 && ny==1){
		labelstr << "-X \"" << Ztype() << "\" -Y \"" << Dtype() << "\" \0";

		// graph does funny things if all of the exponents are the same.  Trying to mess w/that.
		double max=Axis(Z,0), min=Axis(Z,0), dmin=(*this)(0,0,0), dmax=dmin;
		for (int i=0;i<N(Z);i++) {
			if (Axis(Z,i) > max) max = Axis(Z,i);
			if (Axis(Z,i) < min) min = Axis(Z,i);
			if ((*this)(0,0,i) > dmax) dmax = (*this)(0,0,i);
			if ((*this)(0,0,i) < dmin) dmin = (*this)(0,0,i);
		}
		bool z10=0, d10=0;
		if ((int)log10(max) == (int)log10(min)) z10=1;
		if ((int)log10(dmax) == (int)log10(dmin)) d10=1;
		
		if (!combine) for (int i(0);i<N(Z);i++) {
			if (Axis(Z,i) == 0.0) fprintf(gfile, "0.000 %e\n", (*this)(0,0,i));
			if (z10) fprintf(gfile, "%lf ", Axis(Z,i));
			else fprintf(gfile, "%le ", Axis(Z,i));
			if (d10) fprintf(gfile, "%lf\n", (*this)(0,0,i));
			else fprintf(gfile, "%le\n", (*this)(0,0,i));
		}
	} else if (N(Y)==1 && N(X)==2 && keyword("WITHERROR").size()) {
		cout << "Graphing with ERROR!\n";
		labelstr << "-X \"" << Ztype() << "\" -Y \"" << Dtype() << "\" -I e \0";
		if (!combine) {
			for (int i=0;i<nz;i++) {
				if ((*this)(0,0,i) > 1e5  ||  (*this)(0,0,i) < 1e-5)
					fprintf(gfile, "%le %e %e\n", Axis(Z,i), (*this)(0,0,i), (*this)(1,0,i));
				else 
					fprintf(gfile, "%lf %f %f\n", Axis(Z,i), (*this)(0,0,i), (*this)(1,0,i));
			}
			fprintf(gfile, "\n");
		}
	} else if (ny==1){
		labelstr << "-X \"" << Ztype() << "\" -Y \"" << Dtype() << "\" \0";
		if (!combine) for (int j=0;j<nx;j++){
			for (int i=0;i<nz;i++) fprintf(gfile, "%le %e\n", zaxis[i], (*this)(j,0,i));
			fprintf(gfile, "\n");
		}
	}
	
	if (!combine) fprintf(gfile, "\n");
	if (!combine) fclose(gfile);

	command += labels;
	command += " ";
	command += gname;
	command += " ";
	if (dest==oX) command+="-T X&";
	if (dest==oPS){
		command += "-T ps> ";
		command += gname;
		command += ".ps&";
	}
	if (dest==oFIG){
		command += "-T fig> ";
		command += gname;
		command += ".fig&";
	}
	if (dest==oGIF){
		command += "-T ps> ";
		command += gname;
		command += ".ps ; convert ";
		command += gname;
		command += ".ps gif:";
		command += gname;
		command += ".gif&";
	}

//	cout << "Issuing command:  "<<command<<"\n";
	ofstream commandfile("graph.commands", ios::out | ios::app);
	commandfile << command << "\n";
	system (command.c_str());
	
}
	
	
			
void cube::determinetype(ifstream *Pfin, const char fname[])
{
	char lbltest[8];
	for (int i=0;i<7;i++) lbltest[i]=(*Pfin).get();
	lbltest[7]='\0';
	char extension[4];
	strcpy(extension, &fname[strlen(fname)-3]);
	for (int i=0;i<3;i++) extension[i]=tolower(extension[i]);
	cout << "extension = \"" << extension << "\"\n";
	string slbl(lbltest);
	if (strcmp(lbltest, "LBLSIZE") == 0) cubetype=vicar;
	else if (strcmp(extension, "txt") == 0) cubetype=text;
	else if (strcmp(extension, "gif") == 0) cubetype=gif;
	else if (strcmp(extension, "jpg") == 0) cubetype=jpg;
	else if (strcmp(extension, "peg") == 0) cubetype=jpg;
	else if (strcmp(extension, "tif") == 0) cubetype=tif;
	else if (strcmp(extension, "its") == 0) cubetype=fits;
	else if (strcmp(extension, "out") == 0) cubetype=graphout;
	else if (slbl == "Jcube1b") {
		cubetype=Jcube1b;
		hdr["Jcubever"]="1.2";
	}
	else if (slbl == "Jcube1a") {
		cubetype=Jcube1a;
		hdr["Jcubever"]="1.1";
	}
	else if (slbl == "Jcube1 ") {
		cubetype=Jcube1;
		hdr["Jcubever"]="1.0";
	}
	else if (slbl == "Object ") {
		cubetype=isis3;
	}
	else cubetype=isis;
	
	for (int i=6;i>=0;i--) (*Pfin).putback(lbltest[i]);

}
	
int cube::getlabelsize(ifstream *Pfin)
{
	int L;
	
	if(cubetype==vicar){
		string instr;
		for (int i=0;i<1024;i++)
			instr+=Pfin->get();
		
		char cinstr[2000];
		strcpy(cinstr, instr.c_str());
	
		int nlocation=7, n2location=7;
		for (int i=7;instr[i]<'0' || instr[i]>'9';i++) nlocation=i;
		
		sscanf(&cinstr[nlocation], "%d " , &L);
		
		
		for (int i=instr.length()-1; i >= 0; i--) (*Pfin).putback(instr[i]);
		cout << "LBLsize measured to be " << L << "\n";
	}
	
	if(cubetype==isis){
		char tmpheader[1000];
		int i=0,j=0;
		while (strncmp(&tmpheader[i-12], "LABEL_RECORDS", 13)){
			i++;
			tmpheader[i]=(*Pfin).get();
//			cout << i << "=" << tmpheader[i] << '\n';
		}
		while (tmpheader[i+j] != '=') {
			j++;
			tmpheader[i+j]=(*Pfin).get();
//			cout << i+j << "=" << tmpheader[i+j] << '\n';
		}	
		int k=1; 
		tmpheader[i+j+k]=(*Pfin).get();
		while (isspace(tmpheader[i+j+k])) {
			k++;
			tmpheader[i+j+k]=(*Pfin).get();
//			cout << i+j+k << "=" << tmpheader[i+j+k] << '\n';
		}
		while (!isspace(tmpheader[i+j+k])) {
			k++;
			tmpheader[i+j+k]=(*Pfin).get();
//			cout << i+j+k << "=" << tmpheader[i+j+k] << '\n';
		}
		sscanf(&tmpheader[i+j], "=%d", &L);
		for (int m=i+j+k;m>0;m--) { 
//			cout << "Returning " << m << "=" << tmpheader[m] << '\n';
			(*Pfin).putback(tmpheader[m]);}
		L*=512;
	}
		
	
	if(cubetype==isis3){  // added 2015 September 8 for NLDSAR data
		char tmpheader[10000];
		int i=0,j=0;
		while (strncmp(&tmpheader[i-8], "StartByte", 9)){
			i++;
			tmpheader[i]=(*Pfin).get();
//			cout << i << "=" << tmpheader[i] << '\n';
		}
		while (tmpheader[i+j] != '=') {
			j++;
			tmpheader[i+j]=(*Pfin).get();
//			cout << i+j << "=" << tmpheader[i+j] << '\n';
		}	
		int k=1; 
		tmpheader[i+j+k]=(*Pfin).get();
		while (isspace(tmpheader[i+j+k])) {
			k++;
			tmpheader[i+j+k]=(*Pfin).get();
//			cout << i+j+k << "=" << tmpheader[i+j+k] << '\n';
		}
		while (!isspace(tmpheader[i+j+k])) {
			k++;
			tmpheader[i+j+k]=(*Pfin).get();
//			cout << i+j+k << "=" << tmpheader[i+j+k] << '\n';
		}
		sscanf(&tmpheader[i+j], "=%d", &L);
		for (int m=i+j+k;m>0;m--) { 
//			cout << "Returning " << m << "=" << tmpheader[m] << '\n';
			(*Pfin).putback(tmpheader[m]);}
	}
	
		
	if (cubetype==Jcube1){
		char tmpheader[10000];
		int i=0;
		while (strncmp(&tmpheader[i-8], "STARTDATA", 9)){
			i++;
			tmpheader[i]=(*Pfin).get();
		}
		L=i+1;
		for (i=i;i>0;i--) (*Pfin).putback(tmpheader[i]);
	}
	
	cout.flush();
	return L;
}

void cube::readheader(ifstream *Pfin, const char fname[])
{
	int lblsize=getlabelsize(Pfin);
	int ourlabelsize=0;						//  The size of our LBLSIZE when we output the cube
													//  will be different from the incube size, maybe
	if (osuppress < -2) cout << "LBLSIZE = " << lblsize << '\n';
	
	if(cubetype==vicar){
		string instr;
		for (int i=0; i<lblsize-1; i++) {static char c='0'; c = (*Pfin).get(); instr+=c;}

	
		cout << "size=" << lblsize << " -- \"" << instr << "\"\n";
		
		int eqloc;
		while ( (eqloc=instr.find_first_of('=')) != -1){
			int i=eqloc;
			int chars=3;

			while (!isspace(instr[--i]) && i > 0);
			if(i==0)i=-1;
			string lbl=instr.substr(i+1, eqloc-i-1);
			chars+=eqloc-i-1;


			i=eqloc;			
			if (instr[eqloc+1]!=0x27)
				while (!isspace(instr[i]) && ++i < lblsize);
			else {
				i+=2;
				while (instr[i]!=0x27 && ++i < lblsize);
				i++;
			}
			string value=instr.substr(eqloc+1, i-eqloc-1);
			chars+=i-eqloc-1;
			
			hdr[lbl]=value;
//			cout << "Assigning " << lbl << "=" << value << " size " << chars << '\n';
			instr[eqloc] = '|';
			ourlabelsize += chars;
		}
//		cout << "Oldlabelsize with 2 spaces for each:  " << ourlabelsize << '\n';
		ourlabelsize+=5-hdr["LBLSIZE"].length();	// allow more room for changing LBLSIZEs.
		hdr["LBLSIZE"]=lblstr(ourlabelsize);									// making the new LBLSIZE
		if(debug)cout << "NEWLABEL:::" << hdr["LBLSIZE"] << ":::\n";
		if (hdr["INTFMT"]=="'HIGH'") byteorder=bigendian;
		else byteorder=littleendian;
//		cout << "VICAR header read in.\n";
	}
	if(cubetype==isis  ||  cubetype==isis3){ 
		string instr="";
		for (int i=0;i<lblsize;i++) instr += (*Pfin).get();
		
		int eqloc;
		if (osuppress < -2) cout << "In isis readheader, osuppress=" << osuppress << "\n";
		if (osuppress < -2) cout << "lblsize = " << lblsize << "\n";
		while ( (eqloc=instr.find_first_of('=')) != -1){
			int i=eqloc-1;
			
			while (isspace(instr[i])) i--;  // take up space before the =
			int endlbl(i);
			while (!isspace(instr[--i]) && i>0);
			if(i==0)i=-1;
			string lbl=instr.substr(i+1, endlbl-i);
			
			string value;
			i=eqloc+2;
			if (instr[i]=='"') {
				instr[i]='`';
				int newi;
				while ((newi=instr.find_first_of('"'))<=i) instr[newi]='|';
				i=instr.find_first_of('"');
				instr[eqloc+2]='"';
				value=instr.substr(eqloc+2, i-eqloc-1);
				instr[eqloc+2]='|';
				instr[i]='|';
			}
			else if (instr[i]=='('){
				int newi;
				while ((newi=instr.find_first_of(')'))<=i) instr[newi]='|';
				newi=instr.find_first_of(')');
				i=instr.find_first_of(')');;
				value=instr.substr(eqloc+2, i-eqloc-1);
				instr[eqloc+2]='|';
				instr[i]='|';
			}
			else {
				while (instr[i]!='\n' && ++i<lblsize);
				i--;
				value=instr.substr(eqloc+2, i-eqloc-1);
			}
			if (osuppress < -2) cout << "value starts at " << eqloc+1 << " and goes " << i-eqloc << "?\n"; 
			
			hdr[lbl]=value;
			if (osuppress < -2) cout << "Assigning " << lbl << "=" << value << '\n';
			instr[eqloc] = '|';	
			
		}
		hdr["StartByte"]=int2str(lblsize);
		
	}
	
	if (cubetype==text) {
		char instr[1000];
		(*Pfin).getline(instr, 999);
//		(*Pfin) >> nx >> ny >> nz;
		cout << "(" << instr << ")\n";
		sscanf(instr, "%ld %ld %ld", &nx, &ny, &nz);
		cout << "after sscanf, size = ("<<nx<<","<<ny<<","<<nz<<")\n";
	}
	
			
	if (cubetype==Jcube1 || cubetype==Jcube1a || cubetype==Jcube1b) {
		char instr[10000];
		char lblname[100];
		char value[10000];
		int j;
		char c;
		
		(*Pfin).get(instr, 9999);
		(*Pfin).get(instr[j=strlen(instr)]);
		instr[j+1]='\0';
//		cout << "instr " << strlen(instr) << ":  " << instr << "\n";
		
		(*Pfin).get(lblname, 99, ' ');
		(*Pfin).get(value, 9999);
		(*Pfin).get(c);	
		while (strncmp(lblname,"STARTDATA",9)){
//			cout << "lblname:  " << lblname << "\n";
//			cout << "value:    " << value << "\n";
			hdr[lblname]=value;
			if (string(value).find('(') != string::npos &&
				 string(value).find(')') == string::npos) {
				do {
					(*Pfin).get(value, 9999);
					(*Pfin).get(c);
					hdr[lblname] += value;
//					cout << "   +value: " << value << "\n";
				} while (string(value).find(')') == string::npos);
			}
			(*Pfin).get(lblname, 99, ' ');
			(*Pfin).get(value, 9999);
			(*Pfin).get(c);
		}
//		cout << "Header read in successfully!\n";
	}
	
	if (cubetype == fits) {
		FITShdrread(fname);
	}
			
		
				
	begindata = Pfin->tellg();	
		
}


short cube::getshort(ifstream *Pfin)
{
	char byte[2];
	
	if (byteorder==littleendian) for (int i=0; i<2;i++) byte[i]=(*Pfin).get();
	if (byteorder==bigendian)    for (int i=1;i>=0;i--) byte[i]=(*Pfin).get();
	short value;
	void* Pint;
	void* Pchars;
	Pchars=byte;
	Pint=&value;
	
	memcpy(Pint, Pchars, 2);
	
//	cout << "reading " << value << " from " << byte[0] << byte[1] << '\n';
	return value;
}

long cube::get32int(ifstream *Pfin)
{
	char byte[4];
	for (int i=0;i<4;i++) byte[i]=(*Pfin).get();
	
	long value;
	void* Pint;
	void* Pchars;
	Pchars=byte;
	Pint=&value;
	memcpy(Pint, Pchars, 4);
//	cout << "reading " << value << " from " << byte[0] << byte[1] << '\n';
	return value;
}

namespace Jcrap {
	inline void getfloat(ifstream* Pfin, void *pv)	
	/* created 2/22/2k.  Faster, more efficient, new invocation */
	{
		char* byte=(char*)pv;
		for (int i=0;i<4;i++) byte[i]=(*Pfin).get();
	}
	
	inline void getfloatasdouble(ifstream* Pfin, double &thevalue)
	/* created 2012 August 22 for the new double-sized array axis values */
	{
		float floatie_temp;
		char* byte=(char*)((void*)(&floatie_temp));
		for (int i=0;i<4;i++) byte[i]=(*Pfin).get();
		thevalue = double(floatie_temp);
	}
	
	inline void getdouble(ifstream* Pfin, double &thevalue)
	/* created 2012 August 22 for the new double-sized array axis values */
	{
		char* byte=(char*)((void*)(&thevalue));
		for (int i=0;i<8;i++) byte[i]=(*Pfin).get();
	}

	
	inline void putdouble(ofstream* fout, const double& value)
	/* Created 2012 August 22 for new double-sized Axis values */
	{
		const void* Pchars=&value;
		const char* byte=(const char*)Pchars;
		for (int i=0;i<8;i++) (*fout).put(byte[i]);
	}

	
	inline void getbigendint(ifstream* &Pfin, void *pv)
	/* created 3/3/2k for fits files.  Soon to be deprecated by Jcube 2.0 I hope */
	{
		unsigned int inint;
		char *pinint=(char *)&inint;
		
		for (int i=3;i>=0;i--) pinint[i]=(*Pfin).get();
			
		float *d=(float *)pv;
		*d=(float)inint;
	}
	
	float getunsignedchar(ifstream* Pfin)
	{
		unsigned char inchar;
		char *pinchar=(char *)&inchar;
	
		*pinchar=(*Pfin).get();
				
		
		return static_cast<float>(inchar);
	}
	
	
	float getfloat(ifstream *Pfin)
	{
		char byte[4];
		for (int i=0;i<4;i++) byte[i]=(*Pfin).get();
		
		float value;
		void* Pfloat;
		void* Pchars;
		Pchars = byte;
		Pfloat = &value;
		memcpy(Pfloat, Pchars, 4);
		return value;
	}
	
	float getbigendianfloat(ifstream *Pfin)
	{
		char byte[4];
		for (int i=0;i<4;i++) byte[3-i]=(*Pfin).get();
		float value;
		void* Pfloat;
		void* Pchars;
		Pchars = byte;
		Pfloat = &value;
		memcpy(Pfloat, Pchars, 4);
		return value;
	}	
	
	float getbigendianshort(ifstream *Pfin)
	{
//		cout << "Getbigendshort\n";cout.flush();
		char byte[2];
		for (int i=0;i<2;i++) byte[1-i]=(*Pfin).get();
		
		short value;
		void* Pfloat;
		void* Pchars;
		Pchars = byte;
		Pfloat = &value;
		memcpy(Pfloat, Pchars, 2);
//		cout << "Getbigendshort done\n";cout.flush();
		return float(value);
	}	
	
	float getbigendianint(ifstream *Pfin)
	{
//		cout << "Getbigendint\n";cout.flush();
		char byte[4];
		for (int i=0;i<4;i++) byte[3-i]=(*Pfin).get();
		
		int value;
		void* Pfloat;
		void* Pchars;
		Pchars = byte;
		Pfloat = &value;
		memcpy(Pfloat, Pchars, 4);
//		cout << "Getbigendint done\n";cout.flush();
		return float(value);
	}
}

/* Old and inefficient:  deprecated 2/22/2k but left here for backward
compatibility concerns */
float cube::getfloat(ifstream *Pfin)
{
	char byte[4];
	for (int i=0;i<4;i++) byte[i]=(*Pfin).get();
	
	float value;
	void* Pfloat;
	void* Pchars;
	Pchars = byte;
	Pfloat = &value;
	memcpy(Pfloat, Pchars, 4);
	return value;
}


/* Added 2004 June 12 trying to get vicar readins to work */
float cube::getbigendfloat(ifstream *Pfin)
{
	char byte[4];
	for (int i=0;i<4;i++) byte[3-i]=(*Pfin).get();
	
	float value;
	void* Pfloat;
	void* Pchars;
	Pchars = byte;
	Pfloat = &value;
	memcpy(Pfloat, Pchars, 4);
	return value;
}
			
void cube::readdata(ifstream *Pfin, const char fname[])
{
	if (osuppress<-2) cout << "reading data\n";
	if (osuppress<-2) cout << "nx=" << nx << " ny=" << ny << " nz=" << nz << "\n";

	storageSLOW=X; storageMEDIUM=Y; storageFAST=Z;  // default
	if (keyword("Storage_Order").find("XYZ") != string::npos); // stick with default
	else if (keyword("Storage_Order").find("XZY") != string::npos) {	
		storageSLOW=X; storageMEDIUM=Z; storageFAST=Y;
	} else if (keyword("Storage_Order").find("YXZ") != string::npos) {
		storageSLOW=Y; storageMEDIUM=X; storageFAST=Z;
	} else if (keyword("Storage_Order").find("YZX") != string::npos) {// a.k.a. BIL	
		storageSLOW=Y; storageMEDIUM=Z; storageFAST=X;
	} else if (keyword("Storage_Order").find("ZXY") != string::npos) {
		storageSLOW=Z; storageMEDIUM=X; storageFAST=Y;
	} else if (keyword("Storage_Order").find("ZYX")!=string::npos) {// a.k.a. BSQ	
		storageSLOW=Z; storageMEDIUM=Y; storageFAST=X;
	} else if (hdr["AXIS_NAME"]=="(SAMPLE,BAND,LINE)" || hdr["AXES_NAME"]=="(SAMPLE,BAND,LINE)"
		|| hdr["AXIS_NAME"]=="(SAMPLE, BAND, LINE)") {
		storageSLOW=Y; storageMEDIUM=Z; storageFAST=X;
	} else if (hdr["AXIS_NAME"]=="(SAMPLE,LINE,BAND)" || hdr["AXES_NAME"]=="(SAMPLE,LINE,BAND)") {
		storageSLOW=Z; storageMEDIUM=Y; storageFAST=X;
	} else if (keyword("Format").find("BandSequential")!=string::npos ||
		 keyword("ByteOrder").find("Lsb")!=string::npos) {// a.k.a. BSQ	
		storageSLOW=Z; storageMEDIUM=Y; storageFAST=X;
	}
		

	createarrays();
	
	xaxis = new double[nx];
	yaxis = new double[ny];
	zaxis = new double[nz];
	
	if (!osuppress) cout << "Reading  0%"; 
	if (!osuppress) cout.flush();
	
	if(cubetype==vicar){
		for(int x=0;x<nx;x++) xaxis[x]=double(x);		//setting up the axis calibration
		for(int y=0;y<ny;y++) yaxis[y]=double(y);		//labels.
		for(int z=0;z<nz;z++) zaxis[z]=double(z);
		
		cout << "ORG = " << hdr["ORG"] << "\n";
		if (hdr["ORG"] == "'BSQ'") {
			storageSLOW=Z; storageMEDIUM=Y; storageFAST=X;
			for(int z=0;z<nz;z++){
				if (!osuppress) printpercent(z,nz-1);
				for(int y=0;y<ny;y++){
					for(int x=0;x<nx;x++){
						(*this)(x,y,z)=getfloat(Pfin);
					}
				}
			}
		}		
		
		if (hdr["ORG"] == "'BIL'") {
			storageSLOW=Y; storageMEDIUM=Z; storageFAST=X;
			for(int y=0;y<ny;y++){
				if (!osuppress) printpercent(y,ny-1);
				for(int z=0;z<nz;z++){
					for(int x=0;x<nx;x++){
						(*this)(x,y,z)=getfloat(Pfin);
					}
				}
			}
		}
	}

	
	if (cubetype==isis3) {
		
// ISIS3		
		if (keyword("MinimumLatitude").size()) {
			cout << "MinimumLatitude detected -- assigning Y axis labels\n";
			float minlat(-str2float(keyword("MaximumLatitude")));  // opposite to change to south lat convention
			float maxlat(-str2float(keyword("MinimumLatitude")));
			float increment((maxlat-minlat)/float(N(Y)-1));		//setting up the axis calibration
			for(int y=0;y<ny;y++) Axis(Y,y) = minlat + increment*float(y);
		}
		if (keyword("MaximumLongitude").size()) {
			cout << "MaximumLongitude detected -- assigning X axis labels\n";
			float minlon(-str2float(keyword("MinimumLongitude")));  // opposite to change to east long convention
			float maxlon(-str2float(keyword("MaximumLongitude")));
			float increment(fabs(maxlon-minlon)/float(N(X)-1));
			for(long x(0);x<nx;x++) Axis(X,x) = minlon + increment*float(x);
			keyword("cylindrical_map")=string("yes");
		}
		cout << " data size is " << Data.size() << "\n";
	}
		
	if(cubetype==isis){
		for(long x(0);x<nx;x++) xaxis[x]=double(x);		//setting up the axis calibration
		for(long y(0);y<ny;y++) yaxis[y]=double(y);		//labels.
		for(long z(0);z<nz;z++) zaxis[z]=double(z);

// ISIS2		
		if (keyword("MINIMUM_LATITUDE").size()) {
			cout << "MINIMUM_LATITUDE detected -- assigning Y axis labels\n";
			float minlat(-str2float(keyword("MAXIMUM_LATITUDE")));  // opposite to change to south lat convention
			float maxlat(-str2float(keyword("MINIMUM_LATITUDE")));
			float increment((maxlat-minlat)/float(N(Y)-1));		//setting up the axis calibration
			for(int y=0;y<ny;y++) Axis(Y,y) = minlat + increment*float(y);
		}
		if (keyword("WESTERNMOST_LONGITUDE").size()) {
			cout << "WESTERNMOST_LONGITUDE detected -- assigning X axis labels\n";
			float minlon(-str2float(keyword("WESTERNMOST_LONGITUDE")));  // opposite to change to east long convention
			float maxlon(-str2float(keyword("EASTERNMOST_LONGITUDE")));
			float increment((maxlon-minlon)/float(N(X)-1));
			for(long x(0);x<nx;x++) Axis(X,x) = minlon + increment*float(x);
			keyword("cylindrical_map")=string("yes");
		}

		
		if (cubetype==isis) begindata = (str2int(hdr["^QUBE"])-1)*512;
		if (cubetype==isis3)begindata = str2int(hdr["StartByte"])-1;
		if (cubetype==isis) cout << " QUBE=\n" << hdr["^QUBE"] << "\n, ISIS begindata = " << begindata << "\n";
		if (cubetype==isis3)cout << " StartByte=" << hdr["StartByte"] << ", ISIS begindata = " << begindata << "\n";
		Pfin->seekg(begindata,ios::beg);
	
		float (*isisget)(ifstream*s) = &Jcrap::getfloat;  // assume littleendian reals by default
		if (hdr["CORE_ITEM_TYPE"].find("MAC_UNSIGNED_INTEGER") != string::npos) {
			isisget = &getunsignedchar;
			if (!osuppress) cout << "MAC_unsigned\n";
		}
		if (hdr["CORE_ITEM_TYPE"].find("MAC_REAL") != string::npos) {
			isisget = &getbigendianfloat;
			if (!osuppress) cout << "MAC_real\n";
		}
		if (hdr["CORE_ITEM_TYPE"].find("SUN_INTEGER") != string::npos) {
			if (hdr["CORE_ITEM_BYTES"].find("2") != string::npos) {
				isisget = &getbigendianshort;
				if (!osuppress) cout << "SUN_INTEGER2\n";
			}
			if (hdr["CORE_ITEM_BYTES"].find("4") != string::npos) {
				isisget = &getbigendianint;
				if (!osuppress) cout << "SUN_INTEGER4\n";
			}
		}
				
		float (*sideget)(ifstream*s) = &Jcrap::getfloat;
		if (hdr["SAMPLE_SUFFIX_ITEM_TYPE"].find("MAC_UNSIGNED_INTEGER") != string::npos) {
			sideget = &getunsignedchar;
			if (!osuppress) cout << "side = MAC_unsigned\n";
		}
		if (hdr["SAMPLE_SUFFIX_ITEM_TYPE"].find("MAC_REAL") != string::npos) {
			sideget = &getbigendianfloat;
			if (!osuppress) cout << "side = MAC_real\n";
		}
		if (hdr["CORE_ITEM_TYPE"].find("SUN_INTEGER") != string::npos) {
			if (hdr["SAMPLE_SUFFIX_ITEM_BYTES"].find("2") != string::npos) {
				sideget = &getbigendianshort;
				if (!osuppress) cout << "side = SUN_INTEGER2\n";
			}
			if (hdr["SAMPLE_SUFFIX_ITEM_BYTES"].find("4") != string::npos) {
				sideget = &getbigendianint;
				if (!osuppress) cout << "side = SUN_INTEGER4\n";
			}
		}
		
		if (sidecube().N(X)==bottomcube().N(Y)==backcube().N(Z)==0 && keyword("Rotation").size()==0) {
			/*if (osuppress<0)*/cout << "data order " << storageSLOW << storageMEDIUM << storageFAST;
			if (osuppress<0)cout << " reading in automatically (NO SIDEPLANES)\n";
			for(long n100(0);n100<N(storageSLOW);n100++){
				if (osuppress<1) printpercent(n100,N(storageSLOW)-1);
				for(long n10(0);n10<N(storageMEDIUM);n10++){
					if (osuppress<1 && N(storageSLOW)==1) printpercent(n10, N(storageMEDIUM)-1);
					for(long n1(0);n1<N(storageFAST);n1++){
						storageorderaccess(n100,n10,n1)=isisget(Pfin);
					}
				}
			}
		} else if (sidecube().N(X)==bottomcube().N(Y)==backcube().N(Z)==0 
				&& keyword("Rotation").find("90")!=string::npos &&
				storageSLOW==Z && storageMEDIUM==Y && storageFAST==X) {
			if (osuppress<0)cout << " reading in with rotation\n";
			for(long z(0);z<N(Z);z++){
				if (osuppress<1) printpercent(z,N(Z)-1);
				for(long x(0);x<N(X);x++){
					if (osuppress<1 && N(Z)==1) printpercent(x, N(X)-1);
					for(long y(0);y<N(Y);y++){
						(*this)(N(X)-1-x,y,z)=isisget(Pfin);
					}
				}
			}
		} else if(storageSLOW==Y && storageMEDIUM==Z && storageFAST==X) {
			for(long y(0);y<ny+Bottom();y++){
				if (!osuppress) printpercent(y,ny-1);
				for(long z(0);z<nz+Back();z++){
					for(long x(0);x<nx+Side();x++){
//						cout <<" REading in ("<<x<<","<<y<<","<<z<<")\n";
						if ( x<nx && y<ny && z<nz ) {
							(*this)(x,y,z)=isisget(Pfin);
						}else{
							if ((x>nx-1) + (y>ny-1) + (z>nz-1) > 2) continue;
							else if ((x>nx-1) + (y>ny-1) + (z>nz-1) > 1)
								sideget(Pfin); // throw one away
							else if (x>nx-1) sidecube()(x-nx,y,z)=sideget(Pfin);
							else if (y>ny-1) bottomcube()(x,y-ny,z)=sideget(Pfin);
							else if (z>nz-1) backcube()(x,y,z-nz)=sideget(Pfin);
						}
					}
				}
			}
		} else if (storageSLOW==Z && storageMEDIUM==Y && storageFAST==X) {
			for(long z(0);z<nz;z++){
				if (!osuppress) printpercent(z,nz-1);
				for(long y(0);y<ny;y++){
					if (!osuppress && N(Z)==1) printpercent(y, N(Y)-1);
					for(long x(0);x<nx;x++){
						if ( x < (nx-Side()) && z < (nz-Back())) {
							(*this)(x,y,z)=isisget(Pfin);
						}else{
							(*this)(x,y,z)=isisget(Pfin);
							isisget(Pfin);
						}
					}
				}
			}
		} else {
		}
		if (osuppress<-2) cout << "ISIS data read in\n";
	}	
	if (cubetype==text){
		for(int x=0;x<nx;x++) xaxis[x]=double(x);		//setting up the axis calibration
		for(int y=0;y<ny;y++) yaxis[y]=double(y);		//labels.
		for(int z=0;z<nz;z++) zaxis[z]=double(z);

		cout << "Reading in text data\n";
		for(int z=0;z<nz;z++){
			if (!osuppress) printpercent(z,nz-1);
			for(int y=0;y<ny;y++){
				for(int x=0;x<nx;x++){
					(*Pfin) >> (*this)(x,y,z);
//					cout << "Read in :  " << data[x][y][z] << "\n";
				}
			}
		}
	}	
	if (cubetype==graphout){ // reads in x y per line in text, like graph.out files
		string readin;
		int z(0);
		vector<float> x, y;
		while ((*Pfin) >> readin) {
			x.push_back(str2float(readin));
			(*Pfin) >> readin;
			y.push_back(str2float(readin));
			z++;
		}
		*this = cube(1,1,z);
		cubetype = graphout;
		cout << "For graphout file, determined geometry to be 1x1x"<<z<<".\n";
		
		for(int z=0;z<nz;z++){
			Axis(Z,z) = x.at(z);
			(*this)(0,0,z) = y.at(z);
//			cout << "Read in " << Axis(Z,z) << "  :  " << (*this)(0,y,z) << "\n";	
		}
	}	
	if (cubetype==Jcube1){
		for (int x=0;x<nx;x++) Jcrap::getfloatasdouble(Pfin, xaxis[x]);
		for (int y=0;y<ny;y++) Jcrap::getfloatasdouble(Pfin, yaxis[y]);
		for (int z=0;z<nz;z++) Jcrap::getfloatasdouble(Pfin, zaxis[z]);
		
		storageSLOW=Z; storageMEDIUM=Y; storageFAST=X;
		for (int z=0;z<nz;z++){
			if (!osuppress) printpercent(z,nz-1);
			for(int y=0;y<ny;y++){
				for(int x=0;x<nx;x++){
					(*this)(x,y,z)=getfloat(Pfin);
				}
			}
		}
	}
	if (cubetype==Jcube1a){
		for (int x=0;x<nx;x++) Jcrap::getfloatasdouble(Pfin, xaxis[x]);
		for (int y=0;y<ny;y++) Jcrap::getfloatasdouble(Pfin, yaxis[y]);
		for (int z=0;z<nz;z++) Jcrap::getfloatasdouble(Pfin, zaxis[z]);
		

		
		if (memdims==3) {
			for(cube::iterator i=begin();i!=end();i++)
				Jcrap::getfloat(Pfin,(void*)&(*i));
		} else {
			current.resize(3);
			current[X]=-1;
			current[Y]=-1;
			current[Z]=-1;
			begindata = (*Pfin).tellg();
		}
	}
	if (cubetype==Jcube1b){
		for (int x=0;x<nx;x++) Jcrap::getdouble(Pfin, xaxis[x]);
		for (int y=0;y<ny;y++) Jcrap::getdouble(Pfin, yaxis[y]);
		for (int z=0;z<nz;z++) Jcrap::getdouble(Pfin, zaxis[z]);
		

		
		if (memdims==3) {
			for(cube::iterator i=begin();i!=end();i++)
				Jcrap::getfloat(Pfin,(void*)&(*i));
		} else {
			current.resize(3);
			current[X]=-1;
			current[Y]=-1;
			current[Z]=-1;
			begindata = (*Pfin).tellg();
		}
	}
	if (cubetype==fits){
		for(int x=0;x<nx;x++) xaxis[x]=double(x);		//setting up the axis calibration
		for(int y=0;y<ny;y++) yaxis[y]=double(y);		//labels.
		for(int z=0;z<nz;z++) zaxis[z]=double(z);
		byteorder=littleendian;
		
		FITSreaddata(fname);
	}
	if (cubetype==isis3){  // GDAL:  added 2016 August 9 JWB
		GDALDataset *inDataset;
		inDataset = (GDALDataset *) GDALOpen(fname, GA_ReadOnly);
		

		
		if (memdims==3)  {
			for (int z(0);z<N(Z);z++) {
				
				GDALRasterBand *thisBand;
				thisBand = inDataset->GetRasterBand(z+1);
				thisBand->RasterIO( GF_Read, 0, 0, N(X), N(Y), &((*this)(0,0,z)), N(X), N(Y), GDT_Float32, 0, 0);
				
			}
			cout << "isis3 data read in\n";
		} else {
			//cout << "HOLY CRAP -- DON'T GO ANY FARTHER:  GDAL can't read in memdims!=3 cubes\n";
		}
	}
	
	if (!osuppress) printpercent(nz-1, nz-1);
	if (!osuppress) cout << "   --    Dims:  ";	

	return;
}

void putfitshdrline(ofstream *fout,map<string, string, less<string> >::iterator i)
{
	char hdrline[81]="                                                                                ";
	for (int c=0;c<8 && c<i->first.size();c++) hdrline[c] = i->first[c];
	hdrline[8] = '=';
	for (int c=10;c-10 < i->second.size();c++) {
		if (i->second[c-10] == '\t') i->second[c-10] = ' ';
		hdrline[c] = i->second[c-10];
	}
	*fout << hdrline;
}

void cube::putheader(ofstream *fout)
{  
	switch(cubetype){
		case bsqraw:
			break;
		case vicar:
		{
			if (header().find("INTFMT")!=header().end()) hdr.erase(hdr.find("INTFMT"));
			if (byteorder==littleendian) hdr["INTFMT"]="'LOW'";
			else hdr["INTFMT"]="'HIGH'";
			if (hdr.find("LBLSIZE")!=hdr.end()) hdr.erase(hdr.find("LBLSIZE"));
			int recs=0;
			int reallblsize;
			string headertext="";
			
// makeup the string that will become the header (the text part of it, anyway)
			headertext+="LBLSIZE=000 \n";		// placeholder
			for (map<const string, string>::iterator i=hdr.begin();i!=hdr.end();i++){
				headertext+=i->first;
				headertext+="=";
				headertext+=i->second;
				headertext+=" \n"; 
			}
			
			int txtsize;
			txtsize=headertext.size();	
			
			while (txtsize > (reallblsize=(recs*str2int(hdr["RECSIZE"])))) recs++;
			char txtsizestr[500];
			sprintf(txtsizestr, "%i", reallblsize);
			for (int i=0;i<5;i++) {
				if (txtsizestr[i]!='\0')
					headertext[8+i]=txtsizestr[i];
				else
					break;
			}
			
			int dif=reallblsize-txtsize;
			*fout << headertext;
			for (int j=0;j<dif;j++) *fout << " ";
			break;
		}
		case text:
		{
			*fout << nx << " " << ny << " " << nz << "\n";
			break;
		}
		case Jcube1:
		{
			ostringstream sorder;
			sorder << storageSLOW << storageMEDIUM << storageFAST;
			keyword("Storage_Order") = sorder.str();
			cout << "Storage Order:  \"" << sorder.str() << "\"\n";
			
			*fout << "Jcube1 \n";
/*			cout << "test:  " << hdr["Xaxistype"] << "\n";
			cout << "beforecall fout = " << fout << "\n";*/
			printpairendl<pair<const string, string> > P = for_each(hdr.begin(),hdr.end(),
					printpairendl<pair<const string, string> >(fout));
			/*cout << "out of new test\n";*/
			*fout << "STARTDATA -->\n";
			break;
		}
		case Jcube1a:
		{
			ostringstream sorder;
			sorder << storageSLOW << storageMEDIUM << storageFAST;
			keyword("Storage_Order") = sorder.str();
			cout << "Storage Order:  \"" << sorder.str() << "\"\n";
			cout << "SLOW = " << storageSLOW << " MED = " << storageMEDIUM << " FAST = " << storageFAST << "\n";
			
			*fout << "Jcube1a\n";
/*			cout << "test:  " << hdr["Xaxistype"] << "\n";
			cout << "beforecall fout = " << fout << "\n";*/
			printpairendl<pair<const string, string> > P = for_each(hdr.begin(),hdr.end(),
					printpairendl<pair<const string, string> >(fout));
			//cout << "out of new test\n";
			*fout << "STARTDATA -->\n";
			break;
		}	
		case Jcube1b:
		{
			ostringstream sorder;
			sorder << storageSLOW << storageMEDIUM << storageFAST;
			keyword("Storage_Order") = sorder.str();
			cout << "Storage Order:  \"" << sorder.str() << "\"\n";
			cout << "SLOW = " << storageSLOW << " MED = " << storageMEDIUM << " FAST = " << storageFAST << "\n";
			
			*fout << "Jcube1b\n";
/*			cout << "test:  " << hdr["Xaxistype"] << "\n";
			cout << "beforecall fout = " << fout << "\n";*/
			printpairendl<pair<const string, string> > P = for_each(hdr.begin(),hdr.end(),
					printpairendl<pair<const string, string> >(fout));
			//cout << "out of new test\n";
			*fout << "STARTDATA -->\n";
			break;
		}		
		case isis:
		{
/*			ostringstream sorder;
			sorder << storageSLOW << storageMEDIUM << storageFAST;
			keyword("Storage_Order") = sorder.str();
			cout << "Storage Order:  \"" << sorder.str() << "\"\n";
*/			
		    
		    /* ISIS cubes are multispectral image cubes.
		       The order of data is specified by the axes.  
		       For VIMS cubes, SAMPLE goes "across", LINE goes "down", BAND goes "into the page"
		       BAND is the spectral dimension.  At each seperate band there is an image.
		       When writing the file, these dimensions can be written in any order.
		       This should be user-selectable.
		       At this time, it is not user selectable.  The cube will be output BAND, SAMPLE, LINE.
		       
		       Note: we are not worrying about endian-ness (yet) in this cube.  The cube will
		       be stored in the endian-ness of the machine it is running on.
		       This should be user-selectable.

		       Each cube is specified by a number of FILE_RECORDS, consisting of a number of
		       RECORD_BYTES each.  It will be necessary to pad the last record block with spaces
		       if it is not an exact multiple of the RECORD_BYTES.  It will also be necessary to
		       pad the LABEL_RECORDS in the same way.
		       
		       The HISTORY section is where records will be placed that describe how the 
		       cube was generated.  This program will specify the cubes that were used to 
		       derive output cube products in the history label.  Also, the history will include
		       decisions about cube construction, such as co-averaging of pixels and the radius of
		       matching.

		       An ISIS cube is stored in this order: LABEL >> HISTORY >> DATA
		       Note: BOF means Beginning of File.
		       It is necessary to determine how large the LABEL is, in order to determine the
		       offset to the HISTORY record.  In order to determine the location offset of cube
		       DATA, it is necessary to determine how large the HISTORY record is.  
		       The line "LABEL_RECORDS = XXXX"
		       stores the size of the label, both HISTORY and QUBE, in units of RECORDS. 
		       The line "^HISTORY = YYYY" 
		       stores the offset to the cube history from BOF, in units of RECORDS.
		       The line "^QUBE = ZZZZ" 
		       stores the offset to the data from BOF, in units of RECORDS.
		       
		        1. To determine the above, we create the header first.  
		        2. Then, we count the number of bytes (characters) in the header, 
		           not including the HISTORY.  
		        3. If this is not a multiple of RECORD_BYTES, we pad the header with spaces until it 
		            is a multiple of RECORD_BYTES. 
		        4. We then set the HISTORY pointer. (Replace YYYY with the appropriate offset).
		        5. Next, add the history to the header string.  
		        6. If this is not a multiple of RECORD_BYTES, we pad the header again.
		        7. We then set the QUBE pointer. (Replace ZZZZ with the appropriate offset).
		        8. Set the LABEL_RECORDS count.
		        9. Set the FILE_RECORDS count.  The FILE_RECORDS count is not an exact amount.  It is 
		           guaranteed to be greater than or equal to the amount of records in the file.  This is 
		           because we do not pad the data so that it matches the FILE_RECORDS specification.  ISIS 
		           doesn't care.
		       10. Write the finished LABEL into the file.
		    */
		    cout << "An ISIS cube!  How scrumptious!\n";
		    int bands = nz; int lines = ny; int samples = nx;
		    
		    string s = "CCSD3ZF0000100000001NJPL3IF0PDS200000001 = CASSFDU_LABEL\r\n";
		    s += "\r\n";
		    s += "/* File Structure */\r\n";
		    s += "\r\n";
		    s += "RECORD_TYPE = " + hdr["RECORD_TYPE"] + "\r\n";
		    s += "RECORD_BYTES = " + hdr["RECORD_BYTES"] + "\r\n";
		    s += "FILE_RECORDS = WWWWWWWWWW\r\n";                      // Replace WWWWWWWWWW with actual file record count later
		    s += "LABEL_RECORDS = XXXX\r\n";                           // Replace XXXX with actual label record count later
		    s += "FILE_STATE = " + hdr["FILE_STATE"] + "\r\n";
		    s += "\r\n";
		    s += "/* Pointer to ISIS history label */\r\n";
		    s += "\r\n";
		    s += "^HISTORY = YYYY\r\n";                                // Replace YYYY with actual history record offset later
		    s += "OBJECT = HISTORY\r\n";
		    s += "END_OBJECT = HISTORY\r\n";
		    s += "\r\n";
		    s += "/* Qube structure: Standard ISIS Cube of VIMS Data */";
		    s += "\r\n";
		    s += "^QUBE = ZZZZ\r\n";                                   // Replace ZZZZ with cube label record offset later
		    s += "OBJECT = QUBE\r\n";
		    s += "  AXES = 3\r\n";                                     // Fixed at 3
		    s += "  AXIS_NAME = (SAMPLE,LINE,BAND)\r\n";               // Fixed at BAND, LINE, SAMPLE for now.
		    s += "\r\n";
		    s += "  /* Core description. */\r\n";
		    s += "  CORE_ITEMS = (" + int2str(samples) + "," + int2str(lines) + "," + int2str(bands) + ")\r\n";
		    s += "  CORE_ITEM_BYTES = " + hdr["CORE_ITEM_BYTES"] + "\r\n";
		    s += "  CORE_ITEM_TYPE = " + hdr["CORE_ITEM_TYPE"] + "\r\n";
		    s += "  CORE_BASE = " + hdr["CORE_BASE"] + "\r\n";
		    s += "  CORE_MULTIPLIER = " + hdr["CORE_MULTIPLIER"] + "\r\n";
		    s += "  CORE_VALID_MINIMUM = " + hdr["CORE_VALID_MINIMUM"] + "\r\n";
		    s += "  CORE_NULL = " + hdr["CORE_NULL"] + "\r\n";
		    s += "  CORE_LOW_REPR_SATURATION = " + hdr["CORE_LOW_REPR_SATURATION"] + "\r\n";
		    s += "  CORE_LOW_INSTR_SATURATION = " + hdr["CORE_LOW_INSTR_SATURATION"] + "\r\n";
		    s += "  CORE_HIGH_REPR_SATURATION = " + hdr["CORE_HIGH_REPR_SATURATION"] + "\r\n";
		    s += "  CORE_HIGH_INSTR_SATURATION = " + hdr["CORE_HIGH_INSTR_SATURATION"] + "\r\n";
		    s += "  CORE_MINIMUM_DN = " + hdr["CORE_MINIMUM_DN"] + "\r\n";
		    s += "  CORE_NAME = " + hdr["CORE_NAME"] + "\r\n";
		    s += "  CORE_UNIT = " + hdr["CORE_UNIT"] + "\r\n";
		    s += "\r\n";
		    s += "  /* Suffix description. */\r\n";
		    s += "  SUFFIX_ITEMS = (0,0,0)\r\n";                       // We don't store the suffix here
		    s += "  SUFFIX_BYTES = 4\r\n";                             // Fixed at 4 (for now)
		    s += "\r\n";
		    s += "  /* Target Information */\r\n";
		    s += "  TARGET_NAME = " + hdr["TARGET_NAME"] + "\r\n";
		    s += "\r\n";
		    s += "  /* Product Information */\r\n";
		    //
		    // User should change the product name
		    // We kludge this here to be the target plus a date of production
		    //
		    // The time format is like this: 2005-302T08:10:09.000Z
		    //
		   time_t tim=time(NULL);
		   tm *now=localtime(&tim);
		   char product_id[100], product_creation_time[100];
		   sprintf(product_creation_time, "%d-%03dT%02d:%02d:%02d.000Z", 
			   now->tm_year+1900, now->tm_yday, now->tm_mon+1, now->tm_mday, now->tm_min);
		   string starg = hdr["TARGET_NAME"];
		   // Get rid of quotes around the TARGET string
		   starg.replace(starg.find("\""),1,"");
		   starg.replace(starg.find("\""),1,"");
		   sprintf(product_id, "%s_%s", starg.c_str(), product_creation_time);
		   hdr["PRODUCT_ID"] = product_id;
		   hdr["PRODUCT_CREATION_TIME"] = product_creation_time;
		   hdr["PRODUCT_VERSION_TYPE"] = "FINAL";
		   s += "  PRODUCT_ID = \"" + hdr["PRODUCT_ID"] + "\"\r\n";
		   s += "  PRODUCT_CREATION_TIME = \"" + hdr["PRODUCT_CREATION_TIME"] + "\"\r\n";
		   s += "  PRODUCT_VERSION_TYPE = \"" + hdr["PRODUCT_VERSION_TYPE"] + "\"\r\n";
		   s += "\r\n";
		   s += "  /* Spectral axis description */\r\n";
		   s += "  GROUP = BAND_BIN\r\n";
		   //
		   // Think about how to change these based on band selection
		   // 
		   s += "  BAND_BIN_CENTER = " + hdr["BAND_BIN_CENTER"] + "\r\n";
		   s += "  BAND_BIN_UNIT = " + hdr["BAND_BIN_UNIT"] + "\r\n";
		   s += "  BAND_BIN_ORIGINAL_BAND = " + hdr["BAND_BIN_ORIGINAL_BAND"] + "\r\n";
			s += "  END_GROUP = BAND_BIN\r\n";
		   s += "END_OBJECT = QUBE\r\n";
		   s += "END\r\n";
		   // Pad the label
		   int records = pad(s, str2int(hdr["RECORD_BYTES"]));
		   cout << "The header length, which should be the final length, is: " << s.length() << "\n";
		   // Set the label records count
		   int pos = s.find("LABEL_RECORDS = XXXX");
		   string srepl = int2str(records);
		   int nocare = pad(srepl, 4);
		   s.replace(pos+16,4,srepl);
		   //
		   // Replace HISTORY placeholder with a pointer
		   // The HISTORY will come immediately after the LABEL, 
		   // so it starts at LABEL_RECORDS + 1
		   //
		   pos = s.find("^HISTORY = YYYY");
		   srepl = int2str(records+1);
		   nocare = pad(srepl, 4);
		   s.replace(pos+11,4,srepl);
		   // Insert history here
		   s += history();
		   // End the entire cube header
		   s += "END\r\n";
		   // Pad the entire cube label
		   records = pad(s, str2int(hdr["RECORD_BYTES"]));
			//
		   // Replace QUBE placeholder with a pointer
		   // The QUBE data will come immediately after the HISTORY,
		   // so it starts at LABEL_RECORDS+HISTORY_RECORDS+1
		   //
		   pos = s.find("^QUBE = ZZZZ");
		   srepl = int2str(records+1);
		   nocare = pad(srepl, 4);
		   s.replace(pos+8,4,srepl);
		   // Set the file records count
		   int total_bytes = (samples*bands*lines*str2int(hdr["CORE_ITEM_BYTES"])) + 
		     records*str2int(hdr["RECORD_BYTES"]);
		   int file_records_pad = total_bytes % str2int(hdr["RECORD_BYTES"]);
		   int file_records = (total_bytes/str2int(hdr["RECORD_BYTES"]));
		   if (file_records_pad > 0) {
		     file_records++;
		   }
		   pos = s.find("FILE_RECORDS = WWWWWWWWWW");
		   srepl = int2str(file_records);
		   nocare = pad(srepl, 10);
		   s.replace(pos+15,10,srepl);
		   int orignumber = 042002;
		   int mynumber = 015000;
		   printf("The original cube data starts at %d\n", orignumber);
		   printf("The header length in decimal is: %d, the header length according to od is: %d\n",s.length(), mynumber);
		   // Write the header
		   *fout << s;
		}
	}
}


/* Given a string, pad out the string (with spaces) until
   the string length is an integer multiple of the 
   blocksize desired.
   This function is primarily for padding ISIS header labels.
*/
int cube::pad(string &s, int blocksize)
{
  int strlength = 0;
  strlength = s.length();
  int extra = strlength % blocksize;
  int blocks = strlength / blocksize;
  int padlength = 0;
  if (extra > 0) {
    blocks++;
    padlength = (blocks*blocksize) - strlength;
    s.append(padlength, ' ');
  }
  printf("String length: %d, blocksize: %d, padlength: %d, blocks: %d\n",strlength, blocksize, padlength, blocks);
  return blocks;
}

void cube::writedata(ofstream *fout)
{
	cout << "Writing    0% ";
	cout.flush();
	
	switch (cubetype){
		case bsqraw:
			for(int z=0;z<nz;z++){
				printpercent(z,nz-1);
				for(int y=0;y<ny;y++){
					for(int x=0;x<nx;x++){
						putfloat(fout, (*this)(x,y,z));
					}
				}
			}
			printpercent(nz-1,nz-1);
			cout << '\n';
			break;
		case vicar:
			for(int z=0;z<nz;z++){
				printpercent(z,nz-1);
				for(int y=0;y<ny;y++){
					for(int x=0;x<nx;x++){
						putfloat(fout, (*this)(x,y,z));
					}
				}
			}
			printpercent(nz-1,nz-1);
			cout << '\n';
			break;
		case text:
			for (int z=0;z<nz;z++){
				printpercent(z, nz-1);
				for (int y=0;y<ny;y++){
					for (int x=0;x<nx;x++){
						*fout << (*this)(x,y,z) << " ";
					}
					printpercent(nz-1,nz-1);	
					*fout << "\n";
				}
			}
			break;
		case Jcube1:
			for (int x=0;x<nx;x++) putfloat(fout, float(xaxis[x]));
			for (int y=0;y<ny;y++) putfloat(fout, float(yaxis[y]));
			for (int z=0;z<nz;z++) putfloat(fout, float(zaxis[z]));
			for(int z=0;z<nz;z++){
				printpercent(z,nz-1);
				for(int y=0;y<ny;y++){
					for(int x=0;x<nx;x++){
						putfloat(fout, (*this)(x,y,z));
					}
				}
			}
			printpercent(nz-1,nz-1);
			cout << '\n';
			break;
		case Jcube1a:
			for (int x=0;x<nx;x++) putfloat(fout, float(xaxis[x]));
			for (int y=0;y<ny;y++) putfloat(fout, float(yaxis[y]));
			for (int z=0;z<nz;z++) putfloat(fout, float(zaxis[z]));
			
			for (cube::iterator i=begin();i!=end();i++)
				putfloat(fout, *i);
			printpercent(nz-1,nz-1);
			cout << '\n';
			break;
		case Jcube1b:
			for (int x=0;x<nx;x++) Jcrap::putdouble(fout, xaxis[x]);
			for (int y=0;y<ny;y++) Jcrap::putdouble(fout, yaxis[y]);
			for (int z=0;z<nz;z++) Jcrap::putdouble(fout, zaxis[z]);
			
			for (cube::iterator i=begin();i!=end();i++)
				putfloat(fout, *i);
			printpercent(nz-1,nz-1);
			cout << '\n';
			break;
		case isis:
		  //for (cube::iterator i=begin();i!=end();i++)
		  //	putfloat(fout, *i);
			for(int z=0;z<nz;z++){
				printpercent(z,nz-1);
				for(int y=0;y<ny;y++){
					for(int x=0;x<nx;x++){
						putfloat(fout, (*this)(x,y,z));
					}
				}
			}
			printpercent(nz-1,nz-1);
			cout << '\n';
			break;
	}
}	

void cube::putfloat(ofstream* fout, const float& value) const
/* This new putfloat was created 2/22/2k.  It reduces the number of times
the 4 byte floating point value is copied by 2 by using pass-by-reference
and not copying it again to the char[] for no apparent reason */
{
	const void* Pchars=&value;
	const char* byte=(const char*)Pchars;
	for (int i=0;i<4;i++) (*fout).put(byte[i]);
}

void cube::putfloat(fstream* fout, const float& value) const
/* This new putfloat was created 2/22/2k.  It reduces the number of times
the 4 byte floating point value is copied by 2 by using pass-by-reference
and not copying it again to the char[] for no apparent reason */
{
	const void* Pchars=&value;
	const char* byte=(const char*)Pchars;
	for (int i=0;i<4;i++) (*fout).put(byte[i]);
}

inline void cube::putshort(ofstream* &fout, short& value)
/* for fits output, 10.27.2001   pretty much copied from putfloat above */
{
	void* Pchars=&value;
	char* byte=(char*)Pchars;
	if (byteorder==bigendian) for (int i=0;i<2;i++) (*fout).put(byte[i]);
	else for (int i=1;i>=0;i--) (*fout).put(byte[i]);
}

cube cube::binary2cube(const char infilename[], int x, int y, int z, axis A, axis B,
		axis C)
/* for loading starcrash densfiles, though should work for any float binary
file.  Not yet implimented:  ints, doubles, chars.  JB 2004.9.20

Usage:  x, y, and z are the sizes of the cube dimensions, while A, B, and C
are the axes that change fastest (A) next fastest (B) and least fast (C) in the
file layout.  */
{
	FILE *infile;
	infile = fopen(infilename, "r");
	
	cube answer(x,y,z);
	for (long c=0;c<answer.N(C);c++)
		for (long b=0;b<answer.N(B);b++)
			for (long a=0;a<answer.N(A);a++)
				answer(A, a, B, b, C, c) = float(getw(infile));

	return answer;
}
