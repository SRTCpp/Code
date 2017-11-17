#include "Jcube.h"
#include <string.h>

main(int argc, char **argv)
{
	if (argc < 2) cout << "No name\n";
	cube tograph(argv[1]);
	
	char switches[400]="";
	for (int i=2;i<argc;i++) strcat(switches, argv[i]);
	tograph.graph(switches);
}
