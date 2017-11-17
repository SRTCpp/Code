#include "../Jcrap.h"
#include "../Jangle.h"
#include "photontraverse.h"
#include "photongenerators.h"
#include "detectors.h"

#include <list>

class SRTC
{
	public:
		SRTC(atmosphere, photongenerator*, vector<detector*>);
	
		atmosphere& Atmosphere() {return _atmosphere;}
		vector<detector*>& Detectors() {return _detectorlist;}
		
		void run(int=-1);
		
		photongenerator* tophotongenerator() {return _tophotongenerator;}
		
		static int debug, ompdebug;
		double numericalprecision;
		
		static string datapath();
			
	private:
		atmosphere _atmosphere;
			
		unsigned long _iterationlimit;
	
		photongenerator *_tophotongenerator;
		
		vector<detector*> _detectorlist;
};
