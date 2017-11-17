
class Jcrap::pixelnode
{
	public:
		pixelnode(cube&, int, int, methodtype=mFAST);
		pixelnode(pixelnode*, pixelnode*);
		~pixelnode();
		
		bool operator==(pixelnode&);
		bool operator<(pixelnode&);
		
		void assignparent(pixelnode*);
		int children();
		cube& centroid();
		int x();
		int y();
		
		pixelnode& left();
		pixelnode& right();
		
		pair<double, pixelnode*> closest;
		
		double distance(pixelnode&);
		bool unispectral();
		void convert();
		
		double sum;
		
		static cube return1;
		static cube return2;
		static int returncount;
			
	private:
		pixelnode* Parent;
		pixelnode* Left;
		pixelnode* Right;
		cube *Centroid;
		int Children;
		methodtype method;
		
		int xcoord;
		int ycoord;
};
