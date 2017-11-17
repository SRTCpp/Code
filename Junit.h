
class Jcrap::unit : public list<Jcrap::unitelement>
{	
	public:

//		Junit(const Junit&);     // this is unnecessary -- the default one should be fine
		unit(const char* = "");
		unit(double);
		unit(unitelement);
//		~unit();
		
// subdefinition
		enum builtintype{base,dimensionless,system,conglomerate,null,unknown};

// operators!
		friend ostream& operator<<(ostream&, const unit);
		bool operator==(unit);
		unit operator*=(unit);
		unit operator*(unit);
		unit operator/(unit);
		operator double();
		
		
// functions
		unit pow(double);
		unit simplify(void);
		unit simplify(int);
		unit simplify(unit);
		unit convert(unit);
//		unit simplify(const char*);
		unit removedimensionless();
		unit reorder();
		unit diff(unit);
		double order();
		double expadd();
		static void setdefault(int=0);
		builtintype stype();
		friend class unitelement;
		string tostring() const;
		
		
	private:
		static int defaultsimplifyint;
};


class Jcrap::unitelement
{
	public:
		// in addition to default Junitelement(Junitelement) creator
		unitelement(string);	
		unitelement();	
		unitelement(double);
	
		string abbr;		// abbreviation of unit name, "cm"
		string name;  		// name of unit, "centimeters"
		float exp;			// exponent, as 2 for cm^2
		float mult;			// how many of the base units are in one of this unit
								// i.e. 1e-2 for cm
		string type;		// i.e. "length" for cm

		unit subunits;   //  i.e., this has (kg), (m), (s^-2) for Newtons
		
		bool operator==(const unitelement&);
		bool operator< (unitelement r);
		unit operator++(int);				// moves to next in the unitdatabase
		
		const unit::builtintype& stype() { return specialtype; }

	private:	
		unit::builtintype specialtype;			// so I can quit testing unitelement.type
			
		// taking care of database
		static unit unitdatabase;
		static void readunitfile();
		
};
				

ostream& operator<<(ostream& out, const Jcrap::unitelement& j);

ostream& operator<<(ostream& out, const Jcrap::unit j);
