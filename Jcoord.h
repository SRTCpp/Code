class coord {
	public:
		coord(int,int,int);
		coord();
		~coord();
		
		int& A(axis a);
		
	private:
		int x, y, z;
}
