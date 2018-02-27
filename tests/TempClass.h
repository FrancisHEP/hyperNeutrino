class Point
{
	//accessible
	public:
		void setPoint(double x, double y);
		void printPoint();

	//non-accessible
	private:
		double xPos;
		double yPos;
};

void Point::setPoint(double x, double y)
{
	xPos = x;
	yPos = y;
}

void Point::printPoint()
{
	cout<<"x = "<<xPos<<", "<<"y = "<<yPos<<endl;
}
