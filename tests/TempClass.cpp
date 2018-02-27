#include <iostream>
using namespace std;

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

int main(){

	Point M;
	M.setPoint(1.4,1.4);
	M.printPoint();
	return 0;
}
