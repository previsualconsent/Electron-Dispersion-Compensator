// Runge-Kutta.cpp
//--------------------------------------------------
//A Runge-Kutta Method for solving Differential Equations
//of the form y'=f(x,y) ; y(x0)=y0
//--------------------------------------------------

#include <iostream>
#include <iomanip>
#include "gnuplot_i.hpp"
using namespace std;

//Define constants
#define X0 0 
#define Y0 0
#define H 0.2
#define N 5

//Define Functions
double f(double x, double y);
double runge(double x, double y);

void wait_for_key ()
{
#if defined(WIN32) || defined(_WIN32) || defined(__WIN32__) || defined(__TOS_WIN__)  // every keypress registered, also arrow keys
    cout << endl << "Press any key to continue..." << endl;

    FlushConsoleInputBuffer(GetStdHandle(STD_INPUT_HANDLE));
    _getch();
#elif defined(unix) || defined(__unix) || defined(__unix__) || defined(__APPLE__)
    cout << endl << "Press ENTER to continue..." << endl;

    std::cin.clear();
    std::cin.ignore(std::cin.rdbuf()->in_avail());
    std::cin.get();
#endif
    return;
}

//Main Function
int main()
{
double x;
double y;
vector<double> xlist,ylist;
Gnuplot g1;

cout<<"\t*** Euler Method ***"
<<"\n\n";
cout<<" "
<<setw(12)<<"x"<<setw(12)<<"\ty"
<<"\n"
<<"\t------------------------------"
<<"\n";
y=Y0;
for(int i=0;i<=5;i++)
{
xlist.push_back(x);
ylist.push_back(y);
x=X0+(i*H);
y=runge(x,y);

cout<<left<<setw(6)<<i<<"|"
<<setprecision(4)<<left<<setw(8)<<"\t"<<x
<<setprecision(4)<<left<<setw(8)<<"\t"<<y;
cout<<"\n\n";
}
g1.set_style("points").plot_xy(xlist,ylist,"plot");
wait_for_key();

cout<<"\n\n";
return 0;
}

double runge(double x[], double y)
{
double K1 = (H * f(x,y));
double K2 = (H * f((x + 1 / 2 * H), (y + 1 / 2 * K1)));
double K3 = (H * f((x + 1 / 2 * H), (y + 1 / 2 * K2)));
double K4 = (H * f((x + H), (y + K3)));
double runge = (y + ((K1 + 2 * K2 + 2 * K3 + K4)/6));
return runge;
}

double f(double x, double y)
{
double f = x+y;
return f;
}

