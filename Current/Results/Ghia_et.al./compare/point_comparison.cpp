#include <iostream>
#include <fstream>
#include <cmath> 
#include <assert.h> 
#include <string>
#include <sstream>



using namespace std;



double osUy0100[129][2]; 
double ghUy0100[17][2];
//{1,8,9,10,14,23,37,59,65,80,95,110,123,124,125,126,129};
int indexY[17] = 
{1,9,10,11,13,21,30,31,65,104,111,117,122,123,124,125,129};
int main() 
{ 
  ifstream infile;   

  int num = 0; // num must start at 0
  infile.open("Suleman-vx-Re0100-0129x0129.dat");// file containing numbers in 3 columns 
     if(infile.fail()) // checks to see if file opended 
    { 
      cout << "error" << endl; 
      return 1; // no point continuing if the file didn't open...
    } 
       while(!infile.eof()) // reads file to end of *file*, not line
      { 
         infile >> osUy0100[num][1]; // read first column number
         infile >> osUy0100[num][2]; // read second column number
	

         ++num; // go to the next number

         // you can also do it on the same line like this:
         // infile >> exam1[num] >> exam2[num] >> exam3[num]; ++num;
      } 
  infile.close(); 


/*
  infile.open("GhiaEtAl-uy-Re100.dat");// file containing numbers in 3 columns 
     if(infile.fail()) // checks to see if file opended 
    { 
      cout << "error" << endl; 
      return 1; // no point continuing if the file didn't open...
    } 
for (int i = 0; i <=16; i++)
{
         infile >> ghUy0100[num][1]; // read first column number
         infile >> ghUy0100[num][2]; // read second column number
 cout << ghUy0100[num][1] << ghUy0100[num][2] << endl;

         ++num; // go to the next number

         // you can also do it on the same line like this:
         // infile >> exam1[num] >> exam2[num] >> exam3[num]; ++num;
      } 
  infile.close(); 
*/



ofstream myfile;
std:string filename = "comp-vx-100.dat";
double deltaU;
double deltaY;
for (int i = 0; i <=16; i++)
{
deltaY = osUy0100[indexY[i]-1][1];// - ghUy0100[i][1];
deltaU = osUy0100[indexY[i]-1][2];

   myfile.open (filename.c_str(),ios::app);
   myfile << deltaY << " " << deltaU << endl;
   myfile.close();

}


  return 0; // everything went right.
} 
