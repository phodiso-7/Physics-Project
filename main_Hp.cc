#include <cmath>
#include <iostream>
#include <complex>
#include <fstream>
#include <iomanip>

using namespace std;

#include "func.cc"
#include "BrHp.cc"

int main()
{

  dcomp f(double x);
  dcomp j(double x);
  dcomp I1(double x, double y);
  dcomp I2(double x, double y);
  dcomp betaH0(double x);
  dcomp betaH12(double x);
  dcomp betaH1(double x);
  dcomp betaHZ0(double x, double y);
  dcomp betaHZ12(double x, double y);
  dcomp betaHZ1(double x, double y);

  double mh = 125, mH, mHp, dm = 0, alpha = 0, tanAlpha = 0, vSM = 246.22, vd = vSM, vT = 4.1;
  ofstream out1;
  out1.open("BrHpm_mass.dat");
  for (double j = 100; j <= 200; j = j + 1.0)
  {

    mH = j;
    mHp = mH + dm;

    MyAnaBrHp myAnaBrHp(vd, vT, mh, mH, mHp, tanAlpha);

    out1 << mHp << " " << myAnaBrHp.br_Hp_cs << " " << myAnaBrHp.br_Hp_tb << " " << myAnaBrHp.br_Hp_taunu << " " << myAnaBrHp.br_Hp_WZ << endl;
  }
  out1.close();
}
