typedef complex<double> dcomp;

double vSM = 246.22, mP = 0., mG = 0., me = 0.00051099895, mmu = 0.1056583755, mtau = 1.77686, md = 0.00467, mu = 0.00216, ms = 0.0934, mc = 1.27, mb = 4.18, mt = 172.69, mZ = 91.1876;
double Gf = 0.00001166378, g = 0.65742, mW = sqrt(g*g/(4.*Gf*sqrt(2.))), sinW = 0.460682, cosW = sqrt(1-sinW*sinW), tanW = sinW/cosW, gZ = g/cosW, g1 = gZ*sinW;
double Vud = 0.97427, Vus = 0.22534, Vub = 0.00351, Vcd = 0.2252, Vcs = 0.97344, Vcb = 0.0412, Vtd = 0.00867, Vts = 0.0404, Vtb = 0.999146; 
double as = 0.1179, aem = 1./128., aem0 = 1./137.;

dcomp f(double x) 
{
 if (x < 1) return -pow(log(1+sqrt(1-x))-log(1-sqrt(1-x))-dcomp(0,1)*M_PI,2)/4.;
 else return pow(asin(sqrt(1./x)),2);
}

dcomp j(double x) 
{
 if (x < 1) return 0.5*sqrt(1-x)*(log(1+sqrt(1-x))-log(1-sqrt(1-x))-dcomp(0,1)*M_PI);
 else return sqrt(x-1)*asin(sqrt(1./x));
}

dcomp betaH0(double x) 
{
 return -x*(1.-x*f(x));
}

dcomp betaH12(double x) 
{
 return 2*x*(1.+(1-x)*f(x));
}

dcomp betaH1(double x) 
{
 return -(2+3*x+3*x*(2-x)*f(x));
}

dcomp I1(double x, double y)
{
 return 0.5*x*y/(x-y) + 0.5*pow(x*y/(x-y),2)*(f(x)-f(y)) + pow(x/(x-y),2)*y*(j(x)-j(y));
}

dcomp I2(double x, double y)
{
 return -0.5*x*y/(x-y)*(f(x)-f(y));
}

dcomp betaHZ0(double x, double y) 
{
 return I1(x,y);
}

dcomp betaHZ12(double x, double y) 
{
 return I1(x,y)-I2(x,y);
}

dcomp betaHZ1(double x, double y) 
{
 double sinW = 0.460682, cosW = sqrt(1-sinW*sinW), tanW = sinW/cosW;
 return cosW*((1+2./x)*pow(tanW,2)-(5+2./x))*I1(x,y) + 4*cosW*(3-tanW*tanW)*I2(x,y);
}

double beta(double x, double y)
{
 return 1+x*x+y*y-2*x*y-2*x-2*y;
}

double betaf(double x) 
{
 return sqrt(1-4*x);
}

double betaV(double x) 
{
 return (1-4*x+12*x*x)*betaf(x);
}

double betaVp(double x) 
{
 return 3*(1-8*x+20*x*x)/sqrt(4*x-1)*acos((3*x-1)/(2*x*sqrt(x))) -(1-x)*(2-13*x+47*x*x)/(2.*x) -1.5*(1-6*x+4*x*x)*log(x);
}

double betaS(double x) 
{
 return (x-1)*(2-0.5*log(x)) + (1-5*x)/sqrt(4*x-1)*(atan((2*x-1)/sqrt(4*x-1))-atan(1/sqrt(4*x-1)));
}

double G(double x, double y) 
{
 return (2*pow((-1+x),3)-9*(-1+x*x)*y+6*(-1+x)*y*y+6*(1+x-y)*y*sqrt(-beta(x,y))*(atan((-1+x-y)/(sqrt(-beta(x,y))))+atan((-1+x+y)/(sqrt(-beta(x,y)))))-3*(1+x*x+y*y-2*x*y-2*y)*y*log(x))/(12*y);
}

double betaff(double x, double y)
{
 return ((x+y)*(1-x-y)-4*x*y)*sqrt(beta(x,y));
}

double betat(double x, double y)
{
 return pow(x,2)/pow(y,3)*(4*x*y+3*y-4*x)*log(x*(y-1)/(y-x)) + (3*y*y-4*y-3*x*x+1)*log((y-1)/(y-x)) -5./2. + (1-x)/pow(y,2)*(3*pow(y,3)-y*x-2*y*x*x+4*x*x) +x*(4-3*x/2.);
}

double H(double x, double y) 
{
 return (atan((1-x+y)/sqrt(-beta(x,y)))+atan((1-x-y)/sqrt(-beta(x,y))))*(-3*x*x*x+(9*y+7)*x*x-5*(1-y)*(1-y)*x+pow((1-y),3))/(4*x*sqrt(-beta(x,y)))+((-1+x)*(6*y*y+y*(39*x-9)+2*(1-x)*(1-x))-3*y*(y*y+2*y*(3*x-1)-x*(3*x+4)+1)*log(x))/(24*x*y);
}
