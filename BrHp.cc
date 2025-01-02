class MyAnaBrHp
{
  public:
  
  double br_Hp_ud = 0, br_Hp_cs = 0, br_Hp_tb = 0, br_Hp_taunu = 0, br_Hp_WZvirt = 0, br_Hp_WvirtZ = 0, br_Hp_WZ = 0, br_Hp_hW = 0, br_Hp_cascade = 0;
  
  MyAnaBrHp(double vd, double vT, double mh, double mH, double mHp, double tanAlpha);
};

MyAnaBrHp::MyAnaBrHp(double vd, double vT, double mh, double mH, double mHp, double tanAlpha)
{

 double re = pow(me/mHp,2), rmu = pow(mmu/mHp,2), rtau = pow(mtau/mHp,2), rW = pow(mW/mHp,2), rZ = pow(mZ/mHp,2), rh = pow(mh/mHp,2), rH = pow(mH/mHp,2);
 double ru = pow(mu/mHp,2), rc = pow(mc/mHp,2), rt = pow(mt/mHp,2), rd = pow(md/mHp,2), rs = pow(ms/mHp,2), rb = pow(mb/mHp,2);
 
 double cosb = vd/sqrt(vd*vd+4*vT*vT), sinb = -2*vT/(sqrt(vd*vd+4*vT*vT)), cosa = 1./sqrt(1+pow(tanAlpha,2)), sina = tanAlpha*cosa;
 
 //Lagrangian parameters in terms of physical masses and mixing angle
 double lambda1 = (mh*mh*cosa*cosa + mH*mH*sina*sina)/(2*vd*vd);
 double lambda2 = (mh*mh*sina*sina + mH*mH*cosa*cosa - vd*vd*mHp*mHp/(vd*vd+4*vT*vT))/(2*vT*vT);
 double lambda3 = (mh*mh-mH*mH)*sina*cosa/(vd*vT) + 2*mHp*mHp/(vd*vd+4*vT*vT);
 double lambda4 = 2*vT*mHp*mHp/(vd*vd+4*vT*vT);
 
 //into a pair of leptons
 double DW_Hp_enu = 0, DW_Hp_munu = 0, DW_Hp_taunu = 0;
 DW_Hp_enu = pow(mHp,3)*sinb*sinb*betaff(re,0)/(8*M_PI*vd*vd);
 DW_Hp_munu = pow(mHp,3)*sinb*sinb*betaff(rmu,0)/(8*M_PI*vd*vd);
 DW_Hp_taunu = pow(mHp,3)*sinb*sinb*betaff(rtau,0)/(8*M_PI*vd*vd);
 
 double DW_Hp_lnu = DW_Hp_enu + DW_Hp_munu + DW_Hp_taunu;
  
 //into a pair of quarks
 double DW_Hp_ud = 0, DW_Hp_us = 0, DW_Hp_ub = 0, DW_Hp_cd = 0, DW_Hp_cs = 0, DW_Hp_cb = 0, DW_Hp_td = 0, DW_Hp_ts = 0, DW_Hp_tb = 0, DW_Hp_bbW = 0;
 double kQCD_cc = 0.3158, kQCD_bb = 0.5827, kQCD_ss = 0.49398;
 DW_Hp_ud = Vud*Vud*3*pow(mHp,3)*sinb*sinb*betaff(ru,rd)/(8*M_PI*vd*vd);
 DW_Hp_us = Vus*Vus*3*pow(mHp,3)*sinb*sinb*betaff(ru,rs)/(8*M_PI*vd*vd);
 DW_Hp_ub = Vub*Vub*3*pow(mHp,3)*sinb*sinb*betaff(ru,rb)/(8*M_PI*vd*vd);
 
 DW_Hp_cd = Vcd*Vcd*3*pow(mHp,3)*sinb*sinb*betaff(rc,rd)/(8*M_PI*vd*vd)*kQCD_cc; //multiply with kQCD_cc to account for NLO effects
 DW_Hp_cs = Vcs*Vcs*3*pow(mHp,3)*sinb*sinb*betaff(rc,rs)/(8*M_PI*vd*vd)*kQCD_cc;
 DW_Hp_cb = Vcb*Vcb*3*pow(mHp,3)*sinb*sinb*betaff(rc,rb)/(8*M_PI*vd*vd)*kQCD_bb;
 
 if(mHp > mt+md) DW_Hp_td = Vtd*Vtd*3*pow(mHp,3)*sinb*sinb*betaff(rt,rd)/(8*M_PI*vd*vd);
 if(mHp > mt+ms) DW_Hp_ts = Vts*Vts*3*pow(mHp,3)*sinb*sinb*betaff(rt,rs)/(8*M_PI*vd*vd);
 if(mHp > mt+mb) DW_Hp_tb = Vtb*Vtb*3*pow(mHp,3)*sinb*sinb*betaff(rt,rb)/(8*M_PI*vd*vd);
 if(mHp > mW+2*mb && mHp < mt+mb) DW_Hp_bbW = 3*pow(mt,4)*mHp*sinb*sinb*betat(rW,rt)/(128*pow(M_PI,3)*pow(vd,4));
 
 double DW_Hp_tb_tot = DW_Hp_tb + DW_Hp_bbW;

/***********************************/

 //  if(vT == 4.6){
 //  if(mHp == 170) DW_Hp_tb_tot = 0.0000086;
 //  if(mHp == 171) DW_Hp_tb_tot = 0.0000095;
 //  if(mHp == 172) DW_Hp_tb_tot = 0.0000107;
 //  if(mHp == 173) DW_Hp_tb_tot = 0.0000125;
 //  if(mHp == 174) DW_Hp_tb_tot = 0.0000151;
 //  if(mHp == 175) DW_Hp_tb_tot = 0.0000189;
 //  if(mHp == 176) DW_Hp_tb_tot = 0.0000243;
 //  if(mHp == 177) DW_Hp_tb_tot = 0.0000320;
 //  if(mHp == 178) DW_Hp_tb_tot = 0.0000424;
 //  if(mHp == 179) DW_Hp_tb_tot = 0.0000561;}
//  if (vT == 2.9)
//  {
//    if (mHp == 171)DW_Hp_tb_tot = 0.0000041;
//    if (mHp == 172)DW_Hp_tb_tot = 0.0000044;
//    if (mHp == 173)DW_Hp_tb_tot = 0.0000049;
//    if (mHp == 174)DW_Hp_tb_tot = 0.0000056;
//    if (mHp == 175)DW_Hp_tb_tot = 0.0000067;
//    if (mHp == 176)DW_Hp_tb_tot = 0.0000086;
//    if (mHp == 177)DW_Hp_tb_tot = 0.0000115;
//    if (mHp == 178)DW_Hp_tb_tot = 0.0000157;
//    if (mHp == 179)DW_Hp_tb_tot = 0.0000216;
//  }

/***********************************/
 double XvT = 4.6;
 double Xsinb = -2 * XvT / (sqrt(vd * vd + 4 * XvT * XvT));
 double Xcosb = vd / sqrt(vd * vd + 4 * XvT * XvT);
 double XlamHWZ = -g / 2. * (2 * g * XvT * cosW * Xcosb - g1 * vd * sinW * Xsinb);

 if (vT == 4.1)
 {if (mHp == 170)DW_Hp_tb_tot = 0.0000086*pow(sinb/Xsinb,2);
  if (mHp == 171)DW_Hp_tb_tot = 0.0000095*pow(sinb/Xsinb,2);
  if (mHp == 172)DW_Hp_tb_tot = 0.0000107*pow(sinb/Xsinb,2);
  if (mHp == 173)DW_Hp_tb_tot = 0.0000125*pow(sinb/Xsinb,2);
  if (mHp == 174)DW_Hp_tb_tot = 0.0000151*pow(sinb/Xsinb,2);
  if (mHp == 175)DW_Hp_tb_tot = 0.0000189*pow(sinb/Xsinb,2);
  if (mHp == 176)DW_Hp_tb_tot = 0.0000243*pow(sinb/Xsinb,2);
  if (mHp == 177)DW_Hp_tb_tot = 0.0000320*pow(sinb/Xsinb,2);
  if (mHp == 178)DW_Hp_tb_tot = 0.0000424*pow(sinb/Xsinb,2);
  if (mHp == 179)DW_Hp_tb_tot = 0.0000561*pow(sinb/Xsinb,2);}




 
 double DW_Hp_qq = DW_Hp_ud + DW_Hp_us + DW_Hp_ub + DW_Hp_cd + DW_Hp_cs + DW_Hp_cb + DW_Hp_td + DW_Hp_ts + DW_Hp_tb_tot;
 
 double DW_Hp_ff = DW_Hp_qq + DW_Hp_lnu;
 
 //into a pair of gauge bosons (WZ and Wa)
 double DW_Hp_WZ = 0, DW_Hp_WZvirt = 0, DW_Hp_WvirtZ = 0, DW_Hp_Wa = 0;
 double lamHWZ = -g/2.*(2*g*vT*cosW*cosb-g1*vd*sinW*sinb);
 double lamHWa = -g/2.*(2*vT*cosb+vd*sinb);
 if(mHp > mW+mZ) DW_Hp_WZ = pow(lamHWZ,2)*sqrt(beta(rW,rZ))*(2+pow(1-rW-rZ,2)/(4*rW*rZ))/(16*M_PI*mHp);
 if(mHp > mW && mHp < mW+mZ) DW_Hp_WZvirt = 9*gZ*gZ*pow(lamHWZ,2)*H(rW,rZ)*(7./12.-10./9.*pow(sinW,2)+40./27.*pow(sinW,4))/(128*pow(M_PI,3)*mHp);
 if(mHp > mZ && mHp < mW+mZ) DW_Hp_WvirtZ = 9*g*g*pow(lamHWZ,2)*H(rZ,rW)/(256*pow(M_PI,3)*mHp);
 if(mHp > mW) DW_Hp_Wa = 0; //pow(lamHWa,2)*(1-rW)/(8*M_PI*mHp);
  
 double DW_Hp_Wnunu = DW_Hp_WZvirt*0.205, DW_Hp_Wee = DW_Hp_WZvirt*0.03363, DW_Hp_Wmumu = DW_Hp_WZvirt*0.03366, DW_Hp_Wtautau = DW_Hp_WZvirt*0.03367;
 double DW_Hp_Wdd = DW_Hp_WZvirt*0.156, DW_Hp_Wss = DW_Hp_Wdd, DW_Hp_Wbb = DW_Hp_Wdd, DW_Hp_Wuu = DW_Hp_WZvirt*0.116, DW_Hp_Wcc = DW_Hp_Wuu;
 double DW_Hp_Zenu = DW_Hp_WvirtZ/9, DW_Hp_Zmunu = DW_Hp_WvirtZ/9, DW_Hp_Ztaunu = DW_Hp_WvirtZ/9;
 double DW_Hp_Zud = DW_Hp_WvirtZ/3*Vud*Vud, DW_Hp_Zus = DW_Hp_WvirtZ/3*Vus*Vus, DW_Hp_Zub = DW_Hp_WvirtZ/3*Vub*Vub;
 double DW_Hp_Zcd = DW_Hp_WvirtZ/3*Vcd*Vcd, DW_Hp_Zcs = DW_Hp_WvirtZ/3*Vcs*Vcs, DW_Hp_Zcb = DW_Hp_WvirtZ/3*Vcb*Vcb;
  
 double DW_Hp_WZ_tot = DW_Hp_WZ + DW_Hp_WZvirt + DW_Hp_WvirtZ;
 /*******************************/
//  if(vT == 4.6){
//  if(mHp == 171) DW_Hp_WZ_tot = 0.000128;
//  if(mHp == 172) DW_Hp_WZ_tot = 0.000163;
//  if(mHp == 173) DW_Hp_WZ_tot = 0.000206;}
 
//  if(vT == 2.9) {
//  if(mHp == 171) DW_Hp_WZ_tot = 0.000051;
//  if(mHp == 172) DW_Hp_WZ_tot = 0.000065;
//  if(mHp == 173) DW_Hp_WZ_tot = 0.000082;}
 /*******************************/
  if(vT == 4.1){
  if(mHp == 171) DW_Hp_WZ_tot = 0.000128*pow(lamHWZ/XlamHWZ,2);
  if(mHp == 172) DW_Hp_WZ_tot = 0.000163*pow(lamHWZ/XlamHWZ,2);
  if(mHp == 173) DW_Hp_WZ_tot = 0.000206*pow(lamHWZ/XlamHWZ,2);}

 double DW_Hp_VV = DW_Hp_WZ_tot + DW_Hp_Wa;
 
 //into the SM Higgs and W-boson
 double DW_Hp_hW = 0, DW_Hp_hWvirt = 0;
 double lamhHpW = -g/2.*(2*sina*cosb+cosa*sinb);
 if(mHp > mh+mW) DW_Hp_hW = pow(mHp,3)*pow(lamhHpW,2)*pow(beta(rW,rh),1.5)/(16*M_PI*mW*mW);
 if(mHp > mh && mHp < mh+mW) DW_Hp_hWvirt = 9*pow(g,2)*mHp*pow(lamhHpW,2)*G(rh,rW)/(128*pow(M_PI,3));
 double DW_Hp_henu = DW_Hp_hWvirt/9, DW_Hp_hmunu = DW_Hp_hWvirt/9, DW_Hp_htaunu = DW_Hp_hWvirt/9;
 double DW_Hp_hud = DW_Hp_hWvirt/3*Vud*Vud, DW_Hp_hus = DW_Hp_hWvirt/3*Vus*Vus, DW_Hp_hub = DW_Hp_hWvirt/3*Vub*Vub;
 double DW_Hp_hcd = DW_Hp_hWvirt/3*Vcd*Vcd, DW_Hp_hcs = DW_Hp_hWvirt/3*Vcs*Vcs, DW_Hp_hcb = DW_Hp_hWvirt/3*Vcb*Vcb;
 
 double DW_Hp_hW_tot = DW_Hp_hW + DW_Hp_hWvirt;
 
 //into a lighter scalar and a gauge boson (cascade decay)
 double DW_Hp_HW = 0, DW_Hp_HWvirt = 0;
 double lamHHpW = sina*sinb-2*cosa*cosb;
 if(mHp > mH+mW) DW_Hp_HW = g*g*pow(mHp,3)*pow(lamHHpW,2)*pow(beta(rW,rH),1.5)/(64*M_PI*mW*mW);
 if(mHp > mH && mHp < mH+mW) DW_Hp_HWvirt = 9*pow(g,4)*mHp*pow(lamHHpW,2)*G(rH,rW)/(512*pow(M_PI,3));
 double DW_Hp_Henu = DW_Hp_HWvirt/9, DW_Hp_Hmunu = DW_Hp_HWvirt/9, DW_Hp_Htaunu = DW_Hp_HWvirt/9;
 double DW_Hp_Hud = DW_Hp_HWvirt/3*Vud*Vud, DW_Hp_Hus = DW_Hp_HWvirt/3*Vus*Vus, DW_Hp_Hub = DW_Hp_HWvirt/3*Vub*Vub;
 double DW_Hp_Hcd = DW_Hp_HWvirt/3*Vcd*Vcd, DW_Hp_Hcs = DW_Hp_HWvirt/3*Vcs*Vcs, DW_Hp_Hcb = DW_Hp_HWvirt/3*Vcb*Vcb;
 
 double DW_Hp_cascade = DW_Hp_HW + DW_Hp_HWvirt;
  
 //total decay width
 double DW_Hp_tot = DW_Hp_ff + DW_Hp_VV + DW_Hp_hW_tot + DW_Hp_cascade;
 
 //Branching ratios
 br_Hp_cs = DW_Hp_cs/DW_Hp_tot, br_Hp_tb = DW_Hp_tb_tot/DW_Hp_tot, br_Hp_taunu = DW_Hp_taunu/DW_Hp_tot;
 br_Hp_WZvirt = DW_Hp_WZvirt/DW_Hp_tot, br_Hp_WvirtZ = DW_Hp_WvirtZ/DW_Hp_tot, br_Hp_WZ = DW_Hp_WZ_tot/DW_Hp_tot, br_Hp_hW = DW_Hp_hW_tot/DW_Hp_tot, br_Hp_cascade = DW_Hp_cascade/DW_Hp_tot;


}
