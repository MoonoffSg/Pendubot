
nargout_loss_cp_t loss_pendubot(cost_t cost, double m[], double s[][])
{
  nargout_loss_cp_t ret      = new nargout_loss_cp_t();
  gTrig_r           tg       = new gTrig_r();
  gTrig_Full_r      tgF      = new gTrig_Full_r();

  double   cw                = cost.width_c;
  double   b                 = cost.expl;
  
  double   target[];// = new double [D1];
  double   sV[][];//   = new double [s.length][s[0].length];
  double   mdm[][];//  = new double [2][D0];
  double   sdm[][];//  = new double [D0][D0];
  double   Cdm[][];//  = new double [D1 + 2][D0];
  double   mds[][];//  = new double [2][D0 * D0];
  double   sds[][];//  = new double [D0][D0 * D0];
  double   Cds[][];//  = new double [D1 + 2][D0 * D0];
  double   dCdm[][];
  double   dCds[][];// = new double [D1][D1];

  // 1. Some precomputations
  int      D0 = s[0].length,   // state dimension 
           D = D0,
           D1 = D0 + 2 * cost.angle.length;        // state dimension (with sin/cos)
  
  
  double   M[]      = new double [D1];
  for(int i = 0; i < D0 ; i++)
    M[i] = m[i];
    
  double   S[][] = new double [D1][D1];
  for(int i = 0; i < D0 ; i++)
    for(int j = 0 ; j < D0 ; j++)
      S[i][j] = s[i][j];
  
  double   Mdm[][] = new double[D0 + D1 - D0][D0];
  for(int i = 0; i < D0; i++)
    Mdm[i][i] = 1;
  
  double  Sdm[][] = new double [D1 * D1][D0];
  double  Mds[][] = new double [D1][D0 * D0];
  double  Sds[][];// = new double [D1 * D1][D0 * D0];
  Sds = kron(Mdm, Mdm);

  //% 2. Define static penalty as distance from target setpoint
  //ell1 = cost.p(1); ell2 = cost.p(2); C = [ell1 0 ell2 0; 0 ell1 0 ell2];
  //Q = zeros(D1); Q(D+1:D+4,D+1:D+4) = C'*C;
  
  double ell1 = cost.p[0],
         ell2 = cost.p[1];
         
  double C[][] = {
                   {ell1, 0, ell2, 0},
                   {0, ell1, 0, ell2}
                  };
  int indR[] = {D, D + 1, D + 2, D + 3},
      indC[] = {D, D + 1, D + 2, D + 3};
  double Q[][] = new double[D1][D1];
  Q = setMat(Q, indR, indC, matrixMultiply(transMat(C), C));

  // 3. Trigonometric augmentation
  
//  if ((D1-D0) > 0)
  {
    // augment target
    //target = [cost.target(:); gTrig(cost.target(:), 0*s, cost.angle)]
 
    target = new double[cost.target.length + 4];
    sV     = new double[s.length][s[0].length];
    for(int i = 0; i < cost.target.length ; i++)
      target[i] = cost.target[i];
      
    tg = gTrig(cost.target, sV, cost.angle, 1);
    
    target[cost.target.length]      = tg.M[0];
    target[cost.target.length + 1]  = tg.M[1];
    target[cost.target.length + 2]  = tg.M[2];
    target[cost.target.length + 3]  = tg.M[3];
    
    // augment state
    //i = 1:D0; k = D0+1:D1
    int i[] = createInd(0, D0),
        k[] = createInd(D0, D1);
      
    tgF = gTrigF(getVec(M, i), getMat(S, i, i), cost.angle, 1);
    
    M = setVec(M, k, tgF.M);
    S = setMat(S, k, k, tgF.V);
    C = tgF.C;//coeffProdMat(1.0, tgF.C); 
    mdm = tgF.dMdm;
    sdm = tgF.dVdm;
    Cdm = tgF.dCdm;
    mds = tgF.dMdv;
    sds = tgF.dVdv;
    Cds = tgF.dCdv;

    // compute derivatives (for augmentation)
   
    int X[][]  = new int [D1][D1];
    int XT[][] = new int [D1][D1];
//    float I[][] = new float [D1][D1];
    
    for(int i1 = 0 ; i1 < D1 ; i1++)
    for(int j1 = 0 ; j1 < D1 ; j1++)
      X[j1][i1] = j1 + D1 * i1; 
    
    XT = transMat(X);
    
    int ii[]  = new int [D0 * D0];
    for(int i1 = 0 ; i1 < D0 ; i1++)
    for(int j1 = 0 ; j1 < D0 ; j1++)
      ii[j1 + D0 * i1] = X[j1][i1]; 
    
    int kk[]  = new int [(D1 - D0) * (D1 - D0)];
    int indexKK = 0;
    for(int i1 = D0 ; i1 < D1 ; i1++)
    for(int j1 = D0 ; j1 < D1 ; j1++)
      kk[indexKK++] = X[j1][i1]; 
      
    int ik[]  = new int [(D1 - D0) * D0];
    indexKK = 0;
    for(int i1 = D0 ; i1 < D1 ; i1++)
    for(int j1 = 0 ; j1 < D0 ; j1++)
      ik[indexKK++] = X[j1][i1]; 
      
    int ki[]  = new int [(D1 - D0) * D0];
    indexKK = 0;
    for(int i1 = D0 ; i1 < D1 ; i1++)
    for(int j1 = 0 ; j1 < D0 ; j1++)
      ki[indexKK++] = XT[j1][i1]; 
    
  // chainrule
  // Mdm(k,:)  = mdm*Mdm(i,:) + mds*Sdm(ii,:); 
    Mdm = setMat(Mdm, k, Chainrule(mdm, Mdm, i, mds, Sdm, ii));
     
  //Mds(k,:)  = mdm*Mds(i,:) + mds*Sds(ii,:);
    Mds = setMat(Mds, k, Chainrule(mdm, Mds, i, mds, Sds, ii));
   
  //Sdm(kk,:) = sdm*Mdm(i,:) + sds*Sdm(ii,:);
    Sdm = setMat(Sdm, kk, Chainrule(sdm, Mdm, i, sds, Sdm, ii));
  
  //Sds(kk,:) = sdm*Mds(i,:) + sds*Sds(ii,:);
    Sds = setMat(Sds, kk, Chainrule(sdm, Mds, i, sds, Sds, ii));
   
  //dCdm      = Cdm*Mdm(i,:) + Cds*Sdm(ii,:);
    dCdm = Chainrule(Cdm, Mdm, i, Cds, Sdm, ii);
 
  //dCds      = Cdm*Mds(i,:) + Cds*Sds(ii,:);
    dCds = Chainrule(Cdm, Mds, i, Cds, Sds, ii);
   
  //S(i,k) = S(i,i)*C; 
  //S(k,i) = S(i,k)';   
    S = addToSimpl(S, i, k, C);
  
  //SS = kron(eye(length(k)),S(i,i)); CC = kron(C',eye(length(i)));
    double SS[][] = kron(eye(k.length), getMat(S, i, i));
    double CC[][] = kron(transMat(C), eye(i.length));

  //Sdm(ik,:) = SS*dCdm + CC*Sdm(ii,:); 
  //Sdm(ki,:) = Sdm(ik,:);
    Sdm = setMat(Sdm, ik, Chainrule(SS, dCdm, CC, Sdm, ii));
    Sdm = setMat(Sdm, ki, getMat(Sdm, ik));
   
  //Sds(ik,:) = SS*dCds + CC*Sds(ii,:); 
  //Sds(ki,:) = Sds(ik,:);
    Sds = setMat(Sds, ik, Chainrule(SS, dCds, CC, Sds, ii));
    Sds = setMat(Sds, ki, getMat(Sds, ik)); 
  
  }
  
//  //4. Calculate loss!
  //L = 0; dLdm = zeros(1,D0); dLds = zeros(1,D0*D0); S2 = 0;
  double L         = 0,
         dLdm[]    = new double [D0],
         dLds[]    = new double [D0 * D0],
         S2        = 0;
  for(int i_o = 0 ; i_o < 1 ; i_o++) // for i = 1:length(cw) 
  {
    //cost.z = target; cost.W = Q/cw(i)^2;
    cost.z = target;//vectorCoef(target, 1);
    cost.W = coeffProdMat(1.0 / (cw * cw), Q);
   
// [r rdM rdS s2 s2dM s2dS] = lossSat(cost, M, S);
    nargout_lossSat_t dat = new nargout_lossSat_t();
    dat = lossSat(cost, M, S);
    double r = dat.L;
    double rdM[] = dat.dLdm;
    double rdS[][] = dat.dLds;
    double s2 = dat.S;
    double s2dM[] = dat.dSdm;
    double s2dS[][] = dat.dSds;
  
    L = L + r; 
    S2 = S2 + s2;
  
 // dLdm = dLdm + rdM(:)'*Mdm + rdS(:)'*Sdm;
    double rdST[] = new double[rdS.length * rdS[0].length];
    int indK = 0;
    for(int i1 = 0 ; i1 < rdS.length ; i1++)
      for(int j1 = 0 ; j1 < rdS[0].length ; j1++)
        rdST[indK++] = rdS[i1][j1];
    dLdm = vectorSUM(1, dLdm, 1, matrixMultiplyC_T2(rdM, Mdm));
    dLdm = vectorSUM(1, dLdm, 1, matrixMultiplyC_T2(rdST, Sdm));

  // dLds = dLds + rdM(:)'*Mds + rdS(:)'*Sds;
    dLds = vectorSUM(1, dLds, 1, matrixMultiplyC_T2(rdM, Mds));
    dLds = vectorSUM(1, dLds, 1, matrixMultiplyC_T2(rdST, Sds));

  //if (b~=0 || ~isempty(b)) && abs(s2)>1e-12
  
  if(/*(b != 0) && */(abs_db(s2) > 1e-12))
  {
    L = L + b * sqrt_db(s2);
 
    //dLdm = dLdm + b/sqrt(s2) * ( s2dM(:)'*Mdm + s2dS(:)'*Sdm )/2;
    double s2dST[] = new double[s2dS.length * s2dS[0].length];
    double dLdmT[] = new double[dLdm.length];
    double dLdsT[];
    indK = 0;
    for(int i = 0 ; i < rdS.length ; i++)
      for(int j = 0 ; j < rdS[0].length ; j++)
        s2dST[indK++] = s2dS[i][j];
    dLdmT = matrixMultiplyC_T2(s2dM, Mdm);
    dLdmT = vectorSUM(1, dLdmT, 1, matrixMultiplyC_T2(s2dST, Sdm));
    dLdmT = vectorCoef(dLdmT, b/sqrt_db(s2) * 0.5);
    dLdm = vectorSUM(1, dLdm, 1, dLdmT); 

    //dLds = dLds + b/sqrt(s2) * ( s2dM(:)'*Mds + s2dS(:)'*Sds )/2;
    dLdsT = matrixMultiplyC_T2(s2dM, Mds);
    dLdsT = vectorSUM(1, dLdsT, 1, matrixMultiplyC_T2(s2dST, Sds));
    dLdsT = vectorCoef(dLdsT, b/sqrt_db(s2) * 0.5);
    dLds = vectorSUM(1, dLds, 1, dLdsT); 
  }
  
  }
  
  // normalize n = 1
  // n = length(cw); L = L/n; dLdm = dLdm/n; dLds = dLds/n; S2 = S2/n
  ret.L = L;
  ret.S2 = S2;
  ret.dLdm = dLdm;//vectorCoef(dLdm, 1.0);
  ret.dLds = dLds;//vectorCoef(dLds, 1.0);
    
  return ret;
}

//nargout_loss_cp_t loss_pendulum(cost_t cost, double m[], double s[][])
//{
//  nargout_loss_cp_t ret      = new nargout_loss_cp_t();
//  gTrig_r           tg       = new gTrig_r();
//  gTrig_Full_r      tgF      = new gTrig_Full_r();

//  double   cw                = cost.width_c;
//  double   b                 = cost.expl;
  
//  double   target[];// = new double [D1];
//  double   sV[][];//   = new double [s.length][s[0].length];
//  double   C[][];//    = new double [D0][2];
//  double   mdm[][];//  = new double [2][D0];
//  double   sdm[][];//  = new double [D0][D0];
//  double   Cdm[][];//  = new double [D1 + 2][D0];
//  double   mds[][];//  = new double [2][D0 * D0];
//  double   sds[][];//  = new double [D0][D0 * D0];
//  double   Cds[][];//  = new double [D1 + 2][D0 * D0];
//  double   dCdm[][];
//  double   dCds[][];// = new double [D1][D1];

//  // 1. Some precomputations
//  int      D0 = s[0].length,   // state dimension 
//           D1 = D0 + 2;        // state dimension (with sin/cos)
  
  
//  double   M[]      = new double [D1];
//  for(int i = 0; i < D0 ; i++)
//    M[i] = m[i];
    
//  double   S[][] = new double [D1][D1];
//  for(int i = 0; i < D0 ; i++)
//    for(int j = 0 ; j < D0 ; j++)
//      S[i][j] = s[i][j];
  
//  double   Mdm[][] = new double[D0 + D1 - D0][D0];
//  for(int i = 0; i < D0; i++)
//    Mdm[i][i] = 1;
  
//  double  Sdm[][] = new double [D1 * D1][D0];
//  double  Mds[][] = new double [D1][D0 * D0];
//  double  Sds[][];// = new double [D1 * D1][D0 * D0];
//  Sds = kron(Mdm, Mdm);

//  // 2. Define static penalty as distance from target setpoint
  
//  double ell = cost.p;
//  double Q[][] = new double[D1][D1];
  
//  //double t[] = {1.0, ell};
//  //double tt[][] = Outer_Product(t, t);

//  Q[D0][D0]  = ell * ell;
//  Q[D0 + 1][D0 + 1] = ell * ell;

//  // 3. Trigonometric augmentation
  
////  if ((D1-D0) > 0)
//  {
//    // augment target
//    //target = [cost.target(:); gTrig(cost.target(:), 0*s, cost.angle)]
    
//    target = new double[cost.target.length + 2];
//    sV     = new double[s.length][s[0].length];
//    for(int i = 0; i < cost.target.length ; i++)
//      target[i] = cost.target[i];
      
//    tg = gTrig(cost.target, sV, cost.angle, 1);
    
//    target[cost.target.length]      = tg.M[0];
//    target[cost.target.length + 1]  = tg.M[1];

//    // augment state
//    //i = 1:D0; k = D0+1:D1
//    int i[] = createInd(0, D0),
//        k[] = createInd(D0, D1);
      
//    tgF = gTrigF(getVec(M, i), getMat(S, i, i), cost.angle, 1);
    
//    M = setVec(M, k, tgF.M);
//    S = setMat(S, k, k, tgF.V);
//    C = coeffProdMat(1.0, tgF.C); 
//    mdm = tgF.dMdm;
//    sdm = tgF.dVdm;
//    Cdm = tgF.dCdm;
//    mds = tgF.dMdv;
//    sds = tgF.dVdv;
//    Cds = tgF.dCdv;

//    // compute derivatives (for augmentation)
   
//    int X[][]  = new int [D1][D1];
//    int XT[][] = new int [D1][D1];
////    float I[][] = new float [D1][D1];
    
//    for(int i1 = 0 ; i1 < D1 ; i1++)
//    for(int j1 = 0 ; j1 < D1 ; j1++)
//      X[j1][i1] = j1 + D1 * i1; 
    
//    XT = transMat(X);
    
//    int ii[]  = new int [D0 * D0];
//    for(int i1 = 0 ; i1 < D0 ; i1++)
//    for(int j1 = 0 ; j1 < D0 ; j1++)
//      ii[j1 + D0 * i1] = X[j1][i1]; 
    
//    int kk[]  = new int [(D1 - D0) * (D1 - D0)];
//    int indexKK = 0;
//    for(int i1 = D0 ; i1 < D1 ; i1++)
//    for(int j1 = D0 ; j1 < D1 ; j1++)
//      kk[indexKK++] = X[j1][i1]; 
      
//    int ik[]  = new int [(D1 - D0) * D0];
//    indexKK = 0;
//    for(int i1 = D0 ; i1 < D1 ; i1++)
//    for(int j1 = 0 ; j1 < D0 ; j1++)
//      ik[indexKK++] = X[j1][i1]; 
      
//    int ki[]  = new int [(D1 - D0) * D0];
//    indexKK = 0;
//    for(int i1 = D0 ; i1 < D1 ; i1++)
//    for(int j1 = 0 ; j1 < D0 ; j1++)
//      ki[indexKK++] = XT[j1][i1]; 
    
//  // chainrule
//  // Mdm(k,:)  = mdm*Mdm(i,:) + mds*Sdm(ii,:); 
//    Mdm = setMat(Mdm, k, Chainrule(mdm, Mdm, i, mds, Sdm, ii));
     
//  //Mds(k,:)  = mdm*Mds(i,:) + mds*Sds(ii,:);
//    Mds = setMat(Mds, k, Chainrule(mdm, Mds, i, mds, Sds, ii));
   
//  //Sdm(kk,:) = sdm*Mdm(i,:) + sds*Sdm(ii,:);
//    Sdm = setMat(Sdm, kk, Chainrule(sdm, Mdm, i, sds, Sdm, ii));
  
//  //Sds(kk,:) = sdm*Mds(i,:) + sds*Sds(ii,:);
//    Sds = setMat(Sds, kk, Chainrule(sdm, Mds, i, sds, Sds, ii));
   
//  //dCdm      = Cdm*Mdm(i,:) + Cds*Sdm(ii,:);
//    dCdm = Chainrule(Cdm, Mdm, i, Cds, Sdm, ii);
 
//  //dCds      = Cdm*Mds(i,:) + Cds*Sds(ii,:);
//    dCds = Chainrule(Cdm, Mds, i, Cds, Sds, ii);
   
//  //S(i,k) = S(i,i)*C; 
//  //S(k,i) = S(i,k)';   
//    S = addToSimpl(S, i, k, C);
  
//  //SS = kron(eye(length(k)),S(i,i)); CC = kron(C',eye(length(i)));
//    double SS[][] = kron(eye(k.length), getMat(S, i, i));
//    double CC[][] = kron(transMat(C), eye(i.length));

//  //Sdm(ik,:) = SS*dCdm + CC*Sdm(ii,:); 
//  //Sdm(ki,:) = Sdm(ik,:);
//    Sdm = setMat(Sdm, ik, Chainrule(SS, dCdm, CC, Sdm, ii));
//    Sdm = setMat(Sdm, ki, getMat(Sdm, ik));
   
//  //Sds(ik,:) = SS*dCds + CC*Sds(ii,:); 
//  //Sds(ki,:) = Sds(ik,:);
//    Sds = setMat(Sds, ik, Chainrule(SS, dCds, CC, Sds, ii));
//    Sds = setMat(Sds, ki, getMat(Sds, ik)); 
  
//  }
  
////  //4. Calculate loss!
//  //L = 0; dLdm = zeros(1,D0); dLds = zeros(1,D0*D0); S2 = 0;
//  double L         = 0,
//         dLdm[]    = new double [D0],
//         dLds[]    = new double [D0 * D0],
//         S2        = 0;
//  for(int i_o = 0 ; i_o < 1 ; i_o++) // for i = 1:length(cw) 
//  {
//    //cost.z = target; cost.W = Q/cw(i)^2;
//    cost.z = vectorCoef(target, 1);
//    cost.W = coeffProdMat(1.0 / (cw * cw), Q);
   
//// [r rdM rdS s2 s2dM s2dS] = lossSat(cost, M, S);
//    nargout_lossSat_t dat = new nargout_lossSat_t();
//    dat = lossSat(cost, M, S);
//    double r = dat.L;
//    double rdM[] = dat.dLdm;
//    double rdS[][] = dat.dLds;
//    double s2 = dat.S;
//    double s2dM[] = dat.dSdm;
//    double s2dS[][] = dat.dSds;
  
//    L = L + r; 
//    S2 = S2 + s2;
  
// // dLdm = dLdm + rdM(:)'*Mdm + rdS(:)'*Sdm;
//    double rdST[] = new double[rdS.length * rdS[0].length];
//    int indK = 0;
//    for(int i1 = 0 ; i1 < rdS.length ; i1++)
//      for(int j1 = 0 ; j1 < rdS[0].length ; j1++)
//        rdST[indK++] = rdS[i1][j1];
//    dLdm = vectorSUM(1, dLdm, 1, matrixMultiplyC_T2(rdM, Mdm));
//    dLdm = vectorSUM(1, dLdm, 1, matrixMultiplyC_T2(rdST, Sdm));

//  // dLds = dLds + rdM(:)'*Mds + rdS(:)'*Sds;
//    dLds = vectorSUM(1, dLds, 1, matrixMultiplyC_T2(rdM, Mds));
//    dLds = vectorSUM(1, dLds, 1, matrixMultiplyC_T2(rdST, Sds));

//  //if (b~=0 || ~isempty(b)) && abs(s2)>1e-12
  
//  if(/*(b != 0) && */(abs_db(s2) > 1e-12))
//  {
//    L = L + b * sqrt_db(s2);
 
//    //dLdm = dLdm + b/sqrt(s2) * ( s2dM(:)'*Mdm + s2dS(:)'*Sdm )/2;
//    double s2dST[] = new double[s2dS.length * s2dS[0].length];
//    double dLdmT[] = new double[dLdm.length];
//    double dLdsT[];
//    indK = 0;
//    for(int i = 0 ; i < rdS.length ; i++)
//      for(int j = 0 ; j < rdS[0].length ; j++)
//        s2dST[indK++] = s2dS[i][j];
//    dLdmT = matrixMultiplyC_T2(s2dM, Mdm);
//    dLdmT = vectorSUM(1, dLdmT, 1, matrixMultiplyC_T2(s2dST, Sdm));
//    dLdmT = vectorCoef(dLdmT, b/sqrt_db(s2) * 0.5);
//    dLdm = vectorSUM(1, dLdm, 1, dLdmT); 

//    //dLds = dLds + b/sqrt(s2) * ( s2dM(:)'*Mds + s2dS(:)'*Sds )/2;
//    dLdsT = matrixMultiplyC_T2(s2dM, Mds);
//    dLdsT = vectorSUM(1, dLdsT, 1, matrixMultiplyC_T2(s2dST, Sds));
//    dLdsT = vectorCoef(dLdsT, b/sqrt_db(s2) * 0.5);
//    dLds = vectorSUM(1, dLds, 1, dLdsT); 
//  }
  
//  }
  
//  // normalize n = 1
//  // n = length(cw); L = L/n; dLdm = dLdm/n; dLds = dLds/n; S2 = S2/n
//  ret.L = L;
//  ret.S2 = S2;
//  ret.dLdm = vectorCoef(dLdm, 1.0);
//  ret.dLds = vectorCoef(dLds, 1.0);
    
//  return ret;
//}

//nargout_loss_cp_t loss_cp(cost_t cost, double m[], double s[][])
//{
//  nargout_loss_cp_t ret      = new nargout_loss_cp_t();
//  gTrig_r           tg       = new gTrig_r();
//  gTrig_Full_r      tgF      = new gTrig_Full_r();

//  double   cw                = cost.width_c;
//  double   b                 = cost.expl;
  
//  double   target[];// = new double [D1];
//  double   sV[][];//   = new double [s.length][s[0].length];
//  double   C[][];//    = new double [D0][2];
//  double   mdm[][];//  = new double [2][D0];
//  double   sdm[][];//  = new double [D0][D0];
//  double   Cdm[][];//  = new double [D1 + 2][D0];
//  double   mds[][];//  = new double [2][D0 * D0];
//  double   sds[][];//  = new double [D0][D0 * D0];
//  double   Cds[][];//  = new double [D1 + 2][D0 * D0];
//  double   dCdm[][];
//  double   dCds[][];// = new double [D1][D1];

//  // 1. Some precomputations
//  int      D0 = s[0].length,   // state dimension 
//           D1 = D0 + 2;        // state dimension (with sin/cos)
  
  
//  double   M[]      = new double [D1];
//  for(int i = 0; i < D0 ; i++)
//    M[i] = m[i];
    
//  double   S[][] = new double [D1][D1];
//  for(int i = 0; i < D0 ; i++)
//    for(int j = 0 ; j < D0 ; j++)
//      S[i][j] = s[i][j];
  
//  double   Mdm[][] = new double[D0 + D1 - D0][D0];
//  for(int i = 0; i < D0; i++)
//    Mdm[i][i] = 1;
  
//  double  Sdm[][] = new double [D1 * D1][D0];
//  double  Mds[][] = new double [D1][D0 * D0];
//  double  Sds[][];// = new double [D1 * D1][D0 * D0];
//  Sds = kron(Mdm, Mdm);

//  // 2. Define static penalty as distance from target setpoint
  
//  double ell = cost.p;
//  double Q[][] = new double[D1][D1];
  
//  double t[] = {1.0, ell};
//  double tt[][] = Outer_Product(t, t);
//  Q[0][0]    = tt[0][0];
//  Q[0][D0]   = tt[0][1];
//  Q[D0][0]   = tt[1][0];
//  Q[D0][D0]  = tt[1][1];
//  Q[D0 + 1][D0 + 1] = ell * ell;

//  // 3. Trigonometric augmentation
  
////  if ((D1-D0) > 0)
//  {
//    // augment target
//    //target = [cost.target(:); gTrig(cost.target(:), 0*s, cost.angle)]
    
//    target = new double[cost.target.length + 2];
//    sV     = new double[s.length][s[0].length];
//    for(int i = 0; i < cost.target.length ; i++)
//      target[i] = cost.target[i];
      
//    tg = gTrig(cost.target, sV, cost.angle, 1);
    
//    target[cost.target.length]      = tg.M[0];
//    target[cost.target.length + 1]  = tg.M[1];

//    // augment state
//    //i = 1:D0; k = D0+1:D1
//    int i[] = createInd(0, D0),
//        k[] = createInd(D0, D1);
      
//    tgF = gTrigF(getVec(M, i), getMat(S, i, i), cost.angle, 1);
    
//    M = setVec(M, k, tgF.M);
//    S = setMat(S, k, k, tgF.V);
//    C = coeffProdMat(1.0, tgF.C); 
//    mdm = tgF.dMdm;
//    sdm = tgF.dVdm;
//    Cdm = tgF.dCdm;
//    mds = tgF.dMdv;
//    sds = tgF.dVdv;
//    Cds = tgF.dCdv;

//    // compute derivatives (for augmentation)
   
//    int X[][]  = new int [D1][D1];
//    int XT[][] = new int [D1][D1];
////    float I[][] = new float [D1][D1];
    
//    for(int i1 = 0 ; i1 < D1 ; i1++)
//    for(int j1 = 0 ; j1 < D1 ; j1++)
//      X[j1][i1] = j1 + D1 * i1; 
    
//    XT = transMat(X);
    
//    int ii[]  = new int [D0 * D0];
//    for(int i1 = 0 ; i1 < D0 ; i1++)
//    for(int j1 = 0 ; j1 < D0 ; j1++)
//      ii[j1 + D0 * i1] = X[j1][i1]; 
    
//    int kk[]  = new int [(D1 - D0) * (D1 - D0)];
//    int indexKK = 0;
//    for(int i1 = D0 ; i1 < D1 ; i1++)
//    for(int j1 = D0 ; j1 < D1 ; j1++)
//      kk[indexKK++] = X[j1][i1]; 
      
//    int ik[]  = new int [(D1 - D0) * D0];
//    indexKK = 0;
//    for(int i1 = D0 ; i1 < D1 ; i1++)
//    for(int j1 = 0 ; j1 < D0 ; j1++)
//      ik[indexKK++] = X[j1][i1]; 
      
//    int ki[]  = new int [(D1 - D0) * D0];
//    indexKK = 0;
//    for(int i1 = D0 ; i1 < D1 ; i1++)
//    for(int j1 = 0 ; j1 < D0 ; j1++)
//      ki[indexKK++] = XT[j1][i1]; 
    
//  // chainrule
//  // Mdm(k,:)  = mdm*Mdm(i,:) + mds*Sdm(ii,:); 
//    Mdm = setMat(Mdm, k, Chainrule(mdm, Mdm, i, mds, Sdm, ii));
     
//  //Mds(k,:)  = mdm*Mds(i,:) + mds*Sds(ii,:);
//    Mds = setMat(Mds, k, Chainrule(mdm, Mds, i, mds, Sds, ii));
   
//  //Sdm(kk,:) = sdm*Mdm(i,:) + sds*Sdm(ii,:);
//    Sdm = setMat(Sdm, kk, Chainrule(sdm, Mdm, i, sds, Sdm, ii));
  
//  //Sds(kk,:) = sdm*Mds(i,:) + sds*Sds(ii,:);
//    Sds = setMat(Sds, kk, Chainrule(sdm, Mds, i, sds, Sds, ii));
   
//  //dCdm      = Cdm*Mdm(i,:) + Cds*Sdm(ii,:);
//    dCdm = Chainrule(Cdm, Mdm, i, Cds, Sdm, ii);
 
//  //dCds      = Cdm*Mds(i,:) + Cds*Sds(ii,:);
//    dCds = Chainrule(Cdm, Mds, i, Cds, Sds, ii);
   
//  //S(i,k) = S(i,i)*C; 
//  //S(k,i) = S(i,k)';   
//    S = addToSimpl(S, i, k, C);
  
//  //SS = kron(eye(length(k)),S(i,i)); CC = kron(C',eye(length(i)));
//    double SS[][] = kron(eye(k.length), getMat(S, i, i));
//    double CC[][] = kron(transMat(C), eye(i.length));

//  //Sdm(ik,:) = SS*dCdm + CC*Sdm(ii,:); 
//  //Sdm(ki,:) = Sdm(ik,:);
//    Sdm = setMat(Sdm, ik, Chainrule(SS, dCdm, CC, Sdm, ii));
//    Sdm = setMat(Sdm, ki, getMat(Sdm, ik));
   
//  //Sds(ik,:) = SS*dCds + CC*Sds(ii,:); 
//  //Sds(ki,:) = Sds(ik,:);
//    Sds = setMat(Sds, ik, Chainrule(SS, dCds, CC, Sds, ii));
//    Sds = setMat(Sds, ki, getMat(Sds, ik)); 
  
//  }
  
////  //4. Calculate loss!
//  //L = 0; dLdm = zeros(1,D0); dLds = zeros(1,D0*D0); S2 = 0;
//  double L         = 0,
//         dLdm[]    = new double [D0],
//         dLds[]    = new double [D0 * D0],
//         S2        = 0;
//  for(int i_o = 0 ; i_o < 1 ; i_o++) // for i = 1:length(cw) 
//  {
//    //cost.z = target; cost.W = Q/cw(i)^2;
//    cost.z = vectorCoef(target, 1);
//    cost.W = coeffProdMat(1.0 / (cw * cw), Q);
   
//// [r rdM rdS s2 s2dM s2dS] = lossSat(cost, M, S);
//    nargout_lossSat_t dat = new nargout_lossSat_t();
//    dat = lossSat(cost, M, S);
//    double r = dat.L;
//    double rdM[] = dat.dLdm;
//    double rdS[][] = dat.dLds;
//    double s2 = dat.S;
//    double s2dM[] = dat.dSdm;
//    double s2dS[][] = dat.dSds;
  
//    L = L + r; 
//    S2 = S2 + s2;
  
// // dLdm = dLdm + rdM(:)'*Mdm + rdS(:)'*Sdm;
//    double rdST[] = new double[rdS.length * rdS[0].length];
//    int indK = 0;
//    for(int i1 = 0 ; i1 < rdS.length ; i1++)
//      for(int j1 = 0 ; j1 < rdS[0].length ; j1++)
//        rdST[indK++] = rdS[i1][j1];
//    dLdm = vectorSUM(1, dLdm, 1, matrixMultiplyC_T2(rdM, Mdm));
//    dLdm = vectorSUM(1, dLdm, 1, matrixMultiplyC_T2(rdST, Sdm));

//  // dLds = dLds + rdM(:)'*Mds + rdS(:)'*Sds;
//    dLds = vectorSUM(1, dLds, 1, matrixMultiplyC_T2(rdM, Mds));
//    dLds = vectorSUM(1, dLds, 1, matrixMultiplyC_T2(rdST, Sds));

//  //if (b~=0 || ~isempty(b)) && abs(s2)>1e-12
  
//  if(/*(b != 0) && */(abs_db(s2) > 1e-12))
//  {
//    L = L + b * sqrt_db(s2);
 
//    //dLdm = dLdm + b/sqrt(s2) * ( s2dM(:)'*Mdm + s2dS(:)'*Sdm )/2;
//    double s2dST[] = new double[s2dS.length * s2dS[0].length];
//    double dLdmT[] = new double[dLdm.length];
//    double dLdsT[];
//    indK = 0;
//    for(int i = 0 ; i < rdS.length ; i++)
//      for(int j = 0 ; j < rdS[0].length ; j++)
//        s2dST[indK++] = s2dS[i][j];
//    dLdmT = matrixMultiplyC_T2(s2dM, Mdm);
//    dLdmT = vectorSUM(1, dLdmT, 1, matrixMultiplyC_T2(s2dST, Sdm));
//    dLdmT = vectorCoef(dLdmT, b/sqrt_db(s2) * 0.5);
//    dLdm = vectorSUM(1, dLdm, 1, dLdmT); 

//    //dLds = dLds + b/sqrt(s2) * ( s2dM(:)'*Mds + s2dS(:)'*Sds )/2;
//    dLdsT = matrixMultiplyC_T2(s2dM, Mds);
//    dLdsT = vectorSUM(1, dLdsT, 1, matrixMultiplyC_T2(s2dST, Sds));
//    dLdsT = vectorCoef(dLdsT, b/sqrt_db(s2) * 0.5);
//    dLds = vectorSUM(1, dLds, 1, dLdsT); 
//  }
  
//  }
  
//  // normalize n = 1
//  // n = length(cw); L = L/n; dLdm = dLdm/n; dLds = dLds/n; S2 = S2/n
  
//  ret.L = L;
//  ret.S2 = S2;
//  ret.dLdm = vectorCoef(dLdm, 1.0);
//  ret.dLds = vectorCoef(dLds, 1.0);
//  return ret;
//}

//************************************************************************
//               lossSat
//************************************************************************


nargout_lossSat_t lossSat(cost_t cost, double m[], double s[][])
{
  
  nargout_lossSat_t ret      = new nargout_lossSat_t();
  
  //D = length(m);  get state dimension
  int D = m.length;
  
  //set some defaults if necessary
  double W[][] = new double [cost.W.length][cost.W[0].length] ;
  W = cost.W;
  double z[]   = new double [cost.z.length];
  z = cost.z;

  //SW = s*W;
  //iSpW = W/(eye(D)+SW);
  double SW[][]  = new double[D][D],
         SWt[][]  = new double[D][D],
         iSpW[][]  = new double[D][D];
  SW = matrixMultiply(s, W);
  SWt = matrixADD(eye(D), SW);
  iSpW = solveMat(W, SWt);
  
  // 1. Expected cost
  //L = -exp(-(m-z)'*iSpW*(m-z)/2)/sqrt(det(eye(D)+SW)); % in interval [-1,0]
  double L = -0.5 * matrixMultiplyC_V(matrixMultiplyC_T2(vectorSUM(1, m, -1, z), iSpW), vectorSUM(1, m, -1, z));
  L = -exp_db(L) / sqrt_db(detMatrix(SWt));
  
  // 1a. derivatives of expected cost
  //  dLdm = -L*(m-z)'*iSpW;  % wrt input mean
  double dLdm[];
  dLdm = matrixMultiplyC_T2(vectorSUM(-L, m, L, z), iSpW);
  //  dLds = L*(iSpW*(m-z)*(m-z)'-eye(D))*iSpW/2;  % wrt input covariance matrix
  double dLds[][];
  dLds = matrixMultiply(iSpW, Outer_Product(vectorSUM(1, m, -1, z), vectorSUM(1, m, -1, z)));
  dLds = matrixSUB(dLds, eye(D));
  dLds = matrixMultiply(matrixCoef(dLds, L * 0.5), iSpW);
 
  // 2. Variance of cost
  // i2SpW = W/(eye(D)+2*SW);
  double i2SpW[][]  = new double[D][D];
  SWt = matrixADD(eye(D), matrixCoef(SW, 2.0));
  i2SpW = solveMat(W, SWt);
  // r2 = exp(-(m-z)'*i2SpW*(m-z))/sqrt(det(eye(D)+2*SW));
  double r2 = -matrixMultiplyC_V(matrixMultiplyC_T2(vectorSUM(1, m, -1, z), i2SpW), vectorSUM(1, m, -1, z));
  r2 = exp_db(r2) / sqrt_db(detMatrix(SWt));
  double S = r2 - L * L;
  if(S < 1e-12) S = 0; //for numerical reasons
  
  // 2a. derivatives of variance of cost
  // wrt input mean
  // dSdm = -2*r2*(m-z)'*i2SpW-2*L*dLdm;
  double dSdm[];
  dSdm = matrixMultiplyC_T2(vectorCoef(vectorSUM(1, m, -1, z), -2.0 * r2), i2SpW);
  dSdm = vectorSUM(1, dSdm, -2.0 * L, dLdm);
  // wrt input covariance matrix
  // dSds = r2*(2*i2SpW*(m-z)*(m-z)'-eye(D))*i2SpW-2*L*dLds;
  double dSds[][];
  dSds = matrixMultiply(matrixCoef(i2SpW, 2.0), Outer_Product(vectorSUM(1, m, -1, z), vectorSUM(1, m, -1, z)));
  dSds = matrixSUB(dSds, eye(D));
  dSds = matrixMultiply(matrixCoef(dSds, r2), i2SpW);
  dSds = matrixSUB(dSds, matrixCoef(dLds, 2.0 * L));
  
  L = 1 + L;
  ret.L = L;
  ret.dLdm =dLdm;// vectorCoef(dLdm, 1.0);
  ret.dLds = dLds;//coeffProdMat(1.0, dLds);
  ret.S = S;
  ret.dSdm = dSdm;//vectorCoef(dSdm, 1.0);
  ret.dSds = dSds;//coeffProdMat(1.0, dSds);
  
  return ret;
  
}
