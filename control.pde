
tp_p policeInit(policy_t pol,  int nct)
{
    tp_p p = new tp_p();
    
    if(mode_pol == 0)
    {
        //[mm ss cc] = gTrig(mu0, S0, plant.angi)
        //mm = [mu0; mm]; cc = S0*cc; ss = [S0 cc; cc' ss] 
        
        double mu0[] = new double [4];
        double S0t[][] = new double [4][4];
        for(int i = 0; i < S0t.length ; i++)
          S0t[i][i] = 0.1 * 0.1; 
          S0t[2][2] *= 0.01;
          S0t[3][3] *= 0.01;
          
        gTrig_r getgTrig = gTrig(mu0, S0t, pol.angle, pol.e);
        
        double mm[]   = new double [mu0.length + getgTrig.M.length];
        double cc[][] = new double [getgTrig.V.length][getgTrig.V[0].length];
             
        //mm = [mu0; mm]; cc = S0*cc
        for(int i = 0; i < mu0.length ; i++)
          mm[i] = mu0[i];
          
        for(int i = 0; i < getgTrig.M.length ; i++)
          mm[i + mu0.length] = getgTrig.M[i];
        
        cc = matrixMultiply(S0t, getgTrig.C);
        
        //ss = [S0 cc; cc' ss]
        double ss[][] = new double [S0t.length + getgTrig.V.length][S0t[0].length + getgTrig.V[0].length];
     
        //S0 cc
        for(int i = 0; i < S0t.length ; i++)
        for(int j = 0; j < S0t[0].length; j++)
          ss[i][j] = S0t[i][j];
        for(int i = 0; i < cc.length ; i++)
        for(int j = 0; j < cc[0].length; j++)
          ss[i][j + S0t[0].length] = cc[i][j];
        //cc' ss
        for(int i = 0; i < cc[0].length ; i++)
        for(int j = 0; j < cc.length; j++)
          ss[i + S0t.length][j] = cc[j][i];
        for(int i = 0; i < getgTrig.V.length ; i++)
        for(int j = 0; j < getgTrig.V[0].length; j++)
          ss[i + S0t.length][j + S0t[0].length] = getgTrig.V[i][j]; 
    
        //policy.p.inputs = gaussian(mm(poli), ss(poli,poli), nc)';
        double mm_pol[]   = new double [pol.poli.length]; // Угловая скорост / синус / косинус
        double ss_pol[][] = new double [pol.poli.length][pol.poli.length]; // Угловая скорост / синус / косинус
        for(int i = 0; i < pol.poli.length ; i++)
          mm_pol[i] = mm[pol.poli[i]];
        
        for(int i = 0; i < pol.poli.length ; i++)
          for(int j = 0; j < pol.poli.length ; j++)
          ss_pol[i][j] = ss[pol.poli[i]][pol.poli[j]];
           
        p.inputs = transMat(gaussian(mm_pol, ss_pol, nct));
        
        //policy.p.targets = 0.1*randn(nc, length(policy.maxU));
        p.targets = new double [nct];
        for(int i = 0; i < nct; i ++)
         p.targets[i]  = 0.1 * (double)randomGaussian();
         
        //log([1 1 0.7 0.7 0.7 0.7 1 0.01]')
        p.hyp = new double [8];
        p.hyp[0] = log_db(1);
        p.hyp[1] = log_db(1);
        p.hyp[2] = log_db(0.7);
        p.hyp[3] = log_db(0.7);
        p.hyp[4] = log_db(0.7);
        p.hyp[5] = log_db(0.7);
        p.hyp[6] = log_db(1);
        p.hyp[7] = log_db(0.01);
    }
    else
    {
         p.w = new double [pol.poli.length];
         
         p.b = 1.5848954372683464;
         p.w[0] = 0.600854243152309; 
         p.w[1] = 0.09901011125505928;  
         p.w[2] = 0.2537614483114453; 
         p.w[3] = 5.387786993001714; 
         p.w[4] = -1.6444470556774193;  
        
 //i = 0 : 7.954053817258841  
 //i = 1 : 0.25903051327137877  
 //i = 2 : 0.12443211638970815  
 //i = 3 : 0.4226601145726813  
 //i = 4 : 4.917988535997297  
 //i = 5 : -7.580887834370011  
        ////policy.p.w = 1e-2*randn(length(policy.maxU),length(poli));
        //for(int i = 0; i < p.w.length ; i++)
        //  p.w[i] = 0.01 * (double)randomGaussian() ;
       
        //p.b = 0;
    }
  return p;
};

//function [M, S, C, dMdm, dSdm, dCdm, dMds, dSds, dCds, dMdp, dSdp, dCdp] = congp(policy, m, s)
nargout_conpols_t conpols(policy_t _policy, double m[], double s[][])
{
  nargout_conpols_t ret      = new nargout_conpols_t();
  
  if(mode_pol == 0)
  {
    // GP mode
    nargout_gp0_t   gp0      = new nargout_gp0_t();
    policyGP = new dynmodel_t(_policy.p.inputs.length, _policy.p.inputs[0].length, _policy.p.targets.length);
    
    // 1. Extract policy parameters
    policyGP.hyp    = vectorCoef2(_policy.p.hyp, 1.0);
    policyGP.inputs =   _policy.p.inputs ; //coeffProdMat(1.0, _policy.p.inputs);
    policyGP.targets = vectorCoef2(_policy.p.targets, 1.0);
    
    //% fix policy signal and the noise variance 
    //(avoids some potential numerical problems)
    //policy.hyp(end-1,:) = log(1);                  set signal variance to 1
    //policy.hyp(end,:) = log(0.01);                 set noise standard dev to 0.01
    policyGP.hyp[policyGP.hyp.length - 2][0] = log_db(1.0);
    policyGP.hyp[policyGP.hyp.length - 1][0] = log_db(0.01);
    
    //2. Compute predicted control u inv(s)*covariance between input and control
    policyGP = gp2(policyGP);
    gp0 = _gp2(policyGP, m, s); //<>// //<>// //<>//
    
    ret.M =gp0.M;// vectorCoef(gp0.M, 1.0);
    ret.S = gp0.S;//coeffProdMat(1.0, gp0.S);
    ret.C = gp0.V;//coeffProdMat(1.0, gp0.V);
    // GP mode
  }
  else
  {
    double w[] = policy.p.w,//vectorCoef(policy.p.w, 1.0),
           b = policy.p.b;
    int    E = 1,
           D = w.length;
    
    //% 2. Predict control signal
    //M = w*m + b;                                                        % mean
    //S = w*s*w'; S = (S+S')/2;                                     % covariance
    //V = w';                                   % inv(s)*input-output covariance
    
    double M = b + matrixMultiplyC_V(w, m),
           S = matrixMultiplyC_V(matrixMultiplyC_T2(w, s), w),
           V[] = w;//vectorCoef(w, 1.0);
    
    ret.M = new double[1];
    ret.S = new double[1][1];
    ret.M[0] = M;
    ret.S[0][0] = S;
    ret.C = new double[V.length][1];
    for(int i = 0; i < V.length ; i++)
      ret.C[i][0] = V[i];
  }
  
  return ret;
  
}

nargout_conpolsd_t conpolsd(policy_t _policy, double m[], double s[][])
{
  
  nargout_conpolsd_t ret = new nargout_conpolsd_t();
  
  if(mode_pol == 0)
  {
      // GP mode
      nargout_gp2_t getgp2 = new nargout_gp2_t();
      dynmodel_t    policyGP = new dynmodel_t(_policy.p.inputs.length, _policy.p.inputs[0].length, _policy.p.targets.length);
      
      // 1. Extract policy parameters
      policyGP.hyp    =vectorCoef2(_policy.p.hyp, 1.0);
      policyGP.inputs =  _policy.p.inputs;//coeffProdMat(1.0, _policy.p.inputs);
      policyGP.targets = vectorCoef2(_policy.p.targets, 1.0);
      
      //% fix policy signal and the noise variance 
      //(avoids some potential numerical problems)
      //policy.hyp(end-1,:) = log(1);                  set signal variance to 1
      //policy.hyp(end,:) = log(0.01);                 set noise standard dev to 0.01
      policyGP.hyp[policyGP.hyp.length - 2][0] = log(1.0);
      policyGP.hyp[policyGP.hyp.length - 1][0] = log(0.01);
      
      //2. Compute predicted control u inv(s)*covariance between input and control
      policyGP = gp2(policyGP);
      getgp2 = _gp2d(policyGP, m, s);
    
      //d = size(policy.inputs,2);            
      //d2 = size(policy.hyp,1); dimU = size(policy.targets,2);
      //sidx = bsxfun(@plus,(d+1:d2)',(0:dimU-1)*d2)
      int d   = policyGP.inputs[0].length,
          d2  = policyGP.hyp.length,
          sidx[] = new int [d2 - d];        // Считаем, что управляющий сигнал один, поэтому вектор индексов, а не матрица
       for(int i = 0; i < sidx.length ; i++)
        sidx[i] = d  + i;
       
       //dMdh(:,sidx(:)) = 0; dSdh(:,sidx(:)) = 0; dCdh(:,sidx(:)) = 0
       for(int i = 0; i < sidx.length ; i++)
       {
         getgp2.dMdX[0][sidx[i]] = 0;
         getgp2.dSdX[0][sidx[i]] = 0;
       }
       
       for(int i = 0; i < sidx.length ; i++)
       {
         for(int j = 0; j < getgp2.dVdX.length ; j++)
         getgp2.dVdX[j][sidx[i]] = 0;
       }
      
      // % 4. Merge derivatives
      // dMdp = [dMdh dMdi dMdt]; dSdp = [dSdh dSdi dSdt]; dCdp = [dCdh dCdi dCdt];
      ret.M = getgp2.M;    //vectorCoef(getgp2.M, 1.0);
      ret.S = getgp2.S;    //coeffProdMat(1.0, getgp2.S);
      ret.C = getgp2.V;   //coeffProdMat(1.0, getgp2.V);
      ret.dMdm =getgp2.dMdm;// coeffProdMat(1.0, getgp2.dMdm);
      ret.dSdm =getgp2.dSdm; //coeffProdMat(1.0, getgp2.dSdm);
      ret.dCdm =getgp2.dVdm; //coeffProdMat(1.0, getgp2.dVdm);
      ret.dMds =getgp2.dMds; //coeffProdMat(1.0, getgp2.dMds);
      ret.dSds =getgp2.dSds; //coeffProdMat(1.0, getgp2.dSds);
      ret.dCds =getgp2.dVds;// coeffProdMat(1.0, getgp2.dVds);
      //*************
        ret.dMdp = new double [getgp2.dMdX[0].length + getgp2.dMdi[0].length + getgp2.dMdt[0].length];
        int j = 0;
        for(int i = 0; i < getgp2.dMdX[0].length ; i++)
        {
          ret.dMdp[j] = getgp2.dMdX[0][i];
          j++;
        }
        for(int i = 0; i < getgp2.dMdi[0].length ; i++)
        {
          ret.dMdp[j] = getgp2.dMdi[0][i];
          j++;
        }
        for(int i = 0; i < getgp2.dMdt[0].length ; i++)
        {
          ret.dMdp[j] = getgp2.dMdt[0][i];
          j++;
        }
        //*************
        ret.dSdp = new double [getgp2.dSdX[0].length + getgp2.dSdi[0].length + getgp2.dSdt[0].length];
        j = 0;
        for(int i = 0; i < getgp2.dSdX[0].length ; i++)
        {
          ret.dSdp[j] = getgp2.dSdX[0][i];
          j++;
        }
        for(int i = 0; i < getgp2.dSdi[0].length ; i++)
        {
          ret.dSdp[j] = getgp2.dSdi[0][i];
          j++;
        }
        for(int i = 0; i < getgp2.dSdt[0].length ; i++)
        {
          ret.dSdp[j] = getgp2.dSdt[0][i];
          j++;
        }
        //*************
        ret.dCdp = new double [getgp2.dVdX.length][getgp2.dVdX[0].length + getgp2.dVdi[0].length + getgp2.dVdt[0].length];
        j = 0;
        for(int i = 0; i < getgp2.dSdX[0].length ; i++)
        {
          for(int k = 0; k < ret.dCdp.length; k++)
            ret.dCdp[k][j] = getgp2.dVdX[k][i];
          j++;
        }
        for(int i = 0; i < getgp2.dVdi[0].length ; i++)
        {
          for(int k = 0; k < ret.dCdp.length; k++)
            ret.dCdp[k][j] = getgp2.dVdi[k][i];
          j++;
        }
        for(int i = 0; i < getgp2.dVdt[0].length ; i++)
        {
          for(int k = 0; k < ret.dCdp.length; k++)
            ret.dCdp[k][j] = getgp2.dVdt[k][i];
          j++;
        }
      // GP mode  
  }
  else
  {
      //% 1. Extract policy parameters from policy structure
      //w = policy.p.w;                                 % weight matrix
      //b = policy.p.b;                                 % bias/offset
      //[E D] = size(w);    
      
      double w[] = policy.p.w,//vectorCoef(policy.p.w, 1.0),
             b = policy.p.b;
      int    E = 1,
             D = w.length;
      
      //% 2. Predict control signal
      //M = w*m + b;                                                        % mean
      //S = w*s*w'; S = (S+S')/2;                                     % covariance
      //V = w';                                   % inv(s)*input-output covariance
      
      double M = b + matrixMultiplyC_V(w, m),
             S = matrixMultiplyC_V(matrixMultiplyC_T2(w, s), w),
             V[] =w;// vectorCoef(w, 1.0);
      
      ret.M = new double[1];
      ret.S = new double[1][1];
      
      ret.M[0] = M;
      ret.S[0][0] = S;
      ret.C = new double[V.length][1];
      for(int i = 0; i < V.length ; i++)
        ret.C[i][0] = V[i];
      
      //dMdm = w;            dSdm = zeros(E*E,D); dVdm = zeros(D*E,D)
      //dMds = zeros(E,D*D); dSds = kron(w,w);    dVds = zeros(D*E,D*D)
      
      ret.dMdm = new double[1][w.length];
      ret.dMdm[0] = w;//vectorCoef(w, 1.0);
      ret.dSdm = new double[E * E][D]; 
      ret.dCdm = new double[D * E][D];
      ret.dMds = new double[E][D * D];
      ret.dSds = kron(w,w);
      ret.dCds = new double[D * E][D * D];
      //X=reshape(1:D*D,[D D]); XT=X'; dSds=(dSds+dSds(:,XT(:)))/2; % symmetrize
      //X=reshape(1:E*E,[E E]); XT=X'; dSds=(dSds+dSds(XT(:),:))/2
      ret.dSds = symmetrizeB(D, ret.dSds);
      ret.dSds = symmetrizeA(E, ret.dSds);
     
      //wTdw =reshape(permute(reshape(eye(E*D),[E D E D]),[2 1 3 4]),[E*D E*D])
      double wTdw[][] = eye(5);
      
      //dMdp = [eye(E) kron(m',eye(E))]
      ret.dMdp = new double[m.length + 1];
      ret.dMdp[0] = 1;
      for(int i = 0; i < m.length ; i++)
        ret.dMdp[i + 1] = m[i];
      
      //dSdp = [zeros(E*E,E) kron(eye(E),w*s)*wTdw + kron(w*s,eye(E))]
      ret.dSdp = new double[w.length + 1];
      ret.dSdp[0] = 0;
      double dSdp[] = matrixMultiplyC_T2(w, s),
             dSdp_ex[][] = new double[1][w.length + 1];
      for(int i = 0; i < w.length ; i++)
      {
        ret.dSdp[i + 1] = 2 * dSdp[i];
        dSdp_ex[0][i + 1] = ret.dSdp[i + 1];
      }
      //dSdp = (dSdp + dSdp(XT(:),:))/2
      dSdp_ex = symmetrizeA(E, dSdp_ex);
      for(int i = 0; i < (w.length + 1) ; i++)
      {
        ret.dSdp[i] = dSdp_ex[0][i];
      }
      
      //dVdp = [zeros(D*E,E) wTdw]
      ret.dCdp = new double[w.length][w.length + 1];
      for(int i = 0; i < (w.length) ; i++)
      {
        ret.dCdp[i][i + 1] = 1;
      }
  }
  return ret;
  
}

nargout_gSin_t gSin(double m[], double v[][], int i[], double e[])
{
  
  nargout_gSin_t ret = new nargout_gSin_t();
  
  //d = length(m); I = length(i)
  int     d            = m.length,
          I            = i.length;
  
  //mi(1:I,1) = m(i); vi = v(i,i); vii(1:I,1) = diag(vi)
  double  mi[]         = new double [i.length],    
          vi[][]       = new double [i.length][i.length],
          vii[]        = new double [i.length];
  
  for(int j = 0; j < i.length ; j ++)
    mi[j] = m[i[j]];
  
  for(int j = 0; j < i.length ; j ++)
  for(int k = 0; k < i.length ; k ++)
    vi[j][k] = v[i[j]][i[k]];
    
  for(int j = 0; j < i.length ; j ++)
    vii[j] = v[i[j]][i[j]];

  //M = e.*exp(-vii/2).*sin(mi)
  double   M[] = new double [vii.length]; 
  for(int j = 0; j < i.length ; j ++)
    M[j] = e[j] * exp_db(-vii[j] / 2.0) * sin_db(mi[j]);
  
  ret.M = M;//vectorCoef(M, 1.0);
  
  //lq = -bsxfun(@plus,vii,vii')/2; q = exp(lq);
  double lq[][] = new double[vii.length][vii.length],
         q[][];
  lq = bsxfun_plus(vii, vii);    
  lq = coeffProdMat(-0.5 , lq);
  q = EXP_F(lq);
  
  
  //V = (exp(lq+vi)-q).*cos(bsxfun(@minus,mi,mi')) - (exp(lq-vi)-q).*cos(bsxfun(@plus,mi,mi'))
  //V = e*e'.*V/2;                   
  double V[][] = new double[mi.length][mi.length];
  double V_[][] = EXP_F(matrixADD(lq, vi));
  V_ = matrixSUB(V_, q);
  V = COS_F(bsxfun_minus(mi, mi));
  V = matrixDOT(V_, V);
  
  V_= EXP_F(matrixSUB(lq, vi));
  V_ = matrixSUB(V_, q);
  V_ = matrixDOT(V_, COS_F(bsxfun_plus(mi, mi)));
  V = matrixSUB(V, V_);
  
  V = matrixDOT(Outer_Product(e, e), V);
  V = coeffProdMat(0.5, V);
  
  ret.V = V;//coeffProdMat(1.0, V);
  
  //C = zeros(d,I); C(i,:) = diag(e.*exp(-vii/2).*cos(mi))
  double C[][] = new double [d][I],
         C_[][] = new double [I][I];
  for(int j = 0; j < I ; j ++)
    C_[j][j] = e[j] * exp_db(-vii[j] / 2.0) * cos_db(mi[j]);
  
   for(int j = 0; j < I ; j ++)
   for(int k = 0; k < I ; k ++)
     C[i[j]][k] = C_[j][k];
  
  ret.C = C;//coeffProdMat(1.0, C);
  
  //dVdm = zeros(I,I,d); dCdm = zeros(d,I,d); dVdv = zeros(I,I,d,d)
  //dCdv = zeros(d,I,d,d); dMdm = C'
  double dVdm[][][] = new double[I][I][d],
         dCdm[][][] = new double[d][I][d],
         dVdv[][][][] = new double [I][I][d][d],
         dCdv[][][][] = new double[d][I][d][d],
         dMdm[][];
  dMdm = coeffProdMat(1.0, transMat(C));
  ret.dMdm = coeffProdMat(1.0, dMdm);
  
  //U1 = -(exp(lq+vi)-q).*sin(bsxfun(@minus,mi,mi'))
  //U2 = (exp(lq-vi)-q).*sin(bsxfun(@plus,mi,mi'))
  double U1[][] = new double [mi.length][mi.length],
         U2[][] = new double [mi.length][mi.length];
  double U12[][] = EXP_F(matrixADD(lq, vi));
  U12 = matrixSUB(U12, q);
  U1 = SIN_F(bsxfun_minus(mi, mi));
  U1 = matrixDOT(U12, U1);
  U1 = coeffProdMat(-1.0, U1);
  
  U12 = EXP_F(matrixSUB(lq, vi));
  U12 = matrixSUB(U12, q);
  U2 = matrixDOT(U12, SIN_F(bsxfun_plus(mi, mi)));
 
  // for j = 1:I
  //u = zeros(I,1); u(j) = 1/2
  for(int j = 0 ; j < I ; j++)
  {
    double u[] = new double [I];
    u[j] = 0.5;
  
    //dVdm(:,:,i(j)) = e*e'.*(U1.*bsxfun(@minus,u,u') + U2.*bsxfun(@plus,u,u'))
    double dVdm_[][];
    dVdm_ = matrixDOT(U1, bsxfun_minus(u, u));
    dVdm_ = matrixADD(dVdm_, matrixDOT(U2, bsxfun_plus(u, u)));
    dVdm_ = matrixDOT(Outer_Product(e, e), dVdm_);
    for(int j1 = 0 ; j1 < dVdm_.length ; j1++)
    for(int j2 = 0 ; j2 < dVdm_[0].length ; j2++)
      dVdm[j1][j2][i[j]] = dVdm_[j1][j2];
    
    //dVdv(j,j,i(j),i(j)) = exp(-vii(j)) * (1+(2*exp(-vii(j))-1)*cos(2*mi(j)))*e(j)*e(j)/2
    dVdv[j][j][i[j]][i[j]] = exp_db(-vii[j]) * (1 + (2.0 * exp_db(-vii[j]) - 1) 
                                * cos_db(2.0 * mi[j])) * e[j] * e[j] / 2.0;
    
    //for k = [1:j-1 j+1:I]
    for(int k = 0 ; k <= (j - 1) ; k++)
    {
      
      // dVdv(j,k,i(j),i(k)) = (exp(lq(j,k)+vi(j,k)).*cos(mi(j)-mi(k)) + ...
      //                      exp(lq(j,k)-vi(j,k)).*cos(mi(j)+mi(k)))*e(j)*e(k)/2;
      
      dVdv[j][k][i[j]][i[k]] = (exp_db(lq[j][k] + vi[j][k]) * cos_db(mi[j] - mi[k]) +
                               exp_db(lq[j][k] - vi[j][k]) * cos_db(mi[j] + mi[k])) * e[j] * e[k] * 0.5;
      
      //dVdv(j,k,i(j),i(j)) = -V(j,k)/2 
      //dVdv(j,k,i(k),i(k)) = -V(j,k)/2 
      
      dVdv[j][k][i[j]][i[j]] = -V[j][k] / 2.0;
      dVdv[j][k][i[k]][i[k]] = -V[j][k] / 2.0;
      
    }
    for(int k = (j + 1) ; k < I ; k++)
    {
      
      // dVdv(j,k,i(j),i(k)) = (exp(lq(j,k)+vi(j,k)).*cos(mi(j)-mi(k)) + ...
      //                      exp(lq(j,k)-vi(j,k)).*cos(mi(j)+mi(k)))*e(j)*e(k)/2;
      
      dVdv[j][k][i[j]][i[k]] = (exp_db(lq[j][k] + vi[j][k]) * cos_db(mi[j] - mi[k]) +
                               exp_db(lq[j][k] - vi[j][k]) * cos_db(mi[j] + mi[k])) * e[j] * e[k] * 0.5;
      
      //dVdv(j,k,i(j),i(j)) = -V(j,k)/2 
      //dVdv(j,k,i(k),i(k)) = -V(j,k)/2 
      
      dVdv[j][k][i[j]][i[j]] = -V[j][k] / 2.0;
      dVdv[j][k][i[k]][i[k]] = -V[j][k] / 2.0;
      
    }
    
    //dCdm(i(j),j,i(j)) = -M(j)
    //dCdv(i(j),j,i(j),i(j)) = -C(i(j),j)/2
    
    dCdm[i[j]][j][i[j]] = -M[j];
    dCdv[i[j]][j][i[j]][i[j]] = -C[i[j]][j] / 2.0;
    
  }
  
  //dMdv = permute(dCdm,[2 1 3])/2;
   double dMdv[][][] = new double[dCdm[0].length][dCdm.length][dCdm[0][0].length];
     for(int j13 = 0; j13 < dCdm[0][0].length; j13++)
     for(int j11 = 0; j11 < dCdm.length; j11++)
       for(int j12 = 0; j12 < dCdm[0].length; j12++)
         dMdv[j12][j11][j13] = 0.5 * dCdm[j11][j12][j13];
  
  //dMdv = reshape(dMdv,[I d*d])
  ret.dMdv = reshape(dMdv, I, d * d);
  
  // dVdv = reshape(dVdv,[I*I d*d]) 
  ret.dVdv = reshape(dVdv, I * I, d * d);
  
  //dVdm = reshape(dVdm,[I*I d])
  ret.dVdm = reshape(dVdm, I * I, d);
  
  //dCdv = reshape(dCdv,[d*I d*d]) 
  ret.dCdv = reshape(dCdv, d * I, d * d);//----------------
  
  //dCdm = reshape(dCdm,[d*I d]);
  ret.dCdm = reshape(dCdm, d * I, d);

  return ret;
  
}

nargout_gSat_t gSat(double m[], double v[][], int i, double e)
{
  
  nargout_gSat_t ret = new nargout_gSat_t();
  nargout_gSin_t gSinT = new nargout_gSin_t();
  
  //d = length(m); I = length(i); i = i(:)';
  //if nargin < 4; e = ones(1, I); end; e = e(:)';
  int d = m.length,
      I = 1;
      
  //P = [eye(d); 3*eye(d)]
  double P[][] = new double[d + d][d];
  for(int j = 0; j < d ; j++)
  {
    P[j][j] = 1.0;
    P[j + d][j] = 3.0;
  }
  
  //ma = P*m;    madm = P
  double ma[] = matrixMultiplyC(P, m);
  double madm[][] = P;//coeffProdMat(1.0, P);
  
  //va = P*v*P'; vadv = kron(P,P); va = (va+va')/2;
  double va[][] = matrixMultiply(P, v);
    va = matrixMultiply(va, transMat(P));
  double vadv[][] = kron(P, P);
  va = matrixADD(va, transMat(va));
  va = coeffProdMat(0.5, va);
  
  //[M2, S2, C2, Mdma, Sdma, Cdma, Mdva, Sdva, Cdva] = gSin(ma, va, [i d+i], [9*e e]/8);
  
  int ind[] = {i, d + i};
  double ee[] = {9 * e / 8.0, e / 8.0};
  gSinT = gSin(ma, va, ind, ee);
  
  //P = [eye(I) eye(I)] 
  double P2[] = {1, 1};
  
  //Q = [eye(d) 3*eye(d)]
  double Q[][] = new double [d][2 * d];
  for(int j = 0; j < d ; j++)
  {
    Q[j][j] = 1;
    Q[j][d + j] = 3;
  }
  
  //M = P*M2  
  ret.M = matrixMultiplyC_V(P2, gSinT.M);
  
  //S = P*S2*P'; S = (S+S')/2;    
  ret.S = matrixMultiplyC_V(matrixMultiplyC_T2(P2, gSinT.V), P2);
  ret.S = (ret.S + ret.S) / 2.0; // Необходимо будет в случае нескольких углов
  
  //C = Q*C2*P'
  ret.C = matrixMultiplyC(matrixMultiply(Q, gSinT.C), P2);
  
  return ret;
  
}

nargout_gSat_t gSatd(double m[], double v[][], int i, double e)
{
  
  nargout_gSat_t ret = new nargout_gSat_t();
  nargout_gSin_t gSinT = new nargout_gSin_t();
  
  //d = length(m); I = length(i); i = i(:)';
  //if nargin < 4; e = ones(1, I); end; e = e(:)';
  int d = m.length,
      I = 1;
      
  //P = [eye(d); 3*eye(d)]
  double P[][] = new double[d + d][d];
  for(int j = 0; j < d ; j++)
  {
    P[j][j] = 1.0;
    P[j + d][j] = 3.0;
  }
  
  //ma = P*m;    madm = P
  double ma[] = matrixMultiplyC(P, m);
  double madm[][] = P;//coeffProdMat(1.0, P);
  
  //va = P*v*P'; vadv = kron(P,P); va = (va+va')/2;
  double va[][] = matrixMultiply(P, v);
    va = matrixMultiply(va, transMat(P));
  double vadv[][] = kron(P, P);
  va = matrixADD(va, transMat(va));
  va = coeffProdMat(0.5, va);
  
  //[M2, S2, C2, Mdma, Sdma, Cdma, Mdva, Sdva, Cdva] = gSin(ma, va, [i d+i], [9*e e]/8);
  
  int ind[] = {i, d + i};
  double ee[] = {9 * e / 8.0, e / 8.0};
  gSinT = gSin(ma, va, ind, ee);
  
  //P = [eye(I) eye(I)] 
  double P2[] = {1, 1};
  
  //Q = [eye(d) 3*eye(d)]
  double Q[][] = new double [d][2 * d];
  for(int j = 0; j < d ; j++)
  {
    Q[j][j] = 1;
    Q[j][d + j] = 3;
  }
  
  //M = P*M2  
  ret.M = matrixMultiplyC_V(P2, gSinT.M);
  
  //S = P*S2*P'; S = (S+S')/2;    
  ret.S = matrixMultiplyC_V(matrixMultiplyC_T2(P2, gSinT.V), P2);
  ret.S = (ret.S + ret.S) / 2.0; // Необходимо будет в случае нескольких углов
  
  //C = Q*C2*P'
  ret.C = matrixMultiplyC(matrixMultiply(Q, gSinT.C), P2);
  
  //dMdm = P*Mdma*madm
  ret.dMdm = matrixMultiplyC_T2(matrixMultiplyC_T2(P2, gSinT.dMdm), madm);
  
  //dMdv = P*Mdva*vadv
  ret.dMdv = matrixMultiplyC_T2(matrixMultiplyC_T2(P2, gSinT.dMdv), vadv);
  
  //dSdm = kron(P,P)*Sdma*madm
  double kronP[] = {1, 1, 1, 1};
  ret.dSdm = matrixMultiplyC_T2(matrixMultiplyC_T2(kronP, gSinT.dVdm), madm);
  
  //dSdv = kron(P,P)*Sdva*vadv
  ret.dSdv = matrixMultiplyC_T2(matrixMultiplyC_T2(kronP, gSinT.dVdv), vadv);
  
  // dCdm = kron(P,Q)*Cdma*madm
  double P2P[][] = {
                    {1, 1}
                  };
  ret.dCdm = matrixMultiply(matrixMultiply(kron(P2P, Q), gSinT.dCdm), madm);
  
  //dCdv = kron(P,Q)*Cdva*vadv
  ret.dCdv = matrixMultiply(matrixMultiply(kron(P2P, Q), gSinT.dCdv), vadv);
  
  return ret;
  
}

//function [M, S, C, dMdm, dSdm, dCdm, dMds, dSds, dCds, dMdp, dSdp, dCdp] = conCat(con, sat, policy, m, s)

nargout_conpols_t conCat(policy_t _policy, double m[], double s[][])
{
  
  nargout_conpols_t ret = new nargout_conpols_t();
  nargout_conpols_t con = new nargout_conpols_t();
  nargout_gSat_t  sat = new nargout_gSat_t();
  
  //maxU=policy.maxU
  double maxU =_policy.maxU;
  
  //E=length(maxU)
  //D=length(m)
  int   E     = 1,
        D     = m.length;
        
  //F=D+E; j=D+1:F; i=1:D
  int F       = D + E,
      j       = F,
      i[]     = new int [D];
  for(int i1 = 0; i1 < D; i1++) 
     i[i1] = i1;
     
  //M = zeros(F,1); M(i) = m; S = zeros(F); S(i,i) = s
  double M[]   = new double [F],
         S[][] = new double [F][F];
        
  for(int i1 = 0; i1 < m.length; i1++) 
    M[i1] = m[i1];
  
  for(int i1 = 0; i1 < s.length; i1++) 
  for(int i2 = 0; i2 < s[0].length; i2++) 
    S[i1][i2] = s[i1][i2];
    
  //[M(j), S(j,j), Q] = con(policy, m, s);
  con = conpols(_policy, m, s);
  M[j - 1] = con.M[0];
  S[j - 1][j - 1] = con.S[0][0];
  double Q[][] = con.C;//coeffProdMat(1.0, con.C);

  //q = S(i,i)*Q; S(i,j) = q; S(j,i) = q'
  double S_i_i[][] = new double [i.length][i.length];
  for(int i1 = 0; i1 < i.length; i1++) 
  for(int i2 = 0; i2 < i.length; i2++) 
    S_i_i[i1][i2] = S[i1][i2];
  double q[][] = matrixMultiply(S_i_i, Q);
  for(int i1 = 0; i1 < i.length; i1++)
  {
    S[i1][j - 1] = q[i1][0];
    S[j - 1][i1] = q[i1][0];
  }
  
  //[M, S, R] = sat(M, S, j, maxU)
  sat = gSat(M, S, j - 1, maxU);
  
  ret.M = new double[1];
  ret.S = new double[1][1];
  ret.M[0]    = sat.M;
  ret.S[0][0] = sat.S;
  
  //C = [eye(D) Q]*R
  double C_[][] = new double[D][D + 1];
  for(int i1 = 0; i1 < D; i1++) 
   C_[i1][i1] = 1.0;
  for(int i1 = 0; i1 < Q.length; i1++) 
   C_[i1][D] = Q[i1][0];
  
  
  double C2[] = matrixMultiplyC(C_, sat.C);
  ret.C = new double[C2.length][1];
  for(int i1 = 0; i1 < C2.length; i1++) 
    ret.C[i1][0] = C2[i1];

  return ret;
  
}

nargout_conpolsd_t conCatd(policy_t _policy, double m[], double s[][])
{
  
  nargout_conpolsd_t ret = new nargout_conpolsd_t();
  nargout_conpolsd_t con = new nargout_conpolsd_t();
  nargout_gSat_t  sat = new nargout_gSat_t();
  
  //maxU=policy.maxU
  double maxU =_policy.maxU;
  
  //E=length(maxU)
  //D=length(m)
  int   E     = 1,
        D     = m.length;
        
  //F=D+E; j=D+1:F; i=1:D
  int F       = D + E,
      j       = F,
      i[]     = new int [D];
  for(int i1 = 0; i1 < D; i1++) 
     i[i1] = i1;
     
  //M = zeros(F,1); M(i) = m; S = zeros(F); S(i,i) = s
  double M[]   = new double [F],
         S[][] = new double [F][F];
        
  for(int i1 = 0; i1 < m.length; i1++) 
    M[i1] = m[i1];
  
  for(int i1 = 0; i1 < s.length; i1++) 
  for(int i2 = 0; i2 < s[0].length; i2++) 
    S[i1][i2] = s[i1][i2];
  
  //Mdm = zeros(F,D); Sdm = zeros(F*F,D); Mdm(1:D,1:D) = eye(D)
  //Mds = zeros(F,D*D); Sds = kron(Mdm,Mdm)  
  
  double Mdm[][] = new double [F][D],
         Sdm[][] = new double [F * F][D],
         Mds[][] = new double [F][D * D],
         Sds[][];
   
  for(int i1 = 0; i1 < D; i1++) 
    Mdm[i1][i1] = 1.0;
    
  Sds = kron(Mdm, Mdm);

  //X = reshape(1:F*F,[F F]); XT = X'
  
  int X[][] = new int [F][F];
  
  int count = 0;
  for(int i1 = 0; i1 < F; i1++) 
  for(int i2 = 0; i2 < F; i2++) 
    X[i2][i1] = count++;
  
  int XT[][] = transMat(X);

  //I=0*X;I(j,j)=1; jj=X(I==1)'; I=0*X;I(i,j)=1;ij=X(I==1)'; ji=XT(I==1)'
  
  int jj = X[j - 1][j - 1];
  int  ij[] = new int [i.length],
       ji[] = new int [i.length];
       
  for(int i1 = 0; i1 < i.length; i1++) 
  {
    ji[i1] = XT[i[i1]][j - 1];
    ij[i1] = X[i[i1]][j - 1];
  }
  
  //1. Unsquashed controller --------------------------------------------------
  //[M(j), S(j,j), Q, Mdm(j,:), Sdm(jj,:), dQdm, Mds(j,:), Sds(jj,:), dQds, Mdp, Sdp, dQdp] = con(policy, m, s);
  con = conpolsd(_policy, m, s);
  M[j - 1] = con.M[0];
  S[j - 1][j - 1] = con.S[0][0];

  double Q[][] = con.C;//coeffProdMat(1.0, con.C);
  
  for(int i1 = 0; i1 < con.dMdm[0].length; i1++) 
    Mdm[j - 1][i1] = con.dMdm[0][i1];
    
  for(int i1 = 0; i1 < con.dSdm[0].length; i1++) 
    Sdm[jj][i1] = con.dSdm[0][i1];
    
  double dQdm[][] = con.dCdm;//coeffProdMat(1.0, con.dCdm);
  
  for(int i1 = 0; i1 < con.dMds[0].length; i1++) 
    Mds[j - 1][i1] = con.dMds[0][i1];
    
  for(int i1 = 0; i1 < con.dMds[0].length; i1++) 
    Sds[jj][i1] = con.dSds[0][i1];

  double dQds[][] = con.dCds;//coeffProdMat(1.0, con.dCds);
  
  double Mdp[] = con.dMdp;//vectorCoef(con.dMdp, 1.0);
  
  double Sdp[] = con.dSdp;//vectorCoef(con.dSdp, 1.0);
  
  double dQdp[][] = con.dCdp;//coeffProdMat(1.0, con.dCdp); 
 
  //q = S(i,i)*Q; S(i,j) = q; S(j,i) = q'
  double S_i_i[][] = new double [i.length][i.length];
  for(int i1 = 0; i1 < i.length; i1++) 
  for(int i2 = 0; i2 < i.length; i2++) 
    S_i_i[i1][i2] = S[i1][i2];
  double q[][] = matrixMultiply(S_i_i, Q);
  for(int i1 = 0; i1 < i.length; i1++)
  {
    S[i1][j - 1] = q[i1][0];
    S[j - 1][i1] = q[i1][0];
  }

  // update the derivatives
  //SS = kron(eye(E),S(i,i))
  double SS[][] = kron(eye(E), S_i_i);
  
  //QQ = kron(Q',eye(D))
  double QQ[][] = kron(transMat(Q), eye(D));
  
  //Sdm(ij,:) = SS*dQdm;      Sdm(ji,:) = Sdm(ij,:) 
  double SS_dQ[][] = matrixMultiply(SS, dQdm);
  for(int i1 = 0; i1 < ij.length; i1++) 
  for(int i2 = 0; i2 < SS_dQ[0].length; i2++)
  {
    Sdm[ij[i1]][i2] = SS_dQ[i1][i2];
    Sdm[ji[i1]][i2] = SS_dQ[i1][i2];
  }
  
  // Sds(ij,:) = SS*dQds + QQ; Sds(ji,:) = Sds(ij,:)
  double SS_dQ2[][] = matrixMultiply(SS, dQds);
  SS_dQ2 = matrixADD(SS_dQ2, QQ);
  
  for(int i1 = 0; i1 < ij.length; i1++) 
  for(int i2 = 0; i2 < SS_dQ2[0].length; i2++)
  {
    Sds[ij[i1]][i2] = SS_dQ2[i1][i2];
    Sds[ji[i1]][i2] = SS_dQ2[i1][i2];
  }
 
  //2. Apply Saturation -------------------------------------------------------
  //[M, S, R, MdM, SdM, RdM, MdS, SdS, RdS] = sat(M, S, j, maxU);
  
  sat = gSatd(M, S, j - 1, maxU);
  
  ret.M = new double[1];
  ret.S = new double[1][1];
  ret.M[0]    = sat.M;
  ret.S[0][0] = sat.S;

  //apply chain-rule to compute derivatives after concatenation
  //dMdm = MdM*Mdm + MdS*Sdm; dMds = MdM*Mds + MdS*Sds;
  ret.dMdm = Chainrule(sat.dMdm, Mdm, sat.dMdv, Sdm);
  ret.dMds = Chainrule(sat.dMdm, Mds, sat.dMdv, Sds);
  
  //printMat(ret.dMdm);
  //printMat(ret.dMds);
  //dSdm = SdM*Mdm + SdS*Sdm; dSds = SdM*Mds + SdS*Sds
  ret.dSdm = Chainrule(sat.dSdm, Mdm, sat.dSdv, Sdm);
  ret.dSds = Chainrule(sat.dSdm, Mds, sat.dSdv, Sds);
  
  //dRdm = RdM*Mdm + RdS*Sdm; dRds = RdM*Mds + RdS*Sds
  double dRdm[][] = Chainrule(sat.dCdm, Mdm, sat.dCdv, Sdm);
  double dRds[][] = Chainrule(sat.dCdm, Mds, sat.dCdv, Sds);
  
  //dMdp = MdM(:,j)*Mdp + MdS(:,jj)*Sdp;
  //dSdp = SdM(:,j)*Mdp + SdS(:,jj)*Sdp;
  ret.dMdp = vectorADD(vectorCoef(con.dMdp, sat.dMdm[j - 1]), vectorCoef(con.dSdp, sat.dMdv[jj]));
  ret.dSdp = vectorADD(vectorCoef(con.dMdp, sat.dSdm[j - 1]), vectorCoef(con.dSdp, sat.dSdv[jj]));
  
  //dRdp = RdM(:,j)*Mdp + RdS(:,jj)*Sdp;
  double dRdp[][] = matrixADD(Outer_Product(Mat_to_Vec(sat.dCdm, j - 1), con.dMdp), 
                          Outer_Product(Mat_to_Vec(sat.dCdv, jj), con.dSdp));
                          
  //C = [eye(D) Q]*R inv(s)*cov(x,u)
  double C_[][] = new double[D][D + 1];
  for(int i1 = 0; i1 < D; i1++) 
   C_[i1][i1] = 1.0;
  for(int i1 = 0; i1 < Q.length; i1++) 
   C_[i1][D] = Q[i1][0];
  
  
  double C2[] = matrixMultiplyC(C_, sat.C);
  ret.C = new double[C2.length][1];
  for(int i1 = 0; i1 < C2.length; i1++) 
    ret.C[i1][0] = C2[i1];
  
  //update the derivatives
  //RR = kron(R(j,:)',eye(D)); QQ = kron(eye(E),[eye(D) Q]);
  double R[][] = {{sat.C[j - 1]}};
  double RR[][] = kron(R, eye(D));
  double QQ2[][] = kron(eye(E), C_);
  
  //dCdm = QQ*dRdm + RR*dQdm;
  //dCds = QQ*dRds + RR*dQds;
  //dCdp = QQ*dRdp + RR*dQdp;
  ret.dCdm = Chainrule(QQ2, dRdm, RR, dQdm);
  ret.dCds = Chainrule(QQ2, dRds, RR, dQds);
  ret.dCdp = Chainrule(QQ2, dRdp, RR, dQdp);
 
 
  return ret;
  
}

//function [Mnext, Snext] = propagate(m, s, plant, dynmodel, policy)

nargout_propagate_t propagate(double m[], double s[][], dynmodel_t _gpmodel, policy_t _policy)
{
 
 nargout_propagate_t ret = new nargout_propagate_t();
 gTrig_r             gTr = new gTrig_r();
 
 //D0 = length(m);                                        % size of the input mean
 //D1 = D0 + 2*length(angi);          % length after mapping all angles to sin/cos
 //D2 = D1 + length(policy.maxU);          % length after computing control signal
 //D3 = D2 + D0;                                         % length after predicting
 //M = zeros(D3,1); M(1:D0) = m; S = zeros(D3); S(1:D0,1:D0) = s;   % init M and S
 
 int D0 = m.length,
     D1 = D0 + 2 * _policy.angle.length,
     D2 = D1 + 1,
     D3 = D2 + D0;
 
 //M = zeros(D3,1); M(1:D0) = m; S = zeros(D3); S(1:D0,1:D0) = s   % init M and S
 double M[]   = new double [D3],
        S[][] = new double [D3][D3];
 
 for(int i1 = 0; i1 < D0; i1++) 
   M[i1] = m[i1];
    
 for(int i1 = 0; i1 < D0; i1++) 
  for(int i2 = 0; i2 < D0; i2++) 
    S[i1][i2] = s[i1][i2];
  
 //1) Augment state distribution with trigonometric functions
 //i = 1:D0; j = 1:D0; k = D0+1:D1;
 //[M(k), S(k,k) C] = gTrig(M(i), S(i,i), angi);
 
 int i[] = createInd(0, D0),
     j[] = createInd(0, D0),
     k[] = createInd(D0, D1);

 //[M(k), S(k,k) C] = gTrig(M(i), S(i,i), angi);
 
 gTr = gTrig(getVec(M, i), getMat(S, i, i), _policy.angle, 1);
 
 M = setVec(M, k, gTr.M);
 S = setMat(S, k, k, gTr.V);
 double C[][] = gTr.C;//coeffProdMat(1.0, gTr.C); 
 
 //q = S(j,i)*C; S(j,k) = q; S(k,j) = q'
 S = addTo(S, i, j, k, C);
 
 //sn2 = exp(2*dynmodel.hyp(end,:)); sn2(difi) = sn2(difi)/2
 double sn2[] = vectorCoef(Mat_to_VecT(_gpmodel.hyp, _gpmodel.hyp.length - 1), 2.0); 
 sn2 = EXP_F(sn2);
 for(int i1 = 0; i1 < difi.length; i1++)
  sn2[difi[i1]] = sn2[difi[i1]] / 2.0;
   
 //mm=zeros(D1,1); mm(i)=M(i); ss(i,i)=S(i,i)+diag(sn2)
 double mm[]   = new double [D1],
        ss[][] = new double [D1][D1];
 
 mm = setVec(mm, i, M);
 ss = setMat(ss, i, i, getMat(S, i, i));
 for(int i1 = 0; i1 < sn2.length; i1++) 
   ss[i1][i1] += sn2[i1];

 //[mm(k), ss(k,k) C] = gTrig(mm(i), ss(i,i), angi)     % noisy state measurement
 
 gTr = gTrig(getVec(mm, i), getMat(ss, i, i), _policy.angle, 1);
 mm = setVec(mm, k, gTr.M);
 ss = setMat(ss, k, k, gTr.V);
 C = null;
 C = gTr.C;//coeffProdMat(1.0, gTr.C); 
 
 //q = ss(j,i)*C; ss(j,k) = q; ss(k,j) = q'
 ss = addTo(ss, i, j, k, C); 

 //2) Compute distribution of the control signal
 //i = poli; j = 1:D1; k = D1+1:D2;
 i = null; j = null; k = null;
 i = _policy.poli;//vectorCoef(_policy.poli, 1);
 j = createInd(0, D1);
 k = createInd(D1, D2);

 //[M(k) S(k,k) C] = policy.fcn(policy, mm(i), ss(i,i));
 nargout_conpols_t cnCt = new nargout_conpols_t();
 cnCt = conCat(_policy, getVec(mm, i), getMat(ss, i, i));
 
 M = setVec(M, k, cnCt.M);
 S = setMat(S, k, k, cnCt.S);
 C = null;
 C = cnCt.C;//coeffProdMat(1.0, cnCt.C); 
  
 //q = S(j,i)*C; S(j,k) = q; S(k,j) = q'
 S = addTo(S, i, j, k, C);

 //3) Compute dynamics-GP prediction 
 //ii = [dyni D1+1:D2]; j = 1:D2
 int ii[] = new int [_gpmodel.dyni.length + D2 - D1];
 for(int i1 = 0; i1 < _gpmodel.dyni.length; i1++) 
  ii[i1] = _gpmodel.dyni[i1];
 for(int i1 = 0; i1 < (D2 - D1); i1++) 
  ii[i1 + _gpmodel.dyni.length] = i1 + D1;

 //dyn = dynmodel; k = D2+1:D3; i = ii
 j = null;
 j = createInd(0, D2);
 k = null;
 k = createInd(D2, D3);
 i = vectorCoef(ii, 1);
 
 //j = setdiff(j,k) В данном случае не требуется в массиве k нет индексов из j
 //[M(k), S(k,k), C] = dyn.fcn(dyn, M(i), S(i,i));
 // q = S(j,i)*C; S(j,k) = q; S(k,j) = q';

 nargout_gp0_t gp = new nargout_gp0_t();
 
 _gpmodel = gp0(_gpmodel);
 gp = _gp0(_gpmodel, getVec(M, i), getMat(S, i, i));
 
 M = setVec(M, k, gp.M);
 S = setMat(S, k, k, gp.S);
 C = null;
 C = gp.V;//coeffProdMat(1.0, gp.V); 

 //q = S(j,i)*C; S(j,k) = q; S(k,j) = q'
 S = addTo(S, i, j, k, C);
 
 // 4) Compute distribution of the next state
 //P = [zeros(D0,D2) eye(D0)]; P(difi,difi) = eye(length(difi))
 double P[][] = new double [D0][D0 + D2];
 for(int i1 = 0; i1 < D0; i1++) 
 {
   P[i1][i1 + D2] = 1.0;
 }
 
 for(int i1 = 0; i1 < difi.length; i1++) 
 {
   P[difi[i1]][difi[i1]] = 1.0;
 }

 //Mnext = P*M
 ret.Mnext = matrixMultiplyC(P, M); 
 
 //Snext = P*S*P'; Snext = (Snext+Snext')/2
 ret.Snext = matrixMultiply(P, S);
 ret.Snext = matrixMultiply(ret.Snext, transMat(P));
 ret.Snext = matrixADD(ret.Snext, transMat(ret.Snext));
 ret.Snext = coeffProdMat(0.5, ret.Snext);
 
 return ret;
 
}

//function [Mnext, Snext, dMdm, dSdm, dMds, dSds, dMdp, dSdp] = propagated(m, s, plant, dynmodel, policy)

nargout_propagated_t propagated(double m[], double s[][], dynmodel_t _gpmodel, policy_t _policy)
{
 
 nargout_propagated_t ret = new nargout_propagated_t();
 gTrig_Full_r             gTr   = new gTrig_Full_r();
 gTrig_r                  gTr2  = new gTrig_r();
 
 
 //D0 = length(m);                                        % size of the input mean
 //D1 = D0 + 2*length(angi);          % length after mapping all angles to sin/cos
 //D2 = D1 + length(policy.maxU);          % length after computing control signal
 //D3 = D2 + D0;                                         % length after predicting
 //M = zeros(D3,1); M(1:D0) = m; S = zeros(D3); S(1:D0,1:D0) = s;   % init M and S
 
 int D0 = m.length, //<>// //<>// //<>//
     D1 = D0 + 2 * _policy.angle.length,
     D2 = D1 + 1,
     D3 = D2 + D0;
 
 //M = zeros(D3,1); M(1:D0) = m; S = zeros(D3); S(1:D0,1:D0) = s   % init M and S
 double M[]   = new double [D3],
        S[][] = new double [D3][D3];
 
 for(int i1 = 0; i1 < D0; i1++) 
   M[i1] = m[i1];
    
 for(int i1 = 0; i1 < D0; i1++) 
  for(int i2 = 0; i2 < D0; i2++) 
    S[i1][i2] = s[i1][i2];
  

 //Mdm = [eye(D0); zeros(D3-D0,D0)]; Sdm = zeros(D3*D3,D0)
 double Mdm[][] = new double [D3][D0];
 for(int i1 = 0; i1 < D0; i1++) 
   Mdm[i1][i1] = 1.0;
 
 double Sdm[][] = new double [D3 * D3][D0]; 
 
 //Mds = zeros(D3,D0*D0); Sds = kron(Mdm,Mdm)
 double Mds[][] = new double [D3][D0 * D0]; 
 double Sds[][] = kron(Mdm,Mdm);
 
 //X = reshape(1:D3*D3,[D3 D3]); XT = X'; Sds = (Sds + Sds(XT(:),:))/2; 
  Sds = symmetrizeA(D3, Sds);
  
  //X = reshape(1:D0*D0,[D0 D0]); XT = X'; Sds = (Sds + Sds(:,XT(:)))/2
  Sds = symmetrizeB(D0, Sds);

  // 1) Augment state distribution with trigonometric functions
  //i = 1:D0; j = 1:D0; k = D0+1:D1
  int i[] = createInd(0, D0),
      j[] = createInd(0, D0),
      k[] = createInd(D0, D1);

 //[M(k), S(k,k) C] = gTrig(M(i), S(i,i), angi);
 
 gTr = gTrigF(getVec(M, i), getMat(S, i, i), _policy.angle, 1);
 M = setVec(M, k, gTr.M);
 S = setMat(S, k, k, gTr.V);
 double C[][] = gTr.C;//coeffProdMat(1.0, gTr.C); 
 
 //[S Mdm Mds Sdm Sds] = fillIn(S,C,mdm,sdm,Cdm,mds,sds,Cds,Mdm,Sdm,Mds,Sds,[ ],[ ],[ ],i,j,k,D3)
 nargin_fillInSh_t inp = new nargin_fillInSh_t();
 inp.S = S;   //coeffProdMat(1.0, S);
 inp.C = C;   //coeffProdMat(1.0, C);
 inp.mdm =gTr.dMdm; //coeffProdMat(1.0, gTr.dMdm);
 inp.sdm =gTr.dVdm; //coeffProdMat(1.0, gTr.dVdm);
 inp.Cdm =gTr.dCdm; //coeffProdMat(1.0, gTr.dCdm);
 inp.mds =gTr.dMdv; //coeffProdMat(1.0, gTr.dMdv);
 inp.sds = gTr.dVdv;//coeffProdMat(1.0, gTr.dVdv);
 inp.Cds =gTr.dCdv; //coeffProdMat(1.0, gTr.dCdv);
 inp.Mdm =Mdm; //coeffProdMat(1.0, Mdm);
 inp.Sdm =Sdm; //coeffProdMat(1.0, Sdm);
 inp.Mds =Mds; //coeffProdMat(1.0, Mds);
 inp.Sds =Sds; //coeffProdMat(1.0, Sds);
 inp.i = i;//vectorCoef(i, 1); 
 inp.j = j;//vectorCoef(j, 1); 
 inp.k = k;//vectorCoef(k, 1);
 inp.D = D3;
 
 nargout_fillInSh_t fl = fillIn(inp);
 
 S =fl.S; //coeffProdMat(1.0, fl.S);
 Mdm =fl.Mdm;// coeffProdMat(1.0, fl.Mdm); 
 Mds =fl.Mds; //coeffProdMat(1.0, fl.Mds); 
 Sdm =fl.Sdm; //coeffProdMat(1.0, fl.Sdm); 
 Sds =fl.Sds; //coeffProdMat(1.0, fl.Sds);

 //sn2 = exp(2*dynmodel.hyp(end,:)); sn2(difi) = sn2(difi)/2
 double sn2[] = vectorCoef(Mat_to_VecT(_gpmodel.hyp, _gpmodel.hyp.length - 1), 2.0); 
 sn2 = EXP_F(sn2);
 for(int i1 = 0; i1 < difi.length; i1++)
  sn2[difi[i1]] = sn2[difi[i1]] / 2.0;
 

 //mm=zeros(D1,1); mm(i)=M(i); ss(i,i)=S(i,i)+diag(sn2)
 double mm[]   = new double [D1],
        ss[][] = new double [D1][D1];
 
 mm = setVec(mm, i, M);
 ss = setMat(ss, i, i, getMat(S, i, i));
 for(int i1 = 0; i1 < sn2.length; i1++) 
   ss[i1][i1] += sn2[i1];

 //[mm(k), ss(k,k) C] = gTrig(mm(i), ss(i,i), angi)     % noisy state measurement
 
 gTr2 = gTrig(getVec(mm, i), getMat(ss, i, i), _policy.angle, 1);
 mm = setVec(mm, k, gTr2.M);
 ss = setMat(ss, k, k, gTr2.V);
 C = null;
 C = gTr2.C;//coeffProdMat(1.0, gTr2.C); 

 //q = ss(j,i)*C; ss(j,k) = q; ss(k,j) = q'
 ss = addTo(ss, i, j, k, C);
 
 //2) Compute distribution of the control signal
 //i = poli; j = 1:D1; k = D1+1:D2
 i = null; j = null; k = null;
 i = _policy.poli;//vectorCoef(_policy.poli, 1);
 j = createInd(0, D1);
 k = createInd(D1, D2);

 //[M(k) S(k,k) C mdm sdm Cdm mds sds Cds Mdp Sdp Cdp] = policy.fcn(policy, mm(i), ss(i,i))
 nargout_conpolsd_t cnCt = new nargout_conpolsd_t();
 cnCt = conCatd(_policy, getVec(mm, i), getMat(ss, i, i));
 
 M = setVec(M, k, cnCt.M);
 S = setMat(S, k, k, cnCt.S);
 C = null;
 C = cnCt.C;//coeffProdMat(1.0, cnCt.C); 
 
 //[S Mdm Mds Sdm Sds Mdp Sdp] = fillIn(S,C,mdm,sdm,Cdm,mds,sds,Cds,Mdm,Sdm,Mds,Sds,Mdp,Sdp,Cdp,i,j,k,D3)
 nargin_fillInFl_t inp2 = new nargin_fillInFl_t();
 inp2.S = S;  //coeffProdMat(1.0, S);
 inp2.C = C;  //coeffProdMat(1.0, C);
 inp2.mdm =cnCt.dMdm;// coeffProdMat(1.0, cnCt.dMdm);
 inp2.sdm =cnCt.dSdm;// coeffProdMat(1.0, cnCt.dSdm);
 inp2.Cdm = cnCt.dCdm;// coeffProdMat(1.0, cnCt.dCdm);
 inp2.mds =cnCt.dMds;// coeffProdMat(1.0, cnCt.dMds);
 inp2.sds =cnCt.dSds;// coeffProdMat(1.0, cnCt.dSds);
 inp2.Cds =cnCt.dCds;// coeffProdMat(1.0, cnCt.dCds);
 inp2.Mdm =Mdm;// coeffProdMat(1.0, Mdm);
 inp2.Sdm =Sdm;// coeffProdMat(1.0, Sdm);
 inp2.Mds =Mds;// coeffProdMat(1.0, Mds);
 inp2.Sds = Sds;// coeffProdMat(1.0, Sds);
 
 inp2.Mdp = coeffProdMat(1.0, cnCt.dMdp);
 inp2.Sdp = coeffProdMat(1.0, cnCt.dSdp);
 inp2.Cdp = cnCt.dCdp;//coeffProdMat(1.0, cnCt.dCdp);
 
 inp2.i = i;//vectorCoef(i, 1); 
 inp2.j = j;//vectorCoef(j, 1); 
 inp2.k = k;//vectorCoef(k, 1);
 inp2.D = D3;
 
 nargout_fillInFl_t fl2 = fillInFl(inp2);
  
 S =fl2.S;// coeffProdMat(1.0, fl2.S);
 Mdm = fl2.Mdm;//coeffProdMat(1.0, fl2.Mdm); 
 Mds = fl2.Mds;//coeffProdMat(1.0, fl2.Mds); 
 Sdm =fl2.Sdm;// coeffProdMat(1.0, fl2.Sdm); 
 Sds = fl2.Sds;//coeffProdMat(1.0, fl2.Sds);
 double Mdp[][] = fl2.Mdp;//coeffProdMat(1.0, fl2.Mdp);
 double Sdp[][] = fl2.Sdp;//coeffProdMat(1.0, fl2.Sdp);

 //3) Compute distribution of the change in state
 //ii = [dyni D1+1:D2]; j = 1:D2
 int ii[] = new int [_gpmodel.dyni.length + D2 - D1];
 for(int i1 = 0; i1 < _gpmodel.dyni.length; i1++) 
  ii[i1] = _gpmodel.dyni[i1];
 for(int i1 = 0; i1 < (D2 - D1); i1++) 
  ii[i1 + _gpmodel.dyni.length] = i1 + D1;
 
 //dyn = dynmodel; k = D2+1:D3; i = ii
 j = null;
 j = createInd(0, D2);
 k = null;
 k = createInd(D2, D3);
 i = ii;//vectorCoef(ii, 1);

 //[M(k) S(k,k) C mdm sdm Cdm mds sds Cds] = dyn.fcn(dyn, M(i), S(i,i))
 nargout_gp0d_t gpd = new nargout_gp0d_t();
 
 //_gpmodel = gp0(_gpmodel);
 gpd = _gp0d(_gpmodel, getVec(M, i), getMat(S, i, i));
 
 M = setVec(M, k, gpd.M);
 S = setMat(S, k, k, gpd.S);
 C = null;
 C = gpd.V;//coeffProdMat(1.0, gpd.V); 

 //[S Mdm Mds Sdm Sds Mdp Sdp] = fillIn(S,C,mdm,sdm,Cdm,mds,sds,Cds,Mdm,Sdm,Mds,Sds,Mdp,Sdp,[ ],i,j,k,D3)
 nargin_fillInMd_t inp3 = new nargin_fillInMd_t();
 inp3.S =S;// coeffProdMat(1.0, S);
 inp3.C =C;// coeffProdMat(1.0, C);
 inp3.mdm = gpd.dMdm;//coeffProdMat(1.0, gpd.dMdm);
 inp3.sdm = gpd.dSdm;//coeffProdMat(1.0, gpd.dSdm);
 inp3.Cdm =gpd.dVdm;// coeffProdMat(1.0, gpd.dVdm);
 inp3.mds =  gpd.dMds;//coeffProdMat(1.0, gpd.dMds);
 inp3.sds = gpd.dSds;//coeffProdMat(1.0, gpd.dSds);
 inp3.Cds =  gpd.dVds;//coeffProdMat(1.0, gpd.dVds);
 inp3.Mdm =Mdm;// coeffProdMat(1.0, Mdm);
 inp3.Sdm =Sdm;// coeffProdMat(1.0, Sdm);
 inp3.Mds =  Mds;//coeffProdMat(1.0, Mds);
 inp3.Sds = Sds;//coeffProdMat(1.0, Sds);
 
 inp3.Mdp = Mdp;// coeffProdMat(1.0, Mdp);
 inp3.Sdp = Sdp;//coeffProdMat(1.0, Sdp);
  
 inp3.i =i; //vectorCoef(i, 1); 
 inp3.j =j; //vectorCoef(j, 1); 
 inp3.k =k; //vectorCoef(k, 1);
 inp3.D = D3;
 
 nargout_fillInMd_t fl3 = fillInMd(inp3);
 
 S = fl3.S;// coeffProdMat(1.0, fl3.S);
 Mdm =fl3.Mdm;// coeffProdMat(1.0, fl3.Mdm); 
 Mds =fl3.Mds;// coeffProdMat(1.0, fl3.Mds); 
 Sdm = fl3.Sdm;// coeffProdMat(1.0, fl3.Sdm); 
 Sds = fl3.Sds;//coeffProdMat(1.0, fl3.Sds);
 Mdp = fl3.Mdp;//coeffProdMat(1.0, fl3.Mdp);
 Sdp = fl3.Sdp;//coeffProdMat(1.0, fl3.Sdp);

 //4) Compute distribution of the next state
 //P = [zeros(D0,D2) eye(D0)]; P(difi,difi) = eye(length(difi));  P = sparse( P)
 double P[][] = new double [D0][D0 + D2];
 for(int i1 = 0; i1 < D0; i1++) 
 {
   P[i1][i1 + D2] = 1.0;
 }
 
 for(int i1 = 0; i1 < difi.length; i1++) 
 {
   P[difi[i1]][difi[i1]] = 1.0;
 }
 
 //Mnext = P*M; Snext = P*S*P'; Snext = (Snext+Snext')/2

 ret.Mnext = matrixMultiplyC(P, M);
 ret.Snext = matrixMultiply(matrixMultiply(P, S), transMat(P));
 ret.Snext = coeffProdMat(0.5, matrixADD(ret.Snext, transMat(ret.Snext)));
 
 //PP = kron(P,P)
  double PP[][] = kron(P, P);
 
 //dMdm =  P*Mdm; dMds =  P*Mds; dMdp =  P*Mdp
 ret.dMdm = matrixMultiply(P, Mdm);
 ret.dMds = matrixMultiply(P, Mds);
 ret.dMdp = matrixMultiply(P, Mdp);
 
 //dSdm = PP*Sdm; dSds = PP*Sds; dSdp = PP*Sdp
 ret.dSdm = matrixMultiply(PP, Sdm);
 ret.dSds = matrixMultiply(PP, Sds);
 ret.dSdp = matrixMultiply(PP, Sdp);
 
 //sparse_t P_sp = new sparse_t(P);
  
 //ret.Mnext = matrixMultiplyC_SP(P_sp, M);//matrixMultiplyC(P, M);
 //ret.Snext = matrixMultiply_SP(matrixMultiply_SP(P_sp, S), transMat_SP(P_sp));//matrixMultiply(matrixMultiply(P, S), transMat(P));
 //ret.Snext = coeffProdMat(0.5, matrixADD(ret.Snext, transMat(ret.Snext)));
 
 ////PP = kron(P,P)
 //sparse_t PP_sp = new sparse_t(kron(P, P)); 
 ////double PP[][] = kron(P, P);
 
 ////dMdm =  P*Mdm; dMds =  P*Mds; dMdp =  P*Mdp
 //ret.dMdm = matrixMultiply_SP(P_sp, Mdm);//matrixMultiply(P, Mdm);
 //ret.dMds = matrixMultiply_SP(P_sp, Mds);//matrixMultiply(P, Mds);
 //ret.dMdp = matrixMultiply_SP(P_sp, Mdp);//matrixMultiply(P, Mdp);
 
 ////dSdm = PP*Sdm; dSds = PP*Sds; dSdp = PP*Sdp
 //ret.dSdm = matrixMultiply_SP(PP_sp, Sdm);//matrixMultiply(PP, Sdm);
 //ret.dSds = matrixMultiply_SP(PP_sp, Sds);//matrixMultiply(PP, Sds);
 //ret.dSdp = matrixMultiply_SP(PP_sp, Sdp);//matrixMultiply(PP, Sdp);
 
 //X = reshape(1:D0*D0,[D0 D0]); XT = X'     symmetrize dS
 //dSdm = (dSdm + dSdm(XT(:),:))/2
 ret.dSdm = symmetrizeA(D0, ret.dSdm);
 
 //dMds = (dMds + dMds(:,XT(:)))/2
 ret.dMds = symmetrizeB(D0, ret.dMds);
 
 //dSds = (dSds + dSds(XT(:),:))/2
 ret.dSds = symmetrizeA(D0, ret.dSds);
 
 //dSds = (dSds + dSds(:,XT(:)))/2
 ret.dSds = symmetrizeB(D0, ret.dSds);
 
 //dSdp = (dSdp + dSdp(XT(:),:))/2
 ret.dSdp = symmetrizeA(D0, ret.dSdp);

 return ret;
 
}

// Полная версия
//[S Mdm Mds Sdm Sds Mdp Sdp] = fillIn(S,C,mdm,sdm,Cdm,mds,sds,Cds,Mdm,Sdm,Mds,Sds,Mdp,Sdp,Cdp,i,j,k,D3)

//function [S Mdm Mds Sdm Sds Mdp Sdp] = fillIn(S,C,mdm,sdm,Cdm,mds,sds,Cds,Mdm,Sdm,Mds,Sds,Mdp,Sdp,dCdp,i,j,k,D)

nargout_fillInFl_t fillInFl(nargin_fillInFl_t inp)
{
  
  nargout_fillInFl_t ret = new nargout_fillInFl_t();
  
  //X = reshape(1:D*D,[D D]); XT = X'                         vectorized indices
    
  int X[][] = new int [inp.D][inp.D];
  
  int count = 0;
  for(int i1 = 0; i1 < inp.D; i1++) 
  for(int i2 = 0; i2 < inp.D; i2++) 
    X[i2][i1] = count++;
  
  int XT[][] = transMat(X);
  
  //I=0*X; I(i,i)=1; ii=X(I==1)'; 
  int ii[] = new int [inp.i.length * inp.i.length];
  int ct = 0;
  for(int i1 = 0; i1 < inp.i.length; i1++) 
  for(int i2 = 0; i2 < inp.i.length; i2++) 
    ii[ct++] = X[inp.i[i2]][inp.i[i1]];
  
  //I=0*X; I(k,k)=1; kk=X(I==1)'
  int kk[] = new int [inp.k.length * inp.k.length];
  ct = 0;
  for(int i1 = 0; i1 < inp.k.length; i1++) 
  for(int i2 = 0; i2 < inp.k.length; i2++) 
    kk[ct++] = X[inp.k[i2]][inp.k[i1]];
    
  //I=0*X; I(j,i)=1; ji=X(I==1)' 
  int ji[] = new int [inp.j.length * inp.i.length];
  ct = 0;
  for(int i1 = 0; i1 < inp.i.length; i1++) 
  for(int i2 = 0; i2 < inp.j.length; i2++) 
    ji[ct++] = X[inp.j[i2]][inp.i[i1]];
    
  //I=0*X; I(j,k)=1; jk=X(I==1)'; kj=XT(I==1)'
  int jk[] = new int [inp.j.length * inp.k.length];
  int kj[] = new int [inp.k.length * inp.j.length];
  ct = 0;
  for(int i1 = 0; i1 < inp.k.length; i1++) 
  for(int i2 = 0; i2 < inp.j.length; i2++) 
    jk[ct++] = X[inp.j[i2]][inp.k[i1]];
    
  ct = 0;
  for(int i1 = 0; i1 < inp.k.length; i1++) 
  for(int i2 = 0; i2 < inp.j.length; i2++) 
    kj[ct++] = XT[inp.j[i2]][inp.k[i1]];
  
  //Mdm(k,:)  = mdm*Mdm(i,:) + mds*Sdm(ii,:)
  inp.Mdm = setMat(inp.Mdm, inp.k, Chainrule(inp.mdm, inp.Mdm, inp.i, inp.mds, inp.Sdm, ii));
  ret.Mdm = inp.Mdm;//coeffProdMat(1.0, inp.Mdm);
  
  //Mds(k,:)  = mdm*Mds(i,:) + mds*Sds(ii,:)
  inp.Mds = setMat(inp.Mds, inp.k, Chainrule(inp.mdm, inp.Mds, inp.i, inp.mds, inp.Sds, ii));
  ret.Mds = inp.Mds;//coeffProdMat(1.0, inp.Mds);
  
  //Sdm(kk,:) = sdm*Mdm(i,:) + sds*Sdm(ii,:)
  inp.Sdm = setMat(inp.Sdm, kk, Chainrule(inp.sdm, inp.Mdm, inp.i, inp.sds, inp.Sdm, ii));
  ret.Sdm =  inp.Sdm;//coeffProdMat(1.0, inp.Sdm);
  
  //Sds(kk,:) = sdm*Mds(i,:) + sds*Sds(ii,:)
  inp.Sds = setMat(inp.Sds, kk, Chainrule(inp.sdm, inp.Mds, inp.i, inp.sds, inp.Sds, ii));
  ret.Sds =  inp.Sds;//coeffProdMat(1.0, inp.Sds);
  
  //dCdm      = Cdm*Mdm(i,:) + Cds*Sdm(ii,:)
  double dCdm[][] = Chainrule(inp.Cdm, inp.Mdm, inp.i, inp.Cds, inp.Sdm, ii);
  
  //dCds      = Cdm*Mds(i,:) + Cds*Sds(ii,:)
  double dCds[][] = Chainrule(inp.Cdm, inp.Mds, inp.i, inp.Cds, inp.Sds, ii);
  
  //aa = length(k); bb = aa^2; cc = numel(C)
  //mdp = zeros(D,size(Mdp,2)); sdp = zeros(D*D,size(Mdp,2))
  //mdp(k,:)  = reshape(Mdp,aa,[]); Mdp = mdp
  ret.Mdp = new double[inp.D][inp.Mdp[0].length];
  ret.Mdp = setMat(ret.Mdp, inp.k, inp.Mdp);
  
  //sdp(kk,:) = reshape(Sdp,bb,[]); Sdp = sdp
  ret.Sdp = new double[inp.D * inp.D][inp.Mdp[0].length];
  ret.Sdp = setMat(ret.Sdp, kk, inp.Sdp);
  
  //Cdp       = reshape(dCdp,cc,[]); dCdp = Cdp
  double dCdp[][] = coeffProdMat(1.0, inp.Cdp);

  //q = S(j,i)*C; S(j,k) = q; S(k,j) = q'
  double S[][] = addTo(inp.S, inp.i, inp.j, inp.k, inp.C);
  ret.S =S;// coeffProdMat(1.0, S);
  
  //SS = kron(eye(length(k)),S(j,i)); CC = kron(C',eye(length(j)))
  
  double SS[][] = kron(eye(inp.k.length), getMat(inp.S, inp.j, inp.i));
  double CC[][] = kron(transMat(inp.C), eye(inp.j.length));
  
  //Sdm(jk,:) = SS*dCdm + CC*Sdm(ji,:); Sdm(kj,:) = Sdm(jk,:)
  inp.Sdm = setMat(inp.Sdm, jk, Chainrule(SS, dCdm, CC, getMat(inp.Sdm, ji)));
  inp.Sdm = setMat(inp.Sdm, kj, getMat(inp.Sdm, jk));
  ret.Sdm =inp.Sdm;// coeffProdMat(1.0, inp.Sdm);
  
  //Sds(jk,:) = SS*dCds + CC*Sds(ji,:); Sds(kj,:) = Sds(jk,:)
  inp.Sds = setMat(inp.Sds, jk, Chainrule(SS, dCds, CC, getMat(inp.Sds, ji)));
  inp.Sds = setMat(inp.Sds, kj, getMat(inp.Sds, jk));
  ret.Sds = inp.Sds;//coeffProdMat(1.0, inp.Sds);
  
  //Sdp(jk,:) = SS*dCdp + CC*Sdp(ji,:); Sdp(kj,:) = Sdp(jk,:)
  
  ret.Sdp = setMat(ret.Sdp, jk, Chainrule(SS, dCdp, CC, getMat(ret.Sdp, ji)));
  ret.Sdp = setMat(ret.Sdp, kj, getMat(ret.Sdp, jk));
 
  return ret;
  
}
// Короткая версия
//[S Mdm Mds Sdm Sds] = fillIn(S,C,mdm,sdm,Cdm,mds,sds,Cds,Mdm,Sdm,Mds,Sds,[ ],[ ],[ ],i,j,k,D3)
//S,C,mdm,sdm,Cdm,mds,sds,Cds,Mdm,Sdm,Mds,Sds,[ ],[ ],[ ],i,j,k,D3

nargout_fillInSh_t fillIn(nargin_fillInSh_t inp)
{
  
  nargout_fillInSh_t ret = new nargout_fillInSh_t();
  
  //X = reshape(1:D*D,[D D]); XT = X'                         vectorized indices
    
  int X[][] = new int [inp.D][inp.D];
  
  int count = 0;
  for(int i1 = 0; i1 < inp.D; i1++) 
  for(int i2 = 0; i2 < inp.D; i2++) 
    X[i2][i1] = count++;
  
  int XT[][] = transMat(X);
  
  //I=0*X; I(i,i)=1; ii=X(I==1)'; 
  int ii[] = new int [inp.i.length * inp.i.length];
  int ct = 0;
  for(int i1 = 0; i1 < inp.i.length; i1++) 
  for(int i2 = 0; i2 < inp.i.length; i2++) 
    ii[ct++] = X[inp.i[i2]][inp.i[i1]];
  
  //I=0*X; I(k,k)=1; kk=X(I==1)'
  int kk[] = new int [inp.k.length * inp.k.length];
  ct = 0;
  for(int i1 = 0; i1 < inp.k.length; i1++) 
  for(int i2 = 0; i2 < inp.k.length; i2++) 
    kk[ct++] = X[inp.k[i2]][inp.k[i1]];
    
  //I=0*X; I(j,i)=1; ji=X(I==1)' 
  int ji[] = new int [inp.j.length * inp.i.length];
  ct = 0;
  for(int i1 = 0; i1 < inp.i.length; i1++) 
  for(int i2 = 0; i2 < inp.j.length; i2++) 
    ji[ct++] = X[inp.j[i2]][inp.i[i1]];
    
  //I=0*X; I(j,k)=1; jk=X(I==1)'; kj=XT(I==1)'
  int jk[] = new int [inp.j.length * inp.k.length];
  int kj[] = new int [inp.k.length * inp.j.length];
  ct = 0;
  for(int i1 = 0; i1 < inp.k.length; i1++) 
  for(int i2 = 0; i2 < inp.j.length; i2++) 
    jk[ct++] = X[inp.j[i2]][inp.k[i1]];
    
  ct = 0;
  for(int i1 = 0; i1 < inp.k.length; i1++) 
  for(int i2 = 0; i2 < inp.j.length; i2++) 
    kj[ct++] = XT[inp.j[i2]][inp.k[i1]];
  
  //Mdm(k,:)  = mdm*Mdm(i,:) + mds*Sdm(ii,:)
  inp.Mdm = setMat(inp.Mdm, inp.k, Chainrule(inp.mdm, inp.Mdm, inp.i, inp.mds, inp.Sdm, ii));
  ret.Mdm = inp.Mdm;//coeffProdMat(1.0, inp.Mdm);
  
  //Mds(k,:)  = mdm*Mds(i,:) + mds*Sds(ii,:)
  inp.Mds = setMat(inp.Mds, inp.k, Chainrule(inp.mdm, inp.Mds, inp.i, inp.mds, inp.Sds, ii));
  ret.Mds = inp.Mds;//coeffProdMat(1.0, inp.Mds);
  
  //Sdm(kk,:) = sdm*Mdm(i,:) + sds*Sdm(ii,:)
  inp.Sdm = setMat(inp.Sdm, kk, Chainrule(inp.sdm, inp.Mdm, inp.i, inp.sds, inp.Sdm, ii));
  ret.Sdm =  inp.Sdm;//coeffProdMat(1.0, inp.Sdm);
  
  //Sds(kk,:) = sdm*Mds(i,:) + sds*Sds(ii,:)
  inp.Sds = setMat(inp.Sds, kk, Chainrule(inp.sdm, inp.Mds, inp.i, inp.sds, inp.Sds, ii));
  ret.Sds = inp.Sds;//coeffProdMat(1.0, inp.Sds);
  
  //dCdm      = Cdm*Mdm(i,:) + Cds*Sdm(ii,:)
  double dCdm[][] = Chainrule(inp.Cdm, inp.Mdm, inp.i, inp.Cds, inp.Sdm, ii);
  
  //dCds      = Cdm*Mds(i,:) + Cds*Sds(ii,:)
  double dCds[][] = Chainrule(inp.Cdm, inp.Mds, inp.i, inp.Cds, inp.Sds, ii);
  
  //q = S(j,i)*C; S(j,k) = q; S(k,j) = q'
  double S[][] = addTo(inp.S, inp.i, inp.j, inp.k, inp.C);
  ret.S = S;//coeffProdMat(1.0, S);
  
  //SS = kron(eye(length(k)),S(j,i)); CC = kron(C',eye(length(j)))
  
  double SS[][] = kron(eye(inp.k.length), getMat(inp.S, inp.j, inp.i));
  double CC[][] = kron(transMat(inp.C), eye(inp.j.length));
  
  //Sdm(jk,:) = SS*dCdm + CC*Sdm(ji,:); Sdm(kj,:) = Sdm(jk,:)
  inp.Sdm = setMat(inp.Sdm, jk, Chainrule(SS, dCdm, CC, getMat(inp.Sdm, ji)));
  inp.Sdm = setMat(inp.Sdm, kj, getMat(inp.Sdm, jk));
  ret.Sdm = inp.Sdm;//coeffProdMat(1.0, inp.Sdm);
  
  //Sds(jk,:) = SS*dCds + CC*Sds(ji,:); Sds(kj,:) = Sds(jk,:)
  inp.Sds = setMat(inp.Sds, jk, Chainrule(SS, dCds, CC, getMat(inp.Sds, ji)));
  inp.Sds = setMat(inp.Sds, kj, getMat(inp.Sds, jk));
  ret.Sds = inp.Sds;//coeffProdMat(1.0, inp.Sds);
 
  return ret;
  
}

// Средняя версия
//[S Mdm Mds Sdm Sds Mdp Sdp] = fillIn(S,C,mdm,sdm,Cdm,mds,sds,Cds,Mdm,Sdm,Mds,Sds,Mdp,Sdp,[ ],i,j,k,D3)

//function [S Mdm Mds Sdm Sds Mdp Sdp] = fillIn(S,C,mdm,sdm,Cdm,mds,sds,Cds,Mdm,Sdm,Mds,Sds,Mdp,Sdp,[], i,j,k,D)

nargout_fillInMd_t fillInMd(nargin_fillInMd_t inp)
{
  
  nargout_fillInMd_t ret = new nargout_fillInMd_t();

  //X = reshape(1:D*D,[D D]); XT = X'                         vectorized indices
    
  int X[][] = new int [inp.D][inp.D];
  
  int count = 0;
  for(int i1 = 0; i1 < inp.D; i1++) 
  for(int i2 = 0; i2 < inp.D; i2++) 
    X[i2][i1] = count++;
  
  int XT[][] = transMat(X);
  
  //I=0*X; I(i,i)=1; ii=X(I==1)'; 
  int ii[] = new int [inp.i.length * inp.i.length];
  int ct = 0;
  for(int i1 = 0; i1 < inp.i.length; i1++) 
  for(int i2 = 0; i2 < inp.i.length; i2++) 
    ii[ct++] = X[inp.i[i2]][inp.i[i1]];
  
  //I=0*X; I(k,k)=1; kk=X(I==1)'
  int kk[] = new int [inp.k.length * inp.k.length];
  ct = 0;
  for(int i1 = 0; i1 < inp.k.length; i1++) 
  for(int i2 = 0; i2 < inp.k.length; i2++) 
    kk[ct++] = X[inp.k[i2]][inp.k[i1]];
    
  //I=0*X; I(j,i)=1; ji=X(I==1)' 
  int ji[] = new int [inp.j.length * inp.i.length];
  ct = 0;
  for(int i1 = 0; i1 < inp.i.length; i1++) 
  for(int i2 = 0; i2 < inp.j.length; i2++) 
    ji[ct++] = X[inp.j[i2]][inp.i[i1]];
    
  //I=0*X; I(j,k)=1; jk=X(I==1)'; kj=XT(I==1)'
  int jk[] = new int [inp.j.length * inp.k.length];
  int kj[] = new int [inp.k.length * inp.j.length];
  ct = 0;
  for(int i1 = 0; i1 < inp.k.length; i1++) 
  for(int i2 = 0; i2 < inp.j.length; i2++) 
    jk[ct++] = X[inp.j[i2]][inp.k[i1]];
    
  ct = 0;
  for(int i1 = 0; i1 < inp.k.length; i1++) 
  for(int i2 = 0; i2 < inp.j.length; i2++) 
    kj[ct++] = XT[inp.j[i2]][inp.k[i1]];
  
  //Mdm(k,:)  = mdm*Mdm(i,:) + mds*Sdm(ii,:)
  inp.Mdm = setMat(inp.Mdm, inp.k, Chainrule(inp.mdm, inp.Mdm, inp.i, inp.mds, inp.Sdm, ii));
  ret.Mdm = inp.Mdm;//coeffProdMat(1.0, inp.Mdm);
  
  //Mds(k,:)  = mdm*Mds(i,:) + mds*Sds(ii,:)
  inp.Mds = setMat(inp.Mds, inp.k, Chainrule(inp.mdm, inp.Mds, inp.i, inp.mds, inp.Sds, ii));
  ret.Mds = inp.Mds;//coeffProdMat(1.0, inp.Mds);
  
  //Sdm(kk,:) = sdm*Mdm(i,:) + sds*Sdm(ii,:)
  inp.Sdm = setMat(inp.Sdm, kk, Chainrule(inp.sdm, inp.Mdm, inp.i, inp.sds, inp.Sdm, ii));
  ret.Sdm = inp.Sdm;//coeffProdMat(1.0, inp.Sdm);
  
  //Sds(kk,:) = sdm*Mds(i,:) + sds*Sds(ii,:)
  inp.Sds = setMat(inp.Sds, kk, Chainrule(inp.sdm, inp.Mds, inp.i, inp.sds, inp.Sds, ii));
  ret.Sds = inp.Sds;//coeffProdMat(1.0, inp.Sds);
  
  //dCdm      = Cdm*Mdm(i,:) + Cds*Sdm(ii,:)
  double dCdm[][] = Chainrule(inp.Cdm, inp.Mdm, inp.i, inp.Cds, inp.Sdm, ii);
  
  //dCds      = Cdm*Mds(i,:) + Cds*Sds(ii,:)
  double dCds[][] = Chainrule(inp.Cdm, inp.Mds, inp.i, inp.Cds, inp.Sds, ii);
  
  //Mdp(k,:)  = mdm*Mdp(i,:) + mds*Sdp(ii,:)
  inp.Mdp = setMat(inp.Mdp, inp.k, Chainrule(inp.mdm, inp.Mdp, inp.i, inp.mds, inp.Sdp, ii));
  ret.Mdp = inp.Mdp;//coeffProdMat(1.0, inp.Mdp);
  
  //Sdp(kk,:) = sdm*Mdp(i,:) + sds*Sdp(ii,:)
  inp.Sdp = setMat(inp.Sdp, kk, Chainrule(inp.sdm, inp.Mdp, inp.i, inp.sds, inp.Sdp, ii));
  ret.Sdp = inp.Sdp;//coeffProdMat(1.0, inp.Sdp);
  
  //dCdp      = Cdm*Mdp(i,:) + Cds*Sdp(ii,:)
  double dCdp[][] = Chainrule(inp.Cdm, inp.Mdp, inp.i, inp.Cds, inp.Sdp, ii);
 
  //q = S(j,i)*C; S(j,k) = q; S(k,j) = q'
  inp.S = addTo(inp.S, inp.i, inp.j, inp.k, inp.C);
  ret.S = inp.S;//coeffProdMat(1.0, inp.S);
  
  //SS = kron(eye(length(k)),S(j,i)); CC = kron(C',eye(length(j)))
  
  double SS[][] = kron(eye(inp.k.length), getMat(inp.S, inp.j, inp.i));
  double CC[][] = kron(transMat(inp.C), eye(inp.j.length));
  
  //Sdm(jk,:) = SS*dCdm + CC*Sdm(ji,:); Sdm(kj,:) = Sdm(jk,:)
  inp.Sdm = setMat(inp.Sdm, jk, Chainrule(SS, dCdm, CC, getMat(inp.Sdm, ji)));
  inp.Sdm = setMat(inp.Sdm, kj, getMat(inp.Sdm, jk));
  ret.Sdm = inp.Sdm;//coeffProdMat(1.0, inp.Sdm);
  
  //Sds(jk,:) = SS*dCds + CC*Sds(ji,:); Sds(kj,:) = Sds(jk,:)
  inp.Sds = setMat(inp.Sds, jk, Chainrule(SS, dCds, CC, getMat(inp.Sds, ji)));
  inp.Sds = setMat(inp.Sds, kj, getMat(inp.Sds, jk));
  ret.Sds = inp.Sds;//coeffProdMat(1.0, inp.Sds);
  
  //Sdp(jk,:) = SS*dCdp + CC*Sdp(ji,:); Sdp(kj,:) = Sdp(jk,:)
  
  ret.Sdp = setMat(ret.Sdp, jk, Chainrule(SS, dCdp, CC, getMat(ret.Sdp, ji)));
  ret.Sdp = setMat(ret.Sdp, kj, getMat(ret.Sdp, jk));
 
  return ret;
  
}

void learnPolicy()
{
  println("learnPolicy start");
  
  double mM[] = {0, 0, PI, PI};
  double vM[][] = eye(4);
  vM = coeffProdMat(0.01, vM);
  vM[2][2] *= 0.01;
  vM[3][3] *= 0.01;
  
  println("mu0 для learnPolicy: ");
  printVec(mM);
  println("policy.p.inputs: " + policy.p.inputs.length + " x " + policy.p.inputs[0].length);
  println("policy.p.targets: " + policy.p.targets.length);
  println("policy.p.hyp: " + policy.p.hyp.length);
  println("dynmodel.inputs: " + dynmodel.inputs.length + " x " + dynmodel.inputs[0].length);
  println("dynmodel.dyni: " + dynmodel.dyni.length);
  println("H = " + H);
  
  try {
    nargout_minimize_pol_t rt = minimize_pol(mM, vM, dynmodel, policy, cost, H);
    println("minimize_pol - OK, J = " + rt.fX[rt.fX.length-1]);
    policy.p = let_p(rt.pol);
    println("learnPolicy - DONE");
  } catch(Exception e) {
    println("ОШИБКА в minimize_pol: " + e.toString());
    e.printStackTrace();
  }
  //double mM[] = {0, 0, 0, 0};
  //double vM[][] = eye(4);
  //vM = coeffProdMat(0.01, vM);
  //vM[2][2] *= 0.01;
  //vM[3][3] *= 0.01;
  
  //mM[0] = 0;
  //mM[1] = 0;
  //mM[2] = 3.1415;
  //mM[3] = 3.1415;
  //nargout_minimize_pol_t rt = minimize_pol(mM, vM, dynmodel, policy, cost, H);
 
  //println("Stop");
  
  ////**************                          ***********//
  //policy.p = let_p(rt.pol);
  ////policy.p.inputs = coeffProdMat(1.0, rt.pol.inputs);
  ////policy.p.targets = vectorCoef(rt.pol.targets, 1.0);
  ////policy.p.hyp = vectorCoef(rt.pol.hyp, 1.0);
      
}

void applyController()
{
  
   double mM[] = {0, 0, 0, 0};
  double vM[][] = eye(4);
  vM = coeffProdMat(0.01, vM);
  vM[2][2] *= 0.01;
  vM[3][3] *= 0.01;
  
  mM[0] = 0;
  mM[1] = 0;
  mM[2] = 3.1415;
  mM[3] = 3.1415;
  double     mu0[]           = {0, 0, 3.1415, 3.1415}, //  Угловая скорость / Угол
             S0[][]          = {
                                 {0.01, 0, 0, 0},
                                 {0, 0.01, 0, 0},
                                 {0, 0, 0.0001, 0},
                                 {0, 0, 0, 0.0001}
                                };
  get_rollout = rollout(gaussian(mu0, S0), policy, H, plant, cost);

}
