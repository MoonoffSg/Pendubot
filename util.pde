
  class func
  {
    
     int           mode = 0;
     
     //********************************************************************************//
     //                                    Политика                                    //
     int           _H   = 0;
     double        mm0[], 
                   sS0[][];
     
     dynmodel_t    _gpmd; 
     policy_t      _pol; 
     cost_t        _cst;
     nargout_f_t   ret;
     nargout_value_t vl = new nargout_value_t();
     
     //*********************************************************************************//
     //                                      ГП                                         //
     double         xE[][];
     double         yE[]; 
     curb_t         curb;
     //                                                                                 //
     //*********************************************************************************//
     
     nargout_f_t f(double x[])
     {
       ret = new nargout_f_t();
       
       // Функция стоимости для оптимизации политики
       if(mode == 1)
       {
         
         tp_p x_p = rewrap(_pol.p, x);
         vl = value(x_p, mm0, sS0, _gpmd, _pol, _cst, _H);
         
         ret.fx = vl.J;
         
         ret.dfx = unwrap(vl.dJdp);
         //int p_index = 0;
         
         //ret.dfx = new double[vl.dJdp.inputs.length * vl.dJdp.inputs[0].length + vl.dJdp.targets.length + vl.dJdp.hyp.length];
         //for(int i1 = 0; i1 < vl.dJdp.hyp.length ; i1++)
         //  ret.dfx[p_index++] = vl.dJdp.hyp[i1];
         
         //for(int i1 = 0; i1 < vl.dJdp.inputs[0].length ; i1++)
         //for(int j1 = 0; j1 < vl.dJdp.inputs.length ; j1++)
         //  ret.dfx[p_index++] = vl.dJdp.inputs[j1][i1];
         
         //for(int i1 = 0; i1 < vl.dJdp.targets.length ; i1++)
         //  ret.dfx[p_index++] = vl.dJdp.targets[i1];
          
          
       }
       // Функция правдоподобия
       else
       if(mode == 2)
       {
         ret = hypCurb2(x, xE, yE, curb);
         
       }
       return ret;
     }
     
     func(double m0[], double S0[][], dynmodel_t _gpmodel, policy_t _policy, cost_t _cost, int Hh)
     {
 
       mm0 = m0;//vectorCoef(m0, 1.0);
       sS0 = S0;//coeffProdMat(1.0, S0);
       _gpmd = _gpmodel; 
       _pol  = _policy; 
       _cst  = _cost;
       _H    = Hh;
     }
      
     func(double xt[][], double yt[], curb_t curbt) 
     {   
      
       xE = xt;//coeffProdMat(1.0, xt);
       yE = yt;//vectorCoef(yt, 1.0);
       curb = curbt;

     }
  }

gTrig_r  gTrig(double m[], double v[][], int i[], double e)
{
  
  gTrig_r ret = new gTrig_r();
  
  int     d            = m.length,
          I            = i.length;
          
  double  ee[]         = new double [2 * I],
          et[][]       = new double [2 * I][2 * I];
  
  ee = eyeV(ee.length);
  
  double  mi[]         = {m[i[0]], m[i[1]]},    // В состоянии i позиция - угол
          vi[][]       = {
                            {v[i[0]][i[0]], v[i[0]][i[1]]},
                            {v[i[1]][i[0]], v[i[1]][i[1]]}
                         },
          vii[]        = {v[i[0]][i[0]], v[i[1]][i[1]]},
          lq[][]       = { 
                          {0, 0},
                          {0, 0}
                         },
          q[][]           ,
          U1[][]          ,
          U2[][]          ,
          U3[][]          ,
          U4[][]          ;
  
  double  M[]          = new double [2 * I],
          V[][]        = new double [2 * I][2 * I],
          C[][]        = new double [d][2 * I];
                         
  int     Ic[]         = {1, 3},
          Is[]         = {0, 2};      
          
  M[Is[0]] = e * exp_db(-vii[0] * 0.5) * sin_db(mi[0]);
  M[Is[1]] = e * exp_db(-vii[1] * 0.5) * sin_db(mi[1]);
  M[Ic[0]] = e * exp_db(-vii[0] * 0.5) * cos_db(mi[0]);
  M[Ic[1]] = e * exp_db(-vii[1] * 0.5) * cos_db(mi[1]);
  
  lq    = coeffProdMat(-0.5, bsxfun_plus(vii, vii));
  q     = EXP_F(lq);
  
  //U1 = (exp(lq+vi)-q).*sin(bsxfun(@minus,mi,mi'))
  U1    = matrixDOT((matrixSUB(EXP_F(matrixADD(lq, vi)), q)), SIN_F(bsxfun_minus(mi, mi)));
  //U2 = (exp(lq-vi)-q).*sin(bsxfun(@plus,mi,mi'));
  U2    = matrixDOT((matrixSUB(EXP_F(matrixSUB(lq, vi)), q)), SIN_F(bsxfun_plus(mi, mi)));
  //U3 = (exp(lq+vi)-q).*cos(bsxfun(@minus,mi,mi'))
  U3    = matrixDOT((matrixSUB(EXP_F(matrixADD(lq, vi)), q)), COS_F(bsxfun_minus(mi, mi)));
  //U4 = (exp(lq-vi)-q).*cos(bsxfun(@plus,mi,mi'))
  U4    = matrixDOT((matrixSUB(EXP_F(matrixSUB(lq, vi)), q)), COS_F(bsxfun_plus(mi, mi)));
  
  //V(Is,Is) = U3 - U4; V(Ic,Ic) = U3 + U4; V(Is,Ic) = U1 + U2; 
  //V(Ic,Is) = V(Is,Ic)'; V = ee*ee'.*V/2;                               % variance
  // setMat(double RCv[][], int indR[], int indC[], double TRm[][])
  V = setMat(V, Is, Is, matrixSUB(U3, U4));
  V = setMat(V, Ic, Ic, matrixADD(U3, U4));
  V = setMat(V, Is, Ic, matrixADD(U1, U2));
  V = setMat(V, Ic, Is, matrixADD(U1, U2));
  
  et = Outer_Product(ee, ee);
  
  V = coeffProdMat(0.5, matrixDOT(Outer_Product(ee, ee), V));
  
  //C = zeros(d,2*I); C(i,Is) = diag(M(Ic)); C(i,Ic) = diag(-M(Is))
  C = setMat(C, i, Is, diag(getVec(M, Ic)));
  C = setMat(C, i, Ic, diag(vectorCoef(getVec(M, Is), -1)));
  
  ret.M = M;//vectorCoef(M, 1.0);
  ret.V = V;//coeffProdMat(1.0, V);
  ret.C = C;//coeffProdMat(1.0, C);
  
  return ret;
  
}

gTrig_Full_r  gTrigF(double m[], double v[][], int i[], double e)
{
  
  gTrig_Full_r ret = new gTrig_Full_r();
  
  int     d            = m.length,
          I            = i.length;
          
  double  ee[]         = new double [2 * I],
          et[][]       = new double [2 * I][2 * I];
  
  ee = eyeV(ee.length);
  
  double  mi[]         = {m[i[0]], m[i[1]]},    // В состоянии i позиция - угол
          vi[][]       = {
                            {v[i[0]][i[0]], v[i[0]][i[1]]},
                            {v[i[1]][i[0]], v[i[1]][i[1]]}
                         },
          vii[]        = {v[i[0]][i[0]], v[i[1]][i[1]]},
          lq[][]       = { 
                          {0, 0},
                          {0, 0}
                         },
          q[][]           ,
          U1[][]          ,
          U2[][]          ,
          U3[][]          ,
          U4[][]          ;
  
  double  M[]          = new double [2 * I],
          V[][]        = new double [2 * I][2 * I],
          C[][]        = new double [d][2 * I];
                         
  int     Ic[]         = {1, 3},
          Is[]         = {0, 2};      
          

  double  dVdm[][][]   = new double[2 * I][2 * I][d],
          dCdm[][][]   = new double[d][2 * I][d],
          dMdv[][][]   = new double[2 * I][d][d],
          dVdv[][][][] = new double[2 * I][2 * I][d][d], 
          dCdv[][][][] = new double[d][2 * I][d][d], 
          dMdm[][]     = new double[2 * I][d];
  
  M[Is[0]] = e * exp_db(-vii[0] * 0.5) * sin_db(mi[0]);
  M[Is[1]] = e * exp_db(-vii[1] * 0.5) * sin_db(mi[1]);
  M[Ic[0]] = e * exp_db(-vii[0] * 0.5) * cos_db(mi[0]);
  M[Ic[1]] = e * exp_db(-vii[1] * 0.5) * cos_db(mi[1]);
  
  lq    = coeffProdMat(-0.5, bsxfun_plus(vii, vii));
  q     = EXP_F(lq);
  
  //U1 = (exp(lq+vi)-q).*sin(bsxfun(@minus,mi,mi'))
  U1    = matrixDOT((matrixSUB(EXP_F(matrixADD(lq, vi)), q)), SIN_F(bsxfun_minus(mi, mi)));
  //U2 = (exp(lq-vi)-q).*sin(bsxfun(@plus,mi,mi'));
  U2    = matrixDOT((matrixSUB(EXP_F(matrixSUB(lq, vi)), q)), SIN_F(bsxfun_plus(mi, mi)));
  //U3 = (exp(lq+vi)-q).*cos(bsxfun(@minus,mi,mi'))
  U3    = matrixDOT((matrixSUB(EXP_F(matrixADD(lq, vi)), q)), COS_F(bsxfun_minus(mi, mi)));
  //U4 = (exp(lq-vi)-q).*cos(bsxfun(@plus,mi,mi'))
  U4    = matrixDOT((matrixSUB(EXP_F(matrixSUB(lq, vi)), q)), COS_F(bsxfun_plus(mi, mi)));
  
  //V(Is,Is) = U3 - U4; V(Ic,Ic) = U3 + U4; V(Is,Ic) = U1 + U2; 
  //V(Ic,Is) = V(Is,Ic)'; V = ee*ee'.*V/2;                               % variance
  // setMat(double RCv[][], int indR[], int indC[], double TRm[][])
  V = setMat(V, Is, Is, matrixSUB(U3, U4));
  V = setMat(V, Ic, Ic, matrixADD(U3, U4));
  V = setMat(V, Is, Ic, matrixADD(U1, U2));
  V = setMat(V, Ic, Is, transMat(matrixADD(U1, U2)));
  
  et = Outer_Product(ee, ee);
  
  V = coeffProdMat(0.5, matrixDOT(Outer_Product(ee, ee), V));
  
  //C = zeros(d,2*I); C(i,Is) = diag(M(Ic)); C(i,Ic) = diag(-M(Is))
  C = setMat(C, i, Is, diag(getVec(M, Ic)));
  C = setMat(C, i, Ic, diag(vectorCoef(getVec(M, Is), -1)));
  
  dMdm = transMat(C);
  
  for (int j = 0; j < I; j++)
  {
    //u = zeros(I,1); u(j) = 1/2
    double u[] = {0, 0};
    u[j] = 1.0 / 2.0;
    
    //dVdm(Is,Is,i(j)) = e*e'.*(-U1.*bsxfun(@minus,u,u')+U2.*bsxfun(@plus,u,u'))
    dVdm = set3DMass(dVdm, Is, Is, i[j], coeffProdMat(e * e, 
            matrixADD(matrixDOT(coeffProdMat(-1, U1), bsxfun_minus(u, u)),
            matrixDOT(U2, bsxfun_plus(u, u)))));
    
    //dVdm(Ic,Ic,i(j)) = e*e'.*(-U1.*bsxfun(@minus,u,u')-U2.*bsxfun(@plus,u,u'))
    dVdm = set3DMass(dVdm, Ic, Ic, i[j], coeffProdMat(e * e, 
            matrixSUB(matrixDOT(coeffProdMat(-1, U1), bsxfun_minus(u, u)),
            matrixDOT(U2, bsxfun_plus(u, u)))));
    
    //dVdm(Is,Ic,i(j)) = e*e'.*(U3.*bsxfun(@minus,u,u') +U4.*bsxfun(@plus,u,u'))
    double MatTm[][] = coeffProdMat(e * e, 
            matrixADD(matrixDOT(U3, bsxfun_minus(u, u)),
            matrixDOT(U4, bsxfun_plus(u, u))));
    dVdm = set3DMass(dVdm, Is, Ic, i[j], MatTm);
    
    //dVdm(Ic,Is,i(j)) = dVdm(Is,Ic,i(j))'
    dVdm = set3DMass(dVdm, Ic, Is, i[j], transMat(MatTm));
    
    //dVdv(Is(j),Is(j),i(j),i(j)) = exp(-vii(j)) * ...
    //                           (1+(2*exp(-vii(j))-1)*cos(2*mi(j)))*e(j)*e(j)/2
    dVdv[Is[j]][Is[j]][i[j]][i[j]] = exp_db(-vii[j]) * 
                        (1.0 + (2.0 * exp_db(-vii[j]) - 1.0) *
                        cos_db(2.0 * mi[j])) * e * e / 2.0;
    
    //dVdv(Ic(j),Ic(j),i(j),i(j)) = exp(-vii(j)) * ...
    //                           (1-(2*exp(-vii(j))-1)*cos(2*mi(j)))*e(j)*e(j)/2
    dVdv[Ic[j]][Ic[j]][i[j]][i[j]] = exp_db(-vii[j]) * 
                        (1.0 - (2.0 * exp_db(-vii[j]) - 1.0) *
                        cos_db(2.0 * mi[j])) * e * e / 2.0;
   
    //dVdv(Is(j),Ic(j),i(j),i(j)) = exp(-vii(j)) * ...
    //                               (1-2*exp(-vii(j)))*sin(2*mi(j))*e(j)*e(j)/2
    dVdv[Is[j]][Ic[j]][i[j]][i[j]] = exp_db(-vii[j]) * 
                        (1.0 - 2.0 * exp_db(-vii[j])) *
                        sin_db(2.0 * mi[j]) * e * e / 2.0;
    
    //dVdv(Ic(j),Is(j),i(j),i(j)) = dVdv(Is(j),Ic(j),i(j),i(j))
    dVdv[Ic[j]][Is[j]][i[j]][i[j]]  = dVdv[Is[j]][Ic[j]][i[j]][i[j]];
    
    //for k = [1:j-1 j+1:I]
    int k = 1 - j;
    
    //dVdv(Is(j),Is(k),i(j),i(k)) = (exp(lq(j,k)+vi(j,k)).*cos(mi(j)-mi(k)) ...
    //                     + exp(lq(j,k)-vi(j,k)).*cos(mi(j)+mi(k)))*e(j)*e(k)/2
    dVdv[Is[j]][Is[k]][i[j]][i[k]] = (exp_db(lq[j][k] + vi[j][k]) * cos_db(mi[j] - mi[k]) +
                        exp_db(lq[j][k] - vi[j][k]) * cos_db(mi[j] + mi[k])) * e * e / 2.0;
    
    
    //dVdv(Is(j),Is(k),i(j),i(j)) = -V(Is(j),Is(k))/2 
    dVdv[Is[j]][Is[k]][i[j]][i[j]] = -V[Is[j]][Is[k]] / 2.0;
    
    //dVdv(Is(j),Is(k),i(k),i(k)) = -V(Is(j),Is(k))/2
    dVdv[Is[j]][Is[k]][i[k]][i[k]] = -V[Is[j]][Is[k]] / 2.0;
    
    //dVdv(Ic(j),Ic(k),i(j),i(k)) = (exp(lq(j,k)+vi(j,k)).*cos(mi(j)-mi(k)) ...
    //                   - exp(lq(j,k)-vi(j,k)).*cos(mi(j)+mi(k)))*e(j)*e(k)/2
    dVdv[Ic[j]][Ic[k]][i[j]][i[k]] = (exp_db(lq[j][k] + vi[j][k]) * cos_db(mi[j] - mi[k]) -
                        exp_db(lq[j][k] - vi[j][k]) * cos_db(mi[j] + mi[k])) * e * e / 2.0;
     
    //dVdv(Ic(j),Ic(k),i(j),i(j)) = -V(Ic(j),Ic(k))/2
    dVdv[Ic[j]][Ic[k]][i[j]][i[j]] = -V[Ic[j]][Ic[k]] / 2.0;
    
    //dVdv(Ic(j),Ic(k),i(k),i(k)) = -V(Ic(j),Ic(k))/2
    dVdv[Ic[j]][Ic[k]][i[k]][i[k]] = -V[Ic[j]][Ic[k]] / 2.0;
   
    //dVdv(Ic(j),Is(k),i(j),i(k)) = -(exp(lq(j,k)+vi(j,k)).*sin(mi(j)-mi(k)) ...
    //                   + exp(lq(j,k)-vi(j,k)).*sin(mi(j)+mi(k)))*e(j)*e(k)/2
    dVdv[Ic[j]][Is[k]][i[j]][i[k]] = -(exp_db(lq[j][k] + vi[j][k]) * sin_db(mi[j] - mi[k]) +
                        exp_db(lq[j][k] - vi[j][k]) * sin_db(mi[j] + mi[k])) * e * e / 2.0;
   //  printMass4D(dVdv);                   
    //dVdv(Ic(j),Is(k),i(j),i(j)) = -V(Ic(j),Is(k))/2
    dVdv[Ic[j]][Is[k]][i[j]][i[j]] = -V[Ic[j]][Is[k]] / 2.0;
     
    //dVdv(Ic(j),Is(k),i(k),i(k)) = -V(Ic(j),Is(k))/2
    dVdv[Ic[j]][Is[k]][i[k]][i[k]] = -V[Ic[j]][Is[k]] / 2.0;   
     
    //dVdv(Is(j),Ic(k),i(j),i(k)) = (exp(lq(j,k)+vi(j,k)).*sin(mi(j)-mi(k)) ...
    //                   - exp(lq(j,k)-vi(j,k)).*sin(mi(j)+mi(k)))*e(j)*e(k)/2
    dVdv[Is[j]][Ic[k]][i[j]][i[k]] = (exp_db(lq[j][k] + vi[j][k]) * sin_db(mi[j] - mi[k]) -
                        exp_db(lq[j][k] - vi[j][k]) * sin_db(mi[j] + mi[k])) * e * e / 2.0;
     
    //dVdv(Is(j),Ic(k),i(j),i(j)) = -V(Is(j),Ic(k))/2
    dVdv[Is[j]][Ic[k]][i[j]][i[j]] = -V[Is[j]][Ic[k]] / 2.0;
    
    //dVdv(Is(j),Ic(k),i(k),i(k)) = -V(Is(j),Ic(k))/2
    dVdv[Is[j]][Ic[k]][i[k]][i[k]] = -V[Is[j]][Ic[k]] / 2.0;
    
    //dCdm(i(j),Is(j),i(j)) = -M(Is(j)); dCdm(i(j),Ic(j),i(j)) = -M(Ic(j));
    dCdm[i[j]][Is[j]][i[j]] = -M[Is[j]]; 
    dCdm[i[j]][Ic[j]][i[j]] = -M[Ic[j]];
    
    //dCdv(i(j),Is(j),i(j),i(j)) = -C(i(j),Is(j))/2;
    //dCdv(i(j),Ic(j),i(j),i(j)) = -C(i(j),Ic(j))/2;
    dCdv[i[j]][Is[j]][i[j]][i[j]] = -C[i[j]][Is[j]] / 2.0;
    dCdv[i[j]][Ic[j]][i[j]][i[j]] = -C[i[j]][Ic[j]] / 2.0;
  }
  
//  printMass3D(dVdm);
  
  for(int k = 0; k < d ; k++)
  for(int ii = 0; ii < d ; ii++)
  for(int j = 0; j < (2 * I) ; j++)
    dMdv[j][ii][k] = dCdm[ii][j][k] / 2.0;
  
//  printMass3D(dMdv); //<>// //<>// //<>// //<>//
  
  ret.M = M;//vectorCoef(M, 1.0);
  ret.V = V;//coeffProdMat(1.0, V);
  ret.C = C;//coeffProdMat(1.0, C);
  
  ret.dMdv = reshape(dMdv, 2 * I, d * d);
  ret.dVdm = reshape(dVdm, 4 * I * I, d);
  ret.dVdv = reshape(dVdv, 4 * I * I, d * d);
  ret.dCdv = reshape(dCdv, d*2*I, d * d);
  ret.dCdm = reshape(dCdm, d*2*I, d);
  ret.dMdm = dMdm;//coeffProdMat(1.0, dMdm);
  return ret;
  
}

gTrig_r  gTrig(double m[], double v[][], int i, double e)
{
  
  gTrig_r ret = new gTrig_r();
  
  int     d            = m.length,
          I            = 1;
          
  double  ee[]         = new double [2 * I],
          et[][]       = new double [2 * I][2 * I];
  
  ee = eyeV(ee.length);
  
  double  mi           = m[i],    // В состоянии i позиция - угол
          vi           = v[i][i],
          vii          = vi,
          lq           = 0,
          q            = 0,
          U1           = 0,
          U2           = 0,
          U3           = 0,
          U4           = 0;
          
  double  M[]          = new double [2 * I],
          V[][]        = new double [2 * I][2 * I],
          C[][]        = new double [d][2 * I];
                         
  int     Ic           = 2 * I - 1,
          Is           = Ic - 1;      
          
  M[Is] = e * exp_db(-vii * 0.5) * sin_db(mi);
  M[Ic] = e * exp_db(-vii * 0.5) * cos_db(mi);
  lq    = -(vii + vii) * 0.5f;
  q     = exp_db(lq);
  
  U1    = 0;
  U2    = (exp_db(lq - vi) - q) * sin_db(mi + mi);
  U3    = (exp_db(lq + vi) - q);
  U4    = (exp_db(lq - vi) - q) * cos_db(mi + mi);
  
  V[Is][Is] = U3 - U4;
  V[Ic][Ic] = U3 + U4;
  V[Is][Ic] = U1 + U2;
  V[Ic][Is] = V[Is][Ic];
  
  et = Outer_Product(ee, ee);
  
  V[0][0] *= et[0][0] * 0.5f;
  V[0][1] *= et[0][1] * 0.5f;
  V[1][0] *= et[1][0] * 0.5f;
  V[1][1] *= et[1][1] * 0.5f;
  
  C[i][Is] = M[Ic];
  C[i][Ic] = -M[Is];
  
  ret.M =M; //vectorCoef(M, 1.0);
  ret.V =V; //coeffProdMat(1.0, V);
  ret.C = C;//coeffProdMat(1.0, C);
  
  return ret;
  
}

gTrig_Full_r  gTrigF(double m[], double v[][], int i, double e)
{
  
  gTrig_Full_r ret = new gTrig_Full_r();
  
  int     d            = m.length,
          I            = 1;
          
  double  ee[]         = {1, 1},
          et[][]       = {
                          {0, 0},
                          {0, 0}
                         };
  
  double  mi           = m[i],    // В состоянии 1 позиция - угол
          vi           = v[i][i],
          vii          = vi,
          lq           = 0,
          q            = 0,
          U1           = 0,
          U2           = 0,
          U3           = 0,
          U4           = 0;
          
  double  M[]          = {0, 0};
  double  V[][]        = new double[d][d],
          C[][]        = new double[d][2 * I];
                         
  int     Ic           = 2 * (I) - 1,
          Is           = Ic - 1;      

  double  dVdm[][][]   = new double[2 * I][2 * I][d],
          dCdm[][][]   = new double[d][2 * I][d],
          dMdv[][][]   = new double[2 * I][d][d],
          dVdv[][][][] = new double[2 * I][2 * I][d][d], 
          dCdv[][][][] = new double[d][2 * I][d][d], 
          dMdm[][]     = new double[2 * I][d];
  
  M[Is] = e * exp_db(-vii * 0.5) * sin_db(mi);
  M[Ic] = e * exp_db(-vii * 0.5) * cos_db(mi);
  lq    = -(vii + vii) * 0.5f;
  q     = exp_db(lq);
  
  U1    = 0;
  U2    = (exp_db(lq - vi) - q) * sin_db(mi + mi);
  U3    = (exp_db(lq + vi) - q);
  U4    = (exp_db(lq - vi) - q) * cos_db(mi + mi);
  
  V[Is][Is] = U3 - U4;
  V[Ic][Ic] = U3 + U4;
  V[Is][Ic] = U1 + U2;
  V[Ic][Is] = V[Is][Ic];
  
  ee[0] = e;
  ee[1] = e;
  et = Outer_Product(ee, ee);
  
  V[0][0] *= et[0][0] * 0.5f;
  V[0][1] *= et[0][1] * 0.5f;
  V[1][0] *= et[1][0] * 0.5f;
  V[1][1] *= et[1][1] * 0.5f;
  
  C[i][Is] = M[Ic];
  C[i][Ic] = -M[Is];
  
  dMdm = transMat(C);
  
  for (int j = 0; j < I; j++)
  {
    double u = 1.0 / 2.0;
    dVdm[Is][Is][i] = e * e * (-U1 * 0 + U2 * (u + u));
    dVdm[Ic][Ic][i] = e * e *(-U1 * 0 - U2 * (u + u));
    dVdm[Is][Ic][i] = e * e * (U3 * 0 + U4 * (u + u));
    dVdm[Ic][Is][i] = dVdm[Is][Ic][i]; 
    dVdv[Is][Is][i][i] = exp_db(-vii) * (1.0 + (2.0 * exp_db(-vii) - 1.0) * 
                          cos_db(2.0 * mi)) * e * e / 2.0;
    dVdv[Ic][Ic][i][i] = exp_db(-vii) * (1.0 - (2.0 * exp_db(-vii) - 1.0) * 
                          cos_db(2.0 * mi)) * e * e / 2.0;
    dVdv[Is][Ic][i][i] = exp_db(-vii) * (1.0 - 2.0 * exp_db(-vii)) * 
                          sin_db(2.0 * mi) * e * e / 2.0;
    dVdv[Ic][Is][i][i] = dVdv[Is][Ic][i][i];
    
    dCdm[i][Is][i] = -M[Is]; 
    dCdm[i][Ic][i] = -M[Ic];
    
    dCdv[i][Is][i][i] = -C[i][Is] / 2.0;
    dCdv[i][Ic][i][i] = -C[i][Ic] / 2.0;
  }
  
//  printMass3D(dVdm);
  
  for(int k = 0; k < d ; k++)
  for(int ii = 0; ii < d ; ii++)
  for(int j = 0; j < (2 * I) ; j++)
    dMdv[j][ii][k] = dCdm[ii][j][k] / 2.0;
  
//  printMass3D(dMdv); //<>// //<>// //<>// //<>// //<>// //<>//
  
  ret.M = M;//vectorCoef(M, 1.0);
  ret.V = V;//coeffProdMat(1.0, V);
  ret.C = C;//coeffProdMat(1.0, C);
  
  ret.dMdv = reshape(dMdv, 2 * I, d * d);
  ret.dVdm = reshape(dVdm, 4 * I * I, d);
  ret.dVdv = reshape(dVdv, 4 * I * I, d * d);
  ret.dCdv = reshape(dCdv, d*2*I, d * d);
  ret.dCdm = reshape(dCdm, d*2*I, d);
  ret.dMdm = dMdm;//coeffProdMat(1.0, dMdm);
  return ret;
  
}

//===================================================================================
// Генерирует n выборок гауссовой случайной величины p(x) = N(m,S).
// Выборка основана на факторизации Холецкого ковариационной матрицы S.

// Входные аргументы
// M      математическое ожидание гауссовой случайной величины  [D x 1]
// S      ковариационной матрицы                                [D x D]
// n      количество выборок

// Выходные аргументы
// x      матрица с выборками гауссовой случайной величины      [D x n]

double[][] gaussian(double M[], double S[][], int n)
{
  
  int        ml    = M.length;
  double[][] randn = new double[ml][n];
  double[][] gss   = new double[ml][n];
  double[][] _S    = new double[ml][ml];                  
  int        i     = 0,
             j     = 0;
  
  Chol_S(S, ml);
  _S = getSch();
  _S = transMat(_S);
  
  for(i = 0; i < ml; i ++)
    for(j = 0; j < n ; j++)
      randn[i][j]  = (double)randomGaussian() ;
     
      
  gss = matrixMultiply(_S, randn);
  
  for(i = 0; i < ml; i ++)
    for(j = 0; j < n ; j++)
        gss[i][j] += M[i]; 

  return gss;
  
}

double[] gaussian(double M[], double S[][])
{
  
  int        ml    = M.length;
  double[]   randn = new double[ml];
  double[]   gss   = new double[ml];
  double[][] _S    = new double[ml][ml];                  
  int        i     = 0;
  
  Chol_S(S, ml);
  _S = getSch();
  _S = transMat(_S);
  
  for(i = 0; i < ml; i ++)
     randn[i]  = (double)randomGaussian() ;
      
  gss = matrixMultiplyC(_S, randn);
  
  for(i = 0; i < ml; i ++)
     gss[i] += M[i]; 

  return gss;
  
}

//========================================================================================================================
// sq_dist - функция для вычисления матрицы, которая содержит попарных квадратов расстояний между двумя наборами векторов.
// Вектора хранятся в столбцах двух матриц, матрица a (размера D на n) и матрица b (размера D на m). 
// Если указан только один аргумент вторая матрица считается идентичной первой.

double [][] sq_dist(int D, double a[][], int n, double b[][], int m) // a Dxn, b Dxm
{
  
   double[][]  C        = new double[n][m];
   int         i        = 0,
               j        = 0,
               k        = 0;
   double      z        = 0,
               t        = 0;
   
   for(i = 0; i < n; i++) 
     for(j = 0; j < m; j++) 
     {
        z = 0.0;
        for(k = 0; k < D; k++) 
        { 
          t = a[k][i] - b[k][j]; 
          z += t*t; 
        }
        
        C[i][j] = z;
        
      }

  return C;
}

double [][] sq_dist(double a[], int n)
{
  
   double[][]  C        = new double[n][n];
   int         i        = 0,
               j        = 0;
   double      z        = 0,
               t        = 0;
   
   for(i = 0; i < n; i++) 
   {
     for(j = 0; j < n; j++) 
     {
       
          t = a[i] - a[j]; 
          z = t*t; 
       
        
        C[i][j] = z;
        
      }
   }

  return C;
}

//================================================================================================
// STD Среднеквадратическое отклонение. 
// Для векторов Y = STD(X) возвращает среднеквадратическое отклонение. 
// Для матриц Y — вектор-строка, содержащий среднеквадратическое отклонение каждого столбца. 

double [] std(double a[][])
{
  
   int         n        = a.length,
               D        = a[0].length;
   double[]    mean     = new double[D + 1];
   double[]    _var     = new double[D + 1];
   int         i        = 0,
               j        = 0;
   
   
   for(i = 0; i < D; i++) 
   {
     mean[i] = 0;
     for(j = 0; j < n; j++) 
     {
       
       mean[i] += a[j][i]; 
          
     }
     
     mean[i] /= (double)n;
     
   }
  
   for(i = 0; i < D; i++) 
   {
     _var[i] = 0;
     for(j = 0; j < n; j++) 
     {
       
          _var[i] += (a[j][i] - mean[i]) * (a[j][i] - mean[i]); 
          
     }
     
     _var[i] /= (double)(n - 1);
     _var[i] = sqrt_db(_var[i]);
     
   }
   
  return _var;
}
  
 
 //=================================================================================  
 // Функции необходимые для проведения оптимизации. Вспомогательные
 
  p_wp_t p_wp = new p_wp_t();
  
  int wp(p_linesearch_r pp, double SIG, int RHO) // Проверка условий Wolfe-Powell. Инициализация
  {
    
    p_wp.a = RHO * pp.s; 
    p_wp.b = pp.f; 
    p_wp.c = -SIG * pp.s; 
    p_wp.sig = SIG; 
    p_wp.rho = RHO; 
    return 0;
    
  }
  
  int wp(p_linesearch_r pp) // Проверка условий Wolfe-Powell. Основное тело
  {
    
     if (pp.f > (p_wp.a * pp.x + p_wp.b))  
     {
        if(p_wp.a > 0)
            return(-1);
         else 
            return(-2); 
     }
     else
     {
      if (pp.s < -p_wp.c)
            return 0; 
      else if (pp.s > p_wp.c)
                return(1);
           else 
                return (2);
    
      }
  }
  
  double minCubic(double x, double df, double s0, double s1, int extr) // минимизатор кубического аппроксимации
  {
   
    double    INT       = 0.1, // пределы интерполяции
              EXT       = 5.0, // пределы экстраполяции
              A         = -6.0 * df + 3.0 * (s0 + s1) * x,
              B         = 3.0 * df - (2.0 * s0 + s1) * x, 
              discr     = B * B - A * s0 * x,
              z         = 0;
    boolean   notReal   = false;
    
          if(discr < 0)
            notReal    = true;
            
          if (extr == 1) // мы экстраполируем?
          {
            if(notReal) // исправляем плохой z
            {
              z = EXT * x;
            }
            else
            {
              z = -s0 * x * x / (B + sqrt_db(discr));
              if((abs_db(z) > 1e10) || (z < x) || (z > (x * EXT)))
                z = EXT * x;
            }
            z = max_db(z, (1 + INT) * x);      
          }  // в противном случае мы интерполируем
          else
          {
            if(notReal) // исправляем плохой z
            {
              z = x * 0.5;
            }
            else
            {
              z = -s0 * x * x / (B + sqrt_db(discr));
              if((abs_db(z) > 1e10) || (z < 0) || (z > x))
                z = x * 0.5; 
            }
            z = min_db(max_db(z, INT * x), (1.0 - INT) * x); // по крайней мере INT от границ
          }
          
    return z;
    
  }
  
  p_linesearch_r Let(p_linesearch_r pT) // Для копирования структыры типа p_linesearch_r в функции линейного поиска
  {
   
    p_linesearch_r ret = new p_linesearch_r(pT.x, pT.f, pT.df, pT.s);
    
    return ret;
    
  }
    
  nargout_BFGS2_t BFGS2(func getF, double x0[], double fx0, double dfx0[], pt_struct_t p) // Реализация метода BFGS
  {
   
    int                   Dt            = x0.length;
    nargout_BFGS2_t       ret           = new nargout_BFGS2_t();
    nargout_linesearch2_r lineSearchRet = new nargout_linesearch2_r();
    
    int            ok,
                   i,
                   i_t;
    double         x[],
                   dfx[]    = new double[Dt],
                   r[]      = new double[Dt],
                   H[][]    = new double[Dt][Dt],
                   t[]      = new double[Dt],
                   y[]      = new double[Dt],
                   Hy[]     = new double[Dt],
                   ty       = 0,
                   s        = 0,
                   b        = 0,
                   normV;
    double         coffT;
    
    p.H = new double[Dt][Dt];
    ret.fx = new double[1];
    
    if(p.lgth < 0) 
                  i = 0;
      else            
                  i = 1;
    ok = 1;   
    
    x = x0;//vectorCoef(x0, 1.0);
    ret.fx[0] = fx0;
    r = vectorCoef(dfx0, -1.0);
    s = matrixMultiplyC_V(vectorCoef(r, -1), r);
    b = -1.0 / (s - 1); 
    
    H =  eye(Dt);
    
    while (i < abs(p.lgth))
    {
      
      // Вычисление нормы вектора SUM(ABS(V).^p)^(1/p) p = 2
      normV = 0;
      for(i_t = 0; i_t < Dt; i_t++)
        normV += abs_db(r[i_t]) * abs_db(r[i_t]);
        normV = pow_db(normV, 0.5);
    
      b = min_db(b, 1.0 / normV);                    
      b = max_db(b, 1e-7 / normV);
      
      lineSearchRet = LineSearch2(getF, x0, fx0, dfx0, r, s, b, i, p);
      x = lineSearchRet.x;//vectorCoef(lineSearchRet.x, 1.0);
      dfx =lineSearchRet.df;// vectorCoef(lineSearchRet.df, 1.0);
      b = lineSearchRet.a;
      i = lineSearchRet.i;
      fx0 = lineSearchRet.fx;
  
      if (i < 0)
      {
        i = -i; 
        if (ok != 0)
          ok = 0; 
        else 
          break;
      }
      else
      {
        ok = 1; 
        
        t = vectorSUM(1, x, -1, x0);
        y = vectorSUM(1, dfx, -1, dfx0);
        ty =  matrixMultiplyC_V(t, y);               
        Hy = matrixMultiplyC(H, y);
        
        // (ty+y'*Hy)/ty^2*t*t'
        coffT = (ty + matrixMultiplyC_V(y, Hy)) / (ty * ty);
        double t_[] = new double [Dt];
        t_ = vectorCoef(t, coffT);
        double partH[][] = new double[Dt][Dt];
        partH = Outer_Product(t, t_);
        H = matrixADD(H, partH);
        
        // -1/ty*Hy*t'
        t_ = vectorCoef(Hy, -(1.0 / ty));
        partH = Outer_Product(t_, t);
        H = matrixADD(H, partH);
        
        // - 1/ty*t*Hy'
        t_ = vectorCoef(t, -(1.0 / ty));
        partH = Outer_Product(t_, Hy);
        
        // H = H + (ty+y'*Hy)/ty^2*t*t' - 1/ty*Hy*t' - 1/ty*t*Hy'
        H = matrixADD(H, partH);
    
      }
      
      r = matrixMultiplyC(H, dfx);
      r = vectorCoef(r, -1); 
      s = matrixMultiplyC_V(r, dfx);
      x0 = x;//vectorCoef(x, 1.0);
      dfx0 = dfx;//vectorCoef(dfx, 1.0);
      ret.fx = append_db(ret.fx, fx0);
      p.H = H;//coeffProdMat(1.0, H);
    }
    ret.x = x;//vectorCoef(x, 1.0);
    ret.i = i;
    return ret;
    
  }
  
  nargout_linesearch2_r LineSearch2(func getF, double x0[], 
                          double f0, 
                          double df0[], 
                          double d[], 
                          double s, 
                          double a, 
                          int i, 
                          pt_struct_t pp) // Функция реализует линейный поиск минимума
  {
    
    nargout_linesearch2_r ret     = new nargout_linesearch2_r();
    int                   LIMIT   = 0,
                          j       = 0,
                          ok      = 0;
    p_linesearch_r        p0,
                          p1,
                          p2,
                          p3;
    
    nargout_f_t          ret_f;
    double               tmpVector[]    = new double[x0.length];
    
    if (pp.lgth < 0)
         LIMIT = min(pp.MFEPLS, -i - pp.lgth); 
     else 
         LIMIT   = pp.MFEPLS;

    p0 = new p_linesearch_r(0, f0, df0, s);
    p1 = new p_linesearch_r(0, f0, df0, s); // p1 = p0
    p2 = new p_linesearch_r(0, f0, df0, s); // p1 = p0
    p3 = new p_linesearch_r(a, 0, df0, 0);
       
    wp(p0, pp.SIG, 0);
    
    while (true)   
    {
      ok = 0; 
      while((ok == 0) && (j < LIMIT))
      {
        try
        {
          j = j+1; 
          tmpVector = vectorSUM(1, x0, p3.x, d);     
          ret_f = getF.f(tmpVector);  
          p3.f = ret_f.fx;
          p3.df = ret_f.dfx;//vectorCoef(ret_f.dfx, 1.0);
          p3.s = matrixMultiplyC_V(p3.df, d); 
          ok = 1; 
        }
        catch(Throwable e)
        {
          
          p3.x = (p1.x + p3.x) / 2.0; 
          ok = 0; 
          p3.fs = 'N';  
          p3.ss = 'N';
          
        }
       
      }
      
      if ((wp(p3) != 0) || (j >= LIMIT)) break;
      
      p0 = Let(p1);
      p1 = Let(p3);
      p3.x = p0.x + minCubic(p1.x - p0.x, p1.f - p0.f, p0.s, p1.s, 1);       
      
    }
    while (true)   
    {
      
      try {
       
        double isnan = p3.f + p3.s;
        
      }
      catch (Throwable e) 
      {
        p2 = Let(p1);
        break;
      }
      if(abs_db(p3.f + p3.s) > 1e10)
      {
        p2 = Let(p1);
        break;
      }
      
      if (p1.f > p3.f)
      {
        p2 = Let(p3);
      }
      else 
      {
        p2 = Let(p1);
      }
      
      if ((wp(p2) > 1) || (j >= LIMIT)) break; 
      p2.x = p1.x + minCubic(p3.x - p1.x, p3.f - p1.f, p1.s, p3.s, 0);    
      ok = 0; 
      
      while((ok == 0) && (j < LIMIT))
      {
        try
        {
          j = j+1; 
          tmpVector = vectorSUM(1, x0, p2.x, d);       
          ret_f = getF.f(tmpVector);  
          p2.f = ret_f.fx;
          p2.df = ret_f.dfx;//vectorCoef(ret_f.dfx, 1.0);
          p2.s = matrixMultiplyC_V(p2.df, d); 
          ok = 1; 
        }
        catch(Throwable e)
        {
          
          p2.x = (p1.x + p2.x) / 2.0f; 
          ok = 0; 
          p2.fs = 'N';  
          p2.ss = 'N';
          if (LIMIT == j)
          {
            
            p2 = Let(p1);
                    
          }
          
        }
       
      }
      
      if( (wp(p2) > -1) && (p2.s > 0) || (wp(p2) < -1))
      {
        p3 = Let(p2);
      }
      else 
      {
        p1 = Let(p2);
      }
      
    }
    
    //x = x0 +p2.x*d; 
    ret.x = vectorSUM(1, x0, p2.x, d); 
    ret.fx = p2.f; 
    ret.df = p2.df;//vectorCoef(p2.df, 1.0);
    ret.a = p2.x;        
    if (pp.lgth < 0)
      i = i + j;
    else 
      i = i + 1; 
      
      println(pp.S + " " + (i - 1)  + "; Значение функции " + ret.fx);
 
    if ((wp(p2) < 2)) 
    {
      i = -i;
    }
    ret.i = i;
    return ret;
    
  }
  
  
  //==================================================================================
  // Расстояние Махаланобиса по точкам, возведенное в квадрат  (a-b)*Q*(a-b)'
  // Входные аргументы:

  // a матрица, содержащая n векторов-строк [n x D]
  // b матрица, содержащая n векторов-строк [n x D]
  // Q весовая матрица.                     [D x D]

  // Выходные аргументы:
  // K квадратов расстояний по точкам       [n x n]
  
double [][] maha(double a[][], double b[][]) // Матрица Q - единичная                   
{
  int n = a.length,
      D = a[0].length;
      
  double R[][]  = new double [n][n],
         R1[][] = new double [n][n];
  double sm1[]  = new double [n],
         sm2[]  = new double [n];
   
  for(int i = 0; i < n ; i++)
  {
    sm1[i] = 0;
    for(int j = 0; j < D ; j++)
      sm1[i] += a[i][j] * a[i][j];
  }
  
  for(int i = 0; i < n ; i++)
  {
    sm2[i] = 0;
    for(int j = 0; j < D ; j++)
      sm2[i] += b[i][j] * b[i][j];
  }
  
  for(int i = 0; i < n ; i++)
    for(int j = 0; j < n ; j++)
    R[i][j] = sm1[i] + sm2[j];
  
  R1 = matrixMultiply(a, transMat(b));
  for(int i = 0; i < n ; i++)
    for(int j = 0; j < n ; j++)
    R1[i][j] *= -2.0;
  
  R = matrixADD(R, R1);
// K = bsxfun(@plus,sum(a.*a,2),sum(b.*b,2)')-2*a*b';

  return R;

}

double [][] maha(double a[][], double b[][], double Q[][]) // Общий случай                      
{
  int n = a.length,
      D = a[0].length;
      
  double R[][]  = new double [n][n],
         R1[][] = new double [n][n];
  double sm1[]  = new double [n],
         sm2[]  = new double [n];
  double aQ[][] = new double [n][D],
         bQ[][] = new double [n][D];
  
  // aQ = a*Q; K = bsxfun(@plus,sum(aQ.*a,2),sum(b*Q.*b,2)')-2*aQ*b';
  aQ = matrixMultiply(a, Q);
  bQ = matrixMultiply(b, Q);
   
  for(int i = 0; i < n ; i++)
  {
    sm1[i] = 0;
    for(int j = 0; j < D ; j++)
      sm1[i] += aQ[i][j] * a[i][j];
  }
  
  for(int i = 0; i < n ; i++)
  {
    sm2[i] = 0;
    for(int j = 0; j < D ; j++)
      sm2[i] += bQ[i][j] * b[i][j];
  }
  
  for(int i = 0; i < n ; i++)
    for(int j = 0; j < n ; j++)
    R[i][j] = sm1[i] + sm2[j];
  
  R1 = matrixMultiply(aQ, transMat(b));
  for(int i = 0; i < n ; i++)
    for(int j = 0; j < n ; j++)
    R1[i][j] *= -2.0;
  
  R = matrixADD(R, R1);

  return R;

}
