  
// //=================================================================================  
// // Функции основные для проведения оптимизации
 
  nargout_minimize_pol_t minimize_pol(double         _m0[],
                                      double         _S0[][],
                                      dynmodel_t     _gpmodel,
                                      policy_t       _policy,
                                      cost_t         _cost,
                                      int            _H) // Минимизатор политики
  {
    
    nargout_minimize_pol_t ret = new nargout_minimize_pol_t();
    
    //*********                                  ******//
        pt_struct_t opt   = new pt_struct_t(); 
        opt.S             = "               #";
        opt.S2            = "function evaluation #";
        opt.lgth          = 150;
        opt.MFEPLS        = 20;
        opt.MSR           = 100;
        opt.SIG           = 0.5;
        opt.fh            = 1;
    //********                                   ******//
    
    //********                                   ******//
        nargout_f_t fx_dfx;
        func get = new func(_m0, _S0, _gpmodel, _policy, _cost, _H);
        get.mode = 1;
        
        double X[] = unwrap(_policy.p);
      
        fx_dfx = get.f(X);
    //********                                   ******//
    
    println("Начальное значение функции " + fx_dfx.fx);
    
    if(opt.lgth > 0) 
            opt.S = " #"; 
    else
            opt.S = "function evaluation #";
    
    nargout_BFGS2_t getBFGS = new nargout_BFGS2_t();
    getBFGS = BFGS2(get, X, fx_dfx.fx, fx_dfx.dfx, opt);
    
    ret.pol = rewrap(_policy.p, getBFGS.x);
    ret.fX = getBFGS.fx;//vectorCoef(getBFGS.fx, 1.0);
    
    return ret;
    
  }
  
  //<>// //<>// //<>//
