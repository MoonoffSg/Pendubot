
  nargout_minimize_gp_t minimize_gp( double       x[], 
                                     double       gp_inputs[][],
                                     double       gp_targets[],
                                     curb_t       cr) 
  {
    
    nargout_minimize_gp_t ret = new nargout_minimize_gp_t();
    
    //*********                                  ******//
        pt_struct_t opt   = new pt_struct_t(); 
        opt.S             = "               #";
        opt.S2            = "function evaluation #";
        opt.lgth          = 300;
        opt.MFEPLS        = 10;
        opt.MSR           = 100;
        opt.SIG           = 0.5;
    //********                                   ******//
    
    //********                                   ******//
        nargout_f_t  fx_dfx;
        func get = new func(gp_inputs, gp_targets, cr);
        get.mode = 2;
        fx_dfx = get.f(x);
    //********                                   ******//
  
    println("Начальное значение функции " + fx_dfx.fx);
    
    if(opt.lgth > 0) 
            opt.S = " #"; 
    else
            opt.S = "function evaluation #";
    
    nargout_BFGS2_t getBFGS = new nargout_BFGS2_t();
    getBFGS = BFGS2(get, x, fx_dfx.fx, fx_dfx.dfx, opt);
    ret.hyp = getBFGS.x;//vectorCoef(getBFGS.x, 1.0);
    ret.v =getBFGS.fx;// vectorCoef(getBFGS.fx, 1.0);
    return ret;
    
  }
  
