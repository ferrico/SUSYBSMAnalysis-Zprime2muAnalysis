 /***************************************************************************** 
  * Project: RooFit                                                           * 
  *                                                                           * 
  * Copyright (c) 2000-2005, Regents of the University of California          * 
  *                          and Stanford University. All rights reserved.    * 
  *                                                                           * 
  * Redistribution and use in source and binary forms,                        * 
  * with or without modification, are permitted according to the terms        * 
  * listed in LICENSE (http://roofit.sourceforge.net/license.txt)             * 
  *****************************************************************************/ 

 // -- CLASS DESCRIPTION [PDF] -- 
 // Your description goes here... 

 #include "Riostream.h" 

 #include "RooDrellYanC0.h" 
 #include "RooAbsReal.h" 
 #include "RooAbsCategory.h" 

 ClassImp(RooDrellYanC0) 

 RooDrellYanC0::RooDrellYanC0(const char *name, const char *title, 
                        RooAbsReal& _x,
                        RooAbsReal& _M,
                        RooAbsReal& _Gamma,
                        RooAbsReal& _theta) :
   RooAbsPdf(name,title), 
   x("x","x",this,_x),
   M("M","M",this,_M),
   Gamma("Gamma","Gamma",this,_Gamma),
   theta("theta","theta",this,_theta)
 { 
 } 


 RooDrellYanC0::RooDrellYanC0(const RooDrellYanC0& other, const char* name) :  
   RooAbsPdf(other,name), 
   x("x",this,other.x),
   M("M",this,other.M),
   Gamma("Gamma",this,other.Gamma),
   theta("theta",this,other.theta)
 { 
 } 



 Double_t RooDrellYanC0::evaluate() const 
 {
   Double_t num = M*M*Gamma*Gamma;
   Double_t den = (x*x - M*M)*(x*x - M*M) + x*x*x*x*Gamma*Gamma/M/M;
   
   return (num/den) * exp(-theta*x);
 } 
