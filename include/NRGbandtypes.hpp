
#ifndef _HYBFUNC_LORENTZIAN_
#define _HYBFUNC_LORENTZIAN_

double HybFunc_Lorentzian(double omega,
			  void *params);

double HybFunc_Lorentzian_timesEn(double omega,
				  void *params);

#endif

////////

#ifndef _HYBFUNC_CONST_
#define _HYBFUNC_CONST_

double HybFunc_Const(double omega,
		     void *params);

double HybFunc_Const_timesEn(double omega,
			     void *params);

//Campo-Oliveira
double HybFunc_Const_divEn(double omega,
			   void *params);


#endif

////////

#ifndef _HYBFUNC_POWERLAW_
#define _HYBFUNC_POWERLAW_

double HybFunc_PowerLaw(double omega,
			void *params);


double HybFunc_PowerLaw_timesEn(double omega,
				void *params);

//Campo-Oliveira
double HybFunc_PowerLaw_divEn(double omega,
			      void *params);


#endif

////////

#ifndef _HYBFUNC_FROMFILE_
#define _HYBFUNC_FROMFILE_

double HybFunc_FromFile(double omega,
			    void *params);


double HybFunc_FromFile_timesEn(double omega,
				void *params);

//Campo-Oliveira
double HybFunc_FromFile_divEn(double omega,
			      void *params);


#endif

////////

