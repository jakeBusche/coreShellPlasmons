#!/usr/bin/env wolframscript
(* ::Package:: *)

(****************************************************************************)
(*************************** Initialization *********************************)
(****************************************************************************)

BeginPackage["fitFunc`"]

(****************************************************************************)
(*********************** Defining Physical Constants ************************)
(****************************************************************************)

c = 2.998*10^10. (* speed of light (cm/s) *)
cmks = 2.998*10^8 (* speed of light (m/s) *)
eps0 = 8.854187817*10^(-12.) (* permittivity of free space (C^2/(m^2*N)) *)
emks = 1.6021766208*10^(-19.) (* electron charge (C) *)
e = 4.8032042510 * 10^(-10) (* electron charge (statC) *)
nm = 10^(-7.) (* nm / cm *)
hbar = 6.626 * 10^(-27)  (* erg * s *)
hbareV = 6.58211951 * 10^(-16) (* reduced Planck constant (eV*s) *)
me = 9.10938*10^(-28) (* mass of electron (g) *)
meV = 510998.946 (* mass of electron (eV/c^2) *)

(****************************************************************************)
(*************************** Defining Functions *****************************)
(****************************************************************************)
(*

epsDrude : Function that returns a Drude-model dielectric constant at a given
	energy given the plasma energy, damping energy, and high-frequency 
	dielectric constant as inputs

phiBCeq : Function that returns the boundary conditions for the (l,m)^th mode
	of the electrostatic potential for a set of n concentric spherical 
	surfaces. The boundary conditions are returned as a list of symbolic 
	equations.

dBCeq : Function that returns the boundary conditions for the (l,m)^th mode
	of the displacement field for a set of n concentric spherical 
	surfaces. The boundary conditions are returned as a list of symbolic 
	equations.

BCsolns : Function that solves the boundary value problem for the concentric 
	spherical surfaces problem given the symbolic equations from phiBCeq and 
	dBCeq as inputs. The solutions are returned as a list of replacement 
	rules.

alphalm : Function that provides the effective polarizability for the (l,m)^th
	mode of a set of n concentric spherical surfaces. The polarizability is
	returned as a symbolic value. 

resFinder : Function that finds the resonant terms in the polarizability of 
	the (l,m)^th mode of a set of n concentric spherical surfaces and
	separates the 1) prefactors, 2) resonant energies, 3) damping energies,
	and 4) natural energies into four rows of a list. These quantities are
	returned as a matrix in list form.

massFinder : Function that finds the effective masses and oscillator strengths
	of the (l,m)^th mode of a set of n concentric spherical surfaces and 
	places them as the first and second rows of a matrix. These quantities are
	returned as a matrix in list form.

parFinder : Organizes the results of resFinder and massFinder to give a matrix
	with rows containing the 1) Resonance energies, 2) Damping energies, 
	3) natural energies, 4) Complex oscillator strengths, and 5) effective 
	masses of the resonances contained in the (l,m)^th mode of a system of n
	concentric spherical surfaces. These quantities are returned as a matrix 
	in list form.

sigAbsSymb : Calculates the symbolic absorption cross section in cm^2 for the
	(l,m)^th mode of a system of n concentric spherical surfaces. The 
	absorption cross section is returned as a symbolic value that contains
	only numbers and the symbol "En" when given numerical inputs for all
	parameters except the polarizability.

sigAbs : Calculates the numerical absorption cross section in cm^2 at a given
	energy for the (l,m)^th mode of a system of n concentric spherical 
	surfaces. The absorption cross section is returned as a number.

*)

epsDrude[hw0_,hwp0_,hgam0_,epsInf0_] := Module[{hw = hw0,hwp=hwp0,hgam=hgam0,
	epsInf=epsInf0},

	Return[epsInf - hwp^2/(hw^2 + I * hw * hgam)]
]

normDist[x0_,mu0_,sigma0_]:= Module[{x = x0,mu = mu0,sigma = sigma0},
	func = Exp[-(x - mu)^2/(2*sigma^2)];
	Return[func]
] 

phiBCeq[n0_,l0_,m0_] := Module[{n=n0,l=l0,m=m0,Eq,EqList},
	Clear[a,rp,A,B,eps];
	B[1,l,m] = 0;
	A[n+2,l,m] = 0;
	eps[n+2] = eps[n+1];

	For[i=1,i<=n,i++,
		Eq[i,l,m] = A[i,l,m] + B[i,l,m]/a[i]^(2*l+1) - ( A[i+1,l,m] + 
			B[i+1,l,m]/a[i]^(2*l+1) )
	]; 

	Eq[n+1,l,m] = A[n+1,l,m] + B[n+1,l,m]/rp^(2*l+1) - ( A[n+2,l,m] + 
		B[n+2,l,m]/rp^(2*l+1) );
	EqList = {};

	For[i=1,i<=n+1,i++,
		AppendTo[EqList,Eq[i,l,m]]
	];
	Return[EqList]
]

dBCeq[n0_,l0_,m0_] := Module[{n=n0,l=l0,m=m0,Eq,EqList},
	ClearAll[a,rp,Ylmstar,A,B,eps];

	B[1,l,m] = 0;
	A[n+2,l,m] = 0;
	eps[n+2] = eps[n+1];

	For[i=1,i<=n,i++,
		Eq[i,l,m] = eps[i](A[i,l,m]*l*a[i]^(l-1)-B[i,l,m]*(l+1)/a[i]^(l+2)) - 
			eps[i+1](A[i+1,l,m]*l*a[i]^(l-1)-B[i+1,l,m]*(l+1)/a[i]^(l+2))
	];

	Eq[n+1,l,m] = eps[n+1](A[n+1,l,m]*l*rp^(l-1)-B[n+1,l,m]*(l+1)/rp^(l+2)) - 
		eps[n+2](A[n+2,l,m]*l*rp^(l-1)-B[n+2,l,m]*(l+1)/rp^(l+2)) - 
		4 * Pi * Ylmstar/rp^2;
	EqList = {};

	For[i=1,i<=n+1,i++,
		AppendTo[EqList,Eq[i,l,m]]
	];

	Return[EqList]
]

BCsolns[n0_,l0_,m0_] := Module[{n=n0,l=l0,m=m0,
	phiBClist,dBClist,BClist,zeroList,Alist,Blist},
	ClearAll[A,B,solns];

	phiBClist = phiBCeq[n,l,m];
	dBClist = dBCeq[n,l,m];
	BClist = Join[phiBClist,dBClist];
	zeroList = ConstantArray[0,Length[BClist]];

	Alist = {};
	Blist = {};

	For[i = 1,i <= n + 2,i++,
		Which[
			i == 1,
				AppendTo[Alist,A[i,l,m]],
			i == n+2, 
				AppendTo[Blist,B[i,l,m]],
			True,
				AppendTo[Alist,A[i,l,m]];
				AppendTo[Blist,B[i,l,m]];
		]
	];

	solns = Solve[BClist == zeroList,Join[Alist,Blist]];
	Return[solns[[1]]]
]

alphalm[n0_,l0_,m0_] := Module[{n=n0,l=l0,m=m0},
	ClearAll[a,rp,Ylmstar,eps,A,B,coeff];

	BCsolnsAssoc = BCsolns[n,l,m];

	coeff = B[n+2,l,m]/.BCsolnsAssoc;
	coeff = coeff - 4*Pi/(2*l+1) * rp * Ylmstar / eps[n+1];
	coeff = -coeff * eps[n+1] * (2*l+1)/(4*Pi*Ylmstar*a[n]^(l-1))*rp^(l+1);

	Return[FullSimplify[coeff]]
]

resFinder[alpha0_,n0_,l0_,m0_,params0_,epsb0_]:=Module[{alpha = alpha0,
	n = n0,l = l0,m = m0,
	params = params0,epsb = epsb0,
	chilm,
	L0,l1,l2,l3,l4,
	num,den,resEn},
	ClearAll[a,En,aList,epsList,aRules,epsArray,rules,scale];

	Which[
		FreeQ[alpha,a[n]],Return[0],
		True,Null
	];

	aList = Array[a,n];
	epsList = Array[eps,n];
	aRules = Thread[ aList -> params[[All,1]] ];
	epsArray = {};

	Do[
		AppendTo[epsArray,epsDrude[En,params[[i,2]],params[[i,3]],
			params[[i,4]] ]]
	,{i,n}];

	epsRules = Thread[epsList -> epsArray];
	rules = Join[aRules,epsRules,{eps[n+1] -> epsb}];
	chilm = alpha/a[n]^(l+2);
	chilm = chilm/.Dispatch[rules];
	chilm = FullSimplify[chilm];

	num = Numerator[chilm];
	den = Denominator[chilm];
	
	ClearAll[denResList,denList,denFact,numResList,numList,numFact,
		denPref,numPref,denRem,numRem,denListSimp,numListSimp,
		denFactSimp,numFactSimp];

	denResList = En/.Solve[den == 0, En];
	numResList = En/.Solve[num == 0, En];
	numList = {};
	denList = {};
	numListSimp = {};
	denListSimp = {};

	Do[
		Which[
			Re[ denResList[[j]] ] == 0.,
				Null,
			Re[ denResList[[j]] ] =!= 0.,
				AppendTo[denListSimp,En - denResList[[j]]];
		];
	,{j,Length[denResList]}];

	Do[
		Which[
			Re[ denResList[[j]] ] == 0.,
				Null,
			Re[ denResList[[j]] ] =!= 0.,
				AppendTo[numListSimp,En - numResList[[j]]];
		];
	,{j,Length[numResList]}];

	Do[
		AppendTo[denList,En - denResList[[j]]];
	,{j,Length[denResList]}];

	Do[
		AppendTo[numList,En - numResList[[j]]];
	,{j,Length[numResList]}];

	denFact = Apply[Times,denList];
	numFact = Apply[Times,numList];
	denFactSimp = Apply[Times,denListSimp];
	numFactSimp = Apply[Times,numListSimp];
	denPref = PolynomialQuotientRemainder[den,denFact,En][[1]];
	numPref = PolynomialQuotientRemainder[num,numFact,En][[1]];
	denRem = Chop[PolynomialQuotientRemainder[den,denFact,En]][[2]];
	numRem = Chop[PolynomialQuotientRemainder[num,numFact,En]][[2]];

	(*** Uncomment for debugging ***)
	(* Print["---"];
	Print["chilm = "];
	Print[chilm];
	Print["Simplified (polynomial) chilm = "];
	Print[ Expand[num]/Expand[den] ];
	Print["---"]
	Print["Roots of numerator = "];
	Print[numResList];
	Print["Roots of denominator = "];
	Print[denResList];
	Print["Numerator in factored form = "];
	Print[numFact];
	Print["Denominator in factored form = "];
	Print[denFact];
	Print["Numerator prefactor = "];
	Print[numPref];
	Print["Denominator prefactor = "];
	Print[denPref];
	Print["Numerator remainder = "];
	Print[numRem];
	Print["Denominator remainder = "];
	Print[denRem];
	Print["---"];
	Print["chilm in separated fraction form (manual) = "];
	Print[numPref/denPref * Chop[Apart[numFact/denFact,En]]];
	Print["chilm in separated fraction form (manual with simplification) = "];
	Print[numPref/denPref * Chop[Apart[numFactSimp/denFactSimp,En]]];
	Print["chilm in separated fraction form (automatic) = "];
	Print[Apart[chilm,En]]; *)
 
	l1 = {};
	l2 = {};
	l3 = {};
	l4 = {};

	chilmSimp = numPref/denPref*Apart[numFactSimp/denFactSimp,En];
	L0 = Apply[List,chilmSimp];
	L0 = L0[[1]] * Apply[ List,L0[[2]] ];

	Do[
		den = Denominator[Apply[List,L0][[i]]];
		num = Numerator[Apply[List,L0][[i]]];

 		Which[
			den === 1,resEn = Nothing,
			True, resEn = Apply[List,den][[1]];
		];

		Which[
			Re[resEn] > 0.,
				AppendTo[l1,num];
				AppendTo[l2,Re[resEn]];
				AppendTo[l3,Im[resEn]];
				AppendTo[l4,Sqrt[Re[resEn]^2 + Im[resEn]^2]],
			True,Null
		]; 
	,{i,Length[L0]}];

	Return[{l1,l2,l3,l4}];
]

massFinder[alpha0_,n0_,l0_,m0_,params0_,epsb0_]:=Module[{alpha = alpha0,
	n = n0,l = l0,m = m0,
	params = params0,epsb = epsb0,
	L0,L1,p,
	resEn,natEn,A,f,M},

	L0 = {};
	L1 = {};
	p = resFinder[alpha,n,l,m,params,epsb];

	Do[
		resEn = p[[2,i]];
		natEn = p[[4,i]];
		A = p[[1,i]];
		f = 2 * resEn * A / natEn^2;
		M = e^2 * hbareV^2 / ( params[[n,1]]^(2l+1)*natEn^2*Abs[f] );
		AppendTo[L0,f];
		AppendTo[L1,M];
	,{i,Length[ p[[1,All]] ]}];

	Return[{L0,L1}]
]

parFinder[alpha0_,n0_,l0_,m0_,params0_,epsb0_]:=Module[{alpha = alpha0,
	n = n0,l = l0,m = m0,
	params = params0,epsb = epsb0,
	p0,p1},

	p0 = resFinder[alpha,n,l,m,params,epsb];
	p1 = massFinder[alpha,n,l,m,params,epsb];

	Do[
		AppendTo[ p0,p1[[i]] ];
	,{i,Length[ p1[[All,1]] ]}];

	Return[ Part[p0,2;;] ];
]

sigAbsSymb[alpha0_,n0_,l0_,m0_,params0_,epsb0_]:=Module[{alpha = alpha0,
	n = n0,l = l0,m = m0,
	params = params0,epsb = epsb0,
	p,L0,sig,dampEn,natEn,M},
	ClearAll[En];

	p = parFinder[alpha,n,l,m,params,epsb];
	L0 = {};

	Do[
		dampEn = 2 * p[[2,i]];
		natEn = p[[3,i]];
		M = p[[5,i]];
		sig = 4*Pi*hbareV*En/c * En*dampEn/((natEn^2-En^2)^2+En^2*dampEn^2) *
			e^2/M / (10^(-7))^2;
		AppendTo[L0,sig]
	,{i,Length[ p[[1,All]] ]}];

	Return[L0]
]

sigAbs[hw0_,sigAbsSymb0_]:=Module[{hw = hw0,sigSymb = sigAbsSymb0},
	ClearAll[En,sig];

	sig = sigSymb/.Dispatch[En->hw];
	Return[sig]
]

(****************************************************************************)
(************************* Defining Test Functions **************************)
(****************************************************************************)

massFinderTest[alpha0_,n0_,l0_,m0_,params0_,epsb0_]:=Module[{alpha = alpha0,
	n = n0,l = l0,m = m0,
	params = params0,epsb = epsb0,
	L0,L1,p,
	resEn,natEn,A,f,M},

	L0 = {};
	L1 = {};
	p = resFinderTest[alpha,n,l,m,params,epsb];

	Do[
		resEn = p[[2,i]];
		natEn = p[[4,i]];
		A = p[[1,i]];
		f = 2 * resEn * A / natEn^2;
		M = e^2 * hbareV^2 / ( params[[n,1]]^(2l+1)*natEn^2*Abs[f] );
		AppendTo[L0,f];
		AppendTo[L1,M];
	,{i,Length[ p[[1,All]] ]}];

	Return[{L0,L1}]
]

parFinderTest[alpha0_,n0_,l0_,m0_,params0_,epsb0_]:=Module[{alpha = alpha0,
	n = n0,l = l0,m = m0,
	params = params0,epsb = epsb0,
	p0,p1},

	p0 = resFinderTest[alpha,n,l,m,params,epsb];
	p1 = massFinderTest[alpha,n,l,m,params,epsb];

	Do[
		AppendTo[ p0,p1[[i]] ];
	,{i,Length[ p1[[All,1]] ]}];

	Return[ Part[p0,2;;] ];
]

sigAbsSymbTest[alpha0_,n0_,l0_,m0_,params0_,epsb0_]:=Module[{alpha = alpha0,
	n = n0,l = l0,m = m0,
	params = params0,epsb = epsb0,
	p,L0,sig,dampEn,natEn,M},
	ClearAll[En];

	p = parFinderTest[alpha,n,l,m,params,epsb];
	L0 = {};

	Do[
		dampEn = 2 * p[[2,i]];
		natEn = p[[3,i]];
		M = p[[5,i]];
		sig = 4*Pi*hbareV*En/c * En*dampEn/((natEn^2-En^2)^2+En^2*dampEn^2) *
			e^2/M / (10^(-7))^2;
		AppendTo[L0,sig]
	,{i,Length[ p[[1,All]] ]}];

	Return[L0]
]

EndPackage[]
