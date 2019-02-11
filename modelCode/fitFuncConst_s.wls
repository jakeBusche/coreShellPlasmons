#!/usr/bin/env wolframscript
(* ::Package:: *)

(****************************************************************************)
(*************************** Initialization *********************************)
(****************************************************************************)

BeginPackage["fitFuncConst`"]

(****************************************************************************)
(*********************** Defining Physical Constants ************************)
(****************************************************************************)

c = 2.998*10^10. (* speed of light (cm/s) *)
cmks = 2.998*10^8 (* speed of light (m/s) *)
eps0 = 8.854187817*10^(-12.) (* permittivity of free space (C^2/(m^2*N)) *)
emks = 1.6021766208*10^(-19.) (* electron charge (C) *)
e = 4.8032042510 * 10^(-10) (* electron charge (statC) *)
nm2cm = 10^(-7.) (* cm / nm *)
cm2nm = 10^(7.) (* nm / cm *)
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

massFinder : Function that finds the oscillator strengths and effective masses
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

epsDrude[hw0_,hwp0_,hgam0_,epsInf0_] := Module[
    { hw = hw0,hwp = hwp0,hgam = hgam0,
	    epsInf = epsInf0 },

	Return[epsInf - hwp^2/(hw^2 + I * hw * hgam)]
]

phiBCeq[n0_,l0_,m0_,radList0_] := Module[
    { n = n0,l = l0,m = m0,
        radList = radList0,
        Eq,EqList },

	Clear[rp,A,B];
	B[1,l,m] = 0;
	A[n+2,l,m] = 0;

    (* The first n potential BCs are all similar *)
	For[i=1,i<=n,i++,
		Eq[i,l,m] = A[i,l,m] + B[i,l,m]/radList[[i]]^(2*l+1) - ( A[i+1,l,m] + 
			B[i+1,l,m]/radList[[i]]^(2*l+1) )
	]; 

    (* The n+1^th potential BC has a slightly different form *)
	Eq[n+1,l,m] = A[n+1,l,m] + B[n+1,l,m]/rp^(2*l+1) - ( A[n+2,l,m] + 
		B[n+2,l,m]/rp^(2*l+1) );
	EqList = {};

    (* Gather the BCs into a list of equations *)
	For[i=1,i<=n+1,i++,
		AppendTo[EqList,Eq[i,l,m]]
	];
	Return[EqList]
]

dBCeq[n0_,l0_,m0_,radList0_,epsList0_] := Module[
    { n = n0,l = l0,m = m0,
    	radList = radList0,
    	epsList = epsList0,
    	Eq,EqList },

	ClearAll[rp,Ylmstar,A,B];

	B[1,l,m] = 0;
	A[n+2,l,m] = 0;

    (* Length[epsList] must equal Length[radList] + 2 or the computation 
        must cease. *)
	Which[
        Length[epsList] == Length[radList] + 2,
            Null,
        True,
            Print["dBCeq : Number of elements in epsList must" <>
                " be equal to the number of elements in radList plus two."];
            Print["Exiting program."];
            Abort[];
	];

    (* If last two elements of epsList are not equal, overwrite the last one
        with the second-to-last value and alert the user. Otherwise, do
        nothing. *)
	Which[
	    epsList[[n+2]] != epsList[[n+1]], 
	        Print["dBCeq : Last element of epsList was overwritten to"<>
	            " equal second-to-last element. "];
	        epsList[[n+2]] = epsList[[n+1]],
	    True,
	        Null;
	];
	
    (* The first n displacement field BCs are all similar *)
	For[i=1,i<=n,i++,
		Eq[i,l,m] = epsList[[i]](A[i,l,m]*l*radList[[i]]^(l-1) - 
			B[i,l,m]*(l+1)/radList[[i]]^(l+2)) - epsList[[i+1]](
			A[i+1,l,m]*l*radList[[i]]^(l-1) - 
			B[i+1,l,m]*(l+1)/radList[[i]]^(l+2))
	];

    (* The n+1^th displacement field BC has a slightly different form *)
	Eq[n+1,l,m] = epsList[[n+1]](A[n+1,l,m]*l*rp^(l-1) - 
		B[n+1,l,m]*(l+1)/rp^(l+2)) - epsList[[n+2]](A[n+2,l,m]*l*rp^(l-1) - 
		B[n+2,l,m]*(l+1)/rp^(l+2)) - 4 * Pi * Ylmstar/rp^2;
	EqList = {};

    (* Gather the BCs into a list of equations *)
	For[i=1,i<=n+1,i++,
		AppendTo[EqList,Eq[i,l,m]]
	];

	Return[EqList]
]

BCsolns[n0_,l0_,m0_,radList0_,epsList0_] := Module[
    { n = n0,l = l0,m = m0,
    	radList = radList0,
    	epsList = epsList0,
	    phiBClist,dBClist,BClist,zeroList,Alist,Blist},
	ClearAll[A,B,solns];

	phiBClist = phiBCeq[n,l,m,radList];
	dBClist = dBCeq[n,l,m,radList,epsList];

	(* Join the lists of BCs for use in a later Solve[] call *)
	BClist = Join[phiBClist,dBClist];

	(* Build a list of zeros for use in a later Solve[call] *)
	zeroList = ConstantArray[0,Length[BClist]];

	Alist = {};
	Blist = {};

    (* Build lists of the coefficients A[i,l,m], B[i,l,m] for later use *)
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
    
    (* Solve the BCs using the proper list syntax for Solve[] *)
	solns = Solve[BClist == zeroList,Join[Alist,Blist]];
	Return[solns[[1]]]
]

alphalm[n0_,l0_,m0_,radList0_,epsList0_] := Module[
    { n = n0,l = l0,m = m0,
    	radList = radList0,
    	epsList = epsList0 },
	ClearAll[rp,Ylmstar,A,B,pol];

    (* Get the solutions of the BCs *)
	BCsolnsAssoc = BCsolns[n,l,m,radList,epsList];

    (* Get the solution of the Green function outside the particle and the
        source. Subtract off the direct Coulomb potential contribution. Divide
        away leading factors of Ylmstar, rp, etc., to get the expansion
        policient of the Green function (contains spectral response 
        only) times the radius of the particle to the l+2 power, i.e. the
        particle polarizability. *)
	pol = B[n+2,l,m]/.BCsolnsAssoc;
	pol = pol - 4*Pi/(2*l+1) * rp * Ylmstar / epsList[[n+1]];
	pol = -pol * epsList[[n+1]] * (2*l+1)/(4*Pi*Ylmstar*
		radList[[n]]^(l-1))*rp^(l+1);
    
    (* Numerical precision limits require the polarizability to be simplified
        and have hanging small factors removed before use. *)
	Return@Chop@FullSimplify[pol]

]

resFinder[alpha0_,n0_,l0_,m0_,radList0_,epsList0_]:=Module[{alpha = alpha0,
	n = n0,l = l0,m = m0,
	radList = radList0,
	epsList = epsList0,
	chilm,
	L0,l1,l2,l3,l4,
	num,den,resEn},
	ClearAll[En,scale,alphaVar,
	    denResList,denList,denFact,numResList,numList,numFact,
		denPref,numPref,denRem,numRem,denListSimp,numListSimp,
		denFactSimp,numFactSimp];

	alphaVar = Variables@Level[alpha,{-1}];
	Which[Length[alphaVar] > 1,
	        Print["resFinder : More than one variable detected in" <> 
	            " polarizability."];
	        Print["Exiting program."];
	        Abort[],
        True,
            Null;
	];

	(* Replace polarizability variable (from dielectric function) with "En"
	     for ease of use. *)
	alpha = alpha/.alphaVar[[1]] -> En;

    (* Divide radius information from the polarizability to get the unitless
        spectral response *)
	chilm = alpha/Last[radList]^(l+2);
	chilm = FullSimplify[chilm];

    (* Get the roots of the numerator and denominator of chilm. Place them 
        into lists. *)
	num = Numerator[chilm];
	den = Denominator[chilm];
	denResList = En/.Solve[den == 0, En];
	numResList = En/.Solve[num == 0, En];

    (* Represent the numerator and denominator of chilm in factored form. 
        First store the terms in lists then as products. For numListSimp and
        denListSimp, store only terms that have a root with a nonzero real 
        part. *)
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

    (* Find the polynomial quotients of the numerator and denominator and 
        their simplified (all terms with imaginary roots removed) factored 
        forms. Also, find the polynomial remainders (which should be small 
        or zero) for debugging purposes. *)
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

    (* The spectral response with a factored numerator and denominator, 
        roots with nonzero real parts removed, and proper prefactors 
        multiplied in represented as a sum of fractions. *) 
	chilmSimp = numPref/denPref*Apart[numFactSimp/denFactSimp,En]; 

	(* Convert the spectral response from a sum to a list *)
	L0 = Apply[List,chilmSimp];

    (* The Apply[List,chilmSimp] command places an overall prefactors into 
        the first element of the resulting list and the rest of the 
        expression into the second element. Multiply the prefactor back in 
        and convert the rest of the expression into a list. *) 
	L0 = L0[[1]] * Apply[ List,L0[[2]] ];

    (* For each term in chilm, if the term is not the constant contribution
        (see notes for details), record the numerator (oscillator strength),
        real and imaginary roots (natural and damping energy) of the
        denominator, and resonance energy (see notes). Return a matrix of
        these quantities. *)
	Do[
		den = Denominator[Apply[List,L0][[i]]];
		num = Numerator[Apply[List,L0][[i]]];

 		Which[
			den === 1,resEn = Nothing,
			True, resEn = Apply[List,den][[1]];
		];

        (* Mathematica likes to evaluate the terms in chilm such that their 
            roots have positive imaginary parts. Thus, for each term with a
            root of positive real part, there is a root of negative real
            part. Store the numerator of each term and the natural, damping,
            and resonance energies of each. *)
		Which[
			Re[resEn] < 0. && Im[resEn] > 0.,
			    resEn = -resEn;
				AppendTo[l1,num];
				AppendTo[l2,Re[resEn]];
				AppendTo[l3,-2*Im[resEn]];
				AppendTo[l4,Chop[Sqrt[Im[resEn]^2+(resEn-I*Im[resEn])^2]]],
			True,Null
		]; 
	,{i,Length[L0]}];

	Return[{l1,l2,l3,l4}];
]

massFinder[alpha0_,n0_,l0_,m0_,radList0_,epsList0_]:=Module[
    { alpha = alpha0,
	n = n0,l = l0,m = m0,
	radList = radList0,epsList = epsList0,
	L0,L1,p,
	resEn,natEn,A,f,M},

	L0 = {};
	L1 = {};

	(* Get matrix of oscillator parameters *)
	p = resFinder[alpha,n,l,m,radList,epsList];

    (* Define oscillator strengths and effective masses and store them. *)
	Do[
		resEn = p[[2,i]];
		natEn = p[[4,i]];
		f = p[[1,i]];
		M = e^2 * cm2nm^3 / ( Last[radList]^(l+2) * 
			Abs[f]/hbareV * 2*resEn/hbareV );
		AppendTo[L0,f];
		AppendTo[L1,M];
	,{i,Length[ p[[1,All]] ]}];

	Return[{L0,L1}]
]

parFinder[alpha0_,n0_,l0_,m0_,radList0_,epsList0_]:=Module[
    { alpha = alpha0,
	n = n0,l = l0,m = m0,
	radList = radList0,epsList = epsList0,
	p0,p1},

    (* Get the natural, damping, and resonance energies as well as the 
        oscillator strengths and effective masses of the oscillators of the 
        particle. Store them nicely. *)

	p0 = resFinder[alpha,n,l,m,radList,epsList];
	p1 = massFinder[alpha,n,l,m,radList,epsList];

	Do[
		AppendTo[ p0,p1[[i]] ];
	,{i,Length[ p1[[All,1]] ]}];

	Return[ Part[p0,2;;] ];
]

sigAbsSymb[alpha0_,n0_,l0_,m0_,radList0_,epsList0_]:=Module[
    { alpha = alpha0,
	n = n0,l = l0,m = m0,
	radList = radList0,epsList = epsList0,
	p,L0,sig,dampEn,natEn,M},
	ClearAll[En];

    (* Get the oscillator parameters for each oscillator of the particle. *)
	p = parFinder[alpha,n,l,m,radList,epsList];
	L0 = {};

    (* Calculate the absorption cross-section of each of the oscillators using
        a nonradiative dipole oscillator model. The cross-section will be a 
        function of energy only. *)
	Do[
		dampEn = p[[2,i]];
		natEn = p[[3,i]];
		M = p[[5,i]];
		sig = 4*Pi * e^2*cm2nm^3 / (M * c*cm2nm) * En^2 * dampEn * hbareV / (
			(natEn^2 - En^2)^2 + En^2*dampEn^2);
		AppendTo[L0,sig]
	,{i,Length[ p[[1,All]] ]}];

	Return[L0]
]

sigAbs[hw0_,sigAbsSymb0_]:=Module[{hw = hw0,sigSymb = sigAbsSymb0},
	ClearAll[En,sig];

    (* Get the symbolic definition of the absorption cross-section and 
        evaluate it at a given energy. *)
	sig = sigSymb/.Dispatch[En->hw];
	Return[sig]
]

EndPackage[]
