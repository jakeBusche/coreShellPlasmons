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
erg2eV = 6.242*10^11 (* eV / erg *)
eV2erg = 1 / erg2eV (* erg / eV *)
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

(******************** Convenient simplification functions *******************)

epsDrude[hw0_,hwp0_,hgam0_,epsInf0_] := Module[
    { hw = hw0,hwp = hwp0,hgam = hgam0,
	    epsInf = epsInf0 },

	Return[epsInf - hwp^2/(hw^2 + I * hw * hgam)]
]

VariablesList[expr_] := Variables@Level[expr,{-1}]

FlagPrint[flag_,output__] := Which[
    flag == True,
        Print[output],
    True,
        Null
]

(********************** Observable generating functions *********************)

phiBCeq[n0_,l0_,m0_,radList0_] := Module[
    { n = n0,l = l0,m = m0,
        radList = radList0,
        Eq,EqList },

	Clear[rp,A,B];

	(* Two expansion coefficients must be zero from boundary conditions at 
	    infinity and at zero *)
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

    (* Two expansion coefficients must be zero from boundary conditions at 
	    infinity and at zero *)
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

(* Initialize optional arguments *)
Protect[verbose];
Options[BCsolns] = {verbose -> True};

BCsolns[n0_,l0_,m0_,radList0_,epsList0_,precision0_:100,
  opts:OptionsPattern[BCsolns]] := Module[
    { n = n0,l = l0,m = m0,
    	radList = radList0,
    	epsList = epsList0,
    	precision = precision0,
	    phiBClist,dBClist,BClist,zeroList,Alist,Blist,v},
	ClearAll[A,B,solns];

	v = OptionValue[verbose];

	(* Set the precision of the radius and dielectric values to avoid any two
        terms in Solve[] from becoming linearly dependent to the given
        numerical precision *)
	radList = SetPrecision[radList,precision];
	epsList = SetPrecision[epsList,precision];

    (* Get boundary condition equations for the potential and displacement
        field *)
	phiBClist = phiBCeq[n,l,m,radList];
	dBClist = dBCeq[n,l,m,radList,epsList];

	(* Join the lists of BCs for use in a later Solve[] call *)
	BClist = Join[phiBClist,dBClist];

	(* Build a list of zeros for use in a later Solve[] call *)
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
    
    (* Solve the BCs using the proper list syntax for Solve[]. Wrap the 
        Solve[] command in a Do loop and Which statement that tests 
        whether Solve[] will function properly with the radius and dielectric 
        function values evaluated with the given precision. If it does not, 
        produce an error message after n = nTrials failed attempts. *)
    ClearAll[chk,err,nTrials];

    nTrials = 25;

    (* Print update with verbose mode on *)
    FlagPrint[v,"BCsolns : Using boundary conditions to calculate" <>
        " expansion coefficients..."];
 
    Do[    
        chk = Quiet@Check[ Solve[BClist == zeroList,Join[Alist,Blist],
            WorkingPrecision -> precision], err ];

        Which[
            chk =!= err && chk =!= {},
                FlagPrint[v,"BCsolns : Solution satisfying boundary" <>
                    " conditions calculated normally."];
	            solns = chk;
	            Break[],
            chk === err && i =!= nTrials,
                Continue[],
            chk === err && i == nTrials,
                Print["BCsolns : Satisfactory solution not found after " <> 
                    ToString[nTrials] <> " attempts. Try increasing the " <>
                    " precision of the input values."];
        ]
    ,{i,1,nTrials}];

	Return[solns[[1]]]
]

alphalm[n0_,l0_,m0_,radList0_,epsList0_,precision0_:100,
  opts:OptionsPattern[BCsolns]] := Module[
    { n = n0,l = l0,m = m0,
    	radList = radList0,
    	epsList = epsList0,
    	precision = precision0,
    	v },
	ClearAll[rp,Ylmstar,A,B,pol,En];

	v = OptionValue[verbose];

    (* Get the solutions of the BCs *)
	BCsolnsAssoc = BCsolns[n,l,m,radList,epsList,precision,verbose -> v];

    FlagPrint[v,"alphalm : Obtaining polarizability from expansion" <> 
        " coefficient B_{n+2,l,m}..."];    

    (* Get the expansion of the (l,m)^th mode of Green function outside the 
        particle and the source. *)
    pol = B[n+2,l,m] /. BCsolnsAssoc;

    (* Find user-defined name of energy variable in pol. Let En be synonymous
        the energy variable for ease of use. *)
    En = (Select[VariablesList[pol], (# =!= rp && # =!= Ylmstar) &])[[1]];
    
    (* The output forms of the expansion coefficient B[n+2,l,m] are 
        different and depend on the numbers of shells in the particle ( and
        therefore, I assume, on Mathematica's ability to simplify the output
        expression before returning it). Further details of the strategies
        for handling the two cases can be found in the notes. *)

    (*** Uncomment for testing ***)
    (* Print["pol head = ",Head@pol];
    Print["pol terms # = ",Length@pol];
    Which[
        Head@pol === Times,
            Print["pol = ",ScientificForm[pol,3]],
        Head@pol === Plus,
            Print["pol = ",Short[ScientificForm[pol,3],3]];
    ]; *)
    
    Which[
        (* For low number of shells, Mathematica can simplify the expression
            for B[n+2,l,m] automatically and returns a product (see notes for 
            details). *)
        Head@pol === Times,
            ClearAll[polPref,polNum,polDen,polCoul,polCoulRem];

            (
            (* Get rid of symbolic prefactors of the polarizability 
                expression *)
            pol = Cancel[ pol * rp^(l + 1)/Ylmstar ];

            (* Get the (symbolic) numerator and denominator as well as the
                (numeric) prefactor of the polarizability *)
            polPref = Select[ pol, VariablesList[#] === {} & ];
            polDen = (Select[ pol, VariablesList[#] === {En} & ])^(-1);
            polNum = Collect[ 
                Select[ pol, Sort@VariablesList[#] === Sort@{En, rp} & ], 
                rp ];

            (*** Uncomment for testing ***)
            (* Print["Simplified pol = ",ScientificForm[pol,3]];
            Print["pref = ",polPref];
            Print["den = ",polDen];
            Print["num = ",polNum]; *)
            
            (* Simplify the numerator of the polarizability by collecting the
                terms that depend on rp. These should, when divided by the
                denominator and multiplied by the prefactor, give the term 
                of the expansion coefficient that contributes to the Coulomb
                potential of the test charge (see notes for details). *)
            polCoul = polPref * (PolynomialQuotientRemainder[ Select[ polNum,
                Sort@VariablesList[#] === Sort@{En,rp} & ], 
                polDen, En ])[[1]];

            (*** Uncomment for testing: the remainder of the division step 
                that is used to find the Coulomb potential should be 
                zero. ***)
            (* polCoulRem = polPref*(PolynomialQuotientRemainder[ 
                Select[ polNum, VariablesList[#] === {En,rp} & ],
                 polDen, En ])[[2]]; *)

            (* Collect the terms that do not contribute to the Coulomb
                potential of the test charge. These describe the nanoparticle
                response. *)    
            polNum = Plus@@(Select[ polNum, 
                Sort@VariablesList[#] =!= Sort@{En,rp} & ]);

            (*** Uncomment for testing ***)
            (* Print["simp num = ",polNum];
            Print["coul term = ",polCoul];
            Print["coul term reaminder = ",polCoulRem]; *)

            (* Algebra steps to get polarizability (see notes for details) *)
            pol = polPref * polNum / polDen;
            pol = (2*l + 1)/(4 * Pi) * epsList[[n+2]]/radList[[n]]^(2*l+1)* 
                pol;
            pol = radList[[n]]^(l+2) * pol
            ),

        (* For large numbers of shells, Mathematica doesn't attempt to 
            simplify the expression for B[n+2,l,m] and returns a sum (see 
            notes for details) *)
        Head@pol === Plus && Length@pol == 3,
            ClearAll[polCoul];

            (
            
            (* Get the term of B[n+2,l,m] that describes the Coulomb 
                potential of the test charge. Select should return the term 
                as the only member of a Length == 1 list. The forced 
                numerical evaluation with 3 digit precision avoids over-
                pickiness by Select when pattern matching. *)
            polCoul = Select[ List@@pol, 
                N[#,3] === N[ 4*Pi/(2*l+1)/epsList[[n+2]] * rp * Ylmstar,
                    3 ] & ];

            (*** Uncomment for testing ***)
            (* Print["coul term = ",polCoul]; *)

            (* Check to see if polCoul is a single term in a list, as should
	            be the case. If so, subtract it from pol. If not, print an
	            error message and return Null. *) 
            Which[
                Length[ polCoul ] == 1,
                    polCoul = polCoul[[1]];
                    pol = pol - polCoul,
                True,
                    Print["alphalm : Vacuum term did not separate properly."];
            ];
            
            (* Multiply out the symbolic dependence on rp and Ylmstar from 
                the two terms in pol, use Cancel to simplify, and use Together
                to produce a single fraction term. *)
            pol = Together[ Cancel[ rp^(l+1)/Ylmstar * pol ] ];

            (*** Uncomment for testing ***)
            (* Print["pol = ",pol]; *)

            (* Obtain resonant term (n(\omega)/d(\omega)) from B[n+2,l,m] 
	            (see notes for details). *)
            pol = (2*l+1)/(4*Pi)*epsList[[n+2]]/radList[[n]] * pol;

            (* Obtain polarizability from resonant term *)
            pol = radList[[n]]^(l+2) * pol;
            ),

        (* For moderate numbers of shells, Mathematica performs some of its
            simplification procedures but not all of them. The expression for
            B[n+2,l,m] is given as two terms in a sum, one of which contains 
            the frequency dependence, and the other which contains the 
            Coulomb term *)
        Head@pol === Plus && Length@pol == 2,
            ClearAll[pol1,pol2,polConstMag];

            (* The first term of pol contains the frequency-dependent 
                numerator and denominator. The second contains a constant and
                the Coulomb term modulo a factor of Ylmstar/rp^(l+1) *)
            pol1 = pol[[1]];
            pol2 = pol[[2]];

            pol1 = Collect[Cancel[pol1 * rp^(l+1)/Ylmstar],rp];
            pol2 = Collect[Cancel[pol2 * rp^(l+1)/Ylmstar],rp];

            (* When subtracting the Coulomb term from B[n+2,l,m], make sure
                no hanging small factors remain by chopping any terms smaller
                than 10^(-10) times the order of magnitude of the constant 
                term *)
            polConstMag = Floor@Log[10, Abs@Select[pol2,
                VariablesList[#] === {}& ] ];
            pol2 = Chop[pol2 - 4*Pi/(2*l+1)/epsList[[n+2]]*rp^(2*l+1),
                polConstMag * 10^(-10)];

            (* Add the constant and frequency dependent pieces of B[n+2,l,m] 
                and multiply by the proper factors to build the 
                polarizability (see notes for details) *)
            pol = Together[pol1 + pol2];
            pol = pol * (2*l+1)/(4*Pi)*epsList[[n+2]]/radList[[n]]^(2*l+1);
            pol = radList[[n]]^(l+2) * pol;

    ];

    (* Check to see if the algebra above was done improperly. If any 
        variables, besides the energy variable remain, the calculation will
        fail. *)
    Which[
        Length@VariablesList[pol] > 1,
            Print["alphalm : Polarizability was not obtained properly."];
            Print["alphalm : Symbols in polarizability are ",
                VariablesList[pol]];
            Return[Null],
        True,
            FlagPrint[v,"alphalm : Polarizability calculated normally."];
    ];
    
    (* Replace the energy variable with En *)
    pol = pol/.En -> En;

	Return[pol];

]

resFinder[alpha0_,n0_,l0_,m0_,radList0_,epsList0_,precision0_:100,
  opts:OptionsPattern[BCsolns]]:= Module[{alpha = alpha0,
	n = n0,l = l0,m = m0,
	radList = radList0,
	epsList = epsList0,
	precision = precision0,
	chilm,
	L0,l1,l2,l3,l4,
	num,den,resEn,v},
	ClearAll[En,scale,alphaVar,
	    denResList,denList,denFact,numResList,numList,numFact,
		denPref,numPref,denRem,numRem,denListSimp,numListSimp,
		denFactSimp,numFactSimp];

    v = OptionValue[verbose];

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

    (* Set precision of radius and dielectric function values *)
	radList = SetPrecision[radList,precision];
	epsList = SetPrecision[epsList,precision];

    (* Divide radius information from the polarizability to get the unitless
        spectral response *)
	chilm = alpha/Last[radList]^(l+2);

	FlagPrint[v,"resFinder : Calculating mode oscillator strengths and" <>
        " complex resonance frequencies..."];

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
	Print[numPref/denPref * Apart[numFact/denFact,En]];
	Print["chilm in separated fraction form (manual with simplification) = "];
	Print[numPref/denPref * Apart[numFactSimp/denFactSimp,En]];
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

	Which[
		VariablesList[{l1,l2,l3,l4}] === {},
        
            FlagPrint[v,"resFinder : Mode oscillator strengths and "<> 
                        "complex resonance frequencies" <>
	                    " calculated normally."],
        True,
            Print["resFinder : Mode oscillator strengths and " <> 
                "complex resonance frequencies" <>
                "were not obtained properly."];
            Print["resFinder : Symbols in oscillator strengths and " <> 
                "resonance frequencies" <>
                " are ",VariablesList[{l1,l2,l3,l4}]]
	];

    

	Return[{l1,l2,l3,l4}];
]

massFinder[alpha0_,n0_,l0_,m0_,radList0_,epsList0_,precision0_:100,
  opts:OptionsPattern[BCsolns]]:=Module[
    { alpha = alpha0,
	n = n0,l = l0,m = m0,
	radList = radList0,epsList = epsList0,
	precision = precision0,
	L0,L1,p,
	resEn,natEn,A,f,M,v},
    v = OptionValue[verbose];

	L0 = {};
	L1 = {};

	(* Set precision of radius and dielectric function values *)
	radList = SetPrecision[radList,precision];
	epsList = SetPrecision[epsList,precision];

	FlagPrint[v,"massFinder : Calculating mode masses..."];

	(* Get matrix of oscillator parameters *)
	p = resFinder[alpha,n,l,m,radList,epsList,precision,verbose -> False];

    (* Define oscillator strengths and effective masses and store them. *)
	Do[
		resEn = p[[2,i]];
		natEn = p[[4,i]];
		f = p[[1,i]];
		M = e^2 * cm2nm^3 / ( 2* Last[radList]^3 * 
			Abs[f/hbareV] * resEn/hbareV);
		AppendTo[L0,f];
		AppendTo[L1,M];
	,{i,Length[ p[[1,All]] ]}];

	Which[
	    VariablesList[{L1}] === {},
	        FlagPrint[v,"massFinder : Mode masses calculated normally."],
        True,
            Print["massFinder : Mode" <> 
                " masses were not obtained properly."];
            Print["massFinder : Symbols in mode masses" <>
                " are ", VariablesList[{L1}]]
	];

	Return[{L0,L1}]
]

parFinder[alpha0_,n0_,l0_,m0_,radList0_,epsList0_,precision0_:100,
  opts:OptionsPattern[BCsolns]]:=Module[
    { alpha = alpha0,
	n = n0,l = l0,m = m0,
	radList = radList0,
	epsList = epsList0,
	precision = precision0,
	p0,p1,v},
	v = OptionValue[verbose];

	(* Set precision of radius and dielectric function values *)
	radList = SetPrecision[radList,precision];
	epsList = SetPrecision[epsList,precision];

    FlagPrint[v,"parFinder : Calculating mode oscillator parameters..."];

    (* Get the natural, damping, and resonance energies as well as the 
        oscillator strengths and effective masses of the oscillators of the 
        particle. Store them nicely. *)

	p0 = resFinder[alpha,n,l,m,radList,epsList,precision];
	p1 = massFinder[alpha,n,l,m,radList,epsList,precision];

	Do[
		AppendTo[ p0,p1[[i]] ];
	,{i,Length[ p1[[All,1]] ]}];

	Which[
        VariablesList[Part[p0,2;;]] === {},
            FlagPrint[v,"parFinder : Oscillator parameters calculated" <>
                " normally."],
        True,
            Print["parFinder : Oscillator parameters were not obtained" <>
                " properly."];
            Print["parFinder : Variables in oscillator parameters are ",
                VariablesList[Part[p0,2;;]]]
	];

	Return[ Part[p0,2;;] ];
]

sigAbsSymb[alpha0_,n0_,l0_,m0_,radList0_,epsList0_,precision0_,
  opts:OptionsPattern[BCsolns]]:=Module[
    { alpha = alpha0,
	n = n0,l = l0,m = m0,
	radList = radList0,
	epsList = epsList0,
	precision = precision0,
	p,L0,sig,dampEn,natEn,M},
	ClearAll[En];
	v = OptionValue[verbose];

	(* Set precision of radius and dielectric function values *)
	radList = SetPrecision[radList,precision];
	epsList = SetPrecision[epsList,precision];

    FlagPrint[v,"sigAbsSymb : Calculating symbolic absorption " <> 
        "cross-section..."];

    (* Get the oscillator parameters for each oscillator of the particle. *)
	p = parFinder[alpha,n,l,m,radList,epsList,precision,verbose -> v];
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

    (* If the number of terms in the absorption cross-section is not equal to
        the number of surfaces in the model, print an error *)

    FlagPrint[v && Length@VariablesList[L0] == 1 && Length@L0 == n,
        "sigAbsSymb : Symbolic absorption cross-section" <> 
            " calculated normally."];

    FlagPrint[Length@VariablesList[L0] == 1 && Length@L0 =!= n,
        "sigAbsSymb : Number of terms in symbolic absorption " <> 
            "cross-section does not match the chosen number of surfaces." <>
            " A higher calculation precision may be needed."];

    FlagPrint[Length@VariablesList[L0] =!= 1,
        "sigAbsSymb : Symbolic absorption cross-section was not" <> 
                " obtained properly. Variables in symbolic absorption" <> 
                " cross-section are ",VariablesList[L0]];

    Return[L0]

]

sigAbs[hw0_,sigAbsSymb0_]:=Module[{hw = hw0,sigSymb = sigAbsSymb0},
	ClearAll[En,sig];

    (* Get the symbolic definition of the absorption cross-section and 
        evaluate it at a given energy. *)
	sig = sigSymb/.Dispatch[En -> hw];
	Return[sig]
]

EndPackage[]