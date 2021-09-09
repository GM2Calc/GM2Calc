If[Install["../bin/gm2calc.mx"] === $Failed,
   Quit[1];
  ];

errors = 0;
passed = 0;

TestClose[val1_?NumericQ, val2_?NumericQ, eps_:$MachineEpsilon] :=
    If[Abs[val1 - val2] > eps,
       Print["Error: expressions are not equal: ",
             InputForm[val1], " =!= ", InputForm[val2]];
       errors++,
       passed++
      ];

TestClose[val1_, val2_, eps_:$MachineEpsilon] := (
    Print["Non-numeric value found: val1 = ", val1, ", val2 = ", val2, ", eps = ", eps];
    errors++
    );

TestEqual[expr1_, expr2_] :=
    If[expr1 === expr2,
       passed++,
       Print["Error: expressions are not equal: ",
             InputForm[expr1], " =!= ", InputForm[expr2]];
       errors++
    ]

GM2CalcSetSMParameters[
    alphaMZ -> 0.00775531,
    alpha0 -> 0.00729735,
    alphaS -> 0.1184,
    MW -> 80.385,
    MZ -> 91.1876,
    MT -> 173.34,
    mbmb -> 4.18,
    ML -> 1.777,
    MM -> 0.1056583715];

(* default SLHA point *)
point = { MSvmL -> 5.18860573 10^2,
          MSm -> {5.05095249 10^2, 5.25187016 10^2},
          MChi -> {2.01611468 10^2, 4.10040273 10^2, -5.16529941 10^2, 5.45628749 10^2},
          MCha -> {4.0998989 10^2, 5.46057190 10^2},
          MAh -> 1.5 10^3,
          TB -> 40,
          Mu -> 500,
          MassB -> 200,
          MassWB -> 400,
          MassG -> 2000,
          mq2 -> 7000^2 IdentityMatrix[3],
          ml2 -> 500^2 IdentityMatrix[3],
          mu2 -> 7000^2 IdentityMatrix[3],
          md2 -> 7000^2 IdentityMatrix[3],
          me2 -> 500^2 IdentityMatrix[3],
          Au -> 0 IdentityMatrix[3],
          Ad -> 0 IdentityMatrix[3],
          Ae -> 0 IdentityMatrix[3],
          Q -> 1000 };

{myAmu, myDamu} = {amu, Damu} /. GM2CalcAmuSLHAScheme[point];

TestClose[myAmu , 2.339244076106932*^-9];
TestClose[myDamu, 2.334226722205371*^-10];

(* default GM2Calc interface point *)
point = { MAh -> 1500,
          TB -> 10,
          Mu -> 350,
          MassB -> 150,
          MassWB -> 300,
          MassG -> 1000,
          mq2 -> 500^2 IdentityMatrix[3],
          ml2 -> 500^2 IdentityMatrix[3],
          mu2 -> 500^2 IdentityMatrix[3],
          md2 -> 500^2 IdentityMatrix[3],
          me2 -> 500^2 IdentityMatrix[3],
          Au -> 0 IdentityMatrix[3],
          Ad -> 0 IdentityMatrix[3],
          Ae -> 0 IdentityMatrix[3],
          Q -> 454.7 };

{myAmu, myDamu} = {amu, Damu} /. GM2CalcAmuGM2CalcScheme[point];

TestClose[myAmu , 7.964321357104185*^-10];
TestClose[myDamu, 2.309099582759784*^-10];

(* invalid point (GM2Calc interface) *)
point = {
    MAh    -> 1500,
    TB     -> 1000,
    Mu     -> 30000,
    MassB  -> 1000,
    MassWB -> -30000,
    MassG  -> 2000,
    mq2    -> 3000^2 IdentityMatrix[3],
    ml2    ->    5^2 IdentityMatrix[3],
    mu2    -> 3000^2 IdentityMatrix[3],
    md2    -> 3000^2 IdentityMatrix[3],
    me2    -> {{1000^2,0,0}, {0,1000^2,0}, {0,0,3000^2}},
    Au     -> 0 IdentityMatrix[3],
    Ad     -> 0 IdentityMatrix[3],
    Ae     -> 0 IdentityMatrix[3],
    Q      -> 866.360379
};

TestEqual[GM2CalcAmuGM2CalcScheme[point], {}];

(* invalid point w/o tan(beta) resummation (GM2Calc interface) *)
GM2CalcSetFlags[
    loopOrder -> 2,
    tanBetaResummation -> False,
    forceOutput -> False];

point = {
    MAh    -> 1500,
    TB     -> 1000,
    Mu     -> 30000,
    MassB  -> 1000,
    MassWB -> -30000,
    MassG  -> 2000,
    mq2    -> 3000^2 IdentityMatrix[3],
    ml2    -> {{1000^2,0,0}, {0,1000^2,0}, {0,0,3000^2}},
    mu2    -> 3000^2 IdentityMatrix[3],
    md2    -> 3000^2 IdentityMatrix[3],
    me2    -> {{1000^2,0,0}, {0,1000^2,0}, {0,0,3000^2}},
    Au     -> 0 IdentityMatrix[3],
    Ad     -> 0 IdentityMatrix[3],
    Ae     -> 0 IdentityMatrix[3],
    Q      -> 866.360379
};

{myAmu, myDamu} = {amu, Damu} /. GM2CalcAmuGM2CalcScheme[point];

TestEqual[myAmu, Indeterminate];

(* THDM mass basis point *)
GM2CalcSetFlags[loopOrder -> 2, runningCouplings -> False];

GM2CalcSetSMParameters[
    alpha0 -> 0.00729735,
    alphaMZ -> 1/128.94579,
    alphaS -> 0.1184,
    MhSM -> 125.09,
    MW -> 80.385,
    MZ -> 91.1876,
    MT -> 173.34,
    mcmc -> 1.28,
    mu2GeV -> 0.0022,
    mbmb -> 4.18,
    ms2GeV -> 0.096,
    md2GeV -> 0.0047,
    ML -> 1.77684,
    MM -> 0.1056583715,
    ME -> 0.000510998928,
    Mv1 -> 0,
    Mv2 -> 0,
    Mv3 -> 0,
    CKM -> IdentityMatrix[3] ];

point = {
    yukawaType        -> 2,
    Mhh               -> { 125, 400 },
    MAh               -> 420,
    MHp               -> 440,
    sinBetaMinusAlpha -> 0.999,
    lambda6           -> 0.2,
    lambda7           -> 0.1,
    TB                -> 3,
    m122              -> 200^2,
    zetau             -> 0,
    zetad             -> 0,
    zetal             -> 0,
    Piu                -> 0 IdentityMatrix[3],
    Pid                -> 0 IdentityMatrix[3],
    Pil                -> {{0,0,0}, {0,0,0}, {0,0,0}}
};

{myAmu, myDamu} = {amu, Damu} /. GM2CalcAmuTHDMMassBasis[point];

TestClose[myAmu, 1.8787546646721519*^-11];

(* THDM gauge basis point *)
GM2CalcSetFlags[loopOrder -> 2, runningCouplings -> False];

GM2CalcSetSMParameters[
    alpha0 -> 0.00729735,
    alphaMZ -> 1/128.94579,
    alphaS -> 0.1184,
    MhSM -> 125.09,
    MW -> 80.385,
    MZ -> 91.1876,
    MT -> 173.34,
    mcmc -> 1.28,
    mu2GeV -> 0.0022,
    mbmb -> 4.18,
    ms2GeV -> 0.096,
    md2GeV -> 0.0047,
    ML -> 1.77684,
    MM -> 0.1056583715,
    ME -> 0.000510998928,
    Mv1 -> 0,
    Mv2 -> 0,
    Mv3 -> 0,
    CKM -> IdentityMatrix[3] ];

point = {
    yukawaType        -> 2,
    lambda            -> { 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1 },
    TB                -> 3,
    m122              -> 200^2,
    zetau             -> 0,
    zetad             -> 0,
    zetal             -> 0,
    Piu                -> 0 IdentityMatrix[3],
    Pid                -> 0 IdentityMatrix[3],
    Pil                -> {{0,0,0}, {0,0,0}, {0,0,0}}
};

{myAmu, myDamu} = {amu, Damu} /. GM2CalcAmuTHDMGaugeBasis[point];

TestClose[myAmu, 3.2720909955357482*^-11];

Print[];
Print["Passed tests: [", passed, "/", passed + errors,"]"];

Quit[errors];
