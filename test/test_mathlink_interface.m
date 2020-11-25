If[Install["../bin/gm2calc.mx"] === $Failed,
   Quit[1];
  ];

errors = 0;
passed = 0;

TestClose[val1_?NumericQ, val2_?NumericQ, eps_:10^-10] :=
    If[val1 - val2 > eps,
       Print["Error: expressions are not equal: ",
             InputForm[val1], " =!= ", InputForm[val2]];
       errors++,
       passed++
      ];

TestClose[val1_, val2_, eps_:10^-10] := (
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

Print[];
Print["Passed tests: [", passed, "/", passed + errors,"]"];

Quit[errors];
