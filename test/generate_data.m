(* calculate func[args] with precision prec *)
GeneratePoint[func_, prec_, args__] :=
    N[#, prec]& @ { args, func[args] }

(* calculate func[arg] with precision prec on a grid *)
GenerateGridData[func_, narg_, prec_, {min_, max_, step_}] :=
    Module[{vals = Table[x, {x, min, max, step}]},
           GeneratePoint[func, prec, Sequence @@ #]& /@ Tuples[vals, narg]
    ]

(* calculate func with precision prec on a grid close to unity *)
GenerateUnitData[func_, narg_, prec_, n_, frac_, dir_] :=
    Module[{vals = Table[1 + dir frac^k, {k, 1, n, 1}]},
           GeneratePoint[func, prec, Sequence @@ #]& /@ Tuples[vals, narg]
    ]

ExportData[func_, narg_, prec_, zero_] :=
    Module[{data, filename = ToString[func] <> ".txt", step = 1/20},
           Print["Generating data for " <> ToString[func]];
           data = Join[
               GenerateGridData[func, narg, prec, {zero, 5, step}],
               GenerateUnitData[func, narg, prec, prec, 1/10, -1],
               GenerateUnitData[func, narg, prec, prec, 1/10, +1]
           ];
           Print["Writing data to ", filename];
           Export[filename, data, "Table"];
    ]

$MaxExtraPrecision = 10000;
precision = 15;
zero = 0;

ExportData[F1C, 1, precision, 0];
ExportData[F2C, 1, precision, 0];
ExportData[F3C, 1, precision, 0];
ExportData[F4C, 1, precision, 0];

ExportData[F1N, 1, precision, 0];
ExportData[F2N, 1, precision, 0];
ExportData[F3N, 1, precision, 0];
ExportData[F4N, 1, precision, 0];

ExportData[G3, 1, precision, 0];
ExportData[G4, 1, precision, 0];

ExportData[Fa, 2, precision, 1/20];
ExportData[Fb, 2, precision, 1/20];
