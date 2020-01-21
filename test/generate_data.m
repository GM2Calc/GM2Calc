(* calculate func[arg] with precision prec *)
GeneratePoint[func_, prec_, arg_] :=
    N[#, prec]& @ { arg, func[arg] }

(* calculate func[arg] with precision prec on a grid *)
GenerateGridData[func_, prec_, {min_, max_, step_}] :=
    Table[
        GeneratePoint[func, prec, arg],
        {arg, min, max, step}
    ]

(* calculate func with precision prec on a grid close to unity *)
GenerateUnitData[func_, prec_, n_, frac_, dir_] :=
    Table[
        GeneratePoint[func, prec, 1 + dir frac^k],
        {k, 1, n, 1}
    ]

ExportData[func_, prec_] :=
    Module[{data, filename = ToString[func] <> ".txt"},
           Print["Generating data for " <> ToString[func]];
           data = Join[
               GenerateGridData[func, prec, {0, 5, 1/20}],
               GenerateUnitData[func, prec, prec, 1/10, -1],
               GenerateUnitData[func, prec, prec, 1/10, +1]
           ];
           Print["Writing data to ", filename];
           Export[filename, data, "Table"];
    ]

precision = 15

ExportData[F1C, precision];
ExportData[F2C, precision];
ExportData[F3C, precision];
ExportData[F4C, precision];

ExportData[F1N, precision];
ExportData[F2N, precision];
ExportData[F3N, precision];
ExportData[F4N, precision];

ExportData[G3, precision];
ExportData[G4, precision];
