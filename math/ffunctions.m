(* Loop functions for (g-2)_mu *)

Li2[x_] := PolyLog[2, x]

(* arxiv:1003.5820, Eq.(13) *)
F1C[0] := 4

F1C[1] := 1

F1C[x_] := (2/(1 - x)^4)*(2 + 3*x - 6*x^2 + x^3 + 6*x*Log[x])

(* arxiv:1003.5820, Eq.(14) *)
F2C[1] := 1

F2C[x_] := (3/(2*(1 - x)^3))*(-3 + 4*x - x^2 - 2*Log[x])

(* arxiv:1003.5820, Eq.(37) *)
F3C[1] := 1

F3C[x_] := (4/(141*(1 - x)^4))*(
    (1 - x)*(151*x^2 - 335*x + 592)
    + 6*(21*x^3 - 108*x^2 - 93*x + 50)*Log[x]
    - 54*x*(x^2 - 2*x - 2)*Log[x]^2
    - 108*x*(x^2 - 2*x + 12)*Li2[1 - x]
)

(* arxiv:1003.5820, Eq.(38) *)
F4C[1] := 1

F4C[x_] := (-9/(122*(1 - x)^3))*(
    8*(x^2 - 3*x + 2)
    + (11*x^2 - 40*x + 5)*Log[x]
    - 2*(x^2 - 2*x - 2)*Log[x]^2
    - 4*(x^2 - 2*x + 9)*Li2[1 - x]
)

(* arxiv:1003.5820, Eq.(15) *)
F1N[0] := 2

F1N[1] := 1

F1N[x_] := (2/(1 - x)^4)*(1 - 6*x + 3*x^2 + 2*x^3 - 6*x^2*Log[x])

(* arxiv:1003.5820, Eq.(16) *)
F2N[0] := 3

F2N[1] := 1

F2N[x_] := (3/(1 - x)^3)*(1 - x^2 + 2*x*Log[x])

(* arxiv:1003.5820, Eq.(39) *)
F3N[0] := 8/105

F3N[1] := 1

F3N[x_] := (4/(105*(1 - x)^4))*(
    (1 - x)*(-97*x^2 - 529*x + 2)
    + 6*x^2*(13*x + 81)*Log[x]
    + 108*x*(7*x + 4)*Li2[1 - x]
)

(* arxiv:1003.5820, Eq.(40) *)
F4N[0] := (-3*(-9 + Pi^2))/4

F4N[1] := 1

F4N[x_] := (-(9/4/(1 - x)^3))*(
    (x + 3)*(x*Log[x] + x - 1)
    + (6*x + 2)*Li2[1 - x]
)

(* arxiv:1311.1775, Eq.(6.3a) *)
Fa[1, 1] := 1/4

Fa[x_, x_] := (2 + 3 x - 6 x^2 + x^3 + 6 x Log[x])/(2 (-1 + x)^4 x)

Fa[x_, y_] := -(G3[x] - G3[y])/(x - y)

(* arxiv:1311.1775, Eq.(6.3b) *)
Fb[1, 1] := 1/12

Fb[x_, x_] := (-5 + 4 x + x^2 - 2 Log[x] - 4 x Log[x])/(2 (-1 + x)^4)

Fb[x_, y_] := -(G4[x] - G4[y])/(x - y)

(* arxiv:1311.1775, Eq.(6.4a) *)
G3[1] := 1/3

G3[x_] := 1/(2 (x - 1)^3) ((x - 1)(x - 3) + 2 Log[x])

(* arxiv:1311.1775, Eq.(6.4b) *)
G4[1] := 1/6

G4[x_] := 1/(2 (x - 1)^3) ((x - 1)(x + 1) - 2 x Log[x])

(* I[a,b,c] with squared arguments *)
I2abc[a_, a_, a_] := 1/(2 a)

I2abc[a_, a_, c_] := (a - c - c*Log[a/c])/(a - c)^2

I2abc[a_, c_, a_] := I2abc[a, a, c]

I2abc[c_, a_, a_] := I2abc[a, a, c]

I2abc[a_, b_, c_] :=
    (a b Log[a/b] + b c Log[b/c] + c a Log[c/a])/((a - b) (b - c) (a - c))

(* I[a,b,c] with non-squared arguments *)
Iabc[a_, b_, c_] := I2abc[a^2, b^2, c^2]

(* arxiv:hep-ph/0609168, Eq.(70) *)
fPS[z_] := Module[{y = Sqrt[1 - 4z]}, 2z/y (PolyLog[2, 1 - (1-y)/(2z)] - PolyLog[2, 1 - (1+y)/(2z)])]

(* arxiv:hep-ph/0609168, Eq.(71) *)
fS[z_] := (2z - 1) fPS[z] - 2z(2 + Log[z])

(* arxiv:hep-ph/0609168, Eq.(72) *)
fsf[z_] := z/2 (2 + Log[z] - fPS[z])

(* arxiv:1502.04199, Eq.(25) *)
F1[0] := 0

F1[w_] := Module[{x}, w/2 NIntegrate[(2x(1-x)-1)/(w-x(1-x)) Log[w/(x(1-x))], {x,0,1}]]

(* arxiv:1502.04199, Eq.(26) *)
F1t[0] := 0

F1t[w_] := Module[{x}, w/2 NIntegrate[1/(w-x(1-x)) Log[w/(x(1-x))], {x,0,1}]]

(* arxiv:1502.04199, Eq.(27) *)
F2[w_] := Module[{x}, 1/2 NIntegrate[x(x-1)/(w-x(1-x)) Log[w/(x(1-x))], {x,0,1}]]

(* arxiv:1502.04199, Eq.(28) *)
F3[w_] := Module[{x}, 1/2 NIntegrate[(x w(3x(4x-1)+10)-x(1-x))/(w-x(1-x)) Log[w/(x(1-x))], {x,0,1}]]

(* arxiv:1502.04199, Eq.(29) *)
G[wa_, wb_, x_] := Log[(wa*x+wb*(1-x))/(x(1-x))]/(x(1-x)-wa*x-wb*(1-x))

(* Integral of x^n G[wa,wb,x] over {x,0,1} *)
Gn[wa_, wb_, n_] := NIntegrate[x^n G[wa, wb, x], {x,0,1}]
