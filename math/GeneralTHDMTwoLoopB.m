expandAmu = {
    CW2 -> CW^2,
    SW2 -> 1 - CW2,
    xA -> MA0^2 / MZ^2,
    xHp -> MHp^2 / MZ^2,
    xH -> MHH^2 / MZ^2,
    f1 -> 7/2 - 25/(2*CW2) + 4*CW2 - 4*CW2^2,
    f2 -> 2*(17 - 24*CW2 + 56*CW2^2 - 128*CW2^3 + 64*CW2^4),
    f3 -> (25 - 32*CW2 + 4*CW2^2)/(CW2*SW2),
    f4 -> 13/2 - 15*CW2 + 10*CW2^2,
    f5 -> CW2*(5 - 16*CW2 + 8*CW2^2)/SW2,
    f6 -> (7 - 14*CW2 + 4*CW2^2)/(4*CW2*SW2),
    f7 -> 1 - 6*CW2 + 4*CW2^2,
    f8 -> (13 - 20*CW2 + 4*CW2^2)/(CW2*SW2),
    f9 -> 7 - 12*CW2 + 8*CW2^2
}


(* Eq.(72), arxiv:1607.06292 *)
T0[u_, w_, cw2_] :=
    9/cw2^2*(u - w)*(cw2*(u - w)*(u + 2*w) - (u - w)^3 + cw2^2*w)/(cw2^2 + (u - w)^2 - 2*cw2*(u + w))*Phi[u, w, cw2]


(* Eq.(73), arxiv:1607.06292 *)
T1[u_, w_, cw2_] :=
    9/cw2^2*(u - w)*(cw2*w - (u - w)^2)*PolyLog[2, 1 - u/w]


(* Eq.(74), arxiv:1607.06292 *)
T2[u_, w_, cw2_, xH_, xA_, xHp_, sgn_] :=
    Log[u]*(
      + (6*u^2 + cw2*(u - xHp) + 2*cw2^2*(u - xHp))/(2*(u - w))
      + f6*(u - xHp)^2*(3*cw2^2 + 3*cw2*(u - xHp) + (u - xHp)^2)/(cw2*(u - w))
      + sgn*f7*3*u^2*(u - xHp)/((xA - xH)*(u - w))
      - f8*3*u*(u - xHp)^2/(2*(u - w))
      - f9*3*u*(u - xHp)/(2*(u - w))
    )


T2p[u_, w_, cw2_, xH_, xA_, xHp_] :=
    T2[u, w, cw2, xH, xA, xHp, +1]


T2m[u_, w_, cw2_, xH_, xA_, xHp_] :=
    T2[u, w, cw2, xH, xA, xHp, -1]


(* Eq.(75), arxiv:1607.06292 *)
T4[u_, w_, cw2_, xH_, xA_] :=
    (u - w)*Log[u]/4*f5*(xA*(3 + 2*xH) - xA^2 + 3*xH - xH^2 - 3)


(* Eq.(76), arxiv:1607.06292 *)
T5[u_, w_, cw2_] :=
    Log[u]*(
      3/2*u + f6/cw2*((u - w)^3 + 3*cw2*(u - w)^2 + 3*cw2^2*(u - w))
      - 3/2*f8*u*(u - w) - cw2/2 - cw2^2)


(* Eq.(77), arxiv:1607.06292 *)
T6[u_, w_, cw2_] :=
    9/2*((u - w)*(u^2 - 2*u*w + w*(w - cw2))/cw2^2*Log[u/w]*Log[w/cw2]
         + Log[cw2]/cw2*(2*u^2 + u*(cw2 - 4*w) - w*(cw2 - 2*w)))


(* Eq.(78), arxiv:1607.06292 *)
T7[u_, w_, cw2_] :=
    Module[{ s1 = u + w - 1 + Sqrt[1 + (u - w)^2 - 2*(u + w)] },
           f5*(2*(u + w) - (u - w)^2 - 1)*Log[s1/(2*Sqrt[u*w])]*(u + w - 1 - 4*u*w/s1)
    ]


(* Eq.(79), arxiv:1607.06292 *)
T8[u_, w_, cw2_] :=
    Module[{ s2 = u + w - cw2 + Sqrt[(u + w - cw2)^2 - 4*u*w] },
           2*f6*(4*u*w - (u + w - cw2)^2)*Log[s2/(2*Sqrt[u*w])]*((u + w)/cw2 - 4*u*w/(cw2*s2) - 1)
    ]


(* Eq (71), arxiv:1607:06292 *)
amu2LBNonYuk = (AL/(24*Pi*CW2*(1 - CW2)) MM/MZ)^2 (
    + (xA - xH)/(xA - xHp)*T2p[xA, xH, CW2, xH, xA, xHp]
    + T2m[xH, xHp, CW2, xH, xA, xHp]
    + (xA - xH)/(xA - xHp)*T4[xA, xHp, CW2, xH, xA]
    + T4[xH, xA, CW2, xH, xA]
    + T5[xHp, xH, CW2]
    + T5[xHp, xA, CW2]
    + T2p[xHp, xH, CW2, xH, xA, xHp]
    + T2p[xHp, xA, CW2, xH, xA, xHp]
    + T6[xA, xHp, CW2]
    + T6[xH, xHp, CW2]
    + T7[xA, xH, CW2]
    + T7[xHp, xHp, CW2]*(1 - 2*CW2)^2
    + T8[xA, xHp, CW2]
    + T8[xH, xHp, CW2]
    - 16/3*CW2*(1 - CW2)*(1 + 8*CW2 - 8*CW2^2)
    + 8*CW2^2*(1 - CW2)^2/(5*xHp)
    + f2*xHp
    - f3*xHp^2
    + f1*(xA^2 + xH^2)
    + f3*xHp*(xA + xH)
    + f4*(xA + xH)
    - f5*xA*xH
    + T1[xA, xHp, CW2]
    + T1[xH, xHp, CW2]
    + T0[xA, xHp, CW2]
    + T0[xH, xHp, CW2]
) //. expandAmu

Null
