expandAmu = {
    CW2 -> CW^2,
    SW2 -> 1 - CW2,
    xA -> MA0^2 / MZ^2,
    xHp -> MHp^2 / MZ^2,
    xhSM -> MH^2 / MZ^2,
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
T0[u_, w_] :=
    9/CW2^2*(u - w)*(CW2*(u - w)*(u + 2*w) - (u - w)^3 + CW2^2*w)/(CW2^2 + (u - w)^2 - 2*CW2*(u + w))*Phi[u, w, CW2]


(* Eq.(73), arxiv:1607.06292 *)
T1[u_, w_] :=
    9/CW2^2*(u - w)*(CW2*w - (u - w)^2)*PolyLog[2, 1 - u/w]


(* Eq.(74), arxiv:1607.06292 *)
T2[u_, w_, sgn_] :=
    Log[u]*(
      + (6*u^2 + CW2*(u - xHp) + 2*CW2^2*(u - xHp))/(2*(u - w))
      + f6*(u - xHp)^2*(3*CW2^2 + 3*CW2*(u - xHp) + (u - xHp)^2)/(CW2*(u - w))
      + sgn*f7*3*u^2*(u - xHp)/((xA - xH)*(u - w))
      - f8*3*u*(u - xHp)^2/(2*(u - w))
      - f9*3*u*(u - xHp)/(2*(u - w))
    )

T2p[args__] := T2[args, +1]

T2m[args__] := T2[args, -1]


(* Eq.(75), arxiv:1607.06292 *)
T4[u_, w_] :=
    (u - w)*Log[u]/4*f5*(xA*(3 + 2*xH) - xA^2 + 3*xH - xH^2 - 3)


(* Eq.(76), arxiv:1607.06292 *)
T5[u_, w_] :=
    Log[u]*(
      3/2*u + f6/CW2*((u - w)^3 + 3*CW2*(u - w)^2 + 3*CW2^2*(u - w))
      - 3/2*f8*u*(u - w) - CW2/2 - CW2^2)


(* Eq.(77), arxiv:1607.06292 *)
T6[u_, w_] :=
    9/2*((u - w)*(u^2 - 2*u*w + w*(w - CW2))/CW2^2*Log[u/w]*Log[w/CW2]
         + Log[CW2]/CW2*(2*u^2 + u*(CW2 - 4*w) - w*(CW2 - 2*w)))


(* Eq.(78), arxiv:1607.06292 *)
T7[u_, w_] :=
    Module[{ s1 = u + w - 1 + Sqrt[1 + (u - w)^2 - 2*(u + w)] },
           f5*(2*(u + w) - (u - w)^2 - 1)*Log[s1/(2*Sqrt[u*w])]*(u + w - 1 - 4*u*w/s1)
    ]


(* Eq.(79), arxiv:1607.06292 *)
T8[u_, w_] :=
    Module[{ s2 = u + w - CW2 + Sqrt[(u + w - CW2)^2 - 4*u*w] },
           2*f6*(4*u*w - (u + w - CW2)^2)*Log[s2/(2*Sqrt[u*w])]*((u + w)/CW2 - 4*u*w/(CW2*s2) - 1)
    ]


(* Eq (71), arxiv:1607:06292 *)
amu2LBNonYuk = (AL/(24*Pi*CW2*(1 - CW2)) MM/MZ)^2 (
    + (xA - xH)/(xA - xHp)*T2p[xA, xH]
    + T2m[xH, xHp]
    + (xA - xH)/(xA - xHp)*T4[xA, xHp]
    + T4[xH, xA]
    + T5[xHp, xH]
    + T5[xHp, xA]
    + T2p[xHp, xH]
    + T2p[xHp, xA]
    + T6[xA, xHp]
    + T6[xH, xHp]
    + T7[xA, xH]
    + T7[xHp, xHp]*(1 - 2*CW2)^2
    + T8[xA, xHp]
    + T8[xH, xHp]
    - 16/3*CW2*(1 - CW2)*(1 + 8*CW2 - 8*CW2^2)
    + 8*CW2^2*(1 - CW2)^2/(5*xHp)
    + f2*xHp
    - f3*xHp^2
    + f1*(xA^2 + xH^2)
    + f3*xHp*(xA + xH)
    + f4*(xA + xH)
    - f5*xA*xH
    + T1[xA, xHp]
    + T1[xH, xHp]
    + T0[xA, xHp]
    + T0[xH, xHp]
) //. expandAmu


(* Eq.(99), arxiv:1607.06292 *)
b[u_, w_] := AL*Pi/(CW2*(-1 + CW2))*(u + 2*w)


(* Eq.(100), arxiv:1607.06292 *)
Fm0[u_, w_] :=
    1/(AL*Pi) * CW2*(-1 + CW2)/(u + 2*w) * YF1[u, w]


(* Eq.(101), arxiv:1607.06292 *)
Fmp[u_, w_] :=
    (-9*(-1 + CW2))/(AL*Pi) * (T9[u, w]/2 + T10[u, w])


(* Eq.(121), arxiv:1607.06292 *)
YFZ[u_] :=
    Module[{
        z1 = 3*(17 - 48*CW2 + 32*CW2^2), (* Eq.(122) *)
        z2 = 5 - 12*CW2 + 8*CW2^2,       (* Eq.(123) *)
        z3 = 3*(1 - 3*CW2 + 2*CW2^2)     (* Eq.(124) *)
        },
        (
            + z1*u*PolyLog[2, 1 - u]
            + z2/(2*u^2)*(6*(-4 + u)*u + Pi^2*(4 + 3*u) + 6*u*(4 + u)*Log[u]
                          - 6*(4 + 3*u)*PolyLog[2, 1 - u] + 6*u*(2 + u)*Phi[u, 1, 1])
            + z3*u*(6 + Pi^2*(-4 + u)*u + 3*Log[u]*(4 + (-4 + u)*u*Log[u])
                    + 12*(-4 + u)*u*PolyLog[2, 1 - u] + 6*(-2 + u)*Phi[u, 1, 1])
        )
    ]


(* Eq.(125), arxiv:1607.06292 *)
YFW[u_] := (
    - 57/2*CW2 - 4*CW2^3*Pi^2/u^2 + 3*CW2^2*(32 - 3*Pi^2)/(4*u)
    + 3*(16*CW2^3 + 9*CW2^2*u + 12*CW2*u^2 - 19*u^3)*PolyLog[2, 1 - u/CW2]/(2*u^2)
    + 3*CW2*(16*CW2 + 19*u)*(Log[CW2/u])/(2*u)
    + 3*(4*CW2^2 - 50*CW2*u + 19*u^2)*Phi[u, CW2, CW2]/(2*(4*CW2-u)*u)
)

(* Eq.(102), arxiv:1607.06292 *)
YF1[u_, w_] := (
    - 72*CW2*(-1 + CW2)*(u + 2*w)/u - 36*CW2*(-1 + CW2)*(u + 2*w)/u*Log[w]
    + 9*(-8*CW2^2 - 3*u + 2*CW2*(4 + u))*(u + 2*w)/(2*(u-1)*u)*Log[u]
    - 9*(3 - 10*CW2 + 8*CW2^2)*w*(u + 2*w)/((4*w-1)*(u-1))*Phi[w, w, 1]
    + 9*(8*CW2^2 + 3*u - 2*CW2*(4 + u))*w*(u + 2*w)/((4*w-u)*(u-1)*u^2)*Phi[u, w, w]
)


(* Eq.(105), arxiv:1607.06292 *)
YF2[u_] :=
    Module[{
        f0 = 3/4*CW2^2*(-640 + 576*CW2 + 7*Pi^2),       (* Eq.(106) *)
        f1 = 96*CW2^3*(11 - 53*CW2 + 36*CW2^2),         (* Eq.(107) *)
        f2 = -3/4*CW2*(-66*CW2 - 48*CW2^2 + 672*CW2^3), (* Eq.(108) *)
        f3 = -3/4*CW2*(109 - 430*CW2 + 120*CW2^2),      (* Eq.(109) *)
        f4 = 96*CW2^3*(-11 + 9*CW2),                    (* Eq.(110) *)
        f5 = 45/2*CW2^2 + 192*CW2^3,                    (* Eq.(111) *)
        f6 = 3/4*CW2*(157 + 90*CW2),                    (* Eq.(112) *)
        f7 = -3/4*(18 + 61*CW2),                        (* Eq.(113) *)
        f8 = -7 + 61*CW2 - 162*CW2^2 + 96*CW2^3,        (* Eq.(114) *)
        f9 = 1 - 5*CW2 + 10*CW2^2,                      (* Eq.(115) *)
        f10 = -1728*CW2^4*(-1 + CW2),                   (* Eq.(116) *)
        f11 = 3*CW2^3*(-899 + 768*CW2),                 (* Eq.(117) *)
        f12 = 387*CW2^2 - 363*CW2^3,                    (* Eq.(118) *)
        f13 = 9/2*CW2*(57 + 106*CW2),                   (* Eq.(119) *)
        f14 = -15/2*(7 + 45*CW2)                        (* Eq.(120) *)
        },
        (
            + YFW[u]
            + YFZ[u]
            + 8*CW2^3*Pi^2/u^2 + f0/u + 393/8*CW2
            + (f1/u + f2 + f3*u)*Log[CW2]/((4*CW2-1)*(4*CW2-u))
            + (f4/u + f5 + f6*u + f7*u^2)*Log[u]/((u-1)*(4*CW2-u))
            - 3/2*(32*CW2^3/u^2 + 21*CW2^2/u + 15*CW2 - 35*u)*PolyLog[2, 1 - u/CW2]
            + (f8 + f9*u)*9*CW2*(-3 + 4*CW2)/2*Phi[CW2, CW2, 1]/((4*CW2-1)^2*(u-1))
            + (f10/u^2 + f11/u + f12 + f13*u + f14*u^2 + 105/2*u^3)*Phi[u, CW2, CW2]/
              ((4*CW2-u)^2*(u-1))
        )
    ]


(* Eq.(126), arxiv:1607.06292 *)
YF3[u_, w_] :=
    Module[{
        (* Eq.(127) *)
        a1 = -9*CW2*u^3 + 9*CW2*u^2*(3*CW2+w) + 27*CW2^2*u*(w-CW2) + 9*(CW2^4 - 4*CW2^3*w + 3*CW2^2*w^2),
        (* Eq.(128) *)
        a2 = 9*CW2^2*w/2 - 9*u^2*(5*CW2 + w) + u*(36*CW2^2 + 153*CW2*w/4) + 9*u^3,
        (* Eq.(129) *)
        a3 = 9*CW2*u^2 - 9/2*CW2*u*(4*CW2 + w),
        (* Eq.(130) *)
        a4 = -9/2*u^2*w*(2*CW2^2 + 9*CW2*w + 2*w^2) + 9/8*u*w*(32*CW2^3 + 13*CW2^2*w + 35*CW2*w^2) + 9*u^3*w^2,
        (* Eq.(131) *)
        a5 = -9*u^3*(CW2 + w) - 9*u*(3*CW2^3 + 2*CW2*w^2) + 9*u^2*(3*CW2^2 + 4*CW2*w + w^2) + 9/2*CW2^2*(2*CW2^2 - 6*CW2*w + w^2),
        (* Eq.(132) *)
        a6 = -9*u^4*(9*CW2 + w) + u*(81*CW2^3*w - 225*CW2^4) + 9*CW2^4*(w - CW2) - 9/2*u^2*(3*CW2^3 + 37*CW2^2*w) + u^3*(198*CW2^2 + 72*CW2*w) + 9*u^5,
        (* Eq.(133) *)
        a7 = -9*CW2*u^4 + 18*CW2*u^3*(2*CW2 + w) + 36*u*(CW2^4 - 2*CW2^3*w) - 9*CW2*u^2*(6*CW2^2 - CW2*w + w^2) - 9*CW2*(CW2 - 3*w)*(CW2^3 - 2*CW2^2*w + CW2*w^2)
        },
        (
            + 9*u*(2*CW2 - u + w)/w
            + (a1*(Log[u] - Log[CW2]) + 9*CW2^2*(CW2^2 - 4*CW2*w + 3*w^2)*Log[CW2])*(Log[w] - Log[CW2])/(2*w^2*(CW2-w))
            + a2*Log[u]/(w*(4*CW2-u))
            + a3*Log[w]/(w*(CW2-w))
            + a4*Log[CW2]/(w^2*(4*CW2-u)*(CW2-w))
            + a5/(CW2*w^2)*PolyLog[2, 1 - u/CW2]
            + a6/(u*CW2*(4*CW2-u)^2*(CW2-w))*Phi[u, CW2, CW2]
            + a7/(w^2*(CW2-w)*(CW2^2-2*CW2*(u+w)+(u-w)^2))*Phi[u, w, CW2]
        )
    ]


(* Eq.(103), arxiv:1607.06292 *)
T9[u_, w_] := (
    - 2*(CW2^2*w + CW2*(u^2 + u*w - 2*w^2) - (u-w)^3)*Phi[u, w, CW2]/
      ((CW2 - w)*(CW2^2 - 2*CW2*(u+w) + (u-w)^2))
    + 2*CW2^2*(u^2 - 4*u*w + 2*w^2)*Phi[u, w, w]/(w^2*(w-CW2)*(u-4*w))
    - 2*(CW2*u*(u-2*w) + w*(u-w)^2)*PolyLog[2, 1 - u/w]/w^2
)


(* Eq.(104), arxiv:1607.06292 *)
T10[u_, w_] := (
    (u^2 - CW2*w - 2*u*w + w^2)/(2*(CW2-w))*Log[w/u]*Log[w/CW2]
    + CW2*(CW2 + 2*u - 2*w)/(2*(CW2-w))*Log[w/CW2]
    + CW2*u/w*Log[w/u]
    + CW2/w*(w-u)
)


(* Eq (52), arxiv:1607:06292 *)
amu2LBYuk =
    Module[{
        (* Eq.(91), arxiv:1607.06292 *)
        a000 = b[xhSM, xHp]*Fm0[xhSM, xHp],
        (* Eq.(92), arxiv:1607.06292 *)
        a0z0 = -b[xH, 0]*(Fm0[xH, xHp] + Fmp[xH, xHp]),
        (* Eq.(93), arxiv:1607.06292 *)
        a500 = Fm0[xhSM, xHp],
        (* Eq.(94), arxiv:1607.06292 *)
        a5z0 = -1/2 (Fm0[xH, xHp] + Fmp[xH, xHp]),
        (* Eq.(95), arxiv:1607.06292 *)
        a001 = (+ b[xH, 0]*Fm0[xH, xHp]
                - b[xhSM, 0]*Fm0[xhSM, xHp]),
        (* Eq.(96), arxiv:1607.06292 *)
        a0z1 = (
            - (
                + b[xH, xHp]*(
                    + Fm0[xH, xHp]
                    + Fmp[xH, xHp]
                  )
                - YF3[xH, xHp]
                (* SM contributions with opposite sign: *)
                - b[xhSM, xHp]*(
                    + Fm0[xhSM, xHp]
                    + Fmp[xhSM, xHp]
                  )
                + YF3[xhSM, xHp]
            ) + YF2[xH] ),
        (* Eq.(97), arxiv:1607.06292 *)
        a501 = (+ Fm0[xH, xHp]/2
                - Fm0[xhSM, xHp]/2),
        (* Eq.(98), arxiv:1607.06292 *)
        a5z1 = (- Fm0[xH, xHp]
                - Fmp[xH, xHp]
                + Fm0[xhSM, xHp]
                + Fmp[xhSM, xHp])
        },
        (AL/(24*Pi*CW2*(1 - CW2)) MM/MZ)^2 (
            + a000
            + a0z0*(TB - 1/TB)*ZetaL
            + a500*Lambda5
            + a5z0*(TB - 1/TB)*Lambda5*ZetaL
            + (
                + a001*(TB - 1/TB)
                + a0z1*ZetaL
                + a501*(TB - 1/TB)*Lambda5
                + a5z1*Lambda5*ZetaL
            )*aeps
        ) //. expandAmu
    ]


Null
