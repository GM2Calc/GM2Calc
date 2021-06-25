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


T2p[args__] := T2[args, +1]


T2m[args__] := T2[args, -1]


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


(* Eq.(99), arxiv:1607.06292 *)
b[u_, w_, al_, cw2_] := al*Pi/(cw2*(-1 + cw2))*(u + 2*w)


(* Eq.(100), arxiv:1607.06292 *)
Fm0[u_, w_, al_, cw2_] :=
    1/(al*Pi) * cw2*(-1 + cw2)/(u + 2*w) * YF1[u, w, cw2]


(* Eq.(101), arxiv:1607.06292 *)
Fmp[u_, w_, al_, cw2_] :=
    (-9*(-1 + cw2))/(al*Pi) * (T9[u, w, cw2]/2 + T10[u, w, cw2])


(* Eq.(121), arxiv:1607.06292 *)
YFZ[u_, cw2_] :=
    Module[{
        z1 = 3*(17 - 48*cw2 + 32*cw2^2), (* Eq.(122) *)
        z2 = 5 - 12*cw2 + 8*cw2^2,       (* Eq.(123) *)
        z3 = 3*(1 - 3*cw2 + 2*cw2^2)     (* Eq.(124) *)
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
YFW[u_, cw2_] := (
    - 57/2*cw2 - 4*cw2^3*Pi^2/u^2 + 3*cw2^2*(32 - 3*Pi^2)/(4*u)
    + 3*(16*cw2^3 + 9*cw2^2*u + 12*cw2*u^2 - 19*u^3)*PolyLog[2, 1 - u/cw2]/(2*u^2)
    + 3*cw2*(16*cw2 + 19*u)*(Log[cw2/u])/(2*u)
    + 3*(4*cw2^2 - 50*cw2*u + 19*u^2)*Phi[u, cw2, cw2]/(2*(4*cw2-u)*u)
)

(* Eq.(102), arxiv:1607.06292 *)
YF1[u_, w_, cw2_] := (
    - 72*cw2*(-1 + cw2)*(u + 2*w)/u - 36*cw2*(-1 + cw2)*(u + 2*w)/u*Log[w]
    + 9*(-8*cw2^2 - 3*u + 2*cw2*(4 + u))*(u + 2*w)/(2*(u-1)*u)*Log[u]
    - 9*(3 - 10*cw2 + 8*cw2^2)*w*(u + 2*w)/((4*w-1)*(u-1))*Phi[w, w, 1]
    + 9*(8*cw2^2 + 3*u - 2*cw2*(4 + u))*w*(u + 2*w)/((4*w-u)*(u-1)*u^2)*Phi[u, w, w]
)


(* Eq.(105), arxiv:1607.06292 *)
YF2[u_, cw2_] :=
    Module[{
        f0 = 3/4*cw2^2*(-640 + 576*cw2 + 7*Pi^2),       (* Eq.(106) *)
        f1 = 96*cw2^3*(11 - 53*cw2 + 36*cw2^2),         (* Eq.(107) *)
        f2 = -3/4*cw2*(-66*cw2 - 48*cw2^2 + 672*cw2^3), (* Eq.(108) *)
        f3 = -3/4*cw2*(109 - 430*cw2 + 120*cw2^2),      (* Eq.(109) *)
        f4 = 96*cw2^3*(-11 + 9*cw2),                    (* Eq.(110) *)
        f5 = 45/2*cw2^2 + 192*cw2^3,                    (* Eq.(111) *)
        f6 = 3/4*cw2*(157 + 90*cw2),                    (* Eq.(112) *)
        f7 = -3/4*(18 + 61*cw2),                        (* Eq.(113) *)
        f8 = -7 + 61*cw2 - 162*cw2^2 + 96*cw2^3,        (* Eq.(114) *)
        f9 = 1 - 5*cw2 + 10*cw2^2,                      (* Eq.(115) *)
        f10 = -1728*cw2^4*(-1 + cw2),                   (* Eq.(116) *)
        f11 = 3*cw2^3*(-899 + 768*cw2),                 (* Eq.(117) *)
        f12 = 387*cw2^2 - 363*cw2^3,                    (* Eq.(118) *)
        f13 = 9/2*cw2*(57 + 106*cw2),                   (* Eq.(119) *)
        f14 = -15/2*(7 + 45*cw2)                        (* Eq.(120) *)
        },
        (
            + YFW[u, cw2]
            + YFZ[u, cw2]
            + 8*cw2^3*Pi^2/u^2 + f0/u + 393/8*cw2
            + (f1/u + f2 + f3*u)*Log[cw2]/((4*cw2-1)*(4*cw2-u))
            + (f4/u + f5 + f6*u + f7*u^2)*Log[u]/((u-1)*(4*cw2-u))
            - 3/2*(32*cw2^3/u^2 + 21*cw2^2/u + 15*cw2 - 35*u)*PolyLog[2, 1 - u/cw2]
            + (f8 + f9*u)*9*cw2*(-3 + 4*cw2)/2*Phi[cw2, cw2, 1]/((4*cw2-1)^2*(u-1))
            + (f10/u^2 + f11/u + f12 + f13*u + f14*u^2 + 105/2*u^3)*Phi[u, cw2, cw2]/
              ((4*cw2-u)^2*(u-1))
        )
    ]


(* Eq.(126), arxiv:1607.06292 *)
YF3[u_, w_, cw2_] :=
    Module[{
        (* Eq.(127) *)
        a1 = -9*cw2*u^3 + 9*cw2*u^2*(3*cw2+w) + 27*cw2^2*u*(w-cw2) + 9*(cw2^4 - 4*cw2^3*w + 3*cw2^2*w^2),
        (* Eq.(128) *)
        a2 = 9*cw2^2*w/2 - 9*u^2*(5*cw2 + w) + u*(36*cw2^2 + 153*cw2*w/4) + 9*u^3,
        (* Eq.(129) *)
        a3 = 9*cw2*u^2 - 9/2*cw2*u*(4*cw2 + w),
        (* Eq.(130) *)
        a4 = -9/2*u^2*w*(2*cw2^2 + 9*cw2*w + 2*w^2) + 9/8*u*w*(32*cw2^3 + 13*cw2^2*w + 35*cw2*w^2) + 9*u^3*w^2,
        (* Eq.(131) *)
        a5 = -9*u^3*(cw2 + w) - 9*u*(3*cw2^3 + 2*cw2*w^2) + 9*u^2*(3*cw2^2 + 4*cw2*w + w^2) + 9/2*cw2^2*(2*cw2^2 - 6*cw2*w + w^2),
        (* Eq.(132) *)
        a6 = -9*u^4*(9*cw2 + w) + u*(81*cw2^3*w - 225*cw2^4) + 9*cw2^4*(w - cw2) - 9/2*u^2*(3*cw2^3 + 37*cw2^2*w) + u^3*(198*cw2^2 + 72*cw2*w) + 9*u^5,
        (* Eq.(133) *)
        a7 = -9*cw2*u^4 + 18*cw2*u^3*(2*cw2 + w) + 36*u*(cw2^4 - 2*cw2^3*w) - 9*cw2*u^2*(6*cw2^2 - cw2*w + w^2) - 9*cw2*(cw2 - 3*w)*(cw2^3 - 2*cw2^2*w + cw2*w^2)
        },
        (
            + 9*u*(2*cw2 - u + w)/w
            + (a1*(Log[u] - Log[cw2]) + 9*cw2^2*(cw2^2 - 4*cw2*w + 3*w^2)*Log[cw2])*(Log[w] - Log[cw2])/(2*w^2*(cw2-w))
            + a2*Log[u]/(w*(4*cw2-u))
            + a3*Log[w]/(w*(cw2-w))
            + a4*Log[cw2]/(w^2*(4*cw2-u)*(cw2-w))
            + a5/(cw2*w^2)*PolyLog[2, 1 - u/cw2]
            + a6/(u*cw2*(4*cw2-u)^2*(cw2-w))*Phi[u, cw2, cw2]
            + a7/(w^2*(cw2-w)*(cw2^2-2*cw2*(u+w)+(u-w)^2))*Phi[u, w, cw2]
        )
    ]


(* Eq.(103), arxiv:1607.06292 *)
T9[u_, w_, cw2_] := (
    - 2*(cw2^2*w + cw2*(u^2 + u*w - 2*w^2) - (u-w)^3)*Phi[u, w, cw2]/
      ((cw2 - w)*(cw2^2 - 2*cw2*(u+w) + (u-w)^2))
    + 2*cw2^2*(u^2 - 4*u*w + 2*w^2)*Phi[u, w, w]/(w^2*(w-cw2)*(u-4*w))
    - 2*(cw2*u*(u-2*w) + w*(u-w)^2)*PolyLog[2, 1 - u/w]/w^2
)


(* Eq.(104), arxiv:1607.06292 *)
T10[u_, w_, cw2_] := (
    (u^2 - cw2*w - 2*u*w + w^2)/(2*(cw2-w))*Log[w/u]*Log[w/cw2]
    + cw2*(cw2 + 2*u - 2*w)/(2*(cw2-w))*Log[w/cw2]
    + cw2*u/w*Log[w/u]
    + cw2/w*(w-u)
)


(* Eq (52), arxiv:1607:06292 *)
amu2LBYuk =
    Module[{
        (* Eq.(91), arxiv:1607.06292 *)
        a000 = b[xhSM, xHp, AL, CW2]*Fm0[xhSM, xHp, AL, CW2],
        (* Eq.(92), arxiv:1607.06292 *)
        a0z0 = -b[xH, 0, AL, CW2]*(Fm0[xH, xHp, AL, CW2] + Fmp[xH, xHp, AL, CW2]),
        (* Eq.(93), arxiv:1607.06292 *)
        a500 = Fm0[xhSM, xHp, AL, CW2],
        (* Eq.(94), arxiv:1607.06292 *)
        a5z0 = -1/2 (Fm0[xH, xHp, AL, CW2] + Fmp[xH, xHp, AL, CW2]),
        (* Eq.(95), arxiv:1607.06292 *)
        a001 = (+ b[xH, 0, AL, CW2]*Fm0[xH, xHp, AL, CW2]
                - b[xhSM, 0, AL, CW2]*Fm0[xhSM, xHp, AL, CW2]),
        (* Eq.(96), arxiv:1607.06292 *)
        a0z1 = (
            - (
                + b[xH, xHp, AL, CW2]*(
                    + Fm0[xH, xHp, AL, CW2]
                    + Fmp[xH, xHp, AL, CW2]
                  )
                - YF3[xH, xHp, CW2]
                (* SM contributions with opposite sign: *)
                - b[xhSM, xHp, AL, CW2]*(
                    + Fm0[xhSM, xHp, AL, CW2]
                    + Fmp[xhSM, xHp, AL, CW2]
                  )
                + YF3[xhSM, xHp, CW2]
            ) + YF2[xH, CW2] ),
        (* Eq.(97), arxiv:1607.06292 *)
        a501 = (+ Fm0[xH, xHp, AL, CW2]/2
                - Fm0[xhSM, xHp, AL, CW2]/2),
        (* Eq.(98), arxiv:1607.06292 *)
        a5z1 = (- Fm0[xH, xHp, AL, CW2]
                - Fmp[xH, xHp, AL, CW2]
                + Fm0[xhSM, xHp, AL, CW2]
                + Fmp[xhSM, xHp, AL, CW2])
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
