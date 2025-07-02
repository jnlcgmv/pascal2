program TestSun200;
uses
  SysUtils;
var
  T: Double;   // Tiempo en siglos julianos desde J2000
  L, B, R, sun_pos_1, sun_pos_2, sun_pos_3, P2: Double; // Resultados: ecliptic coordinates L,B,R
(*-----------------------------------------------------------------------*)
(* SUN200: ecliptic coordinates L,B,R (in deg and AU) of the             *)
(*         Sun referred to the mean equinox of date                      *)
(*         (T: time in Julian centuries since J2000)                     *)
(*         (   = (JED-2451545.0)/36525             )                     *)
(*-----------------------------------------------------------------------*)
PROCEDURE SUN200(T:DOUBLE;VAR L,B,R:DOUBLE);
  VAR C3,S3:          ARRAY [-1..7] OF DOUBLE;
      C,S:            ARRAY [-8..0] OF DOUBLE;
      M2,M3,M4,M5,M6: DOUBLE;
      D,A,UU:         DOUBLE;
      U,V,DL,DR,DB:   DOUBLE;
      I:              INTEGER;
  FUNCTION FRAC(X:DOUBLE):DOUBLE;
    (* for some compilers TRUNC has to be replaced by LONG_TRUNC *)
    (* or INT (Turbo-Pascal) if the routine is used with T<-24   *)
    BEGIN  X:=X-TRUNC(X); IF (X<0) THEN X:=X+1.0; FRAC:=X  END;
  PROCEDURE ADDTHE(C1,S1,C2,S2:DOUBLE; VAR C,S:DOUBLE);
    BEGIN  C:=C1*C2-S1*S2; S:=S1*C2+C1*S2; END;
  PROCEDURE TERM(I1,I,IT:INTEGER;DLC,DLS,DRC,DRS,DBC,DBS:DOUBLE);
    BEGIN
      WriteLn(Format('C3[I1] = %.18f', [C3[I1]])); 
      WriteLn(Format('S3[I1] = %.18f', [S3[I1]])); 
      WriteLn(Format('C[I] = %.18f', [C[I]])); 
      WriteLn(Format('S[I] = %.18f', [S[I]])); 
      IF IT=0 THEN ADDTHE(C3[I1],S3[I1],C[I],S[I],U,V)
              ELSE BEGIN U:=U*T; V:=V*T END;
      DL:=DL+DLC*U+DLS*V; DR:=DR+DRC*U+DRS*V; DB:=DB+DBC*U+DBS*V;
      WriteLn(Format('U = %.18f', [U]));
      WriteLn(Format('V = %.18f', [V]));
      WriteLn(Format('DLC = %.18f', [DLC]));
      WriteLn(Format('DLS = %.18f', [DLS]));
    END;
  PROCEDURE PERTVEN;  (* Keplerian terms and perturbations by Venus *)
    VAR I: INTEGER;
    BEGIN
      C[0]:=1.0; S[0]:=0.0; C[-1]:=COS(M2); S[-1]:=-SIN(M2);
      FOR I:=-1 DOWNTO -5 DO ADDTHE(C[I],S[I],C[-1],S[-1],C[I-1],S[I-1]);
      TERM(1, 0,0,-0.22,6892.76,-16707.37, -0.54, 0.00, 0.00);
      WriteLn(Format('inc_DR = %.18f', [DR]));
      WriteLn(Format('inc_DL = %.18f', [DL]));
      WriteLn(Format('inc_DB = %.18f', [DB]));
      TERM(1, 0,1,-0.06, -17.35,    42.04, -0.15, 0.00, 0.00);
      TERM(1, 0,2,-0.01,  -0.05,     0.13, -0.02, 0.00, 0.00);
      TERM(2, 0,0, 0.00,  71.98,  -139.57,  0.00, 0.00, 0.00);
      TERM(2, 0,1, 0.00,  -0.36,     0.70,  0.00, 0.00, 0.00);
      TERM(3, 0,0, 0.00,   1.04,    -1.75,  0.00, 0.00, 0.00);
      TERM(0,-1,0, 0.03,  -0.07,    -0.16, -0.07, 0.02,-0.02);
      TERM(1,-1,0, 2.35,  -4.23,    -4.75, -2.64, 0.00, 0.00);
      TERM(1,-2,0,-0.10,   0.06,     0.12,  0.20, 0.02, 0.00);
      TERM(2,-1,0,-0.06,  -0.03,     0.20, -0.01, 0.01,-0.09);
      TERM(2,-2,0,-4.70,   2.90,     8.28, 13.42, 0.01,-0.01);
      TERM(3,-2,0, 1.80,  -1.74,    -1.44, -1.57, 0.04,-0.06);
      TERM(3,-3,0,-0.67,   0.03,     0.11,  2.43, 0.01, 0.00);
      TERM(4,-2,0, 0.03,  -0.03,     0.10,  0.09, 0.01,-0.01);
      TERM(4,-3,0, 1.51,  -0.40,    -0.88, -3.36, 0.18,-0.10);
      TERM(4,-4,0,-0.19,  -0.09,    -0.38,  0.77, 0.00, 0.00);
      TERM(5,-3,0, 0.76,  -0.68,     0.30,  0.37, 0.01, 0.00);
      TERM(5,-4,0,-0.14,  -0.04,    -0.11,  0.43,-0.03, 0.00);
      TERM(5,-5,0,-0.05,  -0.07,    -0.31,  0.21, 0.00, 0.00);
      TERM(6,-4,0, 0.15,  -0.04,    -0.06, -0.21, 0.01, 0.00);
      TERM(6,-5,0,-0.03,  -0.03,    -0.09,  0.09,-0.01, 0.00);
      TERM(6,-6,0, 0.00,  -0.04,    -0.18,  0.02, 0.00, 0.00);
      TERM(7,-5,0,-0.12,  -0.03,    -0.08,  0.31,-0.02,-0.01);
      WriteLn(Format('DR = %.15f', [DR]));
      WriteLn(Format('DL = %.15f', [DL]));
      WriteLn(Format('DB = %.15f', [DB]));
    END;
  PROCEDURE PERTMAR;  (* perturbations by Mars *)
    VAR I: INTEGER;
    BEGIN
      C[-1]:=COS(M4); S[-1]:=-SIN(M4);
      FOR I:=-1 DOWNTO -7 DO ADDTHE(C[I],S[I],C[-1],S[-1],C[I-1],S[I-1]);
      TERM(1,-1,0,-0.22,   0.17,    -0.21, -0.27, 0.00, 0.00);
      TERM(1,-2,0,-1.66,   0.62,     0.16,  0.28, 0.00, 0.00);
      TERM(2,-2,0, 1.96,   0.57,    -1.32,  4.55, 0.00, 0.01);
      TERM(2,-3,0, 0.40,   0.15,    -0.17,  0.46, 0.00, 0.00);
      TERM(2,-4,0, 0.53,   0.26,     0.09, -0.22, 0.00, 0.00);
      TERM(3,-3,0, 0.05,   0.12,    -0.35,  0.15, 0.00, 0.00);
      TERM(3,-4,0,-0.13,  -0.48,     1.06, -0.29, 0.01, 0.00);
      TERM(3,-5,0,-0.04,  -0.20,     0.20, -0.04, 0.00, 0.00);
      TERM(4,-4,0, 0.00,  -0.03,     0.10,  0.04, 0.00, 0.00);
      TERM(4,-5,0, 0.05,  -0.07,     0.20,  0.14, 0.00, 0.00);
      TERM(4,-6,0,-0.10,   0.11,    -0.23, -0.22, 0.00, 0.00);
      TERM(5,-7,0,-0.05,   0.00,     0.01, -0.14, 0.00, 0.00);
      TERM(5,-8,0, 0.05,   0.01,    -0.02,  0.10, 0.00, 0.00);
    END;
  PROCEDURE PERTJUP;  (* perturbations by Jupiter *)
    VAR I: INTEGER;
    BEGIN
      C[-1]:=COS(M5); S[-1]:=-SIN(M5);
      FOR I:=-1 DOWNTO -3 DO ADDTHE(C[I],S[I],C[-1],S[-1],C[I-1],S[I-1]);
      TERM(-1,-1,0,0.01,   0.07,     0.18, -0.02, 0.00,-0.02);
      TERM(0,-1,0,-0.31,   2.58,     0.52,  0.34, 0.02, 0.00);
      TERM(1,-1,0,-7.21,  -0.06,     0.13,-16.27, 0.00,-0.02);
      TERM(1,-2,0,-0.54,  -1.52,     3.09, -1.12, 0.01,-0.17);
      TERM(1,-3,0,-0.03,  -0.21,     0.38, -0.06, 0.00,-0.02);
      TERM(2,-1,0,-0.16,   0.05,    -0.18, -0.31, 0.01, 0.00);
      TERM(2,-2,0, 0.14,  -2.73,     9.23,  0.48, 0.00, 0.00);
      TERM(2,-3,0, 0.07,  -0.55,     1.83,  0.25, 0.01, 0.00);
      TERM(2,-4,0, 0.02,  -0.08,     0.25,  0.06, 0.00, 0.00);
      TERM(3,-2,0, 0.01,  -0.07,     0.16,  0.04, 0.00, 0.00);
      TERM(3,-3,0,-0.16,  -0.03,     0.08, -0.64, 0.00, 0.00);
      TERM(3,-4,0,-0.04,  -0.01,     0.03, -0.17, 0.00, 0.00);
    END;
  PROCEDURE PERTSAT;  (* perturbations by Saturn *)
    BEGIN
      C[-1]:=COS(M6); S[-1]:=-SIN(M6);
      ADDTHE(C[-1],S[-1],C[-1],S[-1],C[-2],S[-2]);
      TERM(0,-1,0, 0.00,   0.32,     0.01,  0.00, 0.00, 0.00);
      TERM(1,-1,0,-0.08,  -0.41,     0.97, -0.18, 0.00,-0.01);
      TERM(1,-2,0, 0.04,   0.10,    -0.23,  0.10, 0.00, 0.00);
      TERM(2,-2,0, 0.04,   0.10,    -0.35,  0.13, 0.00, 0.00);
    END;
  PROCEDURE PERTMOO;  (* difference between the Earth-Moon      *)
    BEGIN             (* barycenter and the center of the Earth *)
      DL := DL +  6.45*SIN(D) - 0.42*SIN(D-A) + 0.18*SIN(D+A)
                              + 0.17*SIN(D-M3) - 0.06*SIN(D+M3);
      DR := DR + 30.76*COS(D) - 3.06*COS(D-A)+ 0.85*COS(D+A)
                              - 0.58*COS(D+M3) + 0.57*COS(D-M3);
      DB := DB + 0.576*SIN(UU);
    END;
  BEGIN  (* SUN200 *)
    DL:=0.0; DR:=0.0; DB:=0.0;
    M2:=2.0 * PI*FRAC(0.1387306+162.5485917*T);
    M3:=2.0 * PI*FRAC(0.9931266+99.9973604*T);
    M4:=2.0 * PI*FRAC(0.0543250+ 53.1666028*T);
    M5:=2.0 * PI*FRAC(0.0551750+ 8.4293972*T);
    M6:=2.0 * PI*FRAC(0.8816500+  3.3938722*T); D :=2.0 * PI*FRAC(0.8274+1236.8531*T);
    A :=2.0 * PI*FRAC(0.3749+1325.5524*T);      UU:=2.0 * PI*FRAC(0.2591+1342.2278*T);
    WriteLn(Format('M2 = %.20f', [M2]));
    WriteLn(Format('M3 = %.20f', [M3]));
    WriteLn(Format('M4 = %.20f', [M4]));
    WriteLn(Format('M5 = %.20f', [M5]));
    WriteLn(Format('M6 = %.20f', [M6]));
    WriteLn(Format('A = %.20f', [A]));
    WriteLn(Format('D = %.20f', [D]));
    WriteLn(Format('UU = %.20f', [UU]));
    C3[0]:=1.0;     S3[0]:=0.0;
    C3[1]:=COS(M3); S3[1]:=SIN(M3);  C3[-1]:=C3[1]; S3[-1]:=-S3[1];
    FOR I:=2 TO 7 DO ADDTHE(C3[I-1],S3[I-1],C3[1],S3[1],C3[I],S3[I]);
    for i := -1 to 7 do
        writeln('C3[', i, '] = ', C3[i]:0:20);
 
    for i := -1 to 7 do
        writeln('S3[', i, '] = ', S3[i]:0:20);
    PERTVEN; PERTMAR; PERTJUP; PERTSAT; PERTMOO;
    DL:=DL + 6.40*SIN(2.0 * PI*(0.6983+0.0561*T))+1.87*SIN(2.0 * PI*(0.5764+0.4174*T))
           + 0.27*SIN(2.0 * PI*(0.4189+0.3306*T))+0.20*SIN(2.0 * PI*(0.3581+2.4814*T));
    L:= 360.0*FRAC(0.7859453 + M3/(2.0 * PI) + ((6191.2+1.1*T)*T+DL)/1296.0E3 );
    R:= 1.0001398 - 0.0000007*T  +  DR*1E-6;
    B:= DB/3600.0;
  END;   (* SUN200 *)
 

begin
  T := 0.244274469541410;  
  SUN200(T, L, B, R);
  P2 := 6.283185307;
  sun_pos_1 := R * 149597870.0 * (COS(L * 2 * PI / 360.0) * COS(B * 2 * PI / 360.0));
  sun_pos_2 := R * 149597870.0 * (COS(23.43929111 * 2 * PI / 360.0) * SIN(L * 2 * PI / 360.0) * COS(B * 2 * PI / 360.0) - SIN(23.43929111 * 2 * PI / 360.0) * SIN(B * 2 * PI / 360.0));
  sun_pos_3 := R * 149597870.0 * (SIN(23.43929111 * 2 * PI / 360.0) * SIN(L * 2 * PI / 360.0) * COS(B * 2 * PI / 360.0) + COS(23.43929111 * 2 * PI / 360.0) * SIN(B * 2 * PI / 360.0));
  WriteLn(Format('Para T=%.15f:', [T]));
  WriteLn(Format('L (longitud eclíptica) = %.15f grados', [L]));
  WriteLn(Format('B (latitud eclíptica)   = %.15f grados', [B]));
  WriteLn(Format('R (distancia)           = %.15f AU', [R]));
  WriteLn(Format('sun_pos_1 (distancia)           = %.15f km', [sun_pos_1]));
  WriteLn(Format('sun_pos_2 (distancia)           = %.15f km', [sun_pos_2]));
  WriteLn(Format('sun_pos_3 (distancia)           = %.15f km', [sun_pos_3]));
end.
