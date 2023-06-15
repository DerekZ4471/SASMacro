data FEV1;
    input patient basefev1 fev1h1 fev1h2 fev1h3 fev1h4 fev1h5 
          fev1h6 fev1h7 fev1h8 drug $;
  datalines;
   201  2.46  2.68 2.76 2.50 2.30 2.14 2.40 2.33 2.20 a
   202  3.50  3.95 3.65 2.93 2.53 3.04 3.37 3.14 2.62 a
   203  1.96  2.28 2.34 2.29 2.43 2.06 2.18 2.28 2.29 a
   204  3.44  4.08 3.87 3.79 3.30 3.80 3.24 2.98 2.91 a
   205  2.80  4.09 3.90 3.54 3.35 3.15 3.23 3.46 3.27 a
   206  2.36  3.79 3.97 3.78 3.69 3.31 2.83 2.72 3.00 a
   207  1.77  3.82 3.44 3.46 3.02 2.98 3.10 2.79 2.88 a
   208  2.64  3.67 3.47 3.19 2.19 2.85 2.68 2.60 2.73 a
   209  2.30  4.12 3.71 3.57 3.49 3.64 3.38 2.28 3.72 a
   210  2.27  2.77 2.77 2.75 2.75 2.71 2.75 2.52 2.60 a
   211  2.44  3.77 3.73 3.67 3.56 3.59 3.35 3.32 3.18 a
   212  2.04  2.00 1.91 1.88 2.09 2.08 1.98 1.70 1.40 a
   214  2.77  3.36 3.42 3.28 3.30 3.31 2.99 3.01 3.08 a
   215  2.96  4.31 4.02 3.38 3.31 3.46 3.49 3.38 3.35 a
   216  3.11  3.88 3.92 3.71 3.59 3.57 3.48 3.42 3.63 a
   217  1.47  1.97 1.90 1.45 1.45 1.24 1.24 1.17 1.27 a
   218  2.73  2.91 2.99 2.87 2.88 2.84 2.67 2.69 2.77 a
   219  3.25  3.59 3.54 3.17 2.92 3.48 3.05 3.27 2.96 a
   220  2.73  2.88 3.06 2.75 2.71 2.83 2.58 2.68 2.42 a
   221  3.30  4.04 3.94 3.84 3.99 3.90 3.89 3.89 2.98 a
   222  2.85  3.38 3.42 3.28 2.94 2.96 3.12 2.98 2.99 a
   223  2.72  4.49 4.35 4.38 4.36 3.77 4.23 3.83 3.89 a
   224  3.68  4.17 4.30 4.16 4.07 3.87 3.87 3.85 3.82 a
   232  2.49  3.73 3.51 3.16 3.26 3.07 2.77 2.92 3.00 a
   201  2.30  3.41 3.48 3.41 3.49 3.33 3.20 3.07 3.15 c
   202  2.91  3.92 4.02 4.04 3.64 3.29 3.10 2.70 2.69 c
   203  2.08  2.52 2.44 2.27 2.23 2.01 2.26 2.34 2.44 c
   204  3.02  4.43 4.30 4.08 4.01 3.62 3.23 2.46 2.97 c
   205  3.26  4.55 4.58 4.44 4.04 4.33 3.87 3.75 3.81 c
   206  2.29  4.25 4.37 4.10 4.20 3.84 3.43 3.79 3.74 c
   207  1.96  3.00 2.80 2.59 2.42 1.61 1.83 1.21 1.50 c
   208  2.70  4.06 3.98 4.06 3.93 3.61 2.91 2.07 2.67 c
   209  2.50  4.37 4.06 3.68 3.64 3.17 3.37 3.20 3.25 c
   210  2.35  2.83 2.79 2.82 2.79 2.80 2.76 2.64 2.69 c
   211  2.34  4.06 3.68 3.59 3.27 2.60 2.72 2.22 2.68 c
   212  2.20  2.82 1.90 2.57 2.30 1.67 1.90 2.07 1.76 c
   214  2.78  3.18 3.13 3.11 2.97 3.06 3.27 3.24 3.33 c
   215  3.43  4.39 4.63 4.19 4.00 4.01 3.66 3.47 3.22 c
   216  3.07  3.90 3.98 4.09 4.03 4.07 3.56 3.83 3.75 c
   217  1.21  2.31 2.19 2.21 2.09 1.75 1.72 1.80 1.36 c
   218  2.60  3.19 3.18 3.15 3.14 3.08 2.96 2.97 2.85 c
   219  2.61  3.54 3.45 3.25 3.01 3.07 2.65 2.47 2.55 c
   220  2.48  2.99 3.02 3.02 2.94 2.69 2.66 2.68 2.70 c
   221  3.73  4.37 4.20 4.17 4.19 4.07 3.86 3.89 3.89 c
   222  2.54  3.26 3.39 3.27 3.20 3.32 3.09 3.25 3.15 c
   223  2.83  4.72 4.97 4.99 4.96 4.95 4.82 4.56 4.49 c
   224  3.47  4.27 4.50 4.34 4.00 4.11 3.93 3.68 3.77 c
   232  2.79  4.10 3.85 4.27 4.01 3.78 3.14 3.94 3.69 c
   201  2.14  2.36 2.36 2.28 2.35 2.31 2.62 2.12 2.42 p
   202  3.37  3.03 3.02 3.19 2.98 3.01 2.75 2.70 2.84 p
   203  1.88  1.99 1.62 1.65 1.68 1.65 1.85 1.96 1.30 p
   204  3.10  3.24 3.37 3.54 3.31 2.81 3.58 3.76 3.05 p
   205  2.91  3.35 3.92 3.69 3.97 3.94 3.63 2.92 3.31 p
   206  2.29  3.04 3.28 3.17 2.99 3.31 3.21 2.98 2.82 p
   207  2.20  2.46 3.22 2.65 3.02 2.25 1.50 2.37 1.94 p
   208  2.70  2.85 2.81 2.96 2.69 2.18 1.91 2.21 1.71 p
   209  2.25  3.45 3.48 3.80 3.60 2.83 3.17 3.22 3.13 p
   210  2.48  2.56 2.52 2.67 2.60 2.68 2.64 2.65 2.61 p
   211  2.12  2.19 2.44 2.41 2.55 2.93 3.08 3.11 3.06 p
   212  2.37  2.14 1.92 1.75 1.58 1.51 1.94 1.84 1.76 p
   214  2.73  2.57 3.08 2.62 2.91 2.71 2.39 2.42 2.73 p
   215  3.15  2.90 2.80 3.17 2.39 3.01 3.22 2.75 3.14 p
   216  2.52  3.02 3.21 3.17 3.13 3.38 3.25 3.29 3.35 p
   217  1.48  1.35 1.15 1.24 1.32 0.95 1.24 1.04 1.16 p
   218  2.52  2.61 2.59 2.77 2.73 2.70 2.72 2.71 2.75 p
   219  2.90  2.91 2.89 3.01 2.74 2.71 2.86 2.95 2.66 p
   220  2.83  2.78 2.89 2.77 2.77 2.69 2.65 2.84 2.80 p
   221  3.50  3.81 3.77 3.78 3.90 3.80 3.78 3.70 3.61 p
   222  2.86  3.06 2.95 3.07 3.10 2.67 2.68 2.94 2.89 p
   223  2.42  2.87 3.08 3.02 3.14 3.67 3.84 3.55 3.75 p
   224  3.66  3.98 3.77 3.65 3.81 3.77 3.89 3.63 3.74 p
   232  2.88  3.04 3.00 3.24 3.37 2.69 2.89 2.89 2.76 p
  ;
  run;
   /*---convert data to univariate form                  ---*/
  data fev1uni; set fev1;
     array f{8} fev1:;
  
     do hour=1 to 8; fev = f{hour}; output; end;
     drop fev1:;
  run;
  data fev1uni; set fev1uni; rename fev=fev1; run;
  proc sort data=fev1uni; by drug hour;
  run;
  