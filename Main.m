% Given Inputs

A = [1 3 3 3 2;
    -5 8 3 5 2;
     0 0 6 0 0;
     1 0 -1 3 0;
     5 -7 1 -5 0];
B = [2; 1; 0; 1; 0];
C = [1, 1, -1, -1, 0];

%///////////////////////////////
%Step 1: Find CC and OO 
CC = ctrb(A,B);
lengthCC = size(CC,1);
CC_rref = rref(CC);

OO = obsv(A,C);
lengthOO = size(OO,1);
OO_rref = rref(OO);

%///////////////////////////////
%Step 2: From C find Vc

% Vc should be the size of n x p(CC) 
Vc = zeros(lengthCC, rank(CC_rref));

% For each LI column, copy it to Vc
for i=1:1:rank(CC_rref)
Vc(:,i) = CC(:,i);
end

%Display Vc
disp("Vc: ")
disp(Vc)

%///////////////////////////////
%Step 2 continued
%V(!o) should be the dependant columns of rref(O)
%V(!o) should have size n x (n-p(OO))
Vo_bar = zeros(lengthOO, lengthOO-rank(OO_rref));

% Copy over the dependent columns from OO_rref to Vo_bar
% For every dependent column, a 1 is added for each column after the rank
j=0;
for i=rank(OO_rref):1:lengthOO
    if(i>rank(OO_rref))
        j=j+1;
        Vo_bar(:,i-rank(OO_rref)) = -1*OO_rref(:,i);
        Vo_bar(rank(OO_rref)+j,i-rank(OO_rref)) = 1;
    end
end

%Display Vobar
disp("Vobar: ")
disp(Vo_bar)

%///////////////////////////////
%Step 3: Finding V c,obar

% Generate [Qc -Qobar]
Qc_negQobar = [Vc -1*Vo_bar];
rref_Qc_negQobar = rref(Qc_negQobar);

% independent_cols = find(any(rref_Qc_negQobar(:, :) == 1, 1))
% for i = 1:1:size(rref_Qc_negQobar,2)
%     if ~ismember(i, independent_cols)
%         dependent_col = -1*rref_Qc_negQobar(:,i)
%     end
% end
% dependent_col = nonzeros(dependent_col);

dependent_col = -1*null(rref_Qc_negQobar);
dependent_col = dependent_col(1:rank(CC_rref), :);

% Find Qc * x 
Qcx = Vc*dependent_col;
Vcobar = Qcx;

disp("Vcobar:");
disp(Vcobar);

%///////////////////////
% Step 4

Vcobar_transpose = transpose(Vcobar);
Vucobar_basis = null(Vcobar_transpose);

% Get [Qcobar_orth -Qc]
Qcobar_negQc = [Vucobar_basis -Vc];
rref_Qcobar_negQc = rref(Qcobar_negQc);

%uw span is needed to get vco
% the size is the m of vcobar x dependant column numbers of rref_Qcobar_negQc
uw_span = zeros(size(Vucobar_basis,2)+1,size(rref_Qcobar_negQc,2)-rank(rref_Qcobar_negQc));
% Copy the dependant columns over
j=0;
for i=rank(Qcobar_negQc)+1:1:size(Qcobar_negQc,2)
    j=j+1;
    uw_span(:,j) = rref_Qcobar_negQc(:,i);
end

% We want it to be for the controllable part only
uw_span = -uw_span;
uw_span = uw_span(1:size(Vucobar_basis,2), :);
disp("uw_span:");
disp(uw_span);

%Finding Vco by multiplying the basis to every dependant vector
for i=1:1:size(uw_span,2)
    Vco(:,i) = Vucobar_basis*uw_span(:,i);
end
disp("Vco:");
disp(Vco);

%///////////////////////
% Step 5
Qcobarorth_negqobar = [Vucobar_basis -1*Vo_bar];
rref_Qcobarorth_negqobar = rref(Qcobarorth_negqobar);

% 
% independent_cols_again = find(any(rref_Qcobarorth_negqobar(:, :) == 1, 1))
% for i = 1:1:size(rref_Qcobarorth_negqobar,2)
%     if ~ismember(i, independent_cols_again)
%         dependent_col = -1*rref_Qcobarorth_negqobar(:,i)
%     end
% end
% dependent_col = dependent_col(1:size(Vucobar_basis,2), :)

% For finding Vcbarobar
dependent_col = null(rref_Qcobarorth_negqobar);
dependent_col = dependent_col(1:size(Vucobar_basis,2), :);

vcbarobar = Vucobar_basis*dependent_col;
disp('v cbar obar:');
disp(vcbarobar);

%//////////////////////////
% Step 6:

U = transpose([Vcobar Vco vcbarobar]);
U_rref = rref(U);
Vcbaro = 2*null(U_rref);
disp('V cbar O:');
disp(Vcbaro);

%/////////////////
% Step 7:

Q = [transpose(Vcobar_transpose) Vco vcbarobar -1*Vcbaro ];
disp("Q:");
disp(Q);

Abar = inv(Q)*A*Q;
disp("A hat:");
disp(Abar);

Bbar = inv(Q)*B;
disp("B hat:");
disp(Bbar);

Cbar = C*Q;
disp(" C hat:");
disp(Cbar);