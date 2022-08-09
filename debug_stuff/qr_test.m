    M_h = options.discretization.M;
    
    Ms = [2 5 10 20 40 60];
    error1 = [];
    error2 = [];
    error3 = [];
    error4 = [];
    error5 = [];
    error6 = [];
    
    for k = 1:numel(Ms)
        B = options.rom.B;
        B = B(:,1:Ms(k));
        themat = (M_h*B)';
        [Q,R,P] = qr(themat,0);
        error1(k) = norm(themat(:,P)-Q*R)
        RR = R(:,P);
        error2(k) = norm(themat-Q*RR)
        H = rank(themat)
        Q1 = Q(:,1:H);
        Q2 = Q(:,H+1:end);
        R1 = RR(1:H,:);
        error3(k) = norm(themat-Q1*R1)
        error4(k) = norm(M_h*B*Q2)
        
        a1 = (R1*R1')\(R1*yM);
        error5(k) = norm(M_h*B*Q1*a1 - yM)
        
        a2 = Q2'*(B'*(Om.*V));
        
        R = Q1*a1 + Q2*a2;
        error6(k) = norm(M_h*B*R - yM)
    end