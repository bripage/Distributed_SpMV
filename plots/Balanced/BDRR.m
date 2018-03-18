function output = BDRR(p,B,phi)
    EXP = load('Mval.txt');
    m = EXP(:,1);
    %//nnz = 8814880;
    d = 64; %size of double in bits
    %B = 60129542144; %network bandwidth in bits
    %nnzProc = nnz/m;
    
    output = (log2(p).*((m./p).*d + phi))./B;
    
end