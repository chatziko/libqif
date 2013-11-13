pi=[0:0.01:1]

//GLeakage Vulnerability Function
function [V]= V(pi,pi2)
    if pi+pi2>1 then
        V=0;
    else
    p= [pi pi2 (1-pi-pi2)]
    //test for all the guesses W
    [m,i]= size (g);
    [n,i]= size (C);
    r= zeros(n,1);
    for w=1:m
    //sum for all the inputs X
    for x=1:n
       r(w)=r(w)+p(x)*g(w,x);
      // continue;
    end
    //continue;
    end
    //take the maximum w
    V= max(r); 
    end
endfunction

//GLeakage Entropy
function [E]= E(pi,pi2) 
    E= -log(V(pi,pi2));
endfunction

//MinEntropy Vulnerability Function
function [Vm]= Vm(pi,pi2)
    if pi+pi2>1 then
        Vm=0;
    else
    p= [pi pi2 (1-pi-pi2)]
    Vm= max(p); 
    end
endfunction

//MinEntropy Entropy
function [Em]= Em(pi,pi2) 
    Em= -log(Vm(pi,pi2));
endfunction

//Shannon Entropy Function
function [Es]= Es(pi,pi2)
    if pi+pi2>1 then
        Es=0;
    else
        p= [pi pi2 (1-pi-pi2)]
        Es=0;
        //sum over x
        [n,i]= size (C);
        for x=1:n
           Es= Es+p(x) * log(p(x));
        end
        Es=-Es;
    end
endfunction

//Function G for Guessing
function [G]= G(p)
    G=0;
    n = size (p);
    for x=1:n
       G = G + x * p(x);
    end
endfunction

//Guessing Entropy
function [Eg]= Eg(pi,pi2)
    if pi+pi2>1 then
        Vm=0;
    else
        p= [pi pi2 (1-pi-pi2)]
        p= gsort(p,'g','d');
        Eg=G(p);
    end
endfunction