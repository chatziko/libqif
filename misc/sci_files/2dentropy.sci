pi=[0:0.01:1]

//GLeakage Vulnerability Function
function [V]= V(pi)
    p= [pi (1-pi)]
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
endfunction

//GLeakage Entropy
function [E]= E(pi) 
    E = -log (V(pi));
endfunction

//MinEntropy Vulnerability Function
function [Vm]= Vm(pi)
    p= [pi (1-pi)]
    //take the maximum pi
    Vm= max(p);
endfunction

//MinEntropy Entropy
function [Em]= Em(pi) 
    Em = -log (Vm(pi));
endfunction

//Shannon Entropy Function
function [Es]= Es(pi)
    p= [pi (1-pi)]
    Es=0;
    //sum over x
    [n,i]= size (C);
    for x=1:n
       Es= Es+p(x) * log(p(x));
    end
    Es=-Es;
endfunction

//Guessing Entropy
//function [Eg]= Eg(pi)
//    p= [pi (1-pi)]
//    p= gsort(p,'g','i');
//    Eg=0;
//    //sum over x
//    [n,i]= size (C);
//    for x=1:n
//       Eg= Eg + x * p(x);
//    end
//endfunction

//Function G for Guessing
function [G]= G(p)
    G=0;
    n = size (p);
    for x=1:n
       G = G + x * p(x);
    end
endfunction

//Guessing Entropy
function [Eg]= Eg(pi)
    p= [pi (1-pi)]
    p= gsort(p,'g','d');
    Eg=G(p);
endfunction
