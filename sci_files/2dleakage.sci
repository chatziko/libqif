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

//GLeakage Conditional Vulnerability Function
function [Vp]= Vp(pi)
    p= [pi (1-pi)]
    //sum for all output Y
    Vp=0;
    [n,i]= size (C);
    for y=1:i
    //test for all the guesses W
    [m,i]= size (g);
    r= zeros(n,1);
    for w=1:m 
    //sum for all the inputs X
    for x=1:n
       r(w)=r(w)+p(x)*g(w,x)*C(x,y);
      // continue;
    end
    //continue;
    end
    
    //take the maximum w and sum
    Vp= Vp+ max(r);
    //continue;
    end
endfunction

//Multiplicative G-Leakage
function [ML]= ML(pi) 
    ML= log (Vp(pi)/V(pi));
endfunction

//MinEntropy Vulnerability Function
function [Vm]= Vm(pi)
    p= [pi (1-pi)]
    //take the maximum pi
    Vm= max(p);
endfunction

//MinEntropy Conditional Vulnerability Function
function [Vpm]= Vpm(pi)
    p= [pi (1-pi)]
    //sum for all output Y
    Vpm=0;
    [n,i]= size (C);
    for y=1:i
    r= zeros(n,1);
    //for all the inputs X
    for x=1:n
       r(x)=p(x)*C(x,y);
    end
    
    //take the maximum w and sum
    Vpm= Vpm+ max(r);
    end
endfunction

//MinEntropy Leakage
function [MLm]= MLm(pi) 
    MLm= log (Vpm(pi)/Vm(pi));
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

//Shannon Conditional Entropy Function
function [Eps]= Eps(pi)
    p= [pi (1-pi)]
    Eps=0;
    //sum over x
    [n,i]= size (C);
    for x=1:n
        //sum over y
        sumY=0;
        for y=1:i
            sumY= SumY + C[x,y] * log (C[x,y]);
        end
       Eps= Eps+p(x) * sumY;
    end
    Eps=-Eps;
endfunction

//Shannon Leakage
function [MLs]= MLs(pi) 
    MLs= Es(pi)-Eps(pi));
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
function [Eg]= Eg(pi)
    p= [pi (1-pi)]
    p= gsort(p,'g','d');
    Eg=G(p);
endfunction

//Guessing Conditional Entropy
function [Epg]= Epg(pi)
    Epg=0;
    //for all y
    [n,i]= size (C);
    for y=1:i
        p= [(pi*C[1,y]) ((1-pi)*C[2,y])]
        p= gsort(p,'g','d');
        Epg= Epg + G(p);
    end    
endfunction

//Guessing Leakage
function [MLg]= MLg(pi) 
    MLg= Eg(pi)-Epg(pi));
endfunction