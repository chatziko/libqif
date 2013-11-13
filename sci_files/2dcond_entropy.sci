pi=[0:0.01:1]

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

//GLeakage Conditional Entropy
function [Ep]= Ep(pi) 
    Ep = -log (Vp(pi));
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

//MinEntropy Conditional Entropy
function [Epm]= Epm(pi) 
    Epm = -log (Vpm(pi));
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

//Function G for Guessing
function [G]= G(p)
    G=0;
    n = size (p);
    for x=1:n
       G = G + x * p(x);
    end
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