pi=[0:0.01:1]

//GLeakage Conditional Vulnerability Function
function [Vp]= Vp(pi,pi2)
    if pi+pi2>1 then
        Vp=0;
    else
    p= [pi pi2 (1-pi-pi2)]
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
end
endfunction

//GLeakage Conditional Entropy
function [Ep]= Ep(pi,pi2) 
    Ep= -log(Vp(pi,pi2));
endfunction

//MinEntropy Conditional Vulnerability Function
function [Vpm]= Vpm(pi,pi2)
    if pi+pi2>1 then
        Vpm=0;
    else
        p= [pi pi2 (1-pi-pi2)]
        //sum for all output Y
        Vp=0;
        [n,i]= size (C);
        for y=1:i
            r= zeros(n,1);
            //for all the inputs X
            for x=1:n
               r(w)=p(x)*C(x,y);
            end
            
            //take the maximum w and sum
            Vpm= Vpm+ max(r);
        end
    end
endfunction

//MinEntropy Conditional Entropy
function [Epm]= Epm(pi,pi2) 
    Epm= -log(Vpm(pi,pi2));
endfunction

//Shannon Conditional Entropy Function
function [Eps]= Eps(pi,pi2)
    if pi+pi2>1 then
        Eps=0;
    else
        p= [pi pi2 (1-pi-pi2)]
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

//Guessing Conditional Entropy
function [Epg]= Epg(pi,pi2)
    Epg=0;
    //for all y
    [n,i]= size (C);
    for y=1:i
        p= [(pi*C[1,y]) (pi2*C[2,y]) ((1-pi-pi2)*C[3,y])]
        p= gsort(p,'g','d');
        Epg= Epg + G(p);
    end    
endfunction