function approximateZero = naive_Newton(fnc, x0, tol)
    maxNumIter = 20;
    itCounter = 0;  % Initialize iteration counter.
    
    syms fp;   % derivative
    fp = diff(fnc);
    
    xcur = x0;
    while ((abs(subs(fnc,xcur))>tol)  &&  (itCounter < maxNumIter) ) 
            xcur = double(xcur - subs(fnc,xcur)/subs(fp,xcur));
            itCounter = itCounter+1;
    end

    if(abs(subs(fnc,xcur))>tol )
            disp(['Warning:  Tolerance not met after ' num2str(maxNumIter) ' iterations.']);
    end

    approximateZero = xcur;
    
