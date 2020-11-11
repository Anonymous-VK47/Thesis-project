function [aprn,eprn,Eprn,Cprn] = data_sep(apr,apw,epr,epw,Epr,Epw,Cpr,Cpw,sn1)

    aprn = zeros(1,3);
    apwn = zeros(1,3);
    eprn = zeros(1,3);
    epwn = zeros(1,3);
    Eprn = zeros(1,3);
    Epwn = zeros(1,3);
    Cprn = zeros(1,3);
    Cpwn = zeros(1,3);

    for x = 1:length(sn1)
        if apr(x)<50
            aprn(1) = aprn(1)+1;
        elseif apr(x)>=50 && apr(x)<70
            aprn(2) = aprn(2)+1;
        else
            aprn(3) = aprn(3)+1;
        end
        if apw(x)<50
            apwn(1) = apwn(1)+1;
        elseif apw(x)>=50 && apw(x)<70
            apwn(2) = apwn(2)+1;
        else
            apwn(3) = apwn(3)+1;
        end
        if epr(x)<50
            eprn(1) = eprn(1)+1;
        elseif epr(x)>=50 && epr(x)<70
            eprn(2) = eprn(2)+1;
        else
            eprn(3) = eprn(3)+1;
        end
        if epw(x)<50
            epwn(1) = epwn(1)+1;
        elseif epw(x)>=50 && epw(x)<70
            epwn(2) = epwn(2)+1;
        else
            epwn(3) = epwn(3)+1;
        end
        if Epr(x)<50
            Eprn(1) = Eprn(1)+1;
        elseif Epr(x)>=50 && Epr(x)<70
            Eprn(2) = Eprn(2)+1;
        else
            Eprn(3) = Eprn(3)+1;
        end
        if Epw(x)<50
            Epwn(1) = Epwn(1)+1;
        elseif Epw(x)>=50 && Epw(x)<70
            Epwn(2) = Epwn(2)+1;
        else
            Epwn(3) = Epwn(3)+1;
        end
        if Cpr(x)<50
            Cprn(1) = Cprn(1)+1;
        elseif Cpr(x)>=50 && Cpr(x)<70
            Cprn(2) = Cprn(2)+1;
        else
            Cprn(3) = Cprn(3)+1;
        end
        if Cpw(x)<50
            Cpwn(1) = Cpwn(1)+1;
        elseif Cpw(x)>=50 && Cpw(x)<70
            Cpwn(2) = Cpwn(2)+1;
        else
            Cpwn(3) = Cpwn(3)+1;
        end
    end
    
    out1 = [aprn;eprn;Eprn;Cprn];

    %{
    for y = 1:3  
        aprn(y) = 100*aprn(y)/length(sn1);
        %apwn(y) = 100*apwn(y)/length(sn1);
        eprn(y) = 100*eprn(y)/length(sn1);
        %epwn(y) = 100*epwn (y)/length(sn1);
        Eprn(y) = 100*Eprn (y)/length(sn1);
        %Epwn(y) = 100*Epwn(y)/length(sn1);
        Cprn(y) = 100*Cprn (y)/length(sn1);
        %Cpwn(y) = 100*Cpwn(y)/length(sn1); 
    end
    out2 = [aprn;eprn;Eprn;Cprn];
    %}
end