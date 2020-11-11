function UpdateStats(n,choice)
    mkfld();
    Student_number = evalin('base','Student_number');
    data2 = fopen('C:\Multibody Dynamics Generator\SN.txt','r');
    log = fscanf(data2,'%s');
    data1 = fopen('C:\Multibody Dynamics Generator\SN.txt','wt');
    %ln = length(log);
    %temp = string(zeros(13,1));

    check = strfind(log,Student_number)+9;

    switch n
        case 1
            % ang_vel attempts
            log(check+1) = string(str2double(log(check+1))+1);
            if choice
                % ang_vel correct
                log(check+2) = string(str2double(log(check+2))+1);
                log(check+4) = string(str2double(log(check+2))/str2double(log(check+1)));
            else
                % ang_vel incorrect
                log(check+3) = string(str2double(log(check+3))+1);
                log(check+5) = string(str2double(log(check+3))/str2double(log(check+1)));
            end
        case 2
            % Energy attempts
            log(check+6) = string(str2double(log(check+6))+1);
            if choice
                % Energy correct
                log(check+7) = string(str2double(log(check+7))+1);
                log(check+9) = string(str2double(log(check+7))/str2double(log(check+6)));
            else
                % Energy incorrect
                log(check+8) = string(str2double(log(check+8))+1);
                log(check+10) = string(str2double(log(check+8))/str2double(log(check+6)));
            end
        case 3
            % EOM attempts
            log(check+11) = string(str2double(log(check+11))+1);
            if choice
                % EOM correct
                log(check+12) = string(str2double(log(check+12))+1);
                log(check+14) = string(str2double(log(check+12))/str2double(log(check+11)));
            else
                % EOM incorrect
                log(check+13) = string(str2double(log(check+13))+1);
                log(check+15) = string(str2double(log(check+13))/str2double(log(check+11)));
            end
        otherwise
            % Control attempts
            log(check+16) = string(str2double(log(check+16))+1);
            if choice
                % Control correct
                log(check+17) = string(str2double(log(check+17))+1);
                log(check+19) = string(str2double(log(check+17))/str2double(log(check+16)));
            else
                % Control incorrect
                log(check+18) = string(str2double(log(check+18))+1);
                log(check+20) = string(str2double(log(check+18))/str2double(log(check+16)));
            end
    end

    fprintf(data1,'%s',log);
    fclose(data1);
end