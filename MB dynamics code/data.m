Student_number = string(input("Insert student number:\n",'s'));
data2 = fopen('SN.txt','r');
log = fscanf(data2,'%s');
data1 = fopen('SN.txt','wt');
ln = length(log);
temp = string(zeros(16,1));
if ln >= 1
    
    check = strfind(log,Student_number)+8;
    
    if check >=1
        % ang_vel attempts
        log(check+1) = string(str2double(log(check+1))+1);
        % ang_vel correct
        log(check+2) = string(str2double(log(check+2))+1);
        % ang_vel incorrect
        log(check+3) = string(str2double(log(check+3))+1);

        % Energy attempts
        log(check+4) = string(str2double(log(check+4))+1);
        % Energy correct
        log(check+5) = string(str2double(log(check+5))+1);
        % Energy incorrect
        log(check+6) = string(str2double(log(check+6))+1);

        % Energy attempts
        log(check+7) = string(str2double(log(check+7))+1);
        % Energy correct
        log(check+8) = string(str2double(log(check+8))+1);
        % Energy incorrect
        log(check+9) = string(str2double(log(check+9))+1);

        % EOM attempts
        log(check+10) = string(str2double(log(check+10))+1);
        % EOM correct
        log(check+11) = string(str2double(log(check+11))+1);
        % EOM incorrect
        log(check+12) = string(str2double(log(check+12))+1);

        % EOM attempts
        log(check+13) = string(str2double(log(check+13))+1);
        % EOM correct
        log(check+14) = string(str2double(log(check+14))+1);
        % EOM incorrect
        log(check+15) = string(str2double(log(check+15))+1);
        fprintf(data1,'%s',log);
    else
        % Student number
        temp(1,:) = Student_number;
        log = [log,(temp)'];
        fprintf(data1,'%s',log);
    end
else
        % Student number
        temp(1,:) = Student_number;
        fprintf(data1,'%s',temp);
end
  
fclose(data1);