function [pass,admin] = signup(r_o_s, Student_number, password, password2, app)
    mkfld()
    data2 = fopen('C:\Multibody Dynamics Generator\SN.txt','r');
    log = fscanf(data2,'%s');
    fclose(data2);
    ln = length(log);
    data1 = fopen('C:\Multibody Dynamics Generator\SN.txt','wt');
    temp = string(zeros(21,1));
    i = 2;
    temp_p = "";
    pass = 0;
    admin = 0;
    admin_check = strcmp(Student_number,"Administrator");
    password_check = 1;

    if ln >= 1
        u_check = strfind(log,Student_number);
        if u_check>=1
            if (r_o_s == 2) && admin_check < 1
                fprintf(data1,'%s',log);
                app.TextArea_2.Value = "Account already exists, try signing in";
            else
                while (log((u_check)-i)~= "~")
                    temp_p = temp_p+log((u_check)-i);
                    i = i+1;
                end
                temp_p = reverse(split(temp_p));
                if strcmp(password,temp_p)
                    fprintf(data1,'%s',log);
                    assignin('base','logged_in',Student_number);
                    app.TextArea_2.Value = "sign in successful";
                    pass = 1;
                    assignin('base','pass',pass);
                else
                    fprintf(data1,'%s',log);
                    app.TextArea_2.Value = "Incorrect password";
                end
            end
        elseif admin_check >=1 && r_o_s == 1
            if strcmp(password,"Admin12*")
                fprintf(data1,'%s',log);
                app.TextArea_2.Value = "sign in successful";
                admin = 1;
                assignin('base','admin',admin);
            else
                fprintf(data1,'%s',log);
                app.TextArea_2.Value = "Incorrect password";
            end

        else
            if (r_o_s == 1)
                fprintf(data1,'%s',log);
                app.TextArea_2.Value = "Account does not exist, try signing up";
            else
                if strcmp(password,password2)
                    p = split(password,"");
                    for i = 1:length(split(password,""))
                        if (p(i) == "~") || (p(i) == "<") || (p(i)== ">")
                            password_check = 0;
                        end
                    end
                    if password_check

                        if ~strcmp(Student_number,"Administrator")
                            % Student number
                            temp(1,:) = "~"+password+"<"+Student_number+">";
                            log = [log,(temp)'];
                            fprintf(data1,'%s',log);
                            app.TextArea_2.Value = "sign up successful, use these credentials to log in";
                        else
                            fprintf(data1,'%s',log);
                            app.TextArea_2.Value = "Enter a valid student number";
                        end
                    else
                        app.TextArea_2.Value = "Invalid password! Do not use '~', '<' or '>' characters";
                    end
                else
                    fprintf(data1,'%s',log);
                    app.TextArea_2.Value = "passwords did not match, please try again";
                end
            end
        end
    else
        if (r_o_s == 1) && admin_check < 1
            app.TextArea_2.Value = "Account does not exist, try signing up";

        elseif admin_check >=1 && r_o_s == 1
            if strcmp(password,"Admin12*")
                fprintf(data1,'%s',log);
                app.TextArea_2.Value = "sign in successful";
                admin = 1;
                assignin('base','admin',admin);
            else
                app.TextArea_2.Value = "Incorrect password";
            end

        else
            if strcmp(password,password2)
                p = split(password,"");
                for i = 1:length(split(password,""))
                    if (p(i)== "~") || (p(i) == "<") || (p(i) == ">")
                            password_check = 0;
                    end
                end
                if password_check
                    if ~strcmp(Student_number,"Administrator")
                        % Student number
                        temp(1,:) = "~"+password+"<"+Student_number+">";
                        log = [log,(temp)'];
                        fprintf(data1,'%s',log);
                        app.TextArea_2.Value = "sign up successful, use these credentials to log in";
                    else
                        fprintf(data1,'%s',log);
                        app.TextArea_2.Value = "Enter a valid student number";
                    end
                else
                    app.TextArea_2.Value = "Invalid password! Do not use '~', '<' or '>' characters";
                end
            else
                app.TextArea_2.Value = "passwords did not match, please try again";
            end
        end
    end

    fclose(data1);
end