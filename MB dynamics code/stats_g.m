function stats_g(m,app,value)
    mkfld();
    data2 = fopen('C:\Multibody Dynamics Generator\SN.txt','r');
    log = fscanf(data2,'%s');
    %data1 = fopen('stats.txt','wt');
    %ln = length(log);
    %temp = string(zeros(1,9));

    sn1 = strfind(log,"<");
    sn2 = strfind(log,">");
    st1 = strfind(log,"~");
    st = string(zeros(21,length(sn1)));
    stats_list = string(zeros(length(sn1),1));
    new = [5,6,10,11,15,16,20,21];
    apr = zeros(1,length(sn1));
    apw = zeros(1,length(sn1));
    epr = zeros(1,length(sn1));
    epw = zeros(1,length(sn1));
    Epr = zeros(1,length(sn1));
    Epw = zeros(1,length(sn1));
    Cpr = zeros(1,length(sn1));
    Cpw = zeros(1,length(sn1));

    for i = 1:length(sn1)
        st(1,i) = join(log(sn1(i)+1:sn2(i)-1));
        for j = 2:21
            if j == new(1)
                if log(sn2(i)+2-1) == "0"
                    apr(i) = 0;
                    st(j,i) = apr(i);
                else
                    apr(i) = 100*str2double(log(sn2(i)+3-1))/str2double(log(sn2(i)+2-1));
                    st(j,i) = apr(i);
                end
            elseif j == new(2)
                if log(sn2(i)+2-1) == "0"
                    apw(i) = 0;
                    st(j,i) = apw(i);
                else
                    apw(i) = 100*str2double(log(sn2(i)+4-1))/str2double(log(sn2(i)+2-1));
                    st(j,i) = apw(i);
                end
            elseif j == new(3)
                if log(sn2(i)+2-1) == "0"
                    epr(i) = 0;
                    st(j,i) = epr(i);
                else
                    epr(i) = 100*str2double(log(sn2(i)+8-1))/str2double(log(sn2(i)+7-1));
                    st(j,i) = epr(i);
                end
            elseif j == new(4)
                if log(sn2(i)+2-1) == "0"
                    epw(i) = 0;
                    st(j,i) = epw(i);
                else
                    epw(i) = 100*str2double(log(sn2(i)+9-1))/str2double(log(sn2(i)+7-1));
                    st(j,i) = epw(i);
                end
            elseif j == new(5)
                if log(sn2(i)+2-1) == "0"
                    Epr(i) = 0;
                    st(j,i) = Epr(i);
                else
                    Epr(i) = 100*str2double(log(sn2(i)+13-1))/str2double(log(sn2(i)+12-1));
                    st(j,i) = Epr(i);
                end
            elseif j == new(6)
                if log(sn2(i)+2-1) == "0"
                    Epw(i) = 0;
                    st(j,i) = Epw(i);
                else
                    Epw(i) = 100*str2double(log(sn2(i)+14-1))/str2double(log(sn2(i)+12-1));
                    st(j,i) = Epw(i);
                end
            elseif j == new(7)
                if log(sn2(i)+2-1) == "0"
                    Cpr(i) = 0;
                    st(j,i) = Cpr(i);
                else
                    Cpr(i) = 100*str2double(log(sn2(i)+18-1))/str2double(log(sn2(i)+17-1));
                    st(j,i) = Cpr(i);
                end
            elseif j == new(8)
                if log(sn2(i)+2-1) == "0"
                    Cpw(i) = 0;
                    st(j,i) = Cpw(i);
                else
                    Cpw(i) = 100*str2double(log(sn2(i)+19-1))/str2double(log(sn2(i)+17-1));
                    st(j,i) = Cpw(i);
                end
            else
                st(j,i) = log(sn2(i)+j-1);
            end
        end
    end

    for i = 1:length(sn1)
        stats_list(i) = sprintf("%s\n\nNumber of angular velocity test attempts: %s\nNumber of correct answers: %s\nNumber of incorrect answers: %s\nPercentage of correct answers: %s\nPercentage of incorrect answers: %s\n\nNumber of Energy derivation test attempts: %s\nNumber of correct answers: %s\nNumber of incorrect answers: %s\nPercentage of correct answers: %s\nPercentage of incorrect answers: %s\n\nNumber of Dynamics derivation test attempts: %s\nNumber of correct answers: %s\nNumber of incorrect answers: %s\nPercentage of correct answers: %s\nPercentage of incorrect answers: %s \n\nNumber of Feedback Linearization test attempts: %s\nNumber of correct answers: %s\nNumber of incorrect answers: %s\nPercentage of correct answers: %s\nPercentage of incorrect answers: %s\n\n\n",st(1,i),st(2,i),st(3,i),st(4,i),st(5,i),st(6,i),st(7,i),st(8,i),st(9,i),st(10,i),st(11,i),st(12,i),st(13,i),st(14,i),st(15,i),st(16,i),st(17,i),st(18,i),st(19,i),st(20,i),st(21,i));
    end

    stats_list1 = join(sort(stats_list));
    stat_text = sprintf("The individual stats are as follows:\n\n\n %s",stats_list1);
    if m == 1
        app.TextArea_2.Value = stat_text;
        drawnow();
    end
    
    if m == 2
        for k = 1:length(stats_list)
            if strfind(stats_list(k),value) >=1
                output = sprintf("The stats of the student are as follows:\n\n\n %s",join(stats_list(k)));
            end
        end
        app.TextArea_2.Value = output;
        drawnow();
    end

    [aprn,eprn,Eprn,Cprn] = data_sep(apr,apw,epr,epw,Epr,Epw,Cpr,Cpw,sn1);
    
    if m == 3
        figure()
        hist([(apr)',(epr)',(Epr)',(Cpr)']);
        
        figure()
        h1 = histogram('Categories',{'0%-50%','50%-70%','70%-100%'},'BinCounts',aprn);
        h1.BarWidth = 0.5;
        ylabel("Number of students")
        xlabel("Categories")
        title("Histogram of Angular velocity test");
        ylim([0 length(sn1)+2])
        
        figure()
        h1 = histogram('Categories',{'0%-50%','50%-70%','70%-100%'},'BinCounts',eprn);
        h1.BarWidth = 0.5;
        ylabel("Number of students")
        xlabel("Categories")
        title("Histogram of Energy derivation test");
        ylim([0 length(sn1)+2])
        
        figure()
        h1 = histogram('Categories',{'0%-50%','50%-70%','70%-100%'},'BinCounts',Eprn);
        h1.BarWidth = 0.5;
        ylabel("Number of students")
        xlabel("Categories")
        title("Histogram of Dynamics derivation test");
        ylim([0 length(sn1)+2])
        
        figure()
        h1 = histogram('Categories',{'0%-50%','50%-70%','70%-100%'},'BinCounts',Cprn);
        h1.BarWidth = 0.5;
        ylabel("Number of students")
        xlabel("Categories")
        title("Histogram of Feedback Linearization test");
        ylim([0 length(sn1)+2])
    end
    
    if m == 4
        figure()
        pie(aprn,{'0%-50%','50%-70%','70%-100%'});
        title("Pie chart of angular velocity test");
        
        figure()
        pie(eprn,{'0%-50%','50%-70%','70%-100%'});
        title("Pie chart of Energy derivation test");
        
        figure()
        pie(Eprn,{'0%-50%','50%-70%','70%-100%'});
        title("Pie chart of Dynamics derivation test");
        
        figure()
        pie(Cprn,{'0%-50%','50%-70%','70%-100%'});
        title("Pie chart of Feedback Linearization test");
    end
    %fclose(data1);
end