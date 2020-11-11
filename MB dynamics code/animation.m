function animation(out, psi_Values)
    global number_of_links;

    %legend('$\phi_1$', '$\dot{\phi}_1$', '$\psi_2$', '$\dot{\psi}_2$', 'FontSize', 14, 'Interpreter', 'latex');
    %grid on;
    %close all;
    
    %% Positions
    
    % This is the base
    pos0=[zeros(length(psi_Values),1),zeros(length(psi_Values),1),zeros(length(psi_Values),1)];
    
    % This is the joint between link 1 and 2
    pos1 = [out(:,1), out(:,2), out(:,3)];
    
    if number_of_links >= 2
        % This is the position of the second joint/mass
        pos2 = [out(:,4), out(:,5), out(:,6)];
    end
    
    if number_of_links == 3
        % This is the position of the tip
        pos3 = [out(:,7), out(:,8), out(:,9)];
    end
    
    %% Set Up figure
    
    h = figure();
    xlabel({'X Position (m)'},'FontSize',14,'FontName','AvantGarde');
    ylabel({'Y Position (m)'},'FontSize',14,'FontName','AvantGarde');
    zlabel({'Z Position (m)'},'FontSize',14,'FontName','AvantGarde');
    if number_of_links == 1
        title({'Single-link robot arm'},'FontWeight','bold','FontSize',20,'FontName','AvantGarde');
    elseif number_of_links == 2
        title({'Two-link robot arm'},'FontWeight','bold','FontSize',20,'FontName','AvantGarde');
    else
        title({'Three-link robot arm'},'FontWeight','bold','FontSize',20,'FontName','AvantGarde');
    end
    grid on
    hold on;
    
    % Represent the links as lines
    h2 = line('Color', 'b', 'LineWidth', 3); 
    if number_of_links >= 2
        h1 = line('Color', 'b', 'LineWidth', 3);
    end
    if number_of_links == 3
       h3 = line('Color', 'b', 'LineWidth', 3);
    end
    
    % Initial Co-ordinates of the base, joint and the mass 
    p0 = plot3(pos0(1,1),pos0(1,2),pos0(1,3),'o','MarkerFaceColor','green','MarkerSize',12); 
    p1 = plot3(pos1(1,1),pos1(1,2),pos1(1,3),'o','MarkerFaceColor','red','MarkerSize',12);
        
    if number_of_links >= 2
        p2 = plot3(pos2(1,1),pos2(1,2),pos2(1,3),'o','MarkerFaceColor','blue','MarkerSize',12);
    end
    
    if number_of_links == 3
        p3 = plot3(pos3(1,1),pos3(1,2),pos3(1,3),'o','MarkerFaceColor','blue','MarkerSize',12);
    end 
 
    %% plot values continuously until all values from simulation have been used
    if ishandle(h)
        for i = 1:length(psi_Values)
            % Update co-ordinates of links
            if number_of_links == 3
                if ishandle(h)
                    h3.XData = [pos2(i,1), pos3(i,1)];
                    h3.YData = [pos2(i,2), pos3(i,2)];
                    h3.ZData = [pos2(i,3), pos3(i,3)];
                end
            end

            % Update co-ordinates of links
            if number_of_links >= 2
                if ishandle(h)
                    h1.XData = [pos1(i,1), pos2(i,1)];
                    h1.YData = [pos1(i,2), pos2(i,2)];
                    h1.ZData = [pos1(i,3), pos2(i,3)];
                end
            end

            % Update co-ordinates of links
            if ishandle(h)
                h2.XData = [pos1(i,1), pos0(i,1)];
                h2.YData = [pos1(i,2), pos0(i,2)];
                h2.ZData = [pos1(i,3), pos0(i,3)]; 

                % Update co-ordinates of mass, joint and base
                p1.XData = pos1(i,1);
                p1.YData = pos1(i,2);
                p1.ZData = pos1(i,3);
            end
            
            % Update co-ordinates of mass, joint and base
            if number_of_links == 3
                if ishandle(h)
                    p3.XData = pos3(i,1);
                    p3.YData = pos3(i,2);
                    p3.ZData = pos3(i,3);
                end
            end

            % Update co-ordinates of mass, joint and base
            if number_of_links >= 2
                if ishandle(h)
                    p2.XData = pos2(i,1);
                    p2.YData = pos2(i,2);
                    p2.ZData = pos2(i,3);
                end
            end

            % Update co-ordinates of mass, joint and base
            if ishandle(h)
                p0.XData = pos0(i,1);
                p0.YData = pos0(i,2);
                p0.ZData = pos0(i,3);
            end

            if ishandle(h)
                drawnow 

                if number_of_links == 1
                    xlim([-0.15 0])
                    ylim([0 0.15])
                    zlim([0 0.15])
                    view(-62,15) 
                elseif number_of_links == 2
                    xlim([-0.25 0.05])
                    ylim([-0.05 0.28])
                    zlim([-0.15 0.25])
                    view(30,35) 
                else
                    xlim([-0.25 0.3])
                    ylim([-0.25 0.3])
                    zlim([-0.25 0.3])
                    view(15,25)
                end
            end
        end
        %pause(0.5);
       % close(h);
    end
end