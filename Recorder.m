classdef Recorder < handle
    %Recorder: stores the state trajectories and plots them
    %This class can store state trajectories, convert them to i.e. voltage
    %or current trajectories and plot them according to the wishes.
    properties
        name
        x_save
        v_limit_reached
        i_limit_reached
        f_limit_reached
        p_g_limit_reached
        q_g_limit_reached
        f_function_save
        h_function_save
        iteration               = 0
        transient_length        = 390
        start_disp              = 1
        m
        n
        node_legend
        node_legend_viol
        line_legend
    end
    methods
        function obj = Recorder(mystate,f_function,h_function, iterations, mygrid)
            x = mystate.getx();
            obj.m = length(mystate.i_re)/2;
            obj.n = length(mystate.v);
            obj.x_save = [x,zeros(length(x),iterations)];
            obj.f_function_save = [f_function,zeros(1,iterations)];
            obj.h_function_save = [h_function,zeros(1,iterations)];
            assert(obj.transient_length >= 1);
            assert(5*obj.n + 4*obj.m + 1 == length(x));
            
            obj.v_limit_reached = false(length(x),iterations);
            obj.i_limit_reached = false(length(mystate.i_re),iterations);
            obj.f_limit_reached = false(length(x),iterations);
            obj.p_g_limit_reached = false(length(x),iterations);
            obj.q_g_limit_reached = false(length(x),iterations);
            
            for i=1:obj.n
                obj.node_legend{i} = ['Node ', num2str(i)];
            end
            
            obj.node_legend_viol = obj.node_legend;
            obj.node_legend_viol{obj.n+1} = ['constraint violation'];
            
            for i=1:obj.m
                obj.line_legend{i} = ['Line ', num2str(i)];
                obj.line_legend{i+obj.m} = ['Line (opp. dir.)', num2str(i)];
            end
            
            obj.line_legend{2*obj.m+1} = ['constraint violation'];
            
            obj.name = mygrid.name;
            
            if strcmp(obj.name, 'wind')
                obj.transient_length = 250;
            else
                obj.transient_length = 390;
            end
        end
        
        function store(obj, mystate, f_function, h_function)
            obj.f_function_save(:,obj.iteration) = f_function;
            obj.h_function_save(:,obj.iteration) = h_function;
            obj.x_save(:,obj.iteration) = mystate.getx();
            obj.start_disp = min(obj.transient_length, obj.iteration);
        end
        
        function count(obj)
            obj.iteration = obj.iteration + 1;
        end
        
        function store_limits(obj, mystate, mygrid)
            obj.v_limit_reached(:,obj.iteration) = logical(Controller.v_limit_reached(mystate, mygrid));
            obj.i_limit_reached(:,obj.iteration) = logical(Controller.i_limit_reached(mystate, mygrid));
            obj.f_limit_reached(:,obj.iteration) = logical(Controller.f_limit_reached(mystate, mygrid));
            obj.p_g_limit_reached(1+2*mygrid.n:3*mygrid.n,obj.iteration) = logical(Controller.S_limit_reached(mystate, mygrid));
            obj.q_g_limit_reached(1+3*mygrid.n:4*mygrid.n,obj.iteration) = logical(Controller.S_limit_reached(mystate, mygrid));
        end
        function plotV(obj)
            figure(1);
            plot(gca, obj.x_save(1:obj.n,1:obj.iteration)');
            ylabel(gca,'Voltage [p.u.]');
            xlabel(gca,'Iterations');
        end
        function plotAll(obj)
            
            %f2.OuterPosition = [0 0 1440 930];
            %f2.WindowStyle = 'docked';
            mp = 2;
            np = 2;
            m = obj.m;
            n = obj.n;
            
            f2 = figure(2);
            clf
            f2.OuterPosition = [0 0 800 1000];
            v = subplot(mp,np,1);
            i = subplot(mp,np,2);
            p_g = subplot(mp,np,3);
            q_g = subplot(mp,np,4);
            box(v,'on')
            box(i,'on')
            box(p_g,'on')
            box(q_g,'on')            
                    
            hold(v,'on')
            plot(v,obj.x_save(1:obj.n,obj.start_disp:obj.iteration)');
            for k=1:obj.n
                x_values = 1:(obj.iteration - (obj.start_disp - 1));                            %shift the index according to the
                x_values = x_values(obj.v_limit_reached(k,obj.start_disp:obj.iteration));       %take out all the points where the limit was not reached
                y_values = obj.x_save(k,obj.start_disp:obj.iteration);                          %trim the voltage values down to the desired length
                y_values = y_values(obj.v_limit_reached(k,obj.start_disp:obj.iteration));      %take out all the points where the limit was not reached
                plot(v,x_values,y_values,'xr','LineWidth',2);
            end
            ylabel(v,'v: Voltage amplitude [p.u.]');
            xlabel(v,'Iterations');
            warning('off','MATLAB:legend:IgnoringExtraEntries');
            legend(v,obj.node_legend_viol);
            hold(v,'off')
            
            i_mag = sqrt(obj.x_save(5*n+1:5*n+2*m,obj.start_disp:obj.iteration).^2 + obj.x_save(5*n+2*m+1:5*n+4*m,obj.start_disp:obj.iteration).^2);        
            hold(i,'on')
            plot(i, i_mag');
            for k=1:size(i_mag,1)
                x_values = 1:(obj.iteration - (obj.start_disp - 1));                            %shift the index according to the
                x_values = x_values(obj.i_limit_reached(k,obj.start_disp:obj.iteration));       %take out all the points where the limit was not reached
                y_values = i_mag(k,:);                          %trim the voltage values down to the desired length
                y_values = y_values(obj.i_limit_reached(k,obj.start_disp:obj.iteration));      %take out all the points where the limit was not reached
                plot(i,x_values,y_values,'xr','LineWidth',2);
            end
            ylabel(i,'i: Line Current Amplitude [p.u.]');
            xlabel(i,'Iterations');
            legend(i,obj.line_legend);
            hold(i,'off')
                   
            hold(p_g,'on')
            plot(p_g, obj.x_save(2*n+1:3*n,obj.start_disp:obj.iteration)');
            p_g.ColorOrderIndex = 1;
            plot(p_g, obj.x_save(4*n+1:5*n,obj.start_disp:obj.iteration)','--');
            for k=(2*n+1):(3*n)
                x_values = 1:(obj.iteration - (obj.start_disp - 1));                            %shift the index according to the
                x_values = x_values(obj.p_g_limit_reached(k,obj.start_disp:obj.iteration));       %take out all the points where the limit was not reached
                y_values = obj.x_save(k,obj.start_disp:obj.iteration);                          %trim the voltage values down to the desired length
                y_values = y_values(obj.p_g_limit_reached(k,obj.start_disp:obj.iteration));      %take out all the points where the limit was not reached
                plot(p_g,x_values,y_values,'xr','LineWidth',2);
            end
            ylabel(p_g,{'p_g: Active Power generated [p.u]';'(dotted line is p_{ref})'});
            xlabel(p_g,'Iterations');
            legend(p_g,obj.node_legend);
            if strcmp(obj.name, 'step_load')
                x_r = 13; y_r = 0.3; w_r = 6; h_r = 0.7;
                rectangle(p_g,'Position', [x_r-w_r/2, y_r-h_r/2, w_r, h_r], ...
                    'EdgeColor', [0.4, 0.1, 0.4], 'LineWidth',1);

                %patch(f,[2200 3500 3500 2200],[-0.04 -0.04 -0.01 -0.01],[1 1 1],'EdgeColor','none');
                x_a = 0.3; y_a = 0.2; w_a = 0.13; h_a = 0.13;
                ax = axes('Units', 'Normalized', ...
                    'Position', [x_a, y_a, w_a, h_a], ...
                    'Box', 'on', ...
                    'LineWidth', 1, ...
                    'Color', [0.95, 0.99, 0.95]);
                
                hold(ax,'on')
                plot(ax, obj.x_save(2*n+1:3*n,obj.start_disp:obj.iteration)');
                ax.ColorOrderIndex = 1;
                plot(ax, obj.x_save(4*n+1:5*n,obj.start_disp:obj.iteration)','--');
                for k=(2*n+1):(3*n)
                    x_values = 1:(obj.iteration - (obj.start_disp - 1));                            %shift the index according to the
                    x_values = x_values(obj.p_g_limit_reached(k,obj.start_disp:obj.iteration));       %take out all the points where the limit was not reached
                    y_values = obj.x_save(k,obj.start_disp:obj.iteration);                          %trim the voltage values down to the desired length
                    y_values = y_values(obj.p_g_limit_reached(k,obj.start_disp:obj.iteration));      %take out all the points where the limit was not reached
                    plot(ax,x_values,y_values,'xr','LineWidth',2);
                end
                title('zoomed')
                axis([x_r-w_r/2, x_r+w_r/2, y_r-h_r/2, y_r+h_r/2]);
                hold(p_g,'off')
            else
                
            end
            
            hold(p_g,'off')
                 
            hold(q_g,'on')
            plot(q_g, obj.x_save(3*n+1:4*n,obj.start_disp:obj.iteration)');
            for k=(3*n+1):(4*obj.n)
                x_values = 1:(obj.iteration - (obj.start_disp - 1));                            %shift the index according to the
                x_values = x_values(obj.q_g_limit_reached(k,obj.start_disp:obj.iteration));       %take out all the points where the limit was not reached
                y_values = obj.x_save(k,obj.start_disp:obj.iteration);                          %trim the voltage values down to the desired length
                y_values = y_values(obj.q_g_limit_reached(k,obj.start_disp:obj.iteration));      %take out all the points where the limit was not reached
                plot(q_g,x_values,y_values,'xr','LineWidth',2);
            end
            ylabel(q_g,'q_g: Reactive Power generated [p.u]');
            xlabel(q_g,'Iterations');
            legend(q_g,obj.node_legend);
            hold(q_g,'off')
            
            f3 = figure(3);
            clf
            f3.OuterPosition = [800 0 800 1000];
            f = subplot(mp,np,2);
            p_ref = subplot(mp,np,1);
            f_function = subplot(mp,np,3);
            h_function = subplot(mp,np,4);
            box(f,'on')
            box(p_ref,'on')
            box(f_function,'on')
            box(h_function,'on')
                     
            hold(f,'on')
            plot(f, obj.x_save(end,obj.start_disp:obj.iteration)');
            for k=(5*n+4*m+1):(5*n+4*m+1)
                x_values = 1:(obj.iteration - (obj.start_disp - 1));                            %shift the index according to the
                x_values = x_values(obj.f_limit_reached(k,obj.start_disp:obj.iteration));       %take out all the points where the limit was not reached
                y_values = obj.x_save(k,obj.start_disp:obj.iteration);                          %trim the voltage values down to the desired length
                y_values = y_values(obj.f_limit_reached(k,obj.start_disp:obj.iteration));      %take out all the points where the limit was not reached
                plot(f,x_values,y_values,'xr','LineWidth',2);
            end
            ylabel(f,'f: Frequency deviation [Hz]');
            xlabel(f,'Iterations');
            
            % ZOOM: Specify the position and the size of the rectangle
            if strcmp(obj.name, 'wind')
                x_r = 1173; y_r = -0.012; w_r = 10; h_r = 0.04;
                rectangle(f,'Position', [x_r-w_r/2, y_r-h_r/2, w_r, h_r], ...
                    'EdgeColor', [0.4, 0.1, 0.4], 'LineWidth',1);

                patch(f,[2200 3500 3500 2200],[-0.04 -0.04 -0.01 -0.01],[1 1 1],'EdgeColor','none');
                x_a = 0.79; y_a = 0.63; w_a = 0.1; h_a = 0.15;
                ax = axes('Units', 'Normalized', ...
                    'Position', [x_a, y_a, w_a, h_a], ...
                    'Box', 'on', ...
                    'LineWidth', 1, ...
                    'Color', [0.95, 0.99, 0.95]);
                hold(ax,'on')
                plot(ax,obj.x_save(end,obj.start_disp:obj.iteration)');
                for k=(5*n+4*m+1):(5*n+4*m+1)
                    x_values = 1:(obj.iteration - (obj.start_disp - 1));                            %shift the index according to the
                    x_values = x_values(obj.f_limit_reached(k,obj.start_disp:obj.iteration));       %take out all the points where the limit was not reached
                    y_values = obj.x_save(k,obj.start_disp:obj.iteration);                          %trim the voltage values down to the desired length
                    y_values = y_values(obj.f_limit_reached(k,obj.start_disp:obj.iteration));      %take out all the points where the limit was not reached
                    plot(ax,x_values,y_values,'xr','LineWidth',2);
                end
                title('zoomed')
                axis([x_r-w_r/2, x_r+w_r/2, y_r-h_r/2, y_r+h_r/2]);
                hold(f,'off')
            else
                
            end

            box on
            plot(p_ref, obj.x_save(4*n+1:5*n,obj.start_disp:obj.iteration)');
            ylabel(p_ref,'p_{ref}: Reference value for prim. frequency controller [p.u]');
            xlabel(p_ref,'Iterations');
            legend(p_ref,obj.node_legend);
            
            box on            
            plot(f_function, obj.f_function_save(:,obj.start_disp:obj.iteration)');
            ylabel(f_function,'f(x): cost function value');
            xlabel(f_function,'Iterations');
            
%             figure(4);
%             new_handle = copyobj(v,figure(4));
%             subplot(1,1,1,new_handle);
%             saveas(figure(4),'volatile_wind1_f','epsc')
            
            box on            
            plot(h_function, obj.h_function_save(:,obj.start_disp:obj.iteration)');
            ylabel(h_function,'h(x): distance to the manifold');
            xlabel(h_function,'Iterations');
            
            %figure(1)
            
            %plot(sum(obj.x_save(3*n+1:4*n,obj.start_disp:obj.iteration)))
            
        end
        
    end
    
end

