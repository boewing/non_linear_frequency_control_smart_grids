classdef Recorder < handle
    %Recorder: stores the state trajectories and plots them
    %This class can store state trajectories, convert them to i.e. voltage
    %or current trajectories and plot them according to the wishes.
    properties
        x_save
        v_limit_reached
        i_limit_reached
        f_limit_reached
        p_g_limit_reached
        q_g_limit_reached
        f_function_save
        h_function_save
        iteration               = 0
        transient_length        = 100
        start_disp              = 1
        m
        n
        node_legend
        line_legend
    end
    methods
        function obj = Recorder(mystate,f_function,h_function, iterations)
            x = mystate.getx();
            obj.m = length(mystate.i)/2;
            obj.n = length(mystate.v);
            obj.x_save = [x,zeros(length(x),iterations)];
            obj.f_function_save = [f_function,zeros(1,iterations)];
            obj.h_function_save = [h_function,zeros(1,iterations)];
            assert(obj.transient_length >= 1);
            assert(5*obj.n + 2*obj.m + 1 == length(x));
            
            obj.v_limit_reached = false(length(x),iterations);
            obj.i_limit_reached = false(length(x),iterations);
            obj.f_limit_reached = false(length(x),iterations);
            obj.p_g_limit_reached = false(length(x),iterations);
            obj.q_g_limit_reached = false(length(x),iterations);
            
            for i=1:obj.n
                obj.node_legend{i} = ['Node ', num2str(i)];
            end
            
            for i=1:obj.m
                obj.line_legend{i} = ['Line ', num2str(i)];
                obj.line_legend{i+obj.m} = ['Line (opp. dir.)', num2str(i)];
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
            obj.p_g_limit_reached(:,obj.iteration) = logical(Controller.S_limit_reached(mystate, mygrid));
            obj.q_g_limit_reached(:,obj.iteration) = logical(Controller.S_limit_reached(mystate, mygrid));
        end
        function plotV(obj)
            figure(1);
            plot(gca, obj.x_save(1:obj.n,1:obj.iteration)');
            ylabel(gca,'Voltage [p.u.]');
            xlabel(gca,'Iterations');
        end
        function plotAll(obj)
            f2 = figure(2);
            clf
            f2.Position = [0 0 1900 1000];
            mp = 2;
            np = 4;
            m = obj.m;
            n = obj.n;
            v = subplot(mp,np,1);
            i = subplot(mp,np,2);
            f = subplot(mp,np,3);
            p_g = subplot(mp,np,5);
            q_g = subplot(mp,np,6);
            p_ref = subplot(mp,np,7);
            f_function = subplot(mp,np,4);
            h_function = subplot(mp,np,8);
            
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
            legend(v,obj.node_legend);
            hold(v,'off')
            
            hold(i,'on')
            plot(i, obj.x_save(5*n+1:5*n+2*m,obj.start_disp:obj.iteration)');
            for k=(5*n+1):(5*n+2*m)
                x_values = 1:(obj.iteration - (obj.start_disp - 1));                            %shift the index according to the
                x_values = x_values(obj.i_limit_reached(k,obj.start_disp:obj.iteration));       %take out all the points where the limit was not reached
                y_values = obj.x_save(k,obj.start_disp:obj.iteration);                          %trim the voltage values down to the desired length
                y_values = y_values(obj.i_limit_reached(k,obj.start_disp:obj.iteration));      %take out all the points where the limit was not reached
                plot(i,x_values,y_values,'xr','LineWidth',2);
            end
            ylabel(i,'i: Current amplitude [p.u.]');
            xlabel(i,'Iterations');
            legend(i,obj.line_legend);
            hold(i,'off')
            
            hold(f,'on')
            plot(f, obj.x_save(end,obj.start_disp:obj.iteration)');
            for k=(5*n+2*m+1):(5*n+2*m+1)
                x_values = 1:(obj.iteration - (obj.start_disp - 1));                            %shift the index according to the
                x_values = x_values(obj.f_limit_reached(k,obj.start_disp:obj.iteration));       %take out all the points where the limit was not reached
                y_values = obj.x_save(k,obj.start_disp:obj.iteration);                          %trim the voltage values down to the desired length
                y_values = y_values(obj.f_limit_reached(k,obj.start_disp:obj.iteration));      %take out all the points where the limit was not reached
                plot(f,x_values,y_values,'xr','LineWidth',2);
            end
            ylabel(f,'f: Frequency deviation [Hz]');
            xlabel(f,'Iterations');
            hold(f,'off')
            
            hold(p_g,'on')
            plot(p_g, obj.x_save(2*n+1:3*n,obj.start_disp:obj.iteration)');
            for k=(2*n+1):(3*n)
                x_values = 1:(obj.iteration - (obj.start_disp - 1));                            %shift the index according to the
                x_values = x_values(obj.p_g_limit_reached(k,obj.start_disp:obj.iteration));       %take out all the points where the limit was not reached
                y_values = obj.x_save(k,obj.start_disp:obj.iteration);                          %trim the voltage values down to the desired length
                y_values = y_values(obj.p_g_limit_reached(k,obj.start_disp:obj.iteration));      %take out all the points where the limit was not reached
                plot(p_g,x_values,y_values,'xr','LineWidth',2);
            end
            ylabel(p_g,'p_g: Active Power generated [p.u]');
            xlabel(p_g,'Iterations');
            legend(p_g,obj.node_legend);
            hold(p_g,'off')
            
            hold(q_g,'on')
            plot(q_g, obj.x_save(3*n+1:4*n,obj.start_disp:obj.iteration)');
            for k=(3*n+1):(4*n)
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
            
            plot(p_ref, obj.x_save(4*n+1:5*n,obj.start_disp:obj.iteration)');
            ylabel(p_ref,'p_{ref}: Reference value for prim. frequency controller [p.u]');
            xlabel(p_ref,'Iterations');
            legend(p_ref,obj.node_legend);
            
            plot(f_function, obj.f_function_save(:,obj.start_disp:obj.iteration)');
            ylabel(f_function,'f(x): cost function value');
            xlabel(f_function,'Iterations');
            
            plot(h_function, obj.h_function_save(:,obj.start_disp:obj.iteration)');
            ylabel(h_function,'h(x): staying in the physical valid space');
            xlabel(h_function,'Iterations');
            
        end
        
    end
    
end

