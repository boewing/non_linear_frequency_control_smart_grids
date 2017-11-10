classdef Recorder < handle
    %Recorder: stores the state trajectories and plots them
    %This class can store state trajectories, convert them to i.e. voltage
    %or current trajectories and plot them according to the wishes.
    properties
        x_save
        iteration
        m
        n
        node_legend
        line_legend
    end
    methods
        function obj = Recorder(mystate, iterations)
            x = mystate.getx();
            obj.m = length(mystate.i)/2;
            obj.n = length(mystate.v);
            obj.x_save = [x,zeros(length(x),iterations)];
            obj.iteration = 1;
            
            for i=1:obj.n
                obj.node_legend{i} = ['Node ', num2str(i)];
            end
            
            for i=1:obj.m
                obj.line_legend{i} = ['Line ', num2str(i)];
                obj.line_legend{i+obj.m} = ['Line (opp. dir.)', num2str(i)];
            end
        end
        
%         function store(obj, mystate, i)
%             x_save(i) = mystate.getx();
%         end
        
        function store(obj, mystate)
            obj.iteration = obj.iteration + 1;
            obj.x_save(:,obj.iteration) = mystate.getx();
        end
        function plotV(obj)
            figure(1);
            plot(gca, obj.x_save(1:obj.n,1:obj.iteration)');
            ylabel(gca,'Voltage [p.u.]');
            xlabel(gca,'Iterations');
        end
        function plotAll(obj)
            f2 = figure(2);
            f2.Position = [0 0 1900 1000];
            mp = 2;
            np = 3;
            m = obj.m;
            n = obj.n;
            v = subplot(mp,np,1);
            i = subplot(mp,np,2);
            f = subplot(mp,np,3);
            p_g = subplot(mp,np,4);
            q_g = subplot(mp,np,5);
            p_ref = subplot(mp,np,6);
            
            plot(v,obj.x_save(1:obj.n,1:obj.iteration)');
            ylabel(v,'v: Voltage amplitude [p.u.]');
            xlabel(v,'Iterations');
            legend(v,obj.node_legend);
            
            plot(i, obj.x_save(5*n+1:5*n+2*m,1:obj.iteration)');
            ylabel(i,'i: Current amplitude [p.u.]');
            xlabel(i,'Iterations');
            legend(i,obj.line_legend);
            
            plot(f, obj.x_save(end,1:obj.iteration)');
            ylabel(f,'f: Frequency deviation [Hz]');
            xlabel(f,'Iterations');
            
            plot(p_g, obj.x_save(2*n+1:3*n,1:obj.iteration)');
            ylabel(p_g,'p_g: Active Power generated [p.u]');
            xlabel(p_g,'Iterations');
            legend(p_g,obj.node_legend);
            
            plot(q_g, obj.x_save(3*n+1:4*n,1:obj.iteration)');
            ylabel(q_g,'q_g: Reactive Power generated [p.u]');
            xlabel(q_g,'Iterations');
            legend(q_g,obj.node_legend);
            
            plot(p_ref, obj.x_save(3*n+1:4*n,1:obj.iteration)');
            ylabel(p_ref,'p_{ref}: Reference value for prim. frequency controller [p.u]');
            xlabel(p_ref,'Iterations');
            legend(p_ref,obj.node_legend);
        end
    end
    
end

