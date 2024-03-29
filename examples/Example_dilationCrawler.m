% 10-Apr-2023 11:50:28
% Auto-generated test script

% Initialize the test suite
% Add test cases here

clf;
w = 10;
s = sRectangle(0,50,0,50);
msh = Mesh(s,'Quads',[5,5]);
msh = msh.generate();

msh = msh.removeElements([11,12,1,21,14,13]);

fem = Fem(msh,'TimeStep',1/250,'TimeHorizon',1.25);

mat = NeoHookean(1.0,0.3);
mat.contact.ContactFriction = .1;

fem = fem.addMaterial(mat);
fem = fem.addContact(SDF);
fem = fem.addGravity();
fem = fem.addDilation([1:5], @(t) -0.5 * sin(w*pi*t) * smoothstep(6*t));
fem = fem.addDilation([11,16:19],@(t) 0.5 * sin(w*pi*t - 1.55*pi) * smoothstep(6*t));

fem.options.Display = @plt;
fem.solver.MaxIteration = 30;
fem.options.isNonlinear = 1;

fem = fem.simulate();

function y = SDF
    y = sLine(150,-15,-1,-1);
 end

 function plt(Fem)
    cla;
    showDilationFem(Fem);
    showContactFem(Fem);
    axis([-15 150 -5 60]);
    axis off;

    % if Fem.solver.Time <= Fem.solver.TimeStep
    %     gif('running.gif','DelayTime',1/60);
    %     disp('making GIF!');
    % else
    %     gif;
    % end
end