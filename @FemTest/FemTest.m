classdef FemTest < matlab.unittest.TestCase

    methods (Test)

        function testSolve(testCase)
            sdf = sRectangle(10, 5);
            msh = Mesh(sdf,'NElem',30);
            msh = msh.generate();

            fem = Fem(msh,'TimeStep',0.1);
            fem = fem.addMaterial( NeoHookean );
            fem = fem.addSupport('left', [1, 1]);
            fem = fem.addSupport('right',[1, 1]);
            fem = fem.addLoad('bottom', [0, 5e-4]);

            fem.options.Display = [];
            fem.solver.isLog = false;

            fem = fem.solve();
            testCase.verifyClass(fem,'Fem');
        end

        function testDisplace(testCase)
            sdf = sRectangle(10, 5);
            msh = Mesh(sdf,'NElem',30);
            msh = msh.generate();

            fem = Fem(msh,'TimeStep',0.1);
            fem = fem.addMaterial( NeoHookean );
            fem = fem.addSupport('left', [1, 1]);
            fem = fem.addDisplace('right',[5, 0]);

            fem.options.Display = [];
            fem.solver.isLog = false;

            fem = fem.solve();
            testCase.verifyNotEmpty(fem.solver.sol.tout);
            testCase.verifyNotEmpty(fem.solver.sol.yout);
        end

        function testGravity(testCase)
            sdf = sRectangle(10, 5);
            msh = Mesh(sdf,'NElem',30);
            msh = msh.generate();

            fem = Fem(msh,'TimeStep',0.1);
            fem = fem.addMaterial( NeoHookean );
            fem = fem.addSupport('left', [1, 1]);
            fem = fem.addGravity();

            fem.options.Display = [];
            fem.solver.isLog = false;

            fem = fem.solve();
            testCase.verifyNotEmpty(fem.solver.sol.tout);
            testCase.verifyNotEmpty(fem.solver.sol.yout);
        end

        function testPressure(testCase)
            sdf = sRectangle(10) - sCircle(3);
            msh = Mesh(sdf,'NElem',80);
            msh = msh.generate();

            fem = Fem(msh,'TimeStep',0.1);
            fem = fem.addMaterial( NeoHookean );
            fem = fem.addSupport('bottom', [1, 1]);
            fem = fem.addSupport('left', [1, 1]);
            fem = fem.addSupport('right', [1, 1]);
            fem = fem.addSupport('top', [1, 1]);
            fem = fem.addPressure('allhole', 1e-1);

            fem.options.Display = [];
            fem.solver.isLog = false;

            fem = fem.solve();
            testCase.verifyNotEmpty(fem.solver.sol.tout);
            testCase.verifyNotEmpty(fem.solver.sol.yout);
        end

        function testEigen(testCase)
            sdf = sRectangle(10, 5);
            msh = Mesh(sdf,'NElem',30);
            msh = msh.generate();

            fem = Fem(msh,'TimeStep',0.1);
            fem = fem.addMaterial( NeoHookean );
            fem = fem.addSupport('left', [1, 1]);
            fem = fem.addSupport('right',[1, 1]);
            fem = fem.addLoad('bottom', [0, 5e-4]);

            fem.options.Display = [];
            fem.solver.isLog = false;

            fem = fem.solve();
            testCase.verifyClass(fem,'Fem');

            fem = fem.eigen();
            testCase.verifyNotEmpty(fem.solver.sol.pod.V);
            testCase.verifyNotEmpty(fem.solver.sol.pod.D);
        end

        function testSimulate(testCase)
            sdf = sRectangle(10, 3);
            msh = Mesh(sdf,'Quads',[10,3]);
            msh = msh.generate();

            fem = Fem(msh,'TimeStep',0.1);
            fem = fem.addMaterial( NeoHookean );
            fem = fem.addSupport('left', [1, 1]);
            fem = fem.addGravity();

            fem.options.Display = [];
            fem.solver.isLog = false;

            fem = fem.simulate();
            testCase.verifyNotEmpty(fem.solver.sol.tout);
            testCase.verifyNotEmpty(fem.solver.sol.yout);
        end

        function testOptCompliance(testCase)
            msh = Mesh(sRectangle(10, 2),'Quads',[50,10]);
            msh = msh.generate();

            fem = Fem(msh,'SpatialFilterRadius',0.3);

            fem = fem.addMaterial(NeoHookean);
            fem = fem.addSupport('se',[1, 1]);
            fem = fem.addSupport('sw',[1, 1]);
            fem = fem.addLoad('bottommid',[0,1]);

            fem.solver.isLog = false;
            fem.options.isNonlinear = false;
            fem.topology.MaxIteration = 15;
            fem.options.Display = [];

            fem = fem.optimize;
            testCase.verifyNotEmpty(fem.topology.sol.x);
        end

        function testOptCompliant(testCase)
            msh = Mesh(sRectangle(10, 10),'Quads',[55,55]);
            msh = msh.generate();

            fem = Fem(msh,'SpatialFilterRadius',0.3);

            fem = fem.addMaterial(NeoHookean);
            fem = fem.addSupport('nw',[1, 1]);
            fem = fem.addSupport('sw',[1, 1]);
            fem = fem.addLoad('leftmid',[-1,0]);
            fem = fem.addOutput('rightmid',[1,0]);

            fem.topology.Type = 'compliant';
            fem.options.isNonlinear = false;
            fem.options.Display = [];
            fem.topology.MaxChange = 0.15;
            fem.topology.MaxIteration = 15;
            fem = fem.optimize;
            testCase.verifyNotEmpty(fem.topology.sol.x);
        end
    end
end
