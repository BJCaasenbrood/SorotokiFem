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
    end
end