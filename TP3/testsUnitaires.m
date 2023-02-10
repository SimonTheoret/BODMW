% Tests unitaires pour AM4
classdef testsUnitaires < matlab.unittest.TestCase
    properties (TestParameter)
        temps = 0;
        position = 0;
    end

    methods (Test)
        %Insérer tests unitaires ici:
        function testFunctionPredicteur(testCase)
            solution = AM4(odefun, testCase.temps, testCase.position, testCase.position, testCase.position, testCase.position, 1);
            vraiesolution = 0;
            test.Case.verifyEqual(solution,vraiesolution)
        end
    end
end


