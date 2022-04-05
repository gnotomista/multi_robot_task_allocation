classdef LTAMulti < LTASingle
    methods (Access = protected)
        function h = gamma_of_h_of_x(obj, i, j)
            switch j
                case 1 % coverage
                    J = norm(obj.optimization_data.p(:,i)-obj.optimization_data.G(:,i))^2;
                case 2 % static formation
                    J = norm(obj.optimization_data.p(:,i)-obj.optimization_data.y(:,i))^2;
            end
            h = -J;
        end
        function dh = dh_dx(obj, i, j)
            switch j
                case 1 % coverage
                    dJ = 2*(obj.optimization_data.p(:,i)-obj.optimization_data.G(:,i));
                case 2 % static formation
                    dJ = 2*(obj.optimization_data.p(:,i)-obj.optimization_data.y(:,i));
            end
            dh = -dJ;
        end
    end
end
