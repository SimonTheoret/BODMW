function [exponent] = getBase2Exponent(value)

    exponent = 0;

    while value >= 2^exponent
        exponent = exponent + 1;
    end
    
    exponent = exponent - 1;

end