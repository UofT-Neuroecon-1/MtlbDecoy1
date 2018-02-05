for obs = 1:numel(Xs)
    for part = 1:numel(Particles{4}.OptimTheta)
        ProbaChoice( Xs{obs} , Particles{4}.model , Particles{4}.OptimTheta{part}, attrSign, param );
    end
end