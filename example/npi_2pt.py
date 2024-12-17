import autowick as w

print("Automatic Wick contraction of $< n(x) pi+(z) pi-(w) n^\dagger(y) >$\n\n")
corr = w.CorrelationFunction([w.Neutron('x'), w.PionPlus('z', 'D'), w.PionMinus('w', 'E'), w.AntiNeutron('y')])
print(corr)

# autogenerate distillation code
print("Distillation code:")
dist = corr.distillation({'x': 't', 'y': 'self.t0', 'z': 't', 'w': 'self.t0'}, 
                         {'z': 'pi+/snk', 'w': 'pi+/src'})
print(dist)