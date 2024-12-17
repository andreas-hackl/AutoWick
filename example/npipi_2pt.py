import autowick as w

print("Automatic Wick contraction of $< p(x) pi+(z) pi-(w) pi-(q)^\dagger pi+(s)^\dagger p^\dagger(y) >$\n\n")
corr = w.CorrelationFunction([w.Proton('x'), w.PionPlus('z', 'D'), w.PionMinus('w', 'E'), w.PionPlus('q', 'F'), w.PionMinus('s', 'G'), w.AntiProton('y')])
print(corr)

# autogenerate distillation code
print("Distillation code:")
dist = corr.distillation({'x': 't', 'y': 'self.t0', 'z': 't', 'w': 't', 'q': 'self.t0', 's': 'self.t0'}, {'z': 'pi+/snk', 'w': 'pi-/snk', 'q': 'pi-/src', 's': 'pi+/src'})
print(dist)
