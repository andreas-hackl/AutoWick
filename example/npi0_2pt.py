import autowick as w

print("Automatic Wick contraction of $< p(x) pi0(z) pi0(w) p(y) >$\n\n")
corr_uu = w.CorrelationFunction([w.Proton('x'), w.Pion0u('z', 'D'), w.Pion0u('w', 'E'), w.AntiProton('y')])
corr_ud = w.CorrelationFunction([w.Proton('x'), w.Pion0u('z', 'D'), w.Pion0d('w', 'E'), w.AntiProton('y')])
corr_du = w.CorrelationFunction([w.Proton('x'), w.Pion0d('z', 'D'), w.Pion0u('w', 'E'), w.AntiProton('y')])
corr_dd = w.CorrelationFunction([w.Proton('x'), w.Pion0d('z', 'D'), w.Pion0d('w', 'E'), w.AntiProton('y')])
corr = corr_uu + corr_dd - corr_ud - corr_du

print(corr)

print(len(corr.abstract_contractions))

# autogenerate distillation code
print("Distillation code:")
dist = corr.distillation({'x': 't', 'y': 'self.t0', 'z': 't', 'w': 'self.t0'}, {'z': 'pi0/snk', 'w': 'pi0/src'})
print(dist)
