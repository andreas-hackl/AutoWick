import numpy as np

tab=4*' '

class Fermion:
    def __init__(self, ftype, pos, spin_idx, color_idx):
        self.ftype = ftype.lower()
        self.cidx = color_idx
        self.sidx = spin_idx
        self.pos = pos

    def __str__(self):
        return '%s_{%s,%s}(%s)'%(self.ftype, self.sidx, self.cidx, self.pos)

    def __eq__(self, other):
        x1 = self.ftype == other.ftype
        x2 = self.cidx  == other.cidx
        x3 = self.sidx  == other.sidx
        x4 = self.pos   == other.pos
        return bool(np.prod([x1,x2,x3,x4]))
    

class AntiFermion(Fermion):
    def __str__(self):
        return '~'+Fermion.__str__(self)
    

class AbstractPropagator:
    def __init__(self, fermion, pos):
        self.fermion = 'l'
        self.pos = pos
    def __str__(self):
        return "L(%s|%s)"%self.pos


class Propagator:

    def __init__(self, fermion, afermion):
        self.quark = fermion.ftype.upper()
        self.pos = (fermion.pos, afermion.pos)
        self.spin = (fermion.sidx, afermion.sidx)
        self.color = (fermion.cidx, afermion.cidx)

    def __str__(self):
        return r"\prop{%s}{%s %s}{%s %s}(%s|%s)"%(self.quark, *self.spin, *self.color, *self.pos)
    
    def short_str(self):
        return r"%s(%s|%s)"%(self.quark, *self.pos)
    
    def abstract(self):
        return AbstractPropagator(self.quark, self.pos)
    
    def path(self):
        return list(self.pos)
    

class Sum:

    def __init__(self):
        self.terms = []

    def __add__(self, other):
        res = Sum()
        if type(other) == Sum:
            res.terms = self.terms + other.terms
            return res
        if type(other) == Product:
            res.terms = self.terms
            res.terms.append(other)
            return res

        raise ValueError()

    def __str__(self):
        s = self.terms[0].__str__()
        for t in self.terms[1:]:
            ts = t.__str__()
            if ts.startswith('-'):
                s += ' ' + ts
            else:
                s += ' + '+ts
        return s
    
    def short_str(self):
        return ' '.join([t.short_string() for t in self.terms])
    
    def abstract(self):
        res = Sum()
        res.terms = [term.abstract() for term in self.terms]
        return res

    

class Product:

    def __init__(self):
        self.operator = []
        self.prefactor = 1

    def add(self, op):
        self.operator.append(op)

    def __str__(self):
        if self.prefactor == 1:
            return ' '.join([o.__str__() for o in self.operator])
        elif self.prefactor == -1:
            return f'{self.prefactor:+d}'[0]+' '+' '.join([o.__str__() for o in self.operator])
        return f'{self.prefactor:+d}'+' '+' '.join([o.__str__() for o in self.operator])
    
    def short_str(self):
        if self.prefactor == 1:
            return ' '.join([o.short_str() for o in self.operator])
        elif self.prefactor == -1:
            return f'{self.prefactor:+d}'[0]+' '+' '.join([o.short_str() for o in self.operator])
        return f'{self.prefactor:+d}'+' '+' '.join([o.short_str() for o in self.operator])

    def __eq__(self, other):
        return self.__str__() == other.__str__()

    def __hash__(self):
        return hash(self.__str__())

    def __mul__(self, other):
        if isinstance(other, (int, float)):
            res = Product()
            for o in self.operator:
                res.add(o)
            res.prefactor = self.prefactor * other
            return res
            
        res = Product()
        for o in self.operator:
            res.add(o)
        for o in other.operator:
            res.add(o)
        res.prefactor = self.prefactor * other.prefactor
        return res

    def __add__(self, other):
        res = Sum()
        res.terms = [self, other]
        return res
    
    def abstract(self):
        res = Product()

        res.prefactor = self.prefactor
        for op in self.operator:
            res.add(op.abstract())

        return res
    

class AbstractMatrix:
    def __init__(self, symb):
        self.symb = symb
    def __str__(self):
        return r"%s"%(self.symb)
    

class SpinMatrix:

    def __init__(self, symb, spin_idx):
        self.symb = symb
        self.spin = spin_idx
    
    def __str__(self):
        return r"\smat{%s}{%s %s}"%(self.symb, *self.spin)
    
    def short_str(self):
        return r"%s"%(self.symb)
    
    def abstract(self):
        return AbstractMatrix(self.symb)
    

class ColorMatrix:

    def __init__(self, symb, color_idx):
        self.symb = symb
        self.color = color_idx

    def __str__(self):
        return r"\cmat{%s}{%s %s}"%(self.symb, *self.color)
    
    def short_str(self):
        return r"%s"%(self.symb)
    
    def abstract(self):
        return  AbstractMatrix(self.symb)
    
class EpsTensor:
    def __init__(self, a, b, c):
        self.a = a
        self.b = b
        self.c = c

class Baryon:
    def __init__(self, ftypes, pos):
        self.ftypes = ftypes
        self.pos = pos

        self.fermions = [
            Fermion(ftypes[0], pos, r"\alpha", r"a"),
            Fermion(ftypes[1], pos, r"\beta", r"b"),
            Fermion(ftypes[2], pos, r"\gamma", r"c"),
        ]

        self.spin = r"\alpha"
        self.eps = EpsTensor('a', 'b', 'c')

        self.Gamma = SpinMatrix(r"\Gamma", (r"\beta", r"\gamma"))


class AntiBaryon:
    def __init__(self, ftypes, pos):
        self.ftypes = ftypes
        self.pos = pos

        self.fermions = [
            AntiFermion(ftypes[2], pos, r"\gamma'", r"c'"),
            AntiFermion(ftypes[1], pos, r"\beta'", r"b'"),
            AntiFermion(ftypes[0], pos, r"\alpha'", r"a'"),
        ]

        self.spin = r"\alpha'"
        self.eps = EpsTensor(r"a'", r"b'", r"c'")

        self.Gamma = SpinMatrix(r"\bar\Gamma", (r"\beta'", r"\gamma'"))


class Meson:
    def __init__(self, ftypes, pos, p, spin, color, spin_mat):
        self.ftypes = ftypes
        self.pos = pos
        self.momentum = p

        self.fermions = [
            AntiFermion(ftypes[0], pos, r"%s"%spin, r"%c"%color),
            Fermion(ftypes[1], pos, r"%s'"%spin, r"%s'"%color),
        ]

        self.Gamma = spin_mat
        self.mom_ins = ColorMatrix(r'T(%s)'%p, (r'%s'%color, r"%s'"%color))


alphabet = {
    'D': (r"\delta", r"d"),
    'E': (r"\varepsilon", r"e"),
    'F': (r"\zeta", r"f"),
    'G': (r"\eta", r"g"),
    'H': (r"\theta", r"h"),
    'I': (r"\iota", r"i"),
    'K': (r"\kappa", r"k"),
    'L': (r"\lambda", r"l"),
}


class PionPlus(Meson):
    def __init__(self, pos, letter, p='0'):
        spin = alphabet[letter][0]
        color = alphabet[letter][1]
        Meson.__init__(self, ['d', 'u'], pos, p, spin, color, SpinMatrix(r'\gamma^5', (r"%s"%spin, r"%s'"%spin)))
        self.letter = letter

    def name(self):
        return "pi_plus"

class PionMinus(Meson):
    def __init__(self, pos, letter, p='0'):
        spin = alphabet[letter][0]
        color = alphabet[letter][1]
        Meson.__init__(self, ['u', 'd'], pos, p, spin, color, SpinMatrix(r'\gamma^5', (r"%s"%spin, r"%s'"%spin)))
        self.letter = letter

    def name(self):
        return "pi_minus"

class Pion0d(Meson):
    def __init__(self, pos, letter, p='0'):
        spin = alphabet[letter][0]
        color = alphabet[letter][1]
        Meson.__init__(self, ['d', 'd'], pos, p, spin, color, SpinMatrix(r'\gamma^5', (r"%s"%spin, r"%s'"%spin)))
        self.letter = letter

    def name(self):
        return "pi_0"

class Pion0u(Meson):
    def __init__(self, pos, letter, p='0'):
        spin = alphabet[letter][0]
        color = alphabet[letter][1]
        Meson.__init__(self, ['u', 'u'], pos, p, spin, color, SpinMatrix(r'\gamma^5', (r"%s"%spin, r"%s'"%spin)))
        self.letter = letter
    
    def name(self):
        return "pi_0"
    

class LocalCurrent:
    def __init__(self, ftypes, pos, p, spin, color, spin_mat):
        self.ftypes = ftypes
        self.pos = pos
        self.momentum = p

        self.fermions = [
            AntiFermion(ftypes[0], pos, r"%s"%spin, r"%c"%color),
            Fermion(ftypes[1], pos, r"%s'"%spin, r"%s'"%color),
        ]

        self.Gamma = spin_mat
        self.mom_ins = ColorMatrix(r'T(%s)'%p, (r'%s'%color, r"%s'"%color))


class LocalCurrentMinus(LocalCurrent):
    def __init__(self, pos, letter, p='0'):
        spin = alphabet[letter][0]
        color = alphabet[letter][1]
        Meson.__init__(self, ['u', 'd'], pos, p, spin, color, SpinMatrix(r'\Gamma', (r"%s"%spin, r"%s'"%spin)))
        self.letter = letter

    def name(self):
        return "J_minus"
    
class LocalCurrentPlus(LocalCurrent):
    def __init__(self, pos, letter, p='0'):
        spin = alphabet[letter][0]
        color = alphabet[letter][1]
        Meson.__init__(self, ['d', 'u'], pos, p, spin, color, SpinMatrix(r'\Gamma', (r"%s"%spin, r"%s'"%spin)))
        self.letter = letter

    def name(self):
        return "J_plus"
    
class LocalCurrent0u(LocalCurrent):
    def __init__(self, pos, letter, p='0'):
        spin = alphabet[letter][0]
        color = alphabet[letter][1]
        Meson.__init__(self, ['u', 'u'], pos, p, spin, color, SpinMatrix(r'\Gamma', (r"%s"%spin, r"%s'"%spin)))
        self.letter = letter

    def name(self):
        return "J_0"


class LocalCurrent0d(LocalCurrent):
    def __init__(self, pos, letter, p='0'):
        spin = alphabet[letter][0]
        color = alphabet[letter][1]
        Meson.__init__(self, ['d', 'd'], pos, p, spin, color, SpinMatrix(r'\Gamma', (r"%s"%spin, r"%s'"%spin)))
        self.letter = letter

    def name(self):
        return "J_0"
    
class Proton(Baryon):
    def __init__(self, pos):
        Baryon.__init__(self, ['u', 'u', 'd'], pos)

    def symb(self):
        return 'p'
    
class Neutron(Baryon):
    def __init__(self, pos):
        Baryon.__init__(self, ['d', 'd', 'u'], pos)
    def symb(self):
        return 'n'

class AntiProton(AntiBaryon):
    def __init__(self, pos):
        AntiBaryon.__init__(self, ['u', 'u', 'd'], pos)
    def symb(self):
        return 'p'
    
class AntiNeutron(AntiBaryon):
    def __init__(self, pos):
        AntiBaryon.__init__(self, ['d', 'd', 'u'], pos)
    def symb(self):
        return 'n'

class AbstractMesonPath:
    def __init__(self, props, mesons):
        self.positions = [p.pos[0] for p in props] + [props[-1].pos[1]]
        self.quarks = [p.quark for p in props]
        self.momenta = [m.momentum for m in mesons]
        self.pos = [self.positions[0], self.positions[-1]]

    def is_loop(self):
        return self.positions[0] == self.positions[-1]
    
    def __str__(self):
        s = ''
        s += f'L(%s|%s) '%(self.positions[0], self.positions[1])
        for i in range(len(self.momenta)):
            s += f'T(%s) gamma_pi '%(self.momenta[i])
            if len(self.quarks) > i+1:
                s += f'L(%s|%s) '%(self.positions[i+1], self.positions[i+2])
        return s



class AbstractMesonLoop:
    def __init__(self, symb, path):
        self.symb = symb
        self.path = path

    def __str__(self):
        return "%s"%(self.symb)
    
    def definition(self):
        return r"%s = %s"%(self.symb, self.path.__str__())

    def distillation(self, pos_t, pion_names):
        header = f"{tab}def {self.symb.replace('S', 'seq')}(self, t): \n"

        pos = self.path.positions
        peramb_str = '[ ' + ', '.join([f'self.perambs[{pos_t[pos[i+1]]}][{pos_t[pos[i]]}]' for i in range(len(pos)-1)]) + ' ]'

        pion_str = '[' + ', '.join([f"self.pions[\"{pion_names[p]}\"]({pos_t[p]})" for p in pos[1:]]) + ' ]'
        
        line = f"{2*tab}return loop({peramb_str}, {pion_str})\n"
        return header + line


class AbstractSeqProp:
    def __init__(self, symb, path):
        self.symb = symb
        self.path = path
    def __str__(self):
        return "%s(%s|%s)"%(self.symb, *self.path.pos)
    def definition(self):
        return r"%s = %s"%(self.symb, self.path.__str__())

    def distillation(self, pos_t, pion_names):
        header = f"{tab}def {self.symb.replace('S', 'seq')}(self, t): \n"

        pos = self.path.positions
        peramb_str = '[ ' + ', '.join([f'self.perambs[{pos_t[pos[i+1]]}][{pos_t[pos[i]]}]' for i in range(len(pos)-1)]) + ' ]'

        pion_str = '[' + ', '.join([f"self.pions[\"{pion_names[p]}\"]({pos_t[p]})" for p in pos[1:-1]]) + ' ]'
        
        line = f"{2*tab}return sequential({peramb_str}, {pion_str})\n"
        return header + line

    

class SequentialProp:
    def __init__(self, symb, props, mesons):
        prod = Product()
        n = len(mesons)
        prod.add(props[0])
        for i in range(n):
            prod.add(mesons[i].mom_ins)
            prod.add(mesons[i].Gamma)
            if len(props) > n:
                prod.add(props[i+1])
        self.prod = prod

        self.symb = symb
        self.color = (props[0].color[0], props[-1].color[1])
        self.spin = (props[0].spin[0], props[-1].spin[1])
        self.pos = (props[0].pos[0], props[-1].pos[1])

        self.props = props
        self.mesons = mesons

    def __str__(self):
        return r"\prop{%s}{%s %s}{%s %s}(%s|%s)"%(self.symb, *self.spin, *self.color, *self.pos)
    
    def short_str(self):
        return r"%s(%s|%s)"%(self.symb, *self.pos)
    
    def definition(self):
        return self.prod.__str__()
    
    def __eq__(self, other):
        if self.pos != other.pos:
            return False
        
        if len(self.props) != len(other.props):
            return False
        for sp, op in zip(self.props, other.props):
            if sp.quark != op.quark:
                return False
            if sp.pos != op.pos:
                return False
        
        return True
    
    def abstract(self, new_symb=''):
        if new_symb == '':
            new_symb = self.symb
            
        path = AbstractMesonPath(self.props, self.mesons)

        if path.is_loop():
            return AbstractMesonLoop(new_symb, path)

        return AbstractSeqProp(new_symb, path)
    
    def path(self):
        mpath = AbstractMesonPath(self.props, self.mesons)
        return mpath.positions
    
def minimal_str(s):
    s = s.replace('(x|y)', '')
    s = s.replace(r'\bar\Gamma', 'G')
    s = s.replace(r'\Gamma', 'G')
    s = s.replace(' ', '')
    return s

class TracelessContribution:
    def __init__(self, prefactor, prod0, prod1, prod2, loops=[]):
        self.prefactor = prefactor
        self.prod0 = prod0
        self.prod1 = prod1
        self.prod2 = prod2
        self.loops = loops

    def __str__(self):
        base_str = f"{self.prefactor} * Traceless[{self.prod0.__str__()}, {self.prod1.__str__()}, {self.prod2.__str__()}]"
        if len(self.loops) == 0:
            return base_str
        return base_str + '*' + '*'.join([f"{loop.__str__()}" for loop in self.loops])

    def hash(self):
        s = f"{self.prefactor} [" + ", ".join([minimal_str(p.__str__()) for p in [self.prod0, self.prod1, self.prod2]]) + ' ]'
        if len(self.loops) > 0:
            ls = "[" + ", ".join([minimal_str(l.__str__()) for l in self.loops]) + ']'
            return (s, ls)
        return s

    def repr(self):
        s = f"Tl[" + ", ".join([minimal_str(p.__str__()) for p in [self.prod0, self.prod1, self.prod2]]) + ' ]'
        if len(self.loops) > 0:
            ls = "[" + ", ".join(sorted([minimal_str(l.__str__()) for l in self.loops])) + ']'
            return (s, ls)
        return s

    def copy(self):
        return TracelessContribution(self.prefactor, self.prod0, self.prod1, self.prod2, self.loops)

    
    def short_str(self):
        s = f"{self.prefactor} L[" + ", ".join([minimal_str(p.__str__()) for p in [self.prod0, self.prod1, self.prod2]]) + ' ]'
        if len(self.loops) > 0:
            ls = " " + " * ".join([minimal_str(l.__str__()) for l in self.loops])
            return s + ls
        return s

    def op_names(self):
        return [minimal_str(p.__str__()) for p in [self.prod0, self.prod1, self.prod2]]

    def equiv(self, other):
        if type(other) != TracelessContribution:
            return False
        return self.repr() == other.repr()
        
    
class TracefullContribution:
    def __init__(self, prefactor, prod0, prod1, prod2, loops=[]):
        self.prefactor = prefactor
        self.prod0 = prod0
        self.prod1 = prod1
        self.prod2 = prod2
        self.loops = loops

    def __str__(self):
        base_str = f"{self.prefactor} * Tracefull[{self.prod0.__str__()}, {self.prod1.__str__()}, {self.prod2.__str__()}]"
        if len(self.loops) == 0:
            return base_str
        return base_str + '*' + '*'.join([f"{loop.__str__()}" for loop in self.loops])

    def hash(self):
        s = f"{self.prefactor} [" + ", ".join([minimal_str(p.__str__()) for p in [self.prod0, self.prod1, self.prod2]]) + ' ]'
        if len(self.loops) > 0:
            ls = "[" + ", ".join([minimal_str(l.__str__()) for l in self.loops]) + ']'
            return (s, ls)
        return s
    
    def short_str(self):
        s = f"{self.prefactor} T[" + ", ".join([minimal_str(p.__str__()) for p in [self.prod0, self.prod1, self.prod2]]) + ' ]'
        if len(self.loops) > 0:
            ls = " " + " * ".join([minimal_str(l.__str__()) for l in self.loops])
            return s + ls
        return s

    def repr(self):
        s = f"Tf[" + ", ".join([minimal_str(p.__str__()) for p in [self.prod0, self.prod1, self.prod2]]) + ' ]'
        if len(self.loops) > 0:
            ls = "[" + ", ".join(sorted([minimal_str(l.__str__()) for l in self.loops])) + ']'
            return (s, ls)
        return s

    def copy(self):
        return TracefullContribution(self.prefactor, self.prod0, self.prod1, self.prod2, self.loops)

    def op_names(self):
        return [minimal_str(p.__str__()) for p in [self.prod0, self.prod1, self.prod2]]

    def equiv(self, other):
        if type(other) != TracefullContribution:
            return False
        return self.repr() == other.repr()