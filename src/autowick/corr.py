import numpy as np
import itertools as it
from autowick.objects import *

tab=4*' '


def permuation_sign(P):
    sign = 1
    for i in range(len(P)-1):
        if P[i] != i:
            sign *= -1
            idx = P.index(min(P[i:]))
            P[i], P[idx] = P[idx], P[i]
    return sign



def product(terms):
    if len(terms) == 1:
        return terms[0]

    n0 = len(terms[0])
    n1 = len(terms[1])

    prod_terms = []
    for i0, i1 in it.product(range(n0), range(n1)):
        prod_terms.append(terms[0][i0]*terms[1][i1])

    return product([prod_terms]+terms[2:])

def contract(flist):
    ftypes = sorted(list(set([f.ftype for f in flist])))
    character = [(-1)**(type(f)==Fermion) * (ftypes.index(f.ftype)+1) for f in flist]
    if sum(character) != 0:
        raise ValueError()
    
    nflist = flist.copy()
    n = len(flist)//2
    
    fermions = {ftype: [] for ftype in ftypes}
    afermions = {ftype: [] for ftype in ftypes}
    for f in flist:
        if type(f) == Fermion:
            fermions[f.ftype].append(f)
        else:
            afermions[f.ftype].append(f)
    
    total_sign = 1
    ordered_list = []
    remainder = nflist
    for ftype in ftypes:
        for fermion, afermion in zip(fermions[ftype], afermions[ftype]):
            # move fermion to the correct position
            i = remainder.index(fermion)
            sign_i = (-1)**(i)
            total_sign *= sign_i
            ordered_list.append(remainder[i])
            remainder.remove(remainder[i])

            # move antifermion to the correct position
            i = remainder.index(afermion)
            sign_i = (-1)**(i)
            total_sign *= sign_i
            ordered_list.append(remainder[i])
            remainder.remove(remainder[i])
            
    ####
    terms = []
    for ftype in ftypes:
        subf = fermions[ftype]
        suba = afermions[ftype]
        n = len(subf)
        subterms = []
        for P in it.permutations(range(n)):
            term = Product()
            term.prefactor = permuation_sign(list(P)) * (-1)**n
            for i, pi in enumerate(P):
                term.add(Propagator(subf[i], suba[pi]))
            
            subterms.append(term)
        terms.append(subterms)
    
    terms = product(terms)
    
    res = terms[0] * total_sign
    for term in terms[1:]:
        res = res + term * total_sign
    return res

def check_loop(path):
    if path[0] != path[-1]:
        return path
    else:
        if path[0] not in ['A', 'B', 'C']:
            return path[:-1]
        return path


def split_hash(h):
    i = h.index('[')
    return int(h[:i]), h[i:]

def cross_hash(h):
    i = h.index('[')
    spre = h[:i]
    s = h[i:]
    ssp = s[1:-1].split(', ')
    ssp = [ssp[1], ssp[0], ssp[2]]
    ns = ', '.join(ssp)
    return spre + '[' + ns + ']'



def in_contr(contr, seq):
    for prod in [contr.prod0, contr.prod1, contr.prod2] + contr.loops:
        for op in prod.operator:
            if type(op) == AbstractSeqProp or type(op) == AbstractMesonLoop :
                if op.definition() == seq.definition():
                    return True
    return False
    

def operation(op1, op2, op_func, unit=0):
    # compare sequential operator
    seqs = []
    for seq in op1.abstract_seqs:
        seqs.append(seq)

    n = len(seqs)
    correspondens = {}
    for o_seq in op2.abstract_seqs:
        done = False
        for i, seq in enumerate(seqs):
            if seq.definition().split("=")[1] == o_seq.definition().split("=")[1]:
                correspondens[o_seq.symb] = seq.symb
                done = True
        if not done:
            seqs.append(o_seq)
            nsym = f"S{n}"
            correspondens[o_seq.symb] = nsym
            n += 1


    for j, (osymb, nsymb) in enumerate(correspondens.items()):
        op2.abstract_seqs[j].symb = nsymb
        
            
    contr_dict_1 = {contr.repr(): contr.copy() for contr in op1.abstract_contractions}
    contr_dict_2 = {contr.repr(): contr.copy() for contr in op2.abstract_contractions}

    tags = list(set(list(contr_dict_1.keys()) + list(contr_dict_2.keys())))

    contrs = []
    for tag in tags:
        if tag in contr_dict_1:
            contr = contr_dict_1[tag].copy()
            if tag in contr_dict_2:
                p1 = contr_dict_1[tag].prefactor
                p2 = contr_dict_2[tag].prefactor

                pre = op_func(p1, p2)
                contr.prefactor = pre
            else:
                p1 = contr_dict_1[tag].prefactor
                pre = op_func(p1, unit)
                contr.prefactor = pre
            if abs(pre) > 1e-10:
                contrs.append(contr)
        elif tag in contr_dict_2:
            contr = contr_dict_2[tag].copy()
            p2 = contr_dict_2[tag].prefactor
            pre = op_func(unit, p2)
            contr.prefactor = pre
            contrs.append(contr)
        else:
            raise
        
    

    used_seqs = []
    for seq in seqs:
        done = False
        for contr in contrs:
            if in_contr(contr, seq) and not done:
                used_seqs.append(seq)
                done = True

    cf = CorrelationFunction()
    cf.abstract_contractions = contrs
    cf.abstract_seqs = used_seqs
    return cf


class CorrelationFunction:

    def __init__(self, particle_list=None):
        if particle_list is not None:
            self.particles = particle_list
            self.baryon = None
            self.antibaryon = None
            self.mesons = []
    
            for p in self.particles:
                if isinstance(p, Baryon):
                    if self.baryon is None:
                        self.baryon = p
                    else:
                        raise ValueError()
                elif isinstance(p, AntiBaryon):
                    if self.antibaryon is None:
                        self.antibaryon = p
                    else:
                        raise ValueError()
                    
                elif isinstance(p, Meson):
                    self.mesons.append(p)

                elif isinstance(p, LocalCurrent):
                    self.mesons.append(p)
    
                else:
                    raise ValueError(f'unknown type: {type(p)}')
                
            self.letters = [m.letter for m in self.mesons]
            if len(set(self.letters)) != len(self.letters):
                raise ValueError(f'at least one letter appears multiple times')
            
    
            cf = self.contract_fermions()
            ndiag = len(cf.terms)
    
            self.seqs = []
    
            self.contractions = []
    
            for i, prod in enumerate(cf.terms):
                contr = self.solve_diag(i, prod)
                self.contractions.append(contr)
    
            self.classified_seqs = self.classify_seqs()
    
    
            self.abstract_contractions = self.abstract()
        else:
            self.abstract_contractions = []
            self.abstract_seqs = []
      

    def contract_fermions(self):
        flist = list(it.chain(*[p.fermions for p in self.particles]))
        return contract(flist)
    
    def solve_diag(self, i, prod):
        sign = prod.prefactor
        in_letters = [op.color[0].upper()[0] for op in prod.operator]
        out_letters = [op.color[1].upper()[0] for op in prod.operator]
        paths = []
        unused_letters = self.letters.copy()
    
        for letter in self.letters:
            if letter in unused_letters:
                
    
                current_letter = letter
                i0 = out_letters.index(current_letter)
                if in_letters[i0] not in unused_letters:
                    finished = False
                    path = [i0]
                    while not finished:
                        inext = in_letters.index(current_letter)
                        
                        unused_letters.remove(current_letter)
                        current_letter = out_letters[inext]

                        path.append(inext)
                        if current_letter not in unused_letters:
                            finished = True

    
                    paths.append(path)
        
        if len(unused_letters) != 0:
            # Use case when there are pion loops

            for letter in self.letters:
                if letter in unused_letters:
                    
                    current_letter = letter
                    i0 = out_letters.index(current_letter)
                    start_letter = in_letters[i0]

                    if start_letter == current_letter:
                        path = [i0, i0]
                        finished=True
                        unused_letters.remove(current_letter)
                    else:
                        finished = False
                        path = [i0]
                        while not finished:
                            inext = in_letters.index(current_letter)
                            unused_letters.remove(current_letter)
                            current_letter = out_letters[inext]

                        
                            

                            path.append(inext)
                            if current_letter == start_letter :
                                path.append(i0)
                                finished = True
                                unused_letters.remove(current_letter)

    
                    paths.append(path)

        

        local_seq = []
        for j, path in enumerate(paths):
            props = [prod.operator[i] for i in check_loop(path)]
            
            mesons = [self.mesons[self.letters.index(out_letters[i])] for i in path[:-1]]
            seq = SequentialProp(f'S{i}_{j}', props, mesons)
            self.seqs.append(seq)
            local_seq.append(seq)

        new_op = {}
        for op in prod.operator + local_seq:
            l0 = op.color[0].upper()[0]
            l1 = op.color[1].upper()[0]
            if l0 not in self.letters and l1 not in self.letters:
                new_op[(l0, l1)] = op 
    
            elif l0 in self.letters and l1 in self.letters:
                new_op[(l0, l1)] = op

        new_op_start = {key[0]: op for key, op in new_op.items()}

        with_loop = len(new_op) > 3

        op0 = Product()
        op1 = Product()
        op2 = Product()
        loops = []

        A = new_op_start['A']
        B = new_op_start['B']
        C = new_op_start['C']
        G0 = self.baryon.Gamma
        G1 = self.antibaryon.Gamma

        if new_op.get(('A', 'A'), None) is not None:
            res_type = TracefullContribution

            if new_op.get(('B', 'B'), None) is not None:
                # + Q[BG, GC] A

                diag_sign = 1

                op0.add(B)
                op0.add(G1)

                op1.add(G0)
                op1.add(C)

                op2.add(A)

            else:
                # - Q[B, GCG] A

                diag_sign = -1

                op0.add(B)

                op1.add(G0)
                op1.add(C)
                op1.add(G1)

                op2.add(A)

        elif new_op.get(('A', 'B'), None) is not None:
            res_type = TracelessContribution

            if new_op.get(('B', 'A'), None) is not None:
                # - Q[AG, GC] B

                diag_sign = -1

                op0.add(A)
                op0.add(G1)

                op1.add(G0)
                op1.add(C)

                op2.add(B)

            else:
                # + Q[AG, B] GC

                diag_sign = 1

                op0.add(A)
                op0.add(G1)

                op1.add(B)

                op2.add(G0)
                op2.add(C)

        elif new_op.get(('A', 'C'), None) is not None:
            res_type = TracelessContribution

            if new_op.get(('B', 'A'), None) is not None:
                # + Q[A, GCG] B

                diag_sign = 1

                op0.add(A)

                op1.add(G0)
                op1.add(C)
                op1.add(G1)

                op2.add(B)

            else:
                # - Q[A, BG] GC

                diag_sign = -1

                op0.add(A)

                op1.add(B)
                op1.add(G1)

                op2.add(G0)
                op2.add(C)

        loops = []
        if with_loop:
            for letter in self.letters:
                if new_op.get((letter, letter), None) is not None:
                    opl = Product()
                    opl.add(new_op[(letter, letter)])
                    loops.append(opl)


        return res_type(sign*diag_sign, op0, op1, op2, loops)
    
    def classify_seqs(self):
        if len(self.seqs) == 0:
            return []
        classified_seq = [[self.seqs[0]]]

        for seq in self.seqs[1:]:
            equivalence_flag = False
            for i, seq_list in enumerate(classified_seq):
                if seq == seq_list[0]:
                    seq_list.append(seq)
                    equivalence_flag = True

            if equivalence_flag == False:
                classified_seq.append([seq])

        return classified_seq
    

    def abstract(self):
          # go to abstract definitions
        
        self.abstract_seqs = [seq_list[0].abstract(new_symb=f'S{i}') for i, seq_list in enumerate(self.classified_seqs)]
        
        abstract_contractions = []

        for contr in self.contractions:
            abstract_prods = []
            for prod in [contr.prod0, contr.prod1, contr.prod2]:
                abstract_prod = Product()
                abstract_prod.prefactor = prod.prefactor
                for op in prod.operator:
                    if isinstance(op, SequentialProp):
                        for i, seq_list in enumerate(self.classified_seqs):
                            if op in seq_list:
                                abstract_prod.add(self.abstract_seqs[i])
                    else:
                        abstract_prod.add(op.abstract())

                abstract_prods.append(abstract_prod)
            
            abstract_loops = []
            for loop in contr.loops:
                if len(loop.operator) != 1:
                    raise ValueError()
                op = loop.operator[0]
                abstract_prod = Product()
                abstract_prod.prefactor = loop.prefactor
                for i, seq_list in enumerate(self.classified_seqs):
                    if op in seq_list:
                        abstract_prod.add(self.abstract_seqs[i])

                abstract_loops.append(abstract_prod)

            
            if isinstance(contr, TracefullContribution):
                abstract_contractions.append(TracefullContribution(contr.prefactor, abstract_prods[0], abstract_prods[1], abstract_prods[2], abstract_loops))
            else:
                abstract_contractions.append(TracelessContribution(contr.prefactor, abstract_prods[0], abstract_prods[1], abstract_prods[2], abstract_loops))

        return abstract_contractions
    

    def __str__(self):
        s = "Definition of sequential propagator: \n\n"
        for seq in self.abstract_seqs:
            s += f"{tab}{seq.definition()}\n"

        s += "\n\nCorrelation Function:\n\n"

        for contr in self.abstract_contractions:
            s += f"{tab}{contr.short_str()}\n" 

        return s

    def __sub__(self, other):
        return operation(self, other, lambda x, y: x - y)

    def __add__(self, other):
        return operation(self, other, lambda x, y: x + y)
        
        

    def distillation(self, pos_t, pion_names):
        s = ''

        for seq in self.abstract_seqs:
            s += seq.distillation(pos_t, pion_names)
            s += '\n'

        s += f"{tab}def diagrams(self, t):\n"
        for seq in self.abstract_seqs:
            s += f"{2*tab}{seq.symb} = self.{seq.symb.replace('S', 'seq')}(t)\n"

        s += f"\n\n{2*tab}def contraction(func, ops, factor):\n"
        s += f"{3*tab}d0, d1 = func(self.eps_src, self.eps_snk[t], ops, self.P_snk, self.P_src)\n"
        s += f"{3*tab}d0 *= factor\n"
        s += f"{3*tab}d1 *= factor\n"
        s += f"{3*tab}return d0, d1\n"

        s += '\n\n'
        s += f'{2*tab}# shortcuts\n\n'
        s += f"{2*tab}L = self.perambs[self.t0][t]\n"

        hashed_contr = {}
        ops = []
        for contr in self.abstract_contractions:
            h = contr.hash()
            ops += contr.op_names()
            if hashed_contr.get(h, None) is None:
                hashed_contr[h] = [contr]
            else:
                hashed_contr[h].append(contr)
        
        ops = list(set(ops))
        
        
        GX_ops = []
        XG_ops = []
        GXG_ops = []
        
        for op in ops:

            if op.startswith('G') and op.endswith('G'):
                GXG_ops.append(op)
            elif op.startswith('G'):
                GX_ops.append(op)
            elif op.endswith('G'):
                XG_ops.append(op)
        
        s += f'{2*tab}# GX\n'
        for op in GX_ops:
            s += f'{2*tab}{op} = self.Gamma_snk @ {op[1:]}\n'
        
        s += f"\n{2*tab}# XG \n"
        for op in XG_ops:
            s += f'{2*tab}{op} = {op[:-1]} @ self.Gamma_src\n'
        
        s += f"\n{2*tab}# GXG \n"
        for op in GXG_ops:
            s += f'{2*tab}{op} = self.Gamma_snk @ {op[1:-1]} @ self.Gamma_src\n'
        
        s += f"\n\n{2*tab}diags = []"
        
        combined_contr = []
        d_cross_contr = {}
        combined_contr_l = []
        d_cross_contr_l = {}
        
        
        
        for h, contrs in hashed_contr.items():
            if isinstance(h, str) and len(contrs) == 2:
                combined_contr.append(split_hash(h))
        
            elif isinstance(h, str) and len(contrs) == 1:
                sp = split_hash(h)
                ch = cross_hash(h)
                chs = split_hash(ch)
                if d_cross_contr.get(chs, None) is None:
                    d_cross_contr[sp] = [contrs[0]]
                else:
                    d_cross_contr[chs].append(contrs[0])
        
            elif isinstance(h, tuple) and len(contrs) == 2:
                sh = split_hash(h[0])
                combined_contr_l.append((sh[0], sh[1], h[1]))
        
            elif isinstance(h, tuple) and len(contrs) == 1:
                sp = split_hash(h[0])
                ht = (sp[0], sp[1], h[1])
                ch = cross_hash(sp[1])
                cht = (sp[0], ch, h[1])
                if d_cross_contr_l.get(cht, None) is None:
                    d_cross_contr_l[ht] = [contrs[0]]
                else:
                    d_cross_contr_l[cht].append(contrs[0])
            else:
                raise RuntimeError("THIS SHOULD NOT HAPPEN! contr isn't either combined nor cross! TODO: Include this case in this code!")
        
        cross_contr = list(d_cross_contr.keys())
        cross_contr_l = list(d_cross_contr_l.keys())
                
        
        
        if len(combined_contr) > 0:
            comb_str = '[' + ', '.join(f'({c[0]}, {c[1]})' for c in combined_contr) + ']'
            s += f"\n\n{2*tab}# combined contractions\n"
            s += f"{2*tab}for pre, ops in {comb_str}:\n"
            s += f"{3*tab}di, dj = contraction(combined_contraction, ops, factor=pre)\n"
            s += f"{3*tab}diags.append(di)\n"
            s += f"{3*tab}diags.append(dj)\n"
            s += f"\n\n"
        
        if len(cross_contr) > 0:
            cross_str = '[' + ', '.join(f'({c[0]}, {c[1]})' for c in cross_contr) + ']'
            s += f"{2*tab}# cross contractions\n"
            s += f"{2*tab}for pre, ops in {cross_str}:\n"
            s += f"{3*tab}di, dj = contraction(cross_contraction, ops, factor=pre)\n"
            s += f"{3*tab}diags.append(di)\n"
            s += f"{3*tab}diags.append(dj)\n"
            s += f"\n\n"
            
            s += f"{2*tab}# LOOP CONTRACTIONS\n"
        
        if len(combined_contr_l) > 0:
            comb_l_str = '[' + ', '.join(f'({c[0]}, {c[1]}, {c[2]})' for c in combined_contr_l) + ']'
            
            s += f"{2*tab}# combined contractions\n"
            s += f"{2*tab}for pre, ops, loops in {comb_l_str}:\n"
            s += f"{3*tab}factor = pre * np.prod(loops)\n"
            s += f"{3*tab}di, dj = contraction(combined_contraction, ops, factor=factor)\n"
            s += f"{3*tab}diags.append(di)\n"
            s += f"{3*tab}diags.append(dj)\n"
            s += f"\n\n"
        
        if len(cross_contr_l) > 0:
            cross_l_str = '[' + ', '.join(f'({c[0]}, {c[1]}, {c[2]})' for c in cross_contr_l) + ']'
            s += f"{2*tab}# cross contractions\n"
            s += f"{2*tab}for pre, ops, loops in {cross_l_str}:\n"
            s += f"{3*tab}factor = pre * np.prod(loops)\n"
            s += f"{3*tab}di, dj = contraction(cross_contraction, ops, factor=factor)\n"
            s += f"{3*tab}diags.append(di)\n"
            s += f"{3*tab}diags.append(dj)\n"
            s += f"\n\n"
        
        s += f"{2*tab}return diags"
        return s
    

    def template_2pt(self, class_name, pos_t, pion_names, pions_src=[], pions_snk=[], factor_str=''):
        s = ''
        s += f"class {class_name}(PerambBaryon):\n"
        s += f"{tab}def __init__(self, eps_src, eps_snk, perambs, mins_src" 
        for pi_min in pions_src:
            s += f", mins_{pi_min}"
        s += ", mins_snk"
        for pi_min in pions_snk:
            s += f", mins_{pi_min}"
        s += ", t0, Gamma_src=None, Gamma_snk=None, P_src=None, P_snk=None, tval=None"
        for keyword in pions_src:
            s += f", Gamma_{keyword}=None"
        for keyword in pions_snk:
            s += f", Gamma_{keyword}=None"
        s += "):\n"
        s += f"{2*tab}self.t0 = t0\n"
        s += f"{2*tab}self.eps_src = np.conjugate(eps_src[self.t0])\n"
        s += f"{2*tab}self.eps_snk = eps_snk\n\n"
        s += f"{2*tab}self.perambs = perambs\n"
        s += f"{2*tab}self.Nt = perambs[0].shape[0]\n\n"
        s += f"{2*tab}self.tval = range(self.Nt) if tval is None else tval\n\n"
        
        for pi_min in pions_src:
            s += f"{2*tab}self.mins_{pi_min} = mins_{pi_min}\n"
        for pi_min in pions_snk:
            s += f"{2*tab}self.mins_{pi_min} = mins_{pi_min}\n"

        s += "\n"

        s += f"{2*tab}g_default = lqcdpy.distillation.mat.pion['g5']\n"
        for keyword in pions_src:
            s += f"{2*tab}self.Gamma_{keyword} = g_default if Gamma_{keyword} is None else Gamma_{keyword}\n"
        for keyword in pions_snk:
            s += f"{2*tab}self.Gamma_{keyword} = g_default if Gamma_{keyword} is None else Gamma_{keyword}\n"

        s += f"{2*tab}self.pions = {{\n"
        for pi_min in pions_src:
            s += f"{3*tab}\"{pi_min}\": lambda t: (self.mins_{pi_min}[t], self.Gamma_{pi_min}),\n"
        for pi_min in pions_snk:
            s += f"{3*tab}\"{pi_min}\": lambda t: (self.mins_{pi_min}[t], self.Gamma_{pi_min}),\n"
        s += f"{2*tab}}}\n\n"
        
        s += f"{2*tab}cg_default = lqcdpy.distillation.mat.nucleon['Cg5']\n\n"
        s += f"{2*tab}self.Gamma_src = cg_default if Gamma_src is None else Gamma_src\n"
        s += f"{2*tab}self.Gamma_snk = cg_default if Gamma_snk is None else Gamma_snk\n\n"

        s += f"{2*tab}self.P_src = np.eye(4) if P_src is None else P_src\n"
        s += f"{2*tab}self.P_snk = np.eye(4) if P_snk is None else P_snk\n\n"


        s += self.distillation(pos_t, pion_names)

        s += f"\n\n\n{tab}def corr(self, t):\n"
        s += f"{2*tab}return {factor_str}sum(self.diagrams(t))\n"

        return s
    
    def template_3pt(self, class_name, pos_t, pion_names, pions_src=[], currents=[], pions_snk=[], factor_str=''):
        s = ''
        s += f"class {class_name}(PerambBaryon):\n"
        s += f"{tab}def __init__(self, eps_src, eps_snk, perambs" 
        for pi_min in pions_src:
            s += f", mins_{pi_min}"
        for current in currents:
            s += f", mins_{current}"
        for pi_min in pions_snk:
            s += f", mins_{pi_min}"
        s += ", t0, tau, Gamma_src=None, Gamma_snk=None, P_src=None, P_snk=None, tval=None"
        for keyword in pions_src:
            s += f", Gamma_{keyword}=None"
        for current in currents:
            s += f", Gamma_{current}=None"
        for keyword in pions_snk:
            s += f", Gamma_{keyword}=None"
        s += "):\n"
        s += f"{2*tab}self.t0 = t0\n"
        s += f"{2*tab}self.tau = tau\n"
        s += f"{2*tab}self.eps_src = np.conjugate(eps_src[self.t0])\n"
        s += f"{2*tab}self.eps_snk = eps_snk\n\n"
        s += f"{2*tab}self.perambs = perambs\n"
        s += f"{2*tab}self.Nt = perambs[0].shape[0]\n\n"
        s += f"{2*tab}self.tval = range(self.Nt) if tval is None else tval\n\n"
        
        for pi_min in pions_src:
            s += f"{2*tab}self.mins_{pi_min} = mins_{pi_min}\n"
        for current in currents:
            s += f"{2*tab}self.mins_{current} = mins_{current}\n"
        for pi_min in pions_snk:
            s += f"{2*tab}self.mins_{pi_min} = mins_{pi_min}\n"

        s += "\n"

        s += f"{2*tab}g_default = lqcdpy.distillation.mat.pion['g5']\n"
        for keyword in pions_src:
            s += f"{2*tab}self.Gamma_{keyword} = g_default if Gamma_{keyword} is None else Gamma_{keyword}\n"
        for current in currents:
            s += f"{2*tab}self.Gamma_{current} = g_default if Gamma_{current} is None else Gamma_{current}\n"
        for keyword in pions_snk:
            s += f"{2*tab}self.Gamma_{keyword} = g_default if Gamma_{keyword} is None else Gamma_{keyword}\n"

        s += f"{2*tab}self.pions = {{\n"
        for pi_min in pions_src:
            s += f"{3*tab}\"{pi_min}\": lambda t: (self.mins_{pi_min}[t], self.Gamma_{pi_min}),\n"
        for current in currents:
            s += f"{3*tab}\"{current}\": lambda t: (self.mins_{current}[t], self.Gamma_{current}),\n"
        for pi_min in pions_snk:
            s += f"{3*tab}\"{pi_min}\": lambda t: (self.mins_{pi_min}[t], self.Gamma_{pi_min}),\n"
        s += f"{2*tab}}}\n\n"
        
        s += f"{2*tab}cg_default = lqcdpy.distillation.mat.nucleon['Cg5']\n\n"
        s += f"{2*tab}self.Gamma_src = cg_default if Gamma_src is None else Gamma_src\n"
        s += f"{2*tab}self.Gamma_snk = cg_default if Gamma_snk is None else Gamma_snk\n\n"

        s += f"{2*tab}self.P_src = np.eye(4) if P_src is None else P_src\n"
        s += f"{2*tab}self.P_snk = np.eye(4) if P_snk is None else P_snk\n\n"

        s += self.distillation(pos_t, pion_names)

        s += f"\n\n\n{tab}def corr(self, t):\n"
        s += f"{2*tab}return {factor_str}sum(self.diagrams(t))\n"

        return s

        
