aa = [
    'A',
    'R',
    'D',
    'N',
    'C',
    'E',
    'Q',
    'G',
    'H',
    'I',
    'L',
    'K',
    'M',
    'F',
    'P',
    'S',
    'T',
    'W',
    'Y',
    'V'
    ]
masses = [
    71.0371137870,
    156.1011110260,
    115.0269430310,
    114.0429274460,
    103.00919,
    129.0425930950,
    128.0585775100,
    57.0214637230,
    137.0589118610,
    113.0840639790,
    113.0840639790,
    128.0949630160,
    131.0404846050,
    147.0684139150,
    97.0527638510,
    87.0320284090,
    101.0476784730,
    186.0793129520,
    163.0633285370,
    99.0684139150
    ]

aa_dict = {aa[i]:masses[i] for i in range(20)}

theoretical_mass = lambda xs: sum([aa_dict[x] for x in xs]) + 18.010564686

no_mod_il = lambda pep: ''.join([p.replace('I','L') for p in pep if p.isalpha()])

def theoretical_mz(sequence,charge):
    aa = ''.join([a for a in sequence if a.isalpha()])
    mods = ''.join([m for m in sequence if not m.isalpha()])
    if len(mods) > 0:
        mods = eval(mods)
    else:
        mods = 0
    return (theoretical_mass(aa) + mods + (int(charge)*1.007276035))/int(charge)

def format_mq_pep(mq_pep_in):
    return mq_pep_in.replace('_','').replace('(ox)','+15.995').replace('C','C+57.021').replace('(ac)','+42.011')
