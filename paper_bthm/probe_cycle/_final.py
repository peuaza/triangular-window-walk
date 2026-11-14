import json
J=json.load(open('cyc_cells.json'))
R=J['rules']; E=J['exp']
NS=['512','1024','2048']
def r3(v): return None if v is None else round(float(v),3)
def rd(v): return None if v is None else int(round(float(v)))
lines=[]
for rule in R:
    C=[R[rule].get(n) for n in NS]
    d={"b0_meet":[r3(c['b0_meet_frac']) if c else None for c in C],
       "b0_cap":[r3(c['b0_cap_frac']) if c else None for c in C],
       "cyc_len":[rd(c['cyc_len_med']) if c else None for c in C],
       "cyc_width":[rd(c['cyc_width_med']) if c else None for c in C],
       "bC_meet":[r3(c['bC_meet_frac']) if c else None for c in C],
       "bC_rescue":[r3(c['bC_rescue_frac']) if c else None for c in C],
       "reads0_med":[rd(c['reads0_med_meet']) if c else None for c in C],
       "readsC_med":[rd(c['readsC_med_all']) if c else None for c in C]}
    lines.append('  "%s":%s'%(rule,json.dumps(d,separators=(',',':'))))
exp={}
for rule in R:
    exp[rule]=[E[rule]['readsC_med_all'],E[rule]['reads0_med_meet']]
txt='{\n "n":[512,1024,2048],\n "cap_steps":[10485860,41943140,167772260],\n'
txt+=' "rules":{\n'+',\n'.join(lines)+'\n },\n'
txt+=' "exp_readsC_reads0":'+json.dumps(exp,separators=(',',':'))+'\n}'
print(txt)
open('_final_compact.json','w').write(txt+'\n')
